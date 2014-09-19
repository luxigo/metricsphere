#include <windows.h>
#include <atlstr.h>
#include "MultiStreamIO.h"

////////////////////////////
// VariableIOClientStream //
////////////////////////////
VariableIOClientStream::~VariableIOClientStream( void )
{
	if( _fp ) close( );
	FreeArray( _data );
}
VariableIOClientStream::VariableIOClientStream( void ) : IOClient( )
{
	_fp = NULL;
}
int VariableIOClientStream::close( void )
{
	if( !_fp )
	{
		fprintf( stderr , "Closing null file pointer!\n" );
		return 0;
	}
	EnterCriticalSection( &lock );
	_endIO = true;
	LeaveCriticalSection( &lock );
	while( server ) Sleep(1);
	if( !_read )
	{
		while( _frontAndBack<_current )
		{
			size_t ioBytes = fwrite( _data + ( _frontAndBack % _bufferMultiplier ) * _blockSize , 1 ,_blockSize , _fp );
			WriteBytes += ioBytes;
			_head += ioBytes;
			_frontAndBack++;
		}
		if( _head<_pseudoHead )
		{
			size_t ioBytes = fwrite( _data + ( _current % _bufferMultiplier ) * _blockSize , 1 ,_pseudoHead-_head , _fp );
			WriteBytes += ioBytes;
			_head += ioBytes;
		}
	}
	int ret = fclose( _fp );
	_fp = NULL;
	return ret;
}
void VariableIOClientStream::Initialize( const char* fileName , int maxIOSize , int bufferMultiplier , bool read )
{
	EnterCriticalSection( &lock );

	if( read ) _fp = fopen( fileName , "rb" );
	else       _fp = fopen( fileName , "wb" );
	_bufferMultiplier = bufferMultiplier;

	_head = _pseudoHead = 0;
	_read = read;
	_blockSize = maxIOSize;
	while( _blockSize<IO_BLOCK_SIZE ) _blockSize += maxIOSize;
	_data = AllocArray< byte >( _blockSize*_bufferMultiplier , 1 , "VariableIOClientStream::_data" );
	if( !_data )
	{
		fprintf(stderr , "Failed to allocate memory for BufferedIOState buffer: %d x %d\n" , _blockSize , _bufferMultiplier );
		exit(0);
	}
	_endIO = false;
	_current      = 0;
	_frontAndBack = 0;
	if( _read )
	{
		int ioBytes = fread( _data , 1 , _blockSize , _fp );
		ReadBytes += ioBytes;
		_head      = ioBytes;
		if( ioBytes<_blockSize ) _endIO = true;
		_frontAndBack++;
	}
	LeaveCriticalSection( &lock );
}

// Assumes that whatever data wants to be read has already been brought into the buffer.
size_t VariableIOClientStream::read( void* buffer , size_t sz )
{
	if( sz > _blockSize ) fprintf( stderr , "Requesting a read that is too large: %d\n" , sz );
	EnterCriticalSection( &lock );
	long long start , end;
	int startBlock , endBlock;
	start = _pseudoHead;
	end   = start + sz;
	if( end>_head && _endIO ) end = _head;
	if( (start/_blockSize) != (end/_blockSize) && !server ) Service();		// If we need data from the next block and there is no server
																			// go ahead and read it.

#if ASSERT_MEMORY_ACCESS
	if( end>_head )	fprintf( stderr , "VariableIOClientStream: Index out of bounds: [%lld , %lld]> %lld\n" , start , end , _head ) , exit(0);
#endif // ASSERT_MEMORY_ACCESS
	startBlock =  start  / _blockSize;
	endBlock   = (end-1) / _blockSize;
	int startIdx = startBlock % _bufferMultiplier;
	int   endIdx =   endBlock % _bufferMultiplier;
	int startOffset =    start % _blockSize;
	int   endOffset =( (end-1) % _blockSize ) + 1;
	if( startIdx == endIdx ) memcpy( buffer , _data + startIdx*_blockSize + startOffset , endOffset-startOffset );
	else
	{
		memcpy( buffer , _data + startIdx*_blockSize + startOffset , _blockSize-startOffset );
		memcpy( (void*)( size_t(buffer) + _blockSize-startOffset ) , _data + endIdx*_blockSize , endOffset );
	}
	_pseudoHead += end-start;
	if( (start/_blockSize) != (end/_blockSize) ) _current++;
	LeaveCriticalSection( &lock );
	return end-start;
}
// Assumes that there is room in the buffer for storing the data.
size_t VariableIOClientStream::write( const void* buffer , size_t sz )
{
	if( sz > _blockSize ) fprintf( stderr , "Requesting a write that is too large: %d\n" , sz );
	EnterCriticalSection( &lock );
	long long start , end;
	int startBlock , endBlock;
	start = _pseudoHead;
	end   = start + sz;
	if( (start/_blockSize) != (end/_blockSize) && !server ) Service();		// If we need data from the next block and there is no server
																			// go ahead and read it.

	startBlock =  start  / _blockSize;
	endBlock   = (end-1) / _blockSize;
	int startIdx = startBlock % _bufferMultiplier;
	int   endIdx =   endBlock % _bufferMultiplier;
	int startOffset = start % _blockSize;
	int   endOffset =( (end-1) % _blockSize ) + 1;
	if( startIdx == endIdx )
	{
		memcpy( _data + startIdx*_blockSize + startOffset , buffer , endOffset-startOffset );
	}
	else
	{
		memcpy( _data + startIdx*_blockSize + startOffset , buffer , _blockSize-startOffset );
		memcpy( _data + endIdx*_blockSize , (void*)( size_t(buffer) + _blockSize-startOffset ) , endOffset );
	}
	_pseudoHead += end-start;
	if( (start/_blockSize) != (end/_blockSize) ) _current++;
	LeaveCriticalSection( &lock );
	return end-start;
}

int VariableIOClientStream::Service( void )
{
	DWORD ioBytes;
	// Try and grab the lock:
	//   [W] _head
	//   [W] _frontAndBack
	//   [W] _endIO 
	//   [R] _data
	EnterCriticalSection(  &lock );

	// If we are reading, see if we can advance the front pointer
	if( _read )
		if( !_endIO && _frontAndBack-_current<_bufferMultiplier )
		{
			ioBytes = fread( _data + ( _frontAndBack % _bufferMultiplier ) * _blockSize , 1 ,_blockSize , _fp );
			if( ioBytes != _blockSize ) _endIO = true;
			ReadBytes += ioBytes;
			_head += ioBytes;
			_frontAndBack++;
			LeaveCriticalSection(&lock);
			if( _endIO ) return COMPLETE;
			else		 return SUCCESS;
		}
	// If we are writing, see if we can flush the back pointer
	if( !_read )
		if( _frontAndBack<_current )
		{
			ioBytes = fwrite( _data + ( _frontAndBack % _bufferMultiplier ) * _blockSize , 1 ,_blockSize , _fp );
			WriteBytes += ioBytes;
			_head += ioBytes;
			_frontAndBack++;
			LeaveCriticalSection(&lock);
			return SUCCESS;
		}
		else if( _endIO && _head<_pseudoHead )
		{
			ioBytes = fwrite( _data + ( _current % _bufferMultiplier ) * _blockSize , 1 ,_pseudoHead-_head , _fp );
			WriteBytes += ioBytes;
			_head += ioBytes;
			LeaveCriticalSection(&lock);
			return COMPLETE;
		}
	if( _endIO ) // If we are done
	{
		LeaveCriticalSection( &lock );
		return COMPLETE;
	}
	else
	{
		LeaveCriticalSection( &lock );
		return NONE;
	}
}
/////////////////////////
// FixedIOClientStream //
/////////////////////////

FixedIOClientStream::FixedIOClientStream(void) : IOClient( )
{
	hFile=NULL;
	off = NullPointer< int >( );
	data = NullPointer< byte >( );
	r=rs=win=b=0;
	blockSize=0;
}
FixedIOClientStream::~FixedIOClientStream( void )
{
	FreeArray( off );
	VFreeArray( data );
}
void FixedIOClientStream::Init( HANDLE hFile , int rowSize , int rows , int bufferMultiplier )
{
	FreeArray( off );
	VFreeArray( data );
	this->hFile = hFile;
	r = rows;
	rs = rowSize;
	b = bufferMultiplier;
	off = AllocArray< int >( b , 1 , "FixedIOClientStream::off" );
	if( !off )
	{
		fprintf( stderr , "Failed to allocate memory for offsets\n" );
		exit(0);
	}
}
void FixedIOClientStream::Reset( bool read , int minWindowSize )
{
	this->read=read;
	win=minWindowSize<r ? minWindowSize : r;
	while( win*rs<IO_BLOCK_SIZE && win<r ) win++;

	blockSize=((win*rs+(BYTES_PER_SECTOR-1))/BYTES_PER_SECTOR+1)*BYTES_PER_SECTOR;
	if(!((win*rs)%BYTES_PER_SECTOR))	blockSize-=BYTES_PER_SECTOR;
	VFreeArray( data );
	data = VAllocArray< byte >( blockSize * b , 1 , "FixedIOClientStream::data" );
	if( !data )
	{
		fprintf( stderr , "Failed to allocate memory for StreamState buffer: %d x %d\n" ,  blockSize , b );
		exit(0);
	}
	if( read )
	{
		DWORD ioBytes;
		MySetFilePointer( hFile , 0 );
		MyReadFile( hFile , data , blockSize , &ioBytes , NULL );
		ReadBytes += ioBytes;
	}
	current = 0;
	back    = 0;
	front   = win;
	off[0]  = 0;
}
void FixedIOClientStream::Unset(void)
{
	read=false;
	win=0;
	blockSize=0;
	VFreeArray( data );
	current=0;
	back=0;
	front=win;
	off[0]=0;
}
Pointer( byte ) FixedIOClientStream::operator[]	(int idx)
{
	Pointer( byte ) rowData;
	EnterCriticalSection(&lock);
#if ASSERT_MEMORY_ACCESS
	if(idx<0 || idx<current || idx>=r || idx>=current+win)
		fprintf(stderr,"StreamState: Index out of bounds: %d\t[%d, %d]\t%d x %d\n",idx,current,current+win,rs,r) , exit(0);
#endif // ASSERT_MEMORY_ACCESS
	int bIndex = (idx/win)%b;		// Which block is it in?
	int wIndex =  idx%win;			// Which row of the block?
	rowData = data + bIndex*blockSize + off[bIndex] + wIndex*rs;
	LeaveCriticalSection(&lock);
	return rowData;
}

bool FixedIOClientStream::Advance( void )
{
	EnterCriticalSection(&lock);
	if(current+1>=back && current+1<front)
	{
		current++;
		LeaveCriticalSection(&lock);
		return true;
	}
	else if(current+1>=r)
	{
		current=r;
		LeaveCriticalSection(&lock);
		return true;
	}
	else
	{
		LeaveCriticalSection(&lock);
		return false;
	}
}
int FixedIOClientStream::Service(void)
{
	DWORD ioBytes;
	int ioRows=win;
	int ioBlockSize=blockSize;
	// Try and grab the lock
	EnterCriticalSection(&lock);
	// First see if we can advance the front pointer
	if(front<r && front+win-back <= win*b)
	{
		LONGLONG locationOnDisk=LONGLONG(front)*rs;
		LONGLONG readStart=(locationOnDisk/BYTES_PER_SECTOR)*BYTES_PER_SECTOR;
		int offset=locationOnDisk-readStart;
		int bIndex=(front/win)%b;
		if(read)
		{
			LeaveCriticalSection(&lock);
			if(front+win<=r)	ioRows=win;
			else				ioRows=r-front;
			ioBlockSize=((offset+ioRows*rs+BYTES_PER_SECTOR-1)/BYTES_PER_SECTOR)*BYTES_PER_SECTOR;
			MySetFilePointer(hFile,readStart);
			MyReadFile( hFile , data + bIndex*blockSize , ioBlockSize , &ioBytes , NULL );
			ReadBytes+=ioBytes;
			EnterCriticalSection(&lock);
		}
		off[bIndex]=offset;
		front+=ioRows;
		LeaveCriticalSection(&lock);
		return SUCCESS;
	}
	// Now try to free up trailing memory
	else if( (back+win<=current || current>=r) && back<r)	// If we won't write out needed data and there is what to write
	{
		LONGLONG locationOnDisk=LONGLONG(back)*rs;
		LONGLONG writeStart=(locationOnDisk/BYTES_PER_SECTOR)*BYTES_PER_SECTOR;
		int offset=locationOnDisk-writeStart;
		int bIndex=(back/win)%b;
		if(!read)											// If we are doing a write, write out the data
		{
			LeaveCriticalSection(&lock);
			if(back+win<=r)	ioRows=win;
			else			ioRows=r-back;
			ioBlockSize=((offset+ioRows*rs+BYTES_PER_SECTOR-1)/BYTES_PER_SECTOR)*BYTES_PER_SECTOR;
			MySetFilePointer(hFile,writeStart);
			MyWriteFile( hFile , data + bIndex*blockSize , ioBlockSize , &ioBytes , NULL );
			WriteBytes+=ioBytes;
			EnterCriticalSection(&lock);
		}
		back+=ioRows;
		if(!read && back<r)
		{
			LONGLONG locationOnDisk=LONGLONG(back)*rs;
			LONGLONG writeStart=(locationOnDisk/BYTES_PER_SECTOR)*BYTES_PER_SECTOR;
			int offset=locationOnDisk-writeStart;		// The number of bytes that need to be copied over from the previous buffer
			int bIndex=(back/win)%b;
			int oldBIndex=(bIndex+b-1)%b;
			if(offset)	memcpy( data + bIndex*blockSize ,  data + oldBIndex*blockSize + ioBlockSize-BYTES_PER_SECTOR , offset );
		}
		LeaveCriticalSection(&lock);
		return SUCCESS;
	}
	// Check if we are done
	else if(back>=r)										// If we have already written out the last row
	{
		LeaveCriticalSection(&lock);
		return COMPLETE;
	}
	else
	{
		LeaveCriticalSection(&lock);
		return NONE;
	}
}
/////////////////////////
// MultiStreamIOClient //
/////////////////////////
MultiStreamIOClient::MultiStreamIOClient( const char* fileName , int rs , int r , int bufferMultiplier , bool writeOnly )
{
	CString s( fileName );
	hFile = CreateFile( s , FILE_READ_DATA | FILE_WRITE_DATA,0,NULL,OPEN_ALWAYS , FILE_FLAG_NO_BUFFERING , NULL );
	if( !hFile )
	{
		fprintf(stderr,"Failed to create file handle\n");
		PrintError();
		exit(0);
	}
	if( writeOnly )
	{
		// Pre-allocate file space
		long long fileSize;
		fileSize=(LONGLONG)(r)*rs;
		fileSize=((fileSize+IOClient::BYTES_PER_SECTOR-1)/IOClient::BYTES_PER_SECTOR)*IOClient::BYTES_PER_SECTOR;
		if(MySetFilePointer(hFile,fileSize)==INVALID_SET_FILE_POINTER)
		{
			PrintError();
			exit(0);
		}
		SetEndOfFile(hFile);
	}
	stream.Init( hFile , rs , r , bufferMultiplier );
	server = NULL;
}
MultiStreamIOClient::MultiStreamIOClient( int rs , int r , int bufferMultiplier , const char* prefix , bool deleteOnClose )
{
	// Create a temporary file
	char* tmp;
	char fullPrefix[512];
	if(prefix)	sprintf(fullPrefix,"%s_%lld_",prefix,unsigned long long(GetCurrentProcessId()));
	else		sprintf(fullPrefix,"scratch_%lld_",unsigned long long(GetCurrentProcessId()));

	tmp = _tempnam( "." , fullPrefix );

	CString s(tmp);
	if(deleteOnClose)	hFile = CreateFile(s,FILE_READ_DATA | FILE_WRITE_DATA,0,NULL,CREATE_ALWAYS, FILE_FLAG_NO_BUFFERING | FILE_FLAG_DELETE_ON_CLOSE, NULL);
	else				hFile = CreateFile(s,FILE_READ_DATA | FILE_WRITE_DATA,0,NULL,CREATE_ALWAYS, FILE_FLAG_NO_BUFFERING, NULL);
	if(!hFile)
	{
		fprintf(stderr,"Failed to create file handle\n");
		PrintError();
		exit(0);
	}

	// Pre-allocate file space
	long long fileSize;
	fileSize=(LONGLONG)(r)*rs;
	fileSize=((fileSize+IOClient::BYTES_PER_SECTOR-1)/IOClient::BYTES_PER_SECTOR)*IOClient::BYTES_PER_SECTOR;
	if(MySetFilePointer(hFile,fileSize)==INVALID_SET_FILE_POINTER)
	{
		PrintError();
		exit(0);
	}
	SetEndOfFile(hFile);

#if NEW_MULTI_STREAM_CODE
	stream.Init( hFile , rs , r , bufferMultiplier );
#else // !NEW_MULTI_STREAM_CODE
	stream.Init(hFile,rs,r,bufferMultiplier);
#endif // NEW_MULTI_STREAM_CODE
	server=NULL;
}
MultiStreamIOClient::MultiStreamIOClient( int rs , int r , int bufferMultiplier , const char* dir , const char* prefix , bool deleteOnClose )
{
	// Create a temporary file
	char fullPrefix[512];
	if(prefix)	sprintf( fullPrefix , "%s_%lld",prefix , unsigned long long( GetCurrentProcessId() ) );
	else		sprintf( fullPrefix , "scratch_%lld"   , unsigned long long( GetCurrentProcessId() ) );
	CString s(_tempnam(dir,fullPrefix));

	if(deleteOnClose)	hFile = CreateFile(s,FILE_READ_DATA | FILE_WRITE_DATA,0,NULL,CREATE_ALWAYS, FILE_FLAG_NO_BUFFERING | FILE_FLAG_DELETE_ON_CLOSE, NULL);
	else				hFile = CreateFile(s,FILE_READ_DATA | FILE_WRITE_DATA,0,NULL,CREATE_ALWAYS, FILE_FLAG_NO_BUFFERING, NULL);

	if(!hFile)
	{
		fprintf(stderr,"Failed to create file handle\n");
		PrintError();
		exit(0);
	}

	// Pre-allocate file space
	long long fileSize;
	fileSize=(LONGLONG)(r)*rs;
	fileSize=((fileSize+IOClient::BYTES_PER_SECTOR-1)/IOClient::BYTES_PER_SECTOR)*IOClient::BYTES_PER_SECTOR;
	if(MySetFilePointer(hFile,fileSize)==INVALID_SET_FILE_POINTER)
	{
		PrintError();
		exit(0);
	}
	SetEndOfFile(hFile);

#if NEW_MULTI_STREAM_CODE
	stream.Init( hFile , rs , r , bufferMultiplier );
#else // !NEW_MULTI_STREAM_CODE
	stream.Init(hFile,rs,r,bufferMultiplier);
#endif // NEW_MULTI_STREAM_CODE
	server=NULL;
}
MultiStreamIOClient::~MultiStreamIOClient(void)
{
	finalize();
}
void MultiStreamIOClient::finalize(void)
{
	if(hFile)	CloseHandle(hFile);
	hFile=NULL;
}
int MultiStreamIOClient::Service		( void )								{ return stream.Service(); }
int		MultiStreamIOClient::rows		( void )						const	{ return stream.r;  }
int		MultiStreamIOClient::rowSize	( void )						const	{ return stream.rs; }
void	MultiStreamIOClient::reset		( bool r , int minWindowSize )			{ stream.Reset( r , minWindowSize ); }
void	MultiStreamIOClient::unset		( void )								{ stream.Unset(); }
Pointer( byte ) MultiStreamIOClient::operator []( int idx ) { return stream[idx]; }
bool	MultiStreamIOClient::Advance	( void )								{ return stream.Advance(); }
void	MultiStreamIOClient::advance	( void )
{
	// BADNESS!!! Server may not be NULL if it was set for an earlier server and not for the newer one.
	if( server )
	{
		while(1)
			if( Advance() )	return;
			else server->SetPending( this );
	}
	else
		if( !Advance() )
		{
			int serviceState = Service();
			while( serviceState==IOClient::SUCCESS ) serviceState = Service();
			if( !Advance() )
			{
//				fprintf( stderr , "Shouldn't happen: %d %d %d\n" , current , front , r );
				fprintf( stderr , "Shouldn't happen: %d %d %d\n" , stream.current , stream.front , stream.r );
				exit(0);
			}
		}
		else Service();	// To make sure that the last row gets written
}
void MultiStreamIOClient::SetServer( MultiStreamIOServer* server ) { IOClient::SetServer( server ); }