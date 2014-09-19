#define SOCKET_DEBUG 0
#define MAJOR_SOCKET_DEBUG 0
#include "Time.h"
#include "GridStream.h"

DWORD MySetFilePointer(HANDLE hFile,LONGLONG distanceToMove)
{
	LARGE_INTEGER li;
	li.QuadPart=distanceToMove;
	DWORD ret=SetFilePointer(hFile,li.LowPart,&li.HighPart,FILE_BEGIN);
	if(ret==INVALID_SET_FILE_POINTER)
	{
		PrintError();
		exit(0);
	}
	return ret;
}
void MyReadFile( HANDLE hFile , Pointer( byte ) lpBuffer , DWORD nNumberOfBytesToRead , LPDWORD lpNumberOfBytesRead , LPOVERLAPPED lpOverlapped )
{
	if( !ReadFile( hFile , lpBuffer , nNumberOfBytesToRead , lpNumberOfBytesRead , lpOverlapped ) )
	{
		fprintf(stderr,"Failed to read file: %d / %d\n",*lpNumberOfBytesRead,nNumberOfBytesToRead);
		PrintError();
		exit(0);
	}
}
void MyWriteFile( HANDLE hFile , Pointer( byte ) lpBuffer , DWORD nNumberOfBytesToWrite , LPDWORD lpNumberOfBytesWritten , LPOVERLAPPED lpOverlapped )
{
	if( !WriteFile( hFile , lpBuffer , nNumberOfBytesToWrite , lpNumberOfBytesWritten , lpOverlapped ) )
	{
		fprintf(stderr,"Failed to write file\n");
		PrintError();
		exit(0);
	}
}
///////////////////////
// NULLStreamingGrid //
///////////////////////
NULLStreamingGrid::NULLStreamingGrid( int rowSize , int rows )
{
	_rows = rows;
	_rowSize = rowSize;
	_data = AllocArray< byte >( rowSize , 1 , "NULLStreamingGrid::_data");
	memset( _data , 0 , rowSize );
}
NULLStreamingGrid::~NULLStreamingGrid( void )
{
	FreeArray( _data );
}

int NULLStreamingGrid::rows( void ) const { return _rows; }
int NULLStreamingGrid::rowSize( void ) const { return _rowSize; }
Pointer( byte ) NULLStreamingGrid::operator [] ( int idx ){ return _data; }

///////////////////////////
// BufferedStreamingGrid //
///////////////////////////
BufferedStreamingGrid::BufferedStreamingGrid(StreamingGrid* sg)
{
	data = NullPointer< byte >( );
	current=0;
	this->sg=sg;
}
BufferedStreamingGrid::~BufferedStreamingGrid(void)
{
	FreeArray( data );
	sg=NULL;
}
int BufferedStreamingGrid::rows		(void) const	{return sg->rows();}
int BufferedStreamingGrid::rowSize	(void) const	{return sg->rowSize();}
//void BufferedStreamingGrid::reset	(bool read,bool write,int minWindowSize)
void BufferedStreamingGrid::reset	(bool read,int minWindowSize)
{
	this->read=read;
	win=minWindowSize<rows() ? minWindowSize : rows();

	FreeArray( data );
	data = AllocArray< byte >( rowSize() * win , 1 , "BufferedStreamingGrid::_data" );
	if(!data)
	{
		fprintf( stderr , "Failed to allocate memory for BufferedStreamingGrid\n" );
		exit( 0 );
	}
	sg->reset(read,1);
	if(read)
	{
		for(int w=0;w<win;w++)
		{
			Pointer( byte ) row = (*sg)[w];
			memcpy( data+w*rowSize() , row , rowSize() );
			sg->advance();
		}
	}
	current=0;
}
void BufferedStreamingGrid::advance(void)
{
	if(read)
	{
		current++;
		if(current+win-1<rows())
		{
			sg->advance();
			memcpy( data + ((current+win-1)%win)*rowSize() , (*sg)[current+win-1] , rowSize() );
		}
	}
	else
	{
		if(current-win+1>=0)
		{
			memcpy( (*sg)[current-win+1] , data+((current-win+1)%win)*rowSize() , rowSize() );
			sg->advance();
		}
		current++;
	}
}
Pointer( byte ) BufferedStreamingGrid::operator[]	(int idx)
{
#if ASSERT_MEMORY_ACCESS
	if(idx<0 || idx<current || idx>=rows() || idx>=current+win)
		fprintf(stderr,"BufferedStreamingGrid: Index out of bounds: %d\t[%d, %d]\t%d x %d\n",idx,current,current+win,rowSize(),rows()) , exit(0);
#endif // ASSERT_MEMORY_ACCESS
	return data+(idx%win)*rowSize();
}
//////////////////////
// MemoryBackedGrid //
//////////////////////
MemoryBackedGrid::MemoryBackedGrid( Pointer( byte ) data , int rs , int r , bool separateColors )
{
	_separateColors = separateColors;
	_del=false;
	_r=r, _rs=rs;
	_data=data;
}
MemoryBackedGrid::MemoryBackedGrid(int rs,int r , bool separateColors )
{
	_separateColors = separateColors;
	_del=true;
	_r=r, _rs=rs;
	_data = AllocArray< byte >( r * rs , 1 , "MemoryBackedGrid::data" );
	if( !_data )	fprintf(stderr,"Failed to allocate memory backed grid of size %d x %d\n" , rs , r ) , exit(0);
}
MemoryBackedGrid::~MemoryBackedGrid(void)
{
	if( _del ) FreeArray( _data );
	_r=_rs=0;
}
int	MemoryBackedGrid::rows(void) const
{
	return _r;
}
int MemoryBackedGrid::rowSize(void) const
{
	return _rs;
}
Pointer( byte ) MemoryBackedGrid::operator[](int idx)
{
#if ASSERT_MEMORY_ACCESS
	if(idx<0 || idx>=_r)	fprintf(stderr,"Index out of bounds: %d [%d, %d]\n",idx,0,_r) , exit(0);
#endif // ASSERT_MEMORY_ACCESS
	return _data + _rs * idx;
}
#if NEW_CONNECTION_BACKED_GRID
//////////////////////////
// DataStreamBackedGrid //
//////////////////////////
DataStreamBackedGrid::DataStreamBackedGrid( DataStream* stream , int rs , int r , bool separateColors )
{
	_separateColors = separateColors;
	_read = true;
	_r = r;
	_rs = rs;
	_stream = stream;
	_data = malloc( rs );
	if( !_data ) fprintfId( stderr , "Failed to allocate\n" ) , exit(0);
}
DataStreamBackedGrid::~DataStreamBackedGrid(void)
{
	free( _data );
	_data = NULL;
	_r = _rs = 0;
	_stream = NULL;
}
int	DataStreamBackedGrid::rows(void) const { return _r; }
int DataStreamBackedGrid::rowSize(void) const { return _rs; }
void* DataStreamBackedGrid::operator[]( int idx )
{
	if( _read && !_readComplete ) _stream->read( _data , _rs ) , _readComplete = true;
#if ASSERT_MEMORY_ACCESS
	if( idx<0 || idx>=_r || idx!=_c) fprintfId( stderr , "Index out of bounds: %d != %d\n" , idx , _c )  , exit(0);
#endif // ASSERT_MEMORY_ACCESS
	return _data;
}
void DataStreamBackedGrid::advance( void )
{
	if( !_read ) _stream->write( _data , _rs );
	else
	{
		if( !_readComplete ) _stream->read( _data , _rs );
		_readComplete = false;
	}
	_c++;
}
void DataStreamBackedGrid::reset( bool read , int minWindowSize )
{
	_read = read;
	_readComplete = !_read;
	_c = 0;
}

#else // !NEW_CONNECTION_BACKED_GRID
//////////////////////
// SocketBackedGrid //
//////////////////////
SocketBackedGrid::SocketBackedGrid( SOCKET sock , int rs , int r , bool bs , bool separateColors )
{
	_separateColors = separateColors;
	_read = true;
	_r = r;
	_rs = rs;
	_sock = sock;
	_bs = bs;
	_data = AllocArray< byte >( rs , 1 , "SocketBackedGrid::data" );
	if( !_data ) fprintfId( stderr , "Failed to allocate\n" ) , exit(0);
}
SocketBackedGrid::~SocketBackedGrid(void)
{
	FreeArray( _data );
	_data = NullPointer< byte >( );
	_r = _rs = 0;
}
int	SocketBackedGrid::rows(void) const { return _r; }
int SocketBackedGrid::rowSize(void) const { return _rs; }
Pointer( byte ) SocketBackedGrid::operator[]( int idx )
{
#if SOCKET_DEBUG
// 	for( int i=0 ; i<_sockCount ; i++ )	printfId( "\tSocket[%d]: %d , port[%d]: %d , bytes[%d]: %d\n" , i , _socks[i] , i , GetLocalSocketPort( _socks[i] ) , i , _rowsizes[i]*Channels*sizeof(Real) ) , fflush( stdout );
	{
		MyWinSock::StdinLock lock;
		if( _read )	printfId( "Reading row[%d] , Socket: %d , bytes: %d\n" , _c , _sock , _rs ) , fflush( stdout );
		else		printfId( "Writing row[%d] , Socket: %d , bytes: %d\n" , _c , _sock , _rs ) , fflush( stdout );
	}
#endif // SOCKET_DEBUG
	if( _read && !_readComplete )
	{
		ReceiveOnSocket( _sock , _data , _rs , _bs , "SocketBackedGrid::operator[]" ) , _readComplete = true;
#if MAJOR_SOCKET_DEBUG
			{
				MyWinSock::StdinLock lock;
				if( idx<MAJOR_SOCKET_DEBUG )
				{
					printfId("Receiving[%d]\t",idx);
					for(int k=0;k<MAJOR_SOCKET_DEBUG;k++) printf("%d ",((int*)_data)[k]);
					printf("\n") , fflush( stdout );
				}
			}
#endif // MAJOR_SOCKET_DEBUG
	}
#if ASSERT_MEMORY_ACCESS
	if( idx<0 || idx>=_r || idx!=_c) fprintfId( stderr , "[SocketBackedGrid] Index out of bounds: %d != %d || [ %d, %d )\n" , idx , _c , 0 , _r ) , exit(0);
#endif // ASSERT_MEMORY_ACCESS
#if SOCKET_DEBUG
	{
		MyWinSock::StdinLock lock;
		if( _read )	printfId( "Done reading socket (%d)\n" , _sock ) , fflush( stdout );
		else		printfId( "Done writing socket (%d)\n" , _sock ) , fflush( stdout );
	}
#endif // SOCKET_DEBUG
	return _data;
}
void SocketBackedGrid::advance( void )
{
#if SOCKET_DEBUG
	{
		MyWinSock::StdinLock lock;
		printfId( "Advancing row[%d]\t(%d)\n" , _c , _sock ) , fflush( stdout );
	}
#endif // SOCKET_DEBUG
	if( !_read )
	{
#if SOCKET_DEBUG
//		printfId( " Sending %d bytes on socket\n", _rs ) , fflush( stdout );
#endif // SOCKET_DEBUG
		SendOnSocket( _sock , _data , _rs , _bs , "SocketBackedGrid::advance" );
#if MAJOR_SOCKET_DEBUG
		{
			MyWinSock::StdinLock lock;
			if( _c<MAJOR_SOCKET_DEBUG )
			{
				printfId("Sending[%d]\t",_c);
				for(int k=0;k<MAJOR_SOCKET_DEBUG;k++) printf("%d ",((int*)_data)[k]);
				printf("\n") , fflush( stdout );
			}
		}
#endif // MAJOR_SOCKET_DEBUG
	}
	else
	{
		if( !_readComplete ) ReceiveOnSocket( _sock , _data , _rs , _bs , "SocketBackedGrid::advance" );
#if MAJOR_SOCKET_DEBUG
		{
			MyWinSock::StdinLock lock;
			if( _c<MAJOR_SOCKET_DEBUG )
			{
				printfId("Receiving[%d]\t",_c);
				for(int k=0;k<MAJOR_SOCKET_DEBUG;k++) printf("%d ",((int*)_data)[k]);
				printf("\n") , fflush( stdout );
			}
		}
#endif // MAJOR_SOCKET_DEBUG
		_readComplete = false;
	}
	_c++;
#if SOCKET_DEBUG
	{
		MyWinSock::StdinLock lock;
		printfId("Done advancing socket\t(%d)\n",_sock) , fflush( stdout );
	}
#endif // SOCKET_DEBUG
}
void SocketBackedGrid::reset( bool read , int minWindowSize )
{
	if(!_read ) EndSendOnSocket( _sock , _bs , "SocketBackedGrid::reset" );
#if SOCKET_DEBUG
	{
		MyWinSock::StdinLock lock;
		if( read )	printfId( "Resetting socket for read \t(%d)\n" , _sock ) , fflush( stdout );
		else		printfId( "Resetting socket for write\t(%d)\n" , _sock ) , fflush( stdout );
	}
#endif // SOCKET_DEBUG
	_read = read;
	_readComplete = !_read;
	_c = 0;
	if( _read ) StartReceiveOnSocket( _sock , _bs , "SocketBackedGrid::reset" );
}
#endif // NEW_CONNECTION_BACKED_GRID

//////////////////////
// SharedMemoryGrid // 
//////////////////////
SharedMemoryGrid::SharedMemoryGrid( HANDLE* startHandles , HANDLE* endHandles , int handleCount , Pointer( byte ) row , int rowSize , int rows )
{
	_read = true;
	_r = rows;
	_rs = rowSize;
	_data = row;
	_hCount = handleCount;
	_startHandles = startHandles;
	_endHandles   =   endHandles;
}
SharedMemoryGrid::~SharedMemoryGrid(void)
{
	_data = NullPointer< byte >( );
	_r = _rs = 0;
}
int	SharedMemoryGrid::rows(void) const { return _r; }
int SharedMemoryGrid::rowSize(void) const { return _rs; }
Pointer( byte ) SharedMemoryGrid::operator[]( int idx )
{
	if( !_dataReady )
//		if( WaitForMultipleObjects( _hCount , _startHandles , true , INFINITE ) != WAIT_OBJECT_0 )
//			fprintfId( stderr , "Wait for multiple failed: " ) , PrintError();
		MyWaitForMultipleObjects( _hCount , _startHandles , 10000 , "SharedMemoryGrid::operator[%d / %d]" , idx , _r );
	_dataReady = true;
#if ASSERT_MEMORY_ACCESS
	if( idx<0 || idx>=_r || idx!=_c ) fprintfId( stderr , "Index out of bounds: %d != %d || [ %d , %d )\n" , idx , _c , 0 , _r ) , exit(0);
#endif // ASSERT_MEMORY_ACCESS
	return _data;
}
void SharedMemoryGrid::advance( void )
{
	for( int i=0 ; i<_hCount ; i++ ) if( !SetEvent( _endHandles[i] ) ) fprintfId( stderr , "SetEvent failed: " ) , PrintError();

	_dataReady = false;
	_c++;
}
void SharedMemoryGrid::reset( bool read , int minWindowSize )
{
	_read = read;
	_dataReady = false;
	_c = 0;
}


////////////////////
// FileBackedGrid //
////////////////////
FileBackedGrid::FileBackedGrid(int rs,int r)
{
	_rSize=_wSize=0;
	_name=_tempnam(".","");
	_fp=fopen(_name,"w+b");
	_read=true;
	_r=r, _rs=rs;
	_data = NullPointer< byte >( );
	_c=_w=0;
}
FileBackedGrid::FileBackedGrid(int rs,int r,FILE* fp)
{
	_rSize=_wSize=0;
	_name=NULL;
	_fp=fp;
	_read=true;
	_r=r, _rs=rs;
	_data = NullPointer< byte >( );
	_c=_w=0;
}
FileBackedGrid::~FileBackedGrid(void)
{
	if(!_read)	for(int i=0;i<_w;i++)	_writeNext();
	FreeArray( _data );
	_r=_rs=_c=_w=0;
	if(_name)
	{
		fclose(_fp);
		remove(_name);
	}
	printf("Read/Write Size: %lld MB/ %lld MB\n",_rSize>>20,_wSize>>20) , fflush( stdout );
}
int	FileBackedGrid::rows(void) const	{	return _r;	}
int FileBackedGrid::rowSize(void) const	{	return _rs;	}
Pointer( byte ) FileBackedGrid::operator[]( int idx )
{
#if ASSERT_MEMORY_ACCESS
	if(idx<0 || idx>=_r || idx<_c || idx>=_c+_w)	fprintf(stderr,"Index out of bounds: %d [0 ,[%d, %d], %d]\n",idx,_c,_c+_w,_r) , exit(0);
#endif // ASSERT_MEMORY_ACCESS
	return _data + _rs * (idx%_w);
}
void FileBackedGrid::reset(bool read,int minWindowSize)
{
	if(!_read)	for(int i=0;i<_w;i++)	_writeNext();
	FreeArray( _data );

	fseek(_fp,0,SEEK_SET);
	_read=read;
	_w=minWindowSize;
	if(_w>_r)	_w=_r;
	_data = AllocArray< byte >( _rs * _w , 1 , "FileBackedGrid::data" );
	if( !_data ) fprintf(stderr,"Failed to allocate in FileBackedGrid::reset\n") , exit(0);
	if(_read)
	{
		_c = -_w;
		for(int i=0;i<_w;i++)	_readNext();
	}
	else	_c = 0;
}
void FileBackedGrid::_readNext(void)
{
	if( _c+_w<_r ) fread( _data + _rs*((_c+_w)%_w) , 1 , _rs , _fp ) , _rSize+=_rs;
	_c++;
}
void FileBackedGrid::_writeNext(void)
{
	_c++;
	if(_c-_w>=0)	fwrite( _data + _rs*((_c+_w)%_w) , 1 , _rs , _fp ) ,	_wSize+=_rs;
}
void FileBackedGrid::advance(void)
{
	if(_read)	_readNext();
	else		_writeNext();
}

//////////////////////////////
// CompressedFileBackedGrid //
//////////////////////////////
CompressedFileBackedGrid::CompressedFileBackedGrid(int rs,int r)
{
	_rSize=_wSize=0;
	_name=_tempnam(".","");
	if(!(_fp=fopen(_name,"wb+")))	fprintf(stderr,"Failed to open: %s\n",_name) , exit(0);
	_read=true;
	_r=r, _rs=rs;
	_c=_w=0;
	_sSize=(_rs+12)*1.002;
	_data = NullPointer< byte >( );
	_scratch = AllocArray< byte >( _sSize , 1 , "CompressedFileBackedGrid::_scratch" );
	if(!_scratch)	fprintf(stderr,"Failed to allocate in CompressedFileBackedGrid::CompressedFileBackedGrid\n") , exit(0);
}
CompressedFileBackedGrid::CompressedFileBackedGrid(int rs,int r,FILE* fp)
{
	_rSize=_wSize=0;
	_name=NULL;
	_fp=fp;
	_read=true;
	_r=r, _rs=rs;
	_c=_w=0;
	_sSize=(_rs+12)*1.002;
	_data = NullPointer< byte >( );
	_scratch = AllocArray< byte >( _sSize , 1 , "CompressedFileBackedGrid::_scratch" );
	if(!_scratch)	fprintf(stderr,"Failed to allocate in CompressedFileBackedGrid::CompressedFileBackedGrid\n") , exit(0);
}
CompressedFileBackedGrid::~CompressedFileBackedGrid(void)
{
	if(!_read)	for(int i=0;i<_w;i++)	_writeNext();
	FreeArray( _data );
	_r=_rs=_c=_w=0;
	if(_name)
	{
		fclose(_fp);
		remove(_name);
		free(_name);
	}
	printf("Read/Write Size: %lld MB / %lld MB\n",_rSize>>20,_wSize>>20) , fflush( stdout );
}
int	CompressedFileBackedGrid::rows(void) const	{	return _r;	}
int CompressedFileBackedGrid::rowSize(void) const	{	return _rs;	}
Pointer( byte ) CompressedFileBackedGrid::operator[](int idx)
{
#if ASSERT_MEMORY_ACCESS
	if(idx<0 || idx>=_r)	fprintf(stderr,"Index out of bounds: %d [0 , %d]\n",idx,_r)								 , exit(0);
	if( _read && (idx<_c || idx>=_c+_w) )	fprintf(stderr,"Read index out of bounds: %d [%d , %d]\n",idx,_c,_c+_w)	 , exit(0);
	if(!_read && (idx>_c || idx<=_c-_w) )	fprintf(stderr,"Write index out of bounds: %d [%d , %d]\n",idx,_c-_w,_c) , exit(0);
#endif // ASSERT_MEMORY_ACCESS
	return _data + _rs*(idx%_w);
}
void CompressedFileBackedGrid::reset(bool read,int minWindowSize)
{
	if(!_read)	for(int i=0;i<_w;i++)	_writeNext();

	FreeArray( _data );

	if(fseek(_fp,0,SEEK_SET))	fprintf(stderr,"fseek failed in CompressedFileBackedGrid::reset\n") , exit(0);
	_read=read;
	_w=minWindowSize;
	if(_w>_r)	_w=_r;
	_data = AllocArray< byte >( _rs * _w , 1 , "CompressedFileBackedGrid::_data" );
	if(!_data)	fprintf(stderr,"Failed to allocate in CompressedFileBackedGrid::reset\n") , exit(0);
	if(_read)	for( _c = -_w ; _c < 0 ; _readNext() )	;
	else		_c = 0;
}
#include <Util/ZLIB/ZLIB.h>
void CompressedFileBackedGrid::_readNext(void)
{
	if(_c+_w<_r)
	{
		uLong destLen=_rs,sourceLen;
		Bytef *dest = (Bytef*)(LONGLONG(_data)+_rs*((_c+_w)%_w));
		const Bytef *source = (Bytef*)PointerAddress( _scratch );
		if(fread(&sourceLen,sizeof(uLong),1,_fp)!=1)	fprintf(stderr,"Read failed\n") , exit(0);
		if(fread(_scratch,1,sourceLen,_fp)!=sourceLen)	fprintf(stderr,"Read failed\n") , exit(0);
		_rSize+=sizeof(uLong)+sourceLen;

		switch(uncompress(dest,&destLen,source,sourceLen))
		{
		case Z_OK:			break;
		case Z_MEM_ERROR:	fprintf(stderr,"uncompress -- Z_MEM_ERROR\n")  , exit(0);
		case Z_BUF_ERROR:	fprintf(stderr,"uncompress -- Z_BUF_ERROR\n")  , exit(0);
		case Z_DATA_ERROR:	fprintf(stderr,"uncompress -- Z_DATA_ERROR\n") , exit(0);
		}
	}
	_c++;
}
void CompressedFileBackedGrid::_writeNext(void)
{
	_c++;
	if(_c-_w>=0 && _c-_w<_r)
	{
		uLong destLen=_sSize,sourceLen=_rs;
		Bytef *dest = (Bytef*)PointerAddress( _scratch );
		Bytef *source = (Bytef*)(LONGLONG(_data)+_rs*((_c+_w)%_w));
		switch(compress(dest,&destLen,source,sourceLen))
		{
		case Z_OK:	break;
		case Z_MEM_ERROR:	fprintf(stderr,"compress -- Z_MEM_ERROR\n")   , exit(0);
		case Z_BUF_ERROR:	fprintf(stderr,"compress -- Z_BUF_ERROR\n")   , exit(0);
		default:			fprintf(stderr,"compress -- unknown error\n") , exit(0);
		}
		if(fwrite(&destLen,sizeof(uLong),1,_fp)!=1) fprintf(stderr,"Write failed\n") , exit(0);
		if(fwrite(_scratch,1,destLen,_fp)!=destLen) fprintf(stderr,"Write failed\n") , exit(0);
		_wSize+=sizeof(uLong)+destLen;
	}
}
void CompressedFileBackedGrid::advance(void)
{
	if(_read)	_readNext();
	else		_writeNext();
}
