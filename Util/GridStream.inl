#define SOCKET_DEBUG 0
#define MAJOR_SOCKET_DEBUG 0
#include "Time.h"

#if NEW_CONNECTION_BACKED_GRID
///////////////////////////////
// MultiDataStreamBackedGrid //
///////////////////////////////
template< int Channels , class Real >
MultiDataStreamBackedGrid< Channels , Real >::MultiDataStreamBackedGrid( DataStream** streams , int* rowSizes , int streamCount , int rows , bool separateColors )
{
	_read = true;
	_separateColors = separateColors;
	_streamCount = streamCount;
	_rows = rows;
	_rowSize = 0;
	for( int i=0 ; i<_streamCount ; i++ ) _rowSize += rowSizes[i];

	_rowSizes = ( int* )malloc( sizeof(int) * _streamCount );
	_subData = (Real**)malloc( sizeof(Real*) *_streamCount );
	_streams = (DataStream**) malloc( sizeof(DataStream*) * _streamCount );
	if( !_rowSizes || !_streams || !_subData ) fprintf( stderr , "Failed to allocate MultiDataStreamBackedGrid 1\n" ) , exit(0);
	for( int i=0 ; i<_streamCount ; i++ )
	{
		_rowSizes[i] = rowSizes[i];
		_streams[i] = streams[i];
		_subData[i] = (Real*)malloc( sizeof(Real) * Channels * rowSizes[i] );
		if( !_subData[i] ) fprintf( stderr , "Failed to allocate MultiDataStreamBackedGrid 2\n" ) , exit(0);
	}

	_data = ( Real* ) malloc( _rowSize*Channels*sizeof(Real) );
	if( !_data ) fprintf( stderr , "Failed to allocate\n" ) , exit(0);
}
template< int Channels , class Real >
MultiDataStreamBackedGrid< Channels , Real >::~MultiDataStreamBackedGrid(void)
{
	for( int i=0 ; i<_streamCount ; i++ ) free( _subData[i] );
	free( _data );
	free( _rowSizes );
	free( _subData );
	free( _streams );

	_rows = _rowSize = 0;
	_data = NULL;
	_rowSizes = NULL;
	_streams = NULL;
	_subData = NULL;
}
template< int Channels , class Real >
int	MultiDataStreamBackedGrid< Channels , Real >::rows( void ) const { return _rows; }
template< int Channels , class Real >
int MultiDataStreamBackedGrid< Channels , Real >::rowSize( void ) const { return _rowSize * Channels * sizeof(Real); }
template< int Channels , class Real >
void* MultiDataStreamBackedGrid< Channels , Real >::operator[]( int idx )
{
	if( _read && !_readComplete )
	{
		int offset = 0;
		for( int i=0 ; i<_streamCount ; i++ )
		{
			_streams[i]->read( _subData[i] , _rowSizes[i]*Channels*sizeof(Real) );
			if( _separateColors ) for( int c=0 ; c<Channels ; c++ ) memcpy( _data + offset + c*_rowSize , _subData[i]+c*_rowSizes[i] , sizeof(Real)*_rowSizes[i] );
			else memcpy( _data+offset*Channels , _subData[i] , _rowSizes[i]*Channels*sizeof(Real) );
			offset += _rowSizes[i];
		}
		_readComplete = true;
	}
#if ASSERT_MEMORY_ACCESS
	if( idx<0 || idx>=_rows || idx!=_current) fprintf( stderr , "Index out of bounds: %d != %d\n" , idx , _current ) , exit(0);
#endif // ASSERT_MEMORY_ACCESS
	return _data;
}
template< int Channels , class Real >
void MultiDataStreamBackedGrid< Channels , Real >::advance( void )
{
	if( !_read )
	{
		int offset = 0;
		for( int i=0 ; i<_streamCount ; i++ )
		{
			if( _separateColors )
				for( int c=0 ; c<Channels ; c++ ) memcpy( _subData[i]+c*_rowSizes[i] , _data + offset + c*_rowSize , sizeof(Real)*_rowSizes[i] );
			else memcpy( _subData[i] , _data+offset*Channels , _rowSizes[i]*Channels*sizeof(Real) );
			_streams[i]->write( _subData[i] , _rowSizes[i]*Channels*sizeof(Real) );
			offset += _rowSizes[i];
		}
	}
	else
	{
		if( !_readComplete )
		{
			int offset = 0;
			for( int i=0 ; i<_streamCount ; i++ )
			{
				if( _separateColors ) for( int c=0 ; c<Channels ; c++ ) memcpy( _subData[i]+c*_rowSizes[i] , _data + offset + c*_rowSize , sizeof(Real)*_rowSizes[i] );
				else memcpy( _subData[i] , _data+offset*Channels , _rowSizes[i]*Channels*sizeof(Real) );
				_streams[i]->read( _subData[i] , _rowSizes[i]*Channels*sizeof(Real) );
				offset += _rowSizes[i];
			}
		}
		_readComplete = false;
	}
	_current++;
}
template< int Channels , class Real >
void MultiDataStreamBackedGrid< Channels , Real >::reset( bool read , int minWindowSize )
{
	_read = read;
	_readComplete = !_read;
	_current = 0;
}
#else // !NEW_CONNECTION_BACKED_GRID
///////////////////////////
// MultiSocketBackedGrid //
///////////////////////////
template< int Channels , class Real >
MultiSocketBackedGrid< Channels , Real >::MultiSocketBackedGrid( SOCKET* socks , int* rowSizes , int sockCount , int rows , bool blockingSend , bool separateColors )
{
	_read = true;
	_separateColors = separateColors;
	_blockingSend = blockingSend;
	_sockCount = sockCount;
	_rows = rows;
	_rowSize = 0;
	for( int i=0 ; i<_sockCount ; i++ ) _rowSize += rowSizes[i];

	_rowSizes = ( int* )malloc( sizeof(int) * _sockCount );
	_socks = ( SOCKET* )malloc( sizeof(SOCKET) * _sockCount );
	_subData = AllocArray< Pointer( Real ) >( _sockCount , 1 , "MultiSocketBackedGrid::subData" );
	if( !_rowSizes || !_socks || !_subData ) fprintf( stderr , "Failed to allocate\n" ) ,  exit(0);
	for( int i=0 ; i<_sockCount ; i++ )
	{
		_rowSizes[i] = rowSizes[i];
		_socks[i] = socks[i];
		_subData[i] = AllocArray< Real >( Channels * rowSizes[i] , 1 , "MultiSocketBackedGrid::subData[]" );
		if( !_subData[i] ) fprintf( stderr , "Failed to allocate\n" ) , exit(0);
	}

	_data = AllocArray< Real >( _rowSize * Channels , 1 , "MultiSocketBackedGrid::data" );
	if( !_data ) fprintf( stderr , "Failed to allocate\n" ) , exit(0);
}
template< int Channels , class Real >
MultiSocketBackedGrid< Channels , Real >::~MultiSocketBackedGrid(void)
{
	for( int i=0 ; i<_sockCount ; i++ ) FreeArray( _subData[i] );
	FreeArray( _data );
	FreeArray( _subData );
	free( _socks );
	free( _rowSizes );

	_rows = _rowSize = 0;
	_rowSizes = NULL;
	_socks = NULL;
}
template< int Channels , class Real >
int	MultiSocketBackedGrid< Channels , Real >::rows( void ) const { return _rows; }
template< int Channels , class Real >
int MultiSocketBackedGrid< Channels , Real >::rowSize( void ) const { return _rowSize * Channels * sizeof(Real); }
template< int Channels , class Real >
Pointer( byte ) MultiSocketBackedGrid< Channels , Real >::operator[]( int idx )
{
#if SOCKET_DEBUG
	{
		MyWinSock::StdinLock lock;
		if( _read ) printfId( "Reading multi-socket[%d]\n" , _current ) , fflush( stdout );
		else		printfId( "Writing multi-socket[%d]\n" , _current ) , fflush( stdout );
		for( int i=0 ; i<_sockCount ; i++ )
			printf( "\tSocket[%d]: %d , port[%d]: %d , bytes[%d]: %d\n" , i , _socks[i] , i , GetLocalSocketPort( _socks[i] ) , i , _rowSizes[i]*Channels*sizeof(Real) ) , fflush( stdout );
	}
#endif // SOCKET_DEBUG
	if( _read && !_readComplete )
	{
		int offset = 0;
		for( int i=0 ; i<_sockCount ; i++ )
		{
			ReceiveOnSocket( _socks[i] , _subData[i] , _rowSizes[i]*Channels*sizeof(Real) , _blockingSend , "MultiSocketBackedGrid<Channels,Real>::[]" );
			if( _separateColors ) for( int c=0 ; c<Channels ; c++ ) memcpy( _data + offset + c*_rowSize , _subData[i]+c*_rowSizes[i] , sizeof(Real)*_rowSizes[i] );
			else memcpy( _data+offset*Channels , _subData[i] , _rowSizes[i]*Channels*sizeof(Real) );
			offset += _rowSizes[i];
#if MAJOR_SOCKET_DEBUG
			{
				MyWinSock::StdinLock lock;
				if( idx<MAJOR_SOCKET_DEBUG )
				{
					printf("Receiving[%d,%d]\t",idx,i);
					for( int k=0 ; k<MAJOR_SOCKET_DEBUG ; k++ ) printf("%d ",((int*)_subData[i])[k]);
					printf("\n") , fflush( stdout );
				}
			}
#endif // MAJOR_SOCKET_DEBUG
		}
		_readComplete = true;
	}
#if ASSERT_MEMORY_ACCESS
	if( idx<0 || idx>=_rows || idx!=_current) fprintf( stderr , "[MultiSocketBackedGrid] Index out of bounds: %d != %d || [ %d , %d )\n" , idx , _current , 0 , _rows ) , exit(0);
#endif // ASSERT_MEMORY_ACCESS
#if SOCKET_DEBUG
	{
		MyWinSock::StdinLock lock;
		if( _read ) printfId( "Done reading multi-socket\n") , fflush( stdout );
		else		printfId( "Done writing multi-socket\n") , fflush( stdout );
	}
#endif // SOCKET_DEBUG
	return ( Pointer( byte ) )_data;
}
template< int Channels , class Real >
void MultiSocketBackedGrid< Channels , Real >::advance( void )
{
#if SOCKET_DEBUG
	{
		MyWinSock::StdinLock lock;
		printfId("Advancing multi-row[%d]\n" , _current) , fflush( stdout );
		for( int i=0 ; i<_sockCount ; i++ )	printf("\tSocket[%d]: %d , port[%d]: %d\n" , i , _socks[i] , i , GetLocalSocketPort( _socks[i] ) ) , fflush( stdout );
	}
#endif // SOCKET_DEBUG

	if( !_read )
	{
		int offset = 0;
		for( int i=0 ; i<_sockCount ; i++ )
		{
			if( _separateColors ) for( int c=0 ; c<Channels ; c++ ) memcpy( _subData[i]+c*_rowSizes[i] , _data + offset + c*_rowSize , sizeof(Real)*_rowSizes[i] );
			else memcpy( _subData[i] , _data+offset*Channels , _rowSizes[i]*Channels*sizeof(Real) );
#if MAJOR_SOCKET_DEBUG
			{
				MyWinSock::StdinLock lock;
				if( _current<MAJOR_SOCKET_DEBUG )
				{
					printfId("Sending[%d,%d]\t",_current,i)
					for(int k=0;k<MAJOR_SOCKET_DEBUG;k++) printf("%d ",((int*)_subData[i])[k]);
					printf("\n") , fflush( stdout );
				}
			}
#endif // MAJOR_SOCKET_DEBUG
			SendOnSocket( _socks[i] , _subData[i] , _rowSizes[i]*Channels*sizeof(Real) , _blockingSend  , "MultiSocketBackedGrid< Channels,Real>::advance" );
			offset += _rowSizes[i];
		}
	}
	else
	{
		if( !_readComplete )
		{
			int offset = 0;
			for( int i=0 ; i<_sockCount ; i++ )
			{
				if( _separateColors ) for( int c=0 ; c<Channels ; c++ ) memcpy( _subData[i]+c*_rowSizes[i] , _data + offset + c*_rowSize , sizeof(Real)*_rowSizes[i] );
				else memcpy( _subData[i] , _data+offset*Channels , _rowSizes[i]*Channels*sizeof(Real) );
				ReceiveOnSocket( _socks[i] , _subData[i] , _rowSizes[i]*Channels*sizeof(Real) , _blockingSend , "MultiSocketBackedGrid<Channels,Real>::advance"  );
				offset += _rowSizes[i];
#if MAJOR_SOCKET_DEBUG
				{
					MyWinSock::StdinLock lock;
					if( _current<MAJOR_SOCKET_DEBUG )
					{
						printfId("Receiving[%d,%d]\t",_current,i);
						for(int k=0;k<MAJOR_SOCKET_DEBUG;k++) printf("%d ",((int*)_subData[i])[k]);
						printf("\n") , fflush( stdout );
					}
				}
#endif // MAJOR_SOCKET_DEBUG
			}
		}
		_readComplete = false;
	}
	_current++;
#if SOCKET_DEBUG
	{
		MyWinSock::StdinLock lock;
		printfId( "Done advancing multi-socket\n" ) , fflush( stdout );
	}
#endif // SOCKET_DEBUG
}
template< int Channels , class Real >
void MultiSocketBackedGrid< Channels , Real >::reset( bool read , int minWindowSize )
{
	if(!_read ) for( int i=0 ; i<_sockCount ; i++ ) EndSendOnSocket( _socks[i] , _blockingSend , "MultiSocketBackedGrid<Channels,Real>::reset" );
#if SOCKET_DEBUG
	{
		MyWinSock::StdinLock lock;
		if( read )	printfId( "Resetting multi-socket for read\n" ) , fflush( stdout );
		else		printfId( "Resetting multi-socket for write\n" ) , fflush( stdout );
		for( int i=0 ; i<_sockCount ; i++ )	printf("\tSocket[%d]: %d , port[%d]: %d\n" , i , _socks[i] , i , GetLocalSocketPort( _socks[i] ) ) , fflush( stdout );
	}
#endif // SOCKET_DEBUG
	_read = read;
	_readComplete = !_read;
	_current = 0;
	if( _read ) for( int i=0 ; i<_sockCount ; i++ ) StartReceiveOnSocket( _socks[i] , _blockingSend , "MultiSocketBackedGrid<Channels,Real >::reset" );
}
#endif // NEW_CONNECTION_BACKED_GRID
