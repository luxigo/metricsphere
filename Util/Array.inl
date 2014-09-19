#include <windows.h>
#define FULL_MISHA_DEBUG    0
#define ASSERT_ARRAY_ACCESS 1
#define ASSERT_VALID_VALUE  0



#include <stdio.h>
#include <emmintrin.h>
#include <vector>

inline bool isfinitef( float fp ){ float f=fp; return ((*(unsigned *)&f)&0x7f800000)!=0x7f800000; }


template< class C >        bool IsValid( const C& c );
#if _DEBUG
template< >         inline bool IsValid< float >( const float& f ) { return isfinitef( f ) &&  ( f==0.f || abs(f)>1e-31f ); }
#else // !_DEBUG
template< >         inline bool IsValid< float >( const float& f ) { return isfinitef( f ); }
#endif // _DEBUG
template< >         inline bool IsValid< __m128 >( const __m128& m )
{
	const __m128* addr = &m;
	if( size_t(addr) & 15 ) return false;
	else                    return true;
}
template< class C > inline bool IsValid( const C& c ){ return true; }


#if FULL_MISHA_DEBUG
class MemoryInfo
{
public:
	void* address;
	char name[512];
};
static std::vector< MemoryInfo > memoryInfo;
#endif // FULL_MISHA_DEBUG

template< class C >
class Array
{
	static int _count;
	void _assertBounds( int idx ) const
	{
		if( idx<min || idx>=max )
		{
			fprintf( stderr , "Array index out-of-bounds: %d <= %d < %d\n" , min , idx , max );
			ASSERT( 0 );
			exit( 0 );
		}
	}
protected:
	bool virt;
	C *data , *_data;
	int min , max;
public:
	int minimum( void ) const { return min; }
	int maximum( void ) const { return max; }

	inline Array( int size , int alignment , bool clear , bool virt=false , const char* name=NULL )
	{
		if( !_count ) printf( "Using Array instead of pointer\n" );
		this->virt = virt;
		if( virt ) _data = data = ( C* ) VirtualAlloc( NULL , size*sizeof(C) , MEM_RESERVE | MEM_COMMIT , PAGE_READWRITE );
#ifdef WIN32
		else       _data = data = ( C* ) _aligned_malloc( size * sizeof( C ) , alignment );
#else // !WIN32
		else       _data = data = ( C* ) memalign( size * sizeof( C ) , alignment );
#endif // WIN32
		if( clear ) memset( data ,  0 , size * sizeof( C ) );
		else        memset( data , -1 , size * sizeof( C ) );
		min = 0;
		max = size;
#if FULL_MISHA_DEBUG
		{
			MyWinSock::StdinLock lock;
			int sz = memoryInfo.size();
			memoryInfo.resize( sz + 1 );
			memoryInfo[sz].address = _data;
			if( name ) strcpy( memoryInfo[sz].name , name );
			else memoryInfo[sz].name[0] = 0;
		}
#endif // FULL_MISHA_DEBUG
		_count++;
	}
	Array( void )
	{
		if( !_count ) printf( "Using Array instead of pointer\n" );
		data = _data = NULL;
		min = max = 0;
		_count++;
	}
	template< class C > static Array FromPointer( C* data , int max )
	{
		Array a;
		a._data = NULL;
		a.data = data;
		a.min = 0;
		a.max = max;
		return a;
	}
	template< class C > static Array FromPointer( C* data , int min , int max )
	{
		Array a;
		a._data = NULL;
		a.data = data;
		a.min = min;
		a.max = max;
		return a;
	}
	inline C& operator[]( int idx )
	{
#if ASSERT_ARRAY_ACCESS
		_assertBounds( idx );
#endif // ASSERT_ARRAY_ACCESS
#if ASSERT_VALID_VALUE
		if( !IsValid( data[idx] ) )
		{
			fprintf( stderr , "Invalid value\n" );
			ASSERT( 0 );
			exit( 0 );
		}
#endif // ASSERT_VALID_VALUE
		return data[idx];
	}
	inline const C& operator[]( int idx ) const
	{
#if ASSERT_ARRAY_ACCESS
		_assertBounds( idx );
#endif // ASSERT_ARRAY_ACCESS
#if ASSERT_VALID_VALUE
		if( !IsValid( data[idx] ) )
		{
			fprintf( stderr , "Invalid value\n" );
			ASSERT( 0 );
			exit( 0 );
		}
#endif // ASSERT_VALID_VALUE
		return data[idx];
	}
	inline Array operator + ( int idx )
	{
#if ASSERT_ARRAY_ACCESS && 0 // Should we really be asserting on the pointer change, maybe just the derefrence
		_assertBounds( idx );
#endif // ASSERT_ARRAY_ACCESS
		Array a;
		a._data = _data;
		a.data = data+idx;
		a.min = min-idx;
		a.max = max-idx;
		return a;
	}
	inline Array& operator += ( int idx  )
	{
#if ASSERT_ARRAY_ACCESS && 0 // Should we really be asserting on the pointer change, maybe just the derefrence
		_assertBounds( idx );
#endif // ASSERT_ARRAY_ACCESS
		min -= idx;
		max -= idx;
		data += idx;
		return (*this);
	}
	Array  operator -  ( int idx ) { return (*this) +  (-idx); }
	Array& operator -= ( int idx ) { return (*this) += (-idx); }
	Array  operator +  ( unsigned int idx ) { return (*this) +  int(idx); }
	Array& operator += ( unsigned int idx ) { return (*this) += int(idx); }
	Array  operator -  ( unsigned int idx ) { return (*this) -  int(idx); }
	Array& operator -= ( unsigned int idx ) { return (*this) -= int(idx); }

	Array( Array& a )
	{
		if( !a )
		{
			data = _data =  NULL;
			min = max = 0;
		}
		else
		{
			_data = a._data;
			data = a.data;
			min = a.min;
			max = a.max;
		}
	}
	template< class D >
	Array( Array< D >& a )
	{
		_data = NULL;
		if( !a )
		{
			data =  NULL;
			min = max = 0;
		}
		else
		{
			int szC = sizeof( C );
			int szD = sizeof( D );
			data = (C*)&a[0];
			min = ( a.minimum() * szD ) / szC;
			max = ( a.maximum() * szD ) / szC;
			if( min*szC!=a.minimum()*szD || max*szC!=a.maximum()*szD )
			{
				fprintf( stderr , "Could not convert array [ %d , %d ] * %d => [ %d , %d ] * %d\n" , a.minimum() , a.maximum() , szD , min , max , szC );
				ASSERT( 0 );
				exit( 0 );
			}
		}
	}

	void Free( void )
	{
		if( _data )
		{
			if( virt ) VirtualFree( _data , 0 , MEM_RELEASE );
#ifdef WIN32
			else       _aligned_free( _data );
#else // !WIN32
			else       free( _data );
#endif // WIN32
#if FULL_MISHA_DEBUG
			{
				MyWinSock::StdinLock lock;
				int idx;
				for( idx=0 ; idx<memoryInfo.size( ) ; idx++ ) if( memoryInfo[idx].address==_data ) break;
				if( idx==memoryInfo.size() )
				{
					fprintf( stderr , "Could not find memory in address table\n" );
					ASSERT( 0 );
				}
				else
				{
					memoryInfo[idx] = memoryInfo[memoryInfo.size()-1];
					memoryInfo.pop_back( );
				}
			}
#endif // FULL_MISHA_DEBUG
		}
		(*this) = Array( );
	}
	C* pointer( void ){ return data; }
	bool operator !( void ) { return data==NULL; }
	operator bool( ) { return data!=NULL; }
};
template< class C > int Array< C >::_count = 0;

template< class C >
class ConstArray
{
	void _assertBounds( int idx ) const
	{
		if( idx<min || idx>=max )
		{
			fprintf( stderr , "Array index out-of-bounds: %d <= %d < %d\n" , min , idx , max );
			ASSERT( 0 );
			exit( 0 );
		}
	}
protected:
	const C *data;
	int min , max;
public:
	int minimum( void ) const { return min; }
	int maximum( void ) const { return max; }

	inline ConstArray( void )
	{
		data = NULL;
		min = max = 0;
	}
	inline ConstArray( Array< C >& a )
	{
		data = a.pointer( );
		min = a.minimum( );
		max = a.maximum( );
	}
	inline ConstArray( ConstArray& a )
	{
		data = a.pointer( );
		min = a.minimum( );
		max = a.maximum( );
	}
	template< class D >
	inline ConstArray( Array< D >& a )
	{
		int szC = sizeof( C );
		int szD = sizeof( D );
		data = ( const C*)a.pointer( );
		min = ( a.minimum() * szD ) / szC;
		max = ( a.maximum() * szD ) / szC;
		if( min*szC!=a.minimum()*szD || max*szC!=a.maximum()*szD )
		{
			fprintf( stderr , "Could not convert array [ %d , %d ] * %d => [ %d , %d ] * %d\n" , a.minimum() , a.maximum() , szD , min , max , szC );
			ASSERT( 0 );
			exit( 0 );
		}
/*
		data = ( const C* )a.pointer( );
		min = a.minimum( );
		max = a.maximum( );
*/
	}
	template< class D >
	inline ConstArray( ConstArray< D >& a )
	{
		int szC = sizeof( C );
		int szD = sizeof( D );
		data = ( const C*)a.pointer( );
		min = ( a.minimum() * szD ) / szC;
		max = ( a.maximum() * szD ) / szC;
		if( min*szC!=a.minimum()*szD || max*szC!=a.maximum()*szD )
		{
			fprintf( stderr , "Could not convert array [ %d , %d ] * %d => [ %d , %d ] * %d\n" , a.minimum() , a.maximum() , szD , min , max , szC );
			ASSERT( 0 );
			exit( 0 );
		}
		/*
		data = ( const C* )a.pointer( );
		min = a.minimum( );
		max = a.maximum( );
		*/
	}
	inline const C& operator[]( int idx ) const
	{
#if ASSERT_ARRAY_ACCESS
		_assertBounds( idx );
#endif // ASSERT_ARRAY_ACCESS
#if ASSERT_VALID_VALUE
		if( !IsValid( data[idx] ) )
		{
			fprintf( stderr , "Invalid value\n" );
			ASSERT( 0 );
			exit( 0 );
		}
#endif // ASSERT_VALID_VALUE
		return data[idx];
	}
	inline ConstArray operator + ( int idx ) const
	{
#if ASSERT_ARRAY_ACCESS
		_assertBounds( idx );
#endif // ASSERT_ARRAY_ACCESS
		ConstArray a;
		a.data = data+idx;
		a.min = min-idx;
		a.max = max-idx;
		return a;
	}
	inline ConstArray& operator += ( int idx  )
	{
#if ASSERT_ARRAY_ACCESS
		_assertBounds( idx );
#endif // ASSERT_ARRAY_ACCESS
		min -= idx;
		max -= idx;
		data += idx;
		return (*this);
	}
	ConstArray  operator -  (          int idx ) const { return (*this) +  (-idx); }
	ConstArray& operator -= (          int idx )       { return (*this) += (-idx); }
	ConstArray  operator +  ( unsigned int idx ) const { return (*this) +  int(idx); }
	ConstArray& operator += ( unsigned int idx )       { return (*this) += int(idx); }
	ConstArray  operator -  ( unsigned int idx ) const { return (*this) -  int(idx); }
	ConstArray& operator -= ( unsigned int idx )       { return (*this) -= int(idx); }

	const C* pointer( void ){ return data; }
	bool operator !( void ) { return data==NULL; }
	operator bool( ) { return data!=NULL; }
};

#if FULL_MISHA_DEBUG
inline void PrintMemoryInfo( void ){ for( int i=0 ; i<memoryInfo.size() ; i++ ) printf( "%d] %s\n" , i , memoryInfo[i].name ); }
#endif // FULL_MISHA_DEBUG
template< class C >
Array< C >& memcpy( Array< C >& destination , const void* source , size_t size )
{
	if( size>destination.maximum()*sizeof(C) )
	{
		fprintf( stderr , "Size of copy exceeds destination maximum: %d > %d\n" , size , destination.maximum()*sizeof( C ) );
		ASSERT( 0 );
		exit( 0 );
	}
	memcpy( &destination[0] , source , size );
	return destination;
}
template< class C , class D >
Array< C >& memcpy( Array< C >& destination , Array< D >& source , size_t size )
{
	if( size>destination.maximum()*sizeof( C ) )
	{
		fprintf( stderr , "Size of copy exceeds destination maximum: %d > %d\n" , size , destination.maximum()*sizeof( C ) );
		ASSERT( 0 );
		exit( 0 );
	}
	if( size>source.maximum()*sizeof( D ) )
	{
		fprintf( stderr , "Size of copy exceeds source maximum: %d > %d\n" , size , source.maximum()*sizeof( D ) );
		ASSERT( 0 );
		exit( 0 );
	}
	memcpy( &destination[0] , &source[0] , size );
	return destination;
}
template< class C , class D >
Array< C >& memcpy( Array< C >& destination , ConstArray< D >& source , size_t size )
{
	if( size>destination.maximum()*sizeof( C ) )
	{
		fprintf( stderr , "Size of copy exceeds destination maximum: %d > %d\n" , size , destination.maximum()*sizeof( C ) );
		ASSERT( 0 );
		exit( 0 );
	}
	if( size>source.maximum()*sizeof( D ) )
	{
		fprintf( stderr , "Size of copy exceeds source maximum: %d > %d\n" , size , source.maximum()*sizeof( D ) );
		ASSERT( 0 );
		exit( 0 );
	}
	memcpy( &destination[0] , &source[0] , size );
	return destination;
}
template< class C , class D >
Array< C >& memcpy( Array< C >& destination , const Array< D >& source , size_t size )
{
	if( size>destination.maximum()*sizeof( C ) )
	{
		fprintf( stderr , "Size of copy exceeds destination maximum: %d > %d\n" , size , destination.maximum()*sizeof( C ) );
		ASSERT( 0 );
		exit( 0 );
	}
	if( size>source.maximum()*sizeof( D ) )
	{
		fprintf( stderr , "Size of copy exceeds source maximum: %d > %d\n" , size , source.maximum()*sizeof( D ) );
		ASSERT( 0 );
		exit( 0 );
	}
	memcpy( &destination[0] , &source[0] , size );
	return destination;
}
template< class D >
void* memcpy( void* destination , Array< D >& source , size_t size )
{
	if( size>source.maximum()*sizeof( D ) )
	{
		fprintf( stderr , "Size of copy exceeds source maximum: %d > %d\n" , size , source.maximum()*sizeof( D ) );
		ASSERT( 0 );
		exit( 0 );
	}
	memcpy( destination , &source[0] , size );
	return destination;
}
template< class D >
void* memcpy( void* destination , ConstArray< D >& source , size_t size )
{
	if( size>source.maximum()*sizeof( D ) )
	{
		fprintf( stderr , "Size of copy exceeds source maximum: %d > %d\n" , size , source.maximum()*sizeof( D ) );
		ASSERT( 0 );
		exit( 0 );
	}
	memcpy( destination , &source[0] , size );
	return destination;
}
template< class C >
Array< C >& memset( Array< C >& destination , int value , size_t size )
{
	if( size>destination.maximum()*sizeof( C ) )
	{
		fprintf( stderr , "Size of copy exceeds destination maximum: %d > %d\n" , size , destination.maximum()*sizeof( C ) );
		ASSERT( 0 );
		exit( 0 );
	}
	memset( &destination[0] , value , size );
	return destination;
}

template< class C >
size_t fread( Array< C > destination , size_t eSize , size_t count , FILE* fp )
{
	if( count*eSize>destination.maximum()*sizeof( C ) )
	{
		fprintf( stderr , "Size of read exceeds source maximum: %d > %d\n" , count*eSize , destination.maximum()*sizeof( C ) );
		ASSERT( 0 );
		exit( 0 );
	}
	return fread( &destination[0] , eSize , count , fp );
}
template< class C >
size_t fwrite( Array< C > source , size_t eSize , size_t count , FILE* fp )
{
	if( count*eSize>source.maximum()*sizeof( C ) )
	{
		fprintf( stderr , "Size of write exceeds source maximum: %d > %d\n" , count*eSize , source.maximum()*sizeof( C ) );
		ASSERT( 0 );
		exit( 0 );
	}
	return fwrite( &source[0] , eSize , count , fp );
}
template< class C >
size_t fwrite( ConstArray< C > source , size_t eSize , size_t count , FILE* fp )
{
	if( count*eSize>source.maximum()*sizeof( C ) )
	{
		fprintf( stderr , "Size of write exceeds source maximum: %d > %d\n" , count*eSize , source.maximum()*sizeof( C ) );
		ASSERT( 0 );
		exit( 0 );
	}
	return fwrite( &source[0] , eSize , count , fp );
}
template< class C >
BOOL ReadFile( HANDLE hFile , Array< C > lpBuffer , DWORD nNumberOfBytesToRead , LPDWORD lpNumberOfBytesRead , LPOVERLAPPED lpOverlapped )
{
	if( lpBuffer.maximum()<nNumberOfBytesToRead )
	{
		fprintf( stderr , "Cannot read from file: %d < %d\n" , lpBuffer.maximum() , nNumberOfBytesToRead );
		ASSERT( 0 );
		exit( 0 );
	}
	return ReadFile( hFile , (LPVOID)PointerAddress< C >( lpBuffer ) , nNumberOfBytesToRead , lpNumberOfBytesRead , lpOverlapped );
}
template< class C >
BOOL WriteFile( HANDLE hFile , Array< C > lpBuffer , DWORD nNumberOfBytesToWrite , LPDWORD lpNumberOfBytesWritten , LPOVERLAPPED lpOverlapped )
{
	if( lpBuffer.maximum()<nNumberOfBytesToWrite )
	{
		fprintf( stderr , "Cannot write to file: %d < %d\n" , lpBuffer.maximum() , nNumberOfBytesToWrite );
		ASSERT( 0 );
		exit( 0 );
	}
	return WriteFile( hFile , (LPCVOID)PointerAddress< C >( lpBuffer ) , nNumberOfBytesToWrite , lpNumberOfBytesWritten , lpOverlapped );
}
template< class C >
BOOL WriteFile( HANDLE hFile , ConstArray< C > lpBuffer , DWORD nNumberOfBytesToWrite , LPDWORD lpNumberOfBytesWritten , LPOVERLAPPED lpOverlapped )
{
	if( lpBuffer.maximum()<nNumberOfBytesToWrite )
	{
		fprintf( stderr , "Cannot write to file: %d < %d\n" , lpBuffer.maximum() , nNumberOfBytesToWrite );
		ASSERT( 0 );
		exit( 0 );
	}
	return WriteFile( hFile , (LPCVOID)PointerAddress< C >( lpBuffer ) , nNumberOfBytesToWrite , lpNumberOfBytesWritten , lpOverlapped );
}
template< class C >
int recv( SOCKET& s , Array< C > destination , int len , int flags )
{
	if( len>destination.maximum()*sizeof( C ) )
	{
		fprintf( stderr , "Size of recv exceeds destination maximum: %d > %d\n" , len , destination.maximum()*sizeof( C ) );
		ASSERT( 0 );
		exit( 0 );
	}
	return recv( s , (char*)&destination[0] , len , flags );
}
template< class C >
int send ( SOCKET s , Array< C > source , int len , int flags )
{
	if( len>source.maximum()*sizeof( C ) )
	{
		fprintf( stderr , "Size of send exceeds source maximum: %d > %d\n" , len , source.maximum()*sizeof( C ) );
		ASSERT( 0 );
		exit( 0 );
	}
	return send( s , (char*)&source[0] , len , flags );
}
template< class C >
int send ( SOCKET s , ConstArray< C > source , int len , int flags )
{
	if( len>source.maximum()*sizeof( C ) )
	{
		fprintf( stderr , "Size of send exceeds source maximum: %d > %d\n" , len , source.maximum()*sizeof( C ) );
		ASSERT( 0 );
		exit( 0 );
	}
	return send( s , (char*)&source[0] , len , flags );
}
