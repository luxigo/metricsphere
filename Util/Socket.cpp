/*
Copyright (c) 2006, Michael Kazhdan and Matthew Bolitho
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution. 

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#include "Socket.h"

int MyWinSock::_wsaCount = 0;
WSADATA MyWinSock::_wsaData;
CRITICAL_SECTION MyWinSock::stdinLock;
CRITICAL_SECTION MyWinSock::stderrLock;
CRITICAL_SECTION MyWinSock::systemLock;
#if STORE_CONNECTION_TABLE
hash_map< int , ConnectionData > MyWinSock::ConnectionTable;
#endif // STORE_CONNECTION_TABLE

void printfId(const char* format,...)
{
	va_list args;
	va_start(args,format);
	printf( "%d] " , GetCurrentThreadId() );
	vprintf(format,args);
	va_end(args);
}
void fprintfId(FILE* fp , const char* format,...)
{
	va_list args;
	va_start( args , format );
	fprintf( fp , "%d] " , GetCurrentThreadId() );
	vfprintf( fp , format , args );
	va_end( args);
}

bool MyWaitForSingleObject( HANDLE hHandle , DWORD milliseconds , const char* deadlockMessage , ... )
{
#if DEBUG_DEADLOCK
	int state;
	while( (state = WaitForSingleObject( hHandle , milliseconds ) ) == WAIT_TIMEOUT )
	{
		MyWinSock::fprintfId( stderr , "Possible single deadlock: " );
		MyWinSock::StderrLock lock;
		va_list args;
		va_start( args , deadlockMessage );
		vfprintf( stderr , deadlockMessage , args );
		va_end( args );
		fprintf( stderr , "\n" );
	}
	if( state !=  WAIT_OBJECT_0 )
	{
		MyWinSock::fprintfId( stderr , "Failed to wait for single event: " ) , PrintError();
		return false;
	}
#else // !DEBUG_DEADLOCK
	if( WaitForSingleObject( hHandle , INFINITE ) != WAIT_OBJECT_0 )
	{
		MyWinSock::StderrLock lock;
		fprintf( stderr , "Failed to wait for single event: " ) , PrintError();
		return false;
	}
#endif // DEADLOCK
	return true;
}
bool MyWaitForMultipleObjects( DWORD nCount , const HANDLE* lpHandles , DWORD milliseconds , const char* deadlockMessage , ... )
{
#if DEBUG_DEADLOCK
	int state;
	while( (state = WaitForMultipleObjects( nCount , lpHandles , TRUE , milliseconds ) ) == WAIT_TIMEOUT )
	{
		MyWinSock::fprintfId( stderr , "Possible multiple deadlock: " );
		MyWinSock::StderrLock lock;
		va_list args;
		va_start( args , deadlockMessage );
		vfprintf( stderr , deadlockMessage , args );
		va_end( args );
		fprintf( stderr , "\n" );
	}
	if( state !=  WAIT_OBJECT_0 )
	{
		MyWinSock::fprintfId( stderr , "Failed to wait for multiple events: " ) , PrintError();
		MyWinSock::StderrLock lock;
		va_list args;
		va_start( args , deadlockMessage );
		vfprintf( stderr , deadlockMessage , args );
		va_end( args );
		fprintf( stderr , "\n" );
		return false;
	}
#else // !DEBUG_DEADLOCK
	if( WaitForMultipleObjects( nCount , lpHandles , TRUE , INFINITE ) != WAIT_OBJECT_0 )
	{
		MyWinSock::StderrLock lock;
		fprintf( stderr , "Failed to wait for multiple events: " ) , PrintError();
		return false;
	}
#endif // DEADLOCK
	return false;
}


void MyWinSock::Load( void )
{
	if(!_wsaCount)
	{
		InitializeCriticalSection( &stdinLock );
		InitializeCriticalSection( &stderrLock );
		InitializeCriticalSection( &systemLock );
		if( WSAStartup( MAKEWORD(2,2), &_wsaData ) ) fprintfId( stderr , "WSAStartup failed: %s\n", LastSocketError() ) , exit(0);
	}
	_wsaCount++;
}
void MyWinSock::UnLoad( void )
{
	_wsaCount--;
	if(!_wsaCount)
	{
		WSACleanup();
		DeleteCriticalSection( &stdinLock );
		DeleteCriticalSection( &stderrLock );
		DeleteCriticalSection( &systemLock );
	}
}
void MyWinSock::fprintfId( FILE* fp , const char* format,...)
{
	EnterCriticalSection( &stdinLock );
	va_list args;
	va_start(args,format);
	fprintf( fp , "%d] " , GetCurrentThreadId() );
	vfprintf( fp , format , args );
	va_end(args);
	LeaveCriticalSection( &stdinLock );
}

void MyWinSock::printfId(const char* format,...)
{
	EnterCriticalSection( &stdinLock );
	va_list args;
	va_start(args,format);
	printf( "%d] " , GetCurrentThreadId() );
	vprintf(format,args);
	va_end(args);
	LeaveCriticalSection( &stdinLock );
}

//////////////////
void StartReceiveOnSocket( SOCKET& s , bool blockingSend , const char* errorMessage , ... )
{
	if( blockingSend )
	{
#if DEBUG_SOCKET
		printfId( " Sending Acknowledgement (%d)\n" , s );
#endif // DEBUG_SOCKET
		int ack;
		if( !SendOnSocket( s , GetPointer( ack ) , sizeof(ack) ) )
		{
			fprintfId( stderr , "Failed to send acknowledgement (%d)\n" , s );
			{
				MyWinSock::StderrLock lock;
				va_list args;
				va_start( args , errorMessage );
				vfprintf( stderr , errorMessage , args );
				va_end( args );
				fprintf( stderr , "\n" );
			}
			exit(0);
		}
#if DEBUG_SOCKET
		printfId( " Done (%d)\n" , s ) , fflush( stdout );
#endif // DEBUG_SOCKET
	}
}

void EndSendOnSocket( SOCKET& s , bool blockingSend , const char* errorMessage , ... )
{
	if( blockingSend )
	{
#if DEBUG_SOCKET
		printfId( " Receiving Acknowledgement (%d)\n" , s ) , fflush( stdout );
#endif // DEBUG_SOCKET
		int ack;
		if( !ReceiveOnSocket( s , GetPointer( ack ) , sizeof(ack) ) )
		{
			fprintfId( stderr , "Failed to receive acknowledgement (%d)\n" , s );
			{
				MyWinSock::StderrLock lock;
				va_list args;
				va_start( args , errorMessage );
				vfprintf( stderr , errorMessage , args );
				va_end( args );
				fprintf( stderr , "\n" );
			}
			exit(0);
		}
#if DEBUG_SOCKET
		printfId( " Done (%d)\n" , s ) , fflush( stdout );
#endif // DEBUG_SOCKET
	}
}

bool StartReceiveOnSocket( SOCKET& s , bool blockingSend )
{
	if( blockingSend )
	{
#if DEBUG_SOCKET
		printfId( " Sending Acknowledgement (%d)\n" , s ) , fflush( stdout );
#endif // DEBUG_SOCKET
		int ack;
		if( !SendOnSocket( s , GetPointer( ack ) , sizeof(ack) ) )
		{
			fprintfId( stderr , "Failed to send acknowledgement (%d)\n" , s );
			return false;
		}
#if DEBUG_SOCKET
		printfId( " Done (%d)\n" , s ) , fflush( stdout );
#endif // DEBUG_SOCKET
	}
	return true;
}

bool EndSendOnSocket( SOCKET& s , bool blockingSend )
{
	if( blockingSend )
	{
#if DEBUG_SOCKET
		printfId( " Receiving Acknowledgement (%d)\n" , s ) , fflush( stdout );
#endif // DEBUG_SOCKET
		int ack;
		if( !ReceiveOnSocket( s , GetPointer( ack ) , sizeof(ack) ) )
		{
			fprintfId( stderr , "Failed to receive acknowledgement (%d)\n" , s );
			return false;
		}
#if DEBUG_SOCKET
		printfId( " Done (%d)\n" , s ) , fflush( stdout );
#endif // DEBUG_SOCKET
	}
	return true;
}

void PrintHostAddress( void )
{
	MyWinSock::SystemLock lock;
	char hostName[512];
	gethostname( hostName , 512 );
	hostent* host = gethostbyname( hostName );
	for( int i=0 ; ; i++ )
		if( host->h_addr_list[i] == NULL ) break;
		else printf( "Address[%d]: %s\n" , i , inet_ntoa( *(struct in_addr*)host->h_addr_list[i] ) ) , fflush( stdout );
	for( int i=0 ; ; i++ )
		if( host->h_aliases[i] == NULL ) break;
		else printf( "Aliases[%d]: %s\n" , i , host->h_aliases[i] ) , fflush( stdout );
}
void GetHostAddress( char* address , char* prefix)
{
	char hostName[512];
	gethostname( hostName , 512 );
	{
		MyWinSock::SystemLock lock;
		hostent* host = gethostbyname( hostName );
		if( !prefix )
		{
			strcpy( address , inet_ntoa(*(struct in_addr*)host->h_addr) );
			return;
		}
		for( int i=0 ; ; i++ )
			if( host->h_addr_list[i] == NULL ) break;
			else if( strstr( inet_ntoa( *(struct in_addr*)host->h_addr_list[i] ) , prefix ) )
			{
				strcpy( address , inet_ntoa(*(struct in_addr*)host->h_addr_list[i]) );
				return;
			}
		strcpy( address , inet_ntoa(*(struct in_addr*)host->h_addr) );
	}
}
int GetLocalSocketPort( SOCKET& s )
{
    struct sockaddr_in local;
	int len=sizeof(local);
	if( getsockname ( s , (struct sockaddr*) &local , &len ) == SOCKET_ERROR )
	{
		fprintfId( stderr , "Error at getsockname(): %s\n" , LastSocketError() );
		return -1;
	}
	return int(ntohs(local.sin_port));
}
const char* GetLocalSocketAddress( SOCKET& s )
{
    struct sockaddr_in local;
	int len=sizeof(local);
	if( getsockname ( s , (struct sockaddr*) &local , &len ) == SOCKET_ERROR )
	{
		fprintfId( stderr , "Error at getsockname(): %s\n" , LastSocketError() );
		return NULL;
	}
	return inet_ntoa( local.sin_addr );
}
int GetPeerSocketPort( SOCKET& s )
{
    struct sockaddr_in peer;
	int len = sizeof( peer );
	if( getpeername ( s , (struct sockaddr*) &peer , &len ) == SOCKET_ERROR )
	{
		fprintfId( stderr , "Error at getpeername(): %s\n" , LastSocketError() );
		return -1;
	}
	return int(ntohs( peer.sin_port) );
}
const char* GetPeerSocketAddress( SOCKET& s )
{
    struct sockaddr_in peer;
	int len=sizeof( peer );
	if( getpeername ( s , (struct sockaddr*) &peer , &len ) == SOCKET_ERROR )
	{
		fprintfId( stderr , "Error at getpeername(): %s\n" , LastSocketError() );
		return NULL;
	}
	return inet_ntoa( peer.sin_addr );
}
SOCKET GetConnectSocket( const char* address , int port , int ms , bool progress )
{
	in_addr addr;
	inet_aton( address , &addr );
	return GetConnectSocket( addr , port , ms , progress );
}
SOCKET GetConnectSocket( in_addr address , int port , int ms , bool progress )
{
    struct sockaddr_in addr_in;
	memset( &addr_in, 0, sizeof(addr_in) );
	addr_in.sin_family = AF_INET;
	addr_in.sin_addr = address;
	addr_in.sin_port= htons ( port );

	SOCKET sock = socket( AF_INET, SOCK_STREAM , 0);
	if ( sock == INVALID_SOCKET )
	{
		fprintfId( stderr , "Error at GetConnectSocket( ... , %d ): %s\n" , port , LastSocketError() );
		return INVALID_SOCKET;
	}
	long long sleepCount = 0;
	while (connect( sock, (const sockaddr*)&addr_in, sizeof(addr_in) ) == SOCKET_ERROR)
	{
		sleepCount++;
		Sleep( 1 );
		if( progress && !(sleepCount%ms) ) printf( "." );
	}
	if( progress ) printf( "\n" ) , fflush( stdout );
	int val = 1;
	setsockopt( sock , IPPROTO_TCP , TCP_NODELAY , (char*)&val , sizeof(val) );
#if STORE_CONNECTION_TABLE
	if( sock!=INVALID_SOCKET )
	{
		MyWinSock::SystemLock lock;
		ConnectionData& data = MyWinSock::ConnectionTable[ sock ];
	    struct sockaddr_in addr;
		int len = sizeof( addr );
		if( getsockname ( sock , (struct sockaddr*) &addr , &len ) == SOCKET_ERROR )
		{
			fprintfId( stderr , "Error at getsockname(): %s\n" , LastSocketError() );
			return NULL;
		}
		memcpy( &data.localAddr , &addr.sin_addr , sizeof( addr.sin_addr ) );
		data.localPort = int( ntohs( addr.sin_port ) );
		if( getpeername ( sock , (struct sockaddr*) &addr , &len ) == SOCKET_ERROR )
		{
			fprintfId( stderr , "Error at getpeername(): %s\n" , LastSocketError() );
			return NULL;
		}
		memcpy( &data.peerAddr , &addr.sin_addr , sizeof( addr.sin_addr ) );
		data.peerPort = int( ntohs( addr.sin_port ) );
	}
#endif // STORE_CONNECTION_TABLE
	return sock;
}
SOCKET AcceptSocket( SOCKET listen )
{
	SOCKET sock = accept( listen , NULL , NULL );
	if ( sock == INVALID_SOCKET )
	{
		fprintfId( stderr , "accept failed: %s\n" , LastSocketError() );
		return INVALID_SOCKET;
	}
	int val = 1;
	setsockopt( sock , IPPROTO_TCP , TCP_NODELAY , (char*)&val , sizeof(val) );
#if STORE_CONNECTION_TABLE
	if( sock!=INVALID_SOCKET )
	{
		MyWinSock::SystemLock lock;
		ConnectionData& data = MyWinSock::ConnectionTable[ sock ];
	    struct sockaddr_in addr;
		int len = sizeof( addr );
		if( getsockname ( sock , (struct sockaddr*) &addr , &len ) == SOCKET_ERROR )
		{
			fprintfId( stderr , "Error at getsockname(): %s\n" , LastSocketError() );
			return NULL;
		}
		memcpy( &data.localAddr , &addr.sin_addr , sizeof( addr.sin_addr ) );
		data.localPort = int( ntohs( addr.sin_port ) );
		if( getpeername ( sock , (struct sockaddr*) &addr , &len ) == SOCKET_ERROR )
		{
			fprintfId( stderr , "Error at getpeername(): %s\n" , LastSocketError() );
			return NULL;
		}
		memcpy( &data.peerAddr , &addr.sin_addr , sizeof( addr.sin_addr ) );
		data.peerPort = int( ntohs( addr.sin_port ) );
	}
#endif // STORE_CONNECTION_TABLE
	return sock;
}

SOCKET GetListenSocket( int& port )
{
	SOCKET listenSocket = socket(AF_INET, SOCK_STREAM, 0);
	if (listenSocket == INVALID_SOCKET)
	{
		fprintfId( stderr , "Error at socket(): %s\n", LastSocketError());
		return INVALID_SOCKET;
	}

    struct sockaddr_in local;
    memset(&local, 0, sizeof(local));
	local.sin_addr.s_addr = htonl(INADDR_ANY);
    local.sin_port = htons(port);
    local.sin_family = AF_INET;

	// Setup the TCP listening socket
	if (bind( listenSocket, (const sockaddr*)&local, sizeof(local) ) == SOCKET_ERROR)
	{
		fprintfId( stderr , "bind failed: %s\n" , LastSocketError());
		closesocket(listenSocket);
		return INVALID_SOCKET;
	}

	if ( listen( listenSocket, SOMAXCONN ) == SOCKET_ERROR )
	{
		fprintfId( stderr , "Error at listen(): %s\n" , LastSocketError() );
		closesocket(listenSocket);
		return INVALID_SOCKET;
	}
	int len=sizeof(local);
	if(getsockname(listenSocket,(struct sockaddr*)&local,&len) == SOCKET_ERROR)
	{
		fprintfId( stderr , "Error at getsockname(): %s\n" , LastSocketError() );
		closesocket(listenSocket);
		return INVALID_SOCKET;
	}
	port=int(ntohs(local.sin_port));
	return listenSocket;
}
void CloseSocket( SOCKET& s )
{
#if STORE_CONNECTION_TABLE
	if( s!=INVALID_SOCKET )
	{
		MyWinSock::SystemLock lock;
		MyWinSock::ConnectionTable.erase( s );
	}
#endif // STORE_CONNECTION_TABLE

	if( s!=INVALID_SOCKET ) closesocket( s );
	s = INVALID_SOCKET;
}

int inet_aton(const char *cp, struct in_addr *inp)
{
    unsigned int a = 0, b = 0, c = 0, d = 0;
    int n = 0, r;
    unsigned long int addr = 0;
    r = sscanf(cp, "%u.%u.%u.%u%n", &a, &b, &c, &d, &n);
    if (r == 0 || n == 0) return 0;
    cp += n;
    if (*cp) return 0;
    if (a > 255 || b > 255 || c > 255 || d > 255) return 0;
    if (inp) {
        addr += a; addr <<= 8;
        addr += b; addr <<= 8;
        addr += c; addr <<= 8;
        addr += d;
        inp->s_addr = htonl(addr);
    }
    return 1;
}

static const char *wstrerror(int err)
{
    switch (err)
	{
        case WSAEINTR: return "Interrupted function call";
        case WSAEACCES: return "Permission denied";
        case WSAEFAULT: return "Bad address";
        case WSAEINVAL: return "Invalid argument";
        case WSAEMFILE: return "Too many open files";
        case WSAEWOULDBLOCK: return "Resource temporarily unavailable";
        case WSAEINPROGRESS: return "Operation now in progress";
        case WSAEALREADY: return "Operation already in progress";
        case WSAENOTSOCK: return "Socket operation on nonsocket";
        case WSAEDESTADDRREQ: return "Destination address required";
        case WSAEMSGSIZE: return "Message too long";
        case WSAEPROTOTYPE: return "Protocol wrong type for socket";
        case WSAENOPROTOOPT: return "Bad protocol option";
        case WSAEPROTONOSUPPORT: return "Protocol not supported";
        case WSAESOCKTNOSUPPORT: return "Socket type not supported";
        case WSAEOPNOTSUPP: return "Operation not supported";
        case WSAEPFNOSUPPORT: return "Protocol family not supported";
        case WSAEAFNOSUPPORT: return "Address family not supported by protocol family"; 
        case WSAEADDRINUSE: return "Address already in use";
        case WSAEADDRNOTAVAIL: return "Cannot assign requested address";
        case WSAENETDOWN: return "Network is down";
        case WSAENETUNREACH: return "Network is unreachable";
        case WSAENETRESET: return "Network dropped connection on reset";
        case WSAECONNABORTED: return "Software caused connection abort";
        case WSAECONNRESET: return "Connection reset by peer";
        case WSAENOBUFS: return "No buffer space available";
        case WSAEISCONN: return "Socket is already connected";
        case WSAENOTCONN: return "Socket is not connected";
        case WSAESHUTDOWN: return "Cannot send after socket shutdown";
        case WSAETIMEDOUT: return "Connection timed out";
        case WSAECONNREFUSED: return "Connection refused";
        case WSAEHOSTDOWN: return "Host is down";
        case WSAEHOSTUNREACH: return "No route to host";
        case WSAEPROCLIM: return "Too many processes";
        case WSASYSNOTREADY: return "Network subsystem is unavailable";
        case WSAVERNOTSUPPORTED: return "Winsock.dll version out of range";
        case WSANOTINITIALISED: return "Successful WSAStartup not yet performed";
        case WSAEDISCON: return "Graceful shutdown in progress";
        case WSAHOST_NOT_FOUND: return "Host not found";
        case WSATRY_AGAIN: return "Nonauthoritative host not found";
        case WSANO_RECOVERY: return "Nonrecoverable name lookup error"; 
        case WSANO_DATA: return "Valid name, no data record of requested type";
        default: return "Unknown error";
    }
}
//static const char* LastSocketError(void) { return wstrerror( WSAGetLastError() ); }
const char* LastSocketError(void) { return wstrerror( WSAGetLastError() ); }
///////////////////////////
// DataStreamConstructor //
///////////////////////////
DataStreamConstructor::DataStreamConstructor( void )
{
	_sock = INVALID_SOCKET;
	_stream = NULL;
	char address[512];
	_myPID = GetCurrentProcessId( );
	GetHostAddress( address );
	inet_aton( address, &_myAddr );
}
void DataStreamConstructor::init( SOCKET sock , bool master , bool cleanUp )
{
	_sock = sock;
	_master = master;
	_cleanUp = cleanUp;
}
DataStream* DataStreamConstructor::getDataStream( void ) { return _stream; }
void DataStreamConstructor::doStep( int sNum )
{
	switch( sNum )
	{
	case 0:
		if( !_master )
		{
			SendOnSocket( _sock , GetPointer(_myAddr) , sizeof(_myAddr) );
			SendOnSocket( _sock , GetPointer(_myPID ) , sizeof(_myPID) );
		}
		break;
	case 1:
		if( _master )
		{
			in_addr addr;
			int pid;
			ReceiveOnSocket( _sock , GetPointer(addr) , sizeof(addr) );
			ReceiveOnSocket( _sock , GetPointer(pid ) , sizeof(pid) );

			if( _myAddr.s_addr == addr.s_addr && _myPID == pid )
			{
				SharedMemoryBuffer::StreamPair sPair;
				SharedMemoryBuffer::StreamPair::CreateSharedBufferPair( sPair );
				SendOnSocket( _sock , GetPointer(sPair.second) , sizeof( sPair.second ) );
				_stream = sPair.first;
				if( _cleanUp ) CloseSocket( _sock );
			}
			else
			{
				SharedMemoryBuffer::SecondStream* sStream = NULL;
				SendOnSocket( _sock , GetPointer(sStream) , sizeof( sStream ) );
				_stream = new Socket( _sock );
			}
		}
		break;
	case 2:
		if( !_master )
		{
			SharedMemoryBuffer::SecondStream* sStream;
			ReceiveOnSocket( _sock , GetPointer(sStream) , sizeof( sStream ) );
			if( sStream )
			{
				if( _cleanUp ) CloseSocket( _sock );
				_stream = sStream;
			}
			else _stream = new Socket( _sock  );
		}
		break;
	}
}


////////////////
// DataStream //
////////////////
DataStream* DataStream::GetDataStream( SOCKET sock , bool master , bool cleanUp )
{
	char address[512];
	in_addr myAddr , addr;
	int pid , myPID = GetCurrentProcessId( );
	GetHostAddress( address );
	inet_aton( address, &myAddr );

	if( master )
	{
		ReceiveOnSocket( sock , GetPointer(addr) , sizeof(addr) );
		ReceiveOnSocket( sock , GetPointer(pid ) , sizeof(pid) );
		if( myAddr.s_addr == addr.s_addr && myPID == pid )
		{
			SharedMemoryBuffer::StreamPair sPair;
			SharedMemoryBuffer::StreamPair::CreateSharedBufferPair( sPair );
			SendOnSocket( sock , GetPointer(sPair.second) , sizeof( sPair.second ) );
			if( cleanUp ) CloseSocket( sock );
			return sPair.first;
		}
		else
		{
			SharedMemoryBuffer::SecondStream* sStream = NULL;
			SendOnSocket( sock , GetPointer(sStream) , sizeof( sStream ) );
			return new Socket( sock );
		}
	}
	else
	{
		SendOnSocket( sock , GetPointer(myAddr) , sizeof(myAddr) );
		SendOnSocket( sock , GetPointer(myPID ) , sizeof(myPID) );
		SharedMemoryBuffer::SecondStream* sStream;
		ReceiveOnSocket( sock , GetPointer(sStream) , sizeof( sStream ) );
		if( sStream )
		{
			if( cleanUp ) CloseSocket( sock );
			return sStream;
		}
		else return new Socket( sock  );
	}
}
////////////
// Socket //
////////////
Socket::Socket( SOCKET sock ) { _sock = sock; _mySocket = false; }
Socket::Socket( const char* address , int port , int ms , bool progress )
{
	_sock = GetConnectSocket( address , port , ms , progress );
	_mySocket = true;
}
Socket::Socket( in_addr address , int port , int ms , bool progress )
{
	_sock = GetConnectSocket( address , port , ms , progress );
	_mySocket = true;
}
Socket::~Socket( void ) { if( _mySocket ) CloseSocket( _sock); }

bool Socket::write( ConstPointer( byte ) buf , int len ) { return SendOnSocket   ( _sock , buf , len ); }
bool Socket::read ( Pointer(       byte ) buf , int len ){ return ReceiveOnSocket( _sock , buf , len ); }

////////////////////////
// SharedMemoryBuffer //
////////////////////////
SharedMemoryBuffer::SharedMemoryBuffer( void )
{
	_buf1 = NullPointer< byte >( );
	_buf2 = NullPointer< byte >( );
	_bufSize1 = _bufSize2 = 0;
	_loadHandle1   = CreateEvent( NULL , false , false , NULL );
	_loadHandle2   = CreateEvent( NULL , false , false , NULL );
	_unloadHandle1 = CreateEvent( NULL , false , true  , NULL );
	_unloadHandle2 = CreateEvent( NULL , false , true  , NULL );
}
SharedMemoryBuffer::~SharedMemoryBuffer( void )
{
	FreeArray( _buf1 ) , _bufSize1 = 0;
	FreeArray( _buf2 ) , _bufSize2 = 0;
	CloseHandle( _loadHandle1 )   , _loadHandle1   = NULL;
	CloseHandle( _loadHandle2 )   , _loadHandle2   = NULL;
	CloseHandle( _unloadHandle1 ) , _unloadHandle1 = NULL;
	CloseHandle( _unloadHandle2 ) , _unloadHandle2 = NULL;
}
SharedMemoryBuffer::FirstStream::FirstStream  ( SharedMemoryBuffer* smb ) { _smb = smb; }
SharedMemoryBuffer::SecondStream::SecondStream( SharedMemoryBuffer* smb ) { _smb = smb; }
SharedMemoryBuffer::FirstStream::~FirstStream  ( void ) { if( _smb ) delete _smb , _smb = NULL; }
SharedMemoryBuffer::SecondStream::~SecondStream( void ) {                          _smb = NULL; }
// The first stream reads on buf1 and writes on buf2
bool SharedMemoryBuffer::FirstStream::read( Pointer( byte ) buf , int len ) 
{
	if( WaitForSingleObject( _smb->_loadHandle1 , INFINITE ) != WAIT_OBJECT_0 ) fprintfId( stderr , "Wait for single failed: " ) , PrintError();
	if( len>_smb->_bufSize1 )
	{
		printf( "Uh oh 1\n" ) , fflush( stdout );
		return false;
	}
	memcpy( buf , _smb->_buf1 , len );
	SetEvent( _smb->_unloadHandle1 );
	return true;
}
bool SharedMemoryBuffer::FirstStream::write( ConstPointer( byte ) buf , int len )
{
	if( WaitForSingleObject( _smb->_unloadHandle2 , INFINITE ) != WAIT_OBJECT_0 ) fprintfId( stderr , "Wait for single failed: " ) , PrintError();
	if( len>_smb->_bufSize2 )
	{
		FreeArray( _smb->_buf2 );
		_smb->_bufSize2 = 0;
		_smb->_buf2 = AllocArray< byte >( len , 1 , "SharedMemoryBuffer::FirstStream::write (_smb->_buf2)" );
		if( !_smb->_buf2 )
		{
			printf( "Uh oh 2\n" ) , fflush( stdout );
			return false;
		}
		_smb->_bufSize2 = len;
	}
	memcpy( _smb->_buf2 , buf , len );
	SetEvent( _smb->_loadHandle2 );
	return true;
}
bool SharedMemoryBuffer::SecondStream::read( Pointer( byte ) buf , int len )
{
	if( WaitForSingleObject( _smb->_loadHandle2 , INFINITE ) != WAIT_OBJECT_0 ) fprintfId( stderr , "Wait for single failed: " ) , PrintError();
	if( len>_smb->_bufSize2 )
	{
		printf( "Uh oh 3\n" ) , fflush( stdout );
		return false;
	}
	memcpy( buf , _smb->_buf2 , len );
	SetEvent( _smb->_unloadHandle2 );
	return true;
}
bool SharedMemoryBuffer::SecondStream::write( ConstPointer( byte ) buf , int len )
{
	if( WaitForSingleObject( _smb->_unloadHandle1 , INFINITE ) != WAIT_OBJECT_0 ) fprintfId( stderr , "Wait for single failed: " ) , PrintError();
	if( len>_smb->_bufSize1 )
	{
		FreeArray( _smb->_buf1 );
		_smb->_bufSize1 = 0;
		_smb->_buf1 = AllocArray< byte >( len , 1 , "SharedMemoryBuffer::SecondStream::write (_smb->_buf1)" );
		if( !_smb->_buf1 )
		{
			printf( "Uh oh 4\n" ) , fflush( stdout );
			return false;
		}
		_smb->_bufSize1 = len;
	}
	memcpy( _smb->_buf1 , buf , len );
	SetEvent( _smb->_loadHandle1 );
	return true;
}
bool SharedMemoryBuffer::StreamPair::CreateSharedBufferPair( StreamPair& pair )
{
	SharedMemoryBuffer* smb = new SharedMemoryBuffer( );
	pair.first  = new SharedMemoryBuffer::FirstStream ( smb );
	pair.second = new SharedMemoryBuffer::SecondStream( smb );
	if( !pair.first || !pair.second )
	{
		fprintf( stderr , "Failed to create shared buffer pair\n" );
		return false;
	}
	return true;
}

SharedMemoryBuffer::StreamPair::StreamPair( void )
{
	first = NULL;
	second = NULL;
}
