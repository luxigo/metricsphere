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

#ifndef SOCKET_INCLUDED
#define SOCKET_INCLUDED

#include <stdarg.h>
#include <WinSock.h>
#include <stdio.h>
#include "Util/Time.h"
#include "Util/Array.h"
#include <hash_map>
using stdext::hash_map;


#define DEBUG_DEADLOCK 0
#define STORE_CONNECTION_TABLE 1

void printfId(const char* format,...);
void fprintfId(FILE* fp , const char* format,...);

class ConnectionData
{
public:
	in_addr localAddr , peerAddr;
	int localPort , peerPort;
};

class MyWinSock
{
	static int _wsaCount;
	static WSADATA _wsaData;
	static CRITICAL_SECTION stdinLock , stderrLock , systemLock;
public:
	static void Load( void );
	static void UnLoad( void );
#if STORE_CONNECTION_TABLE
	static hash_map< int , ConnectionData > ConnectionTable;
#endif // STORE_CONNECTION_TABLE
	class StdinLock
	{
	public:
		StdinLock( void ) { if( _wsaCount>0 ) EnterCriticalSection( &stdinLock ); }
		~StdinLock( void ) { if( _wsaCount>0 ) LeaveCriticalSection( &stdinLock ); }
	};
	class StderrLock
	{
	public:
		StderrLock( void ) { if( _wsaCount>0 ) EnterCriticalSection( &stderrLock ); }
		~StderrLock( void ) { if( _wsaCount>0 ) LeaveCriticalSection( &stderrLock ); }
	};
	class SystemLock
	{
	public:
		SystemLock( void ) { if( _wsaCount>0 ) EnterCriticalSection( &systemLock ); }
		~SystemLock( void ) { if( _wsaCount>0 ) LeaveCriticalSection( &systemLock ); }
	};
	static void printfId( const char* format , ... );
	static void fprintfId( FILE* fp , const char* format , ... );
};

bool MyWaitForSingleObject( HANDLE hHandle , DWORD milliseconds , const char* deadlockMessage , ... );
bool MyWaitForMultipleObjects( DWORD nCount , const HANDLE* lpHandles , DWORD milliseconds , const char* deadlockMessage , ... );

static long long packetsSent = 0 , packetsReceived = 0;
//static const char *LastSocketError( void );
const char *LastSocketError( void );
int inet_aton(const char *cp, struct in_addr *inp);

template<class C>	bool ReceiveOnSocket	( SOCKET& s , Pointer( C ) data , int dataSize );
template<class C>	void ReceiveOnSocket	( SOCKET& s , Pointer( C ) data , int dataSize , const char* errorMessage , ... );
template<class C>	bool ReceiveOnSocket	( SOCKET& s , Pointer( C ) data , int dataSize , bool blockingSend );
template<class C>	void ReceiveOnSocket	( SOCKET& s , Pointer( C ) data , int dataSize , bool blockingSend , const char* errorMessage , ... );

template<class C>	bool SendOnSocket		( SOCKET& s , ConstPointer( C ) data , int dataSize );
template<class C>	void SendOnSocket		( SOCKET& s , ConstPointer( C ) data , int dataSize , const char* errorMessage , ... );
template<class C>	bool SendOnSocket		( SOCKET& s , ConstPointer( C ) data , int dataSize , bool blockingSend );
template<class C>	void SendOnSocket		( SOCKET& s , ConstPointer( C ) data , int dataSize , bool blockingSend , const char* errorMessage , ... );

bool StartReceiveOnSocket( SOCKET& s , bool blockingSend );
bool EndSendOnSocket( SOCKET& s , bool blockingSend );
void StartReceiveOnSocket	( SOCKET& s , bool blockingSend , const char* errorMessage , ... );
void EndSendOnSocket		( SOCKET& s , bool blockingSend , const char* errorMessage , ... );
SOCKET GetListenSocket( int& port );
SOCKET AcceptSocket( SOCKET listen );
SOCKET GetConnectSocket( const char* address , int port , int ms=5 , bool progress=false );
SOCKET GetConnectSocket( in_addr , int port , int ms=5 , bool progress=false );
void CloseSocket( SOCKET& s );
int         GetLocalSocketPort   ( SOCKET& s );
const char* GetLocalSocketAddress( SOCKET& s );
int         GetPeerSocketPort   ( SOCKET& s );
const char* GetPeerSocketAddress( SOCKET& s );
void GetHostAddress( char* address , char* prefix = NULL );
void PrintHostAddress( void );


class DataStream
{
public:
	virtual ~DataStream( void ){ }
	virtual bool write( ConstPointer( byte ) buf , int len ) = 0;
	virtual bool read ( Pointer(      byte ) buf , int len ) = 0;

	static DataStream* GetDataStream( SOCKET sock , bool master , bool cleanUp = true );
};

class DataStreamConstructor
{
	SOCKET _sock;
	bool _master;
	bool _cleanUp;
	DataStream* _stream;
	in_addr _myAddr;
	int _myPID;

public:
	static const int STEPS = 3;
	DataStreamConstructor( void );
	void init( SOCKET sock , bool master , bool cleanUp = true );
	void doStep( int sNum );
	DataStream* getDataStream( void );
};
class Socket : public DataStream
{
	SOCKET _sock;
	bool _mySocket;
public:
	Socket( SOCKET sock );
	Socket( in_addr address , int port , int ms=5 , bool progress=false );
	Socket( const char* address , int port , int ms=5 , bool progress=false );
	~Socket( void );
	bool write( ConstPointer( byte ) buf , int len );
	bool read ( Pointer(      byte ) buf , int len );
};

class SharedMemoryBuffer
{
protected:
	Pointer( byte ) _buf1;
	Pointer( byte ) _buf2;
	int _bufSize1 , _bufSize2;
	HANDLE _loadHandle1 , _unloadHandle1 , _loadHandle2 , _unloadHandle2;
	SharedMemoryBuffer( void );
	~SharedMemoryBuffer( void );
public:
	class FirstStream : public DataStream
	{
		friend class SharedMemoryBuffer;
	protected:
		SharedMemoryBuffer* _smb;
		FirstStream( SharedMemoryBuffer* smb );
	public:
		~FirstStream( void );
		bool write( ConstPointer( byte ) buf , int len );
		bool read ( Pointer(      byte ) buf , int len );
	};
	class SecondStream : public DataStream
	{
		friend class SharedMemoryBuffer;
	protected:
		SharedMemoryBuffer* _smb;
		SecondStream( SharedMemoryBuffer* smb );
	public:
		~SecondStream( void );
		bool write( ConstPointer( byte ) buf , int len );
		bool read (      Pointer( byte ) buf , int len );
	};
	class StreamPair
	{
	public:
		StreamPair( void );
		FirstStream* first;
		SecondStream* second;
		static bool CreateSharedBufferPair( StreamPair& pair );
	};
};

template< class C >
DWORD WINAPI StartThread( LPVOID lpParam )
{
	( (C*) lpParam )->Run();
	return 0;
}

template< class C >
HANDLE SpawnThread( C* c , int threadPriority=THREAD_PRIORITY_NORMAL )
{
	HANDLE threadHandle = NULL;
	DWORD threadID;
	threadHandle = CreateThread
		( 
		NULL,				// default security attributes
		0,					// use default stack size  
		StartThread< C >,	// thread function 
		c,					// argument to thread function 
		0,					// use default creation flags 
		&threadID			// returns the thread identifier
		);
	if( !threadHandle )
	{
		fprintf( stderr , "CreateThread failed: " );
		PrintError();
	}
	SetThreadPriority( threadHandle , threadPriority );
	return threadHandle;
}

#include "Socket.inl"
#endif // SOCKET_INCLUDED
