#include <windows.h>
#include <atlstr.h>
#include "BaseMultiStreamIO.h"
#include "Socket.h"

//////////////
// IOClient //
//////////////
const int IOClient::BYTES_PER_SECTOR = 1<<9;					// 512B Sector sizes
//const int IOClient::IO_BLOCK_SIZE = BYTES_PER_SECTOR<<13;		// 4MB IO Chunks
const int IOClient::IO_BLOCK_SIZE = BYTES_PER_SECTOR<<12;		// 2MB IO Chunks
long long IOClient::ReadBytes  = 0;
long long IOClient::WriteBytes = 0;

IOClient::IOClient( void )
{
	InitializeCriticalSection( &lock );
	server = NULL;
}
IOClient::~IOClient( void )
{
	DeleteCriticalSection( &lock );
}
void IOClient::SetServer( class MultiStreamIOServer* server )
{
	this->server = server;
	if( server ) server->AddClient( this );
}
/////////////////////////
// MultiStreamIOServer //
/////////////////////////
MultiStreamIOServer::MultiStreamIOServer(void)
{
	ioThread = NULL;
	InitializeCriticalSection( &pendingLock );
	InitializeCriticalSection( &clientLock );

	pendingClient = NULL;
	DWORD ioThreadID;
#if 1
	ioThread = CreateThread( 
		NULL,			// default security attributes
		0,				// use default stack size  
		IOThread,		// thread function 
		this,			// argument to thread function 
		0,				// use default creation flags 
		&ioThreadID);	// returns the thread identifier
	if(!ioThread)
	{
		fprintf(stderr,"Failed to create I/O thread\n");
		exit(0);
	}
#endif
}
MultiStreamIOServer::~MultiStreamIOServer(void)
{
	WaitOnIO();
	DeleteCriticalSection( &pendingLock );
	DeleteCriticalSection( &clientLock );
}
bool MultiStreamIOServer::SetPending( IOClient* client )
{
	EnterCriticalSection( &pendingLock );
	if( !pendingClient && client )
	{
		pendingClient = client;
		LeaveCriticalSection( &pendingLock );
		return true;
	}
	else
	{
		LeaveCriticalSection( &pendingLock );
		return false;
	}
}
void MultiStreamIOServer::AddClient( IOClient* client )
{
	EnterCriticalSection( &clientLock );
	clients.push_back( client );
	LeaveCriticalSection( &clientLock );
}
void MultiStreamIOServer::WaitOnIO(void)
{
	if( ioThread )
	{
		while( 1 )
		{
			// As long as the server has work to do, it cannot be terminated.
			EnterCriticalSection( &clientLock );
			if( !clients.size() )
			{
				if( !TerminateThread( ioThread , 0 ) ) fprintf( stderr , "Failed to terminate MultiStreamIOServer thread\n" );
				LeaveCriticalSection( &clientLock );
				break;
			}
			LeaveCriticalSection( &clientLock );
		}

	}
}
int MultiStreamIOServer::clientNum( void )
{
	EnterCriticalSection( &clientLock );
	int sz = clients.size();
	LeaveCriticalSection( &clientLock );
	return sz;
}
DWORD WINAPI MultiStreamIOServer::IOThread( LPVOID lpParam )
{
	MultiStreamIOServer* server = (MultiStreamIOServer*)lpParam;
	std::vector< IOClient* >& clients = server->clients;
	int idx = 0;
	while( 1 )
	{
		EnterCriticalSection( &server->pendingLock );
		if( server->pendingClient )
		{
			if( server->pendingClient->Service() == IOClient::NONE ) Sleep(0);
			else server->pendingClient = NULL;
			LeaveCriticalSection( &server->pendingLock );
		}
		else
		{
			LeaveCriticalSection( &server->pendingLock );
			EnterCriticalSection( &server->clientLock );
			bool ioDone = false;
			for( int i=0 ; i<clients.size() && !ioDone ; i++ )
			{
				idx = (idx+1)%clients.size();
				switch( clients[idx]->Service() )
				{
				case IOClient::COMPLETE:
					clients[idx]->SetServer( NULL );
					clients[idx] = clients[clients.size()-1];
					clients.pop_back();
				case IOClient::SUCCESS:
					ioDone = true;
					break;
				}
			}
			LeaveCriticalSection( &server->clientLock );
			if( !ioDone ) Sleep( 1 );
			else Sleep( 0 );
		}
	}
	return 0;
}
