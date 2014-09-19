#ifndef BASE_MULTI_STREAM_IO
#define BASE_MULTI_STREAM_IO
#include <vector>
#include <windows.h>

#define NEW_MULTI_STREAM_CODE 1
#define NEW_JPEG_IO 1
#define NEW_PNG_IO 1

// BADNESS!!! Apparently bad things can happen if the bufferMultiplier is set to one. Possibly related to the fact that
// rows can overflow buffers? 

class IOClient
{
public:
	static const int BYTES_PER_SECTOR;			// The atomic size of a read/write operation (for some IO applications)
	static const int IO_BLOCK_SIZE;				// A nice size for reads
	static long long ReadBytes , WriteBytes;	// Total bytes read/written
	enum
	{
		NONE,		// If there was no I/O that the client could do
		SUCCESS,	// If the client succeeded in performing the I/O
		COMPLETE	// If the client won't need to do any more I/O
	};

	friend class MultiStreamIOServer;
	CRITICAL_SECTION lock;				// A locking devices so that the server doesn't step on the client's toes
	class MultiStreamIOServer* server;	// The server responsible for processing the job requests (this should only be set by the server)!

	IOClient ( void );
	~IOClient( void );

	void SetServer( class MultiStreamIOServer* server );
	virtual int	Service( void ) = 0;		// The task the server asks the client to do
};

class MultiStreamIOServer
{
	HANDLE ioThread;
	static DWORD WINAPI IOThread( LPVOID lpParam );
	std::vector< IOClient* > clients;
	void WaitOnIO( void );
	CRITICAL_SECTION pendingLock , clientLock;
	IOClient* pendingClient;
public:
	MultiStreamIOServer ( void );
	~MultiStreamIOServer( void );

	int clientNum( void );
	virtual bool SetPending ( IOClient* client );
	virtual void AddClient	( IOClient* client );
};
#endif // BASE_MULTI_STREAM_IO