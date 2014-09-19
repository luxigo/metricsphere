#ifndef MULTI_STREAM_IO
#define MULTI_STREAM_IO
#include <vector>
#include "GridStream.h"

#define NEW_MULTI_STREAM_CODE 1

// BADNESS!!! Apparently bad things can happen if the bufferMultiplier is set to one. Possibly related to the fact that
// rows can overflow buffers? 

class IOClientStream
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
	CRITICAL_SECTION lock;						// A locking devices so that the server doesn't step on the client's toes
	IOClientStream ( void );
	~IOClientStream( void );
	virtual int	Update( void ) = 0;				// The task the server asks the client to do
};
class IOClient
{
public:
	friend class MultiStreamIOServer;
	class MultiStreamIOServer* server;
	IOClientStream* stream;
	int clientIndex;
	IOClient( void );
	void SetServer( class MultiStreamIOServer* server );
};

// When only the  maximum size of the requested read/write is known in advance.
class VariableIOClientStream : public IOClientStream , public IOClient
{
	int		_blockSize;							// The actual size of a read/write operation
	void*	_data;								// The buffer
	int		_bufferMultiplier;					// The number of blocks to store in memory at any given time (must be at least 2)
	FILE*	_fp;								// Where to read/write data
	int		_current;							// The index of the buffer from which data can be requested
	int		_frontAndBack;						// For reading, [ _current , _frontAndBack ] is in the buffer.
												// For writing, [ _frontAndBack , _current ]  is in the buffer
	bool	_read;								// Are we reading data?
	bool	_endIO;								// Have we performed all the requisite I/O?
	long long _pseudoHead;						// Where are we believed to be in the file
	long long _head;							// Where are we actually in the file
public:
	VariableIOClientStream( FILE* fp , int maxIOSize , int bufferMultiplier , bool read );
	~VariableIOClientStream( void );
	int		Update	( void );
	size_t  read (       void* buffer , size_t sz );
	size_t  write( const void* buffer , size_t sz );
};

class FixedIOClientStream : public IOClientStream
{
public:
	int*				off;
	int					blockSize;
	void*				data;
	int					r , rs , win , b;
	HANDLE				hFile;
	bool				read;	// write;
	int					current,front,back;

	FixedIOClientStream ( void );
	~FixedIOClientStream( void );
	void	Init	( HANDLE hFile , int rowSize , int rows , int bufferMultiplier );
	void	Reset	( bool read , int minWindowSize );
	void	Unset	( void );
	int		Update	( void );
	bool	Advance	( void );
	void*	operator[] ( int idx );
};

class MultiStreamIOServer
{
	HANDLE ioThread;
	static DWORD WINAPI IOThread( LPVOID lpParam );
#if NEW_MULTI_STREAM_CODE
	std::vector< IOClientStream* > streams;
	std::vector< IOClient* > clients;
#else // !NEW_MULTI_STREAM_CODE
	std::vector<FixedIOClientStream*> streams;
	std::vector< class MultiStreamIOClient* > clients;
#endif // NEW_MULTI_STREAM_CODE
	void WaitOnIO( void );
public:
	CRITICAL_SECTION	lock;
	int pendingStream;
	MultiStreamIOServer ( void );
	~MultiStreamIOServer( void );

#if NEW_MULTI_STREAM_CODE
	virtual int AddClient	( IOClient* client );
#else // !NEW_MULTI_STREAM_CODE
	virtual int AddClient	( MultiStreamIOClient* client );
#endif // NEW_MULTI_STREAM_CODE
	virtual void StartIO	( void );
	virtual void Reset		( void );
};

#if NEW_MULTI_STREAM_CODE
class MultiStreamIOClient : public StreamingGrid , public IOClient
{
#else // !NEW_MULTI_STREAM_CODE
class MultiStreamIOClient : public StreamingGrid
{
	friend class MultiStreamIOServer;
	int clientIndex;
	FixedIOClientStream stream;
	MultiStreamIOServer* server;
#endif // NEW_MULTI_STREAM_CODE
	HANDLE hFile;
public:
	MultiStreamIOClient	( const char* fileName , int rs , int r , int bufferMultiplier , bool writeOnly=false );
	MultiStreamIOClient	( int rs , int r , int bufferMultiplier , const char* prefix , bool deleteOnClose );
	MultiStreamIOClient	( int rs , int r , int bufferMultiplier , const char* dir , const char* prefix , bool deleteOnClose);
	~MultiStreamIOClient(void);
	void	finalize	(void);
#if !NEW_MULTI_STREAM_CODE
	void	SetServer	(MultiStreamIOServer* manager);
#endif // NEW_MULTI_STREAM_CODE

	int		rows		( void ) const;
	int		rowSize		( void ) const;
	void*	operator[]	( int idx );
	void	advance		( void );
	void	reset		( bool read , int minWindowSize );
	void	unset		( void );
};

#include "MultiStreamIO.inl"
#endif // MULTI_STREAM_IO