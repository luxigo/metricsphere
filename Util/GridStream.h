#ifndef GRID_STREAM_INCLUDED
#define GRID_STREAM_INCLUDED
#define ASSERT_MEMORY_ACCESS 1

#define NEW_CONNECTION_BACKED_GRID 0

#include <Util/Array.h>
#include "BaseMultiStreamIO.h"
#include "Socket.h"
DWORD MySetFilePointer(HANDLE hFile,LONGLONG distanceToMove);
void MyReadFile ( HANDLE hFile , Pointer( byte ) lpBuffer , DWORD nNumberOfBytesToRead  , LPDWORD lpNumberOfBytesRead    , LPOVERLAPPED lpOverlapped );
void MyWriteFile( HANDLE hFile , Pointer( byte ) lpBuffer , DWORD nNumberOfBytesToWrite , LPDWORD lpNumberOfBytesWritten , LPOVERLAPPED lpOverlapped );
class StreamingGrid : public IOClient
{
public:
	virtual			~StreamingGrid	( void )							{;}
	virtual int		rows			( void ) const						= 0;
	virtual int		rowSize			( void ) const						= 0;
	virtual Pointer( byte ) operator[] ( int idx )						= 0;
	virtual void	advance			( void )							{;}
	virtual void	reset			( bool read , int minWindowSize )	{;}
	virtual void	unset			( void )							{;}
	virtual bool	SeparateColors	( void ) const						{ return true; }
	virtual void	SetServer		( MultiStreamIOServer* server )		{;}
	virtual int		Service			( void )							{ return IOClient::NONE; }
	virtual bool	HasAlpha		( void ) const						{ return false; }
};
class NULLStreamingGrid : public StreamingGrid
{
	Pointer( byte ) _data;
	int _rows , _rowSize;
public:
	NULLStreamingGrid	( int rowSize , int rows );
	~NULLStreamingGrid	( void );
	int		rows		( void ) const;
	int		rowSize		( void ) const;
	Pointer( byte ) operator[] ( int idx );
};

class BufferedStreamingGrid : public StreamingGrid
{
	bool read;
	Pointer( byte ) data;
	int current , win;
	StreamingGrid* sg;
public:
	BufferedStreamingGrid	(StreamingGrid* sg);
	~BufferedStreamingGrid	(void);
	int		rows			(void) const;
	int		rowSize			(void) const;
	Pointer( byte )			operator[] ( int idx );
	void	advance			(void);
	void	reset			(bool read,int minWindowSize);
};
#if NEW_CONNECTION_BACKED_GRID
template<int Channels , class Real >
class MultiDataStreamBackedGrid : public StreamingGrid
{
	bool	_read , _readComplete , _blockingSend , _separateColors;
	Real	*_data , **_subData;
	int		_rows , *_rowSizes , _rowSize , _current , _streamCount;
	DataStream**	_streams;
public:
	MultiDataStreamBackedGrid	( DataStream** streams , int *rowSizes , int streamCount , int rows , bool separateColors);
	~MultiDataStreamBackedGrid	( void );
	int		rows		( void ) const;
	int		rowSize		( void ) const;
	Pointer( ( byte ) ) operator[] ( int idx );
	void	advance		( void );
	void	reset		( bool read , int minWindowSize );
	bool	SeparateColors( void ) const { return _separateColors; }
};

class DataStreamBackedGrid : public StreamingGrid
{
	bool	_read , _readComplete , _bs , _separateColors;
	void*	_data;
	int		_r , _rs , _c;
	DataStream*	_stream;
public:
	DataStreamBackedGrid	( DataStream* stream , int rowSize , int rows , bool separateColors );
	~DataStreamBackedGrid	( void );
	int		rows		( void ) const;
	int		rowSize		( void ) const;
	Pointer( ( byte ) ) operator[] ( int idx );
	void	advance		( void );
	void	reset		( bool read , int minWindowSize );
	bool	SeparateColors( void ) const { return _separateColors; }
};
#else // !NEW_CONNECTION_BACKED_GRID
template<int Channels , class Real >
class MultiSocketBackedGrid : public StreamingGrid
{
	bool	_read , _readComplete , _blockingSend , _separateColors;
	Pointer( Real ) _data;
	Pointer( Pointer( Real ) ) _subData;
	int		_rows , *_rowSizes , _rowSize , _current , _sockCount;
	SOCKET*	_socks;
public:
	MultiSocketBackedGrid	( SOCKET* socks , int *rowSizes , int sockCount , int rows , bool blockingSend , bool separateColors);
	~MultiSocketBackedGrid	( void );
	int		rows		( void ) const;
	int		rowSize		( void ) const;
	Pointer( byte ) operator[] ( int idx );
	void	advance		( void );
	void	reset		( bool read , int minWindowSize );
	bool	SeparateColors( void ) const { return _separateColors; }
};

class SocketBackedGrid : public StreamingGrid
{
	bool	_read , _readComplete , _bs , _separateColors;
	Pointer( byte ) _data;
	int		_r , _rs , _c;
	SOCKET	_sock;
public:
	SocketBackedGrid	( SOCKET sock , int rowSize , int rows , bool blockingSend , bool separateColors );
	~SocketBackedGrid	( void );
	int		rows		( void ) const;
	int		rowSize		( void ) const;
	Pointer ( byte ) operator[] ( int idx );
	void	advance		( void );
	void	reset		( bool read , int minWindowSize );
	bool	SeparateColors( void ) const { return _separateColors; }
};
#endif // NEW_CONNECTION_BACKED_GRID
class SharedMemoryGrid : public StreamingGrid
{
	bool    _read , _dataReady;
	Pointer( byte ) _data;
	int     _r , _rs , _c;
	HANDLE *_startHandles , *_endHandles;
	int     _hCount;
public:
	SharedMemoryGrid	( HANDLE* startHandles , HANDLE* endHandles , int handleCount , Pointer( byte ) row , int rowSize , int rows );
	~SharedMemoryGrid	( void );
	int		rows		( void ) const;
	int		rowSize		( void ) const;
	Pointer( byte ) operator[] ( int idx );
	void	advance		( void );
	void	reset		( bool read , int minWindowSize );
	bool	SeparateColors( void ) const { return false; }
};
class MemoryBackedGrid : public StreamingGrid
{
	Pointer( byte ) _data;
	int		_r,_rs;
	bool	_del;
	bool	_separateColors;
public:
	MemoryBackedGrid	( Pointer( byte ) data , int rowSize , int rows , bool separateColors );
	MemoryBackedGrid	( int rowSize , int rows , bool separateColors );
	~MemoryBackedGrid	( void );
	int		rows		( void ) const;
	int		rowSize		( void ) const;
	Pointer( byte ) operator[] ( int idx );
	bool	SeparateColors( void ) const { return _separateColors; }
};
class FileBackedGrid : public StreamingGrid
{
	long long _rSize,_wSize;
	char*	_name;
	FILE*	_fp;
	Pointer( byte ) _data;
	int		_r,_rs,_c,_w;
	bool	_read;
	void	_readNext(void);
	void	_writeNext(void);
public:
	FileBackedGrid		(int rowSize,int rows,FILE* fp);
	FileBackedGrid		(int rowSize,int rows);
	~FileBackedGrid		(void);
	int		rows		(void) const;
	int		rowSize		(void) const;
	Pointer( byte ) operator[] ( int idx );
	void	advance		(void);
	void	reset		(bool read,int minWindowSize);
};
class CompressedFileBackedGrid : public StreamingGrid
{
	int		_sSize;
	long long _rSize,_wSize;
	char*	_name;
	FILE*	_fp;
	Pointer( byte ) _data;
	Pointer( byte ) _scratch;
	int		_r,_rs,_c,_w;
	bool	_read;
	void	_readNext(void);
	void	_writeNext(void);
public:
	CompressedFileBackedGrid	(int rowSize,int rows,FILE* fp);
	CompressedFileBackedGrid	(int rowSize,int rows);
	~CompressedFileBackedGrid	(void);
	int		rows				(void) const;
	int		rowSize				(void) const;
	Pointer( byte ) operator[] ( int idx );
	void	advance				(void);
	void	reset				(bool read,int minWindowSize);
};

#include "GridStream.inl"
#endif // GRID_STREAM_INCLUDED