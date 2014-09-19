#ifndef MEMORY_STREAM_INCLUDED
#define MEMORY_STREAM_INCLUDED

#define USE_VIRTUAL 0
#define NEW_DATA_STREAM_CODE 1

#if USE_VIRTUAL
#define USE_THREADS 1
#endif // USE_VIRTUAL

#include <windows.h>
#if USE_THREADS
#include <tchar.h>
#include <strsafe.h>
#endif // USE_THREADS

template<class Data>
class DataStream
{
public:
	virtual ~DataStream(void){;};
	virtual Data& operator [] (size_t idx) = 0;
	virtual size_t size(void) const = 0;
	virtual void reset(void) = 0;
};
template<class Data>
class MemoryStream : public DataStream<Data>
{
	bool deleteData;
//	Data *data,*_data;
	Pointer( Data ) data;
	size_t _size;
public:
	MemoryStream(size_t size,int byteAlignment=1,bool zeroMemory=true);
//	MemoryStream(Data* data,size_t size);
	MemoryStream( Pointer( Dat ) data , size_t size );
	~MemoryStream(void);
	Data& operator[] (size_t idx);
	size_t size(void) const;
	void reset(void);
};
template<class Data>
class FileMappedStream : public DataStream<Data>
{
	HANDLE hFile;
	HANDLE data;
	size_t _size,_blockSize;
	Data *baseData,*_baseData,*dataPtr;
#if USE_THREADS
	CRITICAL_SECTION lock;
	Data *frontPtr,*backPtr;
	size_t maxReadAhead;

	HANDLE pfThread,wbThread;
	static DWORD WINAPI PrefetchThread(LPVOID lpParam);
	static DWORD WINAPI WritebackThread(LPVOID lpParam);
#endif // USE_THREADS
public:
	FileMappedStream(const char* fileName,size_t blockSize,int byteAlignment=1);
	FileMappedStream(size_t size,size_t blockSize,int byteAlignment=1);
	~FileMappedStream(void);
	Data& operator[] (size_t idx);
	size_t size(void) const;
	void reset(void);
};
#include "DataStream.inl"
#endif // MEMORY_STREAM_INCLUDED