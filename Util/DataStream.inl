//////////////////
// MemoryStream //
//////////////////
template<class Data>
MemoryStream<Data>::MemoryStream(size_t size,int byteAlignment,bool zeroMemory)
{
	data = AllocArray< Data >( size , byteAlignment , "MemoryStream< Data >::data" );
	if( zeroMemory ) memset( data , 0 , sizeof( Data ) * size );
	_size=size;
	deleteData=true;
}
template<class Data>
MemoryStream<Data>::MemoryStream( Pointer( Data ) data,size_t size)
{
	this->data = data;
	_size=size;
	deleteData=false;
}
template<class Data>
MemoryStream<Data>::~MemoryStream(void)
{
	if(deleteData && data)	FreeArray( data );
	data = NullPointer< NULL >( );
	_size=0;
}
template<class Data>
inline Data& MemoryStream<Data>::operator [] (size_t idx)
{
	return *(data+idx);
}
template<class Data>
inline size_t MemoryStream<Data>::size(void) const	{ return _size; }
template<class Data>
inline void MemoryStream<Data>::reset(void){ ; }


//////////////////////
// FileMappedStream //
//////////////////////
#include <cstringt.h>
#include <atlstr.h>

#if USE_THREADS
template<class Data>
DWORD WINAPI FileMappedStream<Data>::PrefetchThread(LPVOID lpParam)
{
	FileMappedStream<Data>* stream = (FileMappedStream<Data>*)lpParam;
	while(1)
	{
		// Try and grab the lock
		EnterCriticalSection(&stream->lock);
		if(	
			// If there is additional memory that can be loaded
			stream->frontPtr<stream->baseData+stream->_size	
			&&
			// If we haven't gotten too far ahead of ourselves
			stream->frontPtr-stream->backPtr<stream->maxReadAhead
			)
		{
			VirtualLock(stream->frontPtr,stream->_blockSize);
			stream->frontPtr+=stream->_blockSize;
			LeaveCriticalSection(&stream->lock);
		}
		else
		{
			LeaveCriticalSection(&stream->lock);
//			SuspendThread(stream->pfThread);
			Sleep(1);
		}
	}
	return 0;
}
template<class Data>
DWORD WINAPI FileMappedStream<Data>::WritebackThread(LPVOID lpParam)
{
	FileMappedStream<Data>* stream = (FileMappedStream<Data>*)lpParam;
	while(1)
	{
		// Try and grab the lock
		EnterCriticalSection(&stream->lock);
		if(
			// If there is data to release
			stream->dataPtr>stream->backPtr	
			)
		{
			VirtualUnlock(stream->backPtr,stream->_blockSize);
			stream->backPtr+=stream->_blockSize;
			LeaveCriticalSection(&stream->lock);
		}
		else
		{
			LeaveCriticalSection(&stream->lock);
//			SuspendThread(stream->wbThread);
			Sleep(1);
		}
	}
	return 0;
}
#endif // USE_THREADS

template<class Data>
FileMappedStream<Data>::FileMappedStream(const char* fileName,size_t blockSize,int byteAlignment)
{
	CString s(fileName);
	hFile = CreateFile(s,FILE_ALL_ACCESS,0,NULL,OPEN_ALWAYS,FILE_ATTRIBUTE_NORMAL, NULL);
	data  = CreateFileMapping(hFile, NULL, PAGE_READWRITE | SEC_COMMIT, 0, 0, NULL);
	_baseData = (Data*)MapViewOfFile(data,FILE_MAP_ALL_ACCESS,0,0,0);
	baseData  = (Data*)(((size_t)(_baseData)+byteAlignment-1) & ~(byteAlignment-1));
	_size=0;
	_blockSize=blockSize;
#if USE_THREADS
	maxReadAhead=blockSize*10;
	frontPtr=backPtr=dataPtr;
	InitializeCriticalSection(&lock);
	DWORD pfThreadID,wbThreadID;
	pfThread=CreateThread( 
		NULL,				// default security attributes
		0,					// use default stack size  
		PrefetchThread,		// thread function 
		this,				// argument to thread function 
		0,					// use default creation flags 
		&pfThreadID);		// returns the thread identifier
	if(!pfThread)
	{
		fprintf(stderr,"Failed to create prefetch thread\n");
		exit(0);
	}
	wbThread=CreateThread( 
		NULL,				// default security attributes
		0,					// use default stack size  
		WritebackThread,	// thread function 
		this,				// argument to thread function 
		0,					// use default creation flags 
		&wbThreadID);		// returns the thread identifier
	if(!wbThread)
	{
		fprintf(stderr,"Failed to create writeback thread\n");
		exit(0);
	}
#else // !USE_THREADS
#if USE_VIRTUAL
	VirtualLock(dataPtr,blockSize);
#endif // USE_VIRTUAL
#endif // USE_THREADS
}

template<class Data>
FileMappedStream<Data>::FileMappedStream(size_t size,size_t blockSize,int byteAlignment)
{
	CString s(tmpnam(NULL));
	hFile = CreateFile(s,FILE_ALL_ACCESS,0,NULL,CREATE_ALWAYS,FILE_ATTRIBUTE_NORMAL | FILE_FLAG_DELETE_ON_CLOSE | FILE_ATTRIBUTE_TEMPORARY,NULL);
	data  = CreateFileMapping(hFile, NULL, PAGE_READWRITE | SEC_COMMIT, 0, size*sizeof(Data)+byteAlignment-1, NULL);
	_baseData = (Data*)MapViewOfFile(data,FILE_MAP_ALL_ACCESS,0,0,size*sizeof(Data)+byteAlignment-1);
	baseData  = (Data*)(((size_t)(_baseData)+byteAlignment-1) & ~(byteAlignment-1));
	_size=size;
	_blockSize=blockSize;
#if USE_THREADS
	maxReadAhead=blockSize*10;
	frontPtr=backPtr=dataPtr;
	InitializeCriticalSection(&lock);
	DWORD pfThreadID,wbThreadID;
	pfThread=CreateThread( 
		NULL,				// default security attributes
		0,					// use default stack size  
		PrefetchThread,		// thread function 
		this,				// argument to thread function 
		0,					// use default creation flags 
		&pfThreadID);		// returns the thread identifier
	if(!pfThread)
	{
		fprintf(stderr,"Failed to create prefetch thread\n");
		exit(0);
	}
	wbThread=CreateThread( 
		NULL,				// default security attributes
		0,					// use default stack size  
		WritebackThread,	// thread function 
		this,				// argument to thread function 
		0,					// use default creation flags 
		&wbThreadID);		// returns the thread identifier
	if(!wbThread)
	{
		fprintf(stderr,"Failed to create writeback thread\n");
		exit(0);
	}
#else // !USE_THREADS
#if USE_VIRTUAL
	VirtualLock(dataPtr,blockSize);
#endif // USE_VIRTUAL
#endif // USE_THREADS
}
template<class Data>
FileMappedStream<Data>::~FileMappedStream(void)
{
#if USE_THREADS
	VirtualUnlock(dataPtr,frontPtr-dataPtr);
#else // !USE_THREADS
#if USE_VIRTUAL
	VirtualUnlock(dataPtr,_blockSize);
#endif // USE_VIRTUAL
#endif // USE_THREADS
	UnmapViewOfFile(_baseData);
	CloseHandle(hFile);
#if USE_THREADS
	TerminateThread(pfThread,0);
	TerminateThread(wbThread,0);
	CloseHandle(pfThread);
	CloseHandle(wbThread);
	DeleteCriticalSection(&lock);
#endif // USE_THREADS
}
template<class Data>
void FileMappedStream<Data>::reset(void)
{
#if USE_THREADS
	EnterCriticalSection(&lock);
	VirtualUnlock(backPtr,frontPtr-backPtr);
	dataPtr = frontPtr = backPtr = baseData;
	LeaveCriticalSection(&lock);
#else // !USE_THREADS
#if USE_VIRTUAL
	VirtualUnlock(dataPtr,_blockSize);
#endif // USE_VIRTUAL
	dataPtr = baseData;
#if USE_VIRTUAL
	VirtualLock(dataPtr,_blockSize);
#endif // USE_VIRTUAL
#endif // USE_THREADS
}
template<class Data>
size_t FileMappedStream<Data>::size(void) const
{
	return _size;
}
template<class Data>
Data& FileMappedStream<Data>::operator [] (size_t idx)
{
#if USE_THREADS
	EnterCriticalSection(&lock);

	dataPtr=baseData+idx;
	// If we are processing faster than we are writing...
	if(frontPtr-backPtr>=maxReadAhead && backPtr<dataPtr)
	{
		VirtualUnlock(backPtr,_blockSize);
		backPtr+=_blockSize;
	}
	// If we are processing faster than we are reading...
	if(dataPtr==frontPtr)
	{
		VirtualLock(dataPtr,_blockSize);
		frontPtr=dataPtr+_blockSize;
		LeaveCriticalSection(&lock);
	}
	else	LeaveCriticalSection(&lock);
#else // !USE_THREADS
#if USE_VIRTUAL
	VirtualUnlock(dataPtr,_blockSize);
#endif // USE_VIRTUAL
	dataPtr=baseData+idx;
#if USE_VIRTUAL
	VirtualLock(dataPtr,_blockSize);
#endif // USE_VIRTUAL
#endif // USE_THREADS
//	ResumeThread(pfThread);
//	ResumeThread(wbThread);
	return *dataPtr;
}

