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

#include <stdio.h>
#include <string.h>
#include <sys/timeb.h>
#ifndef WIN32
#include <sys/time.h>
#else
#include <atlstr.h>
#include <psapi.h>
#endif // WIN32

double Time(void){
#ifdef WIN32
	static bool first=true;
	static __int64 lfreq;

	if (first) {
		first=false;
		LARGE_INTEGER l;
		QueryPerformanceFrequency(&l);
		lfreq=l.QuadPart;
	}

	LARGE_INTEGER l;
	QueryPerformanceCounter(&l);
	return double(l.QuadPart)/double(lfreq);
#else // WIN32
	struct timeval t;
	gettimeofday(&t,NULL);
	return t.tv_sec+(double)t.tv_usec/1000000;
#endif // WIN32
}
#ifdef WIN32
void PrintError(DWORD err)
{
	char error[1024];
	int sz=FormatMessage(FORMAT_MESSAGE_FROM_SYSTEM, NULL, err, 0,(LPTSTR)error, 1024, NULL);
	if(sz>0)
	{
		fprintf(stderr,"\t");
		for(int i=0;i<sz;i++)	printf("%c",error[2*i]);
	}
}
void PrintError(void)	{PrintError(GetLastError());}
void WorkingSetInfo(size_t& current,size_t& peak)
{
    PROCESS_MEMORY_COUNTERS pmc;
	if ( GetProcessMemoryInfo( GetCurrentProcess(), &pmc, sizeof(pmc)) )
	{
#if 1
		current=pmc.WorkingSetSize;
		peak=pmc.PeakWorkingSetSize;
#else
		printf( "\tPageFaultCount: 0x%08X\n", pmc.PageFaultCount ) , fprintf( stdout );
		printf( "\tPeakWorkingSetSize: 0x%08X\n",pmc.PeakWorkingSetSize ) , fprintf( stdout );
		printf( "\tWorkingSetSize: 0x%08X\n", pmc.WorkingSetSize ) , fprintf( stdout );
		printf( "\tQuotaPeakPagedPoolUsage: 0x%08X\n", pmc.QuotaPeakPagedPoolUsage ) , fprintf( stdout );
		printf( "\tQuotaPagedPoolUsage: 0x%08X\n", pmc.QuotaPagedPoolUsage ) , fprintf( stdout );
		printf( "\tQuotaPeakNonPagedPoolUsage: 0x%08X\n", pmc.QuotaPeakNonPagedPoolUsage ) , fprintf( stdout );
		printf( "\tQuotaNonPagedPoolUsage: 0x%08X\n", pmc.QuotaNonPagedPoolUsage ) , fprintf( stdout );
		printf( "\tPagefileUsage: 0x%08X\n", pmc.PagefileUsage ) , fprintf( stdout );
		printf( "\tPeakPagefileUsage: 0x%08X\n", pmc.PeakPagefileUsage ) , fprintf( stdout );
#endif
	}
	else
	{
		fprintf(stderr,"Failed to get process memory info\n");
//		PrintError();
	}
}

#endif // WIN32
