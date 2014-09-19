/*
Copyright (c) 2008, Michael Kazhdan
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
#include <Util/BaseMultiStreamIO.h>

template<class Real>
struct WriteInfo
{
	FILE* fp;
	int width;
};
template<class Real>
struct ReadInfo
{
	FILE* fp;
	int width;
};

template<class Real>
void* InitWrite(char* fileName,int width,int height,int quality , MultiStreamIOServer* ioServer )
{
	WriteInfo<Real>* info = (WriteInfo<Real>*)malloc(sizeof(WriteInfo<Real>));
	info->fp = fopen(fileName,"wb");
	if(!info->fp)	fprintf(stderr,"Failed to open: %s\n",fileName)	,	exit(0);

	fwrite(&width,sizeof(int),1,info->fp);
	fwrite(&height,sizeof(int),1,info->fp);
	info->width=width;
	return info;
}
template<class Real>
void WriteRow(void* pixels,void* v,int j)
{
	WriteInfo<Real>* info = (WriteInfo<Real>*)v;
	fwrite(pixels,sizeof(Real),3*info->width,info->fp);
}
template<class Real>
void FinalizeWrite(void* v)
{
	WriteInfo<Real>* info = (WriteInfo<Real>*)v;
	fclose(info->fp);
	free(info);
}
template<class Real>
void GetImageSize(char* fileName,int& width,int& height)
{
	ReadInfo<Real> info;
	info.fp = fopen( fileName , "rb" );
	if( !info.fp )	fprintf(stderr,"Failed to open: %s\n",fileName)	,	exit(0);
	fread( &width  , sizeof(int) , 1 , info.fp );
	fread( &height , sizeof(int) , 1 , info.fp );
	fclose( info.fp );
}
template<class Real>
void* InitRead(char* fileName,int& width,int& height , MultiStreamIOServer* ioServer )
{
	ReadInfo<Real>* info = (ReadInfo<Real>*)malloc(sizeof(ReadInfo<Real>));
	info->fp=fopen(fileName,"rb");
	if(!info->fp)	fprintf(stderr,"Failed to open: %s\n",fileName)	,	exit(0);
	fread(&width,sizeof(int),1,info->fp);
	fread(&height,sizeof(int),1,info->fp);
	info->width=width;
	return info;
}
template<class Real>
void ReadRow(void* pixels,void* v,int j)
{
	ReadInfo<Real>* info = (ReadInfo<Real>*)v;
	fread(pixels,sizeof(Real),3*info->width,info->fp);
}
template<class Real>
void FinalizeRead(void* v)
{
	ReadInfo<Real>* info = (ReadInfo<Real>*)v;
	fclose(info->fp);
	free(info);
}
