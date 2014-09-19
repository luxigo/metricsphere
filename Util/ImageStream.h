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
#ifndef IMAGE_STREAM_INCLUDED
#define IMAGE_STREAM_INCLUDED


static char JPEGFileExtension[] = "jpg";
static int DefaultOutputTileWidth = 8192;
static int DefaultOutputTileHeight = 8192;
static const char* DefaultOutputTileExtension = JPEGFileExtension;
static double GAMMA = 2.2;

#include "Util/Half/half.h"
#include "Util/GridStream.h"
#include "Util/Color.h"
#include "JPEGStream.h"

template<class Real> bool IsFloatType(void);

extern void*  BMPInitWrite			( char* fn , int width , int height , int quality , MultiStreamIOServer* ioServer );
extern void*  KROInitWriteColor		( char* fn , int width , int height , int quality , MultiStreamIOServer* ioServer );
extern void*  WDPInitWriteColor		( char* fn , int width , int height , int quality , MultiStreamIOServer* ioServer );
extern void*  WDPInitWriteGray		( char* fn , int width , int height , int quality , MultiStreamIOServer* ioServer );
extern void* JPEGInitWriteColor		( char* fn , int width , int height , int quality , MultiStreamIOServer* ioServer );
extern void* JPEGInitWriteGray		( char* fn , int width , int height , int quality , MultiStreamIOServer* ioServer );
extern void*  PNGInitWriteColor		( char* fn , int width , int height , int quality , MultiStreamIOServer* ioServer );
extern void*  PNGInitWriteGray		( char* fn , int width , int height , int quality , MultiStreamIOServer* ioServer );
extern void* TIFFInitWriteColor		( char* fn , int width , int height , int quality , MultiStreamIOServer* ioServer );
extern void* TIFFInitWriteGray		( char* fn , int width , int height , int quality , MultiStreamIOServer* ioServer );
extern void*  PNGInitWriteColorHDR	( char* fn , int width , int height , int quality , MultiStreamIOServer* ioServer );
extern void*  PNGInitWriteGrayHDR	( char* fn , int width , int height , int quality , MultiStreamIOServer* ioServer );
extern void* TIFFInitWriteColorHDR	( char* fn , int width , int height , int quality , MultiStreamIOServer* ioServer );
extern void* TIFFInitWriteGrayHDR	( char* fn , int width , int height , int quality , MultiStreamIOServer* ioServer );

extern void*  KROInitReadColor     ( char* fn , int& width , int& height , MultiStreamIOServer* ioServer );
extern void*  KROInitReadColorAlpha( char* fn , int& width , int& height , MultiStreamIOServer* ioServer );
extern void*  BMPInitReadColor     ( char* fn , int& width , int& height , MultiStreamIOServer* ioServer );
extern void*  WDPInitReadColor     ( char* fn , int& width , int& height , MultiStreamIOServer* ioServer );
extern void*  WDPInitReadGray      ( char* fn , int& width , int& height , MultiStreamIOServer* ioServer );
extern void* JPEGInitReadColor     ( char* fn , int& width , int& height , MultiStreamIOServer* ioServer );
extern void* JPEGInitReadGray      ( char* fn , int& width , int& height , MultiStreamIOServer* ioServer );
extern void*  PNGInitReadColor     ( char* fn , int& width , int& height , MultiStreamIOServer* ioServer );
extern void*  PNGInitReadColorAlpha( char* fn , int& width , int& height , MultiStreamIOServer* ioServer );
extern void*  PNGInitReadGray      ( char* fn , int& width , int& height , MultiStreamIOServer* ioServer );
extern void*  PFMInitReadColor     ( char *fn , int& width , int& height , MultiStreamIOServer* ioServer );
extern void* TIFFInitReadColor     ( char *fn , int& width , int& height , MultiStreamIOServer* ioServer );

extern void  KROGetImageSize( char* fn , int& width , int& height );
extern void  BMPGetImageSize( char* fn , int& width , int& height );
extern void  WDPGetImageSize( char* fn , int& width , int& height );
extern void JPEGGetImageSize( char* fn , int& width , int& height );
extern void  PNGGetImageSize( char* fn , int& width , int& height );
extern void  PFMGetImageSize( char* fn , int& width , int& height );
extern void TIFFGetImageSize( char* fn , int& width , int& height );

extern void  KROWriteRow(void* pixels,void* info,int j);
extern void  WDPWriteColorRow(void* pixels,void* info,int j);
extern void  WDPWriteGrayRow (void* pixels,void* info,int j);
extern void  BMPWriteRow(void* pixels,void* info,int j);
extern void JPEGWriteRow(void* pixels,void* info,int j);
extern void  PNGWriteRow(void* pixels,void* info,int j);
extern void TIFFWriteRow(void* pixels,void* info,int j);
extern void  PNGWriteRowHDR(void* pixels,void* info,int j);
extern void TIFFWriteRowHDR(void* pixels,void* info,int j);

extern void  KROReadColorRowAlpha(void* pixels,void* info,int j);
extern void  KROReadColorRow(void* pixels,void* info,int j);
extern void  BMPReadColorRow(void* pixels,void* info,int j);
extern void  WDPReadColorRow(void* pixels,void* info,int j);
extern void  WDPReadGrayRow (void* pixels,void* info,int j);
extern void JPEGReadRow(void* pixels,void* info,int j);
extern void  PNGReadColorRow(void* pixels,void* info,int j);
extern void  PNGReadGrayRow (void* pixels,void* info,int j);
extern void  PFMReadRow( void* pixels , void* info , int j );
extern void TIFFReadRow( void* pixels , void* info , int j );

extern void  KROFinalizeWrite( void* info );
extern void  WDPFinalizeWrite( void* info );
extern void  BMPFinalizeWrite( void* info );
extern void JPEGFinalizeWrite( void* info );
extern void  PNGFinalizeWrite( void* info );
extern void TIFFFinalizeWrite( void* info );

extern void  KROFinalizeRead( void* info );
extern void  BMPFinalizeRead( void* info );
extern void  WDPFinalizeRead( void* info );
extern void JPEGFinalizeRead( void* info );
extern void  PNGFinalizeRead( void* info );
extern void  PFMFinalizeRead( void* info );
extern void TIFFFinalizeRead( void* info );

template<class Real>	void* InitWrite( char* fn , int width , int height , int quality , MultiStreamIOServer* ioServer );
template<class Real>	void* InitRead ( char* fn , int& width , int& height , MultiStreamIOServer* ioServer );
template<class Real>	void GetImageSize( char* fn , int& width , int& height );
template<class Real>	void WriteRow(void* pixels,void* info);
template<class Real>	void ReadRow(void* pixels,void* info);
template<class Real>	void FinalizeWrite(void* info);
template<class Real>	void FinalizeRead(void* info);


template< class PReal , class BReal , int Channels=3 >
class WriteImageStream : public StreamingGrid
{
	Pointer( PReal ) pixels;
	Pointer( BReal ) buffer;

	int current;
	bool _clamp,_separate , _gammaEncode;
	int _q,_w,_h;
	void Init( char* fileName , MultiStreamIOServer* ioServer );
	double average[Channels];

	void* (*_InitWrite)		(char*,int,int,int,MultiStreamIOServer*);
	void  (*_WriteRow)		(void*,void*,int);
	void  (*_FinalizeWrite)	(void*);
protected:
	void* info;
public:
	WriteImageStream(
		void* (*IW)( char* , int , int , int , MultiStreamIOServer* ) , void (*WR)( void* , void* , int ) , void (*FW)( void* ),
		char* fileName , int width , int height , bool separateColors , bool gammaEncode , int quality , bool clamp ,
		MultiStreamIOServer* ioServer
		);
	~WriteImageStream		(void);
	int		rows			(void) const;
	int		rowSize			(void) const;
	Pointer( byte ) operator[] ( int idx );
	void	advance			(void);
	bool	SeparateColors	(void) const;
};

template< class PReal , class BReal , int Channels=3 >
class ReadImageStream : public StreamingGrid
{
	Pointer( PReal ) pixels;
	Pointer( BReal ) buffer;

	bool _separate , _gammaDecode;
	int current;
	int _w,_h;
	double average[Channels];

	void* (*_InitRead)		( char* , int& , int& , MultiStreamIOServer* );
	void  (*_ReadRow)		( void* , void* , int);
	void  (*_FinalizeRead)	( void* );
protected:
	void* info;
public:
	ReadImageStream	(
		void* (*IR)( char* , int& , int& , MultiStreamIOServer* ) , void (*RR)( void* , void* , int ) , void (*FR)( void* ),
		char* fileName , int &width , int &height , bool separateColors , bool gammaDecode ,
		MultiStreamIOServer* ioServer
		);
	~ReadImageStream		(void);
	int		rows			(void) const;
	int		rowSize			(void) const;
	Pointer( byte ) operator[] ( int idx );
	void	advance			(void);
	bool	SeparateColors	(void)	const;
	int		channels		(void)	const;
};

template<class Real,class OutType>
class WImageStream : public WriteImageStream<OutType,Real>
{
public:
	WImageStream		(char* fileName,int width,int height,bool separateColors , bool gammaEncode ,int quality , MultiStreamIOServer* ioServer );
};
template<class Real>
class KROWImageStream : public WriteImageStream<unsigned char,Real>
{
public:
	KROWImageStream		( char* fileName , int width , int height , bool separateColors , bool gammaEncode , int quality , MultiStreamIOServer* ioServer );
};
template<class Real>
class BMPWImageStream : public WriteImageStream<unsigned char,Real>
{
public:
	BMPWImageStream		(char* fileName,int width,int height,bool separateColors , bool gammaEncode ,int quality , MultiStreamIOServer* ioServer );
};
template<class Real>
class JPEGWImageStream : public WriteImageStream<unsigned char,Real>
{
public:
	JPEGWImageStream	(char* fileName,int width,int height,bool separateColors , bool gammaEncode ,int quality , MultiStreamIOServer* ioServer );
#if NEW_JPEG_IO
	int Service( void );
#endif // NEW_JPEG_IO
};
template<class Real>
class PNGWImageStream : public WriteImageStream<unsigned char,Real>
{
public:
	PNGWImageStream		(char* fileName,int width,int height,bool separateColors , bool gammaEncode ,int quality , MultiStreamIOServer* ioServer );
};
template<class Real>
class TIFFWImageStream : public WriteImageStream<unsigned char,Real>
{
public:
	TIFFWImageStream	(char* fileName,int width,int height,bool separateColors , bool gammaEncode ,int quality , MultiStreamIOServer* ioServer );
};
template<class Real>
class PNGWImageStreamHDR : public WriteImageStream<unsigned __int16,Real>
{
public:
	PNGWImageStreamHDR	(char* fileName,int width,int height,bool separateColors , bool gammaEncode ,int quality , MultiStreamIOServer* ioServer );
};
template<class Real>
class TIFFWImageStreamHDR : public WriteImageStream<unsigned __int16,Real>
{
public:
	TIFFWImageStreamHDR	(char* fileName,int width,int height,bool separateColors , bool gammaEncode ,int quality , MultiStreamIOServer* ioServer );
};

template<class Real>
class WDPWImageStream : public WriteImageStream<half,Real>
{
public:
	WDPWImageStream		(char* fileName,int width,int height,bool separateColors , bool gammaEncode ,int quality , MultiStreamIOServer* ioServer );
};

template<class Real,class InType>
class RImageStream : public ReadImageStream<InType,Real>
{
public:
	RImageStream	( char* fileName , int &width , int &height , bool separateColors , bool gammaDecode , MultiStreamIOServer* ioServer );
};
template<class Real>
class JPEGRImageStream : public ReadImageStream<unsigned char,Real>
{
public:
	JPEGRImageStream	(char* fileName,int &width,int &height,bool separateColors , bool gammaDecode , MultiStreamIOServer* ioServer );
#if NEW_JPEG_IO
	int Service( void );
#endif // NEW_JPEG_IO
};
template<class Real>
class KRORImageStream : public ReadImageStream< unsigned __int16 , Real >
{
public:
	KRORImageStream	( char* fileName , int &width , int &height , bool separateColors , bool gammaDecode , MultiStreamIOServer* ioServer );
};
template<class Real>
class KRORImageStreamAlpha : public ReadImageStream< unsigned __int16 , Real , 4 >
{
public:
	KRORImageStreamAlpha( char* fileName , int &width , int &height , bool separateColors , bool gammaDecode , MultiStreamIOServer* ioServer );
	bool HasAlpha( void ) const { return true; }
};
template<class Real>
class BMPRImageStream : public ReadImageStream< unsigned char , Real >
{
public:
	BMPRImageStream	(char* fileName,int &width,int &height,bool separateColors , bool gammaDecode , MultiStreamIOServer* ioServer );
};
template<class Real>
class PNGRImageStream : public ReadImageStream<unsigned __int16,Real>
{
public:
	PNGRImageStream	(char* fileName,int &width,int &height,bool separateColors , bool gammaDecode , MultiStreamIOServer* ioServer );
};
template<class Real>
class PNGRImageStreamAlpha : public ReadImageStream< unsigned __int16 , Real , 4 >
{
public:
	PNGRImageStreamAlpha (char* fileName,int &width,int &height,bool separateColors , bool gammaDecode , MultiStreamIOServer* ioServer );
	bool HasAlpha( void ) const { return true; }
};
template<class Real>
class WDPRImageStream : public ReadImageStream<half,Real>
{
public:
	WDPRImageStream	(char* fileName,int &width,int &height,bool separateColors , bool gammaDecode , MultiStreamIOServer* ioServer );
};
template<class Real>
class PFMRImageStream : public ReadImageStream< float , Real >
{
public:
	PFMRImageStream	( char* fileName , int &width , int &height , bool separateColors , bool gammaDecode , MultiStreamIOServer* ioServer );
};
template< class Real >
class TIFFRImageStream : public ReadImageStream< unsigned __int16 , Real >
{
public:
	TIFFRImageStream( char* fileName , int &width , int &height , bool separateColors , bool gammaDecode , MultiStreamIOServer* ioServer );
};
template<class Real>
class StreamingInputTile
{
	StreamingGrid* data;
	int w,h,rc;
	int sX,sY;
	char tileName[1024];
	Real** rows;
	double average[3];
public:
	StreamingInputTile(void);
	~StreamingInputTile(void);

	void Init( const char* tName , int rowCount=1 , int startX=0 , int startY=0 , bool separateColors=false , bool gammaDecode=false );
	void init(int idx);
	Color< Real > operator() ( int i , int j );
	Color< Real > operator() ( int i , int j , Real& a );
	int width(void) const;
	int height(void) const;
	int startX(void) const;
	int startY(void) const;
};

template<class Real>
class RGridOfImages : public StreamingGrid
{
	int _c , _r , _w , _h;
	int *_widths , *_heights;
	int _current;
	bool _separate , _gammaDecode;
	StreamingGrid*** _grid;
	char*** _fileNames;
	Pointer( Real ) _buffer;
#if NEW_JPEG_IO
	MultiStreamIOServer* _server;
#endif // NEW_JPEG_IO
public:
	static void GetImageSize( char* fileName , int& width , int& height );
	RGridOfImages		(char* fileName , int& width , int& height , bool separateColors , bool gammaDecode , MultiStreamIOServer* ioServer );
	~RGridOfImages		(void);
	int		rows		(void) const;
	int		rowSize		(void) const;
	Pointer( byte ) operator[] ( int idx );
	void	advance		(void);
	bool	SeparateColors	(void) const;
};
template<class Real>
class WGridOfImages : public StreamingGrid
{
	int _c,_r,_w,_h;
	int *_widths,*_heights;
	int _current,_quality;
	bool _separate , _hdr , _gammaEncode;
	StreamingGrid*** _grid;
	char*** _fileNames;
	Pointer( Real ) _buffer;
#if NEW_JPEG_IO
	MultiStreamIOServer* _server;
#endif // NEW_JPEG_IO
	void _init(char* fileName,const char* header,const char* ext,int width,int height,bool separateColors , bool gammaEncode ,int quality,int tileWidth,int tileHeight , MultiStreamIOServer* ioServer , bool hdr );
public:
//	WGridOfImages		(char* fileName,int width,int height,bool separateColors,int quality,int tileWidth=DEFAULT_TILE_SIZE,int tileHeight=DEFAULT_TILE_SIZE);
//	WGridOfImages		(char* fileName,const char* header,const char* ext,int width,int height,bool separateColors,int quality,int tileWidth=DEFAULT_TILE_SIZE,int tileHeight=DEFAULT_TILE_SIZE);
	WGridOfImages		(char* fileName,int width,int height,bool separateColors , bool gammaEncode ,int quality,int tileWidth,int tileHeight,MultiStreamIOServer* ioServer , bool hdr );
	WGridOfImages		(char* fileName,const char* header,const char* ext,int width,int height,bool separateColors , bool gammaEncode ,int quality,int tileWidth,int tileHeight,MultiStreamIOServer* ioServer , bool hdr );
	~WGridOfImages		(void);
	int		rows		(void) const;
	int		rowSize		(void) const;
	Pointer( byte ) operator[] ( int idx );
	void	advance		(void);
	bool	SeparateColors	(void) const;
};

template<class Real>
class MultiThreadedRGridOfImages : public StreamingGrid
{
	HANDLE *_imageReaderHandles;
	HANDLE *_loadHandles , *_unloadHandles;
	int _c , _r , _w , _h;
	int *_widths , *_heights;
	int _current , _threads;
	bool _separate , _gammaDecode;
	char*** _fileNames;
	Pointer( Real ) _buffer;
	bool _dataReady;

	class ImageReader
	{
		HANDLE _loadHandle , _unloadHandle;
		int _c , _r , _w , _h;
		int *_widths , *_heights;
		bool _separate , _gammaDecode;
		char*** _fileNames;
		Pointer( Real ) _buffer;
#if NEW_JPEG_IO
		MultiStreamIOServer* _server;
#endif // NEW_JPEG_IO
	public:
		void init( int columns , int rows , int* widths , int* heights , char*** fileNames , bool separateColors , bool gammaDecode , Pointer( Real ) buffer , HANDLE loadHandle , HANDLE unloadHandle , MultiStreamIOServer* ioServer );
		void Run( void );
	};
	ImageReader* _imageReaders;
public:
	static void GetImageSize( char* fileName , int& width , int& height );
	MultiThreadedRGridOfImages	( char* fileName , int& width , int& height , bool separateColors , bool gammaDecode , int threads , int threadPriority , MultiStreamIOServer* ioServer );
	~MultiThreadedRGridOfImages ( void );
	int		rows		( void ) const;
	int		rowSize		( void ) const;
	Pointer( byte ) operator[] ( int idx );
	void	advance		( void );
	bool	SeparateColors( void ) const;
};
template<class Real>
class MultiThreadedWGridOfImages : public StreamingGrid
{
	HANDLE *_imageWriterHandles;
	HANDLE *_loadHandles , *_unloadHandles;
	int _c , _r , _w , _h;
	int *_widths , *_heights;
	int _current , _quality , _threads;
	bool _separate , _hdr , _gammaEncode;
	char*** _fileNames;
	Pointer( Real ) _buffer;
	bool _dataReady;

	class ImageWriter
	{
		HANDLE _loadHandle , _unloadHandle;
		int _c , _r , _w , _h;
		int *_widths , *_heights;
		int _quality;
		bool _separate , _hdr , _gammaEncode;
		char*** _fileNames;
		Pointer( Real ) _buffer;
#if NEW_JPEG_IO
		MultiStreamIOServer* _server;
#endif // NEW_JPEG_IO
	public:
		void init( int columns , int rows , int* widths , int* heights , char*** fileNames , bool separateColors , bool gammaEncode , int quality , Pointer( Real ) buffer , HANDLE loadHandle , HANDLE unloadHandle , MultiStreamIOServer* ioServer , bool hdr );
		void Run( void );
	};
	ImageWriter* _imageWriters;

	void _init( char* fileName , const char* header , const char* ext , int width , int height , bool separateColors , bool gammaEncode , int quality , int threads , int threadPriority , int tileWidth , int tileHeight , MultiStreamIOServer* ioServer , bool hdr );
public:
//	MultiThreadedWGridOfImages	( char* fileName , int width , int height , bool separateColors , int quality , int threads , int tileWidth=DEFAULT_TILE_SIZE , int tileHeight=DEFAULT_TILE_SIZE );
//	MultiThreadedWGridOfImages	( char* fileName , const char* header , const char* ext , int width , int height , bool separateColors , int quality , int threads , int tileWidth=DEFAULT_TILE_SIZE , int tileHeight=DEFAULT_TILE_SIZE );
	MultiThreadedWGridOfImages	( char* fileName , int width , int height , bool separateColors  , bool gammaEncode , int quality , int threads , int threadPriority , int tileWidth , int tileHeight , MultiStreamIOServer* ioServer , bool  hdr );
	MultiThreadedWGridOfImages	( char* fileName , const char* header , const char* ext , int width , int height , bool separateColors , bool gammaEncode , int quality , int threads , int threadPriority , int tileWidth , int tileHeight , MultiStreamIOServer* ioServer , bool hdr );
	~MultiThreadedWGridOfImages ( void );
	int		rows		(void) const;
	int		rowSize		(void) const;
	Pointer( byte ) operator[] ( int idx );
	void	advance		(void);
	bool	SeparateColors	(void) const;
};


template<class Real>
class WImageClipper : public StreamingGrid
{
	bool deleteGrid;
	int _w,_h,_ox,_oy,_ow,_oh;
	int current;
	bool _separate , _hdr;
	StreamingGrid* grid;
	Pointer( Real ) buffer;
public:
	WImageClipper		( char* fileName , int width , int height , bool separateColors , bool gammaEncode , int outOffX , int outOffY , int outWidth , int outHeight , int quality , MultiStreamIOServer* server , bool hdr=false , int threads=1 , int threadPriority=THREAD_PRIORITY_NORMAL );
	WImageClipper		( StreamingGrid* outGrid , int width , int height , bool separateColors , int outOffX , int outOffY , int outWidth , int outHeight , bool hdr=false );
	~WImageClipper		(void);
	int		rows		(void)	const;
	int		rowSize		(void) const;
	Pointer( byte ) operator[] ( int idx );
	void	advance		(void);
	bool	SeparateColors	(void) const;
};

template< class Real >
class RChannelExtractor : public StreamingGrid
{
	int _c,_w,_h;
	int current;
	bool _separate;
	StreamingGrid* grid;
	Real* buffer;
public:
	RChannelExtractor	( char* fileName , int& width , int& height , bool separateColors , bool gammaDecode , int channel);
	~RChannelExtractor	(void);
	int		rows		(void)	const;
	int		rowSize		(void) const;
	Pointer( byte ) operator[] ( int idx );
	void	advance		(void);
	bool	SeparateColors	(void) const;
};
template< class Real >
class RImageSampler : public StreamingGrid
{
	int _inW , _inH , _outW , _outH , _inCurrent , _outCurrent;
	bool _separate , _nearest;
	StreamingGrid* _grid;
	Pointer( Real ) _inBuffers[2];
	Pointer( Real ) _outBuffer;
public:
	RImageSampler( char* fileName , int width , int height , bool nearest , bool separateColors , bool gammaDecode , Real logScale , MultiStreamIOServer* ioServer , int threads=0 , int threadPriority = THREAD_PRIORITY_NORMAL );
	~RImageSampler( void );
	int rows( void ) const;
	int rowSize( void ) const;
	Pointer( byte ) operator[] ( int idx );
	void advance( void );
	bool SeparateColors( void ) const;
};
template< class Real >
class WImageSampler : public StreamingGrid
{
	int _inW , _inH , _outW , _outH , _inCurrent , _outCurrent;
	bool _separate , _nearest;
	StreamingGrid* _grid;
	Pointer( Real ) _inBuffers[2];
public:
	WImageSampler( char* fileName , int inWidth , int inHeight , int outWidth , int outHeight , bool nearest , bool separateColors , bool gammaEncode , bool logScale , const int& quality , MultiStreamIOServer* ioServer , int threads=0 , int threadPriority = THREAD_PRIORITY_NORMAL );
	~WImageSampler( void );
	int rows( void ) const;
	int rowSize( void ) const;
	Pointer( byte ) operator[] ( int idx );
	void advance( void );
	bool SeparateColors( void ) const;
};
template< class Real >
class RImageSampler2 : public StreamingGrid
{
	int _inW , _inH , _outW , _outH , _inCurrent , _outCurrent;
	bool _separate;
	StreamingGrid *_nearestGrid  , *_averageGrid;
	Pointer( Real ) _nearestBuffers[2];
	Pointer( Real ) _averageBuffers[2];
	Pointer( Real ) _outBuffer;
public:
	RImageSampler2( char* nearestFileName , char* averageFileName , int width , int height , bool separateColors , bool gammaDecode , MultiStreamIOServer* ioServer , int threads=0 , int threadPriority = THREAD_PRIORITY_NORMAL );
	~RImageSampler2( void );
	int rows( void ) const;
	int rowSize( void ) const;
	Pointer( byte ) operator[] ( int idx );
	void advance( void );
	bool SeparateColors( void ) const;
};
template< class NearestReal , class AverageReal >
class JointRImageSampler
{
	int _inW , _inH , _outW , _outH , _inCurrent , _outCurrent;
	bool _separate , _nearest;
	StreamingGrid *_nearestGrid , *_averageGrid;
	Pointer( NearestReal ) _inNearestBuffers[2];
	Pointer( AverageReal ) _inAverageBuffers[2];
	Pointer( NearestReal ) _outNearestBuffers[2];
	Pointer( AverageReal ) _outAverageBuffers[2];
	struct JointRImageSamplerChild : public StreamingGrid
	{
		bool isNearest;
		JointRImageSampler* parent;
		int rows   ( void ) const { return parent->rows(); }
		int rowSize( void ) const { if( isNearest ) return parent->nearestRowSize() ; else return parent->averageRowSize(); }
		Pointer( byte ) operator[] ( int idx ) { if( isNearest ) return parent->nearestRow( idx ) ; else return parent->averageRow( idx ); }
		void advance( void ) { if( isNearest ) parent->advance(); }
		bool SeparateColors( void ) const { return parent->SeparateColors(); }
	};
public:
	JointRImageSampler( char* nearestFileName , char* averageFileName , int width , int height , bool separateColors , bool gammaDecode , double logScale , MultiStreamIOServer* ioServer , int threads=0 , int threadPriority = THREAD_PRIORITY_NORMAL );
	~JointRImageSampler( void );
	int rows( void ) const;
	int nearestRowSize( void ) const;
	int averageRowSize( void ) const;
	Pointer( byte ) nearestRow( int idx );
	Pointer( byte ) averageRow( int idx );
	void advance( void );
	bool SeparateColors( void ) const;
	JointRImageSamplerChild *nearestChild , *averageChild;
};
template< class Real >
class LogImageReader : public StreamingGrid
{
	StreamingGrid*  _grid;
	Pointer( byte ) _buffer;
	Real _min;
public:
	LogImageReader( char* fileName , int& width , int& height , Real min , bool separateColors , bool gammaDecode , bool getAlpha , MultiStreamIOServer* ioServer , int threads=0 , int threadPriority = THREAD_PRIORITY_NORMAL );
	~LogImageReader( void );
	int rows( void ) const;
	int rowSize( void ) const;
	Pointer( byte ) operator[]( int idx );
	void advance( void ) ;
	bool SeparateColors( void ) const;
	bool HasAlpha( void ) const;
};
template< class Real >
class LogImageWriter : public StreamingGrid
{
	int _index;
	StreamingGrid*  _grid;
	Pointer( byte ) _buffer;
public:
	LogImageWriter( char* fileName , const int& width , const int& height , const bool& separateColors , bool gammaEncode  , const int& quality , MultiStreamIOServer* ioServer , bool hdr=false , int threads=0 , int threadPriority = THREAD_PRIORITY_NORMAL );
	~LogImageWriter( void );
	int rows( void ) const;
	int rowSize( void ) const;
	Pointer( byte ) operator[]( int idx );
	void advance( void ) ;
	bool SeparateColors( void ) const;
	bool HasAlpha( void ) const;
};
template< class Real >
class LogLuminanceImageReader : public StreamingGrid
{
	StreamingGrid*  _grid;
	Pointer( byte ) _buffer;
	Real _min;
public:
	LogLuminanceImageReader( char* fileName , int& width , int& height , Real min , bool separateColors , bool gammaDecode , MultiStreamIOServer* ioServer , int threads=0 , int threadPriority = THREAD_PRIORITY_NORMAL );
	~LogLuminanceImageReader( void );
	int rows( void ) const;
	int rowSize( void ) const;
	Pointer( byte ) operator[]( int idx );
	void advance( void ) ;
	bool SeparateColors( void ) const;
	bool HasAlpha( void ) const;
};
template< class Real >
class LogLuminanceImageWriter : public StreamingGrid
{
	int _index;
	StreamingGrid*  _grid;
	Pointer( byte ) _buffer;
public:
	LogLuminanceImageWriter( char* fileName , const int& width , const int& height , const bool& separateColors , bool gammaEncode  , const int& quality , MultiStreamIOServer* ioServer , bool hdr=false , int threads=0 , int threadPriority = THREAD_PRIORITY_NORMAL );
	~LogLuminanceImageWriter( void );
	int rows( void ) const;
	int rowSize( void ) const;
	Pointer( byte ) operator[]( int idx );
	void advance( void ) ;
	bool SeparateColors( void ) const;
	bool HasAlpha( void ) const;
};

template<class Real>	void GetReadSize( char* fileName , int& width , int& height );
#if NEW_JPEG_IO
template<class Real>	StreamingGrid* GetReadStream ( char* fileName ,       int& width ,       int& height , const bool& separateColors , bool gammaDecode , const bool& getAlpha , MultiStreamIOServer* ioServer ,                  int threads=0 , int threadPriority = THREAD_PRIORITY_NORMAL );
template<class Real>	StreamingGrid* GetWriteStream( char* fileName , const int& width , const int& height , const bool& separateColors , bool gammaEncode , const int& quality   , MultiStreamIOServer* ioServer , bool hdr=false , int threads=0 , int threadPriority = THREAD_PRIORITY_NORMAL );
#else // !NEW_JPEG_IO
template<class Real>	StreamingGrid* GetReadStream ( char* fileName ,      int& width ,       int& height , const bool& separateColors , bool gammaDecode , const bool& getAlpha , MultiStreamIOServer* ioServer , int threads=0 );
template<class Real>	StreamingGrid* GetWriteStream( char* fileName ,const int& width , const int& height , const bool& separateColors , bool gammaEncode , const int& quality   , int threads=0 );
#endif // NEW_JPEG_IO

#include "Stream.inl"
#include "ImageStream.inl"
#endif // GRID_STREAM_INCLUDED