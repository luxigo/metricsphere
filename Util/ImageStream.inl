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
#include "JPEGStream.h"
#include "ChannelConverter.h"


//////////////////////
// WriteImageStream //
//////////////////////
template< class PReal , class BReal , int Channels >
WriteImageStream< PReal , BReal , Channels >::WriteImageStream( void* (*IW)( char* , int , int , int , MultiStreamIOServer* ),void (*WR)( void* , void* , int) , void (*FW)( void* ),					
															    char* fileName , int width , int height , bool separateColors , bool gammaEncode , int quality , bool clamp ,
															    MultiStreamIOServer* ioServer ) 
{
	_InitWrite		= IW;
	_WriteRow		= WR;
	_FinalizeWrite	= FW;
	current=0;
	_q=quality;
	_w=width;
	_h=height;
	_clamp=clamp;
	_gammaEncode = gammaEncode;
	_separate=separateColors;
	info=NULL;
	pixels = NullPointer< PReal >( );
	buffer = NullPointer< BReal >( );
	Init( fileName , ioServer );
}
template< class PReal , class BReal , int Channels >
void WriteImageStream< PReal , BReal , Channels >::Init( char* fileName , MultiStreamIOServer* ioServer )
{
	memset( average , 0 , sizeof( average ) );

	FreeArray( pixels );
	FreeArray( buffer );
	pixels = AllocArray< PReal >( _w * Channels , 1 , "WriteImageStream::pixels" );
	buffer = AllocArray< BReal >( _w * Channels , 1 , "WriteImageStream::buffer" );
	if(!buffer)	fprintf( stderr , "Failed to allocate buffer: %d * %d\n" , _w , Channels ) , exit(0);
	if(!pixels)	fprintf( stderr , "Failed to allocate pixels: %d * %d\n" , _w , Channels ) , exit(0);
	info=_InitWrite(fileName,_w,_h,_q,ioServer);
}

template< class PReal , class BReal , int Channels >
WriteImageStream< PReal , BReal , Channels >::~WriteImageStream(void)
{
	FreeArray( buffer );
	FreeArray( pixels );
	_FinalizeWrite(info);
	for( int c=0 ; c<Channels ; c++ ) average[c] /= _w , average[c] /= _h;
//	printf("Average: %f %f %f\n",average[0],average[1],average[2]) , fflush( stdout );
}
template< class PReal , class BReal , int Channels > int  WriteImageStream< PReal , BReal , Channels >::rows(void) const {return _h;}
template< class PReal , class BReal , int Channels > int  WriteImageStream< PReal , BReal , Channels >::rowSize(void) const {return _w*Channels*sizeof(BReal);}
template< class PReal , class BReal , int Channels > bool WriteImageStream< PReal , BReal , Channels >::SeparateColors(void)	const	{return _separate;}
template< class PReal , class BReal , int Channels >
Pointer( byte ) WriteImageStream< PReal , BReal , Channels >::operator[] (int idx)
{
#if ASSERT_MEMORY_ACCESS
	if(idx!=current)	fprintf(stderr,"Index out of bounds: %d != %d\n",idx,current) , exit(0);
#endif // ASSERT_MEMORY_ACCESS
	return ( Pointer( byte ) )buffer;
}
template< class PReal , class BReal , int Channels >
void WriteImageStream< PReal , BReal , Channels >::advance(void)
{
	if(current<_h)
	{
		if( _gammaEncode )
			for( int x=0 ; x<_w*Channels ; x++ )
				buffer[ x ] = ConvertChannel< double , BReal >( pow( ConvertChannel< BReal , double >( buffer[x] ) , 1.0/GAMMA ) );

		if(!_clamp)
			if(_separate)
			{
				for(int x=0;x<_w;x++)
					for( int c=0 ; c<Channels ; c++ )
					{
						pixels[x*Channels+c]=ConvertChannel< BReal , PReal >( buffer[ x+c*_w] );
						average[c] += pixels[x*Channels+c];
					}
			}
			else
			{
				for(int x=0;x<_w;x++)
					for( int c=0 ; c<Channels ; c++ )
					{
						pixels[x*Channels+c] = ConvertChannel<BReal,PReal>( buffer[x*Channels+c] );
						average[c] += pixels[x*Channels+c];
					}
			}
		else
			for(int x=0;x<_w;x++)
			{
				double _values[ Channels ];

				if( _separate ) for( int c=0 ; c<Channels ; c++ ) _values[c] = ConvertChannel<BReal,double>(buffer[x+c*_w]);
				else            for( int c=0 ; c<Channels ; c++ ) _values[c] = ConvertChannel<BReal,double>(buffer[x*Channels+c]);
				for( int c=0 ; c<Channels ; c++ )
				{
					if		(_values[c]<0)	_values[c]=0;
					else if	(_values[c]>1)	_values[c]=1;
					pixels[x*Channels+c]=ConvertChannel< double , PReal >( _values[c] );
					average[c] += pixels[x*Channels+c];
				}
			}
		_WriteRow( PointerAddress( pixels ) , info , current );
	}
	current++;
}
/////////////////////
// ReadImageStream //
/////////////////////
template< class PReal , class BReal , int Channels >
ReadImageStream< PReal , BReal , Channels >::ReadImageStream( void* (*IR)( char* , int& , int& , MultiStreamIOServer* ) , void (*RR)( void* , void* , int ) , void (*FR)( void* ),					
															  char* fileName , int& width , int& height , bool separateColors , bool gammaDecode ,
															  MultiStreamIOServer* ioServer )
{
	_InitRead		= IR;
	_ReadRow		= RR;
	_FinalizeRead	= FR;

	_gammaDecode	= gammaDecode;
	_separate		= separateColors;
	current=-1;

	info=_InitRead(fileName,width,height , ioServer );
	_w=width;
	_h=height;

	pixels = AllocArray< PReal >( _w * Channels , 1 , "ReadImageStream::pixels" );
	buffer = AllocArray< BReal >( _w * Channels , 1 , "ReadImageStream::buffer" );

	if(!buffer)	fprintf( stderr , "Failed to allocate buffer\n" ) , exit(0);
	if(!pixels)	fprintf( stderr , "Failed to allocate pixels\n" ) , exit(0);

	advance();
}
template< class PReal , class BReal , int Channels >
ReadImageStream< PReal , BReal , Channels >::~ReadImageStream(void)
{
	FreeArray( pixels );
	FreeArray( buffer );
	_FinalizeRead(info);
	for( int c=0 ; c<Channels ; c++ ) average[c] /= _w , average[c] /= _h;
//	printf("Average: %f %f %f\n",average[0],average[1],average[2]) , fflush( stdout );

}
template< class PReal , class BReal , int Channels >	int  ReadImageStream< PReal , BReal , Channels >::rows(void) const {return _h;}
template< class PReal , class BReal , int Channels >	int  ReadImageStream< PReal , BReal , Channels >::rowSize(void) const {	return _w*Channels*sizeof(BReal); }
template< class PReal , class BReal , int Channels >	bool ReadImageStream< PReal , BReal , Channels >::SeparateColors(void) const { return _separate; }
template< class PReal , class BReal , int Channels >	int  ReadImageStream< PReal , BReal , Channels >::channels( void ) const { return Channels; }
template< class PReal , class BReal , int Channels >
Pointer( byte ) ReadImageStream< PReal , BReal , Channels >::operator[] (int idx)
{
#if ASSERT_MEMORY_ACCESS
	if(idx!=current)	fprintf(stderr,"Index out of bounds: %d != %d\n",idx,current) , exit(0);
#endif // ASSERT_MEMORY_ACCESS
	return ( Pointer( byte ) ) buffer;
}
template< class PReal , class BReal , int Channels >
void ReadImageStream< PReal , BReal , Channels >::advance(void)
{
	current++;
	if( current<_h )
	{
		_ReadRow( PointerAddress( pixels ) , info , current );
		if( _separate )
			for(int x=0;x<_w;x++)
			{
				for( int c=0 ; c<Channels ; c++ )
				{
					average[c] += double(pixels[x*Channels+c]);
					buffer[x+c*_w] = ConvertChannel< PReal , BReal >( pixels[x*Channels+c]);
				}
			}
		else
			for(int x=0;x<_w;x++)
			{
				for( int c=0 ; c<Channels ; c++ )
				{
					average[c]+=double(pixels[x*Channels+c]);
					buffer[x*Channels+c] = ConvertChannel<PReal,BReal>(pixels[x*Channels+c]);
				}
			}
		if( _gammaDecode )
			for( int x=0 ; x<_w*Channels ; x++ )
				buffer[ x ] = ConvertChannel< double , BReal >( pow( ConvertChannel< BReal , double >( buffer[x] ) , GAMMA ) );
	}
}

//////////////////
// WImageStream //
//////////////////
template<class Real,class OutType>
WImageStream<Real,OutType>::WImageStream(char* fileName,int width,int height,bool separateColors , bool gammaEncode , int quality , MultiStreamIOServer* ioServer ) :
WriteImageStream( InitWrite<OutType> , WriteRow<OutType> , FinalizeWrite<OutType> , fileName , width , height , separateColors , gammaEncode , quality , false , ioServer )
{
}
/////////////////////
// KROWImageStream //
/////////////////////
template<class Real>
KROWImageStream<Real>::KROWImageStream(char* fileName,int width,int height , bool separateColors , bool gammaEncode , int quality , MultiStreamIOServer* ioServer ) :
WriteImageStream( KROInitWriteColor , KROWriteRow , KROFinalizeWrite , fileName , width , height , separateColors , gammaEncode , quality,true , ioServer )
{
}
/////////////////////
// BMPWImageStream //
/////////////////////
template<class Real>
BMPWImageStream<Real>::BMPWImageStream(char* fileName,int width,int height , bool separateColors , bool gammaEncode , int quality , MultiStreamIOServer* ioServer ) :
WriteImageStream(BMPInitWrite,BMPWriteRow,BMPFinalizeWrite,fileName,width,height,separateColors , gammaEncode , quality,true , ioServer )
{
}
/////////////////////
// JPEGWImageSteam //
/////////////////////
template<class Real>
JPEGWImageStream< Real >::JPEGWImageStream( char* fileName , int width , int height , bool separateColors , bool gammaEncode , int quality , MultiStreamIOServer* ioServer ) :
WriteImageStream(JPEGInitWriteColor,JPEGWriteRow,JPEGFinalizeWrite,fileName,width,height,separateColors, gammaEncode , quality,true , ioServer )
{
}
#if NEW_JPEG_IO
template< class Real >
int JPEGWImageStream< Real >::Service( void )
{
	JPEGWriteInfo* inf = (JPEGWriteInfo*) info;
	return inf->outStream->Service();
}
#endif // NEW_JPEG_IO
////////////////////
// PNGWImageSteam //
////////////////////
template<class Real>
PNGWImageStream<Real>::PNGWImageStream(char* fileName,int width,int height,bool separateColors , bool gammaEncode , int quality , MultiStreamIOServer* ioServer ) :
WriteImageStream(PNGInitWriteColor,PNGWriteRow,PNGFinalizeWrite,fileName,width,height,separateColors , gammaEncode , quality,true , ioServer )
{
}
/////////////////////
// TIFFWImageSteam //
/////////////////////
template<class Real>
TIFFWImageStream<Real>::TIFFWImageStream(char* fileName,int width,int height,bool separateColors , bool gammaEncode , int quality , MultiStreamIOServer* ioServer ) :
WriteImageStream(TIFFInitWriteColor,TIFFWriteRow,TIFFFinalizeWrite,fileName,width,height,separateColors , gammaEncode , quality,true , ioServer )
{
}
///////////////////////
// PNGWImageSteamHDR //
///////////////////////
template<class Real>
PNGWImageStreamHDR<Real>::PNGWImageStreamHDR(char* fileName,int width,int height,bool separateColors , bool gammaEncode , int quality , MultiStreamIOServer* ioServer ) :
WriteImageStream(PNGInitWriteColorHDR,PNGWriteRowHDR,PNGFinalizeWrite,fileName,width,height,separateColors , gammaEncode , quality,true , ioServer )
{
}
////////////////////////
// TIFFWImageSteamHDR //
////////////////////////
template<class Real>
TIFFWImageStreamHDR<Real>::TIFFWImageStreamHDR(char* fileName,int width,int height,bool separateColors , bool gammaEncode , int quality , MultiStreamIOServer* ioServer ) :
WriteImageStream(TIFFInitWriteColorHDR,TIFFWriteRowHDR,TIFFFinalizeWrite,fileName,width,height,separateColors , gammaEncode , quality,true , ioServer )
{
}
/////////////////////
// WDPWImageStream //
/////////////////////
template<class Real>
WDPWImageStream<Real>::WDPWImageStream(char* fileName,int width,int height,bool separateColors , bool gammaEncode , int quality , MultiStreamIOServer* ioServer ) :
WriteImageStream(WDPInitWriteColor,WDPWriteColorRow,WDPFinalizeWrite,fileName,width,height,separateColors , gammaEncode , quality,false , ioServer )
{
}
//////////////////
// RImageStream //
//////////////////
template<class Real,class InType>
RImageStream<Real,InType>::RImageStream( char* fileName , int& width , int& height , bool separateColors , bool gammaDecode , MultiStreamIOServer* ioServer ) :
ReadImageStream( InitRead<InType> , ReadRow<InType> , FinalizeRead<InType> , fileName , width , height , separateColors , gammaDecode , ioServer )
{
}
//////////////////////
// JPEGRImageStream //
//////////////////////
template<class Real>
JPEGRImageStream<Real>::JPEGRImageStream( char* fileName , int& width , int& height , bool separateColors , bool gammaDecode , MultiStreamIOServer* ioServer ) :
ReadImageStream( JPEGInitReadColor , JPEGReadRow , JPEGFinalizeRead , fileName , width , height , separateColors , gammaDecode , ioServer )
{
}
#if NEW_JPEG_IO
template< class Real >
int JPEGRImageStream< Real >::Service( void )
{
	JPEGReadInfo* inf = (JPEGReadInfo*) info;
	return inf->inStream->Service();
}
#endif // NEW_JPEG_IO
/////////////////////
// KRORImageStream //
/////////////////////
template< class Real >
KRORImageStream< Real >::KRORImageStream( char* fileName , int& width , int& height , bool separateColors , bool gammaDecode , MultiStreamIOServer* ioServer ) :
ReadImageStream( KROInitReadColor , KROReadColorRow , KROFinalizeRead , fileName , width , height , separateColors , gammaDecode , ioServer )
{
}
//////////////////////////
// KRORImageStreamAlpha //
//////////////////////////
template< class Real >
KRORImageStreamAlpha< Real >::KRORImageStreamAlpha( char* fileName , int& width , int& height , bool separateColors , bool gammaDecode , MultiStreamIOServer* ioServer ) :
ReadImageStream( KROInitReadColorAlpha , KROReadColorRowAlpha , KROFinalizeRead , fileName , width , height , separateColors , gammaDecode , ioServer )
{
}
/////////////////////
// BMPRImageStream //
/////////////////////
template< class Real >
BMPRImageStream< Real >::BMPRImageStream( char* fileName , int& width , int& height , bool separateColors , bool gammaDecode , MultiStreamIOServer* ioServer ) :
ReadImageStream( BMPInitReadColor , BMPReadColorRow , BMPFinalizeRead , fileName , width , height , separateColors , gammaDecode , ioServer )
{
}
/////////////////////
// PNGRImageStream //
/////////////////////
template< class Real >
PNGRImageStream< Real >::PNGRImageStream( char* fileName , int& width , int& height , bool separateColors , bool gammaDecode , MultiStreamIOServer* ioServer ) :
ReadImageStream( PNGInitReadColor , PNGReadColorRow , PNGFinalizeRead , fileName , width , height , separateColors , gammaDecode , ioServer )
{
}
//////////////////////////
// PNGRImageStreamAlpha //
//////////////////////////
template< class Real >
PNGRImageStreamAlpha< Real >::PNGRImageStreamAlpha( char* fileName , int& width , int& height , bool separateColors , bool gammaDecode , MultiStreamIOServer* ioServer ) :
ReadImageStream( PNGInitReadColorAlpha , PNGReadColorRow , PNGFinalizeRead , fileName , width , height , separateColors , gammaDecode , ioServer )
{
}
/////////////////////
// WDPRImageStream //
/////////////////////
template<class Real>
WDPRImageStream<Real>::WDPRImageStream(char* fileName,int& width,int& height,bool separateColors , bool gammaDecode , MultiStreamIOServer* ioServer ) :
ReadImageStream( WDPInitReadColor , WDPReadColorRow , WDPFinalizeRead , fileName , width , height , separateColors , gammaDecode , ioServer )
{
}
//////////////////////
// PFMRImageStream //
//////////////////////
template<class Real>
PFMRImageStream<Real>::PFMRImageStream( char* fileName , int& width , int& height , bool separateColors , bool gammaDecode , MultiStreamIOServer* ioServer ) :
ReadImageStream( PFMInitReadColor , PFMReadRow , PFMFinalizeRead , fileName , width , height , separateColors , gammaDecode , ioServer )
{
}
//////////////////////
// TIFFRImageStream //
//////////////////////
template< class Real >
TIFFRImageStream< Real >::TIFFRImageStream( char* fileName , int& width , int& height , bool separateColors , bool gammaDecode , MultiStreamIOServer* ioServer ) :
ReadImageStream( TIFFInitReadColor , TIFFReadRow , TIFFFinalizeRead , fileName , width , height , separateColors , gammaDecode , ioServer )
{
}

////////////////////////
// StreamingInputTile //
////////////////////////
template<class Real>
StreamingInputTile<Real>::StreamingInputTile( void )
{
	rc=1;
	data=NULL;
	w=h=0;
	sX=sY=0;
	rows=NULL;
}
template<class Real>
StreamingInputTile<Real>::~StreamingInputTile( void )
{
	if(data)	delete data;
	if(rows)
	{
		for(int i=0;i<rc;i++)	if(rows[i])	delete[] rows[i];
		delete[] rows;
	}
	data=NULL;
	rows=NULL;
	rc=0;
};

template<class Real>
void StreamingInputTile<Real>::Init( const char* tName , int rowCount , int startX , int startY , bool separateColors , bool gammaDecode )
{
	average[0]=average[1]=average[2]=0;
	strcpy(tileName,tName);
	rc=rowCount;
	sX=startX;
	sY=startY;
	if(rowCount<1)	exit(0);
	rows = new Real*[rc];
	data = GetReadStream< Real >( tileName , w , h , separateColors , gammaDecode , false , NULL );
	for(int i=-rc;i<0;i++)	init(i);
}

template<class Real>
void StreamingInputTile<Real>::init( int idx )
{
	if( idx+rc-1==sY )
	{
		if(!data)	exit(0);
		for(int i=0;i<rc;i++)	rows[i] = new Real[w*3];
		for(int i=0;i<rc-1;i++)
		{
			memcpy( rows[i] , (*data)[i] , w*sizeof(Real)*3 );
			data->advance();
		}
	}
	else if(idx==sY+h)
	{
		delete data;
		for(int i=0;i<rc;i++) delete[] rows[i];
		delete[] rows;
		data=NULL;
		rows=NULL;
	}
	if(idx>=sY && data)
	{
		memcpy( rows[(idx+rc-1-sY)%rc] , (*data)[idx+rc-1-sY] , w*3*sizeof(Real) );
		data->advance();
	}
}
template<class Real>
Color<Real> StreamingInputTile<Real>::operator() ( int i , int j )
{
	Color<Real> clr;
	if( this->data->SeparateColors( ) )
	{
		Real* data = &rows[(j-sY)%rc][i-sX];
		clr[0] = data[0*w];
		clr[1] = data[1*w];
		clr[2] = data[2*w];
	}
	else
	{
		Real* data = &rows[(j-sY)%rc][3*(i-sX)];
		clr[0] = data[0];
		clr[1] = data[1];
		clr[2] = data[2];
	}
	return clr;
}
template<class Real>
Color<Real> StreamingInputTile<Real>::operator() ( int i , int j , Real& a )
{
	Color<Real> clr;
	if( this->data->SeparateColors() )
	{
		Real* data=&rows[(j-sY)%rc][i-sX];
		clr[0] = data[0*w];
		clr[1] = data[1*w];
		clr[2] = data[2*w];
	}
	{
		Real* data=&rows[(j-sY)%rc][3*(i-sX)];
		clr[0] = data[0];
		clr[1] = data[1];
		clr[2] = data[2];
	}
	return clr;
}
template<class Real>
int StreamingInputTile<Real>::width(void)	const	{return w;}
template<class Real>
int StreamingInputTile<Real>::height(void)	const	{return h;}
template<class Real>
int StreamingInputTile<Real>::startX(void) const	{return sX;}
template<class Real>
int StreamingInputTile<Real>::startY(void) const	{return sY;}


///////////////////
// RGridOfImages //
///////////////////
template< class Real >
void RGridOfImages< Real >::GetImageSize( char* fileName , int& width , int& height )
{
	int c , r;
	FILE* fp = fopen(fileName,"r");
	char imageName[1024];
	if( fscanf( fp , "Columns: %d " , &c ) != 1) fprintf( stderr , "Failed to read in columns in RGridOfImages\n" ) , exit( 0 );
	if( fscanf( fp , "Rows: %d " , &r )     !=1) fprintf( stderr , "Failed to read in rows in RGridOfImages\n" )    , exit( 0 );
	width = height = 0;
	int* widths  = new int[c];
	int* heights = new int[r];
	for ( int j = 0 ; j < r ; j++ )
		for ( int i = 0 ; i < c ; i++ )
		{
			if( fscanf( fp , " %s " , imageName )!=1 )	fprintf(stderr,"Failed to read in image [%d][%d] in RGridOfImages\n",i,j) , exit(0);
			int ww , hh;
			GetReadSize< Real >( imageName , ww , hh );
			if( !j ) widths[i]  = ww	,	width += ww;
			if( !i ) heights[j] = hh	,	height += hh;
			if( widths[i]!=ww )		fprintf( stderr , "Inconsistent widths in column [%d]: %d != %d\n" , i ,  widths[i] , ww );
			if( heights[j]!=hh )	fprintf( stderr , "Inconsistent heights in row [%d]: %d != %d\n"   , j , heights[j] , hh );
		}
	fclose(fp);
	delete[] widths;
	delete[] heights;

}
template<class Real>
RGridOfImages<Real>::RGridOfImages( char* fileName , int& width , int& height , bool separateColors , bool gammaDecode , MultiStreamIOServer* ioServer)
{
	FILE* fp = fopen(fileName,"r");
	char imageName[1024];
	if(fscanf(fp,"Columns: %d ",&_c) != 1)	fprintf(stderr,"Failed to read in columns in RGridOfImages\n") , exit(0);
	if(fscanf(fp,"Rows: %d ",&_r) !=1)		fprintf(stderr,"Failed to read in rows in RGridOfImages\n")    , exit(0);

	int w,h;
	_gammaDecode = gammaDecode;
	_separate=separateColors;
	_widths  = new int[_c];
	_heights = new int[_r];
	_grid = new StreamingGrid**[_c];
	_fileNames = new char**[_c];
	for ( int i = 0 ; i < _c ; i++ )
	{
		_grid[i] = new StreamingGrid*[_r];
		_fileNames[i] = new char*[_r];
	}
	_w = _h = 0;
	for ( int j = 0 ; j < _r ; j++ )
		for ( int i = 0 ; i < _c ; i++ )
		{
			if(fscanf(fp," %s ",imageName)!=1)	fprintf(stderr,"Failed to read in image [%d][%d] in RGridOfImages\n",i,j) , exit(0);
			_fileNames[i][j]=new char[strlen(imageName)+1];
			strcpy( _fileNames[i][j] , imageName );
			GetReadSize< Real >( _fileNames[i][j] , w , h );
			if( !j ) _widths[i]  = w	,	_w += w;
			if( !i ) _heights[j] = h	,	_h += h;
			if( _widths[i]!=w )  fprintf( stderr , "Inconsistent widths in column [%d]: %d != %d\n" , i , _widths[i]  , w ); 
			if( _heights[j]!=h ) fprintf( stderr , "Inconsistent heights in row   [%d]: %d != %d\n" , j , _heights[j] , h );
		}
	fclose(fp);
	_buffer = AllocArray< Real >( _w * 3 , 1 , "RGridOfImages::_buffer" );
	_current = -1;
	advance();
	width = _w;
	height = _h;
#if NEW_JPEG_IO
	_server = ioServer;
#endif // NEW_JPEG_IO
}
template<class Real>
RGridOfImages<Real>::~RGridOfImages(void)
{
	delete[] _widths;
	delete[] _heights;
	_widths = _heights = NULL;
	for (int i = 0 ; i < _c ; i++ )
	{
		for ( int j = 0 ; j <_r ; j++ )
		{
			if(_grid[i][j])			delete _grid[i][j];
			if(_fileNames[i][j])	delete _fileNames[i][j];
		}
		delete[] _grid[i];
		delete[] _fileNames[i];
	}
	delete[] _grid;
	delete[] _fileNames;
	FreeArray( _buffer );
	_grid = NULL;
	_fileNames = NULL;
}
template<class Real>	int RGridOfImages<Real>::rows(void)		const	{	return _h;	}
template<class Real>	int RGridOfImages<Real>::rowSize(void)	const	{	return _w*3*sizeof(Real);	}
template<class Real>	bool RGridOfImages<Real>::SeparateColors(void)	const	{	return _separate;	}
template<class Real>
Pointer( byte ) RGridOfImages<Real>::operator [] (int idx)
{
#if ASSERT_MEMORY_ACCESS
	if(idx!=_current)	fprintf(stderr,"Index out of bounds: %d != %d\n",idx,_current) , exit(0);
#endif // ASSERT_MEMORY_ACCESS
	return ( Pointer( byte ) )( _buffer );
}
template<class Real>
void RGridOfImages<Real>::advance(void)
{
	_current++;
	if(_current<_h)
	{
		int row;
		int x = 0 , y = _current;
		for ( row = 0 ; y >= _heights[row] ; y -= _heights[row++] ) ;

		for( int col = 0 ; col < _c ; col++ )
		{
			if( y == 0 )
			{
				int w,h;
#if NEW_JPEG_IO
				_grid[col][row] = GetReadStream<Real>( _fileNames[col][row] , w , h , _separate , _gammaDecode , false , _server );
				_grid[col][row]->SetServer( _server );
#else // !NEW_JPEG_IO
				_grid[col][row] = GetReadStream<Real>( _fileNames[col][row] , w , h , _separate , _gammaDecode , false );
#endif // NEW_JPEG_IO
			}
			if( !_grid[col][row] )	fprintf(stderr,"Error: Attempting to read from NULL tile\n") , exit(0);
			Pointer( Real ) subRow = ( Pointer( Real ) )(*_grid[col][row])[y];

			if(_separate)
			{
				memcpy(_buffer+x+0*_w,subRow+0*_widths[col],sizeof(Real)*_widths[col]);
				memcpy(_buffer+x+1*_w,subRow+1*_widths[col],sizeof(Real)*_widths[col]);
				memcpy(_buffer+x+2*_w,subRow+2*_widths[col],sizeof(Real)*_widths[col]);
			}
			else
				memcpy(_buffer+3*x,subRow,sizeof(Real)*_widths[col]*3);

			x+=_widths[col];
			_grid[col][row]->advance();

			if( y >= _grid[col][row]->rows()-1 )
			{
				delete _grid[col][row];
				_grid[col][row]=NULL;
			}
		}
	}
}
///////////////////
// WGridOfImages //
///////////////////
template<class Real>
WGridOfImages<Real>::WGridOfImages(char* fileName,int width,int height,bool separateColors , bool gammaEncode , int quality,int tileWidth,int tileHeight , MultiStreamIOServer* ioServer , bool hdr )
{
//	const char* ext ="jpg";
	char* header = GetFileHeader( fileName );
//	_init( fileName , header , ext , width , height , separateColors , quality , tileWidth , tileHeight );
	_init( fileName , header , DefaultOutputTileExtension , width , height , separateColors , gammaEncode , quality , tileWidth , tileHeight , ioServer , hdr );
	delete[] header;
}
template<class Real>
WGridOfImages<Real>::WGridOfImages(char* fileName,const char* header,const char* ext,int width,int height,bool separateColors , bool gammaEncode , int quality,int tileWidth,int tileHeight , MultiStreamIOServer* ioServer , bool hdr )
{
	_init(fileName,header,ext,width,height,separateColors, gammaEncode , quality,tileWidth,tileHeight , ioServer , hdr );
}
template<class Real>
void WGridOfImages<Real>::_init(char* fileName,const char* header,const char* ext,int width,int height,bool separateColors , bool gammaEncode , int quality,int tileWidth,int tileHeight , MultiStreamIOServer* ioServer , bool hdr )
{
	FILE* fp = fopen(fileName,"w");
	_separate=separateColors;
	_gammaEncode = gammaEncode;
	_w = width;
	_h = height;
	_c = (width + tileWidth -1) / tileWidth;
	_r = (height + tileHeight -1) / tileHeight;
	_quality = quality;
	_hdr = hdr;
	char imageName[1024];
	fprintf(fp,"Columns: %d\n",_c) , fflush( fp );
	fprintf(fp,"Rows: %d\n",_r)    , fflush( fp );

	int w,h;
	_widths  = new int[_c];
	_heights = new int[_r];
	for ( int i = 0 ; i < _c ; i++ )	_widths[i]  = tileWidth;
	for ( int j = 0 ; j < _r ; j++ )	_heights[j] = tileHeight;
	_widths[_c-1]  = width  - (_c-1)*tileWidth;
	_heights[_r-1] = height - (_r-1)*tileHeight;

	_grid = new StreamingGrid**[_c];
	_fileNames = new char**[_c];
	for ( int i = 0 ; i < _c ; i++ )
	{
		_grid[i] = new StreamingGrid*[_r];
		_fileNames[i] = new char*[_r];
	}
	for ( int j = 0 ; j < _r ; j++ )
		for ( int i = 0 ; i < _c ; i++ )
		{
			sprintf( imageName , "%s.%d.%d.%s" , header , i , j , ext );
			_fileNames[i][j]=new char[strlen(imageName)+1];
			strcpy(_fileNames[i][j],imageName);
			fprintf(fp,"%s\n",imageName) , fflush( fp );
			_grid[i][j] = NULL;
		}
	fclose(fp);
	_buffer = AllocArray< Real >( _w * 3 , 1 , "WGridOfImages::_buffer" );
	_current = 0;
#if NEW_JPEG_IO
	_server = ioServer;
#endif // NEW_JPEG_IO
}
template<class Real>
WGridOfImages<Real>::~WGridOfImages(void)
{
	delete[] _widths;
	delete[] _heights;
	_widths = _heights = NULL;
	for (int i = 0 ; i < _c ; i++ )
	{
		for ( int j = 0 ; j <_r ; j++ )
		{
			if(_grid[i][j])			delete _grid[i][j];
			if(_fileNames[i][j])	delete _fileNames[i][j];
		}
		delete[] _grid[i];
	}
	delete[] _grid;
	FreeArray( _buffer );
	_grid = NULL;
	_fileNames = NULL;
}
template<class Real>	int WGridOfImages<Real>::rows(void)		const	{	return _h;	}
template<class Real>	int WGridOfImages<Real>::rowSize(void)	const	{	return _w*3*sizeof(Real);	}
template<class Real>	bool WGridOfImages<Real>::SeparateColors(void)	const	{	return _separate;	}
template<class Real>
Pointer( byte ) WGridOfImages<Real>::operator [] (int idx)
{
#if ASSERT_MEMORY_ACCESS
	if(idx!=_current)	fprintf(stderr,"Index out of bounds: %d != %d\n",idx,_current) , exit(0);
#endif // ASSERT_MEMORY_ACCESS
	return ( Pointer( byte ) )( _buffer );
}

template<class Real>
void WGridOfImages<Real>::advance(void)
{
	if(_current<_h)
	{
		int row;
		int x = 0 , y = _current;
		for ( row = 0 ; y >= _heights[row] ; y -= _heights[row++] ) ;

		for( int col = 0 ; col < _c ; col++ )
		{
			if ( y==0 )
			{
#if NEW_JPEG_IO
				_grid[col][row] = GetWriteStream<Real>( _fileNames[col][row] , _widths[col] , _heights[row] , _separate , _gammaEncode , _quality , _server , _hdr );
#else // !NEW_JPEG_IO
				_grid[col][row] = GetWriteStream<Real>( _fileNames[col][row] , _widths[col] , _heights[row] , _separate , _gammaEncode , _quality );
#endif // NEW_JPEG_IO
			}
			if ( !_grid[col][row])	fprintf(stderr,"Error: Attempting to write to NULL tile\n") , exit(0);
			Pointer( Real ) subRow = ( Pointer( Real ) )(*_grid[col][row])[y];
			if(_separate)
			{
				memcpy(subRow+0*_widths[col],_buffer+x+0*_w,sizeof(Real)*_widths[col]);
				memcpy(subRow+1*_widths[col],_buffer+x+1*_w,sizeof(Real)*_widths[col]);
				memcpy(subRow+2*_widths[col],_buffer+x+2*_w,sizeof(Real)*_widths[col]);
			}
			else
				memcpy(subRow,_buffer+3*x,sizeof(Real)*_widths[col]*3);

			x+=_widths[col];
			_grid[col][row]->advance();

			if ( y>=_grid[col][row]->rows()-1 )
			{
				delete _grid[col][row];
				_grid[col][row] = NULL;
			}
		}
	}
	_current++;
}



/////////////////////////////////////////////
// MultiThreadedRGridOfImages::ImageReader //
/////////////////////////////////////////////
template< class Real >
void MultiThreadedRGridOfImages< Real >::ImageReader::init( int columns , int rows ,
														    int* widths , int* heights ,
															char*** fileNames ,
															bool separateColors ,
															bool gammaDecode ,
															Pointer( Real ) buffer ,
															HANDLE loadHandle , HANDLE unloadHandle ,
															MultiStreamIOServer* ioServer )
{
if( separateColors )
{
	fprintf( stderr , "MultiThreadedRGridOfImages does not support color separation\n" );
	exit( 0 );
}

	_c				= columns;
	_r				= rows;
	_widths			= widths;
	_heights		= heights;
	_separate		= separateColors;
	_gammaDecode	= gammaDecode;
	_fileNames		= fileNames;
	_buffer			= buffer;
	_loadHandle		= loadHandle;
	_unloadHandle	= unloadHandle;
#if NEW_JPEG_IO
	_server = ioServer;
#endif // NEW_JPEG_IO
	_w = _h = 0;
	for( int i=0 ; i<_c ; i++ ) _w += _widths[i];
	for( int i=0 ; i<_r ; i++ ) _h += _heights[i];
}
template< class Real >
void MultiThreadedRGridOfImages< Real >::ImageReader::Run( void )
{
	StreamingGrid** grids = new StreamingGrid*[_c];
	for( int i=0 ; i<_h ; i++ )
	{
//		if( WaitForSingleObject( _loadHandle , INFINITE ) != WAIT_OBJECT_0 ) fprintfId( stderr , "Wait for single failed: " ) , PrintError();
		MyWaitForSingleObject( _loadHandle , 10000 , "MultiThreadedRGridOfImages<Real>::ImageReader::Run" );
		int row , x = 0 , y = i;
		for ( row = 0 ; y >= _heights[row] ; y -= _heights[row++] ) ;
		for( int col = 0 ; col < _c ; col++ )
		{
			if( y == 0 )
			{
				int w,h;
#if NEW_JPEG_IO
				grids[col] = GetReadStream< Real >(_fileNames[col][row] , w , h , _separate , _gammaDecode , false , _server );
#else // !NEW_JPEG_IO
				grids[col] = GetReadStream< Real >(_fileNames[col][row] , w , h , _separate , _gammaDecode , false );
#endif // NEW_JPEG_IO
			}
			if( !grids[col] )	fprintf(stderr,"Error: Attempting to read from NULL tile\n") , exit(0);
			Pointer( Real ) subRow = ( Pointer( Real ) )(*grids[col])[y];
			memcpy( _buffer+3*x , subRow , sizeof(Real)*_widths[col]*3 );

			x += _widths[col];
			grids[col]->advance();

			if( y >= grids[col]->rows()-1 )
			{
				delete grids[col];
				grids[col]=NULL;
			}
		}
		if( !SetEvent( _unloadHandle ) ) fprintfId( stderr , "Set event failed: " ) , PrintError();
	}
	delete[] grids;
}

////////////////////////////////
// MultiThreadedRGridOfImages //
////////////////////////////////
template< class Real >
void MultiThreadedRGridOfImages< Real >::GetImageSize( char* fileName , int& width , int& height )
{
	RGridOfImages< Real >::GetImageSize( fileName , width , height );
}
template<class Real>
MultiThreadedRGridOfImages<Real>::MultiThreadedRGridOfImages( char* fileName , int& width , int& height , bool separateColors , bool gammaDecode , int threads , int threadPriority , MultiStreamIOServer* ioServer )
{
	FILE* fp = fopen( fileName , "r" );
	char imageName[1024];
	if( fscanf( fp , "Columns: %d " , &_c ) !=1 )	fprintf( stderr , "Failed to read in columns in RGridOfImages\n" ) , exit(0);
	if( fscanf( fp , "Rows: %d " , &_r ) !=1 )		fprintf( stderr , "Failed to read in rows in RGridOfImages\n" )    , exit(0);

	int w,h;
	_threads	= threads;
	_gammaDecode= gammaDecode;
	_separate	= separateColors;
	_widths		= new int[_c];
	_heights	= new int[_r];
	_fileNames	= new char**[_c];
	_dataReady	= false;
	if( _threads > _c ) _threads = _c;
	for ( int i = 0 ; i < _c ; i++ ) _fileNames[i] = new char*[_r];

	_w = _h = 0;
	for ( int j = 0 ; j < _r ; j++ )
		for ( int i = 0 ; i < _c ; i++ )
		{
			if( fscanf( fp , " %s " , imageName )!=1 ) fprintf( stderr , "Failed to read in image [%d][%d] in RGridOfImages\n" , i , j ) , exit(0);
			_fileNames[i][j]=new char[strlen(imageName)+1];
			strcpy( _fileNames[i][j] , imageName );
			GetReadSize< Real >( _fileNames[i][j] , w , h );
			if( !j ) _widths[i]  = w	,	_w += w;
			if( !i ) _heights[j] = h	,	_h += h;
			if( _widths[i]!=w )  fprintf( stderr , "Inconsistent widths in column [%d]: %d != %d\n" , i , _widths[i]  , w );
			if( _heights[j]!=h ) fprintf( stderr , "Inconsistent heights in row   [%d]: %d != %d\n" , j , _heights[j] , h );
		}
	fclose(fp);
	_buffer = AllocArray< Real >( _w * 3 , 1 , "MultiThreadedRGridOfImages::_buffer" );

	_imageReaderHandles	= new HANDLE[ _threads ];
	_loadHandles		= new HANDLE[ _threads ];
	_unloadHandles		= new HANDLE[ _threads ];
	_imageReaders		= new ImageReader[ _threads ];
	int offset = 0;
	for( int i=0 ; i<_threads ; i++ )
	{
		_loadHandles[i]		= CreateEvent( NULL , false , false , NULL );
		_unloadHandles[i]	= CreateEvent( NULL , false , true  , NULL );
		int start = (i*_c) / _threads;
		int end = ( (i+1)*_c ) / _threads;
		_imageReaders[i].init( end-start , _r , _widths+start , _heights , _fileNames + start , _separate , _gammaDecode , _buffer+offset*3 , _loadHandles[i] , _unloadHandles[i] , ioServer );
		for( int j=start ; j<end ; j++ ) offset += _widths[j];
		_imageReaderHandles[i] = SpawnThread< ImageReader >( _imageReaders+i , threadPriority );
	}
	_current = -1;
	advance();
	width	= _w;
	height	= _h;
}
template<class Real>
MultiThreadedRGridOfImages<Real>::~MultiThreadedRGridOfImages( void )
{
	delete[] _widths;
	delete[] _heights;
	_widths = _heights = NULL;
	for (int i = 0 ; i < _c ; i++ )
	{
		for ( int j = 0 ; j <_r ; j++ ) if(_fileNames[i][j])	delete _fileNames[i][j];
		delete[] _fileNames[i];
	}
	delete[] _fileNames;
	FreeArray( _buffer );
	_fileNames = NULL;
//	WaitForMultipleObjects( _threads , _imageReaderHandles , TRUE , INFINITE );
	MyWaitForMultipleObjects( _threads , _imageReaderHandles , 10000 , "MultiThreadedRGridOfImages<Real>::~MultiThreadedRGridOfImages" );
	for( int i=0 ; i<_threads ; i++ )
	{
		CloseHandle( _loadHandles[i] );
		CloseHandle( _unloadHandles[i] );
	}
	delete[] _imageReaders;
	delete[] _imageReaderHandles;
	delete[] _loadHandles;
	delete[] _unloadHandles;
}
template<class Real> int MultiThreadedRGridOfImages<Real>::rows(void)	 const	{ return _h;	}
template<class Real> int MultiThreadedRGridOfImages<Real>::rowSize(void) const	{ return _w*3*sizeof(Real);	}
template<class Real> bool MultiThreadedRGridOfImages<Real>::SeparateColors(void) const	{ return _separate;	}
template<class Real>
Pointer( byte ) MultiThreadedRGridOfImages<Real>::operator [] (int idx)
{
	if( !_dataReady )
//		if( WaitForMultipleObjects( _threads , _unloadHandles , TRUE , INFINITE ) != WAIT_OBJECT_0 )
//			fprintfId( stderr , "Wait for multiple failed: " ) , PrintError();
		MyWaitForMultipleObjects( _threads , _unloadHandles , 10000 , "MultiThreadedRGridOfImages<Real>::operator[]" );
	_dataReady = true;
#if ASSERT_MEMORY_ACCESS
	if( idx!=_current ) fprintfId( stderr , "Index out of bounds: %d != %d\n" , idx , _c ) , exit(0);
#endif // ASSERT_MEMORY_ACCESS
	return ( Pointer( byte ) )( _buffer );
}
template<class Real>
void MultiThreadedRGridOfImages<Real>::advance( void )
{
	_current++;
	if( _current<_h )
	{
		if( !_dataReady )
//			if( WaitForMultipleObjects( _threads , _unloadHandles , TRUE , INFINITE ) != WAIT_OBJECT_0 )
//				fprintfId( stderr , "Wait for multiple failed: " ) , PrintError();
			MyWaitForMultipleObjects( _threads , _unloadHandles , 10000 , "MultiThreadedRGridOfImages<Real>::advance" );
		for( int i=0 ; i<_threads ; i++ ) if( !SetEvent( _loadHandles[i] ) ) fprintfId( stderr , "SetEvent failed: " ) , PrintError();
		_dataReady = false;
	}
}

/////////////////////////////////////////////
// MultiThreadedWGridOfImages::ImageWriter //
/////////////////////////////////////////////
template< class Real >
void MultiThreadedWGridOfImages< Real >::ImageWriter::init( int columns , int rows ,
														    int* widths , int* heights ,
															char*** fileNames ,
															bool separateColors , bool gammaEncode , int quality ,
															Pointer( Real ) buffer ,
															HANDLE loadHandle , HANDLE unloadHandle ,
															MultiStreamIOServer* ioServer , bool hdr )
{
if( separateColors )
{
	fprintf( stderr , "MultiThreadedWGridOfImages does not support color separation\n" );
	exit( 0 );
}

	_c				= columns;
	_r				= rows;
	_widths			= widths;
	_heights		= heights;
	_separate		= separateColors;
	_gammaEncode	= gammaEncode;
	_quality		= quality;
	_fileNames		= fileNames;
	_buffer			= buffer;
	_loadHandle		= loadHandle;
	_unloadHandle	= unloadHandle;
	_hdr			= hdr;
#if NEW_JPEG_IO
	_server = ioServer;
#endif // NEW_JPEG_IO
	_w = _h = 0;
	for( int i=0 ; i<_c ; i++ ) _w += _widths[i];
	for( int i=0 ; i<_r ; i++ ) _h += _heights[i];
}
template< class Real >
void MultiThreadedWGridOfImages< Real >::ImageWriter::Run( void )
{
	StreamingGrid** grids = new StreamingGrid*[_c];
	for( int i=0 ; i<_h ; i++ )
	{
#if 1
//		if( WaitForSingleObject( _unloadHandle , INFINITE ) != WAIT_OBJECT_0 ) fprintfId( stderr , "Wait for single failed: " ) , PrintError();
		MyWaitForSingleObject( _unloadHandle , 10000 , "MultiThreadedWGridOfImages<Real>::ImageWriter::Run" );
#else
		while( WaitForSingleObject( _unloadHandle , 10000 ) != WAIT_OBJECT_0 )
		{
			MyWinSock::StdinLock lock;
			fprintfId( stderr , "Deadlock in ImageWriter::Run[%d / %d]: " , i , _h ) , PrintError();
		}
#endif
		int row , x = 0 , y = i;
		for ( row = 0 ; y >= _heights[row] ; y -= _heights[row++] ) ;
		for( int col = 0 ; col < _c ; col++ )
		{

			if ( y==0 )
			{
				MyWinSock::StdinLock lock;
#if NEW_JPEG_IO
				grids[col] = GetWriteStream< Real >( _fileNames[col][row] , _widths[col] , _heights[row] , _separate , _gammaEncode , _quality , _server , _hdr );
#else // !NEW_JPEG_IO
				grids[col] = GetWriteStream< Real >( _fileNames[col][row] , _widths[col] , _heights[row] , _separate , _gammaEncode , _quality );
#endif // NEW_JPEG_IO
			}
			if( !grids[col] )	fprintf(stderr,"Error: Attempting to write to NULL tile\n") , exit(0);
			Pointer( Real ) subRow = ( Pointer( Real ) )(*grids[col])[y];

			memcpy( subRow ,  _buffer+3*x , sizeof(Real)*_widths[col]*3 );

			x += _widths[col];
			grids[col]->advance();

			if( y >= grids[col]->rows()-1 ) delete grids[col] , grids[col]=NULL;
		}
		if( !SetEvent( _loadHandle ) ) fprintfId( stderr , "Set event failed: " ) , PrintError();
	}
	delete[] grids;
}

////////////////////////////////
// MultiThreadedWGridOfImages //
////////////////////////////////
template<class Real>
MultiThreadedWGridOfImages<Real>::MultiThreadedWGridOfImages( char* fileName , int width , int height , bool separateColors , bool gammaEncode , int quality , int threads , int threadPriority , int tileWidth , int tileHeight , MultiStreamIOServer* ioServer , bool hdr )
{
	const char* ext = "jpg";
	char* header = GetFileHeader( fileName );
	_init( fileName , header , ext , width , height , separateColors , gammaEncode , quality , threads , threadPriority , tileWidth , tileHeight , ioServer , hdr );
	delete[] header;
}
template<class Real>
MultiThreadedWGridOfImages<Real>::MultiThreadedWGridOfImages( char* fileName , const char* header , const char* ext , int width , int height , bool separateColors , bool gammaEncode ,  int quality , int threads , int threadPriority , int tileWidth , int  tileHeight , MultiStreamIOServer* ioServer , bool hdr )
{
	_init( fileName , header , ext , width , height , separateColors , gammaEncode , quality , threads , threadPriority , tileWidth , tileHeight , ioServer , hdr );
}
template<class Real>
void MultiThreadedWGridOfImages<Real>::_init( char* fileName , const char* header , const char* ext , int width , int height , bool separateColors , bool gammaEncode , int quality , int threads , int threadPriority , int tileWidth , int tileHeight , MultiStreamIOServer* ioServer , bool hdr )
{
	FILE* fp = fopen(fileName,"w");
	_separate	= separateColors;
	_gammaEncode= gammaEncode;
	_w			= width;
	_h			= height;
	_c			= (width + tileWidth -1) / tileWidth;
	_r			= (height + tileHeight -1) / tileHeight;
	_quality	= quality;
	_threads	= threads;
	_hdr		= hdr;
	_dataReady	= false;
	if( _threads>_c ) _threads = _c;
	char imageName[1024];
	fprintf( fp , "Columns: %d\n" , _c ) , fflush( fp );
	fprintf( fp , "Rows: %d\n" , _r )    , fflush( fp );

	int w , h;
	_widths  = new int[_c];
	_heights = new int[_r];
	for ( int i = 0 ; i < _c ; i++ )	_widths[i]  = tileWidth;
	for ( int j = 0 ; j < _r ; j++ )	_heights[j] = tileHeight;
	_widths[_c-1]  = width  - (_c-1)*tileWidth;
	_heights[_r-1] = height - (_r-1)*tileHeight;

	_fileNames = new char**[_c];
	for ( int i = 0 ; i < _c ; i++ ) _fileNames[i] = new char*[_r];
	for ( int j = 0 ; j < _r ; j++ )
		for ( int i = 0 ; i < _c ; i++ )
		{
			sprintf(imageName,"%s.%d.%d.%s",header,i,j,ext);
			_fileNames[i][j] = new char[strlen(imageName)+1];
			strcpy( _fileNames[i][j] , imageName );
			fprintf( fp , "%s\n" , imageName ) , fflush( fp );
		}
	fclose(fp);
	_buffer = AllocArray< Real >( _w * 3 , 1 , "MultiThreadedWGridOfImages::_buffer" );
	_imageWriterHandles	= new HANDLE[ _threads ];
	_loadHandles		= new HANDLE[ _threads ];
	_unloadHandles		= new HANDLE[ _threads ];
	_imageWriters		= new ImageWriter[ _threads ];
	int offset = 0;
	for( int i=0 ; i<_threads ; i++ )
	{
		_loadHandles[i]		= CreateEvent( NULL , false , true  , NULL );
		_unloadHandles[i]	= CreateEvent( NULL , false , false , NULL );
		int start = (i*_c) / _threads;
		int end = ( (i+1)*_c ) / _threads;
		_imageWriters[i].init( end-start , _r , _widths+start , _heights , _fileNames + start , _separate , _gammaEncode , quality , _buffer+offset*3 , _loadHandles[i] , _unloadHandles[i] , ioServer , hdr );
		for( int j=start ; j<end ; j++ ) offset += _widths[j];
		_imageWriterHandles[i] = SpawnThread< ImageWriter >( _imageWriters+i , threadPriority );
	}
	_current = 0;
}
template<class Real>
MultiThreadedWGridOfImages<Real>::~MultiThreadedWGridOfImages(void)
{
	delete[] _widths;
	delete[] _heights;
	_widths = _heights = NULL;
	for (int i = 0 ; i < _c ; i++ )
	{
		for ( int j = 0 ; j <_r ; j++ ) if( _fileNames[i][j] ) delete _fileNames[i][j];
		delete[] _fileNames[i];
	}
	delete[] _fileNames;
	FreeArray( _buffer );
	_fileNames = NULL;
//	WaitForMultipleObjects( _threads , _imageWriterHandles , TRUE , INFINITE );
//	if( WaitForMultipleObjects( _threads , _imageWriterHandles , TRUE , 10000 ) != WAIT_OBJECT_0 ) fprintfId( stderr , "Dead Lock in ~MultiThreadedWGridOfImages: " ) , PrintError();
	MyWaitForMultipleObjects( _threads , _imageWriterHandles , 10000 , "MultiThreadedWGridOfImages<Real>::~MultiThreadedWGridOfImages" );
	for( int i=0 ; i<_threads ; i++ )
	{
		CloseHandle( _loadHandles[i] );
		CloseHandle( _unloadHandles[i] );
	}
	delete[] _imageWriters;
	delete[] _imageWriterHandles;
	delete[] _loadHandles;
	delete[] _unloadHandles;
}
template<class Real> int MultiThreadedWGridOfImages<Real>::rows(void)	 const	{	return _h;	}
template<class Real> int MultiThreadedWGridOfImages<Real>::rowSize(void) const	{	return _w*3*sizeof(Real);	}
template<class Real> bool MultiThreadedWGridOfImages<Real>::SeparateColors(void)	const	{	return _separate;	}
template<class Real>
Pointer( byte ) MultiThreadedWGridOfImages<Real>::operator [] (int idx)
{
	if( !_dataReady )
	{
#if 1
//		if( WaitForMultipleObjects( _threads , _loadHandles , true , INFINITE ) != WAIT_OBJECT_0 )
//			fprintfId( stderr , "Wait for multiple failed: " ) , PrintError();
		MyWaitForMultipleObjects( _threads , _loadHandles , 10000 , "MultiThreadedWGridOfImages<Real>::operator[]" );
#else
		if( WaitForMultipleObjects( _threads , _loadHandles , TRUE , 10000 ) != WAIT_OBJECT_0 )
			fprintfId( stderr , "Deadlock in MultiThreadedWGridOfImages[]: " ) , PrintError();
#endif
	}
	_dataReady = true;
#if ASSERT_MEMORY_ACCESS
	if( idx!=_current )	fprintf(stderr,"Index out of bounds: %d != %d\n",idx,_current) , exit(0);
#endif // ASSERT_MEMORY_ACCESS
	return ( Pointer( byte ) )( _buffer );
}
template<class Real>
void MultiThreadedWGridOfImages<Real>::advance(void)
{
	if( _current<_h )
	{
		if( !_dataReady )
		{
#if 1
//			if( WaitForMultipleObjects( _threads , _loadHandles , TRUE , INFINITE ) != WAIT_OBJECT_0 )
//				fprintfId( stderr , "Wait for multiple failed: " ) , PrintError();
			MyWaitForMultipleObjects( _threads , _loadHandles , 10000 , "MultiThreadedWGridOfImages<Real>::advance" );
#else
			if( WaitForMultipleObjects( _threads , _loadHandles , TRUE , 10000 ) != WAIT_OBJECT_0 )
				fprintfId( stderr , "Deadlock in MultiThreadedWGridOfImages::advance: " ) , PrintError();
#endif
		}
		for( int i=0 ; i<_threads ; i++ ) if( !SetEvent( _unloadHandles[i] ) ) fprintfId( stderr , "SetEvent failed: " ) , PrintError();
		_dataReady = false;
	}
	_current++;
}

///////////////////
// WImageClipper //
///////////////////
template<class Real>
WImageClipper<Real>::WImageClipper( char* fileName , int width , int height , bool separateColors , bool gammaEncode , int outOffX , int outOffY , int outWidth , int outHeight , int quality , MultiStreamIOServer* ioServer , bool hdr , int threads , int threadPriority )
{
	deleteGrid = true;
	_separate=separateColors;
	_hdr = hdr;
	_w = width;
	_h = height;
	_ox = outOffX;
	_oy = outOffY;
	_ow = outWidth;
	_oh = outHeight;

	if(_ox<0 || _oy<0 || _ox+_ow>_w || _oy+_oh>_h)
	{
		fprintf(stderr,"Clip-region not a subset of the domain in WImageClipper: [%d,%d) x [%d,%d) !< [0,%d) x [0,%d)\n",_ox,_ox+_ow,_oy,_oy+_oh,_w,_h);
		exit(0);
	}

	grid = GetWriteStream<Real>( fileName , _ow , _oh , false , gammaEncode , quality , ioServer , hdr , threads , threadPriority );
	buffer = AllocArray< Real >( _w * 3 , 1 , "WImageClipper::buffer" );
	current = 0;
}
template<class Real>
WImageClipper<Real>::WImageClipper(StreamingGrid* outGrid,int width,int height,bool separateColors , int outOffX,int outOffY,int outWidth,int outHeight , bool hdr )
{
	deleteGrid = false;
	_separate=separateColors;
	_hdr = hdr;
	_w = width;
	_h = height;
	_ox = outOffX;
	_oy = outOffY;
	_ow = outWidth;
	_oh = outHeight;

	if(_ox<0 || _oy<0 || _ox+_ow>_w || _oy+_oh>_h)
	{
		fprintf(stderr,"Clip-region not a subset of the domain in WImageClipper: [%d,%d) x [%d,%d) !< [0,%d) x [0,%d)\n",_ox,_ox+_ow,_oy,_oy+_oh,_w,_h);
		exit(0);
	}

	grid = outGrid;
	buffer = AllocArray< Real >( _w * 3 , 1 , "WImageClipper::buffer" );
	current = 0;
}
template<class Real>
WImageClipper<Real>::~WImageClipper(void)
{
	if( deleteGrid ) delete grid;
	FreeArray( buffer );
}
template<class Real> int WImageClipper<Real>::rows(void)	const { return _h;}
template<class Real> int WImageClipper<Real>::rowSize(void)	const { return sizeof(Real)*_w*3;}
template<class Real> bool WImageClipper<Real>::SeparateColors(void)	const { return _separate; }

template<class Real>
Pointer( byte ) WImageClipper<Real>::operator [] (int idx)
{
#if ASSERT_MEMORY_ACCESS
	if(idx!=current)	fprintf(stderr,"Index out of bounds: %d != %d\n",idx,current) , exit(0);
#endif // ASSERT_MEMORY_ACCESS
	return ( Pointer( byte ) )( buffer );
}
template<class Real>
void WImageClipper<Real>::advance(void)
{
	if(current>=_oy && current<_oy+_oh)
	{
		Pointer( Real ) row = ( Pointer( Real ) )(*grid)[current-_oy];

		if(_separate)
		{
#if 1
			for( int c=0 ; c<3 ;c++ )
				for( int i=0 ; i<_ow ; i++ )
					row[3*i+c] = buffer[_ox + i + c*_w];
#else
			memcpy(row+0*_ow,buffer+_ox+0*_w,sizeof(Real)*_ow);
			memcpy(row+1*_ow,buffer+_ox+1*_w,sizeof(Real)*_ow);
			memcpy(row+2*_ow,buffer+_ox+2*_w,sizeof(Real)*_ow);
#endif
		}
		else
			memcpy(row,buffer+3*_ox,sizeof(Real)*_ow*3);
		grid->advance();
	}
	current++;
}
///////////////////////
// RChannelExtractor //
///////////////////////
template<class Real>
RChannelExtractor<Real>::RChannelExtractor( char* fileName , int& width , int& height , bool separateColors , bool gammaDecode , int channel )
{
	_c = channel;
	grid = GetReadStream<Real>( fileName , _w , _h , separateColors , gammaDecode , false );
	width = _w;
	height = _h;

	buffer = new Real[_w];
	current = -1;
	advance();
}
template<class Real>
RChannelExtractor<Real>::~RChannelExtractor(void)
{
	delete[] buffer;
}
template<class Real>	int RChannelExtractor<Real>::rows(void)		const	{	return _h;	}
template<class Real>	int RChannelExtractor<Real>::rowSize(void)	const	{	return _w*sizeof(Real);	}
template<class Real>	bool RChannelExtractor<Real>::SeparateColors(void)	const	{	return _separate;	}
template<class Real>
Pointer( byte ) RChannelExtractor<Real>::operator [] (int idx)
{
#if ASSERT_MEMORY_ACCESS
	if(idx!=current)	fprintf(stderr,"Index out of bounds: %d != %d\n",idx,current) , exit(0);
#endif // ASSERT_MEMORY_ACCESS
	return ( Pointer( byte ) )( buffer , rowSize() );
}

template<class Real>
void RChannelExtractor<Real>::advance(void)
{
	current++;
	if(current<_h)
	{
		Real* row = (Real*)(*grid)[current];
		if(_separate)	memcpy(buffer,row+_c*_w,sizeof(Real)*_w);
		else	for(int x=0 ; x<_w ; x++)	buffer[x]=row[3*x+_c];
		grid->advance();
	}
}

///////////////////
// RImageSampler //
///////////////////
template< class Real >
RImageSampler< Real >::RImageSampler( char* fileName , int width , int  height , bool nearest , bool separateColors , bool gammaDecode , Real logScale , MultiStreamIOServer* ioServer , int threads , int threadPriority )
{
	if( logScale>0 ) _grid = new LogImageReader< Real >( fileName , _inW , _inH , logScale , separateColors , gammaDecode , false , ioServer , threads , threadPriority );
	else             _grid =      GetReadStream< Real >( fileName , _inW , _inH ,            separateColors , gammaDecode , false , ioServer , threads , threadPriority );
	_outW = width;
	_outH = height;
	_separate = separateColors;
	_nearest = nearest;

	_outBuffer = AllocArray< Real >( 3*_outW , 1 );
	if( !nearest )
	{
		_inBuffers[0] = AllocArray< Real >( 3*_inW , 1 );
		_inBuffers[1] = AllocArray< Real >( 3*_inW , 1 );
	}
	else _inBuffers[0] = _inBuffers[1] = NullPointer< Real >( );
	_inCurrent = 0 , _outCurrent = -1;
	advance();
}
template< class Real >
RImageSampler<Real>::~RImageSampler( void )
{
	FreeArray( _outBuffer );
	FreeArray( _inBuffers[0] );
	FreeArray( _inBuffers[1] );
	delete _grid;
}
template< class Real > int  RImageSampler< Real >::rows( void ) const { return _outH; }
template< class Real > int  RImageSampler< Real >::rowSize(void) const { return _outW * sizeof( Real ) * 3; }
template< class Real > bool RImageSampler< Real >::SeparateColors( void ) const { return _separate; }
template< class Real >
Pointer( byte ) RImageSampler< Real >::operator [] ( int idx )
{
#if ASSERT_MEMORY_ACCESS
	if( idx!=_outCurrent ) fprintf( stderr , "Index out of bounds: %d != %d\n" , idx , _outCurrent ) , exit( 0 );
#endif // ASSERT_MEMORY_ACCESS
	return ( Pointer( byte ) ) _outBuffer;
}
template< class Real >
void RImageSampler< Real >::advance( void )
{
	_outCurrent++;
	double y = double( long long( _outCurrent ) * long long( _inH-1 ) ) / ( _outH-1 );
	int yy;
	if( _nearest )	// Align the input grid to the row we will be reading from
	{
		yy = int( y+0.5 );
		while( _inCurrent<yy && _inCurrent<_inH ) _grid->advance( ) , _inCurrent++;
	}
	else	// Align the input grid to the row we will be reading from and read in the row into the appropriate buffer
	{
		// This needs to be fixed!!!
		yy = int( y );
		while( _inCurrent<yy && _inCurrent<_inH ) _grid->advance( ) , _inCurrent++;
		if( _inCurrent==yy && _inCurrent<_inH )
		{
			Pointer( Real ) row = ( Pointer( Real ) )(*_grid)[_inCurrent];
			memcpy( _inBuffers[_inCurrent&1] , row , sizeof( Real ) * _inW * 3 );
			_grid->advance( );
			_inCurrent++;
		}
		if( _inCurrent==yy+1 && _inCurrent<_inH )
		{
			Pointer( Real ) row = ( Pointer( Real ) )(*_grid)[_inCurrent];
			memcpy( _inBuffers[_inCurrent&1] , row , sizeof( Real ) * _inW * 3 );
			_grid->advance( );
			_inCurrent++;
		}
	}
	if( _outCurrent<_outH )
	{
		if( _nearest )
		{
			Pointer( Real ) row = ( Pointer( Real ) )(*_grid)[_inCurrent];
			for( int i=0 ; i<_outW ; i++ )
			{
//				float x = float( i * ( _inW-1 ) ) / ( _outW-1 );
				double x = float( long long( i ) * long long( _inW-1 ) ) / ( _outW-1 );
				int xx = int( x+0.5 );
				if( _separate ) for( int c=0 ; c<3 ; c++ ) _outBuffer[i+_outW*c] = row[xx+_inW*c];
				else            for( int c=0 ; c<3 ; c++ ) _outBuffer[i*3+c] = row[3*xx+c];
			}
		}
		else
		{
			int y1 = yy;
			int y2 = y1+1;
//			float dy = y-yy;
			double dy = y-yy;
			for( int i=0 ; i<_outW ; i++ )
			{
//				float x = float( i * ( _inW-1 ) ) / ( _outW-1 );
				double x = double( long long( i ) * long long( _inW-1 ) ) / ( _outW-1 );
				int x1 = int( x );
				int x2 = x1 + 1;
//				float dx = x-x1;
				double dx = x-x1;
				if( x2>=_inW ) x2 = x1;
				if( _separate )
					for( int c=0 ; c<3 ; c++ )
						_outBuffer[i+_outW*c] =
						_inBuffers[y1&1][x1+_inW*c] * Real(1.0-dx) * Real(1.0-dy) +
						_inBuffers[y1&1][x2+_inW*c] * Real(    dx) * Real(1.0-dy) +
						_inBuffers[y2&1][x1+_inW*c] * Real(1.0-dx) * Real(    dy) +
						_inBuffers[y2&1][x2+_inW*c] * Real(    dx) * Real(    dy) ;
				else
					for( int c=0 ; c<3 ; c++ )
						_outBuffer[i*3+c] =
						_inBuffers[y1&1][x1*3+c] * Real(1.0-dx) * Real(1.0-dy) +
						_inBuffers[y1&1][x2*3+c] * Real(    dx) * Real(1.0-dy) +
						_inBuffers[y2&1][x1*3+c] * Real(1.0-dx) * Real(    dy) +
						_inBuffers[y2&1][x2*3+c] * Real(    dx) * Real(    dy) ;
			}
		}
	}
}
///////////////////
// WImageSampler //
///////////////////
template< class Real >
WImageSampler< Real >::WImageSampler( char* fileName , int inWidth , int  inHeight , int outWidth , int outHeight , bool nearest , bool separateColors , bool logScale , bool gammaEncode , const int& quality , MultiStreamIOServer* ioServer , int threads , int threadPriority )
{
	if( logScale ) _grid = new LogImageWriter< Real >( fileName , outWidth , outHeight , separateColors , gammaEncode , quality , ioServer , false , threads , threadPriority );
	else           _grid =     GetWriteStream< Real >( fileName , outWidth , outHeight , separateColors , gammaEncode , quality , ioServer ,         threads , threadPriority );
	_inW = inWidth;
	_inH = inHeight;
	_outW = outWidth;
	_outH = outHeight;
	_separate = separateColors;
	_nearest = nearest;

	_inBuffers[0] = AllocArray< Real >( 3* _inW , 1 , "WImageSampler::_inBuffers[0]" );
	_inBuffers[1] = AllocArray< Real >( 3* _inW , 1 , "WImageSampler::_inBuffers[1]" );

	_inCurrent = 0 , _outCurrent = 0;
}
template< class Real >
WImageSampler<Real>::~WImageSampler( void )
{
	FreeArray( _inBuffers[0] );
	FreeArray( _inBuffers[1] );
	delete _grid;
}
template< class Real > int  WImageSampler< Real >::rows( void ) const { return _inH; }
template< class Real > int  WImageSampler< Real >::rowSize(void) const { return _inW * sizeof( Real ) * 3; }
template< class Real > bool WImageSampler< Real >::SeparateColors( void ) const { return _separate; }
template< class Real >
Pointer( byte ) WImageSampler< Real >::operator [] ( int idx )
{
#if ASSERT_MEMORY_ACCESS
	if( idx!=_inCurrent ) fprintf( stderr , "Index out of bounds: %d != %d\n" , idx , _inCurrent ) , exit( 0 );
#endif // ASSERT_MEMORY_ACCESS
	return ( Pointer( byte ) ) _inBuffers[idx&1];
}
template< class Real >
void WImageSampler< Real >::advance( void )
{
	long long start = (_inCurrent-1) * (_outH-1 ) , end = _inCurrent * (_outH-1);

	while( long long( _outCurrent ) * long long( _inH-1 )>=start && long long( _outCurrent ) * long long( _inH-1 )<=end )
	{
		double y = double( long long( _outCurrent ) * long long( _inH-1 ) ) / ( _outH-1 );
		if( _nearest )
		{
			int yy = int( y+0.5 );
			Pointer( Real )  inRow = _inBuffers[yy%2];
			Pointer( Real ) outRow = ( Pointer( Real ) )(*_grid)[_outCurrent];
			for( int i=0 ; i<_outW ; i++ )
			{
				double x = float( long long( i ) * long long( _inW-1 ) ) / ( _outW-1 );
				int xx = int( x+0.5 );
				if( _separate ) for( int c=0 ; c<3 ; c++ ) outRow[i+_outW*c] = inRow[xx+_inW*c];
				else            for( int c=0 ; c<3 ; c++ ) outRow[i*3+c]     = inRow[3*xx+c];
			}
		}
		else
		{
			Pointer( Real ) inRow1 = _inBuffers[(_inCurrent-1)&1];
			Pointer( Real ) inRow2 = _inBuffers[(_inCurrent  )&1];
			Pointer( Real ) outRow = ( Pointer( Real ) )(*_grid)[_outCurrent];
			double dy = y - (_inCurrent-1);
			for( int i=0 ; i<_outW ; i++ )
			{
				double x = float( long long( i ) * long long( _inW-1 ) ) / ( _outW-1 );
				double dx = x - int(x);
				int x1 = int(x);
				int x2 = x1 + 1;
				if( x2>=_inW ) x2 = _inW-1;
				double value[3];
				if( _separate ) for( int c=0 ; c<3 ; c++ ) value[c] = inRow1[x1+_inW*c]*(1.-dx)*(1.-dy) + inRow1[x2+_inW*c]*(dx)*(1.-dy) + inRow2[x1+_inW*c]*(1.-dx)*(dy) + inRow2[x2+_inW*c]*(dx)*(dy);
				else            for( int c=0 ; c<3 ; c++ ) value[c] = inRow1[3*x1+c   ]*(1.-dx)*(1.-dy) + inRow1[3*x2+c   ]*(dx)*(1.-dy) + inRow2[3*x1+c   ]*(1.-dx)*(dy) + inRow2[3*x2+c   ]*(dx)*(dy);
				if( _separate ) for( int c=0 ; c<3 ; c++ ) outRow[i+_outW*c] = value[c];
				else            for( int c=0 ; c<3 ; c++ ) outRow[i*3+c]     = value[c];
			}
		}
		_outCurrent++;
		_grid->advance( );
	}
	_inCurrent++;
}

////////////////////
// RImageSampler2 //
////////////////////
template< class Real >
RImageSampler2< Real >::RImageSampler2( char* nearestFileName , char* averageFileName , int width , int  height , bool separateColors , bool gammaDecode , MultiStreamIOServer* ioServer , int threads , int threadPriority )
{
	int w1 , w2 , h1 , h2;
	_nearestGrid = GetReadStream< Real >( nearestFileName , w1 , h1 , separateColors , gammaDecode , false , ioServer , threads , threadPriority );
	_averageGrid = GetReadStream< Real >( averageFileName , w2 , h2 , separateColors , gammaDecode , false , ioServer , threads , threadPriority );
	
	if( !_nearestGrid ) fprintf( stderr , "Could not read image: %s\n" , nearestFileName ) , exit( 0 );
	if( !_nearestGrid ) fprintf( stderr , "Could not read image: %s\n" , averageFileName ) , exit( 0 );
	if( w1!=w2 || h1!=h2 ) fprintf( stderr , "Dimensions differ: (%d x %d) != (%d x %d)\n" , w1 , h1 , w2 , h2 ) , exit( 0 );
	_inW = w1 , _inH = h1;

	_outW = width;
	_outH = height;
	_separate = separateColors;

	_outBuffer = AllocArray< Real >( 6*_outW );
	_nearestBuffers[0] = AllocArray< Real >( 3*_inW );
	_nearestBuffers[1] = AllocArray< Real >( 3*_inW );
	_averageBuffers[0] = AllocArray< Real >( 3*_inW );
	_averageBuffers[1] = AllocArray< Real >( 3*_inW );

	_inCurrent = 0 , _outCurrent = -1;
	advance();
}
template< class Real >
RImageSampler2< Real >::~RImageSampler2( void )
{
	FreeArray( _outBuffer );
	FreeArray( _nearestBuffers[0] );
	FreeArray( _nearestBuffers[1] );
	FreeArray( _averageBuffers[0] );
	FreeArray( _averageBuffers[1] );
	delete _nearestGrid , _averageGrid;
}
template< class Real > int  RImageSampler2< Real >::rows( void ) const { return _outH; }
template< class Real > int  RImageSampler2< Real >::rowSize(void) const { return _outW * sizeof( Real ) * 6; }
template< class Real > bool RImageSampler2< Real >::SeparateColors( void ) const { return _separate; }
template< class Real >
Pointer( byte ) RImageSampler2< Real >::operator [] ( int idx )
{
#if ASSERT_MEMORY_ACCESS
	if( idx!=_outCurrent ) fprintf( stderr , "Index out of bounds: %d != %d\n" , idx , _outCurrent ) , exit( 0 );
#endif // ASSERT_MEMORY_ACCESS
	return ( Pointer( byte ) ) _outBuffer;
}
template< class Real >
void RImageSampler2< Real >::advance( void )
{
	_outCurrent++;
	Pointer( Real ) row;
	double y = double( long long( _outCurrent ) * long long( _inH-1 ) ) / ( _outH-1 );
	int averageY , nearestY;

	nearestY = int( y+0.5 );
	averageY = int( y );
	while( _inCurrent<averageY && _inCurrent<_inH ) _nearestGrid->advance( ) , _averageGrid->advance( ) , _inCurrent++;
	if( _inCurrent==averageY && _inCurrent<_inH )
	{
		row = ( Pointer( Real ) )(*_nearestGrid)[_inCurrent];
		memcpy( _nearestBuffers[_inCurrent&1] , row , sizeof( Real ) * _inW * 3 );
		_nearestGrid->advance( );

		row = ( Pointer( Real ) )(*_averageGrid)[_inCurrent];
		memcpy( _averageBuffers[_inCurrent&1] , row , sizeof( Real ) * _inW * 3 );
		_averageGrid->advance( );

		_inCurrent++;
	}
	if( _inCurrent==averageY+1 && _inCurrent<_inH )
	{
		row = ( Pointer( Real ) )(*_nearestGrid)[_inCurrent];
		memcpy( _nearestBuffers[_inCurrent&1] , row , sizeof( Real ) * _inW * 3 );
		_nearestGrid->advance( );

		row = ( Pointer( Real ) )(*_averageGrid)[_inCurrent];
		memcpy( _averageBuffers[_inCurrent&1] , row , sizeof( Real ) * _inW * 3 );
		_averageGrid->advance( );

		_inCurrent++;
	}

	if( _outCurrent<_outH )
	{
		Real nearest[3] , average[3] , temp[3] , sum;
		int y1 = averageY;
		int y2 = y1+1;
		double dy = y-averageY;
		for( int i=0 ; i<_outW ; i++ )
		{
			double x = double( long long( i ) * long long( _inW-1 ) ) / ( _outW-1 );
			int nearestX = int( x+0.5 );
			int x1 = int( x );
			int x2 = x1 + 1;
			double dx = x-x1;
			if( x2>=_inW ) x2 = x1;

			sum = 0;
			average[0] = average[1] = average[2] = 0;

			if( _separate ) for( int c=0 ; c<3 ; c++ ) nearest[c] = _nearestBuffers[nearestY&1][nearestX+_inW*c];
			else            for( int c=0 ; c<3 ; c++ ) nearest[c] = _nearestBuffers[nearestY&1][3*nearestX+c];

			if( _separate )
			{
				for( int c=0 ; c<3 ; c++ ) temp[c] = _nearestBuffers[y1&1][x1+_inW*c];
				if( temp[0]==nearest[0] && temp[1]==nearest[1] && temp[2]==nearest[2] )
				{
					double d = (1.0-dx) * (1.0-dy);
					sum += d;
					for( int c=0 ; c<3 ; c++ ) average[c] += _averageBuffers[y1&1][x1+_inW*c] * d;
				}
				for( int c=0 ; c<3 ; c++ ) temp[c] = _nearestBuffers[y1&1][x2+_inW*c];
				if( temp[0]==nearest[0] && temp[1]==nearest[1] && temp[2]==nearest[2] )
				{
					double d = (    dx) * (1.0-dy);
					sum += d;
					for( int c=0 ; c<3 ; c++ ) average[c] += _averageBuffers[y1&1][x2+_inW*c] * d;
				}
				for( int c=0 ; c<3 ; c++ ) temp[c] = _nearestBuffers[y2&1][x1+_inW*c];
				if( temp[0]==nearest[0] && temp[1]==nearest[1] && temp[2]==nearest[2] )
				{
					double d = (1.0-dx) * (    dy);
					sum += d;
					for( int c=0 ; c<3 ; c++ ) average[c] += _averageBuffers[y2&1][x1+_inW*c] * d;
				}
				for( int c=0 ; c<3 ; c++ ) temp[c] = _nearestBuffers[y2&1][x2+_inW*c];
				if( temp[0]==nearest[0] && temp[1]==nearest[1] && temp[2]==nearest[2] )
				{
					double d = (    dx) * (    dy);
					sum += d;
					for( int c=0 ; c<3 ; c++ ) average[c] += _averageBuffers[y2&1][x2+_inW*c] * d;
				}
			}
			else
			{
				for( int c=0 ; c<3 ; c++ ) temp[c] = _nearestBuffers[y1&1][x1*3+c];
				if( temp[0]==nearest[0] && temp[1]==nearest[1] && temp[2]==nearest[2] )
				{
					double d = (1.0-dx) * (1.0-dy);
					sum += d;
					for( int c=0 ; c<3 ; c++ ) average[c] += _averageBuffers[y1&1][x1*3+c] * d;
				}
				for( int c=0 ; c<3 ; c++ ) temp[c] = _nearestBuffers[y1&1][x2*3+c];
				if( temp[0]==nearest[0] && temp[1]==nearest[1] && temp[2]==nearest[2] )
				{
					double d = (    dx) * (1.0-dy);
					sum += d;
					for( int c=0 ; c<3 ; c++ ) average[c] += _averageBuffers[y1&1][x2*3+c] * d;
				}
				for( int c=0 ; c<3 ; c++ ) temp[c] = _nearestBuffers[y2&1][x1*3+c];
				if( temp[0]==nearest[0] && temp[1]==nearest[1] && temp[2]==nearest[2] )
				{
					double d = (1.0-dx) * (    dy);
					sum += d;
					for( int c=0 ; c<3 ; c++ ) average[c] += _averageBuffers[y2&1][x1*3+c] * d;
				}
				for( int c=0 ; c<3 ; c++ ) temp[c] = _nearestBuffers[y2&1][x2*3+c];
				if( temp[0]==nearest[0] && temp[1]==nearest[1] && temp[2]==nearest[2] )
				{
					double d = (    dx) * (    dy);
					sum += d;
					for( int c=0 ; c<3 ; c++ ) average[c] += _averageBuffers[y2&1][x2*3+c] * d;
				}
			}
			average[0] /= sum , average[1] /= sum , average[2] /= sum;

			if( _separate ) for( int c=0 ; c<3 ; c++ ) _outBuffer[i+_outW*c] = nearest[c] , _outBuffer[i+_outW*(3+c)] = average[c];
			else            for( int c=0 ; c<3 ; c++ ) _outBuffer[i*6+c] = nearest[c] , _outBuffer[i*6+(3+c)] = average[c];
		}
	}
}
////////////////////////
// JointRImageSampler //
////////////////////////
template< class NearestReal , class AverageReal >
JointRImageSampler< NearestReal , AverageReal >::JointRImageSampler( char* nearestFileName , char* averageFileName , int width , int  height , bool separateColors , bool gammaDecode , double logScale , MultiStreamIOServer* ioServer , int threads , int threadPriority )
{
	int w , h;
	_nearestGrid = GetReadStream< NearestReal >( nearestFileName , w , h ,  separateColors , false , false , ioServer , threads , threadPriority );
	if( logScale>0 ) _averageGrid = new LogImageReader< AverageReal >( averageFileName , _inW , _inH , logScale , separateColors , gammaDecode , false , ioServer , threads , threadPriority );
	else             _averageGrid =      GetReadStream< AverageReal >( averageFileName , _inW , _inH ,            separateColors , gammaDecode , false , ioServer , threads , threadPriority );
	if( w!=_inW || h!=_inH ) fprintf( stderr , "JointRImageSampler images of different sizes: %d x %d != %d x %d\n" , w , h , _inW , _inH ) , exit( 0 );

	_outW = width;
	_outH = height;
	_separate = separateColors;

	_outNearestBuffers[0] = AllocArray< NearestReal >( 3*_outW , 1 );
	_outNearestBuffers[1] = AllocArray< NearestReal >( 3*_outW , 1 );
	_outAverageBuffers[0] = AllocArray< AverageReal >( 3*_outW , 1 );
	_outAverageBuffers[1] = AllocArray< AverageReal >( 3*_outW , 1 );
	_inNearestBuffers[0]  = AllocArray< NearestReal >( 3*_inW  , 1 );
	_inNearestBuffers[1]  = AllocArray< NearestReal >( 3*_inW  , 1 );
	_inAverageBuffers[0]  = AllocArray< AverageReal >( 3*_inW  , 1 );
	_inAverageBuffers[1]  = AllocArray< AverageReal >( 3*_inW  , 1 );
	_inCurrent = 0 , _outCurrent = -1;
	advance( );
	nearestChild = new JointRImageSamplerChild( );
	averageChild = new JointRImageSamplerChild( );
	nearestChild->isNearest = true;
	averageChild->isNearest = false;
	nearestChild->parent = averageChild->parent = this;
}
template< class NearestReal , class AverageReal >
JointRImageSampler< NearestReal , AverageReal >::~JointRImageSampler( void )
{
	FreeArray( _outNearestBuffers[0] );
	FreeArray( _outNearestBuffers[1] );
	FreeArray( _outAverageBuffers[0] );
	FreeArray( _outAverageBuffers[1] );
	FreeArray( _inNearestBuffers[0] );
	FreeArray( _inNearestBuffers[1] );
	FreeArray( _inAverageBuffers[0] );
	FreeArray( _inAverageBuffers[1] );
	delete _nearestGrid;
	delete _averageGrid;
}
template< class NearestReal , class AverageReal > int  JointRImageSampler< NearestReal , AverageReal >::rows( void ) const { return _outH; }
template< class NearestReal , class AverageReal > int  JointRImageSampler< NearestReal , AverageReal >::nearestRowSize( void ) const { return _outW * sizeof( NearestReal ) * 3; }
template< class NearestReal , class AverageReal > int  JointRImageSampler< NearestReal , AverageReal >::averageRowSize( void ) const { return _outW * sizeof( AverageReal ) * 3; }
template< class NearestReal , class AverageReal > bool JointRImageSampler< NearestReal , AverageReal >::SeparateColors( void ) const { return _separate; }
template< class NearestReal , class AverageReal >
Pointer( byte ) JointRImageSampler< NearestReal , AverageReal >::nearestRow( int idx )
{
#if ASSERT_MEMORY_ACCESS
	if( idx!=_outCurrent && idx!=_outCurrent-1 ) fprintf( stderr , "Index out of bounds: %d <= %d <= %d\n" , _outCurrent-1 , idx , _outCurrent ) , exit( 0 );
#endif // ASSERT_MEMORY_ACCESS
	return ( Pointer( byte ) ) _outNearestBuffers[idx&1];
}
template< class NearestReal , class AverageReal >
Pointer( byte ) JointRImageSampler< NearestReal , AverageReal >::averageRow( int idx )
{
#if ASSERT_MEMORY_ACCESS
	if( idx!=_outCurrent && idx!=_outCurrent-1 ) fprintf( stderr , "Index out of bounds: %d <= %d <= %d\n" , _outCurrent-1 , idx , _outCurrent ) , exit( 0 );
#endif // ASSERT_MEMORY_ACCESS
	return ( Pointer( byte ) ) _outAverageBuffers[idx&1];
}
template< class NearestReal , class AverageReal >
void JointRImageSampler< NearestReal , AverageReal >::advance( void )
{
	if( _outCurrent>=_outH ) return;
	_outCurrent++;
	double y = double( long long( _outCurrent ) * long long( _inH-1 ) ) / ( _outH-1 );

	while( _inCurrent-1<=y && _inCurrent<_inH )
	{
		Pointer( NearestReal ) nRow = ( Pointer( NearestReal ) )(*_nearestGrid)[_inCurrent];
		memcpy( _inNearestBuffers[_inCurrent&1] , nRow , _inW * 3 * sizeof( NearestReal ) );
		Pointer( AverageReal ) aRow = ( Pointer( AverageReal ) )(*_averageGrid)[_inCurrent];
		memcpy( _inAverageBuffers[_inCurrent&1] , aRow , _inW * 3 * sizeof( AverageReal ) );
		_nearestGrid->advance();
		_averageGrid->advance();
		_inCurrent++;
	}
	int yy = int( y+0.5 );
	int y1 = int( y ) , y2 = y1+1;
	double dy = y - y1;
	Pointer( NearestReal ) iNRow  = _inNearestBuffers[yy&1];
	Pointer( NearestReal ) iNRow1 = _inNearestBuffers[y1&1];
	Pointer( NearestReal ) iNRow2 = _inNearestBuffers[y2&1];
	Pointer( AverageReal ) iARow1 = _inAverageBuffers[y1&1];
	Pointer( AverageReal ) iARow2 = _inAverageBuffers[y2&1];
	Pointer( NearestReal ) oNRow  = _outNearestBuffers[_outCurrent&1];
	Pointer( AverageReal ) oARow  = _outAverageBuffers[_outCurrent&1];

	for( int i=0 ; i<_outW ; i++ )
	{
		double x = double( long long( i ) * long long( _inW-1 ) ) / ( _outW-1 );
		int xx =int( x+0.5 );
		int x1 = int( x )  , x2 = x1+1;
		double dx = x - x1;
		NearestReal label[3] , label11[3] , label12[3] , label21[3] , label22[3];
		AverageReal value[3] , value11[3] , value12[3] , value21[3] , value22[3];
		if( _separate ) for( int c=0 ; c<3 ; c++ )
		{
			label  [c] = iNRow [xx+_inW*c];
			label11[c] = iNRow1[x1+_inW*c];
			label12[c] = iNRow2[x1+_inW*c];
			label21[c] = iNRow1[x2+_inW*c];
			label22[c] = iNRow2[x2+_inW*c];
			value  [c] = 0;
			value11[c] = iARow1[x1+_inW*c];
			value12[c] = iARow2[x1+_inW*c];
			value21[c] = iARow1[x2+_inW*c];
			value22[c] = iARow2[x2+_inW*c];
		}
		else for( int c=0 ; c<3 ; c++ )
		{
			label  [c] = iNRow [3*xx+c];
			label11[c] = iNRow1[3*x1+c];
			label12[c] = iNRow2[3*x1+c];
			label21[c] = iNRow1[3*x2+c];
			label22[c] = iNRow2[3*x2+c];
			value  [c] = 0;
			value11[c] = iARow1[3*x1+c];
			value12[c] = iARow2[3*x1+c];
			value21[c] = iARow1[3*x2+c];
			value22[c] = iARow2[3*x2+c];
		}
		double weightSum = 0;
		if( label[0]==label11[0] && label[1]==label11[1] && label[2]==label11[2] )
		{
			double weight = (1.-dx)*(1.-dy);
			weightSum += weight;
			value[0] += value11[0]*weight , value[1] += value11[1]*weight , value[2] += value11[2]*weight;
		}
		if( label[0]==label21[0] && label[1]==label21[1] && label[2]==label21[2] )
		{
			double weight = (dx)*(1.-dy);
			weightSum += weight;
			value[0] += value21[0]*weight , value[1] += value21[1]*weight , value[2] += value21[2]*weight;
		}
		if( label[0]==label12[0] && label[1]==label12[1] && label[2]==label12[2] )
		{
			double weight = (1.-dx)*(dy);
			weightSum += weight;
			value[0] += value12[0]*weight , value[1] += value12[1]*weight , value[2] += value12[2]*weight;
		}
		if( label[0]==label22[0] && label[1]==label22[1] && label[2]==label22[2] )
		{
			double weight = (dx)*(dy);
			weightSum += weight;
			value[0] += value22[0]*weight , value[1] += value22[1]*weight , value[2] += value22[2]*weight;
		}
		weightSum = 1./weightSum;
		if( _separate ) for( int c=0 ; c<3 ; c++ ) oNRow[i+_outW*c] = label[c] , oARow[i+_outW*c] = value[c] * weightSum;
		else            for( int c=0 ; c<3 ; c++ ) oNRow[i*3+c]     = label[c] , oARow[i*3+c]     = value[c] * weightSum;
	}
}


////////////////////
// LogImageReader //
////////////////////
template< class Real >
LogImageReader< Real >::LogImageReader( char* fileName , int& width , int& height , Real min , bool separateColors , bool gammaDecode , bool getAlpha , MultiStreamIOServer* ioServer , int threads , int threadPriority )
{
	_min = min;
	_grid = GetReadStream< Real >( fileName , width , height , separateColors , gammaDecode , getAlpha , ioServer , threads , threadPriority );

	_buffer = AllocArray< byte >( _grid->rowSize() );
}
template< class Real >
LogImageReader<Real>::~LogImageReader( void )
{
	FreeArray( _buffer );
	delete _grid;
}
template< class Real > int  LogImageReader< Real >::rows( void ) const { return _grid->rows(); }
template< class Real > int  LogImageReader< Real >::rowSize(void) const { return _grid->rowSize(); }
template< class Real > bool LogImageReader< Real >::SeparateColors( void ) const { return _grid->SeparateColors(); }
template< class Real > bool LogImageReader< Real >::HasAlpha( void ) const { return _grid->HasAlpha(); }
template< class Real >
Pointer( byte ) LogImageReader< Real >::operator [] ( int idx )
{
	int _sz = _grid->rowSize() / sizeof( Real );
	Pointer( Real ) _in  = ( Pointer( Real ) )(*_grid)[idx];
	Pointer( Real ) _out = ( Pointer( Real ) )_buffer;
	for( int i=0 ; i<_sz ; i++ )
#if 1
		if( _in[i]<_min ) _out[i] = Real( log( double( _min   ) ) );
		else              _out[i] = Real( log( double( _in[i] ) ) );
#else
		if( _in[i]<_min ) _out[i] = log( _min );
		else              _out[i] = log( _in[i] );
#endif
	return _buffer;
}
template< class Real >
void LogImageReader< Real >::advance( void )
{
	_grid->advance();
}
////////////////////
// LogImageWriter //
////////////////////
template< class Real >
LogImageWriter< Real >::LogImageWriter( char* fileName , const int& width , const int& height , const bool& separateColors , bool gammaEncode , const int& quality , MultiStreamIOServer* ioServer , bool hdr , int threads , int threadPriority )
{
	_grid = GetWriteStream< Real >( fileName , width , height , separateColors , gammaEncode , quality , ioServer , hdr , threads , threadPriority );
	_buffer = AllocArray< byte >( _grid->rowSize() );

	_index = 0;
}
template< class Real >
LogImageWriter<Real>::~LogImageWriter( void )
{
	FreeArray( _buffer );
	delete _grid;
}
template< class Real > int  LogImageWriter< Real >::rows( void ) const { return _grid->rows(); }
template< class Real > int  LogImageWriter< Real >::rowSize(void) const { return _grid->rowSize(); }
template< class Real > bool LogImageWriter< Real >::SeparateColors( void ) const { return _grid->SeparateColors(); }
template< class Real > bool LogImageWriter< Real >::HasAlpha( void ) const { return _grid->HasAlpha(); }
template< class Real >
Pointer( byte ) LogImageWriter< Real >::operator [] ( int idx )
{
#if ASSERT_MEMORY_ACCESS
	if( idx!=_index ) fprintf( stderr , "Index out of bounds: %d != %d\n" , idx , _index ) , exit( 0 );
#endif // ASSERT_MEMORY_ACCESS
	return  _buffer;
}
template< class Real >
void LogImageWriter< Real >::advance( void )
{
	Pointer( Real )  in = ( Pointer( Real ) ) _buffer;
	Pointer( Real ) out = ( Pointer( Real ) ) (*_grid)[_index];
	int _sz = _grid->rowSize( ) / sizeof( Real );
#if 1
	for( int i=0 ; i<_sz ; i++ ) out[i] = Real( exp( double( in[i] ) ) );
#else
	for( int i=0 ; i<_sz ; i++ ) out[i] = exp( in[i] );
#endif
	_grid->advance();
	_index++;
}

/////////////////////////////
// LogLuminanceImageReader //
/////////////////////////////
template< class Real >
LogLuminanceImageReader< Real >::LogLuminanceImageReader( char* fileName , int& width , int& height , Real min , bool separateColors , bool gammaDecode , MultiStreamIOServer* ioServer , int threads , int threadPriority )
{
	_min = min;
	_grid = GetReadStream< Real >( fileName , width , height , separateColors , gammaDecode , false , ioServer , threads , threadPriority );

	_buffer = AllocArray< byte >( ( _grid->rowSize() * 4 ) / 3 );
}
template< class Real >
LogLuminanceImageReader<Real>::~LogLuminanceImageReader( void )
{
	FreeArray( _buffer );
	delete _grid;
}
template< class Real > int  LogLuminanceImageReader< Real >::rows          ( void ) const { return _grid->rows(); }
template< class Real > int  LogLuminanceImageReader< Real >::rowSize       ( void ) const { return ( _grid->rowSize() * 4 ) / 3; }
template< class Real > bool LogLuminanceImageReader< Real >::SeparateColors( void ) const { return _grid->SeparateColors(); }
template< class Real > bool LogLuminanceImageReader< Real >::HasAlpha      ( void ) const { return false; }
template< class Real >
Pointer( byte ) LogLuminanceImageReader< Real >::operator [] ( int idx )
{
	int _sz = _grid->rowSize() / sizeof( Real );
	_sz /= 3;
	Pointer( Real ) in  = ( Pointer( Real ) )(*_grid)[idx];
	Pointer( Real ) out = ( Pointer( Real ) )_buffer;
	if( _grid->SeparateColors() )
		for( int i=0 ; i<_sz ; i++ )
		{
			Real r = in[i+0*_sz] , g = in[i+1*_sz] , b = in[i+2*_sz];
			if( r<_min ) r = _min;
			if( g<_min ) g = _min;
			if( b<_min ) b = _min;
			Real luminance = Real( double(r)*0.3 + double(g)*0.59 + double(b)*0.11 );
			r /= luminance , g /= luminance , b /= luminance;
			luminance = Real( log( double( luminance ) ) );
			out[i+0*_sz] = luminance , out[i+1*_sz] = r , out[i+2*_sz] = g , out[i+3*_sz] = b;
		}
	else
		for( int i=0 ; i<_sz ; i++ )
		{
			Real r = in[3*i+0] , g = in[3*i+1] , b = in[3*i+2];
			if( r<_min ) r = _min;
			if( g<_min ) g = _min;
			if( b<_min ) b = _min;
			Real luminance = Real( double(r)*0.3 + double(g)*0.59 + double(b)*0.11 );
			r /= luminance , g /= luminance , b /= luminance;
			luminance = Real( log( double( luminance ) ) );
			out[4*i+0] = luminance , out[4*i+1] = r , out[4*i+2] = g , out[4*i+3] = b;
		}
	return _buffer;
}
template< class Real >
void LogLuminanceImageReader< Real >::advance( void )
{
	_grid->advance();
}
/////////////////////////////
// LogLuminanceImageWriter //
/////////////////////////////
template< class Real >
LogLuminanceImageWriter< Real >::LogLuminanceImageWriter( char* fileName , const int& width , const int& height , const bool& separateColors , bool gammaEncode , const int& quality , MultiStreamIOServer* ioServer , bool hdr , int threads , int threadPriority )
{
	_grid = GetWriteStream< Real >( fileName , width , height , separateColors , gammaEncode , quality , ioServer , hdr , threads , threadPriority );
	_buffer = AllocArray< byte >( ( _grid->rowSize() * 4 ) / 3 );

	_index = 0;
}
template< class Real >
LogLuminanceImageWriter<Real>::~LogLuminanceImageWriter( void )
{
	FreeArray( _buffer );
	delete _grid;
}
template< class Real > int  LogLuminanceImageWriter< Real >::rows          ( void ) const { return _grid->rows(); }
template< class Real > int  LogLuminanceImageWriter< Real >::rowSize       ( void ) const { return ( _grid->rowSize() * 4 ) / 3; }
template< class Real > bool LogLuminanceImageWriter< Real >::SeparateColors( void ) const { return _grid->SeparateColors(); }
template< class Real > bool LogLuminanceImageWriter< Real >::HasAlpha      ( void ) const { return _grid->HasAlpha(); }
template< class Real >
Pointer( byte ) LogLuminanceImageWriter< Real >::operator [] ( int idx )
{
#if ASSERT_MEMORY_ACCESS
	if( idx!=_index ) fprintf( stderr , "Index out of bounds: %d != %d\n" , idx , _index ) , exit( 0 );
#endif // ASSERT_MEMORY_ACCESS
	return  _buffer;
}
template< class Real >
void LogLuminanceImageWriter< Real >::advance( void )
{
	Pointer( Real )  in = ( Pointer( Real ) ) _buffer;
	Pointer( Real ) out = ( Pointer( Real ) ) (*_grid)[_index];
	int _sz = _grid->rowSize( ) / sizeof( Real );
	_sz /= 3;
	if( _grid->SeparateColors() )
		for( int i=0 ; i<_sz ; i++ )
		{
			Real luminance = in[i+0*_sz] , r = in[i+1*_sz] , g = in[i+2*_sz] , b = in[i+3*_sz];
			luminance = Real( exp( double( luminance ) ) );
			out[i+0*_sz] = r*luminance , out[i+1*_sz] = g*luminance , out[i+2*_sz] = b*luminance;
		}
	else
		for( int i=0 ; i<_sz ; i++ )
		{
			Real luminance = in[4*i+0] , r = in[4*i+1] , g = in[4*i+2] , b = in[4*i+3];
			luminance = Real( exp( double( luminance ) ) );
			out[3*i+0] = r*luminance , out[3*i+1] = g*luminance , out[3*i+2] = b*luminance;
		}
	_grid->advance();
	_index++;
}

////////////////////////////////
#include <Util/cmdLineParser.h>
template<class Real>
void GetReadSize( char* fileName , int& width , int& height )
{
	char* ext = GetFileExtension( fileName );
	if( !strcasecmp( ext , "jpeg" ) || !strcasecmp( ext , "jpg" ) )	JPEGGetImageSize( fileName , width , height );
	else if( !strcasecmp( ext , "wdp" ) )							 WDPGetImageSize( fileName , width , height );
	else if( !strcasecmp( ext , "bmp" ) )							 BMPGetImageSize( fileName , width , height );
	else if( !strcasecmp( ext , "kro" ) )							 KROGetImageSize( fileName , width , height );
	else if( !strcasecmp( ext , "png" ) )							 PNGGetImageSize( fileName , width , height );
	else if( !strcasecmp( ext , "tiff" ) || !strcasecmp( ext , "tif" ) )	TIFFGetImageSize( fileName , width , height );
	else if( !strcasecmp( ext , "iGrid" ) )							RGridOfImages< Real >::GetImageSize( fileName , width , height );
	else if( !strcasecmp( ext , "int" ) )							GetImageSize< int     >( fileName , width , height );
	else if( !strcasecmp( ext , "int16" ) )							GetImageSize< __int16 >( fileName , width , height );
	else if( !strcasecmp( ext , "half" ) )							GetImageSize< half    >( fileName , width , height );
	else if( !strcasecmp( ext , "float" ) )							GetImageSize< float   >( fileName , width , height );
	else if( !strcasecmp( ext , "double" ) )						GetImageSize< double  >( fileName , width , height );
	else	fprintf( stderr , "Unsupported input file extension: %s\n" , ext );
	delete[] ext;
}

template<class Real>
StreamingGrid* GetReadStream( char* fileName , int& width , int& height , const bool& separateColors , bool gammaDecode , const bool& getAlpha , MultiStreamIOServer* ioServer , int threads , int threadPriority  )
{
	StreamingGrid* data=NULL;
	char* ext=GetFileExtension(fileName);
	if(!strcasecmp(ext,"jpeg") || !strcasecmp(ext,"jpg") )				data = new JPEGRImageStream		<Real>( fileName , width , height , separateColors , gammaDecode , ioServer );
	else if(!strcasecmp( ext , "wdp" ) )								data = new WDPRImageStream		<Real>( fileName , width , height , separateColors , gammaDecode , ioServer );
	else if(!strcasecmp( ext , "kro" ) )
		if( getAlpha )													data = new KRORImageStreamAlpha <Real>( fileName , width , height , separateColors , gammaDecode , ioServer );
		else															data = new KRORImageStream		<Real>( fileName , width , height , separateColors , gammaDecode , ioServer );
	else if(!strcasecmp( ext , "bmp" ) )								data = new BMPRImageStream		<Real>( fileName , width , height , separateColors , gammaDecode , ioServer );
	else if(!strcasecmp( ext , "png" ) )
		if( getAlpha )													data = new PNGRImageStreamAlpha	<Real>( fileName , width , height , separateColors , gammaDecode , ioServer );
		else															data = new PNGRImageStream		<Real>( fileName , width , height , separateColors , gammaDecode , ioServer );
	else if(!strcasecmp( ext , "pfm" ) )								data = new PFMRImageStream		<Real>( fileName , width , height , separateColors , gammaDecode , ioServer );
	else if(!strcasecmp( ext , "tiff" ) || !strcasecmp( ext , "tif" ) )	data = new TIFFRImageStream		<Real>( fileName , width , height , separateColors , gammaDecode , ioServer );
	else if(!strcasecmp( ext , "iGrid" ) )
		if( !threads )													data = new RGridOfImages	<Real>			( fileName , width , height , separateColors , gammaDecode , ioServer );
		else															data = new MultiThreadedRGridOfImages<Real>	( fileName , width , height , separateColors , gammaDecode , threads , threadPriority , ioServer );
	else if(!strcasecmp(ext,"int") )									data = new RImageStream		<Real,int>		( fileName , width , height , separateColors , gammaDecode , ioServer );
	else if(!strcasecmp(ext,"int16") )									data = new RImageStream		<Real,__int16>	( fileName , width , height , separateColors , gammaDecode , ioServer );
	else if(!strcasecmp(ext,"half") )									data = new RImageStream		<Real,half>		( fileName , width , height , separateColors , gammaDecode , ioServer );
	else if(!strcasecmp(ext,"float") )									data = new RImageStream		<Real,float>	( fileName , width , height , separateColors , gammaDecode , ioServer );
	else if(!strcasecmp(ext,"double") )									data = new RImageStream		<Real,double>	( fileName , width , height , separateColors , gammaDecode , ioServer );
	else	fprintf(stderr,"Unsupported input file extension: %s\n",ext);
	delete[] ext;
	return data;
}
template<class Real>
StreamingGrid* GetWriteStream( char* fileName , const int& width , const int& height , const bool& separateColors , bool gammaEncode , const int& quality , MultiStreamIOServer* ioServer , bool hdr , int threads , int threadPriority )
{
	StreamingGrid* data=NULL;
	char* ext=GetFileExtension(fileName);

	if(!strcasecmp(ext,"jpg") || !strcasecmp(ext,"jpeg"))	data=new JPEGWImageStream	<Real>			( fileName , width , height , separateColors , gammaEncode , quality , ioServer );
	else if(!strcasecmp(ext,"wdp"))							data=new WDPWImageStream	<Real>			( fileName , width , height , separateColors , gammaEncode , quality , ioServer );
	else if(!strcasecmp(ext,"bmp"))							data=new BMPWImageStream	<Real>			( fileName , width , height , separateColors , gammaEncode , quality , ioServer );
	else if(!strcasecmp(ext,"kro"))							data=new KROWImageStream	<Real>			( fileName , width , height , separateColors , gammaEncode , quality , ioServer );
	else if(!strcasecmp(ext,"png"))
		if( hdr )											data=new PNGWImageStreamHDR	<Real>			( fileName , width , height , separateColors , gammaEncode , quality , ioServer );
		else												data=new PNGWImageStream	<Real>			( fileName , width , height , separateColors , gammaEncode , quality , ioServer );
	else if( !strcasecmp(ext,"tiff") || !strcasecmp( ext , "tif" ) )
		if( hdr )											data=new TIFFWImageStreamHDR<Real>			( fileName , width , height , separateColors , gammaEncode , quality , ioServer );
		else												data=new TIFFWImageStream	<Real>			( fileName , width , height , separateColors , gammaEncode , quality , ioServer );
	else if(!strcasecmp(ext,"iGrid"))
		if( !threads )										data=new WGridOfImages		<Real>			( fileName , width , height , separateColors , gammaEncode , quality , DefaultOutputTileWidth , DefaultOutputTileHeight , ioServer , hdr );
		else												data=new MultiThreadedWGridOfImages<Real>	( fileName , width , height , separateColors , gammaEncode , quality , threads , threadPriority , DefaultOutputTileWidth , DefaultOutputTileHeight , ioServer , hdr );
	else if(!strcasecmp(ext,"int"))							data=new WImageStream		<Real,int>		( fileName , width , height , separateColors , gammaEncode , quality , ioServer );
	else if(!strcasecmp(ext,"int16"))						data=new WImageStream		<Real,__int16>	( fileName , width , height , separateColors , gammaEncode , quality , ioServer );
	else if(!strcasecmp(ext,"half"))						data=new WImageStream		<Real,half>		( fileName , width , height , separateColors , gammaEncode , quality , ioServer );
	else if(!strcasecmp(ext,"float"))						data=new WImageStream		<Real,float>	( fileName , width , height , separateColors , gammaEncode , quality , ioServer );
	else if(!strcasecmp(ext,"double"))						data=new WImageStream		<Real,double>	( fileName , width , height , separateColors , gammaEncode , quality , ioServer );
	else	fprintf(stderr,"Unsupported output file extension: %s\n",ext);
	delete[] ext;
	return data;
}
