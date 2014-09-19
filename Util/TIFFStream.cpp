#include <stdlib.h>
#include <Util/TIFF/tiffio.h>
#include <Util/BaseMultiStreamIO.h>

// Code taken from:
// http://www.codeguru.com/cpp/g-m/bitmap/otherformats/article.php/c4933/

struct TIFFWriteInfo
{
	TIFF* tiff;
};
struct TIFFReadInfo
{
	TIFF* tiff;
	tdata_t buff;
	int width;
	uint16 bitsPerSample , samplesPerPixel;
};

void* TIFFInitReadColor( char* fileName , int& width , int& height , MultiStreamIOServer* ioServer )
{
	uint16 temp16;
	uint32 temp32;
	TIFFReadInfo* info = (TIFFReadInfo*)malloc(sizeof(TIFFReadInfo));
	info->tiff = TIFFOpen( fileName , "r" );
	if( !info->tiff ) fprintf( stderr , "TIFFOpen failed\n" ) , exit( 0 );

	TIFFGetField( info->tiff , TIFFTAG_IMAGEWIDTH  , &width  );
	info->width = width;
	TIFFGetField( info->tiff , TIFFTAG_IMAGELENGTH , &height );

	TIFFGetField( info->tiff , TIFFTAG_BITSPERSAMPLE   , &(info->bitsPerSample) );
	if( info->bitsPerSample!=8 && info->bitsPerSample!=16 ) fprintf( stderr , "Unsupported number of bits per sample: %d\n" , info->bitsPerSample ) , exit( 0 );
	TIFFGetField( info->tiff , TIFFTAG_SAMPLESPERPIXEL , &info->samplesPerPixel );
	if( info->samplesPerPixel!=1 && info->samplesPerPixel!=3 && info->samplesPerPixel!=4 ) fprintf( stderr , "Unsupported number of samples per pixel: %d\n" , info->bitsPerSample ) , exit( 0 );
	int scanLineSize1 = TIFFScanlineSize( info->tiff ) , scanLineSize2 = ( info->bitsPerSample / 8 ) * info->samplesPerPixel * info->width;
	if( scanLineSize1!=scanLineSize2 ) fprintf( stderr , "Scanline sizes do not agree: %d!=%d\n" , scanLineSize1 , scanLineSize2 ) , exit( 0 );
	TIFFGetField( info->tiff , TIFFTAG_PLANARCONFIG , &temp16 );
	if( temp16!=PLANARCONFIG_CONTIG ) fprintf( stderr , "Planar configuration must be contiguous: %d != %d\n" , temp16 , PLANARCONFIG_CONTIG );
	TIFFGetField( info->tiff , TIFFTAG_ROWSPERSTRIP , &temp32 );
	if( temp32!=1 ) fprintf( stderr , "Expecting one row per strip: %d\n" , temp32 );
	info->buff = _TIFFmalloc( scanLineSize1 );
	return info;
}
const double CharToInt16 = double((1<<16)-1) / double((1<<8)-1);

void TIFFReadRow( void* pixels , void* v , int j )
{
	TIFFReadInfo* info = (TIFFReadInfo*)v;
	TIFFReadScanline( info->tiff , info->buff , j );
	__int16* _pixels = (__int16*)pixels;

	int w = info->width;
	if( info->bitsPerSample==8 )
	{
		__int8* buffer = (__int8*)info->buff;
		if     ( info->samplesPerPixel==1 ) for( int i=0 ; i<w ; i++ ) _pixels[3*i] = _pixels[3*i+1] = _pixels[3*i+2] = __int16(double(buffer[i])*CharToInt16+0.5);
		else if( info->samplesPerPixel==3 ) for( int i=0 ; i<w ; i++ ) _pixels[3*i] = __int16(double(buffer[3*i])*CharToInt16+0.5) , _pixels[3*i+1] = _int16(double(buffer[3*i+1])*CharToInt16+0.5) , _pixels[3*i+2] = _int16(double(buffer[3*i+2])*CharToInt16+0.5);
		else if( info->samplesPerPixel==4 ) for( int i=0 ; i<w ; i++ ) _pixels[3*i] = __int16(double(buffer[4*i])*CharToInt16+0.5) , _pixels[3*i+1] = _int16(double(buffer[4*i+1])*CharToInt16+0.5) , _pixels[3*i+2] = _int16(double(buffer[4*i+2])*CharToInt16+0.5);
	}
	else if( info->bitsPerSample==16 )
	{
		__int16* buffer = (__int16*)info->buff;
		if     ( info->samplesPerPixel==1 ) for( int i=0 ; i<w ; i++ ) _pixels[3*i] = _pixels[3*i+1] = _pixels[3*i+2] = buffer[i];
		else if( info->samplesPerPixel==3 ) for( int i=0 ; i<w ; i++ ) _pixels[3*i] = buffer[3*i] , _pixels[3*i+1] = buffer[3*i+1] , _pixels[3*i+2] = buffer[3*i+2];
		else if( info->samplesPerPixel==4 ) for( int i=0 ; i<w ; i++ ) _pixels[3*i] = buffer[4*i] , _pixels[3*i+1] = buffer[4*i+1] , _pixels[3*i+2] = buffer[4*i+2];
	}
}
void TIFFFinalizeRead(void* v)
{
	TIFFReadInfo* info = (TIFFReadInfo*)v;
	_TIFFfree( info->buff );
	TIFFClose( info->tiff );
	free(info);
}
void TIFFGetImageSize( char* fn , int& width , int& height )
{
	TIFF* tiff;
	tiff = TIFFOpen( fn , "r" );
	if( !tiff ) fprintf( stderr , "TIFFOpen failed\n" ) , exit( 0 );

	TIFFGetField( tiff , TIFFTAG_IMAGEWIDTH  , &width  );
	TIFFGetField( tiff , TIFFTAG_IMAGELENGTH , &height );
	TIFFClose( tiff );
}

void* TIFFInitWriteColor(char* fileName,int width,int height,int quality , MultiStreamIOServer* ioServer )
{
	TIFFWriteInfo* info = (TIFFWriteInfo*)malloc(sizeof(TIFFWriteInfo));

	info->tiff=TIFFOpen(fileName,"w");
	if (!info->tiff)															fprintf(stderr,"TIFFOpen failed\n")								,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_IMAGEWIDTH, width))					fprintf (stderr, "Can't set ImageWidth tag.\n")					,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_IMAGELENGTH, height))					fprintf (stderr, "Can't set ImageLength tag.\n")				,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_BITSPERSAMPLE, sizeof(char)*8))		fprintf (stderr, "Can't set BitsPerSample tag.\n")				,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_SAMPLESPERPIXEL, 3))					fprintf (stderr, "Can't set SamplesPerPixel tag.\n")			,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_COMPRESSION, COMPRESSION_LZW))		fprintf (stderr, "Can't set Compression tag.\n")				,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_ROWSPERSTRIP, 1))						fprintf (stderr, "Can't set RowsPerStrip tag.\n")				,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG))	fprintf (stderr, "Can't set PlanarConfiguration tag.\n")		,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB))		fprintf (stderr, "Can't set PhotometricInterpretation tag.\n")	,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_RESOLUTIONUNIT, RESUNIT_NONE))		fprintf (stderr, "Can't set ResolutionUnit tag.\n")				,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT))	fprintf (stderr, "Can't set Orientation tag.\n")				,	exit(0);

	return info;
}
void* TIFFInitWriteColorHDR(char* fileName,int width,int height,int quality , MultiStreamIOServer* ioServer )
{
	TIFFWriteInfo* info = (TIFFWriteInfo*)malloc(sizeof(TIFFWriteInfo));

	info->tiff=TIFFOpen(fileName,"w");
	if (!info->tiff)															fprintf(stderr,"TIFFOpen failed\n")								,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_IMAGEWIDTH, width))					fprintf (stderr, "Can't set ImageWidth tag.\n")					,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_IMAGELENGTH, height))					fprintf (stderr, "Can't set ImageLength tag.\n")				,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_BITSPERSAMPLE, sizeof(__int16)*8))	fprintf (stderr, "Can't set BitsPerSample tag.\n")				,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_SAMPLESPERPIXEL, 3))					fprintf (stderr, "Can't set SamplesPerPixel tag.\n")			,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_COMPRESSION, COMPRESSION_LZW))		fprintf (stderr, "Can't set Compression tag.\n")				,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_ROWSPERSTRIP, 1))						fprintf (stderr, "Can't set RowsPerStrip tag.\n")				,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG))	fprintf (stderr, "Can't set PlanarConfiguration tag.\n")		,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB))		fprintf (stderr, "Can't set PhotometricInterpretation tag.\n")	,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_RESOLUTIONUNIT, RESUNIT_NONE))		fprintf (stderr, "Can't set ResolutionUnit tag.\n")				,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT))	fprintf (stderr, "Can't set Orientation tag.\n")				,	exit(0);

	return info;
}
void* TIFFInitWriteGray(char* fileName,int width,int height,int quality , MultiStreamIOServer* ioServer )
{
	TIFFWriteInfo* info = (TIFFWriteInfo*)malloc(sizeof(TIFFWriteInfo));

	info->tiff=TIFFOpen(fileName,"w");
	if (!info->tiff)															fprintf(stderr,"TIFFOpen failed\n")								,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_IMAGEWIDTH, width))					fprintf (stderr, "Can't set ImageWidth tag.\n")					,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_IMAGELENGTH, height))					fprintf (stderr, "Can't set ImageLength tag.\n")				,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_BITSPERSAMPLE, sizeof(char)*8))		fprintf (stderr, "Can't set BitsPerSample tag.\n")				,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_SAMPLESPERPIXEL, 1))					fprintf (stderr, "Can't set SamplesPerPixel tag.\n")			,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_COMPRESSION, COMPRESSION_LZW))		fprintf (stderr, "Can't set Compression tag.\n")				,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_ROWSPERSTRIP, 1))						fprintf (stderr, "Can't set RowsPerStrip tag.\n")				,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG))	fprintf (stderr, "Can't set PlanarConfiguration tag.\n")		,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK))	fprintf (stderr, "Can't set PhotometricInterpretation tag.\n")	,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_RESOLUTIONUNIT, RESUNIT_NONE))		fprintf (stderr, "Can't set ResolutionUnit tag.\n")				,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT))	fprintf (stderr, "Can't set Orientation tag.\n")				,	exit(0);

	return info;
}
void* TIFFInitWriteGrayHDR(char* fileName,int width,int height,int quality , MultiStreamIOServer* ioServer )
{
	TIFFWriteInfo* info = (TIFFWriteInfo*)malloc(sizeof(TIFFWriteInfo));

	info->tiff=TIFFOpen(fileName,"w");
	if (!info->tiff)															fprintf(stderr,"TIFFOpen failed\n")								,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_IMAGEWIDTH, width))					fprintf (stderr, "Can't set ImageWidth tag.\n")					,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_IMAGELENGTH, height))					fprintf (stderr, "Can't set ImageLength tag.\n")				,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_BITSPERSAMPLE, sizeof(__int16)*8))	fprintf (stderr, "Can't set BitsPerSample tag.\n")				,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_SAMPLESPERPIXEL, 1))					fprintf (stderr, "Can't set SamplesPerPixel tag.\n")			,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_COMPRESSION, COMPRESSION_LZW))		fprintf (stderr, "Can't set Compression tag.\n")				,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_ROWSPERSTRIP, 1))						fprintf (stderr, "Can't set RowsPerStrip tag.\n")				,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG))	fprintf (stderr, "Can't set PlanarConfiguration tag.\n")		,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK))	fprintf (stderr, "Can't set PhotometricInterpretation tag.\n")	,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_RESOLUTIONUNIT, RESUNIT_NONE))		fprintf (stderr, "Can't set ResolutionUnit tag.\n")				,	exit(0);
	if (!TIFFSetField(info->tiff, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT))	fprintf (stderr, "Can't set Orientation tag.\n")				,	exit(0);

	return info;
}
void TIFFFinalizeWrite(void* v)
{
	TIFFWriteInfo* info = (TIFFWriteInfo*)v;
	TIFFClose(info->tiff);
	free(info);
}
void TIFFWriteRow(void* pixels,void* v,int j)
{
	unsigned char* _pixels = (unsigned char*)pixels;
	TIFFWriteInfo* info = (TIFFWriteInfo*)v;
	TIFFWriteScanline(info->tiff,pixels,j,0);
}
void TIFFWriteRowHDR(void* pixels,void* v,int j)
{
	unsigned __int16* _pixels = (unsigned __int16*)pixels;
	TIFFWriteInfo* info = (TIFFWriteInfo*)v;
	TIFFWriteScanline(info->tiff,pixels,j,0);
}
