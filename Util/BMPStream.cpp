#include <windows.h>
#include <stdio.h>
#include <Util/BaseMultiStreamIO.h>

/* constants for the biCompression field */
#define BI_RGB        0L
#define BI_RLE8       1L
#define BI_RLE4       2L
#define BI_BITFIELDS  3L

/* Some magic numbers */

#define BMP_BF_TYPE 0x4D42
/* word BM */

#define BMP_BF_OFF_BITS 54
/* 14 for file header + 40 for info header (not sizeof(), but packed size) */

#define BMP_BI_SIZE 40
/* packed size of info header */


struct BMPWriteInfo
{
	FILE* fp;
	int width;
};

struct BMPReadInfo
{
	FILE* fp;
	BYTE* data;
	int width , lineLength;
};

void* BMPInitReadColor( char* fileName , int& width , int& height , MultiStreamIOServer* ioServer )
{
    BITMAPFILEHEADER bmfh;
    BITMAPINFOHEADER bmih;

	BMPReadInfo* info = (BMPReadInfo*)malloc(sizeof(BMPReadInfo));
	info->fp = fopen( fileName , "rb" );
	if(!info->fp)	fprintf(stderr,"Failed to open: %s\n",fileName) , exit(0);


	fread( &bmfh , sizeof( BITMAPFILEHEADER ) , 1 , info->fp );
	fread( &bmih , sizeof( BITMAPINFOHEADER ) , 1 , info->fp );

	if( bmfh.bfType!=BMP_BF_TYPE || bmfh.bfOffBits!=BMP_BF_OFF_BITS ) fprintf( stderr , "Bad bitmap file header\n" ) , exit( 0 );
	if( bmih.biSize!=BMP_BI_SIZE || bmih.biWidth<=0 || bmih.biHeight<=0 || bmih.biPlanes!=1 || bmih.biBitCount!=24 || bmih.biCompression!=BI_RGB ) fprintf( stderr , "Bad bitmap file info\n" ) , exit( 0 );

	info->width = width  = bmih.biWidth;
	height = bmih.biHeight;
	info->lineLength = width * 3;
	if( (info->lineLength % 4) != 0) info->lineLength = (info->lineLength / 4 + 1) * 4;
	if( bmih.biSizeImage!=info->lineLength*height ) fprintf( stderr , "Bad bitmap image size\n" ) , exit( 0 );
	info->data = new BYTE[ info->lineLength ];
	if( !info->data ) fprintf( stderr , "Could not allocate memory for bitmap data\n" ) , exit( 0 );

	fseek( info->fp , (long) bmfh.bfOffBits , SEEK_SET );
	fseek( info->fp , (long) info->lineLength * height , SEEK_CUR );
	return info;
}
void BMPGetImageSize( char* fileName , int& width , int& height )
{
    BITMAPFILEHEADER bmfh;
    BITMAPINFOHEADER bmih;

	FILE* fp = fopen( fileName , "rb" );
	if( !fp ) fprintf( stderr , "Failed to open: %s\n" , fileName ) , exit(0);


	fread( &bmfh , sizeof( BITMAPFILEHEADER ) , 1 , fp );
	fread( &bmih , sizeof( BITMAPINFOHEADER ) , 1 , fp );

	if( bmfh.bfType!=BMP_BF_TYPE || bmfh.bfOffBits!=BMP_BF_OFF_BITS ) fprintf( stderr , "Bad bitmap file header\n" ) , fclose( fp ) , exit( 0 );
	if( bmih.biSize!=BMP_BI_SIZE || bmih.biWidth<=0 || bmih.biHeight<=0 || bmih.biPlanes!=1 || bmih.biBitCount!=24 || bmih.biCompression!=BI_RGB ) fprintf( stderr , "Bad bitmap file info\n" ) , fclose( fp ) , exit( 0 );

	width  = bmih.biWidth;
	height = bmih.biHeight;
	int lineLength = width * 3;
	if( (lineLength % 4) != 0) lineLength = (lineLength / 4 + 1) * 4;
	if( bmih.biSizeImage!=lineLength*height ) fprintf( stderr , "Bad bitmap image size\n" ) , fclose( fp ) , exit( 0 );
	fclose( fp );
}

void* BMPInitWrite(char* fileName,int width,int height,int quality , MultiStreamIOServer* ioServer )
{
	BMPWriteInfo* info=(BMPWriteInfo*)malloc(sizeof(BMPWriteInfo));
	info->fp=fopen(fileName,"wb");
	if(!info->fp)	fprintf(stderr,"Failed to open: %s\n",fileName) , exit(0);
	info->width=width;

	BITMAPFILEHEADER bmfh;
	BITMAPINFOHEADER bmih;
	int lineLength;

	lineLength = width * 3;	/* RGB */
	if ((lineLength % 4) != 0)	lineLength = (lineLength / 4 + 1) * 4;
	/* Write file header */

	bmfh.bfType = BMP_BF_TYPE;
	bmfh.bfSize = BMP_BF_OFF_BITS + lineLength * height;
	bmfh.bfReserved1 = 0;
	bmfh.bfReserved2 = 0;
	bmfh.bfOffBits = BMP_BF_OFF_BITS;

	fwrite(&bmfh,sizeof(BITMAPFILEHEADER),1,info->fp);

	bmih.biSize = BMP_BI_SIZE;
	bmih.biWidth = width;
	bmih.biHeight = -height;
	bmih.biPlanes = 1;
	bmih.biBitCount = 24;		/* RGB */
	bmih.biCompression = BI_RGB;	/* RGB */
	bmih.biSizeImage = lineLength * (DWORD) bmih.biHeight;	/* RGB */
	bmih.biXPelsPerMeter = 2925;
	bmih.biYPelsPerMeter = 2925;
	bmih.biClrUsed = 0;
	bmih.biClrImportant = 0;

	fwrite(&bmih,sizeof(BITMAPINFOHEADER),1,info->fp);

	return info;
}
void BMPWriteRow(void* pixels,void* v,int j)
{
	unsigned char* _pixels = (unsigned char*) pixels;
	BMPWriteInfo* info=(BMPWriteInfo*)v;
	for(int x=0;x<info->width;x++)
	{
		unsigned char tmp=_pixels[x*3];
		_pixels[x*3+0]=_pixels[x*3+2];
		_pixels[x*3+2]=tmp;
	}
	fwrite(pixels,sizeof(unsigned char),info->width*3,info->fp);
	int nbytes=info->width*3;
	while ((nbytes % 4) != 0) {
		putc(0,info->fp);
		nbytes++;
	}
}
void BMPReadColorRow( void* pixels , void* v , int j )
{
	BMPReadInfo* info = (BMPReadInfo*)v;
	fseek( info->fp , -info->lineLength , SEEK_CUR );
    fread( info->data , 1 , info->lineLength , info->fp );
	fseek( info->fp , -info->lineLength , SEEK_CUR );
	if( ferror(info->fp) ) fprintf( stderr , "Error reading bitmap row\n" ) , exit( 0 );
	for( int i=0 ; i<info->width ; i++ )
	{
		BYTE foo = info->data[3*i];
		info->data[3*i] = info->data[3*i+2];
		info->data[3*i+2] = foo;
	}
	unsigned char *p = (unsigned char*) pixels;
	memcpy( pixels , info->data , sizeof( unsigned char ) * 3 * info->width );
}

void BMPFinalizeWrite(void* v)
{
	BMPWriteInfo* info=(BMPWriteInfo*)v;
	fclose(info->fp);
	free(info);
}
void BMPFinalizeRead( void* v  )
{
	BMPReadInfo* info=(BMPReadInfo*)v;
	fclose(info->fp);
	delete[] info->data;
	free(info);
}