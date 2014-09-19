#include <stdio.h>
#include <stdlib.h>
#include <Util/BaseMultiStreamIO.h>
#include <Util/ChannelConverter.h>

// Code courtesy of: http://www.codeguru.com/forum/showthread.php?t=292902
inline void endian_swap(unsigned short& x)
{
	x = (x>>8) | 
		(x<<8);
}

inline void endian_swap(unsigned int& x)
{
	x = (x>>24) | 
		((x<<8) & 0x00FF0000) |
		((x>>8) & 0x0000FF00) |
		(x<<24);
}

// __int64 for MSVC, "long long" for gcc
inline void endian_swap(unsigned __int64& x)
{
	x = (x>>56) | 
		((x<<40) & 0x00FF000000000000) |
		((x<<24) & 0x0000FF0000000000) |
		((x<<8)  & 0x000000FF00000000) |
		((x>>8)  & 0x00000000FF000000) |
		((x>>24) & 0x0000000000FF0000) |
		((x>>40) & 0x000000000000FF00) |
		(x<<56);
}

/*
Header is 20 bytes long :
3 bytes : "KRO" signature in hex 0x4B 0x52 0x4F
1 byte  : 0x01 version
unsigned long : Width
unsigned long : Height
unsigned long : depth = > 8 bits, 16 bits, 32 bits
unsigned long : ncomp => number of compoment, 4 by default, RGB + Alpha
*/

const char KRO_HEADER[] = { 0x4B , 0x52 , 0x4F };
const char KRO_VERSION  = 0x01;
struct KROHeader
{
	char header[3];
	char version;
#if 1
	unsigned int width;
	unsigned int height;
	unsigned int depth;
	unsigned int ncomp;
#else
	unsigned long width;
	unsigned long height;
	unsigned long depth;
	unsigned long ncomp;
#endif
};
struct KROInfo
{
	FILE* fp;
	int width , height;
	int bytesPerChannel , channelsPerPixel;
	void *row;
};

void KROGetImageSize( char* fn , int& width , int& height )
{
	KROInfo info;
	info.fp = fopen( fn , "rb" );
	if( !info.fp ) fprintf( stderr , "Failed to open: %s\n" , fn ) , exit(0);

	KROHeader header;
	fread( &header , 1 , sizeof( KROHeader ) , info.fp );
	endian_swap( header.width ) , endian_swap( header.height ) , endian_swap( header.depth ) , endian_swap( header.ncomp );
	if( header.header[0]!=KRO_HEADER[0] || header.header[1]!=KRO_HEADER[1] || header.header[2]!=KRO_HEADER[2] )
		fprintf( stderr , "Invalid header in %s\n" , fn ) , exit( 0 );

	width  = header.width;
	height = header.height;
	fclose( info.fp );
}
void* KROInitReadColorAlpha( char* fileName , int& width , int& height , MultiStreamIOServer* ioServer )
{
	KROInfo *info = ( KROInfo* )malloc( sizeof( KROInfo ) );
	info->fp = fopen( fileName , "rb" );
	if( !info->fp ) fprintf( stderr , "Failed to open: %s\n" , fileName ) , exit(0);

	KROHeader header;
	fread( &header , 1 , sizeof( KROHeader ) , info->fp );
	endian_swap( header.width ) , endian_swap( header.height ) , endian_swap( header.depth ) , endian_swap( header.ncomp );
	if( header.header[0]!=KRO_HEADER[0] || header.header[1]!=KRO_HEADER[1] || header.header[2]!=KRO_HEADER[2] )
		fprintf( stderr , "Invalid header in %s\n" , fileName ) , exit( 0 );
	if( header.depth!=8 && header.depth!=16 && header.depth!=32 )
		fprintf( stderr , "Only 8 , 16, and 32 bits-per-channel supported in KRO file-format: %d\n" , header.ncomp ) , exit( 0 );
	if( header.ncomp!=4 )
		fprintf( stderr , "RGBA mode not supported for KRO file-format: %d\n" , header.depth ) , exit( 0 );
	info->width  = header.width;
	info->height = header.height;
	info->bytesPerChannel  = header.depth / 8;
	info->channelsPerPixel = header.ncomp;

	info->row = malloc( info->width * info->bytesPerChannel * info->channelsPerPixel );
	if( !info->row ) fprintf( stderr , "Failed to allocate KRO row of size: %d\n" , header.width ) , exit( 0 );

	width  = info->width;
	height = info->height;
	return info;
}
void* KROInitReadColor( char* fileName , int& width , int& height , MultiStreamIOServer* ioServer )
{
	KROInfo *info = ( KROInfo* )malloc( sizeof( KROInfo ) );
	info->fp = fopen( fileName , "rb" );
	if( !info->fp ) fprintf( stderr , "Failed to open: %s\n" , fileName ) , exit(0);

	KROHeader header;
	fread( &header , 1 , sizeof( KROHeader ) , info->fp );
	endian_swap( header.width ) , endian_swap( header.height ) , endian_swap( header.depth ) , endian_swap( header.ncomp );
	if( header.header[0]!=KRO_HEADER[0] || header.header[1]!=KRO_HEADER[1] || header.header[2]!=KRO_HEADER[2] )
		fprintf( stderr , "Invalid header in %s\n" , fileName ) , exit( 0 );
	if( header.depth!=8 && header.depth!=16 && header.depth!=32 )
		fprintf( stderr , "Only 8 , 16, and 32 bits-per-channel supported in KRO file-format: %d\n" , header.ncomp ) , exit( 0 );
	if( header.ncomp!=3 && header.ncomp!=4 )
		fprintf( stderr , "Only RGB and RGBA supported in KRO file-format: %d\n" , header.depth ) , exit( 0 );
	info->width  = header.width;
	info->height = header.height;
	info->bytesPerChannel  = header.depth / 8;
	info->channelsPerPixel = header.ncomp;

	info->row = malloc( info->width * info->bytesPerChannel * info->channelsPerPixel );
	if( !info->row ) fprintf( stderr , "Failed to allocate KRO row of size: %d\n" , header.width ) , exit( 0 );

	width  = info->width;
	height = info->height;
	return info;
}
void* KROInitWriteColor( char* fileName , int width , int height , int quality , MultiStreamIOServer* ioServer )
{
	KROInfo *info = ( KROInfo* )malloc( sizeof( KROInfo ) );
	info->fp = fopen( fileName , "wb" );
	if( !info->fp ) fprintf( stderr , "Failed to open: %s\n" , fileName ) , exit(0);

	KROHeader header;
	header.header[0] = KRO_HEADER[0] , header.header[1] = KRO_HEADER[1] , header.header[2] = KRO_HEADER[2];
	header.version = KRO_VERSION;
	header.width   = width;
	header.height  = height;
	header.depth   = 8;
	header.ncomp   = 3;

	info->row = malloc( info->width * info->bytesPerChannel * info->channelsPerPixel );
	if( !info->row ) fprintf( stderr , "Failed to allocate KRO row of size: %d\n" , header.width ) , exit( 0 );

	endian_swap( header.width ) , endian_swap( header.height ) , endian_swap( header.depth ) , endian_swap( header.ncomp );
	fwrite( &header , 1 , sizeof( KROHeader ) , info->fp );

	info->width  = header.width;
	info->height = header.height;
	info->bytesPerChannel  = header.depth / 8;
	info->channelsPerPixel = header.ncomp;

	return info;
}
void KROReadColorRow( void* pixels , void* v , int j )
{
	KROInfo* info = ( KROInfo* )v;
	unsigned __int16* _pixels = (unsigned __int16*) pixels;
	fread( info->row , 1 , info->bytesPerChannel * info->channelsPerPixel * info->width , info->fp );
	switch( info->bytesPerChannel )
	{
	case 1:
		{
			unsigned char* row = (unsigned char*) info->row;
			for( int x=0 ; x<info->width ; x++ )
			{
				_pixels[x*3+0] = ConvertChannel< unsigned char , unsigned __int16 >( row[x*info->channelsPerPixel+0] );
				_pixels[x*3+1] = ConvertChannel< unsigned char , unsigned __int16 >( row[x*info->channelsPerPixel+1] );
				_pixels[x*3+2] = ConvertChannel< unsigned char , unsigned __int16 >( row[x*info->channelsPerPixel+2] );
			}
		}
		break;
	case 2:
		{
			unsigned __int16* row = (unsigned __int16*) info->row;
			for( int x=0 ; x<info->width ; x++ )
			{
				_pixels[x*3+0] = ConvertChannel< unsigned __int16 , unsigned __int16 >( row[x*info->channelsPerPixel+0] );
				_pixels[x*3+1] = ConvertChannel< unsigned __int16 , unsigned __int16 >( row[x*info->channelsPerPixel+1] );
				_pixels[x*3+2] = ConvertChannel< unsigned __int16 , unsigned __int16 >( row[x*info->channelsPerPixel+2] );
			}
		}
		break;
	case 4:
		{
			float* row = (float*) info->row;
			for( int x=0 ; x<info->width ; x++ )
			{
				_pixels[x*3+0] = ConvertChannel< float , unsigned __int16 >( row[x*info->channelsPerPixel+0] );
				_pixels[x*3+1] = ConvertChannel< float , unsigned __int16 >( row[x*info->channelsPerPixel+1] );
				_pixels[x*3+2] = ConvertChannel< float , unsigned __int16 >( row[x*info->channelsPerPixel+2] );
			}
		}
		break;
	default:
		fprintf( stderr , "Invalid number of bytes per channel: %d\n" , info->bytesPerChannel );
		exit( 0 );
	}
}
void KROReadColorRowAlpha( void* pixels , void* v , int j )
{
	KROInfo* info = ( KROInfo* )v;
	unsigned __int16* _pixels = (unsigned __int16*) pixels;
	fread( info->row , 1 , info->bytesPerChannel * info->channelsPerPixel * info->width , info->fp );
	switch( info->bytesPerChannel )
	{
	case 1:
		{
			unsigned char* row = (unsigned char*) info->row;
			for( int x=0 ; x<info->width ; x++ )
			{
				_pixels[x*4+0] = ConvertChannel< unsigned char , unsigned __int16 >( row[x*info->channelsPerPixel+0] );
				_pixels[x*4+1] = ConvertChannel< unsigned char , unsigned __int16 >( row[x*info->channelsPerPixel+1] );
				_pixels[x*4+2] = ConvertChannel< unsigned char , unsigned __int16 >( row[x*info->channelsPerPixel+2] );
				_pixels[x*4+3] = ConvertChannel< unsigned char , unsigned __int16 >( row[x*info->channelsPerPixel+3] );
			}
		}
		break;
	case 2:
		{
			unsigned __int16* row = (unsigned __int16*) info->row;
			for( int x=0 ; x<info->width ; x++ )
			{
				_pixels[x*4+0] = ConvertChannel< unsigned __int16 , unsigned __int16 >( row[x*info->channelsPerPixel+0] );
				_pixels[x*4+1] = ConvertChannel< unsigned __int16 , unsigned __int16 >( row[x*info->channelsPerPixel+1] );
				_pixels[x*4+2] = ConvertChannel< unsigned __int16 , unsigned __int16 >( row[x*info->channelsPerPixel+2] );
				_pixels[x*4+3] = ConvertChannel< unsigned __int16 , unsigned __int16 >( row[x*info->channelsPerPixel+3] );
			}
		}
		break;
	case 4:
		{
			float* row = (float*) info->row;
			for( int x=0 ; x<info->width ; x++ )
			{
				_pixels[x*4+0] = ConvertChannel< float , unsigned __int16 >( row[x*info->channelsPerPixel+0] );
				_pixels[x*4+1] = ConvertChannel< float , unsigned __int16 >( row[x*info->channelsPerPixel+1] );
				_pixels[x*4+2] = ConvertChannel< float , unsigned __int16 >( row[x*info->channelsPerPixel+2] );
				_pixels[x*4+3] = ConvertChannel< float , unsigned __int16 >( row[x*info->channelsPerPixel+3] );
			}
		}
		break;
	default:
		fprintf( stderr , "Invalid number of bytes per channel: %d\n" , info->bytesPerChannel );
		exit( 0 );
	}
}
void KROWriteRow( void* pixels , void* v , int j )
{
	KROInfo* info = ( KROInfo* )v;
	unsigned __int16* _pixels = (unsigned __int16*) pixels;
	memset( info->row , 0 , sizeof( info->channelsPerPixel * info->bytesPerChannel * info->width ) );
	switch( info->bytesPerChannel )
	{
	case 1:
		{
			unsigned char* row = (unsigned char*) info->row;
			for( int x=0 ; x<info->width ; x++ )
			{
				row[x*info->channelsPerPixel+0] = ConvertChannel< unsigned __int16 , unsigned char >( _pixels[x*3+0] );
				row[x*info->channelsPerPixel+1] = ConvertChannel< unsigned __int16 , unsigned char >( _pixels[x*3+1] );
				row[x*info->channelsPerPixel+2] = ConvertChannel< unsigned __int16 , unsigned char >( _pixels[x*3+2] );
			}
		}
		break;
	case 2:
		{
			unsigned __int16* row = (unsigned __int16*) info->row;
			for( int x=0 ; x<info->width ; x++ )
			{
				row[x*info->channelsPerPixel+0] = ConvertChannel< unsigned __int16 , unsigned __int16 >( _pixels[x*3+0] );
				row[x*info->channelsPerPixel+1] = ConvertChannel< unsigned __int16 , unsigned __int16 >( _pixels[x*3+1] );
				row[x*info->channelsPerPixel+2] = ConvertChannel< unsigned __int16 , unsigned __int16 >( _pixels[x*3+2] );
			}
		}
		break;
	case 4:
		{
			float* row = (float*) info->row;
			for( int x=0 ; x<info->width ; x++ )
			{
				row[x*info->channelsPerPixel+0] = ConvertChannel< unsigned __int16 , float >( _pixels[x*3+0] );
				row[x*info->channelsPerPixel+1] = ConvertChannel< unsigned __int16 , float >( _pixels[x*3+1] );
				row[x*info->channelsPerPixel+2] = ConvertChannel< unsigned __int16 , float >( _pixels[x*3+2] );
			}
		}
		break;
	default:
		fprintf( stderr , "Invalid number of bytes per channel: %d\n" , info->bytesPerChannel );
		exit( 0 );
	}
	fwrite( info->row , 1 , info->bytesPerChannel * info->channelsPerPixel * info->width , info->fp );
}

void KROFinalizeRead( void* v )
{
	KROInfo* info = (KROInfo*)v;
	fclose( info->fp );
	if( info->row ) free( info->row  );
	info->row = NULL;
	free( info );
}
void KROFinalizeWrite( void* v )
{
	KROInfo* info = (KROInfo*)v;
	fclose( info->fp );
	if( info->row ) free( info->row  );
	info->row = NULL;
	free( info );
}
