#include <stdio.h>
#include <stdlib.h>
#include <Util/BaseMultiStreamIO.h>

struct PFMReadInfo
{
	FILE* fp;
	int width , height;
};

void PFMGetImageSize( char* fn , int& width , int& height )
{
	PFMReadInfo info;
	info.fp = fopen( fn , "rb" );
	if( !info.fp ) fprintf( stderr , "Failed to open: %s\n" , fn ) , exit(0);

	char buf[1024];
	fgets( buf , 1024 , info.fp );
	fgets( buf , 1024 , info.fp );
	while (buf[0] == '#') fgets( buf , 1024 , info.fp );
	sscanf(buf, "%d %d" , &info.width , &info.height);
	if( info.width<0 || info.height<0 )
	{
		fprintf(stderr, "Invalid header or image not valid PFM.\n");
		exit( 0 );
	}
	width=info.width;
	height=info.height;

	fclose( info.fp );
}

void* PFMInitReadColor( char* fileName , int& width , int& height , MultiStreamIOServer* ioServer )
{
	PFMReadInfo *info=(PFMReadInfo*)malloc(sizeof(PFMReadInfo));
	info->fp = fopen( fileName , "rb" );
	if( !info->fp ) fprintf( stderr , "Failed to open: %s\n" , fileName ) , exit(0);

	char buf[1024];
	fgets( buf , 1024 , info->fp );
	fgets( buf , 1024 , info->fp );
	while (buf[0] == '#') fgets( buf , 1024 , info->fp );
	sscanf(buf, "%d %d", &info->width, &info->height);
	if( info->width<0 || info->height<0 )
	{
		fprintf(stderr, "Invalid header or image not valid PFM.\n");
		exit( 0 );
	}
	fgets( buf, 1024 , info->fp );
	while ( buf[0] == '#') fgets( buf , 1024 , info->fp );


	width=info->width;
	height=info->height;

	return info;
}
void PFMReadRow( void* pixels , void* v , int j )
{
	float* _pixels = (float*) pixels;
	PFMReadInfo* info = ( PFMReadInfo* )v;
	if( !fread( _pixels , sizeof(float) , info->width * 3 , info->fp ) )
	{
		fprintf( stderr , "Failed to read PFM row of size %d\n" , info->width );
		exit( 0 );
	}
}

void PFMFinalizeRead( void* v )
{
	PFMReadInfo* info = (PFMReadInfo*)v;
	fclose( info->fp );
	free( info );
}