#include <stdlib.h>
#include <Util/PNG/png.h>
#include <Util/MultiStreamIO.h>

struct PNGWriteInfo
{
#if NEW_PNG_IO
	VariableIOClientStream* outStream;
#else // !NEW_PNG_IO
	FILE* fp;
#endif // NEW_PNG_IO
	png_structp png_ptr;
	png_infop info_ptr;
};
struct PNGReadInfo
{
#if NEW_PNG_IO
	VariableIOClientStream* inStream;
#else // !NEW_PNG_IO
	FILE* fp;
#endif // NEW_PNG_IO
	png_structp png_ptr;
	png_infop info_ptr,end_info;
	bool is16bit;
	int width;
	unsigned char* pixels;
	bool hasAlpha;
};
/*
png_set_write_fn( png_strctp png_ptr , png_voidp io_ptr , png_rw_ptr write_data_fn , png_flush_ptr ouput_flush_fn );
typedef void (PNGAPI *png_rw_ptr) PNGARG((png_structp, png_bytep, png_size_t));
typedef void (PNGAPI *png_flush_ptr) PNGARG((png_structp));
png_set_read_fn( png_structp png_ptr , png_voidp io_ptr , png_rw_ptr read_data_fn )
*/

void myPngWriteFunction( png_structp png_ptr , png_bytep buffer , png_size_t byteNum )
{
	VariableIOClientStream* outStream = (VariableIOClientStream*)png_get_io_ptr( png_ptr );
	if( outStream->write( buffer , byteNum ) != byteNum ) fprintf( stderr , "PNG Failed to write %d bytes\n" , byteNum );
}
void myPngReadFunction( png_structp png_ptr , png_bytep buffer , png_size_t byteNum )
{
	VariableIOClientStream* inStream = (VariableIOClientStream*)png_get_io_ptr( png_ptr );
	if( inStream->read( buffer , byteNum ) != byteNum ) fprintf( stderr , "PNG Failed to read %d bytes\n" , byteNum );
}
void myPngFlushFunction( png_structp png_ptr )
{
}
void* PNGInitWriteColor(char* fileName,int width,int height,int quality , MultiStreamIOServer* ioServer )
{
	PNGWriteInfo* info = (PNGWriteInfo*)malloc(sizeof(PNGWriteInfo));
#if NEW_PNG_IO
	info->outStream = new VariableIOClientStream( );
	info->outStream->Initialize( fileName , PNG_ZBUF_SIZE , 2 , false );
	info->outStream->SetServer( ioServer );
#else // !NEW_PNG_IO
	info->fp = fopen(fileName,"wb");
	if(!info->fp)	fprintf(stderr,"Failed to open: %s\n",fileName)	,	exit(0);
#endif // NEW_PNG_IO

	info->png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,0,0,0);
	if(!info->png_ptr)	fprintf(stderr,"Failed to create png write struct\n")	,	exit(0);
	info->info_ptr = png_create_info_struct(info->png_ptr);
	if(!info->info_ptr)	fprintf(stderr,"Failed to create png info struct\n")	,	exit(0);

#if NEW_PNG_IO
	png_set_write_fn( info->png_ptr , info->outStream , myPngWriteFunction , myPngFlushFunction );
#else // !NEW_JOEG_IO
	png_init_io(info->png_ptr, info->fp);
#endif // NEW_PNG_IO


	png_set_compression_level(info->png_ptr,Z_BEST_SPEED);

	png_set_IHDR(info->png_ptr, info->info_ptr, width, height,
		8,PNG_COLOR_TYPE_RGB,
		PNG_INTERLACE_NONE,
		PNG_COMPRESSION_TYPE_DEFAULT,
		PNG_FILTER_TYPE_DEFAULT);
	png_write_info(info->png_ptr, info->info_ptr);

	{
		long int a = 1;
		int swap = (*((unsigned char *) &a) == 1);
		if(swap)	png_set_swap(info->png_ptr);
	}
	return info;
}
void* PNGInitWriteColorHDR(char* fileName,int width,int height,int quality , MultiStreamIOServer* ioServer )
{
	PNGWriteInfo* info = (PNGWriteInfo*)malloc(sizeof(PNGWriteInfo));
#if NEW_PNG_IO
	info->outStream = new VariableIOClientStream( );
	info->outStream->Initialize( fileName , PNG_ZBUF_SIZE , 2 , false );
	info->outStream->SetServer( ioServer );
#else // !NEW_PNG_IO
	info->fp = fopen(fileName,"wb");
	if(!info->fp)	fprintf(stderr,"Failed to open: %s\n",fileName)	,	exit(0);
#endif // NEW_PNG_IO

	info->png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,0,0,0);
	if(!info->png_ptr)	fprintf(stderr,"Failed to create png write struct\n")	,	exit(0);
	info->info_ptr = png_create_info_struct(info->png_ptr);
	if(!info->info_ptr)	fprintf(stderr,"Failed to create png info struct\n")	,	exit(0);

#if NEW_PNG_IO
	png_set_write_fn( info->png_ptr , info->outStream , myPngWriteFunction , myPngFlushFunction );
#else // !NEW_JOEG_IO
	png_init_io(info->png_ptr, info->fp);
#endif // NEW_PNG_IO


	png_set_compression_level(info->png_ptr,Z_BEST_SPEED);

	png_set_IHDR(info->png_ptr, info->info_ptr, width, height,
		16,PNG_COLOR_TYPE_RGB,
		PNG_INTERLACE_NONE,
		PNG_COMPRESSION_TYPE_DEFAULT,
		PNG_FILTER_TYPE_DEFAULT);
	png_write_info(info->png_ptr, info->info_ptr);

	{
		long int a = 1;
		int swap = (*((unsigned char *) &a) == 1);
		if(swap)	png_set_swap(info->png_ptr);
	}
	return info;
}
void* PNGInitWriteGray(char* fileName,int width,int height,int quality , MultiStreamIOServer* ioServer )
{
	PNGWriteInfo* info = (PNGWriteInfo*)malloc(sizeof(PNGWriteInfo));
#if NEW_PNG_IO
	info->outStream = new VariableIOClientStream( );
	info->outStream->Initialize( fileName , PNG_ZBUF_SIZE , 2 , false );
	info->outStream->SetServer( ioServer );
#else // !NEW_PNG_IO
	info->fp = fopen(fileName,"wb");
	if(!info->fp)	fprintf(stderr,"Failed to open: %s\n",fileName)	,	exit(0);
#endif // NEW_PNG_IO

	info->png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,0,0,0);
	if(!info->png_ptr)
	{
		fprintf(stderr,"Failed to create png write struct\n");
		exit(0);
	}
	info->info_ptr = png_create_info_struct(info->png_ptr);
	if(!info->info_ptr)
	{
		fprintf(stderr,"Failed to create png info struct\n");
		exit(0);
	}
#if NEW_PNG_IO
	png_set_write_fn( info->png_ptr , info->outStream , myPngWriteFunction , myPngFlushFunction );
#else // !NEW_JOEG_IO
	png_init_io(info->png_ptr, info->fp);
#endif // NEW_PNG_IO

	png_set_IHDR(info->png_ptr, info->info_ptr, width, height,
		8,PNG_COLOR_TYPE_GRAY,
		PNG_INTERLACE_NONE,
		PNG_COMPRESSION_TYPE_DEFAULT,
		PNG_FILTER_TYPE_DEFAULT);
	png_write_info(info->png_ptr, info->info_ptr);

	{
		long int a = 1;
		int swap = (*((unsigned char *) &a) == 1);
		if(swap)	png_set_swap(info->png_ptr);
	}
	return info;
}
void* PNGInitWriteGrayHDR(char* fileName,int width,int height,int quality , MultiStreamIOServer* ioServer )
{
	PNGWriteInfo* info = (PNGWriteInfo*)malloc(sizeof(PNGWriteInfo));
#if NEW_PNG_IO
	info->outStream = new VariableIOClientStream( );
	info->outStream->Initialize( fileName , PNG_ZBUF_SIZE , 2 , false );
	info->outStream->SetServer( ioServer );
#else // !NEW_PNG_IO
	info->fp = fopen(fileName,"wb");
	if(!info->fp)	fprintf(stderr,"Failed to open: %s\n",fileName)	,	exit(0);
#endif // NEW_PNG_IO

	info->png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,0,0,0);
	if(!info->png_ptr)
	{
		fprintf(stderr,"Failed to create png write struct\n");
		exit(0);
	}
	info->info_ptr = png_create_info_struct(info->png_ptr);
	if(!info->info_ptr)
	{
		fprintf(stderr,"Failed to create png info struct\n");
		exit(0);
	}
#if NEW_PNG_IO
	png_set_write_fn( info->png_ptr , info->outStream , myPngWriteFunction , myPngFlushFunction );
#else // !NEW_JOEG_IO
	png_init_io(info->png_ptr, info->fp);
#endif // NEW_PNG_IO

	png_set_IHDR(info->png_ptr, info->info_ptr, width, height,
		16,PNG_COLOR_TYPE_GRAY,
		PNG_INTERLACE_NONE,
		PNG_COMPRESSION_TYPE_DEFAULT,
		PNG_FILTER_TYPE_DEFAULT);
	png_write_info(info->png_ptr, info->info_ptr);

	{
		long int a = 1;
		int swap = (*((unsigned char *) &a) == 1);
		if(swap)	png_set_swap(info->png_ptr);
	}
	return info;
}
void PNGWriteRow(void* pixels,void* v,int j)
{
	unsigned char* _pixels = (unsigned char*)pixels;
	PNGWriteInfo* info = (PNGWriteInfo*)v;
	png_bytep row_pointer=png_bytep(&_pixels[0]);
	png_write_row(info->png_ptr, row_pointer);
}
void PNGWriteRowHDR(void* pixels,void* v,int j)
{
	unsigned __int16* _pixels = (unsigned __int16*)pixels;
	PNGWriteInfo* info = (PNGWriteInfo*)v;
	png_bytep row_pointer=png_bytep(&_pixels[0]);
	png_write_row(info->png_ptr, row_pointer);
}
void PNGFinalizeWrite(void* v)
{
	PNGWriteInfo* info = (PNGWriteInfo*)v;
	png_write_end(info->png_ptr, NULL);
	png_destroy_write_struct(&info->png_ptr, &info->info_ptr);
#if NEW_PNG_IO
	delete info->outStream;
#else // !NEW_PNG_IO
	fclose(info->fp);
#endif // NEW_PNG_IO
	free(info);
}
void PNGGetImageSize( char* fileName , int& width , int& height )
{
//	png_structp png_ptr;
//	png_infop info_ptr,end_info;
//	bool is16bit;

	FILE* fp;
	png_structp png_ptr;
	png_infop info_ptr , end_info;

	fp = fopen( fileName , "rb" );
	if(!fp)	fprintf(stderr,"Failed to open: %s\n",fileName)	,	exit(0);

	png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING,0,0,0);
	if(!png_ptr)
	{
		fprintf(stderr,"Failed to create PNG read structure\n");
		exit(0);
	}
	info_ptr = png_create_info_struct( png_ptr );
	if(!info_ptr)
	{
		fprintf(stderr,"Failed to create PNG info structure 1\n");
		exit(0);
	}
	end_info = png_create_info_struct(png_ptr);
	if(!end_info)	
	{
		fprintf(stderr,"Failed to create PNG info structure 2\n");
		exit(0);
	}
	png_init_io( png_ptr, fp );
	png_read_info( png_ptr, info_ptr );
	width =  png_get_image_width ( png_ptr , info_ptr );
	height = png_get_image_height( png_ptr , info_ptr );
	png_destroy_read_struct( &png_ptr , &info_ptr , &end_info );
	fclose( fp );
}


void* PNGInitReadColor( char* fileName , int& width , int& height , bool alpha , MultiStreamIOServer* ioServer )
{
	PNGReadInfo* info = (PNGReadInfo*)malloc(sizeof(PNGReadInfo));
#if NEW_PNG_IO
	info->inStream = new VariableIOClientStream( );
	info->inStream->Initialize( fileName , PNG_ZBUF_SIZE , 2 , true );
	info->inStream->SetServer( ioServer );
#else // !NEW_PNG_IO
	info->fp = fopen(fileName,"rb");
	if(!info->fp)	fprintf(stderr,"Failed to open: %s\n",fileName)	,	exit(0);
#endif // NEW_PNG_IO

	info->png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING,0,0,0);
	if(!info->png_ptr)
	{
		fprintf(stderr,"Failed to create PNG read structure\n");
		exit(0);
	}
	info->info_ptr = png_create_info_struct(info->png_ptr);
	if(!info->info_ptr)
	{
		fprintf(stderr,"Failed to create PNG info structure 1\n");
		exit(0);
	}
	info->end_info = png_create_info_struct(info->png_ptr);
	if(!info->end_info)	
	{
		fprintf(stderr,"Failed to create PNG info structure 2\n");
		exit(0);
	}
#if NEW_PNG_IO
	png_set_read_fn( info->png_ptr , info->inStream , myPngReadFunction );
#else // !NEW_JOEG_IO
	png_init_io(info->png_ptr, info->fp);
#endif // NEW_PNG_IO

	png_read_info(info->png_ptr, info->info_ptr);
	info->width=width=png_get_image_width(info->png_ptr,info->info_ptr);
	height=png_get_image_height(info->png_ptr,info->info_ptr);
	int ncomp=png_get_channels(info->png_ptr,info->info_ptr);
	int bit_depth=png_get_bit_depth(info->png_ptr,info->info_ptr);
	int color_type= png_get_color_type(info->png_ptr,info->info_ptr);
	if(width<=0 || height<=0)				exit(0);
	if(ncomp<1 || ncomp>4)					exit(0);
	if(bit_depth!=8 && bit_depth!=16)		exit(0);
	if(color_type==PNG_COLOR_TYPE_PALETTE)	png_set_expand(info->png_ptr)	,	printf("Expanding PNG color pallette\n");
	if (color_type == PNG_COLOR_TYPE_GRAY || color_type == PNG_COLOR_TYPE_GRAY_ALPHA)	png_set_gray_to_rgb(info->png_ptr); // always RGB
	if( alpha )
	{
		if( color_type!=PNG_COLOR_TYPE_RGBA ) return NULL;
		info->hasAlpha = true;
	}
	else
	{
		png_set_strip_alpha( info->png_ptr );
		info->hasAlpha = false;
	}
	//?
	info->is16bit = (bit_depth==16);
	if( alpha )
	{
		if(!info->is16bit)	info->pixels = (unsigned char*) malloc ( width*4 );
		else				info->pixels = NULL;
	}
	else
	{
		if(!info->is16bit)	info->pixels = (unsigned char*) malloc ( width*3 );
		else				info->pixels = NULL;
	}
	{
		long int a = 1;
		int swap = (*((unsigned char *) &a) == 1);
		if(swap)	png_set_swap(info->png_ptr);
	}
	return info;
}
void* PNGInitReadColor     ( char* fileName , int& width , int& height , MultiStreamIOServer* ioServer ) { return PNGInitReadColor( fileName , width , height , false , ioServer ); }
void* PNGInitReadColorAlpha( char* fileName , int& width , int& height , MultiStreamIOServer* ioServer ) { return PNGInitReadColor( fileName , width , height , true  , ioServer ); }
void* PNGInitReadGray(char* fileName,int& width,int& height , MultiStreamIOServer* ioServer )
{
	PNGReadInfo* info = (PNGReadInfo*)malloc(sizeof(PNGReadInfo));

#if NEW_PNG_IO
	info->inStream = new VariableIOClientStream( );
	info->inStream->Initialize( fileName , PNG_ZBUF_SIZE , 2 , true );
	info->inStream->SetServer( ioServer );
#else // !NEW_PNG_IO
	info->fp = fopen(fileName,"rb");
	if(!info->fp)	fprintf(stderr,"Failed to open: %s\n",fileName)	,	exit(0);
#endif // NEW_PNG_IO

	info->png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING,0,0,0);
	if(!info->png_ptr)
	{
		fprintf(stderr,"Failed to create PNG read structure\n");
		exit(0);
	}
	info->info_ptr = png_create_info_struct(info->png_ptr);
	if(!info->info_ptr)
	{
		fprintf(stderr,"Failed to create PNG info structure 1\n");
		exit(0);
	}
	info->end_info = png_create_info_struct(info->png_ptr);
	if(!info->end_info)	
	{
		fprintf(stderr,"Failed to create PNG info structure 2\n");
		exit(0);
	}

#if NEW_PNG_IO
	png_set_read_fn( info->png_ptr , info->inStream , myPngReadFunction );
#else // !NEW_JOEG_IO
	png_init_io(info->png_ptr, info->fp);
#endif // NEW_PNG_IO

	png_read_info(info->png_ptr, info->info_ptr);
	info->width=width=png_get_image_width(info->png_ptr,info->info_ptr);
	height=png_get_image_height(info->png_ptr,info->info_ptr);
	int ncomp=png_get_channels(info->png_ptr,info->info_ptr);
	int bit_depth=png_get_bit_depth(info->png_ptr,info->info_ptr);
	int color_type= png_get_color_type(info->png_ptr,info->info_ptr);
	if(width<=0 || height<=0)				exit(0);
	if(ncomp<1 || ncomp>4)					exit(0);
	if(bit_depth!=8 && bit_depth!=16)		exit(0);
	if(color_type==PNG_COLOR_TYPE_PALETTE)	png_set_expand(info->png_ptr)	,	printf("Expanding PNG color pallette\n");
	if (color_type == PNG_COLOR_TYPE_RGB || color_type == PNG_COLOR_TYPE_RGB_ALPHA)	png_set_rgb_to_gray(info->png_ptr,1,-1,-1); // always RGB
	png_set_strip_alpha(info->png_ptr);
	info->is16bit=(bit_depth==16);
	if(!info->is16bit)	info->pixels = (unsigned char*) malloc ( width*3);
	else				info->pixels = NULL;
	{
		long int a = 1;
		int swap = (*((unsigned char *) &a) == 1);
		if(swap)	png_set_swap(info->png_ptr);
	}
	return info;
}

const double CharToInt16 = double((1<<16)-1) / double((1<<8)-1);

void PNGReadColorRow(void* pixels,void* v,int j)
{
	PNGReadInfo* info = (PNGReadInfo*)v;
	if(info->is16bit)
	{
		png_bytep row=(png_bytep)pixels;
		png_read_row(info->png_ptr,row,NULL);
	}
	else
	{
		png_bytep row=(png_bytep)info->pixels;
		unsigned __int16 *p = (unsigned __int16*) pixels;
		png_read_row(info->png_ptr,row,NULL);
		for(int i=0;i<info->width;i++)
			if( info->hasAlpha )
				for( int c=0 ; c<4 ; c++ ) p[i*4+c] = __int16(double(info->pixels[i*4+c])*CharToInt16+0.5);
			else
				for( int c=0 ; c<3 ; c++ ) p[i*3+c] = __int16(double(info->pixels[i*3+c])*CharToInt16+0.5);

	}
}
void PNGReadGrayRow(void* pixels,void* v,int j)
{
	PNGReadInfo* info = (PNGReadInfo*)v;
	if(info->is16bit)
	{
		png_bytep row=(png_bytep)pixels;
		png_read_row(info->png_ptr,row,NULL);
	}
	else
	{
		png_bytep row=(png_bytep)info->pixels;
		unsigned __int16 *p = (unsigned __int16*) pixels;
		png_read_row(info->png_ptr,row,NULL);
		for(int i=0;i<info->width;i++)	p[i] = __int16(double(info->pixels[i])*CharToInt16+0.5);
	}
}

void PNGFinalizeRead(void* v)
{
	PNGReadInfo* info = (PNGReadInfo*)v;
	png_destroy_read_struct(&info->png_ptr, &info->info_ptr, &info->end_info);
	if(info->pixels)	free(info->pixels);
	info->pixels = NULL;
#if NEW_PNG_IO
	delete info->inStream;
#else // !NEW_PNG_IO
	fclose(info->fp);
#endif // NEW_PNG_IO
	free(info);
}
