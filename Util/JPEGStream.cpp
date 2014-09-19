#include <stdio.h>
#include <stdlib.h>
#include <setjmp.h>
#include <windows.h>

#include "JPEG/jpeglib.h"
#include "JPEG/jerror.h"
#include "JPEG/jinclude.h"
#include "JPEGStream.h"

#if NEW_JPEG_IO

////////////////////////////////////////////////////////////////////
// Code adapted from example.c provided with libjpeg distribution //
////////////////////////////////////////////////////////////////////

////////////////////
// Info for input //
////////////////////

#define INPUT_BUF_SIZE  4096

METHODDEF(void)
init_source (j_decompress_ptr cinfo)
{
	misha_src_ptr src = (misha_src_ptr) cinfo->src;
	src->start_of_file = TRUE;
}
METHODDEF(boolean)
fill_input_buffer (j_decompress_ptr cinfo)
{
	misha_src_ptr src = (misha_src_ptr) cinfo->src;
	size_t nbytes;

	nbytes = src->inStream->read( src->buffer , INPUT_BUF_SIZE );

	if (nbytes <= 0)
	{
		if (src->start_of_file)	/* Treat empty input file as fatal error */
			ERREXIT(cinfo, JERR_INPUT_EMPTY);
		WARNMS(cinfo, JWRN_JPEG_EOF);
		/* Insert a fake EOI marker */
		src->buffer[0] = (JOCTET) 0xFF;
		src->buffer[1] = (JOCTET) JPEG_EOI;
		nbytes = 2;
	}

	src->pub.next_input_byte = src->buffer;
	src->pub.bytes_in_buffer = nbytes;
	src->start_of_file = FALSE;
	return TRUE;
}



METHODDEF(void)
skip_input_data (j_decompress_ptr cinfo, long num_bytes)
{
	misha_src_ptr src = (misha_src_ptr) cinfo->src;

	if (num_bytes > 0)
	{
		while (num_bytes > (long) src->pub.bytes_in_buffer)
		{
			num_bytes -= (long) src->pub.bytes_in_buffer;
			(void) fill_input_buffer(cinfo);
		}
		src->pub.next_input_byte += (size_t) num_bytes;
		src->pub.bytes_in_buffer -= (size_t) num_bytes;
	}
}

METHODDEF(void)
term_source (j_decompress_ptr cinfo)
{
}

GLOBAL(void)
jpeg_stdio_src (j_decompress_ptr cinfo, VariableIOClientStream * inStream )
{
	misha_src_ptr src;
	if (cinfo->src == NULL)
	{

		cinfo->src = (struct jpeg_source_mgr *)	(*cinfo->mem->alloc_small) ((j_common_ptr) cinfo, JPOOL_PERMANENT, SIZEOF(misha_source_mgr));
		src = (misha_src_ptr) cinfo->src;
		src->buffer = (JOCTET *) (*cinfo->mem->alloc_small) ((j_common_ptr) cinfo, JPOOL_PERMANENT,	INPUT_BUF_SIZE * SIZEOF(JOCTET));
	}

	src = (misha_src_ptr) cinfo->src;
	src->pub.init_source = init_source;
	src->pub.fill_input_buffer = fill_input_buffer;
	src->pub.skip_input_data = skip_input_data;
	src->pub.resync_to_restart = jpeg_resync_to_restart; /* use default method */
	src->pub.term_source = term_source;
	src->inStream = inStream;
	src->pub.bytes_in_buffer = 0; /* forces fill_input_buffer on first read */
	src->pub.next_input_byte = NULL; /* until buffer loaded */
}

/////////////////////
// Info for output //
/////////////////////
#define OUTPUT_BUF_SIZE  4096

METHODDEF(void)
init_destination ( j_compress_ptr cinfo )
{
	misha_dest_ptr dest = (misha_dest_ptr) cinfo->dest;

	/* Allocate the output buffer --- it will be released when done with image */
	dest->buffer = (JOCTET *) (*cinfo->mem->alloc_small) ((j_common_ptr) cinfo, JPOOL_IMAGE, OUTPUT_BUF_SIZE * SIZEOF(JOCTET));
	dest->pub.next_output_byte = dest->buffer;
	dest->pub.free_in_buffer = OUTPUT_BUF_SIZE;
}


METHODDEF(boolean)
empty_output_buffer (j_compress_ptr cinfo)
{
	misha_dest_ptr dest = (misha_dest_ptr) cinfo->dest;

	if( dest->outStream->write( dest->buffer , OUTPUT_BUF_SIZE ) != (size_t) OUTPUT_BUF_SIZE ) ERREXIT(cinfo, JERR_FILE_WRITE);
	dest->pub.next_output_byte = dest->buffer;
	dest->pub.free_in_buffer = OUTPUT_BUF_SIZE;

	return TRUE;
}

METHODDEF(void)
term_destination (j_compress_ptr cinfo)
{
	misha_dest_ptr dest = (misha_dest_ptr) cinfo->dest;
	size_t datacount = OUTPUT_BUF_SIZE - dest->pub.free_in_buffer;

	/* Write any data remaining in the buffer */
	if (datacount > 0) if( dest->outStream->write( dest->buffer , datacount ) != datacount ) ERREXIT(cinfo, JERR_FILE_WRITE);
	//  fflush(dest->outfile);
	/* Make sure we wrote the output file OK */
	// if (ferror(dest->outfile))
	//  ERREXIT(cinfo, JERR_FILE_WRITE);
}

GLOBAL(void)
jpeg_stdio_dest (j_compress_ptr cinfo, VariableIOClientStream * outStream )
{
	misha_dest_ptr dest;
	if (cinfo->dest == NULL)/* first time for this JPEG object? */
		cinfo->dest = (struct jpeg_destination_mgr *) (*cinfo->mem->alloc_small) ((j_common_ptr) cinfo, JPOOL_PERMANENT, SIZEOF(misha_destination_mgr));

	dest = (misha_dest_ptr) cinfo->dest;
	dest->pub.init_destination = init_destination;
	dest->pub.empty_output_buffer = empty_output_buffer;
	dest->pub.term_destination = term_destination;
	dest->outStream = outStream;
}
#endif // NEW_JPEG_IO
METHODDEF(void)
my_error_exit (j_common_ptr cinfo)
{
	// cinfo->err really points to a my_error_mgr struct, so coerce pointer
	my_error_ptr myerr = (my_error_ptr) cinfo->err;

	// Always display the message.
	// We could postpone this until after returning, if we chose.
	(*cinfo->err->output_message) (cinfo);

	// Return control to the setjmp point
	longjmp(myerr->setjmp_buffer, 1);
}

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

void* JPEGInitWriteColor(char* fileName,int width,int height,int quality , MultiStreamIOServer* ioServer )
{
	JPEGWriteInfo *info=(JPEGWriteInfo*)malloc(sizeof(JPEGWriteInfo));

#if NEW_JPEG_IO
	info->outStream = new VariableIOClientStream( );
	info->outStream->Initialize( fileName , OUTPUT_BUF_SIZE , 2 , false );
	info->outStream->SetServer( ioServer );
#else // !NEW_JPEG_IO
	info->fp=fopen(fileName,"wb");
	if(!info->fp)	fprintf(stderr,"Failed to open: %s\n",fileName) , exit(0);
#endif // NEW_JPEG_IO

	info->cInfo.err=jpeg_std_error(&info->jErr);
	jpeg_create_compress(&info->cInfo);

#if NEW_JPEG_IO
	jpeg_stdio_dest( &info->cInfo , info->outStream );
#else // !NEW_JOEG_IO
	jpeg_stdio_dest( &info->cInfo , info->fp );
#endif // NEW_JPEG_IO

	info->cInfo.image_width = width;    /* image width and height, in pixels */
	info->cInfo.image_height = height;
	info->cInfo.input_components = 3;           /* # of color components per pixel */
	info->cInfo.in_color_space = JCS_RGB;       /* colorspace of input image */

	jpeg_set_defaults(&info->cInfo);
	jpeg_set_quality(&info->cInfo, quality, TRUE);

	jpeg_start_compress(&info->cInfo, TRUE);
	return info;
}

void* JPEGInitWriteGray(char* fileName,int width,int height,int quality , MultiStreamIOServer* ioServer )
{
	JPEGWriteInfo *info=(JPEGWriteInfo*)malloc(sizeof(JPEGWriteInfo));

#if NEW_JPEG_IO
	info->outStream = new VariableIOClientStream( );
	info->outStream->Initialize( fileName , OUTPUT_BUF_SIZE , 2 , false );
	info->outStream->SetServer( ioServer );
#else // !NEW_JPEG_IO
	info->fp=fopen(fileName,"wb");
	if(!info->fp)	fprintf(stderr,"Failed to open: %s\n",fileName) , exit(0);
#endif // NEW_JPEG_IO

	info->cInfo.err=jpeg_std_error(&info->jErr);
	jpeg_create_compress(&info->cInfo);

#if NEW_JPEG_IO
	jpeg_stdio_dest( &info->cInfo , info->outStream );
#else // !NEW_JOEG_IO
	jpeg_stdio_dest(&info->cInfo,info->fp);
#endif // NEW_JPEG_IO

	info->cInfo.image_width = width;			/* image width and height, in pixels */
	info->cInfo.image_height = height;
	info->cInfo.input_components = 1;			/* # of color components per pixel */
	info->cInfo.in_color_space = JCS_GRAYSCALE;	/* colorspace of input image */

	jpeg_set_defaults(&info->cInfo);
	jpeg_set_quality(&info->cInfo, quality, TRUE);

	jpeg_start_compress(&info->cInfo, TRUE);
	return info;
}
void JPEGWriteRow(void* pixels,void* v,int j)
{
	unsigned char* _pixels = (unsigned char*) pixels;
	JPEGWriteInfo* info=(JPEGWriteInfo*)v;
	JSAMPROW row_pointer[1];
	row_pointer[0] = _pixels;
	(void) jpeg_write_scanlines(&info->cInfo, row_pointer, 1);
}
void JPEGFinalizeWrite(void* v)
{
	JPEGWriteInfo* info=(JPEGWriteInfo*)v;
	jpeg_finish_compress(&info->cInfo);
	jpeg_destroy_compress(&info->cInfo);
#if NEW_JPEG_IO
	delete info->outStream;
#else // !NEW_JPEG_IO
	fclose(info->fp);
#endif // NEW_JPEG_IO
	free(info);
}

void JPEGGetImageSize( char* fileName , int& width , int& height )
{
	FILE* fp;
	struct jpeg_decompress_struct cInfo;
	struct my_error_mgr	jErr;

	fp = fopen( fileName , "rb" );
	if(!fp)	fprintf(stderr,"Failed to open: %s\n",fileName) , exit(0);

	cInfo.err = jpeg_std_error( &jErr.pub );
	jErr.pub.error_exit = my_error_exit;
	if( setjmp( jErr.setjmp_buffer ) )
	{
		jpeg_destroy_decompress( &cInfo );
		fprintf( stderr , "JPEG error occured\n" );
		return;
	}
	jpeg_create_decompress( &cInfo );
	jpeg_stdio_src( &cInfo, fp );	
	(void) jpeg_read_header( &cInfo , TRUE );
	width  = cInfo.image_width;
	height = cInfo.image_height;
	jpeg_destroy_decompress( &cInfo );
	fclose( fp );
}

void* JPEGInitReadColor(char* fileName,int& width,int& height , MultiStreamIOServer* ioServer )
{
	JPEGReadInfo *info=(JPEGReadInfo*)malloc(sizeof(JPEGReadInfo));

#if NEW_JPEG_IO
	info->inStream = new VariableIOClientStream( );
	info->inStream->Initialize( fileName , OUTPUT_BUF_SIZE , 2 , true );
	info->inStream->SetServer( ioServer );
#else // !NEW_JPEG_IO
	info->fp=fopen(fileName,"rb");
	if(!info->fp)	fprintf(stderr,"Failed to open: %s\n",fileName) , exit(0);
#endif // NEW_JPEG_IO

	info->cInfo.err = jpeg_std_error(&info->jErr.pub);
	info->jErr.pub.error_exit = my_error_exit;
	if (setjmp(info->jErr.setjmp_buffer))
	{
		jpeg_destroy_decompress(&info->cInfo);
		fprintf(stderr,"JPEG error occured\n");
		return 0;
	}

	jpeg_create_decompress(&info->cInfo);
#if NEW_JPEG_IO
	jpeg_stdio_src( &info->cInfo , info->inStream );
#else // !NEW_JPEG_IO
	jpeg_stdio_src(&info->cInfo, info->fp);
#endif // NEW_JPEG_IO

	(void) jpeg_read_header(&info->cInfo, TRUE);
	(void) jpeg_start_decompress(&info->cInfo);

	if(info->cInfo.output_components!=3)
	{
		fprintf(stderr,"Only 3 components per pixel supported: %d != 3\n",info->cInfo.output_components);
		return NULL;
	}
	width=info->cInfo.output_width;
	height=info->cInfo.output_height;

	return info;
}
void* JPEGInitReadGray(char* fileName,int& width,int& height , MultiStreamIOServer* ioServer )
{
	JPEGReadInfo *info=(JPEGReadInfo*)malloc(sizeof(JPEGReadInfo));
#if NEW_JPEG_IO
	info->inStream = new VariableIOClientStream( );
	info->inStream->Initialize( fileName , OUTPUT_BUF_SIZE , 2 , true );
	info->inStream->SetServer( ioServer );
#else // !NEW_JPEG_IO
	info->fp=fopen(fileName,"rb");
	if(!info->fp)	fprintf(stderr,"Failed to open: %s\n",fileName) , exit(0);
#endif // NEW_JPEG_IO

	info->cInfo.err = jpeg_std_error(&info->jErr.pub);
	info->jErr.pub.error_exit = my_error_exit;
	if (setjmp(info->jErr.setjmp_buffer))
	{
		jpeg_destroy_decompress(&info->cInfo);
		fprintf(stderr,"JPEG error occured\n");
		return 0;
	}


	jpeg_create_decompress(&info->cInfo);
#if NEW_JPEG_IO
	jpeg_stdio_src( &info->cInfo , info->inStream );
#else // !NEW_JPEG_IO
	jpeg_stdio_src(&info->cInfo, info->fp);
#endif // NEW_JPEG_IO


	(void) jpeg_read_header(&info->cInfo, TRUE);
	(void) jpeg_start_decompress(&info->cInfo);

	if(info->cInfo.output_components!=1)
	{
		fprintf(stderr,"Only one component per pixel supported: %d != 1\n",info->cInfo.output_components);
		return NULL;
	}
	width=info->cInfo.output_width;
	height=info->cInfo.output_height;

	return info;
}
void JPEGReadRow(void* pixels,void* v,int j)
{
	unsigned char* _pixels = (unsigned char*) pixels;
	JPEGReadInfo* info=(JPEGReadInfo*)v;
	if(info->cInfo.output_scanline >= info->cInfo.output_height)
	{
		fprintf(stderr,"Trying to read beyond the end of the jpeg file: %d >= %d\n",info->cInfo.output_scanline,info->cInfo.output_height);
		exit(0);
	}
	JSAMPROW row_pointers[1];
	row_pointers[0]=_pixels;
	int r=jpeg_read_scanlines(&info->cInfo, row_pointers, 1);
}

void JPEGFinalizeRead(void* v)
{
	JPEGReadInfo* info = ( JPEGReadInfo* )v;
	// The JPEG reader will not allow us to close until we've actually read the scan-lines.
	// How annoying!!!
	if( info->cInfo.output_scanline < info->cInfo.output_height )
	{
		unsigned char* pixels = (unsigned char*) malloc( sizeof ( unsigned char) * 3 * info->cInfo.image_width );
		JSAMPROW row_pointers[1];
		row_pointers[0] = pixels;
		while( info->cInfo.output_scanline < info->cInfo.output_height ) jpeg_read_scanlines( &info->cInfo , row_pointers , 1 );
	}
	(void) jpeg_finish_decompress(&info->cInfo);
	jpeg_destroy_decompress(&info->cInfo);
#if NEW_JPEG_IO
	delete info->inStream;
#else // !NEW_JPEG_IO
	fclose(info->fp);
#endif // NEW_JPEG_IO
	free(info);
}
