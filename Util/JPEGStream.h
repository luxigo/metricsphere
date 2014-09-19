#ifndef JPEG_STREAM_INCLUDED
#define JPEG_STREAM_INCLUDED
#include <stdio.h>
#include <stdlib.h>
#include <setjmp.h>
#include <windows.h>

#include "JPEG/jpeglib.h"
#include "JPEG/jerror.h"
#include "JPEG/jinclude.h"

#include <Util/MultiStreamIO.h>
#if NEW_JPEG_IO

typedef struct
{
	struct jpeg_destination_mgr pub;	/* public fields */
	VariableIOClientStream* outStream;	/* target stream */
	JOCTET * buffer;					/* start of buffer */
} misha_destination_mgr;
typedef misha_destination_mgr * misha_dest_ptr;

typedef struct
{
	struct jpeg_source_mgr pub;			/* public fields */
	VariableIOClientStream* inStream; 	/* source stream */
	JOCTET * buffer;					/* start of buffer */
	boolean start_of_file;				/* have we gotten any data yet? */
} misha_source_mgr;

typedef misha_source_mgr * misha_src_ptr;


#endif // NEW_JPEG_IO

////////////////////////////////////////////////////////////////////
// Code adapted from example.c provided with libjpeg distribution //
////////////////////////////////////////////////////////////////////

struct my_error_mgr
{
	struct jpeg_error_mgr pub;    // "public" fields
	jmp_buf setjmp_buffer;        // for return to caller
};
typedef struct my_error_mgr * my_error_ptr;


struct JPEGReadInfo
{
#if NEW_JPEG_IO
	VariableIOClientStream* inStream;
#else // !NEW_JPEG_IO
	FILE* fp;
#endif // NEW_JPEG_IO
	struct jpeg_decompress_struct cInfo;
	struct my_error_mgr	jErr;
};
struct JPEGWriteInfo
{
#if NEW_JPEG_IO
	VariableIOClientStream* outStream;
#else // !NEW_JPEG_IO
	FILE* fp;
#endif // NEW_JPEG_IO
	struct jpeg_compress_struct cInfo;
	struct jpeg_error_mgr jErr;
};
#endif // JPEG_STREAM_INCLUDED
