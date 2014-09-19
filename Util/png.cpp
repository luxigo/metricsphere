#include "PNG/png.h"
#include "image.h"
#include <vector>

//static void read_png(Image& image, FILE* file)
int PNGReadImage(FILE *fp,Image32& img)
{
	png_structp png_ptr =
		png_create_read_struct(PNG_LIBPNG_VER_STRING,
		0, // (png_voidp)user_error_ptr
		0, // user_error_fn
		0  // user_warning_fn
		);
	if(!png_ptr)	return 0;
	png_infop info_ptr = png_create_info_struct(png_ptr);
	if(!info_ptr)	return 0;
	png_infop end_info = png_create_info_struct(png_ptr);
	if(!end_info)	return 0;
	png_init_io(png_ptr, fp);

	// tell the library if we have already read any bytes from header
	// png_set_sig_bytes(png_ptr, 0);
	// callback to handle user chunk data
	// png_set_read_user_chunk_fn(png_ptr,user_chunk_ptr,read_chunk_callback);
	// callback used to control a progress meter
	// png_set_read_status_fn(png_ptr, read_row_callback);
	//
	if (1) {                    // high-level read
		int png_transforms=(PNG_TRANSFORM_STRIP_16  | // 16-bit to 8-bit
			PNG_TRANSFORM_PACKING   | // expand 1, 2, and 4-bit
			0);
		png_read_png(png_ptr, info_ptr, png_transforms, NULL);
		int width=png_get_image_width(png_ptr,info_ptr);
		int height=png_get_image_height(png_ptr,info_ptr);
		int ncomp=png_get_channels(png_ptr,info_ptr);
		int bit_depth=png_get_bit_depth(png_ptr,info_ptr);
		int color_type= png_get_color_type(png_ptr,info_ptr);
		if(width<=0 || height<=0)				return 0;
		if(ncomp<1 || ncomp>4)					return 0;
		if(bit_depth!=8)						return 0;
//		if(color_type==PNG_COLOR_TYPE_PALETTE)	return 0;
		if(color_type==PNG_COLOR_TYPE_PALETTE && 0)
		{
			png_set_palette_to_rgb(png_ptr);
			png_set_expand(png_ptr);
		}
		img.setSize(width,height);
		png_bytep* row_pointers;   // [height]
		row_pointers = png_get_rows(png_ptr, info_ptr);
		for(int y=0;y<img.height();y++)
		{
			unsigned char* buf=(unsigned char*)(row_pointers[y]);
			for(int x=0;x<img.width();x++)
				for(int z=0;z<ncomp;z++)
					if(color_type==PNG_COLOR_TYPE_PALETTE)
					{
						png_color clr=info_ptr->palette[*buf++];
						img(x,y).r=clr.red;
						img(x,y).g=clr.green;
						img(x,y).b=clr.blue;
						img(x,y).a=255;
					}
					else									img(x,y)[z]=*buf++;
		}
	} else {                    // lower-level read, directly into image.
		png_read_info(png_ptr, info_ptr);
		int width=png_get_image_width(png_ptr,info_ptr);
		int height=png_get_image_height(png_ptr,info_ptr);
		int ncomp=png_get_channels(png_ptr,info_ptr);
		int bit_depth=png_get_bit_depth(png_ptr,info_ptr);
		int color_type= png_get_color_type(png_ptr,info_ptr);
		if(width<=0 || height<=0)				return 0;
		if(ncomp<1 || ncomp>4)					return 0;
		if(bit_depth!=8)						return 0;
		if(color_type==PNG_COLOR_TYPE_PALETTE)	return 0;
		img.setSize(width,height);
		if (bit_depth==16) png_set_strip_16(png_ptr);
		if (bit_depth<8) png_set_packing(png_ptr);
		if (color_type == PNG_COLOR_TYPE_GRAY || color_type == PNG_COLOR_TYPE_GRAY_ALPHA)
			png_set_gray_to_rgb(png_ptr); // always RGB
		png_set_filler(png_ptr, 0, PNG_FILLER_AFTER); // always RGBA
		std::vector<png_bytep> row_pointers;
		row_pointers.resize(img.height());
		for(int y=0;y<img.height();y++)	row_pointers[y]=(unsigned char*)&img(0,y);
		png_read_image(png_ptr, &row_pointers[0]);
	}
	png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
	return 1;
}
int PNGWriteImage(Image32& img,FILE* fp)
{
	png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,0,0,0);
	if(!png_ptr)	return 0;
	png_infop info_ptr = png_create_info_struct(png_ptr);
	if(!info_ptr)	return 0;
	png_init_io(png_ptr, fp);
	// turn off compression or set another filter
	// png_set_filter(png_ptr, 0, PNG_FILTER_NONE);
	png_set_IHDR(png_ptr, info_ptr, img.width(), img.height(),
		8,PNG_COLOR_TYPE_RGB_ALPHA,
		PNG_INTERLACE_NONE,
		PNG_COMPRESSION_TYPE_DEFAULT,
		PNG_FILTER_TYPE_DEFAULT);
	if (0) {                    // high-level write
		std::vector<unsigned char> matrix(img.width()*img.height()*4);
		std::vector<png_bytep> row_pointers(img.height());
		for(int y=0;y<img.height();y++)
		{
			row_pointers[y]=&matrix[y*img.width()*4];
			unsigned char* buf=&matrix[y*img.width()*4];
			for(int x=0;x<img.width();x++)
				for(int z=0;z<4;z++)
					*buf++=img(x,y)[z];
		}
		png_set_rows(png_ptr, info_ptr, &row_pointers[0]);
		int png_transforms=0;
		png_write_png(png_ptr, info_ptr, png_transforms, NULL);
	} else {                    // low-level write
		png_write_info(png_ptr, info_ptr);
		// png_set_filler(png_ptr, 0, PNG_FILLER_AFTER);
		//  but no way to provide GRAY data with RGBA fill, so pack each row
		std::vector<unsigned char> buffer(img.width()*4);
		for(int y=0;y<img.height();y++)
		{
			unsigned char* buf=&buffer[0];
			for(int x=0;x<img.width();x++)
				for(int z=0;z<4;z++)
					*buf++=img(x,y)[z];
			png_bytep row_pointer=&buffer[0];
			png_write_row(png_ptr, row_pointer);
		}
	}
	png_write_end(png_ptr, NULL);
	png_destroy_write_struct(&png_ptr, &info_ptr);
	return 1;
}
int PNGReadImage(char* fileName,Image32& img)
{
	FILE* fp=fopen(fileName,"rb");
	if(!fp)	return 0;
	int ret=PNGReadImage(fp,img);
	fclose(fp);
	return ret;
}
int PNGWriteImage(Image32& img,char* fileName)
{
	FILE* fp=fopen(fileName,"wb");
	if(!fp)	return 0;
	int ret=PNGWriteImage(img,fp);
	fclose(fp);
	return ret;
}