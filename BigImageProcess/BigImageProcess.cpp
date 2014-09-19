#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Util/Half/half.h"
#include "Util/MultiStreamIO.h"
#include "Util/CmdLineParser.h"
#include "Util/ImageStream.h"
#include "LinearAlgebra/Vector.h"

#define SEPARATE_COLORS true

cmdLineString In( "in" ) , Out( "out" ) , TileExtension( "tileExt" );
cmdLineIntArray<4> Clip( "clip" );
cmdLineInt Quality( "quality" , 100 ) , Down( "down" ) , TileWidth( "tileWidth" , 8192 ) , TileHeight( "tileHeight" , 8192 );
cmdLineFloatArray<3> DCTermOffset( "dcOffset" );
cmdLineReadable Progress( "progress" ) , HDR( "hdr" ) , GammaCorrection( "gCorrect" ) , Clamp( "clamp" );
cmdLineFloat Gamma( "gamma" , 1.0 ),Brightness( "brighten" , 1.0 );
cmdLineReadable* params[]=
{
	&In , &Out , &Down , &Quality , &Clip , &DCTermOffset , &Gamma , &Brightness , &TileWidth , &TileHeight , &Progress , &TileExtension , &HDR , &Clamp
};

void ShowUsage(char* ex)
{
	printf("Usage %s:\n",ex);
	printf( "\t--%s <input image>\n",In.name);
	printf( "\t[--%s <output image>]\n",Out.name);
	printf( "\t[--%s <2X down sampling passes>]\n",Down.name);
	printf( "\t[--%s <startX , startY , width , height >]\n",Clip.name);
	printf( "\t[--%s <JPEG compression quality>=%d]\n",Quality.name,Quality.value);
	printf( "\t[--%s <default output tile width>=%d]\n" , TileWidth.name , TileWidth.value );
	printf( "\t[--%s <default output tile height>=%d]\n" , TileHeight.name , TileHeight.value );
	printf( "\t[--%s <default output file extension>]\n" , TileExtension.name );
	printf( "\t[--%s <offset for the DC term>]\n",DCTermOffset.name);
	printf( "\t[--%s <gamma correction term>=%f]\n",Gamma.name,Gamma.value);
	printf( "\t[--%s <brightness multiplier>=%f]\n",Brightness.name,Brightness.value);
	printf( "\t\tPixel <- [ (Pixel + Offset) * Brightness ]^Gamma\n");
	printf( "\t[--%s]\n" , HDR.name );
	printf( "\t[--%s]\n" , Clamp.name );
	printf( "\t[--%s]\n" , Progress.name );
}
class MyImage
{
public:
	StreamingGrid* data;
	int width,height;

	MyImage(void)
	{
		data=NULL;
		width=height=0;
	}
	~MyImage(void)
	{
		if(data)	delete data;
		data=NULL;
		width=height=0;
	}
};

double DownStencil2[2][2]=
{
	{1,1},
	{1,1}
};
double DownStencil4[4][4]=
{
	{1,3,3,1},
	{3,9,9,3},
	{3,9,9,3},
	{1,3,3,1}
};
double DownStencil6[6][6]=
{
	{  1,  5, 10, 10,  5,  1},
	{  5, 25, 50, 50, 25,  5},
	{ 10, 50,100,100, 50, 10},
	{ 10, 50,100,100, 50, 10},
	{  5, 25, 50, 50, 25,  5},
	{  1,  5, 10, 10,  5,  1}
};
template< class Real >
void DownSample2(const MyImage& high,MyImage& low,double dcOffset[3],double gamma=1.0,double brightness=1.0)
{
	high.data->reset(true,1);
	low.data->reset(false,1);

	Real *_highRows[2],*highRows[2];
	for(int i=0;i<2;i++)	_highRows[i]=new Real[high.width*3];

	for(int j=0;j<low.height;j++)
	{
		if( Progress.set )
		{
			double r = double(j) / low.height * 100;
			printf("[%.1f%%]\r" , r );
		}

		Pointer( Real ) lowRow = ( Pointer( Real ) )(*low.data)[j];
		// Copy in the next two high rows.
		for(int l=0;l<2;l++)
		{
			int jj=2*j+l;
			if(jj>=0 && jj<high.height)
			{
				Pointer( Real ) highRow = ( Pointer( Real ) )(*high.data)[jj];
				if( SEPARATE_COLORS )
					for(int k=0;k<high.width;k++)
						for(int c=0;c<3;c++)
						{
							highRow[k+c*high.width]=(highRow[k+c*high.width]+dcOffset[c])*brightness;
							highRow[k+c*high.width]=pow(double(highRow[k+c*high.width]),gamma);
						}
				else
					for(int k=0;k<high.width;k++)
						for(int c=0;c<3;c++)
						{
							highRow[k*3+c]=(highRow[k*3+c]+dcOffset[c])*brightness;
							highRow[k*3+c]=pow(double(highRow[k*3+c]),gamma);
						}
				memcpy(_highRows[jj%2],highRow,sizeof(Real)*3*high.width);
			}
			high.data->advance();
		}
		for(int l=0;l<2;l++)
		{
			int jj=2*j+l;
			if(jj<0 || jj>=high.height)	highRows[l]=NULL;
			else						highRows[l]=_highRows[jj%2];
		}

		for(int c=0;c<3;c++)
			for(int i=0;i<low.width;i++)
			{
				double value=0,sum=0;
				for(int l=0;l<2;l++)
				{
					int jj=j*2+l;
					if(jj<0 || jj>=high.height)	continue;
					for(int k=0;k<2;k++)
					{
						int ii=2*i+k;
						if(ii<0 || ii>=high.width)	continue;
						if( SEPARATE_COLORS ) value+=highRows[l][c*high.width+ii]*DownStencil2[k][l];
						else                  value+=highRows[l][c+ii*3]*DownStencil2[k][l];
						sum+=DownStencil2[k][l];
					}
				}
				if( SEPARATE_COLORS ) lowRow[i+c*low.width]=value/sum;
				else                  lowRow[i*3+c]=value/sum;
			}
		low.data->advance();
	}
	for(int i=0;i<2;i++)	delete[] _highRows[i];
}
template<class Real>
void DownSample2(const MyImage& high,MyImage& low)
{
	double offset[3]={0,0,0};
	DownSample2<Real>(high,low,offset);
}

template<class Real>
void DownSample4(const MyImage& high,MyImage& low,double dcOffset[3],double gamma=1.0,double brightness=1.0)
{
	high.data->reset(true,1);
	low.data->reset(false,1);

	Real *_highRows[4],*highRows[4];
	for(int i=0;i<4;i++)	_highRows[i]=new Real[high.width*3];

	// Buffer the first row
	memcpy(_highRows[0],(*high.data)[0],sizeof(Real)*high.width*3);
	if( SEPARATE_COLORS )
		for(int k=0;k<high.width;k++)
			for(int c=0;c<3;c++)
			{
				_highRows[0][k+c*high.width]=(_highRows[0][k+c*high.width]+dcOffset[c])*brightness;
				_highRows[0][k+c*high.width]=pow(double(_highRows[0][k+c*high.width]),gamma);
			}
	else
		for(int k=0;k<high.width;k++)
			for(int c=0;c<3;c++)
			{
				_highRows[0][k*3+c]=(_highRows[0][k*3+c]+dcOffset[c])*brightness;
				_highRows[0][k*3+c]=pow(double(_highRows[0][k*3+c]),gamma);
			}

	high.data->advance();

	for(int j=0;j<low.height;j++)
	{
		if( Progress.set )
		{
			double r = double(j) / low.height * 100;
			printf("[%.1f%%]\r" , r );
		}

		Real* lowRow=(Real*)(*low.data)[j];
		// Copy in the next two high rows.
		for(int l=1;l<3;l++)
		{
			int jj=2*j+l;
			if(jj>=0 && jj<high.height)
			{
				memcpy(_highRows[jj%4],(Real*)(*high.data)[jj],sizeof(Real)*3*high.width);
				if( SEPARATE_COLORS )
					for(int k=0;k<high.width;k++)
						for(int c=0;c<3;c++)
						{
							_highRows[jj%4][k+c*high.width]=(_highRows[jj%4][k+c*high.width]+dcOffset[c])*brightness;
							_highRows[jj%4][k+c*high.width]=pow(double(_highRows[jj%4][k+c*high.width]),gamma);
						}
				else
					for(int k=0;k<high.width;k++)
						for(int c=0;c<3;c++)
						{
							_highRows[jj%4][k*3+c]=(_highRows[jj%4][k*3+c]+dcOffset[c])*brightness;
							_highRows[jj%4][k*3+c]=pow(double(_highRows[jj%4][k*3+c]),gamma);
						}
			}
			high.data->advance();
		}
		for(int l=-1;l<3;l++)
		{
			int jj=2*j+l;
			if(jj<0 || jj>=high.height)	highRows[l+1]=NULL;
			else						highRows[l+1]=_highRows[jj%4];
		}

		for(int c=0;c<3;c++)
			for(int i=0;i<low.width;i++)
			{
				double value=0,sum=0;
				for(int l=-1;l<3;l++)
				{
					int jj=j*2+l;
					if(jj<0 || jj>=high.height)	continue;
					for(int k=-1;k<3;k++)
					{
						int ii=2*i+k;
						if(ii<0 || ii>=high.width)	continue;
						if( SEPARATE_COLORS ) value+=highRows[l+1][c*high.width+ii]*DownStencil4[k+1][l+1];
						else                  value+=highRows[l+1][c+ii*3]*DownStencil4[k+1][l+1];
						sum+=DownStencil4[k+1][l+1];
					}
				}
				if( SEPARATE_COLORS ) lowRow[i+c*low.width]=value/sum;
				else                  lowRow[i*3+c]=value/sum;
			}
		low.data->advance();
	}
	for(int i=0;i<4;i++)	delete[] _highRows[i];
}
template<class Real>
void DownSample4(const MyImage& high,MyImage& low)
{
	double offset[3]={0,0,0};
	DownSample4<Real>(high,low,offset);
}

template<class Real>
void DownSample6(const MyImage& high,MyImage& low,double dcOffset[3],double gamma=1.0,double brightness=1.0)
{
	high.data->reset(true,1);
	low.data->reset(false,1);

	Real *_highRows[6],*highRows[6];
	for(int i=0;i<6;i++)	_highRows[i]=new Real[high.width*3];

	// Buffer the first two row
	for(int j=0;j<2;j++)
	{
		Real* highRow=(Real*)(*high.data)[j];
		if( SEPARATE_COLORS )
			for(int i=0;i<high.width;i++)
				for(int c=0;c<3;c++)
				{
					_highRows[j][i+c*high.width]=(highRow[i+c*high.width]+dcOffset[c])*brightness;
					_highRows[j][i+c*high.width]=pow(double(_highRows[j][i+c*high.width]),gamma);
				}
		else
			for(int i=0;i<high.width;i++)
				for(int c=0;c<3;c++)
				{
					_highRows[j][i*3+c]=(highRow[i*3+c]+dcOffset[c])*brightness;
					_highRows[j][i*3+c]=pow(double(_highRows[j][i*3+c]),gamma);
				}
		high.data->advance();
	}

	for(int j=0;j<low.height;j++)
	{
		if( Progress.set )
		{
			double r = double(j) / low.height * 100;
			printf("[%.1f%%]\r" , r );
		}

		Real* lowRow=(Real*)(*low.data)[j];
		// Copy in the next two high rows.
		for(int l=2;l<4;l++)
		{
			int jj=2*j+l;
			if(jj>=0 && jj<high.height)
			{
				Real* highRow=(Real*)(*high.data)[jj];
				memcpy(_highRows[jj%6],highRow,sizeof(Real)*3*high.width);
				if( SEPARATE_COLORS )
					for(int k=0;k<high.width;k++)
						for(int c=0;c<3;c++)
						{
							_highRows[jj%6][k+c*high.width]=(_highRows[jj%6][k+c*high.width]+dcOffset[c])*brightness;
							_highRows[jj%6][k+c*high.width]=pow(double(_highRows[jj%6][k+c*high.width]),gamma);
						}
				else
					for(int k=0;k<high.width;k++)
						for(int c=0;c<3;c++)
						{
							_highRows[jj%6][k*3+c]=(_highRows[jj%6][k*3+c]+dcOffset[c])*brightness;
							_highRows[jj%6][k*3+c]=pow(double(_highRows[jj%6][k*3+c]),gamma);
						}
			}
			high.data->advance();
		}
		for(int l=-2;l<4;l++)
		{
			int jj=2*j+l;
			if(jj<0 || jj>=high.height)	highRows[l+2]=NULL;
			else						highRows[l+2]=_highRows[jj%6];
		}

		for(int c=0;c<3;c++)
			for(int i=0;i<low.width;i++)
			{
				double value=0,sum=0;
				for(int l=-2;l<4;l++)
				{
					int jj=j*2+l;
					if(jj<0 || jj>=high.height)	continue;
					for(int k=-2;k<4;k++)
					{
						int ii=2*i+k;
						if(ii<0 || ii>=high.width)	continue;
						if( SEPARATE_COLORS ) value+=highRows[l+2][c*high.width+ii]*DownStencil6[k+2][l+2];
						else                  value+=highRows[l+2][c+ii*3]*DownStencil6[k+2][l+2];
						sum+=DownStencil6[k+2][l+2];
					}
				}
				if( SEPARATE_COLORS ) lowRow[i+c*low.width]=value/sum;
				else                  lowRow[i*3+c]=value/sum;
			}
		low.data->advance();
	}
	for(int i=0;i<6;i++)	delete[] _highRows[i];
}
template<class Real>
void DownSample6(const MyImage& high,MyImage& low)
{
	double offset[3]={0,0,0};
	DownSample6<Real>(high,low,offset);
}

template< class Real , bool Remap >
void ClipWindow( const MyImage& in , MyImage& out , int startX , int startY , double dcOffset[3] , double gamma=1.0 , double brightness=1.0 )
{
	for( int j=0 ; j<in.height ; j++ )
	{
		if( Progress.set )
		{
			double r = double(j) / in.height * 100;
			printf("[%.1f%%]\r" , r );
		}

		Pointer( Real ) inRow = ( Pointer( Real ) )(*in.data)[j];
		if(j>=startY && j<startY+out.height)
		{
			Pointer( Real ) outRow = ( Pointer( Real ) )(*out.data)[j-startY];
			if( SEPARATE_COLORS )
			{
				for(int c=0;c<3;c++)
					for(int i=0;i<out.width;i++)
					{
						//inRow[i+startX+c*in.width] /= 256; // Added here because the old floating point output files were in the range [0,256] instead of [0,1]
#if 0
						double temp=inRow[i+startX+c*in.width];
						if(temp<1)	temp=1;
						outRow[i+c*out.width]=(temp+dcOffset[c])*brightness;
						outRow[i+c*out.width]=pow(double(outRow[i+c*out.width]),gamma);
#else
						if( Remap )
						{
							if( Clamp.set )
							{
								Real temp = (inRow[i+startX+c*in.width]+dcOffset[c])*brightness;
								if( temp<0 ) temp = 0;
								if( temp>1 ) temp = 1;
								outRow[i+c*out.width] = temp;
							}
							else outRow[i+c*out.width]=(inRow[i+startX+c*in.width]+dcOffset[c])*brightness;
							outRow[i+c*out.width]=pow(double(outRow[i+c*out.width]),gamma);
						}
						else
							if( Clamp.set )
							{
								Real temp = inRow[i+startX+c*in.width];
								if( temp<0 ) temp = 0;
								if( temp>1>>16 ) temp = 1>>16;
								outRow[i+c*out.width] = temp;
							}
							else outRow[i+c*out.width] = inRow[i+startX+c*in.width];
#endif
					}
			}
			else
			{
				if( Remap )
					for(int i=0;i<out.width;i++)
						for(int c=0;c<3;c++)
						{
						//inRow[i+startX+c*in.width] /= 256; // Added here because the old floating point output files were in the range [0,256] instead of [0,1]
#if 0
							double temp=inRow[(i+startX+in.width)*3 + c];
							if(temp<1)	temp=1;
							outRow[i*3+c]=(temp+dcOffset[c])*brightness;
							outRow[i*3+c]=pow(double(outRow[i*3+c]),gamma);
#else
							outRow[i*3+c]=(inRow[(i+startX)*3+c]+dcOffset[c])*brightness;
							outRow[i*3+c]=pow(double(outRow[i*3+c]),gamma);
						}
				else memcpy( outRow , inRow + startX*3 , sizeof( Real ) * 3 * in.width );
#endif
			}
			out.data->advance();
		}
		in.data->advance();
	}
}
template< class Real >
void ClipWindow(const MyImage& in,MyImage& out,int startX,int startY)
{
	double offset[3] = { 0 , 0 , 0 };
	ClipWindow< Real , false >(in,out,startX,startY,offset);
}
template<class Real>
void ImageInfo(const MyImage& in)
{
	double average[]={0,0,0};
	double lum,brightness=0,variance=0;
	Color<double> min,max;

	for(int j=0;j<in.height;j++)
	{
		if( Progress.set )
		{
			double r = double(j) / in.height * 100;
			printf("[%.1f%%]\r" , r );
		}
		Pointer( Real ) inRow=( Pointer( Real ) )(*in.data)[j];
		for(int c=0;c<3;c++)
			for(int i=0;i<in.width;i++)
			{
				average[c]+=inRow[i+c*in.width];
				variance+=inRow[i+c*in.width]*inRow[i+c*in.width];
				brightness+=log(double(inRow[i+c*in.width]));

				if((!i && !j) || inRow[i+c*in.width]<min[c])	min[c]=inRow[i+c*in.width];
				if((!i && !j) || inRow[i+c*in.width]>max[c])	max[c]=inRow[i+c*in.width];
			}
		in.data->advance();
	}
	variance/=3;
	variance/=in.width;
	variance/=in.height;
	brightness/=3;
	brightness/=in.width;
	brightness/=in.height;
	for(int c=0;c<3;c++)
	{
		average[c]/=in.width;
		average[c]/=in.height;
	}
	lum=(average[0]+average[1]+average[2])/3;
	variance-=lum*lum;
	printf("     Range: [%.5f,%.5f], [%.5f,%.5f], [%.5f,%.5f]\n",min[0],max[0],min[1],max[1],min[2],max[2]);
	printf("   Average: %.5f %.5f %.5f\n",average[0],average[1],average[2]);
	printf(" Luminance: %.5f\t%.5f\n",lum,sqrt(variance));
	printf("Brightness: %.5f\n",exp(brightness));
}
int main( int argc , char* argv[] )
{
	double t = Time();
	MultiStreamIOServer ioServer;
	_putenv( "TMP=" );
	int paramNum=sizeof(params)/sizeof(cmdLineReadable*);
	cmdLineParse(argc-1,&argv[1],paramNum,params,0);

	bool remapColors = DCTermOffset.set || Gamma.set || Brightness.set || Down.set || GammaCorrection.set;

	if( !In.set )
	{
		ShowUsage(argv[0]);
		return EXIT_FAILURE;
	}
	if( Down.set && Clip.set )
	{
		fprintf(stderr,"Only one operation can be implemented at a time\n");
		return EXIT_FAILURE;
	}
	MyImage inImage,outImage;

	if( remapColors ) inImage.data = GetReadStream< float        >( In.value , inImage.width , inImage.height , SEPARATE_COLORS , GammaCorrection.set , false , &ioServer );
	else              inImage.data = GetReadStream< unsigned int >( In.value , inImage.width , inImage.height , SEPARATE_COLORS , GammaCorrection.set , false , &ioServer );

	if( !inImage.data )
	{
		fprintf(stderr,"Failed to read in image from: %s\n",In.value);
		return EXIT_FAILURE;
	}

	if( !Out.set )
	{
		if( remapColors ) ImageInfo< float >( inImage );
		else              ImageInfo< unsigned int >( inImage );
		return EXIT_SUCCESS;
	}

	if( !Clip.set && (!Down.set || !Down.value) )
	{
		Clip.set=true;
		Down.set=false;
		Clip.values[0]=Clip.values[1]=0;
		Clip.values[2]=inImage.width;
		Clip.values[3]=inImage.height;
	}
	if( Down.set )
	{
		outImage.width=inImage.width;
		outImage.height=inImage.height;
		for(int i=0;i<Down.value;i++)
		{
			outImage.width=(outImage.width+1)/2;
			outImage.height=(outImage.height+1)/2;
		}
	}
	else
	{
		int startX,startY,cWidth,cHeight;
		startX=Clip.values[0];
		startY=Clip.values[1];
		cWidth=Clip.values[2];
		cHeight=Clip.values[3];
		if(startX+cWidth>inImage.width)		cWidth=inImage.width-startX;
		if(startY+cHeight>inImage.height)	cHeight=inImage.height-startY;
		if(cWidth<=0 || cHeight<=0)
		{
			fprintf(stderr,"Window width and height must be positive: %d x %d\n",cWidth,cHeight);
			return EXIT_FAILURE;
		}
		outImage.width=cWidth;
		outImage.height=cHeight;
	}

	DefaultOutputTileWidth = TileWidth.value;
	DefaultOutputTileHeight = TileHeight.value;
	if( TileExtension.set ) DefaultOutputTileExtension = TileExtension.value;

	if( remapColors ) outImage.data = GetWriteStream< float        >( Out.value , outImage.width , outImage.height , SEPARATE_COLORS , false , Quality.value , &ioServer , HDR.set );
	else              outImage.data = GetWriteStream< unsigned int >( Out.value , outImage.width , outImage.height , SEPARATE_COLORS , false , Quality.value , &ioServer , HDR.set );
	if( !outImage.data )
	{
		fprintf(stderr,"Failed to write image to: %s\n",Out.value);
		return EXIT_FAILURE;
	}

	double dcOffset[3]={ 0 , 0 , 0 };
	if( DCTermOffset.set )
	{
		dcOffset[0] = DCTermOffset.values[0];
		dcOffset[1] = DCTermOffset.values[1];
		dcOffset[2] = DCTermOffset.values[2];
	}
	if(Down.set && Down.value)
	{
		MyImage tempIn,tempOut;
		tempIn.data=inImage.data;
		tempIn.width=inImage.width;
		tempIn.height=inImage.height;
		inImage.data=NULL;
		for(int i=0;i<Down.value-1;i++)
		{
			tempOut.width=(tempIn.width+1)/2;
			tempOut.height=(tempIn.height+1)/2;
			tempOut.data=new MultiStreamIOClient(sizeof(float)*3*tempOut.width,tempOut.height,2,NULL,true);
			printf("Down-Sampling %d x %d -> %d x %d\n",tempIn.width,tempIn.height,tempOut.width,tempOut.height);
//			if(!i)	DownSample4<float,float>(tempIn,tempOut,dcOffset,Gamma.value,Brightness.value);
//			else	DownSample4<float,float>(tempIn,tempOut);
			if(!i)	DownSample2< float >(tempIn,tempOut,dcOffset,Gamma.value,Brightness.value);
			else	DownSample2< float >(tempIn,tempOut);
			delete tempIn.data;
			tempIn.data=tempOut.data;
			tempIn.width=tempOut.width;
			tempIn.height=tempOut.height;
			tempOut.data=NULL;
		}
		printf("Down-Sampling %d x %d -> %d x %d\n",tempIn.width,tempIn.height,outImage.width,outImage.height);
//		if(Down.value==1)	DownSample4<float,float>(tempIn,outImage,dcOffset,Gamma.value,Brightness.value);
//		else				DownSample4<float,float>(tempIn,outImage);
		if(Down.value==1)	DownSample2< float >(tempIn,outImage,dcOffset,Gamma.value,Brightness.value);
		else				DownSample2< float >(tempIn,outImage);
	}
	else if( Clip.set )
	{
		printf( "Clipping: (%d , %d) x (%d , %d)\n",Clip.values[0],Clip.values[1],Clip.values[0]+outImage.width,Clip.values[1]+outImage.height );
		if( remapColors )
			if( !DCTermOffset.set && !Gamma.set && !Brightness.set ) ClipWindow< float >( inImage , outImage , Clip.values[0] , Clip.values[1] );
			else                                                     ClipWindow< float , true >( inImage , outImage , Clip.values[0] , Clip.values[1] , dcOffset , Gamma.value , Brightness.value );
		else
			ClipWindow< unsigned int >( inImage , outImage , Clip.values[0] , Clip.values[1] );
	}
	if( Progress.set ) printf( "Image Process Time: %f\n" , Time() - t );
	return EXIT_SUCCESS;
}