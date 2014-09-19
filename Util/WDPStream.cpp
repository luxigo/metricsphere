#include <wincodec.h>
#include <wincodecsdk.h>
#include <cstringt.h>
#include <atlstr.h>
#include <windows.h>
#include <guiddef.h>
#include <Util/Half/half.h>
#include <Util/BaseMultiStreamIO.h>

struct WDPWriteInfo
{
	IWICStream *piStream;
	IWICImagingFactory *piFactory;
	IWICBitmapEncoder *piEncoder;
	IWICBitmapFrameEncode *piBitmapFrame;
	int width;
};
struct WDPReadInfo
{
    IWICImagingFactory* pFactory;
	IWICBitmapDecoder *pDecoder;
    IWICBitmapFrameDecode* pBitmapFrame;
	IWICFormatConverter *pConverter;
	WICRect rc;
	int width;
};

void* WDPInitWriteColor(char* fileName,int width,int height,int quality , MultiStreamIOServer* ioServer )
{
	WDPWriteInfo *info=(struct WDPWriteInfo*)malloc(sizeof(struct WDPWriteInfo));
	info->width=width;

	IPropertyBag2 *pPropertybag = NULL;

	WICPixelFormatGUID pixelFormat = GUID_WICPixelFormat48bppRGBHalf;
	GUID containerFormat = GUID_ContainerFormatWmp;

	WICPixelFormatGUID pixelFormat2 = pixelFormat;


	if(!SUCCEEDED(CoInitialize(NULL)))															fprintf(stderr,"CoInitialize failed\n"),				exit(0);
	if(!SUCCEEDED(CoCreateInstance(CLSID_WICImagingFactory, NULL, CLSCTX_INPROC_SERVER, IID_IWICImagingFactory, (LPVOID*)&info->piFactory)))	fprintf(stderr,"CoCreateInstance failed\n"),	exit(0);
	if(!SUCCEEDED(info->piFactory->CreateStream(&info->piStream)))								fprintf(stderr,"CreateStream failed\n"),				exit(0);
	if(!SUCCEEDED(info->piStream->InitializeFromFilename(CString(fileName), GENERIC_WRITE)))	fprintf(stderr,"InitializeFromFilename failed\n"),		exit(0);
	if(!SUCCEEDED(info->piFactory->CreateEncoder(containerFormat, NULL, &info->piEncoder)))		fprintf(stderr,"CreateEncoder failed\n"),				exit(0);
	if(!SUCCEEDED(info->piEncoder->Initialize(info->piStream, WICBitmapEncoderNoCache)))		fprintf(stderr,"Initialize failed\n"),					exit(0);
	if(!SUCCEEDED(info->piEncoder->CreateNewFrame(&info->piBitmapFrame, &pPropertybag)))		fprintf(stderr,"CreateNewFrame failed\n"),				exit(0);
	if(!SUCCEEDED(info->piBitmapFrame->Initialize(pPropertybag)))								fprintf(stderr,"Initialize failed\n"),					exit(0);
	if(!SUCCEEDED(info->piBitmapFrame->SetSize(width,height)))									fprintf(stderr,"SetSize failed\n"),						exit(0);
	if(!SUCCEEDED(info->piBitmapFrame->SetPixelFormat(&pixelFormat2)))							fprintf(stderr,"SetPixelFormat failed\n"),				exit(0);
	if(!IsEqualGUID(pixelFormat2, pixelFormat))													fprintf(stderr,"SetPixelFormat failed after all\n"),	exit(0);
	return info;
}
void* WDPInitWriteGray(char* fileName,int width,int height,int quality , MultiStreamIOServer* ioServer )
{
	WDPWriteInfo *info=(struct WDPWriteInfo*)malloc(sizeof(struct WDPWriteInfo));
	info->width=width;

	IPropertyBag2 *pPropertybag = NULL;

	WICPixelFormatGUID pixelFormat = GUID_WICPixelFormat16bppGrayHalf;
	GUID containerFormat = GUID_ContainerFormatWmp;

	WICPixelFormatGUID pixelFormat2 = pixelFormat;


	if(!SUCCEEDED(CoInitialize(NULL)))															fprintf(stderr,"CoInitialize failed\n"),				exit(0);
	if(!SUCCEEDED(CoCreateInstance(CLSID_WICImagingFactory, NULL, CLSCTX_INPROC_SERVER, IID_IWICImagingFactory, (LPVOID*)&info->piFactory)))	fprintf(stderr,"CoCreateInstance failed\n"),	exit(0);
	if(!SUCCEEDED(info->piFactory->CreateStream(&info->piStream)))								fprintf(stderr,"CreateStream failed\n"),				exit(0);
	if(!SUCCEEDED(info->piStream->InitializeFromFilename(CString(fileName), GENERIC_WRITE)))	fprintf(stderr,"InitializeFromFilename failed\n"),		exit(0);
	if(!SUCCEEDED(info->piFactory->CreateEncoder(containerFormat, NULL, &info->piEncoder)))		fprintf(stderr,"CreateEncoder failed\n"),				exit(0);
	if(!SUCCEEDED(info->piEncoder->Initialize(info->piStream, WICBitmapEncoderNoCache)))		fprintf(stderr,"Initialize failed\n"),					exit(0);
	if(!SUCCEEDED(info->piEncoder->CreateNewFrame(&info->piBitmapFrame, &pPropertybag)))		fprintf(stderr,"CreateNewFrame failed\n"),				exit(0);
	if(!SUCCEEDED(info->piBitmapFrame->Initialize(pPropertybag)))								fprintf(stderr,"Initialize failed\n"),					exit(0);
	if(!SUCCEEDED(info->piBitmapFrame->SetSize(width,height)))									fprintf(stderr,"SetSize failed\n"),						exit(0);
	if(!SUCCEEDED(info->piBitmapFrame->SetPixelFormat(&pixelFormat2)))							fprintf(stderr,"SetPixelFormat failed\n"),				exit(0);
	if(!IsEqualGUID(pixelFormat2, pixelFormat))													fprintf(stderr,"SetPixelFormat failed after all\n"),	exit(0);
	return info;
}
void WDPWriteColorRow(void* pixels,void* v,int j)
{
	WDPWriteInfo* info=(struct WDPWriteInfo*)v;
	if(!SUCCEEDED(info->piBitmapFrame->WritePixels(1, info->width*6, info->width*6, (BYTE*)pixels)))	fprintf(stderr,"WritePixels failed\n"),	exit(0);
}
void WDPWriteGrayRow(void* pixels,void* v,int j)
{
	WDPWriteInfo* info=(struct WDPWriteInfo*)v;
	if(!SUCCEEDED(info->piBitmapFrame->WritePixels(1, info->width*2, info->width*2, (BYTE*)pixels)))	fprintf(stderr,"WritePixels failed\n"),	exit(0);
}
void WDPFinalizeWrite(void* v)
{
	WDPWriteInfo* info=(struct WDPWriteInfo*)v;
	if(!SUCCEEDED(info->piBitmapFrame->Commit()))		fprintf(stderr,"Commit1 failed\n"),		exit(0);
	if(!SUCCEEDED(info->piEncoder->Commit()))			fprintf(stderr,"Commit2 failed\n"),		exit(0);
	info->piFactory->Release();
	info->piEncoder->Release();
	info->piBitmapFrame->Release();
	info->piStream->Release();
}

void WDPGetImageSize( char* fileName , int& width , int& height )
{
	WDPReadInfo info;

	WICPixelFormatGUID pixelFormat;
    unsigned uFrameCount=0;
    unsigned int uiWidth=0, uiHeight=0;
	BOOL convertible;

	if(!SUCCEEDED(CoInitialize(NULL)))															fprintf(stderr,"CoInitialize failed\n"),					exit(0);
	if(!SUCCEEDED(CoCreateInstance(CLSID_WICImagingFactory, NULL, CLSCTX_INPROC_SERVER, IID_IWICImagingFactory, (LPVOID*)&info.pFactory)))	fprintf(stderr,"CoCreateInstance failed\n"),	exit(0);
	if(!SUCCEEDED(info.pFactory->CreateDecoderFromFilename(CString(fileName), NULL, GENERIC_READ,WICDecodeMetadataCacheOnLoad, &info.pDecoder)))	fprintf(stderr,"CreadDecoderFromFilename failed\n"),	exit(0);
	if(!SUCCEEDED(info.pDecoder->GetFrameCount(&uFrameCount)))									fprintf(stderr,"GetFrameCount failed\n"),					exit(0);
	if(uFrameCount!=1)																			fprintf(stderr,"Frame count error:[%s] %d != 1\n",fileName,uFrameCount),	exit(0);    
	if(!SUCCEEDED(info.pDecoder->GetFrame(0, &info.pBitmapFrame)))								fprintf(stderr,"GetFrame failed\n"),						exit(0);
    if(!SUCCEEDED(info.pBitmapFrame->GetSize(&uiWidth, &uiHeight)))								fprintf(stderr,"GetSize failed\n"),							exit(0);
	width = uiWidth;
	height = uiHeight;

	info.pDecoder->Release();
	info.pFactory->Release();
	info.pBitmapFrame->Release();

}

void* WDPInitReadColor(char* fileName,int& width,int& height , MultiStreamIOServer* ioServer )
{
	WDPReadInfo *info=(WDPReadInfo*)malloc(sizeof(WDPReadInfo));

	WICPixelFormatGUID pixelFormat;
    unsigned uFrameCount=0;
    unsigned int uiWidth=0, uiHeight=0;
	BOOL convertible;

	if(!SUCCEEDED(CoInitialize(NULL)))															fprintf(stderr,"CoInitialize failed\n"),					exit(0);
	if(!SUCCEEDED(CoCreateInstance(CLSID_WICImagingFactory, NULL, CLSCTX_INPROC_SERVER, IID_IWICImagingFactory, (LPVOID*)&info->pFactory)))	fprintf(stderr,"CoCreateInstance failed\n"),	exit(0);
	if(!SUCCEEDED(info->pFactory->CreateDecoderFromFilename(CString(fileName), NULL, GENERIC_READ,WICDecodeMetadataCacheOnLoad, &info->pDecoder)))	fprintf(stderr,"CreadDecoderFromFilename failed\n"),	exit(0);
	if(!SUCCEEDED(info->pDecoder->GetFrameCount(&uFrameCount)))									fprintf(stderr,"GetFrameCount failed\n"),					exit(0);
	if(uFrameCount!=1)																			fprintf(stderr,"Frame count error: %d != 1\n",uFrameCount),	exit(0);    
    if(!SUCCEEDED(info->pDecoder->GetFrame(0, &info->pBitmapFrame)))							fprintf(stderr,"GetFrame failed\n"),						exit(0);
    if(!SUCCEEDED(info->pBitmapFrame->GetSize(&uiWidth, &uiHeight)))							fprintf(stderr,"GetSize failed\n"),							exit(0);
	if(!SUCCEEDED(info->pBitmapFrame->GetPixelFormat(&pixelFormat)))							fprintf(stderr,"GetPixelFormat failed\n"),					exit(0);
	if(!SUCCEEDED(info->pFactory->CreateFormatConverter(&info->pConverter)))					fprintf(stderr,"CreateFormatConverter failed\n"),			exit(0);
	if(!SUCCEEDED(info->pConverter->CanConvert(pixelFormat,GUID_WICPixelFormat48bppRGBHalf,&convertible)))	fprintf(stderr,"CanConvert failed\n"),			exit(0);
	if(!convertible)																			fprintf(stderr,"Can't convert  source to destination\n"),	exit(0);
	if(!SUCCEEDED(info->pConverter->Initialize(info->pBitmapFrame,GUID_WICPixelFormat48bppRGBHalf,WICBitmapDitherTypeNone,NULL,0.0,WICBitmapPaletteTypeCustom)))	fprintf(stderr,"Initialize failed\n"),	exit(0);

	info->rc.X = 0;
	info->rc.Y = 0;
	info->rc.Width = uiWidth;
	info->rc.Height = 1;

	width=uiWidth;
	height=uiHeight;
	info->width=width;
	return info;
}
void* WDPInitReadGray(char* fileName,int& width,int& height , MultiStreamIOServer* ioServer )
{
	WDPReadInfo *info=(WDPReadInfo*)malloc(sizeof(WDPReadInfo));

	WICPixelFormatGUID pixelFormat;
    unsigned uFrameCount=0;
    unsigned int uiWidth=0, uiHeight=0;
	BOOL convertible;

	if(!SUCCEEDED(CoInitialize(NULL)))															fprintf(stderr,"CoInitialize failed\n"),					exit(0);
	if(!SUCCEEDED(CoCreateInstance(CLSID_WICImagingFactory, NULL, CLSCTX_INPROC_SERVER, IID_IWICImagingFactory, (LPVOID*)&info->pFactory)))	fprintf(stderr,"CoCreateInstance failed\n"),	exit(0);
	if(!SUCCEEDED(info->pFactory->CreateDecoderFromFilename(CString(fileName), NULL, GENERIC_READ,WICDecodeMetadataCacheOnLoad, &info->pDecoder)))	fprintf(stderr,"CreadDecoderFromFilename failed\n"),	exit(0);
	if(!SUCCEEDED(info->pDecoder->GetFrameCount(&uFrameCount)))									fprintf(stderr,"GetFrameCount failed\n"),					exit(0);
	if(uFrameCount!=1)																			fprintf(stderr,"Frame count error: %d != 1\n",uFrameCount),	exit(0);    
    if(!SUCCEEDED(info->pDecoder->GetFrame(0, &info->pBitmapFrame)))							fprintf(stderr,"GetFrame failed\n"),						exit(0);
    if(!SUCCEEDED(info->pBitmapFrame->GetSize(&uiWidth, &uiHeight)))							fprintf(stderr,"GetSize failed\n"),							exit(0);
	if(!SUCCEEDED(info->pBitmapFrame->GetPixelFormat(&pixelFormat)))							fprintf(stderr,"GetPixelFormat failed\n"),					exit(0);
	if(!SUCCEEDED(info->pFactory->CreateFormatConverter(&info->pConverter)))					fprintf(stderr,"CreateFormatConverter failed\n"),			exit(0);
	if(!SUCCEEDED(info->pConverter->CanConvert(pixelFormat,GUID_WICPixelFormat16bppGrayHalf,&convertible)))	fprintf(stderr,"CanConvert failed\n"),			exit(0);
	if(!convertible)																			fprintf(stderr,"Can't convert  source to destination\n"),	exit(0);
	if(!SUCCEEDED(info->pConverter->Initialize(info->pBitmapFrame,GUID_WICPixelFormat16bppGrayHalf,WICBitmapDitherTypeNone,NULL,0.0,WICBitmapPaletteTypeCustom)))	fprintf(stderr,"Initialize failed\n"),	exit(0);

	info->rc.X = 0;
	info->rc.Y = 0;
	info->rc.Width = uiWidth;
	info->rc.Height = 1;

	width=uiWidth;
	height=uiHeight;
	info->width=width;
	return info;
}
void WDPReadColorRow(void* pixels,void* v,int)
{
	WDPReadInfo* info=(WDPReadInfo*)v;
	if(!SUCCEEDED(info->pConverter->CopyPixels(&info->rc,info->width*6,info->width*6,(BYTE*)pixels)))	fprintf(stderr,"CopyPixels failed\n"),	exit(0);
	info->rc.Y++;
}
void WDPReadGrayRow(void* pixels,void* v,int)
{
	WDPReadInfo* info=(WDPReadInfo*)v;
	if(!SUCCEEDED(info->pConverter->CopyPixels(&info->rc,info->width*2,info->width*2,(BYTE*)pixels)))	fprintf(stderr,"CopyPixels failed\n"),	exit(0);
	info->rc.Y++;
}
void WDPFinalizeRead(void* v)
{
	WDPReadInfo* info=(WDPReadInfo*)v;
	info->pConverter->Release();
	info->pDecoder->Release();
	info->pFactory->Release();
	info->pBitmapFrame->Release();
}
