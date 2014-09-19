#include <winsock2.h>
#include <ws2tcpip.h>
#include "MultigridSolver.h"
#include "Util/MultiStreamIO.h"
#include "Util/Socket.h"
#include "Util/Half/half.h"

#define SHOW_AVERAGE 0
#define MISHA_DEBUG 0

#if MISHA_DEBUG
#include "Util/ImageStream.h"
#endif // MISHA_DEBUG

bool IsDownSamplable(int width,int height,int& newWidth,int& newHeight)
{
	if( (width&3) || (height&3) || (width<20) )	return false;
	newWidth  = width>>1;
	newHeight = height>>1;
	return true;
}
bool IsDownSamplable(int dCount,ClientSocket* cSockets,int sockCount)
{
	for(int i=0;i<sockCount;i++)	if ( (cSockets[i].start>>dCount)&3 || ((cSockets[i].end-cSockets[i].start)>>dCount)<20)	return false;
	return true;
}
template<class Real>
void SetRectangularLaplacianMatrix( SparseMatrix< Real > & lap , int width , int height , double iWeight=0 )
{
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::FullDotProductStencil dotStencil1,d2DotStencil1,dotStencil2,d2DotStencil2;
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(width,dotStencil1,0,0);
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(width,d2DotStencil1,1,1,false);
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(height,dotStencil2,0,0);
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(height,d2DotStencil2,1,1,false);
	SetRectangularLaplacianMatrix( dotStencil1 , d2DotStencil1 , dotStencil2 , d2DotStencil2 , lap , width , height , iWeight );
}
template<class Real>
void SetRectangularLaplacianMatrix( DotProductStencil& dotMajor , DotProductStencil& d2DotMajor , DotProductStencil& dotMinor , DotProductStencil& d2DotMinor ,
								    SparseMatrix< Real > & lap , int width , int height , double iWeight=0 )
{
	lap.Resize( width*height );

	int ii , jj;
	int xStart , xEnd , yStart , yEnd;
	for(int y=0;y<height;y++)
	{
		if		(y<Degree)			jj = y							, yStart = -y		, yEnd = Degree;
		else if	(y>height-1-Degree)	jj = 2*Degree+(y-(height-1))	, yStart = -Degree	, yEnd = height-1-y;
		else						jj = Degree						, yStart = -Degree	, yEnd = Degree;

		for(int x=0;x<width;x++)
		{
			if		(x<Degree)			ii = x						, xStart = -x		, xEnd = Degree;
			else if	(x>width-1-Degree)	ii = 2*Degree+(x-(width-1))	, xStart = -Degree	, xEnd = width-1-x;
			else						ii = Degree					, xStart = -Degree	, xEnd = Degree;

			int idx = x+y*width;
			lap.SetGroupSize( idx , (yEnd-yStart+1)*(xEnd-xStart+1) );
			int _i = 0;
			for(int yy=yStart;yy<=yEnd;yy++)
				for(int xx=xStart;xx<=xEnd;xx++)
				{
					lap.m_ppElements[idx][_i  ].N=(x+xx)+(y+yy)*width;
					lap.m_ppElements[idx][_i++].Value=
						dotMajor.caseTable[ii].values[xx+Degree]*d2DotMinor.caseTable[jj].values[yy+Degree]+
						d2DotMajor.caseTable[ii].values[xx+Degree]*dotMinor.caseTable[jj].values[yy+Degree]+
						(dotMajor.caseTable[ii].values[xx+Degree]*dotMinor.caseTable[jj].values[yy+Degree])*iWeight;
				}
		}
	}
}

template<class Real>
void SetSphericalLaplacianMatrix(SparseMatrix<Real>& lap,int width,int height)
{
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::FullDotProductStencil dotStencil1,d2DotStencil1,dotStencil2,d2DotStencil2;
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(width,dotStencil1,0,0);
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(width,d2DotStencil1,1,1,false);
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(height,dotStencil2,0,0);
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(height,d2DotStencil2,1,1,false);
	SetSphericalLaplacianMatrix( dotStencil1 , d2DotStencil1 , dotStencil2 , d2DotStencil2 , lap , width , height );
}
template<class Real>
void SetSphericalLaplacianMatrix(DotProductStencil& dotMajor,DotProductStencil& d2DotMajor,DotProductStencil& dotMinor,DotProductStencil& d2DotMinor,
								 SparseMatrix<Real>& lap,int width,int height)
{
	lap.Resize(width*height);

	for(int x=0;x<width;x++)
		for(int y=0;y<height;y++)
		{
			int idx1 = x + y*width;
			lap.SetGroupSize( idx1 , (2*Degree+1)*(2*Degree+1) );
			for(int i=-Degree;i<=Degree;i++)
				for(int j=-Degree;j<=Degree;j++)
				{
					int idx2 = (i+Degree)+(j+Degree)*(2*Degree+1);
					int xx=(x+i+width)%width;
					int yy= y+j;
					if(yy<0)
					{
						yy = -yy-1;
						if(xx<width/2)	xx = width/2-1-xx;
						else			xx = width/2+width-1-xx;
					}
					else if(yy>=height)
					{
						yy = height-1-(yy-height);
						xx = width-1-xx;
					}
					lap.m_ppElements[idx1][idx2].N = xx+yy*width;
					lap.m_ppElements[idx1][idx2].Value = 
						dotMajor.caseTable[Degree].values[i+Degree]*d2DotMinor.caseTable[Degree].values[j+Degree]+
						d2DotMajor.caseTable[Degree].values[i+Degree]*dotMinor.caseTable[Degree].values[j+Degree];
				}
		}
}
template<class Real>
void SetSphericalLaplacianMatrix(SparseMatrix<Real>& lap,int width,int height,double iWeight)
{
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::FullDotProductStencil dotStencil1,d2DotStencil1,dotStencil2,d2DotStencil2;
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(width,dotStencil1,0,0);
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(width,d2DotStencil1,1,1,false);
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(height,dotStencil2,0,0);
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(height,d2DotStencil2,1,1,false);
	SetSphericalLaplacianMatrix( dotStencil1 , d2DotStencil1 , dotStencil2 , d2DotStencil2 , lap , width , height , iWeight );
}
template<class Real>
void SetSphericalLaplacianMatrix(DotProductStencil& dotMajor,DotProductStencil& d2DotMajor,DotProductStencil& dotMinor,DotProductStencil& d2DotMinor,
								 SparseMatrix<Real>& lap,int width,int height,double iWeight)
{

	lap.Resize(width*height);

	for(int x=0;x<width;x++)
		for(int y=0;y<height;y++)
		{
			int idx1 = x + y*width;
			lap.SetGroupSize( idx1 , (2*Degree+1)*(2*Degree+1) );
			for(int i=-Degree;i<=Degree;i++)
				for(int j=-Degree;j<=Degree;j++)
				{
					int idx2 = (i+Degree)+(j+Degree)*(2*Degree+1);
					int xx=(x+i+width)%width;
					int yy= y+j;
					if(yy<0)
					{
						yy = -yy-1;
						if(xx<width/2)	xx = width/2-1-xx;
						else			xx = width/2+width-1-xx;
					}
					else if(yy>=height)
					{
						yy = height-1-(yy-height);
						xx = width-1-xx;
					}
					lap.m_ppElements[idx1][idx2].N = xx+yy*width;
					lap.m_ppElements[idx1][idx2].Value = 
						dotMajor.caseTable[Degree].values[i+Degree]*d2DotMinor.caseTable[Degree].values[j+Degree]+
						d2DotMajor.caseTable[Degree].values[i+Degree]*dotMinor.caseTable[Degree].values[j+Degree]+
						(dotMajor.caseTable[Degree].values[i+Degree]*dotMinor.caseTable[Degree].values[j+Degree])*iWeight;
				}
		}
}

DWORD WINAPI AcceptThread(LPVOID lpParam)
{
	AcceptThreadData* params = (AcceptThreadData*)lpParam;
	*params->connectSocket=accept(*params->listenSocket,NULL,NULL);
	if (*params->connectSocket == INVALID_SOCKET)
	{
		fprintf( stderr , "accept failed: %s\n" , LastSocketError() );
		params->success = false;
	}
	else
	{
		int val=1;
		setsockopt(*params->connectSocket,IPPROTO_TCP,TCP_NODELAY,(char*)&val,sizeof(val));
		closesocket(*params->listenSocket);
		*params->listenSocket = INVALID_SOCKET;
		params->success = true;
	}
	return 0;
}
///////////////////////////
// SphericalSynchronizer //
///////////////////////////
template<int Channels>
SphericalSynchronizer<Channels>::SphericalSynchronizer(void)
{
	_cCount = 0;
	_clientSockets = NULL;
}
template<int Channels>
void SphericalSynchronizer<Channels>::init(int width,const int* widths,int cCount,int rPasses,int pPasses,int vCycles)
{
	_width=width;
	_cCount=cCount;
	_rPasses=rPasses;
	_pPasses=pPasses;
	_vCycles=vCycles;

	if(_clientSockets)	delete[] _clientSockets;

	_clientSockets = new ClientSocket[cCount];
	if(!_clientSockets)
	{
		fprintf(stderr,"Failed to allocate SphericalSynchronizer::_clientSockets\n");
		exit(0);
	}
	for(int i=0;i<cCount;i++)
	{
		_clientSockets[i].client = GetListenSocket(_clientSockets[i].port);
		if (_clientSockets[i].client == INVALID_SOCKET)	fprintf(stderr,"Failed to set server in SphericalSynchronizer\n")	,	exit(0);
		if(i)	_clientSockets[i].start = _clientSockets[i-1].end;
		else	_clientSockets[i].start = 0;
		_clientSockets[i].end = _clientSockets[i].start+widths[i];
	}
}
template<int Channels>
SphericalSynchronizer<Channels>::~SphericalSynchronizer(void)
{
	delete[] _clientSockets;
	_clientSockets=NULL;
}
template<int Channels>	int SphericalSynchronizer<Channels>::clients(void)		const	{ return _cCount; }
template<int Channels>	int SphericalSynchronizer<Channels>::port(int cIndex)	const	{ return _clientSockets[cIndex].port; }
template<int Channels>
template<class DataType>
void SphericalSynchronizer<Channels>::Run(void)
{
	for(int i=0;i<_cCount;i++)
	{
		SOCKET temp = accept(_clientSockets[i].client,NULL,NULL);
		if (temp == INVALID_SOCKET)	fprintf( stderr , "accept failed in SphericalSynchronizer::Run(): %s\n" , LastSocketError() )	,	exit(0);
		int val=1;
		setsockopt(temp,IPPROTO_TCP,TCP_NODELAY,(char*)&val,sizeof(val));
		closesocket(_clientSockets[i].client);
		_clientSockets[i].client=temp;
	}

	DataType *buffer1,*buffer2;
	buffer1 = new DataType[_width*Channels*Degree];
	buffer2 = new DataType[_width*Channels*Degree];
	for ( int v = 0 ; v < _vCycles ; v++ )
	{
		//////////////////////////////////////////
		// First synchronize on the restriction //
		//////////////////////////////////////////

		// Synchronize the head of the stream
		for(int k=0;k<_rPasses;k++)
		{
			for(int i=0;i<_cCount;i++)
			{
				int w = _clientSockets[i].end-_clientSockets[i].start;
				if(!ReceiveOnSocket( _clientSockets[i].client , buffer1 , sizeof(DataType)*w*Channels*Degree ))
					fprintf(stderr,"Failed to receive on head restriction (%d %d) %d / %d\n",_rPasses,_pPasses,i,_cCount);
				for(int d=0;d<Degree;d++)
					for( int x = 0 ; x < w ; x++ )
						for( int c = 0; c < Channels ; c++ )
							buffer2[(x+_clientSockets[i].start)*Channels+d*_width*Channels+c] = buffer1[d*w*Channels+w*c+x];
			}
			for(int d=0;d<Degree;d++)
			{
				int dd=Degree-1-d;
				for(int x=0;x<_width/2;x++)
				{
					int xx = _width/2-1-x;
					for(int c=0;c<Channels;c++)	buffer1[xx*Channels+c+dd*_width*Channels] = buffer2[x*Channels+c+d*_width*Channels];
				}
				for(int x=_width/2;x<_width;x++)
				{
					int xx = _width/2+_width-1-x;
					for(int c=0;c<Channels;c++)	buffer1[xx*Channels+c+dd*_width*Channels] = buffer2[x*Channels+c+d*_width*Channels];
				}
			}
			for(int i=0;i<_cCount;i++)
			{
				int w = _clientSockets[i].end-_clientSockets[i].start;
				for(int d=0;d<Degree;d++)
					for( int x = 0 ; x < w ; x++ )
						for( int c = 0 ; c < Channels ; c++ )
							 buffer2[d*w*Channels+w*c+x] = buffer1[(x+_clientSockets[i].start)*Channels+d*_width*Channels+c];
				if(!SendOnSocket( _clientSockets[i].client , buffer2 , sizeof(DataType)*w*Channels*Degree ))
					fprintf(stderr,"Failed to send on head restriction (%d %d) %d / %d\n",_rPasses,_pPasses,i,_cCount);
			}
		}
		// Synchronize the tail of the stream
		for(int k=0;k<_rPasses;k++)
		{
			for(int i=0;i<_cCount;i++)
			{
				int w = _clientSockets[i].end-_clientSockets[i].start;
				if(!ReceiveOnSocket( _clientSockets[i].client , buffer1 , sizeof(DataType)*w*Channels*Degree ))
					fprintf(stderr,"Failed to receive on tail restriction (%d %d) %d / %d\n",_rPasses,_pPasses,i,_cCount);
				for(int d=0;d<Degree;d++)
					for( int x = 0 ; x < w ; x++ )
						for( int c = 0; c < Channels ; c++ )
							buffer2[(x+_clientSockets[i].start)*Channels+d*_width*Channels+c] = buffer1[d*w*Channels+w*c+x];
			}
			for(int d=0;d<Degree;d++)
			{
				int dd=Degree-1-d;
				for(int x=0;x<_width;x++)
				{
					int xx = _width-1-x;
					for(int c=0;c<Channels;c++)	buffer1[xx*Channels+c+dd*_width*Channels] = buffer2[x*Channels+c+d*_width*Channels];
				}
			}
			for(int i=0;i<_cCount;i++)
			{
				int w = _clientSockets[i].end-_clientSockets[i].start;
				for(int d=0;d<Degree;d++)
					for( int x = 0 ; x < w ; x++ )
						for( int c = 0 ; c < Channels ; c++ )
							 buffer2[d*w*Channels+w*c+x] = buffer1[(x+_clientSockets[i].start)*Channels+d*_width*Channels+c];
				if(!SendOnSocket( _clientSockets[i].client , buffer2 , sizeof(DataType)*w*Channels*Degree ))
					fprintf(stderr,"Failed to send on tail restriction (%d %d) %d / %d\n",_rPasses,_pPasses,i,_cCount);
			}
		}

		/////////////////////////////////////////
		// Now synchronize on the prolongation //
		/////////////////////////////////////////

		// Synchronize the head of the stream
		for(int k=0;k<_pPasses;k++)
		{
			for(int i=0;i<_cCount;i++)
			{
				int w = _clientSockets[i].end-_clientSockets[i].start;
				if(!ReceiveOnSocket( _clientSockets[i].client , buffer1 , sizeof(DataType)*w*Channels*Degree ))
					fprintf(stderr,"Failed to receive on head prolongation (%d %d) %d / %d\n",_rPasses,_pPasses,i,_cCount);
				for(int d=0;d<Degree;d++)
					for( int x = 0 ; x < w ; x++ )
						for( int c = 0; c < Channels ; c++ )
							buffer2[(x+_clientSockets[i].start)*Channels+d*_width*Channels+c] = buffer1[d*w*Channels+w*c+x];
			}
			for(int d=0;d<Degree;d++)
			{
				int dd=Degree-1-d;
				for(int x=0;x<_width/2;x++)
				{
					int xx = _width/2-1-x;
					for(int c=0;c<Channels;c++)	buffer1[xx*Channels+c+dd*_width*Channels] = buffer2[x*Channels+c+d*_width*Channels];
				}
				for(int x=_width/2;x<_width;x++)
				{
					int xx = _width/2+_width-1-x;
					for(int c=0;c<Channels;c++)	buffer1[xx*Channels+c+dd*_width*Channels] = buffer2[x*Channels+c+d*_width*Channels];
				}
			}
			for(int i=0;i<_cCount;i++)
			{
				int w = _clientSockets[i].end-_clientSockets[i].start;
				for(int d=0;d<Degree;d++)
					for( int x = 0 ; x < w ; x++ )
						for( int c = 0 ; c < Channels ; c++ )
							 buffer2[d*w*Channels+w*c+x] = buffer1[(x+_clientSockets[i].start)*Channels+d*_width*Channels+c];
				if(!SendOnSocket( _clientSockets[i].client , buffer2 , sizeof(DataType)*w*Channels*Degree ))
					fprintf(stderr,"Failed to send on head prolongation (%d %d) %d / %d\n",_rPasses,_pPasses,i,_cCount);
			}
		}
		// Synchronize the tail of the stream
		for(int k=0;k<_pPasses;k++)
		{
			for(int i=0;i<_cCount;i++)
			{
				int w = _clientSockets[i].end-_clientSockets[i].start;
				if(!ReceiveOnSocket( _clientSockets[i].client , buffer1 , sizeof(DataType)*w*Channels*Degree ))
					fprintf(stderr,"Failed to receive on tail prolongation (%d %d) %d / %d\n",_rPasses,_pPasses,i,_cCount);
				for(int d=0;d<Degree;d++)
					for( int x = 0 ; x < w ; x++ )
						for( int c = 0; c < Channels ; c++ )
							buffer2[(x+_clientSockets[i].start)*Channels+d*_width*Channels+c] = buffer1[d*w*Channels+w*c+x];
			}
			for(int d=0;d<Degree;d++)
			{
				int dd=Degree-1-d;
				for(int x=0;x<_width;x++)
				{
					int xx = _width-1-x;
					for(int c=0;c<Channels;c++)	buffer1[xx*Channels+c+dd*_width*Channels] = buffer2[x*Channels+c+d*_width*Channels];
				}
			}
			for(int i=0;i<_cCount;i++)
			{
				int w = _clientSockets[i].end-_clientSockets[i].start;
				for(int d=0;d<Degree;d++)
					for( int x = 0 ; x < w ; x++ )
						for( int c = 0 ; c < Channels ; c++ )
							 buffer2[d*w*Channels+w*c+x] = buffer1[(x+_clientSockets[i].start)*Channels+d*_width*Channels+c];
				if(!SendOnSocket( _clientSockets[i].client , buffer2 , sizeof(DataType)*w*Channels*Degree ))
					fprintf(stderr,"Failed to send on tail prolongation (%d %d) %d / %d\n",_rPasses,_pPasses,i,_cCount);
			}
		}
	}
	delete[] buffer1;
	delete[] buffer2;
}
template<int Channels>
template<class DataType>
DWORD WINAPI SphericalSynchronizer<Channels>::RunThread(LPVOID lpParam)
{
	SphericalSynchronizer<Channels>* synchronizer = (SphericalSynchronizer*)lpParam;
	synchronizer->Run<DataType>();
	return 0;
}

/////////////////////////////
// SocketedMultigridServer //
/////////////////////////////
#if !USE_MY_WIN_SOCK
template<int Channels>
int SocketedMultigridServer<Channels>::_wsaCount = 0;

template<int Channels>
WSADATA SocketedMultigridServer<Channels>::_wsaData;
#endif // !USE_MY_WIN_SOCK
template<int Channels>
SocketedMultigridServer<Channels>::SocketedMultigridServer(void)
{
	// Initialize Winsock
#if USE_MY_WIN_SOCK
	MyWinSock::Load();
#else // !USE_MY_WIN_SOCK
	if(!_wsaCount)
		if(WSAStartup(MAKEWORD(2,2), &_wsaData))	fprintf(stderr,"WSAStartup failed: %s\n", LastSocketError())	,	exit(0);
	_wsaCount++;
#endif // USE_MY_WIN_SOCK
	clientSockets=NULL;
	_icSynchronizers = _oocSynchronizers = NULL;
	_synchronizerHandles = NULL;
	_syncSockets = NULL;
}
template<int Channels>
SocketedMultigridServer<Channels>::~SocketedMultigridServer(void)
{
	if(clientSockets)	delete[] clientSockets;
	clientSockets=NULL;
#if USE_MY_WIN_SOCK
	MyWinSock::UnLoad();
#else // !USE_MY_WIN_SOCK
	_wsaCount--;
	if(!_wsaCount)	WSACleanup();
#endif // USE_MY_WIN_SOCK
	if(_icSynchronizers)	delete[] _icSynchronizers;
	if(_oocSynchronizers)	delete[] _oocSynchronizers;
	_icSynchronizers = _oocSynchronizers = NULL;
	if(_spherical)
	{
#if FROM_IMAGE
//		WaitForMultipleObjects(_oocDCount+1+_icDCount,_synchronizerHandles,true,INFINITE);
		MyWaitForMultipleObjects( _oocDCount+1+_icDCount , _synchronizerHandles , 10000 , "SocketedMultigridServer<Channels>::~SocketedMultigridServer");
		for(int i=0;i<_oocDCount+1+_icDCount;i++)	CloseHandle(_synchronizerHandles[i]);
#else // !FROM_IMAGE
//		WaitForMultipleObjects( _oocDCount+_icDCount , _synchronizerHandles , true,INFINITE);
		MyWaitForMultipleObjects( _oocDCount+_icDCount , _synchronizerHandles , 10000 , "SocketedMultigridServer<Channels>::~SocketedMultigridServer");
		for(int i=0;i<_oocDCount+_icDCount;i++)	CloseHandle(_synchronizerHandles[i]);
#endif // FROM_IMAGE
		delete[] _synchronizerHandles;
		_synchronizerHandles = NULL;
		if(_syncSockets)
		{
			for(int i=0;i<_icDCount;i++)	if(_syncSockets[i] != INVALID_SOCKET)	closesocket(_syncSockets[i]);
			delete[] _syncSockets;
		}
		_syncSockets = NULL;
	}
}
template<int Channels>
#if FROM_IMAGE
#if SUPPORT_SHARPEN
bool SocketedMultigridServer<Channels>::SetUp(int port,int clientCount,int iters,int inCoreRes,int minMGRes,int vCycles,
											  int quality,int lanes,bool verbose,bool spherical,bool sharpen,double iWeight,double gWeight,bool showProgress,bool noCG)
#else // !SUPPORT_SHARPEN
bool SocketedMultigridServer<Channels>::SetUp(int port,int clientCount,int iters,int inCoreRes,int minMGRes,int vCycles,
											  int quality,int lanes,bool verbose,bool spherical,bool showProgress,bool noCG)
#endif // SUPPORT_SHARPEN
#else // !FROM_IMAGE
bool SocketedMultigridServer<Channels>::SetUp(int port,int clientCount,int width,int height,int iters,int inCoreRes,int minMGRes,double average[Channels],
											  int vCycles,int cWidth,int cHeight,
											  int quality,int lanes,bool verbose,bool spherical,bool showProgress,bool noCG)
#endif // FROM_IMAGE
{
	SOCKET listenSocket = INVALID_SOCKET;
	clientSockets = new ClientSocket[clientCount];

	_spherical	= spherical;
	_verbose	= verbose;
	_noCG		= noCG;
	_cCount		= clientCount;
	_iters		= iters;
	_vCycles	= vCycles;
	_inCoreRes	= inCoreRes;
	_minMGRes	= minMGRes;
#if SUPPORT_SHARPEN
	_sharpen	= sharpen;
	_iWeight	= iWeight;
	_gWeight	= gWeight;
#endif // SUPPORT_SHARPEN
#if !FROM_IMAGE
	for(int c=0;c<Channels;c++)	_average[c]=average[c];
#endif // !FROM_IMAGE

	// Create a SOCKET for connecting to server
	listenSocket = GetListenSocket(port);
	if (listenSocket == INVALID_SOCKET)	return false;

	char hostAddress[512];
	GetHostAddress(hostAddress);
	printf ( "Server Address: %s:%d\n", hostAddress , port );

#if SOCKET_BACKED_GRID
	for(int i=0;i<clientCount*3;i++)
	{
		SOCKET sock = AcceptSocket( listenSocket );
		ClientSocketInfo info;

		ReceiveOnSocket( sock , &info , sizeof( info ) );
		switch( info.type )
		{
		case ClientSocketInfo::DATA:
			clientSockets[info.index].client  = sock;
			break;
		case ClientSocketInfo::X:
			clientSockets[info.index].clientX = sock;
			break;
		case ClientSocketInfo::B:
			clientSockets[info.index].clientB = sock;
			break;
		default:
			fprintf( stderr , "Unrecognized client socket type: %d\n" , info.type );
			return false;
		}
	}
#else // !SOCKET_BACKED_GRID
	for(int i=0;i<clientCount;i++)
	{
		// Accept a client socket
		clientSockets[i].client = accept(listenSocket, NULL, NULL);
		if (clientSockets[i].client == INVALID_SOCKET)
		{
			fprintf( stderr , "accept[%d] failed: %s\n" , i , LastSocketError() );
			return false;
		}
		int val=1;
		setsockopt( clientSockets[i].client , IPPROTO_TCP , TCP_NODELAY , (char*)&val , sizeof(val) );
	}
#endif // SOCKET_BACKED_GRID
    // No longer need the server socket
    closesocket(listenSocket);
	listenSocket = INVALID_SOCKET;

	// Get the client info
	ClientData cd;
#if FROM_IMAGE
	int width = 0 , height;
#endif // FROM_IMAGE
	for( int i=0 ; i<clientCount ; i++ )
	{
		ReceiveOnSocket( clientSockets[i].client,&cd,sizeof(cd) ) ;
#if FROM_IMAGE
		clientSockets[i].index=cd.index;
		if(!i) height = cd.height;
		else if( height!=cd.height)	fprintf(stderr,"Stream heights differ: %d != %d\n",height,cd.height) , exit(0);
		width += cd.width;
		clientSockets[i].start=0;
		clientSockets[i].end=cd.width;
#else // !FROM_IMAGE
		clientSockets[i].start=cd.start;
#endif // FROM_IMAGE
	}
#if FROM_IMAGE
	_cWidth = width;
	_cHeight = height;
	if(!spherical)	_width=width+1,	_height=height+1;
	else			_width=width  , _height=height;
	long long paddedWidth=minMGRes;
	long long paddedHeight=minMGRes;
	int domainW=FiniteElements1D<float,Type,Degree>::DomainSize(_width);
	int domainH=FiniteElements1D<float,Type,Degree>::DomainSize(_height);
	while(paddedWidth<domainW)	paddedWidth*=2;
	while(paddedHeight<domainH)	paddedHeight*=2;
	int blockSize;
	blockSize=paddedWidth/minMGRes;
	blockSize*=2;
	while(paddedWidth-blockSize>domainW)	paddedWidth-=blockSize;

	blockSize=paddedHeight/minMGRes;
	blockSize*=2;
	while(paddedHeight-blockSize>domainH)	paddedHeight-=blockSize;

	paddedWidth=FiniteElements1D<float,Type,Degree>::Dimension(paddedWidth);
	paddedHeight=FiniteElements1D<float,Type,Degree>::Dimension(paddedHeight);

	_width = paddedWidth;
	_height = paddedHeight;
#else // !FROM_IMAGE
	_width=width;
	_height=height;
	_cWidth=cWidth;
	_cHeight=cHeight;
#endif // FROM_IMAGE

	// Sort the clients by start order
	qsort(clientSockets,clientCount,sizeof(ClientSocket),ClientSocket::Sort);

	if(clientSockets[clientCount-1].start>=width)
	{
		fprintf(stderr,"Starting position exceeds image width: %d >= %d\n",clientSockets[clientCount-1].start,width);
		return false;
	}
	// Send the image dimensions to the client
#if FROM_IMAGE
	int start=0;
#endif // FROM_IMAGE
	for(int i=0;i<clientCount;i++)
	{
		ServerData sd;
		sd.width=_width;
		sd.height=_height;
#if FROM_IMAGE
		clientSockets[i].start+=start;
		clientSockets[i].end+=start;
		if(i==clientCount-1) clientSockets[i].end=_width;
		sd.start = clientSockets[i].start;
		sd.end = clientSockets[i].end;
		start = clientSockets[i].end;
//		if( i==clientCount-1 )	sd.end = _width;
//		else					sd.end = clientSockets[i].end;
#else // !FROM_IMAGE
		if( i==clientCount-1 )	sd.end = _width;
		else					sd.end = clientSockets[i+1].start;
		sd.cHeight = _cHeight;
		clientSockets[i].end = sd.end;
#endif // FROM_IMAGE
		if(_cWidth>sd.end)						sd.cEnd=sd.end;
		else if(_cWidth>clientSockets[i].start)	sd.cEnd=_cWidth;
		else									sd.cEnd=clientSockets[i].start;
		sd.iters=_iters;
		sd.vCycles=_vCycles;
		sd.quality=quality;
		sd.lanes=lanes;
		sd.verbose=verbose;
		sd.spherical=spherical;
		sd.progress=showProgress;
#if SUPPORT_SHARPEN
		sd.sharpen=sharpen;
		sd.iWeight=iWeight;
		sd.gWeight=gWeight;
#endif // SUPPORT_SHARPEN
		SendOnSocket ( clientSockets[i].client,&sd,sizeof(sd) );
#if DEBUG_FLAGS
		SendOnSocket ( clientSockets[i].client,&debugFlags,sizeof(debugFlags) );
#endif // DEBUG_FLAGS
	}
	if(spherical)
	{
		// Get the addresses and ports of the clients and pass them on to the neighbors
		for(int i=0;i<clientCount;i++)
		{
			ReceiveOnSocket ( clientSockets[i].client,&clientSockets[i].address,sizeof(clientSockets[i].address) );
			ReceiveOnSocket ( clientSockets[i].client,&clientSockets[i].port,sizeof(clientSockets[i].port) );
		}
		for(int i=0;i<clientCount;i++)
		{
			SendOnSocket ( clientSockets[(i-1+clientCount)%clientCount].client,&clientSockets[i].address,sizeof(clientSockets[i].address) );
			SendOnSocket ( clientSockets[(i-1+clientCount)%clientCount].client,&clientSockets[i].port,sizeof(clientSockets[i].port) );
		}
	}
	else
	{
		// Get the addresses and ports of the clients and pass them on to the neighbors
		for(int i=1;i<clientCount;i++)
		{
			ReceiveOnSocket ( clientSockets[i].client,&clientSockets[i].address,sizeof(clientSockets[i].address) );
			ReceiveOnSocket ( clientSockets[i].client,&clientSockets[i].port,sizeof(clientSockets[i].port) );
		}
		for(int i=1;i<clientCount;i++)
		{
			SendOnSocket ( clientSockets[i-1].client,&clientSockets[i].address,sizeof(clientSockets[i].address) );
			SendOnSocket ( clientSockets[i-1].client,&clientSockets[i].port,sizeof(clientSockets[i].port) );
		}
	}
	// Now everybody is talking to each other...

	// Figure out how many different out-of-core resolutions are needed
	_oocDCount=0;
	while(IsDownSamplable(_oocDCount,clientSockets,clientCount) && long long(width>>_oocDCount)*long long(height>>_oocDCount)>=inCoreRes*inCoreRes)	_oocDCount++;
	for(int i=0;i<clientCount;i++)	SendOnSocket ( clientSockets[i].client,&_oocDCount,sizeof(_oocDCount) );

	// Figure out how many different in-core resolutions are needed
	_icDCount = 1;
	int w = _width>>(_oocDCount-1) , h = _height>>(_oocDCount-1);
	int ww = w , hh = h , ww2 , hh2;
	while(IsDownSamplable(ww,hh,ww2,hh2) && ww2*hh2>=_minMGRes*_minMGRes && ww2>=16 && hh2>=16)
//	while(IsDownSamplable(ww,hh,ww2,hh2) && ww2*hh2>=_minMGRes*_minMGRes)
	{
		_icDCount++;
		ww=ww2;
		hh=hh2;
	}

	if(_spherical)
	{
#if FROM_IMAGE
		_oocSynchronizers = new SphericalSynchronizer<Channels>[2*_oocDCount+1];
#else // !FROM_IMAGE
		_oocSynchronizers = new SphericalSynchronizer<Channels>[2*_oocDCount];
#endif // FROM_IMAGE
		_icSynchronizers = new SphericalSynchronizer<Channels>[2*_icDCount];
		if(!_oocSynchronizers || !_icSynchronizers)
		{
			fprintf(stderr,"Failed to allocate synchronizers\n");
			return false;
		}

		int *oocWidths,icWidths;
		oocWidths = new int[_cCount];

#if FROM_IMAGE
		for(int i=0;i<_cCount;i++)	oocWidths[i] = (clientSockets[i].end-clientSockets[i].start);
		_oocSynchronizers[0].init(_width,oocWidths,_cCount,1,0,1);
		for(int d=0;d<_oocDCount;d++)
		{
			for(int i=0;i<_cCount;i++)	oocWidths[i] = (clientSockets[i].end-clientSockets[i].start)>>d;
			_oocSynchronizers[2*d+1].init(_width>>d,oocWidths,_cCount,_iters*Degree+1,_iters*Degree+1,_vCycles);
			if(_verbose)	_oocSynchronizers[2*d+2].init(_width>>d,oocWidths,_cCount,Degree,Degree,_vCycles);
			else			_oocSynchronizers[2*d+2].init(_width>>d,oocWidths,_cCount,Degree,     0,_vCycles);
		}
#else // !FROM_IMAGE
		for(int d=0;d<_oocDCount;d++)
		{
			for(int i=0;i<_cCount;i++)	oocWidths[i] = (clientSockets[i].end-clientSockets[i].start)>>d;
			_oocSynchronizers[d].init(_width>>d,oocWidths,_cCount,_iters*Degree+1,_iters*Degree+1,_vCycles);
		}
#endif // FROM_IMAGE
		delete[] oocWidths;

		for(int d=0;d<_icDCount;d++)
		{
			icWidths = _width>>(d+_oocDCount-1);
			_icSynchronizers[2*d].init(icWidths,&icWidths,1,_iters*Degree+1,_iters*Degree+1,_vCycles);
			if(_verbose)	_icSynchronizers[2*d+1].init(icWidths,&icWidths,1,Degree,Degree,_vCycles);
			else			_icSynchronizers[2*d+1].init(icWidths,&icWidths,1,Degree,     0,_vCycles);
		}
#if FROM_IMAGE
		_synchronizerHandles = new HANDLE[2*(_oocDCount+_icDCount)+1];
		{
			DWORD synchThreadID;
			if( sharpen )
			{
				_synchronizerHandles[0]=CreateThread(
					NULL,												// default security attributes
					0,													// use default stack size  
					SphericalSynchronizer<Channels>::RunThread<float>,	// thread function 
					&_oocSynchronizers[0],								// argument to thread function 
					0,													// use default creation flags 
					&synchThreadID);									// returns the thread identifier
			}
			else
			{
				_synchronizerHandles[0]=CreateThread(
					NULL,																			// default security attributes
					0,																				// use default stack size  
					SphericalSynchronizer<Channels>::RunThread<ImageData<float,unsigned __int16> >,	// thread function 
					&_oocSynchronizers[0],															// argument to thread function 
					0,																				// use default creation flags 
					&synchThreadID);																// returns the thread identifier
			}
			if(!_synchronizerHandles[0])	fprintf(stderr,"Failed to create sync thread\n")	,	exit(0);
		}
		for(int i=0;i<2*_oocDCount;i++)
#else // !FROM_IMAGE
		_synchronizerHandles = new HANDLE[_oocDCount+_icDCount];
		for(int i=0;i<_oocDCount;i++)
#endif // FROM_IMAGE
		{
			DWORD synchThreadID;
			_synchronizerHandles[i+1]=CreateThread( 
				NULL,												// default security attributes
				0,													// use default stack size  
				SphericalSynchronizer<Channels>::RunThread<float>,	// thread function 
				&_oocSynchronizers[i+1],							// argument to thread function 
				0,													// use default creation flags 
				&synchThreadID);									// returns the thread identifier
			if(!_synchronizerHandles[i+1])	fprintf(stderr,"Failed to create sync thread\n")	,	exit(0);
		}
		for(int i=0;i<2*_icDCount;i++)
		{
			DWORD synchThreadID;
#if FROM_IMAGE
			_synchronizerHandles[i+2*_oocDCount+1]=CreateThread( 
#else // !FROM_IMAGE
			_synchronizerHandles[i+2*_oocDCount]=CreateThread( 
#endif // FROM_IMAGE
				NULL,												// default security attributes
				0,													// use default stack size  
				SphericalSynchronizer<Channels>::RunThread<float>,	// thread function 
				&_icSynchronizers[i],								// argument to thread function 
				0,													// use default creation flags 
				&synchThreadID);									// returns the thread identifier
#if FROM_IMAGE
			if(!_synchronizerHandles[i+2*_oocDCount+1])	fprintf(stderr,"Failed to create sync thread\n")	,	exit(0);
#else // !FROM_IMAGE
			if(!_synchronizerHandles[i+2*_oocDCount])	fprintf(stderr,"Failed to create sync thread\n")	,	exit(0);
#endif // FROM_IMAGE
		}
#if FROM_IMAGE
		int* oocPorts = new int[2*_oocDCount+1];
#else // !FROM_IMAGE
		int* oocPorts = new int[2*_oocDCount];
#endif // FROM_IMAGE
		int* icPorts = new int[2*_icDCount];
		for(int d=0;d<2*_icDCount;d++)	icPorts[d] = _icSynchronizers[d].port(0);
		for(int i=0;i<_cCount;i++)
		{
#if FROM_IMAGE
			for(int d=0;d<2*_oocDCount+1;d++)	oocPorts[d] = _oocSynchronizers[d].port(i);
			SendOnSocket( clientSockets[i].client , oocPorts , sizeof(int)*2*(_oocDCount+1) );
#else // !FROM_IMAGE
			for(int d=0;d<2*_oocDCount;d++)	oocPorts[d] = _oocSynchronizers[d].port(i);
			SendOnSocket( clientSockets[i].client , oocPorts , sizeof(int)*2*_oocDCount );
#endif // FROM_IMAGE
		}
		delete[] oocPorts;

		_syncSockets = new SOCKET[2*_icDCount];
		for(int i=0;i<2*_icDCount;i++)
		{
			_syncSockets[i] = GetConnectSocket(hostAddress,icPorts[i]);
			if(_syncSockets[i] == INVALID_SOCKET)	return false;
		}
		delete[] icPorts;
	}
	return true;
}
template<int Channels>
void SocketedMultigridServer<Channels>::Run(void)
{
	printf("Running server\n");
#if SEND_STENCILS
	DotProductStencil stencils[4];
#endif // SEND_STENCILS
	SolverInfo<Channels>*  solverInfo = new SolverInfo<Channels>[_oocDCount];
	SolverInfo<Channels>* tSolverInfo = new SolverInfo<Channels>[_oocDCount];
	double zeroAverage[Channels];
	for(int c=0;c<Channels;c++)	zeroAverage[c]=0;
	Vector<float> globalB,globalX;
	int gWidth  = _width  >> (_oocDCount-1);
	int gHeight = _height >> (_oocDCount-1);

	globalB.Resize( gWidth * gHeight * Channels );
	globalX.Resize( gWidth * gHeight * Channels );

	if(_verbose)	MultigridSolver<float,ZERO_DERIVATIVE,2,Channels>::verbose = MultigridSolver<float,ZERO_DERIVATIVE,2,Channels>::FULL_VERBOSE;

#if FROM_IMAGE
	memset(_average, 0 , sizeof(_average) );
#endif // FROM_IMAGE

	for(int ii=0;ii<_vCycles;ii++)
	{
		double t;

		t=Time();
		memset(solverInfo,0,sizeof(SolverInfo<Channels>)*_oocDCount);

#if SOCKET_BACKED_GRID
		for( int i=0 ; i<_cCount ; i++ )
		{
			ReceiveOnSocket( clientSockets[i].clientB , true );
			ReceiveOnSocket( clientSockets[i].clientX , true );
		}
		for( int j=0 ; j<gHeight ; j++ )
			for(int i=0;i<_cCount;i++)
			{
				Vector<float> localB,localX;
				int lWidth  = (clientSockets[i].end-clientSockets[i].start)>>(_oocDCount-1);
				int lOffset = clientSockets[i].start>>(_oocDCount-1);
				localX.Resize( lWidth * Channels );
				localB.Resize( lWidth * Channels );
				ReceiveOnSocket( clientSockets[i].clientB , &localB[0] , sizeof(float) * lWidth * Channels , true );
				ReceiveOnSocket( clientSockets[i].clientX , &localX[0] , sizeof(float) * lWidth * Channels , true );
				for(int c=0;c<Channels;c++)
				{
					float *pGlobalB = &globalB[gWidth*j*Channels+gWidth*c+lOffset];
					float *pGlobalX = &globalX[gWidth*j*Channels+gWidth*c+lOffset];
					float *pLocalB  = &localB [lWidth*c];
					float *pLocalX  = &localX [lWidth*c];
					memcpy(pGlobalB,pLocalB,lWidth*sizeof(float));
					memcpy(pGlobalX,pLocalX,lWidth*sizeof(float));
				}
		}
#endif // SOCKET_BACKED_GRID
		for(int i=0;i<_cCount;i++)
		{
			ReceiveOnSocket ( clientSockets[i].client , tSolverInfo , sizeof(SolverInfo<Channels>) * _oocDCount );
			for(int j=0;j<_oocDCount;j++)
			{
				solverInfo[j].bSquareNorm += tSolverInfo[j].bSquareNorm;
				solverInfo[j].rSquareNorm += tSolverInfo[j].rSquareNorm;
				solverInfo[j].xSquareNorm += tSolverInfo[j].xSquareNorm;
				for(int c=0;c<Channels;c++)	solverInfo[j].solutionSum[c] += tSolverInfo[j].solutionSum[c];
			}
#if FROM_IMAGE
			if(ii==0)
			{
				double avg[Channels];
				ReceiveOnSocket ( clientSockets[i].client , avg ,sizeof(avg) );
				for(int c=0 ; c<Channels ; c++)	_average[c]+=avg[c];
#if SHOW_AVERAGE
				printf("%d]\t%f %f %f\n",i,
					(avg[0]/(clientSockets[i].end-clientSockets[i].start))/_cHeight,
					(avg[1]/(clientSockets[i].end-clientSockets[i].start))/_cHeight,
					(avg[2]/(clientSockets[i].end-clientSockets[i].start))/_cHeight);
#endif // SHOW_AVERAGE
			}
#endif // FROM_IMAGE
		}
#if FROM_IMAGE
		if(ii==0)	for(int c=0 ; c<Channels ; c++)	_average[c] /= _cWidth	,	_average[c] /= _cHeight;
#if SHOW_AVERAGE
	printf("Average: %f %f %f\n",_average[0],_average[1],_average[2]);
#endif // SHOW_AVERAGE
#endif // FROM_IMAGE

		// Get the low-res data off of the sockets
#if !SOCKET_BACKED_GRID
		for(int i=0;i<_cCount;i++)
		{
			Vector<float> localB,localX;
			int lWidth  = (clientSockets[i].end-clientSockets[i].start)>>(_oocDCount-1);
			int lHeight = _height>>(_oocDCount-1);
			int lOffset = clientSockets[i].start>>(_oocDCount-1);
			localX.Resize( lWidth*lHeight*Channels );
			localB.Resize( lWidth*lHeight*Channels );
			ReceiveOnSocket (clientSockets[i].client,&localB[0],sizeof(float)*localB.Dimensions() ) ;
			ReceiveOnSocket (clientSockets[i].client,&localX[0],sizeof(float)*localX.Dimensions() ) ;
			for(int y=0;y<gHeight;y++)
				for(int c=0;c<Channels;c++)
				{
					float *pGlobalB = &globalB[gWidth*y*Channels+gWidth*c+lOffset];
					float *pGlobalX = &globalX[gWidth*y*Channels+gWidth*c+lOffset];
					float *pLocalB  = &localB [lWidth*y*Channels+lWidth*c];
					float *pLocalX  = &localX [lWidth*y*Channels+lWidth*c];
					memcpy(pGlobalB,pLocalB,lWidth*sizeof(float));
					memcpy(pGlobalX,pLocalX,lWidth*sizeof(float));
				}
		}
#endif // !SOCKET_BACKED_GRID
#if SEND_STENCILS
		for(int i=0;i<_cCount;i++) ReceiveOnSocket ( clientSockets[i].client , stencils , sizeof( DotProductStencil ) * 4 );
#endif // SEND_STENCILS

		printf("Out-Of-Core Restriction:    %f\n",Time()-t);
		if(_verbose)
			for(int i=_oocDCount-1;i>=0;i--)
			{
				printf("\tError[%d x %d] %g -> %g\t",_width>>(_oocDCount-1-i),_height>>(_oocDCount-1-i),sqrt(solverInfo[i].bSquareNorm),sqrt(solverInfo[i].rSquareNorm));
				for(int c=0;c<Channels;c++)
				{
					if(c==0)	printf("(");
					else		printf(" ");
					printf("%g",solverInfo[i].solutionSum[c]);
					if(c==Channels-1)	printf(")");
					else				printf("");
				}
				printf("\n");
			}


		// Get the base solution
#if SEND_STENCILS
		if(ii==_vCycles-1)	SolveInCore(stencils[0],stencils[1],stencils[2],stencils[3],globalB,globalX,_average);
		else				SolveInCore(stencils[0],stencils[1],stencils[2],stencils[3],globalB,globalX,zeroAverage);
#else // !SEND_STENCILS
		if(ii==_vCycles-1)	SolveInCore(globalB,globalX,_average);
		else				SolveInCore(globalB,globalX,zeroAverage);
#endif // SEND_STENCILS
		// Send the low-res data to the sockets
#if SOCKET_BACKED_GRID
		for( int j=0 ; j<gHeight ; j++ )
			for(int i=0;i<_cCount;i++)
			{
				Vector<float> localB,localX;
				int lWidth  = (clientSockets[i].end-clientSockets[i].start)>>(_oocDCount-1);
				int lOffset = clientSockets[i].start>>(_oocDCount-1);
				localX.Resize( lWidth * Channels );
				localB.Resize( lWidth * Channels );

				for(int c=0;c<Channels;c++)
				{
					float *pGlobalB = &globalB[gWidth*j*Channels+gWidth*c+lOffset];
					float *pGlobalX = &globalX[gWidth*j*Channels+gWidth*c+lOffset];
					float *pLocalB  = &localB [lWidth*c];
					float *pLocalX  = &localX [lWidth*c];
					memcpy(pLocalB,pGlobalB,lWidth*sizeof(float));
					memcpy(pLocalX,pGlobalX,lWidth*sizeof(float));
				}
				SendOnSocket( clientSockets[i].clientX , &localX[0] , sizeof(float) * lWidth * Channels , true );
			}
#else // !SOCKET_BACKED_GRID
		for(int i=0;i<_cCount;i++)
		{
			Vector<float> localB,localX;
			int lWidth  = (clientSockets[i].end-clientSockets[i].start)>>(_oocDCount-1);
			int lHeight = _height>>(_oocDCount-1);
			int lOffset = clientSockets[i].start>>(_oocDCount-1);
			localX.Resize( lWidth*lHeight*Channels );
			localB.Resize( lWidth*lHeight*Channels );

			for(int y=0;y<gHeight;y++)
				for(int c=0;c<Channels;c++)
				{
					float *pGlobalB = &globalB[gWidth*y*Channels+gWidth*c+lOffset];
					float *pGlobalX = &globalX[gWidth*y*Channels+gWidth*c+lOffset];
					float *pLocalB  = &localB [lWidth*y*Channels+lWidth*c];
					float *pLocalX  = &localX [lWidth*y*Channels+lWidth*c];
					memcpy(pLocalB,pGlobalB,lWidth*sizeof(float));
					memcpy(pLocalX,pGlobalX,lWidth*sizeof(float));
				}
			SendOnSocket (clientSockets[i].client,&localB[0],sizeof(float)*localB.Dimensions() ) ;
			SendOnSocket (clientSockets[i].client,&localX[0],sizeof(float)*localX.Dimensions() ) ;
		}
#endif // SOCKET_BACKED_GRID
		t=Time();
		memset(solverInfo,0,sizeof(SolverInfo<Channels>)*_oocDCount);
		for(int i=0;i<_cCount;i++)
		{
			ReceiveOnSocket ( clientSockets[i].client , tSolverInfo , sizeof(SolverInfo<Channels>)*_oocDCount);
			for(int j=0;j<_oocDCount;j++)
			{
				solverInfo[j].bSquareNorm += tSolverInfo[j].bSquareNorm;
				solverInfo[j].rSquareNorm += tSolverInfo[j].rSquareNorm;
				solverInfo[j].xSquareNorm += tSolverInfo[j].xSquareNorm;
				for(int c=0;c<Channels;c++)	solverInfo[j].solutionSum[c] += tSolverInfo[j].solutionSum[c];
			}
		}
		printf("Out-Of-Core Prolongation:    %f\n",Time()-t);
		if(_verbose)
			for(int i=0;i<_oocDCount;i++)
			{
				printf("\tError[%d x %d] %g -> %g\t",_width>>(_oocDCount-1-i),_height>>(_oocDCount-1-i),sqrt(solverInfo[i].bSquareNorm),sqrt(solverInfo[i].rSquareNorm));
				for(int c=0;c<Channels;c++)
				{
					if(c==0)	printf("(");
					else		printf(" ");
					printf("%g",solverInfo[i].solutionSum[c]);
					if(c==Channels-1)	printf(")");
					else				printf("");
				}
				printf("\n");
			}
	}
	delete[] solverInfo;
	delete[] tSolverInfo;
}
template<int Channels>
#if SEND_STENCILS
void SocketedMultigridServer<Channels>::SolveInCore(DotProductStencil& dotMajor,DotProductStencil& d2DotMajor,
													DotProductStencil& dotMinor,DotProductStencil& d2DotMinor,
													Vector<float>& in,Vector<float>& out,double average[Channels])
#else // !SEND_STENCILS
void SocketedMultigridServer<Channels>::SolveInCore(Vector<float>& in,Vector<float>& out,double average[Channels])
#endif // SEND_STENCILS
{
	double t;
	int w = _width>>(_oocDCount-1) , h = _height>>(_oocDCount-1);
	SocketedMultiGridStreamingSolver<Channels>* solvers;

	Vector<float> lowB,lowX;
	StreamingGrid *B,*X;

	solvers=new SocketedMultiGridStreamingSolver<Channels>[_icDCount];

	for(int i=1;i<_icDCount;i++)	solvers[i].parent=&solvers[i-1];
	for(int i=0;i<_icDCount-1;i++)	solvers[i].rChild =solvers[i].pChild =&solvers[i+1];
	for(int i=0;i<_icDCount;i++)
	{
		solvers[i].bSquareNorm=solvers[i].rSquareNorm=solvers[i].xSquareNorm=0;
		for(int c=0;c<Channels;c++)	solvers[i].solutionSum[c]=0;
		solvers[i].setResidual=true;
	}
	SOCKET leftSocket = INVALID_SOCKET , rightSocket = INVALID_SOCKET;
	if(_spherical)
	{
		int port=0;
	    struct sockaddr_in right;
		SOCKET listenSocket=GetListenSocket(port);
		if (listenSocket == INVALID_SOCKET)	fprintf(stderr,"Failed to generate listen socket\n")	,	exit(0);

		char hostName[512];
		gethostname( hostName , 512 );
		{
			MyWinSock::SystemLock lock;
			hostent* host=gethostbyname(hostName);

			rightSocket = socket(AF_INET, SOCK_STREAM, 0);
			if ( rightSocket == INVALID_SOCKET )	fprintf(stderr,"Error at socket(): %s\n", LastSocketError() )	,	exit(0);

			memset(&right, 0, sizeof(right));
			right.sin_family = AF_INET;
			memcpy(&right.sin_addr,host->h_addr,sizeof(struct in_addr));
			right.sin_port= htons ( port );
		}

		// Accept on the left socket
		AcceptThreadData params;
		HANDLE acceptThread = NULL;
		DWORD accepthThreadID;
		params.listenSocket = &listenSocket;
		params.connectSocket = &leftSocket;
		acceptThread=CreateThread( 
			NULL,				// default security attributes
			0,					// use default stack size  
			AcceptThread,		// thread function 
			&params,			// argument to thread function 
			0,					// use default creation flags 
			&accepthThreadID);	// returns the thread identifier
		if(!acceptThread)	fprintf(stderr,"Failed to create accept thread\n")	,	exit(0);

		while (connect( rightSocket, (const sockaddr*)&right, sizeof(right) ) == SOCKET_ERROR)
		{
			fprintf(stderr,"Connecting...\n");
			Sleep(1);
		}
//		WaitForSingleObject(acceptThread,INFINITE);
		MyWaitForSingleObject( acceptThread , 10000 , "SocketedMultigridServer<Channels>::SolveInCore" );
		CloseHandle(acceptThread);
		if(!params.success)	fprintf(stderr,"Failed to accept self\n")	,	exit(0);

	}
	// BADNESS!!! For spherical domains, the server socket has not been set here.
#if SEND_STENCILS
#if SUPPORT_SHARPEN
	if( _sharpen )	solvers[_icDCount-1].Initialize(dotMajor,d2DotMajor,dotMinor,d2DotMinor,_iWeight,0,w,w,h,_iters,leftSocket,_syncSockets,rightSocket,false,_spherical);
	else			solvers[_icDCount-1].Initialize(dotMajor,d2DotMajor,dotMinor,d2DotMinor,         0,w,w,h,_iters,leftSocket,_syncSockets,rightSocket,false,_spherical);
#else // !SUPPORT_SHARPEN
	solvers[_icDCount-1].Initialize(dotMajor,d2DotMajor,dotMinor,d2DotMinor,0,w,w,h,_iters,leftSocket,_syncSockets,rightSocket,false,_spherical);
#endif // SUPPORT_SHARPEN
#else // !SEND_STENCILS
#if SUPPORT_SHARPEN
	if( _sharpen )	solvers[_icDCount-1].Initialize(_iWeight,0,w,w,h,_iters,leftSocket,_syncSockets,rightSocket,false,_spherical);
	else			solvers[_icDCount-1].Initialize(         0,w,w,h,_iters,leftSocket,_syncSockets,rightSocket,false,_spherical);
#else // !SUPPORT_SHARPEN
	solvers[_icDCount-1].Initialize(0,w,w,h,_iters,leftSocket,_syncSockets,rightSocket,false,_spherical);
#endif // SUPPORT_SHARPEN
#endif // SEND_STENCILS

	lowX.Resize(solvers[0].major*solvers[0].minor*Channels);
	lowB.Resize(solvers[0].major*solvers[0].minor*Channels);
	B=new MemoryBackedGrid(&lowB[0],solvers[0].major*Channels*sizeof(float),solvers[0].minor);
	X=new MemoryBackedGrid(&lowX[0],solvers[0].major*Channels*sizeof(float),solvers[0].minor);
	solvers[_icDCount-1].inB=new MemoryBackedGrid(&in[0],w*Channels*sizeof(float),h);
	solvers[_icDCount-1].inX=new MemoryBackedGrid(&out[0],w*Channels*sizeof(float),h);
	solvers[0].outB=B;

	// Data for the interleaved streaming multigrid
	t=Time();
	// Initialize
	solvers[_icDCount-1].InitRestriction();
	solvers[_icDCount-1].SetRestriction();
	solvers[_icDCount-1].SolveRestriction();

	// Clean up
	delete solvers[_icDCount-1].inB;
	solvers[_icDCount-1].inB=NULL;
	delete solvers[_icDCount-1].inX;
	solvers[_icDCount-1].inX=NULL;

	if(_verbose)
	{
		printf("In-Core Restriction: %f\n",Time()-t);
		for(int i=_icDCount-1;i>=0;i--)
		{
			printf("\tError[%d x %d] %g -> %g\t",solvers[i].major,solvers[i].minor,sqrt(solvers[i].bSquareNorm),sqrt(solvers[i].rSquareNorm));
			for(int c=0;c<Channels;c++)
			{
				if(c==0)	printf("(");
				else		printf(" ");
				printf("%g",((solvers[i].solutionSum[c])/solvers[i].major)/solvers[i].minor);
				if(c==Channels-1)	printf(")");
				else				printf("");
			}
			printf("\n");
		}
	}
	solvers[_icDCount-1].UnSetRestriction();
	// Get the base solution
	SparseMatrix<double> lap;
#if SEND_STENCILS
#if SUPPORT_SHARPEN
	if(_spherical)
		if(!_sharpen || _iWeight==0)	SetSphericalLaplacianMatrix(dotMajor,d2DotMajor,dotMinor,d2DotMinor,lap,solvers[0].major,solvers[0].minor);
		else							SetSphericalLaplacianMatrix(dotMajor,d2DotMajor,dotMinor,d2DotMinor,lap,solvers[0].major,solvers[0].minor,_iWeight);
	else
		if(!_sharpen || _iWeight==0)	SetRectangularLaplacianMatrix(dotMajor,d2DotMajor,dotMinor,d2DotMinor,lap,solvers[0].major,solvers[0].minor);
		else							SetRectangularLaplacianMatrix(dotMajor,d2DotMajor,dotMinor,d2DotMinor,lap,solvers[0].major,solvers[0].minor,_iWeight);
#else // !SUPPORT_SHARPEN
	if(_spherical)	SetSphericalLaplacianMatrix(dotMajor,d2DotMajor,dotMinor,d2DotMinor,lap,solvers[0].major,solvers[0].minor);
	else			SetRectangularLaplacianMatrix(dotMajor,d2DotMajor,dotMinor,d2DotMinor,lap,solvers[0].major,solvers[0].minor);
#endif // SUPPORT_SHARPEN
#else // !SEND_STENCILS
#if SUPPORT_SHARPEN
	if(_spherical)
		if(!_sharpen || _iWeight==0)	SetSphericalLaplacianMatrix(lap,solvers[0].major,solvers[0].minor);
		else							SetSphericalLaplacianMatrix(lap,solvers[0].major,solvers[0].minor,_iWeight);
	else
		if(!_sharpen || _iWeight==0)	FiniteElements2D<double,Type,Degree>::LaplacianMatrix(solvers[0].major,solvers[0].minor,lap,true,false);
		else							FiniteElements2D<double,Type,Degree>::LaplacianMatrix(solvers[0].major,solvers[0].minor,_iWeight,lap,true,false);
#else // !SUPPORT_SHARPEN
	if(_spherical)	SetSphericalLaplacianMatrix(lap,solvers[0].major,solvers[0].minor);
	else			FiniteElements2D<double,Type,Degree>::LaplacianMatrix(solvers[0].major,solvers[0].minor,lap,true,false);
#endif // SUPPORT_SHARPEN
#endif // SEND_STENCILS
	{
		Vector<double> myLowX,myLowB;
		myLowB.Resize(solvers[0].major*solvers[0].minor);
		myLowX.Resize(solvers[0].major*solvers[0].minor);
		lowX.Resize(solvers[0].major*solvers[0].minor*Channels);
		for(int c=0;c<Channels;c++)
		{
			for(int i=0;i<solvers[0].major;i++)
				for(int j=0;j<solvers[0].minor;j++)
					myLowB[i+j*solvers[0].major]=lowB[i+j*solvers[0].major*Channels+c*solvers[0].major];
#if SUPPORT_SHARPEN
		if(!_noCG)
			if(!_sharpen || _iWeight==0) MySolveConjugateGradient<Type>(lap,myLowB,6*int(sqrt(lap.groups+1.0)),myLowX);
//			else						 MySolveConjugateGradient<ZERO_VALUE>(lap,myLowB,6*int(sqrt(lap.groups+1.0)),myLowX);
			else						 MySolveConjugateGradient<Type>(lap,myLowB,6*int(sqrt(lap.groups+1.0)),myLowX);
#else // !SUPPORT_SHARPEN
			if(!_noCG)	MySolveConjugateGradient<Type>(lap,myLowB,6*int(sqrt(lap.groups+1.0)),myLowX);
#endif // SUPPORT_SHARPEN

			for(int i=0;i<solvers[0].major;i++)
				for(int j=0;j<solvers[0].minor;j++)
				{
					lowX[i+j*solvers[0].major*Channels+c*solvers[0].major]=myLowX[i+j*solvers[0].major];
					myLowX[i+j*solvers[0].major]=0;
				}
		}
	}
#if MISHA_DEBUG
	{
		StreamingGrid* mishaDebug = GetWriteStream<float>("debug.low.jpeg",solvers[0].major,solvers[0].minor,100,true);
		for(int y=0;y<solvers[0].minor;y++)
		{
			float* row = (float*)(*mishaDebug)[y];
			for(int x=0;x<solvers[0].major;x++)
				for(int c=0;c<Channels;c++)
					row[x+solvers[0].major*c]=lowX[x+y*solvers[0].major*Channels+c*solvers[0].major];
			mishaDebug->advance();
		}
		delete mishaDebug;
	}
#endif // MISHA_DEBUG

	for(int i=0;i<_icDCount;i++)
	{
		solvers[i].bSquareNorm=solvers[i].rSquareNorm=solvers[i].xSquareNorm=0;
		for(int c=0;c<Channels;c++)	solvers[i].solutionSum[c]=0;
		solvers[i].setResidual=_verbose;
	}

	// Solve the prolongation
	t=Time();
	solvers[0].inX=X;
	// Set the child dependencies
	solvers[_icDCount-1].outX=new MemoryBackedGrid(&out[0],w*Channels*sizeof(float),h);
	solvers[_icDCount-1].InitProlongation();
	solvers[0].SetProlongation();
	// Solve
	solvers[0].SolveProlongation();
	if(_verbose)
	{
		printf("In-Core Prolongation:    %f\n",Time()-t);
		for(int i=0;i<_icDCount;i++)
		{
			printf("\tError[%d x %d] %g -> %g\t",solvers[i].major,solvers[i].minor,sqrt(solvers[i].bSquareNorm),sqrt(solvers[i].rSquareNorm));
			for(int c=0;c<Channels;c++)
			{
				if(c==0)	printf("(");
				else		printf(" ");
				printf("%g",((solvers[i].solutionSum[c])/solvers[i].major)/solvers[i].minor);
				if(c==Channels-1)	printf(")");
				else				printf("");
			}
			printf("\n");
		}
	}
	solvers[0].UnSetProlongation();
	delete solvers[_icDCount-1].outX;
	delete[] solvers;
	delete B;
	delete X;
	B=X=NULL;

	int ww=(long long(_cWidth) *w+_width -1)/_width;
	int hh=(long long(_cHeight)*h+_height-1)/_height;
	for(int c=0;c<Channels;c++)
	{
		double newAverage=0;
		int pCount=0;

		for(int j=0;j<hh;j++)
		{
			float* o=&out[j*w*Channels+c*w];
			for(int i=0;i<ww;i++)
			{
				newAverage+=o[i];
				pCount++;
			}
		}

		newAverage/=pCount;
		newAverage=average[c]-newAverage;

#if SUPPORT_SHARPEN
		if(!_sharpen || _iWeight==0)
			for(int j=0;j<h;j++)
			{
				float* o=&out[j*w*Channels+c*w];
				for(int i=0;i<w;i++)	o[i]+=newAverage;
			}
#else // !SUPPORT_SHARPEN
		for(int j=0;j<h;j++)
		{
			float* o=&out[j*w*Channels+c*w];
			for(int i=0;i<w;i++)	o[i]+=newAverage;
		}
#endif // SUPPORT_SHARPEN
	}
#if MISHA_DEBUG
	{
		StreamingGrid* mishaDebug = GetWriteStream<float>("debug.high.jpeg",w,h,true,100);
		for(int y=0;y<h;y++)
		{
			float* row = (float*)(*mishaDebug)[y];
			for(int x=0;x<w;x++)
				for(int c=0;c<Channels;c++)
					row[x+w*c]=out[x+y*w*Channels+c*w];
			mishaDebug->advance();
		}
		delete mishaDebug;
	}
#endif // MISHA_DEBUG

	if(leftSocket!=INVALID_SOCKET)	closesocket(leftSocket);
	if(rightSocket!=INVALID_SOCKET)	closesocket(rightSocket);
}

/////////////////////////////
// SocketedMultigridClient //
/////////////////////////////
#if !USE_MY_WIN_SOCK
template<int Channels>
int SocketedMultigridClient<Channels>::_wsaCount = 0;

template<int Channels>
WSADATA SocketedMultigridClient<Channels>::_wsaData;
#endif // !USE_MY_WIN_SOCK
template<int Channels>
SocketedMultigridClient<Channels>::SocketedMultigridClient(void)
{
	// Initialize Winsock
#if USE_MY_WIN_SOCK
	MyWinSock::Load();
#else // !USE_MY_WIN_SOCK
	if(!_wsaCount)	if(WSAStartup(MAKEWORD(2,2), &_wsaData))	fprintf(stderr,"WSAStartup failed: %s\n", LastSocketError())	,	exit(0);
	_wsaCount++;
#endif
	_serverSocket	= INVALID_SOCKET;
	_leftSocket		= INVALID_SOCKET;
	_rightSocket	= INVALID_SOCKET;
#if SOCKET_BACKED_GRID
	_serverXSocket	= INVALID_SOCKET;
	_serverBSocket	= INVALID_SOCKET;
#endif // SOCKET_BACKED_GRID
	_syncSockets	= NULL;
}
template<int Channels>
SocketedMultigridClient<Channels>::~SocketedMultigridClient(void)
{
	if(_spherical)
		if(_syncSockets)
		{
#if FROM_IMAGE
			for(int i=0;i<_dCount+1;i++)	if(_syncSockets[i] != INVALID_SOCKET)	closesocket(_syncSockets[i]);
#else // !FROM_IMAGE
			for(int i=0;i<_dCount;i++)	if(_syncSockets[i] != INVALID_SOCKET)	closesocket(_syncSockets[i]);
#endif // FROM_IMAGE
			delete[] _syncSockets;
		}
	if( _serverSocket != INVALID_SOCKET )	closesocket( _serverSocket );
	if( _leftSocket   != INVALID_SOCKET )	closesocket( _leftSocket );
	if( _rightSocket  != INVALID_SOCKET )	closesocket( _rightSocket );
#if SOCKET_BACKED_GRID
	if( _serverXSocket != INVALID_SOCKET )	closesocket( _serverXSocket );
	if( _serverBSocket != INVALID_SOCKET )	closesocket( _serverBSocket );
#endif // SOCKET_BACKED_GRID
	_syncSockets	= NULL;
	_serverSocket	= INVALID_SOCKET;
#if SOCKET_BACKED_GRID
	_serverXSocket	= INVALID_SOCKET;
	_serverBSocket	= INVALID_SOCKET;
#endif // SOCKET_BACKED_GRID
	_leftSocket		= INVALID_SOCKET;
	_rightSocket	= INVALID_SOCKET;
#if USE_MY_WIN_SOCK
	MyWinSock::UnLoad();
#else // !USE_MY_WIN_SOCK
	_wsaCount--;
	if(!_wsaCount)	WSACleanup();
#endif // USE_MY_WIN_SOCK
}
template<int Channels> int SocketedMultigridClient<Channels>::start		(void)	const	{ return _start; }
template<int Channels> int SocketedMultigridClient<Channels>::end		(void)	const	{ return _end; }
template<int Channels> int SocketedMultigridClient<Channels>::width		(void)	const	{ return _width; }
template<int Channels> int SocketedMultigridClient<Channels>::height	(void)	const	{ return _height; }
template<int Channels> int SocketedMultigridClient<Channels>::size		(void)	const	{ return _end-_start; }
template<int Channels> int SocketedMultigridClient<Channels>::cEnd		(void)	const	{ return _cEnd; }
template<int Channels> int SocketedMultigridClient<Channels>::cHeight	(void)	const	{ return _cHeight; }
template<int Channels> int SocketedMultigridClient<Channels>::cSize		(void)	const	{ return _cEnd-_start; }
template<int Channels> int SocketedMultigridClient<Channels>::iters		(void)	const	{ return _iters; }
template<int Channels> int SocketedMultigridClient<Channels>::vCycles	(void)	const	{ return _vCycles; }
template<int Channels> int SocketedMultigridClient<Channels>::quality	(void)	const	{ return _quality; }
template<int Channels> bool SocketedMultigridClient<Channels>::sharpen	(void)	const	{ return _sharpen; }

template<int Channels>
#if FROM_IMAGE
bool SocketedMultigridClient<Channels>::SetUp(char* address,int port,int index,int w,int h)
#else // !FROM_IMAGE
bool SocketedMultigridClient<Channels>::SetUp(char* address,int port,int start)
#endif // FROM_IMAGE
{
	SOCKET listenSocket = INVALID_SOCKET;
    struct sockaddr_in local,right;

#if !FROM_IMAGE
	_start=start;
#endif // !FROM_IMAGE
#if SOCKET_BACKED_GRID
	ClientSocketInfo info;
	info.index = index;
	printf( "Connecting" );
	_serverSocket  = GetConnectSocket( address , port , 500 , true );
	printf( "\nConnecting X" );
	_serverXSocket = GetConnectSocket( address , port , 500 , true );
	printf( "\nConnecting B" );
	_serverBSocket = GetConnectSocket( address , port , 500 , true );
	printf( "\n" );
	info.type = ClientSocketInfo::DATA;
	SendOnSocket( _serverSocket  , &info , sizeof(info) );
	info.type = ClientSocketInfo::X;
	SendOnSocket( _serverXSocket , &info , sizeof(info) );
	info.type = ClientSocketInfo::B;
	SendOnSocket( _serverBSocket , &info , sizeof(info) );

#else // !SOCKET_BACKED_GRID
	_serverSocket = socket( AF_INET , SOCK_STREAM , 0 );
	if ( _serverSocket == INVALID_SOCKET )
	{
		fprintf(stderr,"Error at socket(): %s\n", LastSocketError() );
		return false;
	}

    memset(&local, 0, sizeof(local));
	inet_aton(address, &local.sin_addr);
    local.sin_port = htons(port);
    local.sin_family = AF_INET;

	// Connect to server.
	printf("Connecting");
	while (connect( _serverSocket, (const sockaddr*)&local, sizeof(local) ) == SOCKET_ERROR) 
	{
		printf(".");
		Sleep(1000);
	}
	printf("\n");
	{
		int val=1;
		setsockopt( _serverSocket , IPPROTO_TCP , TCP_NODELAY , (char*)&val , sizeof(val) );
	}
#endif // SOCKET_BACKED_GRID

	// Send the start position
	{
		ClientData cd;
#if FROM_IMAGE
		cd.index=index;
		cd.width=w;
		cd.height=h;
#else // !FROM_IMAGE
		cd.start=_start;
#endif // FROM_IMAGE
		SendOnSocket ( _serverSocket , &cd , sizeof(cd) );
	}

	// Get back the other dimension data
	{
		ServerData sd;
		ReceiveOnSocket ( _serverSocket,&sd,sizeof(sd) );
		_width		= sd.width;
		_height		= sd.height;
		_end		= sd.end;
		_cEnd		= sd.cEnd;
#if FROM_IMAGE
		_start		= sd.start;
		_cHeight	= h;
#else // !FROM_IMAGE
		_cHeight	= sd.cHeight;
#endif // FROM_IMAGE
		_iters		= sd.iters;
		_vCycles	= sd.vCycles;
		_verbose	= sd.verbose;
		_lanes		= sd.lanes;
		_spherical	= sd.spherical;
		_progress	= sd.progress;
		_quality	= sd.quality;
#if SUPPORT_SHARPEN
		_sharpen	= sd.sharpen;
		_iWeight	= sd.iWeight;
		_gWeight	= sd.gWeight;
#endif // SUPPORT_SHARPEN
#if DEBUG_FLAGS
		ReceiveOnSocket ( _serverSocket,&debugFlags,sizeof(debugFlags) );
#endif // DEBUG_FLAGS
	}
	// Send information about the left socket
	if(_start!=0 || _spherical)
	{
		int port=0;
		listenSocket=GetListenSocket(port);
		if (listenSocket == INVALID_SOCKET)	return false;

		{
			MyWinSock::SystemLock lock;
			char hostName[512];
			gethostname(hostName,512);
			hostent* host=gethostbyname(hostName);

			SendOnSocket ( _serverSocket,host->h_addr,sizeof(struct in_addr) );
			SendOnSocket ( _serverSocket,&port,sizeof(port) );
		}
	}
	// Get the information about the right socket
	if(_end!=_width || _spherical)
	{
		int port=0;
		_rightSocket = socket(AF_INET, SOCK_STREAM, 0);
		if ( _rightSocket == INVALID_SOCKET )
		{
			fprintf(stderr,"Error at socket(): %s\n", LastSocketError() );
			return false;
		}
		memset(&right, 0, sizeof(right));
		right.sin_family = AF_INET;
		ReceiveOnSocket ( _serverSocket,&right.sin_addr,sizeof(right.sin_addr) );
		ReceiveOnSocket ( _serverSocket,&port,sizeof(port) );
		right.sin_port= htons ( port );
	}

	// Accept on the left socket
	AcceptThreadData params;
	HANDLE acceptThread = NULL;
	DWORD accepthThreadID;
	params.listenSocket = &listenSocket;
	params.connectSocket = &_leftSocket;
	if(_start!=0 || _spherical)
	{
		acceptThread=CreateThread( 
			NULL,				// default security attributes
			0,					// use default stack size  
			AcceptThread,		// thread function 
			&params,			// argument to thread function 
			0,					// use default creation flags 
			&accepthThreadID);	// returns the thread identifier
		if(!acceptThread)
		{
			fprintf(stderr,"Failed to create accept thread\n");
			exit(0);
		}
	}

	// Connect on the right socket
	if(_end!=_width || _spherical)
	{
		while (connect( _rightSocket, (const sockaddr*)&right, sizeof(right) ) == SOCKET_ERROR)
		{
			fprintf(stderr,"Error at connect(): %s\n", LastSocketError() );
			fprintf(stderr,"retrying...\n");
			Sleep(1000);
		}
		int val=1;
		setsockopt(_rightSocket,IPPROTO_TCP,TCP_NODELAY,(char*)&val,sizeof(val));
	}

	if( acceptThread )
	{
		MyWaitForSingleObject( acceptThread , 10000 , "SocketedMultigridClient<Channels>::SetUp" );
		CloseHandle(acceptThread);
		if(!params.success)	return false;
	}

	// Get the depth of the multigrid solver
	ReceiveOnSocket ( _serverSocket,&_dCount,sizeof(_dCount) );

	if(_spherical)
	{
#if FROM_IMAGE
		int* ports = new int[2*_dCount+1];
		_syncSockets = new SOCKET[2*_dCount+1];
		ReceiveOnSocket ( _serverSocket , ports , sizeof(int)*2*(_dCount+1) );
		for(int i=0;i<2*_dCount+1;i++)
#else // !FROM_IMAGE
		int* ports = new int[2*_dCount];
		_syncSockets = new SOCKET[2*_dCount];
		ReceiveOnSocket ( _serverSocket , ports , sizeof(int)*2*_dCount );
		for(int i=0;i<2*_dCount;i++)
#endif // FROM_IMAGE
		{
			_syncSockets[i] = GetConnectSocket(address,ports[i]);
			if(_syncSockets[i] == INVALID_SOCKET)	return false;
		}

		delete[] ports;
	}
	return true;
}
template<int Channels>
#if FROM_IMAGE
template<class PixelType,class LabelType,class StorageType>
void SocketedMultigridClient<Channels>::Run(StreamingGrid* pixels,StreamingGrid* labels,StreamingGrid* out)
#else // !FROM_IMAGE
template<class StorageType>
void SocketedMultigridClient<Channels>::Run(StreamingGrid* laplacian,StreamingGrid* out)
#endif // FROM_IMAGE
{
	double t;
#if SEND_STENCILS
	DotProductStencil stencils[4];
#endif // SEND_STENCILS
	SolverInfo<Channels>* solverInfo = new SolverInfo<Channels>[_dCount];
#if FROM_IMAGE
	SocketedStreamingDivergence<Channels,PixelType,LabelType,StorageType>* sLaplacian = new SocketedStreamingDivergence<Channels,PixelType,LabelType,StorageType>();
#endif // FROM_IMAGE
	SocketedMultiGridStreamingSolver<Channels,StorageType>* solvers;
	SocketedStreamingSolver<Channels>::server.Reset();

#if !SOCKET_BACKED_GRID
	Vector<float> lowB,lowX;
#endif // !SOCKET_BACKED_GRID
	StreamingGrid *B,*X;
	// Also, the template should be StorageType, not Real
	MultiStreamIOClient* in;

	if(_vCycles>1)	in=new MultiStreamIOClient( size()*sizeof(float)*Channels , _height , STREAMING_GRID_BUFFER_MULTIPLIER , NULL , true);
	else			in=NULL;
 
	solvers = new SocketedMultiGridStreamingSolver<Channels,StorageType>[_dCount];
	solvers[_dCount-1].showProgress=_progress;
	for(int i=1;i<_dCount;i++)		solvers[i].parent = &solvers[i-1];
	for(int i=0;i<_dCount-1;i++)	solvers[i].rChild =  solvers[i].pChild = &solvers[i+1];
	for(int i=0;i<_dCount;i++)		solvers[i].laneNum = _lanes;


#if FROM_IMAGE
	sLaplacian->parent = & solvers[_dCount-1];
	solvers[_dCount-1].rChild = sLaplacian;
	sLaplacian->Initialize(pixels,labels,_start,_end,_width,_height,_iters,_leftSocket,_syncSockets,_rightSocket,true,_spherical);
#else // !FROM_IMAGE
	solvers[_dCount-1].Initialize(_start,_end,_width,_height,_iters,_leftSocket,_syncSockets,_rightSocket,true,_spherical);
#endif // FROM_IMAGE

#if SOCKET_BACKED_GRID
#if NON_BLOCKING_SOCKETS
	B = new SocketBackedGrid( _serverBSocket , solvers[0].size() * Channels * sizeof( float ) , solvers[0].minor , false );
	X = new SocketBackedGrid( _serverXSocket , solvers[0].size() * Channels * sizeof( float ) , solvers[0].minor , false );
#else // !NON_BLOCKING_SOCKETS
	B = new SocketBackedGrid( _serverBSocket , solvers[0].size() * Channels * sizeof( float ) , solvers[0].minor , true );
	X = new SocketBackedGrid( _serverXSocket , solvers[0].size() * Channels * sizeof( float ) , solvers[0].minor , true );
#endif // NON_BLOCKING_SOCKETS
//	B = new MemoryBackedGrid( solvers[0].size() * Channels * sizeof( float ) , solvers[0].minor );
//	X = new MemoryBackedGrid( solvers[0].size() * Channels * sizeof( float ) , solvers[0].minor );
#else // !SOCKET_BACKED_GRID
	lowX.Resize(solvers[0].size()*solvers[0].minor*Channels);
	lowB.Resize(solvers[0].size()*solvers[0].minor*Channels);
	X=new MemoryBackedGrid(&lowX[0],solvers[0].size()*Channels*sizeof(float),solvers[0].minor);
	B=new MemoryBackedGrid(&lowB[0],solvers[0].size()*Channels*sizeof(float),solvers[0].minor);
#endif // SOCKET_BACKED_GRID
	for(int ii=0;ii<_vCycles;ii++)
	{
#if FROM_IMAGE
		if(!ii)
		{
			pixels->SetServer(&SocketedStreamingSolver<Channels>::server);
			labels->SetServer(&SocketedStreamingSolver<Channels>::server);
		}
#else // !FROM_IMAGE
		laplacian->SetServer(&SocketedStreamingSolver<Channels>::server);
#endif // FROM_IMAGE
		/////////////////
		// RESTRICTION //
		/////////////////
		for(int i=0;i<_dCount;i++)
		{
			solvers[i].bSquareNorm=solvers[i].rSquareNorm=solvers[i].xSquareNorm=0;
			for ( int c = 0 ; c < Channels ; c++ )	solvers[i].solutionSum[c]=0;
			solvers[i].setResidual=true;
		}
#if !FROM_IMAGE
		solvers[_dCount-1].inB=laplacian;
#endif // !FROM_IMAGE
		if(ii)	solvers[_dCount-1].inX=in;
		else	solvers[_dCount-1].inX=NULL;
		solvers[_dCount-1].outX=NULL;
		solvers[0].inX=NULL;
		solvers[0].outB=B;
		solvers[0].outX=X;
#if FROM_IMAGE
		solvers[_dCount-1].rChild=sLaplacian;
#else // !FROM_IMAGE
		solvers[_dCount-1].rChild=NULL;
#endif // !FROM_IMAGE
		// Data for the interleaved streaming multigrid
		t=Time();
		// Initialize
#if FROM_IMAGE
		if(sLaplacian)	sLaplacian->InitRestriction()	,	sLaplacian->SetRestriction();
		else	solvers[_dCount-1].InitRestriction()	,	solvers[_dCount-1].SetRestriction();
#else // !FROM_IMAGE
		solvers[_dCount-1].InitRestriction();
		solvers[_dCount-1].SetRestriction();
#endif // FROM_IMAGE
		// Solve
		// BADNESS!!! Why do I have to comment this out?
//		if(ii)	in->SetServer(&SocketedStreamingSolver<Channels>::server);
		SocketedStreamingSolver<Channels>::server.StartIO();
#if FROM_IMAGE
		if(sLaplacian)	sLaplacian->SolveRestriction();
		else			solvers[_dCount-1].SolveRestriction();
#else // !FROM_IMAGE
		solvers[_dCount-1].SolveRestriction();
#endif // FROM_IMAGE
		SocketedStreamingSolver<Channels>::server.Reset();

		printf("Out-Of-Core Restriction:    %f\n",Time()-t);
		for (int i = 0 ; i < _dCount ; i++)
		{
			solverInfo[i].bSquareNorm = solvers[i].bSquareNorm;
			solverInfo[i].rSquareNorm = solvers[i].rSquareNorm;
			solverInfo[i].xSquareNorm = solvers[i].xSquareNorm;
			for(int c=0;c<Channels;c++)
			{
				solverInfo[i].solutionSum[c] = solvers[i].solutionSum[c];
				solverInfo[i].solutionSum[c] /= solvers[i].major;
				solverInfo[i].solutionSum[c] /= solvers[i].minor;
			}
		}
		if(_verbose)
			for(int i=_dCount-1;i>=0;i--)
			{
				printf("\tError[%d x %d] %g -> %g\t",solvers[i].size(),solvers[i].minor,sqrt(solvers[i].bSquareNorm),sqrt(solvers[i].rSquareNorm));
				for(int c=0;c<Channels;c++)
				{
					if(c==0)	printf("(");
					else		printf(" ");
					printf("%g",solverInfo[i].solutionSum[c]);
					if(c==Channels-1)	printf(")");
					else				printf("");
				}
				printf("\n");
			}
		SendOnSocket ( _serverSocket , solverInfo , sizeof(SolverInfo<Channels>)*_dCount );

		if(ii==_vCycles-1 && in)	delete in,	in=NULL;
#if FROM_IMAGE
		if(sLaplacian)
		{
			SendOnSocket ( _serverSocket , sLaplacian->average , sizeof(sLaplacian->average) );
			sLaplacian->UnSetRestriction();
			delete sLaplacian;
			sLaplacian=NULL;
		}
		else	solvers[_dCount-1].UnSetRestriction();
#else // !FROM_IMAGE
		solvers[_dCount-1].UnSetRestriction();
#endif // FROM_IMAGE

		t=Time();
#if SOCKET_BACKED_GRID
#else // !SOCKET_BACKED_GRID
		// Send the low res stuff over to the server so that it can perform the in-core solve.
		SendOnSocket (_serverSocket,&lowB[0],sizeof(float)*lowB.Dimensions() ) ;
		SendOnSocket (_serverSocket,&lowX[0],sizeof(float)*lowX.Dimensions() ) ;
#endif // !SOCKET_BACKED_GRID

#if SEND_STENCILS
		stencils[0] = solvers[0].dotMajor;
		stencils[1] = solvers[0].d2DotMajor;
		stencils[2] = solvers[0].dotMinor;
		stencils[3] = solvers[0].d2DotMinor;
		SendOnSocket ( _serverSocket , stencils , sizeof( DotProductStencil ) * 4 );
#endif // SEND_STENCILS

		// Get back the in-core solution.
#if SOCKET_BACKED_GRID
#else // !SOCKET_BACKED_GRID
		ReceiveOnSocket	(_serverSocket,&lowB[0],sizeof(float)*lowB.Dimensions() ) ;
		ReceiveOnSocket	(_serverSocket,&lowX[0],sizeof(float)*lowX.Dimensions() ) ;
#endif // !SOCKET_BACKED_GRID
		printf("In-Core:    %f\n",Time()-t);

		//////////////////
		// PROLONGATION //
		//////////////////
		for(int i=0;i<_dCount;i++)
		{
			solvers[i].bSquareNorm=solvers[i].rSquareNorm=solvers[i].xSquareNorm=0;
			for ( int c = 0 ; c < Channels ; c++ )	solvers[i].solutionSum[c]=0;
			solvers[i].setResidual=_verbose;
		}
		solvers[_dCount-1].inX=NULL;
		solvers[0].outB=NULL;
		solvers[0].outX=NULL;
		if(ii==_vCycles-1)	solvers[_dCount-1].outX=out;
		else				solvers[_dCount-1].outX=in;
		solvers[0].inX=X;

		// Solve the prolongation
		t=Time();
		// Set the child dependencies
		// Initialize
		solvers[_dCount-1].InitProlongation();
		solvers[0].SetProlongation();
		// Solve
		// BADNESS!!! Why do I have to comment this out?
//		if(ii<_vCycles-1)	in->SetServer(&StreamingSolver<Real,Type,Degree,Channels>::server);
		SocketedStreamingSolver<Channels>::server.StartIO();
		solvers[0].SolveProlongation();
		SocketedStreamingSolver<Channels>::server.Reset();

		printf("Out-Of-Core Prolongation:    %f\n",Time()-t);
		for (int i = 0 ; i < _dCount ; i++)
		{
			solverInfo[i].bSquareNorm = solvers[i].bSquareNorm;
			solverInfo[i].rSquareNorm = solvers[i].rSquareNorm;
			solverInfo[i].xSquareNorm = solvers[i].xSquareNorm;
			for(int c=0;c<Channels;c++)
			{
				solverInfo[i].solutionSum[c] = solvers[i].solutionSum[c];
				solverInfo[i].solutionSum[c] /= solvers[i].major;
				solverInfo[i].solutionSum[c] /= solvers[i].minor;
			}
		}
		if(_verbose)
			for(int i=0;i<_dCount;i++)
			{
				printf("\tError[%d x %d] %g -> %g\t",solvers[i].size(),solvers[i].minor,sqrt(solvers[i].bSquareNorm),sqrt(solvers[i].rSquareNorm));
				for(int c=0;c<Channels;c++)
				{
					if(c==0)	printf("(");
					else		printf(" ");
					printf("%g",solverInfo[i].solutionSum[c]);
					if(c==Channels-1)	printf(")");
					else				printf("");
				}
				printf("\n");
			}
		SendOnSocket ( _serverSocket , solverInfo , sizeof(SolverInfo<Channels>)*_dCount );

		solvers[0].UnSetProlongation();
	}
	delete[] solvers;
	delete B;
	delete X;
	B = X = NULL;
	delete[] solverInfo;
}
#if SUPPORT_SHARPEN
template<int Channels>
template<class PixelType,class StorageType>
void SocketedMultigridClient<Channels>::Run(StreamingGrid* pixels,StreamingGrid* out)
{
	double t;
#if SEND_STENCILS
	DotProductStencil stencils[4];
#endif // SEND_STENCILS
	SolverInfo<Channels>* solverInfo = new SolverInfo<Channels>[_dCount];
	SocketedStreamingLaplacian<Channels,PixelType,StorageType>* sLaplacian = new SocketedStreamingLaplacian<Channels,PixelType,StorageType>();
	SocketedMultiGridStreamingSolver<Channels,StorageType>* solvers;
	SocketedStreamingSolver<Channels>::server.Reset();

	Vector<float> lowB,lowX;
	StreamingGrid *B,*X;
	// Also, the template should be StorageType, not Real
	MultiStreamIOClient* in;

	if(_vCycles>1)	in=new MultiStreamIOClient( size()*sizeof(float)*Channels , _height , STREAMING_GRID_BUFFER_MULTIPLIER , NULL , true);
	else			in=NULL;
 
	solvers = new SocketedMultiGridStreamingSolver<Channels,StorageType>[_dCount];
	solvers[_dCount-1].showProgress=_progress;
	for(int i=1;i<_dCount;i++)		solvers[i].parent = &solvers[i-1];
	for(int i=0;i<_dCount-1;i++)	solvers[i].rChild =  solvers[i].pChild = &solvers[i+1];
	for(int i=0;i<_dCount;i++)		solvers[i].laneNum = _lanes;


	sLaplacian->parent = & solvers[_dCount-1];
	solvers[_dCount-1].rChild = sLaplacian;
	sLaplacian->Initialize(pixels,_iWeight,_gWeight,_start,_end,_width,_height,_iters,_leftSocket,_syncSockets,_rightSocket,true,_spherical);

	lowX.Resize(solvers[0].size()*solvers[0].minor*Channels);
	lowB.Resize(solvers[0].size()*solvers[0].minor*Channels);
	B=new MemoryBackedGrid(&lowB[0],solvers[0].size()*Channels*sizeof(float),solvers[0].minor);
	X=new MemoryBackedGrid(&lowX[0],solvers[0].size()*Channels*sizeof(float),solvers[0].minor);
	for(int ii=0;ii<_vCycles;ii++)
	{
		if(!ii)
			pixels->SetServer(&SocketedStreamingSolver<Channels>::server);
		/////////////////
		// RESTRICTION //
		/////////////////
		for(int i=0;i<_dCount;i++)
		{
			solvers[i].bSquareNorm=solvers[i].rSquareNorm=solvers[i].xSquareNorm=0;
			for ( int c = 0 ; c < Channels ; c++ )	solvers[i].solutionSum[c]=0;
			solvers[i].setResidual=true;
		}
		if(ii)	solvers[_dCount-1].inX=in;
		else	solvers[_dCount-1].inX=NULL;
		solvers[_dCount-1].outX=NULL;
		solvers[0].inX=NULL;
		solvers[0].outB=B;
		solvers[0].outX=X;
		solvers[_dCount-1].rChild=sLaplacian;
		// Data for the interleaved streaming multigrid
		t=Time();
		// Initialize
		if(sLaplacian)	sLaplacian->InitRestriction()	,	sLaplacian->SetRestriction();
		else	solvers[_dCount-1].InitRestriction()	,	solvers[_dCount-1].SetRestriction();
		// Solve
		// BADNESS!!! Why do I have to comment this out?
//		if(ii)	in->SetServer(&SocketedStreamingSolver<Channels>::server);
		SocketedStreamingSolver<Channels>::server.StartIO();
		if(sLaplacian)	sLaplacian->SolveRestriction();
		else			solvers[_dCount-1].SolveRestriction();
		SocketedStreamingSolver<Channels>::server.Reset();


		printf("Out-Of-Core Restriction:    %f\n",Time()-t);
		for (int i = 0 ; i < _dCount ; i++)
		{
			solverInfo[i].bSquareNorm = solvers[i].bSquareNorm;
			solverInfo[i].rSquareNorm = solvers[i].rSquareNorm;
			solverInfo[i].xSquareNorm = solvers[i].xSquareNorm;
			for(int c=0;c<Channels;c++)
			{
				solverInfo[i].solutionSum[c] = solvers[i].solutionSum[c];
				solverInfo[i].solutionSum[c] /= solvers[i].major;
				solverInfo[i].solutionSum[c] /= solvers[i].minor;
			}
		}

		if(_verbose)
			for(int i=_dCount-1;i>=0;i--)
			{
				printf("\tError[%d x %d] %g -> %g\t",solvers[i].size(),solvers[i].minor,sqrt(solvers[i].bSquareNorm),sqrt(solvers[i].rSquareNorm));
				for(int c=0;c<Channels;c++)
				{
					if(c==0)	printf("(");
					else		printf(" ");
					printf("%g",solverInfo[i].solutionSum[c]);
					if(c==Channels-1)	printf(")");
					else				printf("");
				}
				printf("\n");
			}
		SendOnSocket ( _serverSocket , solverInfo , sizeof(SolverInfo<Channels>)*_dCount );

		if(ii==_vCycles-1 && in)	delete in,	in=NULL;
		if(sLaplacian)
		{
			SendOnSocket ( _serverSocket , sLaplacian->average , sizeof(sLaplacian->average) );
			sLaplacian->UnSetRestriction();
			delete sLaplacian;
			sLaplacian=NULL;
		}
		else	solvers[_dCount-1].UnSetRestriction();

		t=Time();
		// Send the low res stuff over to the server so that it can perform the in-core solve.
		SendOnSocket (_serverSocket,&lowB[0],sizeof(float)*lowB.Dimensions() ) ;
		SendOnSocket (_serverSocket,&lowX[0],sizeof(float)*lowX.Dimensions() ) ;

#if SEND_STENCILS
		stencils[0] = solvers[0].dotMajor;
		stencils[1] = solvers[0].d2DotMajor;
		stencils[2] = solvers[0].dotMinor;
		stencils[3] = solvers[0].d2DotMinor;
		SendOnSocket ( _serverSocket , stencils , sizeof( DotProductStencil ) * 4 );
#endif // SEND_STENCILS

		// Get back the in-core solution.
		ReceiveOnSocket	(_serverSocket,&lowB[0],sizeof(float)*lowB.Dimensions() ) ;
		ReceiveOnSocket	(_serverSocket,&lowX[0],sizeof(float)*lowX.Dimensions() ) ;
		printf("In-Core:    %f\n",Time()-t);



		//////////////////
		// PROLONGATION //
		//////////////////
		for(int i=0;i<_dCount;i++)
		{
			solvers[i].bSquareNorm=solvers[i].rSquareNorm=solvers[i].xSquareNorm=0;
			for ( int c = 0 ; c < Channels ; c++ )	solvers[i].solutionSum[c]=0;
			solvers[i].setResidual=_verbose;
		}
		solvers[_dCount-1].inX=NULL;
		solvers[0].outB=NULL;
		solvers[0].outX=NULL;
		if(ii==_vCycles-1)	solvers[_dCount-1].outX=out;
		else				solvers[_dCount-1].outX=in;
		solvers[0].inX=X;

		// Solve the prolongation
		t=Time();
		// Set the child dependencies
		// Initialize
		solvers[_dCount-1].InitProlongation();
		solvers[0].SetProlongation();
		// Solve
		// BADNESS!!! Why do I have to comment this out?
//		if(ii<_vCycles-1)	in->SetServer(&StreamingSolver<Real,Type,Degree,Channels>::server);
		SocketedStreamingSolver<Channels>::server.StartIO();
		solvers[0].SolveProlongation();
		SocketedStreamingSolver<Channels>::server.Reset();

		printf("Out-Of-Core Prolongation:    %f\n",Time()-t);
		for (int i = 0 ; i < _dCount ; i++)
		{
			solverInfo[i].bSquareNorm = solvers[i].bSquareNorm;
			solverInfo[i].rSquareNorm = solvers[i].rSquareNorm;
			solverInfo[i].xSquareNorm = solvers[i].xSquareNorm;
			for(int c=0;c<Channels;c++)
			{
				solverInfo[i].solutionSum[c] = solvers[i].solutionSum[c];
				solverInfo[i].solutionSum[c] /= solvers[i].major;
				solverInfo[i].solutionSum[c] /= solvers[i].minor;
			}
		}
		if(_verbose)
			for(int i=0;i<_dCount;i++)
			{
				printf("\tError[%d x %d] %g -> %g\t",solvers[i].size(),solvers[i].minor,sqrt(solvers[i].bSquareNorm),sqrt(solvers[i].rSquareNorm));
				for(int c=0;c<Channels;c++)
				{
					if(c==0)	printf("(");
					else		printf(" ");
					printf("%g",solverInfo[i].solutionSum[c]);
					if(c==Channels-1)	printf(")");
					else				printf("");
				}
				printf("\n");
			}
		SendOnSocket ( _serverSocket , solverInfo , sizeof(SolverInfo<Channels>)*_dCount );

		solvers[0].UnSetProlongation();
	}
	delete[] solvers;
	delete B;
	delete X;
	B=X=NULL;
	delete[] solverInfo;
}
#endif // SUPPORT_SHARPEN