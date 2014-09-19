#ifndef SOCKETED_MULTIGRID_SOLVER_INCLUDED
#define SOCKETED_MULTIGRID_SOLVER_INCLUDED
#include "StreamingSolver.h"
#include "SocketedStreamingSolver.h"
#define NEW_CLIENT_SERVER 1
#if NEW_CLIENT_SERVER
#include "SocketedMultigrid/MultigridServer.h"
#include "SocketedMultigrid/SocketData.h"
#else // !NEW_CLIENT_SERVER

#define SOCKET_BACKED_GRID 0
#define USE_MY_WIN_SOCK 1
#define FROM_IMAGE 1
#define SUPPORT_SHARPEN 1	// Only should be enabled if FROM_IMAGE is
#define SEND_STENCILS 1

class AcceptThreadData
{
public:
	SOCKET *listenSocket,*connectSocket;
	bool success;
};
DWORD WINAPI AcceptThread(LPVOID lpParam);

#if SOCKET_BACKED_GRID
class ClientSocketInfo
{
public:
	enum
	{
		DATA,
		X,
		B
	};
	int index;
	int type;
};
#endif // SOCKET_BACKED_GRID
class ClientData
{
public:
#if FROM_IMAGE
	int index;
	int width;
	int height;
#else // !FROM_IMAGE
	int start;
#endif // FROM_IMAGE
};
class ServerData
{
public:
	int width , height;
	int end , cEnd;
#if FROM_IMAGE
	int start;
#else // !FROM_IMAGE
	int cHeight;
#endif // FROM_IMAGE
	int iters,vCycles,quality,lanes;
	bool verbose,spherical,progress;
#if SUPPORT_SHARPEN
	bool sharpen;
	double gWeight,iWeight;
#endif // SUPPORT_SHARPEN
};

template<int Channels>
class SolverInfo
{
public:
	double bSquareNorm , rSquareNorm , xSquareNorm , solutionSum[Channels];
};

class ClientSocket
{
public:
	SOCKET client;
#if SOCKET_BACKED_GRID
	SOCKET clientX , clientB;
#endif // SOCKET_BACKED_GRID
	int start,end;
	struct in_addr address;
	int port;
#if FROM_IMAGE
	int index;
#endif // FROM_IMAGE

	ClientSocket(void)
	{
		client = INVALID_SOCKET;
		start = end = 0;
		port = 0;
	}
	~ClientSocket(void)
	{
		if(client != INVALID_SOCKET)	closesocket(client);
		client = INVALID_SOCKET;
	}
	static int Sort(const void* v1,const void* v2)
	{
		ClientSocket *s1=(ClientSocket*)v1;
		ClientSocket *s2=(ClientSocket*)v2;
#if FROM_IMAGE
		return s1->index-s2->index;
#else // !FROM_IMAGE
		return s1->start-s2->start;
#endif // FROM_IMAGE
	}
};
template<int Channels>
class SphericalSynchronizer
{
	int _width,_cCount;
	ClientSocket* _clientSockets;
	int _rPasses,_pPasses,_vCycles;
public:

	SphericalSynchronizer(void);
	~SphericalSynchronizer(void);
	void init(int width,const int* widths,int cCount,int rPasses,int pPasses,int vCycles);
	int clients(void)		const;
	int port(int cIndex)	const;
	template<class DataType>	void Run(void);

	template<class DataType>
	static DWORD WINAPI RunThread(LPVOID lpParam);
};
template<int Channels>
class SocketedMultigridServer
{
#if !USE_MY_WIN_SOCK
	static int _wsaCount;
	static WSADATA _wsaData;
#endif // !USE_MY_WIN_SOCK
	bool _verbose,_spherical,_noCG;
	int _cCount,_oocDCount,_icDCount;
	int _width,_height,_cWidth,_cHeight;
	int _iters,_vCycles;
	int _inCoreRes,_minMGRes,*_icPorts;
#if SUPPORT_SHARPEN
	bool _sharpen;
	double _gWeight,_iWeight;
#endif // SUPPORT_SHARPEN

	double _average[Channels];
	ClientSocket* clientSockets;
	SphericalSynchronizer<Channels> *_icSynchronizers,*_oocSynchronizers;
	SOCKET* _syncSockets;
	HANDLE* _synchronizerHandles;
#if SEND_STENCILS
	void SolveInCore(
		DotProductStencil& dotMajor,DotProductStencil& d2DotMajor,DotProductStencil& dotMinor,DotProductStencil& d2DotMinor,
		Vector<float>& in,Vector<float>& out,double average[Channels]);
#else // !SEND_STENCILS
	void SolveInCore(Vector<float>& in,Vector<float>& out,double average[Channels]);
#endif // SEND_STENCILS
public:
	SocketedMultigridServer(void);
	~SocketedMultigridServer(void);
#if FROM_IMAGE
#if SUPPORT_SHARPEN
	bool SetUp(int port,int clientCount,int multiGrid,int inCoreRes,int minMGRes,int iters,
		int quality,int lanes,bool verbose,bool spherical,bool sharpen,double iWeight,double gWeight,bool showProgress,bool noCG=false);
#else // !SUPPORT_SHARPEN
	bool SetUp(int port,int clientCount,int multiGrid,int inCoreRes,int minMGRes,int iters,
		int quality,int lanes,bool verbose,bool spherical,bool showProgress,bool noCG=false);
#endif // SUPPORT_SHARPEN
#else // !FROM_IMAGE
	bool SetUp(int port,int clientCount,int width,int height,int multiGrid,int inCoreRes,int minMGRes,double average[Channels],int iters,int cWidth,int cHeight,
		int quality,int lanes,bool verbose,bool spherical,bool showProgress,bool noCG=false);
#endif // FROM_IMAGE
	void Run(void);
};
template<int Channels>
class SocketedMultigridClient
{
#if !USE_MY_WIN_SOCK
	static int _wsaCount;
	static WSADATA _wsaData;
#endif // !USE_MY_WIN_SOCK
	bool _verbose,_spherical,_progress;
	SOCKET _leftSocket,_serverSocket,_rightSocket,*_syncSockets;
#if SOCKET_BACKED_GRID
	SOCKET _serverXSocket , _serverBSocket;
#endif // SOCKET_BACKED_GRID

	int _start,_end,_width,_height,_cEnd,_cHeight,_dCount,_iters,_vCycles,_quality,_lanes;
#if SUPPORT_SHARPEN
	bool _sharpen;
	double _gWeight,_iWeight;
#endif // SUPPORT_SHARPEN

public:
	SocketedMultigridClient(void);
	~SocketedMultigridClient(void);
#if FROM_IMAGE
	bool SetUp(char* address,int port,int index,int w,int h);
#else // !FROM_IMAGE
	bool SetUp(char* address,int port,int start);
#endif // FROM_IMAGE
#if FROM_IMAGE
	template<class PixelType,class LabelType,class StorageType>
	void Run(StreamingGrid* pixels,StreamingGrid* labels,StreamingGrid* out);
#if SUPPORT_SHARPEN
	template<class PixelType,class StorageType>
	void Run(StreamingGrid* pixels,StreamingGrid* out);
#endif // SUPPORT_SHARPEN
#else // !FROM_IMAGE
	template<class StorageType>
	void Run(StreamingGrid* laplacian,StreamingGrid* out);
#endif // FROM_IMAGE
	int start	(void)	const;
	int end		(void)	const;
	int width	(void)	const;
	int height	(void)	const;
	int size	(void)	const;
	int cEnd	(void)	const;
	int cHeight	(void)	const;
	int cSize	(void)	const;
	int iters	(void)	const;
	int vCycles	(void)	const;
	int quality	(void)	const;
	bool sharpen(void)	const;
};
#include "SocketedMultigridSolver.inl"
#endif // NEW_CLIENT_SERVER
#endif // SOCKETED_MULTIGRID_SOLVER_INCLUDED