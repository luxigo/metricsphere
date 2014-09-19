#ifndef SOCKETED_STREAMING_SOLVER_INCLUDED
#define SOCKETED_STREAMING_SOLVER_INCLUDED

#define STREAMING_GRID_BUFFER_MULTIPLIER 2
#define ALIGNMENT 16

#define NON_BLOCKING_SOCKETS 1
#define SINGLE_OUTPUT 1
#define SAME_SIZE_BUFFERS 1
#define DEBUG_FLAGS 1
#define TIME_IO 1
#define MISHA_DENORMAL_CONTROL 1
#define MISHA_ADJUST_FLOATS 0
#define MISHA_MIN_FLOAT 1e-20
#define PAD_OUT 1
#define MISHA_FIX 1

#define MyModIndex(idx,mod) (ModIndex(idx,mod))


class DebugFlags
{
public:
	bool noGradient,noDivergence,noResidual,noRestriction,noProlongation,noSolver,noInput,noOutput,noTempIO;
	bool noLeftRight , noTopBottom , emptyPackets;
	DebugFlags(void)
	{
		noGradient = noDivergence = noResidual = noRestriction = noProlongation = noSolver = noInput = noOutput = noTempIO = false;
		noLeftRight = noTopBottom = emptyPackets = false;
	}
};
static DebugFlags debugFlags;

#include <emmintrin.h>
#include "LinearAlgebra/SparseMatrix.h"
#include "LaplacianMatrix/LaplacianMatrix1D.h"
#include "LaplacianMatrix/LaplacianMatrix2D.h"
#include "Util/Socket.h"
#include "Util/GridStream.h"
#include "Util/MultiStreamIO.h"
#include "Util/Array.h"

const int Degree=2;
const int Type=ZERO_DERIVATIVE;
const int RealPerWord=sizeof(__m128)/sizeof(float);
const int WordPerDegree = (Degree+Degree+RealPerWord-1)/RealPerWord;


template< int Channels >
class AverageColor
{
	double values[3];
public:
	double& operator[] ( int idx ){ return values[idx]; }
};

class SocketedStreamData
{
protected:
	DataStream *leftStream , *rightStream;
	SOCKET syncXSocket , syncRSocket;
	int _start128,_end128,_size128,_paddedSize128,_padSize128;	// Sizes in terms of __m128's
	int _start,_end,_size,_paddedSize,_padSize;					// Sizes in terms of float's
#if SUPPORT_CYLINDRICAL
	bool _set( int width , int height , int start , int end , DataStream* leftStream , SOCKET syncSocket , DataStream* rightStream , int periodicType );
#else // !SUPPORT_CYLINDRICAL
	bool _set( int width , int height , int start , int end , DataStream* leftStream , SOCKET syncSocket , DataStream* rightStream , bool spherical );
#endif // SUPPORT_CYLINDRICAL
public:
#if SUPPORT_CYLINDRICAL
	bool SetSocketedStreamData( int width , int height , int start , int end ,int iters , DataStream* leftStream , SOCKET syncSocket , DataStream* rightStream , int periodicType );
	bool SetSocketedStreamData( int width , int height , int start , int end ,DataStream* leftStream , SOCKET syncSocket , DataStream* rightStream , int periodicType );
#else // !SUPPORT_CYLINDRICAL
	bool SetSocketedStreamData( int width , int height , int start , int end ,int iters , DataStream* leftStream , SOCKET syncSocket , DataStream* rightStream , bool spherical);
	bool SetSocketedStreamData( int width , int height , int start , int end ,DataStream* leftStream , SOCKET syncSocket , DataStream* rightStream , bool spherical );
#endif // SUPPORT_CYLINDRICAL
	int start	(void)	const;
	int end		(void)	const;
	int size	(void)	const;
};

class TemplateSSE
{
public:
	__declspec (align(16)) __m128 matrixValues[4][5];
	float diagonal[4];
};
typedef FiniteElements2D<double,Type,Degree>::FullMatrixStencil MatrixStencil;
typedef FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::FullDotProductStencil DotProductStencil;

void SetDownSampledStencil( int dim , int iters , DotProductStencil& outDot , DotProductStencil& outD2Dot );
void SetDownSampledStencil( const DotProductStencil& inDot , const DotProductStencil& inD2Dot, int dim , int iters ,
							DotProductStencil& outDot , DotProductStencil& outD2Dot );

template< int Channels , class SyncType >
class SocketedStreamingSolver : public SocketedStreamData
{
private:
	Pointer( SyncType ) syncBuffer;
	Pointer( Pointer( __m128 ) ) RStream;
	Pointer( Pointer( __m128 ) ) BStream;
	Pointer( Pointer( __m128 ) ) XStream;
	Pointer( __m128 ) localXAccum[Degree];
	Pointer( __m128 ) RBuffer[Degree*Channels];
	Pointer( __m128 ) XBuffer[Degree*Channels];

	void freeStreams(void);

	void SyncSolverLeft ( int idx , bool read , bool overlapped=true );
	void SyncSolverRight( int idx , bool read , bool overlapped=true );
	void SyncSolverHead( int idx , bool read );
	void SyncSolverTail( int idx , bool read );
	void SyncResidualHead( int idx , bool read , bool syncPeriodic );
	void SyncResidualTail( int idx , bool read , bool syncPeriodic );
	void Solve			( int idx , int c , int sSolve , int eSolve );
	void SolveInterior	( int idx , int c , int sSolve , int eSolve );
	void SolveReverse			( int idx , int c , int sSolve , int eSolve );
	void SolveInteriorReverse	( int idx , int c , int sSolve , int eSolve );
	void SetResidual			(int idx,int c,int sSolve,int eSolve);
	void SetInteriorResidual	(int idx,int c,int sSolve,int eSolve);
protected:
	bool _deleteServer;
	MultiStreamIOServer *_server;
#if SUPPORT_CYLINDRICAL
	int periodicType;
	bool clearX , clearB;
#else // !SUPPORT_CYLINDRICAL
	bool clearX , clearB , spherical;
#endif // SUPPORT_CYLINDRICAL
	int xSize,bSize,rSize,iters,index;
	Pointer( TemplateSSE ) lapTemplates;
	Pointer( TemplateSSE ) zeroLapTemplates;


	static int OffsetX(int iters);
	static int OffsetB(int iters);
	static int OffsetR(void);
public:
	int laneNum;
	int progressCount;
	int pastIndex;
	double startProgressTime , pastProgressTime;
	bool showProgress;
	int major,minor;
	float laplacianScale,laplacianScaleR;
	bool setResidual,laplacianRescale;
	double bSquareNorm , rSquareNorm , xSquareNorm , solutionSum[Channels];
#if TIME_IO
	double vSync , hSync , rSync;
#endif // TIME_IO

	SocketedStreamingSolver(void);
	~SocketedStreamingSolver(void);

	void Init( int start , int end , int major , int minor , int iters ,
		DataStream* leftStream , SOCKET syncSocket , DataStream *rightStream ,
#if SUPPORT_CYLINDRICAL
		int periodicType , MultiStreamIOServer* server
#else // !SUPPORT_CYLINDRICAL
		bool spherical , MultiStreamIOServer* server
#endif // SUPPORT_CYLINDRICAL
		);
	void Init( const MatrixStencil& lStencil , int start , int end , int major , int minor , int iters ,
		DataStream* leftStream , SOCKET syncSocket , DataStream *rightStream ,
#if SUPPORT_CYLINDRICAL
		int periodicType , MultiStreamIOServer* server
#else // !SUPPORT_CYLINDRICAL
		bool spherical , MultiStreamIOServer* server
#endif // SUPPORT_CYLINDRICAL
		);

	void Set(int rSize=1);
	void Set(int start,int bStart,int xStart,int bEnd,int xEnd,int rSize=1);
	void UnSet(void);
	Pointer( float ) GetRRow( int row , int channel );
	Pointer( float ) GetXRow( int row , int channel );
	Pointer( float ) GetBRow( int row , int channel );
	bool Increment(void);
	template< class StorageType > bool UpdateXInput  ( StreamingGrid* X );
	template< class StorageType > bool UpdateBInput  ( StreamingGrid* B );
	template< class StorageType > bool UpdateXOutput ( StreamingGrid* X );
	template< class StorageType > bool UpdateBOutput ( StreamingGrid* B );

	void Solve(void);
	void SetResidual(void);
};

template<int Channels>
class SocketedMultiGridRestrictionNode
{
public:
	virtual void InitRestriction			(void){;}													// Generic Initializer
	virtual void SetRestriction				(void){;}													// Sets up interleaved dependencies
	virtual void SetRestriction				( Pointer( float ) lB , int c , int idx , int sRestrict , int eRestrict ){;}	// Set parent's row
	virtual void SolveRestriction			(void){;}
	virtual bool IterateRestriction			(void){return false;}
};

template< int Channels , class StorageType , class SyncType >
class SocketedMultiGridStreamingSolver : public SocketedStreamingSolver< Channels , SyncType > , public SocketedMultiGridRestrictionNode<Channels>
{
	typedef typename FiniteElements1D<double,Type,Degree>::FullProlongationStencil ProlongationStencil;
	typedef typename FiniteElements1D<double,Type,Degree>::FullRestrictionStencil RestrictionStencil;
	class ProlongationStencilSSE
	{
	public:
		__declspec (align(16)) __m128 matrixValues[2][4];
	};
	class ProlongationStencilSSE2
	{
	public:
		__declspec (align(16)) __m128 matrixValues[2];
	};
	Pointer( __m128 ) localRAccum;
	Pointer( ProlongationStencilSSE2 ) prolongationStencil;
	ProlongationStencil majorProlongationStencil,minorProlongationStencil;
	RestrictionStencil majorRestrictionStencil,minorRestrictionStencil;
	int prolongationOffset;
	int startProlongation,startRestriction;

	int restrictionBit;
	StreamingGrid *B , *X;	// Defined to be storage type by default
	bool inCore;
public:
	DotProductStencil dotMajor,d2DotMajor,dotMinor,d2DotMinor;

	// There is some messiness in terms of how in/out are defined and used (e.g. do they increment or do they set?)
	StreamingGrid *inX , *outX , *inB , *outB , *outR;
	Pointer( float ) scratchR;

	SocketedMultiGridRestrictionNode *rChild;
	SocketedMultiGridStreamingSolver *parent,*pChild;

	SocketedMultiGridStreamingSolver(void);
	~SocketedMultiGridStreamingSolver(void);

	void Initialize( int start , int end , int major , int minor , int iters ,
		DataStream* leftStream , SOCKET* syncSockets , DataStream* rightStream ,
#if SUPPORT_CYLINDRICAL
		bool memoryMappedFile , int periodicType , MultiStreamIOServer* server
#else // !SUPPORT_CYLINDRICAL
		bool memoryMappedFile , bool spherical , MultiStreamIOServer* server
#endif // SUPPORT_CYLINDRICAL
		);
	void Initialize(
		DotProductStencil& dotMajor , DotProductStencil& d2DotMajor , DotProductStencil& dotMinor , DotProductStencil& d2DotMinor ,
		int start , int end , int major , int minor , int iters , 
		DataStream* leftStream , SOCKET* syncSockets , DataStream* rightStream ,
#if SUPPORT_CYLINDRICAL
		bool memoryMappedFile , int periodicType , MultiStreamIOServer* server
#else // !SUPPORT_CYLINDRICAL
		bool memoryMappedFile , bool spherical , MultiStreamIOServer* server
#endif // SUPPORT_CYLINDRICAL
		);
	void Initialize(
		double iWeight , int start , int end , int major , int minor , int iters ,
		DataStream* leftStream , SOCKET* syncSockets , DataStream* rightStream ,
#if SUPPORT_CYLINDRICAL
		bool memoryMappedFile , int periodicType , MultiStreamIOServer* server
#else // !SUPPORT_CYLINDRICAL
		bool memoryMappedFile , bool spherical , MultiStreamIOServer* server
#endif // SUPPORT_CYLINDRICAL
		);
	void Initialize(
		DotProductStencil& dotMajor , DotProductStencil& d2DotMajor , DotProductStencil& dotMinor , DotProductStencil& d2DotMinor , double iWeight ,
		int start , int end , int major , int minor , int iters,
		DataStream* leftStream , SOCKET* syncSockets , DataStream* rightStream ,
#if SUPPORT_CYLINDRICAL
		bool memoryMappedFile , int periodicType , MultiStreamIOServer* server
#else // !SUPPORT_CYLINDRICAL
		bool memoryMappedFile , bool spherical , MultiStreamIOServer* server
#endif // SUPPORT_CYLINDRICAL
		);
	void InitProlongation	(void);
	void InitRestriction	(void);

	void SetProlongation	(void);
	void SetRestriction		(void);
	void UnSetProlongation	(void);
	void UnSetRestriction	(void);
	bool IterateProlongation(void);
	bool IterateRestriction	(void);
	void SolveProlongation	(void);
	void SolveRestriction	(void);
	void SetRestriction			( Pointer( float ) lowB , int c , int idx , int sRestrict , int eRestrict );
	void SetInteriorRestriction ( Pointer( float ) lowB , int c , int idx , int sRestrict , int eRestrict );
	void SetProlongation		( Pointer( float ) highX , int c , int highIdx , int highStart , int highSize , int highMajor ,int highMinor , double scale );
	void ProlongationUpdate( int idx , double scale );
};

template<class PixelType,class LabelType>
class ImageData
{
public:
	PixelType pixel;
	LabelType label;
};

template<class LabelType,int Channels>
class LabelData
{
public:
	LabelType l[Channels];
	inline bool operator == (const LabelData& ld) const;
	static LabelData BlackLabel( void );
};


template< int Channels , class PixelType , class LabelType , class StorageType , class SyncType >
class SocketedStreamingDivergence : public SocketedStreamData , public SocketedMultiGridRestrictionNode<Channels>
{
	static const int ISize  = 3*Degree+2;
	static const int DXSize = 3*Degree+1;
	static const int DYSize = 3*Degree;

#if SUPPORT_CYLINDRICAL
	int _periodicType;
#else // !SUPPORT_CYLINDRICAL
	bool _spherical;
#endif // SUPPORT_CYLINDRICAL
	int _iWidth,_iHeight;
	int index;
	Pointer( LabelData< LabelType , Channels > ) labels[ISize];
	Pointer( float ) pixels[ISize*Channels];

	Pointer( __m128 ) dx[DXSize*Channels];
	Pointer( __m128 ) dy[DYSize*Channels];
	Pointer( float ) GetPixelRow(int row,int channel);
	Pointer( LabelData< LabelType , Channels > ) GetLabelRow(int row);
	Pointer( float ) GetDXRow(int row,int channel);
	Pointer( float ) GetDYRow(int row,int channel);
	Pointer( ImageData< SyncType , LabelType > ) syncBuffer;


	void SyncImageLeft ( int idx , bool read );
	void SyncImageRight( int idx , bool read );
	void SyncImageHead( int idx , bool read );
	void SyncImageTail( int idx,  bool read );
	StreamingGrid *pixelStream,*labelStream;

	Pointer( __m128 ) localDMajorAccum;
	Pointer( __m128 ) localDMinorAccum;

	typename FiniteElements1D<float,Type,Degree>::DotProduct<Type,Degree>::FullDotProductStencil dotMajorStencil,dotMinorStencil;
	typename FiniteElements1D<float,DERIVATIVE(Type),Degree-1>::DotProduct<Type,Degree>::FullDotProductStencil dDotMajorStencil,dDotMinorStencil;
	typename FiniteElements2D<float,Type,Degree>::FullDivergenceStencil divergenceStencil;

	void _setPartials( int y );
	void _setPartialsX( int y , int start , int end );
	void _setPartialsY( int y , int start , int end );
	void _setPartialsX(int y);
	void _setPartialsY(int y);
public:
	bool blackOut;
	int major,minor;
	AverageColor< Channels > average , meanB;
#if TIME_IO
	double vSync , hSync;
#endif // TIME_IO

	SocketedMultiGridStreamingSolver< Channels , StorageType , SyncType > *parent;
	SocketedStreamingDivergence(void);
	~SocketedStreamingDivergence(void);

	void Initialize( StreamingGrid *pixels , StreamingGrid* labels , int start , int end , int major , int minor , int iters ,
		DataStream* leftStream , SOCKET* syncSockets , DataStream* rightStream ,
#if SUPPORT_CYLINDRICAL
		bool memoryMappedFile , int periodicType , MultiStreamIOServer* server
#else // !SUPPORT_CYLINDRICAL
		bool memoryMappedFile , bool spherical , MultiStreamIOServer* server
#endif // SUPPORT_CYLINDRICAL
		);

	void InitRestriction		(void);
	void SetRestriction			(void);
	void UnSetRestriction		(void);
	void SetRestriction			( Pointer( float ) lB , int c , int idx , int sRestrict , int eRestrict );
	void SetInteriorRestriction	( Pointer( float ) lB , int c , int idx , int sRestrict , int eRestrict );
	void SolveRestriction		(void);
	bool IterateRestriction		(void);
};


template< int Channels , class PixelType , class StorageType , class SyncType >
class SocketedStreamingLaplacian : public SocketedStreamData , public SocketedMultiGridRestrictionNode<Channels>
{
	static const int ISize  = 3*Degree+2;

#if SUPPORT_CYLINDRICAL
	int _periodicType;
#else // !SUPPORT_CYLINDRICAL
	bool _spherical;
#endif // SUPPORT_CYLINDRICAL
	int _iWidth,_iHeight;
	int index;
	Pointer( float ) GetPixelRow(int row,int channel);
	void SyncImageLeft ( int idx , bool read );
	void SyncImageRight( int idx , bool read );
	void SyncImageLeftRight( int idx );
	void SyncImageHead( int idx , bool read );
	void SyncImageTail( int idx , bool read );
	Pointer( SyncType ) syncBuffer;
	StreamingGrid *pixelStream;

	typename FiniteElements1D<float,Type,Degree>::DotProduct<Type,Degree>::FullDotProductStencil dotMajorStencil,dotMinorStencil;
	typename FiniteElements1D<float,Type,Degree>::DotProduct<Type,Degree>::FullDotProductStencil d2DotMajorStencil,d2DotMinorStencil;

	Pointer( __m128 ) localPAccum[Degree];
	Pointer( __m128 ) PStream[ISize*Channels];
	Pointer( TemplateSSE ) lapTemplates;
public:
	int major,minor;
	AverageColor< Channels > average , meanB;
#if TIME_IO
	double vSync , hSync;
#endif // TIME_IO

	SocketedMultiGridStreamingSolver< Channels , StorageType , SyncType > *parent;
	SocketedStreamingLaplacian(void);
	~SocketedStreamingLaplacian(void);
	void Initialize( StreamingGrid *pixels , double iWeight , double gWeight , int start , int end , int major , int minor , int iters , 
		DataStream* leftStream , SOCKET* syncSockets , DataStream* rightStream ,
#if SUPPORT_CYLINDRICAL
		bool memoryMappedFile , int periodicType , MultiStreamIOServer* server
#else // !SUPPORT_CYLINDRICAL
		bool memoryMappedFile , bool spherical , MultiStreamIOServer* server
#endif // SUPPORT_CYLINDRICAL
		);

	void InitRestriction		(void);
	void SetRestriction			(void);
	void UnSetRestriction		(void);
	void SetRestriction			( Pointer( float ) lB , int c , int idx , int sRestrict , int eRestrict );
	void SetInteriorRestriction	( Pointer( float ) lB , int c , int idx , int sRestrict , int eRestrict );
	void SolveRestriction		(void);
	bool IterateRestriction		(void);
};

#include "StreamingSolver128.inl"
#include "SocketedStreamingSolver.inl"
#endif // SOCKETED_STREAMING_SOLVER_INCLUDED
