#ifndef STREAMING_SOLVER_INCLUDED
#define STREAMING_SOLVER_INCLUDED

#define STREAMING_GRID_BUFFER_MULTIPLIER 2
#define ALIGNMENT 16

#define NEW_CODE 1
#define MISHA_FIX 1
#define HIGH_PRECISION 0
#define USE_INTERIOR 1
#define NEW_SOLVE 0		// Enabled, the solver is slightly more accurate but also a little bit slower

#include <emmintrin.h>
#include "LinearAlgebra/SparseMatrix.h"
#include "LaplacianMatrix/LaplacianMatrix1D.h"
#include "LaplacianMatrix/LaplacianMatrix2D.h"
#include "Util/GridStream.h"
#include "Util/MultiStreamIO.h"
#include "Util/Array.h"

static const int Degree = 2;
static const int Type = ZERO_DERIVATIVE;

template< class Real >
class RealType
{
public:
//	typedef __m128 WordClass;
	static const int RealPerWord = sizeof( __m128 ) / sizeof( Real );
};

// This class is responsible for performin Gauss-Seidel relaxations 
// within a single level
template< class Real >
class SphericalStreamingSolver : public RealType< Real >
{
public:
	class PhiStencilData
	{
	public:
		double              phiStencil[5] ,         laplacianD2PhiStencil[5] ,                   phiWeight;
		double            upPhiStencil[8] ,            downEvenPhiStencil[4] ,            downOddPhiStencil[4];
		double upLaplacianD2PhiStencil[8] , downEvenLaplacianD2PhiStencil[4] , downOddLaplacianD2PhiStencil[4];

		void set( int width , const SphericalStencilTable< double >* stencilTable );
	};
	class ThetaStencilData
	{
	public:
		double thetaStencil[5] , laplacianThetaStencil[5] , laplacianD2ThetaStencil[5];
		void set( int i , int N , int samples , const SphericalStencilTable< double >* stencilTable );
	};
	class SphericalLaplacianStencil
	{
	public:
		Real evenMatrixValues[5][8] , oddMatrixValues[5][8];
		double diagonal;
		bool isRegular;
	};
	class SphericalLaplacianStencilSSE
	{
	public:
		__declspec ( align( ALIGNMENT ) ) __m128 matrixValues[4][5];
		Real diagonal[4];
	};
private:
	AdaptiveEquiRectangular< Real > _sphere;
	Pointer( Real ) _scratch0;
	Pointer( Real ) _scratch1;
	Pointer( Real ) _scratch2;
	Pointer( Real ) scratch0;
	Pointer( Real ) scratch1;
	Pointer( Real ) scratch2;
	PhiStencilData* _phiStencilData;
	SphericalLaplacianStencil *laplacianStencils;
	ThetaStencilData* thetaStencils;
	void *_laplacianStencilsSSE , *_zeroCenteredLaplacianStencilsSSE;
	SphericalLaplacianStencilSSE *laplacianStencilsSSE , *zeroCenteredLaplacianStencilsSSE;
	void StencilToStencilSSE( const SphericalLaplacianStencil& inStencil , SphericalLaplacianStencilSSE& outStencil , bool zeroCenter );
	void SetThetaStencil( int row );
	void SetStencil( SphericalLaplacianStencil& stencil , int row , double dWeight , double lWeight );
	void ScaleBRow( int idx , Real scale );
	int _channels;
public:
	int index;
protected:
	int xSize , bSize , rSize , iters , _iters;
	double dWeight , lWeight;
	Pointer( __m128 ) localXPtrs[2*Degree+1];

	Pointer( Real ) localBPtr;
	Pointer( Real ) localRPtr;
	int _columns;
	const SphericalStencilTable< double >* _stencilTable;
public:
	Pointer( __m128 ) XStream;
	Pointer( __m128 ) BStream;
	Pointer( __m128 ) RStream;
	Pointer( Real ) GetXRow( int row , int channel );
	Pointer( Real ) GetBRow( int row , int channel );
	Pointer( Real ) GetRRow( int row , int channel );
	ThetaStencilData&             GetThetaStencil                   ( int row );
	SphericalLaplacianStencil&    GetLaplacianStencil               ( int row );
	SphericalLaplacianStencilSSE& GetLaplacianStencilSSE            ( int row );
	SphericalLaplacianStencilSSE& GetZeroCenteredLaplacianStencilSSE( int row );
	int RowWidth( int row ) const;
	int rows , columns , samples;

	bool setResidual;
	double *bSquareNorm , *rSquareNorm , *xSquareNorm;
	SphericalStreamingSolver ( void );
	~SphericalStreamingSolver( void );

	void Init( int channels , int columns , int rows , int samples , const SphericalStencilTable< double >* stencilTable , double dWeight , double lWeight );
	void SetIterations( int iters , int rSize=1 );
	void SetIterations( int start , int iters , int bStart , int xStart , int bEnd , int xEnd , int rSize=1 );
	void UnSetIterations( void );

	bool Increment( void );
	template< class StorageType > bool UpdateXInput ( StreamingGrid* X , bool overwrite );
	template< class StorageType > bool UpdateBInput ( StreamingGrid* B , bool overwrite );
	template< class StorageType > bool UpdateXOutput( StreamingGrid* X );
	template< class StorageType > bool UpdateBOutput( StreamingGrid* B );

	void Solve( void );

	void Solve   ( int idx , bool reverse = false );
	void SolveSSE( int idx , bool reverse = false );
	void SetResidual   ( int idx );
	void SetResidualSSE( int idx );
	static MultiStreamIOServer server;

	template< class StorageType > void Fullsolve( StreamingGrid* inX , StreamingGrid* inB , StreamingGrid* outX );
	int Channels( void ) const;
};

// This abstract class is responsible for transitioning between levels
template< class Real >
class MultiGridRestrictionNode
{
protected:
	bool _showProgress;
public:
	MultiGridRestrictionNode *parent , *child;
	MultiGridRestrictionNode                ( void ){ parent = child = NULL; _showProgress = false; }
	virtual void InitRestriction			( void ){;}								// Generic Initializer
	virtual void InitProlongation			( void ){;}								// Generic Initializer
	void InitRestriction			( bool showProgress ){ _showProgress=showProgress , InitRestriction ( ); }
	void InitProlongation			( bool showProgress ){ _showProgress=showProgress , InitProlongation( ); }
	virtual void SetRestrictionIterations	( int iters ){;}						// Sets up interleaved dependencies
	virtual void SetProlongationIterations	( int iters ){;}						// Sets up interleaved dependencies
	virtual void SetRestriction				( Pointer( Real )* lB , int idx , int lWidth ){;}	// Set parent's row
	virtual void SolveRestriction			( void ){;}
	virtual void SolveProlongation			( void ){;}
	virtual bool IterateRestriction			( void ){return false;}
	virtual bool IterateProlongation		( void ){return false;}
	virtual int  RowWidth					( int idx ) const { return 0; }
	virtual Pointer( Real ) GetXRow( int idx , int c ) { return NullPointer< Real >( ); }
};

template< class Real , class StorageType=Real >
class MultiGridSphericalStreamingSolver : public SphericalStreamingSolver< Real > , public MultiGridRestrictionNode< Real >
{
	typedef typename FiniteElements1D< double , Type , Degree >::FullRestrictionStencil RestrictionStencil;

	int prolongationOffset;
	int startProlongation , startRestriction;

	int restrictionBit;
	StreamingGrid *B , *X;
	bool inCore;

	bool _clearAverage;
	double* _average , _weightSum;
public:
// There is some messiness in terms of how in/out are defined and used (e.g. do they increment or do they set?)
	StreamingGrid *inX , *outX , *inB , *outB;

	MultiGridSphericalStreamingSolver( void );
	~MultiGridSphericalStreamingSolver( void );
	void Initialize( int channels , int major , int minor , int samples , const SphericalStencilTable< double >* stencilTable , double dWeight , double lWeight , bool memoryMappedFile );
	void InitProlongation	( void );
	void InitRestriction	( void );

	void SetProlongationIterations	( int iters );
	void SetRestrictionIterations	( int iters );
#if MISHA_FIX
	void SetRestrictionIterations	( int iters , bool offSet );
#endif // MISHA_FIX
	void UnSetProlongationIterations( void );
	void UnSetRestrictionIterations	( void );
	bool IterateProlongation( void );
	bool IterateRestriction	( void );
	void SolveProlongation	( void );
	void SolveRestriction	( void );
	Pointer( __m128 ) localRPtrs[Degree+2];
	void SetRestriction( Pointer( Real )* lB , int idx , int lWidth );
	void RestrictionUpdate( int idx );
	/*
	void *_localRAccum;
	__m128* localRAccum;
	*/
//	void SetInteriorRestriction( float* lowB , int idx , int major2 , int minor2 );

	void ProlongationUpdate( int idx );

	static bool IsDownSamplable( const int& hMajor , const int& hMinor );
	static bool IsDownSamplable( const int& hMajor , const int& hMinor , int& lMajor , int& lMinor );

	int  RowWidth( int idx ) const;
	Pointer( Real ) GetXRow( int idx , int c);
};

// This class takes equi-rectangular laplacian constraints and turns them into
// adaptive equi-rectangular constraints that can be fed to the solver and also
// takes the adaptive equi-rectangular constraints and turns them back into 
// equi-rectangular constraints
template< class Real , int Channels >
class StreamingAdaptiveNonAdaptiveConverter : public MultiGridRestrictionNode< Real >
{
	double squareNorm;
	AdaptiveEquiRectangular< Real > _sphere;
	std::vector< int > start;
	std::vector< std::vector< double > > pStencils;
	int columns , rows , index;
	Pointer( Real ) dataRow;
	Pointer( Real ) GetRow( int row , int channel );

public:
	StreamingGrid *outX;

	StreamingAdaptiveNonAdaptiveConverter( void );
	~StreamingAdaptiveNonAdaptiveConverter( void );
	void Initialize( int rows );
	void InitRestriction ( void );
	void InitProlongation( void );
	void SetRestrictionIterations ( int iters );
	void SetProlongationIterations( int iters );
	void SetRestriction( Pointer( Real )* lB , int idx , int lWidth );
	void ProlongationUpdate( int idx );
	void SolveRestriction ( void );
	void SolveProlongation( void );
	bool IterateRestriction ( void );
	bool IterateProlongation( void );
	int  RowWidth( int idx ) const;
};

template< class Real , int PixelChannels , int LabelChannels , class PixelType , class LabelType >
class StreamingNonAdaptiveLaplacian : public MultiGridRestrictionNode< Real > , RealType< Real >
{
	Pointer( PixelType ) _previousPixelRow;
	Pointer( PixelType ) _pixelRow;
	Pointer( LabelType ) _previousLabelRow;
	Pointer( LabelType ) _labelRow;
	Pointer( Real ) _dx;
	Pointer( Real ) _dy;
	Pointer( Real ) _values;
	Pointer( Real ) localDMajor;
	Pointer( Real ) localDMinor;
	Pointer( __m128 ) localDMajorAccum;
	Pointer( __m128 ) localDMinorAccum;
	Pointer( __m128 ) localValueAccum;

	int _currentRow;
	int startRestriction;
	int major , minor , _major , index , dSize , size , samples;

	const LabelType* _unknownIndex;
	bool _readNextRow( void );
	double _vWeight , _gScale;
	double phiWeight , phiStencil[5] , divergenceDPhiStencil[5];
	const SphericalStencilTable< double >* _stencilTable;
public:
	double average[PixelChannels];
	double dXSquareNorm , dYSquareNorm , outputSquareNorm;

	StreamingGrid *pixels , *labels;

	StreamingNonAdaptiveLaplacian( void );
	~StreamingNonAdaptiveLaplacian( void );
	void Initialize( int major , int minor , int samples , const SphericalStencilTable< double >* stencilTable , const LabelType* unknownIndex , double vWeight = 0. , double gWeight = 1. );
	void InitRestriction( void );
	void SetRestrictionIterations( int iters );
	void SetRestriction( Pointer( Real )* lB , int idx , int major2 );
	void SetInteriorRestriction( Pointer( float  )* lB , int idx , int major2 );
	void SetInteriorRestriction( Pointer( double )* lB , int idx , int major2 );
	void SolveRestriction( void );
	bool IterateRestriction ( void );
};
template< class Real , int Channels >
class StreamingNonAdaptiveRestrictionGrabber : public MultiGridRestrictionNode< Real > , RealType< Real >
{
	Pointer( Real ) rows[ Channels ];
public:
	int width , height , index;
	StreamingGrid* out;
	void Initialize( int major , int minor ) { width = major , height = minor , index = 0; }
	void InitRestriction( void ){ if( parent ) parent->InitRestriction(); }
	bool IterateRestriction( void )
	{
		if( _showProgress )
			if     ( index <height ) printf("[%.1f%%]  [%d/%d]         \r" , float(index)/float(height-1)*100 , index , height );
			else if( index==height )printf( "                                              \r" );
		if( index<height )
		{
			Pointer( Real ) row = ( Pointer( Real ) )(*out)[index];
			memset( row , 0 , sizeof( Real ) * Channels * width );
			for( int c=0 ; c<Channels ; c++ ) rows[c] = row + c*width;
			child->SetRestriction( rows , index , width );
			out->advance();
			index++;
			return true;
		}
		return false;
	}
	void SolveRestriction( void )
	{
		while( IterateRestriction( ) ){ ; }
	}
};

#include "StreamingSolver128.inl"
#include "SphericalStreamingSolver.inl"
#endif // STREAMING_SOLVER_INCLUDED
