#ifndef EQUI_RECTANGULAR_INCLUDED
#define EQUI_RECTANGULAR_INCLUDED

#include "LaplacianMatrix/LaplacianMatrix2D.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif // M_PI

template< class Real >
class ThetaStencilTable
{
	struct StencilValues
	{
		Real            thetaWeight   ;
		Real            thetaValues[3];	// For dot-product stencil
		Real divergenceDThetaValues[4];	// For finite differences divergence
		Real   laplacianThetaValues[5];	// For laplacian stencil
		Real laplacianD2ThetaValues[3];	// For laplacian stencil
	};
	std::vector< StencilValues* > _valueTable;
	int _logHeight;
	void _clearTable( void );
public:
	ThetaStencilTable( void );
	~ThetaStencilTable( void );
	bool Init( int height , int samples , bool highPrecision , bool useRestriction );
	Real               thetaWeight (                                   int m , int M ) const;
	bool            setThetaStencil( Real            thetaStencil[5] , int m , int M ) const;
	bool setDivergenceDThetaStencil( Real divergenceDThetaStencil[4] , int m , int M ) const;
	bool   setLaplacianThetaStencil( Real   laplacianThetaStencil[5] , int m , int M , int samples=0 ) const;
	bool setLaplacianD2ThetaStencil( Real laplacianD2ThetaStencil[5] , int m , int M ) const;
	bool read ( FILE* fp );
	bool write( FILE* fp ) const;
};
template< class Real >
class PhiStencilTable
{
	struct StencilValues
	{
		Real            phiWeight   ;
		Real            phiValues[3];	// For dot-product stencil
		Real divergenceDPhiValues[4];	// For finite differences divergence
		Real laplacianD2PhiValues[3];	// For laplacian stencil
	};
	std::vector< StencilValues > _valueTable;
	int _logWidth;

	void _clearTable( void );
public:
	PhiStencilTable( void );
	~PhiStencilTable( void );
	bool Init( int width , bool highPrecision , bool useRestriction );
	Real               phiWeight (                                 int N ) const;
	bool            setPhiStencil( Real            phiStencil[5] , int N ) const;
	bool setDivergenceDPhiStencil( Real divergenceDPhiStencil[4] , int N ) const;
	bool setLaplacianD2PhiStencil( Real laplacianD2PhiStencil[5] , int N ) const;
	bool read ( FILE* fp );
	bool write( FILE* fp ) const;
};
template< class Real >
class SphericalStencilTable : public PhiStencilTable< Real > , public ThetaStencilTable< Real >
{
public:
	bool Init( int M , int N , int samples , bool highPrecision , bool useRestriction )
	{
		return PhiStencilTable::Init( N , highPrecision , useRestriction ) & ThetaStencilTable::Init( M , samples , highPrecision , useRestriction );
	}
	bool read ( const char* fileName )
	{
		FILE* fp = fopen( fileName , "rb" );
		if( !fp ) return false;
		bool ret = PhiStencilTable::read( fp ) && ThetaStencilTable::read( fp );
		fclose( fp );
		return ret;
	};
	bool write( const char* fileName ) const
	{
		FILE* fp = fopen( fileName , "wb" );
		if( !fp ) return false;
		bool ret = PhiStencilTable::write( fp ) && ThetaStencilTable::write( fp );
		fclose( fp );
		return ret;
	}
};

template< class Real , class PolyReal=double > class EquiRectangular
{
	int M , N , Samples;
	const SphericalStencilTable< double >* stencilTable;
	static bool EvaluationStencil( int M , int N , int m , double stencil[3][3] );
	static bool PrimalEvaluationStencil( int M , int N , int m , double stencil[2][2] );
	static bool PCIntegralStencil( int M , int N , int m , double stencil[3][3] );
public:
	static void SetBaseFunctions( int m , int n , int M , int N , PPolynomial< 2 , PolyReal >& thetaBase , PPolynomial< 2 , PolyReal >& phiBase , bool rescale=true );
	static void SetThetaBaseFunction( int m , int M , PPolynomial< 2 , PolyReal >& thetaBase , bool rescale=true );
	static void SetPhiBaseFunction( int n , int N , PPolynomial< 2 , PolyReal >& phiBase , bool rescale=true );
	static void SetDThetaBaseFunction	( int m , int M , PPolynomial< 1 , PolyReal >& thetaBase , bool rescale=true );
	static void SetDPhiBaseFunction	( int n , int N , PPolynomial< 1 , PolyReal >& phiBase , bool rescale=true );

	static double ThetaWeight( int m , int M );
	static double PhiWeight( int N );
	static void            ThetaStencil( int m , int M , double stencil[5] );
	static void              PhiStencil(         int N , double stencil[5] );
	static bool DivergenceDThetaStencil( int m , int M , double stencil[4] );
	static bool DivergenceDPhiStencil  (         int N , double stencil[4] );
	static void   LaplacianThetaStencil( int m , int M , double stencil[5] , int samples );
	static void LaplacianD2ThetaStencil( int m , int M , double stencil[5] );
	static void LaplacianD2PhiStencil  (         int N , double stencil[5] );
	static void            HalfThetaStencil( int m , int M , double stencil[3] );
	static void   HalfLaplacianThetaStencil( int m , int M , double stencil[3] , int samples );
	static void HalfLaplacianD2ThetaStencil( int m , int M , double stencil[3] );
	static void            HalfPhiStencil( int N , double stencil[3] );
	static void HalfLaplacianD2PhiStencil( int N , double stencil[3] );

	static void ThetaEvaluation( int m , int M , double stencil[3] );
	static void       PhiEvaluation( double stencil[3] );
	static void CoarsePhiEvaluation( double evenStencil[3] , double oddStencil[3] );
	static void   FinePhiEvaluation( double stencil[2] );

	static void PrimalThetaEvaluation( int m , int M , double stencil[2] );
	static void       PrimalPhiEvaluation( double stencil[2] );
	static void PrimalCoarsePhiEvaluation( double evenStencil[2] , double oddStencil[3] );
	static void   PrimalFinePhiEvaluation( double stencil[2] );

	static void ThetaPCIntegral( int m , int M , double stencil[3] );
	static void       PhiPCIntegral( double stencil[3] );
	static void CoarsePhiPCIntegral( double evenStencil[2] , double oddStencil[2] );
	static void   FinePhiPCIntegral( double stencil[6] );
	static void AnalyticLaplacianThetaStencil( int m , int M , double stencil[5] );

	static bool LaplacianStencil ( int M , int N , int m , double stencil[5][5] , int samples , const SphericalStencilTable< double >* stencilTable  , bool negate=false );
	static bool DotProductStencil( int M , int N , int m , double stencil[5][5] , const SphericalStencilTable< double >* stencilTable );
	static bool DivergenceDThetaStencil( int M , int N , int m , double stencil[4][5] , int samples );
	static bool DivergenceDPhiStencil  ( int M , int N , int m , double stencil[5][4] , int samples );

	EquiRectangular( int M         , int Samples , const SphericalStencilTable< double >* stencilTable=NULL ) { this->M = M , this->N = 2*M; this->Samples = Samples , this->stencilTable = stencilTable; }
	EquiRectangular( int M , int N , int Samples , const SphericalStencilTable< double >* stencilTable=NULL ) { this->M = M , this->N =   N; this->Samples = Samples , this->stencilTable = stencilTable; }
	EquiRectangular child ( void ) const { return EquiRectangular( M<<1 , N<<1 , Samples , stencilTable ); }
	EquiRectangular parent( void ) const { return EquiRectangular( M>>1 , N>>1 , Samples , stencilTable ); }
	EquiRectangular* Parent( void ) const { if( M>1 && N>1 ) return new EquiRectangular( M>>1 , N>>1 , Samples , stencilTable ) ; else return NULL; }
	int width (void) const { return M; }
	int height(void) const { return N; }
	double area( void ) const { return 4.0*M_PI; }

	static int  Dimension( int M , int N );
	static int  Index( int m , int n , int M , int N );
	static void Index( int idx , int M , int N , int& m , int& n );
	static int  Index( int m , int n , int M , int N , int& mm , int & nn );
	static void Index( int m , int n , int M , int N , double& theta , double& phi );
	int  dimension( void ) const									{ return Dimension( M , N ); }
	int  dimension( int idx ) const									{ return width( ); }
	int  index( int m , int n )	const								{ return Index( m , n , M , N ); }
	int  index( int m , int n , int& mm , int & nn ) const			{ return Index( m , n , M , N , mm , nn ); }
	void index( int idx , int& m , int& n ) const					{ Index( idx , M , N , m , n ); }
	void index( int m , int n , double& theta , double& phi ) const	{ Index( m , n , M , N , theta , phi ); }

	static PolyReal BaseWeight 					( int M , int N , int m , const SphericalStencilTable< double >* stencilTable );
	static PolyReal BaseWeight 					( int M , int N , int m ) { return BaseWeight( M , N , m  , stencilTable ); }
	static void BaseWeights						( int M , int N  , const SphericalStencilTable< double >* stencilTable , Vector< Real >& weights );
	static void BaseWeights						( int M , int N , Vector< Real >& weights ) { BaseWeights( M , N , NULL , weights ); }
	static bool LaplacianMatrix					( int M , int N , int samples , const SphericalStencilTable< double >* stencilTable , SparseMatrix< Real >& laplacianMatrix , double iWeight=0 , double lWeight=1.0 );
	static bool LaplacianMatrix					( int M , int N , int samples , SparseMatrix< Real >& laplacianMatrix , double iWeight=0 , double lWeight=1.0 ){ return LaplacianMatrix( M , N , samples , NULL , laplacianMatrix , iWeight , lWeight ); }
	void elementWeights		( Vector< Real >& weights ) const { BaseWeights( M , N , stencilTable , weights ); }
	void laplacianMatrix	( SparseMatrix< Real >& lMatrix , double iWeight=0 , double lWeight=1.0 ) const	{ LaplacianMatrix( M , N , Samples , stencilTable , lMatrix , iWeight , lWeight ); }
};

#include "EquiRect.inl"
#endif // EQUI_RECTANGULAR_INCLUDED