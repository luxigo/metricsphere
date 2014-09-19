#ifndef ADAPTIVE_EQUI_RECTANGULAR_INCLUDED
#define ADAPTIVE_EQUI_RECTANGULAR_INCLUDED

#include "EquiRect.h"

template< class Real , class PolyReal=double > class AdaptiveEquiRectangular
{
	int _M , _Samples , _Dim;
	const SphericalStencilTable< double >* _stencilTable;
	std::vector< int > rowDimension;
	static void UpStencil( const double inStencil[5] , double stencil[8] );
	static void DownStencil( const double inStencil[5] , double evenStencil[4] , double oddStencil[4] );
	template< int Dim > static void _UpStencil  ( const double inStencil[Dim] , double stencil[Dim+3] );
	template< int Dim > static void _DownStencil( const double inStencil[Dim] , double evenStencil[(Dim>>1)+2] , double oddStencil[(Dim>>1)+2] );
public:
	AdaptiveEquiRectangular( void ){ _stencilTable = NULL; }
	AdaptiveEquiRectangular( int M , int Samples , const SphericalStencilTable< double >* stencilTable=NULL );
	template< class Real2 > AdaptiveEquiRectangular( const AdaptiveEquiRectangular< Real2 >& sphere ) : AdaptiveEquiRectangular( sphere._M , sphere._Samples , sphere._stencilTable ) { ; }
	AdaptiveEquiRectangular child ( void ) const { return AdaptiveEquiRectangular( _M<<1 , _Samples , _stencilTable ); }
	AdaptiveEquiRectangular parent( void ) const { return AdaptiveEquiRectangular( _M>>1 , _Samples , _stencilTable ); }
	AdaptiveEquiRectangular* Parent( void ) const { if( _M>1 ) return new AdaptiveEquiRectangular( _M>>1 , _Samples , _stencilTable ) ; else return NULL; }
	int height(void) const { return _M<<1; }
	int width (void) const { return _M; }
	double area( void ) const { return 4.0*M_PI; }

	int  dimension( void ) const;
	int  dimension( int idx ) const;
	int  index( int m , int n )	const;
	int  index( int m , int n , int& mm , int & nn ) const;
	void index( int idx , int& m , int& n ) const;
	void index( int m , int n , double& theta , double& phi ) const;


	void elementWeights       ( Vector< Real >& weights ) const;
	void laplacianMatrix      ( SparseMatrix< Real >& lMatrix , double iWeight=0 , double lWeight=1.0 ) const;
	void equiRectangularMatrix( SparseMatrix< Real >& zaMatrix ) const;
};
#include "AdaptiveEquiRect.inl"
#endif // ADAPTIVE_EQUI_RECTANGULAR_INCLUDED