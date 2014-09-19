#include <math.h>
#include "FunctionBasis/PPolynomial.h"
////////////////////////
// SphericalLaplacian //
////////////////////////

#ifndef M_PI
#define M_PI		3.14159265358979323846
#endif

template<class Real>
bool SphericalLaplacian<Real>::LaplacianStencil( int thetaDim , int phiDim , int thetaIndex ,
												 typename FiniteElements2D< Real , ZERO_DERIVATIVE , 2 , ZERO_DERIVATIVE , 2 >::FullMatrixStencil::MatrixStencil& s ,
												 int samples , bool negate )
{
	if(thetaIndex<0 || thetaIndex>=thetaDim || thetaDim<6 || phiDim<=3) return false;

	double dotPhi[5] , d2DotPhi[5] , dotTheta[5] , d2DotTheta[5];
	double thetaShift =     M_PI / thetaDim;
	double phiShift   = 2.0*M_PI /   phiDim;
//	const double EPS = 1e-12;
//	const double EPS = 1e-10;
//	const double EPS = 1e-10;
	const double EPS = 0;

	PPolynomial<2> thetaFunction = PPolynomial<2>::GaussianApproximation().scale( thetaShift );
	PPolynomial<2> phiFunction   = PPolynomial<2>::GaussianApproximation().scale( phiShift   );

	PPolynomial<2> dThetaFunction , dPhiFunction;
	PPolynomial<2> northCapThetaFunction , northCapDThetaFunction;
	PPolynomial<2> southCapThetaFunction , southCapDThetaFunction;

double arcArea = ( cos( thetaIndex*thetaShift ) - cos( (thetaIndex+1)*thetaShift ) ) * 2.0 * PI;
	if( thetaIndex<3 )
	{
		northCapThetaFunction  = thetaFunction.shift( 0.5 * thetaShift ) + thetaFunction.shift( -0.5 * thetaShift );
double fScale = northCapThetaFunction.integralSine( 0 , M_PI ) * 2.0 * M_PI / ( arcArea );
//northCapThetaFunction *= sqrt( fScale );
northCapThetaFunction /= fScale;
		northCapDThetaFunction = northCapThetaFunction.derivative();
	}
	if( thetaIndex>=thetaDim-3 )
	{
		southCapThetaFunction  = thetaFunction.shift( (thetaDim+0.5) * thetaShift ) + thetaFunction.shift( (thetaDim-0.5) * thetaShift );
double fScale = southCapThetaFunction.integralSine( 0 , M_PI ) * 2.0 * M_PI / ( arcArea );
//southCapThetaFunction *= sqrt( fScale );
southCapThetaFunction /= fScale;
		southCapDThetaFunction = southCapThetaFunction.derivative();
	}
	thetaFunction = thetaFunction.shift ( (thetaIndex+0.5) * thetaShift );
double fScale = thetaFunction.integralSine( 0 , M_PI ) * phiFunction.Integral() / ( arcArea / phiDim );
//thetaFunction *= sqrt( fScale );
thetaFunction /= fScale;
	dThetaFunction = thetaFunction.derivative();
	dPhiFunction   = phiFunction.derivative();
if(0)
{
	double arcArea = ( cos( thetaIndex*thetaShift ) - cos( (thetaIndex+1)*thetaShift ) ) * 2.0 * PI;
	if		( thetaIndex==0 )			printf( "%2d] %f\n" , thetaIndex , northCapThetaFunction.integralSine( 0 , M_PI ) * 2.0 * M_PI / ( arcArea ) );
	else if	( thetaIndex==thetaDim-1 )	printf( "%2d] %f\n" , thetaIndex , southCapThetaFunction.integralSine( 0 , M_PI ) * 2.0 * M_PI / ( arcArea ) );
	else printf( "%2d] %f\n" , thetaIndex , thetaFunction.integralSine( 0 , M_PI ) * phiFunction.Integral() / ( arcArea / phiDim) );
}

	for( int i=0 ; i<5 ; i++ )
	{
		dotTheta[i]   =   (  thetaFunction *  thetaFunction.shift( (i-2)*thetaShift ) ).integralCosecant	( 0+EPS , M_PI-EPS , samples);
		d2DotTheta[i] = - ( dThetaFunction * dThetaFunction.shift( (i-2)*thetaShift ) ).integralSine		( 0 , M_PI );
		dotPhi[i]     =   (    phiFunction *    phiFunction.shift( (i-2)*phiShift   ) ).Integral			( );
		d2DotPhi[i]   = - (   dPhiFunction *   dPhiFunction.shift( (i-2)*phiShift   ) ).Integral			( );
		for( int j=0 ; j<5 ; j++ ) s.values[i][j] = 0;
	}
	if ( thetaIndex == 0 )
	{
		s.values[2][2] = -( northCapDThetaFunction * northCapDThetaFunction              ).integralSine( 0 , M_PI ) * 2.0*M_PI;
		s.values[3][2] = -( northCapDThetaFunction * dThetaFunction.shift( 1*thetaShift) ).integralSine( 0 , M_PI ) * phiFunction.Integral();
		s.values[4][2] = -( northCapDThetaFunction * dThetaFunction.shift( 2*thetaShift) ).integralSine( 0 , M_PI ) * phiFunction.Integral();
	}
	else if	( thetaIndex == thetaDim-1 )
	{
		s.values[0][2] = -( southCapDThetaFunction * dThetaFunction.shift(-2*thetaShift) ).integralSine( 0 , M_PI ) * phiFunction.Integral();
		s.values[1][2] = -( southCapDThetaFunction * dThetaFunction.shift(-1*thetaShift) ).integralSine( 0 , M_PI ) * phiFunction.Integral();
		s.values[2][2] = -( southCapDThetaFunction * southCapDThetaFunction              ).integralSine( 0 , M_PI ) * 2.0*M_PI;
	}
	else if	( thetaIndex == 1 )
	{
		s.values[1][2] = -( northCapDThetaFunction * dThetaFunction ).integralSine( 0 , M_PI ) * phiFunction.Integral();
		for( int i=2 ; i<5 ; i++ )
			for( int j=0 ; j<5 ; j++ )
				s.values[i][j] = dotTheta[i]*d2DotPhi[j] + d2DotTheta[i]*dotPhi[j];
	}
	else if	( thetaIndex == thetaDim-2 )
	{
		for( int i=0 ; i<3 ; i++ )
			for( int j=0 ; j<5 ; j++ )
				s.values[i][j] = dotTheta[i]*d2DotPhi[j] + d2DotTheta[i]*dotPhi[j];

		s.values[3][2] = -( southCapDThetaFunction * dThetaFunction ).integralSine( 0 , M_PI ) * phiFunction.Integral();
	}
	else if	( thetaIndex == 2 )
	{
		s.values[0][2] = -( northCapDThetaFunction * dThetaFunction ).integralSine( 0 , M_PI ) * phiFunction.Integral();
		for( int i=1 ; i<5 ; i++ )
			for( int j=0 ; j<5 ; j++ )
				s.values[i][j] = dotTheta[i]*d2DotPhi[j] + d2DotTheta[i]*dotPhi[j];
	}
	else if	( thetaIndex == thetaDim-3 )
	{
		for( int i=0 ; i<4 ; i++ )
			for( int j=0 ; j<5 ; j++ )
				s.values[i][j] = dotTheta[i]*d2DotPhi[j] + d2DotTheta[i]*dotPhi[j];
		s.values[4][2] = -( southCapDThetaFunction * dThetaFunction ).integralSine( 0 , M_PI ) * phiFunction.Integral();
	}
	else
		for(int i=0;i<5;i++)
			for(int j=0;j<5;j++)
				s.values[i][j] = dotTheta[i]*d2DotPhi[j] + d2DotTheta[i]*dotPhi[j];

	if( negate )
		for( int i=0 ; i<5 ; i++ )
			for( int j=0 ; j<5 ; j++ )
				s.values[i][j] = -s.values[i][j];
	return true;
}
