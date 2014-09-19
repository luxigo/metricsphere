#ifndef SPHERICAL_LAPLACIAN_INCLUDED
#define SPHERICAL_LAPLACIAN_INCLUDED
#include "LaplacianMatrix/LaplacianMatrix2D.h"

////////////////////////
// SphericalLaplacian //
////////////////////////
template<class Real>
class SphericalLaplacian
{
public:
	static bool LaplacianStencil
		(
		int M , int N , int m ,
		typename FiniteElements2D< Real , ZERO_DERIVATIVE , 2 >::FullMatrixStencil::MatrixStencil& s ,
		int samples , bool negate=false
		);

	/*
	class FullDivergenceStencil
	{
	public:
		class DivergenceStencil
		{
		public:
			Real values1[2*Degree1][2*Degree2+1];
			Real values2[2*Degree1+1][2*Degree2];
		};
		DivergenceStencil caseTable[2*Degree1+1][2*Degree2+2];
	};
	static bool DivergenceStencil(int dim1,int dim2,FullDivergenceStencil& s);

	class FullMatrixStencil
	{
	public:
		class MatrixStencil
		{
		public:
			Real values[2*Degree1+1][2*Degree2+1];
		};
		MatrixStencil caseTable[2*Degree1+1][2*Degree2+1];
	};
	static bool LaplacianStencil(int dim1,int dim2,FullMatrixStencil& s,bool weakForm,bool negate=false);

	class FullProlongationStencil
	{
	public:
		class ProlongationStencil
		{
		public:
			Real values[Degree1+2][Degree2+2];
		};
		ProlongationStencil caseTable[2*Degree1+1][2*Degree2+1];
	};
	static bool ProlongationStencil(int lowD1,int lowD2,FullProlongationStencil &s,int& highD1,int& highD2);
#if NEW_LAPLACIAN_CODE
	class FullRestrictionStencil
	{
	public:
		class RestrictionStencil
		{
		public:
			Real values[(Degree1+3)>>1][(Degree2+3)>>1];
		};
		RestrictionStencil caseTable[2*Degree1+2][2*Degree2+2];
	};
	static bool RestrictionStencil(int highD1,int highD2,FullRestrictionStencil &s,int& lowD1,int& lowD2);
#endif // NEW_LAPLACIAN_CODE
	*/
};
#include "LaplacianMatrix/SphericalLaplacian.inl"
#endif // SPHERICAL_LAPLACIAN_INCLUDED