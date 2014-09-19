#ifndef STREAMING_SOLVER_INCLUDED
#define STREAMING_SOLVER_INCLUDED

#define TEST_COMPRESSION 0
#define NEW_MISHA_CODE 0
// BADNESS!!! There appears to be some kind of badness when running on really small images. Is it possible that
// the multistreamioclient never starts writing?
#define REFLECT_GRADIENT 0
#define ACCUMULATE_PROLONGATION 0	// Maybe someday...
#define MERGE_INTERIOR_STREAMS 0	// Maybe someday...

#define SEPARTE_TEMP_IO_SERVER 0
#define NEW_STREAMING_CODE 1
#define STREAMING_GRID_BUFFER_MULTIPLIER 2
#define MOLLIFY_DC 0

#define USE_SSE_CODE 1
#define STORE_RESIDUAL 0

#define USE_SERVER 1

#if USE_SSE_CODE
#include <emmintrin.h>
#endif // USE_SSE_CODE
#include "LinearAlgebra/SparseMatrix.h"
#include "LaplacianMatrix/LaplacianMatrix1D.h"
#include "LaplacianMatrix/LaplacianMatrix2D.h"
#include "Util/GridStream.h"
#include "Util/MultiStreamIO.h"
#if NEW_MISHA_CODE
#include "misha.h"
#endif // NEW_MISHA_CODE

template<class Real>
class RealType
{
public:
#if USE_SSE_CODE
	typedef __m128 WordClass;
#define ALIGNMENT 16
#else // !USE_SSE_CODE
	typedef Real WordClass;
#define ALIGNMENT 1
#endif // USE_SSE_CODE
	static const int RealPerWord=sizeof(WordClass)/sizeof(Real);
};


template<class Real,int Type,int Degree,int Channels>
class StreamingSolver : public RealType<Real>
{
	WordClass *_localR;
#if MERGE_INTERIOR_STREAMS
	WordClass *_IStream;
#else // !MERGE_INTERIOR_STREAMS
	WordClass *_XStream,*_BStream;
#endif // MERGE_INTERIOR_STREAMS
protected:
	bool symmetric,clearX,clearB;
#if MERGE_INTERIOR_STREAMS
	int iSize,rSize,iters,index,_iters;
#else // !MERGE_INTERIOR_STREAMS
	int xSize,bSize,rSize,iters,index,_iters;
#endif // MERGE_INTERIOR_STREAMS
	WordClass *localXPtrs[2*Degree+1];
	WordClass *_localXAccum[Degree],*localXAccum[Degree];
	Real* localBPtr;
	Real* localRPtr;
	int _major;
public:
#if MERGE_INTERIOR_STREAMS
	Real *IStream,*localR;
#else // !MERGE_INTERIOR_STREAMS
	Real *XStream,*BStream,*localR;
#endif // MERGE_INTERIOR_STREAMS
	Real* GetXRow(int row,int channel);
	Real* GetBRow(int row,int channel);
	int major,minor;
protected:


	// Note: For appropriate indexing, the laplacian stencil needs to be indexed (minor,major)
#if USE_SSE_CODE
#if !NEW_MISHA_CODE
public:
	class LaplacianTemplateSSE
	{
	public:
		__declspec (align(ALIGNMENT)) WordClass matrixValues[4][5];
		float diagonal[4];
	};
#endif // !NEW_MISHA_CODE
protected:
	LaplacianTemplateSSE *lapTemplates2,*_lapTemplates2;
	LaplacianTemplateSSE *lapTemplates3,*_lapTemplates3;	// Like lapTemplates2 but with the center term zeroed out
#else // !USE_SSE_CODE
	typename FiniteElements2D<double,Type,Degree>::FullMatrixStencil laplacianStencil;
#endif // USE_SSE_CODE

#if USE_SSE_CODE
	// Strictly Degree=2
#else // !USE_SSE_CODE
	Real GetLaplacianValue(const double lapValues[][2*Degree+1],int j);
	Real GaussSeidelUpdate(const double lapValues[][2*Degree+1],int j,int iB);
	Real GetInteriorLaplacianValue(const double lapValues[][2*Degree+1],int j);
	Real InteriorGaussSeidelUpdate(const double lapValues[][2*Degree+1],int j,int iB);
#endif // USE_SSE_CODE
public:
	Real laplacianScale,laplacianScaleR;
	bool setResidual,laplacianRescale;
	double bSquareNorm,rSquareNorm,xSquareNorm;
	StreamingSolver(void);
	~StreamingSolver(void);

	void Init(const typename FiniteElements2D<double,Type,Degree>::FullMatrixStencil& lStencil,int major,int minor,bool symmetric);
	void Init(int major,int minor,bool symmetric);
	void SetIterations(int iters,int rSize=1);
	void SetIterations(int start,int iters,int bStart,int xStart,int bEnd,int xEnd,int rSize=1);
	void UnSetIterations(void);

	bool Increment(void);
#if MERGE_INTERIOR_STREAMS
	bool UpdateIInput			(StreamingGrid* I);
	bool UpdateIOutput			(StreamingGrid* I);
#endif // MERGE_INTERIOR_STREAMS
	template<class StorageType>
	bool UpdateXInput			(StreamingGrid* X);
	template<class StorageType>
	bool UpdateBInput			(StreamingGrid* B);
	template<class StorageType>
	bool UpdateXOutput			(StreamingGrid* X);
	template<class StorageType>
	bool UpdateBOutput			(StreamingGrid* B);

	void Solve(void);

	void Solve(int idx);
	void SetResidual(int idx);
	void SolveInterior(int idx);
	void SetInteriorResidual(int idx);
	static MultiStreamIOServer server;
#if SEPARTE_TEMP_IO_SERVER
	static MultiStreamIOServer tempServer;
#endif // SEPARTE_TEMP_IO_SERVER
#if NEW_MISHA_CODE
	class Misha<Channels>* misha[NEW_MISHA_CODE];
#endif // NEW_MISHA_CODE
};

template<class Real>
class MultiGridRestrictionNode
{
public:
	virtual void InitRestriction	(bool symmetric){;}							// Generic Initializer
	virtual void SetRestrictionIterations	(int iters){;}						// Sets up interleaved dependencies
	virtual void SetRestriction		(Real* lB,int idx,int major2,int minor2){;}	// Set parent's row
	virtual void SolveRestriction	(void){;}
	virtual bool IterateRestriction	(void){return false;}
};

template<class Real,int Type,int Degree,int Channels,class StorageType=Real>
class MultiGridStreamingSolver : public StreamingSolver<Real,Type,Degree,Channels>, public MultiGridRestrictionNode<Real>
{
#if MOLLIFY_DC
	double dcTerm[Channels];
#endif // MOLLIFY_DC
	typename FiniteElements1D<double,Type,Degree>::FullProlongationStencil majorProlongationStencil;
	typename FiniteElements1D<double,Type,Degree>::FullProlongationStencil minorProlongationStencil;
#if USE_SSE_CODE
#if !NEW_MISHA_CODE
	class ProlongationStencilSSE
	{
	public:
		__declspec (align(ALIGNMENT)) WordClass matrixValues[2][4];
	};
	class ProlongationStencilSSE2
	{
	public:
		__declspec (align(ALIGNMENT)) WordClass matrixValues[2];
	};
#endif // !NEW_MISHA_CODE
	ProlongationStencilSSE *prolongationStencil2,*_prolongationStencil2;
	ProlongationStencilSSE2 *prolongationStencil3,*_prolongationStencil3;
#endif // USE_SSE_CODE
	int prolongationOffset;
	int startProlongation,startRestriction;

	int restrictionBit;
#if MERGE_INTERIOR_STREAMS
	StreamingGrid *I;
#else // !MERGE_INTERIOR_STREAMS
	StreamingGrid *B,*X;
#endif // MERGE_INTERIOR_STREAMS
public:
// There is some messiness in terms of how in/out are defined and used (e.g. do they increment or do they set?)
	StreamingGrid *inX,*outX,*inB,*outB;

	MultiGridRestrictionNode *rChild;
	MultiGridStreamingSolver *parent,*pChild;


	MultiGridStreamingSolver(void);
	~MultiGridStreamingSolver(void);
	void Initialize(
		typename FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::FullDotProductStencil& dotMajor,
		typename FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::FullDotProductStencil& d2DotMajor,
		typename FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::FullDotProductStencil& dotMinor,
		typename FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::FullDotProductStencil& d2DotMinor,
		int major,int minor,bool symmetric,bool memoryMappedFile=false);
	void Initialize(int major,int minor,bool symmetric,bool memoryMappedFile=false);
	void InitProlongation	(bool symmetric);
	void InitRestriction	(bool symmetric);

	void SetProlongationIterations	(int iters);
	void SetRestrictionIterations	(int iters);
	void UnSetProlongationIterations(void);
	void UnSetRestrictionIterations	(void);
	bool IterateProlongation(void);
	bool IterateRestriction	(void);
	void SolveProlongation	(void);
	void SolveRestriction	(void);
	WordClass *localRPtrs[Degree+2];
	void SetRestriction(Real* lB,int idx,int major2,int minor2);
	WordClass *_localRAccum,*localRAccum;
	void SetInteriorRestriction(Real* lowB,int idx,int major2,int minor2);
#if !USE_SSE_CODE
	void SetInteriorRestrictionDotSum(const WordClass mValues[],const WordClass* xPtrs[],size_t j,WordClass& dotSum) const;
	float InteriorRestrictionUpdate0(const WordClass mValues[],const WordClass* rPtrs[],float& previousDotSum,int j);
	float InteriorRestrictionUpdate1(const WordClass mValues[],const WordClass* rPtrs[],float& previousDotSum,int j);
	void SetRestrictionDotSum(const WordClass& mValues,const WordClass* xPtrs,size_t j,WordClass& dotSum) const;
	float RestrictionUpdate0(const WordClass& mValues,const WordClass* rPtrs,float& previousDotSum,int j);
	float RestrictionUpdate1(const WordClass& mValues,const WordClass* rPtrs,float& previousDotSum,int j);
	void SetRestrictionDotSum(const WordClass mValues[],const WordClass* xPtrs[],size_t j,WordClass& dotSum) const;
	float RestrictionUpdate0(const WordClass mValues[],const WordClass* rPtrs[],float& previousDotSum,int j);
	float RestrictionUpdate1(const WordClass mValues[],const WordClass* rPtrs[],float& previousDotSum,int j);
#endif // !USE_SSE_CODE
	void ProlongationUpdate(int idx);

	static bool IsDownSamplable(const int& hMajor,const int& hMinor,int& lMajor,int& lMinor);
};



template<class Real,int Type,int Degree,int Channels,class PartialType=Real>
class StreamingDivergence : public MultiGridRestrictionNode<Real>, RealType<Real>
{
#if REFLECT_GRADIENT
	static const int ReflectBufferSize=40;
#endif // REFLECT_GRADIENT
	typename FiniteElements1D<Real,Type,Degree>::DotProduct<Type,Degree>::FullDotProductStencil dotMajorStencil,dotMinorStencil;
	typename FiniteElements1D<Real,DERIVATIVE(Type),Degree-1>::DotProduct<Type,Degree>::FullDotProductStencil dDotMajorStencil,dDotMinorStencil;
	typename FiniteElements2D<Real,Type,Degree>::FullDivergenceStencil divergenceStencil;
	int startRestriction;
	int major,minor,_major,index,sMajor,sMinor,dSize;
	WordClass *_localDMajor,*_localDMinor;
	Real *localDMajor,*localDMinor;
	WordClass *_localDMajorAccum,*localDMajorAccum;
	WordClass *_localDMinorAccum,*localDMinorAccum;
	Real* localB;
public:
	StreamingGrid *outB;
	double dXSquareNorm,dYSquareNorm,outputSquareNorm;

	StreamingGrid *dMajor,*dMinor;
	MultiGridRStreamingGridestrictionNode<Real> *rParent,*rChild;

	StreamingDivergence(void);
	~StreamingDivergence(void);
	void InitRestriction(int major,int minor,bool symmetric);
	void SetRestrictionIterations(int iters);
	void SetRestriction(Real* lB,int idx,int major2,int minor2);
	void SetInteriorRestriction(Real* lB,int idx,int major2,int minor2);
	void SolveRestriction(void);
	bool IterateRestriction(void);
};

template<class Real,int Type,int Degree,int Channels,class InputType>
class StreamingLaplacian : public MultiGridRestrictionNode<Real>, RealType<Real>
{
	bool _pad;
	int _w,_h,_ww,_hh,_current;
	InputType *_previousPixelsRow,*_previousLabelsRow;
	Real *_dx,*_dy;

	typename FiniteElements1D<Real,Type,Degree>::DotProduct<Type,Degree>::FullDotProductStencil dotMajorStencil,dotMinorStencil;
	typename FiniteElements1D<Real,DERIVATIVE(Type),Degree-1>::DotProduct<Type,Degree>::FullDotProductStencil dDotMajorStencil,dDotMinorStencil;
	typename FiniteElements2D<Real,Type,Degree>::FullDivergenceStencil divergenceStencil;
	int startRestriction;
	int major,minor,_major,index,sMajor,sMinor,dSize;
	WordClass *_localDMajor,*_localDMinor;
	Real *localDMajor,*localDMinor;
	WordClass *_localDMajorAccum,*localDMajorAccum;
	WordClass *_localDMinorAccum,*localDMinorAccum;
public:
	double average[Channels];
	double dXSquareNorm,dYSquareNorm,outputSquareNorm;

	 *pixels,*labels;
	MultiGridRestrictionNode<Real> *rParent,*rChild;

	StreamingLaplacian(void);
	~StreamingLaplacian(void);
	void InitRestriction(int major,int minor,bool symmetric,bool pad);
	void SetRestrictionIterations(int iters);
	void SetRestriction(Real* lB,int idx,int major2,int minor2);
	void SetInteriorRestriction(Real* lB,int idx,int major2,int minor2);
	void SolveRestriction(void);
	bool IterateRestriction(void);
};

#if USE_SSE_CODE
#include "StreamingSolver128.inl"
#endif // USE_SSE_CODE
#include "StreamingSolver.inl"
#endif // STREAMING_SOLVER_INCLUDED
