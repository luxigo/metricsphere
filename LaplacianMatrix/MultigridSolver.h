#ifndef LAPLACIAN_SOLVER_INCLUDED
#define LAPLACIAN_SOLVER_INCLUDED
#include "StreamingSolver.h"

#define NEW_MULTIGRID_CODE 1
#define USE_MY_CG_SOLVER 1
#define MISHA_TEST 0


#if MISHA_TEST
#include "SocketedStreamingSolver.h"
#endif // MISHA_TEST

template<class Real,int Type,int Degree,int Channels>
class MultigridSolver
{
public:
	static const int VERBOSE,FULL_VERBOSE;
	enum
	{
		CONJUGATE_GRADIENTS,
		JACOBI,
		GAUSS_SEIDEL,
		SOLVER_COUNT
	};

	static int verbose;
	static void SolveInCore(Vector<Real>& in,Vector<Real>& out,int w,int h,int multiGrid,int minMGRes,int solverType,bool symmetric,double average[Channels],int cWidth,int cHeight,int fWidth,int fHeight);

	template<class PartialType>
	static void SolveInCore(StreamingGrid* dX,StreamingGrid* dY,Vector<Real>& out,int w,int h,int multiGrid,int minMGRes,int solverType,bool symmetric,double average[Channels],int cWidth,int cHeight,int fWidth,int fHeight,int iters=1);

	template<class PartialType,class StorageType>
	static void SolveOutOfCore(StreamingGrid* dX,StreamingGrid* dY,StreamingGrid* out,int w,int h,int multiGrid,int inCoreRes,int minMGRes,int solverType,bool symmetric,double average[Channels],int cWidth,int cHeight,int iters=1);
	template<class PartialType,class StorageType>
	static void SolveOutOfCore(StreamingGrid* dX,StreamingGrid* dY,StreamingGrid* out,int w,int h,int multiGrid,int inCoreRes,int minMGRes,int solverType,bool symmetric,double average[Channels],int cWidth,int cHeight,
		std::vector<StreamingGrid*>& inStreams,std::vector<StreamingGrid*>& outStreams,int iters=1);

	template<class InputType,class StorageType>
	static void SolveOutOfCore2(StreamingGrid* pixels,StreamingGrid* labels,StreamingGrid* out,int w,int h,int multiGrid,int inCoreRes,int minMGRes,int solverType,bool symmetric,double average[Channels],int cWidth,int cHeight,bool pad,int iters=1);
	template<class InputType,class StorageType>
	static void SolveOutOfCore2(StreamingGrid* pixels,StreamingGrid* labels,StreamingGrid* out,int w,int h,int multiGrid,int inCoreRes,int minMGRes,int solverType,bool symmetric,double average[Channels],int cWidth,int cHeight,bool pad,
		std::vector<StreamingGrid*>& inStreams,std::vector<StreamingGrid*>& outStreams,int iters=1);

	template<class StorageType>
	static void SolveOutOfCore(StreamingGrid* laplacian,StreamingGrid* out,int w,int h,int multiGrid,int inCoreRes,int minMGRes,int solverType,bool symmetric,double average[Channels],int cWidth,int cHeight,int iters=1);
	template<class StorageType>
	static void SolveOutOfCore(StreamingGrid* laplacian,StreamingGrid* out,int w,int h,int multiGrid,int inCoreRes,int minMGRes,int solverType,bool symmetric,double average[Channels],int cWidth,int cHeight,
		std::vector<StreamingGrid*>& inStreams,std::vector<StreamingGrid*>& outStreams,int iters=1);
};
#include "MultigridSolver.inl"
#endif // LAPLACIAN_SOLVER_INCLUDED