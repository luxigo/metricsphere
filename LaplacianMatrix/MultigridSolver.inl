#include "Util/MultiStreamIO.h"


#if USE_MY_CG_SOLVER
// Since the laplacian matrix is indefinite for ZERO_DERIVATIVE constraints because it kills the DC term,
// we update the matrix multiply to add the DC term back in.
template<int Type>
Vector<double> Multiply(const SparseMatrix<double>& A,const Vector<double>& x);

template< int Type >
Vector< double > Multiply( const SparseMatrix< double >& A , const Vector< double >& x )
{
	return A*x;
}
template<> inline Vector<double> Multiply<ZERO_DERIVATIVE>(const SparseMatrix<double>& A,const Vector<double>& x)
{
	Vector<double> y = A*x;
	double average=0;
	for(int i=0;i<x.Dimensions();i++)	average+=x[i];
	average/=x.Dimensions();
	for(int i=0;i<x.Dimensions();i++)	y[i]+=average;
	return y;
}
template<> inline Vector<double> Multiply<PERIODIC_BOUNDARY>(const SparseMatrix<double>& A,const Vector<double>& x)
{
	Vector<double> y = A*x;
	double average=0;
	for(int i=0;i<x.Dimensions();i++)	average+=x[i];
	average/=x.Dimensions();
	for(int i=0;i<x.Dimensions();i++)	y[i]+=average;
	return y;
}
template<class Type>
Vector<double> Multiply(const SparseMatrix<double>& A,const Vector<double>& x)	{	return A*x;	}

// Conjugate gradient code from J. Shewchuk "An Introduction to the Conjugate Gradient Method Without the Agonizing Pain", V. 1 1/4 (1994), pp. 50
template<int Type>
int MySolveConjugateGradient(const SparseMatrix<double>& A,const Vector<double>& b,const int& iters,Vector<double>& x)
{
	double eps=1e-16;
	Vector<double> r = b - A*x;
	Vector<double> d = r;
	double delta_new = r.Dot(r);
	double delta_0 = delta_new;
	int i;
	for(i=0; i<iters && delta_new>eps*delta_0 ;i++)
	{
		Vector<double> q = Multiply<Type>(A,d);
		double alpha = delta_new / d.Dot(q);
		x = x + d*alpha;
		if( !(i%50) )	r = b - Multiply<Type>(A,x);
		else			r = r - q*alpha;

		double delta_old = delta_new;
		delta_new = r.Dot(r);
		double beta = delta_new / delta_old;
		d = r + d*beta;
	}
	return i;
}
#endif // USE_MY_CG_SOLVER

template<class Real,int Type,int Degree,int Channels>
int MultigridSolver<Real,Type,Degree,Channels>::verbose=0;

template<class Real,int Type,int Degree,int Channels>
const int MultigridSolver<Real,Type,Degree,Channels>::VERBOSE=1;
template<class Real,int Type,int Degree,int Channels>
const int MultigridSolver<Real,Type,Degree,Channels>::FULL_VERBOSE=2;

template<class Real,int Type,int Degree,int Channels>
void MultigridSolver<Real,Type,Degree,Channels>::SolveInCore(Vector<Real>& in,Vector<Real>& out,
															 int w,int h,
															 int multiGrid,int minMGRes,
															 int solverType,bool symmetric,
															 double average[Channels],int cWidth,int cHeight,int fWidth,int fHeight)
{
	double t;
	int dCount  = 1;
	int ww=w,hh=h,ww2,hh2;
	MultiGridStreamingSolver<Real,Type,Degree,Channels>* solvers;

	Vector<Real> lowB,lowX;
	StreamingGrid *B,*X;

	// Figure out how many different solvers we will have
	while(MultiGridStreamingSolver<Real,Type,Degree,Channels>::IsDownSamplable(ww,hh,ww2,hh2) && ww2*hh2>=minMGRes*minMGRes)
	{
		dCount++;
		ww=ww2;
		hh=hh2;
	}
	if(multiGrid<=0)	dCount=0;
	solvers=new MultiGridStreamingSolver<Real,Type,Degree,Channels>[dCount];
	for(int i=1;i<dCount;i++)	solvers[i].parent=&solvers[i-1];
	for(int i=0;i<dCount-1;i++)	solvers[i].rChild =solvers[i].pChild =&solvers[i+1];
	for(int i=0;i<dCount;i++)
	{
		solvers[i].bSquareNorm=solvers[i].rSquareNorm=solvers[i].xSquareNorm=0;
		solvers[i].setResidual=true;
	}
	solvers[dCount-1].Initialize(w,h,symmetric,false);

	lowX.Resize(solvers[0].major*solvers[0].minor*Channels);
	lowB.Resize(solvers[0].major*solvers[0].minor*Channels);
	B=new MemoryBackedGrid(&lowB[0],solvers[0].major*Channels*sizeof(Real),solvers[0].minor);
	X=new MemoryBackedGrid(&lowX[0],solvers[0].major*Channels*sizeof(Real),solvers[0].minor);
	solvers[dCount-1].inB=new MemoryBackedGrid(&in[0],w*Channels*sizeof(Real),h);
	// NEW_CODE_HERE
	solvers[dCount-1].inX=new MemoryBackedGrid(&out[0],w*Channels*sizeof(Real),h);
	solvers[0].outB=B;

	// Data for the interleaved streaming multigrid
	t=Time();
	// Initialize
	solvers[dCount-1].InitRestriction(symmetric);
	solvers[dCount-1].SetRestrictionIterations(multiGrid);
	solvers[dCount-1].SolveRestriction();
	delete solvers[dCount-1].inB;
	solvers[dCount-1].inB=NULL;
	// NEW_CODE_HERE
	delete solvers[dCount-1].inX;
	solvers[dCount-1].inX=NULL;

	if(verbose==FULL_VERBOSE)
	{
		printf("In-Core Restriction: %f\n",Time()-t);
		for(int i=dCount-1;i>=0;i--)	printf("\tError[%d x %d] %g -> %g\n",solvers[i].major,solvers[i].minor,sqrt(solvers[i].bSquareNorm),sqrt(solvers[i].rSquareNorm));
	}
	solvers[dCount-1].UnSetRestrictionIterations();
	// Get the base solution
	SparseMatrix<double> lap;
	if(symmetric)	FiniteElements2D<double,Type,Degree>::LaplacianMatrix(solvers[0].major,solvers[0].minor,lap,true,false);
	else			FiniteElements2D<double,Type,Degree>::LaplacianMatrix(solvers[0].major,solvers[0].minor,lap,false,true);
	{
		Vector<MatrixEntry<double> > diagonal;
		if(solverType!=CONJUGATE_GRADIENTS)	FiniteElements2D<double,Type,Degree>::StripDiagonal(lap,diagonal,true);
		Vector<double> myLowX,myLowB;
		myLowB.Resize(solvers[0].major*solvers[0].minor);
		myLowX.Resize(solvers[0].major*solvers[0].minor);
		lowX.Resize(solvers[0].major*solvers[0].minor*Channels);
		for(int c=0;c<Channels;c++)
		{
			for(int i=0;i<solvers[0].major;i++)
				for(int j=0;j<solvers[0].minor;j++)
					myLowB[i+j*solvers[0].major]=lowB[i+j*solvers[0].major*Channels+c*solvers[0].major];
#if USE_MY_CG_SOLVER
			if(solverType==CONJUGATE_GRADIENTS)	MySolveConjugateGradient<Type>(lap,myLowB,6*int(sqrt(lap.groups+1.0)),myLowX);
#else // !USE_MY_CG_SOLVER
			if(solverType==CONJUGATE_GRADIENTS)	SolveConjugateGradient(lap,myLowB,6*int(sqrt(lap.groups+1.0)),myLowX,1e-16);
#endif // USE_MY_CG_SOLVER
			else
				if(solverType==JACOBI)	SparseMatrix<double>::SolveJacobi(lap,diagonal,myLowB,6*int(sqrt(lap.groups+1.0)),myLowX,false);
				else					SparseMatrix<double>::SolveGaussSeidel(lap,diagonal,myLowB,6*int(sqrt(lap.groups+1.0)),myLowX,true);
			for(int i=0;i<solvers[0].major;i++)
				for(int j=0;j<solvers[0].minor;j++)
				{
					lowX[i+j*solvers[0].major*Channels+c*solvers[0].major]=myLowX[i+j*solvers[0].major];
					myLowX[i+j*solvers[0].major]=0;
				}
		}
	}
	for(int i=0;i<dCount;i++)
	{
		solvers[i].bSquareNorm=solvers[i].rSquareNorm=solvers[i].xSquareNorm=0;
		solvers[i].setResidual=verbose;
	}

	// Solve the prolongation
	t=Time();
	solvers[0].inX=X;
	// Set the child dependencies
	solvers[dCount-1].outX=new MemoryBackedGrid(&out[0],w*Channels*sizeof(Real),h);
	solvers[dCount-1].InitProlongation(symmetric);
	solvers[0].SetProlongationIterations(multiGrid);
	solvers[0].SolveProlongation();
	if(verbose==FULL_VERBOSE)
	{
		printf("In-Core Prolongation:    %f\n",Time()-t);
		for(int i=0;i<dCount;i++)	printf("\tError[%d x %d] %g -> %g\n",solvers[i].major,solvers[i].minor,sqrt(solvers[i].bSquareNorm),sqrt(solvers[i].rSquareNorm));
	}
	solvers[0].UnSetProlongationIterations();
	delete solvers[dCount-1].outX;
	delete[] solvers;
	delete B;
	delete X;
	B=X=NULL;

	ww=(cWidth *w+fWidth -1)/fWidth;
	hh=(cHeight*h+fHeight-1)/fHeight;
	for(int c=0;c<Channels;c++)
	{
		double newAverage=0;
		int pCount=0;

		for(int j=0;j<hh;j++)
		{
			Real* o=&out[j*w*Channels+c*w];
			for(int i=0;i<ww;i++)
			{
				newAverage+=o[i];
				pCount++;
			}
		}
		newAverage/=pCount;
		newAverage=average[c]-newAverage;
		for(int j=0;j<h;j++)
		{
			Real* o=&out[j*w*Channels+c*w];
			for(int i=0;i<w;i++)	o[i]+=newAverage;
		}
	}
}
template<class Real,int Type,int Degree,int Channels>
template<class PartialType>
void MultigridSolver<Real,Type,Degree,Channels>::SolveInCore(StreamingGrid* dX,StreamingGrid* dY,
															 Vector<Real>& out,
															 int w,int h,
															 int multiGrid,int minMGRes,
															 int solverType,bool symmetric,
															 double average[Channels],int cWidth,int cHeight,int fWidth,int fHeight,int iters)
{
	double t;
	int dCount  = 1;
	int ww=w,hh=h,ww2,hh2;

	MultiGridStreamingSolver<Real,Type,Degree,Channels>* solvers;
	Vector<Real> lowB,lowX;
	StreamingGrid *B,*X;

	// Figure out how many different solvers we will have
	while(MultiGridStreamingSolver<Real,Type,Degree,Channels>::IsDownSamplable(ww,hh,ww2,hh2) && ww2*hh2>=minMGRes*minMGRes)
	{
		dCount++;
		ww=ww2;
		hh=hh2;
	}
	if(multiGrid<=0)	dCount=0;

	solvers=new MultiGridStreamingSolver<Real,Type,Degree,Channels>[dCount];
	for(int i=1;i<dCount;i++)	solvers[i].parent=&solvers[i-1];
	for(int i=0;i<dCount-1;i++)	solvers[i].rChild =solvers[i].pChild =&solvers[i+1];

	solvers[dCount-1].Initialize(w,h,symmetric,false);
	lowX.Resize(solvers[0].major*solvers[0].minor*Channels);
	lowB.Resize(solvers[0].major*solvers[0].minor*Channels);
	B=new MemoryBackedGrid(&lowB[0],solvers[0].major*Channels*sizeof(Real),solvers[0].minor);
	X=new MemoryBackedGrid(&lowX[0],solvers[0].major*Channels*sizeof(Real),solvers[0].minor);
	solvers[0].outB=B;

	StreamingDivergence<Real,Type,Degree,Channels,PartialType> sDiv;
	sDiv.dMajor=dX;
	sDiv.dMinor=dY;
	MemoryBackedGrid* in=new MemoryBackedGrid(w*sizeof(Real)*Channels,h);
	for(int ii=0;ii<iters;ii++)
	{
		/////////////////
		// RESTRICTION //
		/////////////////
		sDiv.dXSquareNorm=sDiv.dYSquareNorm=sDiv.outputSquareNorm=0;
		for(int i=0;i<dCount;i++)
		{
			solvers[i].bSquareNorm=solvers[i].rSquareNorm=solvers[i].xSquareNorm=0;
			solvers[i].setResidual=true;
		}

		if(ii)	solvers[dCount-1].inX=in;
		else	solvers[dCount-1].inX=NULL;
		solvers[dCount-1].outX=NULL;
		solvers[0].inX=NULL;
		solvers[0].outB=B;
		// Data for the interleaved streaming multigrid
		t=Time();
		// Set the child dependencies
		sDiv.rParent=&solvers[dCount-1];
		solvers[dCount-1].rChild=&sDiv;

		// Initialize
		sDiv.InitRestriction(w,h,symmetric);
		sDiv.SetRestrictionIterations(multiGrid);

		// Solve
		sDiv.SolveRestriction();
		if(verbose==FULL_VERBOSE)
		{
			printf("In-Core Restriction:    %f\n",Time()-t);
			printf("\tInput Divergence: %g %g -> %g\n",sqrt(sDiv.dXSquareNorm),sqrt(sDiv.dYSquareNorm),sqrt(sDiv.outputSquareNorm));
			for(int i=dCount-1;i>=0;i--)	printf("\tError[%d x %d] %g -> %g\n",solvers[i].major,solvers[i].minor,sqrt(solvers[i].bSquareNorm),sqrt(solvers[i].rSquareNorm));
		}
		solvers[dCount-1].UnSetRestrictionIterations();


		// Get the base solution
		SparseMatrix<double> lap;
		if(symmetric)	FiniteElements2D<double,Type,Degree>::LaplacianMatrix(solvers[0].major,solvers[0].minor,lap,true,false);
		else			FiniteElements2D<double,Type,Degree>::LaplacianMatrix(solvers[0].major,solvers[0].minor,lap,false,true);
		{
			Vector<MatrixEntry<double> > diagonal;
			if(solverType!=CONJUGATE_GRADIENTS)	FiniteElements2D<double,Type,Degree>::StripDiagonal(lap,diagonal,true);
			Vector<double> myLowX,myLowB;
			myLowB.Resize(solvers[0].major*solvers[0].minor);
			myLowX.Resize(solvers[0].major*solvers[0].minor);
			lowX.Resize(solvers[0].major*solvers[0].minor*Channels);
			for(int c=0;c<Channels;c++)
			{
				for(int i=0;i<solvers[0].major;i++)
					for(int j=0;j<solvers[0].minor;j++)
						myLowB[i+j*solvers[0].major]=lowB[i+j*solvers[0].major*Channels+c*solvers[0].major];

				// BADNESS!!! Why doesn't this work in floating point precision?
#if USE_MY_CG_SOLVER
				if(solverType==CONJUGATE_GRADIENTS)	MySolveConjugateGradient<Type>(lap,myLowB,6*int(sqrt(lap.groups+1.0)),myLowX);
#else // !USE_MY_CG_SOLVER
				if(solverType==CONJUGATE_GRADIENTS)	SolveConjugateGradient(lap,myLowB,6*int(sqrt(lap.groups+1.0)),myLowX,1e-16);
#endif // USE_MY_CG_SOLVER
				else
					if(solverType==JACOBI)	SparseMatrix<double>::SolveJacobi(lap,diagonal,myLowB,6*int(sqrt(lap.groups+1.0)),myLowX,false);
					else					SparseMatrix<double>::SolveGaussSeidel(lap,diagonal,myLowB,6*int(sqrt(lap.groups+1.0)),myLowX,true);
				for(int i=0;i<solvers[0].major;i++)
					for(int j=0;j<solvers[0].minor;j++)
					{
						lowX[i+j*solvers[0].major*Channels+c*solvers[0].major]=myLowX[i+j*solvers[0].major];
						myLowX[i+j*solvers[0].major]=0;
					}
			}
		}

		//////////////////
		// PROLONGATION //
		//////////////////
		for(int i=0;i<dCount;i++)
		{
			solvers[i].bSquareNorm=solvers[i].rSquareNorm=solvers[i].xSquareNorm=0;
			solvers[i].setResidual=verbose;
		}
		if(ii==iters-1)	solvers[dCount-1].outX=new MemoryBackedGrid(&out[0],w*sizeof(Real)*Channels,h);
		else			solvers[dCount-1].outX=in;
		solvers[dCount-1].inX=NULL;
		solvers[0].inX=X;
		solvers[0].outB=NULL;
		// Solve the prolongation
		t=Time();
		// Set the child dependencies
		// Initialize
		solvers[dCount-1].InitProlongation(symmetric);
		solvers[0].SetProlongationIterations(multiGrid);
		solvers[0].SolveProlongation();
		if(verbose==FULL_VERBOSE)
		{
			printf("In-Core Prolongation:    %f\n",Time()-t);
			for(int i=0;i<dCount;i++)	printf("\tError[%d x %d] %g -> %g\n",solvers[i].major,solvers[i].minor,sqrt(solvers[i].bSquareNorm),sqrt(solvers[i].rSquareNorm));
		}
		else if(verbose==VERBOSE && ii==iters-1)
		{
			int i=dCount-1;
			printf("\tError[%d x %d] %g -> %g\n",solvers[i].major,solvers[i].minor,sqrt(solvers[i].bSquareNorm),sqrt(solvers[i].rSquareNorm));
		}
		solvers[0].UnSetProlongationIterations();
	}
	delete solvers[dCount-1].outX;
	delete[] solvers;
	delete B;
	delete X;
	delete in;
	B=X=NULL;
	ww=(cWidth *w+fWidth -1)/fWidth;
	hh=(cHeight*h+fHeight-1)/fHeight;
	for(int c=0;c<Channels;c++)
	{
		double newAverage=0;
		int pCount=0;

		for(int j=0;j<hh;j++)
		{
			Real* o=&out[j*w*Channels+c*w];
			for(int i=0;i<ww;i++)
			{
				newAverage+=o[i];
				pCount++;
			}
		}
		newAverage/=pCount;
		newAverage=average[c]-newAverage;
		for(int j=0;j<h;j++)
		{
			Real* o=&out[j*w*Channels+c*w];
			for(int i=0;i<w;i++)	o[i]+=newAverage;
		}
	}
}
template<class Real,int Type,int Degree,int Channels>
template<class PartialType,class StorageType>
void MultigridSolver<Real,Type,Degree,Channels>::SolveOutOfCore(StreamingGrid* dX,StreamingGrid* dY,
																StreamingGrid* out,
																int w,int h,
																int multiGrid,int inCoreRes,int minMGRes,
																int solverType,bool symmetric,
																double average[Channels],int cWidth,int cHeight,int iters)
{
	std::vector<StreamingGrid*> inGrids,outGrids;
	SolveOutOfCore<PartialType,StorageType>(dX,dY,out,w,h,multiGrid,inCoreRes,minMGRes,solverType,symmetric,average,cWidth,cHeight,inGrids,outGrids,iters);
}
template<class Real,int Type,int Degree,int Channels>
template<class PartialType,class StorageType>
void MultigridSolver<Real,Type,Degree,Channels>::SolveOutOfCore(StreamingGrid* dX,StreamingGrid* dY,
																StreamingGrid* out,
																int w,int h,
																int multiGrid,int inCoreRes,int minMGRes,
																int solverType,bool symmetric,
																double average[Channels],int cWidth,int cHeight,
																std::vector<StreamingGrid*>& inGrids,std::vector<StreamingGrid*>& outGrids,int iters)
{
	double t;
	int dCount = 1;
	int ww=w,hh=h,ww2,hh2;
	double zeroAverage[Channels];
	for(int i=0;i<Channels;i++)	zeroAverage[i]=0;

	MultiGridStreamingSolver<Real,Type,Degree,Channels,StorageType>* solvers;
	StreamingSolver<Real,Type,Degree,Channels>::server.Reset();
#if SEPARTE_TEMP_IO_SERVER
	StreamingSolver<Real,Type,Degree,Channels>::tempServer.Reset();
#endif // SEPARTE_TEMP_IO_SERVER

	Vector<Real> lowB,lowX;
	StreamingGrid *B,*X;
#if 1
	// BADNESS!!! Why does enabling this screw up the transition from in-core to out-of-core in the prolongation with multiple iterations?
	// Also, the template should be StorageType, not Real
	MultiStreamIOClient* in;
	if(iters>1)	in=new MultiStreamIOClient(w*sizeof(Real)*Channels,h,STREAMING_GRID_BUFFER_MULTIPLIER);
	else		in=NULL;
#else
	Vector<Real> scratch;
	scratch.Resize(w*h*Channels);
	MemoryBackedGrid* in=new MemoryBackedGrid(&scratch[0],w*sizeof(Real)*Channels,h);
#endif
	// Figure out how many different solvers we will have
	while(MultiGridStreamingSolver<Real,Type,Degree,Channels>::IsDownSamplable(ww,hh,ww2,hh2) && long long(ww2)*hh2>=inCoreRes*inCoreRes)
	{
		dCount++;
		ww=ww2;
		hh=hh2;
	}

	if(multiGrid<=0)	dCount=0;
	solvers=new MultiGridStreamingSolver<Real,Type,Degree,Channels,StorageType>[dCount];
	for(int i=1;i<dCount;i++)	solvers[i].parent = &solvers[i-1];
	for(int i=0;i<dCount-1;i++)	solvers[i].rChild =  solvers[i].pChild = &solvers[i+1];
	solvers[dCount-1].Initialize(w,h,symmetric,true);

	lowX.Resize(solvers[0].major*solvers[0].minor*Channels);
	lowB.Resize(solvers[0].major*solvers[0].minor*Channels);
	B=new MemoryBackedGrid(&lowB[0],solvers[0].major*Channels*sizeof(Real),solvers[0].minor);
	X=new MemoryBackedGrid(&lowX[0],solvers[0].major*Channels*sizeof(Real),solvers[0].minor);

	StreamingDivergence<Real,Type,Degree,Channels,PartialType> sDiv;
	sDiv.dMajor=dX;
	sDiv.dMinor=dY;
	for(int ii=0;ii<iters;ii++)
	{
		dX->SetServer(&StreamingSolver<Real,Type,Degree,Channels>::server);
		dY->SetServer(&StreamingSolver<Real,Type,Degree,Channels>::server);
		/////////////////
		// RESTRICTION //
		/////////////////
		sDiv.dXSquareNorm=sDiv.dYSquareNorm=sDiv.outputSquareNorm=0;
		for(int i=0;i<dCount;i++)
		{
			solvers[i].bSquareNorm=solvers[i].rSquareNorm=solvers[i].xSquareNorm=0;
			solvers[i].setResidual=true;
		}
		if(ii)	solvers[dCount-1].inX=in;
		else	solvers[dCount-1].inX=NULL;
		solvers[dCount-1].outX=NULL;
		solvers[0].inX=NULL;
		solvers[0].outB=B;
		solvers[0].outX=X;
		// Data for the interleaved streaming multigrid
		t=Time();
		// Set the child dependencies
		sDiv.rParent=&solvers[dCount-1];
		solvers[dCount-1].rChild=&sDiv;
		// Initialize
		sDiv.InitRestriction(w,h,symmetric);
		sDiv.SetRestrictionIterations(multiGrid);
		// Solve
		// BADNESS!!! Why do I have to comment this out?
//		if(ii)	in->SetServer(&StreamingSolver<Real,Type,Degree,Channels>::server);
		for(int i=0;i<inGrids.size();i++)	inGrids[i]->SetServer(&StreamingSolver<Real,Type,Degree,Channels>::server);
		StreamingSolver<Real,Type,Degree,Channels>::server.StartIO();
#if SEPARTE_TEMP_IO_SERVER
		StreamingSolver<Real,Type,Degree,Channels>::tempServer.StartIO();
#endif // SEPARTE_TEMP_IO_SERVER
		sDiv.SolveRestriction();
		StreamingSolver<Real,Type,Degree,Channels>::server.Reset();
#if SEPARTE_TEMP_IO_SERVER
		StreamingSolver<Real,Type,Degree,Channels>::tempServer.Reset();
#endif // SEPARTE_TEMP_IO_SERVER
		if(verbose==FULL_VERBOSE)
		{
			printf("Out-Of-Core Restriction:    %f\n",Time()-t);
			printf("\tInput Divergence: %g %g -> %g\n",sqrt(sDiv.dXSquareNorm),sqrt(sDiv.dYSquareNorm),sqrt(sDiv.outputSquareNorm));
			for(int i=dCount-1;i>=0;i--)	printf("\tError[%d x %d] %g -> %g\n",solvers[i].major,solvers[i].minor,sqrt(solvers[i].bSquareNorm),sqrt(solvers[i].rSquareNorm));
		}
		if(ii==iters-1 && in)	delete in,	in=NULL;
		solvers[dCount-1].UnSetRestrictionIterations();
		// Get the base solution
		if(ii==iters-1)		SolveInCore(lowB,lowX,solvers[0].major,solvers[0].minor,multiGrid,minMGRes,solverType,symmetric,average,cWidth,cHeight,w,h);
		else				SolveInCore(lowB,lowX,solvers[0].major,solvers[0].minor,multiGrid,minMGRes,solverType,symmetric,zeroAverage,cWidth,cHeight,w,h);

		//////////////////
		// PROLONGATION //
		//////////////////
		for(int i=0;i<dCount;i++)
		{
			solvers[i].bSquareNorm=solvers[i].rSquareNorm=solvers[i].xSquareNorm=0;
			solvers[i].setResidual=verbose;
		}
		solvers[dCount-1].inX=NULL;
		solvers[0].outB=NULL;
		solvers[0].outX=NULL;
		if(ii==iters-1)	solvers[dCount-1].outX=out;
		else			solvers[dCount-1].outX=in;
		solvers[0].inX=X;

		// Solve the prolongation
		t=Time();
		// Set the child dependencies
		// Initialize
		solvers[dCount-1].InitProlongation(symmetric);
		solvers[0].SetProlongationIterations(multiGrid);
		// Solve
		// BADNESS!!! Why do I have to comment this out?
//		if(ii<iters-1)	in->SetServer(&StreamingSolver<Real,Type,Degree,Channels>::server);
		if(ii==iters-1)	for(int i=0;i<outGrids.size();i++)	outGrids[i]->SetServer(&StreamingSolver<Real,Type,Degree,Channels>::server);
		StreamingSolver<Real,Type,Degree,Channels>::server.StartIO();
#if SEPARTE_TEMP_IO_SERVER
		StreamingSolver<Real,Type,Degree,Channels>::tempServer.StartIO();
#endif // SEPARTE_TEMP_IO_SERVER
		solvers[0].SolveProlongation();
		StreamingSolver<Real,Type,Degree,Channels>::server.Reset();
#if SEPARTE_TEMP_IO_SERVER
		StreamingSolver<Real,Type,Degree,Channels>::tempServer.Reset();
#endif // SEPARTE_TEMP_IO_SERVER
		if(verbose==FULL_VERBOSE)
		{
			printf("Out-Of-Core Prolongation:    %f\n",Time()-t);
			for(int i=0;i<dCount;i++)	printf("\tError[%d x %d] %g -> %g\n",solvers[i].major,solvers[i].minor,sqrt(solvers[i].bSquareNorm),sqrt(solvers[i].rSquareNorm));
		}
		else if(verbose==VERBOSE && ii==iters-1)
		{
			int i=dCount-1;
			printf("\tError[%d x %d] %g -> %g\n",solvers[i].major,solvers[i].minor,sqrt(solvers[i].bSquareNorm),sqrt(solvers[i].rSquareNorm));
		}
		solvers[0].UnSetProlongationIterations();
	}
	delete[] solvers;
	delete B;
	delete X;
	B=X=NULL;
}

///////////////////////////////////////////////////
template<class Real,int Type,int Degree,int Channels>
template<class InputType,class StorageType>
void MultigridSolver<Real,Type,Degree,Channels>::SolveOutOfCore2(StreamingGrid* pixels,StreamingGrid* labels,
																 StreamingGrid* out,
																 int w,int h,
																 int multiGrid,int inCoreRes,int minMGRes,
																 int solverType,bool symmetric,
																 double average[Channels],int cWidth,int cHeight,bool pad,int iters)
{
	std::vector<StreamingGrid*> inGrids,outGrids;
	SolveOutOfCore2<InputType,StorageType>(pixels,labels,out,w,h,multiGrid,inCoreRes,minMGRes,solverType,symmetric,average,cWidth,cHeight,pad,inGrids,outGrids,iters);
}
template<class Real,int Type,int Degree,int Channels>
template<class InputType,class StorageType>
void MultigridSolver<Real,Type,Degree,Channels>::SolveOutOfCore2(StreamingGrid* pixels,StreamingGrid* labels,
																 StreamingGrid* out,
																 int w,int h,
																 int multiGrid,int inCoreRes,int minMGRes,
																 int solverType,bool symmetric,
																 double average[Channels],int cWidth,int cHeight,bool pad,
																 std::vector<StreamingGrid*>& inGrids,std::vector<StreamingGrid*>& outGrids,int iters)
{
	double t;
	int dCount = 1;
	int ww=w,hh=h,ww2,hh2;
	double zeroAverage[Channels];
	for(int i=0;i<Channels;i++)	zeroAverage[i]=0;

	MultiGridStreamingSolver<Real,Type,Degree,Channels,StorageType>* solvers;
	StreamingSolver<Real,Type,Degree,Channels>::server.Reset();
#if SEPARTE_TEMP_IO_SERVER
	StreamingSolver<Real,Type,Degree,Channels>::tempServer.Reset();
#endif // SEPARTE_TEMP_IO_SERVER

	Vector<Real> lowB,lowX;
	StreamingGrid *B,*X;
#if 1
	// BADNESS!!! Why does enabling this screw up the transition from in-core to out-of-core in the prolongation with multiple iterations?
	// Also, the template should be StorageType, not Real
	MultiStreamIOClient* in;
	if(iters>1)	in=new MultiStreamIOClient(w*sizeof(Real)*Channels,h,STREAMING_GRID_BUFFER_MULTIPLIER);
	else		in=NULL;
#else
	Vector<Real> scratch;
	scratch.Resize(w*h*Channels);
	MemoryBackedGrid* in=new MemoryBackedGrid(&scratch[0],w*sizeof(Real)*Channels,h);
#endif
	// Figure out how many different solvers we will have
	while(MultiGridStreamingSolver<Real,Type,Degree,Channels>::IsDownSamplable(ww,hh,ww2,hh2) && long long(ww2)*hh2>=inCoreRes*inCoreRes)
	{
		dCount++;
		ww=ww2;
		hh=hh2;
	}

	if(multiGrid<=0)	dCount=0;
	solvers=new MultiGridStreamingSolver<Real,Type,Degree,Channels,StorageType>[dCount];
	for(int i=1;i<dCount;i++)	solvers[i].parent = &solvers[i-1];
	for(int i=0;i<dCount-1;i++)	solvers[i].rChild =  solvers[i].pChild = &solvers[i+1];
	solvers[dCount-1].Initialize(w,h,symmetric,true);

	lowX.Resize(solvers[0].major*solvers[0].minor*Channels);
	lowB.Resize(solvers[0].major*solvers[0].minor*Channels);
	B=new MemoryBackedGrid(&lowB[0],solvers[0].major*Channels*sizeof(Real),solvers[0].minor);
	X=new MemoryBackedGrid(&lowX[0],solvers[0].major*Channels*sizeof(Real),solvers[0].minor);

	StreamingLaplacian<Real,Type,Degree,Channels,InputType> sLap;
	sLap.pixels=pixels;
	sLap.labels=labels;
	for(int ii=0;ii<iters;ii++)
	{
		pixels->SetServer(&StreamingSolver<Real,Type,Degree,Channels>::server);
		labels->SetServer(&StreamingSolver<Real,Type,Degree,Channels>::server);
		/////////////////
		// RESTRICTION //
		/////////////////
		sLap.dXSquareNorm=sLap.dYSquareNorm=sLap.outputSquareNorm=0;
		for(int i=0;i<dCount;i++)
		{
			solvers[i].bSquareNorm=solvers[i].rSquareNorm=solvers[i].xSquareNorm=0;
			solvers[i].setResidual=true;
		}
		if(ii)	solvers[dCount-1].inX=in;
		else	solvers[dCount-1].inX=NULL;
		solvers[dCount-1].outX=NULL;
		solvers[0].inX=NULL;
		solvers[0].outB=B;
		solvers[0].outX=X;
		// Data for the interleaved streaming multigrid
		t=Time();
		// Set the child dependencies
		sLap.rParent=&solvers[dCount-1];
		solvers[dCount-1].rChild=&sLap;
		// Initialize
		sLap.InitRestriction(w,h,symmetric,pad);
		sLap.SetRestrictionIterations(multiGrid);
		// Solve
		// BADNESS!!! Why do I have to comment this out?
//		if(ii)	in->SetServer(&StreamingSolver<Real,Type,Degree,Channels>::server);
		for(int i=0;i<inGrids.size();i++)	inGrids[i]->SetServer(&StreamingSolver<Real,Type,Degree,Channels>::server);
		StreamingSolver<Real,Type,Degree,Channels>::server.StartIO();
#if SEPARTE_TEMP_IO_SERVER
		StreamingSolver<Real,Type,Degree,Channels>::tempServer.StartIO();
#endif // SEPARTE_TEMP_IO_SERVER
		sLap.SolveRestriction();
		StreamingSolver<Real,Type,Degree,Channels>::server.Reset();
#if SEPARTE_TEMP_IO_SERVER
		StreamingSolver<Real,Type,Degree,Channels>::tempServer.Reset();
#endif // SEPARTE_TEMP_IO_SERVER
		if(verbose==FULL_VERBOSE)
		{
			printf("Out-Of-Core Restriction:    %f\n",Time()-t);
			printf("\tInput Divergence: %g %g -> %g\n",sqrt(sLap.dXSquareNorm),sqrt(sLap.dYSquareNorm),sqrt(sLap.outputSquareNorm));
			for(int i=dCount-1;i>=0;i--)	printf("\tError[%d x %d] %g -> %g\n",solvers[i].major,solvers[i].minor,sqrt(solvers[i].bSquareNorm),sqrt(solvers[i].rSquareNorm));
		}
		if(ii==iters-1 && in)	delete in,	in=NULL;
		solvers[dCount-1].UnSetRestrictionIterations();

		for(int c=0;c<Channels;c++)	average[c]=sLap.average[c];
		// Get the base solution
		if(ii==iters-1)		SolveInCore(lowB,lowX,solvers[0].major,solvers[0].minor,multiGrid,minMGRes,solverType,symmetric,average,cWidth,cHeight,w,h);
		else				SolveInCore(lowB,lowX,solvers[0].major,solvers[0].minor,multiGrid,minMGRes,solverType,symmetric,zeroAverage,cWidth,cHeight,w,h);

		//////////////////
		// PROLONGATION //
		//////////////////
		for(int i=0;i<dCount;i++)
		{
			solvers[i].bSquareNorm=solvers[i].rSquareNorm=solvers[i].xSquareNorm=0;
			solvers[i].setResidual=verbose;
		}
		solvers[dCount-1].inX=NULL;
		solvers[0].outB=NULL;
		solvers[0].outX=NULL;
		if(ii==iters-1)	solvers[dCount-1].outX=out;
		else			solvers[dCount-1].outX=in;
		solvers[0].inX=X;

		// Solve the prolongation
		t=Time();
		// Set the child dependencies
		// Initialize
		solvers[dCount-1].InitProlongation(symmetric);
		solvers[0].SetProlongationIterations(multiGrid);
		// Solve
		// BADNESS!!! Why do I have to comment this out?
//		if(ii<iters-1)	in->SetServer(&StreamingSolver<Real,Type,Degree,Channels>::server);
		if(ii==iters-1)	for(int i=0;i<outGrids.size();i++)	outGrids[i]->SetServer(&StreamingSolver<Real,Type,Degree,Channels>::server);
		StreamingSolver<Real,Type,Degree,Channels>::server.StartIO();
#if SEPARTE_TEMP_IO_SERVER
		StreamingSolver<Real,Type,Degree,Channels>::tempServer.StartIO();
#endif // SEPARTE_TEMP_IO_SERVER
		solvers[0].SolveProlongation();
		StreamingSolver<Real,Type,Degree,Channels>::server.Reset();
#if SEPARTE_TEMP_IO_SERVER
		StreamingSolver<Real,Type,Degree,Channels>::tempServer.Reset();
#endif // SEPARTE_TEMP_IO_SERVER
		if(verbose==FULL_VERBOSE)
		{
			printf("Out-Of-Core Prolongation:    %f\n",Time()-t);
			for(int i=0;i<dCount;i++)	printf("\tError[%d x %d] %g -> %g\n",solvers[i].major,solvers[i].minor,sqrt(solvers[i].bSquareNorm),sqrt(solvers[i].rSquareNorm));
		}
		else if(verbose==VERBOSE && ii==iters-1)
		{
			int i=dCount-1;
			printf("\tError[%d x %d] %g -> %g\n",solvers[i].major,solvers[i].minor,sqrt(solvers[i].bSquareNorm),sqrt(solvers[i].rSquareNorm));
		}
		solvers[0].UnSetProlongationIterations();
	}
	delete[] solvers;
	delete B;
	delete X;
	B=X=NULL;

	for(int c=0;c<Channels;c++)	average[c]=sLap.average[c];
}
///////////////////////////////////////////////////
template<class Real,int Type,int Degree,int Channels>
template<class StorageType>
void MultigridSolver<Real,Type,Degree,Channels>::SolveOutOfCore(StreamingGrid* laplacian,StreamingGrid* out,
																int w,int h,
																int multiGrid,int inCoreRes,int minMGRes,
																int solverType,bool symmetric,
																double average[Channels],int cWidth,int cHeight,int iters)
{
	std::vector<StreamingGrid*> inGrids,outGrids;
	SolveOutOfCore<StorageType>(laplacian,out,w,h,multiGrid,inCoreRes,minMGRes,solverType,symmetric,average,cWidth,cHeight,inGrids,outGrids,iters);
}
template<class Real,int Type,int Degree,int Channels>
template<class StorageType>
void MultigridSolver<Real,Type,Degree,Channels>::SolveOutOfCore(StreamingGrid* laplacian,StreamingGrid* out,
																int w,int h,
																int multiGrid,int inCoreRes,int minMGRes,
																int solverType,bool symmetric,
																double average[Channels],int cWidth,int cHeight,
																std::vector<StreamingGrid*>& inGrids,std::vector<StreamingGrid*>& outGrids,int iters)
{
	double t;
	int dCount = 1;
	int ww=w,hh=h,ww2,hh2;
	double zeroAverage[Channels];
	for(int i=0;i<Channels;i++)	zeroAverage[i]=0;

#if MISHA_TEST
	SocketedMultiGridStreamingSolver<Channels,StorageType>* solvers;
	SocketedStreamingSolver<Channels>::server.Reset();
#else // !MISHA_TEST
	MultiGridStreamingSolver<Real,Type,Degree,Channels,StorageType>* solvers;
	StreamingSolver<Real,Type,Degree,Channels>::server.Reset();
#endif // MISHA_TEST
#if SEPARTE_TEMP_IO_SERVER
#if MISHA_TEST
	SocketedStreamingSolver<Real,Type,Degree,Channels>::tempServer.Reset();
#else // !MISHA_TEST
	StreamingSolver<Real,Type,Degree,Channels>::tempServer.Reset();
#endif // MISHA_TEST
#endif // SEPARTE_TEMP_IO_SERVER

	Vector<Real> lowB,lowX;
	StreamingGrid *B,*X;
	// Also, the template should be StorageType, not Real
	MultiStreamIOClient* in;
	if(iters>1)	in=new MultiStreamIOClient(w*sizeof(Real)*Channels,h,STREAMING_GRID_BUFFER_MULTIPLIER);
	else		in=NULL;
	// Figure out how many different solvers we will have
#if MISHA_TEST && 0
	while(SocketedMultiGridStreamingSolver<Channels>::IsDownSamplable(ww,hh,ww2,hh2) && long long(ww2)*hh2>=inCoreRes*inCoreRes)
#else // !MISHA_TEST
	while(MultiGridStreamingSolver<Real,Type,Degree,Channels>::IsDownSamplable(ww,hh,ww2,hh2) && long long(ww2)*hh2>=inCoreRes*inCoreRes)
#endif // MISHA_TEST
	{
		dCount++;
		ww=ww2;
		hh=hh2;
	}

	if(multiGrid<=0)	dCount=0;
#if MISHA_TEST
	solvers=new SocketedMultiGridStreamingSolver<Channels,StorageType>[dCount];
#else // !MISHA_TEST
	solvers=new MultiGridStreamingSolver<Real,Type,Degree,Channels,StorageType>[dCount];
#endif // MISHA_TEST
	for(int i=1;i<dCount;i++)	solvers[i].parent = &solvers[i-1];
	for(int i=0;i<dCount-1;i++)	solvers[i].rChild =  solvers[i].pChild = &solvers[i+1];
#if MISHA_TEST
	solvers[dCount-1].Initialize(0,w,w,h,NULL,NULL,true);
#else // !MISHA_TEST
	solvers[dCount-1].Initialize(w,h,symmetric,true);
#endif // MISHA_TEST

	lowX.Resize(solvers[0].major*solvers[0].minor*Channels);
	lowB.Resize(solvers[0].major*solvers[0].minor*Channels);
	B=new MemoryBackedGrid(&lowB[0],solvers[0].major*Channels*sizeof(Real),solvers[0].minor);
	X=new MemoryBackedGrid(&lowX[0],solvers[0].major*Channels*sizeof(Real),solvers[0].minor);

	for(int ii=0;ii<iters;ii++)
	{
#if MISHA_TEST
		laplacian->SetServer(&SocketedStreamingSolver<Channels>::server);
#else // !MISHA_TEST
		laplacian->SetServer(&StreamingSolver<Real,Type,Degree,Channels>::server);
#endif // MISHA_TEST
		/////////////////
		// RESTRICTION //
		/////////////////
		for(int i=0;i<dCount;i++)
		{
			solvers[i].bSquareNorm=solvers[i].rSquareNorm=solvers[i].xSquareNorm=0;
			solvers[i].setResidual=true;
		}
		solvers[dCount-1].inB=laplacian;
		if(ii)	solvers[dCount-1].inX=in;
		else	solvers[dCount-1].inX=NULL;
		solvers[dCount-1].outX=NULL;
		solvers[0].inX=NULL;
		solvers[0].outB=B;
		solvers[0].outX=X;
		// Data for the interleaved streaming multigrid
		t=Time();
		// Set the child dependencies
		solvers[dCount-1].rChild=NULL;
		// Initialize
#if MISHA_TEST
		solvers[dCount-1].InitRestriction();
#else // !MISHA_TEST
		solvers[dCount-1].InitRestriction(symmetric);
#endif // MISHA_TEST
		solvers[dCount-1].SetRestrictionIterations(multiGrid);
		// Solve
		// BADNESS!!! Why do I have to comment this out?
#if MISHA_TEST
//		if(ii)	in->SetServer(&SocketedStreamingSolver<Channels>::server);
		for(int i=0;i<inGrids.size();i++)	inGrids[i]->SetServer(&SocketedStreamingSolver<Channels>::server);
		SocketedStreamingSolver<Channels>::server.StartIO();
#else // !MISHA_TEST
//		if(ii)	in->SetServer(&StreamingSolver<Real,Type,Degree,Channels>::server);
		for(int i=0;i<inGrids.size();i++)	inGrids[i]->SetServer(&StreamingSolver<Real,Type,Degree,Channels>::server);
		StreamingSolver<Real,Type,Degree,Channels>::server.StartIO();
#endif // MISHA_TEST
#if SEPARTE_TEMP_IO_SERVER
#if MISHA_TEST
		SocketedStreamingSolver<Channels>::tempServer.StartIO();
#else // !MISHA_TEST
		StreamingSolver<Real,Type,Degree,Channels>::tempServer.StartIO();
#endif // MISHA_TEST
#endif // SEPARTE_TEMP_IO_SERVER
		solvers[dCount-1].SolveRestriction();
#if MISHA_TEST
		SocketedStreamingSolver<Channels>::server.Reset();
#else // !MISHA_TEST
		StreamingSolver<Real,Type,Degree,Channels>::server.Reset();
#endif // MISHA_TEST
#if SEPARTE_TEMP_IO_SERVER
#if MISHA_TEST
		SocketedStreamingSolver<Channels>::tempServer.Reset();
#else // !MISHA_TEST
		StreamingSolver<Real,Type,Degree,Channels>::tempServer.Reset();
#endif // MISHA_TEST
#endif // SEPARTE_TEMP_IO_SERVER
		if(verbose==FULL_VERBOSE)
		{
			printf("Out-Of-Core Restriction:    %f\n",Time()-t);
			for(int i=dCount-1;i>=0;i--)	printf("\tError[%d x %d] %g -> %g\n",solvers[i].major,solvers[i].minor,sqrt(solvers[i].bSquareNorm),sqrt(solvers[i].rSquareNorm));
		}
		if(ii==iters-1 && in)	delete in,	in=NULL;
		solvers[dCount-1].UnSetRestrictionIterations();

		// Get the base solution
		if(ii==iters-1)		SolveInCore(lowB,lowX,solvers[0].major,solvers[0].minor,multiGrid,minMGRes,solverType,symmetric,average,cWidth,cHeight,w,h);
		else				SolveInCore(lowB,lowX,solvers[0].major,solvers[0].minor,multiGrid,minMGRes,solverType,symmetric,zeroAverage,cWidth,cHeight,w,h);

		//////////////////
		// PROLONGATION //
		//////////////////
		for(int i=0;i<dCount;i++)
		{
			solvers[i].bSquareNorm=solvers[i].rSquareNorm=solvers[i].xSquareNorm=0;
			solvers[i].setResidual=verbose;
		}
		solvers[dCount-1].inX=NULL;
		solvers[0].outB=NULL;
		solvers[0].outX=NULL;
		if(ii==iters-1)	solvers[dCount-1].outX=out;
		else			solvers[dCount-1].outX=in;
		solvers[0].inX=X;

		// Solve the prolongation
		t=Time();
		// Set the child dependencies
		// Initialize
#if MISHA_TEST
		solvers[dCount-1].InitProlongation();
#else // !MISHA_TEST
		solvers[dCount-1].InitProlongation(symmetric);
#endif // MISHA_TEST
		solvers[0].SetProlongationIterations(multiGrid);
		// Solve
		// BADNESS!!! Why do I have to comment this out?
//		if(ii<iters-1)	in->SetServer(&StreamingSolver<Real,Type,Degree,Channels>::server);
#if MISHA_TEST
		if(ii==iters-1)	for(int i=0;i<outGrids.size();i++)	outGrids[i]->SetServer(&SocketedStreamingSolver<Channels>::server);
		SocketedStreamingSolver<Channels>::server.StartIO();
#else // !MISHA_TEST
		if(ii==iters-1)	for(int i=0;i<outGrids.size();i++)	outGrids[i]->SetServer(&StreamingSolver<Real,Type,Degree,Channels>::server);
		StreamingSolver<Real,Type,Degree,Channels>::server.StartIO();
#endif // MISHA_TEST
#if SEPARTE_TEMP_IO_SERVER
#if MISHA_TEST
		SocketedStreamingSolver<Channels>::tempServer.StartIO();
#else // !MISHA_TEST
		StreamingSolver<Real,Type,Degree,Channels>::tempServer.StartIO();
#endif // MISHA_TEST
#endif // SEPARTE_TEMP_IO_SERVER
		solvers[0].SolveProlongation();
#if MISHA_TEST
		SocketedStreamingSolver<Channels>::server.Reset();
#else // !MISHA_TEST
		StreamingSolver<Real,Type,Degree,Channels>::server.Reset();
#endif // MISHA_TEST
#if SEPARTE_TEMP_IO_SERVER
#if MISHA_TEST
		SocketedStreamingSolver<Channels>::tempServer.Reset();
#else // !MISHA_TEST
		StreamingSolver<Real,Type,Degree,Channels>::tempServer.Reset();
#endif // MISHA_TEST
#endif // SEPARTE_TEMP_IO_SERVER
		if(verbose==FULL_VERBOSE)
		{
			printf("Out-Of-Core Prolongation:    %f\n",Time()-t);
			for(int i=0;i<dCount;i++)
				printf("\tError[%d x %d] %g -> %g\n",solvers[i].major,solvers[i].minor,sqrt(solvers[i].bSquareNorm),sqrt(solvers[i].rSquareNorm));
		}
		else if(verbose==VERBOSE && ii==iters-1)
		{
			int i=dCount-1;
			printf("\tError[%d x %d] %g -> %g\n",solvers[i].major,solvers[i].minor,sqrt(solvers[i].bSquareNorm),sqrt(solvers[i].rSquareNorm));
		}
		solvers[0].UnSetProlongationIterations();
	}
	delete[] solvers;
	delete B;
	delete X;
	B=X=NULL;
}
