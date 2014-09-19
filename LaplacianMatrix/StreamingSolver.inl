#define MyModIndex(idx,mod) (ModIndex(idx,mod))

/////////////////////
// StreamingSolver //
/////////////////////
template<class Real,int Type,int Degree,int Channels>
MultiStreamIOServer StreamingSolver<Real,Type,Degree,Channels>::server;
#if SEPARTE_TEMP_IO_SERVER
template<class Real,int Type,int Degree,int Channels>
MultiStreamIOServer StreamingSolver<Real,Type,Degree,Channels>::tempServer;
#endif // SEPARTE_TEMP_IO_SERVER
template<class Real,int Type,int Degree,int Channels>
StreamingSolver<Real,Type,Degree,Channels>::StreamingSolver(void)
{
	setResidual=true;
	laplacianRescale=true;
	bSize=rSize=0;
	_localR=NULL;
#if MERGE_INTERIOR_STREAMS
	localR=IStream=NULL;
	_IStream=NULL;
#else // !MERGE_INTERIOR_STREAMS
	localR=XStream=BStream=NULL;
	_XStream=_BStream=NULL;
#endif // MERGE_INTERIOR_STREAMS
#if USE_SSE_CODE
	_lapTemplates2=(LaplacianTemplateSSE*)malloc(sizeof(LaplacianTemplateSSE)*3*(2*Degree+1)+ALIGNMENT-1);
	memset(_lapTemplates2,0,sizeof(LaplacianTemplateSSE)*3*(2*Degree+1)+ALIGNMENT-1);
	lapTemplates2=(LaplacianTemplateSSE*)(((size_t)(_lapTemplates2)+ALIGNMENT-1) & ~(ALIGNMENT-1));

	_lapTemplates3=(LaplacianTemplateSSE*)malloc(sizeof(LaplacianTemplateSSE)*3*(2*Degree+1)+ALIGNMENT-1);
	memset(_lapTemplates3,0,sizeof(LaplacianTemplateSSE)*3*(2*Degree+1)+ALIGNMENT-1);
	lapTemplates3=(LaplacianTemplateSSE*)(((size_t)(_lapTemplates3)+ALIGNMENT-1) & ~(ALIGNMENT-1));
#endif // USE_SSE_CODE
	for(int i=0;i<Degree;i++)	localXAccum[i]=_localXAccum[i]=NULL;
}
template<class Real,int Type,int Degree,int Channels>
StreamingSolver<Real,Type,Degree,Channels>::~StreamingSolver(void)
{
#if MERGE_INTERIOR_STREAMS
	if(_IStream)	free(_IStream);
#else // !MERGE_INTERIOR_STREAMS
	if(_XStream)	free(_XStream);
	if(_BStream)	free(_BStream);
#endif // MERGE_INTERIOR_STREAMS
#if MERGE_INTERIOR_STREAMS
	_IStream=NULL;
	localR=IStream=NULL;
#else // !MERGE_INTERIOR_STREAMS
	_XStream=_BStream=NULL;
	localR=XStream=BStream=NULL;
#endif // MERGE_INTERIOR_STREAMS
	if(_localR)	free(_localR);
	_localR=NULL;
#if USE_SSE_CODE
	free(_lapTemplates2);
	free(_lapTemplates3);
#endif // USE_SSE_CODE
	for(int i=0;i<Degree;i++)
	{
		if(_localXAccum[i])	free(_localXAccum[i]);
		_localXAccum[i]=localXAccum[i]=NULL;
	}
#if NEW_MISHA_CODE
	for(int i=0;i<NEW_MISHA_CODE;i++)
	{
		if(misha[i])	delete misha[i];
		misha[i]=NULL;
	}
#endif // NEW_MISHA_CODE

}
template<class Real,int Type,int Degree,int Channels>
Real* StreamingSolver<Real,Type,Degree,Channels>::GetXRow(int row,int channel)
{
#if MERGE_INTERIOR_STREAMS
	return &IStream[MyModIndex(row,iSize)*2*_major*RealPerWord*Channels+channel*_major*RealPerWord+_major*RealPerWord*Channels];
#else // !MERGE_INTERIOR_STREAMS
	return &XStream[MyModIndex(row,xSize)*  _major*RealPerWord*Channels+channel*_major*RealPerWord];
#endif // MERGE_INTERIOR_STREAMS
}
template<class Real,int Type,int Degree,int Channels>
Real* StreamingSolver<Real,Type,Degree,Channels>::GetBRow(int row,int channel)
{
#if MERGE_INTERIOR_STREAMS
	return &IStream[MyModIndex(row,iSize)*2*_major*RealPerWord*Channels+channel*_major*RealPerWord];
#else // !MERGE_INTERIOR_STREAMS
	return &BStream[MyModIndex(row,bSize)*  _major*RealPerWord*Channels+channel*_major*RealPerWord];
#endif // MERGE_INTERIOR_STREAMS
}
template<class Real,int Type,int Degree,int Channels>
void StreamingSolver<Real,Type,Degree,Channels>::Init(int major,int minor,bool symmetric)
{
	FiniteElements2D<double,Type,Degree>::FullMatrixStencil lStencil;
	if(symmetric)	FiniteElements2D<double,Type,Degree>::LaplacianStencil(major,minor,lStencil,true);
	else			FiniteElements2D<double,Type,Degree>::LaplacianStencil(major,minor,lStencil,false,true);
	Init(lStencil,major,minor,symmetric);
}
template<class Real,int Type,int Degree,int Channels>
void StreamingSolver<Real,Type,Degree,Channels>::Init(const typename FiniteElements2D<double,Type,Degree>::FullMatrixStencil& lStencil,int major,int minor,bool symmetric)
{
	this->major=major;
	this->minor=minor;
	this->symmetric=symmetric;
	_major=(major+RealPerWord-1)/RealPerWord;
#if USE_SSE_CODE
	FiniteElements2D<double,Type,Degree>::FullMatrixStencil laplacianStencil;
#endif // USE_SSE_CODE
	laplacianStencil=lStencil;
	laplacianScale =Real(1.0)/laplacianStencil.caseTable[Degree][Degree].values[Degree][Degree];
	laplacianScaleR=          laplacianStencil.caseTable[Degree][Degree].values[Degree][Degree];
	for(int i=0;i<=2*Degree;i++)
		for(int j=0;j<=2*Degree;j++)
			for(int k=0;k<=2*Degree;k++)
				for(int l=0;l<=2*Degree;l++)
					laplacianStencil.caseTable[i][j].values[k][l]*=laplacianScale;
#if USE_SSE_CODE
	// BADNESS!!! should fix the indexing of lapTemplates so that major index iterates faster
	// Badness!!! Assuming Degree=2
	__declspec (align(16)) float scratch[4];
	if(((major>>2)<<2)!=major)
	{
		fprintf(stderr,"Error: %d is not a multiple of four!!!\n",major);
		exit(0);
	}
	for(int i=0;i<=2*Degree;i++)	// Iterate over the minor index in the mask
	{
		for(int k=0;k<4;k++)
			lapTemplates2[3*i+1].diagonal[k]=1.f/laplacianStencil.caseTable[i][Degree].values[Degree][Degree];
		lapTemplates2[3*i  ].diagonal[0]=1.f/laplacianStencil.caseTable[i][0].values[Degree][Degree];
		lapTemplates2[3*i  ].diagonal[1]=1.f/laplacianStencil.caseTable[i][1].values[Degree][Degree];
		lapTemplates2[3*i  ].diagonal[2]=1.f/laplacianStencil.caseTable[i][Degree].values[Degree][Degree];
		lapTemplates2[3*i  ].diagonal[3]=1.f/laplacianStencil.caseTable[i][Degree].values[Degree][Degree];
		lapTemplates2[3*i+2].diagonal[0]=1.f/laplacianStencil.caseTable[i][Degree].values[Degree][Degree];
		lapTemplates2[3*i+2].diagonal[1]=1.f/laplacianStencil.caseTable[i][Degree].values[Degree][Degree];
		lapTemplates2[3*i+2].diagonal[2]=1.f/laplacianStencil.caseTable[i][2*Degree-1].values[Degree][Degree];
		lapTemplates2[3*i+2].diagonal[3]=1.f/laplacianStencil.caseTable[i][2*Degree  ].values[Degree][Degree];
		for(int k=0;k<4;k++)
			lapTemplates3[3*i+1].diagonal[k]=1.f/laplacianStencil.caseTable[i][Degree].values[Degree][Degree];
		lapTemplates3[3*i  ].diagonal[0]=1.f/laplacianStencil.caseTable[i][0].values[Degree][Degree];
		lapTemplates3[3*i  ].diagonal[1]=1.f/laplacianStencil.caseTable[i][1].values[Degree][Degree];
		lapTemplates3[3*i  ].diagonal[2]=1.f/laplacianStencil.caseTable[i][Degree].values[Degree][Degree];
		lapTemplates3[3*i  ].diagonal[3]=1.f/laplacianStencil.caseTable[i][Degree].values[Degree][Degree];
		lapTemplates3[3*i+2].diagonal[0]=1.f/laplacianStencil.caseTable[i][Degree].values[Degree][Degree];
		lapTemplates3[3*i+2].diagonal[1]=1.f/laplacianStencil.caseTable[i][Degree].values[Degree][Degree];
		lapTemplates3[3*i+2].diagonal[2]=1.f/laplacianStencil.caseTable[i][2*Degree-1].values[Degree][Degree];
		lapTemplates3[3*i+2].diagonal[3]=1.f/laplacianStencil.caseTable[i][2*Degree  ].values[Degree][Degree];
	}
	for(int i=0;i<=2*Degree;i++)	// Iterate over the minor index in the mask
		for(int j=0;j<3;j++)
			for(int k=0;k<=2*Degree;k++)
			{
				int jj;

				if(j==0)	jj=0;
				else		jj=Degree;
				scratch[0]=laplacianStencil.caseTable[i][jj].values[k][2];
				scratch[1]=laplacianStencil.caseTable[i][jj].values[k][3];
				scratch[2]=laplacianStencil.caseTable[i][jj].values[k][4];
				if(j!=2)	scratch[3]=laplacianStencil.caseTable[i][Degree].values[k][1];
				else		scratch[3]=0;
				lapTemplates2[3*i+j].matrixValues[0][k]=_mm_load_ps(scratch);

				if(j==0)	jj=1;
				else		jj=Degree;
				for(int l=0;l<4;l++)	scratch[l]=laplacianStencil.caseTable[i][jj].values[k][l+1];
				lapTemplates2[3*i+j].matrixValues[1][k]=_mm_load_ps(scratch);

				if(j==2)	jj=2*Degree-1;
				else		jj=Degree;
				for(int l=0;l<4;l++)	scratch[l]=laplacianStencil.caseTable[i][jj].values[k][l];
				lapTemplates2[3*i+j].matrixValues[2][k]=_mm_load_ps(scratch);

				if(j==2)	jj=2*Degree;
				else		jj=Degree;
				if(j!=0)	scratch[0]=laplacianStencil.caseTable[i][Degree].values[k][3];
				else		scratch[0]=0;
				scratch[1]=laplacianStencil.caseTable[i][jj].values[k][0];
				scratch[2]=laplacianStencil.caseTable[i][jj].values[k][1];
				scratch[3]=laplacianStencil.caseTable[i][jj].values[k][2];
				lapTemplates2[3*i+j].matrixValues[3][k]=_mm_load_ps(scratch);
			}
	for(int i=0;i<=2*Degree;i++)
		for(int j=0;j<=2*Degree;j++)
			laplacianStencil.caseTable[i][j].values[Degree][Degree]=0;
	for(int i=0;i<=2*Degree;i++)	// Iterate over the minor index in the mask
		for(int j=0;j<3;j++)
			for(int k=0;k<=2*Degree;k++)
			{
				int jj;

				if(j==0)	jj=0;
				else		jj=Degree;
				scratch[0]=laplacianStencil.caseTable[i][jj].values[k][2];
				scratch[1]=laplacianStencil.caseTable[i][jj].values[k][3];
				scratch[2]=laplacianStencil.caseTable[i][jj].values[k][4];
				if(j!=2)	scratch[3]=laplacianStencil.caseTable[i][Degree].values[k][1];
				else		scratch[3]=0;
				lapTemplates3[3*i+j].matrixValues[0][k]=_mm_load_ps(scratch);

				if(j==0)	jj=1;
				else		jj=Degree;
				for(int l=0;l<4;l++)	scratch[l]=laplacianStencil.caseTable[i][jj].values[k][l+1];
				lapTemplates3[3*i+j].matrixValues[1][k]=_mm_load_ps(scratch);

				if(j==2)	jj=2*Degree-1;
				else		jj=Degree;
				for(int l=0;l<4;l++)	scratch[l]=laplacianStencil.caseTable[i][jj].values[k][l];
				lapTemplates3[3*i+j].matrixValues[2][k]=_mm_load_ps(scratch);

				if(j==2)	jj=2*Degree;
				else		jj=Degree;
				if(j!=0)	scratch[0]=laplacianStencil.caseTable[i][Degree].values[k][3];
				else		scratch[0]=0;
				scratch[1]=laplacianStencil.caseTable[i][jj].values[k][0];
				scratch[2]=laplacianStencil.caseTable[i][jj].values[k][1];
				scratch[3]=laplacianStencil.caseTable[i][jj].values[k][2];
				lapTemplates3[3*i+j].matrixValues[3][k]=_mm_load_ps(scratch);
			}
#endif // USE_SSE_CODE
#if NEW_MISHA_CODE
	for(int i=0;i<NEW_MISHA_CODE;i++)	misha[i]=NULL;
	if(major>3000)
	{
		if(!(major&3))
		{
			int subMajor,subSubMajor;
			subMajor=(major>>3)<<2;
			if(!(subMajor&3))
			{
#if NEW_MISHA_CODE == 1
				misha[0]=new Misha<Channels>(major,0,major,lapTemplates2,lapTemplates3);
#elif NEW_MISHA_CODE == 2
				if(!(subMajor&3))
				{
					misha[0]=new Misha<Channels>(major,0,subMajor,lapTemplates2,lapTemplates3);
					misha[1]=new Misha<Channels>(major,subMajor,major,lapTemplates2,lapTemplates3);
				}
#elif NEW_MISHA_CODE == 4
				subSubMajor=(subMajor>>3)<<2;
				misha[0]=new Misha<Channels>(major,0,subSubMajor,lapTemplates2,lapTemplates3);
				misha[1]=new Misha<Channels>(major,subSubMajor,subMajor,lapTemplates2,lapTemplates3);
				subMajor=major-subMajor;
				if(!(subMajor&3))
				{
					subSubMajor=(subMajor>>3)<<2;
					misha[2]=new Misha<Channels>(major,major-subMajor,major-subSubMajor,lapTemplates2,lapTemplates3);
					misha[3]=new Misha<Channels>(major,major-subSubMajor,major,lapTemplates2,lapTemplates3);
				}
#endif
			}
		}
	}
#endif // NEW_MISHA_CODE
}
template<class Real,int Type,int Degree,int Channels>
void StreamingSolver<Real,Type,Degree,Channels>::SetIterations(int iters,int rSize)
{
	SetIterations(-iters*Degree-Degree+1,iters,0,0,0,0,rSize);
}
template<class Real,int Type,int Degree,int Channels>
void StreamingSolver<Real,Type,Degree,Channels>::SetIterations(int start,int iters,int bStart,int xStart,int bEnd,int xEnd,int rSize)
{
	clearX=clearB=true;
	index=start;
	this->iters=iters;
	_iters=iters*Degree;
#if MERGE_INTERIOR_STREAMS
	int bSize=0,xSize=0;
#else // !MERGE_INTERIOR_STREAMS
	bSize=xSize=0;
#endif // MERGE_INTERIOR_STREAMS

	xSize+=(_iters+Degree)>xEnd?(_iters+Degree):xEnd;
	bSize+=_iters>bEnd?_iters:bEnd;
	if(setResidual)
	{
		xSize+=-Degree-1<xStart?Degree+1:-xStart;
		bSize+=       -1<bStart?       1:-bStart;
	}
	else
	{
		xSize+=-Degree<xStart?Degree:-xStart;
		bSize+=0<bStart?0:-bStart;
	}
#if MERGE_INTERIOR_STREAMS
	iSize=xSize>bSize?xSize:bSize;
	if(iSize>minor)
	{
		iSize=minor;
		clearX=clearB=false;
	}
#else // !MERGE_INTERIOR_STREAMS
	if(bSize>minor)
	{
		bSize=minor;
		clearB=false;
	}
	if(xSize>minor)
	{
		xSize=minor;
		clearX=false;
	}
#endif // MERGE_INTERIOR_STREAMS
	if(_localR)	free(_localR);
	localR=NULL;
	_localR=NULL;
#if MERGE_INTERIOR_STREAMS
	if(_IStream)	free(_IStream);
	IStream=NULL;
	_IStream=NULL;
#else // !MERGE_INTERIOR_STREAMS
	if(_XStream)	free(_XStream);
	if(_BStream)	free(_BStream);
	BStream=XStream=NULL;
	_BStream=_XStream=NULL;
#endif // MERGE_INTERIOR_STREAMS

	this->rSize=rSize;
	if(setResidual)
	{
		_localR=(WordClass*)malloc(sizeof(WordClass)*(rSize*_major*Channels+2*Degree)+ALIGNMENT-1);
		localR=(Real*)(((size_t)(_localR+Degree)+ALIGNMENT-1) & ~(ALIGNMENT-1));
		memset(_localR,0,sizeof(WordClass)*(rSize*_major*Channels+2*Degree)+ALIGNMENT-1);
	}
#if MERGE_INTERIOR_STREAMS
	_IStream=(WordClass*)malloc(sizeof(WordClass)*(2*iSize*_major*Channels+2*Degree)+ALIGNMENT-1);
	IStream=(Real*)(((size_t)(_IStream+Degree)+ALIGNMENT-1) & ~(ALIGNMENT-1));
	memset(_IStream,0,sizeof(WordClass)*(2*iSize*_major*Channels+2*Degree)+ALIGNMENT-1);
#else // !MERGE_INTERIOR_STREAMS
	_XStream=(WordClass*)malloc(sizeof(WordClass)*(xSize*_major*Channels+2*Degree)+ALIGNMENT-1);
	XStream=(Real*)(((size_t)(_XStream+Degree)+ALIGNMENT-1) & ~(ALIGNMENT-1));
	memset(_XStream,0,sizeof(WordClass)*(xSize*_major*Channels+2*Degree)+ALIGNMENT-1);
	_BStream=(WordClass*)malloc(sizeof(WordClass)*(bSize*_major*Channels+2*Degree)+ALIGNMENT-1);
	BStream=(Real*)(((size_t)(_BStream+Degree)+ALIGNMENT-1) & ~(ALIGNMENT-1));
	memset(_BStream,0,sizeof(WordClass)*(bSize*_major*Channels+2*Degree)+ALIGNMENT-1);
#endif // MERGE_INTERIOR_STREAMS
	for(int i=0;i<Degree;i++)
	{
		_localXAccum[i]=(WordClass*)malloc(sizeof(WordClass)*(_major+2*Degree)+ALIGNMENT-1);
		localXAccum[i]=(WordClass*)(((size_t)(_localXAccum[i])+ALIGNMENT-1) & ~(ALIGNMENT-1));
	}
}
template<class Real,int Type,int Degree,int Channels>
void StreamingSolver<Real,Type,Degree,Channels>::UnSetIterations(void)
{
	clearX=clearB=true;
	index=0;
	iters=_iters=0;
	rSize=0;
	setResidual=false;
#if MERGE_INTERIOR_STREAMS
	int bSize=0,xSize=0;
#else // !MERGE_INTERIOR_STREAMS
	bSize=xSize=0;
#endif // MERGE_INTERIOR_STREAMS
	if(_localR)	free(_localR);
	localR=NULL;
	_localR=NULL;
#if MERGE_INTERIOR_STREAMS
	if(_IStream)	free(_IStream);
	IStream=NULL;
	_IStream=NULL;
#else // !MERGE_INTERIOR_STREAMS
	if(_XStream)	free(_XStream);
	if(_BStream)	free(_BStream);
	BStream=XStream=NULL;
	_BStream=_XStream=NULL;
#endif // MERGE_INTERIOR_STREAMS
	for(int i=0;i<Degree;i++)
	{
		if(_localXAccum[i])	free(_localXAccum[i]);
		_localXAccum[i]=NULL;
		localXAccum[i]=NULL;
	}
}
#if !USE_SSE_CODE
template<class Real,int Type,int Degree,int Channels>
inline Real StreamingSolver<Real,Type,Degree,Channels>::GetLaplacianValue(const double lapValues[][2*Degree+1],int i)
{
	Real temp=0;
	switch(Degree)
	{
		case 2:
			for(int yy=0;yy<5;yy++)
			{
				if(!localXPtrs[yy])	continue;
				const Real* xPtr=&localXPtrs[yy][i];
				const double* mValues=lapValues[yy];
				temp+=xPtr[0]*mValues[0]+xPtr[1]*mValues[1]+xPtr[2]*mValues[2]+xPtr[3]*mValues[3]+xPtr[4]*mValues[4];
			}
			return temp;
		case 1:
			for(int yy=0;yy<3;yy++)
			{
				if(!localXPtrs[yy])	continue;
				const Real* xPtr=&localXPtrs[yy][i];
				const double* mValues=lapValues[yy];
				temp+=xPtr[0]*mValues[0]+xPtr[1]*mValues[1]+xPtr[2]*mValues[2];
			}
			return temp;
		default:
			for(int yy=0;yy<=2*Degree;yy++)
			{
				if(!localXPtrs[yy])	continue;
				for(int xx=0;xx<=2*Degree;xx++)	temp+=localXPtrs[yy][xx+i]*lapValues[yy][xx];
			}
			return temp;
	}
}
template<class Real,int Type,int Degree,int Channels>
inline Real StreamingSolver<Real,Type,Degree,Channels>::GetInteriorLaplacianValue(const double lapValues[][2*Degree+1],int i)
{
	Real temp=0;
	const Real* xPtr;
	const double* mValues;
	switch(Degree)
	{
		case 2:
			xPtr=&localXPtrs[2][i];
			mValues=lapValues[2];
			temp =xPtr[0]*mValues[0]+xPtr[1]*mValues[1]+xPtr[2]*mValues[2]+xPtr[3]*mValues[3]+xPtr[4]*mValues[4];
			xPtr=&localXAccum[0][i];
			mValues=lapValues[0];
			temp+=xPtr[0]*mValues[0]+xPtr[1]*mValues[1]+xPtr[2]*mValues[2]+xPtr[3]*mValues[3]+xPtr[4]*mValues[4];
			xPtr=&localXAccum[1][i];
			mValues=lapValues[1];
			temp+=xPtr[0]*mValues[0]+xPtr[1]*mValues[1]+xPtr[2]*mValues[2]+xPtr[3]*mValues[3]+xPtr[4]*mValues[4];
			return temp;
		case 1:
			xPtr=&localXPtrs[1][i];
			mValues=lapValues[1];
			temp =xPtr[0]*mValues[0]+xPtr[1]*mValues[1]+xPtr[2]*mValues[2];
			xPtr=&localXAccum[0][i];
			mValues=lapValues[0];
			temp+=xPtr[0]*mValues[0]+xPtr[1]*mValues[1]+xPtr[2]*mValues[2];
			return temp;
		default:
			xPtr=&localXPtrs[Degree][i];
			mValues=lapValues[Degree];
			for(int xx=0;xx<=2*Degree;xx++)	temp+=xPtr[xx]*mValues[xx];
			for(int yy=0;yy<Degree;yy++)
			{
				xPtr=&localXAccum[yy][i];
				mValues=lapValues[yy];
				for(int xx=0;xx<=2*Degree;xx++)	temp+=xPtr[xx]*mValues[xx];
			}
			return temp;
	}
}
template<class Real,int Type,int Degree,int Channels>
inline Real StreamingSolver<Real,Type,Degree,Channels>::GaussSeidelUpdate(const double lapValues[][2*Degree+1],int i,int iB)
{
	return (BStream[iB+i]-        GetLaplacianValue(lapValues,i))/lapValues[Degree][Degree];
}
template<class Real,int Type,int Degree,int Channels>
inline Real StreamingSolver<Real,Type,Degree,Channels>::InteriorGaussSeidelUpdate(const double lapValues[][2*Degree+1],int i,int iB)
{
	return (BStream[iB+i]-GetInteriorLaplacianValue(lapValues,i))/lapValues[Degree][Degree];
}
#endif // !USE_SSE_CODE
template<class Real,int Type,int Degree,int Channels>
bool StreamingSolver<Real,Type,Degree,Channels>::Increment(void)
{
	int idx;
	if(clearX)
	{
#if MERGE_INTERIOR_STREAMS
		if(setResidual)	idx=index-Degree-1+iSize;
		else			idx=index-Degree  +iSize;
#else // !MERGE_INTERIOR_STREAMS
		if(setResidual)	idx=index-Degree-1+xSize;
		else			idx=index-Degree  +xSize;
#endif // MERGE_INTERIOR_STREAMS
		if(idx>=0 && idx<minor)	for(int c=0;c<Channels;c++)	memset(GetXRow(idx,c),0,sizeof(Real)*_major*RealPerWord);
	}
	if(clearB)
	{
#if MERGE_INTERIOR_STREAMS
		if(setResidual)	idx=index-1+iSize;
		else			idx=index  +iSize;
#else // !MERGE_INTERIOR_STREAMS
		if(setResidual)	idx=index-1+bSize;
		else			idx=index  +bSize;
#endif // MERGE_INTERIOR_STREAMS
		if(idx>=0 && idx<minor)	for(int c=0;c<Channels;c++)	memset(GetBRow(idx,c),0,sizeof(Real)*_major*RealPerWord);
	}
	index++;
	// Allow for the writing out of the B vector
	if(setResidual)	return (index-1<minor);
	else			return (index       <minor);
}
#if MERGE_INTERIOR_STREAMS
template<class Real,int Type,int Degree,int Channels>
bool StreamingSolver<Real,Type,Degree,Channels>::UpdateIOutput(StreamingGrid* I)
{
	// Copy the solution
	if(index>=0 && index<minor)
	{
		Real* iPtr=(Real*)(*I)[index];
		for(int c=0;c<Channels;c++)
		{
			int idx;
			idx=MyModIndex(index,iSize)*2*_major*RealPerWord*Channels+c*_major*RealPerWord;
			memcpy(&iPtr[major*c],&IStream[idx],sizeof(Real)*2*major);
		}
		return true;
	}
	return false;
}
template<class Real,int Type,int Degree,int Channels>
bool StreamingSolver<Real,Type,Degree,Channels>::UpdateIInput(StreamingGrid* I)
{
	int idx = index+_iters+Degree-1;
	if(idx>=0 && idx<minor)
	{
		Real* iPtr=(Real*)(*I)[idx];
		for(int c=0;c<Channels;c++)
		{
			int idx1=MyModIndex(index,iSize)*2*_major*RealPerWord*Channels+_major*RealPerWord*Channels+c*_major*RealPerWord;
			memcpy(&iPtr[major*c],GetXRow(index&IStream[idx1],sizeof(Real)*2*major);
		}
	}
	return false;
}
#endif // MERGE_INTERIOR_STREAMS
template<class Real,int Type,int Degree,int Channels>
template<class StorageType>
bool StreamingSolver<Real,Type,Degree,Channels>::UpdateXInput(StreamingGrid* X)
{
	int idx = index+_iters+Degree-1;
	// Read in from the X vector
	if(idx>=0 && idx<minor)
	{
		StorageType* xPtr=(StorageType*)(*X)[idx];
		for(int c=0;c<Channels;c++)
		{
			int idx2=major*c;
			Real* xRow=GetXRow(idx,c);
			for(int jj=0;jj<major;jj++)	xRow[jj]+=Real(xPtr[jj+idx2]);
		}
		return true;
	}
	return false;
}
template<class Real,int Type,int Degree,int Channels>
template<class StorageType>
bool StreamingSolver<Real,Type,Degree,Channels>::UpdateBInput(StreamingGrid* B)
{
	int idx = index+_iters-1;
	// Read in from the B vector
	if(idx>=0 && idx<minor)
	{
		StorageType* bPtr=(StorageType*)(*B)[idx];
		for(int c=0;c<Channels;c++)
		{
			int idx2=major*c;
			Real* bRow=GetBRow(idx,c);
			for(int jj=0;jj<major;jj++)	bRow[jj]=Real(bPtr[jj+idx2]);
		}
		return true;
	}
	return false;
}
template<class Real,int Type,int Degree,int Channels>
template<class StorageType>
bool StreamingSolver<Real,Type,Degree,Channels>::UpdateXOutput(StreamingGrid* X)
{
	// Copy the solution
	if(index>=0 && index<minor)
	{
		StorageType* xPtr=(StorageType*)(*X)[index];
		for(int c=0;c<Channels;c++)
		{
			int idx2=major*c;
			Real* xRow=GetXRow(index,c);
			for(int jj=0;jj<major;jj++)	xPtr[jj+idx2]=StorageType(xRow[jj]);
		}
		return true;
	}
	return false;
}
template<class Real,int Type,int Degree,int Channels>
template<class StorageType>
bool StreamingSolver<Real,Type,Degree,Channels>::UpdateBOutput(StreamingGrid* B)
{
#if STORE_RESIDUAL
	int idx = index-2;
	if(idx>=0 && idx<minor)
	{
		StorageType* bPtr=(StorageType*)(*B)[idx];
		for(int c=0;c<Channels;c++)
		{
			int idx2=major*c;
			Real* localRPtr=&localR[(idx%rSize)*_major*RealPerWord*Channels+c*_major*RealPerWord];
			for(int jj=0;jj<major;jj++)	bPtr[jj+idx2]=StorageType(localRPtr[jj]);
		}
		return true;
	}
#else // !STORE_RESIDUAL
	if(index>=0 && index<minor)
	{
		StorageType* bPtr=(StorageType*)(*B)[index];
		for(int c=0;c<Channels;c++)
		{
			int idx2=major*c;
			Real* bRow=GetBRow(index,c);
			for(int jj=0;jj<major;jj++)	bPtr[jj+idx2]=StorageType(bRow[jj]);
		}
		return true;
	}
#endif // STORE_RESIDUAL
	return false;
}
template<class Real,int Type,int Degree,int Channels>
void StreamingSolver<Real,Type,Degree,Channels>::SolveInterior(int j)
{
	int jj;
	if		(j<Degree)			jj=j;
	else if	(j>minor-1-Degree)	jj=2*Degree+(j-(minor-1));
	else						jj=Degree;
#if USE_SSE_CODE
	jj*=3;
#endif // USE_SSE_CODE
#if NEW_MISHA_CODE
#if NEW_MISHA_CODE == 1
	if(misha[0])
#elif NEW_MISHA_CODE == 2
	if(misha[0] && misha[1])
#elif NEW_MISHA_CODE == 4
	if(misha[0] && misha[1] && misha[2] && misha[3])
#endif
	{
		float* xRows[2*Degree+1];
		float zeroValues[2*Degree+1][Degree][Channels],tempValues1[2*Degree+1][Degree][Channels],tempValues2[2*Degree+1][Degree][Channels];
		for(int y=0;y<2*Degree+1;y++)	xRows[y]=GetXRow(j-Degree+y,0);
		for(int y=0;y<2*Degree+1;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	zeroValues[y][x][c]=0;
#if NEW_MISHA_CODE == 1
		misha[0]->SolveInterior(xRows,GetBRow(j,0),zeroValues,zeroValues);
#elif NEW_MISHA_CODE == 2
		for(int y=0;y<2*Degree+1;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues2[y][x][c]=GetXRow(j-Degree+y,c)[misha[0]->size()+x];
		misha[0]->SolveInterior(xRows,GetBRow(j,0),zeroValues,tempValues2);
		for(int y=0;y<2*Degree+1;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues1[y][x][c]=GetXRow(j-Degree+y,c)[misha[0]->size()-Degree+x];
		misha[1]->SolveInterior(xRows,GetBRow(j,0),tempValues1,zeroValues);
#elif NEW_MISHA_CODE == 4
		if(j&1)
		{
			for(int y=0;y<2*Degree+1;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues2[y][x][c]=GetXRow(j-Degree+y,c)[misha[0]->size()+x];
			misha[0]->SolveInterior(xRows,GetBRow(j,0),zeroValues,tempValues2);
			for(int y=0;y<2*Degree+1;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues1[y][x][c]=GetXRow(j-Degree+y,c)[misha[0]->size()-Degree+x];
			for(int y=0;y<2*Degree+1;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues2[y][x][c]=GetXRow(j-Degree+y,c)[misha[0]->size()+misha[1]->size()+x];
			misha[1]->SolveInterior(xRows,GetBRow(j,0),tempValues1,tempValues2);
			for(int y=0;y<2*Degree+1;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues1[y][x][c]=GetXRow(j-Degree+y,c)[misha[0]->size()+misha[1]->size()-Degree+x];
			for(int y=0;y<2*Degree+1;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues2[y][x][c]=GetXRow(j-Degree+y,c)[misha[0]->size()+misha[1]->size()+misha[2]->size()+x];
			misha[2]->SolveInterior(xRows,GetBRow(j,0),tempValues1,tempValues2);
			for(int y=0;y<2*Degree+1;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues1[y][x][c]=GetXRow(j-Degree+y,c)[misha[0]->size()+misha[1]->size()+misha[2]->size()-Degree+x];
			misha[3]->SolveInterior(xRows,GetBRow(j,0),tempValues1,zeroValues);
		}
		else
		{
			for(int y=0;y<2*Degree+1;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues1[y][x][c]=GetXRow(j-Degree+y,c)[misha[0]->size()+misha[1]->size()+misha[2]->size()-Degree+x];
			misha[3]->SolveInterior(xRows,GetBRow(j,0),tempValues1,zeroValues);
			for(int y=0;y<2*Degree+1;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues1[y][x][c]=GetXRow(j-Degree+y,c)[misha[0]->size()+misha[1]->size()-Degree+x];
			for(int y=0;y<2*Degree+1;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues2[y][x][c]=GetXRow(j-Degree+y,c)[misha[0]->size()+misha[1]->size()+misha[2]->size()+x];
			misha[2]->SolveInterior(xRows,GetBRow(j,0),tempValues1,tempValues2);
			for(int y=0;y<2*Degree+1;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues1[y][x][c]=GetXRow(j-Degree+y,c)[misha[0]->size()-Degree+x];
			for(int y=0;y<2*Degree+1;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues2[y][x][c]=GetXRow(j-Degree+y,c)[misha[0]->size()+misha[1]->size()+x];
			misha[1]->SolveInterior(xRows,GetBRow(j,0),tempValues1,tempValues2);
			for(int y=0;y<2*Degree+1;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues2[y][x][c]=GetXRow(j-Degree+y,c)[misha[0]->size()+x];
			misha[0]->SolveInterior(xRows,GetBRow(j,0),zeroValues,tempValues2);
		}
#endif
	}
	else
#endif // NEW_MISHA_CODE
	for(int c=0;c<Channels;c++)
	{
		localBPtr=GetBRow(j,c);
#if USE_SSE_CODE
		for(int yy=0;yy<=2*Degree;yy++)	localXPtrs[yy]=(WordClass*)GetXRow(j-Degree+yy,c);
		{
			const WordClass* xPtrs[]={localXPtrs[0],localXPtrs[1],localXPtrs[3],localXPtrs[4]};
			WordClass* aPtrs[]={localXAccum[0],localXAccum[1]};
			for(int i=0;i<_major;i++)
			{
				aPtrs[0][i]=_mm_add_ps(xPtrs[0][i],xPtrs[3][i]);
				aPtrs[1][i]=_mm_add_ps(xPtrs[1][i],xPtrs[2][i]);
			}
		}
		float dotSum;
		WordClass dSum;

		__declspec (align(ALIGNMENT)) float scratch[4];
		Real* lX=GetXRow(j,c);
		const WordClass* xPtrs[]={localXAccum[0],localXAccum[1],localXPtrs[2]};

		{
			dotSum=0;
			lX[0]=InteriorGaussSeidelUpdate0(lapTemplates3[jj].matrixValues[0],xPtrs,localBPtr,dotSum,0)*lapTemplates3[jj].diagonal[0];
			WordClass mValues[]=
			{
				lapTemplates3[jj+1].matrixValues[0][0],
				lapTemplates3[jj+1].matrixValues[0][1],
				lapTemplates3[jj+1].matrixValues[0][2]
			};
			int bound=major-4;
			for(int i=4;i<bound;i+=4)	((float*)xPtrs[2])[i]=InteriorGaussSeidelUpdate0(mValues,xPtrs,localBPtr,dotSum,i);
			if(major>=8)	lX[major-4]=InteriorGaussSeidelUpdate0(lapTemplates3[jj+2].matrixValues[0],xPtrs,localBPtr,dotSum,major-4)*lapTemplates3[jj+2].diagonal[0];
		}
		{
			dotSum=0;
			lX[1]=InteriorGaussSeidelUpdate1(lapTemplates3[jj].matrixValues[1],xPtrs,localBPtr,dotSum,1)*lapTemplates3[jj].diagonal[1];
			WordClass mValues[]=
			{
				lapTemplates3[jj+1].matrixValues[1][0],
				lapTemplates3[jj+1].matrixValues[1][1],
				lapTemplates3[jj+1].matrixValues[1][2]
			};
			int bound=major-4;
			for(int i=5;i<bound;i+=4)	((float*)xPtrs[2])[i]=InteriorGaussSeidelUpdate1(mValues,xPtrs,localBPtr,dotSum,i);
			if(major>=8)	lX[major-3]=InteriorGaussSeidelUpdate1(lapTemplates3[jj+2].matrixValues[1],xPtrs,localBPtr,dotSum,major-3)*lapTemplates3[jj+2].diagonal[1];
		}
#if 1
		{
			SetInteriorDotSum(lapTemplates3[jj].matrixValues[2],xPtrs,0,dSum);
			_mm_store_ps(scratch,dSum);
			dotSum=scratch[0]+scratch[1]+scratch[2]+scratch[3];
			WordClass mValues[]=
			{
				lapTemplates3[jj+1].matrixValues[2][0],
				lapTemplates3[jj+1].matrixValues[2][1],
				lapTemplates3[jj+1].matrixValues[2][2]
			};
			if(major<12)	lX[2]=InteriorGaussSeidelUpdate2(lapTemplates3[jj+2].matrixValues[2],xPtrs,localBPtr,dotSum,2)*lapTemplates3[jj].diagonal[2];
			else			lX[2]=InteriorGaussSeidelUpdate2(mValues,xPtrs,localBPtr,dotSum,2)*lapTemplates3[jj].diagonal[2];
			int bound=major-8;
			for(int i=6;i<bound;i+=4)	((float*)xPtrs[2])[i]=InteriorGaussSeidelUpdate2(mValues,xPtrs,localBPtr,dotSum,i);
			if(major>=12)	lX[major-6]=InteriorGaussSeidelUpdate2(lapTemplates3[jj+2].matrixValues[2],xPtrs,localBPtr,dotSum,major-6);
			lX[major-2]=(localBPtr[major-2]-dotSum)*lapTemplates3[jj+2].diagonal[2];
		}
		{
			SetInteriorDotSum(lapTemplates3[jj].matrixValues[3],xPtrs,0,dSum);
			_mm_store_ps(scratch,dSum);
			dotSum=scratch[1]+scratch[2]+scratch[3];
			WordClass mValues[]=
			{
				lapTemplates3[jj+1].matrixValues[3][0],
				lapTemplates3[jj+1].matrixValues[3][1],
				lapTemplates3[jj+1].matrixValues[3][2]
			};
			if(major<12)	lX[3]=InteriorGaussSeidelUpdate3(lapTemplates3[jj+2].matrixValues[3],xPtrs,localBPtr,dotSum,3)*lapTemplates3[jj].diagonal[3];
			else			lX[3]=InteriorGaussSeidelUpdate3(mValues,xPtrs,localBPtr,dotSum,3)*lapTemplates3[jj].diagonal[3];
			int bound=major-8;
			for(int i=7;i<bound;i+=4)	((float*)xPtrs[2])[i]=InteriorGaussSeidelUpdate3(mValues,xPtrs,localBPtr,dotSum,i);
			if(major>=12)	lX[major-5]=InteriorGaussSeidelUpdate3(lapTemplates3[jj+2].matrixValues[3],xPtrs,localBPtr,dotSum,major-5);
			lX[major-1]=(localBPtr[major-1]-dotSum)*lapTemplates3[jj+2].diagonal[3];
		}
#endif
#else // !USE_SSE_CODE
		for(int yy=0;yy<=2*Degree;yy++)	localXPtrs[yy]=GetXRow(j-Degree+yy,c)-Degree;
		for(int d=0;d<Degree;d++)	for(int i=0;i<major+2*Degree;i++)	localXAccum[d][i]=localXPtrs[d][i]+localXPtrs[2*Degree-d][i];
		int idx1=MyModIndex(j,xSize)*_major*RealPerWord*Channels+c*_major*RealPerWord;
		int idx2=MyModIndex(j,bSize)*_major*RealPerWord*Channels+c*_major*RealPerWord;
		for(int i=0;i<major;i++)
		{
			int ii;
			if		(i<Degree)			ii=i;
			else if	(i>major-1-Degree)	ii=2*Degree+(i-(major-1));
			else						ii=Degree;
			XStream[idx1+i]+=InteriorGaussSeidelUpdate(laplacianStencil.caseTable[jj][ii].values,i,idx2);
		}
#endif // USE_SSE_CODE
	}
}
template<class Real,int Type,int Degree,int Channels>
void StreamingSolver<Real,Type,Degree,Channels>::SetInteriorResidual(int j)
{
	int jj;
	if		(j<Degree)			jj=j;
	else if	(j>minor-1-Degree)	jj=2*Degree+(j-(minor-1));
	else						jj=Degree;
#if USE_SSE_CODE
	jj*=3;
#endif // USE_SSE_CODE
	for(int c=0;c<Channels;c++)
	{
		localBPtr=GetBRow(j,c);
		localRPtr=&localR[(j%rSize)*_major*RealPerWord*Channels+c*_major*RealPerWord];
#if USE_SSE_CODE
		for(int yy=0;yy<=2*Degree;yy++)	localXPtrs[yy]=(WordClass*)GetXRow(j-Degree+yy,c);
		{
			const WordClass* xPtrs[]={localXPtrs[0],localXPtrs[1],localXPtrs[3],localXPtrs[4]};
			WordClass* aPtrs[]={localXAccum[0],localXAccum[1]};
			for(int i=0;i<_major;i++)
			{
				aPtrs[0][i]=_mm_add_ps(xPtrs[0][i],xPtrs[3][i]);
				aPtrs[1][i]=_mm_add_ps(xPtrs[1][i],xPtrs[2][i]);
			}
		}
		float dotSum;
		WordClass dSum;

		__declspec (align(ALIGNMENT)) float scratch[4];
		const WordClass* xPtrs[]={localXAccum[0],localXAccum[1],localXPtrs[2]};
		{
			dotSum=0;
			localRPtr[0]=localBPtr[0]-GetInteriorLaplacianValue0(lapTemplates2[jj].matrixValues[0],xPtrs,dotSum,0);
			WordClass mValues[]=
			{
				lapTemplates2[jj+1].matrixValues[0][0],
				lapTemplates2[jj+1].matrixValues[0][1],
				lapTemplates2[jj+1].matrixValues[0][2]
			};
			int bound=major-4;
			for(int i=4;i<bound;i+=4)	localRPtr[i]=localBPtr[i]-GetInteriorLaplacianValue0(mValues,xPtrs,dotSum,i);
			if(major>=8)	localRPtr[major-4]=localBPtr[major-4]-GetInteriorLaplacianValue0(lapTemplates2[jj+2].matrixValues[0],xPtrs,dotSum,major-4);
		}
		{
			dotSum=0;
			localRPtr[1]=localBPtr[1]-GetInteriorLaplacianValue1(lapTemplates2[jj].matrixValues[1],xPtrs,dotSum,1);
			WordClass mValues[]=
			{
				lapTemplates2[jj+1].matrixValues[1][0],
				lapTemplates2[jj+1].matrixValues[1][1],
				lapTemplates2[jj+1].matrixValues[1][2]

			};
			int bound=major-4;
			for(int i=5;i<bound;i+=4)	localRPtr[i]=localBPtr[i]-GetInteriorLaplacianValue1(mValues,xPtrs,dotSum,i);
			if(major>=8)	localRPtr[major-3]=localBPtr[major-3]-GetInteriorLaplacianValue1(lapTemplates2[jj+2].matrixValues[1],xPtrs,dotSum,major-3);
		}
		{
			SetInteriorDotSum(lapTemplates2[jj].matrixValues[2],xPtrs,0,dSum);
			_mm_store_ps(scratch,dSum);
			dotSum=scratch[0]+scratch[1]+scratch[2]+scratch[3];
			WordClass mValues[]=
			{
				lapTemplates2[jj+1].matrixValues[2][0],
				lapTemplates2[jj+1].matrixValues[2][1],
				lapTemplates2[jj+1].matrixValues[2][2]
			};
			if(major<12)	localRPtr[2]=localBPtr[2]-GetInteriorLaplacianValue2(lapTemplates2[jj+2].matrixValues[2],xPtrs,dotSum,2);
			else			localRPtr[2]=localBPtr[2]-GetInteriorLaplacianValue2(mValues,xPtrs,dotSum,2);
			int bound=major-8;
			for(int i=6;i<bound;i+=4)	localRPtr[i]=localBPtr[i]-GetInteriorLaplacianValue2(mValues,xPtrs,dotSum,i);
			if(major>=12)	localRPtr[major-6]=localBPtr[major-6]-GetInteriorLaplacianValue2(lapTemplates2[jj+2].matrixValues[2],xPtrs,dotSum,major-6);
			localRPtr[major-2]=localBPtr[major-2]-dotSum;
		}
		{
			SetInteriorDotSum(lapTemplates2[jj].matrixValues[3],xPtrs,0,dSum);
			_mm_store_ps(scratch,dSum);
			dotSum=scratch[1]+scratch[2]+scratch[3];
			WordClass mValues[]=
			{
				lapTemplates2[jj+1].matrixValues[3][0],
				lapTemplates2[jj+1].matrixValues[3][1],
				lapTemplates2[jj+1].matrixValues[3][2]
			};
			if(major<12)	localRPtr[3]=localBPtr[3]-GetInteriorLaplacianValue3(lapTemplates3[jj+2].matrixValues[3],xPtrs,dotSum,3);
			else			localRPtr[3]=localBPtr[3]-GetInteriorLaplacianValue3(mValues,xPtrs,dotSum,3);
			int bound=major-8;
			for(int i=7;i<bound;i+=4)	localRPtr[i]=localBPtr[i]-GetInteriorLaplacianValue3(mValues,xPtrs,dotSum,i);
			if(major>=12)	localRPtr[major-5]=localBPtr[major-5]-GetInteriorLaplacianValue3(lapTemplates2[jj+2].matrixValues[3],xPtrs,dotSum,major-5);
			localRPtr[major-1]=localBPtr[major-1]-dotSum;
		}
#else // !USE_SSE_CODE
		for(int yy=0;yy<=2*Degree;yy++)	localXPtrs[yy]=GetXRow(j-Degree+yy,c)-Degree;
		for(int d=0;d<Degree;d++)	for(int i=0;i<major+2*Degree;i++)	localXAccum[d][i]=localXPtrs[d][i]+localXPtrs[2*Degree-d][i];
		for(int i=0;i<major;i++)
		{
			int ii;
			if		(i<Degree)			ii=i;
			else if	(i>major-1-Degree)	ii=2*Degree+(i-(major-1));
			else						ii=Degree;
			localRPtr[i]=localBPtr[i]-GetInteriorLaplacianValue(laplacianStencil.caseTable[jj][ii].values,i);
		}
#endif // USE_SSE_CODE
	}
}

template<class Real,int Type,int Degree,int Channels>
void StreamingSolver<Real,Type,Degree,Channels>::Solve(int j)
{
#if NEW_MISHA_CODE
#if NEW_MISHA_CODE == 1
	if(misha[0])
#elif NEW_MISHA_CODE == 2
	if(misha[0] && misha[1])
#elif NEW_MISHA_CODE == 4
	if(misha[0] && misha[1] && misha[2] && misha[3])
#endif
	{
		float* xRows[2*Degree+1];
		float zeroValues[2*Degree+1][Degree][Channels],tempValues1[2*Degree+1][Degree][Channels],tempValues2[2*Degree+1][Degree][Channels];
		for(int y=0;y<2*Degree+1;y++)
			if(j-Degree+y>=0 && j-Degree+y<minor)	xRows[y]=GetXRow(j-Degree+y,0);
			else									xRows[y]=NULL;
		for(int y=0;y<2*Degree+1;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	zeroValues[y][x][c]=tempValues1[y][x][c]=tempValues2[y][x][c]=0;

#if NEW_MISHA_CODE == 1
			misha[0]->Solve(j,minor,xRows,GetBRow(j,0),zeroValues,zeroValues);
#elif NEW_MISHA_CODE == 2
			for(int y=0;y<2*Degree+1;y++)	if(j-Degree+y>=0 && j-Degree+y<minor)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues2[y][x][c]=GetXRow(j-Degree+y,c)[misha[0]->size()+x];
			misha[0]->Solve(j,minor,xRows,GetBRow(j,0),zeroValues,tempValues2);
			for(int y=0;y<2*Degree+1;y++)	if(j-Degree+y>=0 && j-Degree+y<minor)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues1[y][x][c]=GetXRow(j-Degree+y,c)[misha[0]->size()-Degree+x];
			misha[1]->Solve(j,minor,xRows,GetBRow(j,0),tempValues1,zeroValues);
#elif NEW_MISHA_CODE == 4
		if(j&1)
		{
			for(int y=0;y<2*Degree+1;y++)	if(j-Degree+y>=0 && j-Degree+y<minor)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues2[y][x][c]=GetXRow(j-Degree+y,c)[misha[0]->size()+x];
			misha[0]->Solve(j,minor,xRows,GetBRow(j,0),zeroValues,tempValues2);
			for(int y=0;y<2*Degree+1;y++)	if(j-Degree+y>=0 && j-Degree+y<minor)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues1[y][x][c]=GetXRow(j-Degree+y,c)[misha[0]->size()-Degree+x];
			for(int y=0;y<2*Degree+1;y++)	if(j-Degree+y>=0 && j-Degree+y<minor)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues2[y][x][c]=GetXRow(j-Degree+y,c)[misha[0]->size()+misha[1]->size()+x];
			misha[1]->Solve(j,minor,xRows,GetBRow(j,0),tempValues1,tempValues2);
			for(int y=0;y<2*Degree+1;y++)	if(j-Degree+y>=0 && j-Degree+y<minor)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues1[y][x][c]=GetXRow(j-Degree+y,c)[misha[0]->size()+misha[1]->size()-Degree+x];
			for(int y=0;y<2*Degree+1;y++)	if(j-Degree+y>=0 && j-Degree+y<minor)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues2[y][x][c]=GetXRow(j-Degree+y,c)[misha[0]->size()+misha[1]->size()+misha[2]->size()+x];
			misha[2]->Solve(j,minor,xRows,GetBRow(j,0),tempValues1,tempValues2);
			for(int y=0;y<2*Degree+1;y++)	if(j-Degree+y>=0 && j-Degree+y<minor)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues1[y][x][c]=GetXRow(j-Degree+y,c)[misha[0]->size()+misha[1]->size()+misha[2]->size()-Degree+x];
			misha[3]->Solve(j,minor,xRows,GetBRow(j,0),tempValues1,zeroValues);
		}
		else
		{
			for(int y=0;y<2*Degree+1;y++)	if(j-Degree+y>=0 && j-Degree+y<minor)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues1[y][x][c]=GetXRow(j-Degree+y,c)[misha[0]->size()+misha[1]->size()+misha[2]->size()-Degree+x];
			misha[3]->Solve(j,minor,xRows,GetBRow(j,0),tempValues1,zeroValues);
			for(int y=0;y<2*Degree+1;y++)	if(j-Degree+y>=0 && j-Degree+y<minor)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues1[y][x][c]=GetXRow(j-Degree+y,c)[misha[0]->size()+misha[1]->size()-Degree+x];
			for(int y=0;y<2*Degree+1;y++)	if(j-Degree+y>=0 && j-Degree+y<minor)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues2[y][x][c]=GetXRow(j-Degree+y,c)[misha[0]->size()+misha[1]->size()+misha[2]->size()+x];
			misha[2]->Solve(j,minor,xRows,GetBRow(j,0),tempValues1,tempValues2);
			for(int y=0;y<2*Degree+1;y++)	if(j-Degree+y>=0 && j-Degree+y<minor)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues1[y][x][c]=GetXRow(j-Degree+y,c)[misha[0]->size()-Degree+x];
			for(int y=0;y<2*Degree+1;y++)	if(j-Degree+y>=0 && j-Degree+y<minor)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues2[y][x][c]=GetXRow(j-Degree+y,c)[misha[0]->size()+misha[1]->size()+x];
			misha[1]->Solve(j,minor,xRows,GetBRow(j,0),tempValues1,tempValues2);
			for(int y=0;y<2*Degree+1;y++)	if(j-Degree+y>=0 && j-Degree+y<minor)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues2[y][x][c]=GetXRow(j-Degree+y,c)[misha[0]->size()+x];
			misha[0]->Solve(j,minor,xRows,GetBRow(j,0),zeroValues,tempValues2);
		}
#endif
	}
	else
#endif // NEW_MISHA_CODE
	for(int c=0;c<Channels;c++)
	{
#if USE_SSE_CODE
		for(int yy=0;yy<=2*Degree;yy++)
		{
			if(j-Degree+yy>=0 && j-Degree+yy<minor)	localXPtrs[yy]=(WordClass*)GetXRow(j-Degree+yy,c);
			else									localXPtrs[yy]=NULL;
		}
#else // !USE_SSE_CODE
		for(int yy=0;yy<=2*Degree;yy++)
		{
			if(j-Degree+yy>=0 && j-Degree+yy<minor)	localXPtrs[yy]=GetXRow(j-Degree+yy,c)-Degree;
			else									localXPtrs[yy]=NULL;
		}
#endif // USE_SSE_CODE
		int jj;
		if		(j<Degree)			jj=j;
		else if	(j>minor-1-Degree)	jj=2*Degree+(j-(minor-1));
		else						jj=Degree;
#if USE_SSE_CODE
		localBPtr=GetBRow(j,c);
		jj*=3;
		float dotSum;
		WordClass dSum;

		__declspec (align(ALIGNMENT)) float scratch[4];
		Real* lX=GetXRow(j,c);
		const WordClass* xPtrs[]={localXPtrs[0],localXPtrs[1],localXPtrs[2],localXPtrs[3],localXPtrs[4]};
		{
			dotSum=0;
			lX[0]=GaussSeidelUpdate0(lapTemplates3[jj].matrixValues[0],xPtrs,localBPtr,dotSum,0)*lapTemplates2[jj].diagonal[0];
			float d=lapTemplates2[jj+1].diagonal[0];
			WordClass mValues[]=
			{
				lapTemplates3[jj+1].matrixValues[0][0],
				lapTemplates3[jj+1].matrixValues[0][1],
				lapTemplates3[jj+1].matrixValues[0][2],
				lapTemplates3[jj+1].matrixValues[0][3],
				lapTemplates3[jj+1].matrixValues[0][4]
			};
			int bound=major-4;
			for(int i=4;i<bound;i+=4)	lX[i]=GaussSeidelUpdate0(mValues,xPtrs,localBPtr,dotSum,i)*d;
			if(major>=8)	lX[major-4]=GaussSeidelUpdate0(lapTemplates3[jj+2].matrixValues[0],xPtrs,localBPtr,dotSum,major-4)*lapTemplates3[jj+2].diagonal[0];
		}
		{
			dotSum=0;
			lX[1]=GaussSeidelUpdate1(lapTemplates3[jj].matrixValues[1],xPtrs,localBPtr,dotSum,1)*lapTemplates3[jj].diagonal[1];
			WordClass mValues[]=
			{
				lapTemplates3[jj+1].matrixValues[1][0],
				lapTemplates3[jj+1].matrixValues[1][1],
				lapTemplates3[jj+1].matrixValues[1][2],
				lapTemplates3[jj+1].matrixValues[1][3],
				lapTemplates3[jj+1].matrixValues[1][4]
			};
			float d=lapTemplates2[jj+1].diagonal[1];
			int bound=major-4;
			for(int i=5;i<bound;i+=4)	lX[i]=GaussSeidelUpdate1(mValues,xPtrs,localBPtr,dotSum,i)*d;
			if(major>=8)	lX[major-3]=GaussSeidelUpdate1(lapTemplates3[jj+2].matrixValues[1],xPtrs,localBPtr,dotSum,major-3)*lapTemplates3[jj+2].diagonal[1];
		}
		{
			SetDotSum(lapTemplates3[jj].matrixValues[2],xPtrs,0,dSum);
			_mm_store_ps(scratch,dSum);
			dotSum=scratch[0]+scratch[1]+scratch[2]+scratch[3];
			WordClass mValues[]=
			{
				lapTemplates3[jj+1].matrixValues[2][0],
				lapTemplates3[jj+1].matrixValues[2][1],
				lapTemplates3[jj+1].matrixValues[2][2],
				lapTemplates3[jj+1].matrixValues[2][3],
				lapTemplates3[jj+1].matrixValues[2][4]
			};
			float d=lapTemplates2[jj+1].diagonal[2];
			if(major<12)	lX[2]=GaussSeidelUpdate2(lapTemplates3[jj+2].matrixValues[2],xPtrs,localBPtr,dotSum,2)*lapTemplates3[jj].diagonal[2];
			else			lX[2]=GaussSeidelUpdate2(mValues,xPtrs,localBPtr,dotSum,2)*lapTemplates3[jj].diagonal[2];
			int bound=major-8;
			for(int i=6;i<bound;i+=4)	lX[i]=GaussSeidelUpdate2(mValues,xPtrs,localBPtr,dotSum,i)*d;
			if(major>=12)	lX[major-6]=GaussSeidelUpdate2(lapTemplates3[jj+2].matrixValues[2],xPtrs,localBPtr,dotSum,major-6)*lapTemplates3[jj+1].diagonal[2];
			lX[major-2]=(localBPtr[major-2]-dotSum)*lapTemplates3[jj+2].diagonal[2];
		}
		{
			SetDotSum(lapTemplates3[jj].matrixValues[3],xPtrs,0,dSum);
			_mm_store_ps(scratch,dSum);
			dotSum=scratch[1]+scratch[2]+scratch[3];
			WordClass mValues[]=
			{
				lapTemplates3[jj+1].matrixValues[3][0],
				lapTemplates3[jj+1].matrixValues[3][1],
				lapTemplates3[jj+1].matrixValues[3][2],
				lapTemplates3[jj+1].matrixValues[3][3],
				lapTemplates3[jj+1].matrixValues[3][4]
			};
			float d=lapTemplates2[jj+1].diagonal[3];
			if(major<12)	lX[3]=GaussSeidelUpdate3(lapTemplates3[jj+2].matrixValues[3],xPtrs,localBPtr,dotSum,3)*lapTemplates3[jj].diagonal[3];
			else			lX[3]=GaussSeidelUpdate3(mValues,xPtrs,localBPtr,dotSum,3)*lapTemplates3[jj].diagonal[3];
			int bound=major-8;
			for(int i=7;i<bound;i+=4)	lX[i]=GaussSeidelUpdate3(mValues,xPtrs,localBPtr,dotSum,i)*d;
			if(major>=12)	lX[major-5]=GaussSeidelUpdate3(lapTemplates3[jj+2].matrixValues[3],xPtrs,localBPtr,dotSum,major-5)*lapTemplates3[jj+1].diagonal[3];
			lX[major-1]=(localBPtr[major-1]-dotSum)*lapTemplates3[jj+2].diagonal[3];
		}
#else // !USE_SSE_CODE
		Real* lX=GetXRow(j,c);
		int idx1=MyModIndex(j,xSize)*_major*RealPerWord*Channels+c*_major*RealPerWord;
		int idx2=MyModIndex(j,bSize)*_major*RealPerWord*Channels+c*_major*RealPerWord;
		for(int i=0;i<major;i++)
		{
			int ii;
			if		(i<Degree)			ii=i;
			else if	(i>major-1-Degree)	ii=2*Degree+(i-(major-1));
			else						ii=Degree;
			lX[i]+=GaussSeidelUpdate(laplacianStencil.caseTable[jj][ii].values,i,idx2);
		}
#endif USE_SSE_CODE
	}
}

template<class Real,int Type,int Degree,int Channels>
void StreamingSolver<Real,Type,Degree,Channels>::SetResidual(int j)
{
	for(int c=0;c<Channels;c++)
	{
#if USE_SSE_CODE
		for(int yy=0;yy<=2*Degree;yy++)
			if(j-Degree+yy>=0 && j-Degree+yy<minor)	localXPtrs[yy]=(WordClass*)GetXRow(j-Degree+yy,c);
			else									localXPtrs[yy]=NULL;
#else // !USE_SSE_CODE
		for(int yy=0;yy<=2*Degree;yy++)
		{
			if(j-Degree+yy>=0 && j-Degree+yy<minor)	localXPtrs[yy]=GetXRow(j-Degree+yy,c)-Degree;
			else									localXPtrs[yy]=NULL;
		}
#endif // USE_SSE_CODE
		int jj;
		if		(j<Degree)			jj=j;
		else if	(j>minor-1-Degree)	jj=2*Degree+(j-(minor-1));
		else						jj=Degree;
		localBPtr=GetBRow(j,c);
		localRPtr=&localR[(j%rSize)*_major*RealPerWord*Channels+c*_major*RealPerWord];
#if USE_SSE_CODE
		jj*=3;
		float dotSum;
		WordClass dSum;

		__declspec (align(ALIGNMENT)) float scratch[4];
		const WordClass* xPtrs[]={localXPtrs[0],localXPtrs[1],localXPtrs[2],localXPtrs[3],localXPtrs[4]};
		{
			dotSum=0;
			localRPtr[0]=localBPtr[0]-GetLaplacianValue0(lapTemplates2[jj].matrixValues[0],xPtrs,dotSum,0);
			WordClass mValues[]=
			{
				lapTemplates2[jj+1].matrixValues[0][0],
				lapTemplates2[jj+1].matrixValues[0][1],
				lapTemplates2[jj+1].matrixValues[0][2],
				lapTemplates2[jj+1].matrixValues[0][3],
				lapTemplates2[jj+1].matrixValues[0][4]
			};
			int bound=major-4;
			for(int i=4;i<bound;i+=4)	localRPtr[i]=localBPtr[i]-GetLaplacianValue0(mValues,xPtrs,dotSum,i);
			if(major>=8)				localRPtr[bound]=localBPtr[bound]-GetLaplacianValue0(lapTemplates2[jj+2].matrixValues[0],xPtrs,dotSum,bound);
		}
		{
			dotSum=0;
			localRPtr[1]=localBPtr[1]-GetLaplacianValue1(lapTemplates2[jj].matrixValues[1],xPtrs,dotSum,1);
			WordClass mValues[]=
			{
				lapTemplates2[jj+1].matrixValues[1][0],
				lapTemplates2[jj+1].matrixValues[1][1],
				lapTemplates2[jj+1].matrixValues[1][2],
				lapTemplates2[jj+1].matrixValues[1][3],
				lapTemplates2[jj+1].matrixValues[1][4]
			};
			int bound=major-4;
			for(int i=5;i<bound;i+=4)	localRPtr[i]=localBPtr[i]-GetLaplacianValue1(mValues,xPtrs,dotSum,i);
			if(major>=8)	localRPtr[major-3]=localBPtr[major-3]-GetLaplacianValue1(lapTemplates2[jj+2].matrixValues[1],xPtrs,dotSum,major-3);
		}
		{
			SetDotSum(lapTemplates2[jj].matrixValues[2],xPtrs,0,dSum);
			_mm_store_ps(scratch,dSum);
			dotSum=scratch[0]+scratch[1]+scratch[2]+scratch[3];
			WordClass mValues[]=
			{
				lapTemplates2[jj+1].matrixValues[2][0],
				lapTemplates2[jj+1].matrixValues[2][1],
				lapTemplates2[jj+1].matrixValues[2][2],
				lapTemplates2[jj+1].matrixValues[2][3],
				lapTemplates2[jj+1].matrixValues[2][4]
			};
			float d=lapTemplates2[jj+1].diagonal[2];
			if(major<12)	localRPtr[2]=localBPtr[2]-GetLaplacianValue2(lapTemplates2[jj+2].matrixValues[2],xPtrs,dotSum,2);
			else			localRPtr[2]=localBPtr[2]-GetLaplacianValue2(mValues,xPtrs,dotSum,2);
			int bound=major-8;
			for(int i=6;i<bound;i+=4)	localRPtr[i]=localBPtr[i]-GetLaplacianValue2(mValues,xPtrs,dotSum,i);
			if(major>=12)	localRPtr[major-6]=localBPtr[major-6]-GetLaplacianValue2(lapTemplates2[jj+2].matrixValues[2],xPtrs,dotSum,major-6);
			localRPtr[major-2]=localBPtr[major-2]-dotSum;
		}
		{
			SetDotSum(lapTemplates2[jj].matrixValues[3],xPtrs,0,dSum);
			_mm_store_ps(scratch,dSum);
			dotSum=scratch[1]+scratch[2]+scratch[3];
			WordClass mValues[]=
			{
				lapTemplates2[jj+1].matrixValues[3][0],
				lapTemplates2[jj+1].matrixValues[3][1],
				lapTemplates2[jj+1].matrixValues[3][2],
				lapTemplates2[jj+1].matrixValues[3][3],
				lapTemplates2[jj+1].matrixValues[3][4]
			};
			if(major<12)	localRPtr[3]=localBPtr[3]-GetLaplacianValue3(lapTemplates2[jj+2].matrixValues[3],xPtrs,dotSum,3);
			else			localRPtr[3]=localBPtr[3]-GetLaplacianValue3(mValues,xPtrs,dotSum,3);
			int bound=major-8;
			for(int i=7;i<bound;i+=4)	localRPtr[i]=localBPtr[i]-GetLaplacianValue3(mValues,xPtrs,dotSum,i);
			if(major>=12)	localRPtr[major-5]=localBPtr[major-5]-GetLaplacianValue3(lapTemplates2[jj+2].matrixValues[3],xPtrs,dotSum,major-5);
			localRPtr[major-1]=(localBPtr[major-1]-dotSum);
		}
#else // !USE_SSE_CODE
		for(int i=0;i<major;i++)
		{
			int ii;
			if		(i<Degree)			ii=i;
			else if	(i>major-1-Degree)	ii=2*Degree+(i-(major-1));
			else						ii=Degree;
			localRPtr[i]=localBPtr[i]-GetLaplacianValue(laplacianStencil.caseTable[jj][ii].values,i);
		}
#endif // USE_SSE_CODE
	}
}

/////////////////////
template<class Real,int Type,int Degree,int Channels>
void StreamingSolver<Real,Type,Degree,Channels>::Solve(void)
{
	int idx = index+_iters-1;

	// Solve the linear system
	bool skip=false;
	int start,stop;
	if		(index<0)		start=0;
	else if	(index>=minor)	start=minor;
	else					start=index;
	if	(index+_iters-1<0)	stop=0;
	else					stop=index+_iters;
	if(start>=stop)	skip=true;
	if(!skip)
		for(int i=stop-1;i>=start;i-=Degree)	// Iterate over the minor index
		{
			if(i<0 || i>=minor)	continue;
			if(i>=Degree && i<minor-Degree)	SolveInterior(i);
			else							Solve(i);
		}
	idx=index-1;
	if(idx>=0 && idx<minor && setResidual)
	{
		if(idx>=Degree && idx<minor-Degree)
#if NEW_MISHA_CODE
#if NEW_MISHA_CODE == 1
			if(misha[0])
#elif NEW_MISHA_CODE == 2
			if(misha[0] && misha[1])
#elif NEW_MISHA_CODE == 4
			if(misha[0] && misha[1] && misha[2] && misha[3])
#endif
			{
				float* rRow=&localR[(idx%rSize)*_major*RealPerWord*Channels];
				const float* xRows[2*Degree+1];
				float zeroValues[2*Degree+1][Degree][Channels],tempValues1[2*Degree+1][Degree][Channels],tempValues2[2*Degree+1][Degree][Channels];
				for(int y=0;y<2*Degree+1;y++)	xRows[y]=GetXRow(idx-Degree+y,0);
				for(int y=0;y<2*Degree+1;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	zeroValues[y][x][c]=0;
#if NEW_MISHA_CODE == 1
				misha[0]->SetInteriorResidual(rRow,xRows,GetBRow(idx,0),zeroValues,zeroValues);
#elif NEW_MISHA_CODE == 2
				for(int y=0;y<2*Degree+1;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues2[y][x][c]=GetXRow(idx-Degree+y,c)[misha[0]->size()+x];
				misha[0]->SetInteriorResidual(rRow,xRows,GetBRow(idx,0),zeroValues,tempValues2);
				for(int y=0;y<2*Degree+1;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues1[y][x][c]=GetXRow(idx-Degree+y,c)[misha[0]->size()-Degree+x];
				misha[1]->SetInteriorResidual(rRow,xRows,GetBRow(idx,0),tempValues1,zeroValues);
#elif NEW_MISHA_CODE == 4
				for(int y=0;y<2*Degree+1;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues2[y][x][c]=GetXRow(idx-Degree+y,c)[misha[0]->size()+x];
				misha[0]->SetInteriorResidual(rRow,xRows,GetBRow(idx,0),zeroValues,tempValues2);
				for(int y=0;y<2*Degree+1;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues1[y][x][c]=GetXRow(idx-Degree+y,c)[misha[0]->size()-Degree+x];
				for(int y=0;y<2*Degree+1;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues2[y][x][c]=GetXRow(idx-Degree+y,c)[misha[0]->size()+misha[1]->size()+x];
				misha[1]->SetInteriorResidual(rRow,xRows,GetBRow(idx,0),tempValues1,tempValues2);
				for(int y=0;y<2*Degree+1;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues1[y][x][c]=GetXRow(idx-Degree+y,c)[misha[0]->size()+misha[1]->size()-Degree+x];
				for(int y=0;y<2*Degree+1;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues2[y][x][c]=GetXRow(idx-Degree+y,c)[misha[0]->size()+misha[1]->size()+misha[2]->size()+x];
				misha[2]->SetInteriorResidual(rRow,xRows,GetBRow(idx,0),tempValues1,tempValues2);
				for(int y=0;y<2*Degree+1;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues1[y][x][c]=GetXRow(idx-Degree+y,c)[misha[0]->size()+misha[1]->size()+misha[2]->size()-Degree+x];
				misha[3]->SetInteriorResidual(rRow,xRows,GetBRow(idx,0),tempValues1,zeroValues);
#endif
			}
			else
#endif // NEW_MISHA_CODE
			SetInteriorResidual(idx);
		else
#if NEW_MISHA_CODE
#if NEW_MISHA_CODE == 1
			if(misha[0])
#elif NEW_MISHA_CODE == 2
			if(misha[0] && misha[1])
#elif NEW_MISHA_CODE == 4
			if(misha[0] && misha[1] && misha[2] && misha[3])
#endif
			{
				float* rRow=&localR[(idx%rSize)*_major*RealPerWord*Channels];
				const float* xRows[2*Degree+1];
				float zeroValues[2*Degree+1][Degree][Channels],tempValues1[2*Degree+1][Degree][Channels],tempValues2[2*Degree+1][Degree][Channels];
				for(int y=0;y<2*Degree+1;y++)	xRows[y]=GetXRow(idx-Degree+y,0);
				for(int y=0;y<2*Degree+1;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	zeroValues[y][x][c]=0;
#if NEW_MISHA_CODE == 1
				misha[0]->SetResidual(idx,minor,rRow,xRows,GetBRow(idx,0),zeroValues,zeroValues);
#elif NEW_MISHA_CODE == 2
				for(int y=0;y<2*Degree+1;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues2[y][x][c]=GetXRow(idx-Degree+y,c)[misha[0]->size()+x];
				misha[0]->SetResidual(idx,minor,rRow,xRows,GetBRow(idx,0),zeroValues,tempValues2);
				for(int y=0;y<2*Degree+1;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues1[y][x][c]=GetXRow(idx-Degree+y,c)[misha[0]->size()-Degree+x];
				misha[1]->SetResidual(idx,minor,rRow,xRows,GetBRow(idx,0),tempValues1,zeroValues);
#elif NEW_MISHA_CODE == 4
				for(int y=0;y<2*Degree+1;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues2[y][x][c]=GetXRow(idx-Degree+y,c)[misha[0]->size()+x];
				misha[0]->SetResidual(idx,minor,rRow,xRows,GetBRow(idx,0),zeroValues,tempValues2);
				for(int y=0;y<2*Degree+1;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues1[y][x][c]=GetXRow(idx-Degree+y,c)[misha[0]->size()-Degree+x];
				for(int y=0;y<2*Degree+1;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues2[y][x][c]=GetXRow(idx-Degree+y,c)[misha[0]->size()+misha[1]->size()+x];
				misha[1]->SetResidual(idx,minor,rRow,xRows,GetBRow(idx,0),tempValues1,tempValues2);
				for(int y=0;y<2*Degree+1;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues1[y][x][c]=GetXRow(idx-Degree+y,c)[misha[0]->size()+misha[1]->size()-Degree+x];
				for(int y=0;y<2*Degree+1;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues2[y][x][c]=GetXRow(idx-Degree+y,c)[misha[0]->size()+misha[1]->size()+misha[2]->size()+x];
				misha[2]->SetResidual(idx,minor,rRow,xRows,GetBRow(idx,0),tempValues1,tempValues2);
				for(int y=0;y<2*Degree+1;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues1[y][x][c]=GetXRow(idx-Degree+y,c)[misha[0]->size()+misha[1]->size()+misha[2]->size()-Degree+x];
				misha[3]->SetResidual(idx,minor,rRow,xRows,GetBRow(idx,0),tempValues1,zeroValues);
#endif
			}
			else
#endif // NEW_MISHA_CODE

			SetResidual(idx);

		for(int c=0;c<Channels;c++)
		{
			localBPtr=GetBRow(idx,c);
			Real* localXPtr=GetXRow(idx,c);
			localRPtr=&localR[(idx%rSize)*_major*RealPerWord*Channels+c*_major*RealPerWord];
			for(int i=0;i<major;i++)
			{
				rSquareNorm+=double(localRPtr[i])*localRPtr[i];
				bSquareNorm+=double(localBPtr[i])*localBPtr[i];
				xSquareNorm+=double(localXPtr[i])*localXPtr[i];
			}
		}
	}
}
//////////////////////////////
// MultiGridStreamingSolver //
//////////////////////////////
template<class Real,int Type,int Degree,int Channels,class StorageType>
MultiGridStreamingSolver<Real,Type,Degree,Channels,StorageType>::MultiGridStreamingSolver(void)
{
#if MERGE_INTERIOR_STREAMS
	I=NULL;
#else // !MERGE_INTERIOR_STREAMS
	X=B=NULL;
#endif // MERGE_INTERIOR_STREAMS
	inX=outX=inB=outB=NULL;
	parent=NULL;
	rChild=NULL;
	pChild=NULL;
	localRAccum=_localRAccum=NULL;
#if USE_SSE_CODE
	prolongationStencil2=_prolongationStencil2=NULL;
	prolongationStencil3=_prolongationStencil3=NULL;
#endif // USE_SSE_CODE


	int halfD=((Degree+1)>>1);
	int dual=(Degree&1)?0:1;
	int px2;

	// If we update i at depth d:
	// @ depth d+1 i gets projected to: 	[2*i-halfD,2*i+halfD+dual]
	// @ depth d-1 i gets restricted to:	[(i-halfD-dual)/2,(i+halfD)/2]
	// Can start processing the parent once:
	//		i > dual+halfD
	// Can start processing the child once:
	//		i >= halfD/2


	// Restriction
	// Can start processing the parent once:
	//		i-halfD-dual > 0
	//		i > dual+halfD
	startRestriction=dual+halfD-2*Degree+1;

	// Restriction
	// To solve X[i] at depth d-1:
	//		Need to have B[i+iters-1] at depth d-1:
	//		Need to have R[2*(i+iters-1)-halfD,2*(i+iters-1)+halfD+dual]


	// Prolongation
	// Can start processing the child once:
	//		2*i-halfD >= 0
	//		i >= halfD/2
	startProlongation=(halfD+1)>>1;
	px2=2*startProlongation+halfD+dual;
	//		-iters-Degree+1
	prolongationOffset=px2-(1-Degree);
	prolongationOffset++;


	restrictionBit=(halfD+dual+1)&1;
#if MOLLIFY_DC
	for(int i=0;i<Channels;i++)	dcTerm[i]=0;
#endif // MOLLIFY_DC
#if MERGE_INTERIOR_STREAMS
	I=NULL;
#else // !MERGE_INTERIOR_STREAMS
	B=X=NULL;
#endif // MERGE_INTERIOR_STREAMS
}
template<class Real,int Type,int Degree,int Channels,class StorageType>
void MultiGridStreamingSolver<Real,Type,Degree,Channels,StorageType>::Initialize(int major,int minor,bool symmetric,bool memoryMappedFile)
{
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::FullDotProductStencil dotMajor,d2DotMajor,dotMinor,d2DotMinor;

	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(major,dotMajor,0,0);
	if(symmetric)	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(major,d2DotMajor,1,1,false);
	else			FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(major,d2DotMajor,2,0,true);
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(minor,dotMinor,0,0);
	if(symmetric)	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(minor,d2DotMinor,1,1,false);
	else			FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(minor,d2DotMinor,2,0,true);

	Initialize(dotMajor,d2DotMajor,dotMinor,d2DotMinor,major,minor,symmetric,memoryMappedFile);
}
template<class Real,int Type,int Degree,int Channels,class StorageType>
void MultiGridStreamingSolver<Real,Type,Degree,Channels,StorageType>::Initialize(typename FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::FullDotProductStencil& dotMajor,
																				 typename FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::FullDotProductStencil& d2DotMajor,
																				 typename FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::FullDotProductStencil& dotMinor,
																				 typename FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::FullDotProductStencil& d2DotMinor,
																				 int major,int minor,bool symmetric,bool memoryMappedFile)
{
#if MERGE_INTERIOR_STREAMS
	if(I)	delete I;
	I=NULL;
#else // !MERGE_INTERIOR_STREAMS
	if(B)	delete B;
	if(X)	delete X;
	B=X=NULL;
#endif // MERGE_INTERIOR_STREAMS
	if(memoryMappedFile)
	{
#if MERGE_INTERIOR_STREAMS
		I=new MultiStreamIOClient(2*major*Channels*sizeof(Real),minor,STREAMING_GRID_BUFFER_MULTIPLIER);
#else // !MERGE_INTERIOR_STREAMS
#if SEPARTE_TEMP_IO_SERVER
		B=new MultiStreamIOClient(  major*Channels*sizeof(StorageType),minor,STREAMING_GRID_BUFFER_MULTIPLIER,"g:\Temp","poisson_scratch_");
		X=new MultiStreamIOClient(  major*Channels*sizeof(StorageType),minor,STREAMING_GRID_BUFFER_MULTIPLIER,".","poisson_scratch_");
#else // !SEPARTE_TEMP_IO_SERVER
#if TEST_COMPRESSION
		B=new CompressedFileBackedGrid(major*Channels*sizeof(StorageType),minor);
		X=new CompressedFileBackedGrid(major*Channels*sizeof(StorageType),minor);
#else // !TEST_COMPRESSION
		B=new MultiStreamIOClient(  major*Channels*sizeof(StorageType),minor,STREAMING_GRID_BUFFER_MULTIPLIER);
		X=new MultiStreamIOClient(  major*Channels*sizeof(StorageType),minor,STREAMING_GRID_BUFFER_MULTIPLIER);
#endif // TEST_COMPRESSION
#endif // SEPARTE_TEMP_IO_SERVER
#endif // MERGE_INTERIOR_STREAMS
	}
	else
	{
#if MERGE_INTERIOR_STREAMS
		I=new MemoryBackedGrid(2*major*Channels*sizeof(Real),minor);
#else // !MERGE_INTERIOR_STREAMS
printf( "Allocating new memory backed grid B %d x %d in MultiGridStreamingSolver<Real,Type,Degree,Channels,StorageType>::Initialize" , major , minor );
		B=new MemoryBackedGrid( major*Channels*sizeof(Real),minor);
printf( "Allocating new memory backed grid X %d x %d in MultiGridStreamingSolver<Real,Type,Degree,Channels,StorageType>::Initialize" , major , minor );
		X=new MemoryBackedGrid( major*Channels*sizeof(Real),minor);
printf( "Done\n" );
#endif // MERGE_INTERIOR_STREAMS
	}
	FiniteElements2D<double,Type,Degree>::FullMatrixStencil lStencil;
	for(int i=0;i<=2*Degree;i++)
		for(int j=0;j<=2*Degree;j++)
			for(int k=0;k<=2*Degree;k++)
				for(int l=0;l<=2*Degree;l++)
					lStencil.caseTable[j][i].values[l][k]=
					dotMajor.caseTable[i].values[k]*d2DotMinor.caseTable[j].values[l]+
					d2DotMajor.caseTable[i].values[k]*dotMinor.caseTable[j].values[l];
	Init(lStencil,major,minor,symmetric);
	int major2,minor2;
	IsDownSamplable(major,minor,major2,minor2);
	if(parent && IsDownSamplable(major,minor,major2,minor2))
	{
		FiniteElements1D<double,Type,Degree>::ProlongationStencil(major2,majorProlongationStencil,major);
		FiniteElements1D<double,Type,Degree>::ProlongationStencil(minor2,minorProlongationStencil,minor);
		FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::FullDotProductStencil newDotMajor,newD2DotMajor,newDotMinor,newD2DotMinor;
		CombineStencils<double,Type,Degree>(  dotMajor,majorProlongationStencil,major,  newDotMajor);
		CombineStencils<double,Type,Degree>(d2DotMajor,majorProlongationStencil,major,newD2DotMajor);
		CombineStencils<double,Type,Degree>(  dotMinor,minorProlongationStencil,minor,  newDotMinor);
		CombineStencils<double,Type,Degree>(d2DotMinor,minorProlongationStencil,minor,newD2DotMinor);
		parent->Initialize(newDotMajor,newD2DotMajor,newDotMinor,newD2DotMinor,major2,minor2,symmetric,memoryMappedFile);
	}
}
template<class Real,int Type,int Degree,int Channels,class StorageType>
MultiGridStreamingSolver<Real,Type,Degree,Channels,StorageType>::~MultiGridStreamingSolver(void)
{
	if(_localRAccum)	free(_localRAccum);
	_localRAccum=localRAccum=NULL;
#if USE_SSE_CODE
	if(_prolongationStencil2)	free(_prolongationStencil2);
	prolongationStencil2=_prolongationStencil2=NULL;
	if(_prolongationStencil3)	free(_prolongationStencil3);
	prolongationStencil3=_prolongationStencil3=NULL;
#endif // USE_SSE_CODE
#if MERGE_INTERIOR_STREAMS
	if(I)	delete I;
	I=NULL;
#else // !MERGE_INTERIOR_STRAMS
	if(X)	delete X;
	if(B)	delete B;
	X=B=NULL;
#endif // MERGE_INTERIOR_STREAMS
}
template<class Real,int Type,int Degree,int Channels,class StorageType>
void MultiGridStreamingSolver<Real,Type,Degree,Channels,StorageType>::InitProlongation(bool symmetric)
{
#if MOLLIFY_DC
	for(int i=0;i<Channels;i++)	dcTerm[i]/=major*minor;
printf("%d x %d\t",major,minor);
for(int i=0;i<Channels;i++)	printf("%f ",dcTerm[i]);
printf("\n");
#endif // MOLLIFY_DC
	if(parent)	parent->InitProlongation(symmetric);
	if(inX)		inX->reset(true,1);
	if(inB)		inB->reset(true,1);
	if(outX)	outX->reset(false,1);
	if(outB)	outB->reset(false,1);
}
template<class Real,int Type,int Degree,int Channels,class StorageType>
void MultiGridStreamingSolver<Real,Type,Degree,Channels,StorageType>::InitRestriction(bool symmetric)
{
#if MOLLIFY_DC
	for(int i=0;i<Channels;i++)	dcTerm[i]=0;
#endif // MOLLIFY_DC
	int major2,minor2;
	if(parent && IsDownSamplable(major,minor,major2,minor2) )
	{
		_localRAccum=(WordClass*)malloc(sizeof(WordClass)*(_major+2*Degree)+ALIGNMENT-1);
		localRAccum=(WordClass*)(((size_t)(_localRAccum+Degree)+ALIGNMENT-1) & ~(ALIGNMENT-1));
#if USE_SSE_CODE
		_prolongationStencil2=(ProlongationStencilSSE*)malloc(sizeof(ProlongationStencilSSE)*3*(2*Degree+1)+ALIGNMENT-1);
		memset(_prolongationStencil2,0,sizeof(ProlongationStencilSSE)*3*(2*Degree+1)+ALIGNMENT-1);
		prolongationStencil2=(ProlongationStencilSSE*)(((size_t)(_prolongationStencil2)+ALIGNMENT-1) & ~(ALIGNMENT-1));
		__declspec (align(16)) float scratch[4];
		for(int i=0;i<=2*Degree;i++)	// Iterate over the minor index in the mask
			for(int j=0;j<3;j++)
				for(int k=0;k<Degree+2;k++)
				{
					int jj;
					if(j==0)	jj=0;
					else		jj=Degree;
					scratch[0]=majorProlongationStencil.caseTable[jj].values[1]*minorProlongationStencil.caseTable[i].values[k];
					scratch[1]=majorProlongationStencil.caseTable[jj].values[2]*minorProlongationStencil.caseTable[i].values[k];
					scratch[2]=majorProlongationStencil.caseTable[jj].values[3]*minorProlongationStencil.caseTable[i].values[k];
					if(j!=2)	scratch[3]=majorProlongationStencil.caseTable[Degree].values[0]*minorProlongationStencil.caseTable[i].values[k];
					else		scratch[3]=0;
					prolongationStencil2[3*i+j].matrixValues[0][k]=_mm_load_ps(scratch);

					if(j==2)	jj=2*Degree;
					else		jj=Degree;
					if(j!=0)	scratch[0]=majorProlongationStencil.caseTable[Degree].values[3]*minorProlongationStencil.caseTable[i].values[k];
					else		scratch[0]=0;
					scratch[1]=majorProlongationStencil.caseTable[jj].values[0]*minorProlongationStencil.caseTable[i].values[k];
					scratch[2]=majorProlongationStencil.caseTable[jj].values[1]*minorProlongationStencil.caseTable[i].values[k];
					scratch[3]=majorProlongationStencil.caseTable[jj].values[2]*minorProlongationStencil.caseTable[i].values[k];
					prolongationStencil2[3*i+j].matrixValues[1][k]=_mm_load_ps(scratch);
				}
		_prolongationStencil3=(ProlongationStencilSSE2*)malloc(sizeof(ProlongationStencilSSE2)*3+ALIGNMENT-1);
		memset(_prolongationStencil3,0,sizeof(ProlongationStencilSSE2)*3+ALIGNMENT-1);
		prolongationStencil3=(ProlongationStencilSSE2*)(((size_t)(_prolongationStencil3)+ALIGNMENT-1) & ~(ALIGNMENT-1));
		for(int j=0;j<3;j++)
		{
			int jj;
			if(j==0)	jj=0;
			else		jj=Degree;
			scratch[0]=majorProlongationStencil.caseTable[jj].values[1];
			scratch[1]=majorProlongationStencil.caseTable[jj].values[2];
			scratch[2]=majorProlongationStencil.caseTable[jj].values[3];
			if(j!=2)	scratch[3]=majorProlongationStencil.caseTable[Degree].values[0];
			else		scratch[3]=0;
			prolongationStencil3[j].matrixValues[0]=_mm_load_ps(scratch);

			if(j==2)	jj=2*Degree;
			else		jj=Degree;
			if(j!=0)	scratch[0]=majorProlongationStencil.caseTable[Degree].values[3];
			else		scratch[0]=0;
			scratch[1]=majorProlongationStencil.caseTable[jj].values[0];
			scratch[2]=majorProlongationStencil.caseTable[jj].values[1];
			scratch[3]=majorProlongationStencil.caseTable[jj].values[2];
			prolongationStencil3[j].matrixValues[1]=_mm_load_ps(scratch);
		}
#endif // USE_SSE_CODE
	}
	if(parent)	parent->InitRestriction(symmetric);
	if(inX)		inX->reset(true,1);
	if(inB)		inB->reset(true,1);
	if(outX)	outX->reset(false,1);
	if(outB)	outB->reset(false,1);
#if NEW_MISHA_CODE
	for(int i=0;i<NEW_MISHA_CODE;i++)	if(misha[i])	misha[i]->SetStencils(majorProlongationStencil,minorProlongationStencil,prolongationStencil3);
#endif // NEW_MISHA_CODE
}
template<class Real,int Type,int Degree,int Channels,class StorageType>
void MultiGridStreamingSolver<Real,Type,Degree,Channels,StorageType>::SetProlongationIterations(int iters)
{
printf("hi 1\n");
	StreamingSolver<Real,Type,Degree,Channels>::SetIterations(-Degree-iters*Degree+1,iters,0,0,0,prolongationOffset+iters*Degree);
	// Set the iterations for the child so that it trails accordingly
	if(pChild)	pChild->SetProlongationIterations(iters);
#if USE_SERVER
#if MERGE_INTERIOR_STREAMS
	I->SetServer(&server);
#else // !MERGE_INTERIOR_STREAMS
#if SEPARTE_TEMP_IO_SERVER
	B->SetServer(&tempServer);
	X->SetServer(&server);
#else // !SEPARTE_TEMP_IO_SERVER
	B->SetServer(&server);
	X->SetServer(&server);
#endif // SEPARTE_TEMP_IO_SERVER
#endif // MERGE_INTERIOR_STREAMS
#endif // USE_SERVER
#if MERGE_INTERIOR_STREAMS
	I->reset(true,false,iSize);
#else // !MERGE_INTERIOR_STREAMS
//	X->reset(true,false,xSize);
//	B->reset(true,false,bSize);
	X->reset(true,xSize);
	B->reset(true,bSize);
#endif // MERGE_INTERIOR_STREAMS
}
template<class Real,int Type,int Degree,int Channels,class StorageType>
void MultiGridStreamingSolver<Real,Type,Degree,Channels,StorageType>::SetRestrictionIterations(int iters)
{
printf("hi 2\n");
	StreamingSolver<Real,Type,Degree,Channels>::SetIterations(-Degree-iters*Degree+1,iters,0,0,0,0,Degree+2);
	// Set the iterations for the parent so that it trails accordingly
	if(parent)	parent->SetRestrictionIterations(iters);
#if USE_SERVER
#if MERGE_INTERIOR_STREAMS
	I->SetServer(&server);
#else // !MERGE_INTERIOR_STREAMS
#if SEPARTE_TEMP_IO_SERVER
	B->SetServer(&tempServer);
	X->SetServer(&server);
#else // !SEPARTE_TEMP_IO_SERVER
	B->SetServer(&server);
	X->SetServer(&server);
#endif // SEPARTE_TEMP_IO_SERVER
#endif // MERGE_INTERIOR_STREAMS
#endif // USE_SERVER
#if MERGE_INTERIOR_STREAMS
	I->reset(false,true,iSize);
#else // !MERGE_INTERIOR_STREAMS
//	X->reset(false,true,xSize);
//	B->reset(false,true,bSize);
	X->reset(false,xSize);
	B->reset(false,bSize);
#endif // MERGE_INTERIOR_STREAMS
}
template<class Real,int Type,int Degree,int Channels,class StorageType>
void MultiGridStreamingSolver<Real,Type,Degree,Channels,StorageType>::UnSetProlongationIterations(void)
{
	StreamingSolver<Real,Type,Degree,Channels>::UnSetIterations();
	X->unset();
	B->unset();
	if(pChild)	pChild->UnSetProlongationIterations();
}
template<class Real,int Type,int Degree,int Channels,class StorageType>
void MultiGridStreamingSolver<Real,Type,Degree,Channels,StorageType>::UnSetRestrictionIterations(void)
{
	StreamingSolver<Real,Type,Degree,Channels>::UnSetIterations();
	X->unset();
	B->unset();
	if(parent)	parent->UnSetRestrictionIterations();
}
template<class Real,int Type,int Degree,int Channels,class StorageType>
bool MultiGridStreamingSolver<Real,Type,Degree,Channels,StorageType>::IterateProlongation(void)
{
	if(index-1>=minor)	return false;
	int idx=index+_iters+Degree-1;
	if(idx>=0)
	{
#if MERGE_INTERIOR_STREAMS
		if(I)	if(UpdateIInput(I))	I->advance();
#else // !MERGE_INTERIOR_STREAMS
		if(UpdateBInput<StorageType>(B))	B->advance();
#if !STORE_RESIDUAL
		if(UpdateXInput<StorageType>(X))	X->advance();
#endif // !STORE_RESIDUAL
#endif // MERGE_INTERIOR_STREAMS
#if MOLLIFY_DC
		if(idx>=0 && idx<minor)
			for(int c=0;c<Channels;c++)
			{
				Real* lX=GetXRow(idx,c);
				for(int i=0;i<major;i++)	lX[i]-=dcTerm[c];
			}
#endif // MOLLIFY_DC
		if(!parent && inX && idx>=0 && idx<minor)
		{
			Real* inPtr=(Real*)(*inX)[idx];
			for(int c=0;c<Channels;c++)
			{
				int idx1=c*major;
				Real* x=GetXRow(idx,c);
				for(int i=0;i<major;i++)			x[i]=inPtr[idx1+i]*laplacianScaleR;
			}
			inX->advance();
		}

		// Run an interation of the Gauss-Seidel solver
		StreamingSolver<Real,Type,Degree,Channels>::Solve();
#if STORE_RESIDUAL
		if(UpdateXInput<StorageType>(X))	X->advance();
#endif // STORE_RESIDUAL
		if(!pChild && index>=0 && index<minor)
		{
			if(outX)
			{
				Real* outPtr=(Real*)(*outX)[index];
				for(int c=0;c<Channels;c++)
				{
					int idx1=c*major;
					Real* x=GetXRow(index,c);
					for(int i=0;i<major;i++)	outPtr[idx1+i]=x[i]*laplacianScale;
				}
				outX->advance();
			}
			else
			{
//				fprintf(stderr,"Badness: no output stream\n");
//				exit(0);
			}
		}
	}
	if(pChild)
	{
		if(index>=0 && index<minor)		pChild->ProlongationUpdate(index);
		if(index>=startProlongation)	if(pChild->IterateProlongation())	pChild->IterateProlongation();
	}
	return Increment();
}
template<class Real,int Type,int Degree,int Channels,class StorageType>
void MultiGridStreamingSolver<Real,Type,Degree,Channels,StorageType>::SetRestriction(Real* lB,int idx,int major2,int minor2)
{
	if(idx>=Degree && idx<minor2-Degree)
	{
#if NEW_MISHA_CODE
#if NEW_MISHA_CODE == 1
		if(misha[0])
#elif NEW_MISHA_CODE == 2
		if(misha[0] && misha[1])
#elif NEW_MISHA_CODE == 4
		if(misha[0] && misha[1] && misha[2] && misha[3])
#endif
		{
			int halfD=(Degree+1)>>1;
			int startY=2*idx-halfD;
			const float* rRows[Degree+2];
			float zeroValues[Degree+2][Degree][Channels],tempValues1[Degree+2][Degree][Channels],tempValues2[Degree+2][Degree][Channels];
			for(int y=0;y<Degree+2;y++)	rRows[y]=&localR[((startY+y)%rSize)*_major*RealPerWord*Channels];
			for(int y=0;y<Degree+2;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	zeroValues[y][x][c]=0;
#if NEW_MISHA_CODE == 1
			misha[0]->SetInteriorRestriction(lB,rRows,zeroValues,zeroValues);
#elif NEW_MISHA_CODE == 2
			for(int y=0;y<Degree+2;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues2[y][x][c]=rRows[y][_major*RealPerWord*c+misha[0]->size()+x];
			misha[0]->SetInteriorRestriction(lB,rRows,zeroValues,tempValues2);
			for(int y=0;y<Degree+2;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues1[y][x][c]=rRows[y][_major*RealPerWord*c+misha[0]->size()-Degree+x];
			misha[1]->SetInteriorRestriction(lB,rRows,tempValues1,zeroValues);
#else NEW_MISHA_CODE == 4
			for(int y=0;y<Degree+2;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues2[y][x][c]=rRows[y][_major*RealPerWord*c+misha[0]->size()+x];
			misha[0]->SetInteriorRestriction(lB,rRows,zeroValues,tempValues2);
			for(int y=0;y<Degree+2;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues1[y][x][c]=rRows[y][_major*RealPerWord*c+misha[0]->size()-Degree+x];
			for(int y=0;y<Degree+2;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues2[y][x][c]=rRows[y][_major*RealPerWord*c+misha[0]->size()+misha[1]->size()+x];
			misha[1]->SetInteriorRestriction(lB,rRows,tempValues1,tempValues2);
			for(int y=0;y<Degree+2;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues1[y][x][c]=rRows[y][_major*RealPerWord*c+misha[0]->size()+misha[1]->size()-Degree+x];
			for(int y=0;y<Degree+2;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues2[y][x][c]=rRows[y][_major*RealPerWord*c+misha[0]->size()+misha[1]->size()+misha[2]->size()+x];
			misha[2]->SetInteriorRestriction(lB,rRows,tempValues1,tempValues2);
			for(int y=0;y<Degree+2;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues1[y][x][c]=rRows[y][_major*RealPerWord*c+misha[0]->size()+misha[1]->size()+misha[2]->size()-Degree+x];
			misha[3]->SetInteriorRestriction(lB,rRows,tempValues1,zeroValues);
#endif
		}
		else
#endif // NEW_MISHA_CODE
			SetInteriorRestriction(lB,idx,major2,minor2);
		return;
	}
#if NEW_MISHA_CODE
#if NEW_MISHA_CODE == 1
		if(misha[0])
#elif NEW_MISHA_CODE == 2
		if(misha[0] && misha[1])
#elif NEW_MISHA_CODE == 4
		if(misha[0] && misha[1] && misha[2] && misha[3])
#endif
		{
			int halfD=(Degree+1)>>1;
			int startY=2*idx-halfD;
			const float* rRows[Degree+2];
			float zeroValues[Degree+2][Degree][Channels],tempValues1[Degree+2][Degree][Channels],tempValues2[Degree+2][Degree][Channels];
			for(int y=0;y<Degree+2;y++)
				if(startY+y>=0 && startY+y<minor)	rRows[y]=&localR[((startY+y)%rSize)*_major*RealPerWord*Channels];
				else								rRows[y]=NULL;
			for(int y=0;y<Degree+2;y++)	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	zeroValues[y][x][c]=tempValues1[y][x][c]=tempValues2[y][x][c]=0;
#if NEW_MISHA_CODE == 1
			misha[0]->SetRestriction(idx,minor2,lB,rRows,zeroValues,zeroValues);
#elif NEW_MISHA_CODE == 2
			for(int y=0;y<Degree+2;y++)	if(rRows[y])	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues2[y][x][c]=rRows[y][_major*RealPerWord*c+misha[0]->size()+x];
			misha[0]->SetRestriction(idx,minor2,lB,rRows,zeroValues,tempValues2);
			for(int y=0;y<Degree+2;y++)	if(rRows[y])	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues1[y][x][c]=rRows[y][_major*RealPerWord*c+misha[0]->size()-Degree+x];
			misha[1]->SetRestriction(idx,minor2,lB,rRows,tempValues1,zeroValues);
#elif NEW_MISHA_CODE == 4
			for(int y=0;y<Degree+2;y++)	if(rRows[y])	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues2[y][x][c]=rRows[y][_major*RealPerWord*c+misha[0]->size()+x];
			misha[0]->SetRestriction(idx,minor2,lB,rRows,zeroValues,tempValues2);
			for(int y=0;y<Degree+2;y++)	if(rRows[y])	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues1[y][x][c]=rRows[y][_major*RealPerWord*c+misha[0]->size()-Degree+x];
			for(int y=0;y<Degree+2;y++)	if(rRows[y])	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues2[y][x][c]=rRows[y][_major*RealPerWord*c+misha[0]->size()+misha[1]->size()+x];
			misha[1]->SetRestriction(idx,minor2,lB,rRows,tempValues1,tempValues2);
			for(int y=0;y<Degree+2;y++)	if(rRows[y])	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues1[y][x][c]=rRows[y][_major*RealPerWord*c+misha[0]->size()+misha[1]->size()-Degree+x];
			for(int y=0;y<Degree+2;y++)	if(rRows[y])	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues2[y][x][c]=rRows[y][_major*RealPerWord*c+misha[0]->size()+misha[1]->size()+misha[2]->size()+x];
			misha[2]->SetRestriction(idx,minor2,lB,rRows,tempValues1,tempValues2);
			for(int y=0;y<Degree+2;y++)	if(rRows[y])	for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues1[y][x][c]=rRows[y][_major*RealPerWord*c+misha[0]->size()+misha[1]->size()+misha[2]->size()-Degree+x];
			misha[3]->SetRestriction(idx,minor2,lB,rRows,tempValues1,zeroValues);
#endif
		}
		else
#endif // NEW_MISHA_CODE
	{
		int halfD=(Degree+1)>>1;
		int startY=2*idx-halfD;
		for(int c=0;c<Channels;c++)
		{
			Real* myLB=&lB[major2*c];
			for(int yy=0;yy<Degree+2;yy++)
			{
				if(startY+yy>=0 && startY+yy<minor)	localRPtrs[yy]=(WordClass*)(localR + ((startY+yy)%rSize)*_major*RealPerWord*Channels+c*_major*RealPerWord);
				else								localRPtrs[yy]=NULL;
			}

			int jj;
			if		(idx<Degree)			jj=idx;
			else if	(idx>minor2-1-Degree)	jj=2*Degree+(idx-(minor2-1));
			else							jj=Degree;
#if USE_SSE_CODE
			jj*=3;
			float dotSum;
			WordClass dSum;

			__declspec (align(ALIGNMENT)) float scratch[4];
			const WordClass* rPtrs[]={localRPtrs[0],localRPtrs[1],localRPtrs[2],localRPtrs[3]};
			{
				dotSum=0;
				myLB[0]=RestrictionUpdate0(prolongationStencil2[jj].matrixValues[0],rPtrs,dotSum,0);
				WordClass mValues[]=
				{
					prolongationStencil2[jj+1].matrixValues[0][0],
					prolongationStencil2[jj+1].matrixValues[0][1],
					prolongationStencil2[jj+1].matrixValues[0][2],
					prolongationStencil2[jj+1].matrixValues[0][3]
				};
				int bound=major2-2;
				for(int i=2;i<bound;i+=2)	myLB[i]=RestrictionUpdate0(mValues,rPtrs,dotSum,i);
				if(major2>=4)	myLB[major2-2]=RestrictionUpdate0(prolongationStencil2[jj+2].matrixValues[0],rPtrs,dotSum,major2-2);
			}
			{
				SetRestrictionDotSum(prolongationStencil2[jj].matrixValues[1],rPtrs,0,dSum);
				_mm_store_ps(scratch,dSum);
				dotSum=scratch[1]+scratch[2]+scratch[3];
				WordClass mValues[]=
				{
					prolongationStencil2[jj+1].matrixValues[1][0],
					prolongationStencil2[jj+1].matrixValues[1][1],
					prolongationStencil2[jj+1].matrixValues[1][2],
					prolongationStencil2[jj+1].matrixValues[1][3]
				};
				if(major2<6)	myLB[1]=RestrictionUpdate1(prolongationStencil2[jj+2].matrixValues[1],rPtrs,dotSum,1);
				else			myLB[1]=RestrictionUpdate1(mValues,rPtrs,dotSum,1);
				int bound=major2-4;
				for(int i=3;i<bound;i+=2)	myLB[i]=RestrictionUpdate1(mValues,rPtrs,dotSum,i);
				if(major2>=6)				myLB[major2-3]=RestrictionUpdate1(prolongationStencil2[jj+2].matrixValues[1],rPtrs,dotSum,major2-3);
				myLB[major2-1]=dotSum;
			}
#else // !USE_SSE_CODE
			for(int i=0;i<major2;i++)
			{
				int ii;
				if		(i<Degree)			ii=i;
				else if	(i>major2-1-Degree)	ii=2*Degree+(i-(major2-1));
				else						ii=Degree;
				Real temp=0;

				int startX=2*i-halfD;
				for(int yy=0;yy<Degree+2;yy++)
				{
					if(startY+yy<0 || startY+yy>=minor)	continue;
					Real* residualPtr=((Real*)localRPtrs[yy])+startX;
					Real t=0;
					for(int xx=0;xx<Degree+2;xx++)	t+=residualPtr[xx]*majorProlongationStencil.caseTable[ii].values[xx];
					temp+=t*minorProlongationStencil.caseTable[jj].values[yy];
				}
				myLB[i]=temp;
			}
#endif // USE_SSE_CODE
		}
	}
}
template<class Real,int Type,int Degree,int Channels,class StorageType>
void MultiGridStreamingSolver<Real,Type,Degree,Channels,StorageType>::SetInteriorRestriction(Real* lB,int idx,int major2,int minor2)
{
	int halfD=(Degree+1)>>1;
	int startY=2*idx-halfD;
	for(int c=0;c<Channels;c++)
	{
		Real* myLB=&lB[c*major2];
		for(int yy=0;yy<Degree+2;yy++)	localRPtrs[yy]=(WordClass*)(localR + ((startY+yy)%rSize)*_major*RealPerWord*Channels+c*_major*RealPerWord);
#if USE_SSE_CODE
		__declspec (align(16)) float scratch[4];
		__m128 res[Degree+2];
		for(int d=0;d<Degree+2;d++)
		{
			for(int j=0;j<4;j++)	scratch[j]= minorProlongationStencil.caseTable[Degree].values[d];
			res[d]=_mm_load_ps(scratch);
		}
		{
			int d=0;
			for(int i=0;i<_major;i++)	localRAccum[i]=_mm_mul_ps(res[d],_mm_add_ps(localRPtrs[d][i],localRPtrs[Degree+1-d][i]));
		}
		for(int d=1;d<=(Degree>>1);d++)
			for(int i=0;i<_major;i++)	localRAccum[i]=_mm_add_ps(localRAccum[i],_mm_mul_ps(res[d],_mm_add_ps(localRPtrs[d][i],localRPtrs[Degree+1-d][i])));
		if(Degree&1)
		{
			int d=(Degree+1)>>1;
			for(int i=0;i<_major;i++)	localRAccum[i]=_mm_add_ps(localRAccum[i],_mm_mul_ps(res[d],_mm_add_ps(localRPtrs[d][i],localRPtrs[Degree+1-d][i])));
		}
#else // !USE_SSE_CODE
		{
			int d=0;
			for(int i=0;i<major;i++)	localRAccum[i] =(localRPtrs[d][i]+localRPtrs[Degree+1-d][i])*minorProlongationStencil.caseTable[Degree].values[d];
		}
		for(int d=1;d<=(Degree>>1);d++)
			for(int i=0;i<major;i++)	localRAccum[i]+=(localRPtrs[d][i]+localRPtrs[Degree+1-d][i])*minorProlongationStencil.caseTable[Degree].values[d];
		if(Degree&1)
		{
			int d=(Degree+1)>>1;
			for(int i=0;i<major;i++)	localRAccum[i]+=localRPtrs[d][i]*minorProlongationStencil.caseTable[Degree].values[d];
		}
#endif // USE_SSE_CODE
		int jj;
		if		(idx<Degree)			jj=idx;
		else if	(idx>minor-1-Degree)	jj=2*Degree+(idx-(minor-1));
		else							jj=Degree;
#if USE_SSE_CODE
		float dotSum;
		WordClass dSum;

		const WordClass* rPtrs=localRAccum;
		{
			dotSum=0;
			myLB[0]=RestrictionUpdate0(prolongationStencil3[0].matrixValues[0],rPtrs,dotSum,0);
			int bound=major2-2;
			for(int i=2;i<bound;i+=2)	myLB[i]=RestrictionUpdate0(prolongationStencil3[1].matrixValues[0],rPtrs,dotSum,i);
			if(major2>=4)				myLB[major2-2]=RestrictionUpdate0(prolongationStencil3[2].matrixValues[0],rPtrs,dotSum,major2-2);
		}
		{
			SetRestrictionDotSum(prolongationStencil3[0].matrixValues[1],rPtrs,0,dSum);
			_mm_store_ps(scratch,dSum);
			dotSum=scratch[1]+scratch[2]+scratch[3];
			if(major2<6)				myLB[1]=RestrictionUpdate1(prolongationStencil3[2].matrixValues[1],rPtrs,dotSum,1);
			else						myLB[1]=RestrictionUpdate1(prolongationStencil3[1].matrixValues[1],rPtrs,dotSum,1);
			int bound=major2-4;
			for(int i=3;i<bound;i+=2)	myLB[i]=RestrictionUpdate1(prolongationStencil3[1].matrixValues[1],rPtrs,dotSum,i);
			if(major2>=6)				myLB[major2-3]=RestrictionUpdate1(prolongationStencil3[2].matrixValues[1],rPtrs,dotSum,major2-3);
			myLB[major2-1]=dotSum;
		}
#else // !USE_SSE_CODE
		for(int i=0;i<major2;i++)
		{
			int ii;
			if		(i<Degree)			ii=i;
			else if	(i>major2-1-Degree)	ii=2*Degree+(i-(major2-1));
			else						ii=Degree;
			Real temp=0;

			int startX=2*i-halfD;
			for(int yy=0;yy<Degree+2;yy++)
			{
				if(startY+yy<0 || startY+yy>=minor)	continue;
				Real* residualPtr=((Real*)localRPtrs[yy])+startX;
				Real t=0;
				for(int xx=0;xx<Degree+2;xx++)	t+=residualPtr[xx]*majorProlongationStencil.caseTable[ii].values[xx];
				temp+=t*minorProlongationStencil.caseTable[jj].values[yy];
			}
			myLB[i]=temp;
		}
#endif // USE_SSE_CODE
	}
}

template<class Real,int Type,int Degree,int Channels,class StorageType>
bool MultiGridStreamingSolver<Real,Type,Degree,Channels,StorageType>::IterateRestriction(void)
{
	if(index-1>=minor)	return false;
	int idx=index+_iters-1;
	if(idx+Degree>=0 && idx+Degree<minor && inX)
	{
		Real* inPtr=(Real*)(*inX)[idx+Degree];
		for(int c=0;c<Channels;c++)
		{
			int idx1=c*major;
			Real* x=GetXRow(idx+Degree,c);
			for(int i=0;i<major;i++)	x[i]=inPtr[idx1+i]*laplacianScaleR;
		}
		inX->advance();
	}
	if(idx>=0)
	{
		if(idx<minor)
		{
			if(!rChild)
			{
				if(!inB)
				{
					fprintf(stderr,"Badness: no input stream\n");
					exit(0);
				}
				StorageType* inPtr=(StorageType*)(*inB)[idx];
				for(int c=0;c<Channels;c++)
				{
					Real* bRow=GetBRow(idx,c);
					for(int i=0;i<major;i++)	bRow[i]=Real(inPtr[c*major+i]);
				}
				inB->advance();
			}
			else	rChild->SetRestriction(GetBRow(idx,0),idx,major,minor);
		}
		// Run an interation of the Gauss-Seidel solver
		StreamingSolver<Real,Type,Degree,Channels>::Solve();
	}
#if MOLLIFY_DC
	{
		int idx=index-1;
		if(idx>=0 && idx<minor)
			for(int c=0;c<Channels;c++)
			{
				Real* lX=GetXRow(idx,c);
				for(int i=0;i<major;i++)	dcTerm[c]+=lX[i];
			}
	}
#endif // MOLLIFY_DC
#if MERGE_INTERIOR_STREAMS
	if(I)	if(UpdateIOutput(I))	I->advance();
#else // !MERGE_INTERIOR_STREAMS
	// Write out the current solution
	if(UpdateXOutput<StorageType>(X))	X->advance();
	if(UpdateBOutput<StorageType>(B))	B->advance();
#endif // MERGE_INTERIOR_STREAMS
	// Write out the residual
	if(parent)	if(index>=startRestriction && (index&1)==restrictionBit)	parent->IterateRestriction();
	if(!parent && outB && index>=0 && index<minor)
	{
		Real* outPtr=(Real*)(*outB)[index];
		for(int c=0;c<Channels;c++)
		{
			Real* inP=GetBRow(index,c);
			Real* outP=outPtr+c*major;
			for(int i=0;i<major;i++)	outP[i]=inP[i];
		}
		outB->advance();
	}
	if(!parent && outX && index>=0 && index<minor)
	{
		Real* outPtr=(Real*)(*outX)[index];
		for(int c=0;c<Channels;c++)
		{
			Real* inP=GetXRow(index,c);
			Real* outP=outPtr+c*major;
			for(int i=0;i<major;i++)	outP[i]=inP[i]*laplacianScale;
		}
		outX->advance();
	}
	return Increment();
}

template<class Real,int Type,int Degree,int Channels,class StorageType>
void MultiGridStreamingSolver<Real,Type,Degree,Channels,StorageType>::ProlongationUpdate(int idx)
{
#if NEW_MISHA_CODE
#if NEW_MISHA_CODE == 1
	if(misha[0])
#elif NEW_MISHA_CODE == 2
	if(misha[0] && misha[1])
#elif NEW_MISHA_CODE == 4
	if(misha[0] && misha[1] && misha[2] && misha[3])
#endif
	{
		int halfD=(Degree+1)>>1;
		int startY=2*idx-halfD;
		const float* restrictedXRow=parent->GetXRow(idx,0);
		float* xRows[Degree+2];
		float zeroValues[Degree][Channels],tempValues1[Degree][Channels],tempValues2[Degree][Channels];
		for(int y=0;y<Degree+2;y++)	xRows[y]=GetXRow(y+startY,0);
		for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	zeroValues[x][c]=0;
#if NEW_MISHA_CODE == 1
		misha[0]->IncrementProlongation(idx,parent->minor,xRows,restrictedXRow,tempValues1,tempValues2);
#elif NEW_MISHA_CODE == 2
		for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues2[x][c]=restrictedXRow[(_major*RealPerWord*c+misha[0]->size())/2+x];
		misha[0]->IncrementProlongation(idx,parent->minor,xRows,restrictedXRow,tempValues1,tempValues2);
		for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues1[x][c]=restrictedXRow[(_major*RealPerWord*c+misha[0]->size())/2-Degree+x];
		misha[1]->IncrementProlongation(idx,parent->minor,xRows,restrictedXRow,tempValues1,tempValues2);
#elif NEW_MISHA_CODE == 4
		for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues2[x][c]=restrictedXRow[(_major*RealPerWord*c+misha[0]->size())/2+x];
		misha[0]->IncrementProlongation(idx,parent->minor,xRows,restrictedXRow,tempValues1,tempValues2);
		for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues1[x][c]=restrictedXRow[(_major*RealPerWord*c+misha[0]->size())/2-Degree+x];
		for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues2[x][c]=restrictedXRow[(_major*RealPerWord*c+misha[0]->size()+misha[1]->size())/2+x];
		misha[1]->IncrementProlongation(idx,parent->minor,xRows,restrictedXRow,tempValues1,tempValues2);
		for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues1[x][c]=restrictedXRow[(_major*RealPerWord*c+misha[0]->size()+misha[1]->size())/2-Degree+x];
		for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues2[x][c]=restrictedXRow[(_major*RealPerWord*c+misha[0]->size()+misha[1]->size()+misha[2]->size())/2+x];
		misha[2]->IncrementProlongation(idx,parent->minor,xRows,restrictedXRow,tempValues1,tempValues2);
		for(int x=0;x<Degree;x++)	for(int c=0;c<Channels;c++)	tempValues1[x][c]=restrictedXRow[(_major*RealPerWord*c+misha[0]->size()+misha[1]->size()+misha[2]->size())/2-Degree+x];
		misha[3]->IncrementProlongation(idx,parent->minor,xRows,restrictedXRow,tempValues1,tempValues2);
#endif
	}
	else
#endif // NEW_MISHA_CODE
	{
		int major2,minor2;
		major2=parent->major;
		minor2=parent->minor;
		if(idx>=0 && idx<minor2)
		{
			int halfD=(Degree+1)>>1;
			int startY=2*idx-halfD;
			int jj;
			if		(idx<Degree)				jj=idx;
			else if	(idx>minor2-1-Degree)		jj=2*Degree+(idx-(minor2-1));
			else								jj=Degree;
			for(int c=0;c<Channels;c++)
			{
				Real* parentXPtr=parent->GetXRow(idx,c);
				for(int i=0;i<major2;i++)
				{
					int ii;
					if		(i<Degree)					ii=i;
					else if	(i>major2-1-Degree)			ii=2*Degree+(i-(major2-1));
					else								ii=Degree;

					int startX=2*i-halfD;
					for(int yy=0;yy<Degree+2;yy++)
					{
						if(startY+yy<0 || startY+yy>=minor)	continue;
						Real* localXPtr=GetXRow(startY+yy,c)+startX;
						Real value=parentXPtr[i]*minorProlongationStencil.caseTable[jj].values[yy];
						for(int xx=0;xx<Degree+2;xx++)	localXPtr[xx]+=value*majorProlongationStencil.caseTable[ii].values[xx];
					}
				}
			}
		}
	}
}

template<class Real,int Type,int Degree,int Channels,class StorageType>
void MultiGridStreamingSolver<Real,Type,Degree,Channels,StorageType>::SolveProlongation(void)
{
	// Run to completion...
	while( IterateProlongation() ){;}
	// ...and finish up the trailing child
	if(pChild)	pChild->SolveProlongation();
}
template<class Real,int Type,int Degree,int Channels,class StorageType>
void MultiGridStreamingSolver<Real,Type,Degree,Channels,StorageType>::SolveRestriction(void)
{
	// Run to completion...
	while( IterateRestriction() ){;}
	// ...and finish up the trailing parent
	if(parent)	parent->SolveRestriction();
}
template<class Real,int Type,int Degree,int Channels,class StorageType>
bool MultiGridStreamingSolver<Real,Type,Degree,Channels,StorageType>::IsDownSamplable(const int& hMajor,const int&hMinor,int &lMajor,int& lMinor)
{
	return FiniteElements2D<Real,Type,Degree>::IsDownSamplable(hMajor,hMinor,lMajor,lMinor);
}
/////////////////////////
// StreamingDivergence //
/////////////////////////
template<class Real,int Type,int Degree,int Channels,class PartialType>
StreamingDivergence<Real,Type,Degree,Channels,PartialType>::StreamingDivergence(void)
{
	outB=NULL;
	localB=NULL;
	rParent=NULL;
	_localDMajor=_localDMinor=NULL;
	localDMajor=localDMinor=NULL;
	_localDMajorAccum=NULL;
	_localDMinorAccum=NULL;
	localDMajorAccum=NULL;
	localDMinorAccum=NULL;
	dMajor=dMinor=NULL;
	dSize=2*Degree+1;
#if REFLECT_GRADIENT
	dSize+=2*ReflectBufferSize;
#endif // REFLECT_GRADIENT
	// BADNESS!!! Why does this change things? And why is it related to whether or not the verbose flag is set?
	// Could be related to the debug problem. Might be fixed with a memset of localdmajor and localdminor
//	dSize+=12;
	dXSquareNorm=0;
	dYSquareNorm=0;
	outputSquareNorm=0;
}
template<class Real,int Type,int Degree,int Channels,class PartialType>
StreamingDivergence<Real,Type,Degree,Channels,PartialType>::~StreamingDivergence(void)
{
	if(localB)			free(localB);
	if(_localDMajor)	free(_localDMajor);
	if(_localDMinor)	free(_localDMinor);
	_localDMajor=_localDMinor=NULL;
	localDMajor=localDMinor=NULL;
	if(_localDMajorAccum)	free(_localDMajorAccum);
	if(_localDMinorAccum)	free(_localDMinorAccum);
	_localDMajorAccum=_localDMinorAccum=NULL;
	localDMajorAccum=localDMinorAccum=NULL;
}
template<class Real,int Type,int Degree,int Channels,class PartialType>
void StreamingDivergence<Real,Type,Degree,Channels,PartialType>::InitRestriction(int major,int minor,bool symmetric)
{
	index=0;
	this->major=major;
	this->minor=minor;
	sMajor=FiniteElements1D<Real,DERIVATIVE(Type),Degree-1>::Dimension(FiniteElements1D<Real,Type,Degree>::DomainSize(major));
	sMinor=FiniteElements1D<Real,DERIVATIVE(Type),Degree-1>::Dimension(FiniteElements1D<Real,Type,Degree>::DomainSize(minor));

	_major=(sMajor>major)?(sMajor+RealPerWord-1)/RealPerWord:(major+RealPerWord-1)/RealPerWord;

	if(_localDMajor)	free(_localDMajor);
	if(_localDMinor)	free(_localDMinor);
	_localDMajor=(WordClass*)malloc(sizeof(WordClass)*(dSize*_major*Channels*RealPerWord+2*Degree)+ALIGNMENT-1);
	_localDMinor=(WordClass*)malloc(sizeof(WordClass)*(dSize*_major*Channels*RealPerWord+2*Degree)+ALIGNMENT-1);
	memset(_localDMajor,0,sizeof(WordClass)*(dSize*_major*Channels*RealPerWord+2*Degree)+ALIGNMENT-1);
	memset(_localDMinor,0,sizeof(WordClass)*(dSize*_major*Channels*RealPerWord+2*Degree)+ALIGNMENT-1);
	localDMajor=(Real*)(((size_t)(_localDMajor+Degree)+ALIGNMENT-1) & ~(ALIGNMENT-1));
	localDMinor=(Real*)(((size_t)(_localDMinor+Degree)+ALIGNMENT-1) & ~(ALIGNMENT-1));
	if(_localDMajorAccum)	free(_localDMajorAccum);
	if(_localDMinorAccum)	free(_localDMinorAccum);
	_localDMajorAccum=(WordClass*)malloc(sizeof(WordClass)*(_major+2*Degree)+ALIGNMENT-1);
	_localDMinorAccum=(WordClass*)malloc(sizeof(WordClass)*(_major+2*Degree)+ALIGNMENT-1);
#if NEW_STREAMING_CODE
	localDMajorAccum=(WordClass*)(((size_t)(_localDMajorAccum+Degree)+ALIGNMENT-1) & ~(ALIGNMENT-1));
	localDMinorAccum=(WordClass*)(((size_t)(_localDMinorAccum+Degree)+ALIGNMENT-1) & ~(ALIGNMENT-1));
#else // !NEW_STREAMING_CODE
	localDMajorAccum=(WordClass*)(((size_t)(_localDMajorAccum)+ALIGNMENT-1) & ~(ALIGNMENT-1));
	localDMinorAccum=(WordClass*)(((size_t)(_localDMinorAccum)+ALIGNMENT-1) & ~(ALIGNMENT-1));
#endif // NEW_STREAMING_CODE

	FiniteElements1D<Real,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(major,dotMajorStencil,0,0);
	FiniteElements1D<Real,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(minor,dotMinorStencil,0,0);
	FiniteElements1D<Real,DERIVATIVE(Type),Degree-1>::DotProduct<Type,Degree>::DotProductStencil(FiniteElements1D<Real,DERIVATIVE(Type),Degree-1>::Dimension(FiniteElements1D<Real,Type,Degree>::DomainSize(major)),dDotMajorStencil,0,1);
	FiniteElements1D<Real,DERIVATIVE(Type),Degree-1>::DotProduct<Type,Degree>::DotProductStencil(FiniteElements1D<Real,DERIVATIVE(Type),Degree-1>::Dimension(FiniteElements1D<Real,Type,Degree>::DomainSize(minor)),dDotMinorStencil,0,1);
	FiniteElements2D<Real,Type,Degree>::DivergenceStencil(major,minor,divergenceStencil);

	int off1,off2;
	off1=FiniteElements1D<Real,DERIVATIVE(Type),Degree-1>::DotProduct<Type,Degree>::Helper::StopOffset();
	off2=FiniteElements1D<Real,Type,Degree>::DotProduct<Type,Degree>::Helper::StopOffset();

	startRestriction=(off1>off2)?off1:off2;
	startRestriction-=Degree;
	dMajor->reset(true,1);
	dMinor->reset(true,1);
	if(outB)
	{
		outB->reset(false,1);
		if(localB)	free(localB);
		localB=(Real*)malloc(sizeof(Real)*major*Channels);
	}
	if(rParent)	rParent->InitRestriction(symmetric);
}
template<class Real,int Type,int Degree,int Channels,class PartialType>
void StreamingDivergence<Real,Type,Degree,Channels,PartialType>::SetRestrictionIterations(int iters)
{
	if(rParent)	rParent->SetRestrictionIterations(iters);
}
#include "LaplacianMatrix/test.h"
template<class Real,int Type,int Degree,int Channels,class PartialType>
void StreamingDivergence<Real,Type,Degree,Channels,PartialType>::SetInteriorRestriction(Real* lB,int idx,int major2,int minor2)
{
	int dI=idx+FiniteElements1D<Real,DERIVATIVE(Type),Degree-1>::DotProduct<Type,Degree>::Helper::StartOffset();
	int I =idx+FiniteElements1D<Real,Type,Degree>::DotProduct<Type,Degree>::Helper::StartOffset();

	WordClass* majorPtrs[2*Degree+1];
	WordClass* minorPtrs[2*Degree  ];
	for(int c=0;c<Channels;c++)
	{
		Real* myLB=&lB[major2*c];
		for(int xx=0;xx<=2*Degree;xx++)	majorPtrs[xx]=(WordClass*)(&localDMajor[(( I+xx)%dSize)*_major*RealPerWord*Channels+c*_major*RealPerWord]);
		for(int xx=0;xx< 2*Degree;xx++)	minorPtrs[xx]=(WordClass*)(&localDMinor[((dI+xx)%dSize)*_major*RealPerWord*Channels+c*_major*RealPerWord]);
		/*
if(idx==128+testY && c==testC)
{
	printf("DX:\n");
	for(int y=Degree-1;y<=Degree+1;y++)
	{
		for(int x=-Degree+testX;x<=Degree+testX;x++)	printf("\t%.4f",((float*)majorPtrs[y])[128+x]);
		printf("\n");
	}
	printf("DY:\n");
	for(int y=Degree-1;y<=Degree;y++)
	{
		for(int x=-Degree+testX;x<=Degree+testX;x++)	printf("\t%.4f",((float*)minorPtrs[y])[128+x]);
		printf("\n");
	}
}
*/
#if USE_SSE_CODE
		memset(localDMinorAccum,0,sizeof(Real)*_major*RealPerWord);
		__declspec (align(16)) float scratch[4];
		__m128 dot[2*Degree+1],dDot[Degree];
		for(int d=0;d<=Degree;d++)
		{
			for(int j=0;j<4;j++)	scratch[j]= dotMinorStencil.caseTable[Degree].values[d];
			dot[d]=_mm_load_ps(scratch);
		}
		for(int d=0;d<Degree;d++)
		{
			for(int j=0;j<4;j++)	scratch[j]=dDotMinorStencil.caseTable[Degree].values[d];
			dDot[d]=_mm_load_ps(scratch);
		}
		{
			int d=Degree;
			for(int i=0;i<_major;i++)	localDMajorAccum[i]=_mm_mul_ps(dot[d],majorPtrs[d][i]);
		}
		for(int d=0;d<Degree;d++)
			for(int i=0;i<_major;i++)
			{
				localDMajorAccum[i]=_mm_add_ps(localDMajorAccum[i],_mm_mul_ps( dot[d],_mm_add_ps(majorPtrs[d][i],majorPtrs[2*Degree  -d][i])));
				localDMinorAccum[i]=_mm_add_ps(localDMinorAccum[i],_mm_mul_ps(dDot[d],_mm_sub_ps(minorPtrs[d][i],minorPtrs[2*Degree-1-d][i])));
			}
#else // !USE_SSE_CODE
		memset(localDMajorAccum,0,sizeof(Real)*_major*RealPerWord);
		memset(localDMinorAccum,0,sizeof(Real)*_major*RealPerWord);
		for(int d=0;d<Degree;d++)
			for(int i=0;i<_major*RealPerWord;i++)
			{
				((Real*)(localDMajorAccum))[i]+=(majorPtrs[d][i]+majorPtrs[2*Degree  -d][i])* dotMinorStencil.caseTable[Degree].values[d];
				((Real*)(localDMinorAccum))[i]+=(minorPtrs[d][i]-minorPtrs[2*Degree-1-d][i])*dDotMinorStencil.caseTable[Degree].values[d];
			}
			{
				int d=Degree;
				for(int i=0;i<_major*RealPerWord;i++)
					((Real*)(localDMajorAccum))[i]+=majorPtrs[d][i]* dotMinorStencil.caseTable[Degree].values[d];
			}
#endif // USE_SSE_CODE
			for(int j=0;j<major;j++)
			{
				Real temp=0;
				int jj;
				if(j<Degree)				jj=j;
				else if(j>=major-Degree)	jj=2*Degree+(j-(major-1));
				else						jj=Degree;
				int dJ=j+FiniteElements1D<Real,DERIVATIVE(Type),Degree-1>::DotProduct<Type,Degree>::Helper::StartOffset();
				int J =j+FiniteElements1D<Real,Type,Degree>::DotProduct<Type,Degree>::Helper::StartOffset();

				{
					// Partial w.r.t minor index
					const Real* localD=((Real*)(localDMinorAccum))+ J;
					for(int dj=0;dj<=2*Degree;dj++)
					{
						if(J+dj<0 || J+dj>=major)	continue;
						temp+=localD[dj]* dotMajorStencil.caseTable[jj].values[dj];
					}
				}
				{
					// Partial w.r.t major index
					const Real* localD=((Real*)(localDMajorAccum))+dJ;
					for(int dj=0;dj<2*Degree;dj++)
					{
						if(dJ+dj<0 || dJ+dj>=sMajor)	continue;
						temp+=localD[dj]*dDotMajorStencil.caseTable[jj].values[dj];
					}
				}
				myLB[j]=temp;
				outputSquareNorm+=temp*temp;
			}
	}
}

template<class Real,int Type,int Degree,int Channels,class PartialType>
void StreamingDivergence<Real,Type,Degree,Channels,PartialType>::SetRestriction(Real* lB,int idx,int major2,int minor2)
{
	if(idx>Degree && idx<minor-Degree)
	{
		SetInteriorRestriction(lB,idx,major2,minor2);
		return;
	}
	int ii;
	if(idx<Degree)				ii=idx;
	else if(idx>=minor-Degree)	ii=2*Degree+(idx-(minor-1));
	else						ii=Degree;
	int dI=idx+FiniteElements1D<Real,DERIVATIVE(Type),Degree-1>::DotProduct<Type,Degree>::Helper::StartOffset();
	int I =idx+FiniteElements1D<Real,Type,Degree>::DotProduct<Type,Degree>::Helper::StartOffset();
	for(int c=0;c<Channels;c++)
	{
		Real* myLB=&lB[c*major2];
		for(int j=0;j<major;j++)
		{
			Real temp=0;
			int jj;
			if(j<Degree)				jj=j;
			else if(j>=major-Degree)	jj=2*Degree+(j-(major-1));
			else						jj=Degree;
			int dJ=j+FiniteElements1D<Real,DERIVATIVE(Type),Degree-1>::DotProduct<Type,Degree>::Helper::StartOffset();
			int J =j+FiniteElements1D<Real,Type,Degree>::DotProduct<Type,Degree>::Helper::StartOffset();

			// Partial w.r.t minor index
			for(int di=0;di<2*Degree;di++)
			{
				if(dI+di<0 || dI+di>=sMinor)	continue;
				Real* localD=&localDMinor[((dI+di)%dSize)*_major*RealPerWord*Channels+c*_major*RealPerWord+J];
				for(int dj=0;dj<=2*Degree;dj++)
				{
					if(J+dj<0 || J+dj>=major)	continue;
					temp+=localD[dj]*divergenceStencil.caseTable[jj][ii].values2[dj][di];
				}
			}
			// Partial w.r.t major index
			for(int di=0;di<=2*Degree;di++)
			{
				if(I+di<0 || I+di>=minor)	continue;
				Real* localD=&localDMajor[((I+di)%dSize)*_major*RealPerWord*Channels+c*_major*RealPerWord+dJ];
				for(int dj=0;dj<2*Degree;dj++)
				{
					if(dJ+dj<0 || dJ+dj>=sMajor)	continue;
					temp+=localD[dj]*divergenceStencil.caseTable[jj][ii].values1[dj][di];
				}
			}
			myLB[j]=temp;
			outputSquareNorm+=temp*temp;
		}
	}
}
template<class Real,int Type,int Degree,int Channels,class PartialType>
bool StreamingDivergence<Real,Type,Degree,Channels,PartialType>::IterateRestriction(void)
{
	if(index>=minor && index>=sMinor)	return false;
	if(index<minor)
	{
		int columns=dMajor->rowSize()/(Channels*sizeof(PartialType));
		int rows=dMajor->rows();
		if(index<dMajor->rows())
		{
			PartialType* ptr2=(PartialType*)(*dMajor)[index];
			for(int c=0;c<Channels;c++)
			{
				Real* ptr1=&localDMajor[(index%dSize)*_major*RealPerWord*Channels+c*_major*RealPerWord];
				for(int i=0;i<columns;i++)
				{
					ptr1[i]=((Real)ptr2[i*Channels+c]);
					dXSquareNorm+=ptr1[i]*ptr1[i];
				}
#if REFLECT_GRADIENT
				if(columns<sMajor)	ptr1[columns]=0;
				for(int i=1;i<columns && columns+i<sMajor;i++)	ptr1[columns+i]=-((Real)ptr2[(columns-i)*Channels+c])/2;
#endif // REFLECT_GRADIENT
			}
#if REFLECT_GRADIENT
			if(rows-1-index<ReflectBufferSize)
			{
				int iidx=2*rows-1-index;
				// iidx<rows+ReflectBufferSize
				// 2*rows-1-index<rows+ReflectBufferSize
				// rows-1-index<ReflectBufferSize
				if(iidx<minor)
				{
					for(int c=0;c<Channels;c++)
					{
						Real* ptr1=&localDMajor[(iidx%dSize)*_major*RealPerWord*Channels+c*_major*RealPerWord];
						for(int i=0;i<columns;i++)	ptr1[i]=(Real)ptr2[i*Channels+c];
						if(columns<sMajor)	ptr1[columns]=0;
						for(int i=1;i<columns && columns+i<sMajor;i++)	ptr1[columns+i]=-(Real)ptr2[(columns-i)*Channels+c];
					}
				}
			}
#endif // REFLECT_GRADIENT
			dMajor->advance();
		}
#if REFLECT_GRADIENT
		else if(index>=rows+ReflectBufferSize)
#else // !REFLECT_GRADIENT
		else
#endif // REFLECT_GRADIENT
			for(int c=0;c<Channels;c++)
			{
				Real* ptr1=&localDMajor[(index%dSize)*_major*RealPerWord*Channels+c*_major*RealPerWord];
				for(int i=0;i<columns;i++)	ptr1[i]= (Real)0;
			}
	}
	if(index<sMinor)
	{
		int columns=dMinor->rowSize()/(Channels*sizeof(PartialType));
		int rows=dMinor->rows();
		if(index<dMinor->rows())
		{
			PartialType* ptr2=(PartialType*)(*dMinor)[index];
			for(int c=0;c<Channels;c++)
			{
				Real* ptr1=&localDMinor[(index%dSize)*_major*RealPerWord*Channels+c*_major*RealPerWord];
				for(int i=0;i<columns;i++)
				{
					ptr1[i]=(Real)ptr2[i*Channels+c];
					dYSquareNorm+=ptr1[i]*ptr1[i];
				}
#if REFLECT_GRADIENT
				for(int i=0;i<columns && columns+i<major;i++)	ptr1[columns+i]=((Real)ptr2[(columns-i-1)*Channels+c])/2;
#endif // REFLECT_GRADIENT
			}
#if REFLECT_GRADIENT
			if(rows-index<ReflectBufferSize)
			{
				int iidx=2*rows-index;
				// iidx<rows+ReflectBufferSize
				// 2*rows-index<rows+ReflectBufferSize
				// rows-index<ReflectBufferSize
				if(iidx<minor)
				{
					for(int c=0;c<Channels;c++)
					{
						Real* ptr1=&localDMinor[(iidx%dSize)*_major*RealPerWord*Channels+c*_major*RealPerWord];
						for(int i=0;i<columns;i++)	ptr1[i]=-(Real)ptr2[i*Channels+c];
						for(int i=0;i<columns && columns+i<major;i++)	ptr1[columns+i]=-(Real)ptr2[(columns-i-1)*Channels+c];
					}
				}
			}
#endif // REFLECT_GRADIENT
			dMinor->advance();
		}
#if REFLECT_GRADIENT
		else if(index>=rows+ReflectBufferSize || index==rows)
#else // !REFLECT_GRADIENT
		else
#endif // REFLECT_GRADIENT
			for(int c=0;c<Channels;c++)
			{
				Real* ptr1=&localDMinor[(index%dSize)*_major*RealPerWord*Channels+c*_major*RealPerWord];
				for(int i=0;i<columns;i++)	ptr1[i]= (Real)0;
			}
	}
	if(rParent && index>=startRestriction)	rParent->IterateRestriction();
	if(outB && index>=Degree)
	{
		SetRestriction(localB,index-Degree,major,minor);
		PartialType* temp=(PartialType*)(*outB)[index-Degree];
		for(int i=0;i<major*Channels;i++)	temp[i]=PartialType(localB[i]);
		outB->advance();
	}
	index++;
	return true;
}
template<class Real,int Type,int Degree,int Channels,class PartialType>
void StreamingDivergence<Real,Type,Degree,Channels,PartialType>::SolveRestriction(void)
{
	// Run to completion...
	while( IterateRestriction() ){;}
	// Write out the rest of the Laplacian
	if(outB)
	{
		while(index-Degree<minor)
		{
			SetRestriction(localB,index-Degree,major,minor);
			PartialType* temp=(PartialType*)(*outB)[index-Degree];
			for(int i=0;i<major*Channels;i++)	temp[i]=PartialType(localB[i]);
			outB->advance();
			index++;
		}
	}
	// ...and finish up the trailing parent
	if(rParent)	rParent->SolveRestriction();
}
////////////////////////
// StreamingLaplacian //
////////////////////////
template<class Real,int Type,int Degree,int Channels,class InputType>
StreamingLaplacian<Real,Type,Degree,Channels,InputType>::StreamingLaplacian(void)
{
	_localDMajor=_localDMinor=NULL;
	localDMajor=localDMinor=NULL;
	_localDMajorAccum=NULL;
	_localDMinorAccum=NULL;
	localDMajorAccum=NULL;
	localDMinorAccum=NULL;
	pixels=labels=NULL;
	dSize=2*Degree+1;
	dXSquareNorm=0;
	dYSquareNorm=0;
	outputSquareNorm=0;
	_previousPixelsRow=_previousLabelsRow=NULL;
	_dx=_dy=NULL;
}
template<class Real,int Type,int Degree,int Channels,class InputType>
StreamingLaplacian<Real,Type,Degree,Channels,InputType>::~StreamingLaplacian(void)
{
	if(_localDMajor)	free(_localDMajor);
	if(_localDMinor)	free(_localDMinor);
	_localDMajor=_localDMinor=NULL;
	localDMajor=localDMinor=NULL;
	if(_localDMajorAccum)	free(_localDMajorAccum);
	if(_localDMinorAccum)	free(_localDMinorAccum);
	_localDMajorAccum=_localDMinorAccum=NULL;
	localDMajorAccum=localDMinorAccum=NULL;
	if(_previousPixelsRow)	free(_previousPixelsRow);
	if(_previousLabelsRow)	free(_previousLabelsRow);
	_previousPixelsRow=_previousLabelsRow=NULL;
	if(_dx)	free(_dx);
	if(_dy)	free(_dy);
	_dx=_dy=NULL;
}
template<class Real,int Type,int Degree,int Channels,class InputType>
void StreamingLaplacian<Real,Type,Degree,Channels,InputType>::InitRestriction(int major,int minor,bool symmetric,bool pad)
{
	_pad=pad;
	_hh=_h=pixels->rows();
	_ww=_w=pixels->rowSize()/(Channels*sizeof(InputType));
	if(pixels->rowSize()!=_w*Channels*sizeof(InputType))	fprintf(stderr,"Pixel width failure: %d != %d\n",pixels->rowSize(),_w*Channels*sizeof(InputType)),	exit(0);
	if(_h!=labels->rows())									fprintf(stderr,"Label height failure: %d != %d\n",_h,labels->rows()),								exit(0);
	if(labels->rowSize()!=_w*Channels*sizeof(InputType))	fprintf(stderr,"Label width failure: %d != %d\n",labels->rowSize(),_w*Channels*sizeof(InputType)),	exit(0);

	if(_pad)	_hh++	,	_ww++;
	_previousPixelsRow	= (InputType*)malloc(sizeof(InputType)*Channels*_w);
	_previousLabelsRow	= (InputType*)malloc(sizeof(InputType)*Channels*_w);
	_dx					= (Real*)malloc(sizeof(Real)*Channels*(_ww-1));
	_dy					= (Real*)malloc(sizeof(Real)*Channels*(_ww  ));

	_current=0;
	pixels->reset(true,1);
	labels->reset(true,1);

	memcpy(_previousPixelsRow,(*pixels)[_current],sizeof(InputType)*Channels*_w);
	memcpy(_previousLabelsRow,(*labels)[_current],sizeof(InputType)*Channels*_w);
	for(int c=0;c<Channels;c++)	average[c]=0;
	for(int x=0;x<_w;x++)	for(int c=0;c<Channels;c++)	average[c]+=_previousPixelsRow[x+_w*c];
	InputType l1[Channels],l2[Channels];

	// Set the partials with respect to x
	for(int x=0;x<_w-1;x++)
	{
		bool useGradient=true;
		for(int c=0;c<Channels;c++)	l1[c] = _previousLabelsRow[x+_w*c]	,	l2[c] = _previousLabelsRow[x+1+_w*c];
		for(int c=0;c<Channels && useGradient;c++)	useGradient &= (l1[c]==l2[c]);
		if(useGradient)		for(int c=0;c<Channels;c++)	_dx[x*Channels+c] = Real(_previousPixelsRow[x+1+_w*c])-Real(_previousPixelsRow[x+_w*c]);
		else				for(int c=0;c<Channels;c++)	_dx[x*Channels+c] = 0;
	}
	if(_pad)	for(int c=0;c<Channels;c++)	_dx[(_w-1)*Channels+c] = Real(0)-Real(_previousPixelsRow[(_w-1)+_w*c]);
	pixels->advance();
	labels->advance();
	InputType* newPixelsRow = (InputType*)(*pixels)[_current+1];
	InputType* newLabelsRow = (InputType*)(*labels)[_current+1];

	// Set the partials with respect to y
	for(int x=0;x<_w;x++)
	{
		bool useGradient=true;
		for(int c=0;c<Channels;c++)	l1[c] = newLabelsRow[x+_w*c]	,	l2[c] = _previousLabelsRow[x+_w*c];
		for(int c=0;c<Channels && useGradient;c++)	useGradient &= (l1[c]==l2[c]);
		if(useGradient)		for(int c=0;c<Channels;c++)	_dy[x*Channels+c] = Real(newPixelsRow[x+_w*c])-Real(_previousPixelsRow[x+_w*c]);
		else				for(int c=0;c<Channels;c++)	_dy[x*Channels+c] = 0;
	}
	if(_pad)	for(int c=0;c<Channels;c++)	_dy[_w*Channels+c] = Real(0);

	memcpy(_previousPixelsRow,newPixelsRow,sizeof(InputType)*Channels*_w);
	memcpy(_previousLabelsRow,newLabelsRow,sizeof(InputType)*Channels*_w);
	for(int x=0;x<_w;x++)	for(int c=0;c<Channels;c++)	average[c]+=_previousPixelsRow[x+_w*c];
///////////

	index=0;
	this->major=major;
	this->minor=minor;
	sMajor=FiniteElements1D<Real,DERIVATIVE(Type),Degree-1>::Dimension(FiniteElements1D<Real,Type,Degree>::DomainSize(major));
	sMinor=FiniteElements1D<Real,DERIVATIVE(Type),Degree-1>::Dimension(FiniteElements1D<Real,Type,Degree>::DomainSize(minor));

	_major=(sMajor>major)?(sMajor+RealPerWord-1)/RealPerWord:(major+RealPerWord-1)/RealPerWord;

	if(_localDMajor)	free(_localDMajor);
	if(_localDMinor)	free(_localDMinor);
	_localDMajor=(WordClass*)malloc(sizeof(WordClass)*(dSize*_major*Channels*RealPerWord+2*Degree)+ALIGNMENT-1);
	_localDMinor=(WordClass*)malloc(sizeof(WordClass)*(dSize*_major*Channels*RealPerWord+2*Degree)+ALIGNMENT-1);
	memset(_localDMajor,0,sizeof(WordClass)*(dSize*_major*Channels*RealPerWord+2*Degree)+ALIGNMENT-1);
	memset(_localDMinor,0,sizeof(WordClass)*(dSize*_major*Channels*RealPerWord+2*Degree)+ALIGNMENT-1);
	localDMajor=(Real*)(((size_t)(_localDMajor+Degree)+ALIGNMENT-1) & ~(ALIGNMENT-1));
	localDMinor=(Real*)(((size_t)(_localDMinor+Degree)+ALIGNMENT-1) & ~(ALIGNMENT-1));
	if(_localDMajorAccum)	free(_localDMajorAccum);
	if(_localDMinorAccum)	free(_localDMinorAccum);
	_localDMajorAccum=(WordClass*)malloc(sizeof(WordClass)*(_major+2*Degree)+ALIGNMENT-1);
	_localDMinorAccum=(WordClass*)malloc(sizeof(WordClass)*(_major+2*Degree)+ALIGNMENT-1);
#if NEW_STREAMING_CODE
	localDMajorAccum=(WordClass*)(((size_t)(_localDMajorAccum+Degree)+ALIGNMENT-1) & ~(ALIGNMENT-1));
	localDMinorAccum=(WordClass*)(((size_t)(_localDMinorAccum+Degree)+ALIGNMENT-1) & ~(ALIGNMENT-1));
#else // !NEW_STREAMING_CODE
	localDMajorAccum=(WordClass*)(((size_t)(_localDMajorAccum)+ALIGNMENT-1) & ~(ALIGNMENT-1));
	localDMinorAccum=(WordClass*)(((size_t)(_localDMinorAccum)+ALIGNMENT-1) & ~(ALIGNMENT-1));
#endif // NEW_STREAMING_CODE

	FiniteElements1D<Real,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(major,dotMajorStencil,0,0);
	FiniteElements1D<Real,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(minor,dotMinorStencil,0,0);
	FiniteElements1D<Real,DERIVATIVE(Type),Degree-1>::DotProduct<Type,Degree>::DotProductStencil(FiniteElements1D<Real,DERIVATIVE(Type),Degree-1>::Dimension(FiniteElements1D<Real,Type,Degree>::DomainSize(major)),dDotMajorStencil,0,1);
	FiniteElements1D<Real,DERIVATIVE(Type),Degree-1>::DotProduct<Type,Degree>::DotProductStencil(FiniteElements1D<Real,DERIVATIVE(Type),Degree-1>::Dimension(FiniteElements1D<Real,Type,Degree>::DomainSize(minor)),dDotMinorStencil,0,1);
	FiniteElements2D<Real,Type,Degree>::DivergenceStencil(major,minor,divergenceStencil);

	int off1,off2;
	off1=FiniteElements1D<Real,DERIVATIVE(Type),Degree-1>::DotProduct<Type,Degree>::Helper::StopOffset();
	off2=FiniteElements1D<Real,Type,Degree>::DotProduct<Type,Degree>::Helper::StopOffset();

	startRestriction=(off1>off2)?off1:off2;
	startRestriction-=Degree;
	if(rParent)	rParent->InitRestriction(symmetric);
}
template<class Real,int Type,int Degree,int Channels,class InputType>
void StreamingLaplacian<Real,Type,Degree,Channels,InputType>::SetRestrictionIterations(int iters)
{
	if(rParent)	rParent->SetRestrictionIterations(iters);
}
template<class Real,int Type,int Degree,int Channels,class InputType>
void StreamingLaplacian<Real,Type,Degree,Channels,InputType>::SetInteriorRestriction(Real* lB,int idx,int major2,int minor2)
{
	int dI=idx+FiniteElements1D<Real,DERIVATIVE(Type),Degree-1>::DotProduct<Type,Degree>::Helper::StartOffset();
	int I =idx+FiniteElements1D<Real,Type,Degree>::DotProduct<Type,Degree>::Helper::StartOffset();

	WordClass* majorPtrs[2*Degree+1];
	WordClass* minorPtrs[2*Degree  ];
	for(int c=0;c<Channels;c++)
	{
		Real* myLB=&lB[major2*c];
		for(int xx=0;xx<=2*Degree;xx++)	majorPtrs[xx]=(WordClass*)(&localDMajor[(( I+xx)%dSize)*_major*RealPerWord*Channels+c*_major*RealPerWord]);
		for(int xx=0;xx< 2*Degree;xx++)	minorPtrs[xx]=(WordClass*)(&localDMinor[((dI+xx)%dSize)*_major*RealPerWord*Channels+c*_major*RealPerWord]);
#if USE_SSE_CODE
		memset(localDMinorAccum,0,sizeof(Real)*_major*RealPerWord);
		__declspec (align(16)) float scratch[4];
		__m128 dot[2*Degree+1],dDot[Degree];
		for(int d=0;d<=Degree;d++)
		{
			for(int j=0;j<4;j++)	scratch[j]= dotMinorStencil.caseTable[Degree].values[d];
			dot[d]=_mm_load_ps(scratch);
		}
		for(int d=0;d<Degree;d++)
		{
			for(int j=0;j<4;j++)	scratch[j]=dDotMinorStencil.caseTable[Degree].values[d];
			dDot[d]=_mm_load_ps(scratch);
		}
		{
			int d=Degree;
			for(int i=0;i<_major;i++)	localDMajorAccum[i]=_mm_mul_ps(dot[d],majorPtrs[d][i]);
		}
		for(int d=0;d<Degree;d++)
			for(int i=0;i<_major;i++)
			{
				localDMajorAccum[i]=_mm_add_ps(localDMajorAccum[i],_mm_mul_ps( dot[d],_mm_add_ps(majorPtrs[d][i],majorPtrs[2*Degree  -d][i])));
				localDMinorAccum[i]=_mm_add_ps(localDMinorAccum[i],_mm_mul_ps(dDot[d],_mm_sub_ps(minorPtrs[d][i],minorPtrs[2*Degree-1-d][i])));
			}
#else // !USE_SSE_CODE
		memset(localDMajorAccum,0,sizeof(Real)*_major*RealPerWord);
		memset(localDMinorAccum,0,sizeof(Real)*_major*RealPerWord);
		for(int d=0;d<Degree;d++)
			for(int i=0;i<_major*RealPerWord;i++)
			{
				((Real*)(localDMajorAccum))[i]+=(majorPtrs[d][i]+majorPtrs[2*Degree  -d][i])* dotMinorStencil.caseTable[Degree].values[d];
				((Real*)(localDMinorAccum))[i]+=(minorPtrs[d][i]-minorPtrs[2*Degree-1-d][i])*dDotMinorStencil.caseTable[Degree].values[d];
			}
			{
				int d=Degree;
				for(int i=0;i<_major*RealPerWord;i++)
					((Real*)(localDMajorAccum))[i]+=majorPtrs[d][i]* dotMinorStencil.caseTable[Degree].values[d];
			}
#endif // USE_SSE_CODE
			for(int j=0;j<major;j++)
			{
				Real temp=0;
				int jj;
				if(j<Degree)				jj=j;
				else if(j>=major-Degree)	jj=2*Degree+(j-(major-1));
				else						jj=Degree;
				int dJ=j+FiniteElements1D<Real,DERIVATIVE(Type),Degree-1>::DotProduct<Type,Degree>::Helper::StartOffset();
				int J =j+FiniteElements1D<Real,Type,Degree>::DotProduct<Type,Degree>::Helper::StartOffset();

				{
					// Partial w.r.t minor index
					const Real* localD=((Real*)(localDMinorAccum))+ J;
					for(int dj=0;dj<=2*Degree;dj++)
					{
						if(J+dj<0 || J+dj>=major)	continue;
						temp+=localD[dj]* dotMajorStencil.caseTable[jj].values[dj];
					}
				}
				{
					// Partial w.r.t major index
					const Real* localD=((Real*)(localDMajorAccum))+dJ;
					for(int dj=0;dj<2*Degree;dj++)
					{
						if(dJ+dj<0 || dJ+dj>=sMajor)	continue;
						temp+=localD[dj]*dDotMajorStencil.caseTable[jj].values[dj];
					}
				}
				myLB[j]=temp;
				outputSquareNorm+=temp*temp;
			}
	}
}

template<class Real,int Type,int Degree,int Channels,class InputType>
void StreamingLaplacian<Real,Type,Degree,Channels,InputType>::SetRestriction(Real* lB,int idx,int major2,int minor2)
{
	if(idx>Degree && idx<minor-Degree)
	{
		SetInteriorRestriction(lB,idx,major2,minor2);
		return;
	}
	int ii;
	if(idx<Degree)				ii=idx;
	else if(idx>=minor-Degree)	ii=2*Degree+(idx-(minor-1));
	else						ii=Degree;
	int dI=idx+FiniteElements1D<Real,DERIVATIVE(Type),Degree-1>::DotProduct<Type,Degree>::Helper::StartOffset();
	int I =idx+FiniteElements1D<Real,Type,Degree>::DotProduct<Type,Degree>::Helper::StartOffset();
	for(int c=0;c<Channels;c++)
	{
		Real* myLB=&lB[c*major2];
		for(int j=0;j<major;j++)
		{
			Real temp=0;
			int jj;
			if(j<Degree)				jj=j;
			else if(j>=major-Degree)	jj=2*Degree+(j-(major-1));
			else						jj=Degree;
			int dJ=j+FiniteElements1D<Real,DERIVATIVE(Type),Degree-1>::DotProduct<Type,Degree>::Helper::StartOffset();
			int J =j+FiniteElements1D<Real,Type,Degree>::DotProduct<Type,Degree>::Helper::StartOffset();

			// Partial w.r.t minor index
			for(int di=0;di<2*Degree;di++)
			{
				if(dI+di<0 || dI+di>=sMinor)	continue;
				Real* localD=&localDMinor[((dI+di)%dSize)*_major*RealPerWord*Channels+c*_major*RealPerWord+J];
				for(int dj=0;dj<=2*Degree;dj++)
				{
					if(J+dj<0 || J+dj>=major)	continue;
					temp+=localD[dj]*divergenceStencil.caseTable[jj][ii].values2[dj][di];
				}
			}
			// Partial w.r.t major index
			for(int di=0;di<=2*Degree;di++)
			{
				if(I+di<0 || I+di>=minor)	continue;
				Real* localD=&localDMajor[((I+di)%dSize)*_major*RealPerWord*Channels+c*_major*RealPerWord+dJ];
				for(int dj=0;dj<2*Degree;dj++)
				{
					if(dJ+dj<0 || dJ+dj>=sMajor)	continue;
					temp+=localD[dj]*divergenceStencil.caseTable[jj][ii].values1[dj][di];
				}
			}
			myLB[j]=temp;
			outputSquareNorm+=temp*temp;
		}
	}
}
template<class Real,int Type,int Degree,int Channels,class InputType>
bool StreamingLaplacian<Real,Type,Degree,Channels,InputType>::IterateRestriction(void)
{
///////////////

	if(index>=minor && index>=sMinor)	return false;
	if(index<minor)
	{
		if(index<_hh)
			for(int c=0;c<Channels;c++)
			{
				Real* ptr=&localDMajor[(index%dSize)*_major*RealPerWord*Channels+c*_major*RealPerWord];
				for(int i=0;i<_ww-1;i++)
				{
					ptr[i]=_dx[i*Channels+c];
					dXSquareNorm+=ptr[i]*ptr[i];
				}
			}
		else
			for(int c=0;c<Channels;c++)
			{
				Real* ptr=&localDMajor[(index%dSize)*_major*RealPerWord*Channels+c*_major*RealPerWord];
				for(int i=0;i<_ww-1;i++)	ptr[i]= (Real)0;
			}
	}
	if(index<sMinor)
	{
		if(index<_hh-1)
			for(int c=0;c<Channels;c++)
			{
				Real* ptr=&localDMinor[(index%dSize)*_major*RealPerWord*Channels+c*_major*RealPerWord];
				for(int i=0;i<_ww;i++)
				{
					ptr[i]=_dy[i*Channels+c];
					dYSquareNorm+=ptr[i]*ptr[i];
				}
			}
		else
			for(int c=0;c<Channels;c++)
			{
				Real* ptr=&localDMinor[(index%dSize)*_major*RealPerWord*Channels+c*_major*RealPerWord];
				for(int i=0;i<_ww;i++)	ptr[i]= (Real)0;
			}
	}
	_current++;
	if(_current<_hh)
	{
		InputType l1[Channels],l2[Channels];
		// Set the partials with respect to x
		for(int x=0;x<_w-1;x++)
		{
			bool useGradient=true;
			for(int c=0;c<Channels;c++)	l1[c] = _previousLabelsRow[x+_w*c]	,	l2[c] = _previousLabelsRow[x+1+_w*c];
			for(int c=0;c<Channels && useGradient;c++)	useGradient &= (l1[c]==l2[c]);
			if(useGradient)		for(int c=0;c<Channels;c++)	_dx[x*Channels+c] = Real(_previousPixelsRow[x+1+_w*c])-Real(_previousPixelsRow[x+_w*c]);
			else				for(int c=0;c<Channels;c++)	_dx[x*Channels+c] = 0;
		}
		if(_pad)	for(int c=0;c<Channels;c++)	_dx[(_w-1)*Channels+c] = Real(0)-Real(_previousPixelsRow[(_w-1)+_w*c]);

		if(_current<_h)
		{
			pixels->advance();
			labels->advance();
		}
		if(_current+1<_h)
		{
			InputType* newPixelsRow = (InputType*)(*pixels)[_current+1];
			InputType* newLabelsRow = (InputType*)(*labels)[_current+1];

			// Set the partials with respect to y
			for(int x=0;x<_w;x++)
			{
				bool useGradient=true;
				for(int c=0;c<Channels;c++)	l1[c] = newLabelsRow[x+_w*c]	,	l2[c] = _previousLabelsRow[x+_w*c];
				for(int c=0;c<Channels && useGradient;c++)	useGradient &= (l1[c]==l2[c]);
				if(useGradient)		for(int c=0;c<Channels;c++)	_dy[x*Channels+c] = Real(newPixelsRow[x+_w*c])-Real(_previousPixelsRow[x+_w*c]);
				else				for(int c=0;c<Channels;c++)	_dy[x*Channels+c] = 0;
			}
			if(_pad)	for(int c=0;c<Channels;c++)	_dy[_w*Channels+c] = 0;

			memcpy(_previousPixelsRow,newPixelsRow,sizeof(InputType)*Channels*_w);
			memcpy(_previousLabelsRow,newLabelsRow,sizeof(InputType)*Channels*_w);
			for(int x=0;x<_w;x++)	for(int c=0;c<Channels;c++)	average[c]+=_previousPixelsRow[x+_w*c];
		}
		else if(_current+1==_h && _pad)
		{
			// Set the partials with respect to y
			for(int x=0;x<_w;x++)	for(int c=0;c<Channels;c++)	_dy[x*Channels+c] = Real(0)-Real(_previousPixelsRow[x+_w*c]);
			for(int c=0;c<Channels;c++)	_dy[_w*Channels+c] = Real(0);

			memset(_previousPixelsRow,0,sizeof(InputType)*Channels*_w);
			memset(_previousLabelsRow,0,sizeof(InputType)*Channels*_w);
		}
	}
	if(_current==_hh-1)	for(int c=0;c<Channels;c++)	average[c]/=_w	,	average[c]/=_h;


	if(rParent)	if(index>=startRestriction)	rParent->IterateRestriction();
	index++;
	return true;
}
template<class Real,int Type,int Degree,int Channels,class InputType>
void StreamingLaplacian<Real,Type,Degree,Channels,InputType>::SolveRestriction(void)
{
	// Run to completion...
	while( IterateRestriction() ){;}
	// ...and finish up the trailing parent
	if(rParent)	rParent->SolveRestriction();
}
