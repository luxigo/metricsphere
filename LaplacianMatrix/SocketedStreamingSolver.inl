
#define MyFree( buffer ) if( buffer ) free( buffer ) , buffer = NULL

void SetDownSampledStencil( int dim , int iters , DotProductStencil& outDot , DotProductStencil& outD2Dot )
{
	DotProductStencil dot , d2Dot;
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil( dim ,   dot , 0 , 0 );
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil( dim , d2Dot , 1 , 1 , false );
	SetDownSampledStencil( dot , d2Dot , dim , iters , outDot , outD2Dot );

}
void SetDownSampledStencil( const DotProductStencil& inDot , const DotProductStencil& inD2Dot,
						    int dim , int iters ,
							DotProductStencil& outDot , DotProductStencil& outD2Dot )
{

	if( !iters )
	{
		outDot = inDot;
		outD2Dot = inD2Dot;
		return;
	}
	else
	{
		DotProductStencil newDot , newD2Dot;
		FiniteElements1D<double,Type,Degree>::FullProlongationStencil prolongationStencil;
		FiniteElements1D<double,Type,Degree>::ProlongationStencil( dim >>1 , prolongationStencil , dim );
		CombineStencils<double,Type,Degree>(    inDot , prolongationStencil , dim ,   newDot );
		CombineStencils<double,Type,Degree>(  inD2Dot , prolongationStencil , dim , newD2Dot );
		SetDownSampledStencil( newDot , newD2Dot , dim>>1 , iters-1 , outDot , outD2Dot );
	}
}
////////////////////////
// SocketedStreamData //
////////////////////////

#if SUPPORT_CYLINDRICAL
bool SocketedStreamData::_set( int width , int height , int start , int end , DataStream* leftStream , SOCKET syncSocket , DataStream* rightStream , int periodicType )
#else // !SUPPORT_CYLINDRICAL
bool SocketedStreamData::_set(int width,int height,int start,int end,DataStream* leftStream,SOCKET syncSocket,DataStream* rightStream,bool spherical)
#endif // SUPPORT_CYLINDRICAL
{
	if(width&3)														{ fprintf(stderr,"width must be a multiple of four: %d\n",width) ; return false; }
	if(start&3)														{ fprintf(stderr,"start must be a multiple of four: %d\n",start) ; return false; }
	if(end&3)														{ fprintf(stderr,"end must be a multiple of four: %d\n",end)	 ; return false; }
#if SUPPORT_CYLINDRICAL
	if( (start || periodicType!=NO_PERIODIC ) && !leftStream )			{ fprintf(stderr,"left socket invalid for stream\t%dx%d\n",width,height) ; return false; }
	if( (end!=width || periodicType!=NO_PERIODIC ) && !rightStream )	{ fprintf(stderr,"right socket invalid for stream\t%dx%d\n",width,height)  ; return false; }
	if( periodicType==SPHERICAL_PERIODIC && syncSocket==INVALID_SOCKET ){ fprintf(stderr,"sync socket invalid for spherical stream\n") ; return false; }
#else // !SUPPORT_CYLINDRICAL
	if((start || spherical) && !leftStream)							{ fprintf(stderr,"left socket invalid for non-start stream\t%dx%d\n",width,height) ; return false; }
	if((end!=width || spherical) && !rightStream)					{ fprintf(stderr,"right socket invalid for non-end stream\t%dx%d\n",width,height)  ; return false; }
	if(spherical && syncSocket==INVALID_SOCKET)						{ fprintf(stderr,"sync socket invalid for spherical stream\n") ; return false; }
#endif // SUPPORT_CYLINDRICAL

	this->leftStream = leftStream;
	this->rightStream =rightStream;
	syncXSocket = syncRSocket = syncSocket;

	_start128=start>>2;
	_end128=end	>>2;
	_size128=_end128-_start128;

	_start=_start128<<2;
	_end=_end128<<2;
	_size=_size128<<2;

	_paddedSize128=_size128+2*_padSize128;
	_padSize=_padSize128<<2;
	_paddedSize=_paddedSize128<<2;

	return true;
}

// This guys gets called by the streaming solver

bool SocketedStreamData::SetSocketedStreamData(int width,int height,int start,int end,int iters,
											   DataStream* leftStream,SOCKET syncSocket,DataStream* rightStream,
#if SUPPORT_CYLINDRICAL
											   int periodicType )
#else // !SUPPORT_CYLINDRICAL
											   bool spherical)
#endif // SUPPORT_CYLINDRICAL
{
	// It's iters+2 because we need:
	//	1] iters-1 for offsetting the solver
	//	2] 1 for reading in the solver
	//	3] 1 for computing the residual without having to re-sync.
	//	4] 1 for zig-zag since you start one in.
	_padSize128=(iters+2)*WordPerDegree;

#if SUPPORT_CYLINDRICAL
	if( !_set( width , height , start , end , leftStream , syncSocket , rightStream , periodicType ) )	return false;
#else // !SUPPORT_CYLINDRICAL
	if(!_set(width,height,start,end,leftStream,syncSocket,rightStream , spherical))	return false;
#endif // SUPPORT_CYLINDRICAL

	if(end-start<4*_padSize)	{ fprintf(stderr,"stream too narrow: %d < 4*%d\n",end-start,_padSize) ; return false; }
	if(height<2*Degree)			{ fprintf(stderr,"stream too short: %d < 2*%d\n",height,Degree)		  ; return false; }
	return true;
}
// This guys gets called by the streaming laplacian
bool SocketedStreamData::SetSocketedStreamData(int width,int height,int start,int end,
											   DataStream* leftStream,SOCKET syncSocket,DataStream* rightStream,
#if SUPPORT_CYLINDRICAL
											   int periodicType )
#else // !SUPPORT_CYLINDRICAL
											   bool spherical)
#endif // SUPPORT_CYLINDRICAL
{
	_padSize128=WordPerDegree;

#if SUPPORT_CYLINDRICAL
	if( !_set( width , height , start , end , leftStream , syncSocket , rightStream , periodicType ) )	return false;
#else // !SUPPORT_CYLINDRICAL
	if(!_set(width,height,start,end,leftStream,syncSocket,rightStream,spherical))	return false;
#endif // SUPPORT_CYLINDRICAL

	if(end-start<2*_padSize)	{ fprintf(stderr,"stream too narrow: %d < 2*%d\n",end-start,_padSize) ; return false; }
	return true;
}
int SocketedStreamData::start	(void)	const	{ return _start; }
int SocketedStreamData::end		(void)	const	{ return _end; }
int SocketedStreamData::size	(void)	const	{ return _size; }

/////////////////////////////
// SocketedStreamingSolver //
/////////////////////////////
template< int Channels , class SyncType > int SocketedStreamingSolver< Channels , SyncType >::OffsetX(int iters)	{	return (iters+1)*Degree;	}
template< int Channels , class SyncType > int SocketedStreamingSolver< Channels , SyncType >::OffsetB(int iters)	{	return iters*Degree;	}
//template<int Channels>	int SocketedStreamingSolver< Channels , SyncType >::OffsetR(void)	{	return -1;	}
template< int Channels , class SyncType > int SocketedStreamingSolver< Channels , SyncType >::OffsetR(void)		{	return -3;	}

template< int Channels , class SyncType >
SocketedStreamingSolver< Channels , SyncType >::SocketedStreamingSolver(void)
{
	_server = NULL;
	_deleteServer = true;
#if TIME_IO
	vSync = hSync = rSync = 0;
#endif // TIME_IO
	showProgress=false;
	progressCount=0;

	setResidual=true;
	laplacianRescale=true;
	xSize = bSize = rSize = 0;
	laneNum = 1;

	RStream = NullPointer< Pointer( __m128 ) >( );
	XStream = NullPointer< Pointer( __m128 ) >( );
	BStream = NullPointer< Pointer( __m128 ) >( );
	for( int d=0 ; d<Degree*Channels ; d++ ) RBuffer[d] = XBuffer[d] = NullPointer< __m128 >( );

	lapTemplates = AllocArray< TemplateSSE >( 3 * (2*Degree+1) , ALIGNMENT , "SocketedStreamingSolver::lapTemplates" );
	zeroLapTemplates = AllocArray< TemplateSSE >( 3 * (2*Degree+1)  , ALIGNMENT , "SocketedStreamingSolver::zeroLapTemplates" );

	syncBuffer = NullPointer< SyncType >( );
	for( int i=0 ; i<Degree ; i++ ) localXAccum[i] = NullPointer<  __m128 >( );

}
template< int Channels , class SyncType >
SocketedStreamingSolver< Channels , SyncType >::~SocketedStreamingSolver(void)
{
	freeStreams();

	FreeArray( syncBuffer );
	FreeArray( lapTemplates );
	FreeArray( zeroLapTemplates );
	for( int i=0 ; i<Degree ; i++ )
		if( localXAccum[i] )
		{
			localXAccum[i] -= _padSize128;
			FreeArray( localXAccum[i] );
		}
	if( _server && _deleteServer ) delete _server , _server=NULL;

}
template< int Channels , class SyncType >
void SocketedStreamingSolver< Channels , SyncType >::freeStreams(void)
{
	if( XStream )
	{
		for( int i=0 ; i<xSize*Channels ; i++ )
			if( XStream[i] )
			{
				XStream[i] -= _padSize128;
				FreeArray( XStream[i] );
			}
		FreeArray( XStream );
	}
	if( BStream )
	{
		for( int i=0 ; i<bSize*Channels ; i++ )
			if( BStream[i] )
			{
				BStream[i] -= _padSize128;
				FreeArray( BStream[i] );
			}
		FreeArray( BStream );
	}
	if( RStream )
	{
		for( int i=0 ; i<rSize*Channels ; i++ )
			if( RStream[i] )
			{
				RStream[i] -= _padSize128;
				FreeArray( RStream[i] );
			}
		FreeArray( RStream );
	}
	if( RBuffer ) for( int i=0 ; i<Degree*Channels ; i++ ) FreeArray( RBuffer[i] );
	if( XBuffer ) for( int i=0 ; i<Degree*Channels ; i++ ) FreeArray( XBuffer[i] );

	RStream = NullPointer< Pointer( __m128 ) >( );
	XStream = NullPointer< Pointer( __m128 ) >( );
	BStream = NullPointer< Pointer( __m128 ) >( );

	rSize=xSize=bSize=0;
}
template< int Channels , class SyncType >
void SocketedStreamingSolver< Channels , SyncType >::Init( int start , int end , int major , int minor , int iters , 
											 DataStream* leftStream , SOCKET syncSocket , DataStream* rightStream ,
#if SUPPORT_CYLINDRICAL
											 int periodicType , MultiStreamIOServer* server
#else // !SUPPORT_CYLINDRICAL
											 bool spherical , MultiStreamIOServer* server
#endif // SUPPORT_CYLINDRICAL
											 )
{
	MatrixStencil lStencil;
	FiniteElements2D<double,Type,Degree>::LaplacianStencil(major,minor,lStencil,true);
	Init( lStencil , start , end , major , minor , iters ,
		leftStream , syncSocket , rightStream ,
#if SUPPORT_CYLINDRICAL
		periodicType , server
#else // !SUPPORT_CYLINDRICAL
		spherical , server
#endif // SUPPORT_CYLINDRICAL
		);
}
template< int Channels , class SyncType >
void SocketedStreamingSolver< Channels , SyncType >::Init(const MatrixStencil& lStencil,int start,int end,int major,int minor,int iters,
											 DataStream* leftStream,SOCKET syncSocket,DataStream* rightStream,
#if SUPPORT_CYLINDRICAL
											 int periodicType , MultiStreamIOServer* server
#else // !SUPPORT_CYLINDRICAL
											 bool spherical , MultiStreamIOServer* server
#endif // SUPPORT_CYLINDRICAL
											 )
{
#if SUPPORT_CYLINDRICAL
	if( !SetSocketedStreamData( major , minor , start , end , iters , leftStream , syncSocket , rightStream , periodicType ) )	exit(0);
#else // !SUPPORT_CYLINDRICAL
	if(!SetSocketedStreamData(major,minor,start,end,iters,leftStream,syncSocket,rightStream , spherical ) )	exit(0);
#endif // SUPPORT_CYLINDRICAL
	if( _server && _deleteServer ) delete _server , _server = NULL;
	if( server ) _server = server , _deleteServer = false;
	else		 _server = new MultiStreamIOServer() , _deleteServer = true;

#if SUPPORT_CYLINDRICAL
	this->periodicType = periodicType;
#else // !SUPPORT_CYLINDRICAL
	this->spherical=spherical;
#endif // SUPPORT_CYLINDRICAL
	this->major=major;
	this->minor=minor;
	this->iters=iters;
#if SUPPORT_CYLINDRICAL
	if( periodicType==SPHERICAL_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
	if( spherical )
#endif // SUPPORT_CYLINDRICAL
	{
		FreeArray( syncBuffer );
		syncBuffer = AllocArray< SyncType >( _size * Channels * Degree , 1 , "SocketedStreamingSolver::syncBuffer" );
	}

	MatrixStencil laplacianStencil;
	laplacianStencil=lStencil;
	laplacianScale =1.f/laplacianStencil.caseTable[Degree][Degree].values[Degree][Degree];
	laplacianScaleR=    laplacianStencil.caseTable[Degree][Degree].values[Degree][Degree];

	for(int i=0;i<=2*Degree;i++)
		for(int j=0;j<=2*Degree;j++)
			for(int k=0;k<=2*Degree;k++)
				for(int l=0;l<=2*Degree;l++)
					laplacianStencil.caseTable[i][j].values[k][l]*=laplacianScale;
	// BADNESS!!! should fix the indexing of lapTemplates so that major index iterates faster
	__declspec (align(16)) float scratch[4];
	for(int i=0;i<=2*Degree;i++)	// Iterate over the minor index in the mask
	{
		for(int k=0;k<4;k++)
			lapTemplates[3*i+1].diagonal[k]=1.f/laplacianStencil.caseTable[i][Degree].values[Degree][Degree];
		lapTemplates[3*i  ].diagonal[0]=1.f/laplacianStencil.caseTable[i][0].values[Degree][Degree];
		lapTemplates[3*i  ].diagonal[1]=1.f/laplacianStencil.caseTable[i][1].values[Degree][Degree];
		lapTemplates[3*i  ].diagonal[2]=1.f/laplacianStencil.caseTable[i][Degree].values[Degree][Degree];
		lapTemplates[3*i  ].diagonal[3]=1.f/laplacianStencil.caseTable[i][Degree].values[Degree][Degree];
		lapTemplates[3*i+2].diagonal[0]=1.f/laplacianStencil.caseTable[i][Degree].values[Degree][Degree];
		lapTemplates[3*i+2].diagonal[1]=1.f/laplacianStencil.caseTable[i][Degree].values[Degree][Degree];
		lapTemplates[3*i+2].diagonal[2]=1.f/laplacianStencil.caseTable[i][2*Degree-1].values[Degree][Degree];
		lapTemplates[3*i+2].diagonal[3]=1.f/laplacianStencil.caseTable[i][2*Degree  ].values[Degree][Degree];
		for(int k=0;k<4;k++)
			zeroLapTemplates[3*i+1].diagonal[k]=1.f/laplacianStencil.caseTable[i][Degree].values[Degree][Degree];
		zeroLapTemplates[3*i  ].diagonal[0]=1.f/laplacianStencil.caseTable[i][0].values[Degree][Degree];
		zeroLapTemplates[3*i  ].diagonal[1]=1.f/laplacianStencil.caseTable[i][1].values[Degree][Degree];
		zeroLapTemplates[3*i  ].diagonal[2]=1.f/laplacianStencil.caseTable[i][Degree].values[Degree][Degree];
		zeroLapTemplates[3*i  ].diagonal[3]=1.f/laplacianStencil.caseTable[i][Degree].values[Degree][Degree];
		zeroLapTemplates[3*i+2].diagonal[0]=1.f/laplacianStencil.caseTable[i][Degree].values[Degree][Degree];
		zeroLapTemplates[3*i+2].diagonal[1]=1.f/laplacianStencil.caseTable[i][Degree].values[Degree][Degree];
		zeroLapTemplates[3*i+2].diagonal[2]=1.f/laplacianStencil.caseTable[i][2*Degree-1].values[Degree][Degree];
		zeroLapTemplates[3*i+2].diagonal[3]=1.f/laplacianStencil.caseTable[i][2*Degree  ].values[Degree][Degree];
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
				lapTemplates[3*i+j].matrixValues[0][k]=_mm_load_ps(scratch);

				if(j==0)	jj=1;
				else		jj=Degree;
				for(int l=0;l<4;l++)	scratch[l]=laplacianStencil.caseTable[i][jj].values[k][l+1];
				lapTemplates[3*i+j].matrixValues[1][k]=_mm_load_ps(scratch);

				if(j==2)	jj=2*Degree-1;
				else		jj=Degree;
				for(int l=0;l<4;l++)	scratch[l]=laplacianStencil.caseTable[i][jj].values[k][l];
				lapTemplates[3*i+j].matrixValues[2][k]=_mm_load_ps(scratch);

				if(j==2)	jj=2*Degree;
				else		jj=Degree;
				if(j!=0)	scratch[0]=laplacianStencil.caseTable[i][Degree].values[k][3];
				else		scratch[0]=0;
				scratch[1]=laplacianStencil.caseTable[i][jj].values[k][0];
				scratch[2]=laplacianStencil.caseTable[i][jj].values[k][1];
				scratch[3]=laplacianStencil.caseTable[i][jj].values[k][2];
				lapTemplates[3*i+j].matrixValues[3][k]=_mm_load_ps(scratch);
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
				zeroLapTemplates[3*i+j].matrixValues[0][k]=_mm_load_ps(scratch);

				if(j==0)	jj=1;
				else		jj=Degree;
				for(int l=0;l<4;l++)	scratch[l]=laplacianStencil.caseTable[i][jj].values[k][l+1];
				zeroLapTemplates[3*i+j].matrixValues[1][k]=_mm_load_ps(scratch);

				if(j==2)	jj=2*Degree-1;
				else		jj=Degree;
				for(int l=0;l<4;l++)	scratch[l]=laplacianStencil.caseTable[i][jj].values[k][l];
				zeroLapTemplates[3*i+j].matrixValues[2][k]=_mm_load_ps(scratch);

				if(j==2)	jj=2*Degree;
				else		jj=Degree;
				if(j!=0)	scratch[0]=laplacianStencil.caseTable[i][Degree].values[k][3];
				else		scratch[0]=0;
				scratch[1]=laplacianStencil.caseTable[i][jj].values[k][0];
				scratch[2]=laplacianStencil.caseTable[i][jj].values[k][1];
				scratch[3]=laplacianStencil.caseTable[i][jj].values[k][2];
				zeroLapTemplates[3*i+j].matrixValues[3][k]=_mm_load_ps(scratch);
			}
}
template< int Channels , class SyncType >
void SocketedStreamingSolver< Channels , SyncType >::Set(int rSize)
{
	Set(-OffsetX,0,0,0,0,rSize);
}
template< int Channels , class SyncType >
void SocketedStreamingSolver< Channels , SyncType >::Set(int start,int bStart,int xStart,int bEnd,int xEnd,int rSize)
{
	// At index idx, we solve for the X coefficients
	// [idx+Degree-1,idx+2*Degree-1,...,idx+iters*Degree-1]
	// Which requires read access to the X coefficients:
	// [idx-1,idx+(iters+1)*Degree-1]
	// and the B coefficients
	// [idx+Degree-1,idx+iters*Degree-1]

	// After solving at idx, the values of X at coefficient idx+Degree-1 (and below) are finalized.
	// Which means that we can solve for the value of the residual in row idx-1, requiring read access to the X coefficients:
	// [idx-Degree-1,idx+Degree-1]
	// and the B coefficients
	// [idx-1]

	clearX=clearB=true;
	index=start;
	freeStreams();

	xSize += ( OffsetX(iters) > xEnd ? OffsetX(iters) : xEnd ) + 1;
	bSize += ( OffsetB(iters) > bEnd ? OffsetB(iters) : bEnd ) + 1;
	if(setResidual)
	{
		xSize+=-Degree+OffsetR()<xStart?Degree-OffsetR():-xStart;
		bSize+=        OffsetR()<bStart?      -OffsetR():-bStart;
	}
	else
	{
		xSize += -Degree<xStart ? Degree : - xStart;
		bSize +=       0<bStart ?      0 : - bStart;
	}
#if SAME_SIZE_BUFFERS
	bSize = xSize;
#endif // SAME_SIZE_BUFFERS
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

	this->rSize=rSize;
	if(setResidual)
	{
		RStream = AllocArray< Pointer( __m128 ) >( rSize * Channels , 1 , "SocketedStreamingSolver::RStream" );
		for( int i=0 ; i<rSize*Channels ; i++ )
		{
			RStream[ i ]  = AllocArray< __m128 >( _paddedSize128 , ALIGNMENT , "SocketedStreamingSolver::RStream[ i ]" );
			RStream[ i ] += _padSize128;
		}
	}
	XStream = AllocArray< Pointer( __m128 ) >( xSize * Channels , 1  , "SocketedStreamingSolver::XStream" );
	for( int i=0 ; i<xSize*Channels ; i++ )
	{
		XStream[ i ]  = AllocArray< __m128 >( _paddedSize128 , ALIGNMENT , "SocketedStreamingSolver::XStream[ i ]" );
		memset( XStream[i] , 0 , sizeof( __m128 ) * _paddedSize128 );
		XStream[ i ] += _padSize128;
	}

	BStream = AllocArray< Pointer( __m128 ) >( bSize * Channels , 1  , "SocketedStreamingSolver::BStream" );
	for( int i=0 ; i<bSize*Channels ; i++ )
	{
		BStream[ i ] = AllocArray< __m128 >( _paddedSize128 , ALIGNMENT , "SocketedStreamingSolver::BStream[ i ]" );
		BStream[ i ] += _padSize128;
	}

	for( int i=0 ; i<Degree ; i++ )
	{
		localXAccum[i] = AllocArray< __m128 >( _paddedSize128 , ALIGNMENT , "SocketedStreamingSolver::localXAccum[ i ]" );
		localXAccum[i] += _padSize128;
	}

#if SUPPORT_CYLINDRICAL
	if( periodicType==SPHERICAL_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
	if(spherical)
#endif // SUPPORT_CYLINDRICAL
	{
		for( int d=0 ; d<Degree*Channels ; d++ ) XBuffer[d] = AllocArray< __m128 >(  _paddedSize128 , ALIGNMENT , "SocketedStreamingSolver::XBuffer" );
		if( setResidual )
			for( int d=0 ; d<Degree*Channels ; d++ ) RBuffer[d] = AllocArray< __m128 >( _paddedSize128 , ALIGNMENT , "SocketedStreamingSolver::RBuffer" );
	}
}
template< int Channels , class SyncType >
Pointer( float ) SocketedStreamingSolver< Channels , SyncType >::GetXRow( int row , int channel )
{
	if( row<0 )
#if SUPPORT_CYLINDRICAL
		if( periodicType==SPHERICAL_PERIODIC && row >=-Degree)		return ( Pointer( float ) )( XBuffer[(Degree+row)*Channels+channel] + _padSize128 );
#else // !SUPPORT_CYLINDRICAL
		if(spherical && row >=-Degree)		return ( Pointer( float ) )( XBuffer[(Degree+row)*Channels+channel] + _padSize128 );
#endif // SUPPORT_CYLINDRICAL
		else								return NullPointer< float >( );
	else if(row>=minor)
#if SUPPORT_CYLINDRICAL
		if( periodicType==SPHERICAL_PERIODIC && row < minor+Degree)	return ( Pointer( float ) )( XBuffer[(row-minor)*Channels+channel] + _padSize128 );
#else // !SUPPORT_CYLINDRICAL
		if(spherical && row < minor+Degree)	return ( Pointer( float ) )( XBuffer[(row-minor)*Channels+channel] + _padSize128 );
#endif // SUPPORT_CYLINDRICAL
		else								return NullPointer< float >( );
	else									return ( Pointer( float ) )XStream[ MyModIndex( row*Channels+channel , xSize*Channels ) ];
}
template< int Channels , class SyncType >
Pointer( float ) SocketedStreamingSolver< Channels , SyncType >::GetBRow(int row,int channel)
{
	return ( Pointer( float ) )BStream[MyModIndex(row*Channels+channel , bSize*Channels)];
}
template< int Channels , class SyncType >
Pointer( float ) SocketedStreamingSolver< Channels , SyncType >::GetRRow(int row,int channel)
{
	if(row<0)
#if SUPPORT_CYLINDRICAL
		if( periodicType==SPHERICAL_PERIODIC && row >=-Degree)		return ( Pointer( float ) )( RBuffer[(Degree+row)*Channels + channel] + _padSize128 );
#else // !SUPPORT_CYLINDRICAL
		if(spherical && row >=-Degree)		return ( Pointer( float ) )( RBuffer[(Degree+row)*Channels + channel] + _padSize128 );
#endif // SUPPORT_CYLINDRICAL
		else								return NullPointer< float >( );
	else if(row>=minor)
#if SUPPORT_CYLINDRICAL
		if( periodicType==SPHERICAL_PERIODIC && row < minor+Degree)	return ( Pointer( float ) )( RBuffer[(row-minor)*Channels + channel] + _padSize128 );
#else // !SUPPORT_CYLINDRICAL
		if(spherical && row < minor+Degree)	return ( Pointer( float ) )( RBuffer[(row-minor)*Channels + channel] + _padSize128 );
#endif // SUPPORT_CYLINDRICAL
		else								return NullPointer< float >( );
	else									return ( Pointer( float ) )RStream[MyModIndex(row*Channels+channel,rSize*Channels)];
}
template< int Channels , class SyncType >
void SocketedStreamingSolver< Channels , SyncType >::UnSet(void)
{
	clearX=clearB=true;
	index=0;
//	iters=0;
	setResidual=false;

	freeStreams();
	for( int i=0 ; i<Degree ; i++ )
		if( localXAccum[i] )
		{
			localXAccum[i] -= _padSize128;
			FreeArray( localXAccum[i] );
		}
}
template< int Channels , class SyncType >
bool SocketedStreamingSolver< Channels , SyncType >::Increment(void)
{
	int idx;
	if( clearX )
	{
		if(setResidual)	idx=index-Degree+OffsetR()+xSize;
		else			idx=index-Degree          +xSize;
		if( idx>=0 && idx<minor ) for( int c=0 ; c<Channels ; c++ ) memset( GetXRow( idx , c ) , 0 , sizeof( __m128 ) * _size128 );
	}
	if( clearB )
	{
#if SAME_SIZE_BUFFERS
		if(setResidual)	idx=index-Degree+OffsetR()+bSize;
		else			idx=index-Degree          +bSize;
#else // !SAME_SIZE_BUFFERS
		if(setResidual)	idx=index+OffsetR()+bSize;
		else			idx=index          +bSize;
#endif // SAME_SIZE_BUFFERS
		if( idx>=0 && idx<minor ) for( int c=0 ; c<Channels ; c++ ) memset( GetBRow( idx , c ) , 0 , sizeof( __m128 ) * _size128 );
	}
	index++;
	int restrictionIndex  = FiniteElements1D< float , Type , Degree >::FullRestrictionStencil::RestrictionStencil::Start( index+OffsetR() );
	int prolongationIndex = FiniteElements1D< float , Type , Degree >::FullProlongationStencil::ProlongationStencil::Start( index-1 );
	if( restrictionIndex>=minor/2 && index>=minor && prolongationIndex>=minor*2 && (!setResidual || (index+OffsetR()>=minor) ) ) return false;
	else return true;
}
template< int Channels , class SyncType >
template<class StorageType>
bool SocketedStreamingSolver< Channels , SyncType >::UpdateXInput(StreamingGrid* X)
{
	if( !X ) return false;
	int idx = index+OffsetX(iters);
	// Read in from the X vector
	if(idx>=0 && idx<minor)
	{
		Pointer( StorageType ) xPtr = ( Pointer( StorageType ) )(*X)[idx];
		for(int c=0;c<Channels;c++)
		{
			int idx2=_size*c;
			Pointer( float ) xRow = GetXRow(idx,c);
			for(int jj=0;jj<_size;jj++)	xRow[jj]+=float(xPtr[jj+idx2]);
		}
		return true;
	}
	return false;
}
template< int Channels , class SyncType >
template<class StorageType>
bool SocketedStreamingSolver< Channels , SyncType >::UpdateBInput(StreamingGrid* B)
{
	if( !B ) return false;
	int idx = index+OffsetB(iters);
	// Read in from the B vector
	if(idx>=0 && idx<minor)
	{
		Pointer( StorageType ) bPtr = ( Pointer( StorageType ) )(*B)[idx];
		for(int c=0;c<Channels;c++)
		{
			int idx2=_size*c;
			Pointer( float ) bRow = GetBRow( idx , c );
			for(int jj=0;jj<_size;jj++)	bRow[jj]=float(bPtr[jj+idx2]);
		}
		return true;
	}
	return false;
}

template< int Channels , class SyncType >
template<class StorageType>
bool SocketedStreamingSolver< Channels , SyncType >::UpdateXOutput(StreamingGrid* X)
{
	if( !X ) return false;
	// Copy the solution
	if(index>=0 && index<minor)
	{
		Pointer( StorageType ) xPtr = ( Pointer( StorageType ) )(*X)[index];
		for(int c=0;c<Channels;c++)
		{
			int idx2=_size*c;
			Pointer( float ) xRow = GetXRow( index , c );
#if 0	// Sanity check to ensure that the input and output are bounded
			for(int jj=0;jj<_size;jj++)
			{

				xPtr[jj+idx2]=StorageType(xRow[jj]);
				if( !isfinitef(xRow[jj]) || !isfinitef(float(xPtr[jj+idx2])) )
					printf("Infinity creaped in: %f -> %f\n",xRow[jj],float(xPtr[jj+idx2]));
			}
#else
			for(int jj=0;jj<_size;jj++)	xPtr[jj+idx2]=StorageType(xRow[jj]);
#endif
		}
		return true;
	}
	return false;
}
template< int Channels , class SyncType >
template<class StorageType>
bool SocketedStreamingSolver< Channels , SyncType >::UpdateBOutput(StreamingGrid* B)
{
	if( !B ) return false;
	if(index>=0 && index<minor)
	{
		Pointer( StorageType ) bPtr = ( Pointer( StorageType ) )(*B)[index];
		for(int c=0;c<Channels;c++)
		{
			int idx2=_size*c;
			Pointer( float ) bRow = GetBRow( index , c );
			for(int jj=0;jj<_size;jj++)	bPtr[jj+idx2]=StorageType(bRow[jj]);
		}
		return true;
	}
	return false;
}
template< int Channels , class SyncType >
void SocketedStreamingSolver< Channels , SyncType >::SyncSolverHead( int idx , bool read )
{
#if DEBUG_FLAGS
	if( debugFlags.noTopBottom ) return;
#endif // DEBUG_FLAGS
	if( idx<0 || idx>=Degree )
	{
		fprintf( stderr , "Synchronizing solver in non-head row: %d\n" , idx );
		return;
	}
	int r = -idx-1;
	if( syncXSocket!=INVALID_SOCKET )
		if( read )
		{
#if TIME_IO
			double t = Time();
#endif // TIME_IO
			ReceiveOnSocket( syncXSocket , syncBuffer , sizeof( SyncType )*_size*Channels , "SocketedStreamingSolver< Channels , SyncType >::SyncSolverHead" );
#if TIME_IO
			vSync += Time() - t;
#endif // TIME_IO
			for( int c=0 ; c<Channels ; c++ )
			{
				Pointer( float ) xRow = GetXRow( r , c );
				Pointer( SyncType ) sRow = syncBuffer + _size*c;
				for( int i=0 ; i<_size ; i++ ) xRow[i] = float( sRow[i] );
			}
		}
		else
		{
			for( int c=0 ; c<Channels ; c++ )
			{
				Pointer( float ) xRow = GetXRow( idx , c );
				Pointer( SyncType ) sRow = syncBuffer+_size*c;
				for( int i=0 ; i<_size ; i++ ) sRow[i] = SyncType( xRow[i] );
			}
			SendOnSocket( syncXSocket , syncBuffer , sizeof( SyncType )*_size*Channels , "SocketedStreamingSolver< Channels , SyncType >::SyncSolverHead" );
		}
}
template< int Channels , class SyncType >
void SocketedStreamingSolver< Channels , SyncType >::SyncSolverTail( int idx , bool read )
{
#if DEBUG_FLAGS
	if( debugFlags.noTopBottom ) return;
#endif // DEBUG_FLAGS
	if( idx<minor-Degree || idx>=minor )
	{
		fprintf( stderr , "Synchronizing solver in non-tail row: %d\n" , idx );
		return;
	}
	int r = 2 * minor - 1 - idx;
	if( syncXSocket!=INVALID_SOCKET )
		if( read )
		{
#if TIME_IO
			double t = Time();
#endif // TIME_IO
			ReceiveOnSocket( syncXSocket , syncBuffer , sizeof( SyncType )*_size*Channels , "SocketedStreamingSolver< Channels , SyncType >::SyncSolverTail" );
#if TIME_IO
			vSync += Time() - t;
#endif // TIME_IO
			for( int c=0 ; c<Channels ; c++ )
			{
				Pointer( float ) xRow = GetXRow( r , c );
				Pointer( SyncType ) sRow = syncBuffer+_size*c;
				for( int i=0 ; i<_size ; i++ ) xRow[i] = float( sRow[i] );
			}
		}
		else
		{
			for( int c=0 ; c<Channels ; c++ )
			{
				Pointer( float ) xRow = GetXRow( idx , c );
				Pointer( SyncType ) sRow = syncBuffer+_size*c;
				for( int i=0 ; i<_size ; i++ ) sRow[i] = SyncType( xRow[i] );
			}
			SendOnSocket( syncXSocket , syncBuffer , sizeof( SyncType )*_size*Channels , "SocketedStreamingSolver< Channels , SyncType >::SyncSolverTail" );
		}
}

template< int Channels , class SyncType >
void SocketedStreamingSolver< Channels , SyncType >::SyncResidualHead( int idx , bool read , bool syncPeriodic )
{
#if DEBUG_FLAGS
	if( debugFlags.noTopBottom ) return;
#endif // DEBUG_FLAGS
	if( idx<0 || idx>=Degree )
	{
		fprintf( stderr , "Synchronizing residual in non-head row: %d\n" , idx );
		return;
	}
	int r = -idx-1;
	if( syncRSocket!=INVALID_SOCKET )
	{
		if( !syncPeriodic )
			if( read )
			{
#if TIME_IO
				double t = Time();
#endif // TIME_IO
				ReceiveOnSocket( syncRSocket , syncBuffer , sizeof( SyncType )*_size*Channels , "SocketedStreamingSolver< Channels , SyncType >::SyncResidualHead" );
#if TIME_IO
				vSync += Time() - t;
#endif // TIME_IO
				for( int c=0 ; c<Channels ; c++ )
				{
					Pointer( float ) rRow = GetRRow( r , c );
					Pointer( SyncType ) sBuffer = syncBuffer + _size*c;
					for( int i=0 ; i<_size ; i++ ) rRow[i] = float( sBuffer[i] );
				}
			}
			else
			{
				for( int c=0 ; c<Channels ; c++ )
				{
					Pointer( float ) rRow = GetRRow( idx , c );
					Pointer( SyncType ) sBuffer = syncBuffer + _size*c;
					for( int i=0 ; i<_size ; i++ ) sBuffer[i] = SyncType( rRow[i] );
				}
				SendOnSocket( syncRSocket , syncBuffer , sizeof( SyncType )*_size*Channels , "SocketedStreamingSolver< Channels , SyncType >::SyncResidualHead" );
			}
		else
		{
			int size = WordPerDegree*RealPerWord;
			if( read )
			{
#if TIME_IO
				double t;
				t = Time();
#endif // TIME_IO
				leftStream->read( ( Pointer( byte ) )syncBuffer , sizeof( SyncType )*size*Channels );
#if TIME_IO
				hSync += Time() - t;
#endif // TIME_IO
				for( int c=0 ; c<Channels ; c++ )
				{
					Pointer( float ) rRow = GetRRow( r , c ) - size;
					Pointer( SyncType ) sBuffer = syncBuffer + size*c;
					for( int i=0 ; i<size ; i++ ) rRow[i] = float( sBuffer[i] );
				}
#if TIME_IO
				t = Time();
#endif // TIME_IO
				rightStream->read( ( Pointer( byte ) )syncBuffer , sizeof( SyncType )*size*Channels );
#if TIME_IO
				hSync += Time() - t;
#endif // TIME_IO
				for( int c=0 ; c<Channels ; c++ )
				{
					Pointer( float ) rRow = GetRRow( r  , c ) + _size;
					Pointer( SyncType ) sBuffer = syncBuffer + size*c;
					for( int i=0 ; i<size ; i++ ) rRow[i] = float( sBuffer[i] );
				}
			}
			else
			{
				for( int c=0 ; c<Channels ; c++ )
				{
					Pointer( float ) rRow = GetRRow( r , c );
					Pointer( SyncType ) sBuffer = syncBuffer + size*c;
					for( int i=0 ; i<size ; i++ ) sBuffer[i] = SyncType( rRow[i] );
				}
				leftStream->write( ( Pointer( byte ) )syncBuffer , sizeof( SyncType )*size*Channels );
				for( int c=0 ; c<Channels ; c++ )
				{
					Pointer( float ) rRow = GetRRow( r , c ) + _size - size;
					Pointer( SyncType ) sBuffer = syncBuffer + size*c;
					for( int i=0 ; i<size ; i++ ) sBuffer[i] = SyncType( rRow[i] );
				}
				rightStream->write( ( Pointer( byte ) )syncBuffer , sizeof( SyncType )*size*Channels );
			}
		}
	}
}
template< int Channels , class SyncType >
void SocketedStreamingSolver< Channels , SyncType >::SyncResidualTail( int idx , bool read , bool syncPeriodic )
{
#if DEBUG_FLAGS
	if( debugFlags.noTopBottom ) return;
#endif // DEBUG_FLAGS
	if( idx<minor-Degree || idx>=minor )
	{
		fprintf( stderr , "Synchronizing residual in non-tail row: %d\n" , idx );
		return;
	}
	int r = 2 * minor - 1 - idx;
	if( syncRSocket!=INVALID_SOCKET )
	{
		if( !syncPeriodic )
			if( read )
			{
#if TIME_IO
				double t = Time();
#endif // TIME_IO
				ReceiveOnSocket( syncRSocket , syncBuffer , sizeof( SyncType ) * _size * Channels , "SocketedStreamingSolver< Channels , SyncType >::SyncResidualTail" );
#if TIME_IO
				vSync += Time() - t;
#endif // TIME_IO
				for( int c=0 ; c<Channels ; c++ )
				{
					Pointer( float ) rRow = GetRRow( r , c );
					Pointer( SyncType ) sBuffer = syncBuffer + _size*c;
					for( int i=0 ; i<_size ; i++ ) rRow[i] = float( sBuffer[i] );
				}
			}
			else
			{
				for( int c=0 ; c<Channels ; c++ )
				{
					Pointer( float ) rRow = GetRRow( idx , c );
					Pointer( SyncType ) sBuffer = syncBuffer + _size*c;
					for( int i=0 ; i<_size ; i++ ) sBuffer[i] = SyncType( rRow[i] );
				}
				SendOnSocket( syncRSocket , syncBuffer , sizeof(SyncType)*_size*Channels , "SocketedStreamingSolver< Channels , SyncType >::SyncResidualTail" );
			}
		else
		{
			int size = WordPerDegree*RealPerWord;
			if(read)
			{
#if TIME_IO
				double t = Time();
#endif // TIME_IO
				leftStream->read( ( Pointer( byte ) )syncBuffer , sizeof(SyncType) * size * Channels );
#if TIME_IO
				hSync += Time() - t;
#endif // TIME_IO
				for( int c=0 ; c<Channels ; c++ )
				{
					Pointer( float ) rRow = GetRRow( r , c ) - size;
					Pointer( SyncType ) sBuffer = syncBuffer + size*c;
					for( int i=0 ; i<size ; i++ ) rRow[i] = float( sBuffer[i] );
				}
#if TIME_IO
				t = Time();
#endif // TIME_IO
				rightStream->read( ( Pointer( byte ) )syncBuffer , sizeof(SyncType)*size*Channels );
#if TIME_IO
				hSync += Time() - t;
#endif // TIME_IO
				for( int c=0 ; c<Channels ; c++)
				{
					Pointer( float ) rRow = GetRRow( r , c ) + _size;
					Pointer( SyncType ) sBuffer = syncBuffer + size*c;
					for( int i=0 ; i<size ; i++ ) rRow[i] = float( sBuffer[i] );
				}
			}
			else
			{
				for( int c=0 ; c<Channels ; c++ )
				{
					Pointer( float ) rRow = GetRRow( r , c );
					Pointer( SyncType ) sBuffer = syncBuffer + size*c;
					for( int i=0 ; i<size ; i++ ) sBuffer[i] = SyncType( rRow[i] );
				}
				leftStream->write( ( Pointer( byte ) )syncBuffer , sizeof(SyncType)*size*Channels );
				for( int c=0 ; c<Channels ; c++ )
				{
					Pointer( float ) rRow = GetRRow( r , c ) + _size - size;
					Pointer( SyncType ) sBuffer = syncBuffer + size*c;
					for( int i=0 ; i<size ; i++ ) sBuffer[i] = SyncType( rRow[i] );
				}
				rightStream->write( ( Pointer( byte ) )syncBuffer , sizeof(SyncType)*size*Channels );
			}
		}
	}
}

template< int Channels , class SyncType >
void SocketedStreamingSolver< Channels , SyncType >::SyncSolverLeft( int j , bool read , bool overlapped )
{
#if DEBUG_FLAGS
	if( debugFlags.noLeftRight ) return;
#endif // DEBUG_FLAGS
	if( leftStream )
	{
		int xSize = _padSize;
		int offset = _padSize-2*WordPerDegree*RealPerWord;
		int ySize = (iters+1)*Degree+2;
		if(!iters)	ySize+=Degree;
		int size = (xSize+offset) * ySize * Channels + 2 * xSize * Channels;
		int yEnd = Degree+1;
		Pointer( SyncType ) left = AllocArray< SyncType >( size , 1 , "SocketedStreamingSolver::SyncSolverLeft (left)" );
		if( read )
		{
#if TIME_IO
			double t = Time();
#endif // TIME_IO
#if DEBUG_FLAGS
			if( debugFlags.emptyPackets )
			{
				if ( !leftStream->read( ( Pointer( byte ) )left , sizeof( SyncType ) ) ) exit(0);
				FreeArray( left );
				return;
			}
#endif // DEBUG_FLAGS
			if ( !leftStream->read( ( Pointer( byte ) )left , sizeof( SyncType )*size ) )	exit(0);
#if TIME_IO
			hSync += Time() - t;
#endif // TIME_IO
			for( int y=-(ySize-yEnd-1) ; y<=yEnd ; y++ )
			{
				int yy = y + (ySize-yEnd-1);
#if SUPPORT_CYLINDRICAL
				if( (j+y>=0 && j+y<minor) || periodicType==SPHERICAL_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
				if( (j+y>=0 && j+y<minor) || spherical)
#endif // SUPPORT_CYLINDRICAL
					for(int c=0;c<Channels;c++)
					{
						Pointer( float ) tmp = GetXRow( j+y , c );
						if( tmp ) for( int x=0 ; x<xSize ; x++ ) tmp[ -xSize + x ] = float( left[ x * ySize * Channels + yy * Channels + c ] );
					}
			}
			if( overlapped )
				for( int k=0 ; k<iters ; k++ )
				{
					int y = -k*Degree;
					int yy = y + (ySize-yEnd-1);
					if( j+y>=0 && j+y<minor )
						for(int c=0;c<Channels;c++)
						{
							Pointer( float ) tmp = GetXRow( j+y , c );
							if( tmp ) for( int x=0 ; x<offset-WordPerDegree*RealPerWord*k ; x++ ) tmp[ x ] = float( left[ (xSize+x) * ySize * Channels + yy * Channels + c ] );
						}
				}
			for(int c=0;c<Channels;c++)
			{
				Pointer( float ) tmp = GetBRow( j , c );
				if( tmp ) for( int x=0 ; x<xSize ; x++ ) tmp[ -xSize+x ] = float( left[ (xSize+offset) * ySize * Channels + c * xSize + x ] );
			}
#if SUPPORT_CYLINDRICAL
			if( (j+1<minor) || periodicType==SPHERICAL_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if( (j+1<minor) || spherical)
#endif // SUPPORT_CYLINDRICAL
				for(int c=0;c<Channels;c++)
				{
					Pointer( float ) tmp = GetBRow( j+1 , c );
					if( tmp ) for( int x=0 ; x<xSize ; x++ ) tmp[ -xSize+x ] = float( left[ (xSize+offset) * ySize * Channels + Channels * xSize + c * xSize + x ] );
				}
		}
		else
		{
#if DEBUG_FLAGS
			if( debugFlags.emptyPackets )
			{
				if ( !leftStream->write( ( Pointer( byte ) )left , sizeof( SyncType ) ) ) exit(0);
				FreeArray( left );
				return;
			}
#endif // DEBUG_FLAGS
			memset( left , 0 , sizeof( SyncType )*size );
			for( int y=-(ySize-yEnd-1) ; y<=yEnd ; y++ )
			{
				int yy = y + (ySize-yEnd-1);
#if SUPPORT_CYLINDRICAL
				if( (j+y>=0 && j+y<minor) || periodicType==SPHERICAL_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
				if( (j+y>=0 && j+y<minor) || spherical)
#endif // SUPPORT_CYLINDRICAL
					for(int c=0;c<Channels;c++)
					{
						Pointer( float ) tmp = GetXRow( j+y , c );
						if( tmp ) for( int x=0 ; x<xSize+offset ; x++ ) left[ x * ySize * Channels + yy * Channels + c ] = SyncType( tmp[ x-offset ] );
					}
			}
			for(int c=0;c<Channels;c++)
			{
				Pointer( float ) tmp = GetBRow( j , c );
				if( tmp ) for( int x=0 ; x<xSize ; x++ ) left[ (xSize+offset) * ySize * Channels + c * xSize + x ] = SyncType( tmp[ x ] );
			}
#if SUPPORT_CYLINDRICAL
			if( (j+1<minor) || periodicType==SPHERICAL_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if( (j+1<minor) || spherical)
#endif // SUPPORT_CYLINDRICAL
				for(int c=0;c<Channels;c++)
				{
					Pointer( float ) tmp = GetBRow( j+1 , c );
					if( tmp ) for( int x=0 ; x<xSize ; x++ ) left[ (xSize+offset) * ySize * Channels + Channels * xSize + c * xSize + x ] = SyncType( tmp[ x ] );
				}
			if ( !leftStream->write( ( Pointer( byte ) )left , sizeof( SyncType )*size ) )	exit(0);
		}
		FreeArray( left );
	}
}
template< int Channels , class SyncType >
void SocketedStreamingSolver< Channels , SyncType >::SyncSolverRight( int j , bool read , bool overlapped )
{
#if DEBUG_FLAGS
	if( debugFlags.noLeftRight ) return;
#endif // DEBUG_FLAGS
	if( rightStream )
	{
		int xSize = _padSize;
		int ySize = (iters+1)*Degree+2;
		if(!iters)	ySize+=Degree;
		int offset = _padSize-2*WordPerDegree*RealPerWord;
		int size = (xSize+offset) * ySize * Channels + 2 * xSize * Channels;
		int yEnd = Degree+1;
		Pointer( SyncType ) right = AllocArray< SyncType >( size , 1 , "SocketedStreamingSolver::SyncSolverRight (right)" );
		if( read )
		{
#if TIME_IO
			double t = Time();
#endif // TIME_IO
#if DEBUG_FLAGS
			if( debugFlags.emptyPackets )
			{
				if ( !rightStream->read( ( Pointer( byte ) )right , sizeof( SyncType ) ) )	exit(0);
				FreeArray( right );
				return;
			}
#endif // DEBUG_FLAGS
			if ( !rightStream->read( ( Pointer( byte ) )right , sizeof( SyncType )*size ) )	exit(0);
#if TIME_IO
			hSync += Time() - t;
#endif // TIME_IO
			for( int y=-(ySize-yEnd-1) ; y<=yEnd ; y++ )
			{
				int yy = y + (ySize-yEnd-1);
#if SUPPORT_CYLINDRICAL
				if( (j+y>=0 && j+y<minor) || periodicType==SPHERICAL_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
				if( (j+y>=0 && j+y<minor) || spherical)
#endif // SUPPORT_CYLINDRICAL
					for(int c=0;c<Channels;c++)
					{
						Pointer( float ) tmp = GetXRow( j+y , c );
						if( tmp ) for( int x=0 ; x<xSize ; x++ ) tmp[ _size + x ] = float( right[ (x+offset) * ySize * Channels + yy * Channels + c ] );
					}
			}
			if( overlapped )
				for( int k=0 ; k<iters ; k++ )
				{
					int y = -k*Degree;
					int yy = y + (ySize-yEnd-1);
					if( j+y>=0 && j+y<minor )
						for( int c=0 ; c<Channels ; c++ )
						{
							Pointer( float ) tmp = GetXRow( j+y , c );
							if( tmp ) for( int x=WordPerDegree*RealPerWord*k ; x<offset ; x++ ) tmp[ _size - offset + x ] = float( right[ x * ySize * Channels + yy * Channels + c ] );
						}
				}
			for(int c=0;c<Channels;c++)
			{
				Pointer( float ) tmp = GetBRow( j , c );
				if( tmp ) for( int x=0 ; x<xSize ; x++ ) tmp[ _size+x ] = float( right[ (xSize+offset) * ySize * Channels + c * xSize + x ] );
			}
#if SUPPORT_CYLINDRICAL
			if( (j+1<minor) || periodicType==SPHERICAL_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if( (j+1<minor) || spherical)
#endif // SUPPORT_CYLINDRICAL
				for(int c=0;c<Channels;c++)
				{
					Pointer( float ) tmp = GetBRow( j+1 , c );
					if( tmp ) for( int x=0 ; x<xSize ; x++ ) tmp[ _size+x ] = float( right[ (xSize+offset) * ySize * Channels + Channels * xSize + c * xSize + x ] );
				}
		}
		else
		{
#if DEBUG_FLAGS
			if( debugFlags.emptyPackets )
			{
				if ( !rightStream->write( ( Pointer( byte ) )right , sizeof( SyncType ) ) ) exit(0);
				FreeArray( right );
				return;
			}
#endif // DEBUG_FLAGS
			memset(right,0, sizeof( SyncType )*size);
			for(int y=-(ySize-yEnd-1);y<=yEnd;y++)
			{
				int yy = y + (ySize-yEnd-1);
#if SUPPORT_CYLINDRICAL
				if( (j+y>=0 && j+y<minor) || periodicType==SPHERICAL_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
				if( (j+y>=0 && j+y<minor) || spherical)
#endif // SUPPORT_CYLINDRICAL
					for(int c=0;c<Channels;c++)
					{
						Pointer( float ) tmp = GetXRow( j+y , c );
						if( tmp ) for( int x=0 ; x<xSize+offset ; x++ ) right[ x * ySize * Channels + yy * Channels + c ] = SyncType( tmp[ _size-xSize+x ] );
					}
			}
			for(int c=0;c<Channels;c++)
			{
				Pointer( float ) tmp = GetBRow( j , c );
				if( tmp ) for( int x=0 ; x<xSize ; x++ ) right[ (xSize+offset) * ySize * Channels + c * xSize + x ] = SyncType( tmp[ _size - xSize + x ] );
			}
#if SUPPORT_CYLINDRICAL
			if( (j+1<minor) || periodicType==SPHERICAL_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if( (j+1<minor) || spherical)
#endif // SUPPORT_CYLINDRICAL
				for(int c=0;c<Channels;c++)
				{
					Pointer( float ) tmp = GetBRow( j+1 , c );
					if( tmp ) for( int x=0 ; x<xSize ; x++ ) right[ (xSize+offset) * ySize * Channels + Channels * xSize + c * xSize + x ] = SyncType( tmp[ _size - xSize + x ] );
				}
			if ( !rightStream->write( ( Pointer( byte ) )right , sizeof( SyncType )*size ) ) exit(0);
		}
		FreeArray( right );
	}
}
template< int Channels , class SyncType >
void SocketedStreamingSolver< Channels , SyncType >::SolveInterior( int j , int c , int sSolve , int eSolve )
{
	int jj=Degree*3;
	{
		Pointer( float ) localXPtrs[2*Degree+1];
		ConstPointer( float ) localBPtr = GetBRow( j , c );
		for( int y=0 ; y<2*Degree+1 ; y++ ) localXPtrs[y] = GetXRow( j-Degree+y , c );
		{
			Pointer( __m128 ) xPtrs[] = { ( Pointer( __m128 ) )localXPtrs[0] , ( Pointer( __m128 ) )localXPtrs[1] , ( Pointer( __m128 ) )localXPtrs[3] , ( Pointer( __m128 ) )localXPtrs[4] };

			// Perform the accumulation of the vertically symmetric rows and place them into __m128 buffers
			for( int i=-WordPerDegree+((sSolve-_start)>>2) ; i<((eSolve-_start)>>2)+WordPerDegree ; i++ )
			{
				localXAccum[0][i] = _mm_add_ps( xPtrs[0][i] , xPtrs[3][i] );
				localXAccum[1][i] = _mm_add_ps( xPtrs[1][i] , xPtrs[2][i] );
			}
		}

		float dotSum;
		__m128 dSum;
		__declspec (align(16)) float scratch[4];
		ConstPointer( __m128 ) xPtrs[] = { ( Pointer( __m128 ) )localXAccum[0] , ( Pointer( __m128 ) )localXAccum[1] , ( Pointer( __m128 ) )localXPtrs[2] };
		int s,e;

		// Offset 0
		{
			__m128 mValues[]=
			{
				zeroLapTemplates[jj+1].matrixValues[0][0],
				zeroLapTemplates[jj+1].matrixValues[0][1],
				zeroLapTemplates[jj+1].matrixValues[0][2]
			};
#if SUPPORT_CYLINDRICAL
			if( sSolve==0 && periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if(sSolve==0 && !spherical)
#endif // SUPPORT_CYLINDRICAL
			{
				dotSum=0;
				localXPtrs[2][0] = InteriorGaussSeidelUpdate0(zeroLapTemplates[jj].matrixValues[0],xPtrs,localBPtr,dotSum,0)*zeroLapTemplates[jj].diagonal[0];
				s=4;
			}
			else
			{
				SetInteriorDotSum(mValues,xPtrs,-1+((sSolve-_start)>>2),dSum);
				_mm_store_ps(scratch,dSum);
				dotSum=scratch[2]+scratch[3];
				s=sSolve-_start;
			}
#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC )	e=eSolve-_start-4;
#else // SUPPORT_CYLINDRICAL
			if(eSolve==major && !spherical)	e=eSolve-_start-4;
#endif // SUPPORT_CYLINDRICAL
			else							e=eSolve-_start;
			for(int i=s;i<e;i+=4)				localXPtrs[2][i]=InteriorGaussSeidelUpdate0(mValues,xPtrs,localBPtr,dotSum,i);
			int i=eSolve-_start-4;
#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC )		localXPtrs[2][i]=InteriorGaussSeidelUpdate0(zeroLapTemplates[jj+2].matrixValues[0],xPtrs,localBPtr,dotSum,i)*zeroLapTemplates[jj+2].diagonal[0];
#else // !SUPPORT_CYLINDRICAL
			if(eSolve==major && !spherical)		localXPtrs[2][i]=InteriorGaussSeidelUpdate0(zeroLapTemplates[jj+2].matrixValues[0],xPtrs,localBPtr,dotSum,i)*zeroLapTemplates[jj+2].diagonal[0];
#endif // SUPPORT_CYLINDRICAL
		}
		// Offset 1
		{
			__m128 mValues[]=
			{
				zeroLapTemplates[jj+1].matrixValues[1][0],
				zeroLapTemplates[jj+1].matrixValues[1][1],
				zeroLapTemplates[jj+1].matrixValues[1][2]
			};
#if SUPPORT_CYLINDRICAL
			if( sSolve==0 && periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if(sSolve==0 && !spherical)
#endif // SUPPORT_CYLINDRICAL
			{
				dotSum=0;
				localXPtrs[2][1]=InteriorGaussSeidelUpdate1(zeroLapTemplates[jj].matrixValues[1],xPtrs,localBPtr,dotSum,1)*zeroLapTemplates[jj].diagonal[1];
				s=5;
			}
			else
			{
				SetInteriorDotSum(mValues,xPtrs,-1+((sSolve-_start)>>2),dSum);
				_mm_store_ps(scratch,dSum);
				dotSum=scratch[3];
				s=sSolve-_start+1;
			}
#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC)	e=eSolve-_start-4;
#else // !SUPPORT_CYLINDRICAL
			if(eSolve==major && !spherical)	e=eSolve-_start-4;
#endif // SUPPORT_CYLINDRICAL
			else							e=eSolve-_start;
			for(int i=s;i<e;i+=4)				localXPtrs[2][i]=InteriorGaussSeidelUpdate1(mValues,xPtrs,localBPtr,dotSum,i);
			int i=eSolve-_start-3;
#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC )		localXPtrs[2][i]=InteriorGaussSeidelUpdate1(zeroLapTemplates[jj+2].matrixValues[1],xPtrs,localBPtr,dotSum,i)*zeroLapTemplates[jj+2].diagonal[1];
#else // !SUPPORT_CYLINDRICAL
			if(eSolve==major && !spherical)		localXPtrs[2][i]=InteriorGaussSeidelUpdate1(zeroLapTemplates[jj+2].matrixValues[1],xPtrs,localBPtr,dotSum,i)*zeroLapTemplates[jj+2].diagonal[1];
#endif // SUPPORT_CYLINDRICAL
		}
		// Offset 2
		{
			__m128 mValues[]=
			{
				zeroLapTemplates[jj+1].matrixValues[2][0],
				zeroLapTemplates[jj+1].matrixValues[2][1],
				zeroLapTemplates[jj+1].matrixValues[2][2]
			};
#if SUPPORT_CYLINDRICAL
			if( sSolve==0 && periodicType==NO_PERIODIC )	SetInteriorDotSum(zeroLapTemplates[jj].matrixValues[2],xPtrs,0,dSum);
#else // !SUPPORT_CYLINDRICAL
			if(sSolve==0 && !spherical)	SetInteriorDotSum(zeroLapTemplates[jj].matrixValues[2],xPtrs,0,dSum);
#endif // SUPPORT_CYLINDRICAL
			else						SetInteriorDotSum(mValues,xPtrs,(sSolve-_start)>>2,dSum);
			s=(sSolve-_start)+2;
			_mm_store_ps(scratch,dSum);
			dotSum=scratch[0]+scratch[1]+scratch[2]+scratch[3];
#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC )	e=eSolve-_start-8;
#else // !SUPPORT_CYLINDRICAL
			if(eSolve==major && !spherical)	e=eSolve-_start-8;
#endif // SUPPORT_CYLINDRICAL
			else							e=eSolve-_start;
			for(int i=s;i<e;i+=4)	localXPtrs[2][i]=InteriorGaussSeidelUpdate2(mValues,xPtrs,localBPtr,dotSum,i);
#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if(eSolve==major && !spherical)
#endif // SUPPORT_CYLINDRICAL
			{
				int i=eSolve-_start-6;
				localXPtrs[2][i]=InteriorGaussSeidelUpdate2(zeroLapTemplates[jj+2].matrixValues[2],xPtrs,localBPtr,dotSum,i);
				i+=4;
				localXPtrs[2][i]=(localBPtr[i]-dotSum)*zeroLapTemplates[jj+2].diagonal[2];
			}
		}
		// Offset 3
		{
			__m128 mValues[]=
			{
				zeroLapTemplates[jj+1].matrixValues[3][0],
				zeroLapTemplates[jj+1].matrixValues[3][1],
				zeroLapTemplates[jj+1].matrixValues[3][2]
			};
#if SUPPORT_CYLINDRICAL
			if( sSolve==0 && periodicType==NO_PERIODIC )	SetInteriorDotSum(zeroLapTemplates[jj].matrixValues[3],xPtrs,0,dSum);
#else // !SUPPORT_CYLINDRICAL
			if(sSolve==0 && !spherical)	SetInteriorDotSum(zeroLapTemplates[jj].matrixValues[3],xPtrs,0,dSum);
#endif // SUPPORT_CYLINDRICAL
			else						SetInteriorDotSum(mValues,xPtrs,(sSolve-_start)>>2,dSum);
			s=(sSolve-_start)+3;
			_mm_store_ps(scratch,dSum);
			dotSum=scratch[1]+scratch[2]+scratch[3];
#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC )	e=eSolve-_start-8;
#else // !SUPPORT_CYLINDRICAL
			if(eSolve==major && !spherical)	e=eSolve-_start-8;
#endif // SUPPORT_CYLINDRICAL
			else							e=eSolve-_start;
			for(int i=s;i<e;i+=4)	localXPtrs[2][i]=InteriorGaussSeidelUpdate3(mValues,xPtrs,localBPtr,dotSum,i);
#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if(eSolve==major && !spherical)
#endif // SUPPORT_CYLINDRICAL
			{
				int i=eSolve-_start-5;
				localXPtrs[2][i]=InteriorGaussSeidelUpdate3(zeroLapTemplates[jj+2].matrixValues[3],xPtrs,localBPtr,dotSum,i);
				i+=4;
				localXPtrs[2][i]=(localBPtr[i]-dotSum)*zeroLapTemplates[jj+2].diagonal[3];
			}
		}
#if MISHA_ADJUST_FLOATS
		for( int i=sSolve-_start ; i<eSolve-_start ; i++ ) if( localXPtrs[2][i]<MISHA_MIN_FLOAT && localXPtrs[2][i]>-MISHA_MIN_FLOAT ) localXPtrs[2][i] = 0;
#endif // MISHA_ADJUST_FLOATS
	}
}

template< int Channels , class SyncType >
void SocketedStreamingSolver< Channels , SyncType >::Solve( int j , int c , int sSolve , int eSolve )
{
#if DEBUG_FLAGS
	if(debugFlags.noSolver)	return;
#endif // DEBUG_FLAGS
	int jj;
	if		(j<Degree)			jj=j;
	else if	(j>minor-1-Degree)	jj=2*Degree+(j-(minor-1));
	else						jj=Degree;
#if SUPPORT_CYLINDRICAL
	if( periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
	if(!spherical)
#endif // SUPPORT_CYLINDRICAL
	{
		if(sSolve<0)		sSolve=0;
		if(eSolve>major)	eSolve=major;
	}
#if SUPPORT_CYLINDRICAL
	if( jj==Degree || periodicType==SPHERICAL_PERIODIC )	return SolveInterior(j,c,sSolve,eSolve);
#else // !SUPPORT_CYLINDRICAL
	if(jj==Degree || spherical)	return SolveInterior(j,c,sSolve,eSolve);
#endif // SUPPORT_CYLINDRICAL
	jj*=3;

	{
		Pointer( float ) localXPtrs[2*Degree+1];
		ConstPointer( float ) localBPtr = GetBRow(j,c);
		for(int y=0;y<2*Degree+1;y++)
			if(j-Degree+y>=0 && j-Degree+y<minor)	localXPtrs[y] = GetXRow(j-Degree+y,c);
			else									localXPtrs[y] = NullPointer< float >( );
		float dotSum;
		__m128 dSum;
		__declspec (align(16)) float scratch[4];
		ConstPointer( __m128 ) xPtrs[] =
		{
			( ConstPointer( __m128 ) )localXPtrs[0],
			( ConstPointer( __m128 ) )localXPtrs[1],
			( ConstPointer( __m128 ) )localXPtrs[2],
			( ConstPointer( __m128 ) )localXPtrs[3],
			( ConstPointer( __m128 ) )localXPtrs[4]
		};
		int s,e;

		// Offset 0
		{
			__m128 mValues[]=
			{
				zeroLapTemplates[jj+1].matrixValues[0][0],
				zeroLapTemplates[jj+1].matrixValues[0][1],
				zeroLapTemplates[jj+1].matrixValues[0][2],
				zeroLapTemplates[jj+1].matrixValues[0][3],
				zeroLapTemplates[jj+1].matrixValues[0][4]
			};
			float d=zeroLapTemplates[jj+1].diagonal[0];
#if SUPPORT_CYLINDRICAL
			if( sSolve==0 && periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if(sSolve==0 && !spherical)
#endif // SUPPORT_CYLINDRICAL
			{
				dotSum=0;
				localXPtrs[2][0] = GaussSeidelUpdate0(zeroLapTemplates[jj].matrixValues[0],xPtrs,localBPtr,dotSum,0)*zeroLapTemplates[jj].diagonal[0];
				s=4;
			}
			else
			{
				SetDotSum(mValues,xPtrs,-1+((sSolve-_start)>>2),dSum);
				_mm_store_ps(scratch,dSum);
				dotSum=scratch[2]+scratch[3];
				s=sSolve-_start;
			}
#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC )	e=eSolve-_start-4;
#else // !SUPPORT_CYLINDRICAL
			if(eSolve==major && !spherical)	e=eSolve-_start-4;
#endif // SUPPORT_CYLINDRICAL
			else							e=eSolve-_start;
			for(int i=s;i<e;i+=4)			localXPtrs[2][i]=GaussSeidelUpdate0(mValues,xPtrs,localBPtr,dotSum,i)*d;
			int i=eSolve-_start-4;
#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC )		localXPtrs[2][i]=GaussSeidelUpdate0(zeroLapTemplates[jj+2].matrixValues[0],xPtrs,localBPtr,dotSum,i)*zeroLapTemplates[jj+2].diagonal[0];
#else // !SUPPORT_CYLINDRICAL
			if(eSolve==major && !spherical)		localXPtrs[2][i]=GaussSeidelUpdate0(zeroLapTemplates[jj+2].matrixValues[0],xPtrs,localBPtr,dotSum,i)*zeroLapTemplates[jj+2].diagonal[0];
#endif // SUPPORT_CYLINDRICAL
		}
		// Offset 1
		{
			__m128 mValues[]=
			{
				zeroLapTemplates[jj+1].matrixValues[1][0],
				zeroLapTemplates[jj+1].matrixValues[1][1],
				zeroLapTemplates[jj+1].matrixValues[1][2],
				zeroLapTemplates[jj+1].matrixValues[1][3],
				zeroLapTemplates[jj+1].matrixValues[1][4]
			};
			float d=zeroLapTemplates[jj+1].diagonal[1];
#if SUPPORT_CYLINDRICAL
			if( sSolve==0 && periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if(sSolve==0 && !spherical)
#endif // SUPPORT_CYLINDRICAL
			{
				dotSum=0;
				localXPtrs[2][1]=GaussSeidelUpdate1(zeroLapTemplates[jj].matrixValues[1],xPtrs,localBPtr,dotSum,1)*zeroLapTemplates[jj].diagonal[1];
				s=5;
			}
			else
			{
				SetDotSum(mValues,xPtrs,-1+((sSolve-_start)>>2),dSum);
				_mm_store_ps(scratch,dSum);
				dotSum=scratch[3];
				s=sSolve-_start+1;
			}
#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC )	e=eSolve-_start-4;
#else // !SUPPORT_CYLINDRICAL
			if(eSolve==major && !spherical)	e=eSolve-_start-4;
#endif // SUPPORT_CYLINDRICAL
			else							e=eSolve-_start;
			for(int i=s;i<e;i+=4)			localXPtrs[2][i]=GaussSeidelUpdate1(mValues,xPtrs,localBPtr,dotSum,i)*d;
			int i=eSolve-_start-3;
#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC )		localXPtrs[2][i]=GaussSeidelUpdate1(zeroLapTemplates[jj+2].matrixValues[1],xPtrs,localBPtr,dotSum,i)*zeroLapTemplates[jj+2].diagonal[1];
#else // !SUPPORT_CYLINDRICAL
			if(eSolve==major && !spherical)		localXPtrs[2][i]=GaussSeidelUpdate1(zeroLapTemplates[jj+2].matrixValues[1],xPtrs,localBPtr,dotSum,i)*zeroLapTemplates[jj+2].diagonal[1];
#endif // SUPPORT_CYLINDRICAL
		}
		// Offset 2
		{
			__m128 mValues[]=
			{
				zeroLapTemplates[jj+1].matrixValues[2][0],
				zeroLapTemplates[jj+1].matrixValues[2][1],
				zeroLapTemplates[jj+1].matrixValues[2][2],
				zeroLapTemplates[jj+1].matrixValues[2][3],
				zeroLapTemplates[jj+1].matrixValues[2][4]
			};
			float d=zeroLapTemplates[jj+1].diagonal[2];
#if SUPPORT_CYLINDRICAL
			if( sSolve==0 && periodicType==NO_PERIODIC )	SetDotSum(zeroLapTemplates[jj].matrixValues[2],xPtrs,0,dSum);
#else // !SUPPORT_CYLINDRICAL
			if(sSolve==0 && !spherical)	SetDotSum(zeroLapTemplates[jj].matrixValues[2],xPtrs,0,dSum);
#endif // SUPPORT_CYLINDRICAL
			else						SetDotSum(mValues,xPtrs,(sSolve-_start)>>2,dSum);
			s=(sSolve-_start)+2;
			_mm_store_ps(scratch,dSum);
			dotSum=scratch[0]+scratch[1]+scratch[2]+scratch[3];
#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC )	e=eSolve-_start-8;
#else // !SUPPORT_CYLINDRICAL
			if(eSolve==major && !spherical)	e=eSolve-_start-8;
#endif // SUPPORT_CYLINDRICAL
			else							e=eSolve-_start;
			for(int i=s;i<e;i+=4)	localXPtrs[2][i]=GaussSeidelUpdate2(mValues,xPtrs,localBPtr,dotSum,i)*d;
#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if(eSolve==major && !spherical)
#endif // SUPPORT_CYLINDRICAL
			{
				int i=eSolve-_start-6;
				localXPtrs[2][i]=GaussSeidelUpdate2(zeroLapTemplates[jj+2].matrixValues[2],xPtrs,localBPtr,dotSum,i)*d;
				i+=4;
				localXPtrs[2][i]=(localBPtr[i]-dotSum)*zeroLapTemplates[jj+2].diagonal[2];
			}
		}
		// Offset 3
		{
			__m128 mValues[]=
			{
				zeroLapTemplates[jj+1].matrixValues[3][0],
				zeroLapTemplates[jj+1].matrixValues[3][1],
				zeroLapTemplates[jj+1].matrixValues[3][2],
				zeroLapTemplates[jj+1].matrixValues[3][3],
				zeroLapTemplates[jj+1].matrixValues[3][4]
			};
			float d=zeroLapTemplates[jj+1].diagonal[3];
#if SUPPORT_CYLINDRICAL
			if( sSolve==0 && periodicType==NO_PERIODIC )	SetDotSum(zeroLapTemplates[jj].matrixValues[3],xPtrs,0,dSum);
#else // !SUPPORT_CYLINDRICAL
			if(sSolve==0 && !spherical)	SetDotSum(zeroLapTemplates[jj].matrixValues[3],xPtrs,0,dSum);
#endif // SUPPORT_CYLINDRICAL
			else						SetDotSum(mValues,xPtrs,(sSolve-_start)>>2,dSum);
			s=(sSolve-_start)+3;
			_mm_store_ps(scratch,dSum);
			dotSum=scratch[1]+scratch[2]+scratch[3];
#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC )	e=eSolve-_start-8;
#else // !SUPPORT_CYLINDRICAL
			if(eSolve==major && !spherical)	e=eSolve-_start-8;
#endif // SUPPORT_CYLINDRICAL
			else							e=eSolve-_start;
			for(int i=s;i<e;i+=4)	localXPtrs[2][i]=GaussSeidelUpdate3(mValues,xPtrs,localBPtr,dotSum,i)*d;
#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if(eSolve==major && !spherical)
#endif // SUPPORT_CYLINDRICAL
			{
				int i=eSolve-_start-5;
				localXPtrs[2][i]=GaussSeidelUpdate3(zeroLapTemplates[jj+2].matrixValues[3],xPtrs,localBPtr,dotSum,i)*d;
				i+=4;
				localXPtrs[2][i]=(localBPtr[i]-dotSum)*zeroLapTemplates[jj+2].diagonal[3];
			}
		}
#if MISHA_ADJUST_FLOATS
		for( int i=sSolve-_start ; i<eSolve-_start ; i++ ) if( localXPtrs[2][i]<MISHA_MIN_FLOAT && localXPtrs[2][i]>-MISHA_MIN_FLOAT ) localXPtrs[2][i] = 0;
#endif // MISHA_ADJUST_FLOATS
	}
}

template< int Channels , class SyncType >
void SocketedStreamingSolver< Channels , SyncType >::SolveInteriorReverse( int j , int c , int sSolve , int eSolve )
{
	int jj = Degree*3;

	{
		Pointer( float ) localXPtrs[ 2*Degree + 1 ];
		ConstPointer( float ) localBPtr = GetBRow( j , c );
		for( int y=0 ; y<2*Degree+1 ; y++ ) localXPtrs[y] = GetXRow( j-Degree+y , c );
		{
			ConstPointer( __m128 ) xPtrs[] = { ( Pointer( __m128 ) )localXPtrs[0] , ( Pointer( __m128 ) )localXPtrs[1] , ( Pointer( __m128 ) )localXPtrs[3] , ( Pointer( __m128 ) )localXPtrs[4] };

			// Perform the accumulation of the vertically symmetric rows and place them into __m128 buffers
			for(int i=-WordPerDegree+((sSolve-_start)>>2);i<((eSolve-_start)>>2)+WordPerDegree;i++)
			{
				localXAccum[0][i]=_mm_add_ps( xPtrs[0][i] , xPtrs[3][i] );
				localXAccum[1][i]=_mm_add_ps( xPtrs[1][i] , xPtrs[2][i] );
			}
		}
		float dotSum;
		__m128 dSum;
		__declspec (align(16)) float scratch[4];
		ConstPointer( __m128 ) xPtrs[]={ localXAccum[0] , localXAccum[1] , ( Pointer( __m128 ) )localXPtrs[2] };
		int s , e ;

		// Offset 3
		{
			__m128 mValues[]=
			{
				zeroLapTemplates[jj+1].matrixValues[3][0],
				zeroLapTemplates[jj+1].matrixValues[3][1],
				zeroLapTemplates[jj+1].matrixValues[3][2]
			};
#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if( eSolve==major && !spherical )
#endif // SUPPORT_CYLINDRICAL
			{
				dotSum = 0;
				int i = eSolve - _start - 1;
				localXPtrs[2][i] = ReverseInteriorGaussSeidelUpdate3( zeroLapTemplates[jj+2].matrixValues[3] , xPtrs , localBPtr , dotSum , i ) * zeroLapTemplates[jj+2].diagonal[3];
				e = eSolve - _start - 5;
			}
			else
			{
				SetInteriorDotSum( mValues , xPtrs , 1 + ( (eSolve - 1 - _start)>>2 ) , dSum );
				_mm_store_ps( scratch , dSum );
				dotSum = scratch[0] + scratch[1];
				e = eSolve - _start - 1;
			}
#if SUPPORT_CYLINDRICAL
			if( sSolve==0 && periodicType==NO_PERIODIC ) s = 7;
#else // !SUPPORT_CYLINDRICAL
			if( sSolve==0 && !spherical ) s = 7;
#endif // SUPPORT_CYLINDRICAL
			else						  s = sSolve - _start;
			for( int i=e ; i>=s ; i-=4 ) localXPtrs[2][i] = ReverseInteriorGaussSeidelUpdate3( mValues , xPtrs , localBPtr , dotSum , i );
			int i = 3;
#if SUPPORT_CYLINDRICAL
			if( sSolve==0 && periodicType==NO_PERIODIC ) localXPtrs[2][i] = ReverseInteriorGaussSeidelUpdate3( zeroLapTemplates[jj].matrixValues[3] , xPtrs , localBPtr , dotSum , i ) * zeroLapTemplates[jj].diagonal[3];
#else // !SUPPORT_CYLINDRICAL
			if( sSolve==0 && !spherical ) localXPtrs[2][i] = ReverseInteriorGaussSeidelUpdate3( zeroLapTemplates[jj].matrixValues[3] , xPtrs , localBPtr , dotSum , i ) * zeroLapTemplates[jj].diagonal[3];
#endif // SUPPORT_CYLINDRICAL
		}
		// Offset 2
		{
			__m128 mValues[]=
			{
				zeroLapTemplates[jj+1].matrixValues[2][0],
				zeroLapTemplates[jj+1].matrixValues[2][1],
				zeroLapTemplates[jj+1].matrixValues[2][2]
			};
#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if( eSolve==major && !spherical )
#endif // SUPPORT_CYLINDRICAL
			{
				dotSum = 0;
				int i = eSolve - _start - 2;
				localXPtrs[2][i] = ReverseInteriorGaussSeidelUpdate2( zeroLapTemplates[jj+2].matrixValues[2] , xPtrs , localBPtr , dotSum , i ) * zeroLapTemplates[jj+2].diagonal[2];
				e = eSolve - _start - 6;
			}
			else
			{
				SetInteriorDotSum( mValues , xPtrs , 1 + ( (eSolve - 1 - _start)>>2 ) , dSum );
				_mm_store_ps( scratch , dSum );
				dotSum = scratch[0];
				e = eSolve - _start - 2;
			}
#if SUPPORT_CYLINDRICAL
			if( sSolve==0 && periodicType==NO_PERIODIC ) s = 6;
#else // !SUPPORT_CYLINDRICAL
			if( sSolve==0 && !spherical ) s = 6;
#endif // SUPPORT_CYLINDRICAL
			else						  s = sSolve - _start;
			for( int i=e ; i>=s ; i-=4 ) localXPtrs[2][i] = ReverseInteriorGaussSeidelUpdate2( mValues , xPtrs , localBPtr , dotSum , i );
			int i = 2;
#if SUPPORT_CYLINDRICAL
			if( sSolve==0 && periodicType==NO_PERIODIC ) localXPtrs[2][i] = ReverseInteriorGaussSeidelUpdate2( zeroLapTemplates[jj].matrixValues[2] , xPtrs , localBPtr , dotSum , i ) * zeroLapTemplates[jj].diagonal[2];
#else // !SUPPORT_CYLINDRICAL
			if( sSolve==0 && !spherical ) localXPtrs[2][i] = ReverseInteriorGaussSeidelUpdate2( zeroLapTemplates[jj].matrixValues[2] , xPtrs , localBPtr , dotSum , i ) * zeroLapTemplates[jj].diagonal[2];
#endif // SUPPORT_CYLINDRICAL
		}

		// Offset 1
		{
			__m128 mValues[]=
			{
				zeroLapTemplates[jj+1].matrixValues[1][0],
				zeroLapTemplates[jj+1].matrixValues[1][1],
				zeroLapTemplates[jj+1].matrixValues[1][2]
			};
#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if( eSolve==major && !spherical )
#endif // SUPPORT_CYLINDRICAL
			{
				SetInteriorDotSum( zeroLapTemplates[jj+2].matrixValues[1] , xPtrs , (eSolve - 1 - _start )>>2 , dSum );
				_mm_store_ps( scratch , dSum );
				dotSum = scratch[0] + scratch[1] + scratch[2] + scratch[3];
				int i = eSolve - _start - 3;
				localXPtrs[2][i] = ReverseInteriorGaussSeidelUpdate1( mValues , xPtrs , localBPtr , dotSum , i ) * zeroLapTemplates[jj+2].diagonal[1];
				e = eSolve - _start - 7;
			}
			else
			{
				SetInteriorDotSum( mValues , xPtrs , (eSolve - 1 - _start )>>2 , dSum );
				_mm_store_ps( scratch , dSum );
				dotSum = scratch[0] + scratch[1] + scratch[2] + scratch[3];
				e = eSolve - _start - 3;
			}
#if SUPPORT_CYLINDRICAL
			if( sSolve==0 && periodicType==NO_PERIODIC ) s = 9;
#else // !SUPPORT_CYLINDRICAL
			if( sSolve==0 && !spherical ) s = 9;
#endif // SUPPORT_CYLINDRICAL
			else						  s = sSolve - _start;
			for( int i=e ; i>=s ; i-=4 ) localXPtrs[2][i] = ReverseInteriorGaussSeidelUpdate1( mValues , xPtrs , localBPtr , dotSum , i );
#if SUPPORT_CYLINDRICAL
			if( sSolve==0 && periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if( sSolve==0 && !spherical )
#endif // SUPPORT_CYLINDRICAL
			{
				int i = 5;
				localXPtrs[2][i] = ReverseInteriorGaussSeidelUpdate1( zeroLapTemplates[jj].matrixValues[1] , xPtrs , localBPtr , dotSum , i );
				i = 1;
				localXPtrs[2][i] = (localBPtr[i]-dotSum) * zeroLapTemplates[jj].diagonal[1];			}
		}

		// Offset 0
		{
			__m128 mValues[]=
			{
				zeroLapTemplates[jj+1].matrixValues[0][0],
				zeroLapTemplates[jj+1].matrixValues[0][1],
				zeroLapTemplates[jj+1].matrixValues[0][2]
			};
			float d = zeroLapTemplates[jj+1].diagonal[0];

#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if( eSolve==major && !spherical )
#endif // SUPPORT_CYLINDRICAL
			{
				SetInteriorDotSum( zeroLapTemplates[jj+2].matrixValues[0] , xPtrs , (eSolve - 1 - _start)>>2 , dSum );
				_mm_store_ps( scratch , dSum );
				dotSum = scratch[0] + scratch[1] + scratch[2];
				int i = eSolve - _start - 4;
				localXPtrs[2][i] = ReverseInteriorGaussSeidelUpdate0( mValues , xPtrs , localBPtr , dotSum , i ) * zeroLapTemplates[jj+2].diagonal[0];
				e = eSolve - _start - 8;
			}
			else
			{
				SetInteriorDotSum( mValues , xPtrs , ( eSolve - 1 - _start )>>2 , dSum );
				_mm_store_ps( scratch , dSum );
				dotSum = scratch[0] + scratch[1] + scratch[2];
				e = eSolve - _start - 4;
			}
#if SUPPORT_CYLINDRICAL
			if( sSolve==0 && periodicType==NO_PERIODIC ) s = 8;
#else // !SUPPORT_CYLINDRICAL
			if( sSolve==0 && !spherical ) s = 8;
#endif // SUPPORT_CYLINDRICAL
			else						  s = sSolve - _start;
			for( int i=e ; i>=s ; i-=4 ) localXPtrs[2][i] = ReverseInteriorGaussSeidelUpdate0( mValues , xPtrs , localBPtr , dotSum , i );
#if SUPPORT_CYLINDRICAL
			if( sSolve==0 && periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if( sSolve==0 && !spherical )
#endif // SUPPORT_CYLINDRICAL
			{
				int i = 4;
				localXPtrs[2][i] = ReverseInteriorGaussSeidelUpdate0( zeroLapTemplates[jj].matrixValues[0] , xPtrs , localBPtr , dotSum , i );
				i = 0;
				localXPtrs[2][i] = (localBPtr[i]-dotSum) * zeroLapTemplates[jj].diagonal[0];
			}
		}
#if MISHA_ADJUST_FLOATS
		for( int i=sSolve-_start ; i<eSolve-_start ; i++ ) if( localXPtrs[2][i]<MISHA_MIN_FLOAT && localXPtrs[2][i]>-MISHA_MIN_FLOAT ) localXPtrs[2][i] = 0;
#endif // MISHA_ADJUST_FLOATS
	}
}
template< int Channels , class SyncType >
void SocketedStreamingSolver< Channels , SyncType >::SolveReverse( int j , int c , int sSolve , int eSolve )
{
#if DEBUG_FLAGS
	if(debugFlags.noSolver)	return;
#endif // DEBUG_FLAGS
	int jj;
	if		(j<Degree)			jj=j;
	else if	(j>minor-1-Degree)	jj=2*Degree+(j-(minor-1));
	else						jj=Degree;

#if SUPPORT_CYLINDRICAL
	if( periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
	if(!spherical)
#endif // SUPPORT_CYLINDRICAL
	{
		if( sSolve<0 )		sSolve = 0;
		if( eSolve>major )	eSolve = major;
	}
#if SUPPORT_CYLINDRICAL
	if( jj==Degree || periodicType==SPHERICAL_PERIODIC )	return SolveInteriorReverse( j , c , sSolve , eSolve );
#else // !SUPPORT_CYLINDRICAL
	if( jj==Degree || spherical )	return SolveInteriorReverse( j , c , sSolve , eSolve );
#endif // SUPPORT_CYLINDRICAL

	//	if( spherical ) jj=Degree;
	jj*=3;

	{
		Pointer( float ) localXPtrs[ 2*Degree + 1 ];
		ConstPointer( float ) localBPtr = GetBRow( j , c );
		for( int y=0 ; y<2*Degree+1 ; y++ )
#if SUPPORT_CYLINDRICAL
			if( (j-Degree+y>=0 && j-Degree+y<minor) || periodicType==SPHERICAL_PERIODIC )	localXPtrs[y] = GetXRow( j-Degree+y , c );
#else // !SUPPORT_CYLINDRICAL
			if( (j-Degree+y>=0 && j-Degree+y<minor) || spherical )	localXPtrs[y] = GetXRow( j-Degree+y , c );
#endif // SUPPORT_CYLINDRICAL
			else													localXPtrs[y] = NullPointer< float >( );

		float dotSum;
		__m128 dSum;
		__declspec (align(16)) float scratch[4];
		ConstPointer( __m128 ) xPtrs[] = 
		{
			( Pointer( __m128 ) )localXPtrs[0],
			( Pointer( __m128 ) )localXPtrs[1],
			( Pointer( __m128 ) )localXPtrs[2],
			( Pointer( __m128 ) )localXPtrs[3],
			( Pointer( __m128 ) )localXPtrs[4]
		};
		int s , e ;

		// Offset 3
		{
			__m128 mValues[]=
			{
				zeroLapTemplates[jj+1].matrixValues[3][0],
				zeroLapTemplates[jj+1].matrixValues[3][1],
				zeroLapTemplates[jj+1].matrixValues[3][2],
				zeroLapTemplates[jj+1].matrixValues[3][3],
				zeroLapTemplates[jj+1].matrixValues[3][4]
			};
			float d = zeroLapTemplates[jj+1].diagonal[3];

#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if( eSolve==major && !spherical )
#endif // SUPPORT_CYLINDRICAL
			{
				dotSum = 0;
				int i = eSolve - _start - 1;
				localXPtrs[2][i] = ReverseGaussSeidelUpdate3( zeroLapTemplates[jj+2].matrixValues[3] , xPtrs , localBPtr , dotSum , i ) * zeroLapTemplates[jj+2].diagonal[3];
				e = eSolve - _start - 5;
			}
			else
			{
				SetDotSum( mValues , xPtrs , 1 + ( (eSolve - 1 - _start)>>2 ) , dSum );
				_mm_store_ps( scratch , dSum );
				dotSum = scratch[0] + scratch[1];
				e = eSolve - _start - 1;
			}
#if SUPPORT_CYLINDRICAL
			if( sSolve==0 && periodicType==NO_PERIODIC ) s = 7;
#else // !SUPPORT_CYLINDRICAL
			if( sSolve==0 && !spherical ) s = 7;
#endif // SUPPORT_CYLINDRICAL
			else						  s = sSolve - _start;
			for( int i=e ; i>=s ; i-=4 ) localXPtrs[2][i] = ReverseGaussSeidelUpdate3( mValues , xPtrs , localBPtr , dotSum , i ) * d;
			int i = 3;
#if SUPPORT_CYLINDRICAL
			if( sSolve==0 && periodicType==NO_PERIODIC ) localXPtrs[2][i] = ReverseGaussSeidelUpdate3( zeroLapTemplates[jj].matrixValues[3] , xPtrs , localBPtr , dotSum , i ) * zeroLapTemplates[jj].diagonal[3];
#else // !SUPPORT_CYLINDRICAL
			if( sSolve==0 && !spherical ) localXPtrs[2][i] = ReverseGaussSeidelUpdate3( zeroLapTemplates[jj].matrixValues[3] , xPtrs , localBPtr , dotSum , i ) * zeroLapTemplates[jj].diagonal[3];
#endif // SUPPORT_CYLINDRICAL
		}
		// Offset 2
		{
			__m128 mValues[]=
			{
				zeroLapTemplates[jj+1].matrixValues[2][0],
				zeroLapTemplates[jj+1].matrixValues[2][1],
				zeroLapTemplates[jj+1].matrixValues[2][2],
				zeroLapTemplates[jj+1].matrixValues[2][3],
				zeroLapTemplates[jj+1].matrixValues[2][4]
			};
			float d = zeroLapTemplates[jj+1].diagonal[2];

#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if( eSolve==major && !spherical )
#endif // SUPPORT_CYLINDRICAL
			{
				dotSum = 0;
				int i = eSolve - _start - 2;
				localXPtrs[2][i] = ReverseGaussSeidelUpdate2( zeroLapTemplates[jj+2].matrixValues[2] , xPtrs , localBPtr , dotSum , i ) * zeroLapTemplates[jj+2].diagonal[2];
				e = eSolve - _start - 6;
			}
			else
			{
				SetDotSum( mValues , xPtrs , 1 + ( (eSolve - 1 - _start)>>2 ) , dSum );
				_mm_store_ps( scratch , dSum );
				dotSum = scratch[0];
				e = eSolve - _start - 2;
			}
#if SUPPORT_CYLINDRICAL
			if( sSolve==0 && periodicType==NO_PERIODIC ) s = 6;
#else // !SUPPORT_CYLINDRICAL
			if( sSolve==0 && !spherical ) s = 6;
#endif // SUPPORT_CYLINDRICAL
			else						  s = sSolve - _start;
			for( int i=e ; i>=s ; i-=4 ) localXPtrs[2][i] = ReverseGaussSeidelUpdate2( mValues , xPtrs , localBPtr , dotSum , i ) * d;
			int i = 2;
#if SUPPORT_CYLINDRICAL
			if( sSolve==0 && periodicType==NO_PERIODIC ) localXPtrs[2][i] = ReverseGaussSeidelUpdate2( zeroLapTemplates[jj].matrixValues[2] , xPtrs , localBPtr , dotSum , i ) * zeroLapTemplates[jj].diagonal[2];
#else // !SUPPORT_CYLINDRICAL
			if( sSolve==0 && !spherical ) localXPtrs[2][i] = ReverseGaussSeidelUpdate2( zeroLapTemplates[jj].matrixValues[2] , xPtrs , localBPtr , dotSum , i ) * zeroLapTemplates[jj].diagonal[2];
#endif // SUPPORT_CYLINDRICAL
		}

		// Offset 1
		{
			__m128 mValues[]=
			{
				zeroLapTemplates[jj+1].matrixValues[1][0],
				zeroLapTemplates[jj+1].matrixValues[1][1],
				zeroLapTemplates[jj+1].matrixValues[1][2],
				zeroLapTemplates[jj+1].matrixValues[1][3],
				zeroLapTemplates[jj+1].matrixValues[1][4]
			};
			float d = zeroLapTemplates[jj+1].diagonal[1];

#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if( eSolve==major && !spherical )
#endif // SUPPORT_CYLINDRICAL
			{
				SetDotSum( zeroLapTemplates[jj+2].matrixValues[1] , xPtrs , (eSolve - 1 - _start )>>2 , dSum );
				_mm_store_ps( scratch , dSum );
				dotSum = scratch[0] + scratch[1] + scratch[2] + scratch[3];
				int i = eSolve - _start - 3;
				localXPtrs[2][i] = ReverseGaussSeidelUpdate1( mValues , xPtrs , localBPtr , dotSum , i ) * zeroLapTemplates[jj+2].diagonal[1];
				e = eSolve - _start - 7;
			}
			else
			{
				SetDotSum( mValues , xPtrs , (eSolve - 1 - _start )>>2 , dSum );
				_mm_store_ps( scratch , dSum );
				dotSum = scratch[0] + scratch[1] + scratch[2] + scratch[3];
				e = eSolve - _start - 3;
			}
#if SUPPORT_CYLINDRICAL
			if( sSolve==0 && periodicType==NO_PERIODIC ) s = 9;
#else // !SUPPORT_CYLINDRICAL
			if( sSolve==0 && !spherical ) s = 9;
#endif // SUPPORT_CYLINDRICAL
			else						  s = sSolve - _start;
			for( int i=e ; i>=s ; i-=4 ) localXPtrs[2][i] = ReverseGaussSeidelUpdate1( mValues , xPtrs , localBPtr , dotSum , i ) * d;
#if SUPPORT_CYLINDRICAL
			if( sSolve==0 && periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if( sSolve==0 && !spherical )
#endif // SUPPORT_CYLINDRICAL
			{
				int i = 5;
				localXPtrs[2][i] = ReverseGaussSeidelUpdate1( zeroLapTemplates[jj].matrixValues[1] , xPtrs , localBPtr , dotSum , i ) * d;
				i = 1;
				localXPtrs[2][i] = (localBPtr[i]-dotSum) * zeroLapTemplates[jj].diagonal[1];
			}
		}

		// Offset 0
		{
			__m128 mValues[]=
			{
				zeroLapTemplates[jj+1].matrixValues[0][0],
				zeroLapTemplates[jj+1].matrixValues[0][1],
				zeroLapTemplates[jj+1].matrixValues[0][2],
				zeroLapTemplates[jj+1].matrixValues[0][3],
				zeroLapTemplates[jj+1].matrixValues[0][4]
			};
			float d = zeroLapTemplates[jj+1].diagonal[0];

#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if( eSolve==major && !spherical )
#endif // SUPPORT_CYLINDRICAL
			{
				SetDotSum( zeroLapTemplates[jj+2].matrixValues[0] , xPtrs , (eSolve - 1 - _start)>>2 , dSum );
				_mm_store_ps( scratch , dSum );
				dotSum = scratch[0] + scratch[1] + scratch[2];
				int i = eSolve - _start - 4;
				localXPtrs[2][i] = ReverseGaussSeidelUpdate0( mValues , xPtrs , localBPtr , dotSum , i ) * zeroLapTemplates[jj+2].diagonal[0];
				e = eSolve - _start - 8;
			}
			else
			{
				SetDotSum( mValues , xPtrs , ( eSolve - 1 - _start )>>2 , dSum );
				_mm_store_ps( scratch , dSum );
				dotSum = scratch[0] + scratch[1] + scratch[2];
				e = eSolve - _start - 4;
			}
#if SUPPORT_CYLINDRICAL
			if( sSolve==0 && periodicType==NO_PERIODIC ) s = 8;
#else // !SUPPORT_CYLINDRICAL
			if( sSolve==0 && !spherical ) s = 8;
#endif // SUPPORT_CYLINDRICAL
			else						  s = sSolve - _start;
			for( int i=e ; i>=s ; i-=4 ) localXPtrs[2][i] = ReverseGaussSeidelUpdate0( mValues , xPtrs , localBPtr , dotSum , i ) * d;
#if SUPPORT_CYLINDRICAL
			if( sSolve==0 && periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if( sSolve==0 && !spherical )
#endif // SUPPORT_CYLINDRICAL
			{
				int i = 4;
				localXPtrs[2][i] = ReverseGaussSeidelUpdate0( zeroLapTemplates[jj].matrixValues[0] , xPtrs , localBPtr , dotSum , i ) * d;
				i = 0;
				localXPtrs[2][i] = (localBPtr[i]-dotSum) * zeroLapTemplates[jj].diagonal[0];
			}
		}
#if MISHA_ADJUST_FLOATS
		for( int i=sSolve-_start ; i<eSolve-_start ; i++ ) if( localXPtrs[2][i]<MISHA_MIN_FLOAT && localXPtrs[2][i]>-MISHA_MIN_FLOAT ) localXPtrs[2][i] = 0;
#endif // MISHA_ADJUST_FLOATS
	}
}

template< int Channels , class SyncType >
void SocketedStreamingSolver< Channels , SyncType >::SetInteriorResidual(int j,int c,int sSolve,int eSolve)
{
	int jj=Degree*3;
	{
		ConstPointer( float ) localXPtrs[2*Degree+1];
		Pointer( float )       localRPtr = GetRRow(j,c);
		ConstPointer( float ) localBPtr = GetBRow(j,c);
		for(int y=0;y<2*Degree+1;y++)	localXPtrs[y]=GetXRow(j-Degree+y,c);
		{
			ConstPointer( __m128 ) xPtrs[] =
			{
				( ConstPointer( __m128 ) )localXPtrs[0],
				( ConstPointer( __m128 ) )localXPtrs[1],
				( ConstPointer( __m128 ) )localXPtrs[3],
				( ConstPointer( __m128 ) )localXPtrs[4]
			};
			// Perform the accumulation of the vertically symmetric rows and place them into __m128 buffers
			for(int i=-WordPerDegree+((sSolve-_start)>>2);i<((eSolve-_start)>>2)+WordPerDegree;i++)
			{
				localXAccum[0][i]=_mm_add_ps(xPtrs[0][i],xPtrs[3][i]);
				localXAccum[1][i]=_mm_add_ps(xPtrs[1][i],xPtrs[2][i]);
			}
		}


		float dotSum;
		__m128 dSum;
		__declspec (align(16)) float scratch[4];
		ConstPointer( __m128 ) xPtrs[] = { localXAccum[0] , localXAccum[1] , ( ConstPointer( __m128 ) )localXPtrs[2]};
		int s,e;

		// Offset 0
		{
			__m128 mValues[]=
			{
				lapTemplates[jj+1].matrixValues[0][0],
				lapTemplates[jj+1].matrixValues[0][1],
				lapTemplates[jj+1].matrixValues[0][2]
			};
#if SUPPORT_CYLINDRICAL
			if( sSolve==0 && periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if(sSolve==0 && !spherical)
#endif // SUPPORT_CYLINDRICAL
			{
				dotSum=0;
				localRPtr[0]=localBPtr[0]-GetInteriorLaplacianValue0(lapTemplates[jj].matrixValues[0],xPtrs,dotSum,0);
				s=4;
			}
			else
			{
				SetInteriorDotSum(mValues,xPtrs,-1+((sSolve-_start)>>2),dSum);
				s=sSolve-_start;
				_mm_store_ps(scratch,dSum);
				dotSum=scratch[2]+scratch[3];
			}
#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC )	e=eSolve-_start-4;
#else // !SUPPORT_CYLINDRICAL
			if(eSolve==major && !spherical)	e=eSolve-_start-4;
#endif // SUPPORT_CYLINDRICAL
			else							e=eSolve-_start;

			for(int i=s;i<e;i+=4)			localRPtr[i]=localBPtr[i]-GetInteriorLaplacianValue0(mValues,xPtrs,dotSum,i);
			int i=eSolve-_start-4;
#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC )	localRPtr[i]=localBPtr[i]-GetInteriorLaplacianValue0(lapTemplates[jj+2].matrixValues[0],xPtrs,dotSum,i);
#else // !SUPPORT_CYLINDRICAL
			if(eSolve==major && !spherical)	localRPtr[i]=localBPtr[i]-GetInteriorLaplacianValue0(lapTemplates[jj+2].matrixValues[0],xPtrs,dotSum,i);
#endif // SUPPORT_CYLINDRICAL
		}
		// Offset 1
		{
			__m128 mValues[]=
			{
				lapTemplates[jj+1].matrixValues[1][0],
				lapTemplates[jj+1].matrixValues[1][1],
				lapTemplates[jj+1].matrixValues[1][2]
			};
#if SUPPORT_CYLINDRICAL
			if( sSolve==0 && periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if(sSolve==0 && !spherical)
#endif // SUPPORT_CYLINDRICAL
			{
				dotSum=0;
				localRPtr[1]=localBPtr[1]-GetInteriorLaplacianValue1(lapTemplates[jj].matrixValues[1],xPtrs,dotSum,1);
				s=5;
			}
			else
			{
				SetInteriorDotSum(mValues,xPtrs,-1+((sSolve-_start)>>2),dSum);
				s=sSolve-_start+1;
				_mm_store_ps(scratch,dSum);
				dotSum=scratch[3];
			}
#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC )	e=eSolve-_start-4;
#else // !SUPPORT_CYLINDRICAL
			if(eSolve==major && !spherical)	e=eSolve-_start-4;
#endif // SUPPORT_CYLINDRICAL
			else							e=eSolve-_start;
			for(int i=s;i<e;i+=4)			localRPtr[i]=localBPtr[i]-GetInteriorLaplacianValue1(mValues,xPtrs,dotSum,i);
			int i=eSolve-_start-3;
#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC )	localRPtr[i]=localBPtr[i]-GetInteriorLaplacianValue1(lapTemplates[jj+2].matrixValues[1],xPtrs,dotSum,i);
#else // !SUPPORT_CYLINDRICAL
			if(eSolve==major && !spherical)	localRPtr[i]=localBPtr[i]-GetInteriorLaplacianValue1(lapTemplates[jj+2].matrixValues[1],xPtrs,dotSum,i);
#endif // SUPPORT_CYLINDRICAL
		}
		// Offset 2
		{
			__m128 mValues[]=
			{
				lapTemplates[jj+1].matrixValues[2][0],
				lapTemplates[jj+1].matrixValues[2][1],
				lapTemplates[jj+1].matrixValues[2][2]
			};
#if SUPPORT_CYLINDRICAL
			if( sSolve==0 && periodicType==NO_PERIODIC )	SetInteriorDotSum(lapTemplates[jj].matrixValues[2],xPtrs,0,dSum);
#else // !SUPPORT_CYLINDRICAL
			if(sSolve==0 && !spherical)	SetInteriorDotSum(lapTemplates[jj].matrixValues[2],xPtrs,0,dSum);
#endif // SUPPORT_CYLINDRICAL
			else						SetInteriorDotSum(mValues,xPtrs,(sSolve-_start)>>2,dSum);
			s=(sSolve-_start)+2;
			_mm_store_ps(scratch,dSum);
			dotSum=scratch[0]+scratch[1]+scratch[2]+scratch[3];
#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC )	e=eSolve-_start-8;
#else // !SUPPORT_CYLINDRICAL
			if(eSolve==major && !spherical)	e=eSolve-_start-8;
#endif // SUPPORT_CYLINDRICAL
			else							e=eSolve-_start;
			for(int i=s;i<e;i+=4)	localRPtr[i]=localBPtr[i]-GetInteriorLaplacianValue2(mValues,xPtrs,dotSum,i);
#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if(eSolve==major && !spherical)
#endif // SUPPORT_CYLINDRICAL
			{
				int i=eSolve-_start-6;
				localRPtr[i]=localBPtr[i]-GetInteriorLaplacianValue2(lapTemplates[jj+2].matrixValues[2],xPtrs,dotSum,i);
				i+=4;
				localRPtr[i]=localBPtr[i]-dotSum;
			}
		}
		// Offset 3
		{
			__m128 mValues[]=
			{
				lapTemplates[jj+1].matrixValues[3][0],
				lapTemplates[jj+1].matrixValues[3][1],
				lapTemplates[jj+1].matrixValues[3][2]
			};
#if SUPPORT_CYLINDRICAL
			if( sSolve==0 && periodicType==NO_PERIODIC )	SetInteriorDotSum(lapTemplates[jj].matrixValues[3],xPtrs,0,dSum);
#else // !SUPPORT_CYLINDRICAL
			if(sSolve==0 && !spherical)	SetInteriorDotSum(lapTemplates[jj].matrixValues[3],xPtrs,0,dSum);
#endif // SUPPORT_CYLINDRICAL
			else						SetInteriorDotSum(mValues,xPtrs,(sSolve-_start)>>2,dSum);
			s=(sSolve-_start)+3;
			_mm_store_ps(scratch,dSum);
			dotSum=scratch[1]+scratch[2]+scratch[3];
#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC )	e=eSolve-_start-8;
#else // !SUPPORT_CYLINDRICAL
			if(eSolve==major && !spherical)	e=eSolve-_start-8;
#endif // SUPPORT_CYLINDRICAL
			else							e=eSolve-_start;
			for(int i=s;i<e;i+=4)	localRPtr[i]=localBPtr[i]-GetInteriorLaplacianValue3(mValues,xPtrs,dotSum,i);
#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if(eSolve==major && !spherical)
#endif // SUPPORT_CYLINDRICAL
			{
				int i=eSolve-_start-5;
				localRPtr[i]=localBPtr[i]-GetInteriorLaplacianValue3(lapTemplates[jj+2].matrixValues[3],xPtrs,dotSum,i);
				i+=4;
				localRPtr[i]=localBPtr[i]-dotSum;
			}
		}
	}
}
template< int Channels , class SyncType >
void SocketedStreamingSolver< Channels , SyncType >::SetResidual(int j,int c,int sSolve,int eSolve)
{
#if DEBUG_FLAGS
	if(debugFlags.noResidual)	return;
#endif // DEBUG_FLAGS
	int jj;
	if		(j<Degree)			jj=j;
	else if	(j>minor-1-Degree)	jj=2*Degree+(j-(minor-1));
	else						jj=Degree;

#if SUPPORT_CYLINDRICAL
	if( periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
	if(!spherical)
#endif // SUPPORT_CYLINDRICAL
	{
		if(sSolve<0)		sSolve=0;
		if(eSolve>major)	eSolve=major;
	}
#if SUPPORT_CYLINDRICAL
	if( jj==Degree || periodicType==SPHERICAL_PERIODIC )	return SetInteriorResidual(j,c,sSolve,eSolve);
#else // !SUPPORT_CYLINDRICAL
	if(jj==Degree || spherical)	return SetInteriorResidual(j,c,sSolve,eSolve);
#endif // SUPPORT_CYLINDRICAL
	jj*=3;

	{
		ConstPointer( float ) localXPtrs[2*Degree+1];
		Pointer( float )       localRPtr = GetRRow(j,c);
		ConstPointer( float ) localBPtr = GetBRow(j,c);
		for(int y=0;y<2*Degree+1;y++)
			if(j-Degree+y>=0 && j-Degree+y<minor)	localXPtrs[y]=GetXRow(j-Degree+y,c);
			else									localXPtrs[y] = NullPointer< float >( );

		float dotSum;
		__m128 dSum;
		__declspec (align(16)) float scratch[4];
		ConstPointer( __m128 ) xPtrs[] =
		{
			( ConstPointer( __m128 ) )localXPtrs[0],
			( ConstPointer( __m128 ) )localXPtrs[1],
			( ConstPointer( __m128 ) )localXPtrs[2],
			( ConstPointer( __m128 ) )localXPtrs[3],
			( ConstPointer( __m128 ) )localXPtrs[4]
		};
		int s,e;
		// Offset 0
		{
			__m128 mValues[]=
			{
				lapTemplates[jj+1].matrixValues[0][0],
				lapTemplates[jj+1].matrixValues[0][1],
				lapTemplates[jj+1].matrixValues[0][2],
				lapTemplates[jj+1].matrixValues[0][3],
				lapTemplates[jj+1].matrixValues[0][4]
			};
#if SUPPORT_CYLINDRICAL
			if( sSolve==0 && periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if(sSolve==0 && !spherical)
#endif // SUPPORT_CYLINDRICAL
			{
				dotSum=0;
				localRPtr[0]=localBPtr[0]-GetLaplacianValue0(lapTemplates[jj].matrixValues[0],xPtrs,dotSum,0);
				s=4;
			}
			else
			{
				SetDotSum(mValues,xPtrs,-1+((sSolve-_start)>>2),dSum);
				s=sSolve-_start;
				_mm_store_ps(scratch,dSum);
				dotSum=scratch[2]+scratch[3];
			}
#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC)	e=eSolve-_start-4;
#else // !SUPPORT_CYLINDRICAL
			if(eSolve==major && !spherical)	e=eSolve-_start-4;
#endif // SUPPORT_CYLINDRICAL3
			else							e=eSolve-_start;
			for(int i=s;i<e;i+=4)			localRPtr[i]=localBPtr[i]-GetLaplacianValue0(mValues,xPtrs,dotSum,i);
			int i=eSolve-_start-4;
#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC )	localRPtr[i]=localBPtr[i]-GetLaplacianValue0(lapTemplates[jj+2].matrixValues[0],xPtrs,dotSum,i);
#else // !SUPPORT_CYLINDRICAL
			if(eSolve==major && !spherical)	localRPtr[i]=localBPtr[i]-GetLaplacianValue0(lapTemplates[jj+2].matrixValues[0],xPtrs,dotSum,i);
#endif // SUPPORT_CYLINDRICAL
		}
		// Offset 1
		{
			__m128 mValues[]=
			{
				lapTemplates[jj+1].matrixValues[1][0],
				lapTemplates[jj+1].matrixValues[1][1],
				lapTemplates[jj+1].matrixValues[1][2],
				lapTemplates[jj+1].matrixValues[1][3],
				lapTemplates[jj+1].matrixValues[1][4]
			};
#if SUPPORT_CYLINDRICAL
			if( sSolve==0 && periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if(sSolve==0 && !spherical)
#endif // SUPPORT_CYLINDRICAL
			{
				dotSum=0;
				localRPtr[1]=localBPtr[1]-GetLaplacianValue1(lapTemplates[jj].matrixValues[1],xPtrs,dotSum,1);
				s=5;
			}
			else
			{
				SetDotSum(mValues,xPtrs,-1+((sSolve-_start)>>2),dSum);
				s=sSolve-_start+1;
				_mm_store_ps(scratch,dSum);
				dotSum=scratch[3];
			}
#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC )	e=eSolve-_start-4;
#else // !SUPPORT_CYLINDRICAL
			if(eSolve==major && !spherical)	e=eSolve-_start-4;
#endif // SUPPORT_CYLINDRICAL
			else							e=eSolve-_start;
			for(int i=s;i<e;i+=4)			localRPtr[i]=localBPtr[i]-GetLaplacianValue1(mValues,xPtrs,dotSum,i);
			int i=eSolve-_start-3;
#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC )	localRPtr[i]=localBPtr[i]-GetLaplacianValue1(lapTemplates[jj+2].matrixValues[1],xPtrs,dotSum,i);
#else // !SUPPORT_CYLINDRICAL
			if(eSolve==major && !spherical)	localRPtr[i]=localBPtr[i]-GetLaplacianValue1(lapTemplates[jj+2].matrixValues[1],xPtrs,dotSum,i);
#endif // SUPPORT_CYLINDRICAL
		}
		// Offset 2
		{
			__m128 mValues[]=
			{
				lapTemplates[jj+1].matrixValues[2][0],
				lapTemplates[jj+1].matrixValues[2][1],
				lapTemplates[jj+1].matrixValues[2][2],
				lapTemplates[jj+1].matrixValues[2][3],
				lapTemplates[jj+1].matrixValues[2][4]
			};
#if SUPPORT_CYLINDRICAL
			if( sSolve==0 && periodicType==NO_PERIODIC )	SetDotSum(lapTemplates[jj].matrixValues[2],xPtrs,0,dSum);
#else // !SUPPORT_CYLINDRICAL
			if(sSolve==0 && !spherical)	SetDotSum(lapTemplates[jj].matrixValues[2],xPtrs,0,dSum);
#endif // SUPPORT_CYLINDRICAL
			else						SetDotSum(mValues,xPtrs,(sSolve-_start)>>2,dSum);
			s=(sSolve-_start)+2;
			_mm_store_ps(scratch,dSum);
			dotSum=scratch[0]+scratch[1]+scratch[2]+scratch[3];
#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC )	e=eSolve-_start-8;
#else // !SUPPORT_CYLINDRICAL
			if(eSolve==major && !spherical)	e=eSolve-_start-8;
#endif // SUPPORT_CYLINDRICAL
			else							e=eSolve-_start;
			for(int i=s;i<e;i+=4)	localRPtr[i]=localBPtr[i]-GetLaplacianValue2(mValues,xPtrs,dotSum,i);
#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if(eSolve==major && !spherical)
#endif // SUPPORT_CYLINDRICAL
			{
				int i=eSolve-_start-6;
				localRPtr[i]=localBPtr[i]-GetLaplacianValue2(lapTemplates[jj+2].matrixValues[2],xPtrs,dotSum,i);
				i+=4;
				localRPtr[i]=localBPtr[i]-dotSum;
			}
		}
		// Offset 3
		{
			__m128 mValues[]=
			{
				lapTemplates[jj+1].matrixValues[3][0],
				lapTemplates[jj+1].matrixValues[3][1],
				lapTemplates[jj+1].matrixValues[3][2],
				lapTemplates[jj+1].matrixValues[3][3],
				lapTemplates[jj+1].matrixValues[3][4]
			};
#if SUPPORT_CYLINDRICAL
			if( sSolve==0 && periodicType==NO_PERIODIC )	SetDotSum(lapTemplates[jj].matrixValues[3],xPtrs,0,dSum);
#else // !SUPPORT_CYLINDRICAL
			if(sSolve==0 && !spherical)	SetDotSum(lapTemplates[jj].matrixValues[3],xPtrs,0,dSum);
#endif // SUPPORT_CYLINDRICAL
			else						SetDotSum(mValues,xPtrs,(sSolve-_start)>>2,dSum);
			s=(sSolve-_start)+3;
			_mm_store_ps(scratch,dSum);
			dotSum=scratch[1]+scratch[2]+scratch[3];
#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC )	e=eSolve-_start-8;
#else // !SUPPORT_CYLINDRICAL
			if(eSolve==major && !spherical)	e=eSolve-_start-8;
#endif // SUPPORT_CYLINDRICAL
			else							e=eSolve-_start;
			for(int i=s;i<e;i+=4)	localRPtr[i]=localBPtr[i]-GetLaplacianValue3(mValues,xPtrs,dotSum,i);
#if SUPPORT_CYLINDRICAL
			if( eSolve==major && periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if(eSolve==major && !spherical)
#endif // SUPPORT_CYLINDRICAL
			{
				int i=eSolve-_start-5;
				localRPtr[i]=localBPtr[i]-GetLaplacianValue3(lapTemplates[jj+2].matrixValues[3],xPtrs,dotSum,i);
				i+=4;
				localRPtr[i]=localBPtr[i]-dotSum;
			}
		}
	}
}

template< int Channels , class SyncType >
void SocketedStreamingSolver< Channels , SyncType >::Solve( void )
{
	int idx = index+iters*Degree-1;
	int start,stop;

	stop = index+iters*Degree-1;
	start = stop-(iters-1)*Degree;	// = index + Degree - 1;
#if SUPPORT_CYLINDRICAL
	if( periodicType==SPHERICAL_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
	if( spherical )
#endif // SUPPORT_CYLINDRICAL
	{
		if( stop==0 )
		{
			// This is necessary, because when the same resolution row is solved twice, we want to reuse the solution
			for( int d=0 ; d<Degree ; d++ ) SyncSolverHead( d , false ) , SyncSolverHead( d, true );
		}
		if( stop==minor-Degree )
		{
			for( int d=minor-Degree ; d<minor ; d++ ) SyncSolverTail( d , false ) , SyncSolverTail( d , true );
			SyncSolverRight( stop , false ) ,	SyncSolverLeft ( stop , true );
		}
	}
	if( stop==0 ) SyncSolverRight( -1 , false ) , SyncSolverLeft( -1 , true );
	{
		int offset;
		if( stop&1 ) // Solve right to left
		{
			offset = _padSize-2*WordPerDegree*RealPerWord;
			for( int i=stop ; i>=start ; i-=Degree )
			{
				if( i>=0 && i<minor ) for( int c=0 ; c<Channels ; c++ ) SolveReverse( i , c , _end-_padSize-offset , _end+offset );
				offset -= WordPerDegree*RealPerWord;
			}
			SyncSolverRight( stop , false );
			for(int c=0;c<Channels;c++)
			{
				int s128 = _end128-_start128-3*_padSize128;
				for(int l=0;l<laneNum;l++)
				{
					int b128 = _start128+2*_padSize128+((laneNum-l-1)*s128)/laneNum;
					int e128 = _start128+2*_padSize128+((laneNum-l  )*s128)/laneNum;
					int b = b128<<2;
					int e = e128<<2;

					offset=_padSize-2*WordPerDegree*RealPerWord;
					for(int i=stop;i>=start;i-=Degree)
					{
						if(i>=0 && i<minor)	SolveReverse( i , c , b-offset , e-offset );
						offset -= WordPerDegree*RealPerWord;
					}
				}
			}
			SyncSolverLeft( stop , true );
			offset=_padSize-2*WordPerDegree*RealPerWord;
			for(int i=stop;i>=start;i-=Degree)
			{
			if( leftStream )
			{
				if( i>=0 && i<minor ) for( int c=0 ; c<Channels ; c++ ) SolveReverse( i , c , _start+offset , _start+2*_padSize-offset );
			}
			else
			{
				if( i>=0 && i<minor ) for( int c=0 ; c<Channels ; c++ ) SolveReverse( i , c , _start        , _start+2*_padSize-offset );
			}
				offset -= WordPerDegree*RealPerWord;
			}
		}
		else // Solve left to right
		{
			// Solve left
			offset=_padSize-2*WordPerDegree*RealPerWord;
			for( int i=stop ; i>=start ; i-=Degree )
			{
				if(i>=0 && i<minor)	for(int c=0;c<Channels;c++)	Solve( i , c , _start-offset , _start+_padSize+offset );
				offset -= WordPerDegree*RealPerWord;
			}
			// Write out to right buffer
			SyncSolverLeft( stop , false );
			// Solve center
			for(int c=0;c<Channels;c++)
			{
				int s128 = _end128-_start128-3*_padSize128;
				for(int l=0;l<laneNum;l++)
				{
					int b128 = _start128+_padSize128+( l   *s128)/laneNum;
					int e128 = _start128+_padSize128+((l+1)*s128)/laneNum;
					int b = b128<<2;
					int e = e128<<2;

					offset=_padSize-2*WordPerDegree*RealPerWord;
					for(int i=stop;i>=start;i-=Degree)
					{
						if(i>=0 && i<minor)	Solve(i,c,b+offset,e+offset);
						offset -= WordPerDegree*RealPerWord;
					}
				}
			}
			// Read in right buffer
			SyncSolverRight( stop , true );
			// Solve Right
			offset=_padSize-2*WordPerDegree*RealPerWord;
			for( int i=stop ; i>=start ; i-=Degree )
			{
			if( rightStream )
			{
				if( i>=0 && i<minor ) for( int c=0 ; c<Channels ; c++ ) Solve( i , c , _end-2*_padSize+offset , _end-offset );
			}
			else
			{
				if( i>=0 && i<minor ) for( int c=0 ; c<Channels ; c++ ) Solve( i , c , _end-2*_padSize+offset , _end );
			}
				offset -= WordPerDegree*RealPerWord;
			}
		}
	}
#if SUPPORT_CYLINDRICAL
	if( periodicType==SPHERICAL_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
	if( spherical )
#endif // SUPPORT_CYLINDRICAL
	{
		if(start<Degree && stop>=0 )
		{
			int d = start;
			while( d<0 ) d += Degree;
			SyncSolverHead( d , false ) , SyncSolverHead( d , true );
			SyncSolverRight( start , false , false ) , SyncSolverLeft ( start , false , false );
			SyncSolverRight( start , true  , false ) , SyncSolverLeft ( start , true  , false );
		}
		if(stop>minor-Degree-1 && start<minor)
		{
			int d = start;
			while( d<minor-Degree ) d += Degree;
			SyncSolverTail( d , false ) , SyncSolverTail( d , true );
		}
	}
	idx=index+OffsetR();
	if( idx>=0 && idx<minor && setResidual )
	{
		for(int c=0;c<Channels;c++)
		{
			SetResidual(idx,c,_start-WordPerDegree*RealPerWord,_start+_padSize);
			int s128 = _end128-_start128-2*_padSize128;
			for(int l=0;l<laneNum;l++)
			{
				int b128 = _start128+_padSize128+((l  )*s128)/laneNum;
				int e128 = _start128+_padSize128+((l+1)*s128)/laneNum;
				int b = b128<<2;
				int e = e128<<2;
				SetResidual(idx,c,b,e);
			}
			SetResidual(idx,c,_end-_padSize,_end+WordPerDegree*RealPerWord);
		}
		if( idx<Degree )
		{
			SyncResidualHead( idx , false , false) , SyncResidualHead( idx , true ,  false);
			SyncResidualHead( idx , false , true ) , SyncResidualHead( idx , true ,  true );
		}
		if(idx<minor && idx>=minor-Degree)
		{
			SyncResidualTail( idx , false , false) , SyncResidualTail( idx , true , false );
			SyncResidualTail( idx , false , true ) , SyncResidualTail( idx , true ,  true );
		}

		for(int c=0;c<Channels;c++)
		{
			Pointer( float ) localXPtr = GetXRow(idx,c);
			Pointer( float ) localBPtr = GetBRow(idx,c);
			Pointer( float ) localRPtr = GetRRow(idx,c);
			for(int i=0;i<_size;i++)
			{
				xSquareNorm+=double(localXPtr[i])*localXPtr[i];
				rSquareNorm+=double(localRPtr[i])*localRPtr[i];
				bSquareNorm+=double(localBPtr[i])*localBPtr[i];
				solutionSum[c]+=double(localXPtr[i]);
			}
		}
	}
	if( showProgress )
	{
		if( !idx )
		{
			progressCount = 0;
			pastIndex = 0;
			startProgressTime = pastProgressTime = Time();
		}
		if( idx>=0 && idx<minor )
		{
			size_t current,peak;
			WorkingSetInfo(current,peak);
			double t = Time();
			int pCount = (idx*1000) / ( minor-1 );
			if( pCount>progressCount )
			{
				printf("[%.1f%%]  [%d/%d] [%4d/%4d MB] Rows/Second: %f\t(%f)         \r" , float(idx)/float(minor-1)*100 , idx , minor , current>>20 , peak>>20 , float(idx-pastIndex) / ( t-pastProgressTime ) , float(idx) / ( t-startProgressTime ) );
				progressCount = pCount;
				pastProgressTime = t;
				pastIndex = idx;
			}
		}
		else if(idx==minor)		printf( "\n" ) , fflush( stdout );
	}
}

//////////////////////////////////////
// SocketedMultiGridStreamingSolver //
//////////////////////////////////////
template< int Channels , class StorageType , class SyncType >
SocketedMultiGridStreamingSolver< Channels , StorageType , SyncType >::SocketedMultiGridStreamingSolver(void)
{
	X = B = NULL;
	inX = outX = inB = outB = outR = NULL;
	scratchR = NullPointer< float >( );
	parent = NULL;
	rChild = NULL;
	pChild = NULL;
	localRAccum = NullPointer< __m128 >( );
	prolongationStencil = NullPointer< ProlongationStencilSSE2 >( );

	int halfD=((Degree+1)>>1);
	int dual=(Degree&1)?0:1;

	// Restriction
	// To solve X[i] at depth d-1:
	//		Need to have B[i+iters-1] at depth d-1:
	//		Need to have R[2*(i+iters-1)-halfD,2*(i+iters-1)+halfD+dual]

	// If we update i at depth d:
	// @ depth d+1 i gets prolonged to: 	[2*i-halfD,2*i+halfD+dual]
	// @ depth d-1 i gets restricted to:	[(i-halfD-dual+1)/2,(i+halfD)/2]
	// Can start processing the parent once:
	//		i > dual+halfD-1
	// Can start processing the child once:
	//		i >= halfD/2


	// Restriction

	// Can start processing the parent once ((i+OffsetR())+1) gets restricted to a starting index bigger than zero:
	//		(i+OffsetR()+1-halfD-dual+1)/2 >= 1
	//		i+OffsetR()-halfD-dual >= 0
	//		i >= dual+halfD-OffsetR()
	startRestriction = dual + halfD - OffsetR();
	// In the coarser resolution, we only care that the B value has been set from the finer residual
	startRestriction -= 2*Degree;
	// Get the parity of the starting index so we know on what parity to pass control to the coarser resolution
	restrictionBit=startRestriction&1;


	// Prolongation
	// Can start processing the child once i gets prolonged to an index bigger than or equal to one:
	//		2*i-halfD >= 1
	//		i >= (halfD+1)/2
	startProlongation = (halfD+2)>>1;
	// If we prolong at idx=startProlongation, in the finer resolution it will map to a leading index of
	prolongationOffset = 2*startProlongation+halfD+dual;

	B=X=NULL;
}
template< int Channels , class StorageType , class SyncType >
void SocketedMultiGridStreamingSolver< Channels , StorageType , SyncType >::Initialize(int start,int end,int major,int minor,int iters,
																		DataStream* leftStream , SOCKET* syncSockets , DataStream* rightStream ,
#if SUPPORT_CYLINDRICAL
																		bool memoryMappedFile , int periodicType , MultiStreamIOServer* server
#else // !SUPPORT_CYLINDRICAL
																		bool memoryMappedFile , bool spherical , MultiStreamIOServer* server
#endif // SUPPORT_CYLINDRICAL
																		)
{
	DotProductStencil dotMajor,d2DotMajor,dotMinor,d2DotMinor;

	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(major,dotMajor,0,0);
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(major,d2DotMajor,1,1,false);
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(minor,dotMinor,0,0);
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(minor,d2DotMinor,1,1,false);

	Initialize( dotMajor , d2DotMajor , dotMinor , d2DotMinor , start , end , major , minor , iters ,
		leftStream , syncSockets , rightStream ,
#if SUPPORT_CYLINDRICAL
		memoryMappedFile , periodicType , server
#else // !SUPPORT_CYLIDNRICAL
		memoryMappedFile , spherical , server
#endif // SUPPORT_CYLINDRICAL
		);
}
template< int Channels , class StorageType , class SyncType >
void SocketedMultiGridStreamingSolver< Channels , StorageType , SyncType >::Initialize(DotProductStencil& dotMajor,DotProductStencil& d2DotMajor,
																		DotProductStencil& dotMinor,DotProductStencil& d2DotMinor,
																		int start,int end,int major,int minor,int iters,
																		DataStream* leftStream,SOCKET* syncSockets,DataStream* rightStream,
#if SUPPORT_CYLINDRICAL
																		bool memoryMappedFile , int periodicType , MultiStreamIOServer* server
#else // !SUPPORT_CYLIDNRICAL
																		bool memoryMappedFile,bool spherical , MultiStreamIOServer* server
#endif // SUPPORT_CYLINDRICAL
																		)
{
	this->dotMajor		= dotMajor;
	this->dotMinor		= dotMinor;
	this->d2DotMajor	= d2DotMajor;
	this->d2DotMinor	= d2DotMinor;

	MatrixStencil lStencil;
	for(int i=0;i<=2*Degree;i++)
		for(int j=0;j<=2*Degree;j++)
			for(int k=0;k<=2*Degree;k++)
				for(int l=0;l<=2*Degree;l++)
					lStencil.caseTable[j][i].values[l][k]=
					dotMajor.caseTable[i].values[k]*d2DotMinor.caseTable[j].values[l]+
					d2DotMajor.caseTable[i].values[k]*dotMinor.caseTable[j].values[l];
	if(syncSockets)
		Init( lStencil , start , end , major , minor , iters ,
		leftStream , syncSockets[0] , rightStream ,
#if SUPPORT_CYLINDRICAL
		periodicType , server
#else // !SUPPORT_CYLIDNRICAL
		spherical , server
#endif // SUPPORT_CYLINDRICAL
		);
	else
		Init( lStencil , start , end , major , minor , iters ,
		leftStream , INVALID_SOCKET , rightStream ,
#if SUPPORT_CYLINDRICAL
		periodicType , server
#else // !SUPPORT_CYLIDNRICAL
		spherical , server
#endif // SUPPORT_CYLINDRICAL
		);
	if(B)	delete B;
	if(X)	delete X;
	B=X=NULL;
	if(memoryMappedFile)
	{
#if DEBUG_FLAGS
		if(!debugFlags.noTempIO)
#endif // DEBUG_FLAGS
		{
			B = new MultiStreamIOClient( _size * Channels * sizeof(StorageType) , minor , STREAMING_GRID_BUFFER_MULTIPLIER , NULL , true );
			X = new MultiStreamIOClient( _size * Channels * sizeof(StorageType) , minor , STREAMING_GRID_BUFFER_MULTIPLIER , NULL , true);
		}
		inCore = false;
	}
	else
	{
		B = new MemoryBackedGrid( _size * Channels * sizeof( StorageType ) , minor , true );
		X = new MemoryBackedGrid( _size * Channels * sizeof( StorageType ) , minor , true );
		inCore = true;
	}
	FiniteElements1D<double,Type,Degree>::ProlongationStencil(major>>1,majorProlongationStencil,major);
	FiniteElements1D<double,Type,Degree>::ProlongationStencil(minor>>1,minorProlongationStencil,minor);
	int foo;
	FiniteElements1D<double,Type,Degree>::RestrictionStencil(major,majorRestrictionStencil,foo);
	FiniteElements1D<double,Type,Degree>::RestrictionStencil(minor,minorRestrictionStencil,foo);
	if(parent)
	{
		DotProductStencil newDotMajor,newD2DotMajor,newDotMinor,newD2DotMinor;
		CombineStencils<double,Type,Degree>(  dotMajor,majorProlongationStencil,major,  newDotMajor);
		CombineStencils<double,Type,Degree>(d2DotMajor,majorProlongationStencil,major,newD2DotMajor);
		CombineStencils<double,Type,Degree>(  dotMinor,minorProlongationStencil,minor,  newDotMinor);
		CombineStencils<double,Type,Degree>(d2DotMinor,minorProlongationStencil,minor,newD2DotMinor);

		if(syncSockets)
			parent->Initialize( newDotMajor , newD2DotMajor , newDotMinor , newD2DotMinor , start>>1 , end>>1 , major>>1 , minor>>1 , iters ,
			leftStream , syncSockets+1 , rightStream ,
#if SUPPORT_CYLINDRICAL
			memoryMappedFile , periodicType , _server
#else // !SUPPORT_CYLINDRICAL
			memoryMappedFile , spherical , _server
#endif // SUPPORT_CYLINDRICAL
			);
		else
			parent->Initialize( newDotMajor , newD2DotMajor , newDotMinor , newD2DotMinor , start>>1 , end>>1 , major>>1 , minor>>1 , iters ,
			leftStream , NULL , rightStream ,
#if SUPPORT_CYLINDRICAL
			memoryMappedFile , periodicType , _server
#else // !SUPPORT_CYLINDRICAL
			memoryMappedFile , spherical , _server
#endif // SUPPORT_CYLINDRICAL
			);
	}
}
template< int Channels , class StorageType , class SyncType >
void SocketedMultiGridStreamingSolver< Channels , StorageType , SyncType >::Initialize(double iWeight,int start,int end,int major,int minor,int iters,
																		DataStream* leftStream,SOCKET* syncSockets,DataStream* rightStream,
#if SUPPORT_CYLINDRICAL
																		bool memoryMappedFile , int periodicType , MultiStreamIOServer* server
#else // !SUPPORT_CYLINDRICAL
																		bool memoryMappedFile , bool spherical , MultiStreamIOServer* server
#endif // SUPPORT_CYLINDRICAL
																		)
{
	// Warning!!! A certaint amount of care needs to be taken here. It's safe to call this function when iWeight == 0 since
	// the Laplacian stencil is independent of the resolution. However, the dot-product stencil is not!!!
	// Care should be taken:
	// 1] When transitioning from out-of-core to in-core
	// 2] When transitioning from in-core to a conjugate-gradient solver

	DotProductStencil dotMajor,d2DotMajor,dotMinor,d2DotMinor;

	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(major,dotMajor,0,0);
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(major,d2DotMajor,1,1,false);
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(minor,dotMinor,0,0);
	FiniteElements1D<double,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(minor,d2DotMinor,1,1,false);

	Initialize( dotMajor , d2DotMajor , dotMinor , d2DotMinor , iWeight , start , end , major , minor , iters ,
		leftStream , syncSockets , rightStream ,
#if SUPPORT_CYLINDRICAL
		memoryMappedFile , periodicType , server
#else // !SUPPORT_CYLINDRICAL
		memoryMappedFile , spherical , server
#endif // SUPPORT_CYLINDRICAL
		);
}
template< int Channels , class StorageType , class SyncType >
void SocketedMultiGridStreamingSolver< Channels , StorageType , SyncType >::Initialize(DotProductStencil& dotMajor,DotProductStencil& d2DotMajor,
																		DotProductStencil& dotMinor,DotProductStencil& d2DotMinor,
																		double iWeight,
																		int start,int end,int major,int minor,int iters,
																		DataStream* leftStream,SOCKET* syncSockets,DataStream* rightStream,
#if SUPPORT_CYLINDRICAL
																		bool memoryMappedFile , int periodicType , MultiStreamIOServer* server
#else // !SUPPORT_CYLINDRICAL
																		bool memoryMappedFile,bool spherical , MultiStreamIOServer* server
#endif // SUPPORT_CYLINDRICAL
																		)
{
	this->dotMajor		= dotMajor;
	this->dotMinor		= dotMinor;
	this->d2DotMajor	= d2DotMajor;
	this->d2DotMinor	= d2DotMinor;

	MatrixStencil lStencil;
	for(int i=0;i<=2*Degree;i++)
		for(int j=0;j<=2*Degree;j++)
			for(int k=0;k<=2*Degree;k++)
				for(int l=0;l<=2*Degree;l++)
					lStencil.caseTable[j][i].values[l][k]=
					dotMajor.caseTable[i].values[k]*d2DotMinor.caseTable[j].values[l]+
					d2DotMajor.caseTable[i].values[k]*dotMinor.caseTable[j].values[l]+
					(dotMajor.caseTable[i].values[k]*dotMinor.caseTable[i].values[l])*iWeight;
	if(syncSockets)
		Init( lStencil , start , end , major , minor , iters ,
		leftStream , syncSockets[0] , rightStream ,
#if SUPPORT_CYLINDRICAL
		periodicType , server
#else // !SUPPORT_CYLINDRICAL
		spherical , server
#endif // SUPPORT_CYLINDRICAL
		);
	else
		Init( lStencil , start , end , major , minor , iters ,
		leftStream , INVALID_SOCKET , rightStream ,
#if SUPPORT_CYLINDRICAL
		periodicType , server
#else // !SUPPORT_CYLINDRICAL
		spherical , server
#endif // SUPPORT_CYLINDRICAL
		);

	if(B)	delete B;
	if(X)	delete X;
	B=X=NULL;

	if(memoryMappedFile)
	{
#if DEBUG_FLAGS
		if(!debugFlags.noTempIO)
#endif // DEBUG_FLAGS
		{
			B=new MultiStreamIOClient( _size*Channels*sizeof(StorageType) , minor , STREAMING_GRID_BUFFER_MULTIPLIER , NULL , true );
			X=new MultiStreamIOClient( _size*Channels*sizeof(StorageType) , minor , STREAMING_GRID_BUFFER_MULTIPLIER , NULL , true);
		}
		inCore = false;
	}
	else
	{
		B = new MemoryBackedGrid( _size * Channels * sizeof( StorageType ) , minor , true );
		X = new MemoryBackedGrid( _size * Channels * sizeof( StorageType ) , minor , true );
		inCore = true;
	}
	FiniteElements1D<double,Type,Degree>::ProlongationStencil(major>>1,majorProlongationStencil,major);
	FiniteElements1D<double,Type,Degree>::ProlongationStencil(minor>>1,minorProlongationStencil,minor);
	int foo;
	FiniteElements1D<double,Type,Degree>::RestrictionStencil(major,majorRestrictionStencil,foo);
	FiniteElements1D<double,Type,Degree>::RestrictionStencil(minor,minorRestrictionStencil,foo);
	if( parent )
	{
		DotProductStencil newDotMajor,newD2DotMajor,newDotMinor,newD2DotMinor;
		CombineStencils<double,Type,Degree>(  dotMajor,majorProlongationStencil,major,  newDotMajor);
		CombineStencils<double,Type,Degree>(d2DotMajor,majorProlongationStencil,major,newD2DotMajor);
		CombineStencils<double,Type,Degree>(  dotMinor,minorProlongationStencil,minor,  newDotMinor);
		CombineStencils<double,Type,Degree>(d2DotMinor,minorProlongationStencil,minor,newD2DotMinor);


		if(syncSockets)
			parent->Initialize( newDotMajor , newD2DotMajor , newDotMinor , newD2DotMinor , iWeight , start>>1 , end>>1 , major>>1 , minor>>1 , iters ,
			leftStream , syncSockets+1 , rightStream ,
#if SUPPORT_CYLINDRICAL
			memoryMappedFile , periodicType , _server
#else // !SUPPORT_CYLINDRICAL
			memoryMappedFile , spherical , _server
#endif // SUPPORT_CYLINDRICAL
			);
		else
			parent->Initialize( newDotMajor , newD2DotMajor , newDotMinor , newD2DotMinor , iWeight , start>>1 , end>>1 , major>>1 , minor>>1 , iters ,
			leftStream , NULL , rightStream ,
#if SUPPORT_CYLINDRICAL
			memoryMappedFile , periodicType , _server
#else // !SUPPORT_CYLINDRICAL
			memoryMappedFile , spherical , _server
#endif // SUPPORT_CYLINDRICAL
			);
	}
}
template< int Channels , class StorageType , class SyncType >
SocketedMultiGridStreamingSolver< Channels , StorageType , SyncType >::~SocketedMultiGridStreamingSolver(void)
{
	if( localRAccum )
	{
		localRAccum -= _padSize128;
		FreeArray( localRAccum );
	}
	FreeArray( prolongationStencil );
	if(X)	delete X;
	if(B)	delete B;
	X=B=NULL;
	FreeArray( scratchR );
}
template< int Channels , class StorageType , class SyncType >
void SocketedMultiGridStreamingSolver< Channels , StorageType , SyncType >::InitProlongation(void)
{
	if( parent ) parent->InitProlongation();
	if( inX )  inX->reset ( true  , 1 );
	if( inB )  inB->reset ( true  , 1 );
	if( outX ) outX->reset( false , 1 );
	if( outB ) outB->reset( false , 1 );
	if( outR ) outR->reset( false , 1 );
}
template< int Channels , class StorageType , class SyncType >
void SocketedMultiGridStreamingSolver< Channels , StorageType , SyncType >::InitRestriction( void )
{
	{
		localRAccum = AllocArray< __m128 >( _paddedSize128 , ALIGNMENT , "SocketedMultiGridStreamingSolver::localRAccum" );
		localRAccum += _padSize128;
		__declspec (align(16)) float scratch[4];
		prolongationStencil = AllocArray< ProlongationStencilSSE2 >( 3 , ALIGNMENT , "SocketedMultigridStreamingSolver::prolongationStencil" );
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
			prolongationStencil[j].matrixValues[0]=_mm_load_ps(scratch);

			if(j==2)	jj=2*Degree;
			else		jj=Degree;
			if(j!=0)	scratch[0]=majorProlongationStencil.caseTable[Degree].values[3];
			else		scratch[0]=0;
			scratch[1]=majorProlongationStencil.caseTable[jj].values[0];
			scratch[2]=majorProlongationStencil.caseTable[jj].values[1];
			scratch[3]=majorProlongationStencil.caseTable[jj].values[2];
			prolongationStencil[j].matrixValues[1]=_mm_load_ps(scratch);
		}
	}
	if( parent ) parent->InitRestriction();
	if( inX )  inX->reset ( true  , 1 );
	if( inB )  inB->reset ( true  , 1 );
	if( outX ) outX->reset( false , 1 );
	if( outB ) outB->reset( false , 1 );
	if( outR ) outR->reset( false , 1 );
}
template< int Channels , class StorageType , class SyncType >
void SocketedMultiGridStreamingSolver< Channels , StorageType , SyncType >::SetProlongation(void)
{
	// Zig-zagging means that we are not synchronized until the subsequent row
//	SocketedStreamingSolver< Channels , SyncType >::Set(-OffsetX(iters),0,0,0,0);
	SocketedStreamingSolver< Channels , SyncType >::Set(-OffsetX(iters)-1,0,0,0,0);
	// Set the iterations for the child so that it trails accordingly
	if(pChild)	pChild->SetProlongation();
#if DEBUG_FLAGS
	if( !debugFlags.noTempIO )
#endif // DEBUG_FLAGS
	{
		if( B ) B->reset( true , bSize );
		if( X ) X->reset( true , xSize );
		if( B )	B->SetServer( _server );
		if( X ) X->SetServer( _server );
	}
}
template< int Channels , class StorageType , class SyncType >
void SocketedMultiGridStreamingSolver< Channels , StorageType , SyncType >::SetRestriction( void )
{
	// Warning: This is actually somewhat agressive as you only really need something closer Size/2
//	SocketedStreamingSolver< Channels , SyncType >::Set(-OffsetX(iters),0,0,0,0,FiniteElements1D<float,Type,Degree>::FullProlongationStencil::ProlongationStencil::Size-OffsetR());
	SocketedStreamingSolver< Channels , SyncType >::Set(-OffsetX(iters),0,0,0,0,(FiniteElements1D<float,Type,Degree>::FullProlongationStencil::ProlongationStencil::Size+1)/2-OffsetR());
	// Set the iterations for the parent so that it trails accordingly
	if( parent ) parent->SetRestriction();
#if DEBUG_FLAGS
	if( !debugFlags.noTempIO )
#endif // DEBUG_FLAGS
	{
		if( B ) X->reset( false , xSize );
		if( X ) B->reset( false , bSize );
		if( B ) B->SetServer( _server );
		if( X ) X->SetServer( _server );
	}
	FreeArray( scratchR );
	scratchR = AllocArray< float >( _size*2 * Channels , ALIGNMENT , "SocketedMultiGridStreamingSolver::scratchR" );
}
template< int Channels , class StorageType , class SyncType >
void SocketedMultiGridStreamingSolver< Channels , StorageType , SyncType >::UnSetProlongation(void)
{
	SocketedStreamingSolver< Channels , SyncType >::UnSet();
#if DEBUG_FLAGS
	if(!debugFlags.noTempIO)
#endif // DEBUG_FLAGS
	{
		if( X ) X->unset();
		if( B ) B->unset();
	}
	FreeArray( prolongationStencil );
	if( pChild ) pChild->UnSetProlongation();
}
template< int Channels , class StorageType , class SyncType >
void SocketedMultiGridStreamingSolver< Channels , StorageType , SyncType >::UnSetRestriction(void)
{
	SocketedStreamingSolver< Channels , SyncType >::UnSet();
#if DEBUG_FLAGS
	if(!debugFlags.noTempIO)
#endif // DEBUG_FLAGS
	{
		if( X ) X->unset();
		if( B ) B->unset();
	}
	FreeArray( scratchR );
	if( parent ) parent->UnSetRestriction();
}

template< int Channels , class StorageType , class SyncType >
void SocketedMultiGridStreamingSolver< Channels , StorageType , SyncType >::SetRestriction( Pointer( float ) lB , int c , int idx , int sRestrict , int eRestrict )
{
#if DEBUG_FLAGS
	if(debugFlags.noRestriction)	return;
#endif // DEBUG_FLAGS
	int highStart = FiniteElements1D<float,Type,Degree>::FullProlongationStencil::ProlongationStencil::Start(sRestrict);
	int highEnd = FiniteElements1D<float,Type,Degree>::FullProlongationStencil::ProlongationStencil::Start(eRestrict)+FiniteElements1D<float,Type,Degree>::FullProlongationStencil::ProlongationStencil::Size;
	int highStart128 = (highStart-(RealPerWord-1))/RealPerWord;
	int highEnd128   = (highEnd  +(RealPerWord-1))/RealPerWord;
	highStart128 -= _start128;
	highEnd128 -= _start128;

#if SUPPORT_CYLINDRICAL
	if( (idx>=Degree && idx<(minor>>1)-Degree) || periodicType==SPHERICAL_PERIODIC )	SetInteriorRestriction(lB,c,idx,sRestrict,eRestrict);
#else // !SUPPORT_CYLINDRICAL
	if( (idx>=Degree && idx<(minor>>1)-Degree) || spherical)	SetInteriorRestriction(lB,c,idx,sRestrict,eRestrict);
#endif // SUPPORT_CYLINDRICAL
	else
	{
		int startY=FiniteElements1D<float,Type,Degree>::FullProlongationStencil::ProlongationStencil::Start(idx);
		int jj;
		if		(idx<Degree)				jj=idx;
		else if	(idx>(minor>>1)-1-Degree)	jj=2*Degree+(idx-((minor>>1)-1));
		else								jj=Degree;

		{
			ConstPointer( __m128 ) localRPtrs[FiniteElements1D<float,Type,Degree>::FullProlongationStencil::ProlongationStencil::Size];
			for(int y=0;y<FiniteElements1D<float,Type,Degree>::FullProlongationStencil::ProlongationStencil::Size;y++)
				if(startY+y>=0 && startY+y<minor)	localRPtrs[y] = ( Pointer( __m128 ) )GetRRow(startY+y,c);
				else								localRPtrs[y] = NullPointer< __m128 >( );
			__declspec (align(16)) float scratch[4];
			__m128 res[FiniteElements1D<float,Type,Degree>::FullProlongationStencil::ProlongationStencil::Size];
			float dotSum;
			__m128 dSum;
			int s,e;

			for(int d=0;d<FiniteElements1D<float,Type,Degree>::FullProlongationStencil::ProlongationStencil::Size;d++)
			{
				for(int j=0;j<4;j++)	scratch[j]= minorProlongationStencil.caseTable[jj].values[d];
				res[d]=_mm_load_ps(scratch);
			}
			memset( localRAccum+highStart128 , 0 , sizeof(__m128)*(highEnd128-highStart128 ) );
			for(int d=0;d<FiniteElements1D<float,Type,Degree>::FullProlongationStencil::ProlongationStencil::Size;d++)
				if(localRPtrs[d])
					for(int i=highStart128;i<highEnd128;i++) localRAccum[i]=_mm_add_ps(localRAccum[i],_mm_mul_ps(res[d],localRPtrs[d][i]));

			// Offset 0
			{
#if SUPPORT_CYLINDRICAL
				if( sRestrict==0 && periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
				if(sRestrict==0 && !spherical)
#endif // SUPPORT_CYLINDRICAL
				{
					dotSum=0;
					lB[0]=RestrictionUpdate0( prolongationStencil[0].matrixValues[0] , localRAccum , dotSum , 0 );
					s=2;
				}
				else
				{
					SetRestrictionDotSum(prolongationStencil[1].matrixValues[0],localRAccum,highStart128,dSum);
					_mm_store_ps(scratch,dSum);
					dotSum=scratch[3];
					s=sRestrict-_start/2;
				}

#if SUPPORT_CYLINDRICAL
				if( eRestrict==major/2 && periodicType==NO_PERIODIC )	e=(eRestrict-_start/2)-2;
#else // !SUPPORT_CYLINDRICAL
				if(eRestrict==major/2 && !spherical)	e=(eRestrict-_start/2)-2;
#endif // SUPPORT_CYLINDRICAL
				else									e=(eRestrict-_start/2);

				for(int i=s;i<e;i+=2)						lB[i]=RestrictionUpdate0(prolongationStencil[1].matrixValues[0],localRAccum,dotSum,i);
				int i=(eRestrict-_start/2)-2;
#if SUPPORT_CYLINDRICAL
				if( eRestrict==major/2 && periodicType==NO_PERIODIC )	lB[i]=RestrictionUpdate0(prolongationStencil[2].matrixValues[0],localRAccum,dotSum,i);
#else // !SUPPORT_CYLINDRICAL
				if(eRestrict==major/2 && !spherical)	lB[i]=RestrictionUpdate0(prolongationStencil[2].matrixValues[0],localRAccum,dotSum,i);
#endif // SUPPORT_CYLINDRICAL
			}
			// Offset 1
			{
#if SUPPORT_CYLINDRICAL
				if( sRestrict==0 && periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
				if(sRestrict==0 && !spherical)
#endif // SUPPORT_CYLINDRICAL
				{
					SetRestrictionDotSum(prolongationStencil[0].matrixValues[1],localRAccum,0,dSum);
					s=1;
				}
				else
				{
					SetRestrictionDotSum(prolongationStencil[1].matrixValues[1],localRAccum,highStart128+1,dSum);
					s=(sRestrict-_start/2)+1;
				}
				_mm_store_ps(scratch,dSum);
				dotSum=scratch[1]+scratch[2]+scratch[3];
#if SUPPORT_CYLINDRICAL
				if( eRestrict==major/2 && periodicType==NO_PERIODIC )	e=(eRestrict-_start/2)-4;
#else // !SUPPORT_CYLINDRICAL
				if(eRestrict==major/2 && !spherical)	e=(eRestrict-_start/2)-4;
#endif // SUPPORT_CYLINDRICAL
				else									e=(eRestrict-_start/2);
				for(int i=s;i<e;i+=2)	lB[i]=RestrictionUpdate1(prolongationStencil[1].matrixValues[1],localRAccum,dotSum,i);
#if SUPPORT_CYLINDRICAL
				if( eRestrict==major/2 && periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
				if(eRestrict==major/2 && !spherical)
#endif // SUPPORT_CYLINDRICAL
				{
					int i=(eRestrict-_start/2)-3;
					lB[i]=RestrictionUpdate1(prolongationStencil[2].matrixValues[1],localRAccum,dotSum,i);
					i+=2;
					lB[i]=dotSum;
				}
			}
		}
	}
}
template< int Channels , class StorageType , class SyncType >
void SocketedMultiGridStreamingSolver< Channels , StorageType , SyncType >::SetInteriorRestriction( Pointer( float ) lB , int c , int idx , int sRestrict , int eRestrict )
{
	int startY=FiniteElements1D<float,Type,Degree>::FullProlongationStencil::ProlongationStencil::Start(idx);
	int highStart = FiniteElements1D<float,Type,Degree>::FullProlongationStencil::ProlongationStencil::Start(sRestrict);
	int highEnd = FiniteElements1D<float,Type,Degree>::FullProlongationStencil::ProlongationStencil::Start(eRestrict)+FiniteElements1D<float,Type,Degree>::FullProlongationStencil::ProlongationStencil::Size;
	int highStart128 = (highStart-(RealPerWord-1))/RealPerWord;
	int highEnd128   = (highEnd  +(RealPerWord-1))/RealPerWord;
	highStart128 -= _start128;
	highEnd128 -= _start128;
	{
		ConstPointer( __m128 ) localRPtrs[FiniteElements1D<float,Type,Degree>::FullProlongationStencil::ProlongationStencil::Size];
		for(int y=0;y<FiniteElements1D<float,Type,Degree>::FullProlongationStencil::ProlongationStencil::Size;y++)	localRPtrs[y] = ( Pointer( __m128 ) )GetRRow(startY+y,c);

		__declspec (align(16)) float scratch[4];
		__m128 res[FiniteElements1D<float,Type,Degree>::FullProlongationStencil::ProlongationStencil::Size];
		float dotSum;
		__m128 dSum;
		int s,e;
		for(int d=0;d<FiniteElements1D<float,Type,Degree>::FullProlongationStencil::ProlongationStencil::Size;d++)
		{
			for(int j=0;j<4;j++)	scratch[j]= minorProlongationStencil.caseTable[Degree].values[d];
			res[d]=_mm_load_ps(scratch);
		}
		{
			int d=0;
			for(int i=highStart128;i<highEnd128;i++)	localRAccum[i]=_mm_mul_ps(res[d],_mm_add_ps(localRPtrs[d][i],localRPtrs[Degree+1-d][i]));
		}
		for(int d=1;d<=(Degree>>1);d++)
			for(int i=highStart128;i<highEnd128;i++)	localRAccum[i]=_mm_add_ps(localRAccum[i],_mm_mul_ps(res[d],_mm_add_ps(localRPtrs[d][i],localRPtrs[Degree+1-d][i])));
		if(Degree&1)
		{
			int d=(Degree+1)>>1;
			for(int i=highStart128;i<highEnd128;i++)	localRAccum[i]=_mm_add_ps(localRAccum[i],_mm_mul_ps(res[d],_mm_add_ps(localRPtrs[d][i],localRPtrs[Degree+1-d][i])));
		}
		// Offset 0
		{
#if SUPPORT_CYLINDRICAL
			if( sRestrict==0 && periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if(sRestrict==0 && !spherical)
#endif // SUPPORT_CYLINDRICAL
			{
				dotSum=0;
				lB[0] = RestrictionUpdate0(prolongationStencil[0].matrixValues[0],localRAccum,dotSum,0);
				s=2;
			}
			else
			{
				SetRestrictionDotSum(prolongationStencil[1].matrixValues[0],localRAccum,highStart128,dSum);
				_mm_store_ps(scratch,dSum);
				dotSum = scratch[3];
				s = sRestrict-_start/2;
			}
#if SUPPORT_CYLINDRICAL
			if( eRestrict==major/2 && periodicType==NO_PERIODIC )	e=(eRestrict-_start/2)-2;
#else // !SUPPORT_CYLINDRICAL
			if(eRestrict==major/2 && !spherical)	e=(eRestrict-_start/2)-2;
#endif // SUPPORT_CYLINDRICAL
			else									e=(eRestrict-_start/2);

			for(int i=s;i<e;i+=2)					lB[i]=RestrictionUpdate0(prolongationStencil[1].matrixValues[0],localRAccum,dotSum,i);
			int i=(eRestrict-_start/2)-2;
#if SUPPORT_CYLINDRICAL
			if( eRestrict==major/2 && periodicType==NO_PERIODIC )	lB[i]=RestrictionUpdate0(prolongationStencil[2].matrixValues[0],localRAccum,dotSum,i);
#else // !SUPPORT_CYLINDRICAL
			if(eRestrict==major/2 && !spherical)	lB[i]=RestrictionUpdate0(prolongationStencil[2].matrixValues[0],localRAccum,dotSum,i);
#endif // SUPPORT_CYLINDRICAL
		}
		// Offset 1
		{
#if SUPPORT_CYLINDRICAL
			if( sRestrict==0 && periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if(sRestrict==0 && !spherical)
#endif // SUPPORT_CYLINDRICAL
			{
				SetRestrictionDotSum(prolongationStencil[0].matrixValues[1],localRAccum,0,dSum);
				s=1;
			}
			else
			{
				SetRestrictionDotSum(prolongationStencil[1].matrixValues[1],localRAccum,highStart128+1,dSum);
				s=(sRestrict-_start/2)+1;
			}
			_mm_store_ps(scratch,dSum);
			dotSum=scratch[1]+scratch[2]+scratch[3];
#if SUPPORT_CYLINDRICAL
			if( eRestrict==major/2 && periodicType==NO_PERIODIC )	e=(eRestrict-_start/2)-4;
#else // !SUPPORT_CYLINDRICAL
			if(eRestrict==major/2 && !spherical)	e=(eRestrict-_start/2)-4;
#endif // SUPPORT_CYLINDRICAL
			else									e=(eRestrict-_start/2);
			for(int i=s;i<e;i+=2)	lB[i]=RestrictionUpdate1(prolongationStencil[1].matrixValues[1],localRAccum,dotSum,i);
#if SUPPORT_CYLINDRICAL
			if( eRestrict==major/2 && periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if(eRestrict==major/2 && !spherical)
#endif // SUPPORT_CYLINDRICAL
			{
				int i=(eRestrict-_start/2)-3;
				lB[i]=RestrictionUpdate1(prolongationStencil[2].matrixValues[1],localRAccum,dotSum,i);
				i+=2;
				lB[i]=dotSum;
			}
		}
	}
}

template< int Channels , class StorageType , class SyncType >
bool SocketedMultiGridStreamingSolver< Channels , StorageType , SyncType >::IterateRestriction( void )
{
	// index+OffsetR is the last index at which the residual is set after the solve.
	int restrictionIndex = FiniteElements1D< float , Type , Degree >::FullRestrictionStencil::RestrictionStencil::Start( index+OffsetR() );
	if( restrictionIndex>=minor/2 ) return false;
	int idx=index+OffsetB(iters);

#if TIME_IO
	double t = Time();
#endif // TIME_IO
	if( idx+Degree>=0 && idx+Degree<minor )
	{
		if( inX )
		{
			Pointer( SyncType ) inPtr = ( Pointer( SyncType ) )(*inX)[idx+Degree];
			for(int c=0;c<Channels;c++)
			{
				Pointer( SyncType ) _inPtr = inPtr + c*_size;
				Pointer( float ) x = GetXRow( idx+Degree , c );
				for( int i=0 ; i<_size ; i++ ) x[i] = float( _inPtr[i] ) * laplacianScaleR;
			}
			inX->advance();
		}
	}
	if( idx>=0 && idx<minor )
		if( inB )
		{
			Pointer( SyncType ) inPtr = ( Pointer( SyncType ) )(*inB)[idx];
			for( int c=0 ; c<Channels ; c++ )
			{
				Pointer( float ) bRow = GetBRow( idx , c );
				Pointer( SyncType ) _inPtr = inPtr + c*_size;
				for(int i=0;i<_size;i++) bRow[i] = float( _inPtr[i] );
			}
			inB->advance();
		}
#if TIME_IO
	rSync += Time() - t;
#endif // TIME_IO
	if( idx>=0 )
	{
		if( idx<minor )
		{
			if( !inB )
				if( !rChild )
				{
					fprintf(stderr,"Badness: no input stream\n");
					exit(0);
				}
				else
				{
					for(int c=0;c<Channels;c++)
						for(int l=0;l<laneNum;l++)
						{
							int b128 = _start128+((l  )*_size128)/laneNum;
							int e128 = _start128+((l+1)*_size128)/laneNum;
							int b = b128<<2;
							int e = e128<<2;
							rChild->SetRestriction( GetBRow(idx,c) , c , idx , b , e );
						}
				}
		}
		// Run an iteration of the Gauss-Seidel solver
		SocketedStreamingSolver< Channels , SyncType >::Solve();
	}
	// Write out the current solution
#if DEBUG_FLAGS
	if(!debugFlags.noTempIO)
#endif // DEBUG_FLAGS
	{
		if( UpdateXOutput< StorageType >(X) ) X->advance();
		if( UpdateBOutput< StorageType >(B) ) B->advance();
	}
#if 0
	// Write out the residual
	if(parent)	if(index>=startRestriction && (index&1)==restrictionBit)	parent->IterateRestriction();
#endif
#if TIME_IO
	t = Time();
#endif // TIME_IO
	if(!parent && outB && index>=0 && index<minor)
	{
		Pointer( SyncType ) outPtr = ( Pointer( SyncType ) )(*outB)[index];
		for(int c=0;c<Channels;c++)
		{
			Pointer( float ) inP = GetBRow( index , c );
			Pointer( SyncType ) outP = outPtr+c*_size;
			for( int i=0 ; i<_size ; i++ ) outP[i] = SyncType( inP[i] );
		}
		outB->advance();
	}
	if(!parent && outX && index>=0 && index<minor)
	{
		Pointer( SyncType ) outPtr = ( Pointer( SyncType ) )(*outX)[index];
		for(int c=0;c<Channels;c++)
		{
			Pointer( float ) inP = GetXRow(index,c);
			Pointer( SyncType ) outP = outPtr+c*_size;
			for( int i=0 ; i<_size ; i++ ) outP[i] = SyncType( inP[i]*laplacianScale );
		}
		outX->advance();
	}
#if SINGLE_OUTPUT
	if( setResidual && outR && restrictionIndex>=0 && (index&1)==restrictionBit )
#else // !SINGLE_OUTPUT
	if( setResidual && !parent && outR && restrictionIndex>=0 && (index&1)==restrictionBit )
#endif // SINGLE_OUTPUT
	{
		for( int c=0 ; c<Channels;c++ )
		{
			Pointer( float ) outP = scratchR + c*_size/2;
			for( int l=0 ; l<laneNum ; l++ )
			{
				int b128 = (_start128/2) + ( (l  ) * (_size128/2) ) / laneNum;
				int e128 = (_start128/2) + ( (l+1) * (_size128/2) ) / laneNum;
#if MISHA_FIX
				int b = (_start/2) + ( (l  ) * (_size/2) ) / laneNum;;
				int e = (_start/2) + ( (l+1) * (_size/2) ) / laneNum;
#else // !MISHA_FIX
				int b = b128 << 2;
				int e = e128 << 2;
#endif // MISHA_FIX
				SetRestriction( outP , c , restrictionIndex , b , e );
			}
		}
		Pointer( SyncType ) outPtr = ( Pointer( SyncType ) )(*outR)[restrictionIndex];
		for( int i=0 ; i<_size/2*Channels ; i++ ) outPtr[i] = SyncType( scratchR[i] );
//		memcpy( outPtr , scratchR , sizeof(float) * _size/2 * Channels );
		outR->advance();
	}
#if TIME_IO
	rSync += Time() - t;
#endif // TIME_IO
#if 1
	// Write out the residual
	if( parent ) if(index>=startRestriction && (index&1)==restrictionBit)	parent->IterateRestriction();
#endif
	return Increment();
}

template< int Channels , class StorageType , class SyncType >
bool SocketedMultiGridStreamingSolver< Channels , StorageType , SyncType >::IterateProlongation(void)
{
	int prolongationIndex = FiniteElements1D< float , Type , Degree >::FullProlongationStencil::ProlongationStencil::Start( index + Degree - 1  );
	if( prolongationIndex>=minor*2 && index+OffsetR()>=minor ) return false;
	int idx = index + OffsetX( iters );
	if( idx>=0 )
	{
#if DEBUG_FLAGS
		if(!debugFlags.noTempIO)
#endif // DEBUG_FLAGS
		{
			if( UpdateBInput< StorageType>( B ) ) B->advance();
			if( UpdateXInput< StorageType>( X ) ) X->advance();
		}
#if SINGLE_OUTPUT
		if( inX && idx<minor )
#else // !SINGLE_OUTPUT
		if(!parent && inX && idx<minor)
#endif // SINGLE_OUTPUT
		{
#if TIME_IO
			double t = Time();
#endif // TIME_IO
			Pointer( SyncType ) inPtr = ( Pointer( SyncType ) )(*inX)[idx];
			for( int c=0 ; c<Channels ; c++ )
			{
				Pointer( float ) x = GetXRow( idx , c );
				Pointer( SyncType ) _inPtr = inPtr + c*_size;
				for( int i=0 ; i<_size ; i++ )	x[i] += float( _inPtr[i] )*laplacianScaleR;
			}
			inX->advance();
#if TIME_IO
			rSync += Time() - t;
#endif // TIME_IO
		}
		// Run an iteration of the Gauss-Seidel solver
#if SINGLE_OUTPUT
		if( idx<minor && parent && !inX )
#else // !SINGLE_OUTPUT
		if( idx<minor && parent )
#endif // SINGLE_OUTPUT
			for( int c=0 ; c<Channels ; c++ )
				parent->SetProlongation( GetXRow( idx , c ) , c , idx , _start , _size , major , minor , parent->laplacianScale * laplacianScaleR );
		SocketedStreamingSolver< Channels , SyncType >::Solve();
		if(!pChild && index>=0 && index<minor)
		{
			{
#if TIME_IO
				double t = Time();
#endif // TIME_IO
				if( outX )
				{
					Pointer( SyncType ) outPtr = ( Pointer( SyncType ) )(*outX)[index];
					for( int c=0 ; c<Channels ; c++ )
					{
						Pointer( float ) x = GetXRow(index,c);
						Pointer( SyncType ) _outPtr = outPtr + c*_size;
						for( int i=0 ; i<_size ; i++ ) _outPtr[i] = SyncType( x[i]*laplacianScale );
					}
					outX->advance();
				}
#if TIME_IO
				rSync += Time() - t;
#endif // TIME_IO
			}
		}
	}
#if TIME_IO
	double t = Time();
#endif // TIME_IO
	for( int off=0 ; off<2 ; off++ )
#if SINGLE_OUTPUT
		if( outR && prolongationIndex+off>=0 && prolongationIndex+off<2*minor )
#else // !SINGLE_OUTPUT
		if( !pChild && outR && prolongationIndex+off>=0 && prolongationIndex+off<2*minor )
#endif // SINGLE_OUTPUT
		{
			Pointer( SyncType ) outPtr = ( Pointer( SyncType ) )( (*outR)[prolongationIndex+off] );
			memset( scratchR , 0 , _size * 2 * Channels * sizeof( float ) );
			for( int c=0 ; c<Channels ; c++ ) SetProlongation( scratchR + c*_size*2 , c , prolongationIndex+off , _start*2 , _size*2 , major*2 , minor*2 , laplacianScale );
			for( int i=0 ; i<_size*2*Channels ; i++ ) outPtr[i] = SyncType( scratchR[i] );
			outR->advance();
		}
#if TIME_IO
	rSync += Time() - t;
#endif // TIME_IO

	if( pChild && index>=startProlongation )
		if( pChild->IterateProlongation() )
			pChild->IterateProlongation();
	return Increment();
}
template< int Channels , class StorageType , class SyncType >
void SocketedMultiGridStreamingSolver< Channels , StorageType , SyncType >::SetProlongation( Pointer( float ) highX , int c , int highIdx , int highStart , int highSize , int highMajor , int highMinor , double scale )
{
#if DEBUG_FLAGS
	if( debugFlags.noProlongation ) return;
#endif // DEBUG_FLAGS
	int startY = RestrictionStencil::RestrictionStencil::Start( highIdx );
	int jj;
#if SUPPORT_CYLINDRICAL
	if( periodicType!=SPHERICAL_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
	if( !spherical )
#endif // SUPPORT_CYLINDRICAL
		if		( highIdx<Degree )				jj = highIdx;
		else if	( highIdx>highMinor-1-Degree )	jj = 2*Degree+1+(highIdx-(highMinor-1));
		else if ( (highIdx-Degree)&1 )			jj = Degree+1;
		else									jj = Degree;
	else
		if( (highIdx-Degree)&1 )	jj = Degree+1;
		else						jj = Degree;
	ConstPointer( float ) xPtrs[ RestrictionStencil::RestrictionStencil::Size ];
	for( int yy=0 ; yy<RestrictionStencil::RestrictionStencil::Size ; yy++ )
#if SUPPORT_CYLINDRICAL
		if( (yy+startY>=0 && yy+startY<minor) || periodicType==SPHERICAL_PERIODIC )	xPtrs[yy] = GetXRow( yy+startY , c );
#else // !SUPPORT_CYLINDRICAL
		if( (yy+startY>=0 && yy+startY<minor) || spherical )	xPtrs[yy] = GetXRow( yy+startY , c );
#endif // SUPPORT_CYLINDRICAL
		else													xPtrs[yy] = NullPointer< float >( );
		for( int i=0 ; i<highSize ; i++ )
		{
			int ii;
#if SUPPORT_CYLINDRICAL
			if( periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if( !spherical )
#endif // SUPPORT_CYLINDRICAL
				if		( i+highStart<Degree)				ii = i+highStart;
				else if	( i+highStart>highMajor-1-Degree)	ii = 2*Degree+1+(i+highStart-(highMajor-1));
				else if ((i+highStart-Degree)&1)			ii = Degree+1;
				else										ii = Degree;
			else
				if ( (i+highStart-Degree)&1 )	ii = Degree+1;
				else							ii = Degree;

			double value=0;

			int startX = RestrictionStencil::RestrictionStencil::Start(i+highStart) - _start;
			for( int yy=0 ; yy<RestrictionStencil::RestrictionStencil::Size ; yy++ )
			{
				if( !xPtrs[yy] )	continue;
				double tValue=0;

				for( int xx=0 ; xx<RestrictionStencil::RestrictionStencil::Size ; xx++ )
#if SUPPORT_CYLINDRICAL
					if( (xx+startX+_start>=0 && xx+startX+_start<major) || periodicType!=NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
					if( (xx+startX+_start>=0 && xx+startX+_start<major) || spherical )
#endif // SUPPORT_CYLINDRICAL
						tValue += xPtrs[yy][startX+xx]*majorRestrictionStencil.caseTable[ii].values[xx];
				value += tValue*minorRestrictionStencil.caseTable[jj].values[yy];
			}
			highX[i] += value * scale;
		}
}

template< int Channels , class StorageType , class SyncType >
void SocketedMultiGridStreamingSolver< Channels , StorageType , SyncType >::ProlongationUpdate( int j , double scale )
{
#if DEBUG_FLAGS
	if(debugFlags.noProlongation)	return;
#endif // DEBUG_FLAGS
	int startY = RestrictionStencil::RestrictionStencil::Start(j);
	int jj;
#if SUPPORT_CYLINDRICAL
	if( periodicType!=SPHERICAL_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
	if( !spherical )
#endif // SUPPORT_CYLINDRICAL
		if		( j<Degree )			jj = j;
		else if	( j>minor-1-Degree )	jj = 2*Degree+1+(j-(minor-1));
		else if ( (j-Degree)&1 )		jj = Degree+1;
		else							jj = Degree;
	else
		if( (j-Degree)&1 )	jj=Degree+1;
		else				jj=Degree;

	for( int c=0 ; c<Channels ; c++ )
	{
		const float* parentXPtrs[RestrictionStencil::RestrictionStencil::Size];
		float* localXPtr = GetXRow( j , c );
		for( int yy=0 ; yy<RestrictionStencil::RestrictionStencil::Size ; yy++ )
#if SUPPORT_CYLINDRICAL
			if( (yy+startY>=0 && yy+startY<parent->minor) || periodicType==SPHERICAL_PERIODIC )	parentXPtrs[yy]=parent->GetXRow(yy+startY,c);
#else // !SUPPORT_CYLINDRICAL
			if( (yy+startY>=0 && yy+startY<parent->minor) || spherical )	parentXPtrs[yy]=parent->GetXRow(yy+startY,c);
#endif // SUPPORT_CYLINDRICAL
			else															parentXPtrs[yy]=NULL;
		for(int i=0;i<_size;i++)
		{
			int ii;
#if SUPPORT_CYLINDRICAL
			if( periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if(!spherical)
#endif // SUPPORT_CYLINDRICAL
				if		(i+_start<Degree)			ii=i+_start;
				else if	(i+_start>major-1-Degree)	ii=2*Degree+1+(i+_start-(major-1));
				else if ((i+_start-Degree)&1)		ii=Degree+1;
				else								ii=Degree;
			else
				if ((i+_start-Degree)&1)	ii=Degree+1;
				else						ii=Degree;

			double value=0;

			int startX=RestrictionStencil::RestrictionStencil::Start(i+_start)-parent->_start;
			for(int yy=0;yy<RestrictionStencil::RestrictionStencil::Size;yy++)
			{
				if(!parentXPtrs[yy])	continue;
				double tValue=0;

				for(int xx=0;xx<RestrictionStencil::RestrictionStencil::Size;xx++)
#if SUPPORT_CYLINDRICAL
					if( (xx+startX+parent->_start>=0 && xx+startX+parent->_start<parent->major) || periodicType!=NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
					if( (xx+startX+parent->_start>=0 && xx+startX+parent->_start<parent->major) || spherical)
#endif // SUPPORT_CYLINDRICAL
						tValue+=parentXPtrs[yy][startX+xx]*majorRestrictionStencil.caseTable[ii].values[xx];
				value+=tValue*minorRestrictionStencil.caseTable[jj].values[yy];
			}
			localXPtr[i] += value * scale;
		}
	}
}

template< int Channels , class StorageType , class SyncType >
void SocketedMultiGridStreamingSolver< Channels , StorageType , SyncType >::SolveRestriction(void)
{
	// Run to completion...
	while( IterateRestriction() ){;}
	// ...and finish up the trailing parent
	if(parent)	parent->SolveRestriction();

	if( !inCore ) // Why only for out-of-core?
	{
		MultiStreamIOClient *b = (MultiStreamIOClient*) B;
		MultiStreamIOClient *x = (MultiStreamIOClient*) X;
		while( ( b && b->server ) || ( x && x->server ) ) Sleep( 0 );
	}
}
template< int Channels , class StorageType , class SyncType >
void SocketedMultiGridStreamingSolver< Channels , StorageType , SyncType >::SolveProlongation(void)
{
	// Run to completion...
	while( IterateProlongation() ){;}
	// ...and finish up the trailing child
	if(pChild)	pChild->SolveProlongation();
	if( !inCore ) // Why only for out-of-core?
	{
		MultiStreamIOClient *b = (MultiStreamIOClient*) B;
		MultiStreamIOClient *x = (MultiStreamIOClient*) X;
		while( ( b && b->server ) || ( x && x->server ) ) Sleep( 0 );
	}
}

///////////////
// LabelData //
///////////////
template<class LabelType,int Channels>
inline bool LabelData<LabelType,Channels>::operator == (const LabelData& ld) const
{
	bool ret = true;
	for( int c=0 ; c<Channels ; c++ ) ret &= (ld.l[c]==l[c]);
	return ret;
}
template<> inline bool LabelData<unsigned char		,1>::operator == (const LabelData& ld) const { return (ld.l[0] == l[0]); }
template<> inline bool LabelData<unsigned char		,2>::operator == (const LabelData& ld) const { return (ld.l[0] == l[0]) && (ld.l[1] == l[1]); }
template<> inline bool LabelData<unsigned char		,3>::operator == (const LabelData& ld) const { return (ld.l[0] == l[0]) && (ld.l[1] == l[1]) && (ld.l[2] == l[2]); }
template<> inline bool LabelData<unsigned __int16	,1>::operator == (const LabelData& ld) const { return (ld.l[0] == l[0]); }
template<> inline bool LabelData<unsigned __int16	,2>::operator == (const LabelData& ld) const { return (ld.l[0] == l[0]) && (ld.l[1] == l[1]); }
template<> inline bool LabelData<unsigned __int16	,3>::operator == (const LabelData& ld) const { return (ld.l[0] == l[0]) && (ld.l[1] == l[1]) && (ld.l[2] == l[2]); }
template<> inline bool LabelData<unsigned int		,1>::operator == (const LabelData& ld) const { return (ld.l[0] == l[0]); }
template<> inline bool LabelData<unsigned int		,2>::operator == (const LabelData& ld) const { return (ld.l[0] == l[0]) && (ld.l[1] == l[1]); }
template<> inline bool LabelData<unsigned int		,3>::operator == (const LabelData& ld) const { return (ld.l[0] == l[0]) && (ld.l[1] == l[1]) && (ld.l[2] == l[2]); }
template<> inline bool LabelData<unsigned long long	,1>::operator == (const LabelData& ld) const { return (ld.l[0] == l[0]); }
template<> inline bool LabelData<unsigned long long	,2>::operator == (const LabelData& ld) const { return (ld.l[0] == l[0]) && (ld.l[1] == l[1]); }
template<> inline bool LabelData<unsigned long long	,3>::operator == (const LabelData& ld) const { return (ld.l[0] == l[0]) && (ld.l[1] == l[1]) && (ld.l[2] == l[2]); }
template< class LabelType , int Channels >
LabelData< LabelType , Channels > LabelData< LabelType , Channels >::BlackLabel( void )
{
	LabelData ld;
	for( int c=0 ; c<Channels ; c++ ) ld.l[c] = LabelType( -1 );
	return ld;
}

/////////////////////////
// StreamingDivergence //
/////////////////////////
template< int Channels , class PixelType , class LabelType , class StorageType , class SyncType >
SocketedStreamingDivergence< Channels , PixelType , LabelType , StorageType , SyncType >::SocketedStreamingDivergence(void)
{
	blackOut = true;
#if TIME_IO
	vSync = hSync = 0;
#endif // TIME_IO
	for( int i=0 ; i<ISize ; i++ ) labels[i] = NullPointer< LabelData< LabelType , Channels > >( );
	for( int i=0 ; i<ISize*Channels ; i++ )  pixels[i] = NullPointer< float >( );
	for( int i=0 ; i<DXSize*Channels ; i++ ) dx[i] = NullPointer< __m128 >( );
	for( int i=0 ; i<DYSize*Channels ; i++ ) dy[i] = NullPointer< __m128 >( );

	parent = NULL;


	for( int c=0 ; c<Channels ; c++ ) average[c] = meanB[c] = 0;
	syncBuffer = NullPointer< ImageData< SyncType , LabelType > >( );
	localDMajorAccum = NullPointer< __m128 >( );
	localDMinorAccum = NullPointer< __m128 >( );
}
template< int Channels , class PixelType , class LabelType , class StorageType , class SyncType >
SocketedStreamingDivergence< Channels , PixelType , LabelType , StorageType , SyncType >::~SocketedStreamingDivergence(void)
{
	FreeArray( syncBuffer );
	for( int i=0 ; i<ISize ; i++ ) FreeArray( labels[i]  );
	for( int i=0 ; i<ISize*Channels ; i++ ) FreeArray( pixels[i] );
	for( int i=0 ; i<DXSize*Channels ; i++ ) FreeArray( dx[i] );
	for( int i=0 ; i<DYSize*Channels ; i++ ) FreeArray( dy[i] );

	if( localDMajorAccum )
	{
		localDMajorAccum -= _padSize128;
		FreeArray( localDMajorAccum );
	}
	if( localDMinorAccum )
	{
		localDMinorAccum -= _padSize128;
		FreeArray( localDMinorAccum );
	}
}

template< int Channels , class PixelType , class LabelType , class StorageType , class SyncType >
void SocketedStreamingDivergence< Channels , PixelType , LabelType , StorageType , SyncType >::Initialize( StreamingGrid *pxls , StreamingGrid* lbls ,
																					    int start , int end , int major , int minor , int iters ,
																					    DataStream* leftStream,SOCKET* syncSockets,DataStream* rightStream,
#if SUPPORT_CYLINDRICAL
																					    bool memoryMappedFile , int periodicType , MultiStreamIOServer* server
#else // !SUPPORT_CYLINDRICAL
																					    bool memoryMappedFile,bool spherical , MultiStreamIOServer* server
#endif // SUPPORT_CYLINDRICAL
																					  )
{
	if(syncSockets)
	{
#if SUPPORT_CYLINDRICAL
		if( !SetSocketedStreamData( major , minor , start , end , leftStream , syncSockets[0] , rightStream , periodicType ) )	exit(0);
#else // !SUPPORT_CYLINDRICAL
		if(!SetSocketedStreamData(major,minor,start,end,leftStream,syncSockets[0],rightStream,spherical))	exit(0);
#endif // SUPPORT_CYLINDRICAL
	}
	else
	{
#if SUPPORT_CYLINDRICAL
		if( !SetSocketedStreamData( major , minor , start , end , leftStream , INVALID_SOCKET , rightStream , periodicType ) )	exit(0);
#else // !SUPPORT_CYLINDRICAL
		if(!SetSocketedStreamData(major,minor,start,end,leftStream,INVALID_SOCKET,rightStream,spherical))	exit(0);
#endif / SUPPORT_CYLINDRICAL
	}
#if SUPPORT_CYLINDRICAL
	this->_periodicType = periodicType;
#else // !SUPPORT_CYLINDRICAL
	this->_spherical=spherical;
#endif // SUPPORT_CYLINDRICAL
	this->major=major;
	this->minor=minor;
	pixelStream = pxls;
	labelStream = lbls;

	if(parent)
		if(syncSockets)
			parent->Initialize( start , end , major , minor , iters ,
			leftStream , syncSockets+1 , rightStream ,
#if SUPPORT_CYLINDRICAL
			memoryMappedFile , periodicType , server
#else // !SUPPORT_CYLINDRICAL
			memoryMappedFile , spherical , server
#endif // SUPPORT_CYLINDRICAL
			);
		else
			parent->Initialize( start , end , major , minor , iters ,
			leftStream , NULL , rightStream ,
#if SUPPORT_CYLINDRICAL
			memoryMappedFile , periodicType , server
#else // !SUPPORT_CYLINDRICAL
			memoryMappedFile , spherical , server
#endif // SUPPORT_CYLINDRICAL
			);
	for( int i=0 ; i<ISize ; i++ ) labels[i] = AllocArray< LabelData< LabelType , Channels > >( _paddedSize , 1 , "SocketedStreamingDivergence::labels" );
	for( int i=0 ; i<ISize*Channels ; i++ ) pixels[i] = AllocArray< float >( _paddedSize , 1 , "SocketedStreamingDivergence::pixels" );

	for( int i=0 ; i<DXSize*Channels ; i++ ) dx[i] = AllocArray< __m128 >( _paddedSize128 , ALIGNMENT , "SocketedStreamingDivergence::dx" );
	for( int i=0 ; i<DYSize*Channels ; i++ ) dy[i] = AllocArray< __m128 >( _paddedSize128 , ALIGNMENT , "SocketedStreamingDivergence::dy" );

	syncBuffer = AllocArray< ImageData< SyncType , LabelType > >( _size * Channels * Degree , 1 , "SocketedStreamingDivergence::syncBuffer" );
	localDMajorAccum = AllocArray< __m128 >( _paddedSize128 , ALIGNMENT , "SocketedStreamingDivergence::localDMajorAccum" );
	localDMinorAccum = AllocArray< __m128 >( _paddedSize128 , ALIGNMENT , "SocketedStreamingDivergence::localDMinorAccum" );
	localDMajorAccum += _padSize128;
	localDMinorAccum += _padSize128;

}
template< int Channels , class PixelType , class LabelType , class StorageType , class SyncType >
void SocketedStreamingDivergence< Channels , PixelType , LabelType , StorageType , SyncType >::InitRestriction(void)
{
	index = -Degree-1;
	_iHeight=pixelStream->rows();
	_iWidth=pixelStream->rowSize()/(Channels*sizeof(PixelType));

	if(pixelStream->rowSize()!=_iWidth*Channels*sizeof(PixelType))
		fprintf(stderr,"Pixel width failure: %d != %d\n",pixelStream->rowSize(),_iWidth*Channels*sizeof(PixelType)) ,	exit(0);
	if(_iHeight!=labelStream->rows())
		fprintf(stderr,"Label height failure: %d != %d\n",_iHeight,labelStream->rows()) , 								exit(0);
	if(labelStream->rowSize()!=_iWidth*Channels*sizeof(LabelType))
		fprintf(stderr,"Label width failure: %d != %d\n",labelStream->rowSize(),_iWidth*Channels*sizeof(LabelType)) ,	exit(0);

	pixelStream->reset(true,1);
	labelStream->reset(true,1);
	FiniteElements1D<float,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(major,dotMajorStencil,0,0);
	FiniteElements1D<float,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(minor,dotMinorStencil,0,0);
	FiniteElements1D<float,DERIVATIVE(Type),Degree-1>::DotProduct<Type,Degree>::DotProductStencil(FiniteElements1D<float,DERIVATIVE(Type),Degree-1>::Dimension(FiniteElements1D<float,Type,Degree>::DomainSize(major)),dDotMajorStencil,0,1);
	FiniteElements1D<float,DERIVATIVE(Type),Degree-1>::DotProduct<Type,Degree>::DotProductStencil(FiniteElements1D<float,DERIVATIVE(Type),Degree-1>::Dimension(FiniteElements1D<float,Type,Degree>::DomainSize(minor)),dDotMinorStencil,0,1);
	FiniteElements2D<float,Type,Degree>::DivergenceStencil(major,minor,divergenceStencil);

	if(parent)	parent->InitRestriction();
}

template< int Channels , class PixelType , class LabelType , class StorageType , class SyncType >
Pointer( float ) SocketedStreamingDivergence< Channels , PixelType , LabelType , StorageType , SyncType >::GetPixelRow(int row,int channel)
{
	return pixels[MyModIndex(row*Channels+channel,ISize*Channels)] + _padSize;
}
template< int Channels , class PixelType , class LabelType , class StorageType , class SyncType >
Pointer( LabelData< LabelType , Channels > ) SocketedStreamingDivergence< Channels , PixelType , LabelType , StorageType , SyncType >::GetLabelRow( int row ){ return labels[ MyModIndex( row , ISize ) ] + _padSize; }

template< int Channels , class PixelType , class LabelType , class StorageType , class SyncType >
Pointer( float ) SocketedStreamingDivergence< Channels , PixelType , LabelType , StorageType , SyncType >::GetDXRow(int row,int channel)
{
	return ( Pointer( float ) )( dx[MyModIndex(row*Channels + channel , DXSize*Channels)] + _padSize128 );
}
template< int Channels , class PixelType , class LabelType , class StorageType , class SyncType >
Pointer( float ) SocketedStreamingDivergence< Channels , PixelType , LabelType , StorageType , SyncType >::GetDYRow(int row,int channel)
{
	return ( Pointer( float ) )( dy[MyModIndex(row*Channels + channel , DYSize*Channels )] + _padSize128 );
}
template< int Channels , class PixelType , class LabelType , class StorageType , class SyncType >
void SocketedStreamingDivergence< Channels , PixelType , LabelType , StorageType , SyncType >::SyncImageHead( int idx , bool read )
{
#if DEBUG_FLAGS
	if( debugFlags.noTopBottom ) return;
#endif // DEBUG_FLAGS
	if( idx<0 || idx>=Degree )
	{
		fprintf( stderr , "Synchronizing divergence in non-head row: %d\n" , idx );
		return;
	}
	int r = -idx-1;
	if( syncXSocket!=INVALID_SOCKET )
	{
		int size=_size * Channels;
		if( read )
		{
#if TIME_IO
			double t = Time();
#endif // TIME_IO
			ReceiveOnSocket( syncXSocket, syncBuffer , sizeof( ImageData< SyncType , LabelType > ) * size , "SocketedStreamingDivergence< Channels , PixelType , LabelType , StorageType , SyncType >::SyncImageHead" );
#if TIME_IO
			vSync += Time() - t;
#endif // TIME_IO
			Pointer( LabelData< LabelType , Channels > ) labelRow = GetLabelRow( r );
			for(int c=0;c<Channels;c++)
			{
				Pointer( float ) pixelRow = GetPixelRow( r , c );
				for(int x=0;x<_size;x++)
				{
					pixelRow[x] = float( syncBuffer[_size*c+x].pixel );
					labelRow[x].l[c] = syncBuffer[_size*c+x].label;
				}
			}
		}
		else
		{
			Pointer( LabelData< LabelType , Channels > ) labelRow = GetLabelRow( idx );
			for(int c=0;c<Channels;c++)
			{
				Pointer( float ) pixelRow = GetPixelRow( idx , c );
				for(int x=0;x<_size;x++)
				{
					syncBuffer[_size*c+x].pixel = SyncType( pixelRow[x] );
					syncBuffer[_size*c+x].label = labelRow[x].l[c];
				}
			}
			SendOnSocket(syncXSocket,syncBuffer, sizeof( ImageData< SyncType , LabelType > ) * size , "SocketedStreamingDivergence< Channels , PixelType , LabelType , StorageType , SyncType >::SyncImageHead" );
		}
	}
}
template< int Channels , class PixelType , class LabelType , class StorageType , class SyncType >
void SocketedStreamingDivergence< Channels , PixelType , LabelType , StorageType , SyncType >::SyncImageTail( int idx , bool read )
{
#if DEBUG_FLAGS
	if( debugFlags.noTopBottom ) return;
#endif // DEBUG_FLAGS
	if( idx<minor-Degree || idx>=minor )
	{
		fprintf( stderr , "Synchronizing divergence in non-tail row: %d\n" , idx );
		return;
	}
	int r = 2 * minor - 1 - idx;
	if( syncXSocket!=INVALID_SOCKET )
	{
		int size=_size * Channels;
		if( read )
		{
#if TIME_IO
			double t = Time();
#endif // TIME_IO
			ReceiveOnSocket( syncXSocket, syncBuffer , sizeof( ImageData< SyncType , LabelType > ) * size , "SocketedStreamingDivergence< Channels , PixelType , LabelType , StorageType , SyncType >::SyncImageTail" );
#if TIME_IO
			vSync += Time() - t;
#endif // TIME_IO
			Pointer( LabelData< LabelType , Channels > ) labelRow = GetLabelRow( r );
			for(int c=0;c<Channels;c++)
			{
				Pointer( float ) pixelRow = GetPixelRow( r , c );
				for(int x=0;x<_size;x++)
				{
					pixelRow[x] = float( syncBuffer[_size*c+x].pixel );
					labelRow[x].l[c] = syncBuffer[_size*c+x].label;
				}
			}
		}
		else
		{
			Pointer( LabelData< LabelType , Channels > ) labelRow = GetLabelRow( idx );
			for(int c=0;c<Channels;c++)
			{
				Pointer( float ) pixelRow = GetPixelRow( idx , c );
				for(int x=0;x<_size;x++)
				{
					syncBuffer[_size*c+x].pixel = SyncType( pixelRow[x] );
					syncBuffer[_size*c+x].label = labelRow[x].l[c];
				}
			}
			SendOnSocket(syncXSocket,syncBuffer, sizeof( ImageData< SyncType , LabelType > ) * size , "SocketedStreamingDivergence< Channels , PixelType , LabelType , StorageType , SyncType >::SyncImageTail" );
		}
	}
}

template< int Channels , class PixelType , class LabelType , class StorageType , class SyncType >
void SocketedStreamingDivergence< Channels , PixelType , LabelType , StorageType , SyncType >::SyncImageLeft( int j , bool read )
{
#if DEBUG_FLAGS
	if( debugFlags.noLeftRight ) return;
#endif // DEBUG_FLAGS
	if( leftStream )
	{
		int xSize = 2*_padSize;
		int offset = ( sizeof( SyncType )*Channels + sizeof(LabelData<LabelType,Channels>) ) * _padSize;
		Pointer( byte ) left = AllocArray< byte >( ( sizeof( SyncType )*Channels + sizeof(LabelData<LabelType,Channels>) ) * xSize , 1 , "SocketedStreamingDivergence::SyncImageLeft (left)" );
		if(read)
		{
#if TIME_IO
			double t = Time();
#endif // TIME_IO
			if ( !leftStream->read( left , ( sizeof( SyncType )*Channels + sizeof(LabelData<LabelType,Channels>) ) * xSize ) )	exit(0);
#if TIME_IO
			hSync += Time() - t;
#endif // TIME_IO
			for( int c=0 ; c<Channels ; c++ )
			{
				Pointer( float ) pRow = GetPixelRow( j , c ) - _padSize;
				Pointer( SyncType ) lRow = ( Pointer( SyncType ) )( left + ( c * _padSize * sizeof( SyncType ) ) );
				for( int i=0 ; i<_padSize ; i++ ) pRow[i] = float( lRow[i] );
			}
			memcpy(  GetLabelRow(j) - _padSize , left + Channels * _padSize * sizeof( SyncType ) , _padSize * sizeof(LabelData<LabelType,Channels>) );
			for( int c=0 ; c<Channels ; c++ )
			{
				Pointer( float ) pRow = GetPixelRow( j+1 , c ) - _padSize;
				Pointer( SyncType ) lRow = ( Pointer( SyncType ) )( left + offset + c * _padSize * sizeof( SyncType ) );
				for( int i=0 ; i<_padSize ; i++ ) pRow[i] = float( lRow[i] );
			}
			memcpy (  GetLabelRow(j+1) - _padSize , left + offset + Channels * _padSize * sizeof( SyncType ) , _padSize * sizeof(LabelData<LabelType,Channels>) );
		}
		else
		{
			for( int c=0 ; c<Channels ; c++ )
			{
				Pointer( float ) pRow = GetPixelRow( j , c );
				Pointer( SyncType ) lRow = ( Pointer( SyncType ) )( left + ( c * _padSize * sizeof( SyncType ) ) );
				for( int i=0 ; i<_padSize ; i++ ) lRow[i] = SyncType( pRow[i] );
			}
			memcpy ( left + Channels * _padSize * sizeof( SyncType ) , GetLabelRow(j) , _padSize * sizeof(LabelData<LabelType,Channels>) );
			for( int c=0 ; c<Channels ; c++ )
			{
				Pointer( float ) pRow = GetPixelRow( j+1 , c );
				Pointer( SyncType ) lRow = ( Pointer( SyncType ) )( left + offset + c * _padSize * sizeof( SyncType ) );
				for( int i=0 ; i<_padSize ; i++ ) lRow[i] = SyncType( pRow[i] );
			}
			memcpy ( left + offset + Channels * _padSize * sizeof( SyncType ) , GetLabelRow(j+1) , _padSize * sizeof(LabelData<LabelType,Channels>) );
			if ( !leftStream->write( left , ( sizeof( SyncType )*Channels + sizeof(LabelData<LabelType,Channels>) ) * xSize ) )	exit(0);
		}
		FreeArray( left );
	}
}
template< int Channels , class PixelType , class LabelType , class StorageType , class SyncType >
void SocketedStreamingDivergence< Channels , PixelType , LabelType , StorageType , SyncType >::SyncImageRight(int j,bool read)
{
#if DEBUG_FLAGS
	if( debugFlags.noLeftRight ) return;
#endif // DEBUG_FLAGS
	if( rightStream )
	{
		int xSize = 2*_padSize;
		int offset = ( sizeof( SyncType )*Channels + sizeof(LabelData<LabelType,Channels>) ) * _padSize;
		Pointer( byte ) right = AllocArray< byte >( ( sizeof( SyncType )*Channels + sizeof(LabelData<LabelType,Channels>) ) * xSize , 1 , "SocketedStreamingDivergence::SyncImageRight (right)" );
		if( read )
		{
#if TIME_IO
			double t = Time();
#endif // TIME_IO
			if ( !rightStream->read( right , ( sizeof( SyncType )*Channels + sizeof(LabelData<LabelType,Channels>) ) * xSize ) )	exit(0);
#if TIME_IO
			hSync += Time() - t;
#endif // TIME_IO
			for( int c=0 ; c<Channels ; c++ )
			{
				Pointer( float ) pRow = GetPixelRow( j , c ) + _size;
				Pointer( SyncType ) rRow = ( Pointer( SyncType ) )( right + ( c * _padSize * sizeof( SyncType ) ) );
				for( int i=0 ; i<_padSize ; i++ ) pRow[i] = float( rRow[i] );
			}
			memcpy (  GetLabelRow(j) + _size , right + Channels * _padSize * sizeof( SyncType ) , _padSize * sizeof(LabelData<LabelType,Channels>) );
			for( int c=0 ; c<Channels ; c++ )
			{
				Pointer( float ) pRow = GetPixelRow( j+1 , c ) + _size;
				Pointer( SyncType ) rRow = ( Pointer( SyncType ) )( right + offset + c * _padSize * sizeof( SyncType ) );
				for( int i=0 ; i<_padSize ; i++ ) pRow[i] = float( rRow[i] );
			}
			memcpy (  GetLabelRow(j+1) + _size , right + offset + Channels * _padSize * sizeof( SyncType ) , _padSize * sizeof(LabelData<LabelType,Channels>) );
		}
		else
		{
			for( int c=0 ; c<Channels ; c++ )
			{
				Pointer( float ) pRow = GetPixelRow( j , c ) + _size - _padSize;
				Pointer( SyncType ) rRow = ( Pointer( SyncType ) )( right + ( c * _padSize * sizeof( SyncType ) ) );
				for( int i=0 ; i<_padSize ; i++ ) rRow[i] = SyncType( pRow[i] );
			}
			memcpy ( right + Channels * _padSize * sizeof( SyncType ) , GetLabelRow(j) + _size - _padSize , _padSize * sizeof(LabelData<LabelType,Channels>) );
			for( int c=0 ; c<Channels ; c++ )
			{
				Pointer( float ) pRow = GetPixelRow( j+1 , c ) + _size - _padSize;
				Pointer( SyncType ) rRow = ( Pointer( SyncType ) )( right + offset + c * _padSize * sizeof( SyncType ) );
				for( int i=0 ; i<_padSize ; i++ ) rRow[i] = SyncType( pRow[i] );
			}
			memcpy ( right + offset + Channels * _padSize * sizeof( SyncType ) , GetLabelRow(j+1) + _size - _padSize , _padSize * sizeof(LabelData<LabelType,Channels>) );
			if ( !rightStream->write( right , ( sizeof( SyncType )*Channels + sizeof(LabelData<LabelType,Channels>) ) * xSize ) )	exit(0);
		}
		FreeArray( right );
	}
}
template< int Channels , class PixelType , class LabelType , class StorageType , class SyncType >
void SocketedStreamingDivergence< Channels , PixelType , LabelType , StorageType , SyncType >::UnSetRestriction(void)
{
	if(parent)	parent->UnSetRestriction();
}

template< int Channels , class PixelType , class LabelType , class StorageType , class SyncType >
void SocketedStreamingDivergence< Channels , PixelType , LabelType , StorageType , SyncType >::SetRestriction(void)
{
	if(parent)	parent->SetRestriction();
}
template< int Channels , class PixelType , class LabelType , class StorageType , class SyncType >
void SocketedStreamingDivergence< Channels , PixelType , LabelType , StorageType , SyncType >::SetInteriorRestriction( Pointer( float ) lB , int c , int idx , int sRestrict , int eRestrict )
{
	int  off=FiniteElements1D<float,Type,Degree>::DotProduct<Type,Degree>::Helper::StartOffset();
	int dOff=FiniteElements1D<float,DERIVATIVE(Type),Degree-1>::DotProduct<Type,Degree>::Helper::StartOffset();
	int I =idx+ off;
	int dI=idx+dOff;
	Pointer( __m128 ) majorPtrs[2*Degree+1];
	Pointer( __m128 ) minorPtrs[2*Degree  ];

	int highStart = sRestrict-_start-Degree;
	int highEnd = eRestrict-_start+Degree;
	int highStart128 = (highStart-(RealPerWord-1))/RealPerWord;
	int highEnd128   = (highEnd  +(RealPerWord-1))/RealPerWord;
	int myStart = sRestrict-parent->start();
	int myEnd = eRestrict-parent->start();

	{
		for( int xx=0 ; xx<=2*Degree ; xx++ ) majorPtrs[xx] = ( Pointer( __m128 ) )GetDXRow(  I+xx , c );
		for( int xx=0 ; xx< 2*Degree ; xx++ ) minorPtrs[xx] = ( Pointer( __m128 ) )GetDYRow( dI+xx , c );
		__declspec (align(16)) float scratch[4];
		__m128 dot[Degree+1],dDot[Degree];
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
		for(int i=highStart128;i<highEnd128;i++)
		{
			localDMajorAccum[i]=                                                  _mm_mul_ps(dot[Degree],majorPtrs[Degree][i]);
			localDMajorAccum[i]=_mm_add_ps(localDMajorAccum[i],_mm_mul_ps( dot[0],_mm_add_ps(majorPtrs[0][i],majorPtrs[2*Degree  ][i])));
			localDMajorAccum[i]=_mm_add_ps(localDMajorAccum[i],_mm_mul_ps( dot[1],_mm_add_ps(majorPtrs[1][i],majorPtrs[2*Degree-1][i])));
		}
		for(int i=highStart128;i<highEnd128;i++)
		{
			localDMinorAccum[i]=                               _mm_mul_ps(dDot[0],_mm_sub_ps(minorPtrs[0][i],minorPtrs[2*Degree-1][i]));
			localDMinorAccum[i]=_mm_add_ps(localDMinorAccum[i],_mm_mul_ps(dDot[1],_mm_sub_ps(minorPtrs[1][i],minorPtrs[2*Degree-2][i])));
		}
		Pointer( float ) localDX = ( ( Pointer( float ) )localDMinorAccum)+FiniteElements1D<float,Type,Degree>::DotProduct<Type,Degree>::Helper::StartOffset();
		Pointer( float ) localDY = ( ( Pointer( float ) )localDMajorAccum)+FiniteElements1D<float,DERIVATIVE(Type),Degree-1>::DotProduct<Type,Degree>::Helper::StartOffset();

#if SUPPORT_CYLINDRICAL
		if( _periodicType!=SPHERICAL_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
		if(!_spherical)
#endif // SUPPORT_CYLINDRICAL
		{
			if(myStart+_start<Degree)
			{
				for(int j=myStart;j<Degree-_start;j++)
				{
					double temp=0;
					const float* xValues =  dotMajorStencil.caseTable[j].values;
					const float* yValues = dDotMajorStencil.caseTable[j].values;

					// Partial w.r.t minor index
					int J = off+j;
					for(int dj=0;dj<=2*Degree;dj++)
					{
						if( J+dj+_start<0 ||  J+dj+_start>=major)	continue;
						temp+=localDX[j+dj] * xValues[dj];
					}
					// Partial w.r.t major index
					int DJ = dOff + j;
					for(int dj=0;dj<2*Degree;dj++)
					{
						if(DJ+dj+_start<0 || DJ+dj+_start>=major-1)	continue;
						temp+=localDY[j+dj] * yValues[dj];
					}
					lB[j]=temp;
					meanB[c]+=temp;
				}
				myStart = Degree-_start;
			}
			if(myEnd+_start>major-Degree)
			{
				for(int j=major-Degree-_start;j<myEnd;j++)
				{
					double temp=0;
					const float* xValues =  dotMajorStencil.caseTable[2*Degree+(j+_start-(major-1))].values;
					const float* yValues = dDotMajorStencil.caseTable[2*Degree+(j+_start-(major-1))].values;

					// Partial w.r.t minor index
					int J = off+j;
					for(int dj=0;dj<=2*Degree;dj++)
					{
						if( J+dj+_start<0 ||  J+dj+_start>=major)	continue;
						temp+=localDX[j+dj] * xValues[dj];
					}
					// Partial w.r.t major index
					int DJ = dOff + j;
					for(int dj=0;dj<2*Degree;dj++)
					{
						if(DJ+dj+_start<0 || DJ+dj+_start>=major-1)	continue;
						temp+=localDY[j+dj] * yValues[dj];
					}
					lB[j]=temp;
					meanB[c]+=temp;
				}
				myEnd = major-Degree-_start;
			}
		}
		const float* xValues =  dotMajorStencil.caseTable[Degree].values;
		const float* yValues = dDotMajorStencil.caseTable[Degree].values;
		for(int j=myStart;j<myEnd;j++)
		{
			double temp=0;
			// Partial w.r.t minor index
			temp = localDX[j+Degree]*xValues[Degree];
			for(int dj=0;dj<Degree;dj++)	temp+=(localDX[j+dj]+localDX[j+2*Degree-dj])*xValues[dj];

			// Partial w.r.t major index
			for(int dj=0;dj<Degree;dj++)	temp+=(localDY[j+dj]-localDY[j+2*Degree-1-dj])*yValues[dj];

			lB[j]=temp;
			meanB[c]+=temp;
		}
	}
}
template< int Channels , class PixelType , class LabelType , class StorageType , class SyncType >
void SocketedStreamingDivergence< Channels , PixelType , LabelType , StorageType , SyncType >::SetRestriction( Pointer( float ) lB , int c , int idx , int sRestrict , int eRestrict )
{
#if DEBUG_FLAGS
	if(debugFlags.noDivergence)	return;
#endif // DEBUG_FLAGS
#if SUPPORT_CYLINDRICAL
	if( (idx>Degree && idx<minor-Degree) || _periodicType==SPHERICAL_PERIODIC )	 return SetInteriorRestriction(lB,c,idx,sRestrict,eRestrict);
#else // !SUPPORT_CYLINDRICAL
	if( (idx>Degree && idx<minor-Degree) || _spherical)	 return SetInteriorRestriction(lB,c,idx,sRestrict,eRestrict);
#endif // SUPPORT_CYLINDRICAL
	int ii;
#if SUPPORT_CYLINDRICAL
	if( _periodicType==SPHERICAL_PERIODIC )	ii=Degree;
#else // !SUPPORT_CYLINDRICAL
	if(_spherical)					ii=Degree;
#endif // SUPPORT_CYLINDRICAL
	else
		if(idx<Degree)				ii=idx;
		else if(idx>=minor-Degree)	ii=2*Degree+(idx-(minor-1));
		else						ii=Degree;
	int dI=idx+FiniteElements1D<float,DERIVATIVE(Type),Degree-1>::DotProduct<Type,Degree>::Helper::StartOffset();
	int I =idx+FiniteElements1D<float,Type,Degree>::DotProduct<Type,Degree>::Helper::StartOffset();
	int myStart = sRestrict-parent->start();
	int myEnd = eRestrict-parent->start();

	for(int j=myStart;j<myEnd;j++)
	{
		double temp=0;
		int jj;
#if SUPPORT_CYLINDRICAL
		if( _periodicType!=NO_PERIODIC )	jj=Degree;
#else // !SUPPORT_CYLINDRICAL
		if(_spherical)						jj=Degree;
#endif // SUPPORT_CYLINDRICAL
		else
			if(j+_start<Degree)				jj=j+_start;
			else if(j+_start>=major-Degree)	jj=2*Degree+(j+_start-(major-1));
			else							jj=Degree;
		int dJ=j+FiniteElements1D<float,DERIVATIVE(Type),Degree-1>::DotProduct<Type,Degree>::Helper::StartOffset();
		int J =j+FiniteElements1D<float,Type,Degree>::DotProduct<Type,Degree>::Helper::StartOffset();

		// Partial w.r.t minor index
		for(int di=0;di<2*Degree;di++)
		{
#if SUPPORT_CYLINDRICAL
			if( (dI+di<0 || dI+di>=minor-1) && _periodicType!=SPHERICAL_PERIODIC )	continue;
#else // !SUPPORT_CYLINDRICAL
			if( (dI+di<0 || dI+di>=minor-1) && !_spherical)	continue;
#endif // SUPPORT_CYLINDRICAL
			Pointer( float ) localD = GetDYRow( dI+di , c ) + J;
			for(int dj=0;dj<=2*Degree;dj++)
			{
#if SUPPORT_CYLINDRICAL
				if( (J+dj+_start<0 || J+dj+_start>=major) && _periodicType==NO_PERIODIC )	continue;
#else // !SUPPORT_CYLINDRICAL
				if( (J+dj+_start<0 || J+dj+_start>=major) && !_spherical)	continue;
#endif // SUPPORT_CYLINDRICAL
				temp+=localD[dj]*divergenceStencil.caseTable[jj][ii].values2[dj][di];
			}
		}
		// Partial w.r.t major index
		for(int di=0;di<=2*Degree;di++)
		{
#if SUPPORT_CYLINDRICAL
			if( (I+di<0 || I+di>=minor) && _periodicType!=SPHERICAL_PERIODIC )	continue;
#else // !SUPPORT_CYLINDRICAL
			if( (I+di<0 || I+di>=minor) && !_spherical)	continue;
#endif // SUPPORT_CYLINDRICAL
			Pointer( float ) localD = GetDXRow( I+di , c ) + dJ;
			for(int dj=0;dj<2*Degree;dj++)
			{
#if SUPPORT_CYLINDRICAL
				if( (dJ+dj+_start<0 || dJ+dj+_start>=major-1) && _periodicType==NO_PERIODIC )	continue;
#else // !SUPPORT_CYLINDRICAL
				if( (dJ+dj+_start<0 || dJ+dj+_start>=major-1) && !_spherical)	continue;
#endif // SUPPORT_CYLINDRICAL
				temp+=localD[dj]*divergenceStencil.caseTable[jj][ii].values1[dj][di];
			}
		}
		lB[j] = temp;
		meanB[c]+=temp;
	}
}
template< int Channels , class PixelType , class LabelType , class StorageType , class SyncType >
void SocketedStreamingDivergence< Channels , PixelType , LabelType , StorageType , SyncType >::_setPartials( int y )
{
	if( y&1 )
	{
		SyncImageLeft( y , false );
		_setPartialsX( y , -_padSize , _size-_padSize );
		_setPartialsY( y , -_padSize , _size-_padSize );
		SyncImageRight( y , true );
		_setPartialsX( y , _size-_padSize-1 , _size+_padSize );
		_setPartialsY( y , _size-_padSize , _size+_padSize );
	}
	else
	{
		SyncImageRight( y , false );
		_setPartialsX( y , _padSize , _size+_padSize );
		_setPartialsY( y , _padSize , _size+_padSize );
		SyncImageLeft( y , true );
		_setPartialsX( y , -_padSize , _padSize+1 );
		_setPartialsY( y , -_padSize , _padSize );
	}
}
template< int Channels , class PixelType , class LabelType , class StorageType , class SyncType >
void SocketedStreamingDivergence< Channels , PixelType , LabelType , StorageType , SyncType >::_setPartialsX( int y , int start , int end )
{
#if DEBUG_FLAGS
	if(debugFlags.noGradient)	return;
#endif // DEBUG_FLAGS
	Pointer( float ) pixelRow[Channels];
	Pointer( float ) dx[Channels];
	Pointer( LabelData< LabelType , Channels > ) labelRow = GetLabelRow( y );
	for(int c=0;c<Channels;c++)
	{
		pixelRow[c] = GetPixelRow( y , c );
		dx[c] = GetDXRow( y , c );
	}
	if( blackOut )
	{
		LabelData< LabelType , Channels > black = LabelData< LabelType , Channels >::BlackLabel( );
		for( int x=start ; x<end-1 ; x++ )
		{
			bool useGradient=(labelRow[x]==labelRow[x+1]);
			bool b0 = ( labelRow[x  ] == black );
			bool b1 = ( labelRow[x+1] == black );
#if 1
			if( b0 || b1 || useGradient ) for( int c=0 ; c<Channels ; c++ ) dx[c][x] = pixelRow[c][x+1]-pixelRow[c][x];
			else                          for( int c=0 ; c<Channels ; c++ ) dx[c][x] = 0;
#else
			if( b0 && b1 )			for( int c=0 ; c<Channels ; c++ ) dx[c][x] = 0;
			else if( b0 )			for( int c=0 ; c<Channels ; c++ ) dx[c][x] = pixelRow[c][x+1];
			else if( b1 )			for( int c=0 ; c<Channels ; c++ ) dx[c][x] =                 -pixelRow[c][x];
			else if( useGradient )	for( int c=0 ; c<Channels ; c++ ) dx[c][x] = pixelRow[c][x+1]-pixelRow[c][x];
			else					for( int c=0 ; c<Channels ; c++ ) dx[c][x] = 0;
#endif
		}
	}
	else
		for( int x=start ; x<end-1 ; x++ )
		{
			bool useGradient=(labelRow[x]==labelRow[x+1]);
			if( useGradient ) for( int c=0 ; c<Channels ; c++ ) dx[c][x] = pixelRow[c][x+1]-pixelRow[c][x];
			else              for( int c=0 ; c<Channels ; c++ ) dx[c][x] = 0;
		}
}
template< int Channels , class PixelType , class LabelType , class StorageType , class SyncType >
void SocketedStreamingDivergence< Channels , PixelType , LabelType , StorageType , SyncType >::_setPartialsY( int y , int start , int end )
{
#if DEBUG_FLAGS
	if(debugFlags.noGradient)	return;
#endif // DEBUG_FLAGS
	Pointer( float ) dy[Channels];
	Pointer( float ) oldPixelRow[Channels];
	Pointer( float ) newPixelRow[Channels];
	Pointer( LabelData< LabelType , Channels > ) oldLabelRow = GetLabelRow( y-1 );
	Pointer( LabelData< LabelType , Channels > ) newLabelRow = GetLabelRow( y   );
	for( int c=0 ; c<Channels ; c++ )
	{
		oldPixelRow[c] = GetPixelRow( y-1 , c );
		newPixelRow[c] = GetPixelRow( y   , c );
		dy[c] = GetDYRow( y-1 , c );
	}

	if( blackOut )
	{
		LabelData< LabelType , Channels > black = LabelData< LabelType , Channels >::BlackLabel( );
		for( int x=start ; x<end ; x++ )
		{
			bool useGradient=(newLabelRow[x]==oldLabelRow[x]);
			bool b0 = ( oldLabelRow[x]==black );
			bool b1 = ( newLabelRow[x]==black );
#if 1
			if( b0 || b1 || useGradient ) for( int c=0 ; c<Channels ; c++ ) dy[c][x] = newPixelRow[c][x]-oldPixelRow[c][x];
			else                          for( int c=0 ; c<Channels ; c++ ) dy[c][x] = 0;
#else
			if( b0 && b1 )			for( int c=0 ; c<Channels ; c++ ) dy[c][x] = 0;
			else if( b0 )			for( int c=0 ; c<Channels ; c++ ) dy[c][x] = newPixelRow[c][x];
			else if( b1 )			for( int c=0 ; c<Channels ; c++ ) dy[c][x] =                  -oldPixelRow[c][x];
			else if( useGradient )	for( int c=0 ; c<Channels ; c++ ) dy[c][x] = newPixelRow[c][x]-oldPixelRow[c][x];
			else					for( int c=0 ; c<Channels ; c++ ) dy[c][x] = 0;
#endif
		}
	}
	else
		for( int x=start ; x<end ; x++ )
		{
			bool useGradient=(newLabelRow[x]==oldLabelRow[x]);
#if 0
#if SUPPORT_CYLINDRICAL
			if( _periodicType==SPHERICAL_PERIODIC && y==0 )     if( (x>=major/4 && x<major/2) || (x>=3*major/4 && x<major ) ) useGradient = false;
			if( _periodicType==SPHERICAL_PERIODIC && y==minor ) if( (x>=major/2 && x<major) ) useGradient = false;
#else // !SUPPORT_CYLINDRICAL
			if( _spherical && y==0 )     if( (x>=major/4 && x<major/2) || (x>=3*major/4 && x<major ) ) useGradient = false;
			if( _spherical && y==minor ) if( (x>=major/2 && x<major) ) useGradient = false;
#endif // SUPPORT_CYLINDRICAL
#endif
			if( useGradient ) for( int c=0 ; c<Channels ; c++ ) dy[c][x] = newPixelRow[c][x]-oldPixelRow[c][x];
			else              for( int c=0 ; c<Channels ; c++ ) dy[c][x] = 0;
		}
}

template< int Channels , class PixelType , class LabelType , class StorageType , class SyncType >
void SocketedStreamingDivergence< Channels , PixelType , LabelType , StorageType , SyncType >::_setPartialsX( int y )
{
#if DEBUG_FLAGS
	if(debugFlags.noGradient)	return;
#endif // DEBUG_FLAGS
	float *dx[Channels];
	float *pixelRow[Channels];
	LabelData<LabelType,Channels>* labelRow = GetLabelRow(y);
	for(int c=0;c<Channels;c++)
	{
		pixelRow[c] = GetPixelRow(y,c);
		dx[c] = GetDXRow(y,c);
	}
	for(int x=-_padSize;x<_size+_padSize-1;x++)
	{
		bool useGradient=(labelRow[x]==labelRow[x+1]);
		if(useGradient)		for(int c=0;c<Channels;c++)	dx[c][x] = pixelRow[c][x+1]-pixelRow[c][x];
		else				for(int c=0;c<Channels;c++)	dx[c][x] = 0;
	}
}
template< int Channels , class PixelType , class LabelType , class StorageType , class SyncType >
void SocketedStreamingDivergence< Channels , PixelType , LabelType , StorageType , SyncType >::_setPartialsY( int y )
{
#if DEBUG_FLAGS
	if(debugFlags.noGradient)	return;
#endif // DEBUG_FLAGS
	float *dy[Channels];
	float *oldPixelRow[Channels],*newPixelRow[Channels];
	LabelData<LabelType,Channels> *oldLabelRow = GetLabelRow(y-1);
	LabelData<LabelType,Channels> *newLabelRow = GetLabelRow(y  );

	for(int c=0;c<Channels;c++)
	{
		oldPixelRow[c] = GetPixelRow(y-1,c);
		newPixelRow[c] = GetPixelRow(y  ,c);
		dy[c] = GetDYRow( y-1 , c );
	}
	for(int x=-_padSize;x<_size+_padSize;x++)
	{
		bool useGradient=(newLabelRow[x]==oldLabelRow[x]);
		if(useGradient)		for(int c=0;c<Channels;c++)	dy[c][x] = newPixelRow[c][x]-oldPixelRow[c][x];
		else				for(int c=0;c<Channels;c++)	dy[c][x] = 0;
	}
}
template< int Channels , class PixelType , class LabelType , class StorageType , class SyncType >
bool SocketedStreamingDivergence< Channels , PixelType , LabelType , StorageType , SyncType >::IterateRestriction(void)
{
	if( index>=minor ) return false;
	int idx = index+Degree+1;

	// Load the next row of pixel/label values
	if(idx<_iHeight)
	{
		{
			{
				Pointer( LabelType ) labels = ( Pointer( LabelType ) )(*labelStream)[idx];
				Pointer( LabelData< LabelType , Channels > ) labelRow = GetLabelRow(idx);
				if(labelStream->SeparateColors() )
					for( int c = 0 ; c < Channels ; c++ )
						for ( int x = 0 ; x < _iWidth ; x++) labelRow[x].l[c] = labels[_iWidth*c+x];
				else
					for ( int x = 0 ; x < _iWidth ; x++ )
						for ( int c = 0 ; c < Channels ; c++ ) labelRow[x].l[c] = labels[Channels*x+c];
				for ( int x = _iWidth ; x<_size ; x++ )	labelRow[x] = labelRow[_iWidth-1];
			}
			labelStream->advance();
			{
				Pointer( PixelType ) pixels = ( Pointer( PixelType ) )(*pixelStream)[idx];
				for( int c=0 ; c<Channels ; c++ )
				{
					Pointer( float ) pixelRow = GetPixelRow( idx , c );
					if( pixelStream->SeparateColors() )
						for ( int x=0 ; x<_iWidth ; x++ )
							pixelRow[x] = float( pixels[ _iWidth*c + x ] );
					else
						for ( int x=0 ; x<_iWidth ; x++ ) pixelRow[x] = float( pixels[ x*Channels + c ] );
#if PAD_OUT
					for( int x=_iWidth ; x<_size ; x++ ) pixelRow[x] = pixelRow[_iWidth-1];
#endif // PAD_OUT
					// Note: We should probably set the rest of the pixel rows to something consistent
					// so that we don't get unpredictable behavior when we're using padding on a spherical image
					// (not that it is reasonable to pad in the spherical case)
					if( blackOut )
					{
						LabelData< LabelType , Channels > black = LabelData< LabelType , Channels >::BlackLabel( );
						Pointer( LabelData< LabelType , Channels > ) labelRow = GetLabelRow( idx );
						for( int x=0 ; x<_iWidth ; x++ )
						{
							if( labelRow[x]==black ) pixelRow[x] = 0.f;
							average[c] += pixelRow[x];
						}
					}
					else for ( int x = 0 ; x < _iWidth ; x++) average[c] += pixelRow[x];
				}
			}
			pixelStream->advance();
		}
	}
	else if ( idx<minor )
	{
		Pointer( LabelData< LabelType , Channels > ) oldLabels = GetLabelRow( idx-1 );
		Pointer( LabelData< LabelType , Channels > ) newLabels = GetLabelRow( idx   );
		for( int x = 0 ; x < _size ; x++ )	newLabels[x] = oldLabels[x];
#if PAD_OUT
		for( int c=0 ; c<Channels ; c++ )
		{
			Pointer( float ) oldPixels = GetPixelRow( idx-1 , c );
			Pointer( float ) newPixels = GetPixelRow( idx   , c );
			memcpy( newPixels , oldPixels , sizeof( float ) * _size );
		}
#else // !PAD_OUT
		for( int c=0 ; c<Channels ; c++ ) memset( GetPixelRow( idx , c ) , 0 , sizeof(float)*_size );
#endif // PAD_OUT
	}
#if SUPPORT_CYLINDRICAL
	if( _periodicType==SPHERICAL_PERIODIC && idx==Degree-1 ) for( int d=0 ; d<Degree ; d++ ) SyncImageHead( d , false ) , SyncImageHead( d , true );
	if( _periodicType==SPHERICAL_PERIODIC && idx== minor-1 ) for( int d=minor-Degree ; d<minor ; d++ ) SyncImageTail( d , false ) , SyncImageTail( d , true );
#else // !SUPPORT_CYLINDRICAL
	if( _spherical && idx==Degree-1 ) for( int d=0 ; d<Degree ; d++ ) SyncImageHead( d , false ) , SyncImageHead( d , true );
	if( _spherical && idx== minor-1 ) for( int d=minor-Degree ; d<minor ; d++ ) SyncImageTail( d , false ) , SyncImageTail( d , true );
#endif // SUPPORT_CYLINDRICAL
	// If this is the first row, perform the preliminary left/right synchronization
#if SUPPORT_CYLINDRICAL
	if( _periodicType==SPHERICAL_PERIODIC && idx==Degree-1 )
#else // !SUPPORT_CYLINDRICAL
	if( _spherical && idx==Degree-1 )
#endif // SUPPORT_CYLINDRICAL
	{
		SyncImageLeft( -Degree , false ) , SyncImageRight( -Degree , false );
		SyncImageLeft( -Degree , true  ) , SyncImageRight( -Degree , true  );
	}
#if SUPPORT_CYLINDRICAL
	if( _periodicType!=SPHERICAL_PERIODIC && idx==0 )
#else // !SUPPORT_CYLINDRICAL
	if( !_spherical && idx==0 )
#endif // SUPPORT_CYLINDRICAL
	{
		SyncImageLeft( 0 , false ) , SyncImageRight( 0 , false );
		SyncImageLeft( 0 , true  ) , SyncImageRight( 0 , true  );
	}
	idx = index+Degree;
	if( idx>=0 )
	{
#if SUPPORT_CYLINDRICAL
		if( _periodicType==SPHERICAL_PERIODIC && idx==Degree-1)
#else // !SUPPORT_CYLINDRICAL
		if(_spherical && idx==Degree-1)
#endif // SUPPORT_CYLINDRICAL
		{
			for(int y = -Degree ; y <= 0 ; y++) _setPartials(y);
		}
#if SUPPORT_CYLINDRICAL
		if( _periodicType==SPHERICAL_PERIODIC && idx==minor-1 )
#else // !SUPPORT_CYLINDRICAL
		if( _spherical && idx==minor-1 )
#endif // SUPPORT_CYLINDRICAL
		{
			for( int y=minor; y < minor + Degree ; y++ ) _setPartials( y );
		}
		_setPartials( idx );
		if( parent ) parent->IterateRestriction();
	}
	index++;
	return true;
}
template< int Channels , class PixelType , class LabelType , class StorageType , class SyncType >
void SocketedStreamingDivergence< Channels , PixelType , LabelType , StorageType , SyncType >::SolveRestriction(void)
{
	// Run to completion...
	while( IterateRestriction() ){;}
	// ...and finish up the trailing parent
	if(parent)	parent->SolveRestriction();
}
////////////////////////
// StreamingLaplacian //
////////////////////////
template< int Channels , class PixelType , class StorageType , class SyncType >
SocketedStreamingLaplacian< Channels , PixelType , StorageType , SyncType >::SocketedStreamingLaplacian(void)
{
#if TIME_IO
	vSync = hSync = 0;
#endif // TIME_IO
	parent = NULL;
	for( int c=0 ; c<Channels ; c++ ) average[c] = meanB[c] = 0;
	syncBuffer = NullPointer< SyncType >( );
	for( int i=0 ; i<Degree ; i++ ) localPAccum[i] = NullPointer< __m128 >( );
	for( int i=0 ; i<ISize*Channels ; i++ ) PStream[i] = NullPointer< __m128 >( );
	lapTemplates = AllocArray< TemplateSSE >( 3 * ( 2*Degree+1 ) , ALIGNMENT  , "SocketedStreamingLaplacian::lapTemplates" );
}
template< int Channels , class PixelType , class StorageType , class SyncType >
SocketedStreamingLaplacian< Channels , PixelType , StorageType , SyncType >::~SocketedStreamingLaplacian(void)
{
	FreeArray( syncBuffer );
	for( int i=0 ; i<Degree ; i++ )
		if( localPAccum[i] )
		{
			localPAccum[i] -= _padSize128;
			FreeArray( localPAccum[i] );
		}
	for( int i=0 ; i<ISize*Channels ; i++ ) FreeArray( PStream[i] );
	FreeArray( lapTemplates );
}

template< int Channels , class PixelType , class StorageType , class SyncType >
void SocketedStreamingLaplacian< Channels , PixelType , StorageType , SyncType >::Initialize(StreamingGrid *pixels,double iWeight,double gWeight,
																			int start,int end,int major,int minor,int iters,
																			DataStream* leftStream,SOCKET* syncSockets,DataStream* rightStream,
#if SUPPORT_CYLINDRICAL
																			bool memoryMappedFile , int periodicType , MultiStreamIOServer* server
#else // !SUPPORT_CYLINDRICAL
																			bool memoryMappedFile,bool spherical , MultiStreamIOServer* server
#endif // SUPPORT_CYLINDRICAL
																			)
{
	if( syncSockets )
	{
#if SUPPORT_CYLINDRICAL
		if( !SetSocketedStreamData( major , minor , start , end , leftStream , syncSockets[0] , rightStream , periodicType ) ) exit(0);
#else // !SUPPORT_CYLINDRICAL
		if( !SetSocketedStreamData( major , minor , start , end , leftStream , syncSockets[0] , rightStream , spherical ) ) exit(0);
#endif // SUPPORT_CYLINDRICAL
	}
	else
	{
#if SUPPORT_CYLINDRICAL
		if( !SetSocketedStreamData( major , minor , start , end , leftStream , INVALID_SOCKET , rightStream , periodicType ) ) exit(0);
#else // !SUPPORT_CYLINDRICAL
		if( !SetSocketedStreamData( major , minor , start , end , leftStream , INVALID_SOCKET , rightStream , spherical ) ) exit(0);
#endif // SUPPORT_CYLINDRICAL
	}
#if SUPPORT_CYLINDRICAL
	this->_periodicType=periodicType;
#else // !SUPPORT_CYLINDRICAL
	this->_spherical=spherical;
#endif // SUPPORT_CYLINDRICAL
	this->major=major;
	this->minor=minor;
	pixelStream=pixels;

	if(parent)
		if( syncSockets )
			parent->Initialize( iWeight , start , end , major , minor , iters ,
			leftStream , syncSockets+1 , rightStream ,
#if SUPPORT_CYLINDRICAL
			memoryMappedFile , periodicType , server
#else // !SUPPORT_CYLINDRICAL
			memoryMappedFile , spherical , server
#endif // SUPPORT_CYLINDRICAL
			);
		else
			parent->Initialize( iWeight , start , end , major , minor , iters ,
			leftStream , NULL , rightStream ,
#if SUPPORT_CYLINDRICAL
			memoryMappedFile , periodicType , server
#else // !SUPPORT_CYLINDRICAL
			memoryMappedFile , spherical , server
#endif // SUPPORT_CYLINDRICAL
			);
	syncBuffer = AllocArray< SyncType >( _size * Channels * Degree , 1 , "SocketedStreamingLaplacian::syncBuffer" );
	for( int i=0 ; i<Degree ; i++ )
	{
		localPAccum[i]  = AllocArray< __m128 >( _paddedSize128 , ALIGNMENT , "SocketedStreamingLaplacian::localPAccum" );
		localPAccum[i] += _padSize128;
	}
	MatrixStencil lStencil , dStencil , stencil;
	FiniteElements2D<double,Type,Degree>::LaplacianStencil(major,minor,lStencil,true);
	FiniteElements2D<double,Type,Degree>::DotProductStencil(major,minor,dStencil);
	for(int i=0;i<=2*Degree;i++)
		for(int j=0;j<=2*Degree;j++)
			for(int k=0;k<=2*Degree;k++)
				for(int l=0;l<=2*Degree;l++)
					stencil.caseTable[i][j].values[k][l] = lStencil.caseTable[i][j].values[k][l]*gWeight + dStencil.caseTable[i][j].values[k][l]*iWeight;
//					stencil.caseTable[i][j].values[k][l] = lStencil.caseTable[i][j].values[k][l]*gWeight;


	__declspec (align(16)) float scratch[4];
	for( int i=0 ; i<=2*Degree ; i++ )	// Iterate over the minor index in the mask
		for( int j=0 ; j<3 ; j++ )
			for( int k=0 ; k<=2*Degree ; k++ )
			{
				int jj;

				if(j==0)	jj=0;
				else		jj=Degree;
				scratch[0]=stencil.caseTable[i][jj].values[k][2];
				scratch[1]=stencil.caseTable[i][jj].values[k][3];
				scratch[2]=stencil.caseTable[i][jj].values[k][4];
				if(j!=2)	scratch[3]=stencil.caseTable[i][Degree].values[k][1];
				else		scratch[3]=0;
				lapTemplates[3*i+j].matrixValues[0][k] = _mm_load_ps(scratch);

				if(j==0)	jj=1;
				else		jj=Degree;
				for(int l=0;l<4;l++)	scratch[l]=stencil.caseTable[i][jj].values[k][l+1];
				lapTemplates[3*i+j].matrixValues[1][k]=_mm_load_ps(scratch);

				if(j==2)	jj=2*Degree-1;
				else		jj=Degree;
				for(int l=0;l<4;l++)	scratch[l]=stencil.caseTable[i][jj].values[k][l];
				lapTemplates[3*i+j].matrixValues[2][k]=_mm_load_ps(scratch);

				if(j==2)	jj=2*Degree;
				else		jj=Degree;
				if(j!=0)	scratch[0]=stencil.caseTable[i][Degree].values[k][3];
				else		scratch[0]=0;
				scratch[1]=stencil.caseTable[i][jj].values[k][0];
				scratch[2]=stencil.caseTable[i][jj].values[k][1];
				scratch[3]=stencil.caseTable[i][jj].values[k][2];
				lapTemplates[3*i+j].matrixValues[3][k]=_mm_load_ps(scratch);
			}
}
template< int Channels , class PixelType , class StorageType , class SyncType >
void SocketedStreamingLaplacian< Channels , PixelType , StorageType , SyncType >::InitRestriction(void)
{
	index = -Degree-1;
	_iHeight=pixelStream->rows();
	_iWidth=pixelStream->rowSize()/(Channels*sizeof(PixelType));
	if(pixelStream->rowSize()!=_iWidth*Channels*sizeof(PixelType))
		fprintf(stderr,"Pixel width failure: %d != %d\n",pixelStream->rowSize(),_iWidth*Channels*sizeof(PixelType)) ,	exit(0);

	pixelStream->reset(true,1);
	FiniteElements1D<float,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(major,dotMajorStencil,0,0);
	FiniteElements1D<float,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(minor,dotMinorStencil,0,0);
	FiniteElements1D<float,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(major,d2DotMajorStencil,1,1,false);
	FiniteElements1D<float,Type,Degree>::DotProduct<Type,Degree>::DotProductStencil(minor,d2DotMinorStencil,1,1,false);
	if(parent)	parent->InitRestriction();
}

template< int Channels , class PixelType , class StorageType , class SyncType >
Pointer( float ) SocketedStreamingLaplacian< Channels , PixelType , StorageType , SyncType >::GetPixelRow(int row,int channel)
{
//	return ( Pointer( float ) )( PStream[MyModIndex( row*Channels + channel , ISize*Channels )] );
	return ( Pointer( float ) )( PStream[MyModIndex( row*Channels + channel , ISize*Channels )] ) + _padSize;
}
template< int Channels , class PixelType , class StorageType , class SyncType >
void SocketedStreamingLaplacian< Channels , PixelType , StorageType , SyncType >::SyncImageHead( int idx , bool read )
{
#if DEBUG_FLAGS
	if( debugFlags.noTopBottom ) return;
#endif // DEBUG_FLAGS
	if( idx<0 || idx>=Degree )
	{
		fprintf( stderr , "Synchronizing laplacian in non-head row: %d\n" , idx );
		return;
	}
	int r = -idx-1;
	if(syncXSocket!=INVALID_SOCKET)
	{
		int size=_size * Channels;
		if( read )
		{
#if TIME_IO
			double t = Time();
#endif // TIME_IO
			ReceiveOnSocket( syncXSocket , syncBuffer , sizeof( SyncType ) * size , "SocketedStreamingLaplacian< Channels , PixelType , StorageType , SyncType >::SyncImageHead" );
#if TIME_IO
			vSync += Time() - t;
#endif // TIME_IO
			for(int c=0;c<Channels;c++)
			{
				Pointer( float ) pixelRow = GetPixelRow( r , c );
				for( int x=0 ; x<_size ; x++ ) pixelRow[x] = float( syncBuffer[_size*c+x] );
			}
		}
		else
		{
			for(int c=0;c<Channels;c++)
			{
				Pointer( float ) pixelRow = GetPixelRow( idx , c );
				for( int x=0 ; x<_size ; x++ ) syncBuffer[_size*c+x] = SyncType( pixelRow[x] );
			}
			SendOnSocket(syncXSocket,syncBuffer, sizeof( SyncType ) * size , "SocketedStreamingLaplacian< Channels , PixelType , StorageType , SyncType >::SyncImageHead");
		}
	}
}
template< int Channels , class PixelType , class StorageType , class SyncType >
void SocketedStreamingLaplacian< Channels , PixelType , StorageType , SyncType >::SyncImageTail(int idx , bool read )
{
#if DEBUG_FLAGS
	if( debugFlags.noTopBottom ) return;
#endif // DEBUG_FLAGS
	if( idx<minor-Degree || idx>=minor )
	{
		fprintf( stderr , "Synchronizing laplacian in non-tail row: %d\n" , idx );
		return;
	}
	int r = 2 * minor - 1 - idx;
	if( syncXSocket!=INVALID_SOCKET )
	{
		int size = _size * Channels;
		if(read)
		{
#if TIME_IO
			double t = Time();
#endif // TIME_IO
			ReceiveOnSocket( syncXSocket, syncBuffer , sizeof( SyncType ) * size , "SocketedStreamingLaplacian< Channels , PixelType , StorageType , SyncType >::SyncImageTail" );
#if TIME_IO
			vSync += Time() - t;
#endif // TIME_IO
			for( int c=0 ; c<Channels ; c++ )
			{
				Pointer( float ) pixelRow = GetPixelRow( r , c );
				for( int x=0 ; x<_size ; x++ ) pixelRow[x] = float( syncBuffer[_size*c+x] );
			}
		}
		else
		{
			for( int c=0 ; c<Channels ; c++ )
			{
				Pointer( float ) pixelRow = GetPixelRow( idx , c );
				for( int x=0 ; x<_size ; x++ ) syncBuffer[_size*c+x] = SyncType( pixelRow[x] );
			}
			SendOnSocket(syncXSocket,syncBuffer, sizeof( SyncType ) * size , "SocketedStreamingLaplacian< Channels , PixelType , StorageType , SyncType >::SyncImageTail" );
		}
	}
}

template< int Channels , class PixelType , class StorageType , class SyncType >
void SocketedStreamingLaplacian< Channels , PixelType , StorageType , SyncType >::SyncImageLeft(int j,bool read)
{
#if DEBUG_FLAGS
	if( debugFlags.noLeftRight ) return;
#endif // DEBUG_FLAGS
	if( leftStream )
	{
		int xSize = 2*_padSize;
		int offset = sizeof( SyncType ) * Channels * _padSize;
		Pointer( byte ) left = AllocArray< byte >( sizeof( SyncType )*Channels * xSize , 1 , "SocketedStreamingLaplacian::SyncImageLeft (left)" );
		if(read)
		{
#if TIME_IO
			double t = Time();
#endif // TIME_IO
			if ( !leftStream->read( ( Pointer( byte ) )left , sizeof( SyncType ) * Channels * xSize ) )	exit(0);
#if TIME_IO
			hSync += Time() - t;
#endif // TIME_IO
			for( int c=0 ; c<Channels ; c++ )
			{
				Pointer( float ) pRow = GetPixelRow( j , c ) - _padSize;
				Pointer( SyncType ) lRow = ( Pointer( SyncType ) )( left + c * _padSize * sizeof( SyncType ) );
				for( int i=0 ; i<_padSize ; i++ ) pRow[i] = float( lRow[i] );
			}
			for( int c=0 ; c<Channels ; c++ )
			{
				Pointer( float ) pRow = GetPixelRow( j +1 , c ) - _padSize;
				Pointer( SyncType ) lRow = ( Pointer( SyncType ) )( left + offset + c * _padSize * sizeof( SyncType ) );
				for( int i=0 ; i<_padSize ; i++ ) pRow[i] = float( lRow[i] );
			}
		}
		else
		{
			for( int c=0 ; c<Channels ; c++ )
			{
				Pointer( float ) pRow = GetPixelRow( j , c );
				Pointer( SyncType ) lRow = ( Pointer( SyncType ) )( left + c * _padSize * sizeof( SyncType ) );
				for( int i=0 ; i<_padSize ; i++ ) lRow[i] = SyncType( pRow[i] );
			}
			for( int c=0 ; c<Channels ; c++ )
			{
				Pointer( float ) pRow = GetPixelRow( j + 1 , c );
				Pointer( SyncType ) lRow = ( Pointer( SyncType ) )( left + offset + c * _padSize * sizeof( SyncType ) );
				for( int i=0 ; i<_padSize ; i++ ) lRow[i] = SyncType( pRow[i] );
			}
			if ( !leftStream->write( ( Pointer( byte ) )left , sizeof( SyncType )*Channels * xSize ) ) exit(0);
		}
		FreeArray(left);
	}
}
template< int Channels , class PixelType , class StorageType , class SyncType >
void SocketedStreamingLaplacian< Channels , PixelType , StorageType , SyncType >::SyncImageRight(int j,bool read)
{
#if DEBUG_FLAGS
	if( debugFlags.noLeftRight ) return;
#endif // DEBUG_FLAGS
	if( rightStream )
	{
		int xSize = 2*_padSize;
		int offset = sizeof(float) * Channels * _padSize;
		Pointer( byte ) right = AllocArray< byte >( sizeof( SyncType ) * Channels * xSize , 1 , "SocketedStreamingLaplacian::SyncImageRight (right)" );
		if( read )
		{
#if TIME_IO
			double t = Time();
#endif // TIME_IO
			if ( !rightStream->read( ( Pointer( byte ) )right , sizeof( SyncType ) * Channels * xSize ) ) exit(0);
#if TIME_IO
			hSync += Time() - t;
#endif // TIME_IO
			for( int c=0 ; c<Channels ; c++ )
			{
				Pointer( float ) pRow = GetPixelRow( j , c ) + _size;
				Pointer( SyncType ) rRow = ( Pointer( SyncType ) )( right + c * _padSize * sizeof( SyncType ) );
				for( int i=0 ; i<_padSize ; i++ ) pRow[i] = float( rRow[i] );
			}
			for( int c=0 ; c<Channels ; c++ )
			{
				Pointer( float ) pRow = GetPixelRow( j + 1 , c ) + _size;
				Pointer( SyncType ) rRow = ( Pointer( SyncType ) )( right + offset + c * _padSize * sizeof( SyncType ) );
				for( int i=0 ; i<_padSize ; i++ ) pRow[i] = float( rRow[i] );
			}
		}
		else
		{
			for( int c=0 ; c<Channels ; c++ )
			{
				Pointer( float ) pRow = GetPixelRow( j , c ) + _size - _padSize;
				Pointer( SyncType ) rRow = ( Pointer( SyncType ) )( right + c * _padSize * sizeof( SyncType ) );
				for( int i=0 ; i<_padSize ; i++ ) rRow[i] = SyncType( pRow[i] );
			}
			for( int c=0 ; c<Channels ; c++ )
			{
				Pointer( float ) pRow = GetPixelRow( j + 1 , c ) + _size - _padSize;
				Pointer( SyncType ) rRow = ( Pointer( SyncType ) )( right + offset + c * _padSize * sizeof( SyncType ) );
				for( int i=0 ; i<_padSize ; i++ ) rRow[i] = SyncType( pRow[i] );
			}
			if ( !rightStream->write( ( Pointer( byte ) )right , sizeof( SyncType ) * Channels * xSize ) ) exit(0);
		}
		FreeArray( right );
	}
}
template< int Channels , class PixelType , class StorageType , class SyncType >
void SocketedStreamingLaplacian< Channels , PixelType , StorageType , SyncType >::SyncImageLeftRight( int y )
{
	if( y&1 )
	{
		SyncImageLeft ( y , false );
		SyncImageRight( y , true  );
	}
	else
	{
		SyncImageRight( y , false );
		SyncImageLeft ( y , true  );
	}
}
template< int Channels , class PixelType , class StorageType , class SyncType >
void SocketedStreamingLaplacian< Channels , PixelType , StorageType , SyncType >::UnSetRestriction(void)
{
	for( int i=0 ; i<Degree ; i++ )
		if( localPAccum[i] )
		{
			localPAccum[i] -= _padSize128;
			FreeArray( localPAccum[i] );
		}
	for( int i=0 ; i<ISize*Channels ; i++ ) FreeArray( PStream[i] );
	if( parent )	parent->UnSetRestriction();
}

template< int Channels , class PixelType , class StorageType , class SyncType >
void SocketedStreamingLaplacian< Channels , PixelType , StorageType , SyncType >::SetRestriction(void)
{
	for( int i=0 ; i<ISize*Channels ; i++ ) PStream[i] = AllocArray< __m128 >( _paddedSize128 , ALIGNMENT , "SocketedStreamingLaplacian::PStream" );
	if(parent)	parent->SetRestriction();
}
template< int Channels , class PixelType , class StorageType , class SyncType >
void SocketedStreamingLaplacian< Channels , PixelType , StorageType , SyncType >::SetInteriorRestriction( Pointer( float ) lB,int c,int idx,int sRestrict,int eRestrict)
{
	int jj=Degree*3;
	{
		ConstPointer( float ) localPPtrs[2*Degree+1];
		for( int y=0 ; y<2*Degree+1 ; y++ ) localPPtrs[y]=GetPixelRow(idx-Degree+y,c);
		{
			ConstPointer( __m128 ) pPtrs[] =
			{
				(ConstPointer( __m128 ) )localPPtrs[0],
				(ConstPointer( __m128 ) )localPPtrs[1],
				(ConstPointer( __m128 ) )localPPtrs[3],
				(ConstPointer( __m128 ) )localPPtrs[4]
			};

			// Perform the accumulation of the vertically symmetric rows and place them into __m128 buffers
			for(int i=-WordPerDegree+((sRestrict-_start)>>2);i<((eRestrict-_start)>>2)+WordPerDegree;i++)
			{
				localPAccum[0][i]=_mm_add_ps(pPtrs[0][i],pPtrs[3][i]);
				localPAccum[1][i]=_mm_add_ps(pPtrs[1][i],pPtrs[2][i]);
			}
		}


		float dotSum;
		__m128 dSum;
		__declspec (align(16)) float scratch[4];
		ConstPointer( __m128 ) pPtrs[] = { localPAccum[0] , localPAccum[1] , ( ConstPointer( __m128 ) )localPPtrs[2] };
		int s,e;

		// Offset 0
		{
			__m128 mValues[]=
			{
				lapTemplates[jj+1].matrixValues[0][0],
				lapTemplates[jj+1].matrixValues[0][1],
				lapTemplates[jj+1].matrixValues[0][2]
			};
#if SUPPORT_CYLINDRICAL
			if( sRestrict==0 && _periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if(sRestrict==0 && !_spherical)
#endif // SUPPORT_CYLINDRICAL
			{
				dotSum=0;
				lB[0]=GetInteriorLaplacianValue0(lapTemplates[jj].matrixValues[0],pPtrs,dotSum,0);
				s=4;
			}
			else
			{
				SetInteriorDotSum(mValues,pPtrs,-1+((sRestrict-_start)>>2),dSum);
				s=sRestrict-_start;
				_mm_store_ps(scratch,dSum);
				dotSum=scratch[2]+scratch[3];
			}
#if SUPPORT_CYLINDRICAL
			if( eRestrict==major && _periodicType==NO_PERIODIC )	e=eRestrict-_start-4;
#else // !SUPPORT_CYLINDRICAL
			if(eRestrict==major && !_spherical)	e=eRestrict-_start-4;
#endif // SUPPORT_CYLINDRICAL
			else								e=eRestrict-_start;

			for(int i=s;i<e;i+=4)				lB[i]=GetInteriorLaplacianValue0(mValues,pPtrs,dotSum,i);
			int i=eRestrict-_start-4;
#if SUPPORT_CYLINDRICAL
			if( eRestrict==major && _periodicType==NO_PERIODIC )	lB[i]=GetInteriorLaplacianValue0(lapTemplates[jj+2].matrixValues[0],pPtrs,dotSum,i);
#else // !SUPPORT_CYLINDRICAL
			if(eRestrict==major && !_spherical)	lB[i]=GetInteriorLaplacianValue0(lapTemplates[jj+2].matrixValues[0],pPtrs,dotSum,i);
#endif // SUPPORT_CYLINDRICAL
		}
		// Offset 1
		{
			__m128 mValues[]=
			{
				lapTemplates[jj+1].matrixValues[1][0],
				lapTemplates[jj+1].matrixValues[1][1],
				lapTemplates[jj+1].matrixValues[1][2]
			};
#if SUPPORT_CYLINDRICAL
			if( sRestrict==0 && _periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if(sRestrict==0 && !_spherical)
#endif // SUPPORT_CYLINDRICAL
			{
				dotSum=0;
				lB[1]=GetInteriorLaplacianValue1(lapTemplates[jj].matrixValues[1],pPtrs,dotSum,1);
				s=5;
			}
			else
			{
				SetInteriorDotSum(mValues,pPtrs,-1+((sRestrict-_start)>>2),dSum);
				s=sRestrict-_start+1;
				_mm_store_ps(scratch,dSum);
				dotSum=scratch[3];
			}
#if SUPPORT_CYLINDRICAL
			if( eRestrict==major && _periodicType==NO_PERIODIC )	e=eRestrict-_start-4;
#else // !SUPPORT_CYLINDRICAL
			if(eRestrict==major && !_spherical)	e=eRestrict-_start-4;
#endif // SUPPORT_CYLINDRICAL
			else								e=eRestrict-_start;
			for(int i=s;i<e;i+=4)			lB[i]=GetInteriorLaplacianValue1(mValues,pPtrs,dotSum,i);
			int i=eRestrict-_start-3;
#if SUPPORT_CYLINDRICAL
			if( eRestrict==major && _periodicType==NO_PERIODIC )lB[i]=GetInteriorLaplacianValue1(lapTemplates[jj+2].matrixValues[1],pPtrs,dotSum,i);
#else // !SUPPORT_CYLINDRICAL
			if(eRestrict==major && !_spherical)lB[i]=GetInteriorLaplacianValue1(lapTemplates[jj+2].matrixValues[1],pPtrs,dotSum,i);
#endif // SUPPORT_CYLINDRICAL
		}
		// Offset 2
		{
			__m128 mValues[]=
			{
				lapTemplates[jj+1].matrixValues[2][0],
				lapTemplates[jj+1].matrixValues[2][1],
				lapTemplates[jj+1].matrixValues[2][2]
			};
#if SUPPORT_CYLINDRICAL
			if( sRestrict==0 && _periodicType==NO_PERIODIC )	SetInteriorDotSum(lapTemplates[jj].matrixValues[2],pPtrs,0,dSum);
#else // !SUPPORT_CYLINDRICAL
			if(sRestrict==0 && !_spherical)	SetInteriorDotSum(lapTemplates[jj].matrixValues[2],pPtrs,0,dSum);
#endif // SUPPORT_CYLINDRICAL
			else							SetInteriorDotSum(mValues,pPtrs,(sRestrict-_start)>>2,dSum);
			s=(sRestrict-_start)+2;
			_mm_store_ps(scratch,dSum);
			dotSum=scratch[0]+scratch[1]+scratch[2]+scratch[3];
#if SUPPORT_CYLINDRICAL
			if( eRestrict==major && _periodicType==NO_PERIODIC)	e=eRestrict-_start-8;
#else // !SUPPORT_CYLINDRICAL
			if(eRestrict==major && !_spherical)	e=eRestrict-_start-8;
#endif // SUPPORT_CYLINDRICAL
			else								e=eRestrict-_start;
			for(int i=s;i<e;i+=4)	lB[i]=GetInteriorLaplacianValue2(mValues,pPtrs,dotSum,i);
#if SUPPORT_CYLINDRICAL
			if( eRestrict==major && _periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if(eRestrict==major && !_spherical)
#endif // SUPPORT_CYLINDRICAL
			{
				int i=eRestrict-_start-6;
				lB[i]=GetInteriorLaplacianValue2(lapTemplates[jj+2].matrixValues[2],pPtrs,dotSum,i);
				i+=4;
				lB[i]=dotSum;
			}
		}
		// Offset 3
		{
			__m128 mValues[]=
			{
				lapTemplates[jj+1].matrixValues[3][0],
				lapTemplates[jj+1].matrixValues[3][1],
				lapTemplates[jj+1].matrixValues[3][2]
			};
#if SUPPORT_CYLINDRICAL
			if( sRestrict==0 && _periodicType==NO_PERIODIC )	SetInteriorDotSum(lapTemplates[jj].matrixValues[3],pPtrs,0,dSum);
#else // !SUPPORT_CYLINDRICAL
			if(sRestrict==0 && !_spherical)	SetInteriorDotSum(lapTemplates[jj].matrixValues[3],pPtrs,0,dSum);
#endif // SUPPORT_CYLINDRICAL
			else							SetInteriorDotSum(mValues,pPtrs,(sRestrict-_start)>>2,dSum);
			s=(sRestrict-_start)+3;
			_mm_store_ps(scratch,dSum);
			dotSum=scratch[1]+scratch[2]+scratch[3];
#if SUPPORT_CYLINDRICAL
			if( eRestrict==major && _periodicType==NO_PERIODIC )	e=eRestrict-_start-8;
#else // !SUPPORT_CYLINDRICAL
			if(eRestrict==major && !_spherical)	e=eRestrict-_start-8;
#endif // SUPPORT_CYLINDRICAL
			else								e=eRestrict-_start;
			for(int i=s;i<e;i+=4)	lB[i]=GetInteriorLaplacianValue3(mValues,pPtrs,dotSum,i);
#if SUPPORT_CYLINDRICAL
			if( eRestrict==major && _periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if(eRestrict==major && !_spherical)
#endif // SUPPORT_CYLINDRICAL
			{
				int i=eRestrict-_start-5;
				lB[i]=GetInteriorLaplacianValue3(lapTemplates[jj+2].matrixValues[3],pPtrs,dotSum,i);
				i+=4;
				lB[i]=dotSum;
			}
		}
	}
}
template< int Channels , class PixelType , class StorageType , class SyncType >
void SocketedStreamingLaplacian< Channels , PixelType , StorageType , SyncType >::SetRestriction( Pointer( float ) lB,int c,int idx,int sRestrict,int eRestrict)
{
#if DEBUG_FLAGS
	if(debugFlags.noResidual)	return;
#endif // DEBUG_FLAGS
	int jj;
	if		(idx<Degree)			jj=idx;
	else if	(idx>minor-1-Degree)	jj=2*Degree+(idx-(minor-1));
	else							jj=Degree;

#if SUPPORT_CYLINDRICAL
	if( !_periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
	if(!_spherical)
#endif // SUPPORT_CYLINDRICAL
	{
		if(sRestrict<0)		sRestrict=0;
		if(eRestrict>major)	eRestrict=major;
	}
#if SUPPORT_CYLINDRICAL
	if( jj==Degree || _periodicType==SPHERICAL_PERIODIC )	return SetInteriorRestriction(lB,c,idx,sRestrict,eRestrict);
#else // !SUPPORT_CYLINDRICAL
	if(jj==Degree || _spherical)	return SetInteriorRestriction(lB,c,idx,sRestrict,eRestrict);
#endif // SUPPORT_CYLINDRICAL
	jj*=3;

	{
		ConstPointer( float ) localPPtrs[2*Degree+1];
		for(int y=0;y<2*Degree+1;y++)
			if(idx-Degree+y>=0 && idx-Degree+y<minor)	localPPtrs[y] = GetPixelRow(idx-Degree+y,c);
			else										localPPtrs[y] = NullPointer< float >( );

		float dotSum;
		__m128 dSum;
		__declspec (align(16)) float scratch[4];
		ConstPointer( __m128 ) pPtrs[] =
		{
			( ConstPointer( __m128 ) )localPPtrs[0],
			( ConstPointer( __m128 ) )localPPtrs[1],
			( ConstPointer( __m128 ) )localPPtrs[2],
			( ConstPointer( __m128 ) )localPPtrs[3],
			( ConstPointer( __m128 ) )localPPtrs[4]
		};
		int s,e;
		// Offset 0
		{
			__m128 mValues[]=
			{
				lapTemplates[jj+1].matrixValues[0][0],
				lapTemplates[jj+1].matrixValues[0][1],
				lapTemplates[jj+1].matrixValues[0][2],
				lapTemplates[jj+1].matrixValues[0][3],
				lapTemplates[jj+1].matrixValues[0][4]
			};
#if SUPPORT_CYLINDRICAL
			if( sRestrict==0 && _periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if(sRestrict==0 && !_spherical)
#endif // SUPPORT_CYLINDRICAL
			{
				dotSum=0;
				lB[0] = GetLaplacianValue0(lapTemplates[jj].matrixValues[0],pPtrs,dotSum,0);
				s=4;
			}
			else
			{
				SetDotSum(mValues,pPtrs,-1+((sRestrict-_start)>>2),dSum);
				s=sRestrict-_start;
				_mm_store_ps(scratch,dSum);
				dotSum=scratch[2]+scratch[3];
			}
#if SUPPORT_CYLINDRICAL
			if( eRestrict==major && _periodicType==NO_PERIODIC )	e=eRestrict-_start-4;
#else // !SUPPORT_CYLINDRICAL
			if(eRestrict==major && !_spherical)	e=eRestrict-_start-4;
#endif // SUPPORT_CYLINDRICAL
			else								e=eRestrict-_start;
			for(int i=s;i<e;i+=4)				lB[i]=GetLaplacianValue0(mValues,pPtrs,dotSum,i);
			int i=eRestrict-_start-4;
#if SUPPORT_CYLINDRICAL
			if( eRestrict==major && _periodicType==NO_PERIODIC )	lB[i]=GetLaplacianValue0(lapTemplates[jj+2].matrixValues[0],pPtrs,dotSum,i);
#else // !SUPPORT_CYLINDRICAL
			if(eRestrict==major && !_spherical)	lB[i]=GetLaplacianValue0(lapTemplates[jj+2].matrixValues[0],pPtrs,dotSum,i);
#endif // SUPPORT_CYLINDRICAL
		}
		// Offset 1
		{
			__m128 mValues[]=
			{
				lapTemplates[jj+1].matrixValues[1][0],
				lapTemplates[jj+1].matrixValues[1][1],
				lapTemplates[jj+1].matrixValues[1][2],
				lapTemplates[jj+1].matrixValues[1][3],
				lapTemplates[jj+1].matrixValues[1][4]
			};
#if SUPPORT_CYLINDRICAL
			if( sRestrict==0 && _periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if(sRestrict==0 && !_spherical)
#endif // SUPPORT_CYLINDRICAL
			{
				dotSum=0;
				lB[1]=GetLaplacianValue1(lapTemplates[jj].matrixValues[1],pPtrs,dotSum,1);
				s=5;
			}
			else
			{
				SetDotSum(mValues,pPtrs,-1+((sRestrict-_start)>>2),dSum);
				s=sRestrict-_start+1;
				_mm_store_ps(scratch,dSum);
				dotSum=scratch[3];
			}
#if SUPPORT_CYLINDRICAL
			if( eRestrict==major &&  _periodicType==NO_PERIODIC )	e=eRestrict-_start-4;
#else // !SUPPORT_CYLINDRICAL
			if(eRestrict==major && !_spherical)	e=eRestrict-_start-4;
#endif // SUPPORT_CYLINDRICAL3
			else								e=eRestrict-_start;
			for(int i=s;i<e;i+=4)				lB[i]=GetLaplacianValue1(mValues,pPtrs,dotSum,i);
			int i=eRestrict-_start-3;
#if SUPPORT_CYLINDRICAL
			if( eRestrict==major && _periodicType==NO_PERIODIC )	lB[i]=GetLaplacianValue1(lapTemplates[jj+2].matrixValues[1],pPtrs,dotSum,i);
#else // !SUPPORT_CYLINDRICAL
			if(eRestrict==major && !_spherical)	lB[i]=GetLaplacianValue1(lapTemplates[jj+2].matrixValues[1],pPtrs,dotSum,i);
#endif // SUPPORT_CYLINDRICAL
		}
		// Offset 2
		{
			__m128 mValues[]=
			{
				lapTemplates[jj+1].matrixValues[2][0],
				lapTemplates[jj+1].matrixValues[2][1],
				lapTemplates[jj+1].matrixValues[2][2],
				lapTemplates[jj+1].matrixValues[2][3],
				lapTemplates[jj+1].matrixValues[2][4]
			};
#if SUPPORT_CYLINDRICAL
			if( sRestrict==0 && _periodicType==NO_PERIODIC )	SetDotSum(lapTemplates[jj].matrixValues[2],pPtrs,0,dSum);
#else // !SUPPORT_CYLINDRICAL
			if(sRestrict==0 && !_spherical)	SetDotSum(lapTemplates[jj].matrixValues[2],pPtrs,0,dSum);
#endif // SUPPORT_CYLINDRICAL3
			else							SetDotSum(mValues,pPtrs,(sRestrict-_start)>>2,dSum);
			s=(sRestrict-_start)+2;
			_mm_store_ps(scratch,dSum);
			dotSum=scratch[0]+scratch[1]+scratch[2]+scratch[3];
#if SUPPORT_CYLINDRICAL
			if( eRestrict==major && _periodicType==NO_PERIODIC )	e=eRestrict-_start-8;
#else // !SUPPORT_CYLINDRICAL
			if(eRestrict==major && !_spherical)	e=eRestrict-_start-8;
#endif // SUPPORT_CYLINDRICAL
			else								e=eRestrict-_start;
			for(int i=s;i<e;i+=4)	lB[i]=GetLaplacianValue2(mValues,pPtrs,dotSum,i);
#if SUPPORT_CYLINDRICAL
			if( eRestrict==major && _periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if(eRestrict==major && !_spherical)
#endif // SUPPORT_CYLINDRICAL
			{
				int i=eRestrict-_start-6;
				lB[i]=GetLaplacianValue2(lapTemplates[jj+2].matrixValues[2],pPtrs,dotSum,i);
				i+=4;
				lB[i]=dotSum;
			}
		}
		// Offset 3
		{
			__m128 mValues[]=
			{
				lapTemplates[jj+1].matrixValues[3][0],
				lapTemplates[jj+1].matrixValues[3][1],
				lapTemplates[jj+1].matrixValues[3][2],
				lapTemplates[jj+1].matrixValues[3][3],
				lapTemplates[jj+1].matrixValues[3][4]
			};
#if SUPPORT_CYLINDRICAL
			if( sRestrict==0 && _periodicType==NO_PERIODIC )	SetDotSum(lapTemplates[jj].matrixValues[3],pPtrs,0,dSum);
#else // !SUPPORT_CYLINDRICAL
			if(sRestrict==0 && !_spherical)	SetDotSum(lapTemplates[jj].matrixValues[3],pPtrs,0,dSum);
#endif // SUPPORT_CYLINDRICAL
			else							SetDotSum(mValues,pPtrs,(sRestrict-_start)>>2,dSum);
			s=(sRestrict-_start)+3;
			_mm_store_ps(scratch,dSum);
			dotSum=scratch[1]+scratch[2]+scratch[3];
#if SUPPORT_CYLINDRICAL
			if( eRestrict==major && _periodicType==NO_PERIODIC )	e=eRestrict-_start-8;
#else // !SUPPORT_CYLINDRICAL
			if(eRestrict==major && !_spherical)	e=eRestrict-_start-8;
#endif // SUPPORT_CYLINDRICAL
			else								e=eRestrict-_start;
			for(int i=s;i<e;i+=4)	lB[i]=GetLaplacianValue3(mValues,pPtrs,dotSum,i);
#if SUPPORT_CYLINDRICAL
			if( eRestrict==major && _periodicType==NO_PERIODIC )
#else // !SUPPORT_CYLINDRICAL
			if(eRestrict==major && !_spherical)
#endif // SUPPORT_CYLINDRICAL
			{
				int i=eRestrict-_start-5;
				lB[i]=GetLaplacianValue3(lapTemplates[jj+2].matrixValues[3],pPtrs,dotSum,i);
				i+=4;
				lB[i]=dotSum;
			}
		}
	}
}
template< int Channels , class PixelType , class StorageType , class SyncType >
bool SocketedStreamingLaplacian< Channels , PixelType , StorageType , SyncType >::IterateRestriction(void)
{
	if(index>=minor)	return false;

	int idx = index+Degree+1;

	// Load the next row of pixel/label values
	if(idx<_iHeight)
	{
		{
			Pointer( PixelType ) pixels = ( Pointer( PixelType ) )(*pixelStream)[idx];
			for(int c=0;c<Channels;c++)
			{
				Pointer( float ) pixelRow = GetPixelRow( idx , c );
				if( pixelStream->SeparateColors() )
					for ( int x = 0 ; x < _iWidth ; x++)
					{
						pixelRow[x]=float(pixels[_iWidth*c+x]);
						average[c]+=pixelRow[x];
					}
				else
					for ( int x = 0 ; x < _iWidth ; x++)
					{
						pixelRow[x]=float(pixels[x*Channels+c]);
						average[c]+=pixelRow[x];
					}
			}
			pixelStream->advance();
		}
	}
	else if (idx<minor) for(int c=0;c<Channels;c++)	memset( GetPixelRow(idx,c) , 0 , sizeof(float) * _size );

#if SUPPORT_CYLINDRICAL
	if( _periodicType==SPHERICAL_PERIODIC && idx==Degree-1 ) for( int d=0 ; d<Degree ; d++ ) SyncImageHead( d , false ) , SyncImageHead( d , true );
	if( _periodicType==SPHERICAL_PERIODIC && idx== minor-1 ) for( int d=minor-Degree ; d<minor ; d++ ) SyncImageTail( d , false ) , SyncImageTail( d , true );
#else // !SUPPORT_CYLINDRICAL
	if( _spherical && idx==Degree-1 ) for( int d=0 ; d<Degree ; d++ ) SyncImageHead( d , false ) , SyncImageHead( d , true );
	if( _spherical && idx== minor-1 ) for( int d=minor-Degree ; d<minor ; d++ ) SyncImageTail( d , false ) , SyncImageTail( d , true );
#endif // SUPPORT_CYLINDRICAL

	// If this is the first row, perform the preliminary left/right synchronization
#if SUPPORT_CYLINDRICAL
	if( _periodicType==SPHERICAL_PERIODIC && idx==Degree-1 )
#else // !SUPPORT_CYLINDRICAL
	if( _spherical && idx==Degree-1 )
#endif // SUPPORT_CYLINDRICAL
	{
		SyncImageLeft( -Degree , false ) , SyncImageRight( -Degree , false );
		SyncImageLeft( -Degree , true  ) , SyncImageRight( -Degree , true  );
	}
#if SUPPORT_CYLINDRICAL
	if( _periodicType!=SPHERICAL_PERIODIC && idx==0 )
#else // !SUPPORT_CYLINDRICAL
	if( !_spherical && idx==0 )
#endif // SUPPORT_CYLINDRICAL
	{
		SyncImageLeft( 0 , false ) , SyncImageRight( 0 , false );
		SyncImageLeft( 0 , true  ) , SyncImageRight( 0 , true  );
	}

	idx = index+Degree;
	if( idx>=0 )
	{
#if SUPPORT_CYLINDRICAL
		if( _periodicType==SPHERICAL_PERIODIC && idx==Degree-1)
#else // !SUPPORT_CYLINDRICAL
		if( _spherical && idx==Degree-1)
#endif // SUPPORT_CYLINDRICAL
		{
			for(int y = -Degree ; y <= 0 ; y++) SyncImageLeftRight( y );
		}
#if SUPPORT_CYLINDRICAL
		if(_periodicType==SPHERICAL_PERIODIC && idx==minor-1)
#else // !SUPPORT_CYLINDRICAL
		if(_spherical && idx==minor-1)
#endif // SUPPORT_CYLINDRICAL
		{
			for( int y=minor ; y < minor + Degree ; y++ ) SyncImageLeftRight( y );
		}
		SyncImageLeftRight( idx );
		if(parent)	parent->IterateRestriction();
	}
	index++;
	return true;
}
template< int Channels , class PixelType , class StorageType , class SyncType >
void SocketedStreamingLaplacian< Channels , PixelType , StorageType , SyncType >::SolveRestriction(void)
{
	// Run to completion...
	while( IterateRestriction() ){;}
	// ...and finish up the trailing parent
	if(parent)	parent->SolveRestriction();
}
