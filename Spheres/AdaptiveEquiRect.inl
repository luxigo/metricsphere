/////////////////////////////
// AdaptiveEquiRectangular //
/////////////////////////////
template< class Real , class PolyReal >
void AdaptiveEquiRectangular< Real , PolyReal >::elementWeights( Vector< Real >& weights ) const
{
	weights.Resize( _Dim );
	if( _stencilTable )
	{
		int idx = 0;
		for( int m=0 ; m<_M ; m++ )
		{
			int N = rowDimension[ m ];
			Real weight;
			if( m==0 || m==_M-1 )	weight = _stencilTable->thetaWeight( m , _M ) * 2. * M_PI;
			else					weight = _stencilTable->thetaWeight( m , _M ) * _stencilTable->phiWeight( N );
			for( int n=0 ; n<N ; n++ ) weights[ idx++ ] = weight;
		}
		return;
	}
	PolyReal a = M_PI / _M;
	PPolynomial< 2 , PolyReal > thetaBase , phiBase;
	std::vector< PolyReal > phiWeights;
	int base=0;
	
	while( (1<<base)!=(_M<<1) ) base++;
	phiWeights.resize( base+1 );
	EquiRectangular< Real , PolyReal >::SetThetaBaseFunction( 1 , _M , thetaBase , false );
	PolyReal thetaSinWeight = thetaBase.integralSine  ( a , 0. , 3. ) * a;
	PolyReal thetaCosWeight = thetaBase.integralCosine( a , 0. , 3. ) * a;

	for( int i=0 ; i<=base ; i++ )
	{
		if( i )
		{
			EquiRectangular< Real , PolyReal >::SetPhiBaseFunction( 0 , 1<<i , phiBase , false );
			PolyReal b = 2.0 * M_PI / (1<<i);
			phiWeights[i] = phiBase.integral( -2. , 3. ) * b;
		}
		else phiWeights[i] = 2.0 * M_PI;
	}

	int idx = 0;
	for( int m=0 ; m<_M ; m++ )
	{
		int N = rowDimension[ m ];
		Real weight;
		int base=0;
		while( (1<<base)!=N ) base++;
		if( m==0 || m==_M-1 )
		{
			double start = m-1. , end = m+2.;
			if( start<0  ) start =  0;
			if(   end>_M ) end   = _M;
			EquiRectangular< Real , PolyReal >::SetThetaBaseFunction( m , _M , thetaBase , false );
			weight = Real( thetaBase.integralSine( a , start , end ) * a * PolyReal( 2.0 * M_PI ) );
		}
		else weight = Real( ( thetaSinWeight * cos( a * PolyReal( m-1 ) ) + thetaCosWeight * sin( a * PolyReal( m-1 ) ) ) * phiWeights[base] );
		for( int n=0 ; n<N ; n++ ) weights[ idx++ ] = weight;
	}
}
template< int Dim , class Real >
class Stencil
{
public:
	Real value[Dim];
};
template< class Real , class PolyReal >
void AdaptiveEquiRectangular< Real , PolyReal >::laplacianMatrix( SparseMatrix< Real >& lMatrix , double iWeight , double lWeight ) const
{
	int start[5];
	double thetaStencil[5] , laplacianThetaStencil[5] , laplacianD2ThetaStencil[5] , phiStencils[5][5] , laplacianD2PhiStencils[5][5] , phiWeight;
	double            upPhiStencils[5][8] ,            downEvenPhiStencils[5][4] ,            downOddPhiStencils[5][4];
	double upLaplacianD2PhiStencils[5][8] , downEvenLaplacianD2PhiStencils[5][4] , downOddLaplacianD2PhiStencils[5][4];
	double stencils[2][5][8];
	lMatrix.Resize( _Dim );

	for( int i=0 ; i<_M ; i+=_M-1 )
	{
		if( _stencilTable )
		{
			_stencilTable->setThetaStencil           (            thetaStencil , i , _M );
			_stencilTable->setLaplacianD2ThetaStencil( laplacianD2ThetaStencil , i , _M );
		}
		else
		{
			EquiRectangular< double , PolyReal >::ThetaStencil           ( i , _M ,            thetaStencil );
			EquiRectangular< double , PolyReal >::LaplacianD2ThetaStencil( i , _M , laplacianD2ThetaStencil );
		}

		int ii = 0;
		int idx = index( i , 0 );
		int count = 0;
		for( int k=-2 ; k<=2 ; k++ ) if( i+k>=0 && i+k<_M ) count += rowDimension[i+k];
		lMatrix.SetGroupSize( idx , count );

		for( int k=-2 ; k<=2 ; k++ )
		{
			if( i+k<0 || i+k>=_M ) continue;
			else if( k==0 )
			{
				lMatrix[idx][ii  ].N = idx;
				lMatrix[idx][ii++].Value = - laplacianD2ThetaStencil[k+2] * 2.0 * M_PI * lWeight + thetaStencil[k+2] * 2.0 * M_PI * iWeight;
			}
			else
			{
				if( _stencilTable ) phiWeight = _stencilTable->phiWeight( rowDimension[i+k] );
				else phiWeight = EquiRectangular< Real , PolyReal >::PhiWeight( rowDimension[i+k] );
				int start = index( i+k , 0 );
				PolyReal val = - laplacianD2ThetaStencil[k+2] * phiWeight * lWeight + thetaStencil[k+2] * phiWeight * iWeight;
				for( int l=0 ; l<rowDimension[i+k] ; l++ ) lMatrix[idx][ii++] = MatrixEntry< Real >( start+l , Real( val ) );
			}
		}
	}
	for( int i=1 ; i<_M-1 ; i++ )
	{
		if( _stencilTable )
		{
			phiWeight = _stencilTable->phiWeight( rowDimension[i] );
			_stencilTable->setThetaStencil           (            thetaStencil , i , _M );
			_stencilTable->setLaplacianThetaStencil  (   laplacianThetaStencil , i , _M , _Samples );
			_stencilTable->setLaplacianD2ThetaStencil( laplacianD2ThetaStencil , i , _M );
		}
		else
		{
			phiWeight = EquiRectangular< double , PolyReal >::PhiWeight( rowDimension[i] );
			EquiRectangular< double , PolyReal >::ThetaStencil ( i , _M , thetaStencil );
			EquiRectangular< double , PolyReal >::LaplacianThetaStencil  ( i , _M ,   laplacianThetaStencil , _Samples );
			EquiRectangular< double , PolyReal >::LaplacianD2ThetaStencil( i , _M , laplacianD2ThetaStencil );
		}

		int count = 0;
		for( int k=-2 ; k<=2 ; k++ )
			if( i+k>=0 && i+k<_M )
			{
				int dim = rowDimension[i]>rowDimension[i+k] ? rowDimension[i] : rowDimension[i+k];
				if( _stencilTable )
				{
					_stencilTable->setPhiStencil           (            phiStencils[k+2] , dim );
					_stencilTable->setLaplacianD2PhiStencil( laplacianD2PhiStencils[k+2] , dim);
				}
				else
				{
					EquiRectangular< double , PolyReal >::PhiStencil           ( dim ,            phiStencils[k+2] );
					EquiRectangular< double , PolyReal >::LaplacianD2PhiStencil( dim , laplacianD2PhiStencils[k+2] );
				}
				UpStencil  (            phiStencils[k+2] ,            upPhiStencils[k+2] );
				UpStencil  ( laplacianD2PhiStencils[k+2] , upLaplacianD2PhiStencils[k+2] );
				DownStencil(            phiStencils[k+2] ,            downEvenPhiStencils[k+2] ,            downOddPhiStencils[k+2] );
				DownStencil( laplacianD2PhiStencils[k+2] , downEvenLaplacianD2PhiStencils[k+2] , downOddLaplacianD2PhiStencils[k+2] );
				if( i+k==0 || i+k==_M-1 )							count += 1;
				else if(   rowDimension[i]==  rowDimension[i+k] )	count += 5;
				else if(   rowDimension[i]==2*rowDimension[i+k] )	count += 4;
				else if( 2*rowDimension[i]==  rowDimension[i+k] )	count += 8;
				else fprintf( stderr , "Badness in AdaptiveEquiRectangular< Real , PolyReal >::laplacianMatrix\n" ) , exit(0);
				start[k+2] = index( i+k , 0 );
				if( i+k==0 || i+k==_M-1 )							// Pole neighbor
					stencils[0][k+2][0] = - laplacianD2ThetaStencil[2+k] * phiWeight * lWeight + thetaStencil[2+k] * phiWeight * iWeight;
				else if( rowDimension[i]==rowDimension[i+k] )		// Regular neighbor
					for( int l=-2 ; l<=2 ; l++ )
						stencils[0][k+2][l+2] =
						- ( laplacianD2ThetaStencil[2+k] * phiStencils[2+k][2+l] + laplacianThetaStencil[2+k] * laplacianD2PhiStencils[2+k][2+l] ) * lWeight
						+ thetaStencil[2+k] * phiStencils[2+k][2+l] * iWeight;
				else if( rowDimension[i]< rowDimension[i+k] )		// Finer neighbor
					for( int l=-3 ; l<=4 ; l++ )
						stencils[0][k+2][l+3] =
						- ( laplacianD2ThetaStencil[2+k] * upPhiStencils[2+k][l+3] + laplacianThetaStencil[2+k] * upLaplacianD2PhiStencils[2+k][l+3] ) * lWeight
						+  thetaStencil[2+k] * upPhiStencils[2+k][l+3] * iWeight;
				else if( rowDimension[i]> rowDimension[i+k] )		// Coarser neighbor
				{
					for( int l=-2 ; l<=1 ;l++ )
						stencils[0][k+2][l+2] =
						- ( laplacianD2ThetaStencil[2+k] * downEvenPhiStencils[2+k][l+2] + laplacianThetaStencil[2+k] * downEvenLaplacianD2PhiStencils[2+k][l+2] ) * lWeight
						+ thetaStencil[2+k] * downEvenPhiStencils[2+k][l+2] * iWeight;
					for( int l=-1 ; l<=2 ;l++ )
						stencils[1][k+2][l+1] =
						- ( laplacianD2ThetaStencil[2+k] * downOddPhiStencils[2+k][l+1] + laplacianThetaStencil[2+k] * downOddLaplacianD2PhiStencils[2+k][l+1] ) * lWeight
						+ thetaStencil[2+k] * downOddPhiStencils[2+k][l+1] * iWeight;
				}
			}
		for( int j=0 ; j<rowDimension[i] ; j++ )
		{
			int ii = 0;
			int idx = start[2] + j;
			lMatrix.SetGroupSize( idx , count );
			for( int k=-2 ; k<=2 ; k++ )
			{
				if( i+k<0 || i+k>=_M ) continue;
				if( i+k==0 || i+k==_M-1 )
				{
					lMatrix[idx][ii  ].N = start[ k+2 ];
					lMatrix[idx][ii++].Value = stencils[0][k+2][0];
				}
				else if( rowDimension[i]==rowDimension[i+k] )
					for( int l=-2 ; l<=2 ; l++ )
					{
						lMatrix[idx][ii  ].N = start[k+2] + ( (j+l+rowDimension[i+k]) % rowDimension[i+k] );
						lMatrix[idx][ii++].Value = stencils[0][k+2][l+2];
					}
				else if( rowDimension[i]> rowDimension[i+k] )
				{
					if( j&1 )
						for( int l=-1 ; l<=2 ;l++ )
						{
							lMatrix[idx][ii  ].N = start[k+2] + ( (j/2+l+rowDimension[i+k]) % rowDimension[i+k] );
							lMatrix[idx][ii++].Value =	stencils[1][k+2][l+1];
						}
					else
						for( int l=-2 ; l<=1 ;l++ )
						{
							lMatrix[idx][ii  ].N = start[k+2] + ( (j/2+l+rowDimension[i+k]) % rowDimension[i+k] );
							lMatrix[idx][ii++].Value = stencils[0][k+2][l+2];
						}
				}
				else if( rowDimension[i]< rowDimension[i+k] )
					for( int l=-3 ; l<=4 ; l++ )
					{
						lMatrix[idx][ii  ].N = start[k+2] + ( (j*2+l+rowDimension[i+k]) % rowDimension[i+k] );
						lMatrix[idx][ii++].Value = stencils[0][k+2][l+3];
					}
			}
		}
	}
}
template< class Real , class PolyReal >
void AdaptiveEquiRectangular< Real , PolyReal >::equiRectangularMatrix( SparseMatrix< Real >& zaMatrix ) const
{
	// Assuming that N is a power of 2
	int logN = 1;
	std::vector< int > start;
	std::vector< std::vector< int > > pStencils;
	zaMatrix.Resize( dimension() , false );

	while( (1<<logN) < (_M<<1) ) logN++;
	pStencils.resize( logN );
	start.resize( logN );
	for( int i=0 ; i<logN ; i++ )
	{
		if( i==0 )
		{
			pStencils[i].resize( 1 );
			pStencils[i][0] = 1;
			start[i] = 0;
		}
		else
		{
			int dim = 2+pStencils[i-1].size()*2;
			if( dim>=(_M<<1) ) dim = (_M<<1);
			start[i] = 2*start[i-1] - 1;
			if( start[i] <= -(_M<<1) ) start[i] += (_M<<1);
			pStencils[i].resize( dim );
			for( int j=0 ; j<dim ; j++ ) pStencils[i][j] = 0;
			for( int j=0 ; j<pStencils[i-1].size() ; j++ )
			{
				pStencils[i][ (2*j  ) % dim ] += 1 * pStencils[i-1][j];
				pStencils[i][ (2*j+1) % dim ] += 3 * pStencils[i-1][j];
				pStencils[i][ (2*j+2) % dim ] += 3 * pStencils[i-1][j];
				pStencils[i][ (2*j+3) % dim ] += 1 * pStencils[i-1][j];
			}
		}
	}
	for( int i=0 ; i<_M ; i+=_M-1 )
	{
		int idx = index( i , 0 );
		zaMatrix.SetGroupSize( idx , 1 );
		zaMatrix[idx][0].N = EquiRectangular< Real , PolyReal >::Index( i , 0 , _M , _M<<1 );
		zaMatrix[idx][0].Value = 1;
	}
	for( int i=1 ; i<_M-1 ; i++ )
	{
		int idxStart = index( i , 0 );
		int myLogN = 0;
		while( (rowDimension[i]<<myLogN) < (_M<<1) ) myLogN++;
		for( int j=0 ; j<rowDimension[i] ; j++ )
		{
			zaMatrix.SetGroupSize( idxStart+j , pStencils[myLogN].size() );
			for( int k=0 ; k<pStencils[myLogN].size() ; k++ )
			{
				zaMatrix[idxStart+j][k].N = EquiRectangular< Real , PolyReal >::Index( i , start[myLogN]+k+j*(1<<myLogN) , _M , _M<<1 );
				zaMatrix[idxStart+j][k].Value = double(pStencils[myLogN][k]) / ( 1 << (2*myLogN) );
			}
		}
	}
}
template< class Real , class PolyReal >
template< int Dim >
void AdaptiveEquiRectangular< Real , PolyReal >::_UpStencil( const double inStencil[Dim] , double stencil[Dim+3] )
{
	double weights[] = { 0.25 , 0.75 , 0.75 , 0.25 };
	memset( stencil , 0 , sizeof( stencil ) );
	for( int i=0 ; i<Dim ; i++ ) for( int j=0 ; j<4 ; j++ ) stencil[i+j] += inStencil[i] * weights[j];
}

template< class Real , class PolyReal >
template< int Dim >
void AdaptiveEquiRectangular< Real , PolyReal >::_DownStencil( const double inStencil[Dim] , double evenStencil[(Dim>>1)+2] , double oddStencil[(Dim>>1)+2] )
{
	double weights[] = { 0.25 , 0.75 , 0.75 , 0.25 };
	memset( evenStencil , 0 , sizeof( evenStencil ) );
	memset(  oddStencil , 0 , sizeof(  oddStencil ) );
	for( int i=0 ; i<( (Dim>>1)+2 ) ; i++ )
		for( int j=0 ; j<3 ; j++ )
		{
			if( i+j-3>=0 && i+j-3<Dim ) evenStencil[i] += inStencil[i+j-3] * weights[j];
			if( i+j-2>=0 && i+j-2<Dim )  oddStencil[i] += inStencil[i+j-2] * weights[j];
		}
}

template< class Real , class PolyReal >
void AdaptiveEquiRectangular< Real , PolyReal >::UpStencil( const double inStencil[5] , double stencil[8] )
{
	stencil[0] = 0.25 * (                                                  inStencil[0] );
	stencil[1] = 0.25 * (                                 inStencil[0]*3 + inStencil[1] );
	stencil[2] = 0.25 * (                inStencil[0]*3 + inStencil[1]*3 + inStencil[2] );
	stencil[3] = 0.25 * ( inStencil[0] + inStencil[1]*3 + inStencil[2]*3 + inStencil[3] );
	stencil[4] = 0.25 * ( inStencil[1] + inStencil[2]*3 + inStencil[3]*3 + inStencil[4] );
	stencil[5] = 0.25 * ( inStencil[2] + inStencil[3]*3 + inStencil[4]*3                );
	stencil[6] = 0.25 * ( inStencil[3] + inStencil[4]*3                                 );
	stencil[7] = 0.25 * ( inStencil[4]                                                  );
}
template< class Real , class PolyReal >
void AdaptiveEquiRectangular< Real , PolyReal >::DownStencil( const double inStencil[5] , double evenStencil[4] , double oddStencil[4] )
{
	evenStencil[0] = 0.25 * (                                                  inStencil[0] );
	evenStencil[1] = 0.25 * (                inStencil[0]*3 + inStencil[1]*3 + inStencil[2] );
	evenStencil[2] = 0.25 * ( inStencil[1] + inStencil[2]*3 + inStencil[3]*3 + inStencil[4] );
	evenStencil[3] = 0.25 * ( inStencil[3] + inStencil[4]*3                                 );

	oddStencil[0] = 0.25 * (                                 inStencil[0]*3 + inStencil[1] );
	oddStencil[1] = 0.25 * ( inStencil[0] + inStencil[1]*3 + inStencil[2]*3 + inStencil[3] );
	oddStencil[2] = 0.25 * ( inStencil[2] + inStencil[3]*3 + inStencil[4]*3                );
	oddStencil[3] = 0.25 * ( inStencil[4]                                                  );
}
template< class Real , class PolyReal >
AdaptiveEquiRectangular< Real , PolyReal >::AdaptiveEquiRectangular( int M , int Samples , const SphericalStencilTable< double >* stencilTable )
{
	_M = M , _Samples = Samples , _stencilTable = stencilTable;
	rowDimension.resize( M );
	rowDimension[0] = rowDimension[M-1] = 1;
	_Dim = 2;

	for( int j=1 ; j<M-1 ; j++ )
	{
		double theta = M_PI * (j+0.5) / M;
		int cnt = 1;
		//  PI / M / sqrt( 2 ) <= sin( theta ) * 2 * PI / cnt < PI / M * sqrt( 2 )
		//               cnt/2 <= sin( theta ) * sqrt( 2 ) * M < cnt
		while( cnt < sin( theta ) * sqrt( 2.0 ) * _M ) cnt <<= 1;
		rowDimension[j] = cnt;
		_Dim += cnt;
	}
	// sanity check
	for( int j=2 ; j<M-1 ; j++ )
		if( rowDimension[j-1]*2 < rowDimension[j] || rowDimension[j-1] > rowDimension[j]*2 )
			fprintf( stderr , "Badness1 in AdaptiveEquiRectangular constructor\n" ) , exit(0);
	for( int j=2 ; j<M-2 ; j++ )
		if( rowDimension[j-1]!=rowDimension[j] && rowDimension[j]!=rowDimension[j+1] )
			fprintf( stderr , "Badness2 in AdaptiveEquiRectangular constructor\n" ) , exit(0);
}
template< class Real , class PolyReal >
int AdaptiveEquiRectangular< Real , PolyReal >::dimension( void ) const
{
	return _Dim;
}
template< class Real , class PolyReal >
int AdaptiveEquiRectangular< Real , PolyReal >::dimension( int idx ) const
{
	if( idx>=0 && idx<_M )	return rowDimension[idx];
	else return 0;
}
template< class Real , class PolyReal >
int AdaptiveEquiRectangular< Real , PolyReal >::index( int m , int n ) const
{
	int mm , nn;
	return index( m , n , mm , nn );

}
template< class Real , class PolyReal >
int AdaptiveEquiRectangular< Real , PolyReal >::index( int m , int n , int& mm , int & nn ) const
{
	int idx = 0;
	mm = m;
	nn = ( n + rowDimension[mm] ) % rowDimension[mm];
	for( int j=0 ; j<mm ; j++ ) idx += rowDimension[j];
	return idx + nn;
}
template< class Real , class PolyReal >
void AdaptiveEquiRectangular< Real , PolyReal >::index( int idx , int& m , int& n ) const
{
	if( idx<0 || idx>_Dim )
	{
		fprintf( stderr , "Out of bounds index in AdaptiveEquiRectangular< Real , PolyReal >::index\n" );
		return;
	}
	for( int i=0 ; i<_M ; i++ )
	{
		if( idx<rowDimension[i] )
		{
			m = i , n = idx;
			return;
		}
		idx -= rowDimension[i];
	}
}
template< class Real , class PolyReal >
void AdaptiveEquiRectangular< Real , PolyReal >::index( int m , int n , double& theta , double& phi ) const
{
	if		( m<=   0 ) theta = 0		, phi = 0;
	else if ( m>=_M-1 ) theta = M_PI	, phi = 0;
	else
	{
		theta = M_PI * ( m+0.5 ) / _M;
		phi = 2.0 * M_PI * ( n+0.5 ) / rowDimension[m];
	}
}
