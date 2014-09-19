///////////////////////
// ThetaStencilTable //
///////////////////////
template< class Real > ThetaStencilTable< Real >::ThetaStencilTable( void ) { _logHeight = false; }
template< class Real > ThetaStencilTable< Real >::~ThetaStencilTable( void ){ if( _logHeight ) for( int i=0 ; i<_valueTable.size() ; i++ ) delete[] _valueTable[i] ; _logHeight = 0 ; }
template< class Real >
Real ThetaStencilTable< Real >::thetaWeight( int m , int M ) const
{
	int logM=1;
	while( (1<<logM)<M ) logM++;
	if( (1<<logM)!=M || m<0 || m>=M || logM>_logHeight || logM<2 ) return EquiRectangular< Real >::ThetaWeight( m , M );
	else
	{
		if( m>=M/2 ) m = M-1-m;
		return _valueTable[logM][m].thetaWeight;
	}
}
template< class Real >
bool ThetaStencilTable< Real >::setThetaStencil( Real thetaStencil[5] , int m , int M ) const
{
	int logM=1;
	while( (1<<logM)<M ) logM++;
	if( (1<<logM)!=M || m<0 || m>=M ) return false;
	memset( thetaStencil , 0 , sizeof( thetaStencil ) );
	if( logM>_logHeight || logM<2 ) EquiRectangular< Real >::ThetaStencil( m , M , thetaStencil );
	else
	{
		bool flip = false;
		if( m>=M/2 ) m = M-1-m , flip = true;
		for( int d= 0 ; d<=2 ; d++ ) thetaStencil[2+d] = _valueTable[logM][m].thetaValues[d];
		for( int d=-2 ; d< 0 ; d++ ) if( m+d>=0 ) thetaStencil[2+d] = _valueTable[logM][m+d].thetaValues[-d];
		if( flip )
		{
			Real temp;
			for( int d=-2 ; d<0 ; d++ ) temp = thetaStencil[2+d] , thetaStencil[2+d] = thetaStencil[2-d] , thetaStencil[2-d] = temp;
		}
	}
	return true;
}
template< class Real >
bool ThetaStencilTable< Real >::setDivergenceDThetaStencil( Real divergenceDThetaStencil[4] , int m , int M ) const
{
	int logM=1;
	while( (1<<logM)<M ) logM++;
	if( (1<<logM)!=M || m<0 || m>=M ) return false;
	memset( divergenceDThetaStencil , 0 , sizeof( divergenceDThetaStencil ) );
	if( logM>_logHeight || logM<2 || logM>=_valueTable.size() ) EquiRectangular< Real >::DivergenceDThetaStencil( m , M , divergenceDThetaStencil );
	else
	{
		bool flip = false;
		if( m>=M/2 ) m = M-1-m , flip = true;
		for( int d=0 ; d<4 ; d++ ) divergenceDThetaStencil[d] = _valueTable[logM][m].divergenceDThetaValues[d];
		if( flip )
		{
			Real temp;
			for( int d=-2 ; d<0 ; d++ ) temp = divergenceDThetaStencil[2+d] , divergenceDThetaStencil[2+d] = -divergenceDThetaStencil[1-d] , divergenceDThetaStencil[1-d] = -temp;
		}
	}
	return true;
}
template< class Real >
bool ThetaStencilTable< Real >::setLaplacianThetaStencil( Real laplacianThetaStencil[5] , int m , int M , int samples ) const
{
	int logM=1;
	while( (1<<logM)<M ) logM++;
	if( (1<<logM)!=M || m<0 || m>=M ) return false;

	memset( laplacianThetaStencil , 0 , sizeof( laplacianThetaStencil ) );
	if( logM>_logHeight || logM<2 || logM>=_valueTable.size() ) EquiRectangular< Real >::LaplacianThetaStencil( m , M , laplacianThetaStencil , samples );
	else
	{
		bool flip = false;
		if( m>=M/2 ) m = M-1-m , flip = true;
		for( int d=0 ; d<5 ; d++ ) laplacianThetaStencil[d] = _valueTable[logM][m].laplacianThetaValues[d];
		if( flip )
		{
			Real temp;
			for( int d=-2 ; d<0 ; d++ ) temp = laplacianThetaStencil[2+d] , laplacianThetaStencil[2+d] = laplacianThetaStencil[2-d] , laplacianThetaStencil[2-d] = temp;
		}
	}
	return true;
}
template< class Real >
bool ThetaStencilTable< Real >::setLaplacianD2ThetaStencil( Real laplacianD2ThetaStencil[5] , int m , int M ) const
{
	int logM=1;
	while( (1<<logM)<M ) logM++;
	if( (1<<logM)!=M || m<0 || m>=M ) return false;

	memset( laplacianD2ThetaStencil , 0 , sizeof( laplacianD2ThetaStencil ) );
	if( logM>_logHeight || logM<2 || logM>=_valueTable.size() ) EquiRectangular< Real >::LaplacianD2ThetaStencil( m , M , laplacianD2ThetaStencil );
	else
	{
		bool flip = false;
		if( m>=M/2 ) m = M-1-m , flip = true;
		for( int d= 0 ; d<=2 ; d++ ) laplacianD2ThetaStencil[2+d] = _valueTable[logM][m].laplacianD2ThetaValues[d];
		for( int d=-2 ; d< 0 ; d++ ) if( m+d>=0 ) laplacianD2ThetaStencil[2+d] = _valueTable[logM][m+d].laplacianD2ThetaValues[-d];
		if( flip )
		{
			Real temp;
			for( int d=-2 ; d<0 ; d++ ) temp = laplacianD2ThetaStencil[2+d] , laplacianD2ThetaStencil[2+d] =  laplacianD2ThetaStencil[2-d] , laplacianD2ThetaStencil[2-d] =  temp;
		}
	}
	return true;
}
template< class Real >
bool ThetaStencilTable< Real >::Init( int height , int samples , bool highPrecision , bool useRestriction )
{
	int logHeight=1;
	while( (1<<logHeight)<height ) logHeight++;
	if( (1<<logHeight)!=height || height<4 ) return false;

	if( logHeight<=_logHeight) return true;
	_logHeight = logHeight;
	_clearTable();
	_valueTable.resize( logHeight+1 );
	for( int l=0 ; l<_valueTable.size() ; l++ ) _valueTable[l] = NULL;
	int l = logHeight;
	{
		int _height = 1<<l;
		_valueTable[l] = new StencilValues[ _height/2 ];
		for( int i=0 ; i<_height/2 ; i++ )
		{
			if( highPrecision )
			{
				_valueTable[l][i].thetaWeight = EquiRectangular< Real , qfloat >::ThetaWeight( i , _height );
				EquiRectangular< Real , qfloat >::HalfThetaStencil           ( i , _height , _valueTable[l][i].thetaValues );
				EquiRectangular< Real , qfloat >::DivergenceDThetaStencil    ( i , _height , _valueTable[l][i].divergenceDThetaValues );
				EquiRectangular< Real , qfloat >::LaplacianThetaStencil      ( i , _height , _valueTable[l][i].laplacianThetaValues , samples );
				EquiRectangular< Real , qfloat >::HalfLaplacianD2ThetaStencil( i , _height , _valueTable[l][i].laplacianD2ThetaValues );
			}
			else
			{
				_valueTable[l][i].thetaWeight = EquiRectangular< Real >::ThetaWeight( i , _height );
				EquiRectangular< Real >::HalfThetaStencil           ( i , _height , _valueTable[l][i].thetaValues );
				EquiRectangular< Real >::DivergenceDThetaStencil    ( i , _height , _valueTable[l][i].divergenceDThetaValues );
				EquiRectangular< Real >::LaplacianThetaStencil      ( i , _height , _valueTable[l][i].laplacianThetaValues , samples );
				EquiRectangular< Real >::HalfLaplacianD2ThetaStencil( i , _height , _valueTable[l][i].laplacianD2ThetaValues );
			}
		}
	}
	if( !useRestriction )
	{
		for( int l=logHeight-1 ; l>=2 ; l-- )
		{
			int _height = 1<<l;
			_valueTable[l] = new StencilValues[ _height/2 ];
			for( int i=0 ; i<_height/2 ; i++ )
			{
				if( highPrecision )
				{
					_valueTable[l][i].thetaWeight = EquiRectangular< Real , qfloat >::ThetaWeight( i , _height );
					EquiRectangular< Real , qfloat >::HalfThetaStencil           ( i , _height , _valueTable[l][i].thetaValues );
					EquiRectangular< Real , qfloat >::DivergenceDThetaStencil    ( i , _height , _valueTable[l][i].divergenceDThetaValues );
					EquiRectangular< Real , qfloat >::LaplacianThetaStencil      ( i , _height , _valueTable[l][i].laplacianThetaValues , samples );
					EquiRectangular< Real , qfloat >::HalfLaplacianD2ThetaStencil( i , _height , _valueTable[l][i].laplacianD2ThetaValues );
				}
				else
				{
					_valueTable[l][i].thetaWeight = EquiRectangular< Real >::ThetaWeight( i , _height );
					EquiRectangular< Real >::HalfThetaStencil           ( i , _height , _valueTable[l][i].thetaValues );
					EquiRectangular< Real >::DivergenceDThetaStencil    ( i , _height , _valueTable[l][i].divergenceDThetaValues );
					EquiRectangular< Real >::LaplacianThetaStencil      ( i , _height , _valueTable[l][i].laplacianThetaValues , samples );
					EquiRectangular< Real >::HalfLaplacianD2ThetaStencil( i , _height , _valueTable[l][i].laplacianD2ThetaValues );
				}
			}
		}
		return true;
	}
	// Set the coarser stencils by restriction
	for( int l=logHeight-1 ; l>=2 ; l-- )
	{
#if 1
		double pStencil1[] = { 0.25 , 0.50 , 0.25 };
#else
		double pStencil1[] = { 0.50 , 1.00 , 0.50 };
#endif
		double pStencil2[] = { 0.25 , 0.75 , 0.75 , 0.25 };
		double weight1 , weight2;
		Real thetaWeights[4] , thetaStencils[4][5] , divergenceDThetaStencils[4][4] , laplacianThetaStencils[4][5] , laplacianD2ThetaStencils[4][5];
		int _height = 1<<l;
		_valueTable[l] = new StencilValues[ _height/2 ];
		for( int x=0 ; x<_height/2 ; x++ )
		{
			for( int i=0 ; i<4 ; i++ )
			{
				thetaWeights[i] =                        thetaWeight    ( 2*x-1+i , _height*2 );
				setThetaStencil           (            thetaStencils[i] , 2*x-1+i , _height*2 );
				setDivergenceDThetaStencil( divergenceDThetaStencils[i] , 2*x-1+i , _height*2 );
				setLaplacianThetaStencil  (   laplacianThetaStencils[i] , 2*x-1+i , _height*2 );
				setLaplacianD2ThetaStencil( laplacianD2ThetaStencils[i] , 2*x-1+i , _height*2 );
			}
			if( !x ) _valueTable[l][x].thetaWeight = thetaWeights[1] * ( pStencil2[0] + pStencil2[1] ) + thetaWeights[2] * pStencil2[2] + thetaWeights[3] * pStencil2[3];
			else     _valueTable[l][x].thetaWeight = thetaWeights[0] * pStencil2[0] + thetaWeights[1] * pStencil2[1] + thetaWeights[2] * pStencil2[2] + thetaWeights[3] * pStencil2[3];


			for( int i=0 ; i<=2 ; i++ )
			{
				_valueTable[l][x].thetaValues[i] = _valueTable[l][x].laplacianThetaValues[i] = _valueTable[l][x].laplacianD2ThetaValues[i] = 0;
				int y = x+i;
				for( int j=-1 ; j<=2 ; j++ )
				{
					int jj = 2*x+j;
					if( jj<0 ) continue;
					if( !jj ) weight1 = pStencil2[0] + pStencil2[1];
					else      weight1 = pStencil2[1+j];
					for( int k=-1 ; k<=2 ; k++ )
					{
						int kk = 2*y+k;
						if( kk<0 ) continue;
						if( !kk ) weight2 = pStencil2[0] + pStencil2[1];
						else      weight2 = pStencil2[1+k];
						if( kk-jj>=-2 && kk-jj<=2 )
						{
							_valueTable[l][x].thetaValues           [i] +=            thetaStencils[j+1][2+kk-jj] * weight1 * weight2;
							_valueTable[l][x].laplacianD2ThetaValues[i] += laplacianD2ThetaStencils[j+1][2+kk-jj] * weight1 * weight2;
						}
					}
				}
			}
			for( int i=-2 ; i<=2 ; i++ )
			{
				_valueTable[l][x].laplacianThetaValues[i+2] = 0;
				int y = x+i;
				if( y<0 ) continue;
				for( int j=-1 ; j<=2 ; j++ )
				{
					int jj = 2*x+j;
					if( jj<0 ) continue;
					if( !jj ) weight1 = pStencil2[0] + pStencil2[1];
					else      weight1 = pStencil2[1+j];
					for( int k=-1 ; k<=2 ; k++ )
					{
						int kk = 2*y+k;
						if( kk<0 ) continue;
						if( !kk ) weight2 = pStencil2[0] + pStencil2[1];
						else      weight2 = pStencil2[1+k];
						if( kk-jj>=-2 && kk-jj<=2 )
						{
							_valueTable[l][x].laplacianThetaValues[i+2] += laplacianThetaStencils[j+1][2+kk-jj] * weight1 * weight2;
						}
					}
				}
			}
			for( int i=0 ; i<4 ; i++ )
			{
				_valueTable[l][x].divergenceDThetaValues[i] = 0;
				int dy = x + i - 2;
				for( int j=-1 ; j<=2; j++ )
				{
					int jj = 2*x + j;
					if( jj<0 ) continue;
					if( !jj ) weight1 = pStencil2[0] + pStencil2[1];
					else      weight1 = pStencil2[1+j];
					for( int k=-1 ; k<=1 ; k++ )
					{
						int kk = 2*dy + k + 1;
						if( kk<0 ) continue;
						weight2 = pStencil1[1+k];
						if( kk-jj>=-2 && kk-jj<=1 ) _valueTable[l][x].divergenceDThetaValues[i] += divergenceDThetaStencils[j+1][2+kk-jj] * weight1 * weight2;
					}

				}
			}
		}
	}
	return true;
}
template< class Real >
bool ThetaStencilTable< Real >::read( FILE* fp )
{
	_clearTable( );
	fread( &_logHeight , sizeof( int ) , 1 , fp );
	_valueTable.resize( _logHeight+1 );

	for( int l=_logHeight ; l>=2 ; l-- )
	{
		int _height = 1<<l;
		_valueTable[l] = new StencilValues[ _height/2 ];
		fread( &_valueTable[l][0] , sizeof( StencilValues ) , _height/2 , fp );
	}
	return true;
}
template< class Real >
bool ThetaStencilTable< Real >::write( FILE* fp ) const
{
	fwrite( &_logHeight , sizeof( int ) , 1 , fp );
	for( int l=_logHeight ; l>=2 ; l-- )
	{
		int _height = 1<<l;
		fwrite( &_valueTable[l][0] , sizeof( StencilValues ) , _height/2 , fp );
	}
	return true;
}
template< class Real >
void ThetaStencilTable< Real >::_clearTable( void )
{
	for( int l=0 ; l<_valueTable.size() ; l++ ) if( _valueTable[l] ) delete[] _valueTable[l];
	_valueTable.resize( 0 );
}

/////////////////////
// PhiStencilTable //
/////////////////////
template< class Real > PhiStencilTable< Real >::PhiStencilTable( void ) { _logWidth = 0; }
template< class Real > PhiStencilTable< Real >::~PhiStencilTable( void ){ _logWidth = 0; }
template< class Real >
Real PhiStencilTable< Real >::phiWeight( int N ) const
{
	int logN=1;
	while( (1<<logN)<N ) logN++;
	if( (1<<logN)!=N || logN>_logWidth || logN<2 || logN>=_valueTable.size() ) return EquiRectangular< Real >::PhiWeight( N );
	else return _valueTable[logN].phiWeight;
}
template< class Real >
bool PhiStencilTable< Real >::setPhiStencil( Real phiStencil[5] , int N ) const
{
	int logN=1;
	while( (1<<logN)<N ) logN++;
	memset( phiStencil , 0 , sizeof( phiStencil ) );
	if( (1<<logN)!=N || logN>_logWidth || logN<2 || logN>=_valueTable.size() ) EquiRectangular< Real >::PhiStencil( N , phiStencil );
	else for( int d= 0 ; d<=2 ; d++ ) phiStencil[2-d] = phiStencil[2+d] = _valueTable[logN].phiValues[d];
	return true;
}
template< class Real >
bool PhiStencilTable< Real >::setDivergenceDPhiStencil( Real divergenceDPhiStencil[4] , int N ) const
{
	int logN=1;
	while( (1<<logN)<N ) logN++;
	memset( divergenceDPhiStencil , 0 , sizeof( divergenceDPhiStencil ) );
	if( (1<<logN)!=N || logN>_logWidth || logN<2 || logN>=_valueTable.size() ) EquiRectangular< Real >::DivergenceDPhiStencil( N , divergenceDPhiStencil );
	else for( int d=0 ; d<4 ; d++ ) divergenceDPhiStencil[d] = _valueTable[logN].divergenceDPhiValues[d];
	return true;
}
template< class Real >
bool PhiStencilTable< Real >::setLaplacianD2PhiStencil( Real laplacianD2PhiStencil[5] , int N ) const
{
	int logN=1;
	while( (1<<logN)<N ) logN++;
	memset( laplacianD2PhiStencil , 0 , sizeof( laplacianD2PhiStencil ) );
	if( (1<<logN)!=N || logN>_logWidth || logN<2 || logN>=_valueTable.size() ) EquiRectangular< Real >::LaplacianD2PhiStencil( N , laplacianD2PhiStencil );
	else for( int d= 0 ; d<=2 ; d++ ) laplacianD2PhiStencil[2-d] = laplacianD2PhiStencil[2+d] = _valueTable[logN].laplacianD2PhiValues[d];
	return true;
}
template< class Real >
bool PhiStencilTable< Real >::Init( int width , bool highPrecision , bool useRestriction )
{
	int logWidth=1;
	while( (1<<logWidth)<width ) logWidth++;
	if( (1<<logWidth)!=width || width<4 ) return false;
	if( logWidth<=_logWidth ) return true;
	_logWidth = logWidth;

	_clearTable();
	_valueTable.resize( logWidth+1 );
	int l = logWidth;
	{
		int _width = 1<<l;
		if( highPrecision )
		{
			_valueTable[l].phiWeight = EquiRectangular< Real , qfloat >::PhiWeight( _width );
			EquiRectangular< Real , qfloat >::HalfPhiStencil           ( _width , _valueTable[l].phiValues );
			EquiRectangular< Real , qfloat >::DivergenceDPhiStencil    ( _width , _valueTable[l].divergenceDPhiValues );
			EquiRectangular< Real , qfloat >::HalfLaplacianD2PhiStencil( _width , _valueTable[l].laplacianD2PhiValues );
		}
		else
		{
			_valueTable[l].phiWeight = EquiRectangular< Real >::PhiWeight( _width );
			EquiRectangular< Real >::HalfPhiStencil           ( _width , _valueTable[l].phiValues );
			EquiRectangular< Real >::DivergenceDPhiStencil    ( _width , _valueTable[l].divergenceDPhiValues );
			EquiRectangular< Real >::HalfLaplacianD2PhiStencil( _width , _valueTable[l].laplacianD2PhiValues );
		}
	}
	if( !useRestriction )
	{
		for( int l=logWidth-1 ; l>=2 ; l-- )
		{
			int _width = 1<<l;
			if( highPrecision )
			{
				_valueTable[l].phiWeight = EquiRectangular< Real , qfloat >::PhiWeight( _width );
				EquiRectangular< Real , qfloat >::HalfPhiStencil           ( _width , _valueTable[l].phiValues );
				EquiRectangular< Real , qfloat >::DivergenceDPhiStencil    ( _width , _valueTable[l].divergenceDPhiValues );
				EquiRectangular< Real , qfloat >::HalfLaplacianD2PhiStencil( _width , _valueTable[l].laplacianD2PhiValues );
			}
			else
			{
				_valueTable[l].phiWeight = EquiRectangular< Real >::PhiWeight( _width );
				EquiRectangular< Real >::HalfPhiStencil           ( _width , _valueTable[l].phiValues );
				EquiRectangular< Real >::DivergenceDPhiStencil    ( _width , _valueTable[l].divergenceDPhiValues );
				EquiRectangular< Real >::HalfLaplacianD2PhiStencil( _width , _valueTable[l].laplacianD2PhiValues );
			}
		}
		return true;
	}
	// Set the coarser stencils by restriction
	for( int l=logWidth-1 ; l>=2 ; l-- )
	{
#if 1
		double pStencil1[] = { 0.25 , 0.50 , 0.25 };
#else
		double pStencil1[] = { 0.50 , 1.00 , 0.50 };
#endif
		double pStencil2[] = { 0.25 , 0.75 , 0.75 , 0.25 };
		Real phiStencil[5] , divergenceDPhiStencil[4] , laplacianD2PhiStencil[5];
		int _width = 1<<l;
		setPhiStencil           (            phiStencil , _width*2 );
		setDivergenceDPhiStencil( divergenceDPhiStencil , _width*2 );
		setLaplacianD2PhiStencil( laplacianD2PhiStencil , _width*2 );
		_valueTable[l].phiWeight = ( pStencil2[0] + pStencil2[1] + pStencil2[2] + pStencil2[3] ) * phiWeight( _width*2 );
		for( int y=0 ; y<=2 ; y++ )
		{
			_valueTable[l].phiValues[y] = _valueTable[l].laplacianD2PhiValues[y] = 0;
			for( int j=-1 ; j<=2 ; j++ )
			{
				int jj = j;
				for( int k=-1 ; k<=2 ; k++ )
				{
					int kk = 2*y+k;
					if( kk-jj>=-2 && kk-jj<=2 )
					{
						_valueTable[l].phiValues           [y] +=            phiStencil[2+kk-jj] * pStencil2[1+j] * pStencil2[1+k];
						_valueTable[l].laplacianD2PhiValues[y] += laplacianD2PhiStencil[2+kk-jj] * pStencil2[1+j] * pStencil2[1+k];
					}
				}
			}
			for( int i=0 ; i<4 ; i++ )
			{
				_valueTable[l].divergenceDPhiValues[i] = 0;
				int dy = i - 2;
				for( int jj=-1 ; jj<=2; jj++ )
					for( int k=-1 ; k<=1 ; k++ )
					{
						int kk = 2*dy + k + 1;
						if( kk-jj>=-2 && kk-jj<=1 ) _valueTable[l].divergenceDPhiValues[i] += divergenceDPhiStencil[2+kk-jj] * pStencil2[1+jj] * pStencil1[1+k];
					}
			}
		}
	}
	return true;
}
template< class Real >
bool PhiStencilTable< Real >::read( FILE* fp )
{
	_clearTable( );
	fread( &_logWidth , sizeof( int ) , 1 , fp );
	_valueTable.resize( _logWidth+1 );
	fread( &_valueTable[0] , sizeof( StencilValues ) , _logWidth+1 , fp );
	return true;
}
template< class Real >
bool PhiStencilTable< Real >::write( FILE* fp ) const
{
	fwrite( &_logWidth , sizeof( int ) , 1 , fp );
	fwrite( &_valueTable[0] , sizeof( StencilValues ) , _valueTable.size() , fp );
	return true;
}
template< class Real >
void PhiStencilTable< Real >::_clearTable( void )
{
	_valueTable.resize( 0 );
}

/////////////////////
// EquiRectangular //
/////////////////////
template< class Real , class PolyReal >
void EquiRectangular< Real , PolyReal >::SetBaseFunctions( int m , int n , int M , int N , PPolynomial< 2 , PolyReal >& thetaBase , PPolynomial< 2 , PolyReal >& phiBase , bool rescale )
{
	PolyReal thetaScale =     M_PI / M;
	PolyReal phiScale   = 2.0*M_PI / N;

	if( m==0 || m==M-1 )	phiBase = PPolynomial< 2 , PolyReal >::ConstantFunction( PolyReal(N)/PolyReal(2) );
	else					phiBase = PPolynomial< 2 , PolyReal >::GaussianApproximation( ).shift( n+0.5 );

	thetaBase = PPolynomial< 2 , PolyReal >::GaussianApproximation();
	if		( m==  0 )	thetaBase =	thetaBase.shift(   0.5 ) + thetaBase.shift(  -0.5 );
	else if	( m==M-1 )	thetaBase = thetaBase.shift( M+0.5 ) + thetaBase.shift( M-0.5 );
	else				thetaBase = thetaBase.shift( m+0.5 );

	if( rescale ) thetaBase = thetaBase.scale( thetaScale ) , phiBase = phiBase.scale( phiScale );
}

template< class Real , class PolyReal >
void EquiRectangular< Real , PolyReal >::SetThetaBaseFunction( int m , int M , PPolynomial< 2 , PolyReal >& thetaBase , bool rescale )
{
	PolyReal thetaScale = M_PI / M;
	thetaBase = PPolynomial< 2 , PolyReal >::GaussianApproximation();
	if		( m==  0 )	thetaBase =	thetaBase.shift(   0.5 ) + thetaBase.shift(  -0.5 );
	else if	( m==M-1 )	thetaBase = thetaBase.shift( M+0.5 ) + thetaBase.shift( M-0.5 );
	else				thetaBase = thetaBase.shift( m+0.5 );
	if( rescale ) thetaBase = thetaBase.scale( thetaScale );
}
template< class Real , class PolyReal >
void EquiRectangular< Real , PolyReal >::SetPhiBaseFunction( int n , int N , PPolynomial< 2 , PolyReal >& phiBase , bool rescale )
{
	PolyReal phiScale   = 2.0*M_PI / N;
	phiBase = PPolynomial< 2 , PolyReal >::GaussianApproximation().shift( n+0.5 );
	if( rescale ) phiBase = phiBase.scale( phiScale );
}

template< class Real , class PolyReal >
void EquiRectangular< Real , PolyReal >::SetDThetaBaseFunction( int m , int M , PPolynomial< 1 , PolyReal >& thetaBase , bool rescale )
{
	PolyReal thetaScale = M_PI / M;
	thetaBase = PPolynomial< 1 , PolyReal >::GaussianApproximation().shift( m+1. );
	if( rescale ) thetaBase = thetaBase.scale( thetaScale ) / thetaScale;
}
template< class Real , class PolyReal >
void EquiRectangular< Real , PolyReal >::SetDPhiBaseFunction( int n , int N , PPolynomial< 1 , PolyReal >& phiBase , bool rescale )
{
	PolyReal phiScale = 2.0*M_PI / N;
	phiBase = PPolynomial< 1 , PolyReal >::GaussianApproximation().shift( n+1. );
	if( rescale ) phiBase = phiBase.scale( phiScale ) / phiScale;
}

template< class Real , class PolyReal >
double EquiRectangular< Real , PolyReal >::PhiWeight( int N )
{
	return 2. * M_PI / N;
}
template< class Real , class PolyReal >
void EquiRectangular< Real , PolyReal >::PhiStencil( int N , double stencil[5] )
{
	PolyReal b = 2.0 * M_PI / N;
	PPolynomial< 2 , PolyReal > phiBase1 , phiBase2;
	if( N==1 )
	{
		for( int i=-2 ; i<=2 ; i++ ) stencil[i+2] = 2. * M_PI;
		return;
	}

	SetPhiBaseFunction( 0 , N , phiBase1 , false );
	for( int i=-2 ; i<=2 ; i++ )
	{
		SetPhiBaseFunction( i , N , phiBase2 , false );
		stencil[i+2] = ( phiBase1 * phiBase2 ).integral( -1.0 , 2.0 ) * b;
	}
}
template< class Real , class PolyReal >
void EquiRectangular< Real , PolyReal >::HalfPhiStencil( int N , double stencil[3] )
{
	PolyReal b = 2.0 * M_PI / N;
	PPolynomial< 2 , PolyReal > phiBase1 , phiBase2;

	SetPhiBaseFunction( 0 , N , phiBase1 , false );
	for( int i=0 ; i<=2 ; i++ )
	{
		SetPhiBaseFunction( i , N , phiBase2 , false );
		stencil[i] = ( phiBase1 * phiBase2 ).integral( -1.0 , 2.0 ) * b;
	}
}
template< class Real , class PolyReal >
bool EquiRectangular< Real , PolyReal >::DivergenceDPhiStencil( int N , double stencil[4] )
{
	for( int i=0 ; i<4 ; i++ ) stencil[i] = 0;
	if( N<=3 ) return false;

	PolyReal b = 2.0 * M_PI / N;
	PPolynomial< 2 , PolyReal > phiBase;
	PPolynomial< 1 , PolyReal > dPhiBase1 , dPhiBase2;
	SetPhiBaseFunction( 0 , N , phiBase , false );
	dPhiBase1 = phiBase.derivative( );
	for( int j=-2 ; j<=1 ; j++ )
	{
		SetDPhiBaseFunction( j , N , dPhiBase2 , false );
		stencil[j+2] = ( dPhiBase1 * dPhiBase2 ).integral( -1.0 , 2.0 ) / b;
	}
	return true;
}
template< class Real , class PolyReal >
void EquiRectangular< Real , PolyReal >::LaplacianD2PhiStencil( int N , double stencil[5] )
{
	PolyReal b = 2.0 * M_PI / N;
	PPolynomial< 2 , PolyReal > phiBase;
	PPolynomial< 1 , PolyReal > dPhiBase1 , dPhiBase2;

	SetPhiBaseFunction( 0 , N , phiBase , false );
	dPhiBase1 = phiBase.derivative( );
	for( int i=-2 ; i<=2 ; i++ )
	{
		SetPhiBaseFunction( i , N , phiBase , false );
		dPhiBase2 = phiBase.derivative( );
		stencil[i+2] = ( dPhiBase1 * dPhiBase2 ).integral( -1.0 , 2.0 ) / b;
	}
}
template< class Real , class PolyReal >
void EquiRectangular< Real , PolyReal >::HalfLaplacianD2PhiStencil( int N , double stencil[3] )
{
	PolyReal b = 2.0 * M_PI / N;
	PPolynomial< 2 , PolyReal > phiBase;
	PPolynomial< 1 , PolyReal > dPhiBase1 , dPhiBase2;

	SetPhiBaseFunction( 0 , N , phiBase , false );
	dPhiBase1 = phiBase.derivative( );
	for( int i=0 ; i<=2 ; i++ )
	{
		SetPhiBaseFunction( i , N , phiBase , false );
		dPhiBase2 = phiBase.derivative( );
		stencil[i] = ( dPhiBase1 * dPhiBase2 ).integral( -1.0 , 2.0 ) / b;
	}
}
template< class Real , class PolyReal >
double EquiRectangular< Real , PolyReal >::ThetaWeight( int m , int M )
{
	PolyReal a = M_PI / M;
	int start = m-1 , end = m+2;
	PPolynomial< 2 , PolyReal > thetaBase;
	if( start<0 ) start = 0;
	if(   end>M ) end   = M;

	SetThetaBaseFunction( m , M , thetaBase , false );
	return thetaBase.integralSmoothSine( 2 , a , start , end ) * a;
}
template< class Real , class PolyReal >
void EquiRectangular< Real , PolyReal >::ThetaStencil( int m , int M , double stencil[5] )
{
	PolyReal a = M_PI / M;
	int start = m-1 , end = m+2;
	PPolynomial< 2 , PolyReal > thetaBase1 , thetaBase2;
	if( start<0 ) start = 0;
	if(   end>M ) end   = M;

	SetThetaBaseFunction( m , M , thetaBase1 , false );
	for( int i=-2 ; i<=2 ; i++ )
	{
		if( m+i<0 || m+i>=M ) stencil[i+2] = 0;
		else
		{
			int s = start>m+i-1 ? start : m+i-1;
			int e = end  <m+i+2 ? end   : m+i+2;
			SetThetaBaseFunction( m+i , M , thetaBase2 , false );
			stencil[i+2] = ( thetaBase1 * thetaBase2 ).integralSmoothSine( 2 , a , s , e ) * a;
		}
	}
}
template< class Real , class PolyReal >
void EquiRectangular< Real , PolyReal >::HalfThetaStencil( int m , int M , double stencil[3] )
{
	PolyReal a = M_PI / M;
	int start = m-1 , end = m+2;
	PPolynomial< 2 , PolyReal > thetaBase1 , thetaBase2;
	if( start<0 ) start = 0;
	if(   end>M ) end   = M;

	SetThetaBaseFunction( m , M , thetaBase1 , false );
	for( int i=0 ; i<=2 ; i++ )
	{
		if( m+i<0 || m+i>=M ) stencil[i] = 0;
		else
		{
			int s = start>m+i-1 ? start : m+i-1;
			int e = end  <m+i+2 ? end   : m+i+2;
			SetThetaBaseFunction( m+i , M , thetaBase2 , false );
			stencil[i] = ( thetaBase1 * thetaBase2 ).integralSmoothSine( 2 , a , s , e ) * a;
		}
	}
}
template< class Real , class PolyReal >
bool EquiRectangular< Real , PolyReal >::DivergenceDThetaStencil( int m , int M , double stencil[4] )
{
	for( int i=0 ; i<4 ; i++ ) stencil[i] = 0;
	if( m<0 || m>M-1 || M<6 ) return false;

	PolyReal a = M_PI / M;
	int start=m-1 , end=m+2;
	PPolynomial< 2 , PolyReal > thetaBase;
	PPolynomial< 1 , PolyReal > dThetaBase1 , dThetaBase2;

	SetThetaBaseFunction( m , M , thetaBase , false );
	dThetaBase1 = thetaBase.derivative( );
	for( int i=-2 ; i<=1 ; i++ )
	{
		if( i+m<0 || i+m>=M-1 ) continue;
		int s = start>m+i   ? start : m+i  ;
		int e = end  <m+i+2 ? end   : m+i+2;
		SetDThetaBaseFunction( m+i , M , dThetaBase2 , false );
		stencil[i+2] = ( dThetaBase1*dThetaBase2).integralSmoothSine( 1 , a , s , e ) / a;
	}
	return true;
}
template< class Real , class PolyReal >
void EquiRectangular< Real , PolyReal >::LaplacianD2ThetaStencil( int m , int M , double stencil[5] )
{
	PolyReal a = M_PI / M;
	int start = m-1 , end = m+2;
	PPolynomial< 2 , PolyReal > thetaBase;
	PPolynomial< 1 , PolyReal > dThetaBase1 , dThetaBase2;
	if( start<0 ) start = 0;
	if(   end>M ) end   = M;

	SetThetaBaseFunction( m , M , thetaBase , false );
	dThetaBase1 = thetaBase.derivative( );
	for( int i=-2 ; i<=2 ; i++ )
	{
		if( m+i<0 || m+i>=M ) { stencil[i+2]=0 ; continue; }
		int s = start>m+i-1 ? start : m+i-1;
		int e = end  <m+i+2 ? end   : m+i+2;
		SetThetaBaseFunction( m+i , M , thetaBase , false );
		dThetaBase2 = thetaBase.derivative( );
		stencil[i+2] = ( dThetaBase1 * dThetaBase2 ).integralSmoothSine( 1 , a , s , e ) / a;
	}
}
template< class Real , class PolyReal >
void EquiRectangular< Real , PolyReal >::HalfLaplacianD2ThetaStencil( int m , int M , double stencil[3] )
{
	PolyReal a = M_PI / M;
	int start = m-1 , end = m+2;
	PPolynomial< 2 , PolyReal > thetaBase;
	PPolynomial< 1 , PolyReal > dThetaBase1 , dThetaBase2;
	if( start<0 ) start = 0;
	if(   end>M ) end   = M;

	SetThetaBaseFunction( m , M , thetaBase , false );
	dThetaBase1 = thetaBase.derivative( );
	for( int i=0 ; i<=2 ; i++ )
	{
		int s = start>m+i-1 ? start : m+i-1;
		int e = end  <m+i+2 ? end   : m+i+2;
		SetThetaBaseFunction( m+i , M , thetaBase , false );
		dThetaBase2 = thetaBase.derivative( );
		stencil[i] = ( dThetaBase1 * dThetaBase2 ).integralSmoothSine( 1 , a , s , e ) / a;
	}
}
template< class Real , class PolyReal >
void EquiRectangular< Real , PolyReal >::LaplacianThetaStencil( int m , int M , double stencil[5] , int samples )
{
	memset( stencil , 0 , sizeof( stencil ) );

	bool flip;
	if( m>=M/2 )	flip = true , m = M-1-m;
	else			flip = false;
	PolyReal a = M_PI / M;
	int start = m-1 , end = m+2;
	PPolynomial< 2 , PolyReal > thetaBase1 , thetaBase2;
	SetThetaBaseFunction( m , M , thetaBase1 , false );
	if( start<0 ) start = 0;
	if(   end>M ) end   = M;
	PPolynomial< 6 , PolyReal > thetaBase = thetaBase1 * Polynomial< 4 , PolyReal >::Cosecant( PolyReal(start+end)/PolyReal(2.) * a ).scale( PolyReal( 1. ) / a );
	for( int i=-2 ; i<=2 ; i++ )
	{
		if(  i+m<=0 || i+m>=M-1 || m<=0 || m>=M-1 ) continue;
		SetThetaBaseFunction( i+m , M , thetaBase2 , false );
		int s = start>m+i-1 ? start : m+i-1;
		int e = end  <m+i+2 ? end   : m+i+2;
		if( m==1 || i+m==1 )
		{
	// WARNING: It seems that at high resolutions, using a 4th approximation to the Cosecant function leads to instability, so we
			if( samples )
			{
				PPolynomial< 4 , PolyReal > poly = thetaBase1*thetaBase2;
				stencil[i+2] =
					(
					poly.integralCosecant( a , 0 , 1 , samples ) +
					poly.integralCosecant( a , 1 , 2 , samples ) +
					poly.integralCosecant( a , 2 , 3 , samples )
					) * a;
			}
			else
			{
				stencil[i+2]  = Real( ( thetaBase1 * thetaBase2 * Polynomial< 4 , PolyReal >::Cosecant( PolyReal( 1. )                    * a ).scale( PolyReal( 1. ) / a ) ).integral( 0. , 1.  ) * a );
				stencil[i+2] += Real( ( thetaBase1 * thetaBase2 * Polynomial< 4 , PolyReal >::Cosecant( PolyReal( 1. + end )/PolyReal(2.) * a ).scale( PolyReal( 1. ) / a ) ).integral( 1. , end ) * a );
			}
		}
		else
		if( samples )
		{
			PPolynomial< 4 , PolyReal > poly = thetaBase1*thetaBase2;
			stencil[i+2] =
				(
				poly.integralCosecant( a , m-1 , m   , samples ) +
				poly.integralCosecant( a , m   , m+1 , samples ) +
				poly.integralCosecant( a , m+1 , m+2 , samples ) 
				) * a;
		}
		else
		{
			stencil[i+2] = ( thetaBase * thetaBase2 ).integral( start , end ) * a;
		}
	}
	if( flip )
	{
		Real temp;
		temp = stencil[4] , stencil[4] = stencil[0] , stencil[0] = temp;
		temp = stencil[3] ,	stencil[3] = stencil[1] , stencil[1] = temp;
	}
}
template< class Real , class PolyReal >
void EquiRectangular< Real , PolyReal >::HalfLaplacianThetaStencil( int m , int M , double stencil[3] , int samples )
{
	memset( stencil , 0 , sizeof( stencil ) );
	PolyReal a = M_PI / M;
	// NOTE: Assuming that m<M/2, so no flipping required!
	if( m>=M/2 )
	{
		fprintf( stderr , "Error: HalfLaplacianThetaStencil assumes m<M/2: %d %d/%d\n" , m , M );
		exit( 0 );
	}
	int start = m-1 , end = m+2;
	PPolynomial< 2 , PolyReal > thetaBase1 , thetaBase2;
	SetThetaBaseFunction( m , M , thetaBase1 , false );
	if( start<0 ) start = 0;
	if(   end>M ) end   = M;
	PPolynomial< 6 , PolyReal > thetaBase = thetaBase1 * Polynomial< 4 , PolyReal >::Cosecant( PolyReal(start+end)/PolyReal(2.) * a ).scale( PolyReal( 1. ) / a );
	for( int i=0 ; i<=2 ; i++ )
	{
		if(  i+m<=0 || i+m>=M-1 || m<=0 || m>=M-1 ) continue;
		SetThetaBaseFunction( i+m , M , thetaBase2 , false );
		int s = start>m+i-1 ? start : m+i-1;
		int e = end  <m+i+2 ? end   : m+i+2;
		if( m==1 || i+m==1 )
		{
			if( samples )
			{
				PPolynomial< 4 , PolyReal > poly = thetaBase1*thetaBase2;
				stencil[i] =
					(
					poly.integralCosecant( a , 0 , 1 , samples ) +
					poly.integralCosecant( a , 1 , 2 , samples ) +
					poly.integralCosecant( a , 2 , 3 , samples )
					) * a;
			}
			else
			{
				stencil[i]  = Real( ( thetaBase1 * thetaBase2 * Polynomial< 4 , PolyReal >::Cosecant( PolyReal( 1. )                    * a ).scale( PolyReal( 1. ) / a ) ).integral( 0. , 1.  ) * a );
				stencil[i] += Real( ( thetaBase1 * thetaBase2 * Polynomial< 4 , PolyReal >::Cosecant( PolyReal( 1. + end )/PolyReal(2.) * a ).scale( PolyReal( 1. ) / a ) ).integral( 1. , end ) * a );
			}
		}
		else
		{
			if( samples )
			{
				PPolynomial< 4 , PolyReal > poly = thetaBase1*thetaBase2;
				stencil[i] =
					(
					poly.integralCosecant( a , m-1 , m   , samples ) +
					poly.integralCosecant( a , m   , m+1 , samples ) +
					poly.integralCosecant( a , m+1 , m+2 , samples ) 
					) * a;
			}
			else
			{
				stencil[i] = ( thetaBase * thetaBase2 ).integral( start , end ) * a;
			}
		}
	}
}



template< class Real , class PolyReal >
PolyReal EquiRectangular< Real , PolyReal >::BaseWeight( int M , int N , int m , const SphericalStencilTable< double >* stencilTable )
{
	if( stencilTable )
	{
		if( m==0 || m==M-1 ) return stencilTable->thetaWeight( m , M ) * 2. * M_PI;
		else                 return stencilTable->thetaWeight( m , M ) * stencilTable->phiWeight( N );
	}
	else
	{
		PPolynomial< 2 , PolyReal > thetaBase , phiBase;
		PolyReal a =       M_PI / M;
		PolyReal b = 2.0 * M_PI / N;

		double start = m-1. , end = m+2.;
		if( start<0 ) start = 0;
		if(   end>M ) end   = M;
		SetBaseFunctions( m , 0 , M , N , thetaBase , phiBase , false );
		if( m==0 || m==M-1 ) return thetaBase.integralSmoothSine( 2 , a , start , end ) * a * b * PolyReal( N );
		else                 return thetaBase.integralSmoothSine( 2 , a , start , end ) * phiBase.integral( -1.0 , 2.0 ) * a * b;
	}
}
template< class Real , class PolyReal >
void EquiRectangular< Real , PolyReal >::BaseWeights( int M , int N , const SphericalStencilTable< double >* stencilTable , Vector< Real >& weights )
{
	weights.Resize( Dimension( M , N ) );
	if( stencilTable )
	{
		for( int m=0 ; m<M ; m++ )
		{
			double weight;
			if( m==0 || m==M-1 )	weight = stencilTable->thetaWeight( m , M ) * 2. * PI;
			else					weight = stencilTable->thetaWeight( m , M ) * stencilTable->phiWeight( N );
			for( int n=0 ; n<N ; n++ ) weights[ Index( m , n , M , N ) ] = Real( weight );
		}
		return;
	}
	PPolynomial< 2 , PolyReal > thetaBase , phiBase;
	PolyReal a = M_PI / M;
	PolyReal b = 2.0 * M_PI / N;
	EquiRectangular< Real , PolyReal >::SetThetaBaseFunction( 1 , M , thetaBase , false );
	EquiRectangular< Real , PolyReal >::SetPhiBaseFunction  ( 0 , N , phiBase   , false );
	PolyReal thetaSinWeight = thetaBase.integralSmoothSine  ( 2 , a , 0. , 3. ) * a;
	PolyReal thetaCosWeight = thetaBase.integralSmoothCosine( 2 , a , 0. , 3. ) * a;
	PolyReal phiWeight;
	if( stencilTable ) phiWeight = stencilTable->phiWeight( N );
	else phiWeight = phiBase.integral( -1. , 2. ) * b;

	for( int m=0 ; m<M ; m++ )
	{
		Real weight;
		if( m==0 || m==M-1 )
		{
			double start = m-1. , end = m+2.;
			if( start<0 ) start =  0;
			if(   end>M ) end   = M;
			EquiRectangular< Real , PolyReal >::SetThetaBaseFunction( m , M , thetaBase , false );
			weight = Real( thetaBase.integralSmoothSine( 2 , a , start , end ) * a * PolyReal( 2.0 * M_PI ) );
		}
		else weight = Real( ( thetaSinWeight * cos( a * PolyReal( m-1 ) ) + thetaCosWeight * sin( a * PolyReal( m-1 ) ) ) * phiWeight );
		for( int n=0 ; n<N ; n++ ) weights[ Index( m , n , M , N ) ] = weight;
	}
}
template< class Real , class PolyReal >
bool EquiRectangular< Real , PolyReal >::LaplacianMatrix( int M , int N , int samples , const SphericalStencilTable< double >* stencilTable , SparseMatrix< Real >& laplacianMatrix , double iWeight , double lWeight )
{
	double laplacianStencil[5][5];
	double dotProductStencil[5][5];
	laplacianMatrix.Resize( Dimension( M , N ) );
	{
		for( int i=0 ; i < M ; i+=M-1 )
		{
			int ii = 0;
			int idx = Index( i , 0 , M , N );
			LaplacianStencil( M , N , i , laplacianStencil , samples , stencilTable );
			DotProductStencil( M , N , i , dotProductStencil , stencilTable );
			laplacianMatrix.SetGroupSize( idx , 2*N+1 );
			for( int k=-2 ; k<=2 ; k++ )
			{
				if( i+k<0 || i+k>=M ) continue;
				else if( k==0 )
				{
					laplacianMatrix[idx][ii  ].N = idx;
					laplacianMatrix[idx][ii++].Value = laplacianStencil[2][2]*lWeight + dotProductStencil[2][2]*iWeight;
				}
				else
					for( int l=0 ; l<N ; l++ )
					{
						laplacianMatrix[idx][ii  ].N = Index( i+k , l , M , N );
						laplacianMatrix[idx][ii++].Value = laplacianStencil[2+k][2]*lWeight + dotProductStencil[2+k][2]*iWeight;
					}
			}
		}

	}
	for( int i=1 ; i<M-1 ; i++ )
	{
		LaplacianStencil( M , N , i , laplacianStencil , samples , stencilTable );
		DotProductStencil( M , N , i , dotProductStencil , stencilTable );
		for( int j=0; j<N ; j++ )
		{
			int ii = 0;
			int idx = Index( i , j , M , N );
			if		( i==1 || i==M-2 )	laplacianMatrix.SetGroupSize( idx , 3*5+1 );
			else if ( i==2 || i==M-3 )	laplacianMatrix.SetGroupSize( idx , 4*5+1 );
			else						laplacianMatrix.SetGroupSize( idx , 5*5   );

			for( int k=-2 ; k<=2 ; k++ )
			{
				if( i+k<0 || i+k>=M ) continue;
				else if( i+k==0 || i+k==M-1 )
				{
					laplacianMatrix[idx][ii  ].N = Index( i+k , 0 , M , N );
					laplacianMatrix[idx][ii++].Value = laplacianStencil[2+k][2]*lWeight + dotProductStencil[2+k][2]*iWeight;
				}
				else
					for( int l=-2 ; l<=2 ; l++ )
					{
						laplacianMatrix[idx][ii  ].N = Index( i+k , j+l , M , N );
						laplacianMatrix[idx][ii++].Value = laplacianStencil[2+k][2+l]*lWeight + dotProductStencil[2+k][2+l]*iWeight;
					}
			}
		}
	}
	return true;
}

// Hopefully this will make the divergence-matrix work
void _SetDThetaBaseFunction( int m , int M , PPolynomial<1>& thetaBase , bool rescale )
{
	double thetaScale = M_PI / M;
	thetaBase = PPolynomial< 1 >::GaussianApproximation().shift( m+1. );
	if( rescale ) thetaBase = thetaBase.scale( thetaScale ) / thetaScale;
}
void _SetDPhiBaseFunction( int n , int N , PPolynomial<1>& phiBase , bool rescale )
{
	double phiScale = 2.0*M_PI / N;

	phiBase = PPolynomial< 1 >::GaussianApproximation().shift( n+1. );
	if( rescale ) phiBase = phiBase.scale( phiScale ) / phiScale;
}
void PhiStencil( int N , double stencil[5] )
{
	double b = 2.0 * M_PI / N;
	PPolynomial< 2 > thetaBase , phiBase1 , phiBase2;

	EquiRectangular< double >::SetBaseFunctions( 1 , 0 , N , N , thetaBase , phiBase1 , false );
	for( int i=-2 ; i<=2 ; i++ )
	{
		EquiRectangular< double >::SetBaseFunctions( 1 , i , N , N , thetaBase , phiBase2 , false );
		stencil[i+2] = ( phiBase1 * phiBase2 ).integral( -2.0 , 3.0 ) * b;
	}
}
bool DThetaStencil( int M , int N , int m , double stencil[4][5] , int samples )
{
	for( int i=0 ; i<4 ; i++ ) for( int j=0 ; j<5 ; j++ ) stencil[i][j] = 0;
	if( m<0 || m>=M-1 || M<6 || N<=3) return false;

	double a =       M_PI / M;
	double b = 2.0 * M_PI / N;
	double start , end;
	double phiStencil[5];
	PPolynomial< 2 > thetaBase1 , phiBase1 , thetaBase2 , phiBase2 , dThetaBase2;
	PPolynomial< 1 > dThetaBase1;
	EquiRectangular< double >::SetBaseFunctions ( 1 , 0 , M , N , thetaBase1 , phiBase1 , false );
	_SetDThetaBaseFunction( m , M , dThetaBase1 , false );
	PhiStencil( N , phiStencil );
	start = m;
	end   = m+2;
	for( int i=-1 ; i<=2 ; i++ )
	{
		if( i+m <0 || i+m >M-1 ) continue;
		if( i+m==0 || i+m==M-1 )
		{	
			EquiRectangular< double >::SetBaseFunctions( i+m , 0 , M , N , thetaBase2 , phiBase2 , false );
			dThetaBase2 = thetaBase2.derivative( );
			stencil[i+1][2]= - ( dThetaBase1 * dThetaBase2 ).integralSine( a , start , end ) * ( phiBase1 * phiBase2 ).integral( -2.0 , 3.0 ) * b / a;
		}
		else
			for( int j=-2 ; j<=2 ; j++ )
			{
				EquiRectangular< double >::SetBaseFunctions( i+m , j , M , N , thetaBase2 , phiBase2 , false );
				dThetaBase2 = thetaBase2.derivative( );
				stencil[i+1][j+2]= - ( dThetaBase1 * dThetaBase2 ).integralSine( a , start , end ) * phiStencil[j+2] / a;
			}
	}
	return true;
}
void DPhiStencil( int N , double stencil[5] )
{
	double b = 2.0 * M_PI / N;
	PPolynomial< 2 > thetaBase , phiBase , dPhiBase1 , dPhiBase2;

	EquiRectangular< double >::SetBaseFunctions( 1 , 0 , N , N , thetaBase , phiBase , false );
	dPhiBase1 = phiBase.derivative( );
	for( int i=-2 ; i<=2 ; i++ )
	{
		EquiRectangular< double >::SetBaseFunctions( 1 , i , N , N , thetaBase , phiBase , false );
		dPhiBase2 = phiBase.derivative( );
		stencil[i+2] = ( dPhiBase1 * dPhiBase2 ).integral( -2.0 , 3.0 ) / b;
	}
}
void ThetaStencil( int m , int M , double stencil[5] )
{
	double a = M_PI / M;
	double start = m-1. , end = m+2.;
	PPolynomial< 2 > thetaBase1 , thetaBase2 , phiBase;
	if( start<0 ) start = 0;
	if(   end>M ) end   = M;

	EquiRectangular< double >::SetBaseFunctions( m , 0 , M , M , thetaBase1 , phiBase , false );
	for( int i=-2 ; i<=2 ; i++ )
	{
		EquiRectangular< double >::SetBaseFunctions( m+i , 0 , M , M , thetaBase2 , phiBase , false );
		stencil[i+2] = ( thetaBase1 * thetaBase2 ).integralSine( a , start , end ) * a;
	}
}
void ThetaStencil( int m , int M , double stencil[5] , int samples )
{
	for( int i=0 ; i<5 ; i++ ) stencil[i] = 0;

	bool flip;
	if( m>=M/2 )	flip = true , m = M-1-m;
	else			flip = false;
	double a = M_PI / M;
	double start = m-1. , end = m+2.;
	PPolynomial< 2 , double > thetaBase1 , thetaBase2 , phiBase;
	EquiRectangular< double >::SetBaseFunctions( m , 0 , M , M , thetaBase1 , phiBase , false );
	if( start<0 ) start = 0;
	if(   end>M ) end   = M;
	for( int i=-2 ; i<=2 ; i++ )
	{
		if(  i+m<=0 || i+m>=M-1 ) continue;
		EquiRectangular< double >::SetBaseFunctions( i+m , 0 , M , M , thetaBase2 , phiBase , false );
		if( m==1 )
		{
			Polynomial< 10 , double > xCscX = Polynomial< 4 , double >::XCscX( ).scale( 1./a );
			stencil[i+2]  = ( thetaBase2*xCscX ).integral( 0 , 1. ) * a;
			stencil[i+2] += ( thetaBase1 * thetaBase2 ).integralCosecant( a , 1. , end , samples ) * a;
		}
		else if ( i+m==1 )
		{
			Polynomial< 10 , double > xCscX = Polynomial< 4 , double >::XCscX( ).scale( 1./a );
			stencil[i+2]  = ( thetaBase1*xCscX ).integral( 0 , 1. ) * a;
			stencil[i+2] += ( thetaBase1 * thetaBase2 ).integralCosecant( a , 1. , end , samples ) * a;
		}
		else
			stencil[i+2] = ( thetaBase1 * thetaBase2 ).integralCosecant( a , start , end , samples ) * a;
	}
	if( flip )
	{
		double temp;
		temp = stencil[4];
		stencil[4] = stencil[0];
		stencil[0] = temp;
		temp = stencil[3];
		stencil[3] = stencil[1];
		stencil[1] = temp;
	}
}

bool DPhiStencil( int M , int N , int m , double stencil[5][4] , int samples )
{
	for( int i=0 ; i<5 ; i++ ) for( int j=0 ; j<4 ; j++ ) stencil[i][j] = 0;
	if( m<=0 || m>=M-1 || M<6 || N<=3) return false;

	double a =       M_PI / M;
	double b = 2.0 * M_PI / N;
	double thetaStencil[5];
	PPolynomial< 2 > thetaBase , phiBase , dPhiBase2;
	PPolynomial< 1 > dPhiBase1;
	ThetaStencil( m , M , thetaStencil , samples );
	_SetDPhiBaseFunction( 0 , N , dPhiBase1 , false );
	for( int i=-2 ; i<=2 ; i++ )
	{
		if( i+m<=0 || i+m>=M-1 ) continue;
		else
			for( int j=-1 ; j<=2 ; j++ )
			{
				EquiRectangular< double >::SetBaseFunctions( i+m , j , M , N , thetaBase , phiBase , false );
				dPhiBase2 = phiBase.derivative( );
				stencil[i+2][j+1]= - thetaStencil[i+2] * ( dPhiBase1 * dPhiBase2 ).integral( 0.0 , 2.0 ) / b;
			}
	}
	return true;
}

template< class Real , class PolyReal > int EquiRectangular< Real , PolyReal >::Dimension( int M , int N ) { return (M-2)*N+2; }
template< class Real , class PolyReal >
void EquiRectangular< Real , PolyReal >::Index( int idx , int M , int N , int& m , int& n )
{
	if		( idx<=0 ) m = 0 , n = 0;
	else if	( idx>=Dimension( M , N )-1 ) m = M-1 , n = 0;
	else
	{
		// idx = 1 + (m-1)*N + n
		m = ( (idx-1) / N ) + 1;
		n =   (idx-1) % N;
	}
}
template< class Real , class PolyReal >
int EquiRectangular< Real , PolyReal >::Index( int m , int n , int M , int N , int& mm , int& nn )
{
	mm =  m;
	nn = (n+N) % N;
	if( mm<0 )
	{
		mm = -1-mm;
		nn = (N-n) % N;
	}
	else if( mm>M-1 )
	{
		mm = (M-1) - (mm-M);
		nn = (N-n) % N;
	}
	if( mm==  0 )	return 0;
	if( mm==M-1 )	return 1 + ( M-2)*N;
	else			return 1 + (mm-1)*N + nn;
}

template< class Real , class PolyReal >
int EquiRectangular< Real , PolyReal >::Index( int m , int n , int M , int N )
{
	int mm , nn;
	return Index( m , n , M , N , mm , nn );
}
template< class Real , class PolyReal >
void EquiRectangular< Real , PolyReal >::Index( int m , int n , int M , int N , double& theta , double& phi )
{
	int mm , nn;
	Index( m , n , M , N , mm , nn );
	if		( mm==0   )	theta = 0;
	else if	( mm==M-1 ) theta = M_PI;
	else				theta = M_PI * (mm+0.5) / M;
	phi = 2.0 * M_PI * (nn+0.5) / N;
}
template< class Real , class PolyReal >
void EquiRectangular< Real , PolyReal >::PhiEvaluation( double evaluation[3] )
{
	PPolynomial< 2 , PolyReal > thetaBase , phiBase;
	PolyReal p = 0.5;
	for( int i=-1 ; i<=1 ; i++ )
	{
		SetBaseFunctions( 1 , i , 4 , 4 , thetaBase , phiBase , false );
		evaluation[i+1] = Real( phiBase( p ) );
	}
}
template< class Real , class PolyReal >
void EquiRectangular< Real , PolyReal >::PrimalPhiEvaluation( double evaluation[2] )
{
	PPolynomial< 2 , PolyReal > thetaBase , phiBase;
	PolyReal p = 0.;
	for( int i=-1 ; i<=0 ; i++ )
	{
		SetBaseFunctions( 1 , i , 4 , 4 , thetaBase , phiBase , false );
		evaluation[i+1] = Real( phiBase( p ) );
	}
}
template< class Real , class PolyReal >
void EquiRectangular< Real , PolyReal >::CoarsePhiEvaluation( double evenStencil[3] , double oddStencil[3] )
{
	PPolynomial< 2 , PolyReal > thetaBase , phiBase;
	PolyReal pEven = 0.25;
	PolyReal pOdd  = 0.75;
	for( int i=-1 ; i<=1 ; i++ )
	{
		SetBaseFunctions( 1 , i , 4 , 4 , thetaBase , phiBase , false );
		evenStencil[i+1] = Real( phiBase( pEven ) );
		oddStencil[i+1]  = Real( phiBase( pOdd  ) );
	}
}
template< class Real , class PolyReal >
void EquiRectangular< Real , PolyReal >::PrimalCoarsePhiEvaluation( double evenStencil[2] , double oddStencil[3] )
{
	PPolynomial< 2 , PolyReal > thetaBase , phiBase;
	PolyReal pEven = 0.;
	PolyReal pOdd  = 0.5;
	for( int i=-1 ; i<=0 ; i++ )
	{
		SetBaseFunctions( 1 , i , 4 , 4 , thetaBase , phiBase , false );
		evenStencil[i+1] = Real( phiBase( pEven ) );
	}
	for( int i=-1 ; i<=1 ; i++ )
	{
		SetBaseFunctions( 1 , i , 4 , 4 , thetaBase , phiBase , false );
		oddStencil[i+1]  = Real( phiBase( pOdd  ) );
	}
}
template< class Real , class PolyReal >
void EquiRectangular< Real , PolyReal >::FinePhiEvaluation( double stencil[2] )
{
	PPolynomial< 2 , PolyReal > thetaBase , phiBase;
	PolyReal p = 1.;
	for( int i=0 ; i<=1 ; i++ )
	{
		SetBaseFunctions( 1 , i , 4 , 4 , thetaBase , phiBase , false );
		stencil[i] = Real( phiBase( p ) );
	}
}
template< class Real , class PolyReal >
void EquiRectangular< Real , PolyReal >::PrimalFinePhiEvaluation( double stencil[2] )
{
	PPolynomial< 2 , PolyReal > thetaBase , phiBase;
	PolyReal p = 0.;
	for( int i=-1 ; i<=0 ; i++ )
	{
		SetBaseFunctions( 1 , i , 4 , 4 , thetaBase , phiBase , false );
		stencil[i+1] = Real( phiBase( p ) );
	}
}
template< class Real , class PolyReal >
void EquiRectangular< Real , PolyReal >::PhiPCIntegral( double stencil[3] )
{
//	PolyReal b = 2.0 * M_PI / N;
	PolyReal b = 2.0 * M_PI;
	PPolynomial< 2 , PolyReal > thetaBase , phiBase;

	SetBaseFunctions( 1 , 0 , 4 , 4 , thetaBase , phiBase , false );
	for( int i=-1 ; i<=1 ; i++ ) stencil[i+1] = phiBase.integral( i , i+1. ) * b;
}
template< class Real , class PolyReal >
void EquiRectangular< Real , PolyReal >::CoarsePhiPCIntegral( double evenStencil[2] , double oddStencil[2] )
{
	PolyReal b = 2.0 * M_PI;
	PPolynomial< 2 , PolyReal > thetaBase , phiBase;

	SetBaseFunctions( 1 , 0 , 4 , 4 , thetaBase , phiBase , false );
	evenStencil[0] = phiBase.integral( -2. , 0. ) * b;
	evenStencil[1] = phiBase.integral(  0. , 2. ) * b;
	oddStencil[0]  = phiBase.integral( -1. , 1. ) * b;
	oddStencil[1]  = phiBase.integral(  1. , 3. ) * b;
}
template< class Real , class PolyReal >
void EquiRectangular< Real , PolyReal >::FinePhiPCIntegral( double stencil[6] )
{
	PolyReal b = 2.0 * M_PI;
	PPolynomial< 2 , PolyReal > thetaBase , phiBase;

	SetBaseFunctions( 1 , 0 , 4 , 4 , thetaBase , phiBase , false );
	for( int i=0 ; i<6 ; i++ ) stencil[i] = phiBase.integral( -1. + i*.5 , -1. + (i+1)*.5 ) * b;
}
template< class Real , class PolyReal >
void EquiRectangular< Real , PolyReal >::ThetaEvaluation( int m , int M , double evaluation[3] )
{
	PolyReal p;
	PPolynomial< 2 , PolyReal > thetaBase , phiBase;

	if     ( m==0   ) p = 0;
	else if( m==M-1 ) p = M;
	else              p = m+0.5;
	for( int i=-1 ; i<=1 ; i++ )
		if( m+i>=0 && m+i<M )
		{
			SetBaseFunctions( m+i , 0 , M , M , thetaBase , phiBase , false );
			evaluation[i+1] = Real( thetaBase( p ) );
		}
		else evaluation[i+1] = 0;
}
template< class Real , class PolyReal >
void EquiRectangular< Real , PolyReal >::PrimalThetaEvaluation( int m , int M , double evaluation[2] )
{
	PolyReal p = m + 1.;
	PPolynomial< 2 , PolyReal > thetaBase , phiBase;

	for( int i=0 ; i<=1 ; i++ )
		if( m+i>=0 && m+i<M )
		{
			SetBaseFunctions( m+i , 0 , M , M , thetaBase , phiBase , false );
			evaluation[i] = Real( thetaBase( p ) );
		}
		else evaluation[i] = 0;
}
template< class Real , class PolyReal >
void EquiRectangular< Real , PolyReal >::ThetaPCIntegral( int m , int M , double stencil[3] )
{
	PolyReal a = M_PI / M;
	PPolynomial< 2 , PolyReal > thetaBase , phiBase;

	SetBaseFunctions( m , 0 , M , M , thetaBase , phiBase , false );
	for( int i=-1 ; i<=1 ; i++ )
	{
		int start = m+i , end = m+i+1;
		if( m+i>=0 && m+i<M ) stencil[i+1] = thetaBase.integralSmoothSine( 2 , a , start , end ) * a;
		else                  stencil[i+1] = 0;
	}
}

#include <QLIB/QFLOAT.h>
template< class Real , class PolyReal >
bool EquiRectangular< Real , PolyReal >::DivergenceDThetaStencil( int M , int N , int m , double stencil[4][5] , int samples )
{
	for( int i=0 ; i<4 ; i++ ) for( int j=0 ; j<5 ; j++ ) stencil[i][j] = 0;
	double phiStencil[5] , divergenceDThetaStencil[4];
	if( m==0 || m==M-1 ) phiStencil[0] = phiStencil[1] = phiStencil[2] = phiStencil[3] = phiStencil[4] = PhiWeight( N );
	else PhiStencil( N , phiStencil );
	if( !DivergenceDThetaStencil( m , M , divergenceDThetaStencil ) || N<=3 ) return false;
	for( int i=-2 ; i<=1 ; i++ ) for( int j=-2 ; j<=2 ; j++ ) stencil[i+2][j+2] = - divergenceDThetaStencil[i+2] * phiStencil[j+2];
	return true;
}
template< class Real , class PolyReal >
bool EquiRectangular< Real , PolyReal >::DivergenceDPhiStencil( int M , int N , int m , double stencil[5][4] , int samples )
{
	for( int i=0 ; i<5 ; i++ ) for( int j=0 ; j<4 ; j++ ) stencil[i][j] = 0;
	double thetaStencil[5] , divergenceDPhiStencil[4];
	ThetaStencil( m , M , thetaStencil , samples );
	if( m<=0 || m>=M-1 || M<6 ) return false;
	if( !DivergenceDPhiStencil( N , divergenceDPhiStencil ) ) return false;
	for( int i=-2 ; i<=2 ; i++ )
		if( i+m<=0 || i+m>=M-1 ) continue;
		else for( int j=-2 ; j<=1 ; j++ ) stencil[i+2][j+2] = - thetaStencil[i+2] * divergenceDPhiStencil[j+2];
	return true;
}
template< class Real , class PolyReal >
bool EquiRectangular< Real , PolyReal >::DotProductStencil( int M , int N , int m , double stencil[5][5] , const SphericalStencilTable< double >* stencilTable )
{
	if( m<0 || m>=M || M<6 || N<=3 ) return false;
	for( int i=0 ; i<5 ; i++ ) for( int j=0 ; j<5 ; j++ ) stencil[i][j] = 0;
	PolyReal a =       M_PI / M;
	PolyReal b = 2.0 * M_PI / N;
	double thetaStencil[5] , phiStencil[5];
	PPolynomial< 2 , PolyReal > thetaBase1 , phiBase1;
	PPolynomial< 2 , PolyReal > thetaBase2 , phiBase2;
	SetBaseFunctions( m , 0 , M , N , thetaBase1 , phiBase1 , false );
	if( stencilTable ) stencilTable->setThetaStencil( thetaStencil , m , M );
	else ThetaStencil( m , M , thetaStencil );
	PhiStencil( N , phiStencil );

	for( int i=-2 ; i<=2 ; i++ )
	{
		if( i+m<0 || i+m> M-1 ) continue;
		if     ( ( i+m==0 || i+m==M-1 ) && ( m==0 || m==M-1 ) ) stencil[i+2][2] = Real( PolyReal( thetaStencil[i+2] ) * PolyReal( N )  * b );
		else if(   i+m==0 || i+m==M-1   ||   m==0 || m==M-1   )
		{	
			SetBaseFunctions( i+m , 0 , M , N , thetaBase2 , phiBase2 , false );
			stencil[i+2][2] = Real( PolyReal( thetaStencil[i+2] ) * ( phiBase1 * phiBase2 ).integral( -2.0 , 3.0 ) * b );
		}
		else for( int j=-2 ; j<=2 ; j++ ) stencil[i+2][j+2] = thetaStencil[i+2] * phiStencil[j+2];
	}
	return true;
}
template< class Real , class PolyReal >
bool EquiRectangular< Real , PolyReal >::EvaluationStencil( int M , int N , int m , double stencil[3][3] )
{
	if( m<0 || m>=M || M<6 || N<=3 ) return false;
	for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) stencil[i][j] = 0;
	double thetaStencil[3] , phiStencil[3];
	ThetaEvaluation( m , M , thetaStencil );
	PhiEvaluation  ( phiStencil );

	for( int i=-1 ; i<=1 ; i++ )
	{
		if     ( i+m< 0 || i+m> M-1 ) continue;
		else if(   m==0 ||   m==M-1 ) stencil[i+1][1] = thetaStencil[i+1];
		else if( i+m==0 || i+m==M-1 ) stencil[i+1][1] = thetaStencil[i+1];
		else for( int j=-1 ; j<=1 ; j++ ) stencil[i+1][j+1] = thetaStencil[i+1] * phiStencil[j+1];
	}
	return true;
}
template< class Real , class PolyReal >
bool EquiRectangular< Real , PolyReal >::PrimalEvaluationStencil( int M , int N , int m , double stencil[2][2] )
{
	if( m<0 || m>=M-1 ) return false;
	for( int i=0 ; i<2 ; i++ ) for( int j=0 ; j<2 ; j++ ) stencil[i][j] = 0;
	double thetaStencil[2] , phiStencil[2];
	PrimalThetaEvaluation( m , M , thetaStencil );
	PrimalPhiEvaluation  ( phiStencil );

	for( int i=0 ; i<2 ; i++ )
	{
		if( i+m==0 || i+m==M-1 ) stencil[i][1] = thetaStencil[i+1];
		else for( int j=0 ; j<2 ; j++ ) stencil[i][j] = thetaStencil[i] * phiStencil[j];
	}
	return true;
}
template< class Real , class PolyReal >
bool EquiRectangular< Real , PolyReal >::PCIntegralStencil( int M , int N , int m , double stencil[3][3] )
{
	if( m<0 || m>=M || M<6 || N<=3 ) return false;
	for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) stencil[i][j] = 0;
	PolyReal a =       M_PI / M;
	PolyReal b = 2.0 * M_PI / N;
	double thetaStencil[3] , phiStencil[3];
	ThetaPCIntegral( m , M , thetaStencil );
	PhiPCIntegral( phiStencil );

	for( int i=-1 ; i<=1 ; i++ )
		if     ( i+m< 0 || i+m> M-1 ) continue;
		else if(  (m==0 ||   m==M-1) && i==0 ) stencil[1][1] = Real( PolyReal( thetaStencil[1] ) * b * PolyReal( N ) );
		else if(   m==0 ||   m==M-1 ) stencil[i+1][1] = Real( PolyReal( thetaStencil[i+1] ) * b );
		else if( i+m==0 || i+m==M-1 ) stencil[i+1][1] = Real( PolyReal( thetaStencil[i+1] ) * ( PolyReal( phiStencil[0] )+PolyReal( phiStencil[1] )+PolyReal( phiStencil[2] ) ) / PolyReal( N ) );
		else for( int j=-1 ; j<=1 ; j++ ) stencil[i+1][j+1] = Real( PolyReal( thetaStencil[i+1] ) * PolyReal( phiStencil[j+1] ) / PolyReal( N ) );
	return true;
}
template< class Real , class PolyReal >
bool EquiRectangular< Real , PolyReal >::LaplacianStencil( int M , int N , int m , double stencil[5][5] , int samples , const SphericalStencilTable< double >* stencilTable , bool negate )
{
	for( int i=0 ; i<5 ; i++ ) for( int j=0 ; j<5 ; j++ ) stencil[i][j] = 0;
	if( m<0 || m>=M || M<6 || N<=3) return false;

	double laplacianThetaStencil[5] , laplacianD2ThetaStencil[5] , phiStencil[5] , laplacianD2PhiStencil[5];
	PolyReal a =       M_PI / M;
	PolyReal b = 2.0 * M_PI / N;

	PPolynomial< 2 , PolyReal > thetaBase1 , phiBase1 , dPhiBase1;
	PPolynomial< 2 , PolyReal > thetaBase2 , phiBase2 , dPhiBase2;
	if( stencilTable )
	{
		stencilTable->setLaplacianThetaStencil  (   laplacianThetaStencil , m , M , samples );
		stencilTable->setLaplacianD2ThetaStencil( laplacianD2ThetaStencil , m , M );
		stencilTable->setPhiStencil           ( phiStencil            , N );
		stencilTable->setLaplacianD2PhiStencil( laplacianD2PhiStencil , N );
	}
	else
	{
		LaplacianThetaStencil( m , M , laplacianThetaStencil , samples );
		LaplacianD2ThetaStencil( m , M , laplacianD2ThetaStencil );
		PhiStencil( N , phiStencil );
		LaplacianD2PhiStencil( N , laplacianD2PhiStencil );
	}

	SetBaseFunctions( m , 0 , M , N , thetaBase1 , phiBase1 , false );
	dPhiBase1   = phiBase1.derivative( );
	for( int i=-2 ; i<=2 ; i++ )
	{
		if(  i+m< 0 || i+m> M-1 ) continue;
		if( (i+m==0 || i+m==M-1) && (m==0 || m==M-1) )	stencil[i+2][2] = Real( - PolyReal( laplacianD2ThetaStencil[i+2] ) * PolyReal( N ) * b );
		else if( i+m==0 || i+m==M-1 )					stencil[i+2][2] = Real( - PolyReal( laplacianD2ThetaStencil[i+2] ) * phiBase1.integral( -2.0 , 3.0 ) * b );
		else if(   m==0 ||   m==M-1 )
		{	
			SetBaseFunctions( i+m , 0 , M , N , thetaBase2 , phiBase2 , false );
			stencil[i+2][2] = Real( - PolyReal( laplacianD2ThetaStencil[i+2] ) * phiBase2.integral( -2.0 , 3.0 ) * b );
		}
		else for( int j=-2 ; j<=2 ; j++ ) stencil[i+2][j+2] = - laplacianThetaStencil[i+2] * laplacianD2PhiStencil[j+2] - laplacianD2ThetaStencil[i+2] * phiStencil[j+2];
	}
	if( negate ) for( int i=0 ; i<5 ; i++ ) for( int j=0 ; j<5 ; j++ ) stencil[i][j] = -stencil[i][j];

	return true;
}
