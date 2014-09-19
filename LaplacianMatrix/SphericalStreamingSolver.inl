#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif // M_PI

#define MyModIndex( idx , mod ) ( ModIndex( (idx) , (mod) ) )

#define AllocAlignedMemory( _buffer , buffer , BufferType , size , ALIGNMENT )		\
	_buffer = malloc( (size) * sizeof( BufferType ) + ALIGNMENT-1 );				\
	memset( _buffer , 0 , (size) + ALIGNMENT-1 );									\
	buffer  = (BufferType*)( ( (size_t)( _buffer ) + ALIGNMENT-1 ) & ~(ALIGNMENT-1) )

#define MyFree( buffer ) if( buffer ) free( buffer ) , buffer = NULL

//////////////////////////////
// SphericalStreamingSolver //
//////////////////////////////
template< class Real >
MultiStreamIOServer SphericalStreamingSolver< Real >::server;

template< class Real >
SphericalStreamingSolver< Real >::SphericalStreamingSolver( void )
{
	setResidual = true;
	bSize = rSize = 0;
	RStream = NullPointer< __m128 >( );
	XStream = NullPointer< __m128 >( );
	BStream = NullPointer< __m128 >( );
	thetaStencils = NULL;
	laplacianStencils = NULL;
	_laplacianStencilsSSE = _zeroCenteredLaplacianStencilsSSE = NULL;
	_phiStencilData = NULL;
	_scratch0 = NullPointer< Real >( );
	_scratch1 = NullPointer< Real >( );
	_scratch2 = NullPointer< Real >( );
	bSquareNorm = xSquareNorm = rSquareNorm = NULL;
	_stencilTable = NULL;
}

template< class Real >
SphericalStreamingSolver< Real >::~SphericalStreamingSolver( void )
{
	UnSetIterations();
	MyFree( bSquareNorm );
	MyFree( xSquareNorm );
	MyFree( rSquareNorm );
}
template< class Real >
Pointer( Real ) SphericalStreamingSolver< Real >::GetXRow( int row , int channel )
{
	return ( Pointer( Real ) )( XStream + MyModIndex( row , xSize ) * (_columns+2) * _channels + channel * (_columns+2) + 1 );
}
template< class Real >
Pointer( Real ) SphericalStreamingSolver< Real >::GetBRow( int row , int channel )
{
	return ( Pointer( Real ) )( BStream + MyModIndex( row , bSize ) * (_columns+2) * _channels + channel * (_columns+2) + 1 );
}
template< class Real >
Pointer( Real ) SphericalStreamingSolver< Real >::GetRRow( int row , int channel )
{
	return ( Pointer( Real ) )( RStream + MyModIndex( row , rSize ) * (_columns+2) * _channels + channel * (_columns+2) + 1 );
}
template< class Real >
typename SphericalStreamingSolver< Real >::ThetaStencilData& SphericalStreamingSolver< Real >::GetThetaStencil( int row )
{
	return thetaStencils[ MyModIndex( row , bSize+2 ) ];
}
template< class Real >
typename SphericalStreamingSolver< Real >::SphericalLaplacianStencil& SphericalStreamingSolver< Real >::GetLaplacianStencil( int row )
{
	return laplacianStencils[ MyModIndex( row , bSize ) ];
}
template< class Real >
typename SphericalStreamingSolver< Real >::SphericalLaplacianStencilSSE& SphericalStreamingSolver< Real >::GetLaplacianStencilSSE( int row )
{
	return laplacianStencilsSSE[ MyModIndex( row , bSize ) ];
}
template< class Real >
typename SphericalStreamingSolver< Real >::SphericalLaplacianStencilSSE& SphericalStreamingSolver< Real >::GetZeroCenteredLaplacianStencilSSE( int row )
{
	return zeroCenteredLaplacianStencilsSSE[ MyModIndex( row , bSize ) ];
}
void UpStencil( const double inStencil[5] , double stencil[8] )
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
void DownStencil( const double inStencil[5] , double evenStencil[4] , double oddStencil[4] )
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

template< class Real >
void SphericalStreamingSolver< Real >::PhiStencilData::set( int width , const SphericalStencilTable< double >* stencilTable )
{
	if( stencilTable )
	{
		phiWeight = stencilTable->phiWeight   (                         width );
		stencilTable->setPhiStencil           ( phiStencil            , width );
		stencilTable->setLaplacianD2PhiStencil( laplacianD2PhiStencil , width );
	}
	else
	{
#if HIGH_PRECISION
		phiWeight = EquiRectangular< qfloat , qfloat >::PhiWeight( width );
		EquiRectangular< qfloat , qfloat >::PhiStencil ( width ,  phiStencil );
		EquiRectangular< qfloat , qfloat >::LaplacianD2PhiStencil( width , laplacianD2PhiStencil );
#else // !HIGH_PRECISION
		phiWeight = EquiRectangular< double >::PhiWeight( width );
		EquiRectangular< double >::PhiStencil           ( width ,            phiStencil );
		EquiRectangular< double >::LaplacianD2PhiStencil( width , laplacianD2PhiStencil );
#endif // HIGH_PRECISION
	}
	UpStencil  (            phiStencil ,            upPhiStencil );
	UpStencil  ( laplacianD2PhiStencil , upLaplacianD2PhiStencil );
	DownStencil(            phiStencil ,            downEvenPhiStencil ,            downOddPhiStencil );
	DownStencil( laplacianD2PhiStencil , downEvenLaplacianD2PhiStencil , downOddLaplacianD2PhiStencil );
}
template< class Real >
void SphericalStreamingSolver< Real >::ThetaStencilData::set( int row , int rows , int samples , const SphericalStencilTable< double >* stencilTable )
{
	if( stencilTable )
	{
		stencilTable->setThetaStencil           (            thetaStencil , row , rows );
		stencilTable->setLaplacianThetaStencil  (   laplacianThetaStencil , row , rows , samples );
		stencilTable->setLaplacianD2ThetaStencil( laplacianD2ThetaStencil , row , rows );
	}
	else
	{
#if HIGH_PRECISION
		EquiRectangular< qfloat , qfloat >::ThetaStencil           ( row , rows , thetaStencil );
		EquiRectangular< qfloat , qfloat >::LaplacianD2ThetaStencil( row , rows , laplacianD2ThetaStencil );
#else // !HIGH_PRECISION
		EquiRectangular< double >::ThetaStencil           ( row , rows , thetaStencil );
		EquiRectangular< double >::LaplacianD2ThetaStencil( row , rows , laplacianD2ThetaStencil );
#endif // HIGH_PRECISION
		EquiRectangular< double >::LaplacianThetaStencil ( row , rows ,  laplacianD2ThetaStencil , samples );
	}
}

template< class Real >
void SphericalStreamingSolver< Real >::SetThetaStencil( int row )
{

	ThetaStencilData& ts = GetThetaStencil( row );
	ts.set( row , rows , samples , _stencilTable );
}

template< class Real >
void SphericalStreamingSolver< Real >::SetStencil( SphericalLaplacianStencil& stencil , int row , double dWeight , double lWeight )
{
	ThetaStencilData& thetaStencilData = GetThetaStencil( row );
	stencil.isRegular = true;


	int sIndex[5];
	for( int i=-2 ; i<=2 ; i++ )
	{
		int idx = 0;
		while( (1<<idx) <  RowWidth( row+i ) ) idx++;
		sIndex[2+i] = idx;
	}
	// If we are looking at the poles
	if( row==0 || row==rows-1 )
	{
		for( int k=-2 ; k<=2 ; k++ )
			if( row+k<0 || row+k>=rows ) continue;
			else if( k==0 ) stencil.evenMatrixValues[Degree][0] = stencil.oddMatrixValues[Degree][0] = - thetaStencilData.laplacianD2ThetaStencil[k+Degree] * 2.0 * M_PI * lWeight + thetaStencilData.thetaStencil[k+Degree] * 2.0 * M_PI * dWeight;
			else
			{
				int width = RowWidth ( row+k );
				double phiWeight = _phiStencilData[ sIndex[2+k] ].phiWeight;
				double val = - thetaStencilData.laplacianD2ThetaStencil[k+Degree] * phiWeight * lWeight + thetaStencilData.thetaStencil[k+Degree] * phiWeight * dWeight;
				stencil.evenMatrixValues[Degree+k][0] = stencil.oddMatrixValues[Degree+k][0] = val;
			}
		stencil.isRegular = false;
		stencil.diagonal = stencil.evenMatrixValues[Degree][0];
	}
	else
	{

		int width = RowWidth( row );
		double phiWeight = _phiStencilData[ sIndex[2] ].phiWeight;
		for( int k=-2 ; k<=2 ; k++ )
			if( row+k>=0 && row+k<rows )
			{
				int neighborWidth = RowWidth( row+k );
				int idx = sIndex[2]>sIndex[2+k] ? sIndex[2] : sIndex[2+k];
				PhiStencilData& phiStencilData = _phiStencilData[idx];
				if     ( row+k==0 || row+k==rows-1 );
				else if(   width==  neighborWidth )	;
				else if(   width==2*neighborWidth )	;
				else if( 2*width==  neighborWidth )	;
				else fprintf( stderr , "Badness in SphericalStreamingSolver< Real >::SetStencil\n" ) , exit(0);
				if( row+k==0 || row+k==rows-1 )		// Pole neighbor
				{
					stencil.evenMatrixValues[k+2][0] = stencil.oddMatrixValues[k+2][0] =
						- thetaStencilData.laplacianD2ThetaStencil[2+k] * phiWeight * lWeight
						+ thetaStencilData.thetaStencil[2+k] * phiWeight * dWeight;
					stencil.isRegular = false;
				}
				else if( width==neighborWidth )		// Regular neighbor
				{
					for( int l=-2 ; l<=2 ; l++ )
					{
						stencil.evenMatrixValues[k+2][l+2] = stencil.oddMatrixValues[k+2][l+2] =
							- thetaStencilData.laplacianD2ThetaStencil[2+k] * phiStencilData.phiStencil [2+l] * lWeight
							- thetaStencilData.laplacianThetaStencil[2+k] * phiStencilData.laplacianD2PhiStencil[2+l] * lWeight
							+ thetaStencilData.thetaStencil[2+k] * phiStencilData.phiStencil [2+l] * dWeight;
					}
				}
				else if( width<neighborWidth )		// Finer neighbor
				{
					for( int l=-3 ; l<=4 ; l++ )
					{
						stencil.evenMatrixValues[k+2][l+3] = stencil.oddMatrixValues[k+2][l+3] =
							- thetaStencilData.laplacianD2ThetaStencil[2+k] * phiStencilData.upPhiStencil [l+3] * lWeight
							- thetaStencilData.laplacianThetaStencil  [2+k] * phiStencilData.upLaplacianD2PhiStencil[l+3] * lWeight
							+ thetaStencilData.thetaStencil[2+k] * phiStencilData.upPhiStencil [l+3] * dWeight;
					}
					stencil.isRegular = false;
				}
				else if( width>neighborWidth )		// Coarser neighbor
				{
					for( int l=-2 ; l<=1 ;l++ )
						stencil.evenMatrixValues[k+2][l+2] =
							- thetaStencilData.laplacianD2ThetaStencil[2+k] * phiStencilData.downEvenPhiStencil[l+2] * lWeight
							- thetaStencilData.laplacianThetaStencil  [2+k] * phiStencilData.downEvenLaplacianD2PhiStencil[l+2] * lWeight
							+ thetaStencilData.thetaStencil[2+k] * phiStencilData.downEvenPhiStencil [l+2] * dWeight;
					for( int l=-1 ; l<=2 ;l++ )
						stencil.oddMatrixValues [k+2][l+1] =
							- thetaStencilData.laplacianD2ThetaStencil[2+k] * phiStencilData.downOddPhiStencil[l+1] * lWeight
							- thetaStencilData.laplacianThetaStencil  [2+k] * phiStencilData.downOddLaplacianD2PhiStencil[l+1] * lWeight
							+ thetaStencilData.thetaStencil[2+k] * phiStencilData.downOddPhiStencil [l+1] * dWeight;
					stencil.isRegular = false;
				}
				stencil.diagonal = stencil.evenMatrixValues[Degree][Degree];
			}
	}
}
template< class Real >
void SphericalStreamingSolver< Real >::StencilToStencilSSE( const SphericalLaplacianStencil& inStencil , SphericalLaplacianStencilSSE& outStencil , bool zeroCenter )
{
	__declspec ( align(16) ) float scratch[4];

	// Set the diagonal entries
	for( int i=0 ; i<4 ; i++ ) outStencil.diagonal[i] = 1.f / inStencil.evenMatrixValues[Degree][Degree];

	// Now set the matrix entries
	for( int i=0 ; i<=2*Degree ; i++ )
	{
		// OFFSET = 0
		scratch[0] = inStencil.evenMatrixValues[i][2];
		scratch[1] = inStencil.evenMatrixValues[i][3];
		scratch[2] = inStencil.evenMatrixValues[i][4];
		scratch[3] = inStencil.evenMatrixValues[i][1];
		if( zeroCenter && i==Degree ) scratch[0] = 0;
		outStencil.matrixValues[0][i] = _mm_load_ps( scratch );

		// OFFSET = 1
		scratch[0] = inStencil.evenMatrixValues[i][1];
		scratch[1] = inStencil.evenMatrixValues[i][2];
		scratch[2] = inStencil.evenMatrixValues[i][3];
		scratch[3] = inStencil.evenMatrixValues[i][4];
		if( zeroCenter && i==Degree ) scratch[1] = 0;
		outStencil.matrixValues[1][i] = _mm_load_ps( scratch );

		// OFFSET = 2
		scratch[0] = inStencil.evenMatrixValues[i][0];
		scratch[1] = inStencil.evenMatrixValues[i][1];
		scratch[2] = inStencil.evenMatrixValues[i][2];
		scratch[3] = inStencil.evenMatrixValues[i][3];
		if( zeroCenter && i==Degree ) scratch[2] = 0;
		outStencil.matrixValues[2][i] = _mm_load_ps( scratch );

		// OFFSET = 3
		scratch[0] = inStencil.evenMatrixValues[i][3];
		scratch[1] = inStencil.evenMatrixValues[i][0];
		scratch[2] = inStencil.evenMatrixValues[i][1];
		scratch[3] = inStencil.evenMatrixValues[i][2];
		if( zeroCenter && i==Degree ) scratch[3] = 0;
		outStencil.matrixValues[3][i] = _mm_load_ps( scratch );
	}
}

template< class Real >
int SphericalStreamingSolver< Real >::RowWidth( int row ) const
{
	if( row< 0 || row> rows-1 ) return 0;
	return _sphere.dimension( row );
}

template< class Real >
void SphericalStreamingSolver< Real >::Init( int channels , int columns , int rows , int samples , const SphericalStencilTable< double >* stencilTable , double dWeight , double lWeight )
{
	_sphere = AdaptiveEquiRectangular< Real >( rows , samples );
	_stencilTable = stencilTable;
	MyFree( bSquareNorm );
	MyFree( xSquareNorm );
	MyFree( rSquareNorm );
	bSquareNorm = new double[channels];
	xSquareNorm = new double[channels];
	rSquareNorm = new double[channels];

	this->_channels = channels;
	this->samples = samples;
	this->dWeight = dWeight;
	this->lWeight = lWeight;
	this->columns = columns;
	this->rows    = rows;
	_columns = ( columns+RealPerWord-1 ) / RealPerWord;
}
template< class Real >
int SphericalStreamingSolver< Real >::Channels( void ) const { return _channels; }
template< class Real >
void SphericalStreamingSolver< Real >::SetIterations( int iters , int rSize )
{
	SetIterations( -iters * Degree - Degree + 1 , iters , 0 , 0 , 0 , 0 , rSize );
}
template< class Real >
void SphericalStreamingSolver< Real >::SetIterations( int start , int iters , int bStart , int xStart , int bEnd , int xEnd , int rSize )
{
	bool setResidual = this->setResidual;
	UnSetIterations();
	this->setResidual = setResidual;
	index = start;
	this->iters = iters;
	_iters = iters * Degree;

	xSize += (_iters+Degree) > xEnd ? (_iters+Degree) : xEnd;
	bSize += _iters > bEnd ? _iters : bEnd;
	if( setResidual )
	{
		xSize += -Degree-1 < xStart ? Degree+1 : -xStart;
		bSize +=        -1 < bStart ?        1 : -bStart;
	}
	else
	{
		xSize += -Degree < xStart ? Degree : -xStart;
		bSize +=       0 < bStart ?      0 : -bStart;
	}
	this->rSize = rSize;
	if( setResidual ) RStream = AllocArray< __m128 >( rSize * (_columns+2) * _channels , ALIGNMENT ) , memset( RStream , 0 , sizeof( __m128 ) * rSize * (_columns+2) * _channels );
	XStream = AllocArray< __m128 >( xSize * (_columns+2) * _channels , ALIGNMENT ) , memset( XStream , 0 , sizeof( __m128 ) * xSize * (_columns+2) * _channels );
	BStream = AllocArray< __m128 >( bSize * (_columns+2) * _channels , ALIGNMENT ) , memset( BStream , 0 , sizeof( __m128 ) * bSize * (_columns+2) * _channels );

	thetaStencils = new ThetaStencilData[bSize+2];

	laplacianStencils = ( SphericalLaplacianStencil* ) malloc( bSize * sizeof( SphericalLaplacianStencil ) );
	AllocAlignedMemory( _laplacianStencilsSSE             , laplacianStencilsSSE             , SphericalLaplacianStencilSSE , bSize , ALIGNMENT );
	AllocAlignedMemory( _zeroCenteredLaplacianStencilsSSE , zeroCenteredLaplacianStencilsSSE , SphericalLaplacianStencilSSE , bSize , ALIGNMENT );

	int d=0;
	while( (1<<d) <= columns ) d++;
	_phiStencilData = (PhiStencilData*) malloc( sizeof( PhiStencilData ) * d );
	for( int i=0 ; i<d ; i++ ) _phiStencilData[i].set( 1<<i , _stencilTable );


	_scratch0 = AllocArray< Real >( columns+4 ) , scratch0 = _scratch0 + 2;
	_scratch1 = AllocArray< Real >( columns+4 ) , scratch1 = _scratch1 + 2;
	_scratch2 = AllocArray< Real >( columns+4 ) , scratch2 = _scratch2 + 2;
}
template< class Real >
void SphericalStreamingSolver< Real >::UnSetIterations( void )
{
	index = 0;
	iters = _iters = 0;
	rSize = 0;
	setResidual = false;
	bSize = xSize = 0;

	FreeArray( RStream );
	FreeArray( XStream );
	FreeArray( BStream );

	if( thetaStencils ) delete[] thetaStencils , thetaStencils = NULL;

	MyFree( laplacianStencils );
	MyFree( _laplacianStencilsSSE );
	MyFree( _zeroCenteredLaplacianStencilsSSE );
	laplacianStencilsSSE = zeroCenteredLaplacianStencilsSSE = NULL;

	MyFree( _phiStencilData );

	FreeArray( _scratch0 );
	FreeArray( _scratch1 );
	FreeArray( _scratch2 );
}
template< class Real >
bool SphericalStreamingSolver< Real >::Increment( void )
{
	int idx;
	if( setResidual ) idx = index - Degree-1 + xSize;
	else              idx = index - Degree   + xSize;
	if( idx>=0 && idx<rows ) for( int c=0 ; c<_channels ; c++ ) memset( GetXRow( idx , c ) , 0 , sizeof( Real ) * _columns * RealPerWord );
	if( setResidual )idx = index-1 + bSize;
	else             idx = index   + bSize;
	if( idx>=0 && idx<rows ) for( int c=0 ; c<_channels ; c++ ) memset( GetBRow( idx , c ) , 0 , sizeof( Real ) * _columns * RealPerWord );
//	if( setResidual )idx = index-1 + rSize;
//	else             idx = index   + rSize;
//	if( idx>=0 && idx<rows ) for( int c=0 ; c<_channels ; c++ ) memset( GetRRow( idx , c ) , 0 , sizeof( Real ) * _columns * RealPerWord );

	index++;
	// Allow for the writing out of the B vector
	if( setResidual ) return ( index-1<rows);
	else              return ( index  <rows);
}
template< class Real >
template< class StorageType >
bool SphericalStreamingSolver< Real >::UpdateXInput( StreamingGrid* X , bool overwrite )
{
	int idx = index+_iters+Degree-1;
	// Read in from the X vector
	if( idx>=0 && idx<rows )
	{
		Pointer( StorageType ) xPtr = ( Pointer( StorageType ) )( *X )[ idx ];
		for( int c=0 ; c<_channels ; c++ )
		{
			int idx2 = columns*c;
			Pointer( Real ) xRow = GetXRow( idx , c );
			if( overwrite ) for( int jj=0 ; jj<columns ; jj++ ) xRow[jj]  = Real( xPtr[ jj+idx2 ] );
			else			for( int jj=0 ; jj<columns ; jj++ ) xRow[jj] += Real( xPtr[ jj+idx2 ] );
		}
		return true;
	}
	return false;
}
template< class Real >
template< class StorageType >
bool SphericalStreamingSolver< Real >::UpdateBInput( StreamingGrid* B , bool overwrite )
{
	int idx = index+_iters-1;
	// Read in from the B vector
	if( idx>=0 && idx<rows )
	{
		Pointer( StorageType ) bPtr = ( Pointer( StorageType ) )( *B )[ idx ];
		for( int c=0 ; c<_channels ; c++ )
		{
			int idx2 = columns*c;
			Pointer( Real ) bRow = GetBRow( idx , c );
			if( overwrite ) for( int jj=0 ; jj<columns ; jj++ ) bRow[jj]  = Real( bPtr[ jj+idx2 ] );
			else			for( int jj=0 ; jj<columns ; jj++ ) bRow[jj] += Real( bPtr[ jj+idx2 ] );
		}
		return true;
	}
	return false;
}
template< class Real >
template< class StorageType >
bool SphericalStreamingSolver< Real >::UpdateXOutput( StreamingGrid* X )
{
	// Copy the solution
	if( index>=0 && index<rows )
	{
		Pointer( StorageType ) xPtr = ( Pointer( StorageType ) )( *X )[ index ];
		for( int c=0 ; c<_channels ; c++ )
		{
			int idx2 = columns*c;
			Pointer( Real ) xRow = GetXRow( index , c );
			for( int jj=0 ; jj<columns ; jj++ ) xPtr[ jj+idx2 ] = StorageType( xRow[ jj ] );
		}
		return true;
	}
	return false;
}
template< class Real >
template< class StorageType >
bool SphericalStreamingSolver< Real >::UpdateBOutput( StreamingGrid* B )
{
#if MISHA_FIX
	int idx = index + _iters - 1;
#else // !MISHA_FIX
	int idx = index;
#endif // MISHA_FIX
	if( idx>=0 && idx<rows )
	{
		Pointer( StorageType ) bPtr = ( Pointer( StorageType ) )( *B )[ idx ];
		for( int c=0 ; c<_channels ; c++ )
		{
			int idx2 = columns*c;
			Pointer( Real ) bRow = GetBRow( idx , c );
			for( int jj=0 ; jj<columns ; jj++ ) bPtr[ jj+idx2 ] = StorageType( bRow[jj] );
		}
		return true;
	}
	return false;
}
#define SetUpSolveSSE( offset , previous )						\
	const int OFFSET = offset;									\
	__m128 mValues[] =											\
	{															\
		laplacianStencil.matrixValues[OFFSET][0],				\
		laplacianStencil.matrixValues[OFFSET][1],				\
		laplacianStencil.matrixValues[OFFSET][2],				\
		laplacianStencil.matrixValues[OFFSET][3],				\
		laplacianStencil.matrixValues[OFFSET][4]				\
	};															\
	SetFullInteriorDotSum( mValues , xPtrs , previous , dSum );


#define SSESolverUpdate( i , Real )									\
{																	\
	Real val = s2[i-2] + s1[i-1] + s0[i] + s1[i+1] + s2[i+2];		\
	val += (xPtr[(i-2+width)%width] + xPtr[(i+2+width)%width] ) * sValues[0];					\
	val += (xPtr[(i-1+width)%width] + xPtr[(i+1+width)%width] ) * sValues[1];					\
	xPtr[i] = bPtr[i]-val;										\
}
#define FastSSESolverUpdate( i , Real )								\
{																	\
	Real val = s2[i-2] + s1[i-1] + s0[i] + s1[i+1] + s2[i+2];		\
	val += (xPtr[i-2] + xPtr[i+2]) * sValues[0];					\
	val += (xPtr[i-1] + xPtr[i+1]) * sValues[1];					\
	xPtr[i] = bPtr[i]-val;											\
}
template< >
void SphericalStreamingSolver< double >::SolveSSE( int j , bool reverse )
{
	int width = RowWidth( j );
	int _width = width / RealPerWord;

	__m128d scratch128;
	double* scratch = (double*)&scratch128;
	SphericalLaplacianStencil stencil = GetLaplacianStencil( j );

	scratch[0] = stencil.evenMatrixValues[0][2];
	scratch[1] = stencil.evenMatrixValues[1][2];
	const __m128d stencil0a = _mm_load_pd( scratch );
	scratch[0] = stencil.evenMatrixValues[3][2];
	scratch[1] = stencil.evenMatrixValues[4][2];
	const __m128d stencil0b = _mm_load_pd( scratch );

	scratch[0] = stencil.evenMatrixValues[0][3];
	scratch[1] = stencil.evenMatrixValues[1][3];
	const __m128d stencil1a = _mm_load_pd( scratch );
	scratch[0] = stencil.evenMatrixValues[3][3];
	scratch[1] = stencil.evenMatrixValues[4][3];
	const __m128d stencil1b = _mm_load_pd( scratch );

	scratch[0] = stencil.evenMatrixValues[0][4];
	scratch[1] = stencil.evenMatrixValues[1][4];
	const __m128d stencil2a = _mm_load_pd( scratch );
	scratch[0] = stencil.evenMatrixValues[3][4];
	scratch[1] = stencil.evenMatrixValues[4][4];
	const __m128d stencil2b = _mm_load_pd( scratch );

	const double* sValues = stencil.evenMatrixValues[2];
	for( int c=0 ; c<_channels ; c++ )
	{
		ConstPointer( double ) bPtr  = GetBRow( j   , c );
		ConstPointer( double ) xPtr0 = GetXRow( j-2 , c );
		ConstPointer( double ) xPtr1 = GetXRow( j-1 , c );
		Pointer( double )      xPtr  = GetXRow( j   , c );
		ConstPointer( double ) xPtr3 = GetXRow( j+1 , c );
		ConstPointer( double ) xPtr4 = GetXRow( j+2 , c );

		__m128d temp128 , temp128A , temp128B;
		double* temp  = (double*)&temp128;
		double* tempA = (double*)&temp128A;
		double* tempB = (double*)&temp128B;
		for( int i=0 ; i<width ; i++ )
		{
			tempA[0] = xPtr0[i];
			tempA[1] = xPtr1[i];
			tempB[0] = xPtr3[i];
			tempB[1] = xPtr4[i];
			temp128 = _mm_add_pd( _mm_mul_pd( temp128A , stencil0a ) , _mm_mul_pd( temp128B , stencil0b ) ) ; scratch0[i] = temp[0] + temp[1];
			temp128 = _mm_add_pd( _mm_mul_pd( temp128A , stencil1a ) , _mm_mul_pd( temp128B , stencil1b ) ) ; scratch1[i] = temp[0] + temp[1];
			temp128 = _mm_add_pd( _mm_mul_pd( temp128A , stencil2a ) , _mm_mul_pd( temp128B , stencil2b ) ) ; scratch2[i] = temp[0] + temp[1];
		}
		scratch0[-1] = scratch0[width-1];
		scratch1[-1] = scratch1[width-1];
		scratch2[-1] = scratch2[width-1];
		scratch0[-2] = scratch0[width-2];
		scratch1[-2] = scratch1[width-2];
		scratch2[-2] = scratch2[width-2];
		scratch0[width  ] = scratch0[0];
		scratch1[width  ] = scratch1[0];
		scratch2[width  ] = scratch2[0];
		scratch0[width+1] = scratch0[1];
		scratch1[width+1] = scratch1[1];
		scratch2[width+1] = scratch2[1];

		ConstPointer( double ) s0 = scratch0;
		ConstPointer( double ) s1 = scratch1;
		ConstPointer( double ) s2 = scratch2;

		if( reverse )
		{
			SSESolverUpdate( width-1 , double );
			SSESolverUpdate( width-2 , double );
			for( int i=width-3 ; i>=2 ; i-- ) FastSSESolverUpdate( i , double );
			SSESolverUpdate( 1 , double );
			SSESolverUpdate( 0 , double );
		}
		else
		{
			SSESolverUpdate( 0 , double );
			SSESolverUpdate( 1 , double );
			for( int i=2 ; i<width-2 ; i++ ) FastSSESolverUpdate( i , double );
			SSESolverUpdate( width-2 , double );
			SSESolverUpdate( width-1 , double );
		}
	}
}
template< class Real >
void SphericalStreamingSolver< Real >::SolveSSE( int j , bool reverse )
{
	int width = RowWidth( j );
	int _width = width / RealPerWord;
#if NEW_SOLVE
	__m128 scratch128;
	float* scratch = (float*)&scratch128;
	SphericalLaplacianStencil stencil = GetLaplacianStencil( j );

	scratch[0] = stencil.evenMatrixValues[0][2];
	scratch[1] = stencil.evenMatrixValues[1][2];
	scratch[2] = stencil.evenMatrixValues[3][2];
	scratch[3] = stencil.evenMatrixValues[4][2];
	const __m128 stencil0 = _mm_load_ps( scratch );
	scratch[0] = stencil.evenMatrixValues[0][3];
	scratch[1] = stencil.evenMatrixValues[1][3];
	scratch[2] = stencil.evenMatrixValues[3][3];
	scratch[3] = stencil.evenMatrixValues[4][3];
	const __m128 stencil1 = _mm_load_ps( scratch );
	scratch[0] = stencil.evenMatrixValues[0][4];
	scratch[1] = stencil.evenMatrixValues[1][4];
	scratch[2] = stencil.evenMatrixValues[3][4];
	scratch[3] = stencil.evenMatrixValues[4][4];
	const __m128 stencil2 = _mm_load_ps( scratch );
	scratch[0] = stencil.evenMatrixValues[2][0];
	scratch[1] = stencil.evenMatrixValues[2][1];
	scratch[2] = stencil.evenMatrixValues[2][3];
	scratch[3] = stencil.evenMatrixValues[2][4];
	const __m128 stencil3 = _mm_load_ps( scratch );

	const float* sValues = stencil.evenMatrixValues[2];
	for( int c=0 ; c<_channels ; c++ )
	{
		ConstPointer( float ) bPtr  = GetBRow( j   , c );
		ConstPointer( float ) xPtr0 = GetXRow( j-2 , c );
		ConstPointer( float ) xPtr1 = GetXRow( j-1 , c );
		Pointer( float )      xPtr  = GetXRow( j   , c );
		ConstPointer( float ) xPtr3 = GetXRow( j+1 , c );
		ConstPointer( float ) xPtr4 = GetXRow( j+2 , c );

		for( int i=0 ; i<width ; i++ )
		{
			scratch[0] = xPtr0[i];
			scratch[1] = xPtr1[i];
			scratch[2] = xPtr3[i];
			scratch[3] = xPtr4[i];
			__m128 temp128;
			float* temp = (float*)&temp128;

			temp128 = _mm_mul_ps( scratch128 , stencil0 );
			scratch0[i] = temp[0] + temp[1] + temp[2] + temp[3];
			temp128 = _mm_mul_ps( scratch128 , stencil1);
			scratch1[i] = temp[0] + temp[1] + temp[2] + temp[3];
			temp128 = _mm_mul_ps( scratch128 , stencil2 );
			scratch2[i] = temp[0] + temp[1] + temp[2] + temp[3];
		}
		scratch0[-1] = scratch0[width-1];
		scratch1[-1] = scratch1[width-1];
		scratch2[-1] = scratch2[width-1];
		scratch0[-2] = scratch0[width-2];
		scratch1[-2] = scratch1[width-2];
		scratch2[-2] = scratch2[width-2];
		scratch0[width  ] = scratch0[0];
		scratch1[width  ] = scratch1[0];
		scratch2[width  ] = scratch2[0];
		scratch0[width+1] = scratch0[1];
		scratch1[width+1] = scratch1[1];
		scratch2[width+1] = scratch2[1];

		ConstPointer( float ) s0 = scratch0;
		ConstPointer( float ) s1 = scratch1;
		ConstPointer( float ) s2 = scratch2;

		if( reverse )
		{
			SSESolverUpdate( width-1 , float );
			SSESolverUpdate( width-2 , float );
			for( int i=width-3 ; i>=2 ; i-- ) FastSSESolverUpdate( i , float );
			SSESolverUpdate( 1 , float );
			SSESolverUpdate( 0 , float );
		}
		else
		{
			SSESolverUpdate( 0 , float );
			SSESolverUpdate( 1 , float );
			for( int i=2 ; i<width-2 ; i++ ) FastSSESolverUpdate( i , float );
			SSESolverUpdate( width-2 , float );
			SSESolverUpdate( width-1 , float );
		}
	}
#else // !NEW_SOLVE
	SphericalLaplacianStencilSSE& laplacianStencil = GetZeroCenteredLaplacianStencilSSE( j );
	for( int c=0 ; c<_channels ; c++ )
	{
		localBPtr = GetBRow( j , c );
		for( int yy = 0 ; yy <= 2*Degree ; yy++ ) localXPtrs[ yy ] = ( Pointer( __m128 ) ) GetXRow( j-Degree+yy , c );
#if 1
		float dotSum;
		__m128 dSum;
		float* scratch = (float*)&dSum;
		ConstPointer( __m128 ) xPtrs[] = { localXPtrs[0] , localXPtrs[1] , localXPtrs[2] , localXPtrs[3] , localXPtrs[4] };
		Pointer( __m128 ) _xPtrs[] = { localXPtrs[0] , localXPtrs[1] , localXPtrs[2] , localXPtrs[3] , localXPtrs[4] };

		if( reverse )
		{
			// The toroidiality of the row implies that for the second two passes, the front values need to be wrapped to the end
			for( int i=0 ; i<5 ; i++ ) localXPtrs[i][_width] = localXPtrs[i][0];
			{
				SetUpSolveSSE( 3 , 0 );
				dotSum = scratch[1] + scratch[2] + scratch[3];
				for( int i=OFFSET ; i<width ; i+=4 ) ( ( Pointer( float ) )_xPtrs[2] )[ i ] = GaussSeidelUpdate3( SetFullInteriorDotSum , mValues , xPtrs , localBPtr , dotSum , i );
			}
			{
				SetUpSolveSSE( 2 , 0 );
				dotSum = scratch[0] + scratch[1] + scratch[2] + scratch[3];
				for( int i=OFFSET ; i<width ; i+=4 ) ( ( Pointer( float ) )_xPtrs[2] )[ i ] = GaussSeidelUpdate2( SetFullInteriorDotSum , mValues , xPtrs , localBPtr , dotSum , i );
			}

			// The toroidiality of the row implies that for the first two passes the end values need to be wrapped to the front
			for( int i=0 ; i<5 ; i++ ) localXPtrs[i][-1] = localXPtrs[i][_width-1];
			{
				SetUpSolveSSE( 1  , -1 );
				dotSum = scratch[3];
				for( int i=OFFSET ; i<width ; i+=4 ) ( ( Pointer( float ) )_xPtrs[2] )[ i ] = GaussSeidelUpdate1( SetFullInteriorDotSum , mValues , xPtrs , localBPtr , dotSum , i );
			}
			{
				SetUpSolveSSE( 0 , -1 );
				dotSum = scratch[2] + scratch[3];
				for( int i=OFFSET ; i<width ; i+=4 ) ( ( Pointer( float ) )_xPtrs[2] )[ i ] = GaussSeidelUpdate0( SetFullInteriorDotSum , mValues , xPtrs , localBPtr , dotSum , i );
			}
		}
		else
		{
			// The toroidiality of the row implies that for the first two passes the end values need to be wrapped to the front
			for( int i=0 ; i<5 ; i++ ) localXPtrs[i][-1] = localXPtrs[i][_width-1];
			{
				SetUpSolveSSE( 0 , -1 );
				dotSum = scratch[2] + scratch[3];
				for( int i=OFFSET ; i<width ; i+=4 ) ( ( Pointer( float ) )_xPtrs[2] )[ i ] = GaussSeidelUpdate0( SetFullInteriorDotSum , mValues , xPtrs , localBPtr , dotSum , i );
			}
			{
				SetUpSolveSSE( 1  , -1 );
				dotSum = scratch[3];
				for( int i=OFFSET ; i<width ; i+=4 ) ( ( Pointer( float ) )_xPtrs[2] )[ i ] = GaussSeidelUpdate1( SetFullInteriorDotSum , mValues , xPtrs , localBPtr , dotSum , i );
			}

			// The toroidiality of the row implies that for the second two passes, the front values need to be wrapped to the end
			for( int i=0 ; i<5 ; i++ ) localXPtrs[i][_width] = localXPtrs[i][0];
			{
				SetUpSolveSSE( 2 , 0 );
				dotSum = scratch[0] + scratch[1] + scratch[2] + scratch[3];
				for( int i=OFFSET ; i<width ; i+=4 ) ( ( Pointer( float ) )_xPtrs[2] )[ i ] = GaussSeidelUpdate2( SetFullInteriorDotSum , mValues , xPtrs , localBPtr , dotSum , i );
			}
			{
				SetUpSolveSSE( 3 , 0 );
				dotSum = scratch[1] + scratch[2] + scratch[3];
				for( int i=OFFSET ; i<width ; i+=4 ) ( ( Pointer( float ) )_xPtrs[2] )[ i ] = GaussSeidelUpdate3( SetFullInteriorDotSum , mValues , xPtrs , localBPtr , dotSum , i );
			}
		}
#else
		SphericalLaplacianStencil& stencil = GetLaplacianStencil( j );
		float* xPtrs[] = { (float*)localXPtrs[0] , (float*)localXPtrs[1] , (float*)localXPtrs[2] , (float*)localXPtrs[3] , (float*)localXPtrs[4] };

		const int STRIDE = 1;
		for( int off=0 ; off<STRIDE ; off++ )
			if( reverse )
				for( int i=width-1-off ; i>=0 ; i-=STRIDE )
				{
					float val = 0;
					for( int j=-2 ; j<=2 ; j++ )
						for( int k=-2 ; k<=2 ; k++ )
							val += xPtrs[2+j][(i+k+width)%width] * stencil.evenMatrixValues[2+j][2+k];
					xPtrs[2][i] += ( localBPtr[i]-val );
				}
			else
				for( int i=off ; i<width ; i+=STRIDE )
				{
					float val = 0;
					for( int j=-2 ; j<=2 ; j++ )
						for( int k=-2 ; k<=2 ; k++ )
							val += xPtrs[2+j][(i+k+width)%width] * stencil.evenMatrixValues[2+j][2+k];
					xPtrs[2][i] += ( localBPtr[i]-val );
				}
#endif
	}
#endif // NEW_SOLVE
}

#define SetUpResidualSSE( offset , previous )					\
	const int OFFSET = offset;									\
	__m128 mValues[] =											\
	{															\
		laplacianStencil.matrixValues[OFFSET][0],				\
		laplacianStencil.matrixValues[OFFSET][1],				\
		laplacianStencil.matrixValues[OFFSET][2],				\
		laplacianStencil.matrixValues[OFFSET][3],				\
		laplacianStencil.matrixValues[OFFSET][4]				\
	};															\
	float diagonal = laplacianStencil.diagonal[OFFSET];			\
	SetFullInteriorDotSum( mValues , xPtrs , previous , dSum );

#define ResidualSSESolverUpdate( i , Real )							\
{																	\
	Real val = s2[i-2] + s1[i-1] + s0[i] + s1[i+1] + s2[i+2];		\
	val += (xPtr2[(i-2+width)%width] + xPtr2[(i+2+width)%width] ) * sValues[0];					\
	val += (xPtr2[(i-1+width)%width] + xPtr2[(i+1+width)%width] ) * sValues[1];					\
	val += xPtr2[i] * sValues[2];									\
	rPtr[i] = bPtr[i]-val;											\
}
#define FastResidualSSESolverUpdate( i , Real )						\
{																	\
	Real val = s2[i-2] + s1[i-1] + s0[i] + s1[i+1] + s2[i+2];		\
	val += (xPtr2[i-2] + xPtr2[i+2]) * sValues[0];					\
	val += (xPtr2[i-1] + xPtr2[i+1]) * sValues[1];					\
	val += xPtr2[i] * sValues[2];									\
	rPtr[i] = bPtr[i]-val;											\
}
template< >
void SphericalStreamingSolver< double >::SetResidualSSE( int j )
{
	int width = RowWidth( j );
	int _width = width / RealPerWord;

	__m128d scratch128;
	double* scratch = (double*)&scratch128;
	SphericalLaplacianStencil stencil = GetLaplacianStencil( j );

	scratch[0] = stencil.evenMatrixValues[0][2];
	scratch[1] = stencil.evenMatrixValues[1][2];
	const __m128d stencil0a = _mm_load_pd( scratch );
	scratch[0] = stencil.evenMatrixValues[3][2];
	scratch[1] = stencil.evenMatrixValues[4][2];
	const __m128d stencil0b = _mm_load_pd( scratch );

	scratch[0] = stencil.evenMatrixValues[0][3];
	scratch[1] = stencil.evenMatrixValues[1][3];
	const __m128d stencil1a = _mm_load_pd( scratch );
	scratch[0] = stencil.evenMatrixValues[3][3];
	scratch[1] = stencil.evenMatrixValues[4][3];
	const __m128d stencil1b = _mm_load_pd( scratch );

	scratch[0] = stencil.evenMatrixValues[0][4];
	scratch[1] = stencil.evenMatrixValues[1][4];
	const __m128d stencil2a = _mm_load_pd( scratch );
	scratch[0] = stencil.evenMatrixValues[3][4];
	scratch[1] = stencil.evenMatrixValues[4][4];
	const __m128d stencil2b = _mm_load_pd( scratch );

	const double* sValues = stencil.evenMatrixValues[2];
	for( int c=0 ; c<_channels ; c++ )
	{
		ConstPointer( double ) xPtr0 = GetXRow( j-2 , c );
		ConstPointer( double ) xPtr1 = GetXRow( j-1 , c );
		ConstPointer( double ) xPtr2 = GetXRow( j+0 , c );
		ConstPointer( double ) xPtr3 = GetXRow( j+1 , c );
		ConstPointer( double ) xPtr4 = GetXRow( j+2 , c );
		ConstPointer( double ) bPtr  = GetBRow( j   , c );
		Pointer( double )      rPtr  = GetRRow( j   , c );


		__m128d temp128 , temp128A , temp128B;
		double* temp  = (double*)&temp128;
		double* tempA = (double*)&temp128A;
		double* tempB = (double*)&temp128B;
		for( int i=0 ; i<width ; i++ )
		{
			tempA[0] = xPtr0[i];
			tempA[1] = xPtr1[i];
			tempB[0] = xPtr3[i];
			tempB[1] = xPtr4[i];
			temp128 = _mm_add_pd( _mm_mul_pd( temp128A , stencil0a ) , _mm_mul_pd( temp128B , stencil0b ) ) ; scratch0[i] = temp[0] + temp[1];
			temp128 = _mm_add_pd( _mm_mul_pd( temp128A , stencil1a ) , _mm_mul_pd( temp128B , stencil1b ) ) ; scratch1[i] = temp[0] + temp[1];
			temp128 = _mm_add_pd( _mm_mul_pd( temp128A , stencil2a ) , _mm_mul_pd( temp128B , stencil2b ) ) ; scratch2[i] = temp[0] + temp[1];
		}
		scratch0[-1] = scratch0[width-1];
		scratch1[-1] = scratch1[width-1];
		scratch2[-1] = scratch2[width-1];
		scratch0[-2] = scratch0[width-2];
		scratch1[-2] = scratch1[width-2];
		scratch2[-2] = scratch2[width-2];
		scratch0[width  ] = scratch0[0];
		scratch1[width  ] = scratch1[0];
		scratch2[width  ] = scratch2[0];
		scratch0[width+1] = scratch0[1];
		scratch1[width+1] = scratch1[1];
		scratch2[width+1] = scratch2[1];

		ConstPointer( double ) s0 = scratch0;
		ConstPointer( double ) s1 = scratch1;
		ConstPointer( double ) s2 = scratch2;

		ResidualSSESolverUpdate( 0 , double );
		ResidualSSESolverUpdate( 1 , double );
		for( int i=2 ; i<width-2 ; i++ ) FastResidualSSESolverUpdate( i , double );
		ResidualSSESolverUpdate( width-2 , double );
		ResidualSSESolverUpdate( width-1 , double );
	}
}
template< class Real >
void SphericalStreamingSolver< Real >::SetResidualSSE( int j )
{
	int width = RowWidth( j );
	int _width = width / RealPerWord;
	SphericalLaplacianStencilSSE& laplacianStencil = GetLaplacianStencilSSE( j );
	for( int c=0 ; c<_channels ; c++ )
	{
		localBPtr = GetBRow( j , c );
		localRPtr = GetRRow( j , c );

		for( int yy=0 ; yy<=2*Degree ; yy++ ) localXPtrs[yy] = ( Pointer( __m128 ) ) GetXRow( j-Degree+yy , c );
		float dotSum;
		__m128 dSum;
		float* scratch = (float*)&dSum;

		ConstPointer( __m128 ) xPtrs[] = { localXPtrs[0] , localXPtrs[1] , localXPtrs[2] , localXPtrs[3] , localXPtrs[4] };

		// The toroidiality of the row implies that for the first two passes the end values need to be wrapped to the front
		for( int i=0 ; i<5 ; i++ ) localXPtrs[i][-1] = localXPtrs[i][_width-1];

		{
			SetUpResidualSSE( 0 , -1 );
			dotSum = scratch[2] + scratch[3];
			for( int i=OFFSET ; i<width ; i+=4 ) localRPtr[i] = localBPtr[i]-GetLaplacianValue0( SetFullInteriorDotSum , mValues , xPtrs , dotSum , i );
		}
		{
			SetUpResidualSSE( 1 , -1 );
			dotSum = scratch[3];
			for( int i=OFFSET ; i<width ; i+=4 ) localRPtr[i] = localBPtr[i]-GetLaplacianValue1( SetFullInteriorDotSum , mValues , xPtrs , dotSum , i );
		}

		// The toroidiality of the row implies that for the second two passes, the front values need to be wrapped to the end
		for( int i=0 ; i<5 ; i++ ) localXPtrs[i][_width] = localXPtrs[i][0];

		{
			SetUpResidualSSE( 2 , 0 );
			dotSum = scratch[0] + scratch[1] + scratch[2] + scratch[3];
			for( int i=OFFSET ; i<width ; i+=4 ) localRPtr[i] = localBPtr[i]-GetLaplacianValue2( SetFullInteriorDotSum , mValues , xPtrs , dotSum , i );
		}
		{
			SetUpResidualSSE( 3 , 0 );
			dotSum = scratch[1] + scratch[2] + scratch[3];
			for( int i=OFFSET ; i<width ; i+=4 ) localRPtr[i] = localBPtr[i]-GetLaplacianValue3( SetFullInteriorDotSum , mValues , xPtrs , dotSum , i );
		}
	}
}

template< class Real >
void SphericalStreamingSolver< Real >::Solve( int j , bool reverse )
{
	SphericalLaplacianStencil& stencil = GetLaplacianStencil( j );
	int dimensions[5];
	for( int yy=-Degree ; yy<=Degree ; yy++ ) dimensions[Degree+yy] = RowWidth( j+yy );
	double val;
	Pointer( Real ) localXPtrs[2*Degree+1];

	for( int c=0 ; c<_channels ; c++ )
	{
		for( int yy=0 ; yy<=2*Degree ; yy++ )
			if( j-Degree+yy>=0 && j-Degree+yy<rows ) localXPtrs[yy] = GetXRow( j-Degree+yy , c );
			else								     localXPtrs[yy] = NullPointer< Real >( );

		localBPtr = GetBRow( j , c );

		if( j==0 || j==rows-1 ) // Deal with the poles separately
		{
			val = 0;
			for( int k=-Degree ; k<=Degree ; k++ )
				for( int l=0 ; l<dimensions[Degree+k] ; l++ ) val += stencil.evenMatrixValues[Degree+k][0] * localXPtrs[Degree+k][l];
			localXPtrs[ Degree ][0] += localBPtr[0] - val;
		}
		else
			if( reverse )
				for( int i=dimensions[Degree]-1 ; i>=0 ; i-- )
				{
					val = 0;
					for( int k=-Degree ; k<=Degree ; k++ )
					{
						int dim = dimensions[Degree+k];
						Pointer( Real ) xRow = localXPtrs[Degree+k];
						if     ( j+k< 0 || j+k>=rows   ) continue;
						else if( j+k==0 || j+k==rows-1 ) val+= stencil.evenMatrixValues[Degree+k][0] * xRow[0];
						else if( dimensions[Degree]==dim )
							for( int l=-Degree ; l<=Degree ; l++ ) val += stencil.evenMatrixValues[Degree+k][Degree+l] * xRow[ (i+l+dim) % dim ];
						else if( dimensions[Degree]<dim )
							for( int l=-3 ; l<=4 ; l++ ) val += stencil.evenMatrixValues[Degree+k][3+l] * xRow[ (i*2+l+dim) % dim ];
						else if( dimensions[Degree]>dim )
							if( i&1 )
								for( int l=-1 ; l<=2 ; l++ ) val += stencil.oddMatrixValues [Degree+k][1+l] * xRow[ (i/2+l+dim) % dim ];
							else
								for( int l=-2 ; l<=1 ; l++ ) val += stencil.evenMatrixValues[Degree+k][2+l] * xRow[ (i/2+l+dim) % dim ];
					}
					localXPtrs[ Degree ][i] += localBPtr[i] - val;
				}
			else
				for( int i=0 ; i<dimensions[Degree] ; i++ )
				{
					val = 0;
					for( int k=-Degree ; k<=Degree ; k++ )
					{
						int dim = dimensions[Degree+k];
						Pointer( Real ) xRow = localXPtrs[Degree+k];
						if     ( j+k< 0 || j+k>=rows   ) continue;
						else if( j+k==0 || j+k==rows-1 ) val+= stencil.evenMatrixValues[Degree+k][0] * xRow[0];
						else if( dimensions[Degree]==dim )
							for( int l=-Degree ; l<=Degree ; l++ ) val += stencil.evenMatrixValues[Degree+k][Degree+l] * xRow[ (i+l+dim) % dim ];
						else if( dimensions[Degree]<dim )
							for( int l=-3 ; l<=4 ; l++ ) val += stencil.evenMatrixValues[Degree+k][3+l] * xRow[ (i*2+l+dim) % dim ];
						else if( dimensions[Degree]>dim )
							if( i&1 )
								for( int l=-1 ; l<=2 ; l++ ) val += stencil.oddMatrixValues [Degree+k][1+l] * xRow[ (i/2+l+dim) % dim ];
							else
								for( int l=-2 ; l<=1 ; l++ ) val += stencil.evenMatrixValues[Degree+k][2+l] * xRow[ (i/2+l+dim) % dim ];
					}
					localXPtrs[ Degree ][i] += localBPtr[i] - val;
				}
	}
#if 0
if( j<=4 )
{
	double sum = 0;
	printf( "Stencil[%d/%d -- %d]:\n" , j , dimensions[Degree] , rows );
	for( int y=0 ; y<5 ; y++ )
	{
		for( int x=0 ; x<5 ; x++ ) printf( "\t%f" , stencil.evenMatrixValues[y][x] );
		printf( "\n" );
		if( j-2+y<0 ) continue;
		else if( j==0 ) sum = stencil.evenMatrixValues[2][0] + stencil.evenMatrixValues[3][0] * dimensions[Degree+1] + stencil.evenMatrixValues[4][0] * dimensions[Degree+2];
		else if( j-2+y==0 ) sum += stencil.evenMatrixValues[y][0];
		else if( dimensions[Degree]==dimensions[y] ) for( int x=0 ; x<5 ; x++ ) sum += stencil.evenMatrixValues[y][x];
		else if( dimensions[Degree]< dimensions[y] ) for( int x=0 ; x<8 ; x++ ) sum += stencil.evenMatrixValues[y][x];
		else if( dimensions[Degree]> dimensions[y] ) for( int x=0 ; x<4 ; x++ ) sum += stencil.evenMatrixValues[y][x];
	}
	printf( "Sum: %f\n" , sum );
	if( j<=2 ) for( int i=0 ; i<dimensions[Degree]  ; i++ ) printf( "[%d %d] %f %f\n" , j , i , localXPtrs[Degree][i] , localBPtr[i] );
}
if( j==80 ) exit( 0 );
#endif
}

template< class Real >
void SphericalStreamingSolver< Real >::SetResidual( int j )
{
	SphericalLaplacianStencil& stencil = GetLaplacianStencil( j );
	int dimensions[5];
	for( int yy=-Degree ; yy<=Degree ; yy++ ) dimensions[Degree+yy] = RowWidth( j+yy );
	double val;
	ConstPointer( Real ) localXPtrs[2*Degree+1];

	for( int c=0 ; c<_channels ; c++ )
	{
		for( int yy=0 ; yy<=2*Degree ; yy++ )
			if( j-Degree+yy>=0 && j-Degree+yy<rows ) localXPtrs[yy] = GetXRow( j-Degree+yy , c );
			else									 localXPtrs[yy] = NullPointer< Real >( );

		localBPtr = GetBRow( j , c );
		localRPtr = GetRRow( j , c );
		if( j==0 || j==rows-1 ) // Deal with the poles separately
		{
			val = 0;
			for( int k=-Degree ; k<=Degree ; k++ )
				for( int l=0 ; l<dimensions[k+Degree] ; l++ ) val += stencil.evenMatrixValues[Degree+k][0] * localXPtrs[Degree+k][l];
			localRPtr[0] = ( localBPtr[0] - val );
		}
		else
		{
			for( int i=0 ; i<dimensions[Degree] ; i++ )
			{
				val = 0;
				for( int k=-Degree ; k<=Degree ; k++ )
				{
					int dim = dimensions[Degree+k];
					ConstPointer( Real ) xRow = localXPtrs[Degree+k];
					if     ( j+k< 0 || j+k>=rows   ) continue;
					else if( j+k==0 || j+k==rows-1 ) val += stencil.evenMatrixValues[Degree+k][0] * xRow[0];
					else if( dimensions[Degree]==dim )
						for( int l=-Degree ; l<=Degree ; l++ ) val += stencil.evenMatrixValues[Degree+k][Degree+l] * xRow[ (i+l+dim) % dim ];
					else if( dimensions[Degree]<dim )
						for( int l=-3 ; l<=4 ; l++ ) val += stencil.evenMatrixValues[Degree+k][3+l] * xRow[ (i*2+l+dim) % dim ];
					else if( dimensions[Degree]>dim )
						if( i&1 )
							for( int l=-1 ; l<=2 ; l++ ) val += stencil.oddMatrixValues [Degree+k][1+l] * xRow[ (i/2+l+dim) % dim ];
						else
							for( int l=-2 ; l<=1 ; l++ ) val += stencil.evenMatrixValues[Degree+k][2+l] * xRow[ (i/2+l+dim) % dim ];
				}
				localRPtr[i] = ( localBPtr[i] - val );
			}
		}
	}
}
template< >
void SphericalStreamingSolver< float >::ScaleBRow( int idx , float scale )
{
	__declspec ( align( ALIGNMENT ) ) float scratch[4];
	scratch[0] = scratch[1] = scratch[2] = scratch[3] = scale;
	__m128 scratch128 = _mm_load_ps( scratch );
	for( int c=0 ; c<_channels ; c++ )
	{
		Pointer( __m128 ) bPtr = ( Pointer( __m128 ) ) GetBRow( idx , c );
		for( int i=0 ; i<_columns ; i++ ) bPtr[i] = _mm_mul_ps( bPtr[i] , scratch128 );
	}
}
template< >
void SphericalStreamingSolver< double >::ScaleBRow( int idx , double scale )
{
	__declspec ( align( ALIGNMENT ) ) double scratch[2];
	scratch[0] = scratch[1] = scale;
	__m128d scratch128 = _mm_load_pd( scratch );
	for( int c=0 ; c<_channels ; c++ )
	{
		Pointer( __m128d ) bPtr = ( Pointer( __m128d ) ) GetBRow( idx , c );
		for( int i=0 ; i<_columns ; i++ ) bPtr[i] = _mm_mul_pd( bPtr[i] , scratch128 );
	}
}

template< class Real >
void SphericalStreamingSolver< Real >::Solve( void )
{
	if( index+_iters-1>=0 && index+_iters-1<rows )
	{
		SetThetaStencil( index+_iters-1 );
		SphericalLaplacianStencil& stencil = GetLaplacianStencil( index+_iters-1 );
		SetStencil( stencil , index+_iters-1 , dWeight , lWeight );
		for( int i=-Degree ; i<=Degree ; i++ )
			for( int j=0 ; j<8 ; j++ )
			{
				stencil.evenMatrixValues[Degree+i][j] /= stencil.diagonal;
				stencil.oddMatrixValues[Degree+i][j]  /= stencil.diagonal;
			}
		StencilToStencilSSE( stencil , GetZeroCenteredLaplacianStencilSSE( index+_iters-1 ) , true  );
		ScaleBRow( index+_iters-1 , 1.0/stencil.diagonal );
	}
	// Solve the linear system
	bool reverse = true;
	for( int i = index+_iters-1 ; i>=index ; i-=Degree )
	{
		reverse = !reverse;
		if( i<0 || i>=rows ) continue;
#if USE_INTERIOR
		if( GetLaplacianStencil( i ).isRegular ) SolveSSE( i , reverse );
		else									 Solve   ( i , reverse );
#else // !USE_INTERIOR
		Solve( i , reverse );
#endif // USE_INTERIOR
	}

	if( index>=0 && index<rows )
	{
		SphericalLaplacianStencil& stencil = GetLaplacianStencil( index );
		for( int i=-Degree ; i<=Degree ; i++ )
			for( int j=0 ; j<8 ; j++ )
			{
				stencil.evenMatrixValues[Degree+i][j] *= stencil.diagonal;
				stencil.oddMatrixValues[Degree+i][j]  *= stencil.diagonal;
			}
		StencilToStencilSSE( stencil , GetLaplacianStencilSSE( index ) , false  );

		ScaleBRow( index , stencil.diagonal );
	}

	// Set the residual
	int idx = index-1;
	if( idx>=0 && idx<rows && setResidual )
	{
#if USE_INTERIOR
		if( GetLaplacianStencil( idx ).isRegular ) SetResidualSSE( idx );
		else									   SetResidual   ( idx );
#else // !USE_INTERIOR
		SetResidual( idx );
#endif // USE_INTERIOR
		for( int c=0 ; c<_channels ; c++ )
		{
			Pointer( Real ) localBPtr = GetBRow( idx , c );
			Pointer( Real ) localXPtr = GetXRow( idx , c );
			Pointer( Real ) localRPtr = GetRRow( idx , c );
			for( int i=0 ; i<RowWidth( idx ) ; i++ )
			{
				rSquareNorm[c] += double( localRPtr[i] ) * localRPtr[i];
				bSquareNorm[c] += double( localBPtr[i] ) * localBPtr[i];
				xSquareNorm[c] += double( localXPtr[i] ) * localXPtr[i];
			}
		}
	}
}
template< class Real >
template< class StorageType >
void SphericalStreamingSolver< Real >::Fullsolve( StreamingGrid* inX , StreamingGrid* inB , StreamingGrid* outX )
{
	while(1)
	{
		if( index-1>=rows )	return;
		int idx = index+_iters-1;
		if( idx+Degree>=0 && idx+Degree<rows && inX )
		{
			Real* inPtr = (Real*)(*inX)[idx+Degree];
			for( int c=0 ; c<_channels ; c++ )
			{
				int idx1 = c*columns;
				Real* x = GetXRow( idx+Degree , c );
				for( int i=0 ; i<RowWidth( idx+Degree ) ; i++ ) x[i] = inPtr[idx1+i];
			}
			inX->advance();
		}
		if(idx>=0)
		{
			if( idx<rows )
			{
				StorageType* inPtr = ( StorageType* )( *inB )[ idx ];
				for( int c=0 ; c<_channels ; c++ )
				{
					Real* bRow = GetBRow( idx , c );
					for( int i=0 ; i<RowWidth( idx ) ; i++ ) bRow[i] = Real( inPtr[c*columns+i] );
				}
				inB->advance();
			}
			// Run an interation of the Gauss-Seidel solver
			Solve();
		}
		if( outX && index>=0 && index<rows )
		{
			Real* outPtr=( Real* )( *outX )[ index ];
			for( int c=0 ; c<_channels ; c++ )
			{
				Real* inP = GetXRow( index , c );
				Real* outP = outPtr+c*columns;
				for( int i=0 ; i<RowWidth( index ) ; i++ ) outP[i] = inP[i];
			}
			outX->advance();
		}
		if( !Increment() ) return;
	}
}
///////////////////////////////////////
// MultiGridSphericalStreamingSolver //
///////////////////////////////////////
template< class Real , class StorageType >
int MultiGridSphericalStreamingSolver< Real , StorageType >::RowWidth( int idx ) const { return SphericalStreamingSolver< Real >::RowWidth( idx ); }
template< class Real , class StorageType >
Pointer( Real ) MultiGridSphericalStreamingSolver< Real , StorageType >::GetXRow( int idx , int c ) { return SphericalStreamingSolver< Real >::GetXRow( idx , c ); }

template< class Real , class StorageType >
MultiGridSphericalStreamingSolver< Real , StorageType >::MultiGridSphericalStreamingSolver( void ) : MultiGridRestrictionNode()
{
	_average = NULL;
	X = B = NULL;
	inX = outX = inB = outB = NULL;
	parent = NULL;


	/*
	_localRAccum = NULL;
	localRAccum  = NULL;
	_prolongationStencil2 = _prolongationStencil3 = NULL;
	prolongationStencil2 = NULL;
	prolongationStencil3 = NULL;
	*/


	int halfD = (Degree+1)>>1;
	int dual  = (Degree&1) ? 0 : 1;
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
	startRestriction = dual + halfD - 2*Degree + 1;

	// Restriction
	// To solve X[i] at depth d-1:
	//		Need to have B[i+iters-1] at depth d-1:
	//		Need to have R[2*(i+iters-1)-halfD,2*(i+iters-1)+halfD+dual]


	// Prolongation
	// Can start processing the child once:
	//		2*i-halfD >= 0
	//		i >= halfD/2
	startProlongation = (halfD+1)>>1;
	px2 = 2*startProlongation + halfD + dual;
	//		-iters-Degree+1
	prolongationOffset = px2 - (1-Degree);
	prolongationOffset++;
	restrictionBit = (halfD+dual+1)&1;
}
template< class Real , class StorageType >
void MultiGridSphericalStreamingSolver< Real , StorageType >::Initialize( int channels , int columns , int rows , int samples , const SphericalStencilTable< double >* stencilTable , double dWeight , double lWeight , bool memoryMappedFile )
{
	if( _average ) delete[] _average;
	_average = new double[ channels ];
	_clearAverage = (dWeight==0);
	if( X ) delete X , X = NULL;
	if( B ) delete B , B = NULL;

	if( memoryMappedFile )
	{
		inCore = false;
		B = new MultiStreamIOClient( columns * channels * sizeof(StorageType) , rows , STREAMING_GRID_BUFFER_MULTIPLIER , NULL , true );
		X = new MultiStreamIOClient( columns * channels * sizeof(StorageType) , rows , STREAMING_GRID_BUFFER_MULTIPLIER , NULL , true );
	}
	else
	{
		inCore = true;
		B = new MemoryBackedGrid( columns * channels * sizeof(Real) , rows , true );
		X = new MemoryBackedGrid( columns * channels * sizeof(Real) , rows , true );
	}
	Init( channels , columns , rows , samples , stencilTable , dWeight , lWeight );
	if( parent )
		if( IsDownSamplable( columns , rows ) )
//			parent->Initialize( columns>>1 , rows>>1 , samples , dWeight , lWeight , memoryMappedFile );
			( (MultiGridSphericalStreamingSolver*)parent )->Initialize( channels , columns>>1 , rows>>1 , samples , stencilTable , dWeight , lWeight , memoryMappedFile );
		else
			fprintf( stderr , "Is not down-samplable: %d %d\n" , columns , rows ) , exit( 0 );
}
template< class Real , class StorageType >
MultiGridSphericalStreamingSolver< Real , StorageType >::~MultiGridSphericalStreamingSolver( void )
{
	/*
	MyFree( _localRAccum );
	MyFree( _prolongationStencil2 );
	MyFree( _prolongationStencil3 );
	localRAccum          = NULL;
	prolongationStencil2 = NULL;
	prolongationStencil3 = NULL;
	*/

	if( X ) delete X , X = NULL;
	if( B ) delete B , B = NULL;
	if( _average ) delete[] _average;
}
template< class Real , class StorageType >
void MultiGridSphericalStreamingSolver< Real , StorageType >::InitProlongation( void )
{
memset( _average , 0 , sizeof( double ) * Channels() );
_weightSum = 0;
	if( parent ) parent->InitProlongation();
	if( inX )	 inX->reset ( true  , 1 );
	if( inB )	 inB->reset ( true  , 1 );
	if( outX )	 outX->reset( false , 1 );
	if( outB )	 outB->reset( false , 1 );
}
template< class Real , class StorageType >
void MultiGridSphericalStreamingSolver< Real , StorageType >::InitRestriction( void )
{
	memset( _average , 0 , sizeof( double ) * Channels() );
	_weightSum = 0;
	/*
	AllocAlignedMemory( _localRAccum , localRAccum , __m128 , _columns+2 , ALIGNMENT );
	AllocAlignedMemory( _prolongationStencil2 , prolongationStencil2 , ProlongationStencilSSE  , 3 * (2*Degree + 1 ) , ALIGNMENT );
	AllocAlignedMemory( _prolongationStencil3 , prolongationStencil3 , ProlongationStencilSSE2 , 3 , ALIGNMENT );
	*/
#if 0
		__declspec ( align( ALIGNMENT ) ) float scratch[4];
		for( int i=0 ; i<=2*Degree ; i++ ) // Iterate over the minor index in the mask
			for( int j=0 ; j<3 ; j++ )
				for( int k=0 ; k<Degree+2 ; k++ )
				{
					int jj;
					if(j==0)	jj=0;
					else		jj=Degree;
					scratch[0] = majorProlongationStencil.caseTable[jj].values[1]*minorProlongationStencil.caseTable[i].values[k];
					scratch[1] = majorProlongationStencil.caseTable[jj].values[2]*minorProlongationStencil.caseTable[i].values[k];
					scratch[2] = majorProlongationStencil.caseTable[jj].values[3]*minorProlongationStencil.caseTable[i].values[k];
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
#endif
	if( parent ) parent->InitRestriction();
	if( inX )	 inX->reset ( true  , 1 );
	if( inB )	 inB->reset ( true  , 1 );
	if( outX )	 outX->reset( false , 1 );
	if( outB )	 outB->reset( false , 1 );
}
#if MISHA_FIX
template< class Real , class StorageType >
void MultiGridSphericalStreamingSolver< Real , StorageType >::SetRestrictionIterations( int iters ) { SetRestrictionIterations( iters , true ); }
template< class Real , class StorageType >
void MultiGridSphericalStreamingSolver< Real , StorageType >::SetRestrictionIterations( int iters , bool offset )
{
	SphericalStreamingSolver< Real >::SetIterations( -Degree-iters*Degree+1 , iters , 0 , 0 , 0 , 0 , Degree+2 );
	// Set the iterations for the parent so that it trails accordingly
	if( parent ) ( (MultiGridSphericalStreamingSolver*) parent )->SetRestrictionIterations( iters , false );
	X->reset( false , xSize );
	B->reset( false , bSize );
	B->SetServer( &server );
	X->SetServer( &server );
	if( offset ) for( int d=0 ; d<Degree ; d++ ) IterateRestriction( );
}
#else // !MISHA_FIX
template< class Real , class StorageType >
void MultiGridSphericalStreamingSolver< Real , StorageType >::SetRestrictionIterations( int iters )
{
	SphericalStreamingSolver< Real >::SetIterations( -Degree-iters*Degree+1 , iters , 0 , 0 , 0 , 0 , Degree+2 );
	// Set the iterations for the parent so that it trails accordingly
	if( parent ) parent->SetRestrictionIterations( iters );
	X->reset( false , xSize );
	B->reset( false , bSize );
	B->SetServer( &server );
	X->SetServer( &server );
}
#endif // MISHA_FIX
template< class Real , class StorageType >
void MultiGridSphericalStreamingSolver< Real , StorageType >::SetProlongationIterations( int iters )
{
	SphericalStreamingSolver< Real >::SetIterations( -Degree-iters*Degree+1 , iters , 0 , 0 , 0 , prolongationOffset+iters*Degree );
	// Set the iterations for the child so that it trails accordingly
	if( child ) child->SetProlongationIterations( iters );
	X->reset( true , xSize );
	B->reset( true , bSize );
	B->SetServer( &server );
	X->SetServer( &server );
}
template< class Real , class StorageType >
void MultiGridSphericalStreamingSolver< Real , StorageType >::UnSetProlongationIterations( void )
{
	SphericalStreamingSolver< Real >::UnSetIterations( );
	X->unset( );
	B->unset( );
	if( child ) child->UnSetProlongationIterations( );
	MyFree( _scratch0 ) , MyFree( _scratch1 ) , MyFree( _scratch2 );
}
template< class Real , class StorageType >
void MultiGridSphericalStreamingSolver< Real , StorageType >::UnSetRestrictionIterations( void )
{
	SphericalStreamingSolver< Real >::UnSetIterations( );
	X->unset( );
	B->unset( );
	if( parent ) parent->UnSetRestrictionIterations( );
	MyFree( _scratch0 ) , MyFree( _scratch1 ) , MyFree( _scratch2 );
}
template< class Real , class StorageType >
bool MultiGridSphericalStreamingSolver< Real , StorageType >::IterateProlongation( void )
{
	if( _showProgress )
		if     ( index <rows ) printf( "Prolongation: [%.1f%%] [%d/%d]         \r" , float(index)/float(rows-1)*100 , index , rows );
		else if( index==rows ) printf( "                                              \r" );
	if( index-1>=rows ) return false;
	int idx = index+_iters+Degree-1;
	bool xUpdate = false;
	if( idx>=0 )
	{
		if( inCore )
		{
			if( UpdateBInput< Real >( B , true ) ) B->advance();					// -> index+_iters-1
			if( UpdateXInput< Real >( X , true ) ) X->advance() ,  xUpdate = true;	// -> index+_iters-1+Degree
		}
		else
		{
			if( UpdateBInput< StorageType >( B , true ) ) B->advance();
			if( UpdateXInput< StorageType >( X , true ) ) X->advance() , xUpdate = true;
		}
		/*
		if( xUpdate && _clearAverage && idx>=0 && idx<rows )
		{
			int width = RowWidth( idx );
			for( int c=0 ; c<Channels() ; c++ )
			{
				Pointer( Real ) xPtr = GetXRow( idx , c );
				for( int x=0 ; x<width ; x++ ) xPtr[x] -= _average[c];
			}
		}
		*/
		if( !parent && inX && idx>=0 && idx<rows )
		{
			Pointer( Real ) inPtr = ( Pointer( Real ) )( *inX )[ idx ];
			for( int c=0 ; c<Channels() ; c++ ) memcpy( GetXRow( idx , c ) , inPtr + c*columns , sizeof( Real ) * columns );

			inX->advance();
		}

		if( idx<rows && parent ) ProlongationUpdate( idx );

		// Run an iteration of the Gauss-Seidel solver
		SphericalStreamingSolver< Real >::Solve();

		if( outX && index>=0 && index<rows )
		{
			Pointer( Real ) outPtr = ( Pointer( Real ) )( *outX )[ index ];
			for( int c=0 ; c<Channels() ; c++ ) memcpy( outPtr + c*columns , GetXRow( index ,c  ) , sizeof( Real ) * columns );
			outX->advance();
		}
		if( index>=0 && index<rows )
		{
			int width = RowWidth( index );
			double weight = EquiRectangular< double >::BaseWeight( rows , width , index , _stencilTable ) / ( 4. * M_PI );
			_weightSum += weight * width;
			for( int c=0 ; c<Channels() ; c++ )
			{
				Pointer( Real ) xPtr = GetXRow( index , c );
				for( int x=0 ; x<width ; x++ ) _average[c] += xPtr[x] * weight;
			}
		}
	}
#if MISHA_FIX
	if( index>=0 && index<rows && outB )
	{
		Pointer( Real ) outPtr = ( Pointer( Real ) )( *outB )[ index ];
		for( int c=0 ; c<Channels() ; c++ ) memcpy( outPtr + c*columns , GetBRow( index , c ) , sizeof( Real ) * columns );
		outB->advance();
	}
#endif // MISHA_FIX
	if( child && index>=startProlongation ) if( child->IterateProlongation() ) child->IterateProlongation();

	return Increment();
}
template< class Real , class StorageType >
void MultiGridSphericalStreamingSolver< Real , StorageType >::ProlongationUpdate( int j )
{
	MultiGridSphericalStreamingSolver* parent = (MultiGridSphericalStreamingSolver*) this->parent;
	double pStencil[]  = { 1./ 4 , 3./ 4 , 3./ 4 , 1./ 4 };
	double p2Stencil[] = { 1./16 , 3./16 , 6./16 , 10./16 , 12./16 , 12./16 , 10./16 , 6./16 , 3./16 , 1./16 };
	double pStencilEven[4] , pStencilOdd[4];
	double p2StencilEven[10] , p2StencilOdd[10];
	int width = RowWidth( j );

	for( int i=0 ; i<4 ; i++ )
	{
		pStencilEven[i] = pStencil[i] * pStencil[3];
		pStencilOdd[i]  = pStencil[i] * pStencil[2];
	}
	for( int i=0 ;i<10 ; i++ )
	{
		p2StencilEven[i] = p2Stencil[i] * pStencil[3];
		p2StencilOdd[i]  = p2Stencil[i] * pStencil[2];
	}

	for( int c=0 ; c<Channels() ; c++ )
	{
		Pointer( Real ) localXPtr = GetXRow( j , c );
		if( j==0 )
		{
			ConstPointer( Real ) parentXPtr = parent->GetXRow( 0 , c );
			localXPtr[0] += parentXPtr[0] * ( pStencil[0] + pStencil[1] );
		}
		else if ( j==rows-1 )
		{
			ConstPointer( Real ) parentXPtr = parent->GetXRow( rows/2 - 1 , c );
			localXPtr[0] += parentXPtr[0] * ( pStencil[0] + pStencil[1] );
		}
		else
		{
			int startY = RestrictionStencil::RestrictionStencil::Start( j );
			ConstPointer( Real ) parentXPtrs[RestrictionStencil::RestrictionStencil::Size];
			int pWidths[RestrictionStencil::RestrictionStencil::Size]; 
			for( int yy=0 ; yy<RestrictionStencil::RestrictionStencil::Size ; yy++ ) parentXPtrs[yy] = parent->GetXRow( yy+startY , c ) , pWidths[yy] = parent->RowWidth( yy+startY );

			int yOff;
			if( j&1 ) yOff = 2;
			else      yOff = 0;

			for( int yy=0 ; yy<2 ; yy++ )
			{
				double *p , *p2;
				if( (j+yy)&1 ) p = pStencilOdd  , p2 = p2StencilOdd;
				else           p = pStencilEven , p2 = p2StencilEven;
				if( startY+yy==0 || startY+yy==rows/2-1 )
				{
					for( int i=0 ; i<width ; i++ ) localXPtr[i] += parentXPtrs[yy][0] * pStencil[yOff+yy];
				}
				// | | | | | | | | | | | | | | | | |
				else if( width==  pWidths[yy] )
				{
					for( int i=0 ; i<width ; i++ ) localXPtr[i] += parentXPtrs[yy][i] * pStencil[yOff+yy];
				}

				// | |+|+|+|+| | | | | | | | | | | |
				// |   |---|   |   |   |   |   |   |
				else if( width==2*pWidths[yy] )
				{
					{
						int i=0;
						int startX = RestrictionStencil::RestrictionStencil::Start( i );
						int i0 = ( startX+0 + pWidths[yy] ) % pWidths[yy];
						int i1 = ( startX+1 + pWidths[yy] ) % pWidths[yy];
						localXPtr[i] += parentXPtrs[yy][i0]*p[3] + parentXPtrs[yy][i1]*p[1];
					}
					for( int i=1 ; i<width-1 ; i+=2 )
					{
						int startX = RestrictionStencil::RestrictionStencil::Start( i );
						int i0 = startX+0;
						int i1 = startX+1;
						localXPtr[i+0] += parentXPtrs[yy][i0]*p[2] + parentXPtrs[yy][i1]*p[0];
						localXPtr[i+1] += parentXPtrs[yy][i0]*p[3] + parentXPtrs[yy][i1]*p[1];
					}
					{
						int i=width-1;
						int startX = RestrictionStencil::RestrictionStencil::Start( i );
						int i0 = ( startX+0 + pWidths[yy] ) % pWidths[yy];
						int i1 = ( startX+1 + pWidths[yy] ) % pWidths[yy];
						localXPtr[i] += parentXPtrs[yy][i0]*p[2] + parentXPtrs[yy][i1]*p[0];
					}
				}

				// | |*|*|*|*|*|*|*|*|*|*| | | | | |
				// |   |+++|+++|+++|+++|   |   |   |
				// |       |-------|       |       |
				else if( width==4*pWidths[yy] )
				{
					{
						int i=0;
						int startX = RestrictionStencil::RestrictionStencil::Start( RestrictionStencil::RestrictionStencil::Start( i ) );
						int i0 = ( startX+0 + pWidths[yy] ) % pWidths[yy];
						int i1 = ( startX+1 + pWidths[yy] ) % pWidths[yy];
						int i2 = ( startX+2 + pWidths[yy] ) % pWidths[yy];
						localXPtr[i+0] += parentXPtrs[yy][i0]*p2[7] + parentXPtrs[yy][i1]*p2[3];
						localXPtr[i+1] += parentXPtrs[yy][i0]*p2[8] + parentXPtrs[yy][i1]*p2[4] + parentXPtrs[yy][i2]*p2[0];
						localXPtr[i+2] += parentXPtrs[yy][i0]*p2[9] + parentXPtrs[yy][i1]*p2[5] + parentXPtrs[yy][i2]*p2[1];
					}
					for( int i=3 ; i<width-1 ; i+=4 )
					{
						int startX = RestrictionStencil::RestrictionStencil::Start( RestrictionStencil::RestrictionStencil::Start( i ) );
						int i0 = startX+0;
						int i1 = startX+1;
						int i2 = (startX+2) % pWidths[yy]; // for i \in [ width-4 , width )
//						int i2 = startX+2;
						localXPtr[i+0] += parentXPtrs[yy][i0]*p2[6] + parentXPtrs[yy][i1]*p2[2];
						localXPtr[i+1] += parentXPtrs[yy][i0]*p2[7] + parentXPtrs[yy][i1]*p2[3];
						localXPtr[i+2] += parentXPtrs[yy][i0]*p2[8] + parentXPtrs[yy][i1]*p2[4] + parentXPtrs[yy][i2]*p2[0];
						localXPtr[i+3] += parentXPtrs[yy][i0]*p2[9] + parentXPtrs[yy][i1]*p2[5] + parentXPtrs[yy][i2]*p2[1];
					}
					{
						int i = width-1;
						int startX = RestrictionStencil::RestrictionStencil::Start( RestrictionStencil::RestrictionStencil::Start( i ) );
						int i0 = ( startX+0 + pWidths[yy] ) % pWidths[yy];
						int i1 = ( startX+1 + pWidths[yy] ) % pWidths[yy];
						localXPtr[i] += parentXPtrs[yy][i0]*p2[6] + parentXPtrs[yy][i1]*p2[2];
					}
				}
			}
		}
	}
}
template< class Real , class StorageType >
void MultiGridSphericalStreamingSolver< Real , StorageType >::RestrictionUpdate( int j )
{
	double pStencil[]  = { 1./ 4 , 3./ 4 , 3./ 4 , 1./ 4 };
	double p2Stencil[] = { 1./16 , 3./16 , 6./16 , 10./16 , 12./16 , 12./16 , 10./16 , 6./16 , 3./16 , 1./16 };
	int width = RowWidth( j );

	for( int c=0 ; c<_channels ; c++ )
	{
		Real* localBPtr = GetBRow( j , c );
		if( j==0 )
		{
			const Real* childRPtr = child->GetRRow( 0 , c );
			localBPtr[0] += childRPtr[0] * ( pStencil[0] + pStencil[1] );

			for( int yy=1 ; yy<2 ; yy++ )
			childRPtr = child->GetRRow( yy , c );
			for( int i=0 ; i<child->RowWidth( yy ) ; i++ ) localBPtr[0] += childRPtr[i] * pStencil[yy+1];
		}
		else if ( jj==rows-1 );
		{
			const Real* childRPtr = parent->GetXRow( 2*rows - 1 , c );
			localBPtr[0] += childRPtr[0] * ( pStencil[0] + pStencil[1] );

			for( int yy=1 ; yy<2 ; yy++ )
			childRPtr = child->GetRRow( 2*rows-1-yy , c );
			for( int i=0 ; i<child->RowWidth( 2*rows-1-yy ) ; i++ ) localBPtr[0] += childRPtr[i] * pStencil[yy+1];
		}
		else
		{
			int startY = 2*j - 1;
			const Real* childRPtrs[4];
			int cWidths[4]; 
			for( int yy=0 ; yy<4 ; yy++ ) childRPtrs[yy] = child->GetRRow( yy+startY , c ) , cWidths[yy] = child->RowWidth( yy+startY );

			for( int i=0 ; i<width ; i++ )
			{
				double value = 0;

				for( int yy=0 ; yy<4 ; yy++ )
				{
					double tValue;

					// | | | | | | | | | | | | | | | | |
					if     (   width==cWidth[yy] ) tValue = childRPtrs[yy][i];

					// | |+|+|+|+| | | | | | | | | | | |
					// |   |---|   |   |   |   |   |   |
					else if( 2*width==cWidth[yy] )
						for( int xx=-1 ; xx<3 ; xx++ ) tValue += childRPtrs[yy][(2*i+xx+cWidths[yy])%cWidths[yy]] * pStencil[xx+1];

					// | |*|*|*|*|*|*|*|*|*|*| | | | | |
					// |   |+++|+++|+++|+++|   |   |   |
					// |       |-------|       |       |
					else if( 4*width==pWidth[yy] )
						for( int xx=-3 ; xx<7 ; xx++ ) tValue += childRPtrs[yy][(4*i+xx+cWidths[yy])%cWidths[yy]] * pStencil2[xx+3];

					value += tValue * pStencil[yy];
				}
				localBPtr[i] += value;
			}
		}
	}
}
template< class Real , class StorageType >
bool MultiGridSphericalStreamingSolver< Real , StorageType >::IterateRestriction( void )
{
	if( _showProgress )
		if     ( index<rows  ) printf( "Restriction: [%.1f%%] [%d/%d]         \r" , float(index)/float(rows-1)*100 , index , rows );
		else if( index==rows ) printf( "                                              \r" );
//if( rows==4096 )
//printf( "Iterating %d / %d (%d) %f %f %f\n" , index , rows , samples , bSquareNorm[0] + bSquareNorm[1] + bSquareNorm[2] , xSquareNorm[0] + xSquareNorm[1] + xSquareNorm[2] , rSquareNorm[0] + rSquareNorm[1] + rSquareNorm[2] );
	if( index-1>=rows ) return false;
	int idx = index+_iters-1;
	if( idx+Degree>=0 && idx+Degree<rows && inX )
	{
		Pointer( Real ) inPtr = ( Pointer( Real ) )(*inX)[idx+Degree];
		for( int c=0 ; c<Channels() ; c++ ) memcpy( GetXRow( idx+Degree , c) , inPtr + c*columns , sizeof( Real ) * columns );
		inX->advance();
	}
	if( idx>=0 && idx<rows && inB )
	{
		Pointer( Real ) inPtr = ( Pointer( Real ) )(*inB)[idx];
		for( int c=0 ; c<Channels() ; c++ ) memcpy( GetBRow( idx , c ) , inPtr + c*columns , sizeof( Real ) * columns );
		inB->advance();
	}
	if( idx>=0 && idx<rows )
		{
			if( !inB )
				if( !child ) fprintf( stderr , "Badness: no input stream\n" ) , exit(0);
				else
				{
					Pointer( Real )* rows = new Pointer( Real )[ Channels() ];
					for( int c=0 ; c<Channels() ;  c++ ) rows[c] = GetBRow( idx , c );
					child->SetRestriction( rows , idx , RowWidth( idx ) );
					delete[] rows;
				}
#if MISHA_FIX
			if( outB )
			{
				Pointer( Real ) outPtr = ( Pointer( Real ) )( *outB )[ idx ];
				for( int c=0 ; c<Channels() ; c++ ) memcpy( outPtr + c*columns , GetBRow( idx , c ) , sizeof( Real ) * columns );
				outB->advance();
			}
#endif // MISHA_FIX
	}
#if MISHA_FIX
	if( inCore ){ if( UpdateBOutput< Real        >( B ) ) B->advance(); }
	else        { if( UpdateBOutput< StorageType >( B ) ) B->advance(); }
#endif // MISHA_FIX
	// Run an iteration of the Gauss-Seidel solver
	if( idx>=0 ) SphericalStreamingSolver< Real >::Solve();
#if MISHA_FIX
	if( inCore ){ if( UpdateXOutput< Real        >( X ) ) X->advance(); }
	else        { if( UpdateXOutput< StorageType >( X ) ) X->advance(); }
#else // !MISHA_FIX
	// Write out the current solution
	if( inCore )
	{
		if( UpdateXOutput< Real >( X ) ) X->advance();	// index
		if( UpdateBOutput< Real >( B ) ) B->advance();	// index+_iters-1
	}
	else
	{
		if( UpdateXOutput< StorageType >( X ) ) X->advance();
		if( UpdateBOutput< StorageType >( B ) ) B->advance();
	}
#endif // MISHA_FIX

#if !MISHA_FIX
	if( outB && index>=0 && index<rows )
	{
		Pointer( Real ) outPtr = ( Pointer( Real ) )( *outB )[ index ];
		for( int c=0 ; c<Channels() ; c++ ) memcpy( outPtr + c*columns , GetBRow( index , c ) , sizeof( Real ) * columns );
		outB->advance();
	}
#endif // !MISHA_FIX
	if( outX && index>=0 && index<rows )
	{
		Pointer( Real ) outPtr = ( Pointer( Real ) )( *outX )[ index ];
		for( int c=0 ; c<Channels() ; c++ ) memcpy( outPtr + c*columns , GetXRow( index , c ) , sizeof( Real ) * columns );
		outX->advance();
	}
	if( index>=0 && index<rows )
	{
		int width = RowWidth( index );
		double weight = EquiRectangular< double >::BaseWeight( rows , width , index , _stencilTable ) / ( 4. * M_PI );
		_weightSum += weight * width;
		for( int c=0 ; c<Channels() ; c++ )
		{
			Pointer( Real ) xPtr = GetXRow( index , c );
			for( int x=0 ; x<width ; x++ ) _average[c] += xPtr[x] * weight;
		}
	}
	/*
	if( index==rows )
	{
		printf( "Average %d x %d:" , columns , rows );
		for( int c=0 ; c<Channels() ; c++ ) printf( " %f" , _average[c] );
		printf( "\t%f\n" , _weightSum );
	}
	*/
	if( parent ) if(index>=startRestriction && (index&1)==restrictionBit)	parent->IterateRestriction();
	return Increment();
}

template< class Real , class StorageType >
void MultiGridSphericalStreamingSolver< Real , StorageType >::SetRestriction( Pointer( Real )* lB , int idx , int lWidth )
{
	int halfD  = (Degree+1)>>1;
	int startY = 2*idx-halfD;
	double pStencil[]  = { 1./4 , 3./4 , 3./4 , 1./4 };
	double p2Stencil[] = { 1./16 , 3./16 , 6./16 , 10./16 , 12./16 , 12./16 , 10./16 , 6./16 , 3./16 , 1./16 };
	int dimensions[Degree+2];
	Pointer( Real ) localRPtrs[Degree+2];

	for( int yy=0 ; yy<Degree+2 ; yy++ ) dimensions[yy] = RowWidth( startY+yy );
	for( int c=0 ; c<Channels() ; c++ )
	{
		for( int yy=0 ; yy<Degree+2 ; yy++ ) localRPtrs[yy] = ( Pointer( Real ) )GetRRow( startY+yy , c );

		if( idx==0 )
		{
			lB[c][0] += localRPtrs[1][0] * ( pStencil[0] + pStencil[1] );
			for( int yy=2 ; yy<4 ; yy++ ) for( int xx=0 ; xx<dimensions[yy] ; xx++ ) lB[c][0] += localRPtrs[yy][xx] * pStencil[yy];
		}
		else if (idx==rows/2-1 )
		{
			for( int yy=0 ; yy<2 ; yy++ ) for( int xx=0 ; xx<dimensions[yy] ; xx++ ) lB[c][0] += localRPtrs[yy][xx] * pStencil[yy];
			lB[c][0] += localRPtrs[2][0] * ( pStencil[0] + pStencil[1] );
		}
		else
		{
			for( int x=0 ; x<lWidth ; x++ )
			{
				double val = 0;
				for( int yy=0 ; yy<4 ; yy++ )
				{
					int dim = dimensions[yy];
					double tVal = 0;
					if     ( dim==  lWidth ) tVal += localRPtrs[yy][x];
					else if( dim==2*lWidth ) for( int l=-1 ; l<=2 ; l++ ) tVal += localRPtrs[yy][(2*x+l+dim) % dim ] * pStencil[l+1];
					else if( dim==4*lWidth ) for( int l=-3 ; l<=6 ; l++ ) tVal += localRPtrs[yy][(4*x+l+dim) % dim ] * p2Stencil[l+3];
					val += tVal * pStencil[yy];
				}
				lB[c][x] += val;
			}
		}
	}
}
#if 0
template<class Real,int Type,int Degree,int Channels,class StorageType>
void MultiGridSphericalStreamingSolver<Real,Type,Degree,Channels,StorageType>::SetInteriorRestriction(Real* lB,int idx,int major2,int minor2)
{
	int halfD=(Degree+1)>>1;
	int startY=2*idx-halfD;
	for(int c=0;c<Channels;c++)
	{
		Real* myLB=&lB[c*major2];
		for(int yy=0;yy<Degree+2;yy++)	localRPtrs[yy]=(WordClass*)(localR + ((startY+yy)%rSize)*_major*RealPerWord*Channels+c*_major*RealPerWord);
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
		int jj;
		if		(idx<Degree)			jj=idx;
		else if	(idx>minor-1-Degree)	jj=2*Degree+(idx-(minor-1));
		else							jj=Degree;
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
	}
}
#endif
template< class Real , class StorageType >
void MultiGridSphericalStreamingSolver< Real , StorageType >::SolveProlongation( void )
{
	// Run to completion...
	while( IterateProlongation( ) ){;}
	// ...and finish up the trailing child
	if( child ) child->SolveProlongation( );

	while( ( B && B->server ) || ( X && X->server ) ) Sleep( 0 );
}
template< class Real , class StorageType >
void MultiGridSphericalStreamingSolver< Real , StorageType >::SolveRestriction( void )
{
	// Run to completion...
	while( IterateRestriction( ) ){;}
	// ...and finish up the trailing parent
	if(parent)	parent->SolveRestriction();

	while( ( B && B->server ) || ( X && X->server ) ) Sleep( 0 );
}
template< class Real , class StorageType >
bool MultiGridSphericalStreamingSolver< Real , StorageType >::IsDownSamplable( const int& hMajor , const int& hMinor , int &lMajor , int& lMinor )
{
	return FiniteElements2D< Real , ZERO_DERIVATIVE , Degree >::IsDownSamplable( hMajor , hMinor , lMajor , lMinor );
}
template< class Real , class StorageType >
bool MultiGridSphericalStreamingSolver< Real , StorageType >::IsDownSamplable( const int& hMajor , const int& hMinor )
{
	int lMajor , lMinor;
	return IsDownSamplable( hMajor , hMinor , lMajor , lMinor );
}
///////////////////////////////////////////
// StreamingAdaptiveNonAdaptiveConverter //
///////////////////////////////////////////
template< class Real , int Channels >
StreamingAdaptiveNonAdaptiveConverter< Real , Channels >::StreamingAdaptiveNonAdaptiveConverter( void ) : MultiGridRestrictionNode()
{
	outX = NULL;
	dataRow = NullPointer< Real >( );
}
template< class Real , int Channels >
StreamingAdaptiveNonAdaptiveConverter< Real , Channels >::~StreamingAdaptiveNonAdaptiveConverter(void)
{
	FreeArray( dataRow );
}
template< class Real , int Channels >
Pointer( Real ) StreamingAdaptiveNonAdaptiveConverter< Real , Channels >::GetRow( int row , int channel )
{
	return dataRow + channel * columns;
}
template< class Real , int Channels >
void StreamingAdaptiveNonAdaptiveConverter< Real , Channels >::Initialize( int rows )
{
	_sphere = AdaptiveEquiRectangular< Real >( rows , 1000 );
	columns = 2*rows;
	FreeArray( dataRow );
	dataRow = AllocArray< Real >( Channels * columns );

	this->columns = columns;
	this->rows = rows;

	int logColumns = 0;
	while( (1<<logColumns) < columns ) logColumns++;
	pStencils.resize( logColumns );
	start.resize( logColumns );
	for( int i=0 ; i<logColumns ; i++ )
	{
		if( i==0 )
		{
			pStencils[i].resize( 1 );
			pStencils[i][0] = 1.;
			start[i] = 0;
		}
		else
		{
			int dim = 2+pStencils[i-1].size()*2;
			if( dim>=columns ) dim = columns;
			start[i] = 2*start[i-1] - 1;
			if( start[i] <= -columns ) start[i] += columns;
			pStencils[i].resize( dim );
			for( int j=0 ; j<dim ; j++ ) pStencils[i][j] = 0;
			for( int j=0 ; j<pStencils[i-1].size() ; j++ )
			{
				pStencils[i][ (2*j  ) % dim ] += 1. * pStencils[i-1][j];
				pStencils[i][ (2*j+1) % dim ] += 3. * pStencils[i-1][j];
				pStencils[i][ (2*j+2) % dim ] += 3. * pStencils[i-1][j];
				pStencils[i][ (2*j+3) % dim ] += 1. * pStencils[i-1][j];
			}
		}
	}
	for( int i=0 ; i<logColumns ; i++ )
		for( int j=0 ; j<pStencils[i].size() ; j++ )
			pStencils[i][j] /= (1 << (2*i) );
}
template< class Real , int Channels >
int StreamingAdaptiveNonAdaptiveConverter< Real , Channels >::RowWidth( int idx ) const { return columns; }
template< class Real , int Channels >
void StreamingAdaptiveNonAdaptiveConverter< Real , Channels >::InitRestriction( void )
{
	index = 0;
	squareNorm = 0;
	if( parent ) parent->InitRestriction();
}
template< class Real , int Channels >
void StreamingAdaptiveNonAdaptiveConverter< Real , Channels >::InitProlongation( void )
{
	index = 0;
	squareNorm = 0;
	if( parent ) parent->InitProlongation();
}
template< class Real , int Channels >
void StreamingAdaptiveNonAdaptiveConverter< Real , Channels >::SetRestrictionIterations( int iters )
{
	if( parent ) parent->SetRestrictionIterations( iters );
}
template< class Real , int Channels >
void StreamingAdaptiveNonAdaptiveConverter< Real , Channels >::SetProlongationIterations( int iters )
{
	if( child ) child->SetProlongationIterations( iters );
}
template< class Real , int Channels >
void StreamingAdaptiveNonAdaptiveConverter< Real , Channels >::SetRestriction( Pointer( Real )* lB , int idx , int lWidth )
{
	int width = _sphere.dimension( idx );
	int logWidth = 0;
	while( (width<<logWidth) < columns ) logWidth++;
	if( idx>=0 && idx<rows )
		for( int c=0 ; c<Channels ; c++ )
		{
			Pointer( Real ) _row = GetRow( idx , c );
			if( idx==0 || idx==rows-1 )
			{
				double temp = 0;
#if 1
				lB[c][0] = _row[0];
#else
				for( int i=0 ; i<columns ; i++ ) temp += _row[i];
				lB[c][0] += temp / columns;
#endif
			}
			else
			{
				for( int i=0 ; i<width ; i++ )
				{
					double temp = 0;
					for( int j=0 ; j<pStencils[logWidth].size() ; j++ )
						temp += pStencils[logWidth][j] * _row[ (start[logWidth] + j + i * (1<<logWidth ) + columns) % columns ];
					lB[c][i] += temp;
				}
			}

		}
}
template< class Real , int Channels >
void StreamingAdaptiveNonAdaptiveConverter< Real , Channels >::ProlongationUpdate( int idx )
{
	int lWidth = _sphere.dimension( idx );
	int logWidth = 0;
	while( (lWidth<<logWidth) < columns ) logWidth++;
	if( idx>=0 && idx<rows )
		for( int c=0 ; c<Channels ; c++ )
		{
			Pointer( Real ) _row = GetRow( idx , c );
			memset( _row , 0 , sizeof( Real ) * columns );
			Pointer( Real ) pRow = parent->GetXRow( idx , c );
#if 0
			if( idx==0 || idx==rows-1 ) _row[0] += pRow[0];
#else
			if( idx==0 || idx==rows-1 ) for( int i=0 ; i<columns ; i++ ) _row[i] += pRow[0];
#endif
			else
				for( int i=0 ; i<lWidth ; i++ )
					for( int j=0 ; j<pStencils[logWidth].size() ; j++ )
						_row[ (start[logWidth] + j + i * (1<<logWidth ) + columns) % columns ] += pRow[i] * pStencils[logWidth][j];
		}
}
template< class Real , int Channels >
bool StreamingAdaptiveNonAdaptiveConverter< Real , Channels >::IterateRestriction( void )
{
	if( _showProgress )
		if     ( index<rows  ) printf( "Restriction: [%.1f%%] [%d/%d]         \r" , float(index)/float(rows-1)*100 , index , rows );
		else if( index==rows ) printf( "                                              \r" );
	if( index>=rows ) return false;
	if( index>=0 )
	{
		if( child )
		{
			Pointer( Real ) rows[Channels];
			for( int c=0 ; c<Channels ; c++ ) rows[c] = GetRow( index , c );
			child->SetRestriction( rows , index , columns );
			for( int c=0 ; c<Channels ; c++ ) for( int i=0 ; i<columns ; i++ ) squareNorm += rows[c][i] * rows[c][i];
		}
		if( parent ) parent->IterateRestriction();
	}
	index++;
	return true;
}
template< class Real , int Channels >
bool StreamingAdaptiveNonAdaptiveConverter< Real , Channels >::IterateProlongation( void )
{
	// This is a total hack (to address the fact that the multigridrestrictionnode prolongs twice...)
	int idx = index>>1;
	if( _showProgress )
		if     ( idx <rows ) printf( "Prolongation: [%.1f%%] [%d/%d]         \r" , float(idx)/float(rows-1)*100 , idx , rows );
		else if( idx==rows ) printf( "                                              \r" );

	if( idx>=rows ) return false;
	if( idx>=0 && index&1 )
	{
		if( parent )
		{
			ProlongationUpdate( idx );
			if( outX )
			{
				Pointer( Real ) row = ( Pointer( Real ) )(*outX)[idx];
				for( int c=0 ; c<Channels ; c++ )
				{
					Pointer( Real ) _row = GetRow( idx , c );
					if( outX->SeparateColors( ) ) for( int i=0 ; i<columns ; i++ ) row[c*columns +i] = _row[i];
					else                          for( int i=0 ; i<columns ; i++ ) row[i*Channels+c] = _row[i];
				}
				outX->advance( );
			}
		}
		if( child ) child->IterateProlongation();
	}
	index++;
	return true;
}
template< class Real , int Channels >
void StreamingAdaptiveNonAdaptiveConverter< Real , Channels >::SolveRestriction( void )
{
	while( IterateRestriction( ) ){;}
	if( parent ) parent->SolveRestriction();
}
template< class Real , int Channels >
void StreamingAdaptiveNonAdaptiveConverter< Real , Channels >::SolveProlongation( void )
{
	while( IterateProlongation( ) ){;}
	if( child ) child->SolveProlongation();
}


///////////////////////////////////
// StreamingNonAdaptiveLaplacian //
///////////////////////////////////
template< class Real , int PixelChannels , int LabelChannels , class PixelType , class LabelType >
StreamingNonAdaptiveLaplacian< Real , PixelChannels , LabelChannels , PixelType , LabelType >::StreamingNonAdaptiveLaplacian( void )
{
	pixels = labels = NULL;

	dSize = 2*Degree+1;
	size  = 2*Degree+3;

	dXSquareNorm = 0;
	dYSquareNorm = 0;
	outputSquareNorm = 0;

	localDMajorAccum = localDMinorAccum = NullPointer< __m128 >( );
	localDMajor = localDMinor = NullPointer< Real >( );
	_pixelRow = _previousPixelRow = NullPointer< PixelType >( );
	_labelRow = _previousLabelRow = NullPointer< LabelType >( );
	_dx = _dy = NullPointer< Real >( );
	_stencilTable = NULL;
}
template< class Real , int PixelChannels , int LabelChannels , class PixelType , class LabelType >
StreamingNonAdaptiveLaplacian< Real , PixelChannels , LabelChannels , PixelType , LabelType >::~StreamingNonAdaptiveLaplacian(void)
{
	FreeArray( localDMajorAccum );
	FreeArray( localDMinorAccum );
	FreeArray( localDMajor );
	FreeArray( localDMinor );
	FreeArray( _previousPixelRow );
	FreeArray( _previousLabelRow );
	FreeArray( _pixelRow );
	FreeArray( _labelRow );
	FreeArray( _dx );
	FreeArray( _dy );
}
template< class Real , int PixelChannels , int LabelChannels , class PixelType , class LabelType >
bool StreamingNonAdaptiveLaplacian< Real , PixelChannels , LabelChannels , PixelType , LabelType >::_readNextRow( void )
{
	int _channels = PixelChannels;
	if( pixels->HasAlpha() ) _channels++;
	_currentRow++;
	if( _currentRow>=minor ) return false;
	double avg[PixelChannels];
	for( int c=0 ; c<PixelChannels ; c++ ) avg[c] = 0;
	double weight = EquiRectangular< double >::BaseWeight( minor , major , _currentRow , _stencilTable );
	if( _currentRow==0 || _currentRow==minor-1 ) weight /= major;
	Pointer( PixelType ) row = ( Pointer( PixelType ) )(*pixels)[_currentRow];
	for( int c=0 ; c<_channels ; c++ )
		if( pixels->SeparateColors() ) memcpy( _pixelRow+c*major , row+c*major , sizeof(PixelType) * major );
		else                           for( int i=0 ; i<major ; i++ ) _pixelRow[i+c*major] = row[_channels*i+c];
	pixels->advance( );

	if( labels )
	{
		Pointer( LabelType ) row = ( Pointer( LabelType ) )(*labels)[_currentRow];
		for( int c=0 ; c<LabelChannels ; c++ ) 
			if( labels->SeparateColors() ) memcpy( _labelRow+c*major , row+c*major , sizeof(LabelType) * major );
			else                           for( int i=0 ; i<major ; i++ ) _labelRow[i+c*major] = row[LabelChannels*i+c];
		labels->advance( );
		// Force unkown pixel values to black
		if( _unknownIndex )
			for( int i=0 ; i<major ; i++ )
			{
				bool isUnknown=true;
				for( int c=0 ; c<LabelChannels ; c++ ) if( _labelRow[i+c*major]!=_unknownIndex[c] ) isUnknown = false;
				if( isUnknown ) for( int c=0 ; c<PixelChannels ; c++ ) _pixelRow[i+c*major] = 0;
			}
	}

	if( _vWeight )
	{
		for( int c=0 ; c<PixelChannels ; c++ )
		{
			Pointer( Real ) vPtr = _values + (_currentRow%size) * (_major+2) * PixelChannels * RealPerWord + (_major+2) * c * RealPerWord + (Degree+1) * RealPerWord;
			memcpy( vPtr , _pixelRow + c*major  , sizeof( Real ) * major );
			for( int i=0 ; i<RealPerWord ; i++ ) vPtr[-RealPerWord+i] = vPtr[major-RealPerWord+i];
			for( int i=0 ; i<RealPerWord ; i++ ) vPtr[major+i] = vPtr[i];
		}
	}
	for( int i=0 ; i<major ; i++ ) for( int c=0 ; c<PixelChannels ; c++ ) avg[c] += double( _pixelRow[i+major*c] ) * weight;
	for( int c=0 ; c<PixelChannels ; c++ ) average[c] += avg[c];
	if( _currentRow==minor-1 ) for( int c=0 ; c<PixelChannels ; c++ ) average[c] /= 4. * M_PI;
	return true;
}
template< class Real , int PixelChannels , int LabelChannels , class PixelType , class LabelType >
void StreamingNonAdaptiveLaplacian< Real , PixelChannels , LabelChannels , PixelType , LabelType >::Initialize( int major , int minor , int samples , const SphericalStencilTable< double >* stencilTable , const LabelType* unknownIndex , double vWeight , double gScale )
{
	this->major = major;
	this->minor = minor;
	this->samples = samples;
	_unknownIndex = unknownIndex;
	_vWeight = vWeight;
	_gScale  = gScale;
	_major = ( major+RealPerWord-1 ) / RealPerWord;
	_stencilTable = stencilTable;
	if( stencilTable )
	{
		phiWeight  = stencilTable->phiWeight  (                         major );
		stencilTable->setPhiStencil           ( phiStencil            , major );
		stencilTable->setDivergenceDPhiStencil( divergenceDPhiStencil , major );
	}
	else
	{
#if HIGH_PRECISION
		phiWeight = EquiRectangular< qfloat , qfloat >::PhiWeight( major );
		EquiRectangular< qfloat , qfloat >::PhiStencil  ( major ,   phiStencil );
		EquiRectangular< qfloat , qfloat >::DivergenceDPhiStencil( major , divergenceDPhiStencil );
#else // !HIGH_PRECISION
		phiWeight = EquiRectangular< double >::PhiWeight( major );
		EquiRectangular< double >::PhiStencil  ( major ,   phiStencil );
		EquiRectangular< double >::DivergenceDPhiStencil( major , divergenceDPhiStencil );
#endif // HIGH_PRECISION
	}
}
template< class Real , int PixelChannels , int LabelChannels , class PixelType , class LabelType >
void StreamingNonAdaptiveLaplacian< Real , PixelChannels , LabelChannels , PixelType , LabelType >::InitRestriction( void )
{
	int _channels = PixelChannels;
	if( pixels->HasAlpha() ) _channels++;
	FreeArray( localDMajor );
	FreeArray( localDMinor );
	FreeArray( localDMajorAccum );
	FreeArray( localDMinorAccum );

	index = 0;
	for( int c=0 ; c<PixelChannels ; c++) average[c] = 0;

	if( pixels->rows()   !=minor                             )     fprintf( stderr , "Pixel height failure: %d != %d\n" , pixels->rows()    , minor                             ) , exit(0);
	if( pixels->rowSize()!=major*_channels*sizeof(PixelType) )     fprintf( stderr , "Pixel  width failure: %d != %d\n" , pixels->rowSize() , major*_channels*sizeof(PixelType) ) , exit(0);
	if( labels )
	{
		if( labels->rows()!=minor )                                    fprintf( stderr , "Label height failure: %d != %d\n" , labels->rows()    , minor                                 ) , exit(0);
		if( labels->rowSize()!=major*LabelChannels*sizeof(LabelType) ) fprintf( stderr , "Label  width failure: %d != %d\n" , labels->rowSize() , major*LabelChannels*sizeof(LabelType) ) , exit(0);
	}
	_previousPixelRow              = AllocArray< PixelType >(      _channels*major ) , _pixelRow = AllocArray< PixelType >(      _channels*major );
	if( labels ) _previousLabelRow = AllocArray< LabelType >(  LabelChannels*major ) , _labelRow = AllocArray< LabelType >(  LabelChannels*major );
	if( _gScale )
	{
		_dx = AllocArray< Real >( PixelChannels*major );
		_dy = AllocArray< Real >( PixelChannels*major );
	}
	if( _vWeight ) _values = AllocArray< Real >( size*(_major+2)*PixelChannels*RealPerWord + 2*Degree*RealPerWord , ALIGNMENT );

	pixels->reset( true , 1 );
	if( labels ) labels->reset( true , 1 );
	_currentRow = -1;

	_readNextRow( );	// Read in the first row
	memcpy( _previousPixelRow , _pixelRow , sizeof( PixelType ) * _channels * major );
	if( labels ) memcpy( _previousLabelRow , _labelRow , sizeof( LabelType ) * LabelChannels * major );

	_readNextRow( );	// Read in the second row

	// Set the partials with respect to y
	if( _gScale )
		for( int x=0 ; x<major ; x++ )
		{
			bool useGradient = true;
			if( labels ) for( int c=0 ; c<LabelChannels ; c++ ) useGradient &= _labelRow[x+major*c]==_previousLabelRow[x+major*c];
			if( pixels->HasAlpha() && (_pixelRow[x+major*PixelChannels]==0 || _previousPixelRow[x+major*PixelChannels]==0) ) useGradient = true;
			if( _unknownIndex && labels )
			{
				bool isUnknown1=true , isUnknown2=true;
				for( int c=0 ; c<LabelChannels ; c++ ) 
				{
					if( _previousLabelRow[x + c*major]!=_unknownIndex[c] ) isUnknown1 = false;
					if(         _labelRow[x + c*major]!=_unknownIndex[c] ) isUnknown2 = false;
				}
				if( isUnknown1 || isUnknown2 ) useGradient = true;
			}
			if( useGradient ) for( int c=0 ; c<PixelChannels ; c++ ) _dy[x*PixelChannels+c] = Real( _pixelRow[x+major*c] )-Real( _previousPixelRow[x+major*c] );
			else              for( int c=0 ; c<PixelChannels ; c++ ) _dy[x*PixelChannels+c] = 0;
		}

	memcpy( _previousPixelRow , _pixelRow , sizeof(PixelType) * _channels * major );
	if( labels ) memcpy( _previousLabelRow , _labelRow , sizeof(LabelType) * LabelChannels * major );

	index = 0;

	if( _gScale )
	{
		localDMajor = AllocArray< Real >( dSize*(_major+2)*PixelChannels*RealPerWord + 2*Degree*RealPerWord , ALIGNMENT );
		localDMinor = AllocArray< Real >( dSize*(_major+2)*PixelChannels*RealPerWord + 2*Degree*RealPerWord , ALIGNMENT );
		memset( localDMajor , 0 , sizeof(Real) * (dSize*(_major+2)*PixelChannels*RealPerWord + 2*Degree*RealPerWord) );
		memset( localDMinor , 0 , sizeof(Real) * (dSize*(_major+2)*PixelChannels*RealPerWord + 2*Degree*RealPerWord) );
		localDMajorAccum = AllocArray< __m128 >( (_major+2) + 2*Degree , ALIGNMENT );
		localDMinorAccum = AllocArray< __m128 >( (_major+2) + 2*Degree , ALIGNMENT );
	}
	if( _vWeight ) localValueAccum = AllocArray< __m128 >( (_major+2) + 2*Degree , ALIGNMENT );

	int off1 , off2;
	off1 = FiniteElements1D< Real , DERIVATIVE(Type) , Degree-1 >::DotProduct< Type , Degree >::Helper::StopOffset();
	off2 = FiniteElements1D< Real ,            Type  , Degree   >::DotProduct< Type , Degree >::Helper::StopOffset();

	startRestriction  = (off1>off2)?off1:off2;
//	startRestriction -= Degree;
	if( parent ) parent->InitRestriction();
}
template< class Real , int PixelChannels , int LabelChannels , class PixelType , class LabelType >
void StreamingNonAdaptiveLaplacian< Real , PixelChannels , LabelChannels , PixelType , LabelType >::SetRestrictionIterations( int iters )
{
	if( parent ) parent->SetRestrictionIterations( iters );
}
template< class Real , int PixelChannels , int LabelChannels , class PixelType , class LabelType >
bool StreamingNonAdaptiveLaplacian< Real , PixelChannels , LabelChannels , PixelType , LabelType >::IterateRestriction( void )
{
	if( _showProgress )
		if     ( index <minor ) printf( "Restriction: [%.1f%%] [%d/%d]         \r" , float(index)/float(minor-1)*100 , index , minor );
		else if( index==minor ) printf( "                                              \r" );

	int _channels = PixelChannels;
	if( pixels->HasAlpha() ) _channels++;
	if( index>=minor ) return false;
	int previousRow = _currentRow-1;

	// Copy the finite differences from _dx , _dy into localDMajor , localDMinor
	if( _gScale )
	{
		if( previousRow>0 && previousRow<minor-1 )
			for( int c=0 ; c<PixelChannels ; c++ )
			{
				Real* ptrX = &localDMajor[(index%dSize)*(_major+2)*RealPerWord*PixelChannels + c*(_major+2)*RealPerWord + (Degree+1)*RealPerWord];
				for( int i=0 ; i<major ; i++ ) ptrX[i] = _dx[i*PixelChannels+c] , dXSquareNorm += ptrX[i]*ptrX[i];
				for( int i=0 ; i<RealPerWord ; i++ ) ptrX[-RealPerWord+i] = ptrX[major-RealPerWord+i];
				for( int i=0 ; i<RealPerWord ; i++ ) ptrX[major+i] = ptrX[i];
			}
		if( previousRow<minor-1 )
			for( int c=0 ; c<PixelChannels ; c++ )
			{
				Real* ptrY = &localDMinor[(index%dSize)*(_major+2)*RealPerWord*PixelChannels + c*(_major+2)*RealPerWord + (Degree+1)*RealPerWord];
				for( int i=0 ; i<major ; i++ ) ptrY[i] = _dy[i*PixelChannels+c] , dYSquareNorm += ptrY[i]*ptrY[i];
				for( int i=0 ; i<RealPerWord ; i++ ) ptrY[-RealPerWord+i] = ptrY[major-RealPerWord+i];
				for( int i=0 ; i<RealPerWord ; i++ ) ptrY[major+i] = ptrY[i];
			}
	}

	// Set the partials with respect to x
	if( _gScale )
		if( _currentRow>0 && _currentRow<minor-1 )
			for( int x=0 ; x<major ; x++ )
			{
				int x1 = x  , x2 = (x1+1)%major;
				bool useGradient = true;
				if( labels ) for( int c=0 ; c<LabelChannels ; c++ ) useGradient &= _previousLabelRow[x2+major*c]==_previousLabelRow[x1+major*c];
				if( pixels->HasAlpha() && (_previousPixelRow[x2+major*PixelChannels]==0 || _previousPixelRow[x1+major*PixelChannels]==0) ) useGradient = true;
				if( _unknownIndex && labels )
				{
					bool isUnknown1=true , isUnknown2=true;
					for( int c=0 ; c<LabelChannels ; c++ ) 
					{
						if( _previousLabelRow[x1 + c*major]!=_unknownIndex[c] ) isUnknown1 = false;
						if( _previousLabelRow[x2 + c*major]!=_unknownIndex[c] ) isUnknown2 = false;
					}
					if( isUnknown1 || isUnknown2 ) useGradient = true;
				}

				if( useGradient ) for( int c=0 ; c<PixelChannels ; c++ ) _dx[ x*PixelChannels+c ] = Real(_previousPixelRow[x2+major*c])-Real(_previousPixelRow[x1+major*c]);
				else              for( int c=0 ; c<PixelChannels ; c++ ) _dx[ x*PixelChannels+c ] = 0;
			}
		else memset( _dx , 0 , sizeof( Real ) * major * PixelChannels );

	// Set the partial with respect to y
	_readNextRow( );
	if( _currentRow<minor )
	{
		if( _gScale )
			for( int x=0 ; x<major ; x++ )
			{
				bool useGradient=true;
				if( labels ) for( int c=0 ; c<LabelChannels ; c++ ) useGradient &= _labelRow[x+major*c]==_previousLabelRow[x+major*c];
				if( pixels->HasAlpha() && (_pixelRow[x+major*PixelChannels]==0 || _previousPixelRow[x+major*PixelChannels]==0) ) useGradient = true;
				if( _unknownIndex && labels )
				{
					bool isUnknown1=true , isUnknown2=true;
					for( int c=0 ; c<LabelChannels ; c++ ) 
					{
						if( _previousLabelRow[x + c*major]!=_unknownIndex[c] ) isUnknown1 = false;
						if(         _labelRow[x + c*major]!=_unknownIndex[c] ) isUnknown2 = false;
					}
					if( isUnknown1 || isUnknown2 ) useGradient = true	;
				}
				if( useGradient ) for(int c=0 ; c<PixelChannels ; c++ ) _dy[x*PixelChannels+c] = Real(_pixelRow[x+major*c])-Real(_previousPixelRow[x+major*c]);
				else              for(int c=0 ; c<PixelChannels ; c++ ) _dy[x*PixelChannels+c] = 0;
			}

		memcpy( _previousPixelRow , _pixelRow , sizeof(PixelType)*_channels*major );
		if( labels ) memcpy( _previousLabelRow , _labelRow , sizeof(LabelType)*LabelChannels*major );
	}

	if( parent ) if( index>=startRestriction ) parent->IterateRestriction();
	index++;
	return true;
}

#include <Spheres/EquiRect.h>
template< class Real , int PixelChannels , int LabelChannels , class PixelType , class LabelType >
void StreamingNonAdaptiveLaplacian< Real , PixelChannels , LabelChannels , PixelType , LabelType >::SetRestriction( Pointer( Real )* lB , int idx , int major2 )
{
	double thetaStencil[5] , divergenceDThetaStencil[4] , laplacianThetaStencil[5];

	if( _stencilTable )
	{
		_stencilTable->setThetaStencil           (            thetaStencil , idx , minor );
		_stencilTable->setDivergenceDThetaStencil( divergenceDThetaStencil , idx , minor );
		_stencilTable->setLaplacianThetaStencil  (   laplacianThetaStencil , idx , minor , samples );
	}
	else
	{
		EquiRectangular< double >::LaplacianThetaStencil( idx , minor , laplacianThetaStencil , samples );
#if HIGH_PRECISION
		EquiRectangular< qfloat , qfloat >::ThetaStencil           ( idx , minor , thetaStencil );
		EquiRectangular< qfloat , qfloat >::DivergenceDThetaStencil( idx , minor , divergenceDThetaStencil );
#else // !HIGH_PRECISION
		EquiRectangular< double >::ThetaStencil           ( idx , minor , thetaStencil );
		EquiRectangular< double >::DivergenceDThetaStencil( idx , minor , divergenceDThetaStencil );
#endif // HIGH_PRECISION
	}
	
	// WARNING: This is not a safe test that the real type is either float or double...
#if USE_INTERIOR
	if( idx>Degree && idx<minor-Degree-1 && ( sizeof( Real )==sizeof( float ) || sizeof( Real )==sizeof( double ) ) )
	{
		SetInteriorRestriction( lB , idx , major2 );
		return;
	}
#endif // USE_INTERIOR
	int dyStart = idx-2 , yStart = idx-2;
	Pointer( Real ) localDX[2*Degree+1];
	Pointer( Real ) localDY[2*Degree  ];
	Pointer( Real ) localValues[2*Degree+1];
	for( int c=0 ; c<PixelChannels ; c++ )
	{
		Pointer( Real ) myLB = lB[c];

		if( _gScale )
		{
			for( int dy=0 ; dy<=2*Degree ; dy++ ) localDX[dy] = localDMajor + ( ( yStart+dy+dSize)%dSize ) * (_major+2)*RealPerWord*PixelChannels + c*(_major+2)*RealPerWord + (Degree+1)*RealPerWord;
			for( int dy=0 ; dy< 2*Degree ; dy++ ) localDY[dy] = localDMinor + ( (dyStart+dy+dSize)%dSize ) * (_major+2)*RealPerWord*PixelChannels + c*(_major+2)*RealPerWord + (Degree+1)*RealPerWord;

			if( idx==0 || idx==minor-1 )
			{
				Real temp = 0;
				for( int dy=0 ; dy<2*Degree ; dy++ )
				{
					if( dyStart+dy<0 || dyStart+dy>=minor-1 ) continue;
					for( int x=0 ; x<major ; x++ ) temp -= localDY[dy][x] * phiWeight * divergenceDThetaStencil[dy];
				}
				myLB[0] = temp * _gScale;
				outputSquareNorm += temp*temp;
			}
			else
			{
				for( int x=0 ; x<major ; x++ )
				{
					int dxStart = x-2 , xStart = x-2;
					Real temp = 0;

					// Partial w.r.t major index
					for( int dy=0 ; dy<=2*Degree ; dy++ )
					{
						if( yStart+dy<=0 || yStart+dy>=minor-1 ) continue;
						for( int dx=0 ; dx<2*Degree ; dx++ ) temp -= localDX[dy][dxStart+dx] * divergenceDPhiStencil[dx] * laplacianThetaStencil[dy];
					}
					// Partial w.r.t minor index
					for( int dy=0 ; dy<2*Degree ; dy++ )
					{
						if( dyStart+dy<0 || dyStart+dy>=minor-1 ) continue;
						for( int dx=0 ; dx<=2*Degree ; dx++ ) temp -= localDY[dy][xStart+dx] * phiStencil[dx] * divergenceDThetaStencil[dy];
					}
					myLB[x] = temp * _gScale;
					outputSquareNorm += temp*temp;
				}
			}
		}
		else memset( myLB , 0 , sizeof( Real ) * major );
		if( _vWeight )
		{
			for( int dy=0 ; dy<=2*Degree ; dy++ ) localValues[dy] = _values + ( ( yStart+dy+size)%size ) * (_major+2)*RealPerWord*PixelChannels + c*(_major+2)*RealPerWord + (Degree+1)*RealPerWord;
			if( idx==0 || idx==minor-1 )
			{
				Real temp = 0;
				for( int dy=0 ; dy<=2*Degree ; dy++ )
				{
					if( yStart+dy<0 || yStart+dy>minor-1 ) continue;
					else if( yStart+dy==0 || yStart+dy==minor-1 ) temp += localValues[dy][0] * thetaStencil[dy] * 2. * M_PI;
					else for( int x=0 ; x<major ; x++ )           temp += localValues[dy][x] * phiWeight * thetaStencil[dy];
				}
				myLB[0] += temp * _vWeight;
			}
			else
			{
				for( int x=0 ; x<major ; x++ )
				{
					int xStart = x-2;
					Real temp = 0;

					for( int dy=0 ; dy<=2*Degree ; dy++ )
					{
						if( yStart+dy<0 || yStart+dy>minor-1 ) continue;
						else if ( yStart+dy==0 || yStart+dy==minor-1 ) temp += localValues[dy][        0] * thetaStencil[dy] * phiWeight;
						else for( int dx=0 ; dx<=2*Degree ; dx++ )     temp += localValues[dy][xStart+dx] * thetaStencil[dy] * phiStencil[dx];
					}
					myLB[x] += temp * _vWeight;
				}
			}
		}
	}
}
template< class Real , int PixelChannels , int LabelChannels , class PixelType , class LabelType >
void StreamingNonAdaptiveLaplacian< Real , PixelChannels , LabelChannels , PixelType , LabelType >::SetInteriorRestriction( Pointer( float )* lB , int idx , int major2 )
{
	double thetaStencil[5] , divergenceDThetaStencil[4] , laplacianThetaStencil[5];

	if( _stencilTable )
	{
		_stencilTable->setThetaStencil           (            thetaStencil , idx , minor );
		_stencilTable->setDivergenceDThetaStencil( divergenceDThetaStencil , idx , minor );
		_stencilTable->setLaplacianThetaStencil  (   laplacianThetaStencil , idx , minor , samples );
	}
	else
	{
		EquiRectangular< double >::LaplacianThetaStencil( idx , minor , laplacianThetaStencil , samples );
#if HIGH_PRECISION
		EquiRectangular< qfloat , qfloat >::ThetaStencil           ( idx , minor , thetaStencil );
		EquiRectangular< qfloat , qfloat >::DivergenceDThetaStencil( idx , minor , divergenceDThetaStencil );
#else // !HIGH_PRECISION
		EquiRectangular< double >::ThetaStencil           ( idx , minor , thetaStencil );
		EquiRectangular< double >::DivergenceDThetaStencil( idx , minor , divergenceDThetaStencil );
#endif // HIGH_PRECISION
	}

	int dyStart = idx-2 , yStart = idx-2;

	Pointer( __m128 ) localValues[2*Degree+1];
	Pointer( __m128 ) localDX[2*Degree+1];
	Pointer( __m128 ) localDY[2*Degree  ];

	Pointer( __m128 ) localVAccum = localValueAccum+1;
	Pointer( __m128 ) localDXAccum = localDMajorAccum+1;
	Pointer( __m128 ) localDYAccum = localDMinorAccum+1;

	for( int c=0 ; c<PixelChannels ; c++ )
	{
		Pointer( float ) myLB = lB[c];
		if( _gScale )
		{
			for( int dy=0 ; dy<=2*Degree ; dy++ ) localDX[dy] = ( Pointer( __m128 ) )( localDMajor + ( ( yStart+dy+dSize)%dSize ) * (_major+2)*RealPerWord*PixelChannels + c*(_major+2)*RealPerWord + (Degree+1)*RealPerWord );
			for( int dy=0 ; dy< 2*Degree ; dy++ ) localDY[dy] = ( Pointer( __m128 ) )( localDMinor + ( (dyStart+dy+dSize)%dSize ) * (_major+2)*RealPerWord*PixelChannels + c*(_major+2)*RealPerWord + (Degree+1)*RealPerWord );

			__declspec ( align( ALIGNMENT ) ) float scratch[4];
			__m128 dot[2*Degree+1] , divergenceDDot[2*Degree];
			for( int d=0 ; d<=2*Degree ; d++ )
			{
				for( int j=0 ; j<4 ; j++ ) scratch[j] = laplacianThetaStencil[d];
				dot[d] = _mm_load_ps( scratch );
			}
			for( int d=0 ; d< 2*Degree ; d++ )
			{
				for( int j=0 ; j<4 ; j++ ) scratch[j] = divergenceDThetaStencil[d];
				divergenceDDot[d] = _mm_load_ps( scratch );
			}
			for( int x=-1 ; x<_major+1 ; x++ ) localDXAccum[x] = _mm_mul_ps( dot[0] , localDX[0][x] ) , localDYAccum[x] = _mm_mul_ps( divergenceDDot[0] , localDY[0][x] );
			for( int d=1 ; d<=2*Degree ; d++ )
			{
				for( int x=-1 ; x<_major+1 ; x++ ) localDXAccum[x] = _mm_add_ps( localDXAccum[x] , _mm_mul_ps( dot[d] , localDX[d][x] ) );
				if( d<2*Degree ) for( int x=-1 ; x<_major+1 ; x++ ) localDYAccum[x] = _mm_add_ps( localDYAccum[x] , _mm_mul_ps( divergenceDDot[d] , localDY[d][x] ) );
			}

			for( int x=0 ; x<major ; x++ )
			{
				Real temp = 0;
				int dxStart = x-2 , xStart = x-2;
				{
					const Pointer( float ) localDX = ( ( Pointer( float ) )(localDXAccum) ) + dxStart;
					const Pointer( float ) localDY = ( ( Pointer( float ) )(localDYAccum) ) +  xStart;
					temp =
						- ( localDX[0] - localDX[3] ) * divergenceDPhiStencil[0] - (localDX[1] - localDX[2] ) * divergenceDPhiStencil[1]
						- localDY[2] * phiStencil[2] - ( localDY[1] + localDY[3] ) * phiStencil[1] - ( localDY[0] + localDY[4] ) * phiStencil[0];
				}
				myLB[x] = temp * _gScale;
				outputSquareNorm += temp*temp;
			}
		}
		else memset( myLB , 0 , sizeof( float ) * major );

		if( _vWeight )
		{
			for( int dy=0 ; dy<=2*Degree ; dy++ ) localValues[dy] = ( Pointer( __m128 ) )( _values + ( ( yStart+dy+size)%size ) * (_major+2)*RealPerWord*PixelChannels + c*(_major+2)*RealPerWord + (Degree+1)*RealPerWord );

			__declspec ( align( ALIGNMENT ) ) float scratch[4];
			__m128 dot[2*Degree+1];
			for( int d=0 ; d<=2*Degree ; d++ )
			{
				for( int j=0 ; j<4 ; j++ ) scratch[j] = thetaStencil[d];
				dot[d] = _mm_load_ps( scratch );
			}
			for( int x=-1 ; x<_major+1 ; x++ ) localVAccum[x] = _mm_mul_ps( dot[0] , localValues[0][x] );
			for( int d=1 ; d<=2*Degree ; d++ ) for( int x=-1 ; x<_major+1 ; x++ ) localVAccum[x] = _mm_add_ps( localVAccum[x] , _mm_mul_ps( dot[d] , localValues[d][x] ) );

			for( int x=0 ; x<major ; x++ )
			{
				const Pointer( float ) localV = ( ( Pointer( float ) )(localVAccum) ) + x - 2;
				Real temp = localV[2] * phiStencil[2] + ( localV[1] + localV[3] ) * phiStencil[1]  + ( localV[0] + localV[4] ) * phiStencil[0];
				myLB[x] += temp * _vWeight;
			}
		}
	}
}
template< class Real , int PixelChannels , int LabelChannels , class PixelType , class LabelType >
void StreamingNonAdaptiveLaplacian< Real , PixelChannels , LabelChannels , PixelType , LabelType >::SetInteriorRestriction( Pointer( double )* lB , int idx , int major2 )
{
	double thetaStencil[5] , divergenceDThetaStencil[4] , laplacianThetaStencil[5];

	if( _stencilTable )
	{
		_stencilTable->setThetaStencil           (            thetaStencil , idx , minor );
		_stencilTable->setDivergenceDThetaStencil( divergenceDThetaStencil , idx , minor );
		_stencilTable->setLaplacianThetaStencil  (   laplacianThetaStencil , idx , minor , samples );
	}
	else
	{
		EquiRectangular< double >::LaplacianThetaStencil( idx , minor , laplacianThetaStencil , samples );
#if HIGH_PRECISION
		EquiRectangular< qfloat , qfloat >::ThetaStencil           ( idx , minor , thetaStencil );
		EquiRectangular< qfloat , qfloat >::DivergenceDThetaStencil( idx , minor , divergenceDThetaStencil );
#else // !HIGH_PRECISION
		EquiRectangular< double >::ThetaStencil           ( idx , minor , thetaStencil );
		EquiRectangular< double >::DivergenceDThetaStencil( idx , minor , divergenceDThetaStencil );
#endif // HIGH_PRECISION
	}
	int dyStart = idx-2 , yStart = idx-2;

	Pointer( __m128d ) localValues[2*Degree+1];
	Pointer( __m128d ) localDX[2*Degree+1];
	Pointer( __m128d ) localDY[2*Degree  ];
	Pointer( __m128d ) localVAccum;
	Pointer( __m128d ) localDXAccum;
	Pointer( __m128d ) localDYAccum;

	if( _gScale )
	{
		localDXAccum = ( Pointer( __m128d ) )( localDMajorAccum+1 );
		localDYAccum = ( Pointer( __m128d ) )( localDMinorAccum+1 );
	}
	if( _vWeight ) localVAccum  = ( Pointer( __m128d ) )( localValueAccum+1 );

	for( int c=0 ; c<PixelChannels ; c++ )
	{
		Pointer( double ) myLB = lB[c];
		if( _gScale )
		{
			for( int dy=0 ; dy<=2*Degree ; dy++ ) localDX[dy] = ( Pointer( __m128d ) )( localDMajor + ( ( yStart+dy+dSize)%dSize ) * (_major+2)*RealPerWord*PixelChannels + c*(_major+2)*RealPerWord + (Degree+1)*RealPerWord );
			for( int dy=0 ; dy< 2*Degree ; dy++ ) localDY[dy] = ( Pointer( __m128d ) )( localDMinor + ( (dyStart+dy+dSize)%dSize ) * (_major+2)*RealPerWord*PixelChannels + c*(_major+2)*RealPerWord + (Degree+1)*RealPerWord );

			__declspec ( align( ALIGNMENT ) ) double scratch[RealPerWord];
			__m128d dot[2*Degree+1] , divergenceDDot[2*Degree];
			for( int d=0 ; d<=2*Degree ; d++ )
			{
				for( int j=0 ; j<RealPerWord ; j++ ) scratch[j] = laplacianThetaStencil[d];
				dot[d] = _mm_load_pd( scratch );
			}
			for( int d=0 ; d< 2*Degree ; d++ )
			{
				for( int j=0 ; j<RealPerWord ; j++ ) scratch[j] = divergenceDThetaStencil[d];
				divergenceDDot[d] = _mm_load_pd( scratch );
			}
			for( int x=-1 ; x<_major+1 ; x++ ) localDXAccum[x] = _mm_mul_pd( dot[0] , localDX[0][x] ) , localDYAccum[x] = _mm_mul_pd( divergenceDDot[0] , localDY[0][x] );
			for( int d=1 ; d<=2*Degree ; d++ )
			{
				for( int x=-1 ; x<_major+1 ; x++ ) localDXAccum[x] = _mm_add_pd( localDXAccum[x] , _mm_mul_pd( dot[d] , localDX[d][x] ) );
				if( d<2*Degree ) for( int x=-1 ; x<_major+1 ; x++ ) localDYAccum[x] = _mm_add_pd( localDYAccum[x] , _mm_mul_pd( divergenceDDot[d] , localDY[d][x] ) );
			}

			for( int x=0 ; x<major ; x++ )
			{
				Real temp = 0;
				int dxStart = x-2 , xStart = x-2;
				{
					const Pointer( double ) localDX = ( ( Pointer( double ) )(localDXAccum) ) + dxStart;
					const Pointer( double ) localDY = ( ( Pointer( double ) )(localDYAccum) ) +  xStart;
					temp =
						- ( localDX[0] - localDX[3] ) * divergenceDPhiStencil[0] - (localDX[1] - localDX[2] ) * divergenceDPhiStencil[1]
						- localDY[2] * phiStencil[2] - ( localDY[1] + localDY[3] ) * phiStencil[1] - ( localDY[0] + localDY[4] ) * phiStencil[0];
				}
				myLB[x] = temp * _gScale;
				outputSquareNorm += temp*temp;
			}
		}
		else memset( myLB , 0 , sizeof( double ) * major );

		if( _vWeight )
		{
			for( int dy=0 ; dy<=2*Degree ; dy++ ) localValues[dy] = ( Pointer( __m128d ) )( _values + ( ( yStart+dy+size)%size ) * (_major+2)*RealPerWord*PixelChannels + c*(_major+2)*RealPerWord + (Degree+1)*RealPerWord );

			__declspec ( align( ALIGNMENT ) ) double scratch[RealPerWord];
			__m128d dot[2*Degree+1];
			for( int d=0 ; d<=2*Degree ; d++ )
			{
				for( int j=0 ; j<RealPerWord ; j++ ) scratch[j] = thetaStencil[d];
				dot[d] = _mm_load_pd( scratch );
			}
			for( int x=-1 ; x<_major+1 ; x++ ) localVAccum[x] = _mm_mul_pd( dot[0] , localValues[0][x] );
			for( int d=1 ; d<=2*Degree ; d++ ) for( int x=-1 ; x<_major+1 ; x++ ) localVAccum[x] = _mm_add_pd( localVAccum[x] , _mm_mul_pd( dot[d] , localValues[d][x] ) );

			for( int x=0 ; x<major ; x++ )
			{
				const Pointer( double ) localV = ( ( Pointer( double ) )(localVAccum) ) + x - 2;
				Real temp = localV[2] * phiStencil[2] + ( localV[1] + localV[3] ) * phiStencil[1]  + ( localV[0] + localV[4] ) * phiStencil[0];
				myLB[x] += temp * _vWeight;
			}
		}
	}
}

template< class Real , int PixelChannels , int LabelChannels , class PixelType , class LabelType >
void StreamingNonAdaptiveLaplacian< Real , PixelChannels , LabelChannels , PixelType , LabelType >::SolveRestriction( void )
{
	// Run to completion...
	while( IterateRestriction( ) ){;}
	// ...and finish up the trailing parent
	if( parent ) parent->SolveRestriction();
}
