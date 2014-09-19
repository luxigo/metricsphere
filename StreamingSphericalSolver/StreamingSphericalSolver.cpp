#include <stdio.h>
#include <stdlib.h>
#include <Util/CmdLineParser.h>
#include <Util/Time.h>
#include <Spheres/AdaptiveEquiRect.h>
#include <Util/ImageStream.h>
#include <LaplacianMatrix/SphericalStreamingSolver.h>

cmdLineInt Iters( "iters" , 5 ) , MinRes( "minRes" , 64 ) , VCycles( "vCycles" , 1 );
cmdLineString In( "in" ) , Labels( "labels" ) , Out( "out" ) , StencilIO( "stencilIO" ) , TempDir( "tempDir" );
cmdLineFloat IWeight( "iWeight" , 0 ) , GScale( "gScale" , 1 ) , GWeight( "gWeight" , 1 );
cmdLineInt Quality( "quality" , 100 ) , InCoreRes( "inCoreRes" , 1024 );
cmdLineIntArray< 3 > UnknownIndex( "unknownIndex" );
cmdLineReadable Progress( "progress" ) , HighPrecision( "highPrecision" ) , Verbose( "verbose" ) , HDRLabels( "hdrLabels" );

cmdLineReadable* params[]=
{
	&Iters , &MinRes , &VCycles , 
	&In , &Labels , &Out , &StencilIO , &TempDir ,
	&Verbose ,
	&IWeight , &GScale , &GWeight , &InCoreRes , &Quality ,
	&HDRLabels , &UnknownIndex ,
	&Progress , &HighPrecision ,
};
void ShowUsage(char* ex)
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input spherical color composite>\n" , In.name );
	printf( "\t[--%s <input labels]\n" , Labels.name );
	printf( "\t --%s <output stitched spherical images>\n" , Out.name );

	printf( "\t[--%s <stencil input/output file>]\n" , StencilIO.name );

	printf( "\t[--%s <number of v-cycles>=%d]\n" , VCycles.name , VCycles.value );
	printf( "\t[--%s <gauss-seidel iterations>=%d]\n" , Iters.name , Iters.value );
	printf( "\t[--%s]\n" , HighPrecision.name );

	printf( "\t[--%s <minimum multigrid resolution>=%d]\n" , MinRes.name , MinRes.value );
	printf( "\t[--%s <minimum in-core resolution>=%d]\n" , InCoreRes.name , InCoreRes.value );

	printf( "\t[--%s <compression quality>=%d]\n" , Quality.name , Quality.value );

	printf( "\t[--%s <interpolation weight>=%f]\n" , IWeight.name , IWeight.value );
	printf( "\t[--%s <gradient scale>=%f]\n" , GScale.name , GScale.value );
	printf( "\t[--%s <gradient weight>=%f]\n" , GWeight.name , GWeight.value );
	printf( "\t[--%s <label of pixels whose value should be forced to black>]\n" , UnknownIndex.name );
	printf( "\t[--%s <directory for temporary files>]\n" , TempDir.name );
	printf( "\t[--%s]\n" , Verbose.name );
	printf( "\t[--%s]\n" , Progress.name );
	printf( "\t[--%s]\n" , HDRLabels.name );
}

template< class Real >
inline Vector< Real > NoAverageMultiply( const SparseMatrix< Real >& A , const Vector< Real >& x )
{
	Vector< Real > y = A*x;
	double average = 0;
	for( int i=0 ; i<x.Dimensions() ; i++ )	average += x[i];
	average /= x.Dimensions();
	for( int i=0 ; i<x.Dimensions() ; i++ )	y[i] += average;
	return y;
}
template< class Real >
inline Vector< Real > YesAverageMultiply( const SparseMatrix< Real >& A , const Vector< Real >& x )
{
	return  A*x;
}
template< class Real >
int SolveConjugateGradient( const SparseMatrix< Real >& A , const Vector< Real >& b,const int& iters , Vector< Real >& x , bool noAverage=false , bool negativeDefinite=false )
{
	double eps=1e-22;
	Vector< Real > r , d;
	if( negativeDefinite )
		if( noAverage ) r =  NoAverageMultiply( A , x ) - b;
		else            r = YesAverageMultiply( A , x ) - b;
	else
		if( noAverage ) r = b -  NoAverageMultiply( A , x );
		else            r = b - YesAverageMultiply( A , x );
	d = r;

	double delta_new = r.Dot(r);
	double delta_0 = delta_new;
	int i;
	for( i=0 ; i<iters && delta_new>eps*delta_0 ; i++ )
	{
		Vector< Real > q;
		if( noAverage ) q =  NoAverageMultiply( A , d );
		else            q = YesAverageMultiply( A , d );
		if( negativeDefinite ) q *= -1;
		double alpha = delta_new / d.Dot( q );
		x = x + d*alpha;

		if( !(i%50) )
			if( negativeDefinite )
				if( noAverage ) r =  NoAverageMultiply( A , x ) - b;
				else            r = YesAverageMultiply( A , x ) - b;
			else
				if( noAverage ) r = b -  NoAverageMultiply( A , x );
				else            r = b - YesAverageMultiply( A , x );
		else r = r - q*alpha;

		double delta_old = delta_new;
		delta_new = r.Dot( r );
		double beta = delta_new / delta_old;
		d = r + d*beta;
	}
	return i;
}
template< class Real , class IOReal , class LabelType , int Channels >
void StitchImage( int height , int samples , const SphericalStencilTable< double >* stencilTable , StreamingGrid* in , StreamingGrid* labels , StreamingGrid* out , int iters , double iWeight , double gScale , double gWeight , const LabelType* unknownIndex )
{
	int width = height*2;
	StreamingNonAdaptiveLaplacian< Real , Channels , 3 , Real , LabelType > streamingLaplacian;
	StreamingAdaptiveNonAdaptiveConverter< Real , Channels > streamingConverter;

	streamingLaplacian.pixels = in;
	streamingLaplacian.labels = labels;
	streamingLaplacian.Initialize( width , height , samples , stencilTable , unknownIndex , iWeight , gScale*gWeight );
	streamingConverter.Initialize( height );

	const int Dim = Channels;
	int idx , fullDepth = 0 , outOfCoreDepth = 0 , inCoreDepth;
	while( 	!( (width>>fullDepth)&3 ) && !( (height>>fullDepth)&1 ) && (width>>fullDepth)>=MinRes.value && (height>>fullDepth)>=MinRes.value ) fullDepth++;
	int lowHeight = height>>fullDepth;
	int lowWidth = lowHeight * 2;
	fullDepth++;
	while( 	!( (width>>outOfCoreDepth)&3 ) && !( (height>>outOfCoreDepth)&1 ) && (width>>outOfCoreDepth)>=InCoreRes.value && (height>>outOfCoreDepth)>=InCoreRes.value ) outOfCoreDepth++;
	int midHeight = height>>outOfCoreDepth;
	int midWidth = midHeight * 2;
	outOfCoreDepth++;
	// WARNING: We are duplicating the solve on the liminal row
	inCoreDepth = fullDepth - outOfCoreDepth + 1;

	double t = Time();
	MultiStreamIOClient *fullXStream = NULL , *fullBStream = NULL;
	Pointer( Real ) midX = AllocArray< Real >( midWidth * midHeight * Dim );
	Pointer( Real ) midB = AllocArray< Real >( midWidth * midHeight * Dim );
	Pointer( Real ) lowB = AllocArray< Real >( lowWidth * lowHeight * Dim );
	Pointer( Real ) lowX = AllocArray< Real >( lowWidth * lowHeight * Dim );

	memset( lowX , 0 , sizeof( Real ) * lowWidth * lowHeight * Dim );
	memset( lowB , 0 , sizeof( Real ) * lowWidth * lowHeight * Dim );
	if( VCycles.value>1 )
	{
		fullXStream = new MultiStreamIOClient( width * Dim * sizeof( Real ) , height , STREAMING_GRID_BUFFER_MULTIPLIER , NULL , true );
		fullBStream = new MultiStreamIOClient( width * Dim * sizeof( Real ) , height , STREAMING_GRID_BUFFER_MULTIPLIER , NULL , true );
	}

	MemoryBackedGrid midXStream( ( Pointer( byte ) )midX , midWidth * sizeof( Real ) * Dim , midHeight , true );
	MemoryBackedGrid midBStream( ( Pointer( byte ) )midB , midWidth * sizeof( Real ) * Dim , midHeight , true );
	MemoryBackedGrid lowBStream( ( Pointer( byte ) )lowB , lowWidth * sizeof( Real ) * Dim , lowHeight , true );
	MemoryBackedGrid lowXStream( ( Pointer( byte ) )lowX , lowWidth * sizeof( Real ) * Dim , lowHeight , true );

	if( Verbose.set ) printf( "Set up time: %f\n" , Time() - t );

	MultiGridSphericalStreamingSolver< Real , IOReal > *inCoreSolvers , *outOfCoreSolvers;
	inCoreSolvers    = new MultiGridSphericalStreamingSolver< Real , IOReal >[    inCoreDepth ];
	outOfCoreSolvers = new MultiGridSphericalStreamingSolver< Real , IOReal >[ outOfCoreDepth ];

	for( int d=1 ; d<inCoreDepth    ; d++ ) inCoreSolvers   [d].child = &inCoreSolvers   [d-1] , inCoreSolvers   [d-1].parent = &inCoreSolvers   [d];
	for( int d=1 ; d<outOfCoreDepth ; d++ ) outOfCoreSolvers[d].child = &outOfCoreSolvers[d-1] , outOfCoreSolvers[d-1].parent = &outOfCoreSolvers[d];

	streamingLaplacian.parent = &streamingConverter;
	streamingConverter.child = &streamingLaplacian;
	streamingConverter.parent = &outOfCoreSolvers[0];

	outOfCoreSolvers[0].Initialize( Dim ,    width ,    height , samples , stencilTable , iWeight , gWeight , true );
	inCoreSolvers   [0].Initialize( Dim , midWidth , midHeight , samples , stencilTable , iWeight , gWeight , false );

	for( int v=0 ; v<VCycles.value ; v++ )
	{
		double vTime = Time();
		/////////////////
		// Restriction //
		/////////////////

		/////////////////
		// Out-of-Core //
		/////////////////
		t = Time();
		for( int d=0 ; d<outOfCoreDepth ; d++ )
		{
			for( int c=0 ; c<Dim ; c++ ) outOfCoreSolvers[d].bSquareNorm[c] = outOfCoreSolvers[d].rSquareNorm[c] = outOfCoreSolvers[d].xSquareNorm[c] = 0;
			outOfCoreSolvers[d].setResidual = true;
			outOfCoreSolvers[d].inX = outOfCoreSolvers[d].inB = outOfCoreSolvers[d].outX = outOfCoreSolvers[d].outB = NULL;
		}
		if( v )
		{
			outOfCoreSolvers[0].inX = fullXStream;
			outOfCoreSolvers[0].inB = fullBStream;
			outOfCoreSolvers[0].child = NULL;
		}
		else outOfCoreSolvers[0].child = &streamingConverter;

		outOfCoreSolvers[outOfCoreDepth-1].outB = &midBStream;
		outOfCoreSolvers[outOfCoreDepth-1].outX = &midXStream;

		if( !v )
		{
			( ( MultiGridRestrictionNode< Real >* )&streamingLaplacian )->InitRestriction( Progress.set );
			streamingLaplacian.SetRestrictionIterations( Iters.value );
			streamingLaplacian.SolveRestriction( );
		}
		else
		{
			( ( MultiGridRestrictionNode< Real >* ) outOfCoreSolvers )->InitRestriction( Progress.set );
			outOfCoreSolvers[0].SetRestrictionIterations( Iters.value );
			outOfCoreSolvers[0].SolveRestriction( );
		}

		if( Verbose.set )
		{
			printf( "Out-of-Core Restriction Time: %f\n" , Time() - t );
			for( int d=0 ; d<outOfCoreDepth ; d++ )
			{
				double bSquareNorm = 0 , rSquareNorm = 0;
				for( int c=0 ; c<Dim ; c++ ) bSquareNorm += outOfCoreSolvers[d].bSquareNorm[c] , rSquareNorm += outOfCoreSolvers[d].rSquareNorm[c];
				printf( "\tError[%d x %d]:\t%13.10f -> %13.10f\t%g\n" ,	outOfCoreSolvers[d].columns , outOfCoreSolvers[d].rows , sqrt( bSquareNorm ) , sqrt( rSquareNorm ) , sqrt( rSquareNorm/bSquareNorm ) );
			}
		}


		/////////////
		// In-Core //
		/////////////
		t = Time();
		for( int d=0 ; d<inCoreDepth ; d++ )
		{
			for( int c=0 ; c<Dim ; c++ ) inCoreSolvers[d].bSquareNorm[c] = inCoreSolvers[d].rSquareNorm[c] = inCoreSolvers[d].xSquareNorm[c] = 0;
			inCoreSolvers[d].setResidual = true;
			inCoreSolvers[d].inX = inCoreSolvers[d].inB = inCoreSolvers[d].outX = inCoreSolvers[d].outB = NULL;
		}

		inCoreSolvers[0].inX = &midXStream;
		inCoreSolvers[0].inB = &midBStream;

		inCoreSolvers[inCoreDepth-1].outB = &lowBStream;
		inCoreSolvers[inCoreDepth-1].outX = &lowXStream;

		inCoreSolvers[0].InitRestriction( );
		inCoreSolvers[0].SetRestrictionIterations( Iters.value );
		inCoreSolvers[0].SolveRestriction( );

		if( Verbose.set )
		{
			printf( "In-Core Restriction Time: %f\n" , Time() - t );
			for( int d=0 ; d<inCoreDepth ; d++ )
			{
				double bSquareNorm = 0 , rSquareNorm = 0;
				for( int c=0 ; c<Dim ; c++ ) bSquareNorm += inCoreSolvers[d].bSquareNorm[c] , rSquareNorm += inCoreSolvers[d].rSquareNorm[c];
				printf( "\tError[%d x %d]:\t%13.10f -> %13.10f\t%g\n" ,	inCoreSolvers[d].columns , inCoreSolvers[d].rows , sqrt( bSquareNorm ) , sqrt( rSquareNorm ) , sqrt( rSquareNorm/bSquareNorm ) );
			}
		}

		///////////////
		// CG-Solver //
		///////////////
		{
			t = Time();
			AdaptiveEquiRectangular< double > lowSphere( lowHeight , samples , stencilTable );
			SparseMatrix< double > laplacian;
			Vector< double > xx[Dim] , bb[Dim];
			for( int c=0 ; c<Dim ; c++ )
			{
				xx[c].Resize( lowSphere.dimension() );
				bb[c].Resize( lowSphere.dimension() );
			}
			// Get the constraint and solution vectors
			idx = 0;
			for( int j=0 ; j<lowHeight ; j++ )
			{
				for( int i=0 ; i<lowSphere.dimension( j ) ; i++ )
				{
					for( int c=0 ; c<Dim ; c++ )
					{
						bb[c][ idx ] = lowB[ lowWidth*Dim*j + lowWidth*c + i ];
						xx[c][ idx ] = lowX[ lowWidth*Dim*j + lowWidth*c + i ];
					}
					idx++;
				}
			}

			// Get the Laplacian matrix
			lowSphere.laplacianMatrix( laplacian , iWeight , gWeight );
			double bSquareNorm = 0 , rSquareNorm = 0;
			for( int c=0 ; c<Dim ; c++ )
			{
				bSquareNorm += ( bb[c]-laplacian*xx[c] ).SquareNorm();
				SolveConjugateGradient( laplacian , bb[c] , 6*int( sqrt( laplacian.groups + 1.0 ) ) , xx[c] , iWeight==0 );
				rSquareNorm += ( bb[c]-laplacian*xx[c] ).SquareNorm();
			}
			if( Verbose.set )
			{
				printf( "\tError[%d x %d]:\t%13.10f -> %13.10f\t%g\n" ,	lowWidth , lowHeight , sqrt( bSquareNorm ) , sqrt( rSquareNorm ) , sqrt( rSquareNorm/bSquareNorm ) );
				printf( "Conjugate-Gradient Time: %f\n" , Time() - t );
			}

			// Set the solution vectors
			idx = 0;
			for( int j=0 ; j<lowHeight ; j++ )
			{
				for( int i=0 ; i<lowSphere.dimension(j) ; i++ )
				{
					for( int c=0 ; c<Dim ; c++ ) lowX[ lowWidth*Dim*j + lowWidth*c + i ] = xx[c][ idx ];
					idx++;
				}
			}
		}

		//////////////////
		// Prolongation //
		//////////////////

		/////////////
		// In-Core //
		/////////////

		t = Time();
		for( int d=0 ; d<inCoreDepth ; d++ )
		{
			for( int c=0 ; c<Dim ; c++ ) inCoreSolvers[d].bSquareNorm[c] = inCoreSolvers[d].rSquareNorm[c] = inCoreSolvers[d].xSquareNorm[c] = 0;
			inCoreSolvers[d].inX = inCoreSolvers[d].inB = inCoreSolvers[d].outX = inCoreSolvers[d].outB = NULL;
		}
		inCoreSolvers[inCoreDepth-1].inX = &lowXStream;
		inCoreSolvers[0].outX = &midXStream;
		inCoreSolvers[0].outB = &midBStream;

		inCoreSolvers[0].InitProlongation( );
		inCoreSolvers[inCoreDepth-1].SetProlongationIterations( Iters.value );
		inCoreSolvers[inCoreDepth-1].SolveProlongation( );

		if( Verbose.set )
		{
			for( int d=inCoreDepth-1 ; d>=0 ; d-- )
			{
				double bSquareNorm = 0 , rSquareNorm = 0;
				for( int c=0 ; c<Dim ; c++ ) bSquareNorm += inCoreSolvers[d].bSquareNorm[c] , rSquareNorm += inCoreSolvers[d].rSquareNorm[c];
				printf( "\tError[%d x %d]:\t%13.10f -> %13.10f\t%g\n" ,	inCoreSolvers[d].columns , inCoreSolvers[d].rows , sqrt( bSquareNorm ) , sqrt( rSquareNorm ) , sqrt( rSquareNorm/bSquareNorm ) );
			}
			printf( "In-Core Prolongation Time: %f\n" , Time() - t );
		}

		if( v==VCycles.value-1 && iWeight==0 )
		{
			AdaptiveEquiRectangular< double > midSphere( midHeight , samples , stencilTable );
			Vector< double > weights;
			midSphere.elementWeights( weights );
			for( int c=0 ; c<Dim ; c++ )
			{
				double average = 0;
				idx = 0;
				for( int y=0 ; y<midHeight ; y++ ) for( int x=0 ; x<midSphere.dimension( y ) ; x++ ) average += midX[ y*midWidth*Dim + c*midWidth + x ] * weights[idx++];
				average /= 4. * M_PI;
				for( int y=0 ; y<midHeight ; y++ ) for( int x=0 ; x<midSphere.dimension( y ) ; x++ )  midX[ y*midWidth*Dim + c*midWidth + x ] += streamingLaplacian.average[c]-average;
			}
		}

		/////////////////
		// Out-of-Core //
		/////////////////

		t = Time();
		for( int d=0 ; d<outOfCoreDepth ; d++ )
		{
			for( int c=0 ; c<Dim ; c++ ) outOfCoreSolvers[d].bSquareNorm[c] = outOfCoreSolvers[d].rSquareNorm[c] = outOfCoreSolvers[d].xSquareNorm[c] = 0;
			outOfCoreSolvers[d].inX = outOfCoreSolvers[d].inB = outOfCoreSolvers[d].outX = outOfCoreSolvers[d].outB = NULL;
		}
		outOfCoreSolvers[outOfCoreDepth-1].inX = &midXStream;

		if( v==VCycles.value-1 )
		{
			streamingConverter.child = NULL;
			streamingConverter.outX = out;

			outOfCoreSolvers[0].outX = NULL;
			outOfCoreSolvers[0].outB = NULL;
			outOfCoreSolvers[0].child = &streamingConverter;
			( ( MultiGridRestrictionNode< Real >* )&streamingConverter )->InitProlongation( Progress.set );
		}
		else
		{
			outOfCoreSolvers[0].outX = fullXStream;
			outOfCoreSolvers[0].outB = fullBStream;
			outOfCoreSolvers[0].child = NULL;
			( ( MultiGridRestrictionNode< Real >* )outOfCoreSolvers )->InitProlongation( Progress.set );
		}
		outOfCoreSolvers[outOfCoreDepth-1].SetProlongationIterations( Iters.value );
		outOfCoreSolvers[outOfCoreDepth-1].SolveProlongation();

		if( Verbose.set )
		{
			for( int d=outOfCoreDepth-1 ; d>=0 ; d-- )
			{
				double bSquareNorm = 0 , rSquareNorm = 0;
				for( int c=0 ; c<Dim ; c++ ) bSquareNorm += outOfCoreSolvers[d].bSquareNorm[c] , rSquareNorm += outOfCoreSolvers[d].rSquareNorm[c];
				printf( "\tError[%d x %d]:\t%13.10f\t%g -> %13.10f\n" ,	outOfCoreSolvers[d].columns , outOfCoreSolvers[d].rows , sqrt( bSquareNorm ) , sqrt( rSquareNorm ) , sqrt( rSquareNorm/bSquareNorm ) );
			}
			printf( "Out-of-Core Prolongation Time: %f\n" , Time() - t );
			printf( "V-Cycle Time: %f\n" , Time( ) - vTime );
		}
	}
	delete[] inCoreSolvers;
	delete[] outOfCoreSolvers;
	if( fullXStream ) delete fullXStream;
	if( fullBStream ) delete fullBStream;
	FreeArray( midX );
	FreeArray( midB );
	FreeArray( lowB );
	FreeArray( lowX );
}
template< class Real , class IOReal , class LabelType >
int Execute( int height , bool resample )
{
	StreamingGrid *in , *out , *labels=NULL;
	JointRImageSampler< LabelType , Real > *jointIn = NULL;
	int outWidth , outHeight;

	GetReadSize< Real >( In.value , outWidth , outHeight );
	if( resample )
	{
#if 1	// Try to resample using linear interpolation but without blending across the seams.
		if( Labels.set )
		{
			jointIn = new JointRImageSampler< LabelType , Real >( Labels.value , In.value , 2*height , height , true , false , 0 , NULL );
			labels = jointIn->nearestChild;
			in     = jointIn->averageChild;
		}
		else in = new RImageSampler< Real >( In.value , 2*height , height , false , true , false , 0 , NULL );
#else
		in = new RImageSampler< Real >( In.value , 2*height , height , true , true , false , 0 , NULL );
		if( Labels.set ) labels = new RImageSampler< LabelType >( Labels.value , 2*height , height , true  , true , false , 0 , NULL );
#endif
	}
	else
	{
		int w;
		in = GetReadStream< Real >( In.value , w , height , true , false , false , NULL );
		if( Labels.set ) labels = GetReadStream< LabelType >( Labels.value , w , height , true , false , false , NULL );
	}
	if( resample )
		out = new WImageSampler< Real >( Out.value , 2*height , height , outWidth , outHeight , false , false , false , false , Quality.value , NULL );
	else
		out = GetWriteStream< Real >( Out.value , 2*height , height , false , false , Quality.value , NULL );

	SphericalStencilTable< double > stencilTable;
	double t = Time();
	if( StencilIO.set )
	{
		if( !stencilTable.read( StencilIO.value ) ) fprintf( stderr , "Failed to read stencil table: %s\n" , StencilIO.value );
		if( stencilTable.Init( height , 2*height , 8 , true , true ) )
		{
			if( !stencilTable.write( StencilIO.value ) ) fprintf( stderr , "Failed to write stencil table: %s\n" , StencilIO.value );
		}
		else fprintf( stderr , "Failed to set stencil-table: %d x %d\n" , height , 2*height );
	}
	else stencilTable.Init( height , 2*height , 8 , true , true );
	printf( "Spherical Stencil Table set in: %f\n" , Time() - t );
	if( UnknownIndex.set )
	{
		LabelType unknownIndex[] = { LabelType( UnknownIndex.values[0] ) , LabelType( UnknownIndex.values[1] ) , LabelType( UnknownIndex.values[2] ) };
		StitchImage< Real , IOReal , LabelType , 3 >( height , 8 , &stencilTable , in , labels , out , Iters.value , IWeight.value , GScale.value , GWeight.value , unknownIndex );
	}
	else
		StitchImage< Real , IOReal , LabelType , 3 >( height , 8 , &stencilTable , in , labels , out , Iters.value , IWeight.value , GScale.value , GWeight.value , NULL );

	delete in;
	delete out;
	if( labels ) delete labels;
	if( jointIn ) delete jointIn;
	return EXIT_SUCCESS;
}
int main( int argc , char* argv[] )
{

	bool resample = true;
#if 0
	unsigned int control_word_x87  , control_word_sse2;
	__control87_2( _DN_FLUSH , _MCW_DN , &control_word_x87 , &control_word_sse2 );
#endif
	cmdLineParse( argc-1 , &argv[1] , sizeof(params) / sizeof(cmdLineReadable*) , params , 0 );
	char valueString[1024];
	for( int i=0 ; i<sizeof(params) / sizeof(cmdLineReadable*) ; i++ )
		if( params[i]->set )
		{
			params[i]->writeValue( valueString );
			printf( "\t--%s %s\n" , params[i]->name , valueString );
		}

	if( !In.set || !Out.set )
	{
		ShowUsage(argv[0]);
		return EXIT_FAILURE;
	}
	char tmpDir[1024];
	if( TempDir.set) sprintf( tmpDir , "TMP=%s" , TempDir.value );
	else             sprintf( tmpDir , "TMP=" );
	_putenv( tmpDir );

	int w , h;
	GetReadSize< float >( In.value , w , h );
	int height = 1;
	while( 2*height<w || height<h ) height <<= 1;
	if( 2*height==w && height==h ) resample = false;
	if( resample ) printf( "%d x %d -> %d x %d\n" , w , h , 2*height , height );

	IWeight.value = -IWeight.value;
	IWeight.value *= height / 2;
	IWeight.value *= (height*2) / 2;
	IWeight.value /= 4. * M_PI;
	double t = Time( );
	if( HighPrecision.set )
		if( HDRLabels.set ) Execute< double , double , __int16       >( height , resample );
		else                Execute< double , double , unsigned char >( height , resample );
	else
		if( HDRLabels.set ) Execute< float  , half  , __int16       >( height , resample );
		else                Execute< float  , half  , unsigned char >( height , resample );

	size_t current,peak;
	WorkingSetInfo( current , peak );
	printf( "Running Time: %f\n" , Time()-t );
	printf( "Peak working set: %d MB\n" , peak>>20 );
	return EXIT_SUCCESS;
}
