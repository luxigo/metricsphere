////////////////////
// GradientStream //
////////////////////
template<class InReal,class OutReal,int Channels>
GradientStream<InReal,OutReal,Channels>::GradientStream(StreamingGrid* pixels,StreamingGrid* labels,bool pad,double* average)
{
	_pad=pad;
	_average=average;
	_hh=_h=pixels->rows();
	_ww=_w=pixels->rowSize()/(Channels*sizeof(InReal));
	if(pixels->rowSize()!=_w*Channels*sizeof(InReal))	fprintf(stderr,"Pixel width failure: %d != %d\n",pixels->rowSize(),_w*Channels*sizeof(InReal)),	exit(0);
	if(_h!=labels->rows())								fprintf(stderr,"Label height failure: %d != %d\n",_h,labels->rows()),							exit(0);
	if(labels->rowSize()!=_w*Channels*sizeof(InReal))	fprintf(stderr,"Label width failure: %d != %d\n",labels->rowSize(),_w*Channels*sizeof(InReal)),	exit(0);

	if(_pad)	_hh++	,	_ww++;
	_pixels = pixels;
	_labels = labels;
	_previousPixelsRow	= (InReal* )malloc(sizeof(InReal )*Channels*_w);
	_previousLabelsRow	= (InReal* )malloc(sizeof(InReal )*Channels*_w);
	_dx					= (OutReal*)malloc(sizeof(OutReal)*Channels*(_ww-1));
	_dy					= (OutReal*)malloc(sizeof(OutReal)*Channels*(_ww  ));

	for(int c=0;c<Channels;c++)	_average[c]=0;

	_current=0;
	_pixels->reset(true,1);
	_labels->reset(true,1);

	memcpy(_previousPixelsRow,(*_pixels)[_current],sizeof(InReal)*Channels*_w);
	memcpy(_previousLabelsRow,(*_labels)[_current],sizeof(InReal)*Channels*_w);
	for(int x=0;x<_w;x++)	for(int c=0;c<Channels;c++)	_average[c]+=_previousPixelsRow[x+_w*c];

	InReal l1[Channels],l2[Channels];

	// Set the partials with respect to x
	for(int x=0;x<_w-1;x++)
	{
		bool useGradient=true;
		for(int c=0;c<Channels;c++)	l1[c] = _previousLabelsRow[x+_w*c]	,	l2[c] = _previousLabelsRow[x+1+_w*c];
		for(int c=0;c<Channels && useGradient;c++)	useGradient &= (l1[c]==l2[c]);
		if(useGradient)		for(int c=0;c<Channels;c++)	_dx[x*Channels+c] = OutReal(_previousPixelsRow[x+1+_w*c])-OutReal(_previousPixelsRow[x+_w*c]);
		else				for(int c=0;c<Channels;c++)	_dx[x*Channels+c] = 0;
	}
	if(_pad)	for(int c=0;c<Channels;c++)	_dx[(_w-1)*Channels+c] = OutReal(0)-OutReal(_previousPixelsRow[(_w-1)+_w*c]);
	_pixels->advance();
	_labels->advance();
	InReal* newPixelsRow = (InReal*)(*_pixels)[_current+1];
	InReal* newLabelsRow = (InReal*)(*_labels)[_current+1];

	// Set the partials with respect to y
	for(int x=0;x<_w;x++)
	{
		bool useGradient=true;
		for(int c=0;c<Channels;c++)	l1[c] = newLabelsRow[x+_w*c]	,	l2[c] = _previousLabelsRow[x+_w*c];
		for(int c=0;c<Channels && useGradient;c++)	useGradient &= (l1[c]==l2[c]);
		if(useGradient)		for(int c=0;c<Channels;c++)	_dy[x*Channels+c] = OutReal(newPixelsRow[x+_w*c])-OutReal(_previousPixelsRow[x+_w*c]);
		else				for(int c=0;c<Channels;c++)	_dy[x*Channels+c] = 0;
	}
	if(_pad)	for(int c=0;c<Channels;c++)	_dy[_w*Channels+c] = OutReal(0);

	memcpy(_previousPixelsRow,newPixelsRow,sizeof(InReal)*Channels*_w);
	memcpy(_previousLabelsRow,newLabelsRow,sizeof(InReal)*Channels*_w);
	for(int x=0;x<_w;x++)	for(int c=0;c<Channels;c++)	_average[c]+=_previousPixelsRow[x+_w*c];

	_DX=new DX(this);
	_DY=new DY(this);
}
template<class InReal,class OutReal,int Channels>
GradientStream<InReal,OutReal,Channels>::~GradientStream(void)
{
	free(_previousPixelsRow);
	free(_previousLabelsRow);
	free(_dx);
	free(_dy);
	delete _DX;
	delete _DY;
}
template<class InReal,class OutReal,int Channels>
void GradientStream<InReal,OutReal,Channels>::advance(void)
{
	_current++;
	if(_current<_hh)
	{
		InReal l1[Channels],l2[Channels];
		// Set the partials with respect to x
		for(int x=0;x<_w-1;x++)
		{
			bool useGradient=true;
			for(int c=0;c<Channels;c++)	l1[c] = _previousLabelsRow[x+_w*c]	,	l2[c] = _previousLabelsRow[x+1+_w*c];
			for(int c=0;c<Channels && useGradient;c++)	useGradient &= (l1[c]==l2[c]);
			if(useGradient)		for(int c=0;c<Channels;c++)	_dx[x*Channels+c] = OutReal(_previousPixelsRow[x+1+_w*c])-OutReal(_previousPixelsRow[x+_w*c]);
			else				for(int c=0;c<Channels;c++)	_dx[x*Channels+c] = 0;
		}
		if(_pad)	for(int c=0;c<Channels;c++)	_dx[(_w-1)*Channels+c] = OutReal(0)-OutReal(_previousPixelsRow[(_w-1)+_w*c]);

		if(_current<_h)
		{
			_pixels->advance();
			_labels->advance();
		}
		if(_current+1<_h)
		{
			InReal* newPixelsRow = (InReal*)(*_pixels)[_current+1];
			InReal* newLabelsRow = (InReal*)(*_labels)[_current+1];

			// Set the partials with respect to y
			for(int x=0;x<_w;x++)
			{
				bool useGradient=true;
				for(int c=0;c<Channels;c++)	l1[c] = newLabelsRow[x+_w*c]	,	l2[c] = _previousLabelsRow[x+_w*c];
				for(int c=0;c<Channels && useGradient;c++)	useGradient &= (l1[c]==l2[c]);
				if(useGradient)		for(int c=0;c<Channels;c++)	_dy[x*Channels+c] = OutReal(newPixelsRow[x+_w*c])-OutReal(_previousPixelsRow[x+_w*c]);
				else				for(int c=0;c<Channels;c++)	_dy[x*Channels+c] = 0;
			}
			if(_pad)	for(int c=0;c<Channels;c++)	_dy[_w*Channels+c] = OutReal(0);

			memcpy(_previousPixelsRow,newPixelsRow,sizeof(InReal)*Channels*_w);
			memcpy(_previousLabelsRow,newLabelsRow,sizeof(InReal)*Channels*_w);
			for(int x=0;x<_w;x++)	for(int c=0;c<Channels;c++)	_average[c]+=_previousPixelsRow[x+_w*c];
		}
		else if(_current+1==_h && _pad)
		{
			// Set the partials with respect to y
			for(int x=0;x<_w;x++)	for(int c=0;c<Channels;c++)	_dy[x*Channels+c] = OutReal(0)-OutReal(_previousPixelsRow[x+_w*c]);
			for(int c=0;c<Channels;c++)	_dy[_w*Channels+c] = OutReal(0);

			memset(_previousPixelsRow,0,sizeof(InReal)*Channels*_w);
			memset(_previousLabelsRow,0,sizeof(InReal)*Channels*_w);
		}
	}
	if(_current==_hh-1)	for(int c=0;c<Channels;c++)	_average[c]/=_w	,	_average[c]/=_h;
}
template<class InReal,class OutReal,int Channels>
void* GradientStream<InReal,OutReal,Channels>::dX(int idx)
{
	if(idx!=_current || idx>=_hh)	fprintf(stderr,"GradientStream::dX index out of bounds: %d != %d\n",idx,_current),	exit(0);
	return _dx;
}
template<class InReal,class OutReal,int Channels>
void* GradientStream<InReal,OutReal,Channels>::dY(int idx)
{
	if(idx!=_current || idx>=_hh-1)	fprintf(stderr,"GradientStream::dY index out of bounds: %d != %d\n",idx,_current),	exit(0);
	return _dy;
}
template<class InReal,class OutReal,int Channels>
StreamingGrid* GradientStream<InReal,OutReal,Channels>::getDXStream(void)	{return _DX;}
template<class InReal,class OutReal,int Channels>
StreamingGrid* GradientStream<InReal,OutReal,Channels>::getDYStream(void)	{return _DY;}

template<class InReal,class OutReal,int Channels>
int GradientStream<InReal,OutReal,Channels>::width(void)	const	{return _ww;}
template<class InReal,class OutReal,int Channels>
int GradientStream<InReal,OutReal,Channels>::height(void)	const	{return _hh;}
///////////////////////////
// GradientStream::DX/DY //
///////////////////////////
template<class InReal,class OutReal,int Channels>
GradientStream<InReal,OutReal,Channels>::DX::DX(GradientStream<InReal,OutReal,Channels>* gStream)	{	_gStream=gStream;	}
template<class InReal,class OutReal,int Channels>
GradientStream<InReal,OutReal,Channels>::DY::DY(GradientStream<InReal,OutReal,Channels>* gStream)	{	_gStream=gStream;	}
template<class InReal,class OutReal,int Channels>
void GradientStream<InReal,OutReal,Channels>::DX::advance(void)			{							}
template<class InReal,class OutReal,int Channels>
void GradientStream<InReal,OutReal,Channels>::DY::advance(void)			{	_gStream->advance();	}
template<class InReal,class OutReal,int Channels>
void* GradientStream<InReal,OutReal,Channels>::DX::operator[](int idx)	{	return _gStream->dX(idx);	}
template<class InReal,class OutReal,int Channels>
void* GradientStream<InReal,OutReal,Channels>::DY::operator[](int idx)	{	return _gStream->dY(idx);	}
template<class InReal,class OutReal,int Channels>
int GradientStream<InReal,OutReal,Channels>::DX::rows(void)	const	{	return _gStream->height();	}
template<class InReal,class OutReal,int Channels>
int GradientStream<InReal,OutReal,Channels>::DY::rows(void)	const	{	return _gStream->height()-1;}
template<class InReal,class OutReal,int Channels>
int GradientStream<InReal,OutReal,Channels>::DX::rowSize(void)	const	{	return (_gStream->width()-1)*sizeof(OutReal)*Channels;	}
template<class InReal,class OutReal,int Channels>
int GradientStream<InReal,OutReal,Channels>::DY::rowSize(void)	const	{	return (_gStream->width()  )*sizeof(OutReal)*Channels;	}
