#ifndef GRADIENT_STREAM_INCLUDED
#define GRADIENT_STREAM_INCLUDED
#include "ImageStream.h"


template<class InReal,class OutReal,int Channels>
class GradientStream
{
	bool _pad;
	double* _average;
	int _w,_h,_ww,_hh,_current;
	StreamingGrid *_pixels,*_labels;
	InReal *_previousPixelsRow,*_previousLabelsRow;
	OutReal *_dx,*_dy;
	class DX : public StreamingGrid
	{
		GradientStream* _gStream;
	public:
		DX(GradientStream* gStream);
		void advance(void);
		void* operator[]	(int idx);
		int rows			(void)	const;
		int rowSize			(void)	const;
	};
	class DY : public StreamingGrid
	{
		GradientStream* _gStream;
	public:
		DY(GradientStream* gStream);
		void advance(void);
		void* operator[]	(int idx);
		int rows			(void)	const;
		int rowSize			(void)	const;
	};
	DX* _DX;
	DY* _DY;
public:
	GradientStream(StreamingGrid* pixels,StreamingGrid* labels,bool pad,double* average);
	~GradientStream(void);
	int width(void)		const;
	int height(void)	const;
	void* dX(int idx);
	void* dY(int idx);
	void advance(void);
	StreamingGrid* getDXStream(void);
	StreamingGrid* getDYStream(void);
};
#include "GradientStream.inl"
#endif // GRADIENT_STREAM_INCLUDED