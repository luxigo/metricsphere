#ifndef COLOR_INCLUDED
#define COLOR_INCLUDED

template<class Real>
class Color
{
	Real c[3];
public:
	Color(void);
	Color(Real clr);
	Color(Real c0,Real c1,Real c2);
	template<class Real2>
	Color(const Color<Real2>& clr);

	Real luminance(void) const;

	Real& operator[] (int i);
	const Real& operator[] (int i) const;

	Color operator + (const Color& clr) const;
	Color operator - (const Color& clr) const;
	Color operator - (void            ) const;
	Color operator * (const Real& s   ) const;
	Color operator / (const Real& s   ) const;
	// An essential method for getting the conjugate gradient solver to work
	Real  operator * (const Color& clr) const;

	Color& operator += (const Color& clr);
	Color& operator -= (const Color& clr);
	Color& operator *= (const Real& s);
	Color& operator /= (const Real& s);
	bool operator == (const Color<Real>& clr) const;
	bool operator != (const Color<Real>& clr) const;
};

#include "Color.inl"
#endif // COLOR_INCLUDED
