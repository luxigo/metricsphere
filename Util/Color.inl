template<class Real>
Color<Real>::Color(void)					{ c[0]=c[1]=c[2]=0;   }
template<class Real>
Color<Real>::Color(Real clr)				{ c[0]=c[1]=c[2]=clr; }
template<class Real>
Color<Real>::Color(Real c0,Real c1,Real c2)	{ c[0]=c0;	c[1]=c1;	c[2]=c2; }

template<class Real>
template<class Real2>
Color<Real>::Color(const Color<Real2>& clr)
{
	c[0]=Real(clr[0]);
	c[1]=Real(clr[1]);
	c[2]=Real(clr[2]);
}

template<class Real>
Real Color<Real>::luminance(void) const		{return Real(c[0]*0.3+c[1]*0.59+c[2]*0.11);}

template<class Real>
Real& Color<Real>::operator[] (int i) {return c[i];}
template<class Real>
const Real& Color<Real>::operator[] (int i) const {return c[i];}

template<class Real>
Color<Real> Color<Real>::operator + (const Color& clr) const { return Color(c[0]+clr.c[0],c[1]+clr.c[1],c[2]+clr.c[2]); };
template<class Real>
Color<Real> Color<Real>::operator - (const Color& clr) const { return Color(c[0]-clr.c[0],c[1]-clr.c[1],c[2]-clr.c[2]); };
template<class Real>
Color<Real> Color<Real>::operator - (void            ) const { return Color(-c[0],-c[1],-c[2]); };
template<class Real>
Color<Real> Color<Real>::operator * (const Real& s   ) const { return Color(c[0]*s,c[1]*s,c[2]*s); };
template<class Real>
Color<Real> Color<Real>::operator / (const Real& s   ) const { return Color(c[0]/s,c[1]/s,c[2]/s); };
// An essential method for getting the conjugate gradient solver to work
template<class Real>
Real  Color<Real>::operator * (const Color& clr) const { return c[0]*clr[0]+c[1]*clr[1]+c[2]*clr[2]; };

template<class Real>
Color<Real>& Color<Real>::operator += (const Color& clr)
{
	c[0]+=clr.c[0];
	c[1]+=clr.c[1];
	c[2]+=clr.c[2];
	return *this;
};
template<class Real>
Color<Real>& Color<Real>::operator -= (const Color& clr)
{
	c[0]-=clr.c[0];
	c[1]-=clr.c[1];
	c[2]-=clr.c[2];
	return *this;
};
template<class Real>
Color<Real>& Color<Real>::operator *= (const Real& s)
{
	c[0]*=s;
	c[1]*=s;
	c[2]*=s;
	return *this;
};
template<class Real>
Color<Real>& Color<Real>::operator /= (const Real& s)
{
	c[0]/=s;
	c[1]/=s;
	c[2]/=s;
	return *this;
};
template<class Real>
bool Color<Real>::operator == (const Color<Real>& clr)	const
{
	return c[0]==clr[0] && c[1]==clr[1] && c[2]==clr[2];
}
template<class Real>
bool Color<Real>::operator != (const Color<Real>& clr) const
{
	return c[0]!=clr[0] || c[1]!=clr[1] || c[2]!=clr[2];
}
