/*
Copyright (c) 2008, Michael Kazhdan
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution. 

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#include "Util/Half/half.h"

template< class Real >	bool IsFloatType( void );
template< > inline		bool IsFloatType< half >	( void )	{ return true;  }
template< > inline		bool IsFloatType< float >	( void ) 	{ return true;  }
template< > inline		bool IsFloatType< double >	( void )	{ return true;  }
template< class Real >	bool IsFloatType			( void )	{ return false; }

template<class Real>
class ConversionFactor
{
public:
	static const double Factor;
};
template < > const double ConversionFactor<             char >::Factor = double( (long long(1)<< 8)-1 );
template < > const double ConversionFactor<    unsigned char >::Factor = double( (long long(1)<< 8)-1 );
template < > const double ConversionFactor<          __int16 >::Factor = double( (long long(1)<<16)-1 );
template < > const double ConversionFactor< unsigned __int16 >::Factor = double( (long long(1)<<16)-1 );
template < > const double ConversionFactor<              int >::Factor = double( (long long(1)<<32)-1 );
template < > const double ConversionFactor<     unsigned int >::Factor = double( (long long(1)<<32)-1 );
template< class Real > const double ConversionFactor<Real>::Factor	= 1.0;

template< class In , class Out >
class Converter
{
public:
	static const double ConvertFactor;
};
template< class In , class Out >
const double Converter< In , Out >::ConvertFactor = ConversionFactor<Out>::Factor / ConversionFactor<In>::Factor;

template<class In,class Out>
inline Out ConvertChannel( const In& in )
{
//	printf( "%u -> %g -> %u: %d %d %g\n" , in , double(in)*Converter<In,Out>::ConvertFactor , Out( double(in)*Converter<In,Out>::ConvertFactor) , sizeof( In ) , sizeof( Out ) , Converter<In,Out>::ConvertFactor );
	return Out(double(in)*Converter<In,Out>::ConvertFactor);
}
