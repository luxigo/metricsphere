/* Factorial function */

#include "qhead.h"

extern QELT qone[];

int qfac( x, y )
QELT x[], y[];
{
QELT p[NQ], j[NQ];
long i, n;
double dn, dp, dj;

if( x[0] != 0 )
	{
	qinfin( y );
	mtherr( "qfac", DOMAIN );
	return 0;
	}


qtoe( x, (unsigned short *) &dn );
n = dn;

if( n == 0 )
	{
	qmov( qone, y );
	return 0;
	}

if( n == 1 )
	{
	qmov( qone, y );
	return 0;
	}

if( n > 1754 )
	{
	qinfin(y);
	mtherr( "qfac", OVERFLOW );
	return 0;
	}

/* Cheat by using normal arithmetic */
dp = 1.0;
dj = 1.0;
i = 1;
do
	{
	if( i > 17 )
		goto fmore;
	i += 1;
	dj += 1.0;
	dp *= dj;
	}
while( i < n );

etoq( (unsigned short *) &dp, y );
return 0;


fmore:

etoq( (unsigned short *) &dj, j );
etoq( (unsigned short *) &dp, p );

do
	{
	i += 1;
	qadd( qone, j, j );
	qmuli( j, p, p );
	}
while( i < n );

qmov( p, y );
return 0;
}

