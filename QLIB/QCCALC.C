/*	calc.c, complex variable version */
/* Keyboard command interpreter		*/
/* Copyright 1985 by S. L. Moshier	*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "qhead.h"

#ifndef USE_READLINE
#define USE_READLINE 0
#endif
#if USE_READLINE
char *readline();
void add_history();
static char *line_read = (char *)NULL;
#endif

/* Define nonzero for extra functions, making the program bigger.  */
#ifndef MOREFUNS
#define MOREFUNS 1
#endif

/*
typedef struct
	{
	QELT r[NQ];
	QELT i[NQ];
	}qcmplx;
*/

/*
 *#include "config.h"
 */

/* length of command line: */
#define LINLEN 128

#define XON 0x11
#define XOFF 0x13
#define NOTABS 1

#define SALONE 1
#define DECPDP 0
#define INTLOGIN 0
#define INTHELP 1
#ifndef TRUE
#define TRUE 1
#endif

/* initialize printf: */
#define INIPRINTF 0

#if DECPDP
#define TRUE 1
#endif


static char idterp[] = {
"\n\nSteve Moshier's command interpreter V1.3\n"};
#define ISLOWER(c) ((c >= 'a') && (c <= 'z'))
#define ISUPPER(c) ((c >= 'A') && (c <= 'Z'))
#define ISALPHA(c) (ISLOWER(c) || ISUPPER(c))
#define ISDIGIT(c) ((c >= '0') && (c <= '9'))
#define ISATF(c) (((c >= 'a')&&(c <= 'f')) || ((c >= 'A')&&(c <= 'F')))
#define ISXDIGIT(c) (ISDIGIT(c) || ISATF(c))
#define ISOCTAL(c) ((c >= '0') && (c < '8'))
#define ISALNUM(c) (ISALPHA(c) || (ISDIGIT(c))

FILE *fopen();
/* I/O log file: */
static char *savnam = 0;
static FILE *savfil = 0;


#include "qcalc.h"

/* space for values of variables */
#if MOREFUNS
#define NVS 41
#else
#define NVS 34
#endif
static QELT vs[NVS][2*NQ];

/*	the symbol table of temporary variables: */

#define NTEMP 4
struct varent temp[NTEMP] = {
{"T",	OPR | TEMP, &vs[14][0],},
{"T",	OPR | TEMP, &vs[15][0],},
{"T",	OPR | TEMP, &vs[16][0],},
{"\0",	OPR | TEMP, &vs[17][0]},
};

/*	the symbol table of operators		*/
/* EOL is interpreted on null, newline, or ;	*/
struct symbol oprtbl[] = {
{"BOL",		OPR | BOL,	0,},
{"EOL",		OPR | EOL,	0,},
{"-",		OPR | UMINUS,	8,},
/*"~",		OPR | COMP,	8,*/
{",",		OPR | EOE,	1,},
{"=",		OPR | EQU,	2,},
/*"|",		OPR | LOR,	3,*/
/*"^",		OPR | LXOR,	4,*/
/*"&",		OPR | LAND,	5,*/
{"+",		OPR | PLUS,	6,},
{"-",		OPR | MINUS, 6,},
{"*",		OPR | MULT,	7,},
{"/",		OPR | DIV,	7,},
/*"%",		OPR | MOD,	7,*/
{"(",		OPR | LPAREN,	11,},
{")",		OPR | RPAREN,	11,},
{"\0",		ILLEG, 0},
};

#define NOPR 8

/*	the symbol table of indirect variables: */
extern QELT qpi[], qone[];
QELT qcpi[2*NQ];
QELT qci[2*NQ];

struct varent indtbl[] = {
{"a",		VAR | IND,	&vs[32][0],},
{"b",		VAR | IND,	&vs[31][0],},
{"c",		VAR | IND,	&vs[30][0],},
{"d",		VAR | IND,	&vs[29][0],},
{"e",		VAR | IND,	&vs[28][0],},
{"f",		VAR | IND,	&vs[27][0],},
{"g",		VAR | IND,	&vs[26][0],},
{"h",		VAR | IND,	&vs[25][0],},
#if MOREFUNS
{"j",		VAR | IND,	&vs[35][0],},
{"k",		VAR | IND,	&vs[36][0],},
{"l",		VAR | IND,	&vs[37][0],},
{"m",		VAR | IND,	&vs[38][0],},
{"n",		VAR | IND,	&vs[39][0],},
{"o",		VAR | IND,	&vs[40][0],},
#endif
{"p",		VAR | IND,	&vs[24][0],},
{"q",		VAR | IND,	&vs[23][0],},
{"r",		VAR | IND,	&vs[22][0],},
{"s",		VAR | IND,	&vs[13][0],},
{"t",		VAR | IND,	&vs[21][0],},
{"u",		VAR | IND,	&vs[20][0],},
{"v",		VAR | IND,	&vs[19][0],},
{"w",		VAR | IND,	&vs[18][0],},
{"x",		VAR | IND,	&vs[10][0],},
{"y",		VAR | IND,	&vs[11][0],},
{"z",		VAR | IND,	&vs[12][0],},
{"pi",		VAR | IND,	&qcpi[0],},
{"i",		VAR | IND,	&qci[0],},
{"\0",		ILLEG,		0},
};

/*	the symbol table of constants:	*/

#define NCONST 10
struct varent contbl[NCONST] = {
{"C",CONST,&vs[0][0],},
{"C",CONST,&vs[1][0],},
{"C",CONST,&vs[2][0],},
{"C",CONST,&vs[3][0],},
{"C",CONST,&vs[4][0],},
{"C",CONST,&vs[5][0],},
{"C",CONST,&vs[6][0],},
{"C",CONST,&vs[7][0],},
{"C",CONST,&vs[8][0],},
{"\0",CONST,&vs[9][0]},
};

/* the symbol table of string variables: */

static char strngs[4][40];

#define NSTRNG 5
struct strent strtbl[NSTRNG] = {
#if DECPDP
{&strngs[0][0], VAR | STRING, &strngs[0][0],},
{&strngs[1][0], VAR | STRING, &strngs[1][0],},
{&strngs[2][0], VAR | STRING, &strngs[2][0],},
{&strngs[3][0], VAR | STRING, &strngs[3][0],},
#else
{&strngs[0][0], VAR | STRING, &strngs[0][0],},
{&strngs[1][0], VAR | STRING, &strngs[1][0],},
{&strngs[2][0], VAR | STRING, &strngs[2][0],},
{&strngs[3][0], VAR | STRING, &strngs[3][0],},
#endif
{"\0",ILLEG,0,},
};


/* Help messages */
#if INTHELP
static char *intmsg[] = {
"?",
"Unkown symbol",
"Expression ends in illegal operator",
"Precede ( by operator",
")( is illegal",
"Unmatched )",
"Missing )",
"Illegal left hand side",
"Missing symbol",
"Must assign to a variable",
"Divide by zero",
"Missing symbol",
"Missing operator",
"Precede quantity by operator",
"Quantity preceded by )",
"Function syntax",
"Too many function args",
"No more temps",
"Arg list"
};
#endif

/*	the symbol table of functions:	*/
#if SALONE
#ifdef ANSIPROT
int hex(QELT *), cmdh(void), cmdhlp(void), remark(char *), qsys(char *);
int qsave(char *);
int hexinput(QELT *, QELT *, QELT *), cmddig(QELT *);
int cmddm(QELT *), cmdtm(QELT *), cmdem(QELT *);
int take(char *), mxit(void), bits(qcmplx *), cmp(QELT *, QELT *, QELT *);
int intcvts (QELT *, QELT *), todouble (QELT *, QELT *);
struct symbol *parser(void);	/* parser returns pointer to symbol */
int init(void), abmac(void), zgets(char *, int), prhlst(struct symbol *);
int qmovz(QELT *, QELT *);
int zpdtr(QELT *, QELT *, QELT *), zpdtri(QELT *, QELT *, QELT *);
int makec( qcmplx *, qcmplx *, qcmplx * );
int getre( qcmplx *, qcmplx * );
int getim( qcmplx *, qcmplx * );
int qcarg( qcmplx *, qcmplx * );
int qccabs( qcmplx *, qcmplx * );
int qcconj( qcmplx *, qcmplx * );
int qchyp1f1( qcmplx *, qcmplx *, qcmplx *, qcmplx *);
int qcerf( qcmplx *, qcmplx * );
int qcerfw( qcmplx *, qcmplx * );
#else /* not ANSIPROT */
int hex(), cmdh(), cmdhlp();
int cmddm(), cmdtm(), cmdem();
int qtoe(), etoq();
int take(), mxit(), bits();
int cmddig(), qfloor();
int qsave(), qsys();
int makec(), getre(), getim(), qccabs(), qcarg(), qcconj();
int qcsqrt(), qclog(), qcexp(), qcsin(), qccos(), qcatan();
int qcasin(), qcacos(), qctan(), qccot(), qcpow();
int qfac(), qcgamma(), qclgam();
struct symbol *parser();	/* parser returns pointer to symbol */
int qmovz(), shup1(), shdn1(), qtoe113(), e113toq(), qtoe64(), e64toq();
int e24toq(), qtoe24();
int init(), abmac(), zgets(), prhlst(), asctoq();
int qcmov(), qtoasc(), qcneg(), qcadd(), qcsub(), qcmul(), qcdiv();
int qcabs(), qatn2(), shdn1(), shup1(), qmovz();
int qchyp1f1(), qcerf(), qcerfw(), qcpolylog();
#endif /* not ANSIPROT */

static int zpolylog();

struct funent funtbl[] = {
{"h",		OPR | FUNC, cmdh,},
{"help",	OPR | FUNC, cmdhlp,},
{"hex",		OPR | FUNC, hex,},
/*"view",	OPR | FUNC, view,*/
{"cmplx",	OPR | FUNC, makec,},
{"re",		OPR | FUNC, getre,},
{"im",		OPR | FUNC, getim,},
{"arg",		OPR | FUNC, qcarg,},
{"abs",		OPR | FUNC, qccabs,},
{"conj",	OPR | FUNC, qcconj,},
{"acos",	OPR | FUNC, qcacos,},
{"asin",	OPR | FUNC, qcasin,},
{"atan",	OPR | FUNC, qcatan,},
{"cos",		OPR | FUNC, qccos,},
{"cot",		OPR | FUNC, qccot,},
{"exp",		OPR | FUNC, qcexp,},
{"floor",	OPR | FUNC, qfloor,},
{"fac",		OPR | FUNC, qfac,},
{"gamma",	OPR | FUNC, qcgamma,},
{"lgamma",	OPR | FUNC, qclgam,},
{"confhyp",	OPR | FUNC, qchyp1f1,},
{"erf",		OPR | FUNC, qcerf,},
{"erfw",	OPR | FUNC, qcerfw,},
/*
"jv",		OPR | FUNC, qjn,
"yn",		OPR | FUNC, qyn,
"logten",	OPR | FUNC, qlog10,
"expten",	OPR | FUNC, qexp10,
*/
{"log",		OPR | FUNC, qclog,},
/*
"ndtr",		OPR | FUNC, ndtr,
"ndtri",	OPR | FUNC, ndtri,
*/
{"pow",		OPR | FUNC, qcpow,},
{"sin",		OPR | FUNC, qcsin,},
{"sqrt",	OPR | FUNC, qcsqrt,},
{"tan",		OPR | FUNC, qctan,},
{"polylog",	OPR | FUNC, zpolylog,},
{"bits",	OPR | FUNC, bits,},
{"digits",	OPR | FUNC, cmddig,},
{"dm",		OPR | FUNC, cmddm,},
{"tm",		OPR | FUNC, cmdtm,},
{"em",		OPR | FUNC, cmdem,},
{"take",	OPR | FUNC | COMMAN, take,},
{"save",	OPR | FUNC | COMMAN, qsave,},
{"system",	OPR | FUNC | COMMAN, qsys,},
{"exit",	OPR | FUNC, mxit,},
{"\0",		OPR | FUNC,	0},
};

/*	the symbol table of key words */
struct funent keytbl[] = {
{"\0",		ILLEG,	0},
};
#endif /* SALONE */

/* Number of decimals to display */
#if NQ > 12
#define DEFDIS 70
#else
#define DEFDIS 41
#endif
static int ndigits = DEFDIS;

/* Menu stack */
struct funent *menstk[5];
int menptr = 0;

/* Take file stack */
FILE *takstk[10];
int takptr = -1;

/* size of the expression scan list: */
#define NSCAN 20

/* previous token, saved for syntax checking: */
struct symbol *lastok = 0;

/*	variables used by parser: */
static char str[128] = {0};
int uposs = 0;		/* possible unary operator */
double nc = 0;		/*	numeric value of symbol		*/
static QELT qnc[2*NQ] = {0};
char lc[40] = { '\n' };	/*	ASCII string of token	symbol	*/
static char line[LINLEN] = { '\n','\0' };	/* input command line */
static char maclin[LINLEN] = { '\n','\0' };	/* macro command */
char *interl = line;		/* pointer into line */
extern char *interl;
static int maccnt = 0;	/* number of times to execute macro command */
static int comptr = 0;	/* comma stack pointer */
static qcmplx comstk[6];	/* comma argument stack */
static int narptr = 0;	/* pointer to number of args */
static int narstk[5] = {0};	/* stack of number of function args */

/*							main()		*/

/*	Entire program starts here	*/

int main()
{

/*	the scan table:			*/

/*	array of pointers to symbols which have been parsed:	*/
struct symbol *ascsym[NSCAN];

/*	current place in ascsym:			*/
register struct symbol **as;

/*	array of attributes of operators parsed:		*/
int ascopr[NSCAN];

/*	current place in ascopr:			*/
register int *ao;

#if LARGEMEM
/*	array of precedence levels of operators:		*/
long asclev[NSCAN];
/*	current place in asclev:			*/
long *al;
long symval;	/* value of symbol just parsed */
#else
int asclev[NSCAN];
int *al;
int symval;
#endif

qcmplx acc;	/* the accumulator, for arithmetic */
int accflg;	/* flags accumulator in use	*/
qcmplx val;	/* value to be combined into accumulator */
register struct symbol *psym;	/* pointer to symbol just parsed */
struct varent *pvar;	/* pointer to an indirect variable symbol */
struct funent *pfun;	/* pointer to a function symbol */
int att;	/* attributes of symbol just parsed */
int i;		/* counter	*/
int offset;	/* parenthesis level */
int lhsflg;	/* kluge to detect illegal assignments */
int errcod;	/* for syntax error printout */


/* Perform general initialization */

init();

menstk[0] = &funtbl[0];
menptr = 0;
cmdhlp();		/* print out list of symbols */


/*	Return here to get next command line to execute	*/
getcmd:

/* initialize registers and mutable symbols */

accflg = 0;	/* Accumulator not in use				*/
qclear( acc.r );	/* Clear the accumulator			*/
qclear( acc.i );
offset = 0;	/* Parenthesis level zero				*/
comptr = 0;	/* Start of comma stack					*/
narptr = -1;	/* Start of function arg counter stack	*/

psym = (struct symbol *)&contbl[0];
for( i=0; i<NCONST; i++ )
	{
	psym->attrib = CONST;	/* clearing the busy bit */
	++psym;
	}
psym = (struct symbol *)&temp[0];
for( i=0; i<NTEMP; i++ )
	{
	psym->attrib = VAR | TEMP;	/* clearing the busy bit */
	++psym;
	}

psym = (struct symbol *)&strtbl[0];
for( i=0; i<NSTRNG; i++ )
	{
	psym->attrib = STRING | VAR;
	++psym;
	}

/*	List of scanned symbols is empty:	*/
as = &ascsym[0];
*as = 0;
--as;
/*	First item in scan list is Beginning of Line operator	*/
ao = &ascopr[0];
*ao = oprtbl[0].attrib & 0xf;	/* BOL */
/*	value of first item: */
al = &asclev[0];
*al = oprtbl[0].sym;

lhsflg = 0;		/* illegal left hand side flag */
psym = &oprtbl[0];	/* pointer to current token */
qmov( qpi, qcpi );
qclear( &qcpi[NQ] );
qclear( qci );
qmov( qone, &qci[NQ] );

/*	get next token from input string	*/

gettok:
lastok = psym;		/* last token = current token */
psym = parser();	/* get a new current token */
/*printf( "%s attrib %7o value %7o\n", psym->spel, psym->attrib & 0xffff,
		psym->sym );*/

/* Examine attributes of the symbol returned by the parser	*/
att = psym->attrib;
if( att & ILLEG )
	{
	errcod = 1;
	goto synerr;
	}

/*	Push functions onto scan list without analyzing further */
if( att & FUNC )
	{
	/* A command is a function whose argument is
	 * a pointer to the rest of the input line.
	 * A second argument is also passed: the address
	 * of the last token parsed.
	 */
	if( att & COMMAN )
		{
		pfun = (struct funent *)psym;
		( *(pfun->fun))( interl, lastok );
		abmac();	/* scrub the input line */
		goto getcmd;	/* and ask for more input */
		}
	++narptr;	/* offset to number of args */
	narstk[narptr] = 0;
	i = lastok->attrib & 0xffff; /* attrib=short, i=int */
	if( ((i & OPR) == 0)
			|| (i == (OPR | RPAREN))
			|| (i == (OPR | FUNC)) )
		{
		errcod = 15;
		goto synerr;
		}

	++lhsflg;
	++as;
	*as = psym;
	++ao;
	*ao = FUNC;
	++al;
	*al = offset + UMINUS;
	goto gettok;
	}

/* deal with operators */
if( att & OPR )
	{
	att &= 0xf;
	/* expression cannot end with an operator other than
	 * (, ), BOL, or a function
	 */
	if( (att == RPAREN) || (att == EOL) || (att == EOE))
		{
		i = lastok->attrib & 0xffff; /* attrib=short, i=int */
		if( (i & OPR) 
			&& (i != (OPR | RPAREN))
			&& (i != (OPR | LPAREN))
			&& (i != (OPR | FUNC))
			&& (i != (OPR | BOL)) )
				{
				errcod = 2;
				goto synerr;
				}
		}
	++lhsflg;	/* any operator but ( and = is not a legal lhs */

/*	operator processing, continued */

	switch( att )
		{
 	case EOE:
		lhsflg = 0;
		break; 
	case LPAREN:
		/* ( must be preceded by an operator of some sort. */
		if( ((lastok->attrib & OPR) == 0) )
			{
			errcod = 3;
			goto synerr;
			}
		/* also, a preceding ) is illegal */
		if( (unsigned short )lastok->attrib
		   == (unsigned short)(OPR|RPAREN))
			{
			errcod = 4;
			goto synerr;
			}
		/* Begin looking for illegal left hand sides: */
		lhsflg = 0;
		offset += RPAREN;	/* new parenthesis level */
		goto gettok;
	case RPAREN:
		offset -= RPAREN;	/* parenthesis level */
		if( offset < 0 )
			{
			errcod = 5;	/* parenthesis error */
			goto synerr;
			}
		goto gettok;
	case EOL:
		if( offset != 0 )
			{
			errcod = 6;	/* parenthesis error */
			goto synerr;
			}
		break;
	case EQU:
		if( --lhsflg )	/* was incremented before switch{} */
			{
			errcod = 7;
			goto synerr;
			}
	case UMINUS:
	case COMP:
		goto pshopr;	/* evaluate right to left */
	default:	;
		}


/*	evaluate expression whenever precedence is not increasing	*/

symval = psym->sym + offset;

while( symval <= *al )
	{
	/* if just starting, must fill accumulator with last
	 * thing on the line
	 */
	if( (accflg == 0) && (as >= ascsym) && (((*as)->attrib & FUNC) == 0 ))
		{
		pvar = (struct varent *)*as;
		qcmov( (qcmplx *) pvar->value, &acc );
		--as;
		accflg = 1;
		}

/* handle beginning of line type cases, where the symbol
 * list ascsym[] may be empty.
 */
	switch( *ao )
		{
	case BOL:	
		qtoasc( acc.r, str, ndigits );
		printf( "%s\n", str ); /* This is the answer */
		qtoasc( acc.i, str, ndigits );
		if( acc.i[0] == 0 )
			str[0] = '+';
		printf( "%si\n", str ); /* imaginary part */
		if( savfil )
			fprintf( savfil, "%s\n", str );
		goto getcmd;	/* all finished */
	case UMINUS:
		qcneg( &acc );
		goto nochg;
/*
	case COMP:
		acc = ~acc;
		goto nochg;
*/
	default:	;
		}
/* Now it is illegal for symbol list to be empty,
 * because we are going to need a symbol below.
 */
	if( as < &ascsym[0] )
		{
		errcod = 8;
		goto synerr;
		}
/* get attributes and value of current symbol */
	att = (*as)->attrib;
	pvar = (struct varent *)*as;
	if( att & FUNC )
		{
		qclear( val.r );
		qclear( val.i );
		}
	else
		qcmov( (qcmplx *) pvar->value, &val );

/* Expression evaluation, continued. */

	switch( *ao )
		{
	case FUNC:
		pfun = (struct funent *)*as;
	/* Call the function with appropriate number of args */
	i = narstk[ narptr ];
	--narptr;
	switch(i)
			{
			case 0:
			( *(pfun->fun) )(&acc, &acc);
			break;
			case 1:
			( *(pfun->fun) )(&acc,&comstk[comptr-1],&acc);
			break;
			case 2:
			( *(pfun->fun) )(&acc, &comstk[comptr-2],
				&comstk[comptr-1],&acc);
			break;
			case 3:
			( *(pfun->fun) )(&acc, &comstk[comptr-3],
				&comstk[comptr-2], &comstk[comptr-1],&acc);
			break;
			case 4:
			( *(pfun->fun) )(&acc, &comstk[comptr-4],
					 &comstk[comptr-3],
				&comstk[comptr-2], &comstk[comptr-1],&acc);
			break;
			default:
			errcod = 16;
			goto synerr;
			}
		comptr -= i;
		accflg = 1;	/* in case at end of line */
		break;
	case EQU:
		if( ( att & TEMP) || ((att & VAR) == 0) || (att & STRING) )
			{
			errcod = 9;
			goto synerr;	/* can only assign to a variable */
			}
		pvar = (struct varent *)*as;
		qcmov( &acc, (qcmplx *) pvar->value );
		break;
	case PLUS:
		qcadd( &acc, &val, &acc );	break;
	case MINUS:
		qcsub( &acc, &val, &acc );	break;
	case MULT:
		qcmul( &acc, &val, &acc );	break;
	case DIV:
		if( (acc.r[1] == 0) && (acc.i[1] == 0) )
			{
/* divzer: */
			errcod = 10;
			goto synerr;
			}
		qcdiv( &acc, &val, &acc );	break;
/*
	case MOD:
		if( acc == 0 )
			goto divzer;
		acc = val % acc;	break;
	case LOR:
		acc |= val;		break;
	case LXOR:
		acc ^= val;		break;
	case LAND:
		acc &= val;		break;
*/
	case EOE:
		if( narptr < 0 )
			{
			errcod = 18;
			goto synerr;
			}
		narstk[narptr] += 1;
		qcmov( &acc, &comstk[comptr++] );
/*	printf( "\ncomptr: %d narptr: %d %d\n", comptr, narptr, acc );*/
		qcmov( &val, &acc );
		break;
		}


/*	expression evaluation, continued		*/

/* Pop evaluated tokens from scan list:		*/
	/* make temporary variable not busy	*/
	if( att & TEMP )
		(*as)->attrib &= ~BUSY;
	if( as < &ascsym[0] )	/* can this happen? */
		{
		errcod = 11;
		goto synerr;
		}
	--as;
nochg:
	--ao;
	--al;
	if( ao < &ascopr[0] )	/* can this happen? */
		{
		errcod = 12;
		goto synerr;
		}
/* If precedence level will now increase, then			*/
/* save accumulator in a temporary location			*/
	if( symval > *al )
		{
		/* find a free temp location */
		pvar = &temp[0];
		for( i=0; i<NTEMP; i++ )
			{
			if( (pvar->attrib & BUSY) == 0)
				goto temfnd;
			++pvar;
			}
		errcod = 17;
		printf( "no more temps\n" );
		pvar = &temp[0];
		goto synerr;

	temfnd:
		pvar->attrib |= BUSY;
		qcmov( &acc, (qcmplx *) pvar->value );
		/*printf( "temp %d\n", acc );*/
		accflg = 0;
		++as;	/* push the temp onto the scan list */
		*as = (struct symbol *)pvar;
		}
	}	/* End of evaluation loop */


/*	Push operator onto scan list when precedence increases	*/

pshopr:
	++ao;
	*ao = psym->attrib & 0xf;
	++al;
	*al = psym->sym + offset;
	goto gettok;
	}	/* end of OPR processing */


/* Token was not an operator.  Push symbol onto scan list.	*/
if( (lastok->attrib & OPR) == 0 )
	{
	errcod = 13;
	goto synerr;	/* quantities must be preceded by an operator */
	}
	/* ...but not by ) */
if( (unsigned short )lastok->attrib == (unsigned short)(OPR | RPAREN) )
	{
	errcod = 14;
	goto synerr;
	}
++as;
*as = psym;
goto gettok;

synerr:

#if INTHELP
printf( "%s ", intmsg[errcod] );
#endif
printf( " error %d\n", errcod );
if( savfil )
	fprintf( savfil, " error %d\n", errcod );
abmac();	/* flush the command line */
goto getcmd;
exit(0);
}	/* end of program */

/*						parser()	*/

/* Get token from input string and identify it.		*/


static char number[40] = {0};

struct symbol *parser( )
{
register struct symbol *psym;
register char *pline;
struct varent *pvar;
struct strent *pstr;
char *cp, *plc, *pn;
int i;
/* reference for old Whitesmiths compiler: */
/*
 *extern FILE *stdout;
 */

pline = interl;		/* get current location in command string	*/


/*	If at beginning of string, must ask for more input	*/
if( pline == line )
	{

	if( maccnt > 0 )
		{
		--maccnt;
		cp = maclin;
		plc = pline;
		while( (*plc++ = *cp++) != 0 )
			;
		goto mstart;
		}
	if( takptr < 0 )
		{	/* no take file active: prompt keyboard input */
		printf("* ");
		if( savfil )
			fprintf( savfil, "* " );
		}
/* 	Various ways of typing in a command line. */

/*
 * Old Whitesmiths call to print "*" immediately
 * use RT11 .GTLIN to get command string
 * from command file or terminal
 */

/*
 *	fflush(stdout);
 *	gtlin(line);
 */

 
#if USE_READLINE
	if (takptr < 0)
	  {
	    if (line_read)
	      {
		free (line_read);
		line_read = (char *)NULL;
	      }
	    /* Get a line from the user. */
	    line_read = readline ("");
	    /* If the line has any text in it, save it on the history. */
	    if (line_read && *line_read)
	      add_history (line_read);
	    /* Copy to local buffer. */
	    strcpy(line,line_read);
	  }
	else
#endif
	zgets( line, TRUE );	/* keyboard input for other systems: */


mstart:
	uposs = 1;	/* unary operators possible at start of line */
	}

ignore:
/* Skip over spaces */
while( *pline == ' ' )
	++pline;

/* unary minus after operator */
if( uposs && (*pline == '-') )
	{
	psym = &oprtbl[2];	/* UMINUS */
	++pline;
	goto pdon3;
	}
	/* COMP */
/*
if( uposs && (*pline == '~') )
	{
	psym = &oprtbl[3];
	++pline;
	goto pdon3;
	}
*/
if( uposs && (*pline == '+') )	/* ignore leading plus sign */
	{
	++pline;
	goto ignore;
	}

/* end of null terminated input */
if( (*pline == '\n') || (*pline == '\0') || (*pline == '\r') )
	{
	pline = line;
	goto endlin;
	}
if( *pline == ';' )
	{
	++pline;
endlin:
	psym = &oprtbl[1];	/* EOL */
	goto pdon2;
	}


/*						parser()	*/


/* Test for numeric input */
if( (ISDIGIT(*pline)) || (*pline == '.') )
	{
	nc = 0.0;	/* initialize numeric input to zero */
	qclear( qnc );
	pn = pline; /* save place */
	qclear( &qnc[NQ] );
	if( *pline == '0' )
		{ /* leading "0" may mean octal or hex radix */
		++pline;
		if( *pline == '.' )
			goto decimal; /* 0.ddd */
		/* leading "0x" means hexadecimal radix */
		if( (*pline == 'x') || (*pline == 'X') )
			{
			/* Copy the 0x.  */
			pline = pn;
			pn = number;
			*pn++ = *pline++;
			*pn++ = *pline++;
			while( ISXDIGIT(*pline) || *pline == '.' )
				*pn++ = *pline++;
			if (*pline == 'p' || *pline == 'P')
				{
				*pn++ = *pline++;
				if( (*pline == '-') || (*pline == '+') )
					*pn++ = *pline++;
				while( ISDIGIT(*pline) )
					*pn++ = *pline++;
				}
			else
				{
				printf( "0x requires p exponent field.\n" );
				pstr = &strtbl[NSTRNG-1];
				pstr->attrib |= ILLEG;
				psym = (struct symbol *)pstr;
				goto pdon0;
				}
			goto numcvt;
			}
		else
			{
			while( ISOCTAL( *pline ) )
				{
				i = (*pline++ & 0xff) - 060;
				nc = ( nc * 8.0) + i;
				etoq( (unsigned short *) &nc, qnc );
				}
			goto numdon;
			}
		}
	else
		{
		/* no leading "0" means decimal radix */
/******/
decimal:
		pn = number;
		while( (ISDIGIT(*pline)) || (*pline == '.') )
			*pn++ = *pline++;
/* get possible exponent field */
		if( (*pline == 'e') || (*pline == 'E') )
			*pn++ = *pline++;
		else
			goto numcvt;
		if( (*pline == '-') || (*pline == '+') )
			*pn++ = *pline++;
		while( ISDIGIT(*pline) )
			*pn++ = *pline++;
numcvt:
		*pn++ = ' ';
		*pn++ = 0;
		asctoq( number, qnc );
/*		sscanf( number, "%le", &nc );*/
		}
/* output the number	*/
numdon:
	/* search the symbol table of constants 	*/
	pvar = &contbl[0];
	for( i=0; i<NCONST; i++ )
		{
		if( (pvar->attrib & BUSY) == 0 )
			goto confnd;
		if( qcmp( pvar->value, qnc) == 0 )
			{
			psym = (struct symbol *)pvar;
			goto pdon2;
			}
		++pvar;
		}
	printf( "no room for constant\n" );
	psym = (struct symbol *)&contbl[0];
	goto pdon2;

confnd:
	pvar->spel= contbl[0].spel;
	pvar->attrib = CONST | BUSY;
	qmov( qnc, pvar->value );
	qclear( &(pvar->value[NQ]) );
	psym = (struct symbol *)pvar;
	goto pdon2;
	}

/* check for operators */
psym = &oprtbl[3];
for( i=0; i<NOPR; i++ )
	{
	if( *pline == *(psym->spel) )
		goto pdon1;
	++psym;
	}

/* if quoted, it is a string variable */
if( *pline == '"' )
	{
	/* find an empty slot for the string */
	pstr = strtbl;	/* string table	*/
	for( i=0; i<NSTRNG-1; i++ ) 
		{
		if( (pstr->attrib & BUSY) == 0 )
			goto fndstr;
		++pstr;
		}
	printf( "No room for string\n" );
	pstr->attrib |= ILLEG;
	psym = (struct symbol *)pstr;
	goto pdon0;

fndstr:
	pstr->attrib |= BUSY;
	plc = (char *)(pstr->string);
	++pline;
	for( i=0; i<39; i++ )
		{
		*plc++ = *pline;
		if( (*pline == '\n') || (*pline == '\0') || (*pline == '\r') )
			{
illstr:
			pstr = &strtbl[NSTRNG-1];
			pstr->attrib |= ILLEG;
			printf( "Missing string terminator\n" );
			psym = (struct symbol *)pstr;
			goto pdon0;
			}
		if( *pline++ == '"' )
			goto finstr;
		}

	goto illstr;	/* no terminator found */

finstr:
	*(--plc) = '\0';
	pstr->attrib |= BUSY;
	psym = (struct symbol *)pstr;
	goto pdon2;
	}
/* If none of the above, search function and symbol tables:	*/

/* copy character string to array lc[] */
plc = &lc[0];
while( ISALPHA(*pline) )
	{
	/* convert to lower case characters */
	if( ISUPPER( *pline ) )
		*pline += 040;
	*plc++ = *pline++;
	}
*plc = 0;	/* Null terminate the output string */

/*						parser()	*/

psym = (struct symbol *)menstk[menptr];	/* function table	*/
plc = &lc[0];
cp = psym->spel;
do
	{
	if( strcmp( plc, cp ) == 0 )
		goto pdon3;	/* following unary minus is possible */
	++psym;
	cp = psym->spel;
	}
while( *cp != '\0' );

psym = (struct symbol *)&indtbl[0];	/* indirect symbol table */
plc = &lc[0];
cp = psym->spel;
do
	{
	if( strcmp( plc, cp ) == 0 )
		goto pdon2;
	++psym;
	cp = psym->spel;
	}
while( *cp != '\0' );

pdon0:
pline = line;	/* scrub line if illegal symbol */
goto pdon2;

pdon1:
++pline;
if( (psym->attrib & 0xf) == RPAREN )
pdon2:	uposs = 0;
else
pdon3:	uposs = 1;

interl = pline;
return( psym );
}		/* end of parser */

/*	exit from current menu */

int cmdex()
{

if( menptr == 0 )
	{
	printf( "Main menu is active.\n" );
	}
else
	--menptr;

cmdh();
return(0);
}


/*			gets()		*/

int zgets( gline, echo )
char *gline;
int echo;
{
register char *pline;
register int i;


scrub:
pline = gline;
getsl:
	if( (pline - gline) >= LINLEN )
		{
		printf( "\nLine too long\n *" );
		goto scrub;
		}
	if( takptr < 0 )
		{	/* get character from keyboard */
#if DECPDP
		gtlin( gline );
		return(0);
#else
		*pline = getchar();
#endif
		}
	else
		{	/* get a character from take file */
		i = fgetc( takstk[takptr] );
		if( i == -1 )
			{	/* end of take file */
			if( takptr >= 0 )
				{	/* close file and bump take stack */
				fclose( takstk[takptr] );
				takptr -= 1;
				}
			if( takptr < 0 )	/* no more take files:   */
				printf( "*" ); /* prompt keyboard input */
			goto scrub;	/* start a new input line */
			}
		*pline = i;
		}

	*pline &= 0x7f;
	/* xon or xoff characters need filtering out. */
	if ( *pline == XON || *pline == XOFF )
		goto getsl;

	/*	control U or control C	*/
	if( (*pline == 025) || (*pline == 03) )
		{
		printf( "\n" );
		goto scrub;
		}

	/*  Backspace or rubout */
	if( (*pline == 010) || (*pline == 0177) )
		{
		pline -= 1;
		if( pline >= gline )
			{
			if ( echo )
				printf( "\010\040\010" );
			goto getsl;
			}
		else
			goto scrub;
		}
	if ( echo )
		printf( "%c", *pline );
	if( (*pline != '\n') && (*pline != '\r') )
		{
		++pline;
		goto getsl;
		}
	*pline = 0;
	if ( echo )
		printf( "%c", '\n' );	/* \r already echoed */
	if( savfil )
		fprintf( savfil, "%s\n", gline );
return 0;
}


/*		help function  */
int cmdhlp()
{

printf( "%s", idterp );
printf( "\nFunctions:\n" );
prhlst( (struct symbol *) &funtbl[0] );
printf( "\nVariables:\n" );
prhlst( (struct symbol *) &indtbl[0] );
printf( "\nOperators:\n" );
prhlst( (struct symbol *) &oprtbl[2] );
printf("\n");
return(0);
}


int cmdh()
{

prhlst( (struct symbol *) menstk[menptr] );
printf( "\n" );
return(0);
}

/* print keyword spellings */

int prhlst(ps)
struct symbol *ps;
{
register int j, k;
/* int m; */

j = 0;
while( *(ps->spel) != '\0' )
	{
#if NOTABS
	k = strlen( ps->spel ) + 1;
	j += k;
	if( j > 72 )
		{
		printf( "\n" );
		j = k;
		}
	if( ps->attrib & MENU )
		{
		printf( "%s/ ", ps->spel );
		j += 1;
		}
	else
		{
		printf( "%s ", ps->spel );
		}
#else
	k = strlen( ps->spel )  -  1;
/* size of a tab field is 2**3 chars */
	m = ((k >> 3) + 1) << 3;
	j += m;
	if( j > 72 )
		{
		printf( "\n" );
		j = m;
		}
	if( ps->attrib & MENU )
		printf( "%s/\t", ps->spel );
	else
		printf( "%s\t", ps->spel );
#endif
	++ps;
	}
return 0;
}


#if SALONE
int init(){return 0;}
#endif


/*	macro commands */

/*	define macro */
int cmddm(arg)
QELT *arg;
{

zgets( maclin, TRUE );
return(0);
}

/*	type (i.e., display) macro */
int cmdtm(arg)
QELT *arg;
{

printf( "%s\n", maclin );
return(0);
}

/*	execute macro # times */
int cmdem( arg )
QELT *arg;
{
double dn;
int n;

qtoe( arg, (unsigned short *) &dn );
n = dn;
if( n <= 0 )
	n = 1;
maccnt = n;
return( n );
}


/* open a take file */

int take( fname )
char *fname;
{
FILE *f;

while( *fname == ' ' )
	fname += 1;
f = fopen( fname, "r" );

if( f == 0 )
	{
	printf( "Can't open take file %s\n", fname );
	takptr = -1;	/* terminate all take file input */
	return(-1);
	}
takptr += 1;
takstk[ takptr ]  =  f;
printf( "Running %s\n", fname );
return(0);
}


/*	abort macro execution */
int abmac()
{

maccnt = 0;
interl = line;
return 0;
}


/* display integer part in hex, octal, and decimal
 */

int hex(qx)
QELT *qx;
{
long z;
double x;
double fabs();

qtoe( qx, (unsigned short *) &x );
if( fabs(x) >= 2.147483648e9 )
	{
	printf( "hex: too large\n" );
	return(x);
	}

z = x;
printf( "0%lo  0x%lx  %ld.\n", z, z, z );
return 0;
}

int rbits();

int bits( x )
qcmplx *x;
{

rbits( x->r );
rbits( x->i );
return 0;
}


int rbits( u )
QELT u[];
{
QELT x[NQ+1];
unsigned short e113[8];
unsigned short d[4];
int i, j;

qmovz( u, x );
j = 0;
for( i=0; i<NQ; i++ )
	{
#if WORDSIZE == 32
	printf( "0x%08x,", x[i] );
	if( savfil )
		fprintf( savfil, "0x%08x,", x[i] );
	if( ++j > 6 )
#else
	printf( "0x%04x,", x[i] & 0xffff );
	if( savfil )
		fprintf( savfil, "0x%04x,", x[i] & 0xffff );
	if( ++j > 7 )
#endif
		{
		j = 0;
		printf( "\n" );
		if( savfil )
		  fprintf( savfil, "\n" );
		}
	}
printf( "\n" );
if( savfil )
	  fprintf( savfil, "\n" );
qtoe( u, d );
printf( "qtoe: " );
if( savfil )
	fprintf( savfil, "qtoe: " );
#if WORDSIZE == 32
for( i=0; i<2; i++ )
  {
	printf( "%08x ", d[i] );
	if( savfil )
		fprintf( savfil, "0x%08x,", d[i] );
  }
#else
for( i=0; i<4; i++ )
  {
	printf( "%04x ", d[i] & 0xffff );
	if( savfil )
		fprintf( savfil, "0x%04x,", d[i] & 0xffff );
  }
#endif
printf( "\n printf: %.16e\n", *(double *)d );
if( savfil )
	fprintf( savfil, "\n printf: %.16e\n", *(double *)d );

etoq( d, x );
printf( "etoq:\n" );
if( savfil )
	fprintf( savfil, "etoq:\n" );
j = 0;
for( i=0; i<NQ; i++ )
	{
#if WORDSIZE == 32
	printf( "0x%08x,", x[i] );
	if( savfil )
		fprintf( savfil, "0x%08x,", x[i] );
	if( ++j > 6 )
#else
	printf( "0x%04x,", x[i] & 0xffff );
	if( savfil )
		fprintf( savfil, "0x%04x,", x[i] & 0xffff );
	if( ++j > 7 )
#endif
		{
		j = 0;
		printf( "\n" );
		if( savfil )
		  fprintf( savfil, "\n" );
		}
	}
printf( "\n" );

qtoe113( u, e113 );
printf( "qtoe113: " );
if( savfil )
	fprintf( savfil, "qtoe113: " );
for( i=0; i<8; i++ )
  {
	printf( "%04x,", e113[i] & 0xffff );
	if( savfil )
		fprintf( savfil, "%04x,", e113[i] & 0xffff );
  }
e113toq( e113, x );
printf( "\ne113toq:\n" );
if( savfil )
	fprintf( savfil, "\ne113toq:\n" );
j = 0;
for( i=0; i<NQ; i++ )
	{
#if WORDSIZE == 32
	printf( "0x%08x,", x[i] );
	if( savfil )
		fprintf( savfil, "0x%08x,", x[i] );
	if( ++j > 6 )
#else
	printf( "0x%04x,", x[i] & 0xffff );
	if( savfil )
		fprintf( savfil, "0x%04x,", x[i] & 0xffff );
	if( ++j > 7 )
#endif
		{
		j = 0;
		printf( "\n" );
		if( savfil )
		  fprintf( savfil, "\n" );
		}
	}
printf( "\n" );
if( savfil )
	  fprintf( savfil, "\n" );
return(0);
}


/* Exit to monitor. */
int mxit()
{

if( savfil )
	fclose( savfil );
exit(0);
}


int cmddig( x )
QELT x[];
{
double dx;

qtoe( x, (unsigned short *) &dx );
ndigits = dx;
if( ndigits <= 0 )
	ndigits = DEFDIS;
return(0);
}



int qsave(x)
char *x;
{

if( savfil )
	fclose( savfil );
while( *x == ' ' )
	x += 1;
if( *x == '\0' )
	savnam = "calc.sav";
else
	savnam	= x;
savfil = fopen( savnam, "w" );
if( savfil == 0 )
	printf( "Error opening %s\n", savnam );
return 0;
}



int qsys(x)
char *x;
{

system( x+1 );
cmdh();
return 0;
}

/* z = x + iy ) */
int makec( x, y, z )
qcmplx *x, *y, *z;
{
if (x == (qcmplx *)0 || y == (qcmplx *)0 || z == (qcmplx *)0)
  {
    printf("Usage: c = cmplx(r1,r2) creates r1 + r2 i\n");
    return 0;
  }
qmov( x->r, z->r );
qmov( y->r, z->i );
return 0;
}

/* Re(z) */
int getre( x, y )
qcmplx *x, *y;
{
qmov( x->r, y->r );
qclear( y->i );
return 0;
}

/* Im(z) */
int getim( x, y )
qcmplx *x, *y;
{
qmov( x->i, y->r );
qclear( y->i );
return 0;
}

/* Arg z */
int qcarg( x, y )
qcmplx *x, *y;
{
#if ANSIC
qatn2( x->i, x->r, y->r );
#else
qatn2( x->r, x->i, y->r );
#endif
qclear( y->i );
return 0;
}

/* |z| to complex data structure */
int qccabs( x, y )
qcmplx *x, *y;
{
qcabs( x, y->r );
qclear( y->i );
return 0;
}

int qcconj( x, y )
qcmplx *x, *y;
{
qcmov( x, y );
if( y->i[1] != 0 )
	qneg( y->i );
return 0;
}

static int
zpolylog( qn, x, y )
qcmplx *qn, *x, *y;
{
  long n;
  QELT t[NQ];

  qifrac (qn->r, &n, t);
  qcpolylog ((int) n, x, y);
  return 0;
}
