/* testmain.c, (c) 1996, Samuel Ferguson.  This file contains routines
for developing and testing interval arithmetic software.  */

#include "system_headers.h"
#include "sphere.h"
#include "interval.h"
#include "i_sphere.h"
#include "i_bounds.h"
#include "second_partials.h"
#include "macros.h"
#include "i_taylor.h"

/* #pragma fenv_access	*/
#define GET_INTERVAL			0
#define CONSTANTS					0
#define DELTA_PARTIALS		0
#define A_PARTIALS				0
#define SOL_PARTIALS			0
#define BVOL_PARTIALS			0
#define DIH_PARTIALS			0
#define GMA_PARTIALS			0
#define TOMSP_PARTIALS		0
#define VORVOL_PARTIALS		0
#define VOR_PARTIALS			0
#define OCTA_VOR_PARTIALS	0
#define ROUGH							0
#define MINMAX_A					0
#define GMA								0
#define VOR								0
#define DIH								0
#define MAXMIN_U					0
#define QR_CRAD						0
#define CRAD							0
#define SCRAD							0
#define CRAD3LEN					0

#define OUTPUT( var )		num2dec(&form, var, &dec);			\
												dec2str(&form, &dec, buf );			\
												puts( buf );


extern double i_pi_const[2];
extern double i_pi_2_const[2];
extern double i_doct_const[2];
extern double i_two_pi_5_const[2];

/* External prototypes */

void i_init( void );
static double daan_atan (double x);
double sam_atan (double x);


/* Prototypes */

void get_constants( void );
void query_rounding( void );
void set_rounding( void );
void test0( void );
void checksqrt( void );
void checkatan( void );
void getbits( double x, char c[64] );
void printbits( double x );
void printintbits( int n );
void get_trash( void );
void get_more_trash( void );
void get_second_trash( void );
void get_taylor_trash( void );
void get_comp_trash( void );
void get_new_trash( void );

/* Routines */

/* Determine useful constants . . . to avoid division,
etc. */
void get_constants( void )
{
/* This is pretty annoying.  Metrowerks is propagating
constants in a way which isn't completely general:  e.g.,
through changes in rounding direction.  Silly silly.  Goes
away if I turn off optimization, though. */
	/*					  					123456789012345678901234567890  */
	char pi_str[255] =   "3.14159265358979323846264338328";
	char doct_str[255] = "0.7209029495174650928";
	char two_pi_5_str[255] = "1.2566370614359172953850574";
	char sphcut_str[255] = "0.679673818908243874192785026784";
	/*	decimal dec;	*/ /* Can't figure out how to use str2dec, etc.  No docs. */
	double temp[2], x, y;
	
	char buf[256];
	decimal dec;
	decform form;
	
	/* form.style  = FLOATDECIMAL; */
	/* form.style = FIXEDDECIMAL; */
	form.style = FLOATDECIMAL;
	form.digits = 21;

	
	/*	str2dec( pi_str, */
	/* pi */
	ROUND_DOWN;
	temp[0] = atof( pi_str );
	ROUND_UP;
	temp[1] = atof( pi_str );
	ROUND_NEAR;
	printf("Pi = \n[%.20g, \n %.20g]\n\n", temp[0], temp[1]);
	OUTPUT( temp[0] );
	OUTPUT( temp[1] );
	
	/* pi/2 */
	ROUND_DOWN;
	temp[0] = temp[0]*0.5;
	ROUND_UP;
	temp[1] = temp[1]*0.5;
	ROUND_NEAR;
	printf("Pi/2 = \n[%.20g, \n %.20g]\n\n", temp[0], temp[1]);
	OUTPUT( temp[0] );
	OUTPUT( temp[1] );
	
	/* delta_oct */
	ROUND_DOWN;
	temp[0] = atof( doct_str );
	ROUND_UP;
	temp[1] = atof( doct_str );
	ROUND_NEAR;
	printf("delta_oct = \n[%.20g, \n %.20g]\n\n", temp[0], temp[1]);
	OUTPUT( temp[0] );
	OUTPUT( temp[1] );
	
	/* 2 Pi/5 */
	ROUND_DOWN;
	temp[0] = atof( two_pi_5_str );
	ROUND_UP;
	temp[1] = atof( two_pi_5_str );
	ROUND_NEAR;
	printf("2 Pi/5 = \n[%.20g, \n %.20g]\n\n", temp[0], temp[1]);
	OUTPUT( temp[0] );
	OUTPUT( temp[1] );

	/* 2 atan( 0.25 sqrt(2) ) */
	ROUND_DOWN;
	temp[0] = atof( sphcut_str );
	ROUND_UP;
	temp[1] = atof( sphcut_str );
	ROUND_NEAR;
	printf("2 atan( 0.25 sqrt(2) ) = \n[%.20g, \n %.20g]\n\n", temp[0], temp[1]);
	OUTPUT( temp[0] );
	OUTPUT( temp[1] );
	
	/* sqrt( 2 ) */
	ROUND_DOWN;
	temp[0] = sqrt( 2.0 );
	ROUND_UP;
	temp[1] = sqrt( 2.0 );
	ROUND_NEAR;
	printf("Sqrt[2] = \n[%.20g, \n %.20g]\n\n", temp[0], temp[1]);
	OUTPUT( temp[0] );
	OUTPUT( temp[1] );
	/*
	printbits( temp[0] );
	printbits( temp[1] );
	*/
	
	/* 2 sqrt( 2 ) */
	ROUND_DOWN;
	temp[0] = 2.0*sqrt( 2.0 );
	ROUND_UP;
	temp[1] = 2.0*sqrt( 2.0 );
	ROUND_NEAR;
	printf("2 Sqrt[2] = \n[%.20g, \n %.20g]\n\n", temp[0], temp[1]);
	OUTPUT( temp[0] );
	OUTPUT( temp[1] );
	
	/* 1/3 */
	ROUND_DOWN;
	x = 1.0;
	y = 3.0;
	ROUND_DOWN;
	temp[0] = x/y;
	ROUND_UP;
	temp[1] = x/y;
	ROUND_NEAR;
	printf("1/3 = \n[%.20g, \n %.20g]\n\n", temp[0], temp[1]);
	OUTPUT( temp[0] );
	OUTPUT( temp[1] );
	/*
	printbits( temp[0] );
	printbits( temp[1] );
	*/
	
	/*
	printf("1/3 = \n[%.25g, \n %.25g]\n\n", temp[0], temp[1]);
	printf("1/3 = \n[%.17f, \n %.17f]\n\n", temp[0], temp[1]);
	printf("1/3 = \n[%.25f, \n %.25f]\n\n", temp[0], temp[1]);
	printf("1/3 = \n[%.9g, \n %.9g]\n\n", temp[0], temp[1]);
	*/
	/*
	num2dec(&form, temp[0], &dec);
	dec2str(&form, &dec, buf );
	printf("result: \n");
	puts( buf );
	printf("\n");

	num2dec(&form, temp[1], &dec);
	dec2str(&form, &dec, buf );
	printf("result: \n");
	puts( buf );
	printf("\n");
	*/
	/*
	ROUND_NEAR;
	return;
	*/
	
	/* 2/3 */
	ROUND_DOWN;
	x = 2.0;
	y = 3.0;
	ROUND_DOWN;
	temp[0] = x/y;
	ROUND_UP;
	temp[1] = x/y;
	ROUND_NEAR;
	printf("2/3 = \n[%.20g, \n %.20g]\n\n", temp[0], temp[1]);
	OUTPUT( temp[0] );
	OUTPUT( temp[1] );
	/*
	printbits( temp[0] );
	printbits( temp[1] );
	*/
	
	/* 1/6 */
	ROUND_DOWN;
	x = 1.0;
	y = 6.0;
	ROUND_DOWN;
	temp[0] = x/y;
	ROUND_UP;
	temp[1] = x/y;
	ROUND_NEAR;
	printf("1/6 = \n[%.20g, \n %.20g]\n\n", temp[0], temp[1]);
	OUTPUT( temp[0] );
	OUTPUT( temp[1] );
	/*
	printbits( temp[0] );
	printbits( temp[1] );
	*/
	
	/* 1/12 */
	ROUND_DOWN;
	x = 1.0;
	y = 12.0;
	ROUND_DOWN;
	temp[0] = x/y;
	ROUND_UP;
	temp[1] = x/y;
	ROUND_NEAR;
	printf("1/12 = \n[%.20g, \n %.20g]\n\n", temp[0], temp[1]);
	OUTPUT( temp[0] );
	OUTPUT( temp[1] );
	/*
	printbits( temp[0] );
	printbits( temp[1] );
	*/
	
	/* 1/36 */
	ROUND_DOWN;
	x = 1.0;
	y = 36.0;
	ROUND_DOWN;
	temp[0] = x/y;
	ROUND_UP;
	temp[1] = x/y;
	ROUND_NEAR;
	printf("1/36 = \n[%.20g, \n %.20g]\n\n", temp[0], temp[1]);
	OUTPUT( temp[0] );
	OUTPUT( temp[1] );
	/*
	printbits( temp[0] );
	printbits( temp[1] );
	*/

	/* 1/48 */
	ROUND_DOWN;
	x = 1.0;
	y = 48.0;
	ROUND_DOWN;
	temp[0] = x/y;
	ROUND_UP;
	temp[1] = x/y;
	ROUND_NEAR;
	printf("1/48 = \n[%.20g, \n %.20g]\n\n", temp[0], temp[1]);
	OUTPUT( temp[0] );
	OUTPUT( temp[1] );
	/*
	printbits( temp[0] );
	printbits( temp[1] );
	*/

	/* 1/sqrt(2) */
	ROUND_DOWN;
	x = 2.0;
	ROUND_DOWN;
	y = sqrt( x );
	ROUND_UP;
	temp[1] = 1/y;
	y = sqrt( x );
	ROUND_DOWN;
	temp[0] = 1/y;
	ROUND_NEAR;
	printf("1/sqrt(2) = \n[%.20g, \n %.20g]\n\n", 
		temp[0], temp[1]);
	OUTPUT( temp[0] );
	OUTPUT( temp[1] );
	/*
	printbits( temp[0] );
	printbits( temp[1] );
	*/

	/* 1/(2*sqrt(2)) */
	ROUND_DOWN;
	x = 2.0;
	ROUND_DOWN;
	y = sqrt( x );
	ROUND_UP;
	temp[1] = 0.5/y;
	y = sqrt( x );
	ROUND_DOWN;
	temp[0] = 0.5/y;
	ROUND_NEAR;
	printf("1/(2*sqrt(2)) = \n[%.20g, \n %.20g]\n\n", 
		temp[0], temp[1]);
	OUTPUT( temp[0] );
	OUTPUT( temp[1] );
	/*
	printbits( temp[0] );
	printbits( temp[1] );
	*/

	/* 2.51 */
	ROUND_DOWN;
	x = 251.0;
	ROUND_DOWN;
	y = 100.0;
	ROUND_UP;
	temp[1] = x/y;
	ROUND_DOWN;
	temp[0] = x/y;
	ROUND_NEAR;
	printf("2.51 = \n[%.20g, \n %.20g]\n\n", 
		temp[0], temp[1]);
	OUTPUT( temp[0] );
	OUTPUT( temp[1] );
	/*
	printbits( temp[0] );
	printbits( temp[1] );
	*/
	
	/* 4/2.51 */
	ROUND_DOWN;
	y = 251.0;
	ROUND_DOWN;
	x = 400.0;
	ROUND_UP;
	temp[1] = x/y;
	ROUND_DOWN;
	temp[0] = x/y;
	ROUND_NEAR;
	printf("4/2.51 = \n[%.20g, \n %.20g]\n\n", 
		temp[0], temp[1]);
	OUTPUT( temp[0] );
	OUTPUT( temp[1] );
	/*
	printbits( temp[0] );
	printbits( temp[1] );
	*/

	/* 2^(-49) */
	ROUND_NEAR;
	x = ldexp(1.0,-49);
	printf("2^-49 = \n%.20g\n\n", 
		x);
	OUTPUT( x );
	printbits( ATANERR );


	/*
	ROUND_DOWN;
	ROUND_DOWN;
	temp[0] = 2.0*atan( 0.25*sqrt(2.0) );
	ROUND_UP;
	temp[1] = 2.0*atan( 0.25*sqrt(2.0) );
	ROUND_NEAR;
	printf("sphcut = \n[%.20g, \n %.20g]\n\n", temp[0], temp[1]);
	OUTPUT( temp[0] );
	OUTPUT( temp[1] );
	printbits( temp[0] );
	printbits( temp[1] );
	*/
/*							123456789012345678901234567890 */
/* From Mathematica,
	-1.04 pt = -0.0575885914952024245192767955091\
	80892114057249225535003474089391396344229756\
45081820025893370263472960366628049101874842327\
4129683669400747065065701663378\
74794313620470512490692867959089\
950674902225358794018948809,
	-0.52 pt = -0.0287942957476012122596383977545\
	90446057028624612767501737044695698172114878\
22540910012946685131736480183314024550937421163\
7064841834700373532532850831689\
373971568102352562453464339795449\
753374511126793970094744044
*/
	ROUND_NEAR;
}


void query_rounding( void )
{
#if GENERATINGPOWERPC
	int query;
	
	query = fegetround();
	switch( query ) {
		case FE_TONEAREST:
			printf("Rounding to nearest.\n");
			break;
		case FE_TOWARDZERO:
			printf("Rounding toward zero.\n");
			break;
		case FE_UPWARD:
			printf("Rounding upward.\n");
			break;
		case FE_DOWNWARD:
			printf("Rounding downward.\n");
			break;
		default:
			printf("No rounding mode recognized.\n");
			break;
		}
#endif
}


void set_rounding( void )
{
	int mode;
	
	printf("Rounding modes are 0: nearest, 1: toward zero,");
	printf("  2:  upward, 3:  downward\n");
	printf("Enter new rounding mode:  ");
	scanf("%d", &mode);
	printf("Got %d\n", mode);
	switch( mode ) {
		case 0:
			ROUND_NEAR;
			break;
		case 1:
			ROUND_ZERO;
			break;
		case 2:
			ROUND_UP;
			break;
		case 3:
			ROUND_DOWN;
			break;
		default:
			printf("Sorry, that's an unsupported mode.\n");
			break;
		}
}

void test0( void )
{
/*
	double x;
	char str[255];
*/
	
	printf("Hello.  Welcome to IEEE testing.\n\n");
	query_rounding();
	while( 1 ) {
		set_rounding();
		query_rounding();
		printf("1234567890123456789012456789012345678901234567890\n");
		/*
		printf("Enter a number:   ");
		scanf("%s", str);
		sscanf(str, "%lf", &x);
		printf("You gave me:      %s\n", str);
		printf("Interpreted as:   %.22g\n", x);
		printf("\n");
		*/
		checksqrt();
		}
}


void main( void )
{
/*
	double x;
	char str[255];
*/
	int n;
	time_t tstop;
	char *charptr;
	
	printf("Hello.  Welcome to IEEE testing.\n\n");
	time(&tstop);
	charptr = ctime( &tstop );
	printf("Time:    \t");
	puts( charptr );

	/*
	while( 1 ) {
		printf("Enter an integer:  ");
		scanf("%d", &n);
		printf("Got %d\n", n);
		printf("123456789012345678901234567890123456789012345678901234567890\n");
		printintbits( n );
	}
	*/
	n = 1;
	while( n ) {
		/*
		printf("Enter a floating-point (double) value:  ");
		scanf("%lf", &x);
		printbits( x );
		*/
		/*	set_rounding();	*/
		/*	checksqrt();	*/
		/*	checkatan();	*/
#if DA_SYSTEM == 1
/*
		printf("Bitwise OR of floating-point exception macros: %#x\n", 
			fetestexcept( FE_ALL_EXCEPT ) );
		printf("Rounding direction macros: %#x\n\n", 
			fegetround() );
*/
#elif DA_SYTEM == 3
		printf("Bitwise OR of floating-point exception macros: %#x\n", 
			fpgetmask() );
		printf("Rounding direction macros: %#x\n\n", 
			fpgetround() );
#endif
		get_new_trash();
		n = 0;
		/*
		t_print_const();
		get_taylor_trash();
		get_second_trash();
		get_more_trash();
		get_comp_trash();
		get_constants();
		get_new_trash();
		n = 0;
		*/
	}
}


void checksqrt( void )
{
	double x, y, diff, y2;
	
	printf("Enter x:  ");
	scanf("%lf", &x);
	y = sqrt( x );
	ROUND_NEAR;	/* testing. */
	y2 = y*y;
	printf("sqrt( %.21g ) = %.21g\n", x, y);
	printf("y2 = %.21g\n", y2);
	diff = y*y - x;
	printf("sqrt(x)^2 - x = %.21g\n", diff);
	printf("(no MAC) sqrt(x)^2 - x = %.21g\n", y2 - x);
}


void checkatan( void )
{
	double x, y, diff, y2;
	
	printf("Enter x:  ");
	scanf("%lf", &x);
	ROUND_NEAR;	/* testing. */
	y = atan( x );
	printf("(to nearest) %.21g\n", y);
	ROUND_UP;	/* testing. */
	y = atan( x );
	printf("(up)         %.21g\n", y);
	ROUND_DOWN;	/* testing. */
	y = atan( x );
	printf("(down)       %.21g\n", y);
	ROUND_NEAR;	/* testing. */
}


void getbits( double x, char c[64] )
{
	int i, j, m, n[2], sh;
	int *ptr;
	double *dptr;
	
	dptr = &x;
	ptr = (int *) dptr;
	for( j=0; j<2; j++ ) {
		n[j] = ptr[j];
		}
	/*
	for( j=0; j<2; j++ )
		printintbits( n[j] );
	*/
	for( j=0; j<2; j++ ) {
		m = 1;
		sh = n[j];
		for( i=0; i<32; i++ ) {
			c[32*(j+1)-1-i] = sh & m;
			sh = sh >> 1;
			}
		}
}


void printbits( double x )
{
	char c[64];
	int i;
	
	getbits( x, c );
	printf("Printing binary expansion of %.21g\n", x);
	for( i=0; i<12; i++ )
		printf("%d", c[i] );
	printf("  (sign bit and exponent (biased))\n");
	for( i=0; i<6; i++ )
		printf("1234567890");
	printf("\n");
	for( i=12; i<64; i++ )
		printf("%d", c[i] );
	printf("  (mantissa)\n");
}


void printintbits( int n )
{
	char c[32];
	int i, m;
	unsigned int un;
	
	un = n;
	m = 1;
	for( i=0; i<32; i++ ) {
		c[31-i] = un & m;
		un = un >> 1;
		}
	for( i=0; i<32; i++ )
		printf("%d", c[i] );
	printf("\n");
}


void get_trash( void )
{
	double y[12], x[12], out[2], sy[7], dout, diff, eps;
	double delta[2], sqrtdelta[2], oosqrtdelta[2];
	double temp[2];
	int i;
	
	i_init();
	
	ROUND_NEAR;
	
	printf("Enter edge lengths:  ");
	for( i=1; i<7; i++ )
		scanf("%lf", sy + i);
	printf("Got the following:  \n");
	for( i=1; i<7; i++ )
		printf("%f\t", sy[i]);
	printf("\n\n");
	printf("Enter epsilon (interval width):  ");
	scanf("%lf", &eps);
	
	ROUND_UP;
	for( i=0; i<6; i++ ) {
		y[2*i] = sy[i+1];
		y[2*i+1] = sy[i+1] + eps;
		}
	ROUND_NEAR;

	for( i=0; i<12; i+=2 )
		printf("%f\t", y[i]);
	printf("\n");
	for( i=1; i<12; i+=2 )
		printf("%f\t", y[i]);
	printf("\n\n");
	
	i_bigdelta( x, delta );
	i_sqrt( delta, sqrtdelta );
	temp[0] = 1.0;
	temp[1] = 1.0;
	i_div( temp, sqrtdelta, oosqrtdelta );
	
	i_bvol(  y, oosqrtdelta, out );
	printf("i_bvol gives (%.21g, %.21g)\n", out[0], out[1]);
	printf("out[1]-out[0] = %.21g\n", out[1]-out[0]);
	dout = bvol( sy );
	printf("bvol gives %.21g\n", dout);
	diff = fabs(2.0*dout - out[0] - out[1]);
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
	
	i_afunc(  y,  out );
	printf("i_afunc gives (%.21g, %.21g)\n", out[0], out[1]);
	printf("out[1]-out[0] = %.21g\n", out[1]-out[0]);
	dout = afunc( sy );
	printf("afunc gives %.21g\n", dout);
	diff = fabs(2.0*dout - out[0] - out[1]);
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
	
	/*
	i_solid(  y, sqrtdelta, out );
	printf("i_solid gives (%.21g, %.21g)\n", out[0], out[1]);
	printf("out[1]-out[0] = %.21g\n", out[1]-out[0]);
	*/
	dout = solid( sy );
	printf("solid gives %.21g\n", dout);
	diff = fabs(2.0*dout - out[0] - out[1]);
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
	
	i_tvol(  y,  out );
	printf("i_tvol gives (%.21g, %.21g)\n", out[0], out[1]);
	printf("out[1]-out[0] = %.21g\n", out[1]-out[0]);
	/*
	dout = tvol( sy );
	printf("tvol gives %.21g\n", dout);
	*/
	dout = tetvolume( sy );
	printf("tetvolume gives %.21g\n", dout);
	diff = fabs(2.0*dout - out[0] - out[1]);
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
	
	i_gma(  y, sqrtdelta, out );
	printf("i_gma gives (%.21g, %.21g)\n", out[0], out[1]);
	printf("out[1]-out[0] = %.21g\n", out[1]-out[0]);
	dout = gma( sy );
	printf("gma gives %.21g\n", dout);
	diff = fabs(2.0*dout - out[0] - out[1]);
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
	
	i_dih(  x,  out );
	printf("i_dih gives (%.21g, %.21g)\n", out[0], out[1]);
	printf("out[1]-out[0] = %.21g\n", out[1]-out[0]);
	dout = dih( sy );
	printf("dih gives %.21g\n", dout);
	diff = fabs(2.0*dout - out[0] - out[1]);
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
	
	/*
	i_dbvol(  y,  out );
	i_bigdelta(  x,  out );
	i_istet(  y );
	i_tomsrho(  x,  out );
	*/
	i_crad(  y,  out );
	printf("i_crad gives (%.21g, %.21g)\n", out[0], out[1]);
	printf("out[1]-out[0] = %.21g\n", out[1]-out[0]);
	dout = circumradius( sy );
	printf("crad gives %.21g\n", dout);
	diff = fabs(2.0*dout - out[0] - out[1]);
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
	/*
	i_tomsu(  x,  out );
	i_tomsv(  x,  out );
	i_auxPfun(  x,  out );
	i_tomsP(  x,  out );
	i_voronoivol(  y,  out );
	*/
	/*
	i_vor(  y, x, sqrtdelta, out );
	printf("i_vor gives (%.21g, %.21g)\n", out[0], out[1]);
	printf("out[1]-out[0] = %.21g\n", out[1]-out[0]);
	dout = vor( sy );
	printf("vor gives %.21g\n", dout);
	diff = fabs(2.0*dout - out[0] - out[1]);
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
	*/
	i_crad3len(  y, out );
	printf("i_crad3len gives (%.21g, %.21g)\n", out[0], out[1]);
	printf("out[1]-out[0] = %.21g\n", out[1]-out[0]);
	dout = crad3len( sy );
	printf("crad3len gives %.21g\n", dout);
	diff = fabs(2.0*dout - out[0] - out[1]);
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
	
	ROUND_NEAR;
}


void get_more_trash( void )
{
	double y[12], out[2], sy[7], dout, diff, eps, temp, temp2;
	double x[12], xp[12], delta[2], sx[7], delta_part[12];
	double part[12], sum, data[2], sqrtdelta[2], delta32[2];
	double sol_part[12], sol_xpart[12], vor_xpart[12];
	double a_part[12], bvol_part[12], bvol_xpart[12];
	double ck_part[12];
	double dih_part[12], gma_part[12], gma_xpart[12];
	double tomsP_part[12];
	double rel_part[12], vorvol_part[12], vor_part[12];
	double slope, correps, max, out2[2], out3[2];
	double ooy[12], oosqrtdelta[2];
	double ab[4], xy[4], al[2], al_ab[4];
	double vpart[12], spart[12];
	int i, j, delta_partials[6], count, c126, c135, finfo[2];
	int interval;
	

	i_init();
	
	ROUND_NEAR;
	eps = 0.0;
	
	printf("Enter data type for cells:  1 for intervals, 0 otherwise:  ");
	scanf("%d", &interval);
	
	if( interval == 0 ) {
		printf("Enter edge lengths:\n");
		
		for( i=1; i<7; i++ )
			scanf("%lf", sy + i);
		printf("Got the following:  \n");
		for( i=1; i<7; i++ )
			printf("%f\t", sy[i]);
		printf("\n\n");
		printf("Enter epsilon (interval width):  ");
		scanf("%lf", &eps);
		
		ROUND_UP;
		for( i=0; i<6; i++ ) {
			y[2*i] = sy[i+1];
			y[2*i+1] = sy[i+1] + eps;
			}
		ROUND_NEAR;
		
		for( i=1; i<7; i++ )
			sx[i] = sy[i]*sy[i];
	} else {
	printf("Enter edge lengths (as intervals):\n");
	for( i=0; i<12; i++ )
		scanf("%lf", y + i );
	for( i=0; i<6; i++ )
		sy[i+1] = y[2*i];
	}

/*
	for( i=0; i<12; i+=2 )
		printf("%f\t", y[i]);
	printf("\n");
	for( i=1; i<12; i+=2 )
		printf("%f\t", y[i]);
	printf("\n\n");
*/	

	/* Compute square of edge lengths */
	ROUND_UP;
	for( i=1; i<12; i+=2 )
		x[i] = y[i]*y[i];
	ROUND_DOWN;
	for( i=0; i<12; i+=2 )
		x[i] = y[i]*y[i];

	/* Compute ooy = 1/y */
	ROUND_UP;
	for( i=1; i<12; i+=2 )
		ooy[i] = 1.0/y[i-1];
	ROUND_DOWN;
	for( i=0; i<12; i+=2 )
		ooy[i] = 1.0/y[i+1];
	
/*	
	for( i=0; i<12; i+=2 )
		printf("%f\t", x[i]);
	printf("\n");
	for( i=1; i<12; i+=2 )
		printf("%f\t", x[i]);
	printf("\n\n");
*/

#if CONSTANTS
	ROUND_NEAR;
	printf("Checking out interval constants:\n");
	printf("i_pi_const = [%.21g, %.21g]\n", 
		i_pi_const[0], i_pi_const[1]);
	printf("i_pi_2_const = [%.21g, %.21g]\n", 
		i_pi_2_const[0], i_pi_2_const[1]);
	printf("i_doct_const = [%.21g, %.21g]\n", 
		i_doct_const[0], i_doct_const[1]);
	printf("i_two_pi_5_const = [%.21g, %.21g]\n", 
		i_two_pi_5_const[0], i_two_pi_5_const[1]);
	printf("\n");
#endif

#if DELTA_PARTIALS
	count = 0;
	for( i=0; i<6; i++ ) {
		delta_partials[i] = s_delta_partial_sign( i, x );
		if( delta_partials[i] != 0 )
			count++;
	}
	printf("delta_partials count = %d\n", count);
	
	for( i=0; i<6; i++ )
		delta_partial( i, x, part + 2*i );
	delta_partial( 3, x, out );
	printf("delta_4 is (%f, %f)\n", out[0], out[1]);
#endif
	
	s_delta_partials( x, delta_part );
	
#if DELTA_PARTIALS
	sum = 0.0;
	for( i=0; i<12; i++ ) {
		temp = delta_part[i] - part[i];
		printf("diff = %g\n", fabs(temp));
		sum += temp*temp;
	}
	printf("norm on delta_partials = %g\n", sqrt(sum));
	for( i=0; i<12; i+=2 )
		printf("part = (%f, %f)\n", part[i], part[i+1]);
	
	for( i=0; i<12; i+=2 )
		printf("delta_part = (%f, %f)\n", delta_part[i], delta_part[i+1]);
	ROUND_NEAR;
	max = -100.0;
	for( i=0; i<12; i+=2 ) {
		diff = delta_part[i+1] - delta_part[i];
		if( diff < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
		if( diff > max )
			max = diff;
	}
	printf("maxwidth = %g\n", max);
	printf("\n");
#endif

#if DELTA_PARTIALS
	i_delta_partials( x, part );
	
	sum = 0.0;
	for( i=0; i<12; i++ ) {
		temp = delta_part[i] - part[i];
		printf("diff = %g\n", fabs(temp));
		sum += temp*temp;
	}
	printf("norm on i_delta_partials = %g\n", sqrt(sum));
	for( i=0; i<12; i+=2 )
		printf("part = (%f, %f)\n", part[i], part[i+1]);
	
	i_delta_partials( x, delta_part );
	for( i=0; i<12; i+=2 )
		printf("i_delta_part = (%f, %f)\n", delta_part[i], delta_part[i+1]);
	ROUND_NEAR;
	max = -100.0;
	for( i=0; i<12; i+=2 ) {
		diff = delta_part[i+1] - delta_part[i];
		if( diff < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
		if( diff > max )
			max = diff;
	}
	printf("maxwidth = %g\n", max);
	printf("\n");
#endif

	/* First find delta_min */
	for( i=0; i<12; i++ )
		xp[i] = x[i];		/* copy x, then modify copy */
	j = 0;
	for( i=0; i<6; i++ ) {
		data[0] = delta_part[j];
		data[1] = delta_part[j+1];
		if( data[0] > 0 )	/* positive sign */
			xp[j+1] = xp[j];
		if( data[1] < 0 )	/* negative sign */
			xp[j] = xp[j+1];
		j += 2;
	}
	delta[0] = rough_min_delta( xp );
	
	/* Now find delta_max */
	for( i=0; i<12; i++ )
		xp[i] = x[i];		/* copy x, then modify copy */
	j = 0;
	for( i=0; i<6; i++ ) {
		data[0] = delta_part[j];
		data[1] = delta_part[j+1];
		if( data[0] > 0 )	/* positive sign */
			xp[j] = xp[j+1];
		if( data[1] < 0 )	/* negative sign */
			xp[j+1] = xp[j];
		j += 2;
	}
	delta[1] = rough_max_delta( xp );
	
	i_sqrt( delta, sqrtdelta );
	ROUND_DOWN;
	oosqrtdelta[0] = 1.0/sqrtdelta[1];
	ROUND_UP;
	oosqrtdelta[1] = 1.0/sqrtdelta[0];
	
	a_partials( y, a_part );
	
#if A_PARTIALS
	for( i=0; i<12; i+=2 )
		printf("a_part = (%f, %f)\n", a_part[i], a_part[i+1]);
	printf("\n");
#endif

	s_solid_xpars( y, ooy, delta, sqrtdelta, delta_part, sol_xpart );

#if SOL_PARTIALS
	ROUND_NEAR;
	for( i=0; i<12; i+=2 )
		printf("sol_xpart = (%f, %f)\n", sol_xpart[i], sol_xpart[i+1]);
	max = -100.0;
	for( i=0; i<12; i+=2 ) {
		diff = sol_xpart[i+1] - sol_xpart[i];
		if( diff < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
		if( diff > max )
			max = diff;
	}
	printf("maxwidth = %g\n", max);
	printf("\n");
#endif

	s_solid_partials( y, delta, sqrtdelta, delta_part, sol_part );

	/* y to x pars */
	for( i=0; i<12; i+=2 ) {
		i_div( sol_part + i, y + i, ck_part + i );
		ck_part[i] *= 0.5;
		ck_part[i+1] *= 0.5;
	}

#if SOL_PARTIALS
	ROUND_NEAR;
	for( i=0; i<12; i+=2 )
		printf("ck_part = (%f, %f)\n", ck_part[i], ck_part[i+1]);

	ROUND_NEAR;
	for( i=0; i<12; i+=2 )
		printf("sol_part = (%f, %f)\n", sol_part[i], sol_part[i+1]);
	max = -100.0;
	for( i=0; i<12; i+=2 ) {
		diff = sol_part[i+1] - sol_part[i];
		if( diff < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
		if( diff > max )
			max = diff;
	}
	printf("maxwidth = %g\n", max);
	printf("\n");
#endif

	s_bvol_xpars( y, ooy, delta, sqrtdelta, delta_part, 
		sol_xpart, bvol_xpart );

#if BVOL_PARTIALS
	ROUND_NEAR;
	for( i=0; i<12; i+=2 )
		printf("bvol_xpart = (%f, %f)\n", bvol_xpart[i], 
			bvol_xpart[i+1]);
	max = -100.0;
	for( i=0; i<12; i+=2 ) {
		diff = bvol_xpart[i+1] - bvol_xpart[i];
		if( diff < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
		if( diff > max )
			max = diff;
	}
	printf("maxwidth = %g\n", max);
	printf("\n");
#endif

	s_bvol_partials( y, delta, sqrtdelta, delta_part, sol_part, bvol_part );

	/* y to x pars */
	for( i=0; i<12; i+=2 ) {
		i_div( bvol_part + i, y + i, ck_part + i );
		ck_part[i] *= 0.5;
		ck_part[i+1] *= 0.5;
	}
	
#if BVOL_PARTIALS
	ROUND_NEAR;
	for( i=0; i<12; i+=2 )
		printf("ck_part = (%f, %f)\n", ck_part[i], ck_part[i+1]);

	for( i=0; i<12; i+=2 )
		printf("bvol_part = (%f, %f)\n", bvol_part[i], bvol_part[i+1]);
	ROUND_NEAR;
	max = -100.0;
	for( i=0; i<12; i+=2 ) {
		diff = bvol_part[i+1] - bvol_part[i];
		if( diff < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
		if( diff > max )
			max = diff;
	}
	printf("maxwidth = %g\n", max);
	printf("\n");
#endif

	s_dih_partials( x, y, dih_part );

#if DIH_PARTIALS
	for( i=0; i<12; i+=2 )
		printf("dih_part = (%f, %f)\n", dih_part[i], dih_part[i+1]);
	ROUND_NEAR;
	max = -100.0;
	for( i=0; i<12; i+=2 ) {
		diff = dih_part[i+1] - dih_part[i];
		if( diff < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
		if( diff > max )
			max = diff;
	}
	printf("maxwidth = %g\n", max);
	printf("\n");
#endif

	s_gma_xpars( oosqrtdelta, delta_part, 
		bvol_xpart, gma_xpart );

#if GMA_PARTIALS
	ROUND_NEAR;
	for( i=0; i<12; i+=2 )
		printf("gma_xpart = (%f, %f)\n", gma_xpart[i], 
			gma_xpart[i+1]);
	max = -100.0;
	for( i=0; i<12; i+=2 ) {
		diff = gma_xpart[i+1] - gma_xpart[i];
		if( diff < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
		if( diff > max )
			max = diff;
	}
	printf("maxwidth = %g\n", max);
	printf("\n");
#endif

	s_gma_partials( y, sqrtdelta, delta_part, bvol_part, 
		gma_part );

	/* y to x pars */
	for( i=0; i<12; i+=2 ) {
		i_div( gma_part + i, y + i, ck_part + i );
		ck_part[i] *= 0.5;
		ck_part[i+1] *= 0.5;
	}
	
#if GMA_PARTIALS
	ROUND_NEAR;
	for( i=0; i<12; i+=2 )
		printf("ck_part = (%f, %f)\n", ck_part[i], ck_part[i+1]);

	for( i=0; i<12; i+=2 )
		printf("gma_part = (%f, %f)\n", gma_part[i], gma_part[i+1]);
	ROUND_NEAR;
	max = -100.0;
	for( i=0; i<12; i+=2 ) {
		diff = gma_part[i+1] - gma_part[i];
		if( diff < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
		if( diff > max )
			max = diff;
	}
	printf("maxwidth = %g\n", max);
	printf("\n");
#endif

	/*
	slope = 0.4003278004833681;
	correps = 0.0385;
	for( i=0; i<12; i+=2 ) {
		j = i+1;
		ROUND_DOWN;
		rel_part[i] = gma_part[i] + slope*sol_part[i] + 
			correps*dih_part[i];
		ROUND_UP;
		rel_part[j] = gma_part[j] + slope*sol_part[j] + 
			correps*dih_part[j];
	}
	for( i=0; i<12; i+=2 )
		printf("rel_part = (%f, %f)\n", rel_part[i], rel_part[i+1]);
	printf("\n");
	*/
	
	i_mult( delta, sqrtdelta, delta32 );
	tomsP_partials( x, delta, delta32, delta_part, tomsP_part );

#if TOMSP_PARTIALS
	for( i=0; i<12; i+=2 )
		printf("tomsP_part = (%f, %f)\n", tomsP_part[i], tomsP_part[i+1]);
	ROUND_NEAR;
	max = -100.0;
	for( i=0; i<12; i+=2 ) {
		diff = tomsP_part[i+1] - tomsP_part[i];
		if( diff < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
		if( diff > max )
			max = diff;
	}
	printf("maxwidth = %g\n", max);
	printf("\n");
#endif

	vorvol_partials( y, x, delta, sqrtdelta, delta_part, vorvol_part );

#if VORVOL_PARTIALS
	for( i=0; i<12; i+=2 )
		printf("vorvol_part = (%f, %f)\n", vorvol_part[i], vorvol_part[i+1]);
	ROUND_NEAR;
	max = -100.0;
	for( i=0; i<12; i+=2 ) {
		diff = vorvol_part[i+1] - vorvol_part[i];
		if( diff < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
		if( diff > max )
			max = diff;
	}
	printf("maxwidth = %g\n", max);
	printf("\n");
#endif

	vor_partials( y, x, delta, sqrtdelta, delta_part, 
		sol_part, vor_part );

#if VOR_PARTIALS
	for( i=0; i<12; i+=2 )
		printf("vor_part = (%f, %f)\n", 
			vor_part[i], vor_part[i+1]);
	ROUND_NEAR;
	max = -100.0;
	for( i=0; i<12; i+=2 ) {
		diff = vor_part[i+1] - vor_part[i];
		if( diff < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
		if( diff > max )
			max = diff;
	}
	printf("maxwidth = %g\n", max);
	printf("\n");
#endif

	/* y to x pars */
	for( i=0; i<12; i+=2 ) {
		i_div( vor_part + i, y + i, ck_part + i );
		ck_part[i] *= 0.5;
		ck_part[i+1] *= 0.5;
	}

	vor_xpars( x, delta, sqrtdelta, delta_part, sol_xpart, 
		vor_xpart );

#if VOR_PARTIALS
	ROUND_NEAR;
	for( i=0; i<12; i+=2 )
		printf("ck_part = (%f, %f)\n", ck_part[i], ck_part[i+1]);

	for( i=0; i<12; i+=2 )
		printf("vor_xpart = (%f, %f)\n", 
			vor_xpart[i], vor_xpart[i+1]);
	ROUND_NEAR;
	max = -100.0;
	for( i=0; i<12; i+=2 ) {
		diff = vor_xpart[i+1] - vor_xpart[i];
		if( diff < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
		if( diff > max )
			max = diff;
	}
	printf("maxwidth = %g\n", max);
	printf("\n");
#endif

	octa_vor_partials( y, x, delta, sqrtdelta, delta_part, 
		sol_part, vor_part );

#if OCTA_VOR_PARTIALS
	ROUND_NEAR;
	for( i=0; i<12; i+=2 )
		printf("octa_vor_part = [%.15g, %.15g]\n", vor_part[i], vor_part[i+1]);
	ROUND_NEAR;
	max = -100.0;
	for( i=0; i<12; i+=2 ) {
		diff = vor_part[i+1] - vor_part[i];
		if( diff < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
		if( diff > max )
			max = diff;
	}
	printf("maxwidth = %g\n", max);
	printf("\n");
#endif

	/* y to x pars */
	for( i=0; i<12; i+=2 ) {
		i_div( vor_part + i, y + i, ck_part + i );
		ck_part[i] *= 0.5;
		ck_part[i+1] *= 0.5;
	}

	octa_vor_xpars( y, ooy, x, delta, sqrtdelta, delta_part, 
		sol_xpart, vor_part );

#if OCTA_VOR_PARTIALS
	ROUND_NEAR;
	for( i=0; i<12; i+=2 )
		printf("ck_part = (%f, %f)\n", ck_part[i], ck_part[i+1]);

	ROUND_NEAR;
	for( i=0; i<12; i+=2 )
		printf("octa_vor_xpart = [%.15g, %.15g]\n", vor_part[i], vor_part[i+1]);
	ROUND_NEAR;
	max = -100.0;
	for( i=0; i<12; i+=2 ) {
		diff = vor_part[i+1] - vor_part[i];
		if( diff < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
		if( diff > max )
			max = diff;
	}
	printf("maxwidth = %g\n", max);
	printf("\n");
#endif

	ROUND_NEAR;
	/*	return;	*/

#if ROUGH
	temp = rough_min_delta(  x );
	printf("rough_min_delta gives %.21g\n", temp);
	printf("s_min_delta gives %.21g\n", delta[0]);
	printf("diff = %g\n", delta[0] - temp);
	ROUND_NEAR;
	dout = bigdelta( sx );
	printf("bigdelta gives %.21g\n", dout);
	diff = dout - temp;
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));

	temp = rough_max_delta(  x );
	printf("rough_max_delta gives %.21g\n", temp);
	printf("s_max_delta gives %.21g\n", delta[1]);
	printf("diff = %g\n", temp - delta[1]);
	ROUND_NEAR;
	dout = bigdelta( sx );
	printf("bigdelta gives %.21g\n", dout);
	diff = temp - dout;
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
	
	temp = rough_max_bvol(  y );
	printf("rough_max__bvol gives %.21g\n", temp);
	temp2 = s_max_bvol( y, sqrtdelta );
	printf("s_max_bvol gives %.21g\n", temp2);
	printf("diff = %g\n", temp - temp2);
	ROUND_NEAR;
	dout = bvol( sy );
	printf("bvol gives %.21g\n", dout);
	diff = temp - dout;
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
#endif

#if MINMAX_A
	temp = min_a(  y );
	printf("min_a gives %.21g\n", temp);
	ROUND_NEAR;
	dout = afunc( sy );
	printf("afunc gives %.21g\n", dout);
	diff = dout - temp;
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
	
	temp = max_a(  y );
	printf("max_a gives %.21g\n", temp);
	ROUND_NEAR;
	dout = afunc( sy );
	printf("afunc gives %.21g\n", dout);
	diff = temp - dout;
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
#endif

#if ROUGH
	temp = rough_max_solid(  y );
	printf("rough_max_solid gives %.21g\n", temp);
	temp2 = max_solid( y, sqrtdelta );
	printf("max_solid gives %.21g\n", temp2);
	printf("diff = %g\n", temp - temp2);
	ROUND_NEAR;
	dout = solid( sy );
	printf("solid gives %.21g\n", dout);
	diff = temp - dout;
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
	
	i_solid(  y, sqrtdelta, out );
	printf("i_solid gives (%.21g, %.21g)\n", out[0], out[1]);
	printf("out[1]-out[0] = %.21g\n", out[1]-out[0]);
	ROUND_NEAR;
	dout = solid( sy );
	printf("solid gives %.21g\n", dout);
	diff = fabs(2.0*dout - out[0] - out[1]);
	printf("diff = %.21g\n\n", diff);
#endif	

#if ROUGH
	temp = rough_min_tvol(  y );
	printf("rough_min_tvol gives %.21g\n", temp);
	ROUND_NEAR;
	dout = tetvolume( sy );
	printf("tetvolume gives %.21g\n", dout);
	diff = dout - temp;
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
	
	temp = rough_max_tvol(  y );
	printf("rough_max_tvol gives %.21g\n", temp);
	ROUND_NEAR;
	dout = tetvolume( sy );
	printf("tetvolume gives %.21g\n", dout);
	diff = temp - dout;
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
#endif

#if ROUGH
	temp = rough_min_gma(  y );
	printf("rough_min_gma gives %.21g\n", temp);
	temp2 = s_min_gma( y, sqrtdelta );
	printf("s_min_gma gives %.21g\n", temp2);
	printf("diff = %g\n", temp2 - temp);
	ROUND_NEAR;
	dout = gma( sy );
	printf("gma gives %.21g\n", dout);
	diff = dout - temp2;
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));

	temp = rough_max_gma(  y );
	printf("rough_max_gma gives %.21g\n", temp);
	temp2 = s_max_gma( y, sqrtdelta );
	printf("s_max_gma gives %.21g\n", temp2);
	printf("diff = %g\n", temp - temp2);
	ROUND_NEAR;
	dout = gma( sy );
	printf("gma gives %.21g\n", dout);
	diff = temp2 - dout;
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
#endif

#if GMA
	temp = s_min_gma( y, sqrtdelta );
	temp2 = s_max_gma( y, sqrtdelta );
	ROUND_NEAR;
	printf("s_gma gives [%.21g, %.21g]\n", temp, temp2);
	diff = temp2 - temp;
	printf("diff = %.21g\n\n", diff);
#endif

#if DIH
	temp = s_max_dih(  x );
	printf("s_max_dih gives %.21g\n", temp);
	ROUND_NEAR;
	dout = dih( sy );
	printf("dih gives %.21g\n", dout);
	diff = temp - dout;
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
	
	temp = s_min_dih(  x );
	printf("s_min_dih gives %.21g\n", temp);
	ROUND_NEAR;
	dout = dih( sy );
	printf("dih gives %.21g\n", dout);
	diff = dout - temp;
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
#endif

#if MAXMIN_U
	temp = s_min_u(  y );
	printf("s_min_u gives %.21g\n", temp);
	ROUND_NEAR;
	dout = tomsu( sy );
	printf("tomsu gives %.21g\n", dout);
	diff = dout - temp;
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));

	temp = s_max_u(  y );
	printf("s_max_u gives %.21g\n", temp);
	ROUND_NEAR;
	dout = tomsu( sy );
	printf("tomsu gives %.21g\n", dout);
	diff = temp - dout;
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
#endif

#if QR_CRAD
	temp = min_qr_crad( y );
	printf("min_qr_crad gives %.21g\n", temp);
	ROUND_NEAR;
	dout = circumradius( sy );
	printf("circumradius gives %.21g\n", dout);
	diff = dout - temp;
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));

	temp = max_qr_crad(  y );
	printf("max_qr_crad gives %.21g\n", temp);
	ROUND_NEAR;
	dout = circumradius( sy );
	printf("circumradius gives %.21g\n", dout);
	diff = temp - dout;
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
#endif

	/*
	i_dbvol(  y,  out );
	i_bigdelta(  x,  out );
	i_istet(  y );
	i_tomsrho(  x,  out );
	*/
	
#if CRAD
	i_crad(  y,  out );
	printf("i_crad gives (%.21g, %.21g)\n", out[0], out[1]);
	printf("out[1]-out[0] = %.21g\n", out[1]-out[0]);
	ROUND_NEAR;
	dout = circumradius( sy );
	printf("crad gives %.21g\n", dout);
	diff = fabs(2.0*dout - out[0] - out[1]);
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
#endif

#if SCRAD
	c126 = s_crad2_6(  x,  out3 );
	printf("c126 = %d\n", c126 );
	i_sqrt( out3, out2 );
	printf("s_crad2_6 gives (%.21g, %.21g)\n", out3[0], out3[1]);
	printf("i_sqrt gives  (%.21g, %.21g)\n", out2[0], out2[1]);
	printf("out[1]-out2[1] = %.21g\n", out[1]-out2[1]);
	printf("out2[0]-out[0] = %.21g\n", out2[0]-out[0]);
	ROUND_NEAR;
	dout = circumradius( sy );
	printf("crad gives %.21g\n", dout);
	diff = fabs(2.0*dout - out2[0] - out2[1]);
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));

	s_crad2_1(  x,  out3, finfo );
	c126 = finfo[0];
	c135 = finfo[1];
	printf("c126 = %d\t", c126 );
	printf("c135 = %d\n", c135 );
	i_sqrt( out3, out2 );
	printf("s_crad2_1 gives (%.21g, %.21g)\n", out3[0], out3[1]);
	printf("i_sqrt gives  (%.21g, %.21g)\n", out2[0], out2[1]);
	printf("out[1]-out2[1] = %.21g\n", out[1]-out2[1]);
	printf("out2[0]-out[0] = %.21g\n", out2[0]-out[0]);
	ROUND_NEAR;
	dout = circumradius( sy );
	printf("crad gives %.21g\n", dout);
	diff = fabs(2.0*dout - out2[0] - out2[1]);
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
#endif
	/*
	i_tomsu(  x,  out );
	i_tomsv(  x,  out );
	i_auxPfun(  x,  out );
	i_tomsP(  x,  out );
	i_voronoivol(  y,  out );
	*/
	i_vor(  y, x, sqrtdelta, out );
	ROUND_NEAR;
	printf("i_vor gives (%.21g, %.21g)\n", out[0], out[1]);
	printf("diff = %.21g\n", out[1]-out[0]);
	ROUND_NEAR;

	best_i_vor_1(  y, x, sqrtdelta, out );
	ROUND_NEAR;
	printf("best_i_vor_1 gives (%.21g, %.21g)\n", out[0], out[1]);
	printf("diff = %.21g\n", out[1]-out[0]);
	ROUND_NEAR;
	
/*
	dout = vor( sy );
	printf("vor gives %.21g\n", dout);
	diff = fabs(2.0*dout - out[0] - out[1]);
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
*/

#if VOR
	s_density_vor_6(  y, x, out2 );
	ROUND_NEAR;
	printf("s_density_vor_6 gives (%.21g, %.21g)\n", 
		out2[0], out2[1]);
	printf("out2[1]-out2[0] = %.21g\n", out2[1]-out2[0]);
/*
	ROUND_NEAR;
	printf("out[1]-out2[1] = %.21g\n", out[1]-out2[1]);
	printf("out2[0]-out[0] = %.21g\n", out2[0]-out[0]);
	dout = vor( sy );
	printf("vor gives %.21g\n", dout);
	diff = fabs(2.0*dout - out2[0] - out2[1]);
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
*/
#endif

	s_density_vor_1(  y, x, out2 );
	ROUND_NEAR;
	printf("s_density_vor_1 gives (%.21g, %.21g)\n", 
		out2[0], out2[1]);
	printf("diff = %.21g\n\n", out2[1]-out2[0]);

/*
	ROUND_NEAR;
	printf("out[1]-out2[1] = %.21g\n", out[1]-out2[1]);
	printf("out2[0]-out[0] = %.21g\n", out2[0]-out[0]);
	dout = vor( sy );
	printf("vor gives %.21g\n", dout);
	diff = fabs(2.0*dout - out2[0] - out2[1]);
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
*/

/*
	printf("scoring_system gives:\t%d\n", scoring_system( x ) );
	i_score(  y, x, sqrtdelta, out );
	printf("i_score gives (%.21g, %.21g)\n", out[0], out[1]);
	printf("out[1]-out[0] = %.21g\n", out[1]-out[0]);
	ROUND_NEAR;
	dout = score( sy );
	printf("score gives %.21g\n", dout);
	diff = fabs(2.0*dout - out[0] - out[1]);
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
*/

	printf("octa_scoring gives:\t%d\n", octa_scoring( x ) );
	s_octa_vor(  y, x, out );
	printf("s_octa_vor gives (%.21g, %.21g)\n", out[0], out[1]);
	printf("out[1]-out[0] = %.21g\n", out[1]-out[0]);
	ROUND_NEAR;
	dout = score( sy );
	printf("score gives %.21g\n", dout);
	diff = fabs(2.0*dout - out[0] - out[1]);
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));

/*
	out2[0] = max_solid( y, sqrtdelta );
	out2[1] = out[1] + 0.396*out2[0] - 0.5*0.5383016645753285;
	printf("sph = %.21g,  sc = %.21g\n", out2[0], out[1]);
	printf("relation test:  %.21g\n\n", out2[1]);
*/

#if CRAD3LEN
	i_crad3len(  y, out );
	printf("i_crad3len gives (%.21g, %.21g)\n", out[0], out[1]);
	printf("out[1]-out[0] = %.21g\n", out[1]-out[0]);
	ROUND_NEAR;
	dout = crad3len( sy );
	printf("crad3len gives %.21g\n", dout);
	diff = fabs(2.0*dout - out[0] - out[1]);
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
	
	s_crad3x(  x, out );
	printf("s_crad3x gives (%.21g, %.21g)\n", out[0], out[1]);
	printf("out[1]-out[0] = %.21g\n", out[1]-out[0]);
	ROUND_NEAR;
	dout = crad3len( sy );
	printf("crad3len gives %.21g\n", dout);
	diff = fabs(2.0*dout - out[0] - out[1]);
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
#endif
	
	/*
	i_voronoivol(x, sqrtdelta, out );
	ROUND_NEAR;
	printf("i_voronoivol gives (%.21g, %.21g)\n", out[0], out[1]);
	printf("diff = %.21g\n", out[1]-out[0]);
	*/
	i_vor_trunc( y, x, sqrtdelta, out );
	ROUND_NEAR;
	printf("i_vor_trunc gives (%.21g, %.21g)\n", out[0], out[1]);
	printf("diff = %.21g\n", out[1]-out[0]);
	/*
	alt_vor_trunc( y, x, sqrtdelta, out2, out );
	ROUND_NEAR;
	printf("alt_vor_trunc gives (%.21g, %.21g)\n", 
		out[0], out[1]);
	printf("diff = %.21g\n", out[1]-out[0]);
	
	best_i_vor_6(  y, x, sqrtdelta, out );
	ROUND_NEAR;
	printf("best_i_vor_6 gives (%.21g, %.21g)\n", out[0], out[1]);
	printf("diff = %.21g\n", out[1]-out[0]);
	ROUND_NEAR;
	
	i_crad3len(  y, ab + 2 );
	ROUND_NEAR;
	printf("i_crad3len gives (%.21g, %.21g)\n", ab[2], ab[3]);
	i_tomsu( x, out );
	ROUND_NEAR;
	printf("i_tomsu32 gives (%.21g, %.21g)\n", 
		pow( out[0], 1.5 ), pow( out[1], 1.5 ) );
	*/
	
	crad3len_pars( y, x, part );
	ROUND_NEAR;
	printf("crad3len_partials = \n");
	for( i=0; i<6; i+=2 ) {
		printf("[%.18g, %.18g]\n", part[i], part[i+1]);
	}

	crad3len2_xpars( x, part );
	ROUND_NEAR;
	printf("crad3len2_xpartials = \n");
	for( i=0; i<6; i+=2 ) {
		printf("[%.18g, %.18g]\n", part[i], part[i+1]);
	}
	
	/*
	s_crad3x(  x + 2, ab + 2 );
	ROUND_NEAR;
	printf("s_crad3x gives (%.21g, %.21g)\n", ab[2], ab[3]);
	s_crad3x2(  x + 2, xy + 2 );
	ROUND_NEAR;
	printf("s_crad3x2 gives (%.21g, %.21g)\n", xy[2], xy[3]);

	ab[0] = 0.5*y[2];
	ab[1] = 0.5*y[3];
	xy[0] = 0.25*x[2];
	xy[1] = 0.25*x[3];
	
	printf("a = [%.18g, %.18g]\n", ab[0], ab[1] );
	printf("b = [%.18g, %.18g]\n", ab[2], ab[3] );
	printf("x = [%.18g, %.18g]\n", xy[0], xy[1] );
	printf("y = [%.18g, %.18g]\n", xy[2], xy[3] );
	printf("\n");
	*/
	/*
	old_dih_rog( xy, al );
	ROUND_NEAR;
	printf("old_alpha =  [%.21g, %.21g]\n", al[0], al[1] );
	alpha_ab( xy, al_ab );
	printf("\n");
	
	i_dih_rog( xy, al );
	ROUND_NEAR;
	printf("alpha =  [%.21g, %.21g]\n", al[0], al[1] );
	alpha_ab( xy, al_ab );
	printf("\n");
	
	ROUND_NEAR;
	printf("alpha_a =  [%.21g, %.21g]\n", 
		al_ab[0], al_ab[1] );
	printf("alpha_b =  [%.21g, %.21g]\n", 
		al_ab[2], al_ab[3] );
	if( al_ab[0] > al_ab[1] || al_ab[2] > al_ab[3] )
		printf("Wups.  Swapped interval.\n");
	printf("\n");
	wedgevol_ab( ab, xy, al, al_ab, part );
	ROUND_NEAR;
	printf("wedgevol_a =  [%.21g, %.21g]\n", 
		part[0], part[1] );
	printf("wedgevol_b =  [%.21g, %.21g]\n", 
		part[2], part[3] );
	if( part[0] > part[1] || part[2] > part[3] )
		printf("Wups.  Swapped interval.\n");
	printf("\n");
	wedgesol_ab( ab, al, al_ab, part );
	ROUND_NEAR;
	printf("wedgesol_a =  [%.21g, %.21g]\n", 
		part[0], part[1] );
	printf("wedgesol_b =  [%.21g, %.21g]\n", 
		part[2], part[3] );
	if( part[0] > part[1] || part[2] > part[3] )
		printf("Wups.  Swapped interval.\n");
	printf("\n");
	rogvol_ab( ab, xy, part );
	ROUND_NEAR;
	printf("rogvol_a =  [%.21g, %.21g]\n", 
		part[0], part[1] );
	printf("rogvol_b =  [%.21g, %.21g]\n", 
		part[2], part[3] );
	if( part[0] > part[1] || part[2] > part[3] )
		printf("Wups.  Swapped interval.\n");
	printf("\n");
	rogsol_ab( ab, xy, part );
	ROUND_NEAR;
	printf("rogsol_a =  [%.21g, %.21g]\n", 
		part[0], part[1] );
	printf("rogsol_b =  [%.21g, %.21g]\n", 
		part[2], part[3] );
	if( part[0] > part[1] || part[2] > part[3] )
		printf("Wups.  Swapped interval.\n");
	printf("\n");
	*/
	/* printf("Starting trunc_pars . . .\n"); */
	trunc_pars( y, x, sol_part, part );
	
	ROUND_NEAR;
	for( i=0; i<12; i+=2 ) {
		printf("trunc_part = (%.18f, %.18f)\n", 
			part[i], part[i+1]);
		diff = part[i+1] - part[i];
		if( diff < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
	}
	printf("\n");
	
	/*
	approx_trunc_pars( y, part );
	ROUND_NEAR;
	for( i=0; i<12; i+=2 ) {
		printf("approx_part = (%.18f, %.18f)\n", 
			part[i], part[i+1]);
		diff = part[i+1] - part[i];
		if( diff < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
	}
	printf("\n");
	*/
	/*
	approx_rog_wed_pars( y, vpart, spart );
	ROUND_NEAR;
	for( i=0; i<12; i+=2 ) {
		printf("approx_rog_wed_vpart = (%.18f, %.18f)\n", 
			vpart[i], vpart[i+1]);
		diff = vpart[i+1] - vpart[i];
		if( diff < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
	}
	printf("\n");
	for( i=0; i<12; i+=2 ) {
		printf("approx_rog_wed_spart = (%.18f, %.18f)\n", 
			spart[i], spart[i+1]);
		diff = spart[i+1] - spart[i];
		if( diff < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
	}
	printf("\n");
	*/
	/*
	approx_fullwed_pars( y, vpart, spart );
	ROUND_NEAR;
	for( i=0; i<12; i+=2 ) {
		printf("approx_fullwed_vpart = (%.18f, %.18f)\n", 
			vpart[i], vpart[i+1]);
		diff = vpart[i+1] - vpart[i];
		if( diff < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
	}
	printf("\n");
	for( i=0; i<12; i+=2 ) {
		printf("approx_fullwed_spart = (%.18f, %.18f)\n", 
			spart[i], spart[i+1]);
		diff = spart[i+1] - spart[i];
		if( diff < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
	}
	printf("\n");
	*/
	/* test_ab_pars( ab ); */
	
	ROUND_NEAR;
}
/* 2.029102 2.0341248 2.011234 2.023309 2.0123 2.83982134 */


void get_second_trash( void )
{
	double y[12], sy[7], eps, ubv_ij[6][12];
	double x[12], sx[7], temp[2];
	int i, j;
	int interval;
	/*
	double y[12], out[2], sy[7], dout, diff, eps, temp, temp2;
	double x[12], xp[12], delta[2], sx[7], delta_part[12];
	double part[12], sum, data[2], sqrtdelta[2], delta32[2];
	double sol_part[12], a_part[12], bvol_part[12];
	double dih_part[12], gma_part[12], tomsP_part[12];
	double rel_part[12], vorvol_part[12], vor_part[12];
	double slope, correps, max, out2[2], out3[2];
	int i, j, delta_partials[6], count, c126, c135, finfo[2];
	int interval;
	*/
	

	i_init();
	
	ROUND_NEAR;
	eps = 0.0;
	
	printf("Enter data type for cells:  1 for intervals, 0 otherwise:  ");
	scanf("%d", &interval);
	
	if( interval == 0 ) {
		printf("Enter edge lengths:\n");
		
		for( i=1; i<7; i++ )
			scanf("%lf", sy + i);
		printf("Got the following:  \n");
		for( i=1; i<7; i++ )
			printf("%f\t", sy[i]);
		printf("\n\n");
		printf("Enter epsilon (interval width):  ");
		scanf("%lf", &eps);
		
		ROUND_UP;
		for( i=0; i<6; i++ ) {
			y[2*i] = sy[i+1];
			y[2*i+1] = sy[i+1] + eps;
			}
		ROUND_NEAR;
		
		for( i=1; i<7; i++ )
			sx[i] = sy[i]*sy[i];
	} else {
	printf("Enter edge lengths (as intervals):\n");
	for( i=0; i<12; i++ )
		scanf("%lf", y + i );
	for( i=0; i<6; i++ )
		sy[i+1] = y[2*i];
	}

	/* Compute square of edge lengths */
	ROUND_UP;
	for( i=1; i<12; i+=2 )
		x[i] = y[i]*y[i];
	ROUND_DOWN;
	for( i=0; i<12; i+=2 )
		x[i] = y[i]*y[i];
	
	sp_gmavolsecparbds( x, ubv_ij );
	/* sp_vorsecparbds( x, ubv_ij ); */
	/* sp_secparbounds( x, ubv_ij ); */
	/* sp_solsecparbds( y, x, ubv_ij ); */
	/* sp_dihsecparbds( x, ubv_ij ); */
	sp_findmaxmin( ubv_ij, temp );
	
	ROUND_NEAR;
	printf("ubv_ij = \n");
	for( i=0; i<6; i++ ) {
		for( j=2*i; j<12; j+=2 ) {
			printf("[%.18g,\t%.18g]\n", 
				ubv_ij[i][j], ubv_ij[i][j+1]);
		}
		printf("\n");
	}
	printf("sp_findmaxmin = [%.18g,\t%.18g]\n",
		temp[0], temp[1]);

	printf("p4mat = {\n");
	for( i=0; i<6; i++ ) {
		printf("{");
		for( j=2*i; j<12; j+=2 ) {
			printf("%.18g", 
				0.5*(ubv_ij[i][j] + ubv_ij[i][j+1]));
			if( j != 10 )
				printf(",");
			if( j == 4 )
				printf("\n");
		}
		printf("}");
		if( i != 5 )
			printf(",\n");
	}
	printf("};\n");
}


void get_taylor_trash( void )
{
	double y[12], out[2], sy[7], dout, diff, eps, temp, temp2;
	double x[12], xp[12], delta[2], sx[7], delta_part[12];
	double part[12], sum, data[2], sqrtdelta[2], delta32[2];
	double sol_part[12], a_part[12], bvol_part[12];
	double dih_part[12], gma_part[12], tomsP_part[12];
	double rel_part[12], vorvol_part[12], vor_part[12];
	double slope, correps, max, out2[2], out3[2];
	double xt[12], yt[12], xh[12], yh[12], yh2, xh2;
	double sol_val[2], vorvol_val[2], ooyt[12];
	double oosqrtdelta[2];
	int i, j, delta_partials[6], count, c126, c135, finfo[2];
	int interval;	

	i_init();
	
	ROUND_NEAR;
	eps = 0.0;

	interval = 0;
	/*
	printf("Enter data type for cells:  1 for intervals, 0 otherwise:  ");
	scanf("%d", &interval);
	*/
	if( interval == 0 ) {
		printf("Enter edge lengths:\n");
		
		for( i=1; i<7; i++ )
			scanf("%lf", sy + i);
		/*
		printf("Got the following:  \n");
		for( i=1; i<7; i++ )
			printf("%f\t", sy[i]);
		printf("\n\n");
		*/
		printf("Enter epsilon (interval width):  ");
		scanf("%lf", &eps);
		
		ROUND_UP;
		for( i=0; i<6; i++ ) {
			y[2*i] = sy[i+1];
			y[2*i+1] = sy[i+1] + eps;
			}
		ROUND_NEAR;
		
		for( i=1; i<7; i++ )
			sx[i] = sy[i]*sy[i];
	} else {
	printf("Enter edge lengths (as intervals):\n");
	for( i=0; i<12; i++ )
		scanf("%lf", y + i );
	for( i=0; i<6; i++ )
		sy[i+1] = y[2*i];
	}

/*
	for( i=0; i<12; i+=2 )
		printf("%f\t", y[i]);
	printf("\n");
	for( i=1; i<12; i+=2 )
		printf("%f\t", y[i]);
	printf("\n\n");
*/	

	/* Compute square of edge lengths */
	ROUND_UP;
	for( i=1; i<12; i+=2 )
		x[i] = y[i]*y[i];
	ROUND_DOWN;
	for( i=0; i<12; i+=2 )
		x[i] = y[i]*y[i];

	/* Compute t and h for both x and y */
	ROUND_DOWN;
	for( i=0; i<12; i+=2 ) {
		xt[i] = 0.5*(x[i+1] + x[i]);
		yt[i] = 0.5*(y[i+1] + y[i]);
	}
	ROUND_UP;
	for( i=0; i<12; i+=2 ) {
		xt[i+1] = 0.5*(x[i+1] + x[i]);
		yt[i+1] = 0.5*(y[i+1] + y[i]);

		xh[i+1] = 0.5*(x[i+1] - x[i]);
		xh[i] = -xh[i+1];
		/* printf("xh[%d] = %.18g\n", i, xh[i]); */
		
		yh[i+1] = 0.5*(y[i+1] - y[i]);
		yh[i] = -yh[i+1];
		/* printf("yh[%d] = %.18g\n", i,yh[i]); */
	}
	/* Compute ooyt = 1/yt */
	ROUND_UP;
	for( i=1; i<12; i+=2 )
		ooyt[i] = 1.0/yt[i-1];
	ROUND_DOWN;
	for( i=0; i<12; i+=2 )
		ooyt[i] = 1.0/yt[i+1];

	/* Compute xh2 and yh2 */
	xh2 = 0.0;
	yh2 = 0.0;
	for( i=1; i<12; i+=2 ) {
		xh2 += xh[i]*xh[i];
		yh2 += yh[i]*yh[i];
	}
	xh2 *= 0.5;
	yh2 *= 0.5;
	for( i=0; i<10; i+=2 ) {
		for( j=i+2; j<12; j+=2 ) {
			xh2 += xh[i]*xh[j];
			yh2 += yh[i]*yh[j];
		}
	}
	
	/*
	printf("y = \n");
	for( i=0; i<12; i+=2 ) {
		printf("[%.18g, \t%.18g]\n", y[i], y[i+1]);
	}
	printf("x = \n");
	for( i=0; i<12; i+=2 ) {
		printf("[%.18g, \t%.18g]\n", x[i], x[i+1]);
	}
	printf("yh = \n");
	for( i=0; i<12; i+=2 ) {
		printf("[%.18g, \t%.18g]\n", yh[i], yh[i+1]);
	}
	printf("yh2 = %.18g\n", yh2);
	printf("xh = \n");
	for( i=0; i<12; i+=2 ) {
		printf("[%.18g, \t%.18g]\n", xh[i], xh[i+1]);
	}
	printf("xh2 = %.18g\n", xh2);
	*/
	
	/* Do preliminary computations */
	s_delta_partials( xt, delta_part );
	delta[0] = rough_min_delta( xt );
	delta[1] = rough_max_delta( xt );
	i_sqrt( delta, sqrtdelta );
	out[0] = 1.0;
	out[1] = 1.0;
	i_div( out, sqrtdelta, oosqrtdelta );

	/*
	ROUND_NEAR;
	printf("t_sqrtdelta gives (%.18g, %.18g)\n", 
		sqrtdelta[0], sqrtdelta[1]);
	*/

	/*
	s_solid_partials( yt, delta, sqrtdelta, delta_part, part );
	for( i=0; i<12; i+=2 )
		printf("sol_part = (%f, %f)\n", part[i], part[i+1]);
	ROUND_NEAR;
	max = -100.0;
	for( i=0; i<12; i+=2 ) {
		diff = part[i+1] - part[i];
		if( diff < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
		if( diff > max )
			max = diff;
	}
	printf("maxwidth = %g\n", max);
	printf("\n");

	vorvol_xpars( xt, delta, sqrtdelta, delta_part, part );
	for( i=0; i<12; i+=2 )
		printf("vorvol_part = (%f, %f)\n", part[i], part[i+1]);
	ROUND_NEAR;
	max = -100.0;
	for( i=0; i<12; i+=2 ) {
		diff = part[i+1] - part[i];
		if( diff < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
		if( diff > max )
			max = diff;
	}
	printf("maxwidth = %g\n", max);
	printf("\n");

	s_dih_xpars( xt, part );
	for( i=0; i<12; i+=2 )
		printf("dih_part = (%f, %f)\n", part[i], part[i+1]);
	ROUND_NEAR;
	max = -100.0;
	for( i=0; i<12; i+=2 ) {
		diff = part[i+1] - part[i];
		if( diff < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
		if( diff > max )
			max = diff;
	}
	printf("maxwidth = %g\n", max);
	printf("\n");

	t_vorvol( xt, delta, sqrtdelta, delta_part, xh, xh2,
		vorvol_val );
	ROUND_NEAR;
	printf("t_vorvol gives (%.18g, %.18g)\n", 
		vorvol_val[0], vorvol_val[1]);
	printf("diff = %.21g\n", vorvol_val[1]-vorvol_val[0]);
	*/
	
	/*
	i_dih(  xt, out );
	ROUND_NEAR;
	printf("i_dih gives (%.18g, %.18g)\n", out[0], out[1]);
	*/
	t_solid(  yt, ooyt, delta, sqrtdelta, delta_part, xh, xh2, 
		sol_val );
	ROUND_NEAR;
	printf("t_solid gives (%.18g, %.18g)\n", 
		sol_val[0], sol_val[1]);
	printf("diff = %.21g\n", sol_val[1]-sol_val[0]);


	t_gma( yt, ooyt, delta, sqrtdelta, delta_part, xh, 
		xh2, out );
	ROUND_NEAR;
	printf("t_gma gives (%.18g, %.18g)\n", out[0], out[1]);
	printf("diff = %.21g\n", out[1]-out[0]);

	/*
	dout = solid( sy );
	printf("solid gives %.21g\n", dout);
	diff = 0.5*(sol_val[1]+sol_val[0]) - dout;
	printf("diff = %.18g\n\n", diff);
	*/

	t_vorvol( xt, delta, sqrtdelta, delta_part, xh, xh2,
		vorvol_val );
	
	
	t_vor( xt, yt, ooyt, delta, sqrtdelta, delta_part, xh, 
		xh2, out );
	ROUND_NEAR;
	printf("t_vor gives (%.18g, %.18g)\n", out[0], out[1]);
	printf("diff = %.21g\n", out[1]-out[0]);
	/*
	dout = vor( sy );
	printf("vor gives %.21g\n", dout);
	diff = 0.5*(out[1]+out[0]) - dout;
	printf("diff = %.18g\n\n", diff);
	*/
	
	t_dih(  xt, xh, xh2, out );
	ROUND_NEAR;
	printf("t_dih gives (%.18g, %.18g)\n", out[0], out[1]);
	printf("diff = %.21g\n", out[1]-out[0]);
	/*
	dout = dih( sy );
	printf("dih gives %.21g\n", dout);
	diff = 0.5*(out[1]+out[0]) - dout;
	printf("diff = %.18g\n\n", diff);
	*/

/*
	i_gma( yt, sqrtdelta, out );
	ROUND_NEAR;
	printf("i_gma gives (%.18g, %.18g)\n", out[0], out[1]);
	printf("diff = %.21g\n", out[1]-out[0]);

	i_bvol( yt, sqrtdelta, out );
	ROUND_NEAR;
	printf("i_bvol gives (%.18g, %.18g)\n", out[0], out[1]);
	printf("diff = %.21g\n", out[1]-out[0]);

	temp = s_min_bvol( yt, sqrtdelta );
	temp2 = s_max_bvol( yt, sqrtdelta );
	ROUND_NEAR;
	printf("s_bvol gives [%.21g, %.21g]\n", temp, temp2);
	diff = temp2 - temp;
	printf("diff = %.21g\n\n", diff);
*/

	printf("\n");
	/* Start here. */

/*	
	for( i=0; i<12; i+=2 )
		printf("%f\t", x[i]);
	printf("\n");
	for( i=1; i<12; i+=2 )
		printf("%f\t", x[i]);
	printf("\n\n");
*/

#if CONSTANTS
	ROUND_NEAR;
	printf("Checking out interval constants:\n");
	printf("i_pi_const = [%.21g, %.21g]\n", 
		i_pi_const[0], i_pi_const[1]);
	printf("i_pi_2_const = [%.21g, %.21g]\n", 
		i_pi_2_const[0], i_pi_2_const[1]);
	printf("i_doct_const = [%.21g, %.21g]\n", 
		i_doct_const[0], i_doct_const[1]);
	printf("i_two_pi_5_const = [%.21g, %.21g]\n", 
		i_two_pi_5_const[0], i_two_pi_5_const[1]);
	printf("\n");
#endif

#if DELTA_PARTIALS
	count = 0;
	for( i=0; i<6; i++ ) {
		delta_partials[i] = s_delta_partial_sign( i, x );
		if( delta_partials[i] != 0 )
			count++;
	}
	printf("delta_partials count = %d\n", count);
	
	for( i=0; i<6; i++ )
		delta_partial( i, x, part + 2*i );
	delta_partial( 3, x, out );
	printf("delta_4 is (%f, %f)\n", out[0], out[1]);
#endif
	
	s_delta_partials( x, delta_part );
	
#if DELTA_PARTIALS
	sum = 0.0;
	for( i=0; i<12; i++ ) {
		temp = delta_part[i] - part[i];
		printf("diff = %g\n", fabs(temp));
		sum += temp*temp;
	}
	printf("norm on delta_partials = %g\n", sqrt(sum));
	for( i=0; i<12; i+=2 )
		printf("part = (%f, %f)\n", part[i], part[i+1]);
	
	for( i=0; i<12; i+=2 )
		printf("delta_part = (%f, %f)\n", delta_part[i], delta_part[i+1]);
	ROUND_NEAR;
	max = -100.0;
	for( i=0; i<12; i+=2 ) {
		diff = delta_part[i+1] - delta_part[i];
		if( diff < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
		if( diff > max )
			max = diff;
	}
	printf("maxwidth = %g\n", max);
	printf("\n");
#endif
	
	/* First find delta_min */
	for( i=0; i<12; i++ )
		xp[i] = x[i];		/* copy x, then modify copy */
	j = 0;
	for( i=0; i<6; i++ ) {
		data[0] = delta_part[j];
		data[1] = delta_part[j+1];
		if( data[0] > 0 )	/* positive sign */
			xp[j+1] = xp[j];
		if( data[1] < 0 )	/* negative sign */
			xp[j] = xp[j+1];
		j += 2;
	}
	delta[0] = rough_min_delta( xp );
	
	/* Now find delta_max */
	for( i=0; i<12; i++ )
		xp[i] = x[i];		/* copy x, then modify copy */
	j = 0;
	for( i=0; i<6; i++ ) {
		data[0] = delta_part[j];
		data[1] = delta_part[j+1];
		if( data[0] > 0 )	/* positive sign */
			xp[j] = xp[j+1];
		if( data[1] < 0 )	/* negative sign */
			xp[j+1] = xp[j];
		j += 2;
	}
	delta[1] = rough_max_delta( xp );
	
	i_sqrt( delta, sqrtdelta );

	/*
	ROUND_NEAR;
	printf("i_sqrtdelta gives (%.18g, %.18g)\n", 
		sqrtdelta[0], sqrtdelta[1]);
	*/
	
	a_partials( y, a_part );
	
#if A_PARTIALS
	for( i=0; i<12; i+=2 )
		printf("a_part = (%f, %f)\n", a_part[i], a_part[i+1]);
	printf("\n");
#endif

	s_solid_partials( y, delta, sqrtdelta, delta_part, sol_part );

#if SOL_PARTIALS
	for( i=0; i<12; i+=2 )
		printf("sol_part = (%f, %f)\n", sol_part[i], sol_part[i+1]);
	ROUND_NEAR;
	max = -100.0;
	for( i=0; i<12; i+=2 ) {
		diff = sol_part[i+1] - sol_part[i];
		if( diff < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
		if( diff > max )
			max = diff;
	}
	printf("maxwidth = %g\n", max);
	printf("\n");
#endif

	s_bvol_partials( y, delta, sqrtdelta, delta_part, sol_part, bvol_part );
	
#if BVOL_PARTIALS
	for( i=0; i<12; i+=2 )
		printf("bvol_part = (%f, %f)\n", bvol_part[i], bvol_part[i+1]);
	ROUND_NEAR;
	max = -100.0;
	for( i=0; i<12; i+=2 ) {
		diff = bvol_part[i+1] - bvol_part[i];
		if( diff < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
		if( diff > max )
			max = diff;
	}
	printf("maxwidth = %g\n", max);
	printf("\n");
#endif

	s_dih_partials( x, y, dih_part );

#if DIH_PARTIALS
	for( i=0; i<12; i+=2 )
		printf("dih_part = (%f, %f)\n", dih_part[i], dih_part[i+1]);
	ROUND_NEAR;
	max = -100.0;
	for( i=0; i<12; i+=2 ) {
		diff = dih_part[i+1] - dih_part[i];
		if( diff < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
		if( diff > max )
			max = diff;
	}
	printf("maxwidth = %g\n", max);
	printf("\n");
#endif

	s_gma_partials( y, sqrtdelta, delta_part, bvol_part, gma_part );

#if GMA_PARTIALS
	for( i=0; i<12; i+=2 )
		printf("gma_part = (%f, %f)\n", gma_part[i], gma_part[i+1]);
	ROUND_NEAR;
	max = -100.0;
	for( i=0; i<12; i+=2 ) {
		diff = gma_part[i+1] - gma_part[i];
		if( diff < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
		if( diff > max )
			max = diff;
	}
	printf("maxwidth = %g\n", max);
	printf("\n");
#endif

	/*
	slope = 0.4003278004833681;
	correps = 0.0385;
	for( i=0; i<12; i+=2 ) {
		j = i+1;
		ROUND_DOWN;
		rel_part[i] = gma_part[i] + slope*sol_part[i] + 
			correps*dih_part[i];
		ROUND_UP;
		rel_part[j] = gma_part[j] + slope*sol_part[j] + 
			correps*dih_part[j];
	}
	for( i=0; i<12; i+=2 )
		printf("rel_part = (%f, %f)\n", rel_part[i], rel_part[i+1]);
	printf("\n");
	*/
	
	i_mult( delta, sqrtdelta, delta32 );
	tomsP_partials( x, delta, delta32, delta_part, tomsP_part );

#if TOMSP_PARTIALS
	for( i=0; i<12; i+=2 )
		printf("tomsP_part = (%f, %f)\n", tomsP_part[i], tomsP_part[i+1]);
	ROUND_NEAR;
	max = -100.0;
	for( i=0; i<12; i+=2 ) {
		diff = tomsP_part[i+1] - tomsP_part[i];
		if( diff < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
		if( diff > max )
			max = diff;
	}
	printf("maxwidth = %g\n", max);
	printf("\n");
#endif

	vorvol_partials( y, x, delta, sqrtdelta, delta_part, vorvol_part );

#if VORVOL_PARTIALS
	for( i=0; i<12; i+=2 )
		printf("vorvol_part = (%f, %f)\n", vorvol_part[i], vorvol_part[i+1]);
	ROUND_NEAR;
	max = -100.0;
	for( i=0; i<12; i+=2 ) {
		diff = vorvol_part[i+1] - vorvol_part[i];
		if( diff < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
		if( diff > max )
			max = diff;
	}
	printf("maxwidth = %g\n", max);
	printf("\n");
#endif

	vor_partials( y, x, delta, sqrtdelta, delta_part, sol_part, vor_part );

#if VOR_PARTIALS
	for( i=0; i<12; i+=2 )
		printf("vor_part = (%f, %f)\n", vor_part[i], vor_part[i+1]);
	ROUND_NEAR;
	max = -100.0;
	for( i=0; i<12; i+=2 ) {
		diff = vor_part[i+1] - vor_part[i];
		if( diff < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
		if( diff > max )
			max = diff;
	}
	printf("maxwidth = %g\n", max);
	printf("\n");
#endif

	octa_vor_partials( y, x, delta, sqrtdelta, delta_part, sol_part, vor_part );

#if OCTA_VOR_PARTIALS
	ROUND_NEAR;
	for( i=0; i<12; i+=2 )
		printf("octa_vor_part = [%.15g, %.15g]\n", vor_part[i], vor_part[i+1]);
	ROUND_NEAR;
	max = -100.0;
	for( i=0; i<12; i+=2 ) {
		diff = vor_part[i+1] - vor_part[i];
		if( diff < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
		if( diff > max )
			max = diff;
	}
	printf("maxwidth = %g\n", max);
	printf("\n");


	printf("\n");
#endif

	ROUND_NEAR;
	/*	return;	*/

#if ROUGH
	temp = rough_min_delta(  x );
	printf("rough_min_delta gives %.21g\n", temp);
	printf("s_min_delta gives %.21g\n", delta[0]);
	printf("diff = %g\n", delta[0] - temp);
	ROUND_NEAR;
	dout = bigdelta( sx );
	printf("bigdelta gives %.21g\n", dout);
	diff = dout - temp;
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));

	temp = rough_max_delta(  x );
	printf("rough_max_delta gives %.21g\n", temp);
	printf("s_max_delta gives %.21g\n", delta[1]);
	printf("diff = %g\n", temp - delta[1]);
	ROUND_NEAR;
	dout = bigdelta( sx );
	printf("bigdelta gives %.21g\n", dout);
	diff = temp - dout;
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
	
	temp = rough_max_bvol(  y );
	printf("rough_max__bvol gives %.21g\n", temp);
	temp2 = s_max_bvol( y, sqrtdelta );
	printf("s_max_bvol gives %.21g\n", temp2);
	printf("diff = %g\n", temp - temp2);
	ROUND_NEAR;
	dout = bvol( sy );
	printf("bvol gives %.21g\n", dout);
	diff = temp - dout;
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
#endif

#if MINMAX_A
	temp = min_a(  y );
	printf("min_a gives %.21g\n", temp);
	ROUND_NEAR;
	dout = afunc( sy );
	printf("afunc gives %.21g\n", dout);
	diff = dout - temp;
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
	
	temp = max_a(  y );
	printf("max_a gives %.21g\n", temp);
	ROUND_NEAR;
	dout = afunc( sy );
	printf("afunc gives %.21g\n", dout);
	diff = temp - dout;
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
#endif

#if ROUGH
	temp = rough_max_solid(  y );
	printf("rough_max_solid gives %.21g\n", temp);
	temp2 = max_solid( y, sqrtdelta );
	printf("max_solid gives %.21g\n", temp2);
	printf("diff = %g\n", temp - temp2);
	ROUND_NEAR;
	dout = solid( sy );
	printf("solid gives %.21g\n", dout);
	diff = temp - dout;
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
#endif

	i_solid(  y, sqrtdelta, out );
	ROUND_NEAR;
	printf("i_solid gives (%.21g, %.21g)\n", out[0], out[1]);
	printf("diff = %.21g\n", out[1]-out[0]);
	/*
	dout = solid( sy );
	printf("solid gives %.21g\n", dout);
	diff = fabs(2.0*dout - out[0] - out[1]);
	printf("diff = %.21g\n\n", diff);
	*/

#if ROUGH
	temp = rough_min_tvol(  y );
	printf("rough_min_tvol gives %.21g\n", temp);
	ROUND_NEAR;
	dout = tetvolume( sy );
	printf("tetvolume gives %.21g\n", dout);
	diff = dout - temp;
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
	
	temp = rough_max_tvol(  y );
	printf("rough_max_tvol gives %.21g\n", temp);
	ROUND_NEAR;
	dout = tetvolume( sy );
	printf("tetvolume gives %.21g\n", dout);
	diff = temp - dout;
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
#endif

#if ROUGH
	temp = rough_min_gma(  y );
	printf("rough_min_gma gives %.21g\n", temp);
	temp2 = s_min_gma( y, sqrtdelta );
	printf("s_min_gma gives %.21g\n", temp2);
	printf("diff = %g\n", temp2 - temp);
	ROUND_NEAR;
	dout = gma( sy );
	printf("gma gives %.21g\n", dout);
	diff = dout - temp2;
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));

	temp = rough_max_gma(  y );
	printf("rough_max_gma gives %.21g\n", temp);
	temp2 = s_max_gma( y, sqrtdelta );
	printf("s_max_gma gives %.21g\n", temp2);
	printf("diff = %g\n", temp - temp2);
	ROUND_NEAR;
	dout = gma( sy );
	printf("gma gives %.21g\n", dout);
	diff = temp2 - dout;
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
#endif

	temp = s_min_gma( y, sqrtdelta );
	temp2 = s_max_gma( y, sqrtdelta );
	ROUND_NEAR;
	printf("s_gma gives [%.21g, %.21g]\n", temp, temp2);
	diff = temp2 - temp;
	printf("diff = %.21g\n", diff);

#if DIH
	temp = s_max_dih(  x );
	printf("s_max_dih gives %.21g\n", temp);
	ROUND_NEAR;
	dout = dih( sy );
	printf("dih gives %.21g\n", dout);
	diff = temp - dout;
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
	
	temp = s_min_dih(  x );
	printf("s_min_dih gives %.21g\n", temp);
	ROUND_NEAR;
	dout = dih( sy );
	printf("dih gives %.21g\n", dout);
	diff = dout - temp;
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
#endif

#if MAXMIN_U
	temp = s_min_u(  y );
	printf("s_min_u gives %.21g\n", temp);
	ROUND_NEAR;
	dout = tomsu( sy );
	printf("tomsu gives %.21g\n", dout);
	diff = dout - temp;
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));

	temp = s_max_u(  y );
	printf("s_max_u gives %.21g\n", temp);
	ROUND_NEAR;
	dout = tomsu( sy );
	printf("tomsu gives %.21g\n", dout);
	diff = temp - dout;
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
#endif

#if QR_CRAD
	temp = min_qr_crad( y );
	printf("min_qr_crad gives %.21g\n", temp);
	ROUND_NEAR;
	dout = circumradius( sy );
	printf("circumradius gives %.21g\n", dout);
	diff = dout - temp;
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));

	temp = max_qr_crad(  y );
	printf("max_qr_crad gives %.21g\n", temp);
	ROUND_NEAR;
	dout = circumradius( sy );
	printf("circumradius gives %.21g\n", dout);
	diff = temp - dout;
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
#endif

	/*
	i_dbvol(  y,  out );
	i_bigdelta(  x,  out );
	i_istet(  y );
	i_tomsrho(  x,  out );
	*/
	
#if CRAD
	i_crad(  y,  out );
	printf("i_crad gives (%.21g, %.21g)\n", out[0], out[1]);
	printf("out[1]-out[0] = %.21g\n", out[1]-out[0]);
	ROUND_NEAR;
	dout = circumradius( sy );
	printf("crad gives %.21g\n", dout);
	diff = fabs(2.0*dout - out[0] - out[1]);
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
#endif

#if SCRAD
	c126 = s_crad2_6(  x,  out3 );
	printf("c126 = %d\n", c126 );
	i_sqrt( out3, out2 );
	printf("s_crad2_6 gives (%.21g, %.21g)\n", out3[0], out3[1]);
	printf("i_sqrt gives  (%.21g, %.21g)\n", out2[0], out2[1]);
	printf("out[1]-out2[1] = %.21g\n", out[1]-out2[1]);
	printf("out2[0]-out[0] = %.21g\n", out2[0]-out[0]);
	ROUND_NEAR;
	dout = circumradius( sy );
	printf("crad gives %.21g\n", dout);
	diff = fabs(2.0*dout - out2[0] - out2[1]);
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));

	s_crad2_1(  x,  out3, finfo );
	c126 = finfo[0];
	c135 = finfo[1];
	printf("c126 = %d\t", c126 );
	printf("c135 = %d\n", c135 );
	i_sqrt( out3, out2 );
	printf("s_crad2_1 gives (%.21g, %.21g)\n", out3[0], out3[1]);
	printf("i_sqrt gives  (%.21g, %.21g)\n", out2[0], out2[1]);
	printf("out[1]-out2[1] = %.21g\n", out[1]-out2[1]);
	printf("out2[0]-out[0] = %.21g\n", out2[0]-out[0]);
	ROUND_NEAR;
	dout = circumradius( sy );
	printf("crad gives %.21g\n", dout);
	diff = fabs(2.0*dout - out2[0] - out2[1]);
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
#endif
	/*
	i_tomsu(  x,  out );
	i_tomsv(  x,  out );
	i_auxPfun(  x,  out );
	i_tomsP(  x,  out );
	i_voronoivol(  y,  out );
	*/
	
	/*
	i_vor(  y, x, sqrtdelta, out );
	ROUND_NEAR;
	printf("i_vor gives (%.21g, %.21g)\n", out[0], out[1]);
	printf("diff = %.21g\n", out[1]-out[0]);
	ROUND_NEAR;
	*/

	best_i_vor_1(  y, x, sqrtdelta, out );
	ROUND_NEAR;
	printf("best_i_vor_1 gives (%.21g, %.21g)\n", out[0], out[1]);
	printf("diff = %.21g\n", out[1]-out[0]);
	ROUND_NEAR;

	s_dih(  x, out );
	ROUND_NEAR;
	printf("s_dih gives (%.21g, %.21g)\n", out[0], out[1]);
	printf("diff = %.21g\n", out[1]-out[0]);
	ROUND_NEAR;
	
/*
	dout = vor( sy );
	printf("vor gives %.21g\n", dout);
	diff = fabs(2.0*dout - out[0] - out[1]);
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
*/

#if VOR
	s_density_vor_6(  y, x, out2 );
	ROUND_NEAR;
	printf("s_density_vor_6 gives (%.21g, %.21g)\n", 
		out2[0], out2[1]);
	printf("out2[1]-out2[0] = %.21g\n", out2[1]-out2[0]);
/*
	ROUND_NEAR;
	printf("out[1]-out2[1] = %.21g\n", out[1]-out2[1]);
	printf("out2[0]-out[0] = %.21g\n", out2[0]-out[0]);
	dout = vor( sy );
	printf("vor gives %.21g\n", dout);
	diff = fabs(2.0*dout - out2[0] - out2[1]);
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
*/
#endif

/*
	s_density_vor_1(  y, x, out2 );
	ROUND_NEAR;
	printf("s_density_vor_1 gives (%.21g, %.21g)\n", 
		out2[0], out2[1]);
	printf("diff = %.21g\n\n", out2[1]-out2[0]);
*/

/*
	ROUND_NEAR;
	printf("out[1]-out2[1] = %.21g\n", out[1]-out2[1]);
	printf("out2[0]-out[0] = %.21g\n", out2[0]-out[0]);
	dout = vor( sy );
	printf("vor gives %.21g\n", dout);
	diff = fabs(2.0*dout - out2[0] - out2[1]);
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
*/

/*
	printf("scoring_system gives:\t%d\n", scoring_system( x ) );
	i_score(  y, x, sqrtdelta, out );
	printf("i_score gives (%.21g, %.21g)\n", out[0], out[1]);
	printf("out[1]-out[0] = %.21g\n", out[1]-out[0]);
	ROUND_NEAR;
	dout = score( sy );
	printf("score gives %.21g\n", dout);
	diff = fabs(2.0*dout - out[0] - out[1]);
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
*/

/*
	printf("octa_scoring gives:\t%d\n", octa_scoring( x ) );
	s_octa_vor(  y, x, out );
	printf("s_octa_vor gives (%.21g, %.21g)\n", out[0], out[1]);
	printf("out[1]-out[0] = %.21g\n", out[1]-out[0]);
	ROUND_NEAR;
	dout = score( sy );
	printf("score gives %.21g\n", dout);
	diff = fabs(2.0*dout - out[0] - out[1]);
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
*/

/*
	out2[0] = max_solid( y, sqrtdelta );
	out2[1] = out[1] + 0.396*out2[0] - 0.5*0.5383016645753285;
	printf("sph = %.21g,  sc = %.21g\n", out2[0], out[1]);
	printf("relation test:  %.21g\n\n", out2[1]);
*/

#if CRAD3LEN
	i_crad3len(  y, out );
	printf("i_crad3len gives (%.21g, %.21g)\n", out[0], out[1]);
	printf("out[1]-out[0] = %.21g\n", out[1]-out[0]);
	ROUND_NEAR;
	dout = crad3len( sy );
	printf("crad3len gives %.21g\n", dout);
	diff = fabs(2.0*dout - out[0] - out[1]);
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
	
	s_crad3x(  x, out );
	printf("s_crad3x gives (%.21g, %.21g)\n", out[0], out[1]);
	printf("out[1]-out[0] = %.21g\n", out[1]-out[0]);
	ROUND_NEAR;
	dout = crad3len( sy );
	printf("crad3len gives %.21g\n", dout);
	diff = fabs(2.0*dout - out[0] - out[1]);
	printf("diff = %.21g\n", diff);
	printf("relative diff:  %.21g\n\n", diff/fabs(dout));
#endif
	
	/*
	i_voronoivol(x, sqrtdelta, out );
	ROUND_NEAR;
	printf("i_voronoivol gives (%.21g, %.21g)\n", out[0], out[1]);
	printf("diff = %.21g\n", out[1]-out[0]);
	i_vor_trunc( y, x, sqrtdelta, out );
	ROUND_NEAR;
	printf("i_vor_trunc gives (%.21g, %.21g)\n", out[0], out[1]);
	printf("diff = %.21g\n", out[1]-out[0]);
	best_i_vor_6(  y, x, sqrtdelta, out );
	ROUND_NEAR;
	printf("best_i_vor_6 gives (%.21g, %.21g)\n", out[0], out[1]);
	printf("diff = %.21g\n", out[1]-out[0]);
	ROUND_NEAR;
	*/
	
	ROUND_NEAR;
	printf("\n");
}


void get_comp_trash( void )
{
	double y[12], x[12], sy[7], diff, eps, temp, temp2;
	double out[2], out_i[12], outvals[10];
	int i, j, rel[5];
	

	i_init();
	
	ROUND_NEAR;
	eps = 0.0;
	
	printf("Enter edge lengths:\n");
	
	for( i=1; i<7; i++ )
		scanf("%lf", sy + i);
	/*
	printf("Got the following:  \n");
	for( i=1; i<7; i++ )
		printf("%f\t", sy[i]);
	printf("\n\n");
	*/
	printf("Enter epsilon (interval width):  ");
	scanf("%lf", &eps);
	
	ROUND_UP;
	for( i=0; i<6; i++ ) {
		y[2*i] = sy[i+1];
		y[2*i+1] = sy[i+1] + eps;
		}

	ROUND_NEAR;
	printf("Enter individual values vector (5):  \n");
	for( i=0; i<5; i++ )
		scanf("%d", rel + i);
	printf("Got the following:  \n");
	for( i=0; i<5; i++ )
		printf("%d\t", rel[i]);
	printf("\n");
	printf("Enter relation vector (5):  \n");
	for( i=0; i<5; i++ )
		scanf("%lf", sy + i);
	printf("Got the following:  \n");
	for( i=0; i<5; i++ )
		printf("%g\t", sy[i]);
	printf("\n");

	/* Compute square of edge lengths */
	ROUND_UP;
	for( i=1; i<12; i+=2 )
		x[i] = y[i]*y[i];
	ROUND_DOWN;
	for( i=0; i<12; i+=2 )
		x[i] = y[i]*y[i];
	
	t_composite( rel, sy, y, x, out, out_i, outvals );
	ROUND_NEAR;
	printf("comp = [%.18g, \t%.18g]\n", out[0], out[1] );
	printf("comp_i = \n");
	for( i=0; i<12; i+=2 ) {
		printf("[%.18g, \t%.18g]\n", out_i[i], out_i[i+1] );
	}
	printf("\n");
	for( i=0; i<5; i++ ) {
		if( rel[i] ) {
			printf("outval[%d] = [%.18g, \t%.18g]\n",
				i, outvals[2*i], outvals[2*i+1] );
		}
	}
	ROUND_NEAR;
}


void get_new_trash( void )
{
	double y[12], x[12], out[2], sy[7], dout, diff, eps;
	double delta[2], sqrtdelta[2], oosqrtdelta[2];
	double temp[2], partials[12];
	int i;
	
	i_init();
	
	while( 1 ) {
		ROUND_NEAR;
		
		printf("Enter edge lengths:  ");
		for( i=1; i<7; i++ )
			scanf("%lf", sy + i);
		printf("Got the following:  \n");
		for( i=1; i<7; i++ )
			printf("%f\t", sy[i]);
		printf("\n\n");
		printf("Enter epsilon (interval width):  ");
		scanf("%lf", &eps);
		
		ROUND_UP;
		for( i=0; i<6; i++ ) {
			y[2*i] = sy[i+1];
			y[2*i+1] = sy[i+1] + eps;
			}
		ROUND_NEAR;

		for( i=0; i<12; i+=2 )
			printf("%f\t", y[i]);
		printf("\n");
		for( i=1; i<12; i+=2 )
			printf("%f\t", y[i]);
		printf("\n\n");

		/* Compute square of edge lengths */
		ROUND_UP;
		for( i=1; i<12; i+=2 )
			x[i] = y[i]*y[i];
		ROUND_DOWN;
		for( i=0; i<12; i+=2 )
			x[i] = y[i]*y[i];
		
		i_bigdelta( x, delta );
		ROUND_NEAR;
		printf("i_bigdelta gives      (%.21g, %.21g)\n", delta[0], delta[1]);
		printf("out[1]-out[0] = %.21g\n", delta[1]-delta[0]);
		i_bigdelta_best( x, delta );
		ROUND_NEAR;
		printf("i_bigdelta_best gives (%.21g, %.21g)\n", delta[0], delta[1]);
		printf("out[1]-out[0] = %.21g\n", delta[1]-delta[0]);

		i_sqrt( delta, sqrtdelta );
		temp[0] = 1.0;
		temp[1] = 1.0;
		i_div( temp, sqrtdelta, oosqrtdelta );
		
		
		i_dih_old(  x,  out );
		ROUND_NEAR;
		printf("i_dih_old gives  (%.21g, %.21g)\n", out[0], out[1]);
		printf("out[1]-out[0] = %.21g\n", out[1]-out[0]);
		
		i_dih(  x,  out );
		ROUND_NEAR;
		printf("i_dih gives      (%.21g, %.21g)\n", out[0], out[1]);
		printf("out[1]-out[0] = %.21g\n", out[1]-out[0]);
		
		i_dih_alt(  x, delta, out );
		ROUND_NEAR;
		printf("i_dih_alt gives  (%.21g, %.21g)\n", out[0], out[1]);
		printf("out[1]-out[0] = %.21g\n", out[1]-out[0]);
		
		/*
		s_dih_partials( x, y, partials );
		ROUND_NEAR;
		printf("s_dih_partials:  \n");
		for( i=0; i<12; i+=2 )
			printf("[%20.17g, %20.17g]\n", partials[i], partials[i+1]);
		*/
		
		i_dih_partials( x, y, partials );
		/*
		ROUND_NEAR;
		printf("i_dih_partials:  \n");
		for( i=0; i<12; i+=2 )
			printf("[%20.17g, %20.17g]\n", partials[i], partials[i+1]);
		*/

		i_dih_best(  x, partials, out );
		ROUND_NEAR;
		printf("i_dih_best gives (%.21g, %.21g)\n", out[0], out[1]);
		printf("out[1]-out[0] = %.21g\n", out[1]-out[0]);
		
		ROUND_NEAR;
	}
}

/* 2.1 2.2 2.3 2.01 2.02 2.03 */
/* 2.123 2.29102 2.341248 2.011234 2.023309 2.03982134 */
/* 2.0123 2.029102 2.0341248 2.011234 2.023309 2.73982134 */
/* 2.0123 2.029102 2.0341248 2.011234 2.023309 2.63982134 */
/* 2.423 2.439102 2.0541248 2.001234 2.013309 2.74982134 */
/* 2.0 2.0 2.0 2.0 2.0 2.828427124746189847 */
/* 2.0 2.0 2.0 2.0 2.0 2.828427124710598761 */
/* 2.0 2.0 2.0 2.0 2.0 2.828427124746190291 */
/* 2.828427124746189847 2.0 2.0 2.0 2.0 2.0 */

/* 2.029102 2.0341248 2.011234 2.023309 2.0123 2.83982134 */

/* 2.73982134 2.029102 2.0341248 2.011234 2.023309 2.0123 */
/* 2.63982134 2.029102 2.0341248 2.011234 2.023309 2.0123 */
/* 2.74982134 2.439102 2.0541248 2.001234 2.013309 2.423 */
/* 2.74982134 2.039102 2.4541248 2.001234 2.313309 2.0423 */
/* 2.74982134 2.439102 2.4541248 2.001234 2.313309 2.423 */
/* 2.83982134 2.029102 2.0341248 2.011234 2.023309 2.0123 */

/* 2.1 2.2 2.3 2.01 2.02 2.03
2.123 2.29102 2.341248 2.011234 2.023309 2.03982134 */
/*
2.000000000 2.000000000 2.000000000 2.001992188 2.000000000 2.001992188 
2.000000000 2.000000000 2.255000000 2.256992188 2.801062294 2.802306150 
*/
/*
2.0 2.0 2.0 2.0 2.255000000 2.801062294
*/
/* 2.0 2.0 2.079189 2.007949 2.145492 2.51 */
/*
2.000000    2.000001
2.000000    2.000000
2.079189    2.079190
2.007949    2.007950
2.145492    2.145493
2.510000    2.510000
*/
/* 2 2 2 2 2.3084 2.787982 */
/*
2.828427124710598761    2.828427124746189847
2.000000000000000000    2.000000002424120904
2.000000000000000000    2.000000002424120904
2.000000000000000000    2.000000002424120904
2.000000000000000000    2.000000002424120904
2.000000000000000000    2.000000002424120904
*/
/* 2.7578125 2.125 2.0625 2.03125 2.015625 2.0078125 */

/*
From: "Tom Hales" <hales>
Date: Thu, 20 Feb 1997 12:42:31 -0500 (EST)
Message-Id: <199702201742.MAA11198@orthogonal.math.lsa.umich.edu>
To: samf@math.lsa.umich.edu
Subject: truncated Vor
Cc: hales@math.lsa.umich.edu
Content-Length: 221


If y = {2.1,2.2,2.3,2.8,2.4,2.5},
then the Voronoi score, truncated at Sqrt[2], gives me -0.4062055021818615.

The analytic Voronoi score gives me -0.4112805590580618.

Is this consistent with your functions?

Best,
Tom

2.1 2.2 2.3 2.8 2.4 2.5
*/