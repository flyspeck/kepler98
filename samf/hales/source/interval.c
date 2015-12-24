/* interval.c, by Samuel Ferguson, (c) 1997. */
/* This is a collection of interval arithmetic routines, 
designed for verification of sphere-packing conjectures.  */

#include "system_headers.h"
#include "macros.h"
#include "interval.h"

/* #pragma fenv_access	*/	/* Tell the compiler not to do anything too clever.  */

/* Need this for i_acos().  We require that OHSQRT2 = 1/Sqrt[2] - epsilon. */

#define OHSQRT2		0.707106781
#define ONE_3		0.33333333333333333333333
#define ONE_5		0.2
#define AQ_CUT		9.765625e-4
/*					123456789012345678901234567890	*/

/* Global variables */

double i_pi_const[2];
double i_pi_2_const[2];
double i_doct_const[2];
double i_two_pi_5_const[2];
double i_sqrt2[2];

/* Local prototypes */
void euclidean( double x0, double y0, double out[2] );


/* Routines */
void i_init( void )
{
	char pi_str[255] =   "3.14159265358979323846264338328";
	char doct_str[255] = "0.7209029495174650928";
	char two_pi_5_str[255] = "1.2566370614359172953850574";
	/*	decimal dec;	*/ /* Can't figure out how to use str2dec, etc.  No docs. */
	
	/*	str2dec( pi_str, */
	ROUND_DOWN;
	i_pi_const[0] = atof( pi_str );
	i_pi_2_const[0] = i_pi_const[0]*0.5;
	i_doct_const[0] = atof( doct_str );
	i_two_pi_5_const[0] = atof( two_pi_5_str );
	i_sqrt2[0] = sqrt( 2.0 );
	ROUND_UP;
	i_pi_const[1] = atof( pi_str );
	i_pi_2_const[1] = i_pi_const[1]*0.5;
	i_doct_const[1] = atof( doct_str );
	i_two_pi_5_const[1] = atof( two_pi_5_str );
	i_sqrt2[1] = sqrt( 2.0 );
	/*
	ROUND_NEAR;
	printf("i_pi_const:  %g\n", i_pi_const[1] - i_pi_const[0]);
	printf("i_doct_const:  %g\n", i_doct_const[1] - i_doct_const[0]);
	printf("i_two_pi_5_const:  %g\n", i_two_pi_5_const[1] - i_two_pi_5_const[0]);
	printf("Bitwise OR of floating-point exception macros: %#x\n", 
		fetestexcept( FE_ALL_EXCEPT ) );
	*/
}


void i_add( double x[2], double y[2], double out[2] )
{	
	ROUND_DOWN;
	out[0] = x[0] + y[0];
	ROUND_UP;
	out[1] = x[1] + y[1];
}


void i_sub( double x[2], double y[2], double out[2] )
{	
	ROUND_DOWN;
	out[0] = x[0] - y[1];
	ROUND_UP;
	out[1] = x[1] - y[0];
}


void i_negate( double x[2], double out[2] )
{
	out[0] = -x[1];
	out[1] = -x[0];
}


void i_smult( double x, double y[2], double out[2] )
{
	if( x >= 0.0 ) {
		ROUND_DOWN;
		out[0] = x*y[0];
		ROUND_UP;
		out[1] = x*y[1];
		}
	else {
		ROUND_DOWN;
		out[0] = x*y[1];
		ROUND_UP;
		out[1] = x*y[0];
		}
}


void i_mult( double x[2], double y[2], double out[2] )
{	
	if( x[0] >= 0.0 ) {
		if( y[0] >= 0.0 ) {
			ROUND_DOWN;
			out[0] = x[0]*y[0];
			ROUND_UP;
			out[1] = x[1]*y[1];
			}
		else {
			ROUND_DOWN;
			out[0] = y[0]*x[1];
			ROUND_UP;
			if( y[1] >= 0.0 ) {
				out[1] = x[1]*y[1];
				}
			else {
				out[1] = x[0]*y[1];
				}
			}
		}
	else if( y[0] >= 0.0 ) {
	/* assume x[0] < 0.0 */
		ROUND_DOWN;
		out[0] = x[0]*y[1];
		ROUND_UP;
		if( x[1] >= 0.0 ) {
			out[1] = x[1]*y[1];
			}
		else {
			out[1] = y[0]*x[1];
			}
		}
	else if( x[1] < 0.0 ) {
	/* assume x[0], y[0] < 0.0 */
		ROUND_UP;
		out[1] = x[0]*y[0];
		ROUND_DOWN;
		if( y[1] >= 0.0 ) {
			out[0] = x[0]*y[1];
			}
		else {
			out[0] = x[1]*y[1];
			}
		}
	else if( y[1] < 0.0 ) {
	/* assume x[1] >= 0.0 */
		ROUND_UP;
		out[1] = x[0]*y[0];
		ROUND_DOWN;
		out[0] = y[0]*x[1];
		}
	else {
	/* assume x[0], y[0] < 0.0 and x[1], y[1] > 0.0 */
		ROUND_DOWN;
		out[0] = MIN( x[0]*y[1], y[0]*x[1] );
		ROUND_UP;
		out[1] = MAX( x[0]*y[0], x[1]*y[1] );
		}
}


void i_div( double x[2], double y[2], double out[2] )
{
	if( y[0] > 0.0 ) {
		if( x[0] > 0.0 ) {
			ROUND_DOWN;
			out[0] = x[0]/y[1];
			ROUND_UP;
			out[1] = x[1]/y[0];
			}
		else if( x[1] < 0.0 ) {
			ROUND_DOWN;
			out[0] = x[0]/y[0];
			ROUND_UP;
			out[1] = x[1]/y[1];
			}
		else {
			ROUND_DOWN;
			out[0] = x[0]/y[0];
			ROUND_UP;
			out[1] = x[1]/y[0];
			}
		}
	else if( y[1] < 0.0 ) {
		if( x[0] > 0.0 ) {
			ROUND_DOWN;
			out[0] = x[1]/y[0];
			ROUND_UP;
			out[1] = x[0]/y[1];
			}
		else if( x[1] < 0.0 ) {
			ROUND_DOWN;
			out[0] = x[1]/y[1];
			ROUND_UP;
			out[1] = x[0]/y[0];
			}
		else {
			ROUND_DOWN;
			out[0] = x[1]/y[0];
			ROUND_UP;
			out[1] = x[0]/y[0];
			}
		}
	else {	/* [y[0],y[1]] contains zero */
		set_infinity( out );
		}
}


void i_sqrt( double x[2], double out[2] )
{
	ROUND_DOWN;
	if( x[0] < 0.0 )
		out[0] = 0.0;
	else
		out[0] = sqrt( x[0] );
	ROUND_UP;
	out[1] = sqrt( x[1] );
}


void i_atan( double x[2], double out[2] )
{
	ROUND_DOWN;
	out[0] = atan( x[0] ) - ATANERR;
	ROUND_UP;
	out[1] = atan( x[1] ) + ATANERR;
}


/* Assume x >= 0 */
/* Error analysis seems pretty simple. */

double max_atanquot( double x )
{
	double temp, x2, y;
	
	ROUND_UP;
	if( x > AQ_CUT ) {
		y = (atan( x ) + ATANERR)/x;
	}
	else {	/* Use truncated series */
		temp = ONE_3*x*x;
		x2 = x*x;
		y = ONE_5*x2*x2 + ATANERR;
		y -= temp;
		y += 1.0;
	}
	return( y );
}


double min_atanquot( double x )
{
	double temp, x2, y;
	
	ROUND_DOWN;
	if( x > AQ_CUT ) {
		y = (atan( x ) - ATANERR)/x;
	}
	else {	/* Use truncated series */
		temp = ONE_3*x*x;
		x2 = x*x;
		y = ONE_5*x2*x2 - ATANERR;
		y -= temp;
		y += 1.0;
	}
	return( y );
}


/* Damn, but this is *ugly*.  Would probably make
sense to just rewrite the functions so they don't
need to call acos(), use atan() instead.  */

void i_acos( double x[2], double out[2] )
{
	double xv, y, t;
	
	xv = x[0];
	if( xv >= 1.0 )
		out[1] = 0.0;
	else if( xv > OHSQRT2 ) {
#if GOT_FMADD
		ROUND_DOWN;
		t = xv*xv - 1.0;
		ROUND_UP;
		y = sqrt( -t );
#else
		ROUND_DOWN;
		t = xv*xv;
		ROUND_UP;
		y = sqrt( 1.0 - t );
#endif
		out[1] = atan( y/xv ) + ATANERR;
		}
	else if( xv > - OHSQRT2 ) {
		if( xv >= 0.0 ) {
#if GOT_FMADD
			ROUND_UP;
			t = xv*xv - 1.0;
			ROUND_DOWN;
			y = sqrt( -t );
#else
			ROUND_UP;
			t = xv*xv;
			ROUND_DOWN;
			y = sqrt( 1.0 - t );
#endif
			t = atan( xv/y );
			ROUND_UP;
			out[1] = PI_2_HI - t + ATANERR;
			}
		else {
#if GOT_FMADD
			ROUND_DOWN;
			t = xv*xv - 1.0;
			ROUND_UP;
			y = sqrt( -t );
#else
			ROUND_DOWN;
			t = xv*xv;
			ROUND_UP;
			y = sqrt( 1.0 - t );
#endif
			t = atan( (-xv)/y );
			out[1] = PI_2_HI + t + ATANERR;
			}
		}
	else if ( xv > -1.0 ) {
#if GOT_FMADD
		ROUND_UP;
		t = xv*xv - 1.0;
		ROUND_DOWN;
		y = sqrt( -t );
#else
		ROUND_UP;
		t = xv*xv;
		ROUND_DOWN;
		y = sqrt( 1.0 - t );
#endif
		t = atan( y/(-xv) );
		ROUND_UP;
		out[1] = PI_HI - t + ATANERR;
		}
	else if( ISNAN( xv ) )
		out[1] = xv;
	else
		out[1] = PI_HI;
	
	/* Now do it again for the upper bound, switch UP and DOWN.  */
	xv = x[1];
	if( xv >= 1.0 )
		out[0] = 0.0;
	else if( xv > OHSQRT2 ) {
#if GOT_FMADD
		ROUND_UP;
		t = xv*xv - 1.0;
		ROUND_DOWN;
		y = sqrt( -t );
#else
		ROUND_UP;
		t = xv*xv;
		ROUND_DOWN;
		y = sqrt( 1.0 - t );
#endif
		out[0] = atan( y/xv ) - ATANERR;
		}
	else if( xv > - OHSQRT2 ) {
		if( xv >= 0.0 ) {
#if GOT_FMADD
			ROUND_DOWN;
			t = xv*xv - 1.0;
			ROUND_UP;
			y = sqrt( -t );
#else
			ROUND_DOWN;
			t = xv*xv;
			ROUND_UP;
			y = sqrt( 1.0 - t );
#endif
			t = atan( xv/y );
			ROUND_DOWN;
			out[0] = PI_2_LO - t - ATANERR;
			}
		else {
#if GOT_FMADD
			ROUND_UP;
			t = xv*xv - 1.0;
			ROUND_DOWN;
			y = sqrt( -t );
#else
			ROUND_UP;
			t = xv*xv;
			ROUND_DOWN;
			y = sqrt( 1.0 - t );
#endif
			t = atan( (-xv)/y );
			out[0] = PI_2_LO + t - ATANERR;
			}
		}
	else if ( xv > -1.0 ) {
#if GOT_FMADD
		ROUND_DOWN;
		t = xv*xv - 1.0;
		ROUND_UP;
		y = sqrt( -t );
#else
		ROUND_DOWN;
		t = xv*xv;
		ROUND_UP;
		y = sqrt( 1.0 - t );
#endif
		t = atan( y/(-xv) );
		ROUND_DOWN;
		out[0] = PI_LO - t - ATANERR;
		}
	else if( ISNAN( xv ) )
		out[0] = xv;
	else
		out[0] = PI_LO;
}


/* We assume that len[0] > 0  */
void i_squarelen( double len[2], double out[2] )
{
	ROUND_DOWN;
	out[0] = len[0]*len[0];
	ROUND_UP;
	out[1] = len[1]*len[1];
}


void set_infinity( double out[2] )
{
#if DA_SYSTEM == 1
	out[0] = -INFINITY;
	out[1] = INFINITY;
#elif DA_SYSTEM == 2
	out[0] = -HUGE_VAL;
	out[1] = HUGE_VAL;
#elif DA_SYSTEM == 3
	out[0] = -__infinity;
	out[1] = __infinity;
#else
	out[0] = -HUGE_VAL;
	out[1] = HUGE_VAL;
#endif
}


void i_old_recognize( double num, double out[2] )
{
	int i, n, done;
	double x, y, z, bound[2], m;
	
	for( i=0; i<2; i++ ) {
		ROUND_NEAR;
		x = fabs( num );
		done = 0;
		n = 0;
		while( !done ) {
			z = pow( 10.0, (double) n );
			m = floor( z*x );
			/*
			printf("z*x = %40.30f\n", z*x);
			printf("m   = %40.30f\n", m);
			*/
			if( (z*x - m)/m < 1.0e-12 ) {
				done = 1;
			} else {
				n++;
				if( n > 12 ) {
					printf("Went too far in i_recognize.\n");
					done = 1;
				}
			}
		}
		y = m;
		if( i==0 ) {
			ROUND_DOWN;
		} else {
			ROUND_UP;
		}
		bound[i] = y/z;
	}
	
	if( num > 0.0 ) {
		out[0] = bound[0];
		out[1] = bound[1];
	} else {
		out[0] = -bound[1];
		out[1] = -bound[0];
	}
	ROUND_NEAR;
	if( fabs( out[0] - num ) > 1.0e-14 || 
			fabs( out[1] - num ) > 1.0e-14 ) {
		printf("Problem with bound in i_recognize.\n");
		printf("%.18f in [%.18f, %.18f]?\n", 
			num, bound[0], bound[1]);
	}
}


void i_recognize( double num, double out[2] )
{
	double coeffs[2];
	
	if( fabs( num ) < 1e-14 ) {
		out[0] = -fabs( num );
		out[1] = fabs( num );
	} else {
		euclidean( 1.0, num, coeffs );
		ROUND_DOWN;
		out[0] = coeffs[0]/(-coeffs[1]);
		euclidean( 1.0, num, coeffs );
		ROUND_UP;
		out[1] = coeffs[0]/(-coeffs[1]);
	}
}


void euclidean( double x0, double y0, double out[2] )
{
	int done;
	double x, xp, y, yp, a;
	double p, pp, q, qp, r, rp, s, sp;
	
	ROUND_NEAR;
	xp = fabs( x0 );
	yp = fabs( y0 );
	if( yp < xp ) {
		x = yp;
		y = xp;

		p = 0.0;
		q = 1.0;
		r = 1.0;
		s = 0.0;
	} else {
		x = xp;
		y = yp;
		
		p = 1.0;
		q = 0.0;
		r = 0.0;
		s = 1.0;
	}
	
	if( x < 1.0e-8 ) {
		done = 1;
	} else {
		done = 0;
	}
	
	while( !done ) {
		a = floor( y/x );
		
		xp = y - a*x;
		yp = x;
		
		x = xp;
		y = yp;
		
		pp = floor( r - a*p + 0.5 );
		qp = floor( s - a*q + 0.5 );
		rp = p;
		sp = q;
		
		p = pp;
		q = qp;
		r = rp;
		s = sp;
		if( x < 1.0e-8 )
			done = 1;
	}
	if( y0 < 0.0 && x0 > 0.0 )
		q = -q;
	x = fabs( y0/x0 + p/q );
	if( x > 1.0e-12 ) {
		printf("Euclidean algorithm failed to converge\n");
		printf("on (%.18f, %.18f), giving (%f, %f),\n",
			x0, y0, p, q );
		printf("with a relative error of %g\n", x );
	}
	out[0] = p;
	out[1] = q;
}



