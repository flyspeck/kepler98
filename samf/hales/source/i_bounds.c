/* i_bounds.c,  by Samuel Ferguson, (c) 1997,
containing routines for computing bounds on 
functions contained in i_sphere.c  */

#include "system_headers.h"
#include "i_sphere.h"
#include "interval.h"
#include "i_bounds.h"
#include "macros.h"


/* #pragma fenv_access	*/	/* Tell the compiler not to do 
anything too clever.  */


#define SQUARE( x )		( x )*( x )

/* Need this for acos().  We require that 
OHSQRT2 = 1/Sqrt[2] - epsilon. */
#define OHSQRT2		0.707106781


/* External variables */

extern double i_pi_const[2];
extern double i_pi_2_const[2];
extern double i_doct_const[2];
extern double i_two_pi_5_const[2];


/* Routines */

/* rough_min_delta() gives a correct, but rough, lower
bound on bigdelta(), using the sign trick.  (No information
about the derivatives is used.)  */

double rough_min_delta( double x[12] )
{
	double pterms, nterms, p1, p2, p3, p4;

	ROUND_UP;
	p1 = x[1]*x[7]*(x[1] + x[7]);
	p2 = x[3]*x[9]*(x[3] + x[9]);
	p3 = x[5]*x[11]*(x[5] + x[11]);
	p4 = x[5]*(x[3]*x[7] + x[1]*x[9]) + x[11]*(x[1]*x[3] + 
		x[7]*x[9]);
	nterms = p1 + p2 + p3 + p4;
	
	ROUND_DOWN;
	p1 = x[0]*x[6]*(x[2] + x[4] + x[8] + x[10]);
	p2 = x[2]*x[8]*(x[0] + x[4] + x[6] + x[10]);
	p3 = x[4]*x[10]*(x[0] + x[2] + x[6] + x[8]);
	pterms = p1 + p2 + p3;
	
	p1 = pterms - nterms;
	return( p1 );
}


double rough_max_delta( double x[12] )
{
	double pterms, nterms, p1, p2, p3, p4;

	ROUND_DOWN;
	p1 = x[0]*x[6]*(x[0] + x[6]);
	p2 = x[2]*x[8]*(x[2] + x[8]);
	p3 = x[4]*x[10]*(x[4] + x[10]);
	p4 = x[4]*(x[2]*x[6] + x[0]*x[8]) + x[10]*(x[0]*x[2] + 
		x[6]*x[8]);
	nterms = p1 + p2 + p3 + p4;
	
	ROUND_UP;
	p1 = x[1]*x[7]*(x[3] + x[5] + x[9] + x[11]);
	p2 = x[3]*x[9]*(x[1] + x[5] + x[7] + x[11]);
	p3 = x[5]*x[11]*(x[1] + x[3] + x[7] + x[9]);
	pterms = p1 + p2 + p3;
	
	p1 = pterms - nterms;
	return( p1 );
}


/* min_a( y ) = a( y[0], y[2], y[4], y[7], y[9], y[11] ) */

double min_a( double y[12] )
{
	double y2[12];
	double p1, p2, p3, y123;
	int i;

	ROUND_DOWN;
	for( i=0; i<6; i+=2 )
		y2[i] = y[i]*y[i];
	ROUND_UP;
	for( i=7; i<12; i+=2 )
		y2[i] = y[i]*y[i];
	ROUND_DOWN;
	p1 = y2[2] + y2[4] - y2[7];
	p2 = y2[0] + y2[4] - y2[9];
	p3 = y2[0] + y2[2] - y2[11];
	y123 = y[0]*y[2]*y[4];
	p1 = y123 + 0.5*(y[0]*p1 + y[2]*p2 + y[4]*p3);
	if( p1 < 0.0 )
		p1 = 0.0;
	return( p1 );
}


/* max_a( y ) = a( y[1], y[3], y[5], y[6], y[8], y[10] ) */

double max_a( double y[12] )
{
	double y2[12];
	double p1, p2, p3, y123;
	int i;

	ROUND_DOWN;
	for( i=1; i<6; i+=2 )
		y2[i] = y[i]*y[i];
	ROUND_UP;
	for( i=6; i<12; i+=2 )
		y2[i] = y[i]*y[i];
	ROUND_UP;
	p1 = y2[3] + y2[5] - y2[6];
	p2 = y2[1] + y2[5] - y2[8];
	p3 = y2[1] + y2[3] - y2[10];
	y123 = y[1]*y[3]*y[5];
	p1 = y123 + 0.5*(y[1]*p1 + y[3]*p2 + y[5]*p3);
	return( p1 );
}


/* s_min_u() will be correct for small simplices (the sign of the
derivatives depends on x_j < (2 Sqrt[2])^2, which is the case only
for small simplices).  Otherwise, use min_u().  */

double s_min_u( double x[6] )
{
	double pterms, nterms, val;
	
	ROUND_UP;
	nterms = x[0]*x[0] + x[2]*x[2] + x[4]*x[4];
	ROUND_DOWN;
	pterms = 2.0*(x[0]*(x[4] + x[2]) + x[2]*x[4]);
	val = pterms - nterms;
	if( val < 0.0 )
		val = 0.0;
	return( val );
}


/* Same caveat applies for s_max_u() as for s_min_u().  */

double s_max_u( double x[6] )
{
	double pterms, nterms, val;
	
	ROUND_DOWN;
	nterms = x[1]*x[1] + x[3]*x[3] + x[5]*x[5];
	ROUND_UP;
	pterms = 2.0*(x[1]*(x[5] + x[3]) + x[3]*x[5]);
	val = pterms - nterms;
	return( val );
}


/* Caveat:  no long ( > 2Sqrt[2] ) edges.  */
/* delta4( x1, x2_, x3_, x4^, x5_, x6_ ) */

double s_min_delta4( double x[12] )
{
	double temp, p1, p2, p3, val;
	
	ROUND_DOWN;
	temp = x[2] + x[4] - 2.0*x[7] + x[8] + x[10];
	p1 = x[0]*(-x[0] + temp);
	p2 = x[1]*(-x[1] + temp);
	val = MIN( p1, p2 );
	ROUND_UP;
	p2 = x[2]*x[4] + x[8]*x[10];	/* negative terms */
	ROUND_DOWN;
	p1 = x[2]*x[8] + x[4]*x[10];	/* positive terms */
	p3 = p1 - p2;
	val += p3;
	return( val );
}


/* Caveat:  no long ( > 2Sqrt[2] ) edges.  */
/* delta4( x1, x2^, x3^, x4_, x5^, x6^ ) */

double s_max_delta4( double x[12] )
{
	double max1, max2, temp2, p1, p2, p3, val, nterms, pterms;
	
	ROUND_DOWN;
	max1 = 0.5*(x[3] + x[5] - 2.0*x[6] + x[9] + x[11]);
	/* max1 is lower bound on potential maximum */
	ROUND_UP;
	temp2 = x[3] + x[5] - 2.0*x[6] + x[9] + x[11];
	max2 = 0.5*temp2;
	/* max2 is upper bound on potential maximum */
	if( max1 > x[1] || max2 < x[0] ) {
		p1 = x[0]*(-x[0] + temp2);
		p2 = x[1]*(-x[1] + temp2);
		val = MAX( p1, p2 );
		}
	else {
		ROUND_DOWN;
		nterms = max1*(max1 + 2.0*x[6]);
		ROUND_UP;
		pterms = max2*(x[3] + x[5] + x[9] + x[11]);
		val = pterms - nterms;
		}
	ROUND_DOWN;
	p2 = x[3]*x[5] + x[9]*x[11];
	ROUND_UP;
	p1 = x[3]*x[9] + x[5]*x[11];
	p3 = p1 - p2;
	val += p3;
	return( val );
}


void i_delta4( double x[12], double out[2] )
{
	double temp, p1[2], p2[2], p3[2];
	
	ROUND_DOWN;
	temp = -x[1] + x[2] + x[4] - 2.0*x[7] + x[8] + x[10];
	if( temp > 0.0 ) {
		p1[0] = x[0]*temp;
	} else {
		p1[0] = x[1]*temp;
	}
	p2[0] = x[2]*x[4] + x[8]*x[10];
	p3[0] = x[2]*x[8] + x[4]*x[10];
	ROUND_UP;
	temp = -x[0] + x[3] + x[5] - 2.0*x[6] + x[9] + x[11];
	if( temp > 0.0 ) {
		p1[1] = x[1]*temp;
	} else {
		p1[1] = x[0]*temp;
	}
	p2[1] = x[3]*x[5] + x[9]*x[11];
	p3[1] = x[3]*x[9] + x[5]*x[11];
	out[1] = p1[1] - p2[0] + p3[1];
	ROUND_DOWN;
	out[0] = p1[0] - p2[1] + p3[0];
}


double rough_min_tvol( double y[12] )
{
	double x[12], mindelta, temp;
	int i;
	
	ROUND_UP;
	for( i=1; i<12; i+=2 )
		x[i] = y[i]*y[i];
	ROUND_DOWN;
	for( i=0; i<12; i+=2 )
		x[i] = y[i]*y[i];
	
	mindelta = rough_min_delta( x );
	/*	ROUND_DOWN;	*/
	temp = sqrt( mindelta );
	temp = temp*ONE_12_LO;
	return( temp );
}


double rough_max_tvol( double y[12] )
{
	double x[12], maxdelta, temp;
	int i;
	
	ROUND_UP;
	for( i=1; i<12; i+=2 )
		x[i] = y[i]*y[i];
	ROUND_DOWN;
	for( i=0; i<12; i+=2 )
		x[i] = y[i]*y[i];
	
	maxdelta = rough_max_delta( x );
	/*	ROUND_UP;	*/
	temp = sqrt( maxdelta );
	temp = temp*ONE_12_HI;
	return( temp );
}


double rough_max_solid( double y[12] )
{
	double x[12], max, temp, mina;
	int i;
	
	ROUND_DOWN;
	for( i=0; i<12; i+=2 )
		x[i] = y[i]*y[i];
	ROUND_UP;
	for( i=1; i<12; i+=2 )
		x[i] = y[i]*y[i];
	
	max = rough_max_delta( x );
	/*	ROUND_UP;	*/
	temp = sqrt( max );
	mina = 2.0*min_a( y );
	ROUND_UP;
	max = temp/mina;
	temp = 2.0*(atan( max ) + ATANERR);
	return( temp );
}


double max_solid( double y[12], double sqrtdelta[2] )
{
	double max, temp, mina;
	
	temp = sqrtdelta[1];
	mina = 2.0*min_a( y );
	ROUND_UP;
	max = temp/mina;
	temp = 2.0*(atan( max ) + ATANERR);
	return( temp );
}


double rough_min_bvol( double y[12] )
{
	double x[12], yp[12], min, temp, maxa[4];
	int i;
	
	ROUND_DOWN;
	for( i=0; i<12; i+=2 )
		x[i] = y[i]*y[i];
	ROUND_UP;
	for( i=1; i<12; i+=2 )
		x[i] = y[i]*y[i];
	
	min = rough_min_delta( x );
	/*	ROUND_DOWN;	*/
	min = sqrt( min );
	
	maxa[0] = 2.0*max_a( y );
	for( i=0; i<2; i++ ) {
		yp[0 + i] = y[0 + i];
		yp[2 + i] = y[8 + i];
		yp[4 + i] = y[10 + i];
		yp[6 + i] = y[6 + i];
		yp[8 + i] = y[2 + i];
		yp[10 + i] = y[4 + i];
		}
	maxa[1] = 2.0*max_a( yp );
	for( i=0; i<2; i++ ) {
		yp[0 + i] = y[2 + i];
		yp[2 + i] = y[6 + i];
		/*	yp[4 + i] = y[10 + i];	*/
		yp[6 + i] = y[8 + i];
		yp[8 + i] = y[0 + i];
		/*	yp[10 + i] = y[4 + i];	*/
		}
	maxa[2] = 2.0*max_a( yp );
	for( i=0; i<2; i++ ) {
		yp[0 + i] = y[4 + i];
		/*	yp[2 + i] = y[6 + i];	*/
		yp[4 + i] = y[8 + i];
		yp[6 + i] = y[10 + i];
		/*	yp[8 + i] = y[0 + i];	*/
		yp[10 + i] = y[2 + i];
		}
	maxa[3] = 2.0*max_a( yp );
	ROUND_DOWN;
	temp = 0.0;
	for( i=0; i<4; i++ ) {
		temp += atan( min/maxa[i] ) - ATANERR;
	}
	return( temp*TWO_3_LO );
}


double rough_max_bvol( double y[12] )
{
	double x[12], yp[12], max, temp, mina[4];
	int i;
	
	ROUND_DOWN;
	for( i=0; i<12; i+=2 )
		x[i] = y[i]*y[i];
	ROUND_UP;
	for( i=1; i<12; i+=2 )
		x[i] = y[i]*y[i];
	
	max = rough_max_delta( x );
	/*	ROUND_UP;	*/
	max = sqrt( max );
	
	mina[0] = 2.0*min_a( y );
	for( i=0; i<2; i++ ) {
		yp[0 + i] = y[0 + i];
		yp[2 + i] = y[8 + i];
		yp[4 + i] = y[10 + i];
		yp[6 + i] = y[6 + i];
		yp[8 + i] = y[2 + i];
		yp[10 + i] = y[4 + i];
		}
	mina[1] = 2.0*min_a( yp );
	for( i=0; i<2; i++ ) {
		yp[0 + i] = y[2 + i];
		yp[2 + i] = y[6 + i];
		/*	yp[4 + i] = y[10 + i];	*/
		yp[6 + i] = y[8 + i];
		yp[8 + i] = y[0 + i];
		/*	yp[10 + i] = y[4 + i];	*/
		}
	mina[2] = 2.0*min_a( yp );
	for( i=0; i<2; i++ ) {
		yp[0 + i] = y[4 + i];
		/*	yp[2 + i] = y[6 + i];	*/
		yp[4 + i] = y[8 + i];
		yp[6 + i] = y[10 + i];
		/*	yp[8 + i] = y[0 + i];	*/
		yp[10 + i] = y[2 + i];
		}
	mina[3] = 2.0*min_a( yp );
	ROUND_UP;
	temp = 0.0;
	for( i=0; i<4; i++ ) {
		temp += atan( max/mina[i] ) + ATANERR;
	}
	return( temp*TWO_3_HI );
}


double s_min_bvol( double y[12], double sqrtdelta[2] )
{
	double yp[12], min, temp, maxa[4];
	int i;
	
	min = sqrtdelta[0];
	
	maxa[0] = 2.0*max_a( y );
	for( i=0; i<2; i++ ) {
		yp[0 + i] = y[0 + i];
		yp[2 + i] = y[8 + i];
		yp[4 + i] = y[10 + i];
		yp[6 + i] = y[6 + i];
		yp[8 + i] = y[2 + i];
		yp[10 + i] = y[4 + i];
		}
	maxa[1] = 2.0*max_a( yp );
	for( i=0; i<2; i++ ) {
		yp[0 + i] = y[2 + i];
		yp[2 + i] = y[6 + i];
		/*	yp[4 + i] = y[10 + i];	*/
		yp[6 + i] = y[8 + i];
		yp[8 + i] = y[0 + i];
		/*	yp[10 + i] = y[4 + i];	*/
		}
	maxa[2] = 2.0*max_a( yp );
	for( i=0; i<2; i++ ) {
		yp[0 + i] = y[4 + i];
		/*	yp[2 + i] = y[6 + i];	*/
		yp[4 + i] = y[8 + i];
		yp[6 + i] = y[10 + i];
		/*	yp[8 + i] = y[0 + i];	*/
		yp[10 + i] = y[2 + i];
		}
	maxa[3] = 2.0*max_a( yp );
	ROUND_DOWN;
	temp = 0.0;
	for( i=0; i<4; i++ ) {
		temp += atan( min/maxa[i] ) - ATANERR;
	}
	return( temp*TWO_3_LO );
}


double s_max_bvol( double y[12], double sqrtdelta[2] )
{
	double yp[12], max, temp, mina[4];
	int i;
	
	max = sqrtdelta[1];
	
	mina[0] = 2.0*min_a( y );
	for( i=0; i<2; i++ ) {
		yp[0 + i] = y[0 + i];
		yp[2 + i] = y[8 + i];
		yp[4 + i] = y[10 + i];
		yp[6 + i] = y[6 + i];
		yp[8 + i] = y[2 + i];
		yp[10 + i] = y[4 + i];
		}
	mina[1] = 2.0*min_a( yp );
	for( i=0; i<2; i++ ) {
		yp[0 + i] = y[2 + i];
		yp[2 + i] = y[6 + i];
		/*	yp[4 + i] = y[10 + i];	*/
		yp[6 + i] = y[8 + i];
		yp[8 + i] = y[0 + i];
		/*	yp[10 + i] = y[4 + i];	*/
		}
	mina[2] = 2.0*min_a( yp );
	for( i=0; i<2; i++ ) {
		yp[0 + i] = y[4 + i];
		/*	yp[2 + i] = y[6 + i];	*/
		yp[4 + i] = y[8 + i];
		yp[6 + i] = y[10 + i];
		/*	yp[8 + i] = y[0 + i];	*/
		yp[10 + i] = y[2 + i];
		}
	mina[3] = 2.0*min_a( yp );
	ROUND_UP;
	temp = 0.0;
	for( i=0; i<4; i++ ) {
		temp += atan( max/mina[i] ) + ATANERR;
	}
	return( temp*TWO_3_HI );
}


double s_max_dih( double x[12] )
{
	double xv[6], mindx4, mincosdih, u1, u2;
		
	mindx4 = s_min_delta4( x );
	/*	ROUND_DOWN;	*/
	if( mindx4 > 0.0 ) {
		xv[1] = x[1];
		xv[3] = x[3];
		xv[5] = x[11];
		u1 = s_max_u( xv );
		xv[3] = x[5];
		xv[5] = x[9];
		u2 = s_max_u( xv );
		/*	ROUND_UP;	*/
		u1 = sqrt( u1*u2 );
		}
	else {
		xv[0] = x[0];
		xv[2] = x[2];
		xv[4] = x[10];
		u1 = s_min_u( xv );
		xv[2] = x[4];
		xv[4] = x[8];
		u2 = s_min_u( xv );
		/*	ROUND_DOWN;	*/
		u1 = sqrt( u1*u2 );
		}
	ROUND_DOWN;
	mincosdih = mindx4/u1;
	u1 = max_acos( mincosdih );
	return( u1 );
}


double s_min_dih( double x[12] )
{
	double xv[6], maxdx4, maxcosdih, u1, u2;
	
	maxdx4 = s_max_delta4( x );
	/*	ROUND_UP;	*/
	if( maxdx4 > 0.0 ) {
		xv[0] = x[0];
		xv[2] = x[2];
		xv[4] = x[10];
		u1 = s_min_u( xv );
		xv[2] = x[4];
		xv[4] = x[8];
		u2 = s_min_u( xv );
		/*	ROUND_DOWN;	*/
		u1 = sqrt( u1*u2 );
		}
	else {
		xv[1] = x[1];
		xv[3] = x[3];
		xv[5] = x[11];
		u1 = s_max_u( xv );
		xv[3] = x[5];
		xv[5] = x[9];
		u2 = s_max_u( xv );
		/*	ROUND_UP;	*/
		u1 = sqrt( u1*u2 );
		}
	ROUND_UP;
	maxcosdih = maxdx4/u1;
	u1 = min_acos( maxcosdih );
	return( u1 );
}


/* s_dih computes the dihedral angle associated with
edge 1, and assumes a simplex with no long edge (greater
than 2 Sqrt[2]) */

void s_dih( double x[12], double out[2] )
{
	out[0] = s_min_dih( x );
	out[1] = s_max_dih( x );
}


/* As auxiliary functions, max/min_acos() are somewhat
misplaced.  However, this seemed a likely spot. */

double max_acos( double xv )
{
	double y, t;
	
	if( xv >= 1.0 )
		return( 0.0 );
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
		return( atan( y/xv ) + ATANERR);
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
			t = atan( xv/y ) - ATANERR;
			ROUND_UP;
			return( PI_2_HI - t );
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
			t = atan( (-xv)/y ) + ATANERR;
			return( PI_2_HI + t );
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
		t = atan( y/(-xv) ) - ATANERR;
		ROUND_UP;
		return( PI_HI - t );
		}
	else if( ISNAN( xv ) )
		return( xv );
	else
		return( PI_HI );
}


double min_acos( double xv )
{
	double y, t;
	
	if( xv >= 1.0 )
		return( 0.0 );
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
		return( atan( y/xv ) - ATANERR);
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
			t = atan( xv/y ) + ATANERR;
			ROUND_DOWN;
			return( PI_2_LO - t );
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
			t = atan( (-xv)/y ) - ATANERR;
			return( PI_2_LO + t );
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
		t = atan( y/(-xv) ) + ATANERR;
		ROUND_DOWN;
		return( PI_LO - t );
		}
	else if( ISNAN( xv ) )
		return( xv );
	else
		return( PI_LO );
}


/* max/min_qr_crad_test checks the circumradius of all tetrahedra
in the given cell.  It assumes that all tetrahedra are
quasi-regular (otherwise, the simplification upon which this
algorithm depends fails), and returns 0 if the test fails,
1 otherwise.  */

/* max_crad < 1.41?  */

int max_qr_crad_test( double y[12] )
{
	double x[6], x2[6], top, bot, ans;
	double pterms, nterms, p1, p2, p3, p4;
	int i, j;
	
	ROUND_UP;
	j = 1;
	for( i=0; i<6; i++ ) {
		x[i] = y[j]*y[j];
		j += 2;
		}
	/*	i_tomsrho( x, top );	*/
	pterms = 2.0*(x[0]*x[1]*x[3]*x[4] + x[0]*x[2]*x[3]*x[5] + 
		x[1]*x[2]*x[4]*x[5]);
	ROUND_DOWN;
	for( i=0; i<6; i++ )
		x2[i] = x[i]*x[i];
	nterms = x2[0]*x2[3] + x2[1]*x2[4] + x2[2]*x2[5];
	ROUND_UP;
	top = pterms - nterms;
	/*	i_bigdelta( x, bot );	*/
	ROUND_DOWN;
	p1 = x[0]*x[3]*(x[1] + x[2] + x[4] + x[5]);
	p2 = x[1]*x[4]*(x[0] + x[2] + x[3] + x[5]);
	p3 = x[2]*x[5]*(x[0] + x[1] + x[3] + x[4]);
	pterms = p1 + p2 + p3;
	ROUND_UP;
	p1 = x[0]*x[3]*(x[0] + x[3]);
	p2 = x[1]*x[4]*(x[1] + x[4]);
	p3 = x[2]*x[5]*(x[2] + x[5]);
	p4 = x[2]*(x[1]*x[3] + x[0]*x[4]) + x[5]*(x[0]*x[1] + 
		x[3]*x[4]);
	nterms = p1 + p2 + p3 + p4;
	ROUND_DOWN;
	bot = pterms - nterms;
	/*	i_div( top, bot, temp );	*/
	ROUND_UP;
	p1 = top/bot;
	ans = 0.5*sqrt( p1 );
	if( ans < 1.41 )
		return( 1 );
	else
		return( 0 );
}


/* min_crad > 1.41?  */

int min_qr_crad_test( double y[12] )
{
	double x[6], x2[6], top, bot, ans;
	double pterms, nterms, p1, p2, p3, p4;
	int i, j;
	
	/*	ROUND_UP;	*/
	ROUND_DOWN;
	j = 0;
	for( i=0; i<6; i++ ) {
		x[i] = y[j]*y[j];
		j += 2;
		}
	pterms = 2.0*(x[0]*x[1]*x[3]*x[4] + x[0]*x[2]*x[3]*x[5] + 
		x[1]*x[2]*x[4]*x[5]);
	/*	ROUND_DOWN;	*/
	ROUND_UP;
	for( i=0; i<6; i++ )
		x2[i] = x[i]*x[i];
	nterms = x2[0]*x2[3] + x2[1]*x2[4] + x2[2]*x2[5];
	/*	ROUND_UP;	*/
	ROUND_DOWN;
	top = pterms - nterms;
	/*	ROUND_DOWN;	*/
	ROUND_UP;
	p1 = x[0]*x[3]*(x[1] + x[2] + x[4] + x[5]);
	p2 = x[1]*x[4]*(x[0] + x[2] + x[3] + x[5]);
	p3 = x[2]*x[5]*(x[0] + x[1] + x[3] + x[4]);
	pterms = p1 + p2 + p3;
	/*	ROUND_UP;	*/
	ROUND_DOWN;
	p1 = x[0]*x[3]*(x[0] + x[3]);
	p2 = x[1]*x[4]*(x[1] + x[4]);
	p3 = x[2]*x[5]*(x[2] + x[5]);
	p4 = x[2]*(x[1]*x[3] + x[0]*x[4]) + x[5]*(x[0]*x[1] + 
		x[3]*x[4]);
	nterms = p1 + p2 + p3 + p4;
	/*	ROUND_DOWN;	*/
	ROUND_UP;
	bot = pterms - nterms;
	/*	ROUND_UP;	*/
	ROUND_DOWN;
	p1 = top/bot;
	ans = 0.5*sqrt( p1 );
	if( ans > 1.41 )
		return( 1 );
	else
		return( 0 );
}


double max_qr_crad( double y[12] )
{
	double x[6], x2[6], top, bot, ans;
	double pterms, nterms, p1, p2, p3, p4;
	int i, j;
	
	ROUND_UP;
	j = 1;
	for( i=0; i<6; i++ ) {
		x[i] = y[j]*y[j];
		j += 2;
		}
	/*	i_tomsrho( x, top );	*/
	pterms = 2.0*(x[0]*x[1]*x[3]*x[4] + x[0]*x[2]*x[3]*x[5] + 
		x[1]*x[2]*x[4]*x[5]);
	ROUND_DOWN;
	for( i=0; i<6; i++ )
		x2[i] = x[i]*x[i];
	nterms = x2[0]*x2[3] + x2[1]*x2[4] + x2[2]*x2[5];
	ROUND_UP;
	top = pterms - nterms;
	/*	i_bigdelta( x, bot );	*/
	ROUND_DOWN;
	p1 = x[0]*x[3]*(x[1] + x[2] + x[4] + x[5]);
	p2 = x[1]*x[4]*(x[0] + x[2] + x[3] + x[5]);
	p3 = x[2]*x[5]*(x[0] + x[1] + x[3] + x[4]);
	pterms = p1 + p2 + p3;
	ROUND_UP;
	p1 = x[0]*x[3]*(x[0] + x[3]);
	p2 = x[1]*x[4]*(x[1] + x[4]);
	p3 = x[2]*x[5]*(x[2] + x[5]);
	p4 = x[2]*(x[1]*x[3] + x[0]*x[4]) + x[5]*(x[0]*x[1] + 
		x[3]*x[4]);
	nterms = p1 + p2 + p3 + p4;
	ROUND_DOWN;
	bot = pterms - nterms;
	/*	i_div( top, bot, temp );	*/
	ROUND_UP;
	p1 = top/bot;
	ans = 0.5*sqrt( p1 );
	return( ans );
}


double min_qr_crad( double y[12] )
{
	double x[6], x2[6], top, bot, ans;
	double pterms, nterms, p1, p2, p3, p4;
	int i, j;
	
	/*	ROUND_UP;	*/
	ROUND_DOWN;
	j = 0;
	for( i=0; i<6; i++ ) {
		x[i] = y[j]*y[j];
		j += 2;
		}
	pterms = 2.0*(x[0]*x[1]*x[3]*x[4] + x[0]*x[2]*x[3]*x[5] + 
		x[1]*x[2]*x[4]*x[5]);
	/*	ROUND_DOWN;	*/
	ROUND_UP;
	for( i=0; i<6; i++ )
		x2[i] = x[i]*x[i];
	nterms = x2[0]*x2[3] + x2[1]*x2[4] + x2[2]*x2[5];
	/*	ROUND_UP;	*/
	ROUND_DOWN;
	top = pterms - nterms;
	/*	ROUND_DOWN;	*/
	ROUND_UP;
	p1 = x[0]*x[3]*(x[1] + x[2] + x[4] + x[5]);
	p2 = x[1]*x[4]*(x[0] + x[2] + x[3] + x[5]);
	p3 = x[2]*x[5]*(x[0] + x[1] + x[3] + x[4]);
	pterms = p1 + p2 + p3;
	/*	ROUND_UP;	*/
	ROUND_DOWN;
	p1 = x[0]*x[3]*(x[0] + x[3]);
	p2 = x[1]*x[4]*(x[1] + x[4]);
	p3 = x[2]*x[5]*(x[2] + x[5]);
	p4 = x[2]*(x[1]*x[3] + x[0]*x[4]) + x[5]*(x[0]*x[1] + 
		x[3]*x[4]);
	nterms = p1 + p2 + p3 + p4;
	/*	ROUND_DOWN;	*/
	ROUND_UP;
	bot = pterms - nterms;
	/*	ROUND_UP;	*/
	ROUND_DOWN;
	p1 = top/bot;
	ans = 0.5*sqrt( p1 );
	return( ans );
}


double rough_min_gma( double y[12] )
{
	double mtvol, mbvol, temp;
	
	mtvol = rough_max_tvol( y );
	/*	ROUND_UP;	*/
	temp = DOCT_HI*mtvol;
	mbvol = rough_min_bvol( y );
	/*	ROUND_DOWN;	*/
	return( mbvol - temp );
}


double rough_max_gma( double y[12] )
{
	double mtvol, mbvol, temp;
	
	mtvol = rough_min_tvol( y );
	/*	ROUND_DOWN;	*/
	temp = DOCT_LO*mtvol;
	mbvol = rough_max_bvol( y );
	/*	ROUND_UP;	*/
	return( mbvol - temp );
}


double s_min_gma( double y[12], double sqrtdelta[2] )
{
	double mtvol, mbvol, temp;
	
	ROUND_UP;
	mtvol = sqrtdelta[1]*ONE_12_HI;
	temp = DOCT_HI*mtvol;
	mbvol = s_min_bvol( y, sqrtdelta );
	/*	ROUND_DOWN;	*/
	return( mbvol - temp );
}


double s_max_gma( double y[12], double sqrtdelta[2] )
{
	double mtvol, mbvol, temp;
	
	ROUND_DOWN;
	mtvol = sqrtdelta[0]*ONE_12_LO;
	temp = DOCT_LO*mtvol;
	mbvol = s_max_bvol( y, sqrtdelta );
	/*	ROUND_UP;	*/
	return( mbvol - temp );
}


/* indexing is shifted in the normal C fashion */
void delta_partial( int n, double x[12], double out[2] )
{
	double pterms[2], nterms[2], xp[12];
	int indexing[6][6] = {{6, 8, 4, 0, 2, 10},
						{8, 10, 0, 2, 4, 6},
						{10, 0, 8, 4, 6, 2},
						{0, 2, 4, 6, 8, 10},
						{2, 10, 6, 8, 4, 0},
						{4, 0, 2, 10, 6, 8}};
	int i, j, k;
	
	k = 0;
	for( i=0; i<6; i++ ) {
		j = indexing[n][i];
		xp[k] = x[j];
		xp[k + 1] = x[j + 1];
		k += 2;
	}
	ROUND_DOWN;
	pterms[0] = xp[0]*(xp[2] + xp[4] + xp[8] + xp[10]) +
		xp[2]*xp[8] + xp[4]*xp[10];
	nterms[0] = xp[0]*(xp[0] + 2.0*xp[6]) + xp[2]*xp[4] + 
		xp[8]*xp[10];
	ROUND_UP;
	pterms[1] = xp[1]*(xp[3] + xp[5] + xp[9] + xp[11]) +
		xp[3]*xp[9] + xp[5]*xp[11];
	nterms[1] = xp[1]*(xp[1] + 2.0*xp[7]) + xp[3]*xp[5] + 
		xp[9]*xp[11];
	i_sub( pterms, nterms, out );
}


void i_delta_partials( double x[12], double part[12] )
{
	double pterms[2], nterms[2], xp[12];
	int indexing[6][6] = {{6, 8, 4, 0, 2, 10},
						{8, 10, 0, 2, 4, 6},
						{10, 0, 8, 4, 6, 2},
						{0, 2, 4, 6, 8, 10},
						{2, 10, 6, 8, 4, 0},
						{4, 0, 2, 10, 6, 8}};
	int i, j, k, n, m;
	
	ROUND_DOWN;
	m = 0;
	for( n=0; n<6; n++ ) {
		k = 0;
		for( i=0; i<6; i++ ) {
			j = indexing[n][i];
			xp[k] = x[j];
			xp[k + 1] = x[j + 1];
			k += 2;
		}
		pterms[0] = xp[0]*(xp[2] + xp[4] + xp[8] + xp[10]) +
			xp[2]*xp[8] + xp[4]*xp[10];
		nterms[0] = xp[0]*(xp[0] + 2.0*xp[6]) + xp[2]*xp[4] + 
			xp[8]*xp[10];
		ROUND_UP;
		pterms[1] = xp[1]*(xp[3] + xp[5] + xp[9] + xp[11]) +
			xp[3]*xp[9] + xp[5]*xp[11];
		nterms[1] = xp[1]*(xp[1] + 2.0*xp[7]) + xp[3]*xp[5] + 
			xp[9]*xp[11];
		part[m + 1] = pterms[1] - nterms[0];
		ROUND_DOWN;
		part[m] = pterms[0] - nterms[1];
		m += 2;
	}
}


/* Caveat:  no long ( > 2Sqrt[2] ) edges.  */
void s_delta_partials( double x[12], double part[12] )
{
	double pterms, nterms, xp[12];
	double max1, max2, temp, temp2, p1, p2, p3, val;
	int i, j, k, n, ind;
	int indexing[6][6] = {{6, 8, 4, 0, 2, 10},
						{8, 10, 0, 2, 4, 6},
						{10, 0, 8, 4, 6, 2},
						{0, 2, 4, 6, 8, 10},
						{2, 10, 6, 8, 4, 0},
						{4, 0, 2, 10, 6, 8}};
	
	ind = 0;
	for( n=0; n<6; n++ ) {
		k = 0;
		for( i=0; i<6; i++ ) {
			j = indexing[n][i];
			xp[k] = x[j];
			xp[k + 1] = x[j + 1];
			k += 2;
		}
		ROUND_DOWN;
		temp = xp[2] + xp[4] - 2.0*xp[7] + xp[8] + xp[10];
		p1 = xp[0]*(-xp[0] + temp);
		p2 = xp[1]*(-xp[1] + temp);
		val = MIN( p1, p2 );
		ROUND_UP;
		p2 = xp[2]*xp[4] + xp[8]*xp[10];
		ROUND_DOWN;
		p1 = xp[2]*xp[8] + xp[4]*xp[10];
		p3 = p1 - p2;
		val += p3;	/* lower bound */
		part[ind] = val;
		
		max1 = 0.5*(xp[3] + xp[5] - 2.0*xp[6] + xp[9] + xp[11]);
		/* max1 is lower bound on potential maximum */
		ROUND_UP;
		temp2 = xp[3] + xp[5] - 2.0*xp[6] + xp[9] + xp[11];
		max2 = 0.5*temp2;
		/* max2 is upper bound on potential maximum */
		if( max1 > xp[1] || max2 < xp[0] ) {
			p1 = xp[0]*(-xp[0] + temp2);
			p2 = xp[1]*(-xp[1] + temp2);
			val = MAX( p1, p2 );
			}
		else {
			ROUND_DOWN;
			nterms = max1*(max1 + 2.0*xp[6]);
			ROUND_UP;
			pterms = max2*(xp[3] + xp[5] + xp[9] + xp[11]);
			val = pterms - nterms;
			}
		ROUND_DOWN;
		p2 = xp[3]*xp[5] + xp[9]*xp[11];
		ROUND_UP;
		p1 = xp[3]*xp[9] + xp[5]*xp[11];
		p3 = p1 - p2;
		val += p3;	/* upper bound */
		ind++;
		part[ind] = val;
		
		ind++;
	}
}


/* Caveat:  no long ( > 2Sqrt[2] ) edges.  */
int s_delta_partial_sign( int n, double x[12] )
{
	double pterms, nterms, xp[12];
	double max1, max2, temp, temp2, p1, p2, p3, val;
	int i, j, k;
	int indexing[6][6] = {{6, 8, 4, 0, 2, 10},
						{8, 10, 0, 2, 4, 6},
						{10, 0, 8, 4, 6, 2},
						{0, 2, 4, 6, 8, 10},
						{2, 10, 6, 8, 4, 0},
						{4, 0, 2, 10, 6, 8}};
	
	k = 0;
	for( i=0; i<6; i++ ) {
		j = indexing[n][i];
		xp[k] = x[j];
		xp[k + 1] = x[j + 1];
		k += 2;
	}
	ROUND_DOWN;
	temp = xp[2] + xp[4] - 2.0*xp[7] + xp[8] + xp[10];
	p1 = xp[0]*(-xp[0] + temp);
	p2 = xp[1]*(-xp[1] + temp);
	val = MIN( p1, p2 );
	ROUND_UP;
	p2 = xp[2]*xp[4] + xp[8]*xp[10];
	ROUND_DOWN;
	p1 = xp[2]*xp[8] + xp[4]*xp[10];
	p3 = p1 - p2;
	val += p3;	/* lower bound */
	
	if( val > 0.0 )	 /* sign of derivative is positive*/
		return( 1 );
	
	max1 = 0.5*(xp[3] + xp[5] - 2.0*xp[6] + xp[9] + xp[11]);
	/* max1 is lower bound on potential maximum */
	ROUND_UP;
	temp2 = xp[3] + xp[5] - 2.0*xp[6] + xp[9] + xp[11];
	max2 = 0.5*temp2;
	/* max2 is upper bound on potential maximum */
	if( max1 > xp[1] || max2 < xp[0] ) {
		p1 = xp[0]*(-xp[0] + temp2);
		p2 = xp[1]*(-xp[1] + temp2);
		val = MAX( p1, p2 );
		}
	else {
		ROUND_DOWN;
		nterms = max1*(max1 + 2.0*xp[6]);
		ROUND_UP;
		pterms = max2*(xp[3] + xp[5] + xp[9] + xp[11]);
		val = pterms - nterms;
		}
	ROUND_DOWN;
	p2 = xp[3]*xp[5] + xp[9]*xp[11];
	ROUND_UP;
	p1 = xp[3]*xp[9] + xp[5]*xp[11];
	p3 = p1 - p2;
	val += p3;	/* upper bound */
	if( val < 0.0 )
		return( -1 );	/* sign is negative */
	else
		return( 0 );	/* sign indeterminate */
}


void a_partials( double y[12], double part[12] )
{
	double y2[12], pterms[2], nterms[2], temp[2];
	int i;

	ROUND_DOWN;
	for( i=0; i<12; i+=2 )
		y2[i] = y[i]*y[i];
	ROUND_UP;
	for( i=1; i<12; i+=2 )
		y2[i] = y[i]*y[i];
	
	temp[1] = y[1]*y[3] + y[1]*y[5] + y[3]*y[5];
	ROUND_DOWN;
	temp[0] = y[0]*y[2] + y[0]*y[4] + y[2]*y[4];
	/* a_y1 = y2*y3 + (y2^2 + y3^2 - y4^2)/2 + y1*y2 + y1*y3 */
	pterms[1] = temp[1] + 0.5*(y2[3] + y2[5]);
	nterms[1] = 0.5*y2[7];
	ROUND_DOWN;
	pterms[0] = temp[0] + 0.5*(y2[2] + y2[4]);
	nterms[0] = 0.5*y2[6];
	I_SUB( pterms, nterms, part );
	/*	ROUND_UP;	*/
	/* a_y2 = y1*y3 + y1*y2 + (y1^2 + y3^2 - y5^2)/2 + y2*y3 */
	pterms[1] = temp[1] + 0.5*(y2[1] + y2[5]);
	nterms[1] = 0.5*y2[9];
	ROUND_DOWN;
	pterms[0] = temp[0] + 0.5*(y2[0] + y2[4]);
	nterms[0] = 0.5*y2[8];
	i_sub( pterms, nterms, part + 2 );
	/*	ROUND_UP;	*/
	/* a_y3 = y1*y2 + y1*y3 + y2*y3 + (y1^2 + y2^2 - y6^2)/2 */
	pterms[1] = temp[1] + 0.5*(y2[1] + y2[3]);
	nterms[1] = 0.5*y2[11];
	ROUND_DOWN;
	pterms[0] = temp[0] + 0.5*(y2[0] + y2[2]);
	nterms[0] = 0.5*y2[10];
	i_sub( pterms, nterms, part + 4 );
	/*	ROUND_UP;	*/
	/* a_y4 = -y1*y4 */
	part[6] = -(y[1]*y[7]);
	ROUND_DOWN;
	part[7] = -(y[0]*y[6]);
	/* a_y5 = -y2*y5 */
	ROUND_UP;
	part[8] = -(y[3]*y[9]);
	ROUND_DOWN;
	part[9] = -(y[2]*y[8]);
	/* a_y6 = -y3*y6 */
	ROUND_UP;
	part[10] = -(y[5]*y[11]);
	ROUND_DOWN;
	part[11] = -(y[4]*y[10]);
}


void s_solid_partials( double y[12], double delta[2], 
	double sqrtdelta[2], double delta_part[12], 
	double sol_part[12] )
{
	double temp[2], temp2[2], coeff[2], afval[2], afval2[2];
	double a_part[12];
	int i;
	
	i_afunc( y, afval );
	I_SQUARELEN( afval, afval2 );
	temp[0] = 4.0*afval2[0];
	temp[1] = 4.0*afval2[1];
	I_ADD( temp, delta, temp2 );
	I_MULT( sqrtdelta, temp2, temp );
	afval2[0] = 4.0;
	afval2[1] = 4.0;
	i_div( afval2, temp, coeff );
	a_partials( y, a_part );
	for( i=0; i<12; i+=2 ) {
		i_mult( afval, y + i, temp );
		i_mult( temp, delta_part + i, temp2 );
		i_mult( delta, a_part + i, temp );
		I_SUB( temp2, temp, afval2 );
		i_mult( coeff, afval2, sol_part + i );
	}
}


void s_solid_xpars( double y[12], double ooy[12],
	double delta[2], double sqrtdelta[2], 
	double delta_part[12], double sol_part[12] )
{
	double temp[2], temp2[2], coeff[2], afval[2], afval2[2];
	double a_part[12];
	int i;
	
	i_afunc( y, afval );
	I_SQUARELEN( afval, afval2 );
	temp[0] = 4.0*afval2[0];
	temp[1] = 4.0*afval2[1];
	I_ADD( temp, delta, temp2 );
	I_MULT( sqrtdelta, temp2, temp );
	afval2[0] = 2.0;
	afval2[1] = 2.0;
	i_div( afval2, temp, coeff );
	a_partials( y, a_part );
	for( i=0; i<12; i+=2 ) {
		i_mult( delta, ooy + i, temp );
		i_mult( temp, a_part + i, temp2 );
		i_mult( afval, delta_part + i, temp );
		I_SUB( temp, temp2, afval2 );
		i_mult( coeff, afval2, sol_part + i );
	}
}


void s_bvol_partials( double y[12], double delta[2], 
	double sqrtdelta[2], double delta_part[12], 
	double sol_part[12], double bvol_part[12] )
{
	double temp[2], temp2[2], coeff[2], afval[2], afval2[2];
	double a_part[12], yp[12], dp[12], foo;
	int i, j, k, n, ind;
	int forward[3][6] =	{{8, 10, 0, 2, 4, 6},
						{8, 6, 4, 2, 0, 10},
						{10, 2, 6, 4, 8, 0}};
						
	int back[3][6] = 	{{4, 6, 8, 10, 0, 2},
						{8, 6, 4, 2, 0, 10},
						{10, 2, 6, 4, 8, 0}};
	
	for( i=0; i<12; i++ )
		bvol_part[i] = sol_part[i];
	for( n=0; n<3; n++ ) {
		/* swap entries around in y and delta_part */
		j = 0;
		for( i=0; i<6; i++ ) {
			ind = forward[n][i];
			yp[j] = y[ind];
			dp[j] = delta_part[ind];
			j++;
			ind++;
			yp[j] = y[ind];
			dp[j] = delta_part[ind];
			j++;
		}
		i_afunc( yp, afval );
		i_squarelen( afval, afval2 );
		temp[0] = 4.0*afval2[0];
		temp[1] = 4.0*afval2[1];
		i_add( temp, delta, temp2 );
		i_mult( sqrtdelta, temp2, temp );
		afval2[0] = 4.0;
		afval2[1] = 4.0;
		i_div( afval2, temp, coeff );
		a_partials( yp, a_part );
		for( k=0; k<12; k+=2 ) {
			i_mult( afval, yp + k, temp );
			i_mult( temp, dp + k, temp2 );
			i_mult( delta, a_part + k, temp );
			i_sub( temp2, temp, afval2 );
			i_mult( coeff, afval2, yp + k );
		}
		/* now swap back */
		j = 0;
		for( i=0; i<6; i++ ) {
			ind = back[n][i];
			bvol_part[j] += yp[ind];
			j++;
			ind++;
			bvol_part[j] += yp[ind];
			j++;
		}
	}
	ROUND_UP;
	for( i=1; i<12; i+=2 ) {
		foo = bvol_part[i];
		if( foo > 0.0 )
			foo *= ONE_3_HI;
		else
			foo *= ONE_3_LO;
		bvol_part[i] = foo;
	}
	ROUND_DOWN;
	for( i=0; i<12; i+=2 ) {
		foo = bvol_part[i];
		if( foo > 0.0 )
			foo *= ONE_3_LO;
		else
			foo *= ONE_3_HI;
		bvol_part[i] = foo;
	}
}


void s_bvol_xpars( double y[12], double ooy[12],
	double delta[2], double sqrtdelta[2], 
	double delta_part[12], double sol_xpart[12], 
	double bvol_part[12] )
{
	double temp[2], temp2[2], coeff[2], afval[2], afval2[2];
	double a_part[12], yp[12], ooyp[12], dp[12], foo;
	int i, j, k, n, ind;
	int forward[3][6] =	{{8, 10, 0, 2, 4, 6},
						{8, 6, 4, 2, 0, 10},
						{10, 2, 6, 4, 8, 0}};
						
	int back[3][6] = 	{{4, 6, 8, 10, 0, 2},
						{8, 6, 4, 2, 0, 10},
						{10, 2, 6, 4, 8, 0}};
	
	for( i=0; i<12; i++ )
		bvol_part[i] = sol_xpart[i];
	for( n=0; n<3; n++ ) {
		/* swap entries around in y and delta_part */
		j = 0;
		for( i=0; i<6; i++ ) {
			ind = forward[n][i];
			yp[j] = y[ind];
			dp[j] = delta_part[ind];
			ooyp[j] = ooy[ind];
			j++;
			ind++;
			yp[j] = y[ind];
			dp[j] = delta_part[ind];
			ooyp[j] = ooy[ind];
			j++;
		}
		i_afunc( yp, afval );
		i_squarelen( afval, afval2 );
		temp[0] = 4.0*afval2[0];
		temp[1] = 4.0*afval2[1];
		i_add( temp, delta, temp2 );
		i_mult( sqrtdelta, temp2, temp );
		afval2[0] = 2.0;
		afval2[1] = 2.0;
		i_div( afval2, temp, coeff );
		a_partials( yp, a_part );
		for( k=0; k<12; k+=2 ) {
			i_mult( delta, ooyp + k, temp );
			i_mult( temp, a_part + k, temp2 );
			i_mult( afval, dp + k, temp );
			i_sub( temp, temp2, afval2 );
			i_mult( coeff, afval2, yp + k );
		}
		/* now swap back */
		j = 0;
		for( i=0; i<6; i++ ) {
			ind = back[n][i];
			bvol_part[j] += yp[ind];
			j++;
			ind++;
			bvol_part[j] += yp[ind];
			j++;
		}
	}
	ROUND_UP;
	for( i=1; i<12; i+=2 ) {
		foo = bvol_part[i];
		if( foo > 0.0 )
			foo *= ONE_3_HI;
		else
			foo *= ONE_3_LO;
		bvol_part[i] = foo;
	}
	ROUND_DOWN;
	for( i=0; i<12; i+=2 ) {
		foo = bvol_part[i];
		if( foo > 0.0 )
			foo *= ONE_3_LO;
		else
			foo *= ONE_3_HI;
		bvol_part[i] = foo;
	}
}


/* This actually gives 0.5*uu_partials, which is convenient. */
void s_uu_partials( double x[12], double u2[2],
	double u3[2], double uu_part[12] )
{
	double temp, temp2;
	
	ROUND_DOWN;
	temp = x[2] + x[10] - x[1];
	temp2 = x[4] + x[8] - x[1];
	uu_part[0] = temp*u3[0] + temp2*u2[0];
	ROUND_UP;
	temp = x[3] + x[11] - x[0];
	temp2 = x[5] + x[9] - x[0];
	uu_part[1] = temp*u3[1] + temp2*u2[1];

	/*	ROUND_UP;	*/
	temp = x[1] + x[11] - x[2];
	uu_part[3] = temp*u3[1];
	ROUND_DOWN;
	temp = x[0] + x[10] - x[3];
	uu_part[2] = temp*u3[0];
	
	/*	ROUND_DOWN;	*/
	temp = x[0] + x[8] - x[5];
	uu_part[4] = temp*u2[0];
	ROUND_UP;
	temp = x[1] + x[9] - x[4];
	uu_part[5] = temp*u2[1];
	
	uu_part[6] = 0.0;
	uu_part[7] = 0.0;
	
	/*	ROUND_UP;	*/
	temp = x[1] + x[5] - x[8];
	uu_part[9] = temp*u2[1];
	ROUND_DOWN;
	temp = x[0] + x[4] - x[9];
	uu_part[8] = temp*u2[0];
	
	/*	ROUND_DOWN;	*/
	temp = x[0] + x[2] - x[11];
	uu_part[10] = temp*u3[0];
	ROUND_UP;
	temp = x[1] + x[3] - x[10];
	uu_part[11] = temp*u3[1];
}


/* This actually gives 0.5*uu_partials, which is convenient. */
void i_uu_partials( double x[12], double u2[2],
	double u3[2], double uu_part[12] )
{
	double temp[2], temp2[2], temp3[2], temp4[2];
	
	ROUND_DOWN;
	temp[0] = x[2] + x[10] - x[1];
	temp2[0] = x[4] + x[8] - x[1];
	ROUND_UP;
	temp[1] = x[3] + x[11] - x[0];
	temp2[1] = x[5] + x[9] - x[0];
	i_mult( temp, u3, temp3 );
	i_mult( temp2, u2, temp4 );
	I_ADD( temp3, temp4, uu_part );

	/* ROUND_UP; */
	temp[1] = x[1] + x[11] - x[2];
	temp2[1] = x[1] + x[9] - x[4];
	temp3[1] = x[1] + x[5] - x[8];
	temp4[1] = x[1] + x[3] - x[10];
	ROUND_DOWN;
	temp[0] = x[0] + x[10] - x[3];
	temp2[0] = x[0] + x[8] - x[5];
	temp3[0] = x[0] + x[4] - x[9];
	temp4[0] = x[0] + x[2] - x[11];
	i_mult( temp, u3, uu_part + 2 );
	i_mult( temp2, u2, uu_part + 4 );
	i_mult( temp3, u2, uu_part + 8 );
	i_mult( temp4, u3, uu_part + 10 );
	
	uu_part[6] = 0.0;
	uu_part[7] = 0.0;
}


void delta4_partials( double x[12], double part[12] )
{
	ROUND_DOWN;
	part[0] = x[2] + x[4] + x[8] + x[10] - 2.0*x[1] - 2.0*x[7];
	ROUND_UP;
	part[1] = x[3] + x[5] + x[9] + x[11] - 2.0*x[0] - 2.0*x[6];
	
	/*	ROUND_UP;	*/
	part[3] = x[1] + x[9] - x[4];
	ROUND_DOWN;
	part[2] = x[0] + x[8] - x[5];
	
	/*	ROUND_DOWN;	*/
	part[4] = x[0] + x[10] - x[3];
	ROUND_UP;
	part[5] = x[1] + x[11] - x[2];
	
	part[6] = -2.0*x[1];
	part[7] = -2.0*x[0];
	
	/*	ROUND_UP;	*/
	part[9] = x[1] + x[3] - x[10];
	ROUND_DOWN;
	part[8] = x[0] + x[2] - x[11];
	
	/*	ROUND_DOWN;	*/
	part[10] = x[0] + x[4] - x[9];
	ROUND_UP;
	part[11] = x[1] + x[5] - x[8];
}


void s_dih_partials( double x[12], double y[12], 
	double dih_part[12] )
{
	double u2[2], u3[2], uu[2], xv[6], uu_part[12];
	double delta4_part[12], delta4[2], temp[2], temp2[2];
	int i;

	xv[0] = x[0];
	xv[1] = x[1];
	xv[2] = x[2];
	xv[3] = x[3];
	xv[4] = x[10];
	xv[5] = x[11];
	i_tomsu( xv, u2 );

	xv[2] = x[4];
	xv[3] = x[5];
	xv[4] = x[8];
	xv[5] = x[9];
	i_tomsu( xv, u3 );

	ROUND_DOWN;
	uu[0] = u2[0]*u3[0];
	ROUND_UP;
	uu[1] = u2[1]*u3[1];
		
	s_uu_partials( x, u2, u3, uu_part );
	delta4_partials( x, delta4_part );
	delta4[0] = s_min_delta4( x );
	delta4[1] = s_max_delta4( x );
	
	i_mult( delta4, delta4, u2 );
	i_sub( uu, u2, u3 );
	i_sqrt( u3, u2 );
	/*	i_mult( uu, u2, u3 );	*/
	ROUND_DOWN;
	u3[0] = uu[0]*u2[0];
	ROUND_UP;
	u3[1] = uu[1]*u2[1];
	u2[0] = 1.0;
	u2[1] = 1.0;
	i_div( u2, u3, xv );
	
	for( i=0; i<12; i+=2 ) {
		i_mult( delta4_part + i, uu, u2 );
		i_mult( uu_part + i, delta4, u3 );
		i_sub( u3, u2, temp );
		i_mult( xv, temp, temp2 );
		/* Here temp2 is the x_i partial */
		i_mult( y + i, temp2, temp );
		dih_part[i] = 2.0*temp[0];
		dih_part[i+1] = 2.0*temp[1];
		/* Now we have the y_i partial */
	}
}


void i_dih_partials( double x[12], double y[12], 
	double dih_part[12] )
{
	double u2[2], u3[2], uu[2], xv[6], uu_part[12];
	double delta4_part[12], delta4[2], temp[2], temp2[2];
	double delta[2];
	int i;

	xv[0] = x[0];
	xv[1] = x[1];
	xv[2] = x[2];
	xv[3] = x[3];
	xv[4] = x[10];
	xv[5] = x[11];
	i_tomsu( xv, u2 );

	xv[2] = x[4];
	xv[3] = x[5];
	xv[4] = x[8];
	xv[5] = x[9];
	i_tomsu( xv, u3 );

	ROUND_DOWN;
	uu[0] = u2[0]*u3[0];
	ROUND_UP;
	uu[1] = u2[1]*u3[1];
		
	i_uu_partials( x, u2, u3, uu_part );
	delta4_partials( x, delta4_part );
	i_delta4( x, delta4 );
	
	i_bigdelta_best( x, delta );
	if( delta[0] < 0.0 )
		delta[0] = 0.0;
	ROUND_DOWN;
	u3[0] = 4.0*x[0]*delta[0];
	ROUND_UP;
	u3[1] = 4.0*x[1]*delta[1];

	i_sqrt( u3, u2 );
	/*	i_mult( uu, u2, u3 );	*/
	ROUND_DOWN;
	u3[0] = uu[0]*u2[0];
	ROUND_UP;
	u3[1] = uu[1]*u2[1];
	u2[0] = 1.0;
	u2[1] = 1.0;
	i_div( u2, u3, xv );
	
	for( i=0; i<12; i+=2 ) {
		i_mult( delta4_part + i, uu, u2 );
		i_mult( uu_part + i, delta4, u3 );
		i_sub( u3, u2, temp );
		i_mult( xv, temp, temp2 );
		/* Here temp2 is the x_i partial */
		i_mult( y + i, temp2, temp );
		dih_part[i] = 2.0*temp[0];
		dih_part[i+1] = 2.0*temp[1];
		/* Now we have the y_i partial */
	}
}


void i_dih_xpars( double x[12], double dih_part[12] )
{
	double u2[2], u3[2], uu[2], xv[6], uu_part[12];
	double delta4_part[12], delta4[2], temp[2];
	double delta[2];
	int i;

	xv[0] = x[0];
	xv[1] = x[1];
	xv[2] = x[2];
	xv[3] = x[3];
	xv[4] = x[10];
	xv[5] = x[11];
	i_tomsu( xv, u2 );

	xv[2] = x[4];
	xv[3] = x[5];
	xv[4] = x[8];
	xv[5] = x[9];
	i_tomsu( xv, u3 );

	ROUND_DOWN;
	uu[0] = u2[0]*u3[0];
	ROUND_UP;
	uu[1] = u2[1]*u3[1];

	i_uu_partials( x, u2, u3, uu_part );
	delta4_partials( x, delta4_part );
	i_delta4( x, delta4 );
	
	i_bigdelta_best( x, delta );
	if( delta[0] < 0.0 )
		delta[0] = 0.0;
	ROUND_DOWN;
	u3[0] = 4.0*x[0]*delta[0];
	ROUND_UP;
	u3[1] = 4.0*x[1]*delta[1];
/*	
	ROUND_NEAR;
	printf("u3 = [%.18f, %.18f]\n", u3[0], u3[1]);
*/	
	i_sqrt( u3, u2 );
	/*	i_mult( uu, u2, u3 );	*/
	ROUND_DOWN;
	u3[0] = uu[0]*u2[0];
	ROUND_UP;
	u3[1] = uu[1]*u2[1];
	u2[0] = 1.0;
	u2[1] = 1.0;
	i_div( u2, u3, xv );
	
	for( i=0; i<12; i+=2 ) {
		i_mult( delta4_part + i, uu, u2 );
		i_mult( uu_part + i, delta4, u3 );
		i_sub( u3, u2, temp );
		i_mult( xv, temp, dih_part + i );
	}
}


void old_dih_xpars( double x[12], double dih_part[12] )
{
	double u2[2], u3[2], uu[2], xv[6], uu_part[12];
	double delta4_part[12], delta4[2], temp[2];
	int i;

	xv[0] = x[0];
	xv[1] = x[1];
	xv[2] = x[2];
	xv[3] = x[3];
	xv[4] = x[10];
	xv[5] = x[11];
	i_tomsu( xv, u2 );

	xv[2] = x[4];
	xv[3] = x[5];
	xv[4] = x[8];
	xv[5] = x[9];
	i_tomsu( xv, u3 );

	ROUND_DOWN;
	uu[0] = u2[0]*u3[0];
	ROUND_UP;
	uu[1] = u2[1]*u3[1];

	i_uu_partials( x, u2, u3, uu_part );
	delta4_partials( x, delta4_part );
	i_delta4( x, delta4 );
	
	i_mult( delta4, delta4, u2 );
	i_sub( uu, u2, u3 );
/*	
	ROUND_NEAR;
	printf("u3 = [%.18f, %.18f]\n", u3[0], u3[1]);
*/	
	i_sqrt( u3, u2 );
	/*	i_mult( uu, u2, u3 );	*/
	ROUND_DOWN;
	u3[0] = uu[0]*u2[0];
	ROUND_UP;
	u3[1] = uu[1]*u2[1];
	u2[0] = 1.0;
	u2[1] = 1.0;
	i_div( u2, u3, xv );
	
	for( i=0; i<12; i+=2 ) {
		i_mult( delta4_part + i, uu, u2 );
		i_mult( uu_part + i, delta4, u3 );
		i_sub( u3, u2, temp );
		i_mult( xv, temp, dih_part + i );
	}
}


void s_dih_xpars( double x[12], double dih_part[12] )
{
	double u2[2], u3[2], uu[2], xv[6], uu_part[12];
	double delta4_part[12], delta4[2], temp[2];
	int i;

	xv[0] = x[0];
	xv[1] = x[1];
	xv[2] = x[2];
	xv[3] = x[3];
	xv[4] = x[10];
	xv[5] = x[11];
	i_tomsu( xv, u2 );

	xv[2] = x[4];
	xv[3] = x[5];
	xv[4] = x[8];
	xv[5] = x[9];
	i_tomsu( xv, u3 );

	ROUND_DOWN;
	uu[0] = u2[0]*u3[0];
	ROUND_UP;
	uu[1] = u2[1]*u3[1];
		
	s_uu_partials( x, u2, u3, uu_part );
	delta4_partials( x, delta4_part );
	delta4[0] = s_min_delta4( x );
	delta4[1] = s_max_delta4( x );
	
	i_mult( delta4, delta4, u2 );
	i_sub( uu, u2, u3 );
	i_sqrt( u3, u2 );
	/*	i_mult( uu, u2, u3 );	*/
	ROUND_DOWN;
	u3[0] = uu[0]*u2[0];
	ROUND_UP;
	u3[1] = uu[1]*u2[1];
	u2[0] = 1.0;
	u2[1] = 1.0;
	i_div( u2, u3, xv );
	
	for( i=0; i<12; i+=2 ) {
		i_mult( delta4_part + i, uu, u2 );
		i_mult( uu_part + i, delta4, u3 );
		i_sub( u3, u2, temp );
		i_mult( xv, temp, dih_part + i );
	}
}


void s_gma_partials( double y[12], double sqrtdelta[2], 
	double delta_part[12], double bvol_part[12], 
	double part[12] )
{
	double temp[2], temp2[2], coeff[2];
	int i;
	
	ROUND_DOWN;
	temp[0] = 12.0*sqrtdelta[0];
	ROUND_UP;
	temp[1] = 12.0*sqrtdelta[1];
	i_div( i_doct_const, temp, coeff );
	
	for( i=0; i<12; i+=2 ) {
		i_mult( coeff, delta_part + i, temp );
		i_mult( temp, y + i, temp2 );
		i_sub( bvol_part + i, temp2, part + i );
	}
}


void s_gma_xpars( double oosqrtdelta[2], double delta_part[12],
	double bvol_xpart[12], double part[12] )
{
	double temp[2], coeff[2];
	int i;
	
	ROUND_DOWN;
	coeff[0] = 0.5*ONE_12_LO*i_doct_const[0]*oosqrtdelta[0];
	ROUND_UP;
	coeff[1] = 0.5*ONE_12_HI*i_doct_const[1]*oosqrtdelta[1];
	
	for( i=0; i<12; i+=2 ) {
		i_mult( coeff, delta_part + i, temp );
		i_sub( bvol_xpart + i, temp, part + i );
	}
}


void auxp_partials( double x[6], double part[6] )
{
	ROUND_DOWN;
	part[0] = 2.0*(x[2] - x[1]) + x[4];
	part[2] = 2.0*(x[0] - x[3]) + x[4];
	part[4] = x[0] + x[2];
	ROUND_UP;
	part[1] = 2.0*(x[3] - x[0]) + x[5];
	part[3] = 2.0*(x[1] - x[2]) + x[5];
	part[5] = x[1] + x[3];
}


void tomsu_partials( double x[6], double part[6] )
{
	ROUND_DOWN;
	part[0] = 2.0*(x[2] + x[4] - x[1]);
	part[2] = 2.0*(x[0] + x[4] - x[3]);
	part[4] = 2.0*(x[0] + x[2] - x[5]);
	ROUND_UP;
	part[1] = 2.0*(x[3] + x[5] - x[0]);
	part[3] = 2.0*(x[1] + x[5] - x[2]);
	part[5] = 2.0*(x[1] + x[3] - x[4]);
}


void tomsv_partials( double x[12], double part[12] )
{
	double nterms[12], temp[6];
	
	ROUND_DOWN;
	nterms[0] = 2.0*(x[0]*x[6] + x[2]*x[10]);
	nterms[2] = 2.0*(x[2]*x[8] + x[0]*x[10]);
	nterms[4] = x[10]*x[10];
	nterms[6] = x[0]*x[0];
	nterms[8] = x[2]*x[2];
	nterms[10] = 2.0*(x[4]*x[10] + x[0]*x[2]);
	temp[0] = x[2] + x[10];
	temp[2] = x[0] + x[10];
	temp[4] = x[0] + x[2];
	
	ROUND_UP;
	nterms[1] = 2.0*(x[1]*x[7] + x[3]*x[11]);
	nterms[3] = 2.0*(x[3]*x[9] + x[1]*x[11]);
	nterms[5] = x[11]*x[11];
	nterms[7] = x[1]*x[1];
	nterms[9] = x[3]*x[3];
	nterms[11] = 2.0*(x[5]*x[11] + x[1]*x[3]);
	temp[1] = x[3] + x[11];
	temp[3] = x[1] + x[11];
	temp[5] = x[1] + x[3];
	
	/*	ROUND_UP;	*/
	part[1] = x[7]*temp[1] + x[3]*x[9] + 
		x[5]*x[11] - nterms[0];
	part[3] = x[1]*x[7] + x[9]*temp[3] + 
		x[5]*x[11] - nterms[2];
	part[5] = x[11]*temp[5] - nterms[4];
	part[7] = x[1]*temp[1] - nterms[6];
	part[9] = x[3]*temp[3] - nterms[8];
	part[11] = x[1]*x[7] + x[3]*x[9] + 
		x[5]*temp[5] - nterms[10];

	ROUND_DOWN;
	part[0] = x[6]*temp[0] + x[2]*x[8] + 
		x[4]*x[10] - nterms[1];
	part[2] = x[0]*x[6] + x[8]*temp[2] + 
		x[4]*x[10] - nterms[3];
	part[4] = x[10]*temp[4] - nterms[5];
	part[6] = x[0]*temp[0] - nterms[7];
	part[8] = x[2]*temp[2] - nterms[9];
	part[10] = x[0]*x[6] + x[2]*x[8] + 
		x[4]*temp[4] - nterms[11];
}


void i_auxp( double x[6], double out[2] )
{
	double temp[2], temp2[2];
	
	ROUND_DOWN;
	temp[0] = x[0] - x[3];
	temp2[0] = temp[0]*temp[0];
	ROUND_UP;	
	temp[1] = x[1] - x[2];
	temp2[1] = temp[1]*temp[1];
	/*	ROUND_UP;	*/
	out[1] = x[1]*x[5] + x[3]*x[5] - temp2[0];
	ROUND_DOWN;
	out[0] = x[0]*x[4] + x[2]*x[4] - temp2[1];
}


/* This gives the x_i partials of tomsP */
void tomsP_partials( double x[12], double delta[2],
	double delta32[2], double deltapart[12], 
	double part[12] )
{
	double xv[6], uval[2], vval[2], pval[2], rden[2];
	double auxp_par[12], tomsu_par[12], tomsv_par[12];
	double p1[2], p2[2], p3[2], lhs[2];
	int i;
	
	for( i=0; i<4; i++ )
		xv[i] = x[i];
	xv[4] = x[10];
	xv[5] = x[11];

	i_tomsu( xv, uval );
	i_auxp( xv, pval );
	i_tomsv( x, vval );
	
	auxp_partials( xv, part );
	for( i=0; i<4; i++ )
		auxp_par[i] = part[i];
	auxp_par[10] = part[4];
	auxp_par[11] = part[5];
	for( i=4; i<10; i++ )
		auxp_par[i] = 0.0;
	
	tomsu_partials( xv, part );
	for( i=0; i<4; i++ )
		tomsu_par[i] = part[i];
	tomsu_par[10] = part[4];
	tomsu_par[11] = part[5];
	for( i=4; i<10; i++ )
		tomsu_par[i] = 0.0;
	
	tomsv_partials( x, tomsv_par );
	
	i_squarelen( uval, p2 );
	ROUND_DOWN;
	p3[0] = 48.0*p2[0];
	ROUND_UP;
	p3[1] = 48.0*p2[1];
	i_mult( delta32, p3, p2 );
	ROUND_DOWN;
	rden[0] = 1.0/p2[1];
	ROUND_UP;
	rden[1] = 1.0/p2[0];
	
	for( i=0; i<12; i+=2 ) {
		i_mult( tomsv_par + i, pval, p1 );
		i_mult( auxp_par + i, vval, p2 );
		i_add( p1, p2, p3 );
		i_mult( p3, uval, p1 );
		i_mult( p1, delta, lhs );
		i_mult( tomsu_par + i, delta, p1 );
		i_mult( deltapart + i, uval, p2 );
		p3[0] = 0.5*p2[0];
		p3[1] = 0.5*p2[1];
		i_add( p1, p3, p2 );
		i_mult( p2, pval, p1 );
		i_mult( p1, vval, p2 );
		i_sub( lhs, p2, p1 );
		i_mult( p1, rden, part + i );
	}
}


void vorvol_partials( double y[12], double x[12], 
	double delta[2], double sqrtdelta[2],
	double delta_part[12], double part[12] )
{
	double xp[12], dp[12], psum[12], delta32[2];
	int i, j, k, n, ind;
	int forward[3][6] =	{{4, 0, 2, 10, 6, 8},
						{2, 4, 0, 8, 10, 6}};
	/*	
	int back[3][6] = 	{{2, 4, 0, 8, 10, 6},
						{4, 0, 2, 10, 6, 8}};
	*/
	
	ROUND_DOWN;
	delta32[0] = delta[0]*sqrtdelta[0];
	ROUND_UP;
	delta32[1] = delta[1]*sqrtdelta[1];
	
	tomsP_partials( x, delta, delta32, delta_part, psum );
	for( n=0; n<2; n++ ) {
		/* swap entries around in x and delta_part */
		j = 0;
		for( i=0; i<6; i++ ) {
			ind = forward[n][i];
			xp[j] = x[ind];
			dp[j] = delta_part[ind];
			j++;
			ind++;
			xp[j] = x[ind];
			dp[j] = delta_part[ind];
			j++;
		}
		
		tomsP_partials( xp, delta, delta32, dp, part );

		/* now swap back */
		if( n == 0 )
			k = 1;
		else
			k = 0;
		j = 0;
		for( i=0; i<6; i++ ) {
			ind = forward[k][i];
			psum[j] += part[ind];
			j++;
			ind++;
			psum[j] += part[ind];
			j++;
		}
	}
	for( i=0; i<12; i+=2 ) {
		i_mult( psum + i, y + i, xp );
		part[i] = 2.0*xp[0];
		part[i+1] = 2.0*xp[1];
	}
}


void vorvol_xpars( double x[12], double delta[2], 
	double sqrtdelta[2], double delta_part[12], 
	double part[12] )
{
	double xp[12], dp[12], psum[12], delta32[2];
	int i, j, k, n, ind;
	int forward[3][6] =	{{4, 0, 2, 10, 6, 8},
						{2, 4, 0, 8, 10, 6}};
	/*	
	int back[3][6] = 	{{2, 4, 0, 8, 10, 6},
						{4, 0, 2, 10, 6, 8}};
	*/
	
	ROUND_DOWN;
	delta32[0] = delta[0]*sqrtdelta[0];
	ROUND_UP;
	delta32[1] = delta[1]*sqrtdelta[1];
	
	tomsP_partials( x, delta, delta32, delta_part, psum );
	for( n=0; n<2; n++ ) {
		/* swap entries around in x and delta_part */
		j = 0;
		for( i=0; i<6; i++ ) {
			ind = forward[n][i];
			xp[j] = x[ind];
			dp[j] = delta_part[ind];
			j++;
			ind++;
			xp[j] = x[ind];
			dp[j] = delta_part[ind];
			j++;
		}
		
		tomsP_partials( xp, delta, delta32, dp, part );

		/* now swap back */
		if( n == 0 )
			k = 1;
		else
			k = 0;
		j = 0;
		for( i=0; i<6; i++ ) {
			ind = forward[k][i];
			psum[j] += part[ind];
			j++;
			ind++;
			psum[j] += part[ind];
			j++;
		}
	}
	for( i=0; i<12; i++ ) {
		part[i] = psum[i];
	}
}


void vor_partials( double y[12], double x[12], 
	double delta[2], double sqrtdelta[2], 
	double delta_part[12], double sol_part[12], 
	double part[12] )
{
	double temp[2], temp2[2], vorvol_part[12], foo;
	int i;
	
	vorvol_partials( y, x, delta, sqrtdelta, delta_part, 
		vorvol_part );
	
	for( i=0; i<12; i+=2 ) {
		i_mult( i_doct_const, vorvol_part + i, temp );
		ROUND_DOWN;
		temp2[0] = 3.0*temp[0];
		ROUND_UP;
		temp2[1] = 3.0*temp[1];
		i_sub( sol_part + i, temp2, temp );
		ROUND_DOWN;
		foo = 4.0*temp[0];
		if( foo > 0.0 )
			foo *= ONE_3_LO;
		else
			foo *= ONE_3_HI;
		part[i] = foo;
		ROUND_UP;
		foo = 4.0*temp[1];
		if( foo > 0.0 )
			foo *= ONE_3_HI;
		else
			foo *= ONE_3_LO;
		part[i+1] = foo;
	}
}


void vor_xpars( double x[12], double delta[2], 
	double sqrtdelta[2], double delta_part[12], 
	double sol_xpart[12], double part[12] )
{
	double one3[2], t1[2], t2[2], t3[2], vorvol_part[12];
	int i;
	
	vorvol_xpars( x, delta, sqrtdelta, delta_part, 
		vorvol_part );
	
	one3[0] = ONE_3_LO;
	one3[1] = ONE_3_HI;
	
	for( i=0; i<12; i+=2 ) {
		i_mult( i_doct_const, vorvol_part + i, t1 );
		i_mult( one3, sol_xpart + i, t2 );
		I_SUB( t2, t1, t3 );
		part[i] = 4.0*t3[0];
		part[i+1] = 4.0*t3[1];
	}
}


void octa_vor_partials( double y[12], double x[12], 
	double delta[2], double sqrtdelta[2], 
	double delta_part[12], double sol_part[12], 
	double part[12] )
{
	double alt_y[12], alt_x[12], alt_delta_part[12];
	double alt_sol_part[12], vor_part[12];
	double alt_vor_part[12];
	int i;
	
	vor_partials( y, x, delta, sqrtdelta, delta_part,
		sol_part, vor_part );
	
	/* Figure out delta_partials, etc. (how to recycle) */
	/* (1 5 6 4 2 3)  (becomes (0 8 10 6 2 4))			*/
	
	/* Do lots of permutations: */
	for( i=0; i<2; i++ ) {
		alt_y[i     ] 	= y[i];
		alt_y[i +  8] 	= y[i + 2];
		alt_y[i + 10] 	= y[i + 4];
		alt_y[i +  6] 	= y[i + 6];
		alt_y[i +  2] 	= y[i + 8];
		alt_y[i +  4] 	= y[i + 10];
	}
	for( i=0; i<2; i++ ) {
		alt_x[i     ] 	= x[i];
		alt_x[i +  8] 	= x[i + 2];
		alt_x[i + 10] 	= x[i + 4];
		alt_x[i +  6] 	= x[i + 6];
		alt_x[i +  2] 	= x[i + 8];
		alt_x[i +  4] 	= x[i + 10];
	}
	for( i=0; i<2; i++ ) {
		alt_delta_part[i     ] 	= delta_part[i];
		alt_delta_part[i +  8] 	= delta_part[i + 2];
		alt_delta_part[i + 10] 	= delta_part[i + 4];
		alt_delta_part[i +  6] 	= delta_part[i + 6];
		alt_delta_part[i +  2] 	= delta_part[i + 8];
		alt_delta_part[i +  4] 	= delta_part[i + 10];
	}
	
	/* Now need to compute alt_sol_part */
	
	s_solid_partials( alt_y, delta, sqrtdelta, alt_delta_part, 
		alt_sol_part );
	
	vor_partials( alt_y, alt_x, delta, sqrtdelta, 
		alt_delta_part, alt_sol_part, alt_vor_part );
	
	/* Now permute back: */
	for( i=0; i<2; i++ ) {
		alt_x[i     ] 	= alt_vor_part[i];
		alt_x[i +  8] 	= alt_vor_part[i + 2];
		alt_x[i + 10] 	= alt_vor_part[i + 4];
		alt_x[i +  6] 	= alt_vor_part[i + 6];
		alt_x[i +  2] 	= alt_vor_part[i + 8];
		alt_x[i +  4] 	= alt_vor_part[i + 10];
	}
	
	ROUND_DOWN;
	for( i=0; i<12; i+=2 ) {
		part[i] = 0.5*(vor_part[i] + alt_x[i]);
	}
	ROUND_UP;
	for( i=1; i<12; i+=2 ) {
		part[i] = 0.5*(vor_part[i] + alt_x[i]);
	}
}


void octa_vor_xpars( double y[12], double ooy[12],
	double x[12], double delta[2], double sqrtdelta[2], 
	double delta_part[12], double sol_xpart[12], 
	double part[12] )
{
	double alt_y[12], alt_x[12], alt_delta_part[12];
	double alt_sol_part[12], vor_part[12];
	double alt_vor_part[12], alt_ooy[12];
	int i;
	
	vor_xpars( x, delta, sqrtdelta, delta_part,
		sol_xpart, vor_part );
	
	/* Figure out delta_partials, etc. (how to recycle) */
	/* (1 5 6 4 2 3)  (becomes (0 8 10 6 2 4))			*/
	
	/* Do lots of permutations: */
	for( i=0; i<2; i++ ) {
		alt_y[i     ] 	= y[i];
		alt_y[i +  8] 	= y[i + 2];
		alt_y[i + 10] 	= y[i + 4];
		alt_y[i +  6] 	= y[i + 6];
		alt_y[i +  2] 	= y[i + 8];
		alt_y[i +  4] 	= y[i + 10];

		alt_ooy[i     ] 	= ooy[i];
		alt_ooy[i +  8] 	= ooy[i + 2];
		alt_ooy[i + 10] 	= ooy[i + 4];
		alt_ooy[i +  6] 	= ooy[i + 6];
		alt_ooy[i +  2] 	= ooy[i + 8];
		alt_ooy[i +  4] 	= ooy[i + 10];
	}
	for( i=0; i<2; i++ ) {
		alt_x[i     ] 	= x[i];
		alt_x[i +  8] 	= x[i + 2];
		alt_x[i + 10] 	= x[i + 4];
		alt_x[i +  6] 	= x[i + 6];
		alt_x[i +  2] 	= x[i + 8];
		alt_x[i +  4] 	= x[i + 10];
	}
	for( i=0; i<2; i++ ) {
		alt_delta_part[i     ] 	= delta_part[i];
		alt_delta_part[i +  8] 	= delta_part[i + 2];
		alt_delta_part[i + 10] 	= delta_part[i + 4];
		alt_delta_part[i +  6] 	= delta_part[i + 6];
		alt_delta_part[i +  2] 	= delta_part[i + 8];
		alt_delta_part[i +  4] 	= delta_part[i + 10];
	}
	
	/* Now need to compute alt_sol_part */
	
	s_solid_xpars( alt_y, alt_ooy, delta, sqrtdelta, 
		alt_delta_part, alt_sol_part );
	
	vor_xpars( alt_x, delta, sqrtdelta, 
		alt_delta_part, alt_sol_part, alt_vor_part );
	
	/* Now permute back: */
	for( i=0; i<2; i++ ) {
		alt_x[i     ] 	= alt_vor_part[i];
		alt_x[i +  8] 	= alt_vor_part[i + 2];
		alt_x[i + 10] 	= alt_vor_part[i + 4];
		alt_x[i +  6] 	= alt_vor_part[i + 6];
		alt_x[i +  2] 	= alt_vor_part[i + 8];
		alt_x[i +  4] 	= alt_vor_part[i + 10];
	}
	
	ROUND_DOWN;
	for( i=0; i<12; i+=2 ) {
		part[i] = 0.5*(vor_part[i] + alt_x[i]);
	}
	ROUND_UP;
	for( i=1; i<12; i+=2 ) {
		part[i] = 0.5*(vor_part[i] + alt_x[i]);
	}
}


/* If the face is acute (satisfied if no edge is longer
than 2 Sqrt[2]), then crad3len is monotonic in the edge
lengths (or x_i). */
/* crad3len = Sqrt[x1 x2 x3/tomsu[x1, x2, x3]] */
void s_crad3x( double x[6], double out[2] )
{
	double pterms[2], nterms[2], top[2], bot[2];

	ROUND_DOWN;
	top[0] = x[0]*x[2]*x[4];
	nterms[0] = x[0]*x[0] + x[2]*x[2] + x[4]*x[4];
	pterms[1] = 2.0*(x[1]*x[3] + x[1]*x[5] + x[3]*x[5]);
	ROUND_UP;
	top[1] = x[1]*x[3]*x[5];
	nterms[1] = x[1]*x[1] + x[3]*x[3] + x[5]*x[5];
	pterms[0] = 2.0*(x[0]*x[2] + x[0]*x[4] + x[2]*x[4]);
	bot[0] = pterms[0] - nterms[0];
	ROUND_DOWN;
	bot[1] = pterms[1] - nterms[1];
	pterms[0] = top[0]/bot[0];
	ROUND_UP;
	pterms[1] = top[1]/bot[1];
	i_sqrt( pterms, out );
}


/* Caveat:  no long ( > 2Sqrt[2] ) edges.  */
/* crad3x2 = crad3x^2 */
void s_crad3x2( double x[6], double out[2] )
{
	double pterms[2], nterms[2], top[2], bot[2];

	ROUND_DOWN;
	top[0] = x[0]*x[2]*x[4];
	nterms[0] = x[0]*x[0] + x[2]*x[2] + x[4]*x[4];
	pterms[1] = 2.0*(x[1]*x[3] + x[1]*x[5] + x[3]*x[5]);
	ROUND_UP;
	top[1] = x[1]*x[3]*x[5];
	nterms[1] = x[1]*x[1] + x[3]*x[3] + x[5]*x[5];
	pterms[0] = 2.0*(x[0]*x[2] + x[0]*x[4] + x[2]*x[4]);
	bot[0] = pterms[0] - nterms[0];
	ROUND_DOWN;
	bot[1] = pterms[1] - nterms[1];
	out[0] = top[0]/bot[0];
	ROUND_UP;
	out[1] = top[1]/bot[1];
}


/* Caveat:  no long ( > 2Sqrt[2] ) edges.  */
/* Formula:  x1 = x2 + x6 + 
	( sqrt(x2*x6*(8-x2)*(8-x6)) - x2*x6 )/4 */
/* No useful partials, so use crude bounds . . . */
void face_edge( double x2[2], double x6[2], double out[2] )
{
	double temp[2], garg[2];
	
	ROUND_DOWN;
	temp[0] = x2[0]*x6[0];
	garg[0] = sqrt( temp[0]*(8.0 - x2[1])*(8.0 - x6[1]) );
	ROUND_UP;
	temp[1] = x2[1]*x6[1];
	garg[1] = sqrt( temp[1]*(8.0 - x2[0])*(8.0 - x6[0]) );
	out[1] = x2[1] + x6[1] + 0.25*( garg[1] - temp[0] );
	ROUND_DOWN;
	out[0] = x2[0] + x6[0] + 0.25*( garg[0] - temp[1] );
}


/* Caveat:  no long ( > 2Sqrt[2] ) edges.  */
/* Formula:  d x6/ d x2 = 
	1 - 1/4 (x1 + (x2-4) sqrt( x1 (8-x1)/((x2 (8-x2)) ) ) ) */

void face_x6_x2( double x1[2], double x2[2], double out[2] )
{
	double top[2], bot[2], temp[2];
	
	ROUND_DOWN;
	top[0] = x1[0]*(8.0 - x1[1]);
	bot[0] = x2[0]*(8.0 - x2[1]);
	ROUND_UP;
	top[1] = x1[1]*(8.0 - x1[0]);
	bot[1] = x2[1]*(8.0 - x2[0]);
	temp[1] = top[1]/bot[0];
	ROUND_DOWN;
	temp[0] = top[0]/bot[1];
	bot[0] = sqrt( temp[0] );
	top[0] = x1[0] + (x2[0] - 4.0)*bot[0];
	ROUND_UP;
	bot[1] = sqrt( temp[1] );
	top[1] = x1[1] + (x2[1] - 4.0)*bot[1];
	out[1] = 1.0 - 0.25*top[0];
	ROUND_DOWN;
	out[0] = 1.0 - 0.25*top[1];
}


/* Formula:  y1_par = (x1+x2-x3)*(x1-x2+x3)*y2*y3/ufun^(3/2) */
void crad3len_pars( double y[6], double x[6],
	double part[6] )
{
	double xdiff[6], pterms[2], nterms[2], uval, sqrtuval;
	double uval32[2], foo[2], intprod[6], prod[6];
	int i;
	
	ROUND_DOWN;
	/* -x1 + x2 + x3 */
	xdiff[0] = -x[1] + x[2] + x[4];
	/* x1 - x2 + x3 */
	xdiff[2] = x[0] - x[3] + x[4];
	/* x1 + x2 - x3 */
	xdiff[4] = x[0] + x[2] - x[5];
	ROUND_UP;
	xdiff[1] = -x[0] + x[3] + x[5];
	xdiff[3] = x[1] - x[2] + x[5];
	xdiff[5] = x[1] + x[3] - x[4];
	
	for( i=0; i<6; i+=2 ) {
		i_mult( y + i, xdiff + i, intprod + i );
	}
	i_mult( intprod + 2, intprod + 4, prod );
	i_mult( intprod, intprod + 4, prod + 2 );
	i_mult( intprod, intprod + 2, prod + 4 );
	
	ROUND_DOWN;
	/* use pterms & nterms to compute tomsu */
	pterms[0] = 2.0*(x[0]*(x[4] + x[2]) + x[2]*x[4]);
	nterms[0] = x[0]*x[0] + x[2]*x[2] + x[4]*x[4];

	ROUND_UP;
	pterms[1] = 2.0*(x[1]*(x[5] + x[3]) + x[3]*x[5]);
	nterms[1] = x[1]*x[1] + x[3]*x[3] + x[5]*x[5];
	uval = pterms[1] - nterms[0];
	sqrtuval = sqrt( uval );
	uval32[1] = uval*sqrtuval;
	
	ROUND_DOWN;
	uval = pterms[0] - nterms[1];
	if( uval < 0.0 )
		uval = 0.0;
	sqrtuval = sqrt( uval );
	uval32[0] = uval*sqrtuval;
	foo[0] = 1.0/uval32[1];
	ROUND_UP;
	foo[1] = 1.0/uval32[0];
	for( i=0; i<6; i+=2 ) {
		i_mult( foo, prod + i, part + i );
	}
}


/* Formula:  y1_par = (x1+x2-x3)*(x1-x2+x3)*y2*y3/ufun^(3/2) */
/* Note that due to constraints on edge lengths, we are
guaranteed that -x1 + x2 + x6 will be positive (for example).*/
void s_crad3len_pars( double y[6], double x[6],
	double part[6] )
{
	double xdiff[3], pterms[2], nterms[2], uval, sqrtuval;
	double uval32[2], foo, intprod[3], prod[6];
	int i;
	
	ROUND_DOWN;
	/* -x1 + x2 + x3 */
	xdiff[0] = -x[1] + x[2] + x[4];
	/* x1 - x2 + x3 */
	xdiff[1] = x[0] - x[3] + x[4];
	/* x1 + x2 - x3 */
	xdiff[2] = x[0] + x[2] - x[5];
	
	/* use pterms & nterms to compute tomsu */
	pterms[0] = 2.0*(x[0]*(x[4] + x[2]) + x[2]*x[4]);
	nterms[0] = x[0]*x[0] + x[2]*x[2] + x[4]*x[4];
	
	intprod[0] = xdiff[0]*y[0];
	intprod[1] = xdiff[1]*y[2];
	intprod[2] = xdiff[2]*y[4];
	
	prod[0] = intprod[1]*intprod[2];
	prod[2] = intprod[0]*intprod[2];
	prod[4] = intprod[0]*intprod[1];
	
	ROUND_UP;
	xdiff[0] = -x[0] + x[3] + x[5];
	xdiff[1] = x[1] - x[2] + x[5];
	xdiff[2] = x[1] + x[3] - x[4];
	
	pterms[1] = 2.0*(x[1]*(x[5] + x[3]) + x[3]*x[5]);
	nterms[1] = x[1]*x[1] + x[3]*x[3] + x[5]*x[5];
	uval = pterms[1] - nterms[0];
	sqrtuval = sqrt( uval );
	uval32[1] = uval*sqrtuval;

	intprod[0] = xdiff[0]*y[1];
	intprod[1] = xdiff[1]*y[3];
	intprod[2] = xdiff[2]*y[5];
	
	prod[1] = intprod[1]*intprod[2];
	prod[3] = intprod[0]*intprod[2];
	prod[5] = intprod[0]*intprod[1];
	
	ROUND_DOWN;
	uval = pterms[0] - nterms[1];
	if( uval < 0.0 )
		uval = 0.0;
	sqrtuval = sqrt( uval );
	uval32[0] = uval*sqrtuval;
	foo = 1.0/uval32[1];
	for( i=0; i<6; i+=2 ) {
		part[i] = foo*prod[i];
	}
	ROUND_UP;
	foo = 1.0/uval32[0];
	for( i=1; i<6; i+=2 ) {
		part[i] = foo*prod[i];
	}
}


/* Formula:  x1_par = x2*(x1-x2+x3)*x3*(x1+x2-x3)/ufun^2 */
void crad3len2_xpars( double x[6], double part[6] )
{
	double xdiff[6], pterms[2], nterms[2], uval2[2];
	double foo, intprod[6], prod[6];
	int i;
	
	ROUND_DOWN;
	/* -x1 + x2 + x3 */
	xdiff[0] = -x[1] + x[2] + x[4];
	/* x1 - x2 + x3 */
	xdiff[2] = x[0] - x[3] + x[4];
	/* x1 + x2 - x3 */
	xdiff[4] = x[0] + x[2] - x[5];
	ROUND_UP;
	xdiff[1] = -x[0] + x[3] + x[5];
	xdiff[3] = x[1] - x[2] + x[5];
	xdiff[5] = x[1] + x[3] - x[4];

	for( i=0; i<6; i+=2 ) {
		i_mult( x + i, xdiff + i, intprod + i );
	}
	i_mult( intprod + 2, intprod + 4, prod );
	i_mult( intprod, intprod + 4, prod + 2 );
	i_mult( intprod, intprod + 2, prod + 4 );
	
	/* use pterms & nterms to compute tomsu */
	pterms[0] = 2.0*(x[0]*(x[4] + x[2]) + x[2]*x[4]);
	nterms[0] = x[0]*x[0] + x[2]*x[2] + x[4]*x[4];
	
	ROUND_UP;
	pterms[1] = 2.0*(x[1]*(x[5] + x[3]) + x[3]*x[5]);
	nterms[1] = x[1]*x[1] + x[3]*x[3] + x[5]*x[5];
	foo = pterms[1] - nterms[0];
	uval2[1] = foo*foo;

	ROUND_DOWN;
	foo = pterms[0] - nterms[1];
	if( foo < 0.0 )
		foo = 0.0;
	uval2[0] = foo*foo;
	pterms[0] = 1.0/uval2[1];
	ROUND_UP;
	pterms[1] = 1.0/uval2[0];
	for( i=0; i<6; i+=2 ) {
		i_mult( pterms, prod + i, part + i );
	}
}


/* Formula:  x1_par = x2*(x1-x2+x3)*x3*(x1+x2-x3)/ufun^2 */
/* Note that due to constraints on edge lengths, we are
guaranteed that -x1 + x2 + x6 will be positive (for example).*/
void s_crad3len2_xpars( double x[6], double part[6] )
{
	double xdiff[3], pterms[2], nterms[2], uval2[2];
	double foo, intprod[3], prod[6];
	int i;
	
	ROUND_DOWN;
	/* -x1 + x2 + x3 */
	xdiff[0] = -x[1] + x[2] + x[4];
	/* x1 - x2 + x3 */
	xdiff[1] = x[0] - x[3] + x[4];
	/* x1 + x2 - x3 */
	xdiff[2] = x[0] + x[2] - x[5];
	
	/* use pterms & nterms to compute tomsu */
	pterms[0] = 2.0*(x[0]*(x[4] + x[2]) + x[2]*x[4]);
	nterms[0] = x[0]*x[0] + x[2]*x[2] + x[4]*x[4];
	
	intprod[0] = xdiff[0]*x[0];
	intprod[1] = xdiff[1]*x[2];
	intprod[2] = xdiff[2]*x[4];
	
	prod[0] = intprod[1]*intprod[2];
	prod[2] = intprod[0]*intprod[2];
	prod[4] = intprod[0]*intprod[1];
	
	ROUND_UP;
	xdiff[0] = -x[0] + x[3] + x[5];
	xdiff[1] = x[1] - x[2] + x[5];
	xdiff[2] = x[1] + x[3] - x[4];
	
	pterms[1] = 2.0*(x[1]*(x[5] + x[3]) + x[3]*x[5]);
	nterms[1] = x[1]*x[1] + x[3]*x[3] + x[5]*x[5];
	foo = pterms[1] - nterms[0];
	uval2[1] = foo*foo;

	intprod[0] = xdiff[0]*x[1];
	intprod[1] = xdiff[1]*x[3];
	intprod[2] = xdiff[2]*x[5];
	
	prod[1] = intprod[1]*intprod[2];
	prod[3] = intprod[0]*intprod[2];
	prod[5] = intprod[0]*intprod[1];
	
	ROUND_DOWN;
	foo = pterms[0] - nterms[1];
	if( foo < 0.0 )
		foo = 0.0;
	uval2[0] = foo*foo;
	foo = 1.0/uval2[1];
	for( i=0; i<6; i+=2 ) {
		part[i] = foo*prod[i];
	}
	ROUND_UP;
	foo = 1.0/uval2[0];
	for( i=1; i<6; i+=2 ) {
		part[i] = foo*prod[i];
	}
}


/* Convention: ab[4], xy[4], where x = a*a, y = b*b. */

/* alpha_a = a/(2-a^2)*sqrt( (2-b^2)/(b^2-a^2) ) */
/* increasing in a, decreasing in b */
/* Note that usually y > x, but that this need not
always be true (for an interval).  Gotta be careful
about that. */
/* alpha_b = -b/(sqrt( b^2 - a^2 )*sqrt( 2 - b^2 ) ) */
void alpha_ab( double xy[4], double part[4] )
{
	double num[2], den[2], foo;
	double p1, p2, bden[2], temp[2];
	
	ROUND_DOWN;
	num[0] = xy[0]*(2.0 - xy[3]);
	foo = 2.0 - xy[1];
	den[0] = foo*foo*(xy[2] - xy[1]);
	
	p1 = xy[2] - xy[1];
	p2 = 2.0 - xy[3];
	bden[0] = p1*p2;
	
	ROUND_UP;
	num[1] = xy[1]*(2.0 - xy[2]);
	foo = 2.0 - xy[0];
	den[1] = foo*foo*(xy[3] - xy[0]);
	foo = num[1]/den[0];
	part[1] = sqrt( foo );
	
	p1 = xy[3] - xy[0];
	p2 = 2.0 - xy[2];
	bden[1] = p1*p2;
	temp[1] = sqrt( xy[3]/bden[0] );
	
	ROUND_DOWN;
	foo = num[0]/den[1];
	part[0] = sqrt( foo );

	temp[0] = sqrt( xy[2]/bden[1] );
	part[2] = - temp[1];
	part[3] = - temp[0];
}


/* wedgevol_a = ((2-3*a^2)*alpha +a*(2-a^2)*alpha_a)/6 */
/* wedgevol_b = a*(2-a^2)*alpha_b/6 */
/* coefficients are decreasing in a */
/* c1 < 0, c2 > 0 */
void wedgevol_ab( double ab[4], double xy[4], double al[2], 
	double al_ab[4], double part[4] )
{
	double c1[2], c2[2];
	
	ROUND_DOWN;
	c1[0] = ONE_6_HI*(2.0 - 3.0*xy[1]);
	c2[0] = ONE_6_LO*ab[1]*(2.0 - xy[1]);
	part[0] = c1[0]*al[1] + c2[0]*al_ab[0];
	ROUND_UP;
	c1[1] = ONE_6_LO*(2.0 - 3.0*xy[0]);
	c2[1] = ONE_6_HI*ab[0]*(2.0 - xy[0]);
	part[1] = c1[1]*al[0] + c2[1]*al_ab[1];
	part[3] = c2[0]*al_ab[3];
	ROUND_DOWN;
	part[2] = c2[1]*al_ab[2];
}


/* ws_a = al_a*(1 - a/sqrt(2)) - al/sqrt(2) */
/* ws_b = al_b*(1 - a/sqrt(2)) */
/* note (1 - a/sqrt(2) ) > 0 */
/* also, al_a > 0, al_b < 0 */
void wedgesol_ab( double ab[4], double al[2],
	double al_ab[4], double part[4] )
{
	double c1[2];
	
	ROUND_DOWN;
	c1[0] = 1.0 - ab[1]*ONE_SQRT2_HI;
	ROUND_UP;
	c1[1] = 1.0 - ab[0]*ONE_SQRT2_LO;
	
	part[1] = al_ab[1]*c1[1] - al[0]*ONE_SQRT2_LO;
	part[3] = al_ab[3]*c1[0];
	ROUND_DOWN;
	part[0] = al_ab[0]*c1[0] - al[1]*ONE_SQRT2_HI;
	part[2] = al_ab[2]*c1[1];
}


/* rv_a = sqrt( (2-b^2)/(b^2-a^2) )*(b^2 - 2*a^2)/6 */
/* rv_b = a*b*( 2 + a^2 - 2*b^2 )/sqrt( 
			(b^2 - a^2)*(2 - b^2) )/6 */
/* in reduced terms, these become */
/* rv_a = sqrt( (2-y)/(y-x) )*(y-2*x)/6 */
/* rv_b =  a*b*( 2+x-2*y )/sqrt( (y-x)*(2-y) )/6 */
/* Be careful about y > x (not always true).  */
void rogvol_ab( double ab[4], double xy[4], double part[4] )
{
	double two_y[2], y_x[2], y_2x[2], twopx_2y[2];
	double foo, quot[2], den[2], atb[2], temp[2];
	
	ROUND_DOWN;
	two_y[0] = 2.0 - xy[3];
	y_x[0] = xy[2] - xy[1];
	y_2x[0] = xy[2] - 2.0*xy[1];
	twopx_2y[0] = 2.0 + xy[0] - 2.0*xy[3];
	foo = y_x[0]*two_y[0];
	den[0] = 6.0*sqrt( foo );
	atb[0] = ab[0]*ab[2];
	
	ROUND_UP;
	two_y[1] = 2.0 - xy[2];
	y_x[1] = xy[3] - xy[0];
	y_2x[1] = xy[3] - 2.0*xy[0];
	twopx_2y[1] = 2.0 + xy[1] - 2.0*xy[2];
	foo = y_x[1]*two_y[1];
	den[1] = 6.0*sqrt( foo );
	atb[1] = ab[1]*ab[3];
	foo = two_y[1]/y_x[0];
	quot[1] = ONE_6_HI*sqrt( foo );
	temp[1] = atb[1]/den[0];
	
	ROUND_DOWN;
	foo = two_y[0]/y_x[1];
	quot[0] = ONE_6_LO*sqrt( foo );
	temp[0] = atb[0]/den[1];
	
	i_mult( quot, y_2x, part );
	i_mult( temp, twopx_2y, part + 2 );
}


/* rs_a=-sqrt( (z-y)/(y-x) )/(a+c) */
/* rs_b=(a*c - y)/(b*sqrt( (z-y)*(y-x) )) */
/* old rs_b=(2*a + c*x - y*(a+c))/(b*(a+c)*sqrt( (2-y)*(y-x) )) */
/* is the old formula for rs_b correct?  yes, but more */
/* complicated than it needs to be */
void rogsol_ab( double ab[4], double xy[4], double part[4] )
{
	int i;
	double abc[6], xyz[6];
	double z_y[2], y_x[2], ac_y[2], apc[2];
	double den[2], num[2];
	double temp;
	
	for( i=0; i<4; i++ ) {
		abc[i] = ab[i];
		xyz[i] = xy[i];
	}
	abc[4] = SQRT2_LO;
	abc[5] = SQRT2_HI;
	xyz[4] = 2.0;
	xyz[5] = 2.0;
	
	ROUND_DOWN;
	temp = xyz[4] - xyz[3];
	if( temp < 0.0 )
		temp = 0.0;
	z_y[0] = temp;
	temp = xyz[2] - xyz[1];
	if( temp < 0.0 )
		temp = 0.0;
	y_x[0] = temp;
	ac_y[0] = abc[0]*abc[4] - xyz[3];
	apc[0] = abc[0] + abc[4];
	den[0] = abc[2]*sqrt( z_y[0]*y_x[0] );
	
	ROUND_UP;
	z_y[1] = xyz[5] - xyz[2];
	y_x[1] = xyz[3] - xyz[0];
	ac_y[1] = abc[1]*abc[5] - xyz[2];
	apc[1] = abc[1] + abc[5];
	den[1] = abc[3]*sqrt( z_y[1]*y_x[1] );
	num[1] = sqrt( z_y[1]/y_x[0] )/apc[0];
	part[0] = -num[1];
	i_div( ac_y, den, part + 2 );
	
	ROUND_DOWN;
	num[0] = sqrt( z_y[0]/y_x[1] )/apc[1];
	part[1] = -num[0];
}


void old_rogsol_ab( double ab[4], double xy[4], double part[4] )
{
	double two_y[2], y_x[2];
	double sqrt2_y[2], sqrty_x[2], apc[2];
	double dena[2], denb[2], numb[2], pterms[2], nterms[2];
	double c[2], p1;
	
	c[0] = SQRT2_LO;
	c[1] = SQRT2_HI;
	
	ROUND_DOWN;
	two_y[0] = 2.0 - xy[3];
	y_x[0] = xy[2] - xy[1];
	sqrt2_y[0] = sqrt( two_y[0] );
	sqrty_x[0] = sqrt( y_x[0] );
	apc[0] = ab[0] + c[0];
	dena[0] = apc[0]*sqrty_x[0];
	denb[0] = ab[2]*dena[0]*sqrt2_y[0];
	pterms[0] = 2.0*ab[0] + xy[0]*c[0];
	nterms[0] = xy[2]*(ab[0] + c[0]);
	
	ROUND_UP;
	two_y[1] = 2.0 - xy[2];
	y_x[1] = xy[3] - xy[0];
	sqrt2_y[1] = sqrt( two_y[1] );
	sqrty_x[1] = sqrt( y_x[1] );
	apc[1] = ab[1] + c[1];
	dena[1] = apc[1]*sqrty_x[1];
	denb[1] = ab[3]*dena[1]*sqrt2_y[1];
	pterms[1] = 2.0*ab[1] + xy[1]*c[1];
	nterms[1] = xy[3]*(ab[1] + c[1]);
	numb[1] = pterms[1] - nterms[0];
	
	apc[1] = sqrt2_y[1]/dena[0];
	ROUND_DOWN;
	apc[0] = sqrt2_y[0]/dena[1];
	numb[0] = pterms[0] - nterms[1];
	i_div( numb, denb, part + 2 );
	
	part[0] = -apc[1];
	part[1] = -apc[0];
	
	p1 = two_y[1];
}


void wedge_pars( double y[12], double x[12],
	double vol_pars[12], double sol_pars[12] )
{
	double al[2], al_pars[12], foo;
	double wv_al[2], wv_coe[2];
	double ws_al[2], ws_coe[2];
	int i;
	
	/* printf("Starting wedge_pars.\n"); */
	s_dih( x, al );
	s_dih_partials( x, y, al_pars );
	
	ROUND_DOWN;
	foo = -ONE_2SQRT2_HI;
	ws_coe[0] = 1.0 + foo*y[1];
	ws_al[0] = foo*al[1];
	wv_coe[0] = ONE_48_LO*y[1]*(8.0 - x[1]);
	wv_al[0] = ONE_48_HI*(8.0 + (-3.0)*x[1])*al[1];
	ROUND_UP;
	foo = -ONE_2SQRT2_LO;
	ws_coe[1] = 1.0 + foo*y[0];
	ws_al[1] = foo*al[0];
	wv_coe[1] = ONE_48_HI*y[0]*(8.0 - x[0]);
	wv_al[1] = ONE_48_LO*(8.0 + (-3.0)*x[0])*al[0];
	
	/* Do sol_pars */
	for( i=0; i<12; i+=2 ) {
		i_mult( ws_coe, al_pars + i, sol_pars + i);
	}
	ROUND_DOWN;
	sol_pars[0] += ws_al[0];
	ROUND_UP;
	sol_pars[1] += ws_al[1];
	
	/* Do vol_pars */
	for( i=0; i<12; i+=2 ) {
		i_mult( wv_coe, al_pars + i, vol_pars + i);
	}
	ROUND_DOWN;
	vol_pars[0] += wv_al[0];
	ROUND_UP;
	vol_pars[1] += wv_al[1];
}


/* Takes sph_pars as an argument. */
/* What a mess.  This is going to be fun to debug. */
void trunc_pars( double y[12], double x[12], 
	double sph_pars[12], double part[12] )
{
	double xy[4], ab[4], facey[6], facex[6];
	double al[2], al_ab[4];
	double vol_pars[12], sol_pars[12];
	double bpars[6];
	double temp[4], temp2[4], voltemp[4], soltemp[4];
	double vbtemp[2], sbtemp[2];
	double wed_vpars[12], wed_spars[12];
	double wv_temp[12], ws_temp[12];
	double xp[12], yp[12], foo;
	int i, j, k, kp, m, mp;
	
	int b_list[2][3] = {{2, 4, 6}, {0, 4, 8}};
	int swap_list[2][6] = {{2, 0, 4, 8, 6, 10},
												{4, 2, 0, 10, 8, 6}};
	
	for( i=0; i<12; i++ ) {
		vol_pars[i] = 0.0;
		sol_pars[i] = 0.0;
		wed_vpars[i] = 0.0;
		wed_spars[i] = 0.0;
	}

	for( i=0; i<2; i++ ) {
		/* Get correct face */
		for( j=0; j<3; j++ ) {
			k = 2*j;
			kp = k + 1;
			m = b_list[i][j];
			mp = m + 1;
			facey[k] = y[ m ];
			facey[kp] = y[ mp ];
			facex[k] = x[ m ];
			facex[kp] = x[ mp ];
		}
		s_crad3x2(  facex, xy + 2 );
		ROUND_DOWN;
		ab[2] = sqrt( xy[2] );
		ROUND_UP;
		ab[3] = sqrt( xy[3] );
		s_crad3len_pars( facey, facex, bpars );
		/* crad3len_pars are guaranteed to be positive */

		for( j=0; j<2; j++ ) {
			vbtemp[j] = 0.0;
			sbtemp[j] = 0.0;
		}
		for( j=0; j<2; j++ ) {
			m = b_list[i][j];
			mp = m + 1;
			ab[0] = 0.5*y[m];
			ab[1] = 0.5*y[mp];
			xy[0] = 0.25*x[m];
			xy[1] = 0.25*x[mp];
			
			i_dih_rog( xy, al );	/* cache this? */
			alpha_ab( xy, al_ab );
			wedgevol_ab( ab, xy, al, al_ab, temp );
			wedgesol_ab( ab, al, al_ab, temp2 );
			rogvol_ab( ab, xy, voltemp );
			rogsol_ab( ab, xy, soltemp );
			ROUND_DOWN;
			/* Compute rog - wedge */
			for( k=0; k<4; k+=2 ) {
				kp = k + 1;
				voltemp[k] -= temp[kp];
				soltemp[k] -= temp2[kp];
			}
			/* Accumulate b_pars: */
			vbtemp[0] += voltemp[2];
			sbtemp[0] += soltemp[2];
			ROUND_UP;
			for( k=1; k<4; k+=2 ) {
				kp = k - 1;
				voltemp[k] -= temp[kp];
				soltemp[k] -= temp2[kp];
			}
			/* Accumulate b_pars: */
			vbtemp[1] += voltemp[3];
			sbtemp[1] += soltemp[3];
			/* Change a_pars to y_pars and remap */
			ROUND_DOWN;
			vol_pars[m] += 0.5*voltemp[0];
			sol_pars[m] += 0.5*soltemp[0];
			ROUND_UP;
			vol_pars[mp] += 0.5*voltemp[1];
			sol_pars[mp] += 0.5*soltemp[1];
		}
		/* Change b_pars to y_pars and remap */
		for( j=0; j<3; j++ ) {
			k = 2*j;
			m = b_list[i][j];
			mp = m + 1;
			i_mult( vbtemp, bpars + k, temp );
			i_mult( sbtemp, bpars + k, temp2 );
			ROUND_DOWN;
			vol_pars[m] += temp[0];
			sol_pars[m] += temp2[0];
			ROUND_UP;
			vol_pars[mp] += temp[1];
			sol_pars[mp] += temp2[1];
		}
	}
	
	ROUND_NEAR;
	/*
	for( i=0; i<12; i+=2 ) {
		printf("vol_pars = (%.18f, %.18f)\n", 
			vol_pars[i], vol_pars[i+1]);
		foo = vol_pars[i+1] - vol_pars[i];
		if( foo < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
	}
	for( i=0; i<12; i+=2 ) {
		printf("sol_pars = (%.18f, %.18f)\n", 
			sol_pars[i], sol_pars[i+1]);
		foo = sol_pars[i+1] - sol_pars[i];
		if( foo < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
	}
	*/
	
	
	/* That should do it for rogers and rogers-wedges, but
	we still need to do the full wedges.  */
	
	wedge_pars( y, x, wed_vpars, wed_spars );
	
	/* edge 2 (2) */	/*	(2 1 3 5 4 6) -> (2 0 4 8 6 10)	*/
	/* edge 3 (4) */	/*	(3 2 1 6 5 4) -> (4 2 0 10 8 6)	*/
	for( j=0; j<2; j++ ) {
		for( i=0; i<6; i++ ) {
			k = 2*i;
			kp = k + 1;
			m = swap_list[j][i];
			mp = m + 1;
			xp[ k ] = x[ m ];
			yp[ k ] = y[ m ];
			xp[ kp ] = x[ mp ];
			yp[ kp ] = y[ mp ];
		}
		wedge_pars( yp, xp, wv_temp, ws_temp );
		
		/* Now swap back and accumulate. */
		ROUND_DOWN;
		for( i=0; i<6; i++ ) {
			k = 2*i;
			m = swap_list[j][i];
			wed_vpars[m] += wv_temp[k];
			wed_spars[m] += ws_temp[k];
		}
		ROUND_UP;
		for( i=0; i<6; i++ ) {
			k = 2*i + 1;
			m = swap_list[j][i] + 1;
			wed_vpars[m] += wv_temp[k];
			wed_spars[m] += ws_temp[k];
		}
	}
	/*
	ROUND_NEAR;
	for( i=0; i<12; i+=2 ) {
		printf("wed_vpars = (%.18f, %.18f)\n", 
			wed_vpars[i], wed_vpars[i+1]);
		foo = wed_vpars[i+1] - wed_vpars[i];
		if( foo < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
	}
	for( i=0; i<12; i+=2 ) {
		printf("wed_spars = (%.18f, %.18f)\n", 
			wed_spars[i], wed_spars[i+1]);
		foo = wed_spars[i+1] - wed_spars[i];
		if( foo < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
	}
	ROUND_UP;
	*/
	/* Now combine everything. */
	ROUND_DOWN;
	for( i=0; i<12; i+=2 ) {
		wv_temp[i] = wed_vpars[i] + vol_pars[i];
		ws_temp[i] = wed_spars[i] + sol_pars[i];
	}
	ROUND_UP;
	for( i=1; i<12; i+=2 ) {
		wv_temp[i] = wed_vpars[i] + vol_pars[i];
		ws_temp[i] = wed_spars[i] + sol_pars[i];
	}
	/*
	ROUND_NEAR;
	for( i=0; i<12; i+=2 ) {
		printf("wv_temp = (%.18f, %.18f)\n", 
			wv_temp[i], wv_temp[i+1]);
		foo = wv_temp[i+1] - wv_temp[i];
		if( foo < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
	}
	for( i=0; i<12; i+=2 ) {
		printf("ws_temp = (%.18f, %.18f)\n", 
			ws_temp[i], ws_temp[i+1]);
		foo = ws_temp[i+1] - ws_temp[i];
		if( foo < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
	}
	ROUND_UP;
	*/
	/* Clean this up later . . . */
	temp[1] = TWOSQRT2_HI*DOCT_HI*ONE_3_HI;
	ROUND_DOWN;
	temp[0] = TWOSQRT2_LO*DOCT_LO*ONE_3_LO;
	temp2[0] = ONE_3_LO - temp[1];
	ROUND_UP;
	temp2[1] = ONE_3_HI - temp[0];
	/* Note temp > 0, temp2 < 0 */
	/* Fix this. */
	
	/* Generate constant multiples . . . */
	for( i=1; i<12; i+=2 ) {
		foo = ws_temp[i];
		if( foo > 0.0 )
			sol_pars[i] = foo*temp[1];
		else
			sol_pars[i] = foo*temp[0];
	}
	for( i=1; i<12; i+=2 ) {
		foo = wv_temp[i];
		if( foo > 0.0 )
			vol_pars[i] = foo*DOCT_HI;
		else
			vol_pars[i] = foo*DOCT_LO;
	}
	ROUND_DOWN;
	for( i=0; i<12; i+=2 ) {
		foo = ws_temp[i];
		if( foo > 0.0 )
			sol_pars[i] = foo*temp[0];
		else
			sol_pars[i] = foo*temp[1];
	}
	for( i=0; i<12; i+=2 ) {
		foo = wv_temp[i];
		if( foo > 0.0 )
			vol_pars[i] = foo*DOCT_LO;
		else
			vol_pars[i] = foo*DOCT_HI;
	}
	/*
	printf("Constant multiples . . .\n");
	ROUND_NEAR;
	for( i=0; i<12; i+=2 ) {
		printf("vol_pars = (%.18f, %.18f)\n", 
			vol_pars[i], vol_pars[i+1]);
		foo = vol_pars[i+1] - vol_pars[i];
		if( foo < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
	}
	for( i=0; i<12; i+=2 ) {
		printf("sol_pars = (%.18f, %.18f)\n", 
			sol_pars[i], sol_pars[i+1]);
		foo = sol_pars[i+1] - sol_pars[i];
		if( foo < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
	}
	ROUND_DOWN;
	*/
	for( i=0; i<12; i+=2 ) {
		j = i + 1;
		wv_temp[i] = sol_pars[i] - vol_pars[j];
	}
	ROUND_UP;
	for( i=1; i<12; i+=2 ) {
		j = i - 1;
		wv_temp[i] = sol_pars[i] - vol_pars[j];
	}
	/*
	printf("Differences . . .\n");
	ROUND_NEAR;
	for( i=0; i<12; i+=2 ) {
		printf("wv_temp = (%.18f, %.18f)\n", 
			wv_temp[i], wv_temp[i+1]);
		foo = wv_temp[i+1] - wv_temp[i];
		if( foo < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
	}
	ROUND_UP;
	*/
	for( i=1; i<12; i+=2 ) {
		j = i - 1;
		foo = sph_pars[j];
		if( foo > 0.0 )
			sol_pars[i] = foo*temp2[1];
		else
			sol_pars[i] = foo*temp2[0];
	}
	ROUND_DOWN;
	for( i=0; i<12; i+=2 ) {
		j = i + 1;
		foo = sph_pars[j];
		if( foo > 0.0 )
			sol_pars[i] = foo*temp2[0];
		else
			sol_pars[i] = foo*temp2[1];
	}
	/*
	printf("sph_pars . . .\n");
	ROUND_NEAR;
	for( i=0; i<12; i+=2 ) {
		printf("sph_pars = (%.18f, %.18f)\n", 
			sph_pars[i], sph_pars[i+1]);
		foo = sph_pars[i+1] - sph_pars[i];
		if( foo < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
	}
	ROUND_DOWN;
	*/
	/*
	printf("sph_pars mult . . .\n");
	ROUND_NEAR;
	printf("temp = [%.18f, %.18f]\n", temp[0], temp[1]);
	foo = temp[1] - temp[0];
	if( foo < 0.0 )
		printf("Uh oh.  Backwards interval.\n");
	printf("temp2 = [%.18f, %.18f]\n", temp2[0], temp2[1]);
	foo = temp2[1] - temp2[0];
	if( foo < 0.0 )
		printf("Uh oh.  Backwards interval.\n");
	for( i=0; i<12; i+=2 ) {
		printf("sol_pars = (%.18f, %.18f)\n", 
			sol_pars[i], sol_pars[i+1]);
		foo = sol_pars[i+1] - sol_pars[i];
		if( foo < 0.0 )
			printf("Uh oh.  Backwards interval.\n");
	}
	ROUND_DOWN;
	*/
	for( i=0; i<12; i+=2 ) {
		part[i] = 4.0*(sol_pars[i] + wv_temp[i]);
	}
	ROUND_UP;
	for( i=1; i<12; i+=2 ) {
		part[i] = 4.0*(sol_pars[i] + wv_temp[i]);
	}
}


void approx_trunc_pars( double y[12], double part[12] )
{
	double x[12], yp[12], xp[12], ym[12], xm[12], delta[2];
	double sqrtdelta[2], sqrtdeltap[2], sqrtdeltam[2];
	double vval[2], vvalp[2], vvalm[2];
	double eps, val, valp, valm, deriv, deriv2;
	int i, j, jp;
	
	eps = 1.0e-7;
	ROUND_DOWN;
	for( i=0; i<12; i+=2 ) {
		x[i] = y[i]*y[i];
	}
	ROUND_UP;
	for( i=1; i<12; i+=2 ) {
		x[i] = y[i]*y[i];
	}
	i_bigdelta( x, delta );
	i_sqrt( delta, sqrtdelta );
	i_vor_trunc( y, x, sqrtdelta, vval );
	ROUND_NEAR;
	val = 0.5*( vval[0] + vval[1] );
	for( i=0; i<12; i+=2 ) {
		for( j=0; j<12; j++ ) {
			yp[j] = y[j];
			xp[j] = x[j];
			ym[j] = y[j];
			xm[j] = x[j];
		}
		j = i;
		jp = i + 1;
		ROUND_DOWN;
		yp[j] = y[j] + eps;
		xp[j] = yp[j]*yp[j];
		ym[j] = y[j] - eps;
		xm[j] = ym[j]*ym[j];
		ROUND_UP;
		yp[jp] = y[jp] + eps;
		xp[jp] = yp[jp]*yp[jp];
		ym[jp] = y[jp] - eps;
		xm[jp] = ym[jp]*ym[jp];
		
		i_bigdelta( xp, delta );
		i_sqrt( delta, sqrtdeltap );
		i_vor_trunc( yp, xp, sqrtdeltap, vvalp );
		
		i_bigdelta( xm, delta );
		i_sqrt( delta, sqrtdeltam );
		i_vor_trunc( ym, xm, sqrtdeltam, vvalm );
		
		ROUND_NEAR;
		valp = 0.5*( vvalp[0] + vvalp[1] );
		valm = 0.5*( vvalm[0] + vvalm[1] );
		deriv = (valp - val)/eps;
		deriv2 = (valp - valm)/(2.0*eps);
		/* printf("estimate: %g\n", fabs(deriv - deriv2) ); */
		part[j] = deriv2;
		part[jp] = deriv2;
	}
}


void approx_rog_wed_pars( double y[12], double vpart[12],
	double spart[12] )
{
	double x[12], yp[12], xp[12], ym[12], xm[12];
	double vvalp[2], vvalm[2];
	double svalp[2], svalm[2];
	double eps, valp, valm, deriv2;
	double salp, salm, sderiv2;
	int i, j, jp;
	
	eps = 1.0e-7;
	ROUND_DOWN;
	for( i=0; i<12; i+=2 ) {
		x[i] = y[i]*y[i];
	}
	ROUND_UP;
	for( i=1; i<12; i+=2 ) {
		x[i] = y[i]*y[i];
	}
	for( i=0; i<12; i+=2 ) {
		for( j=0; j<12; j++ ) {
			yp[j] = y[j];
			xp[j] = x[j];
			ym[j] = y[j];
			xm[j] = x[j];
		}
		j = i;
		jp = i + 1;
		ROUND_DOWN;
		yp[j] = y[j] + eps;
		xp[j] = yp[j]*yp[j];
		ym[j] = y[j] - eps;
		xm[j] = ym[j]*ym[j];
		ROUND_UP;
		yp[jp] = y[jp] + eps;
		xp[jp] = yp[jp]*yp[jp];
		ym[jp] = y[jp] - eps;
		xm[jp] = ym[jp]*ym[jp];
		
		test_rog_wed( yp, xp, vvalp, svalp );
		
		test_rog_wed( ym, xm, vvalm, svalm );
		
		ROUND_NEAR;
		valp = 0.5*( vvalp[0] + vvalp[1] );
		salp = 0.5*( svalp[0] + svalp[1] );
		valm = 0.5*( vvalm[0] + vvalm[1] );
		salm = 0.5*( svalm[0] + svalm[1] );
		deriv2 = (valp - valm)/(2.0*eps);
		sderiv2 = (salp - salm)/(2.0*eps);
		vpart[j] = deriv2;
		vpart[jp] = deriv2;
		spart[j] = sderiv2;
		spart[jp] = sderiv2;
	}
}


void approx_fullwed_pars( double y[12], double vpart[12],
	double spart[12] )
{
	double x[12], yp[12], xp[12], ym[12], xm[12];
	double vvalp[2], vvalm[2];
	double svalp[2], svalm[2];
	double eps, valp, valm, deriv2;
	double salp, salm, sderiv2;
	int i, j, jp;
	
	eps = 1.0e-7;
	ROUND_DOWN;
	for( i=0; i<12; i+=2 ) {
		x[i] = y[i]*y[i];
	}
	ROUND_UP;
	for( i=1; i<12; i+=2 ) {
		x[i] = y[i]*y[i];
	}
	for( i=0; i<12; i+=2 ) {
		for( j=0; j<12; j++ ) {
			yp[j] = y[j];
			xp[j] = x[j];
			ym[j] = y[j];
			xm[j] = x[j];
		}
		j = i;
		jp = i + 1;
		ROUND_DOWN;
		yp[j] = y[j] + eps;
		xp[j] = yp[j]*yp[j];
		ym[j] = y[j] - eps;
		xm[j] = ym[j]*ym[j];
		ROUND_UP;
		yp[jp] = y[jp] + eps;
		xp[jp] = yp[jp]*yp[jp];
		ym[jp] = y[jp] - eps;
		xm[jp] = ym[jp]*ym[jp];
		
		test_fullwed( yp, xp, vvalp, svalp );
		
		test_fullwed( ym, xm, vvalm, svalm );
		
		ROUND_NEAR;
		valp = 0.5*( vvalp[0] + vvalp[1] );
		salp = 0.5*( svalp[0] + svalp[1] );
		valm = 0.5*( vvalm[0] + vvalm[1] );
		salm = 0.5*( svalm[0] + svalm[1] );
		deriv2 = (valp - valm)/(2.0*eps);
		sderiv2 = (salp - salm)/(2.0*eps);
		vpart[j] = deriv2;
		vpart[jp] = deriv2;
		spart[j] = sderiv2;
		spart[jp] = sderiv2;
	}
}


void test_ab_pars( double ab[4] )
{
	double xyz[6], abp[4], abm[4], xyzp[6], xyzm[6];
	double eps, temp2[2];
	double ivalp[2], ivalm[2];
	double valp, valm, diff;
	double alp[2], alm[2];
	int i, ip;
	
	xyz[4] = 2.0;
	xyz[5] = 2.0;
	xyzp[4] = 2.0;
	xyzp[5] = 2.0;
	xyzm[4] = 2.0;
	xyzm[5] = 2.0;
	
	ROUND_DOWN;
	xyz[0] = ab[0]*ab[0];
	xyz[2] = ab[2]*ab[2];
	ROUND_UP;
	xyz[1] = ab[1]*ab[1];
	xyz[3] = ab[3]*ab[3];

	eps = 1.0e-7;
	
	for( i=0; i<4; i+=2 ) {
		for( ip=0; ip<4; ip++ ) {
			abp[ip] = ab[ip];
			abm[ip] = ab[ip];
			xyzp[ip] = xyz[ip];
			xyzm[ip] = xyz[ip];
		}
		ip = i + 1;
		ROUND_DOWN;
		abp[i] = ab[i] + eps;
		abm[i] = ab[i] - eps;
		xyzp[i] = abp[i]*abp[i];
		xyzm[i] = abm[i]*abm[i];
		ROUND_UP;
		abp[ip] = ab[ip] + eps;
		abm[ip] = ab[ip] - eps;
		xyzp[ip] = abp[ip]*abp[ip];
		xyzm[ip] = abm[ip]*abm[ip];

		i_dih_rog( xyzp, alp );
		i_dih_rog( xyzm, alm );
		ROUND_NEAR;
		valp = 0.5*(alp[0] + alp[1]);
		valm = 0.5*(alm[0] + alm[1]);
		diff = (valp - valm)/(2.0*eps);
		printf("dih_rog_part = %.18g\n", diff);
		
		i_rogersvol2( xyzp, temp2 );
		I_SQRT( temp2, ivalp );
		i_rogersvol2( xyzm, temp2 );
		I_SQRT( temp2, ivalm );
		ROUND_NEAR;
		valp = 0.5*(ivalp[0] + ivalp[1]);
		valm = 0.5*(ivalm[0] + ivalm[1]);
		diff = (valp - valm)/(2.0*eps);
		printf("rogvol_part = %.18g\n", diff);
		
		i_rog_sph( abp, ivalp );
		i_rog_sph( abm, ivalm );
		ROUND_NEAR;
		valp = 0.5*(ivalp[0] + ivalp[1]);
		valm = 0.5*(ivalm[0] + ivalm[1]);
		diff = (valp - valm)/(2.0*eps);
		printf("rogsol_part = %.18g\n", diff);
		
		i_wedge_vol( alp, abp, ivalp );
		i_wedge_vol( alm, abm, ivalm );
		ROUND_NEAR;
		valp = 0.5*(ivalp[0] + ivalp[1]);
		valm = 0.5*(ivalm[0] + ivalm[1]);
		diff = (valp - valm)/(2.0*eps);
		printf("wedgevol_part = %.18g\n", diff);

		i_wedge_sph( alp, abp, ivalp );
		i_wedge_sph( alm, abm, ivalm );
		ROUND_NEAR;
		valp = 0.5*(ivalp[0] + ivalp[1]);
		valm = 0.5*(ivalm[0] + ivalm[1]);
		diff = (valp - valm)/(2.0*eps);
		printf("wedgesol_part = %.18g\n", diff);
	}
}


void i_sign_dihpars( double x[12], double delta[2],
	double delta_i[12], double out[12] )
{
	double delta14[2], p1[2], p2[2], p3[2];
	
	ROUND_DOWN;
	delta14[0] = -2.0*x[1] + x[2] + x[4] - 2.0*x[7] + 
							x[8] + x[10];
	ROUND_UP;
	delta14[1] = -2.0*x[0] + x[3] + x[5] - 2.0*x[6] + 
							x[9] + x[11];
	
	i_mult( delta_i, delta_i + 6, p1 );
	i_mult( delta, delta14, p2 );
	ROUND_DOWN;
	p3[0] = p1[0] - 2.0*p2[1];
	ROUND_UP;
	p3[1] = p1[1] - 2.0*p2[0];
	i_mult( x, p3, p1 );
	i_mult( delta, delta_i + 6, p2 );
	ROUND_DOWN;
	out[0] = p1[0] + p2[0];
	ROUND_UP;
	out[1] = p1[1] + p2[1];
	
	out[2]  = -delta_i[5];
	out[3]  = -delta_i[4];
	out[4]  = -delta_i[3];
	out[5]  = -delta_i[2];
	out[6]  = 1.0;
	out[7]  = 1.0;
	out[8]  = -delta_i[11];
	out[9]  = -delta_i[10];
	out[10] = -delta_i[9];
	out[11] = -delta_i[8];
}


/*
{-2*x[1]*x[4]^2 + 2*x[2]*x[4]*x[5] + 2*x[3]*x[4]*x[6], 
  2*x[1]*x[4]*x[5] - 2*x[2]*x[5]^2 + 2*x[3]*x[5]*x[6], 
  2*x[1]*x[4]*x[6] + 2*x[2]*x[5]*x[6] - 2*x[3]*x[6]^2, 
  -2*x[1]^2*x[4] + 2*x[1]*x[2]*x[5] + 2*x[1]*x[3]*x[6], 
  2*x[1]*x[2]*x[4] - 2*x[2]^2*x[5] + 2*x[2]*x[3]*x[6], 
  2*x[1]*x[3]*x[4] + 2*x[2]*x[3]*x[5] - 2*x[3]^2*x[6]}
*/


void i_tomsrho_xpars( double x[12], double out[12] )
{
	int i, ip;
	double nterms[12], pterms[12], diff[12];
	
	ROUND_DOWN;
	pterms[0] = x[6]*(x[2]*x[8] + x[4]*x[10]);
	pterms[2] = x[8]*(x[0]*x[6] + x[4]*x[10]);
	pterms[4] = x[10]*(x[0]*x[6] + x[2]*x[8]);
	pterms[6] = x[0]*(x[2]*x[8] + x[4]*x[10]);
	pterms[8] = x[2]*(x[0]*x[6] + x[4]*x[10]);
	pterms[10] = x[4]*(x[0]*x[6] + x[2]*x[8]);
	
	nterms[0] = x[0]*x[6]*x[6];
	nterms[2] = x[2]*x[8]*x[8];
	nterms[4] = x[4]*x[10]*x[10];
	nterms[6] = x[0]*x[0]*x[6];
	nterms[8] = x[2]*x[2]*x[8];
	nterms[10] = x[4]*x[4]*x[10];
	ROUND_UP;
	pterms[0+1] = x[6+1]*(x[2+1]*x[8+1] + x[4+1]*x[10+1]);
	pterms[2+1] = x[8+1]*(x[0+1]*x[6+1] + x[4+1]*x[10+1]);
	pterms[4+1] = x[10+1]*(x[0+1]*x[6+1] + x[2+1]*x[8+1]);
	pterms[6+1] = x[0+1]*(x[2+1]*x[8+1] + x[4+1]*x[10+1]);
	pterms[8+1] = x[2+1]*(x[0+1]*x[6+1] + x[4+1]*x[10+1]);
	pterms[10+1] = x[4+1]*(x[0+1]*x[6+1] + x[2+1]*x[8+1]);
	
	nterms[0+1] = x[0+1]*x[6+1]*x[6+1];
	nterms[2+1] = x[2+1]*x[8+1]*x[8+1];
	nterms[4+1] = x[4+1]*x[10+1]*x[10+1];
	nterms[6+1] = x[0+1]*x[0+1]*x[6+1];
	nterms[8+1] = x[2+1]*x[2+1]*x[8+1];
	nterms[10+1] = x[4+1]*x[4+1]*x[10+1];
	
	for( i=0; i<12; i+=2 ) {
		ip = i + 1;
		diff[ip] = pterms[ip] - nterms[i];
	}
	ROUND_DOWN;
	for( i=0; i<12; i+=2 ) {
		ip = i + 1;
		diff[i] = pterms[i] - nterms[ip];
	}
	for( i=0; i<12; i++ ) {
		out[i] = 2.0*diff[i];
	}
}


void i_sign_crad2_xpars( double x[12], double delta[2], 
	double delta_partials[12], double out[12] )
{
	int i;
	double rho_xpars[12], rho_val[2];
	double t1[2], t2[2];
	
	i_tomsrho_xpars( x, rho_xpars );
	i_tomsrho( x, rho_val );
	
	/*
	ROUND_NEAR;
	printf("rho_xpars = \n");
	for( i=0; i<12; i+=2 ) {
		printf("[%.18f, %.18f]\n", rho_xpars[i], rho_xpars[i+1]);
	}
	*/
	
	for( i=0; i<12; i+=2 ) {
		i_mult( delta, rho_xpars + i, t1 );
		i_mult( rho_val, delta_partials + i, t2 );
		i_sub( t1, t2, out + i );
	}
	/*
	ROUND_NEAR;
	printf("signs (out) = \n");
	for( i=0; i<12; i+=2 ) {
		printf("[%.18f, %.18f]\n", out[i], out[i+1]);
	}
	*/
}

