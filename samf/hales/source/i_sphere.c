/* i_sphere.c, by Samuel Ferguson, (c) 1997. */

/* Sam's funky program for doing stuff with tetrahedrons 
and sphere packings */

/* This presentation of the routines for sphere packing 
computations is taken from the presentation in sphere.m, 
the collection of Mathematica routines.  Note that I'm 
now using C indexing, not Mathematica-type indexing.  
Be warned.  */

#include "system_headers.h"
#include "i_sphere.h"
#include "interval.h"
#include "i_bounds.h"
#include "macros.h"

#define SQUARE( x )		( x )*( x )
#define DEBUGTRUNC			0

/* Useful definitions */
#define TWOPI				6.283185307179586476925286766559
#define TWOSQRT2		2.8284271247461900976033774484193961571

/* These are constants which appear in the modified definition 
for small simplices. */
#define CIRCUMCUT		1.39			/* originally 2.78/2 */
#define EDGELENCUT	2.12			/* originally 2.12   */
#define LONGCUT			2.8284272		/* something less than 2 Sqrt[2] */

/* External variables */

extern double i_doct_const[2];
extern double i_sqrt2[2];

/* Routines */

/* Fix this. */
void i_bvol( double y[12], double sqrtdelta[2], 
	double out[2] )
{
	double yp[12], afvs[8], tp[2], sum[2], temp[2], frac[2];
	int i;
	
	ROUND_DOWN;
	tp[0] = 0.5*sqrtdelta[0];
	ROUND_UP;
	tp[1] = 0.5*sqrtdelta[1];
	
	i_afunc( y, afvs );
	for( i=0; i<2; i++ ) {
		yp[0 + i] = y[0 + i];
		yp[2 + i] = y[8 + i];
		yp[4 + i] = y[10 + i];
		yp[6 + i] = y[6 + i];
		yp[8 + i] = y[2 + i];
		yp[10 + i] = y[4 + i];
		}
	i_afunc( yp, afvs + 2 );
	for( i=0; i<2; i++ ) {
		yp[0 + i] = y[2 + i];
		yp[2 + i] = y[6 + i];
		/*	yp[4 + i] = y[10 + i];	*/
		yp[6 + i] = y[8 + i];
		yp[8 + i] = y[0 + i];
		/*	yp[10 + i] = y[4 + i];	*/
		}
	i_afunc( yp, afvs + 4 );
	for( i=0; i<2; i++ ) {
		yp[0 + i] = y[4 + i];
		/*	yp[2 + i] = y[6 + i];	*/
		yp[4 + i] = y[8 + i];
		yp[6 + i] = y[10 + i];
		/*	yp[8 + i] = y[0 + i];	*/
		yp[10 + i] = y[2 + i];
		}
	i_afunc( yp, afvs + 6 );
	
	sum[0] = 0.0;
	sum[1] = 0.0;
	for( i=0; i<8; i+=2 ) {
		i_div( tp, afvs + i, frac );
		i_atan( frac, temp );
		I_ADD( sum, temp, frac );
		sum[0] = frac[0];
		sum[1] = frac[1];
		}
		
	ROUND_DOWN;
	out[0] = sum[0]*TWO_3_LO;
	ROUND_UP;
	out[1] = sum[1]*TWO_3_HI;
}


/* min_a( y ) = a( y[0], y[2], y[4], y[7], y[9], y[11] ) */
/* max_a( y ) = a( y[1], y[3], y[5], y[6], y[8], y[10] ) */
void i_afunc( double y[12], double out[2] )
{
	double y2[12], temp1[2], temp2[2], temp3[2];
	double p1[2], p2[2], p3[2], y123[2];
	int i;

	ROUND_DOWN;
	for( i=0; i<12; i+=2 )
		y2[i] = y[i]*y[i];
	ROUND_UP;
	for( i=1; i<12; i+=2 )
		y2[i] = y[i]*y[i];
	ROUND_DOWN;
	p1[0] = y2[2] + y2[4] - y2[7];
	p2[0] = y2[0] + y2[4] - y2[9];
	p3[0] = y2[0] + y2[2] - y2[11];
	ROUND_UP;
	p1[1] = y2[3] + y2[5] - y2[6];
	p2[1] = y2[1] + y2[5] - y2[8];
	p3[1] = y2[1] + y2[3] - y2[10];
	ROUND_DOWN;
	y123[0] = y[0]*y[2]*y[4];
	ROUND_UP;
	y123[1] = y[1]*y[3]*y[5];
	I_MULT( y, p1, temp1 );
	I_MULT( y + 2, p2, temp2 );
	I_MULT( y + 4, p3, temp3 );
	ROUND_DOWN;
	out[0] = y123[0] + 0.5*(temp1[0] + temp2[0] + temp3[0]);
	if( out[0] < 0.0 )
		out[0] = 0.0;
	ROUND_UP;
	out[1] = y123[1] + 0.5*(temp1[1] + temp2[1] + temp3[1]);
}


void i_solid( double y[12], double sqrtdelta[2], 
	double out[2] )
{
	double max, temp, mina;
	
	temp = sqrtdelta[0];
	mina = 2.0*max_a( y );
	ROUND_DOWN;
	max = temp/mina;
	out[0] = 2.0*(atan( max ) - ATANERR);

	temp = sqrtdelta[1];
	mina = 2.0*min_a( y );
	ROUND_UP;
	max = temp/mina;
	out[1] = 2.0*(atan( max ) + ATANERR);
}


void i_tvol( double y[12], double out[2] )
{
	double x[12], temp1[2], temp2[2];
	int i;

	ROUND_DOWN;
	for( i=0; i<12; i+=2 )
		x[i] = y[i]*y[i];
	ROUND_UP;
	for( i=1; i<12; i+=2 )
		x[i] = y[i]*y[i];
	
	i_bigdelta( x, temp1 );
	I_SQRT( temp1, temp2 );
	
	ROUND_DOWN;
	out[0] = temp2[0]*ONE_12_LO;
	ROUND_UP;
	out[1] = temp2[1]*ONE_12_HI;
}


void i_gma( double y[12], double sqrtdelta[2], 
	double out[2] )
{
	double temp1[2], temp2[2];
	
	i_bvol( y, sqrtdelta, temp1 );
	/* ROUND_UP; */
	temp2[1] = sqrtdelta[1]*ONE_12_HI*i_doct_const[1];
	ROUND_DOWN;
	temp2[0] = sqrtdelta[0]*ONE_12_LO*i_doct_const[0];
	out[0] = temp1[0] - temp2[1];
	ROUND_UP;
	out[1] = temp1[1] - temp2[0];
}


/* i_dih computes the dihedral angle associated with
edge 1 */

void i_dih( double x[12], double out[2] )
{
	double num[2], den[2], pterms[2], nterms[2];
	
	ROUND_DOWN;
	pterms[0] = x[0]*(x[2] + x[4] + x[8] + x[10]) +
		x[2]*x[8] + x[4]*x[10];
	nterms[0] = x[0]*(x[0] + 2.0*x[6]) + x[2]*x[4] + x[8]*x[10];
	ROUND_UP;
	pterms[1] = x[1]*(x[3] + x[5] + x[9] + x[11]) +
		x[3]*x[9] + x[5]*x[11];
	nterms[1] = x[1]*(x[1] + 2.0*x[7]) + x[3]*x[5] + x[9]*x[11];
	num[1] = nterms[1] - pterms[0];
	ROUND_DOWN;
	num[0] = nterms[0] - pterms[1];
	
	i_bigdelta( x, nterms );
	i_mult( nterms, x, pterms );
	i_sqrt( pterms, nterms );
	den[0] = 2.0*nterms[0];
	den[1] = 2.0*nterms[1];
/*
	ROUND_NEAR;
	printf("num = [%.18f, %.18f]\n", num[0], num[1]);
	printf("den = [%.18f, %.18f]\n", den[0], den[1]);
*/
	if( den[0] > 0.0 ) {
		i_div( num, den, nterms );
		i_atan( nterms, pterms );
		ROUND_DOWN;
		out[0] = PI_2_LO + pterms[0];
		ROUND_UP;
		out[1] = PI_2_HI + pterms[1];
	} else {
		if( num[0] > 0.0 ) {
			ROUND_DOWN;
			pterms[0] = atan( num[0]/den[1] ) - ATANERR;
			out[0] = PI_2_LO + pterms[0];
			out[1] = PI_HI;
		} else if( num[1] < 0.0 ) {
			ROUND_UP;
			pterms[1] = atan( num[1]/den[1] ) + ATANERR;
			out[0] = 0.0;
			out[1] = PI_2_HI + pterms[1];
		} else { /* mixed, punt */
			out[0] = 0.0;
			out[1] = PI_HI;			
		}
	}
}


void i_dih_alt( double x[12], double delta[2], double out[2] )
{
	double num[2], den[2], pterms[2], nterms[2];
	
	ROUND_DOWN;
	pterms[0] = x[0]*(x[2] + x[4] + x[8] + x[10]) +
		x[2]*x[8] + x[4]*x[10];
	nterms[0] = x[0]*(x[0] + 2.0*x[6]) + x[2]*x[4] + x[8]*x[10];
	ROUND_UP;
	pterms[1] = x[1]*(x[3] + x[5] + x[9] + x[11]) +
		x[3]*x[9] + x[5]*x[11];
	nterms[1] = x[1]*(x[1] + 2.0*x[7]) + x[3]*x[5] + x[9]*x[11];
	num[1] = nterms[1] - pterms[0];
	ROUND_DOWN;
	num[0] = nterms[0] - pterms[1];
	
	i_mult( delta, x, pterms );
	i_sqrt( pterms, nterms );
	den[0] = 2.0*nterms[0];
	den[1] = 2.0*nterms[1];
	if( den[0] > 0.0 ) {
		i_div( num, den, nterms );
		i_atan( nterms, pterms );
		ROUND_DOWN;
		out[0] = PI_2_LO + pterms[0];
		ROUND_UP;
		out[1] = PI_2_HI + pterms[1];
	} else {
		if( num[0] > 0.0 ) {
			ROUND_DOWN;
			pterms[0] = atan( num[0]/den[1] ) - ATANERR;
			out[0] = PI_2_LO + pterms[0];
			out[1] = PI_HI;
		} else if( num[1] < 0.0 ) {
			ROUND_UP;
			pterms[1] = atan( num[1]/den[1] ) + ATANERR;
			out[0] = 0.0;
			out[1] = PI_2_HI + pterms[1];
		} else { /* mixed, punt */
			out[0] = 0.0;
			out[1] = PI_HI;			
		}
	}
}


void i_dih_old( double x[12], double out[2] )
{
	double xv[6], u1[2], u2[2];
	double temp[2], num[2], den[2], pterms[2], nterms[2];
	int i;
	
	ROUND_DOWN;
	pterms[0] = x[0]*(x[2] + x[4] + x[8] + x[10]) +
		x[2]*x[8] + x[4]*x[10];
	nterms[0] = x[0]*(x[0] + 2.0*x[6]) + x[2]*x[4] + x[8]*x[10];
	ROUND_UP;
	pterms[1] = x[1]*(x[3] + x[5] + x[9] + x[11]) +
		x[3]*x[9] + x[5]*x[11];
	nterms[1] = x[1]*(x[1] + 2.0*x[7]) + x[3]*x[5] + x[9]*x[11];
	I_SUB( pterms, nterms, num );
	
	for( i=0; i<2; i++ ) {
		xv[0 + i] = x[0 + i];
		xv[2 + i] = x[2 + i];
		xv[4 + i] = x[10 + i];
		}
	i_tomsu( xv, u1 );
	
	for( i=0; i<2; i++ ) {
		/*	xv[0 + i] = x[0 + i];	*/
		xv[2 + i] = x[4 + i];
		xv[4 + i] = x[8 + i];
		}
	i_tomsu( xv, u2 );
	
	I_MULT( u1, u2, temp );
	I_SQRT( temp, den );
	i_div( num, den, temp );
	i_acos( temp, out );
}


void i_dih_best( double x[12], double dih_part[12], 
	double out[2] )
{
	int i, j;
	double xp[12], data[2];
	
	/* First find dih_max */
	for( i=0; i<12; i++ )
		xp[i] = x[i];		/* copy x, then modify copy */
	j = 0;
	for( i=0; i<6; i++ ) {
		data[0] = dih_part[j];
		data[1] = dih_part[j+1];
		if( data[0] > 0 )	/* positive sign */
			xp[j] = xp[j+1];
		if( data[1] < 0 )	/* negative sign */
			xp[j+1] = xp[j];
		j += 2;
	}
	i_dih( xp, data );
	out[1] = data[1];
	
	/* Next find dih_min */
	for( i=0; i<12; i++ )
		xp[i] = x[i];		/* copy x, then modify copy */
	j = 0;
	for( i=0; i<6; i++ ) {
		data[0] = dih_part[j];
		data[1] = dih_part[j+1];
		if( data[0] > 0 )	/* positive sign */
			xp[j+1] = xp[j];
		if( data[1] < 0 )	/* negative sign */
			xp[j] = xp[j+1];
		j += 2;
	}
	i_dih( xp, data );
	out[0] = data[0];
}


void i_dbvol( double y[12], double out[2] )
{
	double tp[2], temp[2], frac[2];
	
	i_tvol( y, temp );
	ROUND_DOWN;
	tp[0] = 6.0*temp[0];
	ROUND_UP;
	tp[1] = 6.0*temp[1];
	
	i_afunc( y, temp );
	i_div( tp, temp, frac );
	i_atan( frac, temp );
	
	ROUND_DOWN;
	out[0] = temp[0]*TWO_3_LO;
	ROUND_UP;
	out[1] = temp[1]*TWO_3_HI;
}


void i_dbvol_alt( double y[12], double sqrtdelta[2],
	double out[2] )
{
	double max, temp, mina;
	
	temp = sqrtdelta[0];
	mina = 2.0*max_a( y );
	ROUND_DOWN;
	max = temp/mina;
	out[0] = TWO_3_LO*atan( max ) - ATANERR;

	temp = sqrtdelta[1];
	mina = 2.0*min_a( y );
	ROUND_UP;
	max = temp/mina;
	out[1] = TWO_3_HI*atan( max ) + ATANERR;
}


void i_bigdelta( double x[12], double out[2] )
{
	double pterms[2], nterms[2], p1, p2, p3, p4;

	ROUND_DOWN;
	p1 = x[0]*x[6]*(x[2] + x[4] + x[8] + x[10]);
	p2 = x[2]*x[8]*(x[0] + x[4] + x[6] + x[10]);
	p3 = x[4]*x[10]*(x[0] + x[2] + x[6] + x[8]);
	pterms[0] = p1 + p2 + p3;
	p1 = x[0]*x[6]*(x[0] + x[6]);
	p2 = x[2]*x[8]*(x[2] + x[8]);
	p3 = x[4]*x[10]*(x[4] + x[10]);
	p4 = x[4]*(x[2]*x[6] + x[0]*x[8]) + x[10]*(x[0]*x[2] + x[6]*x[8]);
	nterms[0] = p1 + p2 + p3 + p4;
	ROUND_UP;
	p1 = x[1]*x[7]*(x[3] + x[5] + x[9] + x[11]);
	p2 = x[3]*x[9]*(x[1] + x[5] + x[7] + x[11]);
	p3 = x[5]*x[11]*(x[1] + x[3] + x[7] + x[9]);
	pterms[1] = p1 + p2 + p3;
	p1 = x[1]*x[7]*(x[1] + x[7]);
	p2 = x[3]*x[9]*(x[3] + x[9]);
	p3 = x[5]*x[11]*(x[5] + x[11]);
	p4 = x[5]*(x[3]*x[7] + x[1]*x[9]) + x[11]*(x[1]*x[3] + x[7]*x[9]);
	nterms[1] = p1 + p2 + p3 + p4;
	I_SUB( pterms, nterms, out );
	/* Only care about simplices for which delta > 0 */
	if( out[0] < 0.0 )
		out[0] = 0.0;
	if( out[1] < 0.0 )
		out[1] = 0.0;
}


void i_bigdelta_best( double x[12], double out[2] )
{
	int i, j;
	double delta_part[12], xp[12], data[2];
	
	/* Generate best possible bounds for bigdelta */
	i_delta_partials(  x, delta_part );

	/* First find delta_max */
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
	out[1] = rough_max_delta( xp );
	
	/* Next find delta_min */
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
	out[0] = rough_min_delta( xp );
}


void i_bigdelta_partials_best( double x[12], double delta[2],
	double delta_part[12] )
{
	int i, j;
	double xp[12], data[2];
	
	/* Generate best possible bounds for bigdelta */
	i_delta_partials(  x, delta_part );

	/* First find delta_max */
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
	
	/* Next find delta_min */
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
}


int i_istet( double y[12] )
{
	double x[12], vol[2];
	int i;
	
	ROUND_DOWN;
	for( i=0; i<12; i+=2 )
		x[i] = y[i]*y[i];
	ROUND_UP;
	for( i=1; i<12; i+=2 )
		x[i] = y[i]*y[i];
	
	i_bigdelta( x, vol );
	/* If a subset is valid, then return true? */
	if( vol[1] > 0.0 )
		return( 1 );
	else
		return( 0 );
}


void i_tomsrho( double x[12], double out[2] )
{
	double pterms[2], nterms[2];
	double x2[12];
	int i;

	ROUND_DOWN;
	pterms[0] = 2.0*(x[0]*x[2]*x[6]*x[8] + x[0]*x[4]*x[6]*x[10] + x[2]*x[4]*x[8]*x[10]);
	for( i=0; i<12; i+=2 )
		x2[i] = x[i]*x[i];
	nterms[0] = x2[0]*x2[6] + x2[2]*x2[8] + x2[4]*x2[10];
	ROUND_UP;
	pterms[1] = 2.0*(x[1]*x[3]*x[7]*x[9] + x[1]*x[5]*x[7]*x[11] + x[3]*x[5]*x[9]*x[11]);
	for( i=1; i<12; i+=2 )
		x2[i] = x[i]*x[i];
	nterms[1] = x2[1]*x2[7] + x2[3]*x2[9] + x2[5]*x2[11];
	I_SUB( pterms, nterms, out );
	/* If delta > 0, rho > 0, only care about delta > 0 */
	if( out[0] < 0.0 )
		out[0] = 0.0;
	if( out[1] < 0.0 )
		out[1] = 0.0;
}


void i_crad( double y[12], double out[2] )
{
	double x[12], top[2], bot[2], temp[2];
	int i;
	
	ROUND_DOWN;
	for( i=0; i<12; i+=2 )
		x[i] = y[i]*y[i];
	ROUND_UP;
	for( i=1; i<12; i+=2 )
		x[i] = y[i]*y[i];
	i_tomsrho( x, top );
	i_bigdelta( x, bot );
	i_div( top, bot, temp );
	I_SQRT( temp, top );
	i_smult( 0.5, top, out );
}


void i_crad2( double x[12], double out[2] )
{
	double top[2], bot[2], temp[2];
	
	i_tomsrho( x, top );
	i_bigdelta_best( x, bot );
	i_div( top, bot, temp );
	out[0] = 0.25*temp[0];
	out[1] = 0.25*temp[1];
}


void i_crad2_best( double x[12], double out[2] )
{
	int i, j;
	double delta[2], delta_xpars[12], signs[12];
	double xp[12], data[2];
	
	i_bigdelta_partials_best( x, delta, delta_xpars );
	i_sign_crad2_xpars( x, delta, delta_xpars, signs );
	
	/*
	ROUND_NEAR;
	printf("delta = [%.18f, %.18f]\n", delta[0], delta[1]);
	printf("delta_i = \n");
	for( i=0; i<12; i+=2 ) {
		printf("[%.18f, %.18f]\n", 
			delta_xpars[i], delta_xpars[i+1]);
	}
	printf("i_sign = \n");
	for( i=0; i<12; i+=2 ) {
		printf("[%.18f, %.18f]\n", signs[i], signs[i+1]);
	}
	*/
	
	/* First find crad2_max */
	for( i=0; i<12; i++ )
		xp[i] = x[i];		/* copy x, then modify copy */
	j = 0;
	for( i=0; i<6; i++ ) {
		data[0] = signs[j];
		data[1] = signs[j+1];
		if( data[0] > 0 )	/* positive sign */
			xp[j] = xp[j+1];
		if( data[1] < 0 )	/* negative sign */
			xp[j+1] = xp[j];
		j += 2;
	}
	i_crad2( xp, data );
	out[1] = data[1];
	
	/* Next find crad2_min */
	for( i=0; i<12; i++ )
		xp[i] = x[i];		/* copy x, then modify copy */
	j = 0;
	for( i=0; i<6; i++ ) {
		data[0] = signs[j];
		data[1] = signs[j+1];
		if( data[0] > 0 )	/* positive sign */
			xp[j+1] = xp[j];
		if( data[1] < 0 )	/* negative sign */
			xp[j] = xp[j+1];
		j += 2;
	}
	i_crad2( xp, data );
	out[0] = data[0];
}


void i_tomsu( double x[6], double out[2] )
{
	double pterms[2], nterms[2];
	
	ROUND_DOWN;
	pterms[0] = 2.0*(x[0]*(x[4] + x[2]) + x[2]*x[4]);
	nterms[0] = x[0]*x[0] + x[2]*x[2] + x[4]*x[4];
	ROUND_UP;
	pterms[1] = 2.0*(x[1]*(x[5] + x[3]) + x[3]*x[5]);
	nterms[1] = x[1]*x[1] + x[3]*x[3] + x[5]*x[5];
	out[1] = pterms[1] - nterms[0];
	ROUND_DOWN;
	out[0] = pterms[0] - nterms[1];
	if( out[0] < 0.0 )
		out[0] = 0.0;
}


void i_tomsv( double x[12], double out[2] )
{
	double pterms[2], nterms[2];
	double t1, t2, t3;
	
	ROUND_DOWN;
	t1 = x[0]*x[6];
	t2 = x[2]*x[8];
	t3 = x[4]*x[10];
	pterms[0] = t1*(x[2] + x[10]) + t2*(x[0] + x[10]) + t3*(x[0] + x[2]);
	nterms[0] = t1*x[0] + t2*x[2] + t3*x[10] + 2.0*x[0]*x[2]*x[10];
	ROUND_UP;
	t1 = x[1]*x[7];
	t2 = x[3]*x[9];
	t3 = x[5]*x[11];
	pterms[1] = t1*(x[3] + x[11]) + t2*(x[1] + x[11]) + t3*(x[1] + x[3]);
	nterms[1] = t1*x[1] + t2*x[3] + t3*x[11] + 2.0*x[1]*x[3]*x[11];
	out[1] = pterms[1] - nterms[0];
	ROUND_DOWN;
	out[0] = pterms[0] - nterms[1];
}


void i_auxPfun( double x[6], double out[2] )
{
	double p1[2], p2[2], temp[2];
	double pterms[2], nterms[2];
	
	ROUND_DOWN;
	pterms[0] = 2.0*x[0]*x[2] + x[0]*x[4] + x[2]*x[4];
	nterms[0] = x[0]*x[0] + x[2]*x[2];
	ROUND_UP;
	pterms[1] = 2.0*x[1]*x[3] + x[1]*x[5] + x[3]*x[5];
	nterms[1] = x[1]*x[1] + x[3]*x[3];
	I_SUB( pterms, nterms, p1 );
	
	i_tomsu( x, temp );
	i_smult( 48.0, temp, p2 );
	i_div( p1, p2, out );
}


void i_tomsP( double x[12], double out[2] )
{
	double xv[6], p1[2], p2[2], p3[2], temp[2];
	int i;
	
	for( i=0; i<4; i++ )
		xv[i] = x[i];
	xv[4] = x[10];
	xv[5] = x[11];
	
	i_auxPfun( xv, p1 );
	i_tomsv( x, p2 );
	i_bigdelta( x, temp );
	I_SQRT( temp, p3 );
	I_MULT( p1, p2, temp );
	i_div( temp, p3, out );
}


/* This is tomsP without sqrt( bigdelta ) */
void i_tomsPwod( double x[12], double out[2] )
{
	double xv[6], p1[2], p2[2];
	int i;
	
	for( i=0; i<2; i++ ) {
		xv[0 + i] = x[0 + i];
		xv[2 + i] = x[2 + i];
		xv[4 + i] = x[10 + i];
		}
	i_auxPfun( xv, p1 );
	i_tomsv( x, p2 );
	I_MULT( p1, p2, out );
}

/*
tomschi[{x1_,x2_,x3_,x4_,x5_,x6_}]:=
	x1 x4 (x5 + x6) + x2 x5 (x4 + x6) + 
		x3 x6 (x4 + x5) - 2 x4 x5 x6 - 
		x1 x4^2 - x2 x5^2 - x3 x6^2
*/

void i_tomschi( double x[12], double out[2] )
{
	double pterms[2], nterms[2];
	
	ROUND_DOWN;
	pterms[0] = x[0]*x[6]*(x[8] + x[10]) +
		x[2]*x[8]*(x[6] + x[10]) + x[4]*x[10]*(x[6] + x[8]);
	nterms[0] = 2.0*x[6]*x[8]*x[10] + x[0]*x[6]*x[6] +
		x[2]*x[8]*x[8] + x[4]*x[10]*x[10];
	ROUND_UP;
	pterms[1] = x[1]*x[7]*(x[9] + x[11]) +
		x[3]*x[9]*(x[7] + x[11]) + x[5]*x[11]*(x[7] + x[9]);
	nterms[1] = 2.0*x[7]*x[9]*x[11] + x[1]*x[7]*x[7] +
		x[3]*x[9]*x[9] + x[5]*x[11]*x[11];
	out[1] = pterms[1] - nterms[0];
	ROUND_DOWN;
	out[0] = pterms[0] - nterms[1];
}


/* This version makes no assumptions about whether
the faces of the tetrahedron are acute or not.  */

double max_tomschi( double x[12] )
{
	double pterms, nterms;
	
	ROUND_DOWN;
	nterms = 2.0*x[6]*x[8]*x[10] + x[0]*x[6]*x[6] +
		x[2]*x[8]*x[8] + x[4]*x[10]*x[10];
	ROUND_UP;
	pterms = x[1]*x[7]*(x[9] + x[11]) +
		x[3]*x[9]*(x[7] + x[11]) + x[5]*x[11]*(x[7] + x[9]);
	return( pterms - nterms );
}


/* s_max_tomschi() depends on the assumption that all faces
are acute.  */

double s_max_tomschi( double x[12] )
{
	double pterms, nterms;
	
	ROUND_DOWN;
	nterms = 2.0*x[6]*x[8]*x[10] + x[1]*x[6]*x[6] +
		x[3]*x[8]*x[8] + x[5]*x[10]*x[10];
	ROUND_UP;
	pterms = x[1]*x[7]*(x[9] + x[11]) +
		x[3]*x[9]*(x[7] + x[11]) + x[5]*x[11]*(x[7] + x[9]);
	return( pterms - nterms );
}

/* s_min_tomschi() depends on the assumption that all faces
are acute.  */

double s_min_tomschi( double x[12] )
{
	double pterms, nterms, min, temp;
	int i, j, k, ind[3];
	
	min = 1.0e10;
	
	for( i=0; i<2; i++ ) {
		for( j=0; j<2; j++ ) {
			for( k=0; k<2; k++ ) {
				ROUND_UP;
				ind[0] = 6 + i;
				ind[1] = 8 + j;
				ind[2] = 10 + k;
				nterms = 2.0*x[ ind[0] ]*x[ ind[1] ]*x[ ind[2] ] + 
					x[0]*x[ ind[0] ]*x[ ind[0] ] +
					x[2]*x[ ind[1] ]*x[ ind[1] ] + 
					x[4]*x[ ind[2] ]*x[ ind[2] ];
				ROUND_DOWN;
				pterms = x[0]*x[ ind[0] ]*(x[ ind[1] ] + x[ ind[2] ]) +
					x[2]*x[ ind[1] ]*(x[ ind[0] ] + x[ ind[2] ]) + 
					x[4]*x[ ind[2] ]*(x[ ind[0] ] + x[ ind[1] ]);
				temp =  pterms - nterms;
				if( temp < min )
					min = temp;
			}
		}
	}
	return( min );
}


/* s_crad2_6() depends on the assumption that all faces
are acute, and that the long edge is x6.  
s_crad2_6() returns the square of the circumradius of the
simplex in out[], and returns in addition the
orientation of the face.  */

int s_crad2_6( double x[12], double out[2] )
{
	double xp[12], xpp[12], top[2], bot[2];
	int c456, c126, i;
	
	for( i=0; i<2; i++ ) {
		xp[i] = x[6 + i];
		xp[2 + i] = x[8 + i];
		xp[4 + i] = x[4 + i];
		xp[6 + i] = x[i];
		xp[8 + i] = x[2 + i];
		xp[10 + i] = x[10 + i];
	}
	
	if( s_min_tomschi( x ) > 0.0 )
		c456 = 1;
	else if( s_max_tomschi( x ) < 0.0 )
		c456 = -1;
	else
		c456 = 0; /* no information */
	
	if( s_min_tomschi( xp ) > 0.0 )
		c126 = 1;
	else if( s_max_tomschi( xp ) < 0.0 )
		c126 = -1;
	else
		c126 = 0; /* no information */
	
	/* min_crad comes from xp, max_crad comes from xpp */
	
	for( i=0; i<12; i++ ) {
		xp[i] = x[i];
		xpp[i] = x[i];
	}
	
	xp[11] = xp[10];
	xpp[10] = xpp[11];
	
	if( c456 > 0 ) {
		xp[1] = xp[0];
		xp[3] = xp[2];
		xpp[0] = xpp[1];
		xpp[2] = xpp[3];
	}
	else if( c456 < 0 ) {
		xp[0] = xp[1];
		xp[2] = xp[3];
		xpp[1] = xpp[0];
		xpp[3] = xpp[2];
	}
	if( c126 > 0 ) {
		xp[7] = xp[6];
		xp[9] = xp[8];
		xpp[6] = xpp[7];
		xpp[8] = xpp[9];
	}
	else if( c126 < 0 ) {
		xp[6] = xp[7];
		xp[8] = xp[9];
		xpp[7] = xpp[6];
		xpp[9] = xpp[8];
	}
	
	i = c456*c126;
	if( i > 0 ) {
		xp[5] = xp[4];
		xpp[4] = xpp[5];
	}
	else if( i < 0 ) {
		xp[4] = xp[5];
		xpp[5] = xpp[4];
	}
	
	i_tomsrho( xp, top );
	i_bigdelta( xp, bot );
	ROUND_DOWN;
	out[0] = 0.25*top[0]/bot[1];

	i_tomsrho( xpp, top );
	i_bigdelta( xpp, bot );
	ROUND_UP;
	out[1] = 0.25*top[1]/bot[0];
	return( c126 );
}


/* s_crad2_1() depends on the assumption that all faces
are acute, and that the long edge is x1.  
s_crad2_1() returns the square of the circumradius of the
simplex.  */

void s_crad2_1( double x[12], double out[2], int faceinfo[2] )
{
	double xp[12], xpp[12], top[2], bot[2];
	int c135, c126, i;
	
	for( i=0; i<2; i++ ) {
		xp[i] = x[6 + i];
		xp[2 + i] = x[8 + i];
		xp[4 + i] = x[4 + i];
		xp[6 + i] = x[i];
		xp[8 + i] = x[2 + i];
		xp[10 + i] = x[10 + i];
	}
	
	if( s_min_tomschi( xp ) > 0.0 )
		c126 = 1;
	else if( s_max_tomschi( xp ) < 0.0 )
		c126 = -1;
	else
		c126 = 0; /* no information */

	for( i=0; i<2; i++ ) {
		/* xp[i] = x[6 + i]; */
		xp[2 + i] = x[10 + i];
		xp[4 + i] = x[2 + i];
		/* xp[6 + i] = x[i]; */
		xp[8 + i] = x[4 + i];
		xp[10 + i] = x[8 + i];
	}

	if( s_min_tomschi( xp ) > 0.0 )
		c135 = 1;
	else if( s_max_tomschi( xp ) < 0.0 )
		c135 = -1;
	else
		c135 = 0; /* no information */
	
	faceinfo[0] = c126;
	faceinfo[1] = c135;
	
	/* min_crad comes from xp, max_crad comes from xpp */
	
	for( i=0; i<12; i++ ) {
		xp[i] = x[i];
		xpp[i] = x[i];
	}
	
	xp[1] = xp[0];
	xpp[0] = xpp[1];
	
	if( c135 > 0 ) {	/* determines 2 (2) and 6 (10) */
		xp[11] = xp[10];
		xp[3] = xp[2];
		xpp[10] = xpp[11];
		xpp[2] = xpp[3];
	}
	else if( c135 < 0 ) {
		xp[10] = xp[11];
		xp[2] = xp[3];
		xpp[11] = xpp[10];
		xpp[3] = xpp[2];
	}
	if( c126 > 0 ) {	/* determines 3 (4) and 5 (8) */
		xp[5] = xp[4];
		xp[9] = xp[8];
		xpp[4] = xpp[5];
		xpp[8] = xpp[9];
	}
	else if( c126 < 0 ) {
		xp[4] = xp[5];
		xp[8] = xp[9];
		xpp[5] = xpp[4];
		xpp[9] = xpp[8];
	}
	
	i = c135*c126;
	if( i > 0 ) {			/* determines 4 (6) */
		xp[7] = xp[6];
		xpp[6] = xpp[7];
	}
	else if( i < 0 ) {
		xp[6] = xp[7];
		xpp[7] = xpp[6];
	}
	
	i_tomsrho( xp, top );
	i_bigdelta( xp, bot );
	ROUND_DOWN;
	out[0] = 0.25*top[0]/bot[1];

	i_tomsrho( xpp, top );
	i_bigdelta( xpp, bot );
	ROUND_UP;
	out[1] = 0.25*top[1]/bot[0];
}


void i_rogersvol2( double xyz[6], double out[2] )
{
	double x, y, z, temp, min, max, xp[2], yp[2];
	double xpar, ypar, yzero[2], t1, t2;
	int gotone, gotx, goty;
	
	/* Check for NaNs in s_density_vor() */

	/* First compute partials, see if we can reduce the
	workload.  */
	
	gotx = 0;
	goty = 0;
	ROUND_DOWN;
	xpar = xyz[2] - 2.0*xyz[1];
	ypar = xyz[0] + xyz[4] - 2.0*xyz[3];
	ROUND_UP;
	if( xpar > 0.0 ) {
		xp[0] = xyz[0];
		xp[1] = xyz[1];
		gotx = 1;
	} else {
		xpar = xyz[3] - 2.0*xyz[0];
		if( xpar < 0.0 ) {
			xp[0] = xyz[1];
			xp[1] = xyz[0];
			gotx = 1;
		}
	}
	if( ypar > 0.0 ) {
		yp[0] = xyz[2];
		yp[1] = xyz[3];
		goty = 1;
	} else {
		ypar = xyz[1] + xyz[5] - 2.0*xyz[2];
		if( ypar < 0.0 ) {
			yp[0] = xyz[3];
			yp[1] = xyz[2];
			goty = 1;
		}
	}
	if( gotx && goty ) {
		t1 = yp[1] - xp[1];
		t2 = xyz[5] - yp[1];
		max = xp[1]*t1*t2;
		if( max < 0.0 )
			max = 0.0;
		out[1] = max*ONE_36_HI;
		ROUND_DOWN;
		t1 = yp[0] - xp[0];
		if( t1 < 0.0 )
			t1 = 0.0;
		t2 = xyz[4] - yp[0];
		if( t2 < 0.0 )
			t2 = 0.0;
		min = xp[0]*t1*t2;
		out[0] = min*ONE_36_LO;
		return;
	} else if( gotx ) {	/* x_partial bounded away from zero */
		/* Compute the minimum:  */
		
		ROUND_DOWN;
		z = xyz[4];
		
		x = xp[0];
		y = xyz[2];
		t1 = y - x;
		if( t1 < 0.0 )
			t1 = 0.0;
		t2 = z - y;
		if( t2 < 0.0 )
			t2 = 0.0;
		min = x*t1*t2;
		
		y = xyz[3];
		t1 = y - x;
		if( t1 < 0.0 )
			t1 = 0.0;
		t2 = z - y;
		if( t2 < 0.0 )
			t2 = 0.0;
		temp = x*t1*t2;
		if( temp < min )
			min = temp;
		
		out[0] = min*ONE_36_LO;
		
		/* Compute the maximum:  */
		z = xyz[5];
			/* Check for boundary maximum.  */
		gotone = 0;
		
		x = xp[1];
		yzero[1] = 0.5*(x + z);
		ROUND_DOWN;
		yzero[0] = 0.5*(x + z);
		ROUND_UP;
		if( xyz[2] <= yzero[1] && yzero[0] <= xyz[3] ) {
			max = x*(yzero[1] - x)*(z - yzero[0]);
			gotone = 1;
		}
			/* Check vertices (if none of the boundary stuff came
		through).  */
		
		if( !gotone ) {
			x = xp[1];
			y = xyz[2];
			max = x*(y - x)*(z - y);
			
			y = xyz[3];
			temp = x*(y - x)*(z - y);
			if( temp > max )
				max = temp;
		}
		out[1] = max*ONE_36_HI;
		return;
	} else if( goty ) {	/* y_partial bounded away from zero */
		/* Compute the minimum:  */
		
		ROUND_DOWN;
		z = xyz[4];
		
		x = xyz[0];
		y = yp[0];
		t1 = y - x;
		if( t1 < 0.0 )
			t1 = 0.0;
		t2 = z - y;
		if( t2 < 0.0 )
			t2 = 0.0;
		min = x*t1*t2;
				
		x = xyz[1];
		t1 = y - x;
		if( t1 < 0.0 )
			t1 = 0.0;
		temp = x*t1*t2;
		if( temp < min )
			min = temp;
		
		out[0] = min*ONE_36_LO;
		/* Compute the maximum:  */
		z = xyz[5];
			/* Check for boundary maximum.  */
		gotone = 0;
		
		y = yp[1];
		x = 0.5*y;
		if( xyz[0] <= x && x <= xyz[1] ) {
			max = x*(y - x)*(z - y);
			gotone = 1;
		}
				
			/* Check vertices (if none of the boundary stuff came
		through.  */
		
		if( !gotone ) {
			x = xyz[0];
			y = yp[1];
			max = x*(y - x)*(z - y);
						
			x = xyz[1];
			temp = x*(y - x)*(z - y);
			if( temp > max )
				max = temp;
		}
		if( max < 0.0 )
			max = 0.0;
		out[1] = max*ONE_36_HI;
		return;
	} else {
	/* If we don't have any useful bounds on the partials, do
	it the hard way:  */
	
	i_raw_rogersvol2( xyz, out );
	}
}

/* i_raw_rogersvol2 computes bounds on the rogers volume,
assuming that no partial was bounded away from zero.  */

void i_raw_rogersvol2( double xyz[6], double out[2] )
{
	double x, y, z, temp, min, max, xp[2], yp[2], t1, t2;
	int gotone;
	
	/* Check for NaNs in s_density_vor() */
	
	/* First compute the minimum-- use the fact that we are
	quadratic with negative leading coefficient in each
	variable.  */
	
	ROUND_DOWN;
	z = xyz[4];
	
	x = xyz[0];
	y = xyz[2];
	t1 = y - x;
	if( t1 < 0.0 )
		t1 = 0.0;
	t2 = z - y;
	if( t2 < 0.0 )
		t2 = 0.0;
	min = x*t1*t2;
	
	y = xyz[3];
	t1 = y - x;
	if( t1 < 0.0 )
		t1 = 0.0;
	t2 = z - y;
	if( t2 < 0.0 )
		t2 = 0.0;
	temp = x*t1*t2;
	if( temp < min )
		min = temp;
	
	x = xyz[1];
	t1 = y - x;
	if( t1 < 0.0 )
		t1 = 0.0;
	temp = x*t1*t2;
	if( temp < min )
		min = temp;
	
	y = xyz[2];
	t1 = y - x;
	if( t1 < 0.0 )
		t1 = 0.0;
	t2 = z - y;
	if( t2 < 0.0 )
		t2 = 0.0;
	temp = x*t1*t2;
	if( temp < min )
		min = temp;
	
	out[0] = min*ONE_36_LO;

	/* Now compute the maximum.  */
	
	max = 0.0;
	z = xyz[5];
	xp[0] = z*ONE_3_LO;
	ROUND_UP;
	xp[1] = z*ONE_3_HI;
	yp[0] = 2.0*xp[0];
	yp[1] = 2.0*xp[1];
	/* Check for global maximum.  */
	if( xyz[0] <= xp[1] && xp[0] <= xyz[1] &&
			xyz[2] <= yp[1] && yp[0] <= xyz[3] )
		max = xp[1]*(yp[1] - xp[0])*(z - yp[0]);
	else {
		/* Check for boundary maximum.  */
		gotone = 0;
		
		y = xyz[2];
		x = 0.5*y;
		if( xyz[0] <= x && x <= xyz[1] ) {
			temp = x*(y - x)*(z - y);
			gotone = 1;
		}
		if( temp > max )
			max = temp;
		
		y = xyz[3];
		x = 0.5*y;
		if( xyz[0] <= x && x <= xyz[1] ) {
			temp = x*(y - x)*(z - y);
			gotone = 1;
		}
		if( temp > max )
			max = temp;

		x = xyz[1];
		yp[1] = 0.5*(x + z);
		ROUND_DOWN;
		yp[0] = 0.5*(x + z);
		ROUND_UP;
		if( xyz[2] <= yp[1] && yp[0] <= xyz[3] ) {
			temp = x*(yp[1] - x)*(z - yp[0]);
			gotone = 1;
		}
		if( temp > max )
			max = temp;
		
		x = xyz[0];
		yp[1] = 0.5*(x + z);
		ROUND_DOWN;
		yp[0] = 0.5*(x + z);
		ROUND_UP;
		if( xyz[2] <= yp[1] && yp[0] <= xyz[3] ) {
			temp = x*(yp[1] - x)*(z - yp[0]);
			gotone = 1;
		}
		if( temp > max )
			max = temp;
		
		/* Check vertices (if none of the boundary stuff came
		through.  */
		
		if( !gotone ) {
			x = xyz[0];
			y = xyz[2];
			temp = x*(y - x)*(z - y);
			if( temp > max )
				max = temp;
			
			y = xyz[3];
			temp = x*(y - x)*(z - y);
			if( temp > max )
				max = temp;
			
			x = xyz[1];
			temp = x*(y - x)*(z - y);
			if( temp > max )
				max = temp;
			
			y = xyz[2];
			temp = x*(y - x)*(z - y);
			if( temp > max )
				max = temp;
		}
	}
	if( max < 0.0 )
		max = 0.0;
	out[1] = max*ONE_36_HI;
}

/*  rogers_density( a, b, c ) = 
	4*atan( sqrt( (b^2 - a^2) (c^2 - b^2) )/( 
		ab + b^2 + ac + bc ) )/( a sqrt( (b^2 - a^2) (c^2 - b^2) ) )
		
	and is monotonic decreasing in 1 < a < b < c.  */

void i_rogers_density( double abc[6], double out[2] )
{
	double a, b, c, top, bot, den, x, t1, t2;
	/* int i; */
	
	/* Check for NaNs in s_density_vor() */
	
	a = abc[1];
	b = abc[3];
	c = abc[5];
	
	/* Use atan(x)/x, which is decreasing for x >= 0.  Call
	this function atanquot(x), for arctan quotient.  With this
	formulation, rogers_density = 4 atanquot(x)/den, where
	x = sqrt( ( (b-a)*(c-b) )/( (a+b)*(b+c) ) ) and
	den = a*(a+b)*(b+c)  */
	
	/* Make x as large as possible */
	ROUND_DOWN;
	bot = (a + b)*(b + c);
	ROUND_UP;
	t1 = b - a;
	t2 = c - b;
	if( t2 < 0.0 )
		t2 = 0.0;
	top = t1*t2;
	x = sqrt( top/bot );
	den = a*(a + b)*(b + c);
	ROUND_DOWN;
	top = min_atanquot( x );
	/* ROUND_DOWN; */
	out[0] = 4.0*top/den;

	a = abc[0];
	b = abc[2];
	c = abc[4];
	
	/* Adjust values of b and c, since some estimates may
	be bad . . . */
	if( b < a )
		b = a;
	/*
	if( c < b )
		c = b;
	*/
	
	/* Make x as small as possible */
	ROUND_UP;
	bot = (a + b)*(b + c);
	ROUND_DOWN;
	t1 = b - a;
	t2 = c - b;
	if( t2 < 0.0 )
		t2 = 0.0;
	top = t1*t2;
	x = sqrt( top/bot );
	den = a*(a + b)*(b + c);
	ROUND_UP;
	top = max_atanquot( x );
	/* ROUND_UP; */
	out[1] = 4.0*top/den;
}


void old_i_rogers_density( double abc[6], double out[2] )
{
	double a, b, c, apb, bpc, spart, aspart, den;
	/* int i; */
	
	/* Check for NaNs in s_density_vor() */
	
	a = abc[1];
	b = abc[3];
	c = abc[5];
		
	ROUND_DOWN;
	spart = (b - a)*(c - b);
	ROUND_UP;
	apb = a + b;
	bpc = b + c;
	aspart = a*sqrt( apb*(b - a)*bpc*(c - b) );
	den = apb*bpc;
	ROUND_DOWN;
	apb = sqrt( spart/den );
	out[0] = 4.0*(atan( apb ) - ATANERR)/aspart;

	a = abc[0];
	b = abc[2];
	c = abc[4];
	
	/* Note that this method ignores the case where a > b, which
	can happen with octahedra (this arose since I did this
	originally for flat quad clusters, where that can't
	happen).  */
		
	if( b < c ) {
		ROUND_UP;
		spart = (b - a)*(c - b);
		ROUND_DOWN;
		apb = a + b;
		bpc = b + c;
		aspart = a*sqrt( apb*(b - a)*bpc*(c - b) );
		den = apb*bpc;
		ROUND_UP;
		apb = sqrt( spart/den );
		out[1] = 4.0*(atan( apb ) + ATANERR)/aspart;
	} else {	/* This sometimes happens . . . */
		/* c = b; */
		/*	ROUND_DOWN;	*/
		den = a*b*(a + b);
		ROUND_UP;
		out[1] = 2.0/den;
	}
}


void i_voronoivol( double x[12], double sqrtdelta[2], 
	double out[2] )
{
	double xp[12], p1[2], p2[2], p3[2], temp[2];
	int i;
	
	i_tomsPwod( x, p1 );
	
	for( i=0; i<2; i++ ) {
		xp[0 + i] = x[2 + i];
		xp[2 + i] = x[4 + i];
		xp[4 + i] = x[0 + i];
		xp[6 + i] = x[8 + i];
		xp[8 + i] = x[10 + i];
		xp[10 + i] = x[6 + i];
		}
	
	i_tomsPwod( xp, p2 );

	for( i=0; i<2; i++ ) {
		xp[0 + i] = x[4 + i];
		xp[2 + i] = x[0 + i];
		xp[4 + i] = x[2 + i];
		xp[6 + i] = x[10 + i];
		xp[8 + i] = x[6 + i];
		xp[10 + i] = x[8 + i];
		}
	
	i_tomsPwod( xp, p3 );
	
	I_ADD( p1, p2, temp );
	I_ADD( temp, p3, p1 );
	i_div( p1, sqrtdelta, out );
}


/* This version of vor() uses Tom's new formulation in
terms of the density of a Rogers simplex.  Hopefully, it
will give good bounds.  */

/* At the moment, we require that each face be acute, and
that the long edge be edge 6.  */

void s_density_vor_6( double y[12], double x[12], double out[2] )
{
	double rho[2], rho2[2], etas[6], eta2s[6], abc[6], xyz[6];
	double facex[6], sum[2], rvol[2], den[2], temp[2];
	double temp2[2], psum[2], maxeta, maxeta2, max;
	int chi126, i;
	
	chi126 = s_crad2_6( x, rho2 );
	xyz[4] = rho2[0];
	xyz[5] = rho2[1];
	I_SQRT( rho2, rho );
	abc[4] = rho[0];
	abc[5] = rho[1];
	
	max = rho2[0] + rho2[1];  /* keep track of NaNs */
	
	/* We compute the bounds on the circumradii of each face
	together, in order to find an acceptable lower bound on
	the circumradius of the simplex. */
	
	/* (2 3 4) (-> (2 4 6)) */
	for( i=0; i<2; i++ ) {
		facex[0 + i] = x[2 + i];
		facex[2 + i] = x[4 + i];
		facex[4 + i] = x[6 + i];
	}
	s_crad3x2( facex, eta2s );
	I_SQRT( eta2s, etas );
	
	maxeta = etas[0];
	maxeta2 = eta2s[0];
	
	max += eta2s[0] + eta2s[1];	/* keep track of NaNs */
	
	/* (1 3 5) (-> (0 4 8)) */
	for( i=0; i<2; i++ ) {
		facex[0 + i] = x[0 + i];
		/*	facex[2 + i] = x[4 + i];	*/
		facex[4 + i] = x[8 + i];
	}
	s_crad3x2( facex, eta2s + 2 );
	i_sqrt( eta2s + 2, etas + 2 );
	if( maxeta < etas[2] )
		maxeta = etas[2];
	if( maxeta2 < eta2s[2] )
		maxeta2 = eta2s[2];
		
	max += eta2s[2] + eta2s[3];	/* keep track of NaNs */
	
	/* (1 2 6) (-> (0 2 10)) */
	for( i=0; i<2; i++ ) {
		/*	facex[0 + i] = x[0 + i];	*/
		facex[2 + i] = x[2 + i];
		facex[4 + i] = x[10 + i];
	}
	s_crad3x2( facex, eta2s + 4 );
	i_sqrt( eta2s + 4, etas + 4 );
	if( maxeta < etas[4] )
		maxeta = etas[4];
	if( maxeta2 < eta2s[4] )
		maxeta2 = eta2s[4];
	
	max += eta2s[4] + eta2s[5];	/* keep track of NaNs */
	if( ISNAN( max ) ) {
		out[0] = max;
		out[1] = max;
	}
	
	/* ensure we have legitimate values for rho[0], rho2[0] */
	if( abc[4] < maxeta )
		abc[4] = maxeta;
	if( xyz[4] < maxeta2 )
		xyz[4] = maxeta2;
	
	/* first face */
	xyz[2] = eta2s[0];
	xyz[3] = eta2s[1];
	abc[2] = etas[0];
	abc[3] = etas[1];

	abc[0] = 0.5*y[2];
	abc[1] = 0.5*y[3];
	xyz[0] = 0.25*x[2];
	xyz[1] = 0.25*x[3];
	
	i_rogersvol2( xyz, temp );
	I_SQRT( temp, rvol );
	
	i_rogers_density( abc, den );
	
	I_SUB( den, i_doct_const, temp );
	I_MULT( rvol, temp, sum );
	
	abc[0] = 0.5*y[4];
	abc[1] = 0.5*y[5];
	xyz[0] = 0.25*x[4];
	xyz[1] = 0.25*x[5];

	i_rogersvol2( xyz, temp );
	I_SQRT( temp, rvol );
	
	i_rogers_density( abc, den );
	
	I_SUB( den, i_doct_const, temp );
	I_MULT( rvol, temp, temp2 );
	I_ADD( temp2, sum, temp );
	sum[0] = temp[0];
	sum[1] = temp[1];

	/* second face */
	xyz[2] = eta2s[2];
	xyz[3] = eta2s[3];
	abc[2] = etas[2];
	abc[3] = etas[3];
	
	/*
	abc[0] = 0.5*y[4];
	abc[1] = 0.5*y[5];
	xyz[0] = 0.25*x[4];
	xyz[1] = 0.25*x[5];
	*/
	
	i_rogersvol2( xyz, temp );
	I_SQRT( temp, rvol );
	
	i_rogers_density( abc, den );
	
	I_SUB( den, i_doct_const, temp );
	I_MULT( rvol, temp, temp2 );
	I_ADD( temp2, sum, temp );
	sum[0] = temp[0];
	sum[1] = temp[1];
	
	abc[0] = 0.5*y[0];
	abc[1] = 0.5*y[1];
	xyz[0] = 0.25*x[0];
	xyz[1] = 0.25*x[1];

	i_rogersvol2( xyz, temp );
	I_SQRT( temp, rvol );
	
	i_rogers_density( abc, den );
	
	I_SUB( den, i_doct_const, temp );
	I_MULT( rvol, temp, temp2 );
	I_ADD( temp2, sum, temp );
	sum[0] = temp[0];
	sum[1] = temp[1];

	/* (1 2 6) is the only face (under the given constraints) 
	for which chi can be negative.  */
	
	/* third face */
	xyz[2] = eta2s[4];
	xyz[3] = eta2s[5];
	abc[2] = etas[4];
	abc[3] = etas[5];
	
	/*
	abc[0] = 0.5*y[0];
	abc[1] = 0.5*y[1];
	xyz[0] = 0.25*x[0];
	xyz[1] = 0.25*x[1];
	*/
	
	i_rogersvol2( xyz, temp );
	I_SQRT( temp, rvol );
	
	i_rogers_density( abc, den );
	
	I_SUB( den, i_doct_const, temp );
	I_MULT( rvol, temp, psum );
	
	abc[0] = 0.5*y[2];
	abc[1] = 0.5*y[3];
	xyz[0] = 0.25*x[2];
	xyz[1] = 0.25*x[3];

	i_rogersvol2( xyz, temp );
	I_SQRT( temp, rvol );
	
	i_rogers_density( abc, den );
	
	I_SUB( den, i_doct_const, temp );
	I_MULT( rvol, temp, temp2 );
	I_ADD( temp2, psum, temp );
	psum[0] = temp[0];
	psum[1] = temp[1];
	if( chi126 == 1 )	{	/* chi > 0 */
		I_ADD( sum, psum, temp );
	} else if( chi126 == -1 ) { /* chi < 0 */
		I_SUB( sum, psum, temp );
	} else {
		temp[0] = -1.0;
		temp[1] = 1.0;
		I_MULT( temp, psum, temp2 );
		I_ADD( sum, temp2, temp );
	}
	sum[0] = temp[0];
	sum[1] = temp[1];
	
	out[0] = 4.0*sum[0];
	out[1] = 4.0*sum[1];
}


/* This version of vor() uses Tom's new formulation in
terms of the density of a Rogers simplex.  Hopefully, it
will give good bounds.  */

/* At the moment, we require that each face be acute, and
that the long edge be edge 1.  */

void s_density_vor_1( double y[12], double x[12], double out[2] )
{
	double rho[2], rho2[2], etas[6], eta2s[6], abc[6], xyz[6];
	double facex[6], sum[2], rvol[2], den[2], temp[2];
	double temp2[2], psum[2], maxeta, maxeta2, max;
	int chi126, chi135, faceinfo[2], i;
	
	s_crad2_1( x, rho2, faceinfo );
	chi126 = faceinfo[0];
	chi135 = faceinfo[1];
	
	xyz[4] = rho2[0];
	xyz[5] = rho2[1];
	I_SQRT( rho2, rho );
	abc[4] = rho[0];
	abc[5] = rho[1];
	
	max = rho2[0] + rho2[1];  /* keep track of NaNs */
	
	/* We compute the bounds on the circumradii of each face
	together, in order to find an acceptable lower bound on
	the circumradius of the simplex. */
	
	/* (2 3 4) (-> (2 4 6)) */
	for( i=0; i<2; i++ ) {
		facex[0 + i] = x[2 + i];
		facex[2 + i] = x[4 + i];
		facex[4 + i] = x[6 + i];
	}
	s_crad3x2( facex, eta2s );
	I_SQRT( eta2s, etas );
	
	/* We want to find the maximum lower bound on the etas, 
	that is, find the largest hard lower bound on the
	circumradius of a face.  This is a lower bound on the
	circumradius of the simplex.  */
	maxeta = etas[0];
	maxeta2 = eta2s[0];
	
	max += eta2s[0] + eta2s[1];	/* keep track of NaNs */
	
	/* (1 3 5) (-> (0 4 8)) */
	for( i=0; i<2; i++ ) {
		facex[0 + i] = x[0 + i];
		/*	facex[2 + i] = x[4 + i];	*/
		facex[4 + i] = x[8 + i];
	}
	s_crad3x2( facex, eta2s + 2 );
	i_sqrt( eta2s + 2, etas + 2 );
	if( maxeta < etas[2] )
		maxeta = etas[2];
	if( maxeta2 < eta2s[2] )
		maxeta2 = eta2s[2];
		
	max += eta2s[2] + eta2s[3];	/* keep track of NaNs */
	
	/* (1 2 6) (-> (0 2 10)) */
	for( i=0; i<2; i++ ) {
		/*	facex[0 + i] = x[0 + i];	*/
		facex[2 + i] = x[2 + i];
		facex[4 + i] = x[10 + i];
	}
	s_crad3x2( facex, eta2s + 4 );
	i_sqrt( eta2s + 4, etas + 4 );
	if( maxeta < etas[4] )
		maxeta = etas[4];
	if( maxeta2 < eta2s[4] )
		maxeta2 = eta2s[4];
	
	max += eta2s[4] + eta2s[5];	/* keep track of NaNs */
	if( ISNAN( max ) ) {
		out[0] = max;
		out[1] = max;
	}
	
	/* ensure we have legitimate values for rho[0], rho2[0] */
	if( abc[4] < maxeta )
		abc[4] = maxeta;
	if( xyz[4] < maxeta2 )
		xyz[4] = maxeta2;
	
	/* first face */		/* 2 3 4  (2 4 6) */
	xyz[2] = eta2s[0];
	xyz[3] = eta2s[1];
	abc[2] = etas[0];
	abc[3] = etas[1];

	abc[0] = 0.5*y[2];
	abc[1] = 0.5*y[3];
	xyz[0] = 0.25*x[2];
	xyz[1] = 0.25*x[3];
	
	i_rogersvol2( xyz, temp );
	I_SQRT( temp, rvol );
	/* printf("i_rogersvol = [%g, %g]\n", rvol[0], rvol[1]); */
	
	i_rogers_density( abc, den );
	/* printf("i_rogers_density = [%g, %g]\n", den[0], den[1]); */
	
	I_SUB( den, i_doct_const, temp );
	I_MULT( rvol, temp, sum );
	/* printf("sum = [%g, %g]\n", sum[0], sum[1]); */
	
	abc[0] = 0.5*y[4];
	abc[1] = 0.5*y[5];
	xyz[0] = 0.25*x[4];
	xyz[1] = 0.25*x[5];

	i_rogersvol2( xyz, temp );
	I_SQRT( temp, rvol );
	
	i_rogers_density( abc, den );
	
	I_SUB( den, i_doct_const, temp );
	I_MULT( rvol, temp, temp2 );
	I_ADD( temp2, sum, temp );
	sum[0] = temp[0];
	sum[1] = temp[1];
	
	/* Note that chi(1 3 5) can be negative.  */

	/* second face */		/* 1 3 5  (0 4 8) */
	xyz[2] = eta2s[2];
	xyz[3] = eta2s[3];
	abc[2] = etas[2];
	abc[3] = etas[3];
	
	/*
	abc[0] = 0.5*y[4];
	abc[1] = 0.5*y[5];
	xyz[0] = 0.25*x[4];
	xyz[1] = 0.25*x[5];
	*/
	
	i_rogersvol2( xyz, temp );
	I_SQRT( temp, rvol );
	
	i_rogers_density( abc, den );
	
	I_SUB( den, i_doct_const, temp );
	I_MULT( rvol, temp, psum );
	
	abc[0] = 0.5*y[0];
	abc[1] = 0.5*y[1];
	xyz[0] = 0.25*x[0];
	xyz[1] = 0.25*x[1];

	i_rogersvol2( xyz, temp );
	I_SQRT( temp, rvol );
	
	i_rogers_density( abc, den );
	
	I_SUB( den, i_doct_const, temp );
	I_MULT( rvol, temp, temp2 );
	I_ADD( temp2, psum, temp );
	psum[0] = temp[0];
	psum[1] = temp[1];
	if( chi135 == 1 )	{	/* chi > 0 */
		I_ADD( sum, psum, temp );
	} else if( chi135 == -1 ) { /* chi < 0 */
		I_SUB( sum, psum, temp );
	} else {
		temp[0] = -1.0;
		temp[1] = 1.0;
		I_MULT( temp, psum, temp2 );
		I_ADD( sum, temp2, temp );
	}
	sum[0] = temp[0];
	sum[1] = temp[1];

	/* Note that chi(1 2 6) can be negative.  */
	
	/* third face */		/* 1 2 6  (0 2 10) */
	xyz[2] = eta2s[4];
	xyz[3] = eta2s[5];
	abc[2] = etas[4];
	abc[3] = etas[5];
	
	/*
	abc[0] = 0.5*y[0];
	abc[1] = 0.5*y[1];
	xyz[0] = 0.25*x[0];
	xyz[1] = 0.25*x[1];
	*/
	
	i_rogersvol2( xyz, temp );
	I_SQRT( temp, rvol );
	
	i_rogers_density( abc, den );
	
	I_SUB( den, i_doct_const, temp );
	I_MULT( rvol, temp, psum );
	
	abc[0] = 0.5*y[2];
	abc[1] = 0.5*y[3];
	xyz[0] = 0.25*x[2];
	xyz[1] = 0.25*x[3];

	i_rogersvol2( xyz, temp );
	I_SQRT( temp, rvol );
	
	i_rogers_density( abc, den );
	
	I_SUB( den, i_doct_const, temp );
	I_MULT( rvol, temp, temp2 );
	I_ADD( temp2, psum, temp );
	psum[0] = temp[0];
	psum[1] = temp[1];
	if( chi126 == 1 )	{	/* chi > 0 */
		I_ADD( sum, psum, temp );
	} else if( chi126 == -1 ) { /* chi < 0 */
		I_SUB( sum, psum, temp );
	} else {
		temp[0] = -1.0;
		temp[1] = 1.0;
		I_MULT( temp, psum, temp2 );
		I_ADD( sum, temp2, temp );
	}
	sum[0] = temp[0];
	sum[1] = temp[1];
	
	out[0] = 4.0*sum[0];
	out[1] = 4.0*sum[1];
}


/* This version of vor() uses the "analytic continuation"
formula, and appears to be quite useless.  Might be helpful
eventually, who knows?  */

void i_vor( double y[12], double x[12], double sqrtdelta[2], 
	double out[2] )
{
	double temp1[2], temp2[2], temp3[2];
	
	i_dbvol_alt( y, sqrtdelta, temp1 );
	i_voronoivol( x, sqrtdelta, temp2 );
	I_MULT( temp2, i_doct_const, temp3 );
	I_SUB( temp1, temp3, temp2 );
	i_smult( 4.0, temp2, out );
}


void i_vor_alt( double x[12], double sqrtdelta[2], 
	double sol[2], double out[2] )
{
	double temp1[2], temp2[2], temp3[2];
	
	ROUND_DOWN;
	temp1[0] = ONE_3_LO*sol[0];
	ROUND_UP;
	temp1[1] = ONE_3_HI*sol[1];
	
	i_voronoivol( x, sqrtdelta, temp2 );
	I_MULT( temp2, i_doct_const, temp3 );
	I_SUB( temp1, temp3, temp2 );
	out[0] = 4.0*temp2[0];
	out[1] = 4.0*temp2[1];
}


/* This version of vor() combines both the rogers_density
version and the analytic continuation version, in an effort
to produce the best possible bounds (the rogers_density 
version seems to have trouble when the circumradius punches
through one of the faces, and that is making it hard to
verify some of the bounds).  Of course, doing both
computations is something of a pain, but it may pay off
in the end.  */

void best_i_vor_1( double y[12], double x[12], 
	double sqrtdelta[2], double out[2] )
{
	double ivor[2];
	
	/* Generally, we expect s_density_vor to give much better
	bounds. */
	
	i_vor( y, x, sqrtdelta, ivor );
	s_density_vor_1( y, x, out );
	
	/* Now compute intersection of intervals */
	
	if( ivor[0] > out[0] )
		out[0] = ivor[0];
	else if( ISNAN( out[0] ) )
		out[0] = ivor[0];
	/* (if one is a NaN, use the other one, which might not be) */
	if( ivor[1] < out[1] )
		out[1] = ivor[1];
	else if( ISNAN( out[1] ) )
		out[1] = ivor[1];
	
}


void best_i_vor_6( double y[12], double x[12], 
	double sqrtdelta[2], double out[2] )
{
	double ivor[2];
	
	/* Generally, we expect s_density_vor to give much better
	bounds. */
	
	i_vor( y, x, sqrtdelta, ivor );
	s_density_vor_6( y, x, out );
	
	/* Now compute intersection of intervals */
	
	if( ivor[0] > out[0] )
		out[0] = ivor[0];
	else if( ISNAN( out[0] ) )
		out[0] = ivor[0];
	/* (if one is a NaN, use the other one, which might not be) */
	if( ivor[1] < out[1] )
		out[1] = ivor[1];
	else if( ISNAN( out[1] ) )
		out[1] = ivor[1];
	
}


void i_crad3len( double y[6], double out[2] )
{
	double x[6], xc[2], yc[2], temp[2], temp2[2], ac[2];
	double bc[2];
	int i;
	
	ROUND_DOWN;
	for( i=0; i<6; i+=2 )
		x[i] = y[i]*y[i];
	ROUND_UP;
	for( i=1; i<6; i+=2 )
		x[i] = y[i]*y[i];
	ROUND_DOWN;
	temp[0] = x[0] + x[2] - x[5];
	ROUND_UP;
	temp[1] = x[1] + x[3] - x[4];
	temp2[0] = 2.0*y[0];
	temp2[1] = 2.0*y[1];
	i_div( temp, temp2, xc );
	I_MULT( xc, xc, temp );
	I_SUB( x + 2, temp, temp2 );
	I_SQRT( temp2, yc );
	ac[0] = 0.25*x[0];
	ac[1] = 0.25*x[1];
	I_MULT( y, xc, temp );
	I_SUB( x + 2, temp, temp2 );
	i_div( temp2, yc, bc );
	bc[0] *= 0.5;
	bc[1] *= 0.5;
	I_MULT( bc, bc, temp );
	I_ADD( ac, temp, temp2 );
	I_SQRT( temp2, out );
}


/* crad3len = Sqrt[x1 x2 x3/tomsu[x1, x2, x3]] */
void i_crad3x2( double x[6], double out[2] )
{
	double pterms[2], nterms[2], top[2], bot[2];

	ROUND_DOWN;
	top[0] = x[0]*x[2]*x[4];
	nterms[0] = x[0]*x[0] + x[2]*x[2] + x[4]*x[4];
	pterms[0] = 2.0*(x[0]*x[2] + x[0]*x[4] + x[2]*x[4]);
	ROUND_UP;
	top[1] = x[1]*x[3]*x[5];
	nterms[1] = x[1]*x[1] + x[3]*x[3] + x[5]*x[5];
	pterms[1] = 2.0*(x[1]*x[3] + x[1]*x[5] + x[3]*x[5]);
	bot[1] = pterms[1] - nterms[0];
	ROUND_DOWN;
	bot[0] = pterms[0] - nterms[1];
	out[0] = top[0]/bot[1];
	ROUND_UP;
	out[1] = top[1]/bot[0];
}


int i_isqrtet( double y[12] )
{
	int ok, i;
	
	ok = 1;
	i = 0;
	/* if any subset is qr, call the whole thing qr?  */
	while( (ok == 1 && i < 6) ) {
		if( !( y[i+1] >= 2.0 && y[i] <= 2.51 ) )
			ok = 0;
		i += 2;
		}
	return( ok );
}


/* Here is my philosophy for a complex scoring system:  the
problem at hand is characterizing cells which contain 
tetrahedra which may potentially be scored in different ways.
If a test on a cell comes up true (that is, true for all
tetrahedra within the cell), it should return "1".  If a
test comes up false (that is, false for all tetrahedra
within the cell), then it should return "0".  If the answer
is ambiguous/bivalent/indeterminate, that is, if there
are tetrahedra for which the test comes up true, and 
tetrahedra for which the test comes up false, or if the
characterization is indeterminate, the test should return
"-1".  */

/* scoring_system() answers the question of whether a given 
flat quarter should be gma (compression) scored or not.  
We assume that the long edge is 6.  */

int scoring_system( double x[12] )
{
	int i;
	double face[6], cf1[2], cf2[2];
		
	/* Face (1 2 6) -> (0 2 10) */
	face[0] = x[0];
	face[1] = x[1];
	face[2] = x[2];
	face[3] = x[3];
	face[4] = x[10];
	face[5] = x[11];
	
	s_crad3x2( face, cf1 );

	/* Face (4 5 6) -> (6 8 10) */
	face[0] = x[6];
	face[1] = x[7];
	face[2] = x[8];
	face[3] = x[9];
	
	s_crad3x2( face, cf2 );
	
	/* Check for NaNs */
	i = ISNAN( cf1[0] + cf1[1] + cf2[0] + cf2[1] );
	
	if( cf1[1] <= 2.0 && cf2[1] <= 2.0 )
		return( 1 );
	else if( (cf1[0] <= 2.0 && cf2[0] <= 2.0) || i )
		return( -1 );
	else
		return( 0 );
}


int old_scoring_system( double y[12], double x[12] )
{
	int i, j, k, opp, opp_low, big, pt_a, pt_b, pt_c;
	int pt_i, pt_ii, is_small;
	int adj_edges[12][2] = {{2, 10}, {4, 8},
												{0, 10}, {4, 6},
												{0, 8}, {2, 6},
												{2, 4}, {8, 10},
												{0, 4}, {6, 10},
												{0, 2}, {6, 8}};
	double face[6], cf1[2], cf2[2], cf3[2], cf4[2];
	
	j = -1;
	i = 0;
	while( j < 0 && i < 12 ) {
		if( y[i] >= 2.51 )
			j = i;
		i += 2;
		}
	if( j == -1 ) {
		printf("There appears to be no long edge.");
		return( -1 );
		}
	if( j < 6 )
		opp = j + 7;
	else
		opp = j - 5;
	/* opp is the index of the upper end of the edge opposite 
		the long edge. */
	opp_low = opp - 1;
	big = j;  /* big is the index of the (lower end of the) 
		long edge.  */
	
	/* Here's case a:  */
	if( y[opp] < 2.06 )
		pt_a = 1;
	else if( y[opp_low] >= 2.06 ) {
		return( 0 );
		/* pt_a = 0; */
		}
	else
		pt_a = -1;
	
	/* Here's case b:  */
	pt_b = -10;
	i = 0;
	while( pt_b < 0 && i < 12 ) {
		if( y[adj_edges[i][0] + 1] < 2.08 && 
			y[adj_edges[i][1] + 1] < 2.08 )
			pt_b = 1;
		else if( y[adj_edges[i][0]] < 2.08 && 
			y[adj_edges[i][1]] < 2.08 )
			pt_b = -1;
		i++;
		}
	if( pt_b == -10 )	{ /* fell through all cases */
		return( 0 );
		/* pt_b = 0; */
		}
	
	/* Here's case c:  */
		/* This case has two parts, i) and ii) */
	pt_i = 1;
	j = 1;
	for( i=0; i<12; i+=2 ) {
		if( i != big && i != opp_low && pt_i ) {
			if( y[i] >= 2.12 )
				pt_i = 0;	/* a whole edge is over */
			else if( y[i+1] >= 2.12 )
				j = 0;		/* part of an edge is over */
		}
	}
	if( (pt_i > 0) && (j < 1) )
		pt_i = -1;  /* some satisfy part i */
	
	j = big + 1;
	face[0] = x[big];
	face[1] = x[j];
	k = adj_edges[big][0];
	face[2] = x[k];
	face[3] = x[k + 1];
	k = adj_edges[big][1];
	face[4] = x[k];
	face[5] = x[k + 1];
	
	s_crad3x2( face, cf1 );

	/*
	face[0] = x[big];
	face[1] = x[j];
	*/
	k = adj_edges[j][0];
	face[2] = x[k];
	face[3] = x[k + 1];
	k = adj_edges[j][1];
	face[4] = x[k];
	face[5] = x[k + 1];
	
	s_crad3x2( face, cf2 );
	
	/* Check for NaNs */
	i = ISNAN( cf1[0] + cf1[1] + cf2[0] + cf2[1] );
	
	/* 1.39^2 = 1.9321 */
	if( cf1[1] < 1.9321 && cf2[1] < 1.9321 )
		pt_ii = 1;
	else if( (cf1[0] < 1.9321 && cf2[0] < 1.9321) || i )
		pt_ii = -1;
	else
		pt_ii = 0;
	
	if( pt_i > 0 || pt_ii > 0 )
		pt_c = 1;
	else if( pt_i < 0 || pt_ii < 0 )
		pt_c = -1;
	else {
		return( 0 );
		/* pt_c = 0; */
		}
	
	face[0] = x[opp_low];
	face[1] = x[opp];
	k = adj_edges[opp_low][0];
	face[2] = x[k];
	face[3] = x[k + 1];
	k = adj_edges[opp_low][1];
	face[4] = x[k];
	face[5] = x[k + 1];
	
	s_crad3x2( face, cf3 );

	/*
	face[0] = x[opp_low];
	face[1] = x[opp];
	*/
	k = adj_edges[opp][0];
	face[2] = x[k];
	face[3] = x[k + 1];
	k = adj_edges[opp][1];
	face[4] = x[k];
	face[5] = x[k + 1];
	
	s_crad3x2( face, cf4 );
	
	/* Check for NaNs */
	j = ISNAN( cf3[0] + cf3[1] + cf4[0] + cf4[1] );
	k = ( i || j );

	/* Check to see if all faces are small */
	if( cf1[1] < 2.0 && cf2[1] < 2.0 && cf3[1] < 2.0 && 
		cf4[1] < 2.0 ) {
		is_small = 1;
	}
	else if( k || (cf1[0] < 2.0 && cf2[0] < 2.0 && 
		cf3[0] < 2.0 && cf4[0] < 2.0) ) {
		is_small = -1;
	} else {
		/*
		is_small = 0;
		*/
		return( 0 );
	}

	/* Now combine all of the results */
	if( pt_a > 0 && pt_b > 0 && pt_c > 0 && is_small )
		return( 1 );
	else
		return( -1 );
}


/* octa_scoring() answers the question of whether a
given upright quarter should be gma (compression) 
scored or not.  We assume that the long edge is 1.  */

int octa_scoring( double x[12] )
{
	int i;
	double face[6], cf1[2], cf2[2];
		
	/* Face (1 2 6) -> (0 2 10) */
	face[0] = x[0];
	face[1] = x[1];
	face[2] = x[2];
	face[3] = x[3];
	face[4] = x[10];
	face[5] = x[11];
	
	s_crad3x2( face, cf1 );

	/* Face (1 3 5) -> (0 4 8) */
	face[2] = x[4];
	face[3] = x[5];
	face[4] = x[8];
	face[5] = x[9];
	
	s_crad3x2( face, cf2 );
	
	/* Check for NaNs */
	i = ISNAN( cf1[0] + cf1[1] + cf2[0] + cf2[1] );
	
	if( cf1[1] <= 2.0 && cf2[1] <= 2.0 )
		return( 1 );
	else if( (cf1[0] <= 2.0 && cf2[0] <= 2.0) || i )
		return( -1 );
	else
		return( 0 );
}


int old_octa_scoring( double y[12], double x[12] )
{
	int i, j, k, opp, opp_low, big, pt_c;
	int pt_i, pt_ii, is_small;
	int adj_edges[12][2] = {{2, 10}, {4, 8},
												{0, 10}, {4, 6},
												{0, 8}, {2, 6},
												{2, 4}, {8, 10},
												{0, 4}, {6, 10},
												{0, 2}, {6, 8}};
	double face[6], cf1[2], cf2[2], cf3[2], cf4[2];
	
	j = -1;
	i = 0;
	while( j < 0 && i < 12 ) {
		if( y[i] >= 2.51 )
			j = i;
		i += 2;
		}
	if( j == -1 ) {
		printf("There appears to be no long edge.");
		return( -1 );
		}
	if( j < 6 )
		opp = j + 7;
	else
		opp = j - 5;
	/* opp is the index of the upper end of the edge opposite 
		the long edge. */
	opp_low = opp - 1;
	big = j;  /* big is the index of the (lower end of the) 
		long edge.  */
	
	/* Cases (a) and (b) do not enter in the octahedral scoring
	system.  */
		
	/* Here's case c:  */
		/* This case has two parts, i) and ii) */
	pt_i = 1;
	j = 1;
	for( i=0; i<12; i+=2 ) {
		if( i != big && i != opp_low && pt_i ) {
			if( y[i] >= 2.12 )
				pt_i = 0;	/* a whole edge is over */
			else if( y[i+1] >= 2.12 )
				j = 0;		/* part of an edge is over */
		}
	}
	if( (pt_i > 0) && (j < 1) )
		pt_i = -1;  /* some satisfy part i */
	
	j = big + 1;
	face[0] = x[big];
	face[1] = x[j];
	k = adj_edges[big][0];
	face[2] = x[k];
	face[3] = x[k + 1];
	k = adj_edges[big][1];
	face[4] = x[k];
	face[5] = x[k + 1];
	
	s_crad3x2( face, cf1 );

	/*
	face[0] = x[big];
	face[1] = x[j];
	*/
	k = adj_edges[j][0];
	face[2] = x[k];
	face[3] = x[k + 1];
	k = adj_edges[j][1];
	face[4] = x[k];
	face[5] = x[k + 1];
	
	s_crad3x2( face, cf2 );
	
	/* Check for NaNs */
	i = ISNAN( cf1[0] + cf1[1] + cf2[0] + cf2[1] );
	
	/* 1.39^2 = 1.9321 */
	if( cf1[1] < 1.9321 && cf2[1] < 1.9321 )
		pt_ii = 1;
	else if( (cf1[0] < 1.9321 && cf2[0] < 1.9321) || i )
		pt_ii = -1;
	else
		pt_ii = 0;
	
	if( pt_i > 0 || pt_ii > 0 )
		pt_c = 1;
	else if( pt_i < 0 || pt_ii < 0 )
		pt_c = -1;
	else {
		return( 0 );
		/* pt_c = 0; */
		}
	
	face[0] = x[opp_low];
	face[1] = x[opp];
	k = adj_edges[opp_low][0];
	face[2] = x[k];
	face[3] = x[k + 1];
	k = adj_edges[opp_low][1];
	face[4] = x[k];
	face[5] = x[k + 1];
	
	s_crad3x2( face, cf3 );

	/*
	face[0] = x[opp_low];
	face[1] = x[opp];
	*/
	k = adj_edges[opp][0];
	face[2] = x[k];
	face[3] = x[k + 1];
	k = adj_edges[opp][1];
	face[4] = x[k];
	face[5] = x[k + 1];
	
	s_crad3x2( face, cf4 );
	
	/* Check for NaNs */
	j = ISNAN( cf3[0] + cf3[1] + cf4[0] + cf4[1] );
	k = ( i || j );

	/* Check to see if all faces are small */
	if( cf1[1] < 2.0 && cf2[1] < 2.0 && cf3[1] < 2.0 && 
		cf4[1] < 2.0 ) {
		is_small = 1;
	}
	else if( k || (cf1[0] < 2.0 && cf2[0] < 2.0 && 
		cf3[0] < 2.0 && cf4[0] < 2.0) ) {
		is_small = -1;
	} else {
		/*
		is_small = 0;
		*/
		return( 0 );
	}

	/* Now combine all of the results */
	if( pt_c > 0 && is_small )
		return( 1 );
	else
		return( -1 );
}


void i_score( double y[12], double x[12],
	double sqrtdelta[2], double out[2] )
{
	int sys;
	double gsc[2], vsc[2], temp;
	
	sys = scoring_system(  x );
	if( sys > 0 ) {
		out[0] = s_min_gma( y, sqrtdelta );
		out[1] = s_max_gma( y, sqrtdelta );
		}
	else if( sys == 0 )
		i_vor( y, x, sqrtdelta, out );
	else {
		gsc[0] = s_min_gma( y, sqrtdelta );
		gsc[1] = s_max_gma( y, sqrtdelta );
		i_vor( y, x, sqrtdelta, vsc );
		temp = gsc[0] + gsc[1] + vsc[0] + vsc[1];
		if( ISNAN( temp ) ) {
			out[0] = temp;
			out[1] = temp;
			return;
		}
		if( gsc[1] > vsc[1] )
			out[1] = gsc[1];
		else
			out[1] = vsc[1];
		if( gsc[0] < vsc[0] )
			out[0] = gsc[0];
		else
			out[0] = vsc[0];
		}
}


void s_octa_vor( double y[12], double x[12], double out[2] )
{
	double xp[12], yp[12], sc1[2], sc2[2], sum[2];
	int i, j, k;
	
	for( i=0; i<2; i++ ) {
		j = i;
		k = i;
		xp[j] = x[k];
		yp[j] = y[k];
		j = 2 + i;
		k = 8 + i;
		xp[j] = x[k];
		yp[j] = y[k];
		j = 4 + i;
		k = 10 + i;
		xp[j] = x[k];
		yp[j] = y[k];
		j = 6 + i;
		k = j;
		xp[j] = x[k];
		yp[j] = y[k];
		j = 8 + i;
		k = 2 + i;
		xp[j] = x[k];
		yp[j] = y[k];
		j = 10 + i;
		k = 4 + i;
		xp[j] = x[k];
		yp[j] = y[k];
	}
	s_density_vor_1( y, x, sc1 );
	s_density_vor_1( yp, xp, sc2 );
	I_ADD( sc1, sc2, sum );
	out[0] = 0.5*sum[0];
	out[1] = 0.5*sum[1];
}


void best_octa_vor( double y[12], double x[12], 
	double sqrtdelta[2], double out[2] )
{
	double xp[12], yp[12], sc1[2], sc2[2], sum[2];
	int i, j, k;
	
	for( i=0; i<2; i++ ) {
		j = i;
		k = i;
		xp[j] = x[k];
		yp[j] = y[k];
		j = 2 + i;
		k = 8 + i;
		xp[j] = x[k];
		yp[j] = y[k];
		j = 4 + i;
		k = 10 + i;
		xp[j] = x[k];
		yp[j] = y[k];
		j = 6 + i;
		k = j;
		xp[j] = x[k];
		yp[j] = y[k];
		j = 8 + i;
		k = 2 + i;
		xp[j] = x[k];
		yp[j] = y[k];
		j = 10 + i;
		k = 4 + i;
		xp[j] = x[k];
		yp[j] = y[k];
	}
	best_i_vor_1( y, x, sqrtdelta, sc1 );
	best_i_vor_1( yp, xp, sqrtdelta, sc2 );
	I_ADD( sc1, sc2, sum );
	out[0] = 0.5*sum[0];
	out[1] = 0.5*sum[1];
}


void i_octavor( double y[12], double x[12], 
	double sqrtdelta[2], double sol[2], double out[2] )
{
	double xp[12], yp[12], sc1[2], sc2[2], sum[2];
	int i, j, k;
	
	for( i=0; i<2; i++ ) {
		j = i;
		k = i;
		xp[j] = x[k];
		yp[j] = y[k];
		j = 2 + i;
		k = 8 + i;
		xp[j] = x[k];
		yp[j] = y[k];
		j = 4 + i;
		k = 10 + i;
		xp[j] = x[k];
		yp[j] = y[k];
		j = 6 + i;
		k = j;
		xp[j] = x[k];
		yp[j] = y[k];
		j = 8 + i;
		k = 2 + i;
		xp[j] = x[k];
		yp[j] = y[k];
		j = 10 + i;
		k = 4 + i;
		xp[j] = x[k];
		yp[j] = y[k];
	}
	i_vor_alt( x, sqrtdelta, sol, sc1 );
	i_vor( yp, xp, sqrtdelta, sc2 );
	I_ADD( sc1, sc2, sum );
	out[0] = 0.5*sum[0];
	out[1] = 0.5*sum[1];
}


/* Solid angle of a Roger's simplex arising in Voronoi
truncation:
	sph = 2*atan( sqrt( ((b-a)*(c-b))/((a+b)*(b+c)) ) ),
where c = Sqrt[2].  */

void i_rog_sph( double ab[4], double out[2] )
{
	double a[2], b[2], bp[2], bpsign[2], num[2], den[2];
	int i, gotbpart;
	
	for( i=0; i<2; i++ ) {
		a[i] = ab[i];
		b[i] = ab[i + 2];
	}
		
	/* Decreasing in a, but the b partial depends on
	the sign of a c - b^2, where c = Sqrt[2].  
	In addition, the function is increasing in c.  */
	
	ROUND_DOWN;
	bpsign[0] = a[0]*SQRT2_LO - b[1]*b[1];
	ROUND_UP;
	bpsign[1] = a[1]*SQRT2_HI - b[0]*b[0];
	
	gotbpart = 0;
	if( bpsign[0] > 0.0 ) {	/* partial is positive */
		bp[0] = b[0];
		bp[1] = b[1];
		gotbpart = 1;
	} 
	else if( bpsign[1] < 0.0 ) { /* partial is negative */
		bp[0] = b[1];
		bp[1] = b[0];
		gotbpart = 1;
	}
	if( gotbpart ) {	/* if we got a bounded partial */
		ROUND_DOWN;
		num[0] = (bp[0] - a[1])*(SQRT2_LO - bp[0]);
		den[1] = (a[0] + bp[1])*(bp[1] + SQRT2_HI);
		ROUND_UP;
		num[1] = (bp[1] - a[0])*(SQRT2_HI - bp[1]);
		den[0] = (a[1] + bp[0])*(bp[0] + SQRT2_LO);
	} else {	/* no bounded partial for b */
		ROUND_DOWN;
		num[0] = (b[0] - a[1])*(SQRT2_LO - b[1]);
		den[1] = (a[0] + b[0])*(b[0] + SQRT2_HI);
		ROUND_UP;
		num[1] = (b[1] - a[0])*(SQRT2_HI - b[0]);
		den[0] = (a[1] + b[1])*(b[1] + SQRT2_LO);
	}
	bp[1] = num[1]/den[1];
	ROUND_DOWN;
	bp[0] = num[0]/den[0];
	I_SQRT( bp, a );
	I_ATAN( a, b );
	out[0] = 2.0*b[0];
	out[1] = 2.0*b[1];
#if DEBUGTRUNC
	ROUND_NEAR;
	printf("sph = [%.18g, %.18g]\n", out[0], out[1]);
#endif
}


/* vol = alpha*a*(2 - a^2)/6, which is decreasing in a.  */
void i_wedge_vol( double alpha[2], double a[2], double out[2] )
{
	double a2[2];
	
	ROUND_UP;
	a2[1] = a[1]*a[1];
	ROUND_DOWN;
	a2[0] = a[0]*a[0];
	out[0] = alpha[0]*a[1]*(2.0 - a2[1])*ONE_6_LO;
	ROUND_UP;
	out[1] = alpha[1]*a[0]*(2.0 - a2[0])*ONE_6_HI;
}


/* In the obtuse case, it is possible that alpha will be
negative.  Need to deal with this. */

void obtuse_wedge_vol( double alpha[2], double a[2], 
	double out[2] )
{
	double a2[2], temp[2];
	
	ROUND_UP;
	a2[1] = a[1]*a[1];
	ROUND_DOWN;
	a2[0] = a[0]*a[0];
	temp[0] = a[1]*(2.0 - a2[1])*ONE_6_LO;
	if( alpha[0] >= 0.0 ) {
		out[0] = alpha[0]*temp[0];
	}
	ROUND_UP;
	temp[1] = a[0]*(2.0 - a2[0])*ONE_6_HI;
	if( alpha[1] >= 0.0 ) {
		out[1] = alpha[1]*temp[1];
	} else {
		out[1] = alpha[0]*temp[0];
		ROUND_DOWN;
		out[0] = alpha[1]*temp[1];
	}
}


/* sph = (1 - a Sqrt[2]/2)*alpha, also decreasing in a.  */

void i_wedge_sph( double alpha[2], double a[2], double out[2] )
{
	double temp[2];
	
#if GOT_FMADD
	temp[0] = 0.5*a[0];
	temp[1] = 0.5*a[1];
	ROUND_UP;
	out[1] = (1.0 - temp[0]*SQRT2_LO)*alpha[1];
	ROUND_DOWN;
	out[0] = (1.0 - temp[1]*SQRT2_HI)*alpha[0];
#else
	ROUND_DOWN;
	temp[1] = a[0]*SQRT2_LO;
	ROUND_UP;
	temp[0] = a[1]*SQRT2_HI;
	out[1] = (1.0 - 0.5*temp[1])*alpha[1];
	ROUND_DOWN;
	out[0] = (1.0 - 0.5*temp[0])*alpha[0];
#endif
}


/* In the obtuse case, it is possible that alpha will be
negative.  Need to deal with this. */

void obtuse_wedge_sph( double alpha[2], double a[2], 
	double out[2] )
{
	double temp[2], temp2[2];
	
	temp[0] = 0.5*a[0];
	temp[1] = 0.5*a[1];
	ROUND_DOWN;
	temp2[0] = 1.0 + (-temp[1])*SQRT2_HI;
	if( alpha[0] >= 0.0 ) {
		out[0] = temp2[0]*alpha[0];
	}
	ROUND_UP;
	temp2[1] = 1.0 + (-temp[0])*SQRT2_LO;
	if( alpha[1] >= 0.0 ) {
		out[1] = temp2[1]*alpha[1];
	} else {
		out[1] = temp2[0]*alpha[0];
		ROUND_DOWN;
		out[0] = temp2[1]*alpha[1];
	}
}


/* i_dih_rog = atan( sqrt( (2 - b^2)/(b^2 - a^2) ) ).
Note that dih is increasing in a, decreasing in b.  Use
x = a^2, y = b^2.  This gives
	atan( sqrt( (2 - y)/(y - x) ) ) 
Assume that x < y.  */

void i_dih_rog( double xy[4], double out[2] )
{
	double num[2], den[2], temp;
	
	ROUND_DOWN;
	num[0] = 2.0 - xy[3];
	den[0] = xy[2] - xy[1];
	ROUND_UP;
	num[1] = 2.0 - xy[2];
	den[1] = xy[3] - xy[0];
	
	temp = num[1]/den[0];
	out[1] = atan( sqrt( temp ) ) + ATANERR;

	ROUND_DOWN;
	temp = num[0]/den[1];
	out[0] = atan( sqrt( temp ) ) - ATANERR;
}


/* i_dih_rog = acos( sqrt( (b^2 - a^2)/(2 - a^2) ) ).
Note that dih is increasing in a, decreasing in b.  Use
x = a^2, y = b^2.  This gives
	acos( sqrt( (y - x)/(2 - x) ) ) 
Assume that x < y.  */

void old_dih_rog( double xy[4], double out[2] )
{
	double num[2], den[2], temp, temp2[2];
	
#if DEBUGTRUNC
	ROUND_NEAR;
	printf("x = [%.18g, %.18g]\n", xy[0], xy[1]);
	printf("y = [%.18g, %.18g]\n", xy[2], xy[3]);
#endif

	ROUND_DOWN;
	num[1] = xy[2] - xy[1];
	den[0] = 2.0 - xy[0];
	ROUND_UP;
	num[0] = xy[3] - xy[0];
	den[1] = 2.0 - xy[1];
	temp = num[0]/den[0];
	temp2[0] = sqrt( temp );
	if( temp2[0] > 1.0 )
		temp2[0] = 1.0;
	ROUND_DOWN;
	temp = num[1]/den[1];
	temp2[1] = sqrt( temp );
	out[0] = acos( temp2[0] ) - ATANERR;
	ROUND_DOWN;
	if( temp2[1] > 1.0 )
		temp2[1] = 1.0;
	out[1] = acos( temp2[1] ) + ATANERR;
#if DEBUGTRUNC
	ROUND_NEAR;
	printf("cos([]) = [%.18g, %.18g]\n", temp2[0], temp2[1]);
#endif
}


/* Voronoi truncation.  This is pure Voronoi scoring, where
the Voronoi cell is intersected with a sphere of radius
Sqrt[2], to account for possible truncation.  I assume
that the long edge has length >= 2 Sqrt[2], since
otherwise two quarters would be formed instead.  In addition,
I assume that no side clipping is possible.  */

void i_vor_trunc( double y[12], double x[12], 
	double sqrtdelta[2], double out[2] )
{
	double ab[4], xyz[6], be[6], al[4], al3[4], ga[6];
	double facex[6], vol[2], sol[2];
	double eta[2], eta2[2];
	double temp[2], temp2[2];
	double xp[12], foo;
	int i;
	
	
	/* (2 3 4) (-> (2 4 6)) */
	for( i=0; i<2; i++ ) {
		facex[i] = x[2 + i];
		facex[2 + i] = x[4 + i];
		facex[4 + i] = x[6 + i];
	}
	s_crad3x2( facex, eta2 );
	
	/* Watch out for NaN's */
	if( ISNAN( eta2[0] ) || ISNAN( eta2[1] ) ) {
		out[0] = eta2[0] + eta2[1];
		out[1] = out[0];
		return;
	}
	
	vol[0] = 0.0;
	vol[1] = 0.0;
	sol[0] = 0.0;
	sol[1] = 0.0;
	xyz[4] = 2.0;
	xyz[5] = 2.0;
	
	/* Check to see if Roger's appear:  */
	if( eta2[0] < 2.0 ) {
		I_SQRT( eta2, eta );
		
		/* first face */
		xyz[2] = eta2[0];
		xyz[3] = eta2[1];
		ab[2] = eta[0];
		ab[3] = eta[1];
		
		/* Now get two Roger's simplices: */
		
		/* Edge 2: (2) */
		ab[0] = 0.5*y[2];
		ab[1] = 0.5*y[3];
		xyz[0] = 0.25*x[2];
		xyz[1] = 0.25*x[3];
		
		i_rogersvol2( xyz, temp2 );
		I_SQRT( temp2, temp );
		i_rog_sph( ab, temp2 );
		
		ROUND_DOWN;
		vol[0] += temp[0];
		sol[0] += temp2[0];
		ROUND_UP;
		vol[1] += temp[1];
		sol[1] += temp2[1];
		
		i_dih_rog( xyz, al + 2 );
		
		/* Edge 3: (4) */
		ab[0] = 0.5*y[4];
		ab[1] = 0.5*y[5];
		xyz[0] = 0.25*x[4];
		xyz[1] = 0.25*x[5];

		i_rogersvol2( xyz, temp2 );
		I_SQRT( temp2, temp );
		i_rog_sph( ab, temp2 );
		
		ROUND_DOWN;
		vol[0] += temp[0];
		sol[0] += temp2[0];
		ROUND_UP;
		vol[1] += temp[1];
		sol[1] += temp2[1];

		i_dih_rog( xyz, al3 + 2 );
		
	} else {
		al[2] = 0.0;
		al[3] = 0.0;
		al3[2] = 0.0;
		al3[3] = 0.0;
	}
	
	/* (1 3 5) (-> (0 4 8)) */
	for( i=0; i<2; i++ ) {
		facex[i] = x[i];
		/*	facex[2 + i] = x[4 + i];	*/
		facex[4 + i] = x[8 + i];
	}
	s_crad3x2( facex, eta2 );
	
	/* Watch out for NaN's */
	if( ISNAN( eta2[0] ) || ISNAN( eta2[1] ) ) {
		out[0] = eta2[0] + eta2[1];
		out[1] = out[0];
		return;
	}
	
	/* Check to see if Roger's appear:  */
	if( eta2[0] < 2.0 ) {
		I_SQRT( eta2, eta );
		
		/* second face */
		xyz[2] = eta2[0];
		xyz[3] = eta2[1];
		ab[2] = eta[0];
		ab[3] = eta[1];
		
		/* Now get two Roger's simplices: */
		
		/* Edge 1: (0) */
		ab[0] = 0.5*y[0];
		ab[1] = 0.5*y[1];
		xyz[0] = 0.25*x[0];
		xyz[1] = 0.25*x[1];
		
		i_rogersvol2( xyz, temp2 );
		I_SQRT( temp2, temp );
		i_rog_sph( ab, temp2 );
		
		ROUND_DOWN;
		vol[0] += temp[0];
		sol[0] += temp2[0];
		ROUND_UP;
		vol[1] += temp[1];
		sol[1] += temp2[1];
		
		i_dih_rog( xyz, al );
		
		/* Edge 3: (4) */
		ab[0] = 0.5*y[4];
		ab[1] = 0.5*y[5];
		xyz[0] = 0.25*x[4];
		xyz[1] = 0.25*x[5];

		i_rogersvol2( xyz, temp2 );
		I_SQRT( temp2, temp );
		i_rog_sph( ab, temp2 );
		
		ROUND_DOWN;
		vol[0] += temp[0];
		sol[0] += temp2[0];
		ROUND_UP;
		vol[1] += temp[1];
		sol[1] += temp2[1];

		i_dih_rog( xyz, al3 );
		
	} else {
		al[0] = 0.0;
		al[1] = 0.0;
		al3[0] = 0.0;
		al3[1] = 0.0;
	}

#if DEBUGTRUNC
	printf("spherical angle subtotal (from rogers)\n is [%.20g, %.20g]\n",
	sol[0], sol[1]);
	printf("volume subtotal (from rogers)\n is [%.20g, %.20g]\n",
	vol[0], vol[1]);
#endif
	/* Compute dihedral angles for the spherical wedges. */
	
	s_dih( x, be );
	
	/* edge 2 (2) */	/*	(2 1 3 5 4 6) -> (2 0 4 8 6 10)	*/
	for( i=0; i<2; i++ ) {
		xp[i] = x[2 + i];
		xp[2 + i] = x[i];
		xp[4 + i] = x[4 + i];
		xp[6 + i] = x[8 + i];
		xp[8 + i] = x[6 + i];
		xp[10 + i] = x[10 + i];
	}
	s_dih( xp, be + 2 );
	
	/* edge 3 (4) */	/*	(3 2 1 6 5 4) -> (4 2 0 10 8 6)	*/
	for( i=0; i<2; i++ ) {
		xp[i] = x[4 + i];
		xp[2 + i] = x[2 + i];
		xp[4 + i] = x[i];
		xp[6 + i] = x[10 + i];
		xp[8 + i] = x[8 + i];
		xp[10 + i] = x[6 + i];
	}
	s_dih( xp, be + 4 );
	
	I_SUB( be, al, ga );
	i_sub( be + 2, al + 2, ga + 2 );
	i_sub( be + 4, al3, temp );
	i_sub( temp, al3 + 2, ga + 4 );
	
#if DEBUGTRUNC	
	ROUND_NEAR;
	printf(" be = \n");
	for( i=0; i<6; i+=2 ) {
		printf("\t[%.20g, %.20g]\n", be[i], be[i+1]);
		if( be[i+1] - be[i] < 0.0 )
			printf("BADNESS, diff = %.20g\n", be[i+1] - be[i]);
	}
	printf(" al = \n");
	for( i=0; i<4; i+=2 ) {
		printf("\t[%.20g, %.20g]\n", al[i], al[i+1]);
		if( al[i+1] - al[i] < 0.0 )
			printf("BADNESS, diff = %.20g\n", al[i+1] - al[i]);
	}
	printf(" al3 = \n");
	for( i=0; i<4; i+=2 ) {
		printf("\t[%.20g, %.20g]\n", al3[i], al3[i+1]);
		if( al3[i+1] - al3[i] < 0.0 )
			printf("BADNESS, diff = %.20g\n", al3[i+1] - al3[i]);
	}
	printf(" ga = \n");
	for( i=0; i<6; i+=2 ) {
		printf("\t[%.20g, %.20g]\n", ga[i], ga[i+1]);
		if( ga[i+1] - ga[i] < 0.0 )
			printf("BADNESS, diff = %.20g\n", ga[i+1] - ga[i]);
	}
#endif
	
	
	/* Now we have the dihedral angles, so we can compute
	the spherical angle and volume for the spherical wedges,
	adding them to the total. */
	
#if DEBUGTRUNC
	printf("sol is [%.20g, %.20g]\n",
		sol[0], sol[1]);
	if( sol[1] - sol[0] < 0.0 )
		printf("BADNESS, diff = %.20g\n", sol[1] - sol[0]);
#endif

	for( i=0; i<6; i+=2 ) {
		ab[0] = 0.5*y[i];
		ab[1] = 0.5*y[i+1];
		i_wedge_vol( ga + i, ab, temp );
		i_wedge_sph( ga + i, ab, temp2 );
#if DEBUGTRUNC
		ROUND_NEAR;
		printf("i_wedge_vol = [%.20g, %.20g]\n", temp[0], temp[1]);
		printf("i_wedge_sph = [%.20g, %.20g]\n", 
			temp2[0], temp2[1]);
#endif
		ROUND_DOWN;
		vol[0] += temp[0];
		sol[0] += temp2[0];
		ROUND_UP;
		vol[1] += temp[1];
		sol[1] += temp2[1];
#if DEBUGTRUNC
		printf("sol is [%.20g, %.20g]\n",
			sol[0], sol[1]);
		if( sol[1] - sol[0] < 0.0 )
			printf("BADNESS, diff = %.20g\n", sol[1] - sol[0]);
#endif
	}
	
#if DEBUGTRUNC
	printf("sol is [%.20g, %.20g]\n",
		sol[0], sol[1]);
	if( sol[1] - sol[0] < 0.0 )
		printf("BADNESS, diff = %.20g\n", sol[1] - sol[0]);
#endif

	/* Find the spherical angle of the central spherical
	section, and compute its volume. */
	
	i_solid( y, sqrtdelta, eta2 );
	I_SUB( eta2, sol, eta );

#if DEBUGTRUNC
	printf("central spherical angle is [%.20g, %.20g]\n",
		eta[0], eta[1]);
	if( eta[1] - eta[0] < 0.0 )
		printf("BADNESS, diff = %.20g\n", eta[1] - eta[0]);
#endif
	
	ROUND_DOWN;
	foo = TWOSQRT2_LO*eta[0];
	temp[0] = foo*ONE_3_LO;
	vol[0] += temp[0];
	
	ROUND_UP;
	foo = TWOSQRT2_HI*eta[1];
	temp[1] = foo*ONE_3_HI;
	vol[1] += temp[1];
	
#if DEBUGTRUNC
	printf("central chunk volume is [%.20g, %.20g]\n", 
		temp[0], temp[1]);
	if( temp[1] - temp[0] < 0.0 )
		printf("BADNESS, diff = %.20g\n", temp[1] - temp[0]);
	printf("truncated voronoi volume is [%.20g, %.20g]\n", 
		vol[0], vol[1]);
	if( vol[1] - vol[0] < 0.0 )
		printf("BADNESS, diff = %.20g\n", vol[1] - vol[0]);
	printf("spherical angle is [%.20g, %.20g]\n",
		eta2[0], eta2[1]);
	if( eta2[1] - eta2[0] < 0.0 )
		printf("BADNESS, diff = %.20g\n", eta2[1] - eta2[0]);
#endif

	/* Spherical angle of tet is now eta2 */
	
	ROUND_DOWN;
	temp[0] = DOCT_LO*vol[0];
	temp2[0] = ONE_3_LO*eta2[0];
	ROUND_UP;
	temp[1] = DOCT_HI*vol[1];
	temp2[1] = ONE_3_HI*eta2[1];
	
	eta[1] = temp2[1] - temp[0];
	ROUND_DOWN;
	eta[0] = temp2[0] - temp[1];
	
	out[0] = 4.0*eta[0];
	out[1] = 4.0*eta[1];
}


/* Here's the cleaned-up version. */
void alt_vor_trunc( double y[12], double x[12], 
	double sph_in[2], double out[2] )
{
	double ab[4], xyz[6], be[6], al[8], ga[6];
	double facex[6], vol[2], sol[2];
	double eta[2], eta2[2];
	double temp[2], temp2[2];
	double xp[12], foo;
	int i, j, k, ell;
	
	for( i=0; i<2; i++ ) {
		vol[i] = 0.0;
		sol[i] = 0.0;
	}
	
	xyz[4] = 2.0;
	xyz[5] = 2.0;
	
	for( k=0; k<4; k+=2 ) {
		switch( k ) {
			case 0:		/* first face */
				/* (2 3 4) (-> (2 4 6)) */
				for( i=0; i<2; i++ ) {
					facex[i] = x[2 + i];
					facex[2 + i] = x[4 + i];
					facex[4 + i] = x[6 + i];
				}
				s_crad3x2( facex, eta2 );
				break;
			case 2:		/* second face */
				/* (1 3 5) (-> (0 4 8)) */
				for( i=0; i<2; i++ ) {
					facex[i] = x[i];
					/*	facex[2 + i] = x[4 + i];	*/
					facex[4 + i] = x[8 + i];
				}
				s_crad3x2( facex, eta2 );
				break;
			default:
				break;
		}	
		/* Watch out for NaN's */
		if( ISNAN( eta2[0] ) || ISNAN( eta2[1] ) ) {
			out[0] = eta2[0] + eta2[1];
			out[1] = out[0];
			return;
		}
		/* Check to see if Roger's appear:  */
		if( eta2[0] < 2.0 ) {
			I_SQRT( eta2, eta );
			
			xyz[2] = eta2[0];
			xyz[3] = eta2[1];
			ab[2] = eta[0];
			ab[3] = eta[1];
			
			for( j=0; j<2; j++ ) {
			/* Now get two Roger's simplices: */
				ell = j + k;
				switch( ell ) {
					case 0:		/* Edge 2: (2) */
						ab[0] = 0.5*y[2];
						ab[1] = 0.5*y[3];
						xyz[0] = 0.25*x[2];
						xyz[1] = 0.25*x[3];
						break;
					case 2:		/* Edge 1: (0) */
						ab[0] = 0.5*y[0];
						ab[1] = 0.5*y[1];
						xyz[0] = 0.25*x[0];
						xyz[1] = 0.25*x[1];
						break;
					default:	/* Edge 3: (4) */
						ab[0] = 0.5*y[4];
						ab[1] = 0.5*y[5];
						xyz[0] = 0.25*x[4];
						xyz[1] = 0.25*x[5];
						break;
				}
				i_rogersvol2( xyz, temp2 );
				I_SQRT( temp2, temp );
				i_rog_sph( ab, temp2 );
				i_dih_rog( xyz, al + 2*ell );
				
				ROUND_DOWN;
				vol[0] += temp[0];
				sol[0] += temp2[0];
				
				ROUND_UP;
				vol[1] += temp[1];
				sol[1] += temp2[1];
			}
		}
	}
#if DEBUGTRUNC
	printf("vol = [%.20g, %.20g]\n", vol[0], vol[1]);
	printf("sol = [%.20g, %.20g]\n", sol[0], sol[1]);
#endif
	
	/* Compute dihedral angles for the spherical wedges. */
	
	s_dih( x, be );
	
	/* edge 2 (2) */	/*	(2 1 3 5 4 6) -> (2 0 4 8 6 10)	*/
	for( i=0; i<2; i++ ) {
		xp[i] = x[2 + i];
		xp[2 + i] = x[i];
		xp[4 + i] = x[4 + i];
		xp[6 + i] = x[8 + i];
		xp[8 + i] = x[6 + i];
		xp[10 + i] = x[10 + i];
	}
	s_dih( xp, be + 2 );
	
	/* edge 3 (4) */	/*	(3 2 1 6 5 4) -> (4 2 0 10 8 6)	*/
	for( i=0; i<2; i++ ) {
		xp[i] = x[4 + i];
		xp[2 + i] = x[2 + i];
		xp[4 + i] = x[i];
		xp[6 + i] = x[10 + i];
		xp[8 + i] = x[8 + i];
		xp[10 + i] = x[6 + i];
	}
	s_dih( xp, be + 4 );
	
	ROUND_DOWN;
	ga[4] = be[4] - al[3] - al[7];
	if( ga[4] < 0.0 )	/* Ensure ga >= 0 */
		ga[4] = 0.0;
	ga[0] = be[0] - al[5];
	if( ga[0] < 0.0 )
		ga[0] = 0.0;
	ga[2] = be[2] - al[1];
	if( ga[2] < 0.0 )
		ga[2] = 0.0;
	ROUND_UP;
	ga[5] = be[5] - al[2] - al[6];
	ga[1] = be[1] - al[4];
	ga[3] = be[3] - al[0];

	/* Now we have the dihedral angles, so we can compute
	the spherical angle and volume for the spherical wedges. */
	
	for( i=0; i<6; i+=2 ) {
		ab[0] = 0.5*y[i];
		ab[1] = 0.5*y[i+1];
		i_wedge_vol( ga + i, ab, temp );
		i_wedge_sph( ga + i, ab, temp2 );
		ROUND_DOWN;
		vol[0] += temp[0];
		sol[0] += temp2[0];
		ROUND_UP;
		vol[1] += temp[1];
		sol[1] += temp2[1];
	}
	
	/* Compute the solid angle of the tetrahedron. */
	
	/* i_solid( y, sqrtdelta, sol_out ); */

	/* Find the spherical angle of the central spherical
	section, and compute its volume. */
	
	I_SUB( sph_in, sol, temp );

#if DEBUGTRUNC
	printf("central spherical angle is [%.20g, %.20g]\n",
		temp[0], temp[1]);
	if( temp[1] - temp[0] < 0.0 )
		printf("BADNESS, diff = %.20g\n", temp[1] - temp[0]);
#endif
	
	ROUND_DOWN;
	foo = TWOSQRT2_LO*temp[0];
	vol[0] += foo*ONE_3_LO;
	
	ROUND_UP;
	foo = TWOSQRT2_HI*temp[1];
	vol[1] += foo*ONE_3_HI;

	/* Now make the final computation. */
	
	/* Formula:  sc = -doct*( sum( fullwedvol ) - 
			sum( wedvol ) + sum( rogvol ) ) + 2Sqrt[2]/3 doct*(
			sum( fullwedsol ) - sum( wedsol ) + sum( rogsol ) ) +
			(1-2Sqrt[2] doct)/3 sol
	*/

	ROUND_DOWN;
	temp[0] = DOCT_LO*vol[0];
	temp2[0] = ONE_3_LO*sph_in[0];
	ROUND_UP;
	temp[1] = DOCT_HI*vol[1];
	temp2[1] = ONE_3_HI*sph_in[1];
	
	eta[1] = temp2[1] - temp[0];
	ROUND_DOWN;
	eta[0] = temp2[0] - temp[1];
	
	out[0] = 4.0*eta[0];
	out[1] = 4.0*eta[1];
}


/* Here's the cleaned-up version. */
void obtuse_vor_trunc( double y[12], double x[12], 
	double sph_in[2], double out[2] )
{
	double ab[4], xyz[6], be[6], al[8], ga[6];
	double facex[6], vol[2], sol[2];
	double eta[2], eta2[2];
	double temp[2], temp2[2];
	double xp[12], foo;
	int i, j, k, ell;
	
	for( i=0; i<2; i++ ) {
		vol[i] = 0.0;
		sol[i] = 0.0;
	}
	
	xyz[4] = 2.0;
	xyz[5] = 2.0;
	
	for( k=0; k<4; k+=2 ) {
		switch( k ) {
			case 0:		/* first face */
				/* (2 3 4) (-> (2 4 6)) */
				for( i=0; i<2; i++ ) {
					facex[i] = x[2 + i];
					facex[2 + i] = x[4 + i];
					facex[4 + i] = x[6 + i];
				}
				s_crad3x2( facex, eta2 );
				break;
			case 2:		/* second face */
				/* (1 3 5) (-> (0 4 8)) */
				for( i=0; i<2; i++ ) {
					facex[i] = x[i];
					/*	facex[2 + i] = x[4 + i];	*/
					facex[4 + i] = x[8 + i];
				}
				s_crad3x2( facex, eta2 );
				break;
			default:
				break;
		}	
		/* Watch out for NaN's */
		if( ISNAN( eta2[0] ) || ISNAN( eta2[1] ) ) {
			out[0] = eta2[0] + eta2[1];
			out[1] = out[0];
			return;
		}
		/* Check to see if Roger's appear:  */
		if( eta2[0] < 2.0 ) {
			I_SQRT( eta2, eta );
			
			xyz[2] = eta2[0];
			xyz[3] = eta2[1];
			ab[2] = eta[0];
			ab[3] = eta[1];
			
			for( j=0; j<2; j++ ) {
			/* Now get two Roger's simplices: */
				ell = j + k;
				switch( ell ) {
					case 0:		/* Edge 2: (2) */
						ab[0] = 0.5*y[2];
						ab[1] = 0.5*y[3];
						xyz[0] = 0.25*x[2];
						xyz[1] = 0.25*x[3];
						break;
					case 2:		/* Edge 1: (0) */
						ab[0] = 0.5*y[0];
						ab[1] = 0.5*y[1];
						xyz[0] = 0.25*x[0];
						xyz[1] = 0.25*x[1];
						break;
					default:	/* Edge 3: (4) */
						ab[0] = 0.5*y[4];
						ab[1] = 0.5*y[5];
						xyz[0] = 0.25*x[4];
						xyz[1] = 0.25*x[5];
						break;
				}
				i_rogersvol2( xyz, temp2 );
				I_SQRT( temp2, temp );
				i_rog_sph( ab, temp2 );
				i_dih_rog( xyz, al + 2*ell );
				
				ROUND_DOWN;
				vol[0] += temp[0];
				sol[0] += temp2[0];
				
				ROUND_UP;
				vol[1] += temp[1];
				sol[1] += temp2[1];
			}
		}
	}
#if DEBUGTRUNC
	printf("vol = [%.20g, %.20g]\n", vol[0], vol[1]);
	printf("sol = [%.20g, %.20g]\n", sol[0], sol[1]);
#endif
	
	/* Compute dihedral angles for the spherical wedges. */
	
	s_dih( x, be );
	
	/* edge 2 (2) */	/*	(2 1 3 5 4 6) -> (2 0 4 8 6 10)	*/
	for( i=0; i<2; i++ ) {
		xp[i] = x[2 + i];
		xp[2 + i] = x[i];
		xp[4 + i] = x[4 + i];
		xp[6 + i] = x[8 + i];
		xp[8 + i] = x[6 + i];
		xp[10 + i] = x[10 + i];
	}
	s_dih( xp, be + 2 );
	
	/* edge 3 (4) */	/*	(3 2 1 6 5 4) -> (4 2 0 10 8 6)	*/
	for( i=0; i<2; i++ ) {
		xp[i] = x[4 + i];
		xp[2 + i] = x[2 + i];
		xp[4 + i] = x[i];
		xp[6 + i] = x[10 + i];
		xp[8 + i] = x[8 + i];
		xp[10 + i] = x[6 + i];
	}
	s_dih( xp, be + 4 );
	
	ROUND_DOWN;
	ga[4] = be[4] - al[3] - al[7];
	ga[0] = be[0] - al[5];
	ga[2] = be[2] - al[1];
	ROUND_UP;
	ga[5] = be[5] - al[2] - al[6];
	ga[1] = be[1] - al[4];
	ga[3] = be[3] - al[0];

	/* Now we have the dihedral angles, so we can compute
	the spherical angle and volume for the spherical wedges. */
	
	for( i=0; i<6; i+=2 ) {
		ab[0] = 0.5*y[i];
		ab[1] = 0.5*y[i+1];
		obtuse_wedge_vol( ga + i, ab, temp );
		obtuse_wedge_sph( ga + i, ab, temp2 );
		ROUND_DOWN;
		vol[0] += temp[0];
		sol[0] += temp2[0];
		ROUND_UP;
		vol[1] += temp[1];
		sol[1] += temp2[1];
	}
	
	/* Compute the solid angle of the tetrahedron. */
	
	/* i_solid( y, sqrtdelta, sol_out ); */

	/* Find the spherical angle of the central spherical
	section, and compute its volume. */
	
	I_SUB( sph_in, sol, temp );

#if DEBUGTRUNC
	printf("central spherical angle is [%.20g, %.20g]\n",
		temp[0], temp[1]);
	if( temp[1] - temp[0] < 0.0 )
		printf("BADNESS, diff = %.20g\n", temp[1] - temp[0]);
#endif
	
	ROUND_DOWN;
	foo = TWOSQRT2_LO*temp[0];
	vol[0] += foo*ONE_3_LO;
	
	ROUND_UP;
	foo = TWOSQRT2_HI*temp[1];
	vol[1] += foo*ONE_3_HI;

	/* Now make the final computation. */
	
	/* Formula:  sc = -doct*( sum( fullwedvol ) - 
			sum( wedvol ) + sum( rogvol ) ) + 2Sqrt[2]/3 doct*(
			sum( fullwedsol ) - sum( wedsol ) + sum( rogsol ) ) +
			(1-2Sqrt[2] doct)/3 sol
	*/

	ROUND_DOWN;
	temp[0] = DOCT_LO*vol[0];
	temp2[0] = ONE_3_LO*sph_in[0];
	ROUND_UP;
	temp[1] = DOCT_HI*vol[1];
	temp2[1] = ONE_3_HI*sph_in[1];
	
	eta[1] = temp2[1] - temp[0];
	ROUND_DOWN;
	eta[0] = temp2[0] - temp[1];
	
	out[0] = 4.0*eta[0];
	out[1] = 4.0*eta[1];
}


void test_rog_wed( double y[12], double x[12], 
	double vol[2], double sol[2] )
{
	double ab[4], xyz[6], al[8];
	double facex[6];
	double eta[2], eta2[2];
	double temp[2], temp2[2];
	int i, j, k, ell;
	
	for( i=0; i<2; i++ ) {
		vol[i] = 0.0;
		sol[i] = 0.0;
	}
	
	xyz[4] = 2.0;
	xyz[5] = 2.0;
	
	for( k=0; k<4; k+=2 ) {
		switch( k ) {
			case 0:		/* first face */
				/* (2 3 4) (-> (2 4 6)) */
				for( i=0; i<2; i++ ) {
					facex[i] = x[2 + i];
					facex[2 + i] = x[4 + i];
					facex[4 + i] = x[6 + i];
				}
				s_crad3x2( facex, eta2 );
				break;
			case 2:		/* second face */
				/* (1 3 5) (-> (0 4 8)) */
				for( i=0; i<2; i++ ) {
					facex[i] = x[i];
					/*	facex[2 + i] = x[4 + i];	*/
					facex[4 + i] = x[8 + i];
				}
				s_crad3x2( facex, eta2 );
				break;
			default:
				break;
		}	
		/* Check to see if Roger's appear:  */
		if( eta2[0] < 2.0 ) {
			I_SQRT( eta2, eta );
			
			xyz[2] = eta2[0];
			xyz[3] = eta2[1];
			ab[2] = eta[0];
			ab[3] = eta[1];
			
			for( j=0; j<2; j++ ) {
			/* Now get two Roger's simplices: */
				ell = j + k;
				switch( ell ) {
					case 0:		/* Edge 2: (2) */
						ab[0] = 0.5*y[2];
						ab[1] = 0.5*y[3];
						xyz[0] = 0.25*x[2];
						xyz[1] = 0.25*x[3];
						break;
					case 2:		/* Edge 1: (0) */
						ab[0] = 0.5*y[0];
						ab[1] = 0.5*y[1];
						xyz[0] = 0.25*x[0];
						xyz[1] = 0.25*x[1];
						break;
					default:	/* Edge 3: (4) */
						ab[0] = 0.5*y[4];
						ab[1] = 0.5*y[5];
						xyz[0] = 0.25*x[4];
						xyz[1] = 0.25*x[5];
						break;
				}
				i_rogersvol2( xyz, temp2 );
				I_SQRT( temp2, temp );
				i_rog_sph( ab, temp2 );
				i_dih_rog( xyz, al + 2*ell );
				
				ROUND_DOWN;
				vol[0] += temp[0];
				sol[0] += temp2[0];
				
				ROUND_UP;
				vol[1] += temp[1];
				sol[1] += temp2[1];
			}
		}
	}
}


void test_fullwed( double y[12], double x[12], 
	double vol[2], double sol[2] )
{
	double ab[4], be[6];
	double temp[2], temp2[2];
	double xp[12];
	int i;
	
	for( i=0; i<2; i++ ) {
		vol[i] = 0.0;
		sol[i] = 0.0;
	}
	
	/* Compute dihedral angles for the spherical wedges. */
	
	s_dih( x, be );
	
	/* edge 2 (2) */	/*	(2 1 3 5 4 6) -> (2 0 4 8 6 10)	*/
	for( i=0; i<2; i++ ) {
		xp[i] = x[2 + i];
		xp[2 + i] = x[i];
		xp[4 + i] = x[4 + i];
		xp[6 + i] = x[8 + i];
		xp[8 + i] = x[6 + i];
		xp[10 + i] = x[10 + i];
	}
	s_dih( xp, be + 2 );
	
	/* edge 3 (4) */	/*	(3 2 1 6 5 4) -> (4 2 0 10 8 6)	*/
	for( i=0; i<2; i++ ) {
		xp[i] = x[4 + i];
		xp[2 + i] = x[2 + i];
		xp[4 + i] = x[i];
		xp[6 + i] = x[10 + i];
		xp[8 + i] = x[8 + i];
		xp[10 + i] = x[6 + i];
	}
	s_dih( xp, be + 4 );
	
	/* Now we have the dihedral angles, so we can compute
	the spherical angle and volume for the spherical wedges. */
	
	for( i=0; i<6; i+=2 ) {
		ab[0] = 0.5*y[i];
		ab[1] = 0.5*y[i+1];
		i_wedge_vol( be + i, ab, temp );
		i_wedge_sph( be + i, ab, temp2 );
		ROUND_DOWN;
		vol[0] += temp[0];
		sol[0] += temp2[0];
		ROUND_UP;
		vol[1] += temp[1];
		sol[1] += temp2[1];
	}
}

