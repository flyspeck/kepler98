/* second_partials.c, by Samuel Ferguson, (c) 1997. */
/* Routines for computing second partials of VorVolAnalytic, 
etc. */

/* 970113:  Change routines so we can do everything in terms
of x.  This should significantly improve the performance
(when used appropriately) of the taylor expansion routines. */


#include "system_headers.h"
#include "i_sphere.h"
#include "interval.h"
#include "i_bounds.h"
#include "second_partials.h"
#include "macros.h"

#define DEBUG		0
#define DEBUG2	0

#define NOLONG	0		/* some edge lengths may exceed 2sqrt(2) */


/* To simplify all routines, I will assume that the index of
all second partials lies in the upper triangle:  that is,
i <= j */

/*
afunction[x__]:=
x[1]*(x[2] + x[6] - x[1])
*/

void sp_afunction( double x[12], double apars[12], 
	double out[2] )
{
	double xp[12], xpp[12];
	int i;
	
	for( i=0; i<12; i++ ) {
		xp[i] = x[i];
		xpp[i] = x[i];
	}
	
	if( apars[0] >= 0.0 ) {
		xp[1] = x[0];		/* for min */
		xpp[0] = x[1];	/* for max */
	} 
	else if( apars[1] <= 0.0 ) {
		xp[0] = x[1];		/* for min */
		xpp[1] = x[0];	/* for max */
	}
	
	ROUND_DOWN;
	out[0] = xp[0]*(xp[2] + xp[10] - xp[1]);
	ROUND_UP;
	out[1] = xpp[1]*(xpp[3] + xpp[11] - xpp[0]);
}


/*
newafun[x__]:=
	-x[[1]]^2 + 2 x[[1]] x[[2]] - x[[2]]^2 + 
	x[[1]] x[[6]] + x[[2]] x[[6]]
*/
void sp_newafun( double x[12], double apars[12], 
	double out[2] )
{
	double xp[12], xpp[12], pterms[2], nterms[2];
	int i,ip;
	
	for( i=0; i<12; i++ ) {
		xp[i] = x[i];
		xpp[i] = x[i];
	}
	
	for( i=0; i<4; i+=2 ) {
		ip = i + 1;
		if( apars[i] >= 0.0 ) {
			xp[ip] = x[i];		/* for min */
			xpp[i] = x[ip];	/* for max */
		} 
		else if( apars[ip] <= 0.0 ) {
			xp[i] = x[ip];		/* for min */
			xpp[ip] = x[i];	/* for max */
		}
	}
	ROUND_DOWN;
	pterms[0] = 2.0*xp[0]*xp[2] + xp[0]*xp[10] + xp[2]*xp[10];
	nterms[0] = xpp[0]*xpp[0] + xpp[2]*xpp[2];
	ROUND_UP;
	pterms[1] = 2.0*xpp[1]*xpp[3] + xpp[1]*xpp[11] +
								xpp[3]*xpp[11];
	nterms[1] = xp[1]*xp[1] + xp[3]*xp[3];
	out[1] = pterms[1] - nterms[0];
	ROUND_DOWN;
	out[0] = pterms[0] - nterms[1];
}


/*
bfunction[x__]:=
	tomschi[{x[4],x[5],x[3],x[1],x[2],x[6]}]

bfunction[xtet] :=
		-(x[1]^2*x[4]) - x[2]^2*x[5] - 2*x[1]*x[2]*x[6] + 
		(x[1] + x[2])*x[3]*x[6] - x[3]*x[6]^2 + 
		x[2]*x[5]*(x[1] + x[6]) + x[1]*x[4]*(x[2] + x[6])
*/

void sp_bfunction( double x[12], double bpars[12], 
	double out[2] )
{
	double pterms[2], nterms[2], xp[12], xpp[12];
	int i, j;
	
	for( i=0; i<12; i++ ) {
		xp[i] = x[i];
		xpp[i] = x[i];
	}
	
	for( i=0; i<12; i+=2 ) {
		j = i + 1;
		if( bpars[i] >= 0.0 ) {
			xp[j] = x[i];		/* for min */
			xpp[i] = x[j];	/* for max */
		} 
		else if( bpars[j] <= 0.0 ) {
			xp[i] = x[j];		/* for min */
			xpp[j] = x[i];	/* for max */
		}
	}
	
	ROUND_DOWN;
	pterms[0] = xp[4]*xp[10]*(xp[0] + xp[2]) + 
		xp[2]*xp[8]*(xp[0] + xp[10]) + 
		xp[0]*xp[6]*(xp[2] + xp[10]);
	nterms[0] = xpp[0]*xpp[0]*xpp[6] + xpp[2]*xpp[2]*xpp[8] + 
		2.0*xpp[0]*xpp[2]*xpp[10] + xpp[4]*xpp[10]*xpp[10];
	ROUND_UP;
	pterms[1] = xpp[5]*xpp[11]*(xpp[1] + xpp[3]) + 
		xpp[3]*xpp[9]*(xpp[1] + xpp[11]) + 
		xpp[1]*xpp[7]*(xpp[3] + xpp[11]);
	nterms[1] = xp[1]*xp[1]*xp[7] + xp[3]*xp[3]*xp[9] + 
		2.0*xp[1]*xp[3]*xp[11] + xp[5]*xp[11]*xp[11];
	out[1] = pterms[1] - nterms[0];
	ROUND_DOWN;
	out[0] = pterms[0] - nterms[1];
}


/* cfunction = tomsu */
void sp_cfunction( double x[12], double cpars[12],
	double out[2] )
{
	double pterms[2], nterms[2], xp[12], xpp[12];
	int i, j;
	
	for( i=0; i<4; i++ ) {
		xp[i] = x[i];
		xpp[i] = x[i];
	}
	for( i=10; i<12; i++ ) {
		xp[i] = x[i];
		xpp[i] = x[i];
	}
	
	for( i=0; i<4; i+=2 ) {
		j = i + 1;
		if( cpars[i] >= 0.0 ) {
			xp[j] = x[i];		/* for min */
			xpp[i] = x[j];	/* for max */
		} 
		else if( cpars[j] <= 0.0 ) {
			xp[i] = x[j];		/* for min */
			xpp[j] = x[i];	/* for max */
		}
	}
	for( i=10; i<12; i+=2 ) {
		j = i + 1;
		if( cpars[i] >= 0.0 ) {
			xp[j] = x[i];		/* for min */
			xpp[i] = x[j];	/* for max */
		} 
		else if( cpars[j] <= 0.0 ) {
			xp[i] = x[j];		/* for min */
			xpp[j] = x[i];	/* for max */
		}
	}

	ROUND_DOWN;
	pterms[0] = 2.0*(xp[0]*(xp[10] + xp[2]) + xp[2]*xp[10]);
	nterms[0] = xpp[0]*xpp[0] + xpp[2]*xpp[2] + xpp[10]*xpp[10];
	ROUND_UP;
	pterms[1] = 2.0*(xpp[1]*(xpp[11] + xpp[3]) + xpp[3]*xpp[11]);
	nterms[1] = xp[1]*xp[1] + xp[3]*xp[3] + xp[11]*xp[11];
	out[1] = pterms[1] - nterms[0];
	ROUND_DOWN;
	out[0] = pterms[0] - nterms[1];
	if( out[0] < 0.0 )
		out[0] = 0.0;
}


void sp_dfunction( double x[12], double delta_part[12],
	double out_delta[2], double out_sqrtdelta[2] )
{
	double xp[12], data[2];
	
	int i, j;
	
	/* Generate best possible bounds for bigdelta */
		/* s_delta_partials(  x, delta_part ); */

	/* First find delta_min */
	for( i=0; i<12; i++ )
		xp[i] = x[i];		/* copy x, then modify copy */
	j = 0;
	for( i=0; i<6; i++ ) {
		data[0] = delta_part[j];
		data[1] = delta_part[j+1];
		/*	printf("partials data: %lf\t%lf\n", data[0], data[1]);	*/
		if( data[0] > 0 )	/* positive sign */
			xp[j+1] = xp[j];
		if( data[1] < 0 )	/* negative sign */
			xp[j] = xp[j+1];
		j += 2;
	}
	out_delta[0] = rough_min_delta( xp );
	
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
	out_delta[1] = rough_max_delta( xp );
	
	i_sqrt( out_delta, out_sqrtdelta );
}


void sp_ufunction( double aval[2], double bval[2], 
	double out[2] )
{
	I_MULT( aval, bval, out );
}


/* (u_i*v - u*v_i)/v^2 */
void sp_ubv_i( double uval[2], double vval[2],
	double uval_i[12], double vval_i[12], double out[12] )
{
	double den[2], oden[2], num[2], temp[2], temp2[2];
	int i;
	
	temp[0] = 1.0;
	temp[1] = 1.0;
	i_mult( vval, vval, den );
	i_div( temp, den, oden );
	
	for( i=0; i<12; i+=2 ) {
		i_mult( uval_i + i, vval, temp );
		i_mult( uval, vval_i + i, temp2 );
		I_SUB( temp, temp2, num );
		i_mult( num, oden, out + i );
	}
}


/* u' = a'*b + a*b' */
void sp_u_i( double aval[2], double bval[2],
	double aval_i[12], double bval_i[12], double out[12] )
{
	double t1[2], t2[2];
	int i;
	
	for( i=0; i<12; i+=2 ) {
		i_mult( aval_i + i, bval, t1 );
		i_mult( aval, bval_i + i, t2 );
		i_add( t1, t2, out + i );
	}
}


/* u_ij = a_ij*b + a_i*b_j + a_j*b_i + a*b_ij */
void sp_u_ij( double aval[2], double bval[2],
	double aval_i[12], double bval_i[12],
	double aval_ij[6][12], double bval_ij[6][12], 
	double out[6][12] )
{
	double t1[2], t2[2], t3[2];
	int i, j, m, n;
	
	for( i=0; i<6; i++ ) {
		m = 2*i;
		for( j=i; j<6; j++ ) {
			n = 2*j;
			i_mult( aval_ij[i] + n, bval, t1 );
			i_mult( aval_i + m, bval_i + n, t2 );
			I_ADD( t1, t2, t3 );
			i_mult( aval_i + n, bval_i + m, t1 );
			I_ADD( t1, t3, t2 );
			i_mult( aval, bval_ij[i] + n, t1 );
			i_add( t2, t1, out[i] + n );
		}
	}
}


/* ( v*(u_ij*v + u_i*v_j - u_j*v_i - u*v_ij) + 
		2*v_j*(u*v_i - u_i*v) )/v^3   */
void sp_ubv_ij( double u[2], double v[2],
	double u_i[12], double v_i[12],
	double u_ij[6][12], double v_ij[6][12], 
	double out[6][12] )
{
	double t1[2], t2[2], t3[2], t4[2], num[2], oden[2];
	int i, j, m, n;
	
	/* This is a serious mess.  But, it can't really be
	helped.  Kinda hard to check.  */
	
	i_mult( v, v, t1 );							/* v*v			-> t1 */
	i_mult( v, t1, t2 );						/* v*t1			-> t2 */
	t3[0] = 1.0;
	t3[1] = 1.0;
	i_div( t3, t2, oden );
	
	for( i=0; i<6; i++ ) {
		m = 2*i;
		for( j=i; j<6; j++ ) {
			n = 2*j;
			i_mult( u_ij[i] + n, v, t1 );		/* u_ij*v 	-> t1 */
			i_mult( u_i + m, v_i + n, t2 );	/* u_i*v_j 	-> t2 */
			I_ADD( t1, t2, t3 );						/* t1 + t2 	-> t3 */
			i_mult( u_i + n, v_i + m, t1 );	/* u_j*v_i	-> t1 */
			I_SUB( t3, t1, t2 );						/* t3 - t1	-> t2 */
			i_mult( u, v_ij[i] + n, t1 );		/* u*v_ij		-> t1 */
			I_SUB( t2, t1, t3 );						/* t2 - t1	-> t3 */
			i_mult( v, t3, t4 );						/* v*t3			-> t4 */
			
			i_mult( u, v_i + m, t1 );				/* u*v_i		-> t1 */
			i_mult( u_i + m, v, t2 );				/* u_i*v		-> t2 */
			I_SUB( t1, t2, t3 );						/* t1 - t2	-> t3 */
			i_mult( v_i + n, t3, t2 );			/* v_j*t3		-> t2 */
			t1[0] = 2.0*t2[0];							/* 2*t2			-> t1 */
			t1[1] = 2.0*t2[1];
			
			I_ADD( t4, t1, num );						/* t4 + t1	-> num */
			
			i_mult( num, oden, out[i] + n );/* num*oden	-> out_ij */
		}
	}
}


/*
apars = 
{-2*x[1] + x[2] + x[6], x[1], 0, 0, 0, x[1]}
*/

void sp_apars( double x[12], double out[12] )
{
	int i, n;
	
	for( i=0; i<6; i++ ) {
		n = 2*i;
		switch( i ) {
			case 0:
				ROUND_DOWN;
				out[n] = -2.0*x[1] + x[2] + x[10];
				ROUND_UP;
				out[n+1] = -2.0*x[0] + x[3] + x[11];
				break;
			case 1:
				out[n] = x[0];
				out[n+1] = x[1];
				break;
			case 5:
				out[n] = x[0];
				out[n+1] = x[1];
				break;
			default:
				out[n] = 0.0;
				out[n+1] = 0.0;
				break;
		}
	}
}


/*
newapars = 
{-2 x[1] + 2 x[2] + x[6], 2 x[1] - 2 x[2] + x[6], 0, 0, 0, 
  x[1] + x[2]}
*/
void sp_newapars( double x[12], double out[12] )
{
	int i, n;
	
	for( i=0; i<6; i++ ) {
		n = 2*i;
		switch( i ) {
			case 0:
				ROUND_DOWN;
				out[n] = 2.0*(x[2] - x[1]) + x[10];
				ROUND_UP;
				out[n+1] = 2.0*(x[3] - x[0]) + x[11];
				break;
			case 1:
				ROUND_DOWN;
				out[n] = 2.0*(x[0] - x[3]) + x[10];
				ROUND_UP;
				out[n+1] = 2.0*(x[1] - x[2]) + x[11];
				break;
			case 5:
				ROUND_DOWN;
				out[n] = x[0] + x[2];
				ROUND_UP;
				out[n+1] = x[1] + x[3];
				break;
			default:
				out[n] = 0.0;
				out[n+1] = 0.0;
				break;
		}
	}
}


/*
bpars = 
	{-2*x[1]*x[4] + x[2]*x[5] - 2*x[2]*x[6] + x[3]*x[6] + 
   x[4]*(x[2] + x[6]), x[1]*x[4] - 2*x[2]*x[5] - 
   2*x[1]*x[6] + x[3]*x[6] + x[5]*(x[1] + x[6]), 
  (x[1] + x[2])*x[6] - x[6]^2, 
  -x[1]^2 + x[1]*(x[2] + x[6]), 
  -x[2]^2 + x[2]*(x[1] + x[6]), 
  -2*x[1]*x[2] + (x[1] + x[2])*x[3] + x[1]*x[4] + 
   x[2]*x[5] - 2*x[3]*x[6]}
*/
/*
bsecondpars = 	(- x + x + x)		(these appear to be useless)
								(x - + + x x)
								(+ + 0 0 0 x)
								(x + 0 0 0 +)
								(+ x 0 0 0 +)
								(x x x + + -)
								 0 2 4 6 8 10
*/

void sp_bpars( double x[12], double bsecpars[6][12], 
	double out[12] )
{
	double pterms[2], nterms[2], xp[12], xpp[12];
	int i, j, m, n;
	
	for( m=0; m<6; m++ ) {
		n = 2*m;
		
		for( i=0; i<12; i++ ) {
			xp[i] = x[i];
			xpp[i] = x[i];
		}
		for( i=0; i<12; i+=2 ) {
			j = i + 1;
			if( bsecpars[m][i] >= 0.0 ) {
				xp[j] = x[i];		/* for min */
				xpp[i] = x[j];	/* for max */
			} 
			else if( bsecpars[m][j] <= 0.0 ) {
				xp[i] = x[j];		/* for min */
				xpp[j] = x[i];	/* for max */
			}
		}
		switch( m ) {
			case 0:	/* -2*x[1]*x[4] + x[2]*x[5] - 2*x[2]*x[6] + 
							x[3]*x[6] + x[4]*(x[2] + x[6]) */
				ROUND_DOWN;
				pterms[0] = xp[2]*xp[8] + xp[4]*xp[10] + 
					xp[6]*(xp[2] + xp[10]);
				nterms[0] = 2.0*(xpp[0]*xpp[6] + xpp[2]*xpp[10]);
				ROUND_UP;
				pterms[1] = xpp[3]*xpp[9] + xpp[5]*xpp[11] + 
					xpp[7]*(xpp[3] + xpp[11]);
				nterms[1] = 2.0*(xp[1]*xp[7] + xp[3]*xp[11]);
				break;
			case 1:	/* x[1]*x[4] - 2*x[2]*x[5] - 2*x[1]*x[6] + 
							x[3]*x[6] + x[5]*(x[1] + x[6]) */
				ROUND_DOWN;
				pterms[0] = xp[0]*xp[6] + xp[4]*xp[10] + 
					xp[8]*(xp[0] + xp[10]);
				nterms[0] = 2.0*(xpp[2]*xpp[8] + xpp[0]*xpp[10]);
				ROUND_UP;
				pterms[1] = xpp[1]*xpp[7] + xpp[5]*xpp[11] + 
					xpp[9]*(xpp[1] + xpp[11]);
				nterms[1] = 2.0*(xp[3]*xp[9] + xp[1]*xp[11]);
				break;
			case 2:	/* (x[1] + x[2])*x[6] - x[6]^2 */
				ROUND_DOWN;
				pterms[0] = xp[10]*(xp[0] + xp[2]);
				nterms[0] = xpp[10]*xpp[10];
				ROUND_UP;
				pterms[1] = xpp[11]*(xpp[1] + xpp[3]);
				nterms[1] = xp[11]*xp[11];
				break;
			case 3:	/* -x[1]^2 + x[1]*(x[2] + x[6]) */
				ROUND_DOWN;
				pterms[0] = xp[0]*(xp[2] + xp[10]);
				nterms[0] = xpp[0]*xpp[0];
				ROUND_UP;
				pterms[1] = xpp[1]*(xpp[3] + xpp[11]);
				nterms[1] = xp[1]*xp[1];
				break;
			case 4:	/* -x[2]^2 + x[2]*(x[1] + x[6]) */
				ROUND_DOWN;
				pterms[0] = xp[2]*(xp[0] + xp[10]);
				nterms[0] = xpp[2]*xpp[2];
				ROUND_UP;
				pterms[1] = xpp[3]*(xpp[1] + xpp[11]);
				nterms[1] = xp[3]*xp[3];
				break;
			case 5:	/* -2*x[1]*x[2] + (x[1] + x[2])*x[3] + 
							x[1]*x[4] + x[2]*x[5] - 2*x[3]*x[6] */
				ROUND_DOWN;
				pterms[0] = xp[4]*(xp[0] + xp[2]) + xp[0]*xp[6] + 
					xp[2]*xp[8];
				nterms[0] = 2.0*(xpp[0]*xpp[2] + xpp[4]*xpp[10]);
				ROUND_UP;
				pterms[1] = xpp[5]*(xpp[1] + xpp[3]) + xpp[1]*xpp[7] +
					xpp[3]*xpp[9];
				nterms[1] = 2.0*(xp[1]*xp[3] + xp[5]*xp[11]);
				break;
		}
		out[n+1] = pterms[1] - nterms[0];
		ROUND_DOWN;
		out[n] = pterms[0] - nterms[1];
	}
}


/*
cpars = 
	{-2*x[1] + 2*x[2] + 2*x[6], 2*x[1] - 2*x[2] + 2*x[6], 0, 
  0, 0, 2*x[1] + 2*x[2] - 2*x[6]}
*/

void sp_cpars( double x[12], double out[12] )
{
	double temp[2];
	int m, n;
	
	for( m=0; m<6; m++ ) {
		n = 2*m;
		switch( m ) {
			case 0:
				ROUND_DOWN;
				temp[0] = -x[1] + x[2] + x[10];
				ROUND_UP;
				temp[1] = -x[0] + x[3] + x[11];
				out[n] = 2.0*temp[0];
				out[n+1] = 2.0*temp[1];
				break;
			case 1:
				ROUND_DOWN;
				temp[0] = x[0] - x[3] + x[10];
				ROUND_UP;
				temp[1] = x[1] - x[2] + x[11];
				out[n] = 2.0*temp[0];
				out[n+1] = 2.0*temp[1];
				break;
			case 5:
				ROUND_DOWN;
				temp[0] = x[0] + x[2] - x[11];
				ROUND_UP;
				temp[1] = x[1] + x[3] - x[10];
				out[n] = 2.0*temp[0];
				out[n+1] = 2.0*temp[1];
				break;
			default:
				out[n] = 0.0;
				out[n+1] = 0.0;
				break;
		}
	}
}


/* d_i = delta_i/(2*sqrt(delta)) */
void sp_dpars( double delta_pars[12],
	double sqrt_delta[2], double out[12] )
{
	double o2sqdelta[2];
	int j;
	
	ROUND_DOWN;
	o2sqdelta[0] = 0.5/sqrt_delta[1];
	ROUND_UP;
	o2sqdelta[1] = 0.5/sqrt_delta[0];
	
	for( j=0; j<12; j+=2 ) {
		i_mult( delta_pars + j, o2sqdelta, out + j );
	}
}


/*
deltapars = 

{-2*x[1]*x[4] + x[2]*x[4] + x[3]*x[4] - x[4]^2 + 
   x[2]*x[5] - x[3]*x[5] + x[4]*x[5] - x[2]*x[6] + 
   x[3]*x[6] + x[4]*x[6], 
  x[1]*x[4] - x[3]*x[4] + x[1]*x[5] - 2*x[2]*x[5] + 
   x[3]*x[5] + x[4]*x[5] - x[5]^2 - x[1]*x[6] + 
   x[3]*x[6] + x[5]*x[6], 
  x[1]*x[4] - x[2]*x[4] - x[1]*x[5] + x[2]*x[5] + 
   x[1]*x[6] + x[2]*x[6] - 2*x[3]*x[6] + x[4]*x[6] + 
   x[5]*x[6] - x[6]^2, -x[1]^2 + x[1]*x[2] + x[1]*x[3] - 
   x[2]*x[3] - 2*x[1]*x[4] + x[1]*x[5] + x[2]*x[5] + 
   x[1]*x[6] + x[3]*x[6] - x[5]*x[6], 
  x[1]*x[2] - x[2]^2 - x[1]*x[3] + x[2]*x[3] + x[1]*x[4] + 
   x[2]*x[4] - 2*x[2]*x[5] + x[2]*x[6] + x[3]*x[6] - 
   x[4]*x[6], -(x[1]*x[2]) + x[1]*x[3] + x[2]*x[3] - 
   x[3]^2 + x[1]*x[4] + x[3]*x[4] + x[2]*x[5] + 
   x[3]*x[5] - x[4]*x[5] - 2*x[3]*x[6]}
*/

/* Caveat:  no long ( > 2Sqrt[2] ) edges.  
	Get this from i_bounds.c. */

void sp_deltapars( double x[12], double delta_ij[6][12],
	double part[12] )
{
	double pterms[2], nterms[2], xp[12], xpp[12];
	double xrp[12], xrpp[12], temp[2];
	int i, j, jp, k, kp, n, ind;
	int indexing[6][6] = {{6, 8, 4, 0, 2, 10},
						{8, 10, 0, 2, 4, 6},
						{10, 0, 8, 4, 6, 2},
						{0, 2, 4, 6, 8, 10},
						{2, 10, 6, 8, 4, 0},
						{4, 0, 2, 10, 6, 8}};
		
	ind = 0;
	for( n=0; n<6; n++ ) {
		for( i=0; i<12; i++ ) {
			xrp[i] = x[i];
			xrpp[i] = x[i];
		}
		for( i=0; i<12; i+=2 ) {
			j = i + 1;
			if( delta_ij[n][i] >= 0.0 ) {
				xrp[j] = x[i];		/* for min */
				xrpp[i] = x[j];	/* for max */
			} 
			else if( delta_ij[n][j] <= 0.0 ) {
				xrp[i] = x[j];		/* for min */
				xrpp[j] = x[i];	/* for max */
			}
		}
		/* Now switch everything around. */
		k = 0;
		for( i=0; i<6; i++ ) {
			j = indexing[n][i];
			xp[k] = xrp[j];
			xpp[k] = xrpp[j];
			kp = k + 1;
			jp = j + 1;
			xp[kp] = xrp[jp];
			xpp[kp] = xrpp[jp];
			k += 2;
		}
/* -x[1]^2 + x[1]*x[2] + x[1]*x[3] - 
   x[2]*x[3] - 2*x[1]*x[4] + x[1]*x[5] + x[2]*x[5] + 
   x[1]*x[6] + x[3]*x[6] - x[5]*x[6] */
		ROUND_DOWN;
		pterms[0] = xp[0]*(xp[2] + xp[4] + xp[8] + xp[10]) +
			xp[2]*xp[8] + xp[4]*xp[10];
		nterms[0] = xpp[0]*(xpp[0] + 2.0*xpp[6]) + xpp[2]*xpp[4] + 
			xpp[8]*xpp[10];
		ROUND_UP;
		pterms[1] = xpp[1]*(xpp[3] + xpp[5] + xpp[9] + xpp[11]) +
			xpp[3]*xpp[9] + xpp[5]*xpp[11];
		nterms[1] = xp[1]*(xp[1] + 2.0*xp[7]) + xp[3]*xp[5] + 
			xp[9]*xp[11];
		temp[1] = pterms[1] - nterms[0];
		ROUND_DOWN;
		temp[0] = pterms[0] - nterms[1];
		part[ind] = temp[0];
		ind++;
		part[ind] = temp[1];
		ind++;
	}
}


/*
asecondpars = 
{{-2, 1, 0, 0, 0, 1}, {1, 0, 0, 0, 0, 0}, 
  {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, 
  {0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0}}
*/

void sp_asecondpars( double bigout[6][12] )
{
	double out[2];
	int i, j, n;
	
	for( i=0; i<6; i++ ) {
		for( j=i; j<6; j++ ) {
			n = 2*j;
			if( i==0 ) {
				if( j==0 ) {
					out[0] = -2.0;
					out[1] = -2.0;
				} 
				else if( j==1 || j==5 ) {
					out[0] = 1.0;
					out[1] = 1.0;
				}
				else {
					out[0] = 0.0;
					out[1] = 0.0;
				}
			}
			else {
				out[0] = 0.0;
				out[1] = 0.0;
			}
			bigout[i][n] = out[0];
			bigout[i][n+1] = out[1];
		}
	}
}


/*
newasecpars = 
{{-2, 2, 0, 0, 0, 1}, {2, -2, 0, 0, 0, 1}, 
  {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, 
  {0, 0, 0, 0, 0, 0}, {1, 1, 0, 0, 0, 0}}
*/
void sp_newasecpars( double bigout[6][12] )
{
	double out[2];
	int i, j, n;
	
	for( i=0; i<6; i++ ) {
		for( j=i; j<6; j++ ) {
			n = 2*j;
			if( i==0 ) {
				if( j==0 ) {
					out[0] = -2.0;
					out[1] = -2.0;
				} 
				else if( j==1 ) {
					out[0] = 2.0;
					out[1] = 2.0;
				}
				else if( j==5 ) {
					out[0] = 1.0;
					out[1] = 1.0;
				}
				else {
					out[0] = 0.0;
					out[1] = 0.0;
				}
			}
			else if( i==1 ) {
				if( j==1 ) {
					out[0] = -2.0;
					out[1] = -2.0;
				} else if( j==5 ) {
					out[0] = 1.0;
					out[1] = 1.0;
				}
				else {
					out[0] = 0.0;
					out[1] = 0.0;
				}
			}
			else {
				out[0] = 0.0;
				out[1] = 0.0;
			}
			bigout[i][n] = out[0];
			bigout[i][n+1] = out[1];
		}
	}
}


/*
bsecondpars = 

		{{-2*x[4], 
		x[4] + x[5] - 2*x[6], 
		x[6], 
		-2*x[1] + x[2] + x[6], 
		x[2], 
		-2*x[2] + x[3] + x[4]},

		{x[4] + x[5] - 2*x[6], 
		-2*x[5], 
		x[6], 
		x[1], 
		x[1] - 2*x[2] + x[6], 
		-2*x[1] + x[3] + x[5]},

		{x[6], 
		x[6], 
		0, 
		0, 
		0, 
		x[1] + x[2] - 2*x[6]},

		{-2*x[1] + x[2] + x[6], 
		x[1], 
		0, 
		0, 
		0,
		x[1]},

		{x[2], 
		x[1] - 2*x[2] + x[6], 
		0, 
		0, 
		0, 
		x[2]},

		{-2*x[2] + x[3] + x[4], 
		-2*x[1] + x[3] + x[5], 
		x[1] + x[2] - 2*x[6], 
		x[1], x[2], 
		-2*x[3]}}
*/

void sp_bsecondpars( double x[12], double bigout[6][12] )
{
	double out[2];
	int i, j, n;
	
	for( i=0; i<6; i++ ) {
		for( j=i; j<6; j++ ) {
			n = 2*j;
			switch( i ) {
		/*	{-2*x[4], 
				x[4] + x[5] - 2*x[6], 
				x[6], 
			   -2*x[1] + x[2] + x[6], 
			   x[2], 
			   -2*x[2] + x[3] + x[4]} */
				case 0:
					switch( j ) {
						case 0:
							out[0] = -2.0*x[7];
							out[1] = -2.0*x[6];
							break;
						case 1:
							ROUND_DOWN;
							out[0] = x[6] + x[8] - 2.0*x[11];
							ROUND_UP;
							out[1] = x[7] + x[9] - 2.0*x[10];
							break;
						case 2:
							out[0] = x[10];
							out[1] = x[11];
							break;
						case 3:
							ROUND_DOWN;
							out[0] = -2.0*x[1] + x[2] + x[10];
							ROUND_UP;
							out[1] = -2.0*x[0] + x[3] + x[11];
							break;
						case 4:
							out[0] = x[2];
							out[1] = x[3];
							break;
						case 5:
							ROUND_DOWN;
							out[0] = -2.0*x[3] + x[4] + x[6];
							ROUND_UP;
							out[1] = -2.0*x[2] + x[5] + x[7];
							break;
					}
					break;
				case 1:
		/* 	{x[4] + x[5] - 2*x[6], 
			  -2*x[5], 
			  x[6], 
			  x[1], 
			  x[1] - 2*x[2] + x[6], 
			  -2*x[1] + x[3] + x[5]} */
					switch( j ) {
						case 1:
							out[0] = -2.0*x[9];
							out[1] = -2.0*x[8];
							break;
						case 2:
							out[0] = x[10];
							out[1] = x[11];
							break;
						case 3:
							out[0] = x[0];
							out[1] = x[1];
							break;
						case 4:
							ROUND_DOWN;
							out[0] = x[0] - 2.0*x[3] + x[10];
							ROUND_UP;
							out[1] = x[1] - 2.0*x[2] + x[11];
							break;
						case 5:
							ROUND_DOWN;
							out[0] = -2.0*x[1] + x[4] + x[8];
							ROUND_UP;
							out[1] = -2.0*x[0] + x[5] + x[9];
							break;
					}
					break;
				case 2:
		/*  {x[6], 
			  x[6], 
			  0, 
			  0, 
			  0, 
			  x[1] + x[2] - 2*x[6]} */
					switch( j ) {
						case 5:
							ROUND_DOWN;
							out[0] = x[0] + x[2] - 2.0*x[11];
							ROUND_UP;
							out[1] = x[1] + x[3] - 2.0*x[10];
							break;
						default:
							out[0] = 0.0;
							out[1] = 0.0;
							break;
					}
					break;
				case 3:
		/*  {-2*x[1] + x[2] + x[6], 
			  x[1], 
			  0, 
			  0, 
			  0,
			  x[1]} */
					switch( j ) {
						case 5:
							out[0] = x[0];
							out[1] = x[1];
							break;
						default:
							out[0] = 0.0;
							out[1] = 0.0;
							break;
					}
					break;
				case 4:
		/*  {x[2], 
			  x[1] - 2*x[2] + x[6], 
			  0, 
			  0, 
			  0, 
			  x[2]} */
					switch( j ) {
						case 5:
							out[0] = x[2];
							out[1] = x[3];
							break;
						default:
							out[0] = 0.0;
							out[1] = 0.0;
							break;
					}
					break;
				case 5:
		/*  {-2*x[2] + x[3] + x[4], 
			  -2*x[1] + x[3] + x[5], 
			   x[1] + x[2] - 2*x[6], 
			   x[1], x[2], 
			   -2*x[3]} */
				  out[0] = -2.0*x[5];
				  out[1] = -2.0*x[4];
					break;
				default:
					printf("This can't happen.\n");
					break;
			}
		  bigout[i][n] = out[0];
		  bigout[i][n+1] = out[1];
		}
	}
}


/*
csecondpars = 
{{-2, 2, 0, 0, 0, 2}, {2, -2, 0, 0, 0, 2}, 
  {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, 
  {0, 0, 0, 0, 0, 0}, {2, 2, 0, 0, 0, -2}}
*/

void sp_csecondpars( double bigout[6][12] )
{
	double out[2];
	int i, j, n;
	
	/* There's gotta be a better way to fill in these
	matrices.  */
	for( i=0; i<6; i++ ) {
		for( j=i; j<6; j++ ) {
			n = 2*j;
			if( i==0 ) {
				if( j==0 ) {
					out[0] = -2.0;
					out[1] = -2.0;
				}
				else if( j==1 || j==5 ) {
					out[0] = 2.0;
					out[1] = 2.0;
				}
				else {
					out[0] = 0.0;
					out[1] = 0.0;
				}
			}
			else if( i==j ) {
				if( i==1 ) {
					out[0] = -2.0;
					out[1] = -2.0;
				}
				else if( i==5 ) {
					out[0] = -2.0;
					out[1] = -2.0;
				}
				else {
					out[0] = 0.0;
					out[1] = 0.0;
				}
			}
			else if( i==1 && j==5 ) {
					out[0] = 2.0;
					out[1] = 2.0;
			}
			else {
				out[0] = 0.0;
				out[1] = 0.0;
			}
			bigout[i][n] = out[0];
			bigout[i][n+1] = out[0];
		}
	}
}


/* d_ij = (2*delta_ij*delta - 
	delta_i*delta_j)/(4 delta^(3/2)) */
void sp_dsecondpars( double delta[2], double sqrtdelta[2], 
	double delta_pars[12], double delta_ij[6][12],
	double out[6][12] )
{
	double delta_i[2], delta_j[2];
	double num[2], oden[2], den[2], temp[2];
	int i, j, n, m;
	
	ROUND_DOWN;
	temp[0] = 4.0*delta[0]*sqrtdelta[0];
	ROUND_UP;
	temp[1] = 4.0*delta[1]*sqrtdelta[1];
	oden[1] = 1.0/temp[0];
	ROUND_DOWN;
	oden[0] = 1.0/temp[1];
	

	for( n=0; n<6; n++ ) {
		i = 2*n;
		for( m=n; m<6; m++ ) {
			j = 2*m;
			
			delta_i[0] = delta_pars[i];
			delta_i[1] = delta_pars[i+1];
			delta_j[0] = delta_pars[j];
			delta_j[1] = delta_pars[j+1];
			
			i_mult( delta_ij[n] + j, delta, temp );
			den[0] = 2.0*temp[0];
			den[1] = 2.0*temp[1];
			
			i_mult( delta_i, delta_j, temp );
			I_SUB( den, temp, num );
			
			i_mult( num, oden, out[n] + j );
		}
	}
}


/*
deltasecondpars = 
		{{-2*x[4], 
		x[4] + x[5] - x[6], 
		x[4] - x[5] + x[6], 
		-2*x[1] + x[2] + x[3] - 2*x[4] + x[5] + x[6], 
		x[2] - x[3] + x[4], 
		-x[2] + x[3] + x[4]}, 

		{x[4] + x[5] - x[6], 
		-2*x[5], 
		-x[4] + x[5] + x[6], 
		x[1] - x[3] + x[5], 
		x[1] - 2*x[2] + x[3] + x[4] - 2*x[5] + x[6], 
		-x[1] + x[3] + x[5]},

		{x[4] - x[5] + x[6], 
		-x[4] + x[5] + x[6], 
		-2*x[6], 
		x[1] - x[2] + x[6],
		-x[1] + x[2] + x[6], 
		x[1] + x[2] - 2*x[3] + x[4] + x[5] - 2*x[6]},

		{-2*x[1] + x[2] + x[3] - 2*x[4] + x[5] + x[6], 
		x[1] - x[3] + x[5], 
		x[1] - x[2] + x[6], 
		-2*x[1], 
		x[1] + x[2] - x[6], 
		x[1] + x[3] - x[5]}, 

		{x[2] - x[3] + x[4], 
		x[1] - 2*x[2] + x[3] + x[4] - 2*x[5] + x[6], 
		-x[1] + x[2] + x[6], 
		x[1] + x[2] - x[6], 
		-2*x[2], 
		x[2] + x[3] - x[4]}, 

		{-x[2] + x[3] + x[4], 
		-x[1] + x[3] + x[5], 
		x[1] + x[2] - 2*x[3] + x[4] + x[5] - 2*x[6], 
		x[1] + x[3] - x[5], 
		x[2] + x[3] - x[4], 
		-2*x[3]}}
*/

void sp_deltasecondpars( double x[12], double bigout[6][12] )
{
	double out[2];
	int i, j, n;
	
	for( i=0; i<6; i++ ) {
		for( j=i; j<6; j++ ) {
			n = 2*j;
			switch( i ) {
				case 0:
		/*	{-2*x[4], 
				x[4] + x[5] - x[6], 
				x[4] - x[5] + x[6], 
				-2*x[1] + x[2] + x[3] - 2*x[4] + x[5] + x[6], 
				x[2] - x[3] + x[4], 
				-x[2] + x[3] + x[4]} */
					switch( j ) {
						case 0:
							out[0] = -2.0*x[7];
							out[1] = -2.0*x[6];
							break;
						case 1:
							ROUND_DOWN;
							out[0] = x[6] + x[8] - x[11];
							ROUND_UP;
							out[1] = x[7] + x[9] - x[10];
							break;
						case 2:
							ROUND_DOWN;
							out[0] = x[6] - x[9] + x[10];
							ROUND_UP;
							out[1] = x[7] - x[8] + x[11];
							break;
						case 3:
							ROUND_DOWN;
							out[0] = -2.0*x[1] + x[2] + x[4] - 2.0*x[7] + 
								x[8] + x[10];
							ROUND_UP;
							out[1] = -2.0*x[0] + x[3] + x[5] - 2.0*x[6] + 
								x[9] + x[11];
							break;
						case 4:
							ROUND_DOWN;
							out[0] = x[2] - x[5] + x[6];
							ROUND_UP;
							out[1] = x[3] - x[4] + x[7];
							break;
						case 5:
							ROUND_DOWN;
							out[0] = -x[3] + x[4] + x[6];
							ROUND_UP;
							out[1] = -x[2] + x[5] + x[7];
							break;
					}
					break;
				case 1:
		/*	{x[4] + x[5] - x[6], 
				-2*x[5], 
				-x[4] + x[5] + x[6], 
				x[1] - x[3] + x[5], 
				x[1] - 2*x[2] + x[3] + x[4] - 2*x[5] + x[6], 
				-x[1] + x[3] + x[5]} */
					switch( j ) {
						case 1:
							out[0] = -2.0*x[9];
							out[1] = -2.0*x[8];
							break;
						case 2:
							ROUND_DOWN;
							out[0] = -x[7] + x[8] + x[10];
							ROUND_UP;
							out[1] = -x[6] + x[9] + x[11];
							break;
						case 3:
							ROUND_DOWN;
							out[0] = x[0] - x[5] + x[8];
							ROUND_UP;
							out[1] = x[1] - x[4] + x[9];
							break;
						case 4:
							ROUND_DOWN;
							out[0] = x[0] - 2.0*x[3] + x[4] + x[6] - 
								2.0*x[9] + x[10];
							ROUND_UP;
							out[1] = x[1] - 2.0*x[2] + x[5] + x[7] - 
								2.0*x[8] + x[11];
							break;
						case 5:
							ROUND_DOWN;
							out[0] = -x[1] + x[4] + x[8];
							ROUND_UP;
							out[1] = -x[0] + x[5] + x[9];
							break;
					}
					break;
				case 2:
		/*	{x[4] - x[5] + x[6], 
				-x[4] + x[5] + x[6], 
				-2*x[6], 
				x[1] - x[2] + x[6],
				-x[1] + x[2] + x[6], 
				x[1] + x[2] - 2*x[3] + x[4] + x[5] - 2*x[6]} */
					switch( j ) {
						case 2:
							out[0] = -2.0*x[11];
							out[1] = -2.0*x[10];
							break;
						case 3:
							ROUND_DOWN;
							out[0] = x[0] - x[3] + x[10];
							ROUND_UP;
							out[1] = x[1] - x[2] + x[11];
							break;
						case 4:
							ROUND_DOWN;
							out[0] = -x[1] + x[2] + x[10];
							ROUND_UP;
							out[1] = -x[0] + x[3] + x[11];
							break;
						case 5:
							ROUND_DOWN;
							out[0] = x[0] + x[2] - 2.0*x[5] + x[6] + 
								x[8] - 2.0*x[11];
							ROUND_UP;
							out[1] = x[1] + x[3] - 2.0*x[4] + x[7] + 
								x[9] - 2.0*x[10];
							break;
					}
					break;
				case 3:
		/*	{-2*x[1] + x[2] + x[3] - 2*x[4] + x[5] + x[6], 
				x[1] - x[3] + x[5], 
				x[1] - x[2] + x[6], 
				-2*x[1], 
				x[1] + x[2] - x[6], 
				x[1] + x[3] - x[5]} */
					switch( j ) {
						case 3:
							out[0] = -2.0*x[1];
							out[1] = -2.0*x[0];
							break;
						case 4:
							ROUND_DOWN;
							out[0] = x[0] + x[2] - x[11];
							ROUND_UP;
							out[1] = x[1] + x[3] - x[10];
							break;
						case 5:
							ROUND_DOWN;
							out[0] = x[0] + x[4] - x[9];
							ROUND_UP;
							out[1] = x[1] + x[5] - x[8];
							break;
					}
					break;
				case 4:
		/*	{x[2] - x[3] + x[4], 
				x[1] - 2*x[2] + x[3] + x[4] - 2*x[5] + x[6], 
				-x[1] + x[2] + x[6], 
				x[1] + x[2] - x[6], 
				-2*x[2], 
				x[2] + x[3] - x[4]} */
					switch( j ) {
						case 4:
							out[0] = -2.0*x[3];
							out[1] = -2.0*x[2];
							break;
						case 5:
							ROUND_DOWN;
							out[0] = x[2] + x[4] - x[7];
							ROUND_UP;
							out[1] = x[3] + x[5] - x[6];
							break;
					}
					break;
				case 5:
		/*	{-x[2] + x[3] + x[4], 
				-x[1] + x[3] + x[5], 
				x[1] + x[2] - 2*x[3] + x[4] + x[5] - 2*x[6], 
				x[1] + x[3] - x[5], 
				x[2] + x[3] - x[4], 
				-2*x[3]} */
					out[0] = -2.0*x[5];
					out[1] = -2.0*x[4];
					break;
				default:
					printf("This can't happen.\n");
					break;
			}
			bigout[i][n] = out[0];
			bigout[i][n+1] = out[1];
		}
	}
}


/* Here is the master routine.  Hopefully.  */
void sp_secparbounds( double x[12], double out[6][12] )
{
	double u[2], v[2], a[2], b[2], c[2], d[2];
	double u_i[12], v_i[12], a_i[12], b_i[12], c_i[12], d_i[12];
	double a_ij[6][12], b_ij[6][12], c_ij[6][12], d_ij[6][12];
	double u_ij[6][12], v_ij[6][12];
	double delta[2], delta_i[12], delta_ij[6][12];
	
#if DEBUG || DEBUG2
	int i, j;
	double delta_part[12];
#endif


	sp_deltasecondpars( x, delta_ij );
	/* sp_asecondpars( a_ij ); */
	sp_newasecpars( a_ij );
	sp_bsecondpars( x, b_ij );
	sp_csecondpars( c_ij );
	
	/* sp_apars( x, a_i ); */
	sp_newapars( x, a_i );
	sp_bpars( x, b_ij, b_i );
	sp_cpars( x, c_i );
	sp_deltapars( x, delta_ij, delta_i );

	/* sp_afunction( x, a_i, a ); */
	sp_newafun( x, a_i, a );
	sp_bfunction( x, b_i, b );
	sp_cfunction( x, c_i, c );
	sp_dfunction( x, delta_i, delta, d );
	
	sp_dpars( delta_i, d, d_i );
	sp_dsecondpars( delta, d, delta_i, delta_ij, d_ij );

	i_mult( a, b, u );
	i_mult( c, d, v );
	sp_u_i( a, b, a_i, b_i, u_i );
	sp_u_i( c, d, c_i, d_i, v_i );
	sp_u_ij( a, b, a_i, b_i, a_ij, b_ij, u_ij );
	sp_u_ij( c, d, c_i, d_i, c_ij, d_ij, v_ij );
	sp_ubv_ij( u, v, u_i, v_i, u_ij, v_ij, out );
	/* This omits a factor of 1/48 from the definition.  */

#if DEBUG
	ROUND_NEAR;
	printf("u_ij = \n");
	for( i=0; i<6; i++ ) {
		for( j=2*i; j<12; j+=2 ) {
			printf("[%.18g,\t%.18g]\n", 
				out[i][j], out[i][j+1]);
		}
		printf("\n");
	}
#endif

#if DEBUG2
	ROUND_NEAR;
	printf("a = [%.18g,\t%.18g]\n", a[0], a[1]);
	printf("b = [%.18g,\t%.18g]\n", b[0], b[1]);
	printf("c = [%.18g,\t%.18g]\n", c[0], c[1]);
	printf("d = [%.18g,\t%.18g]\n", d[0], d[1]);
	printf("u = [%.18g,\t%.18g]\n", u[0], u[1]);
	printf("v = [%.18g,\t%.18g]\n", v[0], v[1]);
	printf("a_i = \n");
	for( i=0; i<12; i+=2 ) {
		printf("[%.18g,\t%.18g]\n", a_i[i], a_i[i+1]);
	}
	printf("\n");
	printf("b_i = \n");
	for( i=0; i<12; i+=2 ) {
		printf("[%.18g,\t%.18g]\n", b_i[i], b_i[i+1]);
	}
	printf("\n");
	printf("c_i = \n");
	for( i=0; i<12; i+=2 ) {
		printf("[%.18g,\t%.18g]\n", c_i[i], c_i[i+1]);
	}
	printf("\n");
	printf("d_i = \n");
	for( i=0; i<12; i+=2 ) {
		printf("[%.18g,\t%.18g]\n", d_i[i], d_i[i+1]);
	}
	printf("\n");
	printf("delta_i = \n");
	for( i=0; i<12; i+=2 ) {
		printf("[%.18g,\t%.18g]\n", delta_i[i], delta_i[i+1]);
	}
	printf("\n");
	s_delta_partials( x, delta_part );
	printf("delta_part = \n");
	for( i=0; i<12; i+=2 ) {
		printf("[%.18g,\t%.18g]\n", delta_part[i], delta_part[i+1]);
	}
	printf("\n");
	printf("u_i = \n");
	for( i=0; i<12; i+=2 ) {
		printf("[%.18g,\t%.18g]\n", u_i[i], u_i[i+1]);
	}
	printf("\n");
	printf("v_i = \n");
	for( i=0; i<12; i+=2 ) {
		printf("[%.18g,\t%.18g]\n", v_i[i], v_i[i+1]);
	}
	printf("\n");

	printf("a_ij = \n");
	for( i=0; i<6; i++ ) {
		for( j=2*i; j<12; j+=2 ) {
			printf("[%.18g,\t%.18g]\n", 
				a_ij[i][j], a_ij[i][j+1]);
		}
		printf("\n");
	}
	printf("b_ij = \n");
	for( i=0; i<6; i++ ) {
		for( j=2*i; j<12; j+=2 ) {
			printf("[%.18g,\t%.18g]\n", 
				b_ij[i][j], b_ij[i][j+1]);
		}
		printf("\n");
	}
	printf("c_ij = \n");
	for( i=0; i<6; i++ ) {
		for( j=2*i; j<12; j+=2 ) {
			printf("[%.18g,\t%.18g]\n", 
				c_ij[i][j], c_ij[i][j+1]);
		}
		printf("\n");
	}
	printf("d_ij = \n");
	for( i=0; i<6; i++ ) {
		for( j=2*i; j<12; j+=2 ) {
			printf("[%.18g,\t%.18g]\n", 
				d_ij[i][j], d_ij[i][j+1]);
		}
		printf("\n");
	}
	printf("delta_ij = \n");
	for( i=0; i<6; i++ ) {
		for( j=2*i; j<12; j+=2 ) {
			printf("[%.18g,\t%.18g]\n", 
				delta_ij[i][j], delta_ij[i][j+1]);
		}
		printf("\n");
	}
	printf("u_ij = \n");
	for( i=0; i<6; i++ ) {
		for( j=2*i; j<12; j+=2 ) {
			printf("[%.18g,\t%.18g]\n", 
				u_ij[i][j], u_ij[i][j+1]);
		}
		printf("\n");
	}
	printf("v_ij = \n");
	for( i=0; i<6; i++ ) {
		for( j=2*i; j<12; j+=2 ) {
			printf("[%.18g,\t%.18g]\n", 
				v_ij[i][j], v_ij[i][j+1]);
		}
		printf("\n");
	}

#endif
}


/* Need to make sure that the involutions don't pick up
the bottom half of the matrix. */
void sp_vorsecparbds( double x[12], double out[6][12] )
{
	double xp[12], secpar[6][12];
	/* Involutions */
	int inv1[6] = 	{0, 2, 1, 3, 5, 4};
	int tinv1[6] = 	{0, 4, 2, 6, 10, 8};
	int inv2[6] = 	{2, 1, 0, 5, 4, 3};
	int tinv2[6] = 	{4, 2, 0, 10, 8, 6};
	int i, j, m, mp, n, np, a, b, ap;
	
	sp_secparbounds( x, out );
	
	/*
	ROUND_NEAR;
	printf("p1mat = {\n");
	for( i=0; i<6; i++ ) {
		printf("{");
		for( j=2*i; j<12; j+=2 ) {
			printf("%.18g", 
				0.5*(out[i][j] + out[i][j+1]));
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
	*/
	
	/* Case i: */
	for( i=0; i<6; i++ ) {
		m = 2*i;
		mp = m + 1;
		n = tinv1[i];
		np = n + 1;
		xp[m] = x[n];
		xp[mp] = x[np];
	}
	sp_secparbounds( xp, secpar );

	ROUND_DOWN;
	for( i=0; i<6; i++ ) {
		m = 2*i;
		a = inv1[i];
		for( j=i; j<6; j++ ) {
			n = 2*j;
			b = inv1[j];
			ap = a;
			if( b < a ) {	/* swap */
				ap = b;
				b = a;
			}
			out[i][n] += secpar[ap][2*b];
		}
	}
	ROUND_UP;
	for( i=0; i<6; i++ ) {
		m = 2*i + 1;
		a = inv1[i];
		for( j=i; j<6; j++ ) {
			n = 2*j + 1;
			b = inv1[j];
			ap = a;
			if( b < a ) {	/* swap */
				ap = b;
				b = a;
			}
			out[i][n] += secpar[ap][2*b + 1];
		}
	}

	
	/* Case ii: */
	for( i=0; i<6; i++ ) {
		m = 2*i;
		mp = m + 1;
		n = tinv2[i];
		np = n + 1;
		xp[m] = x[n];
		xp[mp] = x[np];
	}
	sp_secparbounds( xp, secpar );
	ROUND_DOWN;
	for( i=0; i<6; i++ ) {
		m = 2*i;
		a = inv2[i];
		for( j=i; j<6; j++ ) {
			n = 2*j;
			b = inv2[j];
			ap = a;
			if( b < a ) {	/* swap */
				ap = b;
				b = a;
			}
			out[i][n] += secpar[ap][2*b];
		}
	}
	ROUND_UP;
	for( i=0; i<6; i++ ) {
		m = 2*i + 1;
		a = inv2[i];
		for( j=i; j<6; j++ ) {
			n = 2*j + 1;
			b = inv2[j];
			ap = a;
			if( b < a ) {	/* swap */
				ap = b;
				b = a;
			}
			out[i][n] += secpar[ap][2*b + 1];
		}
	}
}


void sp_findmaxmin( double ubv_ij[6][12], double out[2] )
{
	double max, min, temp;
	int i, j;

	for( i=0; i<6; i++ ) {
		for( j=2*i; j<12; j+=2 ) {
			temp = ubv_ij[i][j];
			if( min > temp || ISNAN( temp ) )
				min = temp;
			temp = ubv_ij[i][j+1];
			if( max < temp || ISNAN( temp ) )
				max = temp;
		}
	}
	out[0] = min;
	out[1] = max;
}


/*
solasecpars = 
		{{y[2] + y[3], 
		y[1] + y[2] + y[3], 
		y[1] + y[2] + y[3], 
		-y[4], 
		0, 
		0},
		
		{y[1] + y[2] + y[3], 
		y[1] + y[3], 
		y[1] + y[2] + y[3], 
		0, 
		-y[5], 
		0},
		
		{y[1] + y[2] + y[3], 
		y[1] + y[2] + y[3], 
		y[1] + y[2], 
		0, 
		0, 
		-y[6]}, 
		
		{-y[4], 
		0, 
		0, 
		-y[1], 
		0, 
		0}, 
		
		{0, 
		-y[5], 
		0, 
		0, 
		-y[2], 
		0}, 
		
		{0, 
		0, 
		-y[6], 
		0, 
		0, 
		-y[3]}}
*/

void sp_solasecpars( double y[12], double a_ij[6][12] )
{
	int i, j, m, n;
	double out[2], y123[2];
	
	ROUND_DOWN;
	y123[0] = y[0] + y[2] + y[4];
	ROUND_UP;
	y123[1] = y[1] + y[3] + y[5];
	
	for( i=0; i<6; i++ ) {
		m = 2*i;
		for( j=0; j<6; j++ ) {
			n = 2*j;
			switch( i ) {
				case 0:
/*		{y[2] + y[3], 
			y[1] + y[2] + y[3], 
			y[1] + y[2] + y[3], 
			-y[4], 
			0, 
			0} */
					switch( j ) {
						case 0:
							ROUND_DOWN;
							out[0] = y[2] + y[4];
							ROUND_UP;
							out[1] = y[3] + y[5];
							break;
						case 1:
							out[0] = y123[0];
							out[1] = y123[1];
							break;
						case 2:
							out[0] = y123[0];
							out[1] = y123[1];
							break;
						case 3:
							out[0] = -y[7];
							out[1] = -y[6];
							break;
						default:
							out[0] = 0.0;
							out[1] = 0.0;
							break;
					}
					break;
				case 1:
/*			{y[1] + y[2] + y[3], 
				y[1] + y[3], 
				y[1] + y[2] + y[3], 
				0, 
				-y[5], 
				0} */
					switch( j ) {
						case 1:
							ROUND_DOWN;
							out[0] = y[0] + y[4];
							ROUND_UP;
							out[1] = y[1] + y[5];
							break;
						case 2:
							out[0] = y123[0];
							out[1] = y123[1];
							break;
						case 4:
							out[0] = -y[9];
							out[1] = -y[8];
							break;
						default:
							out[0] = 0.0;
							out[1] = 0.0;
							break;
					}
					break;
				case 2:
/*			{y[1] + y[2] + y[3], 
				y[1] + y[2] + y[3], 
				y[1] + y[2], 
				0, 
				0, 
				-y[6]} */
					switch( j ) {
						case 2:
							ROUND_DOWN;
							out[0] = y[0] + y[2];
							ROUND_UP;
							out[1] = y[1] + y[3];
							break;
						case 5:
							out[0] = -y[11];
							out[1] = -y[10];
							break;
						default:
							out[0] = 0.0;
							out[1] = 0.0;
							break;
					}
					break;
				case 3:
/*			{-y[4], 
				0, 
				0, 
				-y[1], 
				0, 
				0} */
					switch( j ) {
						case 3:
							out[0] = -y[1];
							out[1] = -y[0];
							break;
						default:
							out[0] = 0.0;
							out[1] = 0.0;
							break;
					}
					break;
				case 4:
/*			{0, 
				-y[5], 
				0, 
				0, 
				-y[2], 
				0} */
					switch( j ) {
						case 4:
							out[0] = -y[3];
							out[1] = -y[2];
							break;
						default:
							out[0] = 0.0;
							out[1] = 0.0;
							break;
					}
					break;
				case 5:
/*			{0, 
				0, 
				-y[6], 
				0, 
				0, 
				-y[3]} */
					out[0] = -y[5];
					out[1] = -y[4];
					break;
			}
			a_ij[i][n] = out[0];
			a_ij[i][n+1] = out[1];
		}
	}
}


/* d_(y_i) = 2*d_(x_i)*y_i */

void sp_sold_i( double y[12], double d_xi[12], 
	double out[12] )
{
	double temp[2];
	int i;
	
	for( i=0; i<12; i+=2 ) {
		i_mult( d_xi + i, y + i, temp );
		out[i] = 2.0*temp[0];
		out[i+1] = 2.0*temp[1];
	}
}


/* d_(y_i y_j) = 2*( 2*d_xij*y_i*y_j + d_xi*delta(i,j) )  */
void sp_sold_ij( double y[12], double d_xi[12], 
	double d_xij[6][12], double out[6][12] )
{
	double t1[2], t2[2];
	int i, j, m, n, mp, np;
	
	for( i=0; i<6; i++ ) {
		m = 2*i;
		mp = m + 1;
		for( j=0; j<6; j++ ) {
			n = 2*j;
			np = n + 1;
			ROUND_DOWN;
			t1[0] = 4.0*y[m]*y[n];
			ROUND_UP;
			t1[1] = 4.0*y[mp]*y[np];
			i_mult( d_xij[i] + n, t1, t2 );
			if( i==j ) {
				ROUND_DOWN;
				t1[0] = t2[0] + 2.0*d_xi[m];
				ROUND_UP;
				t1[1] = t2[1] + 2.0*d_xi[mp];
				out[i][n] = t1[0];
				out[i][np] = t1[1];
			} else {
				out[i][n] = t2[0];
				out[i][np] = t2[1];
			}
		}
	}
}


/* sol_ij = ( c_ij*(1 + c^2/4) - 
						c*c_i*c_j/2)/(1 + c^2/4)^2	*/
void sp_sol_ij( double c[2], double c_i[12],
	double c_ij[6][12], double sol_ij[6][12] )
{
	double opc2[2], oopc2s[2], ohc[2];
	double t1[2], t2[2], t3[2], temp, temp2;
	int i, m, n;
	
	i_mult( c, c, t1 );
	ROUND_DOWN;
	temp = 1.0 + 0.25*t1[0];
	opc2[0] = temp;
	temp2 = temp*temp;
	ROUND_UP;
	oopc2s[1] = 1.0/temp2;
	temp = 1.0 + 0.25*t1[1];
	opc2[1] = temp;
	temp2 = temp*temp;
	ROUND_DOWN;
	oopc2s[0] = 1.0/temp2;
	
	ohc[0] = 0.5*c[0];
	ohc[1] = 0.5*c[1];
	
	for( i=0; i<6; i++ ) {
		m = 2*i;
		for( n=0; n<12; n+=2 ) {
			i_mult( c_ij[i] + n, opc2, t1 );
			i_mult( ohc, c_i + n, t2 );
			i_mult( t2, c_i + m, t3 );
			I_SUB( t1, t3, t2 );
			i_mult( t2, oopc2s, sol_ij[i] + n );
		}
	}
}


/* Another main routine.  Modified to give sol_xij,
instead of sol_yij.  */
void sp_solsecparbds( double y[12], double x[12],
	double out[6][12] )
{
	double a[2], c[2], d[2], a_i[12], c_i[12], d_i[12];
	double a_yi[12], a_yij[6][12];
	double a_ij[6][12], c_ij[6][12], d_ij[6][12];
	double delta[2], delta_part[12], delta_ij[6][12];
	double oo2y[12];
	int i;
#if DEBUG
	int j;
#endif

	ROUND_DOWN;
	for( i=0; i<12; i+=2 ) {
		oo2y[i] = 0.5/y[i+1];
	}
	ROUND_UP;
	for( i=1; i<12; i+=2 ) {
		oo2y[i] = 0.5/y[i-1];
	}

	sp_deltasecondpars( x, delta_ij );
	sp_deltapars( x, delta_ij, delta_part );
	i_afunc( y, a );
	a_partials( y, a_yi );
	sp_ytoxpars( oo2y, a_yi, a_i );
	sp_dfunction( x, delta_part, delta, d );
	sp_dpars( delta_part, d, d_i );
	sp_solasecpars( y, a_yij );
	sp_ytoxsecpars( oo2y, a_i, a_yij, a_ij );
	sp_dsecondpars( delta, d, delta_part, delta_ij, d_ij );
	i_div( d, a, c );
	sp_ubv_i( d, a, d_i, a_i, c_i );
	sp_ubv_ij( d, a, d_i, a_i, d_ij, a_ij, c_ij );
	sp_sol_ij( c, c_i, c_ij, out );
	
/*  This is the old routine.
	sp_deltasecondpars( x, delta_ij );
	sp_deltapars( x, delta_ij, delta_part );
	i_afunc( y, a );
	a_partials( y, a_i );
	sp_dfunction( x, delta_part, delta, d );
	sp_dpars( delta_part, d, d_xi );
	sp_sold_i( y, d_xi, d_i );
	sp_solasecpars( y, a_ij );
	sp_dsecondpars( delta, d, delta_part, delta_ij, d_xij );
	sp_sold_ij( y, d_xi, d_xij, d_ij );
	i_div( d, a, c );
	sp_ubv_i( d, a, d_i, a_i, c_i );
	sp_ubv_ij( d, a, d_i, a_i, d_ij, a_ij, c_ij );
	sp_sol_ij( c, c_i, c_ij, out );
*/

#if DEBUG
	ROUND_NEAR;
	printf("a = [%.18g,\t%.18g]\n", a[0], a[1]);
	printf("c = [%.18g,\t%.18g]\n", c[0], c[1]);
	printf("d = [%.18g,\t%.18g]\n", d[0], d[1]);
	printf("a_i = \n");
	for( i=0; i<12; i+=2 ) {
		printf("[%.18g,\t%.18g]\n", a_i[i], a_i[i+1]);
	}
	printf("\n");
	printf("c_i = \n");
	for( i=0; i<12; i+=2 ) {
		printf("[%.18g,\t%.18g]\n", c_i[i], c_i[i+1]);
	}
	printf("\n");
	printf("d_i = \n");
	for( i=0; i<12; i+=2 ) {
		printf("[%.18g,\t%.18g]\n", d_i[i], d_i[i+1]);
	}

	printf("a_ij = \n");
	for( i=0; i<6; i++ ) {
		for( j=2*i; j<12; j+=2 ) {
			printf("[%.18g,\t%.18g]\n", 
				a_ij[i][j], a_ij[i][j+1]);
		}
		printf("\n");
	}
	printf("c_ij = \n");
	for( i=0; i<6; i++ ) {
		for( j=2*i; j<12; j+=2 ) {
			printf("[%.18g,\t%.18g]\n", 
				c_ij[i][j], c_ij[i][j+1]);
		}
		printf("\n");
	}
	printf("d_ij = \n");
	for( i=0; i<6; i++ ) {
		for( j=2*i; j<12; j+=2 ) {
			printf("[%.18g,\t%.18g]\n", 
				d_ij[i][j], d_ij[i][j+1]);
		}
		printf("\n");
	}
#endif
}


void sp_solsecyparbds( double y[12], double x[12],
	double out[6][12] )
{
	double a[2], c[2], d[2], a_i[12], c_i[12], d_xi[12];
	double a_ij[6][12], c_ij[6][12], d_ij[6][12];
	double delta[2], delta_part[12], delta_ij[6][12];
	double d_i[12], d_xij[6][12];
	
	sp_deltasecondpars( x, delta_ij );
	sp_deltapars( x, delta_ij, delta_part );
	i_afunc( y, a );
	a_partials( y, a_i );
	sp_dfunction( x, delta_part, delta, d );
	sp_dpars( delta_part, d, d_xi );
	sp_sold_i( y, d_xi, d_i );
	sp_solasecpars( y, a_ij );
	sp_dsecondpars( delta, d, delta_part, delta_ij, d_xij );
	sp_sold_ij( y, d_xi, d_xij, d_ij );
	i_div( d, a, c );
	sp_ubv_i( d, a, d_i, a_i, c_i );
	sp_ubv_ij( d, a, d_i, a_i, d_ij, a_ij, c_ij );
	sp_sol_ij( c, c_i, c_ij, out );
}


/*
dihafun = 
-x[1]^2 + x[1]*x[2] + x[1]*x[3] - x[2]*x[3] - 
  2*x[1]*x[4] + x[1]*x[5] + x[2]*x[5] + x[1]*x[6] + 
  x[3]*x[6] - x[5]*x[6]
*/
/* Caveat:  no long ( > 2Sqrt[2] ) edges.  */
/* min: delta4( x1, x2_, x3_, x4^, x5_, x6_ ) */
/* max: delta4( x1, x2^, x3^, x4_, x5^, x6^ ) */

void sp_dih_a( double x[12], double out[2] )
{
	double max1, max2, temp, p1, p2, p3, val, nterms, pterms;
	
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
	out[0] = val;
	
	max1 = 0.5*(x[3] + x[5] - 2.0*x[6] + x[9] + x[11]);
	/* max1 is lower bound on potential maximum */
	ROUND_UP;
	temp = x[3] + x[5] - 2.0*x[6] + x[9] + x[11];
	max2 = 0.5*temp;
	/* max2 is upper bound on potential maximum */
	if( max1 > x[1] || max2 < x[0] ) {
		p1 = x[0]*(-x[0] + temp);
		p2 = x[1]*(-x[1] + temp);
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
	out[1] = val;
}


/*
dihcfun[x__]:=tomsu[{x[[1]],x[[3]],x[[5]]}]
*/
void sp_dih_c( double x[12], double cpars[12],
	double out[2] )
{
	double pterms[2], nterms[2], xp[12], xpp[12];
	int i, j;
	
	for( i=0; i<2; i++ ) {
		xp[i] = x[i];
		xpp[i] = x[i];
	}
	for( i=4; i<6; i++ ) {
		xp[i] = x[i];
		xpp[i] = x[i];
	}
	for( i=8; i<10; i++ ) {
		xp[i] = x[i];
		xpp[i] = x[i];
	}
	
	for( i=0; i<2; i+=2 ) {
		j = i + 1;
		if( cpars[i] >= 0.0 ) {
			xp[j] = x[i];		/* for min */
			xpp[i] = x[j];	/* for max */
		} 
		else if( cpars[j] <= 0.0 ) {
			xp[i] = x[j];		/* for min */
			xpp[j] = x[i];	/* for max */
		}
	}
	for( i=4; i<6; i+=2 ) {
		j = i + 1;
		if( cpars[i] >= 0.0 ) {
			xp[j] = x[i];		/* for min */
			xpp[i] = x[j];	/* for max */
		} 
		else if( cpars[j] <= 0.0 ) {
			xp[i] = x[j];		/* for min */
			xpp[j] = x[i];	/* for max */
		}
	}
	for( i=8; i<10; i+=2 ) {
		j = i + 1;
		if( cpars[i] >= 0.0 ) {
			xp[j] = x[i];		/* for min */
			xpp[i] = x[j];	/* for max */
		} 
		else if( cpars[j] <= 0.0 ) {
			xp[i] = x[j];		/* for min */
			xpp[j] = x[i];	/* for max */
		}
	}

	ROUND_DOWN;
	pterms[0] = 2.0*(xp[0]*(xp[8] + xp[4]) + xp[4]*xp[8]);
	nterms[0] = xpp[0]*xpp[0] + xpp[4]*xpp[4] + xpp[8]*xpp[8];
	ROUND_UP;
	pterms[1] = 2.0*(xpp[1]*(xpp[9] + xpp[5]) + xpp[5]*xpp[9]);
	nterms[1] = xp[1]*xp[1] + xp[5]*xp[5] + xp[9]*xp[9];
	out[1] = pterms[1] - nterms[0];
	ROUND_DOWN;
	out[0] = pterms[0] - nterms[1];
	if( out[0] < 0.0 )
		out[0] = 0.0;
}


/*
dihapars = 
{-2*x[1] + x[2] + x[3] - 2*x[4] + x[5] + x[6], 
  x[1] - x[3] + x[5], x[1] - x[2] + x[6], -2*x[1], 
  x[1] + x[2] - x[6], x[1] + x[3] - x[5]}
*/
void sp_dih_apars( double x[12], double out[12] )
{
	double temp[2];
	int m, n;
	
	for( m=0; m<6; m++ ) {
		n = 2*m;
		switch( m ) {
			case 0:
				ROUND_DOWN;
				temp[0] = -2.0*x[1] + x[2] + x[4] - 2.0*x[7] + 
					x[8] + x[10];
				ROUND_UP;
				temp[1] = -2.0*x[0] + x[3] + x[5] - 2.0*x[6] + 
					x[9] + x[11];
				break;
			case 1:
				ROUND_DOWN;
				temp[0] = x[0] - x[5] + x[8];
				ROUND_UP;
				temp[1] = x[1] - x[4] + x[9];
				break;
			case 2:
				ROUND_DOWN;
				temp[0] = x[0] - x[3] + x[10];
				ROUND_UP;
				temp[1] = x[1] - x[2] + x[11];
				break;
			case 3:
				ROUND_DOWN;
				temp[0] = -2.0*x[1];
				ROUND_UP;
				temp[1] = -2.0*x[0];
				break;
			case 4:
				ROUND_DOWN;
				temp[0] = x[0] + x[2] - x[11];
				ROUND_UP;
				temp[1] = x[1] + x[3] - x[10];
				break;
			case 5:
				ROUND_DOWN;
				temp[0] = x[0] + x[4] - x[9];
				ROUND_UP;
				temp[1] = x[1] + x[5] - x[8];
				break;
			default:
				break;
		}
		out[n] = temp[0];
		out[n+1] = temp[1];
	}
}


/*
dih_cpars = 
{-2*x[1] + 2*x[3] + 2*x[5], 0, 2*x[1] - 2*x[3] + 2*x[5], 
  0, 2*x[1] + 2*x[3] - 2*x[5], 0}
*/
void sp_dih_cpars( double x[12], double out[12] )
{
	double temp[2];
	int m, n;
	
	for( m=0; m<6; m++ ) {
		n = 2*m;
		switch( m ) {
			case 0:
				ROUND_DOWN;
				temp[0] = -x[1] + x[4] + x[8];
				ROUND_UP;
				temp[1] = -x[0] + x[5] + x[9];
				out[n] = 2.0*temp[0];
				out[n+1] = 2.0*temp[1];
				break;
			case 2:
				ROUND_DOWN;
				temp[0] = x[0] - x[5] + x[8];
				ROUND_UP;
				temp[1] = x[1] - x[4] + x[9];
				out[n] = 2.0*temp[0];
				out[n+1] = 2.0*temp[1];
				break;
			case 4:
				ROUND_DOWN;
				temp[0] = x[0] + x[4] - x[9];
				ROUND_UP;
				temp[1] = x[1] + x[5] - x[8];
				out[n] = 2.0*temp[0];
				out[n+1] = 2.0*temp[1];
				break;
			default:
				out[n] = 0.0;
				out[n+1] = 0.0;
				break;
		}
	}
}


/*
dihasecpars = 
	{{-2, 1, 1, -2, 1, 1}, 
	{1, 0, -1, 0, 1, 0}, 
	{1, -1, 0, 0, 0, 1}, 
	{-2, 0, 0, 0, 0, 0}, 
	{1, 1, 0, 0, 0, -1}, 
	{1, 0, 1, 0, -1, 0}}
*/
void sp_dihasecpars( double a_ij[6][12] )
{
	int i, j, m, n;
	double out[2];
	
	for( i=0; i<6; i++ ) {
		m = 2*i;
		for( j=0; j<6; j++ ) {
			n = 2*j;
			switch( i ) {
				case 0:
/*			{-2, 1, 1, -2, 1, 1} 	*/
					if( j==0 || j==3 ) {
						out[0] = -2.0;
						out[1] = -2.0;
					} else {
						out[0] = 1.0;
						out[1] = 1.0;
					}
					break;
				case 1:
/*			{1, 0, -1, 0, 1, 0} 	*/
					if( j==2 ) {
						out[0] = -1.0;
						out[1] = -1.0;
					} else if( j==4 ) {
						out[0] = 1.0;
						out[1] = 1.0;
					} else {
						out[0] = 0.0;
						out[1] = 0.0;
					}
					break;
				case 2:
/*			{1, -1, 0, 0, 0, 1} 	*/
					if( j==5 ) {
						out[0] = 1.0;
						out[1] = 1.0;
					} else {
						out[0] = 0.0;
						out[1] = 0.0;
					}
					break;
				case 3:
/*			{-2, 0, 0, 0, 0, 0} 	*/
					out[0] = 0.0;
					out[1] = 0.0;
					break;
				case 4:
/*			{1, 1, 0, 0, 0, -1} 	*/
					if( j==5 ) {
						out[0] = -1.0;
						out[1] = -1.0;
					} else {
						out[0] = 0.0;
						out[0] = 0.0;
					}
					break;
				case 5:
/*			{1, 0, 1, 0, -1, 0} 	*/
					out[0] = 0.0;
					out[1] = 0.0;
					break;
			}
			a_ij[i][n] = out[0];
			a_ij[i][n+1] = out[1];
		}
	}
}


/*
dihcsecpars = 
	{{-2, 0, 2, 0, 2, 0}, 
	{0, 0, 0, 0, 0, 0}, 
	{2, 0, -2, 0, 2, 0}, 
	{0, 0, 0, 0, 0, 0}, 
	{2, 0, 2, 0, -2, 0}, 
	{0, 0, 0, 0, 0, 0}}
*/	/* These all could be optimized a bit more. */
void sp_dihcsecpars( double c_ij[6][12] )
{
	int i, j, m, n;
	double out[2];
	
	for( i=0; i<6; i++ ) {
		m = 2*i;
		for( j=0; j<6; j++ ) {
			n = 2*j;
			switch( i ) {
				case 0:
/*			{-2, 0, 2, 0, 2, 0} 	*/
					if( j==0 ) {
						out[0] = -2.0;
						out[1] = -2.0;
					} else if( j==2 || j==4 ) {
						out[0] = 2.0;
						out[1] = 2.0;
					} else {
						out[0] = 0.0;
						out[1] = 0.0;
					}
					break;
				case 2:
/*			{2, 0, -2, 0, 2, 0} 	*/
					if( j==2 ) {
						out[0] = -2.0;
						out[1] = -2.0;
					} else if( j==4 ) {
						out[0] = 2.0;
						out[1] = 2.0;
					} else {
						out[0] = 0.0;
						out[1] = 0.0;
					}
					break;
				case 4:
/*			{2, 0, 2, 0, -2, 0} 	*/
					if( j==4 ) {
						out[0] = -2.0;
						out[1] = -2.0;
					} else {
						out[0] = 0.0;
						out[1] = 0.0;
					}
					break;
				default:
/*			{0, 0, 0, 0, 0, 0} 		*/
					out[0] = 0.0;
					out[1] = 0.0;
					break;
			}
			c_ij[i][n] = out[0];
			c_ij[i][n+1] = out[1];
		}
	}
}


/* Another main routine.  */
void sp_dihsecparbds( double x[12], double out[6][12] )
{
	double a[2], b[2], c[2], d[2], v[2];
	double a_i[12], b_i[12], c_i[12], d_i[12];
	double v_i[12], ubv_i[12];
	double a_ij[6][12], b_ij[6][12], c_ij[6][12];
	double ubv_ij[6][12], v_ij[6][12];
	double t1[2], t2[2], t3[2], temp_i[12];
	double  p[2], q[2], p_i[12], q_i[12], p_ij[6][12];
	int i, j, m;

	/* a = dih_a, b = dih_c, c = sp_cfun */
	sp_dihasecpars( a_ij );
	sp_dihcsecpars( b_ij );
	sp_csecondpars( c_ij );
	sp_dih_apars( x, a_i );
	sp_cpars( x, c_i );
	sp_dih_cpars( x, b_i );
	/* sp_dih_a() is an improved version of i_delta4(),
		but is restricted to cases where all edge lengths
		are less than 2Sqrt[2] */
#if NOLONG
	sp_dih_a( x, a );
#else
	i_delta4( x, a );
#endif
	sp_dih_c( x, b_i, b );
	sp_cfunction( x, c_i, c );
	
	I_MULT( b, c, p );
	I_SQRT( p, v );	/* v = sqrt( b*c ) */
	
	/* d = -v/sqrt( v^2 - a^2 ) */
	/* I_MULT( v, v, t1 ); */
	I_MULT( a, a, t2 );
	I_SUB( p, t2, q );
	I_SQRT( q, t1 );
	i_div( v, t1, t3 );
	d[0] = -t3[1];
	d[1] = -t3[0];
	
	/* d_i = 1/(2*d)*(b*c/(b*c - a^2))_i */
	/* p = b*c, q = b*c - a^2 */
	sp_u_i( b, c, b_i, c_i, p_i );
	t1[0] = -2.0*a[1];
	t1[1] = -2.0*a[0];
	for( i=0; i<12; i+=2 ) {
		i_mult( t1, a_i + i, temp_i + i );
	}
	ROUND_DOWN;
	for( i=0; i<12; i+=2 ) {
		q_i[i] = temp_i[i] + p_i[i];
	}
	ROUND_UP;
	for( i=1; i<12; i+=2 ) {
		q_i[i] = temp_i[i] + p_i[i];
	}
	sp_ubv_i( p, q, p_i, q_i, temp_i );
	t1[0] = 0.5;
	t1[1] = 0.5;
	i_div( t1, d, t2 );
	for( i=0; i<12; i+=2 ) {
		i_mult( t2, temp_i + i, d_i + i );
	}
	
	/* q = 1/(2*sqrt(p)) */
	sp_u_i( b, c, b_i, c_i, temp_i );
	t1[0] = 0.5;
	t1[1] = 0.5;
	i_div( t1, v, q );
	/* v_i = 1/(2*sqrt(b*c))*(bc)_i */
	/* v_i = q*p_i */
	for( i=0; i<12; i+=2 ) {
		i_mult( temp_i + i, q, v_i + i );
	}
	/* q_i = -2*(q^3)*p_i */
	I_MULT( q, q, t1 );
	I_MULT( q, t1, t2 );
	t1[0] = -2.0*t2[1];
	t1[1] = -2.0*t2[0];
	for( i=0; i<12; i+=2 ) {
		i_mult( t1, p_i + i, q_i + i );
	}
	/* v_ij = (q*p_i)_j */
	sp_u_ij( b, c, b_i, c_i, b_ij, c_ij, p_ij );
	for( i=0; i<6; i++ ) {
		m = 2*i;
		for( j=m; j<12; j+=2 ) {
			i_mult( q_i + j, p_i + m, t1 );
			i_mult( q, p_ij[i] + j, t3 );
			i_add( t1, t3, v_ij[i] + j );
		}
	}
	sp_ubv_ij( a, v, a_i, v_i, a_ij, v_ij, ubv_ij );
	sp_ubv_i( a, v, a_i, v_i, ubv_i );
	/* dih_ij = d_j*(a/v)_i + d*(a/v)_ij */
	for( i=0; i<6; i++ ) {
		m = 2*i;
		for( j=m; j<12; j+=2 ) {
			i_mult( d_i + j, ubv_i + m, t1 );
			i_mult( d, ubv_ij[i] + j, t2 );
			i_add( t1, t2, out[i] + j );
		}
	}
	
	/* Start here. */
	
#if DEBUG
	ROUND_NEAR;
	printf("a = [%.18g,\t%.18g]\n", a[0], a[1]);
	printf("b = [%.18g,\t%.18g]\n", b[0], b[1]);
	printf("c = [%.18g,\t%.18g]\n", c[0], c[1]);
	printf("d = [%.18g,\t%.18g]\n", d[0], d[1]);
	printf("a_i = \n");
	for( i=0; i<12; i+=2 ) {
		printf("[%.18g,\t%.18g]\n", a_i[i], a_i[i+1]);
	}
	printf("\n");
	printf("b_i = \n");
	for( i=0; i<12; i+=2 ) {
		printf("[%.18g,\t%.18g]\n", b_i[i], b_i[i+1]);
	}
	printf("\n");
	printf("c_i = \n");
	for( i=0; i<12; i+=2 ) {
		printf("[%.18g,\t%.18g]\n", c_i[i], c_i[i+1]);
	}
	printf("\n");
	printf("d_i = \n");
	for( i=0; i<12; i+=2 ) {
		printf("[%.18g,\t%.18g]\n", d_i[i], d_i[i+1]);
	}
	/* Compute dih_i */
	for( i=0; i<12; i+=2 ) {
		i_mult( d, ubv_i + i, d_i + i );
	}
	ROUND_NEAR;
	printf("\n");
	printf("dih_i = \n");
	for( i=0; i<12; i+=2 ) {
		printf("[%.18g,\t%.18g]\n", d_i[i], d_i[i+1]);
	}

	printf("a_ij = \n");
	for( i=0; i<6; i++ ) {
		for( j=2*i; j<12; j+=2 ) {
			printf("[%.18g,\t%.18g]\n", 
				a_ij[i][j], a_ij[i][j+1]);
		}
		printf("\n");
	}
	printf("b_ij = \n");
	for( i=0; i<6; i++ ) {
		for( j=2*i; j<12; j+=2 ) {
			printf("[%.18g,\t%.18g]\n", 
				b_ij[i][j], b_ij[i][j+1]);
		}
		printf("\n");
	}
	printf("c_ij = \n");
	for( i=0; i<6; i++ ) {
		for( j=2*i; j<12; j+=2 ) {
			printf("[%.18g,\t%.18g]\n", 
				c_ij[i][j], c_ij[i][j+1]);
		}
		printf("\n");
	}
#endif
}


/* Convert ypars to xpars.  Here oo2y is 1/(2*y).  */
void sp_ytoxpars( double oo2y[12], double ypars[12],
	double xpars[12] )
{
	int i;
	
	for( i=0; i<12; i+=2 ) {
		i_mult( oo2y + i, ypars + i, xpars + i );
	}
}


/* Convert u_(y_i y_j) to u_(x_i x_j).  */
void sp_ytoxsecpars( double oo2y[12], double xpars[12],
	double ysecpars[6][12], double xsecpars[6][12] )
{
	double t1[2], t2[2];
	int i, j, m;
	
	for( i=0; i<6; i++ ) {
		m = 2*i;
		for( j=m; j<12; j+=2 ) {
			if( j==m ) {
				t2[0] = 2.0*xpars[j];
				t2[1] = 2.0*xpars[j+1];
				i_sub( ysecpars[i] + j, t2, t1 );
			} else {
				t1[0] = ysecpars[i][j];
				t1[1] = ysecpars[i][j+1];
			}
			i_mult( oo2y + j, t1, t2 );
			i_mult( oo2y + m, t2, xsecpars[i] + j );
		}
	}
}


/* Convert xpars to ypars.  */
void sp_xtoypars( double y[12], double xpars[12],
	double ypars[12] )
{
	int i;
	double pars[12];
	
	for( i=0; i<12; i+=2 ) {
		i_mult( y + i, xpars + i, pars + i );
	}
	for( i=0; i<12; i++ ) {
		ypars[i] = 2.0*pars[i];
	}
}


/* Convert u_(x_i x_j) to u_(y_i y_j).  */
void sp_xtoysecpars( double y[12], double xpars[12],
	double xsecpars[6][12], double ysecpars[6][12] )
{
	int i, j, k;
	double temp[2];
	
	for( i=0; i<6; i++ ) {
		k = 2*i;
		for( j=k; j<12; j+=2 ) {
			ROUND_DOWN;
			temp[0] = 4.0*y[j]*y[k];
			ROUND_UP;
			temp[1] = 4.0*y[j+1]*y[k+1];
			i_mult( temp, xsecpars[i] + j, ysecpars[i] + j );
		}
	}
	for( i=0; i<6; i++ ) {
		j = 2*i;
		ROUND_DOWN;
		ysecpars[i][j] += 2.0*xpars[j];
		ROUND_UP;
		ysecpars[i][j+1] += 2.0*xpars[j+1];
	}
}


/* Here is the master routine.  Hopefully.  */
void sp_gmavolsecparbds( double x[12], double out[6][12] )
{
	double d[2], d_i[12], delta[2];
	double delta_i[12], delta_ij[6][12];
	
#if DEBUG || DEBUG2
	int i, j;
#endif

	sp_deltasecondpars( x, delta_ij );
	sp_deltapars( x, delta_ij, delta_i );
	sp_dfunction( x, delta_i, delta, d );
	sp_dpars( delta_i, d, d_i );
	sp_dsecondpars( delta, d, delta_i, delta_ij, out );
	/* This omits a factor of 1/12 from the definition.  */

#if DEBUG
	ROUND_NEAR;
	printf("d_ij = \n");
	for( i=0; i<6; i++ ) {
		for( j=2*i; j<12; j+=2 ) {
			printf("[%.18g,\t%.18g]\n", 
				out[i][j], out[i][j+1]);
		}
		printf("\n");
	}
#endif

#if DEBUG2
	ROUND_NEAR;
	printf("d = [%.18g,\t%.18g]\n", d[0], d[1]);
	printf("\n");
	printf("d_i = \n");
	for( i=0; i<12; i+=2 ) {
		printf("[%.18g,\t%.18g]\n", d_i[i], d_i[i+1]);
	}
	printf("\n");
	printf("delta_i = \n");
	for( i=0; i<12; i+=2 ) {
		printf("[%.18g,\t%.18g]\n", delta_i[i], delta_i[i+1]);
	}
	printf("d_ij = \n");
	for( i=0; i<6; i++ ) {
		for( j=2*i; j<12; j+=2 ) {
			printf("[%.18g,\t%.18g]\n", 
				d_ij[i][j], d_ij[i][j+1]);
		}
		printf("\n");
	}
	printf("delta_ij = \n");
	for( i=0; i<6; i++ ) {
		for( j=2*i; j<12; j+=2 ) {
			printf("[%.18g,\t%.18g]\n", 
				delta_ij[i][j], delta_ij[i][j+1]);
		}
		printf("\n");
	}
#endif
}


/* b = 4*x1*delta */
void dih_bfun( double x[12], double delta[2], double out[2] )
{
	ROUND_DOWN;
	out[0] = 4.0*x[0]*delta[0];
	ROUND_UP;
	out[1] = 4.0*x[1]*delta[1];
}


/* b_i = 4*x1*delta_i + 4*delta*D(i,1) */
void dih_bfun_i( double x[12], double delta[2],
	double delta_i[12], double out[12] )
{
	double temp[2];
	int i;
	
	temp[0] = 4.0*x[0];
	temp[1] = 4.0*x[1];
	for( i=0; i<12; i+=2 ) {
		i_mult( temp, delta_i + i, out + i );
	}
	ROUND_DOWN;
	out[0] += 4.0*delta[0];
	ROUND_UP;
	out[1] += 4.0*delta[1];
}


/* b_ij = 4*x1*delta_ij + 4*delta_i*D(j,1) 
	+ 4*delta_j*D(i,1) */
void dih_bfun_ij( double x[12], double delta_i[12],
	double delta_ij[6][12], double out[6][12] )
{
	double temp[2];
	int i, j;
	
	temp[0] = 4.0*x[0];
	temp[1] = 4.0*x[1];
	for( i=0; i<6; i++ ) {
		for( j=2*i; j<12; j+=2 ) {
			i_mult( temp, delta_ij[i] + j, out[i] + j );
		}
	}
	ROUND_DOWN;
	for( j=0; j<12; j+=2 ) {
		out[0][j] += 4.0*delta_i[j];
	}
	out[0][0] += 4.0*delta_i[0];
	ROUND_UP;
	for( j=1; j<12; j+=2 ) {
		out[0][j] += 4.0*delta_i[j];
	}
	out[0][1] += 4.0*delta_i[1];
}


/* a = delta4 */
void dih_afun( double delta_i[12], double out[2] )
{
	out[0] = delta_i[6];
	out[1] = delta_i[7];
}


void dih_afun_i( double delta_ij[6][12], double out[12] )
{
	int i, j;
	
	for( i=0; i<3; i++ ) {
		j = 2*i;
		out[j  ] = delta_ij[i][6];
		out[j+1] = delta_ij[i][7];
	}
	for( j=6; j<12; j++ ) {
		out[j] = delta_ij[3][j];
	}
}


void dih_afun_ij( double out[6][12] )
{
	double secpars[6][6] = {{-2,1,1,-2,1,1},
												{1,0,-1,0,1,0},
												{1,-1,0,0,0,1},
												{-2,0,0,0,0,0},
												{1,1,0,0,0,-1},
												{1,0,1,0,-1,0}};
	double temp;
	int i, j, k, kp;
	
	for( i=0; i<6; i++ ) {
		for( j=i; j<6; j++ ) {
			k = 2*j;
			kp = k + 1;
			temp = secpars[i][j];
			out[i][k] = temp;
			out[i][kp] = temp;
		}
	}
}


/* c = 1/(a*a + b) */
void dih_cfun( double a[2], double b[2], double out[2] )
{
	double num[2], den[2];
	
	i_mult( a, a, num );
	I_ADD( num, b, den );
	num[0] = 1.0;
	num[1] = 1.0;
	i_div( num, den, out );
}


/* c_i = -(2*a*a_i+b_i)/(a*a + b)^2 */
void dih_cfun_i( double a[2], double b[2], double a_i[12],
	double b_i[12], double out[12] )
{
	double aapb[2], num[2], den[2], temp[2], twoa[2];
	double temp2[2];
	int i;
	
	twoa[0] = 2.0*a[0];
	twoa[1] = 2.0*a[1];
	i_mult( a, a, temp );
	I_ADD( temp, b, aapb );
	i_mult( aapb, aapb, temp2 );
	temp[0] = 1.0;
	temp[1] = 1.0;
	i_div( temp, temp2, den );
	for( i=0; i<12; i+=2 ) {
		i_mult( twoa, a_i + i, temp );
		i_add( temp, b_i + i, temp2 );
		num[0] = -temp2[1];
		num[1] = -temp2[0];
		i_mult( num, den, out + i );
	}
}


/* b = delta */
void sol_bfun( double delta[2], double out[2] )
{
	out[0] = delta[0];
	out[1] = delta[1];
}


void sol_bfun_i( double y[12], double delta_i[12], 
	double out[12] )
{
	sp_xtoypars( y, delta_i, out );
}


void sol_bfun_ij( double y[12], double delta_i[12],
	double delta_ij[6][12], double out[6][12] )
{
	sp_xtoysecpars( y, delta_i, delta_ij, out );
}


/* a = a */
void sol_afun( double y[12], double out[2] )
{
	i_afunc( y, out );
}


void sol_afun_i( double y[12], double x[12], double part[12] )
{
	double pterms[2], nterms[2], temp[2];

	temp[1] = y[1]*y[3] + y[1]*y[5] + y[3]*y[5];
	ROUND_DOWN;
	temp[0] = y[0]*y[2] + y[0]*y[4] + y[2]*y[4];
	/* a_y1 = y2*y3 + (y2^2 + y3^2 - y4^2)/2 + y1*y2 + y1*y3 */
	pterms[1] = temp[1] + 0.5*(x[3] + x[5]);
	nterms[1] = 0.5*x[7];
	ROUND_DOWN;
	pterms[0] = temp[0] + 0.5*(x[2] + x[4]);
	nterms[0] = 0.5*x[6];
	I_SUB( pterms, nterms, part );
	/*	ROUND_UP;	*/
	/* a_y2 = y1*y3 + y1*y2 + (y1^2 + y3^2 - y5^2)/2 + y2*y3 */
	pterms[1] = temp[1] + 0.5*(x[1] + x[5]);
	nterms[1] = 0.5*x[9];
	ROUND_DOWN;
	pterms[0] = temp[0] + 0.5*(x[0] + x[4]);
	nterms[0] = 0.5*x[8];
	i_sub( pterms, nterms, part + 2 );
	/*	ROUND_UP;	*/
	/* a_y3 = y1*y2 + y1*y3 + y2*y3 + (y1^2 + y2^2 - y6^2)/2 */
	pterms[1] = temp[1] + 0.5*(x[1] + x[3]);
	nterms[1] = 0.5*x[11];
	ROUND_DOWN;
	pterms[0] = temp[0] + 0.5*(x[0] + x[2]);
	nterms[0] = 0.5*x[10];
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


/* {{y[2]+y[3],y[1]+y[2]+y[3],y[1]+y[2]+y[3],-y[4],0,0},{y[1]+y[2]+y[3],
    y[1]+y[3],y[1]+y[2]+y[3],0,-y[5],0},{y[1]+y[2]+y[3],y[1]+y[2]+y[3],
    y[1]+y[2],0,0,-y[6]},{-y[4],0,0,-y[1],0,0},{0,-y[5],0,0,-y[2],0},{0,
    0,-y[6],0,0,-y[3]}} */
void sol_afun_ij( double y[12], double out[6][12] )
{
	sp_solasecpars( y, out );
}


/* c = 1/(a*a + 0.25*delta) */
void sol_cfun( double a[2], double b[2], double out[2] )
{
	double num[2], den[2];
	
	i_mult( a, a, num );
	ROUND_DOWN;
	den[0] = num[0] + 0.25*b[0];
	ROUND_UP;
	den[1] = num[1] + 0.25*b[1];
	num[0] = 1.0;
	num[1] = 1.0;
	i_div( num, den, out );
}


/* c_i = -(2*a*a_i+0.25*delta_i)/(a*a+0.25*delta)^2 */
void sol_cfun_i( double a[2], double b[2], double a_i[12],
	double b_i[12], double out[12] )
{
	double aapb[2], num[2], den[2], temp[2], twoa[2];
	double temp2[2];
	int i;
	
	twoa[0] = 2.0*a[0];
	twoa[1] = 2.0*a[1];
	i_mult( a, a, temp );
	ROUND_DOWN;
	aapb[0] = temp[0] + 0.25*b[0];
	ROUND_UP;
	aapb[1] = temp[1] + 0.25*b[1];
	i_mult( aapb, aapb, temp2 );
	temp[0] = 1.0;
	temp[1] = 1.0;
	i_div( temp, temp2, den );
	for( i=0; i<12; i+=2 ) {
		i_mult( twoa, a_i + i, temp );
		ROUND_DOWN;
		temp2[0] = temp[0] + 0.25*b_i[i];
		ROUND_UP;
		temp2[1] = temp[1] + 0.25*b_i[i+1];
		num[0] = -temp2[1];
		num[1] = -temp2[0];
		i_mult( num, den, out + i );
	}
}


/* b^(3/2)[c(b^(1/2)/a)_i] = 
	b*c*(0.5*a*b_i - b*a_i) */
void clear_pars( double a[2], double b[2], double c[2],
	double a_i[12], double b_i[12], double out[12] )
{
	int i;
	double prod[2], p1[2], p2[2], sum[2];
	
	i_mult( b, c, prod );
	for( i=0; i<12; i+=2 ) {
		i_mult( a, b_i + i, p1 );
		i_mult( b, a_i + i, p2 );
		ROUND_DOWN;
		sum[0] = 0.5*p1[0] - p2[1];
		ROUND_UP;
		sum[1] = 0.5*p1[1] - p2[0];
		i_mult( prod, sum, out + i );
	}
}


/* b^(3/2)[(c(b^(1/2)/a)_i))j] = 
	b*c_j*(0.5*a*b_i - b*a_i) +
	c*((0.5*a_j*b_i-b_j*a_i+0.5*a*b_ij-b*a_ij)*b -
		b_j*(0.5*a*b_i-b*a_i)) */
void clear_secpars( double a[2], double b[2], double c[2],
	double a_i[12], double b_i[12], double c_i[12],
	 double a_ij[6][12], double b_ij[6][12], 
	 double out[6][12] )
{
	int i, j, k;
	double p1[2], p2[2], p3[2], p4[2], p5[2];
	double clear_i[12];
	
	for( i=0; i<12; i+=2 ) {
		i_mult( a, b_i + i, p1 );
		i_mult( b, a_i + i, p2 );
		ROUND_DOWN;
		p4[0] = 0.5*p1[0] - p2[1];
		ROUND_UP;
		p4[1] = 0.5*p1[1] - p2[0];
		i_mult( b, p4, clear_i + i );
	}
	
	for( i=0; i<6; i++ ) {
		k = 2*i;
		for( j=k; j<12; j+=2 ) {
			i_mult( a_i + j, b_i + k, p1 );
			i_mult( b_i + j, a_i + k, p2 );
			i_mult( a, b_ij[i] + j, p3 );
			i_mult( b, a_ij[i] + j, p4 );
			ROUND_DOWN;
			p5[0] = 0.5*p1[0] - p2[1] + 0.5*p3[0] - p4[1];
			ROUND_UP;
			p5[1] = 0.5*p1[1] - p2[0] + 0.5*p3[1] - p4[0];
			i_mult( p5, b, p1 ); /* save p1 */
			i_mult( a, b_i + k, p2 );
			i_mult( b, a_i + k, p3 );
			ROUND_DOWN;
			p4[0] = 0.5*p2[0] - p3[1];
			ROUND_UP;
			p4[1] = 0.5*p2[1] - p3[0];
			i_mult( b_i + j, p4, p2 );
			ROUND_DOWN;
			p3[0] = p1[0] - 0.5*p2[1];
			ROUND_UP;
			p3[1] = p1[1] - 0.5*p2[0];	/* done with p1 */
			i_mult( c, p3, p4 );
			i_mult( c_i + j, clear_i + k, p5 );
			ROUND_DOWN;
			out[i][j] = p4[0] + p5[0];
			ROUND_UP;
			out[i][j+1] = p4[1] + p5[1];
		}
	}
}


void clear_dih_i( double y[12], double x[12], 
	double clear_i[12] )
{
	int i;
	double a[2], b[2], c[2], a_i[12], b_i[12];
	double delta[2], delta_i[12], delta_ij[6][12];
	double part[12], temp[2];
	
	i_bigdelta_partials_best( x, delta, delta_i );
	if( delta[0] < 0.0 )
		delta[0] = 0.0;
	sp_deltasecondpars( x, delta_ij );
	dih_bfun( x, delta, b );
	dih_bfun_i( x, delta, delta_i, b_i );
	dih_afun( delta_i, a );
	dih_afun_i( delta_ij, a_i );
	dih_cfun( a, b, c );
	clear_pars( a, b, c, a_i, b_i, part );
	
	ROUND_DOWN;
	delta[0] = 8.0*y[0]*x[0];
	ROUND_UP;
	delta[1] = 8.0*y[1]*x[1];
	temp[1] = 1.0/delta[0];
	ROUND_DOWN;
	temp[0] = 1.0/delta[1];
	for( i=0; i<12; i+=2 ) {
		i_mult( temp, part + i, clear_i + i );
	}
}


void clear_dih_ij( double y[12], double x[12], 
	double clear_i[12], double clear_ij[6][12] )
{
	int i, j;
	double a[2], b[2], c[2], a_i[12], b_i[12];
	double c_i[12], a_ij[6][12], b_ij[6][12];
	double delta[2], delta_i[12], delta_ij[6][12];
	double part[12], secpart[6][12], temp[2];
	
	i_bigdelta_partials_best( x, delta, delta_i );
	if( delta[0] < 0.0 )
		delta[0] = 0.0;
	sp_deltasecondpars( x, delta_ij );
	dih_bfun( x, delta, b );
	dih_bfun_i( x, delta, delta_i, b_i );
	dih_bfun_ij( x, delta_i, delta_ij, b_ij );
	dih_afun( delta_i, a );
	dih_afun_i( delta_ij, a_i );
	dih_afun_ij( a_ij );
	dih_cfun( a, b, c );
	dih_cfun_i( a, b, a_i, b_i, c_i );
	clear_pars( a, b, c, a_i, b_i, part );
	clear_secpars( a, b, c, a_i, b_i, c_i, a_ij,
		b_ij, secpart );
	
	ROUND_DOWN;
	delta[0] = 8.0*y[0]*x[0];
	ROUND_UP;
	delta[1] = 8.0*y[1]*x[1];
	temp[1] = 1.0/delta[0];
	ROUND_DOWN;
	temp[0] = 1.0/delta[1];
	for( i=0; i<12; i+=2 ) {
		i_mult( temp, part + i, clear_i + i );
	}
	for( i=0; i<6; i++ ) {
		for( j=2*i; j<12; j+=2 ) {
			i_mult( temp, secpart[i] + j, clear_ij[i] + j );
		}
	}
}


void clear_sol_i( double y[12], double x[12], 
	double clear_i[12] )
{
	double a[2], b[2], c[2], a_i[12], b_i[12];
	double delta[2], delta_i[12], delta_ij[6][12];
	
	i_bigdelta_partials_best( x, delta, delta_i );
	if( delta[0] < 0.0 )
		delta[0] = 0.0;
	sp_deltasecondpars( x, delta_ij );
	sol_bfun( delta, b );
	sol_bfun_i( y, delta_i, b_i );
	sol_afun( y, a );
	sol_afun_i( y, x, a_i );
	sol_cfun( a, b, c );
	clear_pars( a, b, c, a_i, b_i, clear_i );
}


void clear_sol_ij( double y[12], double x[12], 
	double clear_i[12], double clear_ij[6][12] )
{
	double a[2], b[2], c[2], a_i[12], b_i[12];
	double c_i[12], a_ij[6][12], b_ij[6][12];
	double delta[2], delta_i[12], delta_ij[6][12];
	
	i_bigdelta_partials_best( x, delta, delta_i );
	if( delta[0] < 0.0 )
		delta[0] = 0.0;
	sp_deltasecondpars( x, delta_ij );
	sol_bfun( delta, b );
	sol_bfun_i( y, delta_i, b_i );
	sol_bfun_ij( y, delta_i, delta_ij, b_ij );
	sol_afun( y, a );
	sol_afun_i( y, x, a_i );
	sol_afun_ij( y, a_ij );
	sol_cfun( a, b, c );
	sol_cfun_i( a, b, a_i, b_i, c_i );
	clear_pars( a, b, c, a_i, b_i, clear_i );
	clear_secpars( a, b, c, a_i, b_i, c_i, a_ij,
		b_ij, clear_ij );
}


