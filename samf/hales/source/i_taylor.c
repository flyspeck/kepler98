/* i_taylor.c  (c) 1997, Samuel Ferguson.  This file contains 
routines which implement sphere-packing functions using 
second-order Taylor expansions. */


#include "system_headers.h"
#include "i_sphere.h"
#include "interval.h"
#include "i_bounds.h"
#include "macros.h"
#include "i_taylor.h"


#define DEBUG	0

/* Bounds on the second partials: */

/* Dihedral angle second partials */
/* [-0.23841250069568981,	0.16915087490595893] */
#define DIH_VAL_LO		-0.2384125007
#define	DIH_VAL_HI		 0.169150875

/* Solid angle second partial bounds */
/* [-0.10400745566587083,   0.13847858049517325] */
#define SOL_VAL_LO		-0.1040074557
#define	SOL_VAL_HI		 0.1384785805

/* Voronoi volume second partial bounds */
/* [-0.23738923824968214,   0.19941810082761346] */
#define VOR_VAL_LO		-0.2373892383
#define	VOR_VAL_HI		 0.1994181009

/* Voronoi second partial bounds (computed) */
/* vor_c = [-0.713720996103926, 0.869176515630031] */
#define VOR_C_LO			-0.7137209962
#define	VOR_C_HI		 	 0.8691765157

/* Gamma volume second partial bounds */
/* [-0.13621002209281341,	0.10165389226429553] */
#define GMA_VAL_LO		-0.1362100221
#define	GMA_VAL_HI		 0.1016538923

/* Gamma second partial bounds (computed) */
/* gma_c = [-0.211959198389001, 0.282832314019063] */
#define GMA_C_LO			-0.2119591984
#define	GMA_C_HI		 	 0.2828323141

/*
vor_c = [-0.713720996103926, 0.869176515630031]
gma_c = [-0.211959198389001, 0.282832314019063]
												 12345678901234567890
*/

extern double i_doct_const[2];


/* Routines */

void t_dih( double xt[12], double h[12], double h2, 
	double out[2] )
{
	double val[2], val_i[12], t1[2], temp;
	int i;
	
	i_dih( xt, val );
	s_dih_xpars( xt, val_i );
	for( i=0; i<12; i+=2 ) {
		i_mult( h + i, val_i + i, t1 );
		ROUND_DOWN;
		val[0] += t1[0];
		ROUND_UP;
		val[1] += t1[1];
	}
	ROUND_DOWN;
	temp = h2*DIH_VAL_LO;
	out[0] = temp + val[0];
	ROUND_UP;
	temp = h2*DIH_VAL_HI;
	out[1] = temp + val[1];
}


void t_dih_base( double xt[12], double out[2], 
	double out_i[12] )
{	
	i_dih( xt, out );
	s_dih_xpars( xt, out_i );
}


void t_solid( double yt[12], double ooyt[12], 
	double delta[2], double sqrtdelta[2], 
	double delta_part[12], double h[12], 
	double h2, double out[2] )
{
	double val[2], val_i[12], t1[2], temp;
	int i;
	
	i_solid( yt, sqrtdelta, val );
	s_solid_xpars( yt, ooyt, delta, sqrtdelta, 
		delta_part, val_i );
	for( i=0; i<12; i+=2 ) {
		i_mult( h + i, val_i + i, t1 );
		ROUND_DOWN;
		val[0] += t1[0];
		ROUND_UP;
		val[1] += t1[1];
	}
	ROUND_DOWN;
	temp = h2*SOL_VAL_LO;
	out[0] = temp + val[0];
	ROUND_UP;
	temp = h2*SOL_VAL_HI;
	out[1] = temp + val[1];
}


void t_solid_base( double yt[12], double ooyt[12], 
	double delta[2], double sqrtdelta[2], double delta_part[12],
	double out[2], double out_i[12] )
{
	i_solid( yt, sqrtdelta, out );
	s_solid_xpars( yt, ooyt, delta, sqrtdelta, 
		delta_part, out_i );
}


void t_vorvol( double xt[12], double delta[2],
	double sqrtdelta[2], double delta_part[12], 
	double h[12], double h2, double out[2] )
{
	double val[2], val_i[12], t1[2], temp;
	int i;
	
	i_voronoivol( xt, sqrtdelta, val );
	vorvol_xpars( xt, delta, sqrtdelta, delta_part, val_i );
	for( i=0; i<12; i+=2 ) {
		i_mult( h + i, val_i + i, t1 );
#if DEBUG
		ROUND_NEAR;
		printf("t1 = [%.18g, %.18g]\n", t1[0], t1[1]);
#endif
		ROUND_DOWN;
		val[0] += t1[0];
		ROUND_UP;
		val[1] += t1[1];
	}
	ROUND_DOWN;
	temp = h2*VOR_VAL_LO;
	out[0] = temp + val[0];
	ROUND_UP;
	temp = h2*VOR_VAL_HI;
	out[1] = temp + val[1];
#if DEBUG
	ROUND_NEAR;
	printf("h2 = %.18g\n", h2);
	printf("val in vorvol =\n");
	printf("[%.18g, %.18g]\n", val[0], val[1]);
	printf("h2 term in vorvol =\n");
	printf("[%.18g, %.18g]\n", 
		h2*VOR_VAL_LO, h2*VOR_VAL_HI);
#endif
}


void t_vor( double xt[12], double yt[12], double ooyt[12], 
	double delta[2], double sqrtdelta[2], double delta_part[12], 
	double h[12], double h2, double out[2] )
{
	double val[2], c_val[2];
	double solid_i[12], val_i[12];
	double t1[2];
	int i;
	
	c_val[0] = VOR_C_LO;
	c_val[1] = VOR_C_HI;

	i_vor( yt, xt, sqrtdelta, val );
	s_solid_xpars( yt, ooyt, delta, sqrtdelta, delta_part, solid_i );
	vor_xpars( xt, delta, sqrtdelta, delta_part, 
		solid_i, val_i );
	
	for( i=0; i<12; i+=2 ) {
		i_mult( h + i, val_i + i, t1 );
		ROUND_DOWN;
		val[0] += t1[0];
		ROUND_UP;
		val[1] += t1[1];
	}
	
	i_smult( h2, c_val, t1 );
	I_ADD( val, t1, out );
}


void t_vor_base( double xt[12], double yt[12],
	double delta[2],double sqrtdelta[2], double delta_part[12], 
	double solid_i[12], double out[2], double out_i[12] )
{	
	i_vor( yt, xt, sqrtdelta, out );
	vor_xpars( xt, delta, sqrtdelta, delta_part, 
		solid_i, out_i );
}


void t_gma( double yt[12], double ooyt[12],
	double delta[2], double sqrtdelta[2], 
	double delta_part[12], double h[12], double h2, 
	double out[2] )
{
	double val[2], c_val[2];
	double solid_i[12], bvol_i[12], val_i[12];
	double t1[2], oosqrtdelta[2];
	int i;
	
	/* compute 1/sqrtdelta */
	ROUND_DOWN;
	oosqrtdelta[0] = 1.0/sqrtdelta[1];
	ROUND_UP;
	oosqrtdelta[1] = 1.0/sqrtdelta[0];
	
	c_val[0] = GMA_C_LO;
	c_val[1] = GMA_C_HI;
	
	i_gma( yt, sqrtdelta, val );
	s_solid_xpars( yt, ooyt, delta, sqrtdelta, 
		delta_part, solid_i );
	s_bvol_xpars( yt, ooyt, delta, sqrtdelta, 
		delta_part, solid_i, bvol_i );
	s_gma_xpars( oosqrtdelta, delta_part, bvol_i, val_i );
	
	for( i=0; i<12; i+=2 ) {
		i_mult( h + i, val_i + i, t1 );
		ROUND_DOWN;
		val[0] += t1[0];
		ROUND_UP;
		val[1] += t1[1];
	}
	
	i_smult( h2, c_val, t1 );
	I_ADD( val, t1, out );
}


void t_gma_base( double yt[12], double ooyt[12],
	double delta[2], double sqrtdelta[2], 
	double delta_part[12], double solid_i[12],
	double out[2], double out_i[12] )
{
	double oosqrtdelta[2], bvol_i[12];

	/* compute  1/sqrtdelta */
	ROUND_DOWN;
	oosqrtdelta[0] = 1.0/sqrtdelta[1];
	ROUND_UP;
	oosqrtdelta[1] = 1.0/sqrtdelta[0];
	
	i_gma( yt, sqrtdelta, out );
	s_bvol_xpars( yt, ooyt, delta, sqrtdelta, 
		delta_part, solid_i, bvol_i );
	s_gma_xpars( oosqrtdelta, delta_part, bvol_i, out_i );
}


/* Relation constants are assumed to be known
exactly, so they are not expressed as intervals. */

/* Order of relation constants:
	sol, gma, vor, octavor, dih1, dih2, dih3	
	 0		1		 2			3		 		4		  5		  6	*/

/* 	rel describes which values we want computed
independently.
		relconst gives the relation constants.  */

/* Dependencies:
void t_dih_base( double xt[12], double out[2], 
	double out_i[12] )
void t_solid_base( double yt[12], double ooyt[12], 
	double delta[2], double sqrtdelta[2], double delta_part[12],
	double out[2], double out_i[12] )
void t_vor_base( double xt[12], double yt[12],
	double delta[2],double sqrtdelta[2], double delta_part[12], 
	double solid_i[12], double out[2], double out_i[12] )
void t_gma_base( double yt[12], double ooyt[12],
	double delta[2], double sqrtdelta[2], 
	double delta_part[12], double solid_i[12],
	double out[2], double out_i[12] )
*/
void t_composite( int rel[7], double relconst[7], 
	double y[12], double x[12], double out[2], 
	double out_i[12], double outvals[14] )
{
	double yt[12], ooyt[12], xt[12], xtp[12];
	double delta[2], sqrtdelta[2], oosqrtdelta[2];
	double delta_part[12], h[12], h2, hsum[2];
	double val[2], val_i[12], c_val[2];
	double solid_val[2], solid_i[12];
	double bvol_i[12];
	double gma_val[2], gma_i[12];
	double vor_val[2], vor_i[12];
	double octa_val[2], octa_i[12];
	double dih_val[2], dih_i[12], dih_ip[12];
	double temp, t1[2], local[2];
	double localcval[2];
	int i, ip, j, k;
	
	/* Compute t for both x and y, and h for x */
	ROUND_DOWN;
	for( i=0; i<12; i+=2 ) {
		ip = i + 1;
		xt[i] = 0.5*(x[ip] + x[i]);
		yt[i] = 0.5*(y[ip] + y[i]);
	}
	ROUND_UP;
	for( i=0; i<12; i+=2 ) {
		ip = i + 1;
		xt[ip] = 0.5*(x[ip] + x[i]);
		h[ip] = 0.5*(x[ip] - x[i]);
		h[i] = -h[ip];
		yt[ip] = 0.5*(y[ip] + y[i]);
	}
	
	/* Need ooyt for solid, gma, and octavor */
	/* Compute ooyt = 1/yt and hsum */
	temp = 0.0;
	ROUND_UP;
	for( i=1; i<12; i+=2 ) {
		ooyt[i] = 1.0/yt[i-1];
		temp += h[i];
	}
	hsum[1] = temp;
	ROUND_DOWN;
	temp = 0.0;
	for( i=0; i<12; i+=2 ) {
		ooyt[i] = 1.0/yt[i+1];
		temp += h[i];
	}
	hsum[0] = temp;

	/* Compute h2 */
	ROUND_UP;
	h2 = 0.0;
	for( i=1; i<12; i+=2 ) {
		h2 += h[i]*h[i];
	}
	h2 *= 0.5;
	for( i=0; i<10; i+=2 ) {
		for( j=i+2; j<12; j+=2 ) {
			h2 += h[i]*h[j];
		}
	}
	
	/* Do preliminary computations */
	s_delta_partials( xt, delta_part );
	i_bigdelta( xt, delta );
	I_SQRT( delta, sqrtdelta );
	
	ROUND_DOWN;
	oosqrtdelta[0] = 1.0/sqrtdelta[1];
	ROUND_UP;
	oosqrtdelta[1] = 1.0/sqrtdelta[0];

	/* zero out val, val_i, and c_val */
	val[0] = 0.0;
	val[1] = 0.0;
	c_val[0] = 0.0;
	c_val[1] = 0.0;
	for( i=0; i<12; i++ ) {
		val_i[i] = 0.0;
	}
	/* zero out outvals */
	for( i=0; i<10; i++ ) {
		outvals[i] = 0.0;
	}
	
	/* solid (need this for most things) */
	i_solid( yt, sqrtdelta, solid_val );
	s_solid_xpars( yt, ooyt, delta, sqrtdelta, 
		delta_part, solid_i );

	if( rel[0] ) {  /* Compute local bounds . . . */
		local[0] = solid_val[0];
		local[1] = solid_val[1];
		for( i=0; i<12; i+=2 ) {
			i_mult( h + i, solid_i + i, t1 );
			ROUND_DOWN;
			local[0] += t1[0];
			ROUND_UP;
			local[1] += t1[1];
		}
		localcval[0] = SOL_VAL_LO;
		localcval[1] = SOL_VAL_HI;
		
		I_SMULT( h2, localcval, t1 );
		i_add( local, t1, outvals );
	}
	
	temp = relconst[0];
	if( temp != 0.0 ) {
		I_SMULT( temp, solid_val, t1 );
		ROUND_DOWN;
		val[0] += t1[0];
		c_val[0] += temp*SOL_VAL_LO;
		ROUND_UP;
		val[1] += t1[1];
		c_val[1] += temp*SOL_VAL_HI;
		for( i=0; i<12; i+=2 ) {
			i_smult( temp, solid_i + i, t1 );
			ROUND_DOWN;
			val_i[i] += t1[0];
			ROUND_UP;
			val_i[i+1] += t1[1];
		}
	}

	/* gamma */
	temp = relconst[1];
	if( rel[1] || (temp != 0.0) ) {
		i_gma( yt, sqrtdelta, gma_val );
		s_bvol_xpars( yt, ooyt, delta, sqrtdelta, 
			delta_part, solid_i, bvol_i );
		s_gma_xpars( oosqrtdelta, delta_part, bvol_i, gma_i );
		
		if( rel[1] ) {  /* Compute local bounds . . . */
			local[0] = gma_val[0];
			local[1] = gma_val[1];
			for( i=0; i<12; i+=2 ) {
				i_mult( h + i, gma_i + i, t1 );
				ROUND_DOWN;
				local[0] += t1[0];
				ROUND_UP;
				local[1] += t1[1];
			}
			localcval[0] = GMA_C_LO;
			localcval[1] = GMA_C_HI;
			
			I_SMULT( h2, localcval, t1 );
			i_add( local, t1, outvals + 2 );
		}
		
		if( temp != 0.0 ) {
			I_SMULT( temp, gma_val, t1 );
			ROUND_DOWN;
			val[0] += t1[0];
			c_val[0] += temp*GMA_C_LO;
			ROUND_UP;
			val[1] += t1[1];
			c_val[1] += temp*GMA_C_HI;
			for( i=0; i<12; i+=2 ) {
				i_smult( temp, gma_i + i, t1 );
				ROUND_DOWN;
				val_i[i] += t1[0];
				ROUND_UP;
				val_i[i+1] += t1[1];
			}
		}
	}
	
	/* voronoi */
	temp = relconst[2];
	if( rel[2] || (temp != 0.0) ) {
		i_vor_alt( xt, sqrtdelta, solid_val, vor_val );
		vor_xpars( xt, delta, sqrtdelta, delta_part, 
			solid_i, vor_i );
		
		if( rel[2] ) {  /* Compute local bounds . . . */
			local[0] = vor_val[0];
			local[1] = vor_val[1];
			for( i=0; i<12; i+=2 ) {
				i_mult( h + i, vor_i + i, t1 );
				ROUND_DOWN;
				local[0] += t1[0];
				ROUND_UP;
				local[1] += t1[1];
			}
			localcval[0] = VOR_C_LO;
			localcval[1] = VOR_C_HI;
			
			I_SMULT( h2, localcval, t1 );
			i_add( local, t1, outvals + 4 );
		}

		if( temp != 0.0 ) {
			I_SMULT( temp, vor_val, t1 );
			ROUND_DOWN;
			val[0] += t1[0];
			c_val[0] += temp*VOR_C_LO;
			ROUND_UP;
			val[1] += t1[1];
			c_val[1] += temp*VOR_C_HI;
			for( i=0; i<12; i+=2 ) {
				i_smult( temp, vor_i + i, t1 );
				ROUND_DOWN;
				val_i[i] += t1[0];
				ROUND_UP;
				val_i[i+1] += t1[1];
			}
		}
	}

	/* octahedral voronoi (where y1 is long) */
	temp = relconst[3];
	if( rel[3] || (temp != 0.0) ) {
		i_octavor( yt, xt, sqrtdelta, solid_val, octa_val );
		octa_vor_xpars( yt, ooyt, xt, delta, sqrtdelta, 
			delta_part, solid_i, octa_i );
		
		if( rel[3] ) {  /* Compute local bounds . . . */
			local[0] = octa_val[0];
			local[1] = octa_val[1];
			for( i=0; i<12; i+=2 ) {
				i_mult( h + i, octa_i + i, t1 );
				ROUND_DOWN;
				local[0] += t1[0];
				ROUND_UP;
				local[1] += t1[1];
			}
			localcval[0] = VOR_C_LO;
			localcval[1] = VOR_C_HI;
			
			I_SMULT( h2, localcval, t1 );
			i_add( local, t1, outvals + 6 );
		}

		if( temp != 0.0 ) {
			I_SMULT( temp, octa_val, t1 );
			ROUND_DOWN;
			val[0] += t1[0];
			c_val[0] += temp*VOR_C_LO;
			ROUND_UP;
			val[1] += t1[1];
			c_val[1] += temp*VOR_C_HI;
			for( i=0; i<12; i+=2 ) {
				i_smult( temp, octa_i + i, t1 );
				ROUND_DOWN;
				val_i[i] += t1[0];
				ROUND_UP;
				val_i[i+1] += t1[1];
			}
		}
	}

	/* dihedral (dih1, dih2, dih3) */
	for( j=4; j<7; j++ ) {
		ip = 2*j;
		temp = relconst[j];
		if( rel[j] || (temp != 0.0) ) {
			switch( j ) {
				case 4:
					i_dih( xt, dih_val );
					s_dih_xpars( xt, dih_i );
					break;
					/* permute xt for dih2, dih3 */
				case 5:
					/* dih2:  (1 2 3 4 5 6) -> (2 1 3 5 4 6) */
					for( k=0; k<2; k++ ) {
						xtp[     k] = xt[2  + k];
						xtp[2  + k] = xt[     k];
						xtp[4  + k] = xt[4  + k];
						xtp[6  + k] = xt[8  + k];
						xtp[8  + k] = xt[6  + k];
						xtp[10 + k] = xt[10 + k];
					}
					i_dih( xtp, dih_val );
					s_dih_xpars( xtp, dih_ip );
					/* permute dih_i for dih2, dih3 */
					for( k=0; k<2; k++ ) {
						dih_i[     k] = dih_ip[2  + k];
						dih_i[2  + k] = dih_ip[     k];
						dih_i[4  + k] = dih_ip[4  + k];
						dih_i[6  + k] = dih_ip[8  + k];
						dih_i[8  + k] = dih_ip[6  + k];
						dih_i[10 + k] = dih_ip[10 + k];
					}
					break;
				case 6:
					/* dih3:  (1 2 3 4 5 6) -> (3 2 1 6 5 4) */
					for( k=0; k<2; k++ ) {
						xtp[     k] = xt[4  + k];
						xtp[2  + k] = xt[2  + k];
						xtp[4  + k] = xt[     k];
						xtp[6  + k] = xt[10 + k];
						xtp[8  + k] = xt[8  + k];
						xtp[10 + k] = xt[6  + k];
					}
					i_dih( xtp, dih_val );
					s_dih_xpars( xtp, dih_ip );
					/* permute dih_i for dih2, dih3 */
					for( k=0; k<2; k++ ) {
						dih_i[     k] = dih_ip[4  + k];
						dih_i[2  + k] = dih_ip[2  + k];
						dih_i[4  + k] = dih_ip[     k];
						dih_i[6  + k] = dih_ip[10 + k];
						dih_i[8  + k] = dih_ip[8  + k];
						dih_i[10 + k] = dih_ip[6  + k];
					}
					break;
			}

			if( rel[j] ) {  /* Compute local bounds . . . */
				local[0] = dih_val[0];
				local[1] = dih_val[1];
				for( i=0; i<12; i+=2 ) {
					i_mult( h + i, dih_i + i, t1 );
					ROUND_DOWN;
					local[0] += t1[0];
					ROUND_UP;
					local[1] += t1[1];
				}
				localcval[0] = DIH_VAL_LO;
				localcval[1] = DIH_VAL_HI;
				
				I_SMULT( h2, localcval, t1 );
				i_add( local, t1, outvals + ip );
			}

			if( temp != 0.0 ) {
				I_SMULT( temp, dih_val, t1 );
				ROUND_DOWN;
				val[0] += t1[0];
				c_val[0] += temp*DIH_VAL_LO;
				ROUND_UP;
				val[1] += t1[1];
				c_val[1] += temp*DIH_VAL_HI;
				for( i=0; i<12; i+=2 ) {
					i_smult( temp, dih_i + i, t1 );
					ROUND_DOWN;
					val_i[i] += t1[0];
					ROUND_UP;
					val_i[i+1] += t1[1];
				}
			}
		}
	}
	
	/* Now compute the bounds . . . */
	for( i=0; i<12; i+=2 ) {
		i_mult( h + i, val_i + i, t1 );
		ROUND_DOWN;
		val[0] += t1[0];
		ROUND_UP;
		val[1] += t1[1];
	}
	
	I_SMULT( h2, c_val, t1 );
	I_ADD( val, t1, out );
	
	i_mult( hsum, c_val, t1 );
	ROUND_DOWN;
	for( i=0; i<12; i+=2 ) {
		out_i[i] = t1[0] + val_i[i];
	}
	ROUND_UP;
	for( i=1; i<12; i+=2 ) {
		out_i[i] = t1[1] + val_i[i];
	}
}


/* Here relconst is composed of intervals, not scalars. */
void fat_composite( int rel[7], double relconst[14], 
	double y[12], double x[12], double out[2], 
	double out_i[12], double outvals[14] )
{
	double yt[12], ooyt[12], xt[12], xtp[12];
	double delta[2], sqrtdelta[2], oosqrtdelta[2];
	double delta_part[12], h[12], h2, hsum[2];
	double val[2], val_i[12], c_val[2];
	double solid_val[2], solid_i[12];
	double bvol_i[12];
	double gma_val[2], gma_i[12];
	double vor_val[2], vor_i[12];
	double octa_val[2], octa_i[12];
	double dih_val[2], dih_i[12], dih_ip[12];
	double stemp, *ftemp, t1[2], t2[2], local[2];
	double localcval[2];
	int i, ip, j, k, need_fun[7];
	
	for( j=0; j<7; j++ ) {
		i = 2*j;
		ip = i + 1;
		if( relconst[i] != 0.0 || relconst[ip] != 0.0 ) {
			need_fun[j] = 1;
		} else {
			need_fun[j] = 0;
		}
	}
	
	/* Compute t for both x and y, and h for x */
	ROUND_DOWN;
	for( i=0; i<12; i+=2 ) {
		ip = i + 1;
		xt[i] = 0.5*(x[ip] + x[i]);
		yt[i] = 0.5*(y[ip] + y[i]);
	}
	ROUND_UP;
	for( i=0; i<12; i+=2 ) {
		ip = i + 1;
		xt[ip] = 0.5*(x[ip] + x[i]);
		h[ip] = 0.5*(x[ip] - x[i]);
		h[i] = -h[ip];
		yt[ip] = 0.5*(y[ip] + y[i]);
	}
	
	/* Need ooyt for solid, gma, and octavor */
	/* Compute ooyt = 1/yt and hsum */
	stemp = 0.0;
	ROUND_UP;
	for( i=1; i<12; i+=2 ) {
		ooyt[i] = 1.0/yt[i-1];
		stemp += h[i];
	}
	hsum[1] = stemp;
	ROUND_DOWN;
	stemp = 0.0;
	for( i=0; i<12; i+=2 ) {
		ooyt[i] = 1.0/yt[i+1];
		stemp += h[i];
	}
	hsum[0] = stemp;

	/* Compute h2 */
	ROUND_UP;
	h2 = 0.0;
	for( i=1; i<12; i+=2 ) {
		h2 += h[i]*h[i];
	}
	h2 *= 0.5;
	for( i=0; i<10; i+=2 ) {
		for( j=i+2; j<12; j+=2 ) {
			h2 += h[i]*h[j];
		}
	}
	
	/* Do preliminary computations */
	s_delta_partials( xt, delta_part );
	i_bigdelta( xt, delta );
	I_SQRT( delta, sqrtdelta );
	
	ROUND_DOWN;
	oosqrtdelta[0] = 1.0/sqrtdelta[1];
	ROUND_UP;
	oosqrtdelta[1] = 1.0/sqrtdelta[0];

	/* zero out val, val_i, and c_val */
	val[0] = 0.0;
	val[1] = 0.0;
	c_val[0] = 0.0;
	c_val[1] = 0.0;
	for( i=0; i<12; i++ ) {
		val_i[i] = 0.0;
	}
	/* zero out outvals */
	for( i=0; i<10; i++ ) {
		outvals[i] = 0.0;
	}
	
	/* solid (need this for most things) */
	i_solid( yt, sqrtdelta, solid_val );
	s_solid_xpars( yt, ooyt, delta, sqrtdelta, 
		delta_part, solid_i );

	localcval[0] = SOL_VAL_LO;
	localcval[1] = SOL_VAL_HI;
		
	if( rel[0] ) {  /* Compute local bounds . . . */
		local[0] = solid_val[0];
		local[1] = solid_val[1];
		for( i=0; i<12; i+=2 ) {
			i_mult( h + i, solid_i + i, t1 );
			ROUND_DOWN;
			local[0] += t1[0];
			ROUND_UP;
			local[1] += t1[1];
		}
		I_SMULT( h2, localcval, t1 );
		i_add( local, t1, outvals );
	}
	
	ftemp = relconst;
	if( need_fun[0] ) {
		i_mult( ftemp, solid_val, t1 );
		i_mult( ftemp, localcval, c_val );
		ROUND_DOWN;
		val[0] += t1[0];
		ROUND_UP;
		val[1] += t1[1];
		for( i=0; i<12; i+=2 ) {
			i_mult( ftemp, solid_i + i, t1 );
			ROUND_DOWN;
			val_i[i] += t1[0];
			ROUND_UP;
			val_i[i+1] += t1[1];
		}
	}

	/* gamma */
	ftemp = relconst + 2;
	if( rel[1] || need_fun[1] ) {
		i_gma( yt, sqrtdelta, gma_val );
		s_bvol_xpars( yt, ooyt, delta, sqrtdelta, 
			delta_part, solid_i, bvol_i );
		s_gma_xpars( oosqrtdelta, delta_part, bvol_i, gma_i );
		
		localcval[0] = GMA_C_LO;
		localcval[1] = GMA_C_HI;
			
		if( rel[1] ) {  /* Compute local bounds . . . */
			local[0] = gma_val[0];
			local[1] = gma_val[1];
			for( i=0; i<12; i+=2 ) {
				i_mult( h + i, gma_i + i, t1 );
				ROUND_DOWN;
				local[0] += t1[0];
				ROUND_UP;
				local[1] += t1[1];
			}
			I_SMULT( h2, localcval, t1 );
			i_add( local, t1, outvals + 2 );
		}
		
		if( need_fun[1] ) {
			i_mult( ftemp, gma_val, t1 );
			i_mult( ftemp, localcval, t2 );
			ROUND_DOWN;
			val[0] += t1[0];
			c_val[0] += t2[0];
			ROUND_UP;
			val[1] += t1[1];
			c_val[1] += t2[1];
			for( i=0; i<12; i+=2 ) {
				i_mult( ftemp, gma_i + i, t1 );
				ROUND_DOWN;
				val_i[i] += t1[0];
				ROUND_UP;
				val_i[i+1] += t1[1];
			}
		}
	}
	
	/* voronoi */
	ftemp = relconst + 4;
	if( rel[2] || need_fun[2] ) {
		i_vor_alt( xt, sqrtdelta, solid_val, vor_val );
		vor_xpars( xt, delta, sqrtdelta, delta_part, 
			solid_i, vor_i );
		
		localcval[0] = VOR_C_LO;
		localcval[1] = VOR_C_HI;
			
		if( rel[2] ) {  /* Compute local bounds . . . */
			local[0] = vor_val[0];
			local[1] = vor_val[1];
			for( i=0; i<12; i+=2 ) {
				i_mult( h + i, vor_i + i, t1 );
				ROUND_DOWN;
				local[0] += t1[0];
				ROUND_UP;
				local[1] += t1[1];
			}
			I_SMULT( h2, localcval, t1 );
			i_add( local, t1, outvals + 4 );
		}

		if( need_fun[2] ) {
			i_mult( ftemp, vor_val, t1 );
			i_mult( ftemp, localcval, t2 );
			ROUND_DOWN;
			val[0] += t1[0];
			c_val[0] += t2[0];
			ROUND_UP;
			val[1] += t1[1];
			c_val[1] += t2[1];
			for( i=0; i<12; i+=2 ) {
				i_mult( ftemp, vor_i + i, t1 );
				ROUND_DOWN;
				val_i[i] += t1[0];
				ROUND_UP;
				val_i[i+1] += t1[1];
			}
		}
	}

	/* octahedral voronoi (where y1 is long) */
	ftemp = relconst + 6;
	if( rel[3] || need_fun[3] ) {
		i_octavor( yt, xt, sqrtdelta, solid_val, octa_val );
		octa_vor_xpars( yt, ooyt, xt, delta, sqrtdelta, 
			delta_part, solid_i, octa_i );
		
		localcval[0] = VOR_C_LO;
		localcval[1] = VOR_C_HI;
			
		if( rel[3] ) {  /* Compute local bounds . . . */
			local[0] = octa_val[0];
			local[1] = octa_val[1];
			for( i=0; i<12; i+=2 ) {
				i_mult( h + i, octa_i + i, t1 );
				ROUND_DOWN;
				local[0] += t1[0];
				ROUND_UP;
				local[1] += t1[1];
			}
			I_SMULT( h2, localcval, t1 );
			i_add( local, t1, outvals + 6 );
		}

		if( need_fun[3] ) {
			i_mult( ftemp, octa_val, t1 );
			i_mult( ftemp, localcval, t2 );
			ROUND_DOWN;
			val[0] += t1[0];
			c_val[0] += t2[0];
			ROUND_UP;
			val[1] += t1[1];
			c_val[1] += t2[1];
			for( i=0; i<12; i+=2 ) {
				i_mult( ftemp, octa_i + i, t1 );
				ROUND_DOWN;
				val_i[i] += t1[0];
				ROUND_UP;
				val_i[i+1] += t1[1];
			}
		}
	}

	/* dihedral (dih1, dih2, dih3) */
	for( j=4; j<7; j++ ) {
		ip = 2*j;
		ftemp = relconst + ip;
		if( rel[j] || need_fun[j] ) {
			switch( j ) {
				case 4:
					i_dih( xt, dih_val );
					s_dih_xpars( xt, dih_i );
					break;
					/* permute xt for dih2, dih3 */
				case 5:
					/* dih2:  (1 2 3 4 5 6) -> (2 1 3 5 4 6) */
					for( k=0; k<2; k++ ) {
						xtp[     k] = xt[2  + k];
						xtp[2  + k] = xt[     k];
						xtp[4  + k] = xt[4  + k];
						xtp[6  + k] = xt[8  + k];
						xtp[8  + k] = xt[6  + k];
						xtp[10 + k] = xt[10 + k];
					}
					i_dih( xtp, dih_val );
					s_dih_xpars( xtp, dih_ip );
					/* permute dih_i for dih2, dih3 */
					for( k=0; k<2; k++ ) {
						dih_i[     k] = dih_ip[2  + k];
						dih_i[2  + k] = dih_ip[     k];
						dih_i[4  + k] = dih_ip[4  + k];
						dih_i[6  + k] = dih_ip[8  + k];
						dih_i[8  + k] = dih_ip[6  + k];
						dih_i[10 + k] = dih_ip[10 + k];
					}
					break;
				case 6:
					/* dih3:  (1 2 3 4 5 6) -> (3 2 1 6 5 4) */
					for( k=0; k<2; k++ ) {
						xtp[     k] = xt[4  + k];
						xtp[2  + k] = xt[2  + k];
						xtp[4  + k] = xt[     k];
						xtp[6  + k] = xt[10 + k];
						xtp[8  + k] = xt[8  + k];
						xtp[10 + k] = xt[6  + k];
					}
					i_dih( xtp, dih_val );
					s_dih_xpars( xtp, dih_ip );
					/* permute dih_i for dih2, dih3 */
					for( k=0; k<2; k++ ) {
						dih_i[     k] = dih_ip[4  + k];
						dih_i[2  + k] = dih_ip[2  + k];
						dih_i[4  + k] = dih_ip[     k];
						dih_i[6  + k] = dih_ip[10 + k];
						dih_i[8  + k] = dih_ip[8  + k];
						dih_i[10 + k] = dih_ip[6  + k];
					}
					break;
			}

			localcval[0] = DIH_VAL_LO;
			localcval[1] = DIH_VAL_HI;
				
			if( rel[j] ) {  /* Compute local bounds . . . */
				local[0] = dih_val[0];
				local[1] = dih_val[1];
				for( i=0; i<12; i+=2 ) {
					i_mult( h + i, dih_i + i, t1 );
					ROUND_DOWN;
					local[0] += t1[0];
					ROUND_UP;
					local[1] += t1[1];
				}
				I_SMULT( h2, localcval, t1 );
				i_add( local, t1, outvals + ip );
			}

			if( need_fun[j] ) {
				i_mult( ftemp, dih_val, t1 );
				i_mult( ftemp, localcval, t2 );
				ROUND_DOWN;
				val[0] += t1[0];
				c_val[0] += t2[0];
				ROUND_UP;
				val[1] += t1[1];
				c_val[1] += t2[1];
				for( i=0; i<12; i+=2 ) {
					i_mult( ftemp, dih_i + i, t1 );
					ROUND_DOWN;
					val_i[i] += t1[0];
					ROUND_UP;
					val_i[i+1] += t1[1];
				}
			}
		}
	}
	
	/* Now compute the bounds . . . */
	for( i=0; i<12; i+=2 ) {
		i_mult( h + i, val_i + i, t1 );
		ROUND_DOWN;
		val[0] += t1[0];
		ROUND_UP;
		val[1] += t1[1];
	}
	
	I_SMULT( h2, c_val, t1 );
	I_ADD( val, t1, out );
	
	i_mult( hsum, c_val, t1 );
	ROUND_DOWN;
	for( i=0; i<12; i+=2 ) {
		out_i[i] = t1[0] + val_i[i];
	}
	ROUND_UP;
	for( i=1; i<12; i+=2 ) {
		out_i[i] = t1[1] + val_i[i];
	}
}


void t_print_const( void )
{
	double t1[2], t2[2], t3[2], sol[2], one3[2];
	
	i_init();
	
	one3[0] = ONE_3_LO;
	one3[1] = ONE_3_HI;
	sol[0] = SOL_VAL_LO;
	sol[1] = SOL_VAL_HI;
	
	t1[0] = VOR_VAL_LO;
	t1[1] = VOR_VAL_HI;
	i_mult( t1, i_doct_const, t2 );
	i_mult( sol, one3, t3 );
	i_sub( t3, t2, t1 );
	i_smult( 4.0, t1, t2 );
	ROUND_NEAR;
	printf("vor_c = [%.18g, %.18g]\n", t2[0], t2[1]);
	
	t1[0] = GMA_VAL_LO;
	t1[1] = GMA_VAL_HI;
	i_mult( t1, i_doct_const, t2 );
	i_mult( sol, one3, t3 );
	i_smult( 4.0, t3, t1 );
	i_sub( t1, t2, t3 );
	ROUND_NEAR;
	printf("gma_c = [%.18g, %.18g]\n", t3[0], t3[1]);	
}

/* 2.73982134 2.029102 2.0341248 2.011234 2.023309 2.0123 */
