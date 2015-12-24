/* i_appendix.c, by Samuel Ferguson, (c) 1998. 			*/
/* Auxiliary routines for doing calculations in the 	*/
/* appendix of Sphere Packings IV. 							*/

#include "system_headers.h"
#include "sphere.h"
#include "interval.h"
#include "i_sphere.h"
#include "i_bounds.h"
#include "macros.h"

#define TRASH		0
#define DEBUG		0

/* External variables */

/* Global variables */

int rog_dne;
int rog_inside;

double t0_val[2];
double oot0_val[2];
double phi0_val[2];
double zetapt_val[2];
double sol_coeff[2];
double fake_anc_const[8];
double eta0_const[6];
double crown_const[2];

/* Prototypes */

void appendix_init( void );
void appendix_get_trash( void );
void phi_fun( double h[2], double t[2], double out[2] );
void a_fun( double h[2], double out[2] );
void b_fun( double y[2], double out[2] );
void b0_fun( double y[2], double out[2] );
void v0_fun( double y[12], double x[12], double out[2] );
void v1_fun( double y[12], double x[12], double out[2] );
void rog_vor( double abc[6], double xyz[6], double out[2] );
void rog_sol( double abc[6], double out[2] );
void rog_dih( double xyz[6], double out[2] );
void eta0_fun( double h[2], double out[2] );
void old_eta0_fun( double h[2], double out[2] );
void new_crown_fun( double h[2], double eta[2], double out[2] );
void crown_fun( double h[2], double out[2] );
void anc_den( double abc[6], double xyz[6], double out[2] );
void anc_coeff( double y2[2], double out[2] );
void anc_coeff_y2( double y2[2], double out[2] );
void anc_fun( double y[6], double x[6], double out[2] );
void old_anc_fun( double y[6], double x[6], double out[2] );
void older_anc_fun( double y[6], double x[6], double out[2] );
void max_anc_best( double y[6], double x[6], 
	double anc_2and6[2], double out[2] );
void kappa_fun( double y[12], double x[12], 
	double dih_val[2], double out[2] );
void dih_sign26( double x[12], double out[4] );
void dih_sign35( double x[12], double out[4] );
void set_fake_anc_const( void );
void fake_anc_fun( double y[6], double out[2] );
void rog_vol_ab( double abc[6], double xyz[6], double part[4] );
void rog_sol_ab( double abc[6], double xyz[6], double part[4] );
void rog_dih_ab( double abc[6], double xyz[6], double part[4] );
void rog_vor_ab( double abc[6], double xyz[6], double part[4] );
void anc_y6par( double y[6], double x[6], double out[2] );
void old_anc_y6par( double y[6], double x[6], double out[2] );
void anc_2and6pars( double y[6], double x[6], double out[4] );
void fake_kappa_fun( double y[12], double dih_val[2], 
	double out[2] );
void fake_kappa_fun2( double y[12], double dih_val[2], 
	double out[2] );
void fake_anchor( double y[6], double out[2] );
void fake_anchor2( double y[6], double out[2] );
void fake_anchor_y6( double y[6], double out[2] );
void fake_anchor2_y6( double y[6], double out[2] );
void cos_arc_fun( double y[6], double x[6], double out[2] );
void cos2_beta_psi( double cos2psi[2], double cos2theta[2],
	double out[2] );
void cos2_dih3( double xp[12], double sign[2], double out[2] );
void cos2_dih2( double xp[12], double sign[2], double out[2] );
void cos2_rog_dih( double xyz[6], double out[2] );
void set_eta0_const( void );
void tconst_fun( int i, double out[2] );
void sconst_fun( int index, double out[2] );
void tomDfun( int n, int k, double out[2] );
void tomZfun( int n, int k, double out[2] );


/* Code */


void appendix_init( void )
{
	int i;
	double temp[2], temp2[2], tet[12];
	
	set_fake_anc_const();
	set_eta0_const();
		
	/* initialize local constants */
	t0_val[0] = 0.5*TWO51_LO;
	t0_val[1] = 0.5*TWO51_HI;
	ROUND_DOWN;
	oot0_val[0] = 1.0/t0_val[1];
	ROUND_UP;
	oot0_val[1] = 1.0/t0_val[0];
	phi_fun( t0_val, t0_val, phi0_val );
	/* compute zetapt_val here */
	tet[0] = 2.0;
	tet[1] = 2.0;
	ROUND_DOWN;
	temp[0] = sqrt(tet[0])/5.0;
	ROUND_UP;
	temp[1] = sqrt(tet[1])/5.0;
	I_ATAN( temp, tet );
	temp2[0] = 0.5;
	temp2[1] = 0.5;
	i_div( temp2, tet, temp );
	for( i=0; i<12; i++ ) {
		tet[i] = 2.0;
	}
	/* compute pt */
	temp2[0] = rough_min_gma( tet );
	temp2[1] = rough_max_gma( tet );
	i_mult( temp, temp2, zetapt_val );
	
	/* sol_coeff */
	ROUND_DOWN;
	sol_coeff[0] = phi0_val[0];
	ROUND_UP;
	sol_coeff[1] = phi0_val[1];
}


void appendix_get_trash( void )
{
	double sy[6], sx[6], y[12], x[12], eps, val[4];
	double abc[6], xyz[6], part[4], dih_val[2];
	double dih_part[12];
	int i, j;
	int interval;
	
	appendix_init();

	ROUND_NEAR;
	j = 0;
	printf("Enter data type for cells:  1 for intervals, 0 otherwise:  ");
	scanf("%d", &interval);
	
	while( 1 ) {
		if( interval == 0 ) {
			printf("Enter edge lengths:\n");
			
			for( i=0; i<6; i++ )
				scanf("%lf", sy + i);
			printf("Got the following:  \n");
			for( i=0; i<6; i++ )
				printf("%f\t", sy[i]);
			printf("\n\n");
			printf("Enter epsilon (interval width):  ");
			scanf("%lf", &eps);
			
			ROUND_UP;
			for( i=0; i<6; i++ ) {
				y[2*i] = sy[i];
				y[2*i+1] = sy[i] + eps;
				}
			ROUND_NEAR;
		} else {
		printf("Enter edge lengths (as intervals):\n");
		for( i=0; i<12; i++ )
			scanf("%lf", y + i );
		sy[0] = y[0];
		sy[1] = y[2];
		}
		
		ROUND_DOWN;
		for( i=0; i<12; i+=2 ) {
			x[i] = y[i]*y[i];
		}
		ROUND_UP;
		for( i=1; i<12; i+=2 ) {
			x[i] = y[i]*y[i];
		}
		ROUND_NEAR;
		sx[0] = sy[0]*sy[0];
		sx[1] = sy[1]*sy[1];
		printf("x = \n");
		for( i=0; i<12; i+=2 ) {
			printf("[%.16g, %.16g]\n", x[i], x[i+1]);
		}
		printf("\n");
		
		i_dih_partials( x, y, dih_part );
		i_dih_best( x, dih_part, dih_val );

		kappa_fun( y, x, dih_val, val );
		ROUND_NEAR;
		printf("kappa_val = [%.18g, %.18g]\n", val[0], val[1]);
		printf("\n\n");

		dih_sign26( x, val );
		ROUND_NEAR;
		printf("x2_sign = [%.18g, %.18g]\n", val[0], val[1]);
		printf("x6_sign = [%.18g, %.18g]\n", val[2], val[3]);
		printf("\n\n");

		ROUND_NEAR;
		printf("fake_anc_const = \n");
		for( i=0; i<8; i+=2 ) {
			printf("[%.20f, %.20f]\n", 
				fake_anc_const[i], fake_anc_const[i+1]);
		}

		anc_fun( y, x, val );
		ROUND_NEAR;
		printf("anc_val      = [%.18g, %.18g]\n", val[0], val[1]);
		printf("\n\n");

		old_anc_fun( y, x, val );
		ROUND_NEAR;
		printf("old_anc_val  = [%.18g, %.18g]\n", val[0], val[1]);
		printf("\n\n");

		fake_anc_fun( y, val );
		ROUND_NEAR;
		printf("fake_anc_val = [%.18g, %.18g]\n", val[0], val[1]);
		printf("\n\n");
		
		abc[0] = 0.5*y[0];
		abc[1] = 0.5*y[1];
		xyz[0] = 0.25*x[0];
		xyz[1] = 0.25*x[1];
		
		eta0_fun( abc, val );
		old_eta0_fun( abc, dih_val );
		ROUND_DOWN;
		part[0] = val[0] - dih_val[0];
		ROUND_UP;
		part[1] = dih_val[1] - val[1];
		ROUND_NEAR;
		printf("eta0_fun     = [%.18g, %.18g]\n", val[0], val[1]);
		printf("old_eta0_fun = [%.18g, %.18g]\n", dih_val[0], dih_val[1]);
		printf("eta0_diff = [%.18g, %.18g]\n",
			part[0], part[1]);
		printf("\n");

		eta0_fun( abc, dih_val );
		new_crown_fun( abc, dih_val, val );
		crown_fun( abc, dih_val );
		ROUND_DOWN;
		part[1] = val[0] - dih_val[0];
		ROUND_UP;
		part[0] = val[1] - dih_val[1];
		ROUND_NEAR;
		printf("new_crown_fun = [%.18g, %.18g]\n", val[0], val[1]);
		printf("crown_fun     = [%.18g, %.18g]\n", dih_val[0], dih_val[1]);
		printf("crown_diff = [%.18g, %.18g]\n",
			part[0], part[1]);
		printf("\n");
		
		anc_fun( y, x, val );
		old_anc_fun( y, x, dih_val );
		ROUND_DOWN;
		part[1] = val[0] - dih_val[0];
		ROUND_UP;
		part[0] = val[1] - dih_val[1];
		ROUND_NEAR;
		printf("anc_fun     = [%.18g, %.18g]\n", val[0], val[1]);
		printf("old_anc_fun = [%.18g, %.18g]\n", dih_val[0], dih_val[1]);
		printf("anc_diff = [%.18g, %.18g]\n",
			part[0], part[1]);
		printf("\n");
		
		s_crad3x( x, abc + 2 );
		s_crad3x2( x, xyz + 2 );
		abc[4] = SQRT2_LO;
		abc[5] = SQRT2_HI;
		xyz[4] = 2.0;
		xyz[5] = 2.0;
		ROUND_NEAR;
		printf("a = [%.18g, %.18g]\n", abc[0], abc[1]);		
		printf("b = [%.18g, %.18g]\n", abc[2], abc[3]);		
		printf("c = [%.18g, %.18g]\n\n", abc[4], abc[5]);		
		
		rog_sol_ab( abc, xyz, part );
		ROUND_NEAR;
		printf("rs_a     = [%.18g, %.18g]\n", part[0], part[1]);
		printf("rs_b     = [%.18g, %.18g]\n", part[2], part[3]);
		old_rogsol_ab( abc, xyz, part );
		ROUND_NEAR;
		printf("old_rs_a = [%.18g, %.18g]\n", part[0], part[1]);
		printf("old_rs_b = [%.18g, %.18g]\n", part[2], part[3]);
		printf("\n\n");

		rog_vol_ab( abc, xyz, part );
		ROUND_NEAR;
		printf("rv_a     = [%.18g, %.18g]\n", part[0], part[1]);
		printf("rv_b     = [%.18g, %.18g]\n", part[2], part[3]);
		rogvol_ab( abc, xyz, part );
		ROUND_NEAR;
		printf("old_rv_a = [%.18g, %.18g]\n", part[0], part[1]);
		printf("old_rv_b = [%.18g, %.18g]\n", part[2], part[3]);
		printf("\n\n");

		rog_dih_ab( abc, xyz, part );
		ROUND_NEAR;
		printf("rd_a     = [%.18g, %.18g]\n", part[0], part[1]);
		printf("rd_b     = [%.18g, %.18g]\n", part[2], part[3]);

		rog_vor_ab( abc, xyz, part );
		ROUND_NEAR;
		printf("rvor_a     = [%.18g, %.18g]\n", part[0], part[1]);
		printf("rvor_b     = [%.18g, %.18g]\n\n", part[2], part[3]);

		anc_fun( y, x, val );
		ROUND_NEAR;
		printf("anc_val      = [%.18g, %.18g]\n", val[0], val[1]);
		printf("\n\n");

		printf("rog_inside = %d\n", rog_inside);

		anc_2and6pars( y, x, val );
		old_anc_y6par( y, x, dih_val );
		ROUND_DOWN;
		part[1] = val[2] - dih_val[0];
		ROUND_UP;
		part[0] = val[3] - dih_val[1];
		ROUND_NEAR;
		printf("anc_y2     = [%.18g, %.18g]\n", val[0], val[1]);
		printf("anc_y6     = [%.18g, %.18g]\n", val[2], val[3]);
		printf("old_anc_y6 = [%.18g, %.18g]\n", dih_val[0], dih_val[1]);
		printf("ancy6_diff = [%.18g, %.18g]\n",
			part[0], part[1]);
		printf("\n");
	}
}


/* phi is always negative, for t between t0 and Sqrt[2], anyway */
void phi_fun( double h[2], double t[2], double out[2] )
{
	
	ROUND_DOWN;
	out[0] = TWO_3_HI*(2.0-DOCT_HI*h[1]*t[1]*(h[1] + t[1]));
	ROUND_UP;
	out[1] = TWO_3_LO*(2.0-DOCT_LO*h[0]*t[0]*(h[0] + t[0]));
}


void a_fun( double h[2], double out[2] )
{
	double p1[2], p2[2], temp[2];
	
	phi_fun( h, t0_val, temp );
	ROUND_DOWN;
	p1[0] = 1.0 - h[1]*oot0_val[1];
	p2[0] = temp[0] - phi0_val[1];
	ROUND_UP;
	p1[1] = 1.0 - h[0]*oot0_val[0];
	p2[1] = temp[1] - phi0_val[0];
	i_mult( p1, p2, out );
}


void b_fun( double y[2], double out[2] )
{
	double temp[2], aval[2];
	
	temp[0] = 0.5*y[0];
	temp[1] = 0.5*y[1];
	a_fun( temp, aval );
	/* ROUND_UP; */
	out[1] = aval[1] + phi0_val[1];
	ROUND_DOWN;
	out[0] = aval[0] + phi0_val[0];
}


void b0_fun( double y[2], double out[2] )
{
	double temp[2], aval[2];
	
	temp[0] = 0.5*y[0];
	temp[1] = 0.5*y[1];
	a_fun( temp, aval );
	/* ROUND_UP; */
	out[1] = aval[1] + phi0_val[1] - zetapt_val[0];
	ROUND_DOWN;
	out[0] = aval[0] + phi0_val[0] - zetapt_val[1];
}


void v0_fun( double y[12], double x[12], double out[2] )
{
	double delta4[2], delta6[2], u135[2], temp[6];
	double by1[2], by2[2], by3[2];
	double pt1[2], pt2[2], pt3[2];
	
	temp[0] = x[0];
	temp[1] = x[1];
	temp[2] = x[4];
	temp[3] = x[5];
	temp[4] = x[8];
	temp[5] = x[9];
	i_tomsu( temp, u135 );
	ROUND_DOWN;
	temp[0] = u135[0]*y[2];
	ROUND_UP;
	temp[1] = u135[1]*y[3];
	
	b_fun( y, by1 );
	b_fun( y + 2, by2 );
	b_fun( y + 4, by3 );
	
	i_mult( temp, by2, pt2 );

	delta_partial( 3, x, delta4 );
	delta_partial( 5, x, delta6 );
	
	i_mult( delta4, y + 4, temp );
	i_mult( temp, by3, pt3 );
	
	i_mult( delta6, y, temp );
	i_mult( temp, by1, pt1 );
	ROUND_DOWN;
	out[0] = pt2[0] - pt1[1] - pt3[1];
	ROUND_UP;
	out[1] = pt2[1] - pt1[0] - pt3[0];
}


void v1_fun( double y[12], double x[12], double out[2] )
{
	double delta4[2], delta6[2], u135[2], temp[6];
	double by1[2], by2[2], by3[2];
	double pt1[2], pt2[2], pt3[2];
	
	temp[0] = x[0];
	temp[1] = x[1];
	temp[2] = x[4];
	temp[3] = x[5];
	temp[4] = x[8];
	temp[5] = x[9];
	i_tomsu( temp, u135 );
	ROUND_DOWN;
	temp[0] = u135[0]*y[2];
	ROUND_UP;
	temp[1] = u135[1]*y[3];
	
	b0_fun( y, by1 );
	b0_fun( y + 2, by2 );
	b0_fun( y + 4, by3 );
	
	i_mult( temp, by2, pt2 );

	delta_partial( 3, x, delta4 );
	delta_partial( 5, x, delta6 );
	
	i_mult( delta4, y + 4, temp );
	i_mult( temp, by3, pt3 );
	
	i_mult( delta6, y, temp );
	i_mult( temp, by1, pt1 );
	ROUND_DOWN;
	out[0] = pt2[0] - pt1[1] - pt3[1];
	ROUND_UP;
	out[1] = pt2[1] - pt1[0] - pt3[0];
}


void rog_vor( double abc[6], double xyz[6], double out[2] )
{
	double temp[2], rvol[2], den[2], sum[2], delta_oct[2];
	
	delta_oct[0] = DOCT_LO;
	delta_oct[1] = DOCT_HI;
	
	i_rogersvol2( xyz, temp );
	I_SQRT( temp, rvol );
#if DEBUG
	ROUND_NEAR;
	printf("rogvol = [%0.18g\t\t%0.18g]\n", rvol[0], rvol[1]);
#endif
	
	i_rogers_density( abc, den );
#if DEBUG
	ROUND_NEAR;
	printf("rogden = [%0.18g\t\t%0.18g]\n", den[0], den[1]);
#endif
	
	I_SUB( den, delta_oct, temp );
	i_mult( rvol, temp, sum );
	out[0] = 4.0*sum[0];
	out[1] = 4.0*sum[1];
}


/* Solid angle of a Rogers simplex arising in Voronoi
truncation:
	sph = 2*atan( sqrt( ((b-a)*(c-b))/((a+b)*(b+c)) ) ).  */

void rog_sol( double abc[6], double out[2] )
{
	double a[2], b[2], c[2], bp[2], bpsign[2], num[2], den[2];
	double t1, t2;
	int i, gotbpart;
	
	for( i=0; i<2; i++ ) {
		a[i] = abc[i];
		b[i] = abc[i + 2];
		c[i] = abc[i + 4];
	}
		
	/* Decreasing in a, but the b partial depends on
	the sign of a c - b^2, where c = Sqrt[2].  
	In addition, the function is increasing in c.  */
	
	ROUND_DOWN;
	bpsign[0] = a[0]*c[0] - b[1]*b[1];
	ROUND_UP;
	bpsign[1] = a[1]*c[1] - b[0]*b[0];
	
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
		t1 = bp[0] - a[1];
		if( t1 < 0.0 )
			t1 = 0.0;
		t2 = c[0] - bp[0];
		if( t2 < 0.0 )
			t2 = 0.0;
		num[0] = t1*t2;
		den[1] = (a[0] + bp[1])*(bp[1] + c[1]);
		ROUND_UP;
		t1 = bp[1] - a[0];
		t2 = c[1] - bp[1];
		num[1] = t1*t2;
		den[0] = (a[1] + bp[0])*(bp[0] + c[0]);
	} else {	/* no bounded partial for b */
		ROUND_DOWN;
		t1 = b[0] - a[1];
		if( t1 < 0.0 )
			t1 = 0.0;
		t2 = c[0] - b[1];
		if( t2 < 0.0 )
			t2 = 0.0;
		num[0] = t1*t2;
		den[1] = (a[0] + b[0])*(b[0] + c[1]);
		ROUND_UP;
		t1 = b[1] - a[0];
		t2 = c[1] - b[0];
		num[1] = t1*t2;
		den[0] = (a[1] + b[1])*(b[1] + c[0]);
	}
	bp[1] = num[1]/den[1];
	ROUND_DOWN;
	bp[0] = num[0]/den[0];
	I_SQRT( bp, a );
	I_ATAN( a, b );
	out[0] = 2.0*b[0];
	if( out[0] < 0.0 )
		out[0] = 0.0;
	out[1] = 2.0*b[1];
#if DEBUGTRUNC
	ROUND_NEAR;
	printf("sph = [%.18g, %.18g]\n", out[0], out[1]);
#endif
}


/* rog_dih = atan( sqrt( (c^2 - b^2)/(b^2 - a^2) ) ).
Note that dih is increasing in a, decreasing in b.  Use
x = a^2, y = b^2, z = c^2.  This gives
	atan( sqrt( (z - y)/(y - x) ) ) 
*/

void rog_dih( double xyz[6], double out[2] )
{
	double num[2], den[2], temp;
	
	ROUND_DOWN;
	num[0] = xyz[4] - xyz[3];
	if( num[0] < 0.0 ) {
		num[0] = 0.0;
	}
	den[0] = xyz[2] - xyz[1];
	if( den[0] < 0.0 ) {
		den[0] = 0.0;
	}
	ROUND_UP;
	num[1] = xyz[5] - xyz[2];
	den[1] = xyz[3] - xyz[0];
	
	if( num[1] > 0.0 ) {
		temp = num[1]/den[0];
		out[1] = atan( sqrt( temp ) ) + ATANERR;
	}
	else
		out[1] = 0.0;

	ROUND_DOWN;
	if( num[0] > 0.0 ) {
		temp = num[0]/den[1];
		out[0] = atan( sqrt( temp ) ) - ATANERR;
	} else
		out[0] = 0.0;
}


/* increasing in h */
void eta0_fun( double hval[2], double out[2] )
{
	double h2, den[2], hval2[2], temp;
	
	ROUND_UP;
	hval2[1] = hval[1]*hval[1];
	ROUND_DOWN;
	hval2[0] = hval[0]*hval[0];
	
	h2 = hval2[1];
	den[0] = eta0_const[0] + 
		h2*(eta0_const[2] + h2*eta0_const[4]);
	ROUND_UP;
	h2 = hval2[0];
	den[1] = eta0_const[1] + 
		h2*(eta0_const[3] + h2*eta0_const[5]);
	
	temp = hval2[1]/den[0];
	out[1] = sqrt( temp );
	ROUND_DOWN;
	temp = h2/den[1];
	out[0] = sqrt( temp );
}


void old_eta0_fun( double h[2], double out[2] )
{
	double temp[2], x[6];
	
	temp[0] = TWO51_LO;
	temp[1] = TWO51_HI;
	ROUND_DOWN;
	x[0] = 4.0*h[0]*h[0];
	x[2] = 4.0;
	x[4] = temp[0]*temp[0];
	ROUND_UP;
	x[1] = 4.0*h[1]*h[1];
	x[3] = 4.0;
	x[5] = temp[1]*temp[1];
	
	s_crad3x2( x, temp );
	I_SQRT( temp, out );
}


/* This differs from Tom's crown function:  this one
is not multiplied by 2*Pi */
/* Note that the crown is decreasing in h */
void new_crown_fun( double h[2], double eta[2], double out[2] )
{
	double p1, p2, p3, p4[2], p5[2];
	
	/* This is rather confusing . . . */

	ROUND_DOWN;
	p5[0] = 2.0*h[0]*eta[0]*(eta[0] + h[0]);
	ROUND_UP;
	p5[1] = 2.0*h[1]*eta[1]*(eta[1] + h[1]);
	
	p4[0] = crown_const[1] - p5[1];	/* h_hi */
	ROUND_UP;
	p4[1] = crown_const[0] - p5[0];	/* h_lo */
	
	if( p4[0] > 0.0 ) {
		ROUND_UP;
		p3 = 3.0*eta[1];
		ROUND_DOWN;
		p2 = eta[1] - h[1];
		p1 = DOCT_LO*p2/p3;		/* h_hi */
	
		out[0] = p1*p4[0];
	} else {
		ROUND_DOWN;
		p3 = 3.0*eta[1];
		ROUND_UP;
		p2 = eta[1] - h[1];
		p1 = DOCT_HI*p2/p3;		/* h_hi */
		
		ROUND_DOWN;
		out[0] = p1*p4[0];
	}
	if( p4[1] > 0.0 ) {
		ROUND_DOWN;
		p3 = 3.0*eta[0];
		ROUND_UP;
		p2 = eta[0] - h[0];
		p1 = DOCT_HI*p2/p3;		/* h_lo */
	
		out[1] = p1*p4[1];
	} else {
		ROUND_UP;
		p3 = 3.0*eta[0];
		ROUND_DOWN;
		p2 = eta[0] - h[0];
		p1 = DOCT_LO*p2/p3;		/* h_lo */
		
		ROUND_UP;
		out[1] = p1*p4[1];
	}
}


/* This differs from Tom's crown function:  this one
is not multiplied by 2*Pi */
void crown_fun( double h[2], double out[2] )
{
	double temp2[2], eta[2], phi_val[2], diff[2], temp[2];
	
	eta0_fun( h, eta );
	phi_fun( h, eta, phi_val );
	ROUND_DOWN;
	diff[0] = phi_val[0] - phi0_val[1];
	temp[0] = h[0]/eta[1];
	ROUND_UP;
	diff[1] = phi_val[1] - phi0_val[0];
	temp[1] = h[1]/eta[0];
	
	temp2[1] = 1.0 - temp[0];
	ROUND_DOWN;
	temp2[0] = 1.0 - temp[1];
	
	i_mult( diff, temp2, out );
}


void anc_den( double abc[6], double xyz[6], double out[2] )
{
	double temp[2], rvol[2], den[2], sum[2], delta_oct[2];
	
	delta_oct[0] = 4.0*DOCT_LO;
	delta_oct[1] = 4.0*DOCT_HI;
	
	i_rogersvol2( xyz, temp );
	I_SQRT( temp, rvol );
#if DEBUG
	ROUND_NEAR;
	printf("rogvol = [%0.18g\t\t%0.18g]\n", rvol[0], rvol[1]);
#endif
	
	i_rogers_density( abc, den );
#if DEBUG
	ROUND_NEAR;
	printf("rogden = [%0.18g\t\t%0.18g]\n", den[0], den[1]);
#endif

	temp[0] = -phi0_val[1];
	temp[1] = -phi0_val[0];
	ROUND_DOWN;
	sum[0] = (4.0 + 3.0*temp[0])*den[0];
	ROUND_UP;
	sum[1] = (4.0 + 3.0*temp[1])*den[1];

	I_SUB( sum, delta_oct, temp );
	i_mult( rvol, temp, out );
}


void anc_coeff( double y2[2], double out[2] )
{
	double p1, temp;
	
	ROUND_DOWN;
	p1 = ONE_6_LO*(TWO51_LO + 0.5*y2[1]);
	temp = TWO51_LO - y2[1];
	if( temp < 0.0 )
		temp = 0.0;
	out[0] = DOCT_LO*temp*temp*p1;
	ROUND_UP;
	p1 = ONE_6_HI*(TWO51_HI + 0.5*y2[0]);
	temp = TWO51_HI - y2[0];
	out[1] = DOCT_HI*temp*temp*p1;
}


void anc_coeff_y2( double y2[2], double out[2] )
{	
	ROUND_DOWN;
	out[0] = 0.25*DOCT_LO*(y2[0] - TWO51_HI)*(y2[0] + TWO51_LO);
	ROUND_UP;
	out[1] = 0.25*DOCT_HI*(y2[1] - TWO51_LO)*(y2[1] + TWO51_HI);
}


/* Use density to get better bounds. */
void anc_fun( double y[6], double x[6], double out[2] )
{
	int i;
	double rog1[6], rog2[6], rog12[6], rog22[6];
	double eta0_val[2], h[2], temp[2], crad[2], crown_val[2];
	double rsc1[2], rd1[2];
	double rsc2[2], rd2[2];
	double p1[2], p2[2], p4[2], diff[2];
	
	h[0] = 0.5*y[0];
	h[1] = 0.5*y[1];
	rog1[0] = h[0];
	rog1[1] = h[1];
	rog2[0] = 0.5*y[2];
	rog2[1] = 0.5*y[3];
	
	eta0_fun( h, eta0_val );
	new_crown_fun( h, eta0_val, crown_val );
	rog1[4] = eta0_val[0];
	rog1[5] = eta0_val[1];
	rog2[4] = eta0_val[0];
	rog2[5] = eta0_val[1];
	
	h[0] = 0.5*y[2];
	h[1] = 0.5*y[3];
	
	i_crad3x2( x, temp );
	I_SQRT( temp, crad );
	rog1[2] = crad[0];
	rog1[3] = crad[1];
	rog2[2] = crad[0];
	rog2[3] = crad[1];
	
	if( rog1[2] < rog1[0] ) {	/* b_lo < a_lo */
		rog1[2] = rog1[0];
	}
	if( rog2[2] < rog2[0] ) {	/* b_lo < a_lo */
		rog2[2] = rog2[0];
	}
	/* don't make any other adjustments to a, b, c */
#if DEBUG
	ROUND_NEAR;
	printf("rog1 = \n");
	for( i=0; i<6; i+=2 )
		printf("[%0.18g\t%0.18g]\n", rog1[i], rog1[i+1]);
	printf("rog2 = \n");
	for( i=0; i<6; i+=2 )
		printf("[%0.18g\t%0.18g]\n", rog2[i], rog2[i+1]);
#endif

	ROUND_DOWN;
	for( i=0; i<6; i+=2 ) {
		rog12[i] = rog1[i]*rog1[i];
		rog22[i] = rog2[i]*rog2[i];
	}
	ROUND_UP;
	for( i=1; i<6; i+=2 ) {
		rog12[i] = rog1[i]*rog1[i];
		rog22[i] = rog2[i]*rog2[i];
	}
	
	if( rog1[2] > rog1[5] ) {	/* b_lo > c_hi */
		out[0] = 0.0;
		out[1] = 0.0;
		rog_dne = 1;
		return;
	} else {
		if( rog1[0] > rog1[3] ) {	/* a_lo > b_hi */
			for( i=0; i<2; i++ ) {
				rsc1[i] = 0.0;
				rd1[i] = 0.0;
			}
		} else {
			anc_den( rog1, rog12, rsc1 );
			rog_dih( rog12, rd1 );
		}
		if( rog2[0] > rog2[3] ) {	/* a_lo > b_hi */
			for( i=0; i<2; i++ ) {
				rsc2[i] = 0.0;
				rd2[i] = 0.0;
			}
		} else {
			anc_den( rog2, rog22, rsc2 );
			rog_dih( rog22, rd2 );
		}
	}

	/* check to see if we include a discontinuity */
	rog_inside = 1;
	if( rog1[3] >= rog1[4] ) {	/* b_hi >= c_lo */
		rog_inside = 0;
	} else {
		if( rog1[1] >= rog1[2] ) {	/* a_hi >= b_lo */
			rog_inside = 0;
		}
		if( rog2[1] >= rog2[2] ) {	/* a_hi >= b_lo */
			rog_inside = 0;
		}
	}
	
#if DEBUG
	ROUND_NEAR;
	printf("rsc1 = [%0.18g\t\t%0.18g]\n", rsc1[0], rsc1[1]);
	printf("rsc2 = [%0.18g\t\t%0.18g]\n", rsc2[0], rsc2[1]);
	printf("rd1 = [%0.18g\t\t%0.18g]\n", rd1[0], rd1[1]);
	printf("rd2 = [%0.18g\t\t%0.18g]\n", rd2[0], rd2[1]);
#endif

	h[0] = y[2];
	h[1] = y[3];
	anc_coeff( h, diff );

	ROUND_DOWN;
	p1[0] = rd1[0]*(-crown_val[1]);
	p2[0] = rsc1[0] + rsc2[0];
	
	ROUND_UP;
	p1[1] = rd1[1]*(-crown_val[0]);
	p2[1] = rsc1[1] + rsc2[1];
	
	i_mult( diff, rd2, p4 );
	
	ROUND_DOWN;
	out[0] = p1[0] + p2[0] - p4[1];
	ROUND_UP;
	out[1] = p1[1] + p2[1] - p4[0];
}


/* Use density to get better bounds. */
void old_anc_fun( double y[6], double x[6], double out[2] )
{
	int i;
	double rog1[6], rog2[6], rog12[6], rog22[6];
	double eta0_val[2], h[2], temp[2], crad[2], crown_val[2];
	double rsc1[2], rd1[2];
	double rsc2[2], rd2[2];
	double p1[2], p2[2], p4[2], diff[2], phi_val[2];
	
	h[0] = 0.5*y[0];
	h[1] = 0.5*y[1];
	rog1[0] = h[0];
	rog1[1] = h[1];
	rog2[0] = 0.5*y[2];
	rog2[1] = 0.5*y[3];
	
	crown_fun( h, crown_val );
	eta0_fun( h, eta0_val );
	rog1[4] = eta0_val[0];
	rog1[5] = eta0_val[1];
	rog2[4] = eta0_val[0];
	rog2[5] = eta0_val[1];
	
	h[0] = 0.5*y[2];
	h[1] = 0.5*y[3];
	phi_fun( h, t0_val, phi_val );
	
	i_crad3x2( x, temp );
	I_SQRT( temp, crad );
	rog1[2] = crad[0];
	rog1[3] = crad[1];
	rog2[2] = crad[0];
	rog2[3] = crad[1];
	
	if( rog1[2] < rog1[0] ) {	/* b_lo < a_lo */
		rog1[2] = rog1[0];
	}
	if( rog2[2] < rog2[0] ) {	/* b_lo < a_lo */
		rog2[2] = rog2[0];
	}
	/* don't make any other adjustments to a, b, c */
#if DEBUG
	ROUND_NEAR;
	printf("rog1 = \n");
	for( i=0; i<6; i+=2 )
		printf("[%0.18g\t%0.18g]\n", rog1[i], rog1[i+1]);
	printf("rog2 = \n");
	for( i=0; i<6; i+=2 )
		printf("[%0.18g\t%0.18g]\n", rog2[i], rog2[i+1]);
#endif

	ROUND_DOWN;
	for( i=0; i<6; i+=2 ) {
		rog12[i] = rog1[i]*rog1[i];
		rog22[i] = rog2[i]*rog2[i];
	}
	ROUND_UP;
	for( i=1; i<6; i+=2 ) {
		rog12[i] = rog1[i]*rog1[i];
		rog22[i] = rog2[i]*rog2[i];
	}
	
	if( rog1[2] > rog1[5] ) {	/* b_lo > c_hi */
		out[0] = 0.0;
		out[1] = 0.0;
		rog_dne = 1;
		return;
	} else {
		if( rog1[0] > rog1[3] ) {	/* a_lo > b_hi */
			for( i=0; i<2; i++ ) {
				rsc1[i] = 0.0;
				rd1[i] = 0.0;
			}
		} else {
			anc_den( rog1, rog12, rsc1 );
			rog_dih( rog12, rd1 );
		}
		if( rog2[0] > rog2[3] ) {	/* a_lo > b_hi */
			for( i=0; i<2; i++ ) {
				rsc2[i] = 0.0;
				rd2[i] = 0.0;
			}
		} else {
			anc_den( rog2, rog22, rsc2 );
			rog_dih( rog22, rd2 );
		}
	}

	/* check to see if we include a discontinuity */
	rog_inside = 1;
	if( rog1[3] >= rog1[4] ) {	/* b_hi >= c_lo */
		rog_inside = 0;
	} else {
		if( rog1[1] >= rog1[2] ) {	/* a_hi >= b_lo */
			rog_inside = 0;
		}
		if( rog2[1] >= rog2[2] ) {	/* a_hi >= b_lo */
			rog_inside = 0;
		}
	}
	
#if DEBUG
	ROUND_NEAR;
	printf("rsc1 = [%0.18g\t\t%0.18g]\n", rsc1[0], rsc1[1]);
	printf("rsc2 = [%0.18g\t\t%0.18g]\n", rsc2[0], rsc2[1]);
	printf("rd1 = [%0.18g\t\t%0.18g]\n", rd1[0], rd1[1]);
	printf("rd2 = [%0.18g\t\t%0.18g]\n", rd2[0], rd2[1]);
#endif

	ROUND_DOWN;
	p1[0] = rd1[0]*(-crown_val[1]);
	p2[0] = rsc1[0] + rsc2[0];
	temp[0] = 0.25*y[2]*FB251_LO;
	diff[0] = phi_val[0] - phi0_val[1];
	ROUND_UP;
	p1[1] = rd1[1]*(-crown_val[0]);
	p2[1] = rsc1[1] + rsc2[1];
	temp[1] = 0.25*y[3]*FB251_HI;
	diff[1] = phi_val[1] - phi0_val[0];
	ROUND_DOWN;
	h[0] = (1.0 - temp[1])*rd2[0];
	ROUND_UP;
	h[1] = (1.0 - temp[0])*rd2[1];
	i_mult( diff, h, p4 );
	
	ROUND_DOWN;
	out[0] = p1[0] + p2[0] - p4[1];
	ROUND_UP;
	out[1] = p1[1] + p2[1] - p4[0];
}


void older_anc_fun( double y[6], double x[6], double out[2] )
{
	int i;
	double rog1[6], rog2[6], rog12[6], rog22[6];
	double eta0_val[2], h[2], temp[2], crad[2], crown_val[2];
	double rs1[2], rv1[2], rd1[2];
	double rs2[2], rv2[2], rd2[2];
	double p1[2], p2[2], p3[2], p4[2], diff[2], phi_val[2];
	
	h[0] = 0.5*y[0];
	h[1] = 0.5*y[1];
	rog1[0] = h[0];
	rog1[1] = h[1];
	rog2[0] = 0.5*y[2];
	rog2[1] = 0.5*y[3];
	
	crown_fun( h, crown_val );
	eta0_fun( h, eta0_val );
	rog1[4] = eta0_val[0];
	rog1[5] = eta0_val[1];
	rog2[4] = eta0_val[0];
	rog2[5] = eta0_val[1];
	
	h[0] = 0.5*y[2];
	h[1] = 0.5*y[3];
	phi_fun( h, t0_val, phi_val );
	
	i_crad3x2( x, temp );
	I_SQRT( temp, crad );
	rog1[2] = crad[0];
	rog1[3] = crad[1];
	rog2[2] = crad[0];
	rog2[3] = crad[1];
	
	if( rog1[2] < rog1[0] ) {	/* b_lo < a_lo */
		rog1[2] = rog1[0];
	}
	if( rog2[2] < rog2[0] ) {	/* b_lo < a_lo */
		rog2[2] = rog2[0];
	}
	/* don't make any other adjustments to a, b, c */
#if DEBUG
	ROUND_NEAR;
	printf("rog1 = \n");
	for( i=0; i<6; i+=2 )
		printf("[%0.18g\t%0.18g]\n", rog1[i], rog1[i+1]);
	printf("rog2 = \n");
	for( i=0; i<6; i+=2 )
		printf("[%0.18g\t%0.18g]\n", rog2[i], rog2[i+1]);
#endif

	ROUND_DOWN;
	for( i=0; i<6; i+=2 ) {
		rog12[i] = rog1[i]*rog1[i];
		rog22[i] = rog2[i]*rog2[i];
	}
	ROUND_UP;
	for( i=1; i<6; i+=2 ) {
		rog12[i] = rog1[i]*rog1[i];
		rog22[i] = rog2[i]*rog2[i];
	}
	
	if( rog1[2] > rog1[5] ) {	/* b_lo > c_hi */
		out[0] = 0.0;
		out[1] = 0.0;
		rog_dne = 1;
		return;
	} else {
		if( rog1[0] > rog1[3] ) {	/* a_lo > b_hi */
			for( i=0; i<2; i++ ) {
				rs1[i] = 0.0;
				rv1[i] = 0.0;
				rd1[i] = 0.0;
			}
		} else {
			rog_sol( rog1, rs1 );
			rog_vor( rog1, rog12, rv1 );
			rog_dih( rog12, rd1 );
		}
		if( rog2[0] > rog2[3] ) {	/* a_lo > b_hi */
			for( i=0; i<2; i++ ) {
				rs2[i] = 0.0;
				rv2[i] = 0.0;
				rd2[i] = 0.0;
			}
		} else {
			rog_sol( rog2, rs2 );
			rog_vor( rog2, rog22, rv2 );
			rog_dih( rog22, rd2 );
		}
	}

	/* check to see if we include a discontinuity */
	rog_inside = 1;
	if( rog1[3] >= rog1[4] ) {	/* b_hi >= c_lo */
		rog_inside = 0;
	} else {
		if( rog1[1] >= rog1[2] ) {	/* a_hi >= b_lo */
			rog_inside = 0;
		}
		if( rog2[1] >= rog2[2] ) {	/* a_hi >= b_lo */
			rog_inside = 0;
		}
	}
	
#if DEBUG
	ROUND_NEAR;
	printf("rs1 = [%0.18g\t\t%0.18g]\n", rs1[0], rs1[1]);
	printf("rs2 = [%0.18g\t\t%0.18g]\n", rs2[0], rs2[1]);
	printf("rv1 = [%0.18g\t\t%0.18g]\n", rv1[0], rv1[1]);
	printf("rv2 = [%0.18g\t\t%0.18g]\n", rv2[0], rv2[1]);
	printf("rd1 = [%0.18g\t\t%0.18g]\n", rd1[0], rd1[1]);
	printf("rd2 = [%0.18g\t\t%0.18g]\n", rd2[0], rd2[1]);
#endif

	ROUND_DOWN;
	p1[0] = rd1[0]*(-crown_val[1]);
	p2[0] = (rs1[0] + rs2[0])*(-phi0_val[1]);
	p3[0] = rv1[0] + rv2[0];
	temp[0] = 0.25*y[2]*FB251_LO;
	diff[0] = phi_val[0] - phi0_val[1];
	ROUND_UP;
	p1[1] = rd1[1]*(-crown_val[0]);
	p2[1] = (rs1[1] + rs2[1])*(-phi0_val[0]);
	p3[1] = rv1[1] + rv2[1];
	temp[1] = 0.25*y[3]*FB251_HI;
	diff[1] = phi_val[1] - phi0_val[0];
	ROUND_DOWN;
	h[0] = (1.0 - temp[1])*rd2[0];
	ROUND_UP;
	h[1] = (1.0 - temp[0])*rd2[1];
	i_mult( diff, h, p4 );
	
	ROUND_DOWN;
	out[0] = p1[0] + p2[0] + p3[0] - p4[1];
	ROUND_UP;
	out[1] = p1[1] + p2[1] + p3[1] - p4[0];
}


/* Note:  this only gives real bounds on max_anc, it does
not give real bounds on min_anc (but we don't need those). */
void max_anc_best( double y[6], double x[6], 
	double anc_2and6[2], double out[2] )
{
	int i;
	double yp[6], xp[6];
	
	/* Use y2 and y6 partials to improve bounds on anc_fun */

	/* First find anc_max */
	for( i=0; i<6; i++ ) {
		yp[i] = y[i];	/* copy x and y, then modify copies */
		xp[i] = x[i];
	}
	if( anc_2and6[2] > 0 ) {	/* positive sign */
		yp[4] = yp[5];
		xp[4] = xp[5];
	}
	if( anc_2and6[3] < 0 ) {	/* negative sign */
		yp[5] = yp[4];
		xp[5] = xp[4];
	}
	if( anc_2and6[0] > 0 ) {	/* positive sign */
		yp[2] = yp[3];
		xp[2] = xp[3];
	}
	if( anc_2and6[1] < 0 ) {	/* negative sign */
		yp[3] = yp[2];
		xp[3] = xp[2];
	}
	
	anc_fun( yp, xp, out );
}


void kappa_fun( double y[12], double x[12], 
	double dih_val[2], double out[2] )
{
	double h[2], crown_val[2], anc1[2], anc2[2];
	double face[6], face2[6];
	int outside_state, inside_state;
	
	outside_state = 0;	/* assume rog1, rog2 do exist */
	inside_state = 1;  /* assume interval does not include 
								 discontinuity */
	/* (y1, y2, y6) */
	face[0] = y[0];
	face[1] = y[1];
	face[2] = y[2];
	face[3] = y[3];
	face[4] = y[10];
	face[5] = y[11];
	face2[0] = x[0];
	face2[1] = x[1];
	face2[2] = x[2];
	face2[3] = x[3];
	face2[4] = x[10];
	face2[5] = x[11];
	rog_dne = 0;
	anc_fun( face, face2, anc1 );
	if( rog_dne ) {
		outside_state = 1;	/* does not exist */
	}
	if( !rog_inside ) {
		inside_state = 0;
	}
	
	/* (y1, y3, y5) */
	face[2] = y[4];
	face[3] = y[5];
	face[4] = y[8];
	face[5] = y[9];
	face2[2] = x[4];
	face2[3] = x[5];
	face2[4] = x[8];
	face2[5] = x[9];
	rog_dne = 0;
	anc_fun( face, face2, anc2 );
	if( rog_dne ) {
		outside_state += 2;	/* does not exist */
	}
	if( !rog_inside ) {
		inside_state = 0;
	}
	
	h[0] = 0.5*y[0];
	h[1] = 0.5*y[1];
	crown_fun( h, crown_val );
	
	rog_dne = outside_state;
	
	if( inside_state ) {
		rog_inside = 1;	/* can use partials */
	}
	
	ROUND_DOWN;
	out[0] = crown_val[0]*dih_val[1] + anc1[0] + anc2[0];
	ROUND_UP;
	out[1] = crown_val[1]*dih_val[0] + anc1[1] + anc2[1];
}


/*
1	2	3	4	5	6
0	2	4	6	8	10

x2:
pterms = x2*x4 + x1*x5 + x6*(2*x3 + x6)
nterms = (x4 + x5)*x6 + x1*(x4 + x6) + x2*(x5 + x6)

x6:
pterms = x2^2 + x1*x3 + 2*x2*x5 + x4*x6
nterms = x1*(x2 + x4) + x3*x6 + x2*(x3 + x4 + x6)

Recall that the dihedral coefficient is negative.
*/

void dih_sign26( double x[12], double out[4] )
{
	double pterms[4], nterms[4];
	
	ROUND_DOWN;
	pterms[0] = x[2]*x[6] + x[0]*x[8] + x[10]*(2.0*x[4] + x[10]);
	pterms[2] = x[2]*x[2] + x[0]*x[4] + 2.0*x[2]*x[8] + x[6]*x[10];
	nterms[0] = x[10]*(x[6] + x[8]) + x[0]*(x[6] + x[10]) + x[2]*(x[8] + x[10]);
	nterms[2] = x[0]*(x[2] + x[6]) + x[4]*x[10] + x[2]*(x[4] + x[6] + x[10]);
	ROUND_UP;
	pterms[1] = x[3]*x[7] + x[1]*x[9] + x[11]*(2.0*x[5] + x[11]);
	pterms[3] = x[3]*x[3] + x[1]*x[5] + 2.0*x[3]*x[9] + x[7]*x[11];
	nterms[1] = x[11]*(x[7] + x[9]) + x[1]*(x[7] + x[11]) + x[3]*(x[9] + x[11]);
	nterms[3] = x[1]*(x[3] + x[7]) + x[5]*x[11] + x[3]*(x[5] + x[7] + x[11]);
	
	out[1] = pterms[1] - nterms[0];
	out[3] = pterms[3] - nterms[2];
	ROUND_DOWN;
	out[0] = pterms[0] - nterms[1];
	out[2] = pterms[2] - nterms[3];
}


void dih_sign35( double x[12], double out[4] )
{
	int i;
	double xp[12], dihpar[4];
	/* 1	2	3	4	5	6
		1	5	6	4	2	3
		
		0	2	4	6	8	10
		0	8	10	6	2	4
	*/
	for( i=0; i<2; i++ ) {
		xp[     i] = x[     i	];
		xp[2  + i] = x[8  + i	];
		xp[4  + i] = x[10 + i	];
		xp[6  + i] = x[6  + i	];
		xp[8  + i] = x[2  + i	];
		xp[10 + i] = x[4  + i	];
	}
	
	dih_sign26( xp, dihpar );
	
	out[0] = dihpar[2];
	out[1] = dihpar[3];
	out[2] = dihpar[0];
	out[3] = dihpar[1];
}


void set_fake_anc_const( void )
{
	double x, y;
	
	/* -0.022153 */
	x =  -22153.0;
	y = 1000000.0;
	ROUND_DOWN;
	fake_anc_const[0] = x/y;
	ROUND_UP;
	x =  -22153.0;
	y = 1000000.0;
	fake_anc_const[1] = x/y;
	
	/* 0.0389 */
	x =   389.0;
	y = 10000.0;
	ROUND_DOWN;
	fake_anc_const[2] = x/y;
	ROUND_UP;
	x =   389.0;
	y = 10000.0;
	fake_anc_const[3] = x/y;
	
	/* -0.015 */
	x =  -15.0;
	y = 1000.0;
	ROUND_DOWN;
	fake_anc_const[4] = x/y;
	ROUND_UP;
	x =  -15.0;
	y = 1000.0;
	fake_anc_const[5] = x/y;
	
	/* -0.015 */
	fake_anc_const[6] = fake_anc_const[4];
	fake_anc_const[7] = fake_anc_const[5];
}


void fake_anc_fun( double y[6], double out[2] )
{
	ROUND_DOWN;
	out[0] = fake_anc_const[0] + fake_anc_const[2]*y[0] +
				fake_anc_const[4]*y[3] +
				fake_anc_const[6]*y[5];
	ROUND_UP;
	out[1] = fake_anc_const[1] + fake_anc_const[3]*y[1] +
				fake_anc_const[5]*y[2] +
				fake_anc_const[7]*y[4];
}


/* rv_a = sqrt( (c^2-b^2)/(b^2-a^2) )*(b^2 - 2*a^2)/6 */
/* rv_b = a*b*( c^2 + a^2 - 2*b^2 )/sqrt( 
			(b^2 - a^2)*(c^2 - b^2) )/6 */
/* in reduced terms, these become */
/* rv_a = sqrt( (z-y)/(y-x) )*(y-2*x)/6 */
/* rv_b =  a*b*( z+x-2*y )/sqrt( (y-x)*(z-y) )/6 */
void rog_vol_ab( double abc[6], double xyz[6], double part[4] )
{
	double z_y[2], y_x[2], y_2x[2], zpx_2y[2];
	double foo, quot[2], den[2], atb[2], temp[2];
	
	ROUND_DOWN;
	foo = xyz[4] - xyz[3];
	if( foo < 0.0 )
		foo = 0.0;
	z_y[0] = foo;
	foo = xyz[2] - xyz[1];
	if( foo < 0.0 )
		foo = 0.0;
	y_x[0] = foo;
	y_2x[0] = xyz[2] - 2.0*xyz[1];
	zpx_2y[0] = xyz[4] + xyz[0] - 2.0*xyz[3];
	foo = y_x[0]*z_y[0];
	den[0] = 6.0*sqrt( foo );
	atb[0] = abc[0]*abc[2];
	
	ROUND_UP;
	foo = xyz[5] - xyz[2];
	if( foo < 0.0 )
		foo = 0.0;
	z_y[1] = foo;
	y_x[1] = xyz[3] - xyz[0];
	y_2x[1] = xyz[3] - 2.0*xyz[0];
	zpx_2y[1] = xyz[5] + xyz[1] - 2.0*xyz[2];
	foo = y_x[1]*z_y[1];
	den[1] = 6.0*sqrt( foo );
	atb[1] = abc[1]*abc[3];
	foo = z_y[1]/y_x[0];
	quot[1] = ONE_6_HI*sqrt( foo );
	temp[1] = atb[1]/den[0];
	
	ROUND_DOWN;
	foo = z_y[0]/y_x[1];
	quot[0] = ONE_6_LO*sqrt( foo );
	temp[0] = atb[0]/den[1];
	
	i_mult( quot, y_2x, part );
	i_mult( temp, zpx_2y, part + 2 );
}


/* rs_a=-sqrt( (z-y)/(y-x) )/(a+c) */
/* rs_b=(a*c - y)/(b*sqrt( (z-y)*(y-x) )) */
void rog_sol_ab( double abc[6], double xyz[6], double part[4] )
{
	double z_y[2], y_x[2], ac_y[2], apc[2];
	double den[2], num[2];
	double temp;
	
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


/* rogdih_a = a/(z-x)*sqrt((z-y)/(y-x)) */
/* rogdih_b = -b/(z-y)*sqrt((z-y)/(y-x)) */
void rog_dih_ab( double abc[6], double xyz[6], double part[4] )
{
	double z_y[2], y_x[2], z_x[2];
	double quot, temp;
	
	ROUND_DOWN;
	temp = xyz[4] - xyz[3];
	if( temp < 0.0 )
		temp = 0.0;
	z_y[0] = temp;
	temp = xyz[2] - xyz[1];
	if( temp < 0.0 )
		temp = 0.0;
	y_x[0] = temp;
	temp = xyz[4] - xyz[1];
	if( temp < 0.0 )
		temp = 0.0;
	z_x[0] = temp;
	
	ROUND_UP;
	temp = xyz[5] - xyz[2];
	if( temp < 0.0 )
		temp = 0.0;
	z_y[1] = temp;
	y_x[1] = xyz[3] - xyz[0];
	temp = xyz[5] - xyz[0];
	if( temp < 0.0 )
		temp = 0.0;
	z_x[1] = temp;
	quot = sqrt( z_y[1]/y_x[0] );
	part[1] = abc[1]*quot/z_x[0];
	temp = abc[3]*quot/z_y[0];
	part[2] = -temp;
	
	ROUND_DOWN;
	quot = sqrt( z_y[0]/y_x[1] );
	part[0] = abc[0]*quot/z_x[1];
	temp = abc[2]*quot/z_y[1];
	part[3] = -temp;
}


void rog_vor_ab( double abc[6], double xyz[6], double part[4] )
{
	double rv_ab[4], rs_ab[4], delta_oct[2], one_3[2];
	double p1[4], p2[4];
	
	delta_oct[0] = DOCT_LO;
	delta_oct[1] = DOCT_HI;
	one_3[0] = ONE_3_LO;
	one_3[1] = ONE_3_HI;
	
	rog_vol_ab( abc, xyz, rv_ab );
	rog_sol_ab( abc, xyz, rs_ab );
	
	i_mult( rs_ab, one_3, p1 );
	i_mult( rs_ab + 2, one_3, p1 + 2 );
	i_mult( rv_ab, delta_oct, p2 );
	i_mult( rv_ab + 2, delta_oct, p2 + 2 );
	ROUND_DOWN;
	part[0] = 4.0*(p1[0] - p2[1]);
	part[2] = 4.0*(p1[2] - p2[3]);
	ROUND_UP;
	part[1] = 4.0*(p1[1] - p2[0]);
	part[3] = 4.0*(p1[3] - p2[2]);
}


void anc_y6par( double y[6], double x[6], double out[2] )
{
	int i, ab1_dne, ab2_dne;
	double rog1[6], rog2[6], rog12[6], rog22[6];
	double eta0_val[2], h[2], temp[2], crad[2], crown_val[2];
	double rs1[2], rv1[2], rd1[2];
	double rs2[2], rv2[2], rd2[2];
	double p1[2], p2[2], p3[2], p4[2], diff[2];
	double ab_pars[4], part[2], crad3_pars[6];
	
	h[0] = 0.5*y[0];
	h[1] = 0.5*y[1];
	rog1[0] = h[0];
	rog1[1] = h[1];
	rog2[0] = 0.5*y[2];
	rog2[1] = 0.5*y[3];
	
	eta0_fun( h, eta0_val );
	new_crown_fun( h, eta0_val, crown_val );
	rog1[4] = eta0_val[0];
	rog1[5] = eta0_val[1];
	rog2[4] = eta0_val[0];
	rog2[5] = eta0_val[1];
	
	i_crad3x2( x, temp );
	I_SQRT( temp, crad );
	rog1[2] = crad[0];
	rog1[3] = crad[1];
	rog2[2] = crad[0];
	rog2[3] = crad[1];
	
	if( rog1[2] < rog1[0] ) {	/* b_lo < a_lo */
		rog1[2] = rog1[0];
	}
	if( rog2[2] < rog2[0] ) {	/* b_lo < a_lo */
		rog2[2] = rog2[0];
	}
	/* don't make any other adjustments to a, b, c */
#if DEBUG
	ROUND_NEAR;
	printf("rog1 = \n");
	for( i=0; i<6; i+=2 )
		printf("[%0.18g\t%0.18g]\n", rog1[i], rog1[i+1]);
	printf("rog2 = \n");
	for( i=0; i<6; i+=2 )
		printf("[%0.18g\t%0.18g]\n", rog2[i], rog2[i+1]);
#endif

	ROUND_DOWN;
	for( i=0; i<6; i+=2 ) {
		rog12[i] = rog1[i]*rog1[i];
		rog22[i] = rog2[i]*rog2[i];
	}
	ROUND_UP;
	for( i=1; i<6; i+=2 ) {
		rog12[i] = rog1[i]*rog1[i];
		rog22[i] = rog2[i]*rog2[i];
	}
	
	ab1_dne = 0;
	ab2_dne = 0;
	/* check to see if we include a discontinuity */
	if( rog1[3] >= rog1[4] ) {	/* b_hi >= c_lo */
		if( rog1[2] >= rog1[5] ) {	/* b_lo >= c_hi */
			out[0] = 0.0;
			out[1] = 0.0;
			return;
		}
		else {
			set_infinity( out );
			return;
		}
	} else {	/* no bc discontinuity */
		if( rog1[1] >= rog1[2] ) {	/* a_hi >= b_lo */
			if( rog1[0] >= rog1[3] ) {	/* a_lo >= b_hi */
				ab1_dne = 1;
			}
			else {
				set_infinity( out );
				return;
			}
		}
		if( rog2[1] >= rog2[2] ) {	/* a_hi >= b_lo */
			if( rog2[0] >= rog2[3] ) {	/* a_lo >= b_hi */
				ab2_dne = 1;
			}
			else {
				set_infinity( out );
				return;
			}
		}
	}
	
#if DEBUG
	printf("ab1_dne = %d\n", ab1_dne);
	printf("ab2_dne = %d\n", ab2_dne);
#endif

	if( !ab1_dne ) {
		rog_sol_ab( rog1, rog12, ab_pars );
		rs1[0] = ab_pars[2];
		rs1[1] = ab_pars[3];
		rog_vor_ab( rog1, rog12, ab_pars );
		rv1[0] = ab_pars[2];
		rv1[1] = ab_pars[3];
		rog_dih_ab( rog1, rog12, ab_pars );
		rd1[0] = ab_pars[2];
		rd1[1] = ab_pars[3];
	} else {
		for( i=0; i<2; i++ ) {
			rs1[i] = 0.0;
			rv1[i] = 0.0;
			rd1[i] = 0.0;
		}
	}
	if( !ab2_dne ) {
		rog_sol_ab( rog2, rog22, ab_pars );
		rs2[0] = ab_pars[2];
		rs2[1] = ab_pars[3];
		rog_vor_ab( rog2, rog22, ab_pars );
		rv2[0] = ab_pars[2];
		rv2[1] = ab_pars[3];
		rog_dih_ab( rog2, rog22, ab_pars );
		rd2[0] = ab_pars[2];
		rd2[1] = ab_pars[3];
	} else {
		for( i=0; i<2; i++ ) {
			rs2[i] = 0.0;
			rv2[i] = 0.0;
			rd2[i] = 0.0;
		}
	}
	
	i_mult( rd1, crown_val, p1 );
	I_ADD( rs1, rs2, temp );
	i_mult( temp, phi0_val, p2 );
	I_ADD( rv1, rv2, p3 );
	
	h[0] = y[2];
	h[1] = y[3];
	anc_coeff( h, diff );

	i_mult( diff, rd2, p4 );
	
	ROUND_DOWN;
	part[0] = p3[0] - p1[1] - p2[1] - p4[1];
	ROUND_UP;
	part[1] = p3[1] - p1[0] - p2[0] - p4[0];
	/* this is the b partial */
	
	/* now use chain rule */
	rog1[0] = y[4];
	rog1[1] = y[5];
	rog1[2] = y[2];
	rog1[3] = y[3];
	rog1[4] = y[0];
	rog1[5] = y[1];
	
	rog2[0] = x[4];
	rog2[1] = x[5];
	rog2[2] = x[2];
	rog2[3] = x[3];
	rog2[4] = x[0];
	rog2[5] = x[1];
	
	crad3len_pars( rog1, rog2, crad3_pars );
	i_mult( crad3_pars, part, out );
}


void old_anc_y6par( double y[6], double x[6], double out[2] )
{
	int i, ab1_dne, ab2_dne;
	double rog1[6], rog2[6], rog12[6], rog22[6];
	double eta0_val[2], h[2], temp[2], crad[2], crown_val[2];
	double rs1[2], rv1[2], rd1[2];
	double rs2[2], rv2[2], rd2[2];
	double p1[2], p2[2], p3[2], p4[2], diff[2], phi_val[2];
	double ab_pars[4], part[2], crad3_pars[6];
	
	h[0] = 0.5*y[0];
	h[1] = 0.5*y[1];
	rog1[0] = h[0];
	rog1[1] = h[1];
	rog2[0] = 0.5*y[2];
	rog2[1] = 0.5*y[3];
	
	crown_fun( h, crown_val );
	eta0_fun( h, eta0_val );
	rog1[4] = eta0_val[0];
	rog1[5] = eta0_val[1];
	rog2[4] = eta0_val[0];
	rog2[5] = eta0_val[1];
	
	h[0] = 0.5*y[2];
	h[1] = 0.5*y[3];
	phi_fun( h, t0_val, phi_val );
	
	i_crad3x2( x, temp );
	I_SQRT( temp, crad );
	rog1[2] = crad[0];
	rog1[3] = crad[1];
	rog2[2] = crad[0];
	rog2[3] = crad[1];
	
	if( rog1[2] < rog1[0] ) {	/* b_lo < a_lo */
		rog1[2] = rog1[0];
	}
	if( rog2[2] < rog2[0] ) {	/* b_lo < a_lo */
		rog2[2] = rog2[0];
	}
	/* don't make any other adjustments to a, b, c */
#if DEBUG
	ROUND_NEAR;
	printf("rog1 = \n");
	for( i=0; i<6; i+=2 )
		printf("[%0.18g\t%0.18g]\n", rog1[i], rog1[i+1]);
	printf("rog2 = \n");
	for( i=0; i<6; i+=2 )
		printf("[%0.18g\t%0.18g]\n", rog2[i], rog2[i+1]);
#endif

	ROUND_DOWN;
	for( i=0; i<6; i+=2 ) {
		rog12[i] = rog1[i]*rog1[i];
		rog22[i] = rog2[i]*rog2[i];
	}
	ROUND_UP;
	for( i=1; i<6; i+=2 ) {
		rog12[i] = rog1[i]*rog1[i];
		rog22[i] = rog2[i]*rog2[i];
	}
	
	ab1_dne = 0;
	ab2_dne = 0;
	/* check to see if we include a discontinuity */
	if( rog1[3] >= rog1[4] ) {	/* b_hi >= c_lo */
		if( rog1[2] >= rog1[5] ) {	/* b_lo >= c_hi */
			out[0] = 0.0;
			out[1] = 0.0;
			return;
		}
		else {
			set_infinity( out );
			return;
		}
	} else {	/* no bc discontinuity */
		if( rog1[1] >= rog1[2] ) {	/* a_hi >= b_lo */
			if( rog1[0] >= rog1[3] ) {	/* a_lo >= b_hi */
				ab1_dne = 1;
			}
			else {
				set_infinity( out );
				return;
			}
		}
		if( rog2[1] >= rog2[2] ) {	/* a_hi >= b_lo */
			if( rog2[0] >= rog2[3] ) {	/* a_lo >= b_hi */
				ab2_dne = 1;
			}
			else {
				set_infinity( out );
				return;
			}
		}
	}
	
#if DEBUG
	printf("ab1_dne = %d\n", ab1_dne);
	printf("ab2_dne = %d\n", ab2_dne);
#endif

	if( !ab1_dne ) {
		rog_sol_ab( rog1, rog12, ab_pars );
		rs1[0] = ab_pars[2];
		rs1[1] = ab_pars[3];
		rog_vor_ab( rog1, rog12, ab_pars );
		rv1[0] = ab_pars[2];
		rv1[1] = ab_pars[3];
		rog_dih_ab( rog1, rog12, ab_pars );
		rd1[0] = ab_pars[2];
		rd1[1] = ab_pars[3];
	} else {
		for( i=0; i<2; i++ ) {
			rs1[i] = 0.0;
			rv1[i] = 0.0;
			rd1[i] = 0.0;
		}
	}
	if( !ab2_dne ) {
		rog_sol_ab( rog2, rog22, ab_pars );
		rs2[0] = ab_pars[2];
		rs2[1] = ab_pars[3];
		rog_vor_ab( rog2, rog22, ab_pars );
		rv2[0] = ab_pars[2];
		rv2[1] = ab_pars[3];
		rog_dih_ab( rog2, rog22, ab_pars );
		rd2[0] = ab_pars[2];
		rd2[1] = ab_pars[3];
	} else {
		for( i=0; i<2; i++ ) {
			rs2[i] = 0.0;
			rv2[i] = 0.0;
			rd2[i] = 0.0;
		}
	}
	
	i_mult( rd1, crown_val, p1 );
	I_ADD( rs1, rs2, temp );
	i_mult( temp, phi0_val, p2 );
	I_ADD( rv1, rv2, p3 );
	
	ROUND_DOWN;
	temp[0] = 0.25*y[2]*FB251_LO;
	diff[0] = phi_val[0] - phi0_val[1];
	ROUND_UP;
	temp[1] = 0.25*y[3]*FB251_HI;
	diff[1] = phi_val[1] - phi0_val[0];
	part[1] = 1.0 - temp[0];
	ROUND_DOWN;
	part[0] = 1.0 - temp[1];
	i_mult( part, rd2, h );
	i_mult( diff, h, p4 );
	
	ROUND_DOWN;
	part[0] = p3[0] - p1[1] - p2[1] - p4[1];
	ROUND_UP;
	part[1] = p3[1] - p1[0] - p2[0] - p4[0];
	/* this is the b partial */
	
	/* now use chain rule */
	rog1[0] = y[4];
	rog1[1] = y[5];
	rog1[2] = y[2];
	rog1[3] = y[3];
	rog1[4] = y[0];
	rog1[5] = y[1];
	
	rog2[0] = x[4];
	rog2[1] = x[5];
	rog2[2] = x[2];
	rog2[3] = x[3];
	rog2[4] = x[0];
	rog2[5] = x[1];
	
	crad3len_pars( rog1, rog2, crad3_pars );
	i_mult( crad3_pars, part, out );
}


void anc_2and6pars( double y[6], double x[6], double out[4] )
{
	int i, ab1_dne, ab2_dne;
	double rog1[6], rog2[6], rog12[6], rog22[6];
	double eta0_val[2], h[2], temp[2], crad[2], crown_val[2];
	double rs1[2], rv1[2], rd1[2];
	double rd2[2];
	double rs2ab[4], rv2ab[4], rd2ab[4];
	double p1[2], p2[2], p3[2], p4[2], diff[2];
	double ab_pars[4], part[2], crad3_pars[6];
	
	h[0] = 0.5*y[0];
	h[1] = 0.5*y[1];
	rog1[0] = h[0];
	rog1[1] = h[1];
	rog2[0] = 0.5*y[2];
	rog2[1] = 0.5*y[3];
	
	eta0_fun( h, eta0_val );
	new_crown_fun( h, eta0_val, crown_val );
	rog1[4] = eta0_val[0];
	rog1[5] = eta0_val[1];
	rog2[4] = eta0_val[0];
	rog2[5] = eta0_val[1];
	
	i_crad3x2( x, temp );
	I_SQRT( temp, crad );
	rog1[2] = crad[0];
	rog1[3] = crad[1];
	rog2[2] = crad[0];
	rog2[3] = crad[1];
	
	if( rog1[2] < rog1[0] ) {	/* b_lo < a_lo */
		rog1[2] = rog1[0];
	}
	if( rog2[2] < rog2[0] ) {	/* b_lo < a_lo */
		rog2[2] = rog2[0];
	}
	/* don't make any other adjustments to a, b, c */
#if DEBUG
	ROUND_NEAR;
	printf("rog1 = \n");
	for( i=0; i<6; i+=2 )
		printf("[%0.18g\t%0.18g]\n", rog1[i], rog1[i+1]);
	printf("rog2 = \n");
	for( i=0; i<6; i+=2 )
		printf("[%0.18g\t%0.18g]\n", rog2[i], rog2[i+1]);
#endif

	ROUND_DOWN;
	for( i=0; i<6; i+=2 ) {
		rog12[i] = rog1[i]*rog1[i];
		rog22[i] = rog2[i]*rog2[i];
	}
	ROUND_UP;
	for( i=1; i<6; i+=2 ) {
		rog12[i] = rog1[i]*rog1[i];
		rog22[i] = rog2[i]*rog2[i];
	}
	
	ab1_dne = 0;
	ab2_dne = 0;
	/* check to see if we include a discontinuity */
	if( rog1[3] >= rog1[4] ) {	/* b_hi >= c_lo */
		if( rog1[2] >= rog1[5] ) {	/* b_lo >= c_hi */
			out[0] = 0.0;
			out[1] = 0.0;
			return;
		}
		else {
			set_infinity( out );
			return;
		}
	} else {	/* no bc discontinuity */
		if( rog1[1] >= rog1[2] ) {	/* a_hi >= b_lo */
			if( rog1[0] >= rog1[3] ) {	/* a_lo >= b_hi */
				ab1_dne = 1;
			}
			else {
				set_infinity( out );
				return;
			}
		}
		if( rog2[1] >= rog2[2] ) {	/* a_hi >= b_lo */
			if( rog2[0] >= rog2[3] ) {	/* a_lo >= b_hi */
				ab2_dne = 1;
			}
			else {
				set_infinity( out );
				return;
			}
		}
	}
	
#if DEBUG
	printf("ab1_dne = %d\n", ab1_dne);
	printf("ab2_dne = %d\n", ab2_dne);
#endif

	if( !ab1_dne ) {
		rog_sol_ab( rog1, rog12, ab_pars );
		rs1[0] = ab_pars[2];
		rs1[1] = ab_pars[3];
		rog_vor_ab( rog1, rog12, ab_pars );
		rv1[0] = ab_pars[2];
		rv1[1] = ab_pars[3];
		rog_dih_ab( rog1, rog12, ab_pars );
		rd1[0] = ab_pars[2];
		rd1[1] = ab_pars[3];
	} else {
		for( i=0; i<2; i++ ) {
			rs1[i] = 0.0;
			rv1[i] = 0.0;
			rd1[i] = 0.0;
		}
	}
	if( !ab2_dne ) {
		rog_sol_ab( rog2, rog22, rs2ab );
		rog_vor_ab( rog2, rog22, rv2ab );
		rog_dih_ab( rog2, rog22, rd2ab );
		rog_dih( rog22, rd2 );
	} else {
		for( i=0; i<4; i++ ) {
			rs2ab[i] = 0.0;
			rv2ab[i] = 0.0;
			rd2ab[i] = 0.0;
		}
		rd2[0] = 0.0;
		rd2[1] = 0.0;
	}
	
	i_mult( rd1, crown_val, p1 );
	ROUND_DOWN;
	temp[0] = rs1[0] + rs2ab[2];
	p3[0] = rv1[0] + rv2ab[2];
	ROUND_UP;
	temp[1] = rs1[1] + rs2ab[3];
	p3[1] = rv1[1] + rv2ab[3];
	i_mult( temp, phi0_val, p2 );
	
	h[0] = y[2];
	h[1] = y[3];
	anc_coeff( h, diff );

	i_mult( diff, rd2ab + 2, p4 );
	
	ROUND_DOWN;
	part[0] = p3[0] - p1[1] - p2[1] - p4[1];
	ROUND_UP;
	part[1] = p3[1] - p1[0] - p2[0] - p4[0];
	/* this is the b partial */
	
	/* now use chain rule */
	
	crad3len_pars( y, x, crad3_pars );
	i_mult( crad3_pars + 4, part, out + 2 );
	/* Takes care of the y6 partial, now finish y2 partial */
	
	i_mult( crad3_pars + 2, part, p1 );
	i_mult( diff, rd2ab, p2 );
	i_mult( phi0_val, rs2ab, p3 );
	
	anc_coeff_y2( h, diff );
	i_mult( diff, rd2, p4 );
	ROUND_DOWN;
	out[0] = p1[0] - p4[1] + 0.5*(rv2ab[0] - p2[1] - p3[1]);
	ROUND_UP;
	out[1] = p1[1] - p4[1] + 0.5*(rv2ab[1] - p2[0] - p3[0]);
}


/* Since anc_fun < fake_anc_fun (on the appropriate domain of
(664200787)), we have kappa_fun < fake_kappa_fun on the
domain.  Note that only the upper bound of fake_kappa_fun
actually says anything useful.  */

void fake_kappa_fun( double y[12], double dih_val[2], 
	double out[2] )
{
	double h[2], crown_val[2], anc1[2], anc2[2];
	double face[6];
	
	/* (y1, y2, y6) */
	face[0] = y[0];
	face[1] = y[1];
	face[2] = y[2];
	face[3] = y[3];
	face[4] = y[10];
	face[5] = y[11];
	fake_anchor( face, anc1 );
	
	/* (y1, y3, y5) */
	face[2] = y[4];
	face[3] = y[5];
	face[4] = y[8];
	face[5] = y[9];
	fake_anchor( face, anc2 );
	
	h[0] = 0.5*y[0];
	h[1] = 0.5*y[1];
	crown_fun( h, crown_val );
	
	ROUND_DOWN;
	out[0] = crown_val[0]*dih_val[1] + anc1[0] + anc2[0];
	ROUND_UP;
	out[1] = crown_val[1]*dih_val[0] + anc1[1] + anc2[1];
}


void fake_kappa_fun2( double y[12], double dih_val[2], 
	double out[2] )
{
	double h[2], crown_val[2], anc1[2], anc2[2];
	double face[6];
	
	/* (y1, y2, y6) */
	face[0] = y[0];
	face[1] = y[1];
	face[2] = y[2];
	face[3] = y[3];
	face[4] = y[10];
	face[5] = y[11];
	fake_anchor2( face, anc1 );
	
	/* (y1, y3, y5) */
	face[2] = y[4];
	face[3] = y[5];
	face[4] = y[8];
	face[5] = y[9];
	fake_anchor2( face, anc2 );
	
	h[0] = 0.5*y[0];
	h[1] = 0.5*y[1];
	crown_fun( h, crown_val );
	
	ROUND_DOWN;
	out[0] = crown_val[0]*dih_val[1] + anc1[0] + anc2[0];
	ROUND_UP;
	out[1] = crown_val[1]*dih_val[0] + anc1[1] + anc2[1];
}


/* increasing in y1 */
/* decreasing in x  */
void fake_anchor( double y[6], double out[2] )
{
	double x_hi, x_lo, y1_hi, y1_lo, ypart, xpart, val;
	
	y1_lo = y[0];
	y1_hi = y[1];
	ROUND_UP;
	x_hi = y[3] + y[5];
	ROUND_DOWN;
	x_lo = y[2] + y[4];
	if( x_hi < 4.58 ) {
		val = 0.1833230667013778 - 0.02783887375001181*y1_lo;
		ypart = 0.471224150208274927 + y1_lo*val;
		xpart = -0.3239044460786886 + 0.0346916048615081*x_hi;
		out[0] = ypart + x_hi*xpart;
		ROUND_UP;
		val = 0.1833230667013778 - 0.02783887375001181*y1_hi;
		ypart = 0.471224150208274927 + y1_hi*val;
		xpart = -0.3239044460786886 + 0.0346916048615081*x_lo;
		out[1] = ypart + x_lo*xpart;
		return;
	} else if( x_lo > 4.58 ) {
		out[0] = 0.0;
		out[1] = 0.0;
		return;
	} else {	/* mixed case */
		val = 0.1833230667013778 - 0.02783887375001181*y1_lo;
		ypart = 0.471224150208274927 + y1_lo*val;
		xpart = -0.3239044460786886 + 0.0346916048615081*x_lo;
		val = ypart + x_hi*xpart;
		if( val > 0.0 ) {
			val = 0.0;
		}
		out[0] = val;
		ROUND_UP;
		val = 0.1833230667013778 - 0.02783887375001181*y1_hi;
		ypart = 0.471224150208274927 + y1_hi*val;
		xpart = -0.3239044460786886 + 0.0346916048615081*x_hi;
		val = ypart + x_lo*xpart;
		if( val < 0.0 ) {
			val = 0.0;
		}
		out[1] = val;
	}
}


/* increasing in y1 */
/* decreasing in x  */
void fake_anchor2( double y[6], double out[2] )
{
	double x_hi, x_lo, y1_hi, y1_lo, ypart, xpart, val;
	
	y1_lo = y[0];
	y1_hi = y[1];
	ROUND_UP;
	x_hi = y[3] + y[5];
	ROUND_DOWN;
	x_lo = y[2] + y[4];
	val = 0.1833230667013778 - 0.02783887375001181*y1_lo;
	ypart = 0.153841877737499288 + y1_lo*val;
	xpart = -0.167484482694134 + 0.01541738104479281*x_hi;
	out[0] = ypart + x_hi*xpart;
	ROUND_UP;
	val = 0.1833230667013778 - 0.02783887375001181*y1_hi;
	ypart = 0.153841877737499288 + y1_hi*val;
	xpart = -0.167484482694134 + 0.01541738104479281*x_lo;
	out[1] = ypart + x_lo*xpart;
}


/*
-0.323904446078688623`+0.0693832097230162059` (y2+y6)
*/

void fake_anchor_y6( double y[6], double out[2] )
{
	double x_hi, x_lo;
	
	ROUND_UP;
	x_hi = y[3] + y[5];
	ROUND_DOWN;
	x_lo = y[2] + y[4];
	if( x_hi < 4.58 ) {
		out[0] = -0.3239044460786886 + 0.0693832097230162059*x_lo;
		ROUND_UP;
		out[1] = -0.3239044460786886 + 0.0693832097230162059*x_hi;
		return;
	} else if( x_lo > 4.58 ) {
		out[0] = 0.0;
		out[1] = 0.0;
		return;
	} else {	/* mixed case */
		set_infinity( out );
		return;
	}
}

/*
-0.16748448269413398`+0.0308347620895856211` (y2+y6)
*/

void fake_anchor2_y6( double y[6], double out[2] )
{
	double x_hi, x_lo;
	
	ROUND_DOWN;
	x_lo = y[2] + y[4];
	out[0] = -0.16748448269413398 + 0.0308347620895856211*x_lo;
	ROUND_UP;
	x_hi = y[3] + y[5];
	out[1] = -0.16748448269413398 + 0.0308347620895856211*x_hi;
}


void cos_arc_fun( double y[6], double x[6], double out[2] )
{
	double num[2], den[2];
	
	ROUND_DOWN;
	num[0] = x[0] + x[2] - x[5];
	den[0] = 2.0*y[0]*y[2];
	
	ROUND_UP;
	num[1] = x[1] + x[3] - x[4];
	den[1] = 2.0*y[1]*y[3];
	
	i_div( num, den, out );
}


void cos2_beta_psi( double cos2psi[2], double cos2theta[2],
	double out[2] )
{
	double num[2], den[2];
	
	ROUND_DOWN;
	num[0] = cos2psi[0] - cos2theta[1];
	den[0] = 1.0 - cos2theta[1];
	ROUND_UP;
	num[1] = cos2psi[1] - cos2theta[0];
	den[1] = 1.0 - cos2theta[0];
	i_div( num, den, out );
}


void cos2_dih3( double xp[12], double sign[2], double out[2] )
{
	int i;
	double x[12];
	double xv[6], u1[2], u2[2];
	double num[2], den[2], pterms[2], nterms[2];
	
	/* (1 2 3 4 5  6) -> (3 2 1  6 5 4), so
		(0 2 4 6 8 10) -> (4 2 0 10 8 6) 	*/
	for( i=0; i<2; i++ ) {
		x[i     ] = xp[i +  4];
		x[i +  2] = xp[i +  2];
		x[i +  4] = xp[i     ];
		x[i +  6] = xp[i + 10];
		x[i +  8] = xp[i +  8];
		x[i + 10] = xp[i +  6];
	}
	
	ROUND_DOWN;
	pterms[0] = x[0]*(x[2] + x[4] + x[8] + x[10]) +
		x[2]*x[8] + x[4]*x[10];
	nterms[0] = x[0]*(x[0] + 2.0*x[6]) + x[2]*x[4] + x[8]*x[10];
	ROUND_UP;
	pterms[1] = x[1]*(x[3] + x[5] + x[9] + x[11]) +
		x[3]*x[9] + x[5]*x[11];
	nterms[1] = x[1]*(x[1] + 2.0*x[7]) + x[3]*x[5] + x[9]*x[11];
	I_SUB( pterms, nterms, sign );
	i_mult( sign, sign, num );
	
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
	
	I_MULT( u1, u2, den );
	i_div( num, den, out );
}


void cos2_dih2( double xp[12], double sign[2], double out[2] )
{
	int i;
	double x[12];
	double xv[6], u1[2], u2[2];
	double num[2], den[2], pterms[2], nterms[2];
	
	/* (1 2 3 4 5  6) -> (2 1 3 5 4  6), so
		(0 2 4 6 8 10) -> (2 0 4 8 6 10) 	*/
	for( i=0; i<2; i++ ) {
		x[i     ] = xp[i +  2];
		x[i +  2] = xp[i     ];
		x[i +  4] = xp[i +  4];
		x[i +  6] = xp[i +  8];
		x[i +  8] = xp[i +  6];
		x[i + 10] = xp[i + 10];
	}
	
	ROUND_DOWN;
	pterms[0] = x[0]*(x[2] + x[4] + x[8] + x[10]) +
		x[2]*x[8] + x[4]*x[10];
	nterms[0] = x[0]*(x[0] + 2.0*x[6]) + x[2]*x[4] + x[8]*x[10];
	ROUND_UP;
	pterms[1] = x[1]*(x[3] + x[5] + x[9] + x[11]) +
		x[3]*x[9] + x[5]*x[11];
	nterms[1] = x[1]*(x[1] + 2.0*x[7]) + x[3]*x[5] + x[9]*x[11];
	I_SUB( pterms, nterms, sign );
	i_mult( sign, sign, num );
	
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
	
	I_MULT( u1, u2, den );
	i_div( num, den, out );
}


/* (cos(rog_dih))^2 = (y-x)/(z-x) */
void cos2_rog_dih( double xyz[6], double out[2] )
{
	double temp, num[2], den[2];
	
	ROUND_DOWN;
	temp = xyz[2] - xyz[1];
	if( temp < 0 ) {
		temp = 0.0;
	}
	num[0] = temp;
	temp = xyz[4] - xyz[1];
	if( temp < 0 ) {
		temp = 0.0;
	}
	den[0] = temp;
	ROUND_UP;
	temp = xyz[3] - xyz[0];
	if( temp < 0 ) {
		temp = 0.0;
	}
	num[1] = temp;
	temp = xyz[5] - xyz[0];
	if( temp < 0 ) {
		temp = 0.0;
	}
	den[1] = temp;
	out[1] = num[1]/den[0];
	ROUND_DOWN;
	out[0] = num[0]/den[1];
}


void set_eta0_const( void )
{
	double eta0rat[6] = 
		{529046001.0, 10080160000.0,
		103001.0, 126002.0,
		10000.0, 63001.0};
	double crownrat[2] = 
		{15813251.0, 2000000.0};
	double temp[6];
	int i;
		
	ROUND_DOWN;
	for( i=0; i<6; i+=2 ) {
		temp[i] = eta0rat[i]/eta0rat[i+1];
	}
	crown_const[0] = crownrat[0]/crownrat[1];
	ROUND_UP;
	for( i=0; i<6; i+=2 ) {
		temp[i+1] = eta0rat[i]/eta0rat[i+1];
	}
	crown_const[1] = crownrat[0]/crownrat[1];
	
	eta0_const[0] = - temp[1];
	eta0_const[1] = - temp[0];
	eta0_const[2] = temp[2];
	eta0_const[3] = temp[3];
	eta0_const[4] = - temp[5];
	eta0_const[5] = - temp[4];
	
	for( i=0; i<6; i+=2 ) {
		if( eta0_const[i] == eta0_const[i+1] ) {
			printf("Problem with eta0_const[%d]\n", i/2+1);
		}
	}
	if( crown_const[0] == crown_const[1] ) {
		printf("Problem with crown_const\n");
	}
/*	
	ROUND_NEAR;
	printf("eta0rat = \n");
	for( i=0; i<6; i++ ) {
		printf("%30.18f\t", eta0rat[i]);
		if( i%2 )
			printf("\n");
	}
	for( i=0; i<6; i+=2 ) {
		printf("eta0_const[%d] = [%.18f, %.18f]\n",
			i/2+1, eta0_const[i], eta0_const[i+1]);
	}
	printf("crownrat = \n");
	printf("%30.18f\t", crownrat[0]);
	printf("%30.18f\n", crownrat[1]);
	printf("crown_const = [%.18f, %.18f]\n",
		crown_const[0], crown_const[1]);
*/
}


void tconst_fun( int index, double out[2] )
{
	int i;
	double ttable[10] = 
		{0.0, 0.0, 0.0, 0.1317, 0.27113, 0.41056,
		 0.549999, 0.6045, 0.6978, 0.7891};
	double temp[2], pt[2], temp2[2], zeta[2];
	double tet[12];
	
	/* compute zeta here */
	tet[0] = 2.0;
	tet[1] = 2.0;
	ROUND_DOWN;
	temp[0] = sqrt(tet[0])/5.0;
	ROUND_UP;
	temp[1] = sqrt(tet[1])/5.0;
	I_ATAN( temp, tet );
	temp2[0] = 0.5;
	temp2[1] = 0.5;
	i_div( temp2, tet, zeta );
	for( i=0; i<12; i++ ) {
		tet[i] = 2.0;
	}
	/* compute pt */
	pt[0] = rough_min_gma( tet );
	pt[1] = rough_max_gma( tet );
	
	ROUND_DOWN;
	temp[0] = 4.0*PI_LO*zeta[0] - 8.0;
	ROUND_UP;
	temp[1] = 4.0*PI_HI*zeta[1] - 8.0;
	
	i_mult( temp, pt, temp2 );
	if( index < 10 ) {
		if( index < 4 ) {
			printf("t_%d is not defined.\n", index);
			out[0] = 0.0;
			out[1] = 0.0;
		} else {
			i_recognize( ttable[index-1], out );
		}
	} else {
		out[0] = temp2[0];
		out[1] = temp2[1];
	}
}


void sconst_fun( int index, double out[2] )
{
	double stable[9] = 
		{0.0, 0.0, 0.0, 0.0, -0.05704, -0.11408,
		 -0.17112, -0.22816, -0.1972};

	if( index < 10 ) {
		if( index < 5 ) {
			printf("s_%d is not defined.\n", index);
			out[0] = 0.0;
			out[1] = 0.0;
		} else {
			i_recognize( stable[index-1], out );
		}
	} else {
		out[0] = 0.0;
		out[1] = 0.0;
	}
	
}


void tomDfun( int n, int k, double out[2] )
{
	double temp[2], temp2[2], temp3[2], val;
	
	if( k > n ) {
		printf("k > n in tomDfun.\n");
		out[0] = 0.0;
		out[1] = 0.0;
		return;
	}
	if( n + k < 4 ) {
		printf("n + k < 4 in tomDfun.\n");
		out[0] = 0.0;
		out[1] = 0.0;
		return;
	}
	
	val = -0.06585;
	i_recognize( val, temp3 );
	temp2[0] = (double) k;
	temp2[1] = (double) k;
	i_mult( temp2, temp3, temp );
	
	tconst_fun( n + k, temp2 );
	i_add( temp, temp2, out );
}


void tomZfun( int n, int k, double out[2] )
{
	double temp[2], temp2[2], temp3[2], epsilon[2], val;
	
	val = 0.00005;
	i_recognize( val, epsilon );
	if( n == 3 && k == 1 ) {
		out[0] = epsilon[0];
		out[1] = epsilon[1];
	} else {
		sconst_fun( n + k, temp );
		temp2[0] = (double) k;
		temp2[1] = (double) k;
		i_mult( epsilon, temp2, temp3 );
		i_sub( temp, temp3, out );
	}
}


/*
2.123 2.29102 2.341248 2.011234 2.023309 2.03982134
2.0123 2.1765 2.0987 2.1567 2.0324 2.2187
2.59415 2.0 2.0 2.59415 2.0 2.0
*/

/*
2.695999999999999730 2.828427124746190291
2.449999999999999734 2.510000000000000231
2.0  2.510000000000000231
2.769999999999999574 2.769999999999999574
2.0  2.510000000000000231
2.449999999999999734 2.510000000000000231

2.695999999999999730 2.696129323364009789
2.464999999999999858 2.465058593749999805
2.0 2.000498046875000213
2.769999999999999574 2.769999999999999574
2.0 2.000498046875000213
2.464999999999999858 2.465058593749999805
*/
