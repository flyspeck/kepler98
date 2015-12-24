/* i_voronoi.c,  by Samuel Ferguson, (c) 1998,
containing routines for computing bounds on 
truncated voronoi functions.  */

#include "system_headers.h"
#include "i_sphere.h"
#include "interval.h"
#include "i_bounds.h"
#include "i_voronoi.h"
#include "i_appendix.h"
#include "second_partials.h"
#include "macros.h"


/* #pragma fenv_access	*/	/* Tell the compiler not to do 
anything too clever.  */

#define ACUTEFACES	1

#define DEBUG				0


/* External variables */

extern double i_doct_const[2];
extern double t0_val[2];
extern double phi0_val[2];
extern double zetapt_val[2];
extern double sol_coeff[2];

/* Global variables */

double super_coeff[2];

/* Routines */



/*
quo =

1/6*((a*(-a^2 + b^2)*Sqrt[-b + c]*Sqrt[b + c])/
     (Sqrt[-a + b]*Sqrt[a + b]) - 
    4*c^3*ArcTan[(Sqrt[-a + b]*Sqrt[-b + c])/
       (Sqrt[a + b]*Sqrt[b + c])] + 
    (-a + c)*(-a^2 - a*c + 2*c^2)*
     ArcTan[(Sqrt[-b + c]*Sqrt[b + c])/
       (Sqrt[-a + b]*Sqrt[a + b])])

which becomes

( a*sqrt((y-x)*(z-y)) - 
4*c*z*atan( sqrt( (y-x)*(z-y)/((a+b)*(b+c)) ) ) +
(c-a)*(c-a)*(a+2*c)*atan( sqrt( (z-y)/(y-x) ) ) )/6

quoa = 

(Sqrt[-b + c]*(b + c)*(-3*a^2 + b^2 + 2*c^2) + 
    3*Sqrt[-a + b]*Sqrt[a + b]*Sqrt[b + c]*(a^2 - c^2)*
     ArcTan[(Sqrt[-b + c]*Sqrt[b + c])/
       (Sqrt[-a + b]*Sqrt[a + b])])/
  (6*Sqrt[-a + b]*Sqrt[a + b]*Sqrt[b + c])

which becomes

sqrt((x-y)/(y-x))*(-3*x+y+2*z) - 
3*(z-x)*atan( sqrt((z-y)/(y-x)) )

quob = 

-((a*(-b + c)^(3/2)*(b + c)^(3/2))/
    (3*b*Sqrt[-a + b]*Sqrt[a + b]))

which becomes

-(a*(z-y)*sqrt((z-y)/(y-x))/(3*b))

quoaa = 

-((a*(-Sqrt[-b + c]*(b + c)*(3*a^2 - 4*b^2 + c^2) + 
        3*Sqrt[-a + b]*Sqrt[a + b]*(a^2 - b^2)*
         Sqrt[b + c]*ArcTan[(Sqrt[-b + c]*Sqrt[b + c])/
           (Sqrt[-a + b]*Sqrt[a + b])]))/
    (3*(-a + b)^(3/2)*(a + b)^(3/2)*Sqrt[b + c]))

which becomes

a*((3*x-4*y+z)*sqrt((z-y)/(y-x))/(3*(y-x)) +
atan(sqrt((z-y)/(y-x))))

quoab =

-((b*(-b + c)^(3/2)*(b + c)^(3/2))/
    (3*(-a + b)^(3/2)*(a + b)^(3/2)))
    
which becomes

-b*((z-y)/(y-x))*sqrt((z-y)/(y-x))/3

quobb =

-((a*Sqrt[-b + c]*Sqrt[b + c]*
      (a^2*(2*b^2 + c^2) - b^2*(b^2 + 2*c^2)))/
    (3*b^2*(-a + b)^(3/2)*(a + b)^(3/2)))

which becomes

a*(sqrt( (z-y)/(y-x) )*(y*(y+2*z)-x*(2*y+z)))/(3*y*(y-x))

or

a*(sqrt( (z-y)/(y-x) )*(2*y*(x-z)+y*y-x*z))/(3*y*(y-x))


Note that we get good monotonicity results:

quo >= 0

quo_a >= 0
quo_b <= 0

quo_aa >= 0
quo_ab <= 0
quo_bb >= 0

quo_aab <= 0

*/


/* secpars order:  aa, ab, bb */
void i_quoin_all( double abc[6], double xyz[6],
	double out_quo[2], double out_quo_pars[4], 
	double out_quo_secpars[6] )
{
	int i;
	double temp, num[2], den[2], p1[2], p2[2], p3[2];
	double y_x[2], z_y[2], z_x[2], sqrty_xz_y[2];
	double sqrtz_ydy_x[2], quot[2], atansqrt[2];
	double sum[2];
	
	/* all zero if b_hi <= a_lo or c_hi <= b_lo */
	if( abc[3] <= abc[0] || abc[5] <= abc[2] ) {	
		for( i=0; i<2; i++ ) {
			out_quo[i] = 0.0;
			out_quo_pars[i] = 0.0;
			out_quo_secpars[i] = 0.0;
		}
		for( i=2; i<4; i++ ) {
			out_quo_pars[i] = 0.0;
			out_quo_secpars[i] = 0.0;
		}
		for( i=4; i<6; i++ ) {
			out_quo_secpars[i] = 0.0;
		}
		return;
	}
	ROUND_DOWN;
	temp = xyz[2] - xyz[1];
	if( temp < 0.0 )
		temp = 0.0;
	y_x[0] = temp;
	temp = xyz[4] - xyz[3];
	if( temp < 0.0 )
		temp = 0.0;
	z_y[0] = temp;
	temp = xyz[4] - xyz[1];
	if( temp < 0.0 )
		temp = 0.0;
	z_x[0] = temp;
	ROUND_UP;
	y_x[1] = xyz[3] - xyz[0];
	z_y[1] = xyz[5] - xyz[2];
	z_x[1] = xyz[5] - xyz[0];
	sqrty_xz_y[1] = sqrt( y_x[1]*z_y[1] );
	sqrtz_ydy_x[1] = sqrt( z_y[1]/y_x[0] );
	ROUND_DOWN;
	sqrty_xz_y[0] = sqrt( y_x[0]*z_y[0] );
	sqrtz_ydy_x[0] = sqrt( z_y[0]/y_x[1] );
	i_atan( sqrtz_ydy_x, atansqrt );
	if( atansqrt[0] < 0.0 )
		atansqrt[0] = 0.0;
	
	/* Begin with quo */
	ROUND_DOWN;
	den[0] = (abc[0] + abc[2])*(abc[2] + abc[4]);
	ROUND_UP;
	den[1] = (abc[1] + abc[3])*(abc[3] + abc[5]);
	quot[1] = sqrty_xz_y[1]/den[0];
	ROUND_DOWN;
	quot[0] = sqrty_xz_y[0]/den[1];
	i_atan( quot, num );
	if( num[0] < 0.0 )
		num[0] = 0.0;
	ROUND_DOWN;
	p2[0] = 4.0*abc[4]*xyz[4]*num[0];
	p1[0] = abc[0]*sqrty_xz_y[0];
	temp = abc[4] - abc[1];
	if( temp < 0.0 )
		temp = 0.0;
	p3[0] = temp*temp*(abc[0] + 2.0*abc[4])*atansqrt[0];
	ROUND_UP;
	p2[1] = 4.0*abc[5]*xyz[5]*num[1];
	p1[1] = abc[1]*sqrty_xz_y[1];
	temp = abc[5] - abc[0];
	p3[1] = temp*temp*(abc[1] + 2.0*abc[5])*atansqrt[1];
	out_quo[1] = (p1[1] - p2[0] + p3[1])*ONE_6_HI;
	ROUND_DOWN;
	temp = (p1[0] - p2[1] + p3[0])*ONE_6_LO;
	if( temp < 0.0 )
		temp = 0.0;
	out_quo[0] = temp;
	
	/* Now do quo_a */
	den[0] = 3.0*xyz[0];
	ROUND_UP;
	den[1] = 3.0*xyz[1];
	num[1] = xyz[3] + 2.0*xyz[5] - den[0];
	ROUND_DOWN;
	temp = xyz[2] + 2.0*xyz[4] - den[1];
	if( temp < 0.0 )
		temp = 0.0;
	num[0] = temp;
	p1[0] = sqrtz_ydy_x[0]*num[0];
	p2[0] = 3.0*z_x[0]*atansqrt[0];
	ROUND_UP;
	p1[1] = sqrtz_ydy_x[1]*num[1];
	p2[1] = 3.0*z_x[1]*atansqrt[1];
	temp = p1[1] - p2[0];
	if( temp > 0.0 )
		out_quo_pars[1] = temp*ONE_6_HI;
	else
		out_quo_pars[1] = temp*ONE_6_LO;
	ROUND_DOWN;
	temp = p1[0] - p2[1];
	if( temp > 0.0 )
		out_quo_pars[0] = temp*ONE_6_HI;
	else
		out_quo_pars[0] = temp*ONE_6_LO;
	
	/* Now do quo_b */
	num[0] = abc[0]*z_y[0]*sqrtz_ydy_x[0];
	den[0] = 3.0*abc[2];
	ROUND_UP;
	num[1] = abc[1]*z_y[1]*sqrtz_ydy_x[1];
	den[1] = 3.0*abc[3];
	p1[1] = num[1]/den[0];
	out_quo_pars[2] = -p1[1];
	ROUND_DOWN;
	p1[0] = num[0]/den[1];
	out_quo_pars[3] = -p1[0];
	
	/* Now do quo_aa */
	p3[0] = 3.0*xyz[0] - 4.0*xyz[3] + xyz[4];
	p2[0] = 3.0*y_x[0];
	ROUND_UP;
	p3[1] = 3.0*xyz[1] - 4.0*xyz[2] + xyz[5];
	p2[1] = 3.0*y_x[1];
	num[1] = sqrtz_ydy_x[1]/p2[0];
	ROUND_DOWN;
	num[0] = sqrtz_ydy_x[0]/p2[1];
	i_mult( num, p3, p1 );
	ROUND_DOWN;
	sum[0] = p1[0] + atansqrt[0];
	ROUND_UP;
	sum[1] = p1[1] + atansqrt[1];
	i_mult( abc, sum, out_quo_secpars );
	
	/* Now do quo_ab */
	ROUND_DOWN;
	p1[0] = abc[2]*z_y[0]*sqrtz_ydy_x[0]/p2[1];
	ROUND_UP;
	p1[1] = abc[3]*z_y[1]*sqrtz_ydy_x[1]/p2[0];
	out_quo_secpars[2] = -p1[1];
	out_quo_secpars[3] = -p1[0];
	
	/* Now do quo_bb */
	p1[1] = xyz[3]*(xyz[3] + 2.0*xyz[5]);
	p2[1] = xyz[1]*(2.0*xyz[3] + xyz[5]);
	den[1] = 3.0*xyz[3]*y_x[1];
	ROUND_DOWN;
	p1[0] = xyz[2]*(xyz[2] + 2.0*xyz[4]);
	p2[0] = xyz[0]*(2.0*xyz[2] + xyz[4]);
	den[0] = 3.0*xyz[2]*y_x[0];
	temp = p1[0] - p2[1];
	if( temp < 0.0 )
		temp = 0.0;
	p3[0] = temp;
	num[0] = abc[0]*sqrtz_ydy_x[0]*p3[0];
	out_quo_secpars[4] = num[0]/den[1];
	ROUND_UP;
	p3[1] = p1[1] - p2[0];
	num[1] = abc[1]*sqrtz_ydy_x[1]*p3[1];
	out_quo_secpars[5] = num[1]/den[0];
}


void i_quoin( double abc[6], double xyz[6],
	double out_quo[2] )
{
	int i;
	double temp, num[2], den[2], p1[2], p2[2], p3[2];
	double y_x[2], z_y[2], z_x[2], sqrty_xz_y[2];
	double sqrtz_ydy_x[2], quot[2], atansqrt[2];
	
	/* all zero if b_hi <= a_lo or c_hi <= b_lo */
	if( abc[3] <= abc[0] || abc[5] <= abc[2] ) {	
		for( i=0; i<2; i++ ) {
			out_quo[i] = 0.0;
		}
		return;
	}
	ROUND_DOWN;
	temp = xyz[2] - xyz[1];
	if( temp < 0.0 )
		temp = 0.0;
	y_x[0] = temp;
	temp = xyz[4] - xyz[3];
	if( temp < 0.0 )
		temp = 0.0;
	z_y[0] = temp;
	temp = xyz[4] - xyz[1];
	if( temp < 0.0 )
		temp = 0.0;
	z_x[0] = temp;
	ROUND_UP;
	y_x[1] = xyz[3] - xyz[0];
	z_y[1] = xyz[5] - xyz[2];
	z_x[1] = xyz[5] - xyz[0];
	sqrty_xz_y[1] = sqrt( y_x[1]*z_y[1] );
	sqrtz_ydy_x[1] = sqrt( z_y[1]/y_x[0] );
	ROUND_DOWN;
	sqrty_xz_y[0] = sqrt( y_x[0]*z_y[0] );
	sqrtz_ydy_x[0] = sqrt( z_y[0]/y_x[1] );
	i_atan( sqrtz_ydy_x, atansqrt );
	if( atansqrt[0] < 0.0 )
		atansqrt[0] = 0.0;
	
	/* Begin with quo */
	ROUND_DOWN;
	den[0] = (abc[0] + abc[2])*(abc[2] + abc[4]);
	ROUND_UP;
	den[1] = (abc[1] + abc[3])*(abc[3] + abc[5]);
	quot[1] = sqrty_xz_y[1]/den[0];
	ROUND_DOWN;
	quot[0] = sqrty_xz_y[0]/den[1];
	i_atan( quot, num );
	if( num[0] < 0.0 )
		num[0] = 0.0;
	ROUND_DOWN;
	p2[0] = 4.0*abc[4]*xyz[4]*num[0];
	p1[0] = abc[0]*sqrty_xz_y[0];
	temp = abc[4] - abc[1];
	if( temp < 0.0 )
		temp = 0.0;
	p3[0] = temp*temp*(abc[0] + 2.0*abc[4])*atansqrt[0];
	ROUND_UP;
	p2[1] = 4.0*abc[5]*xyz[5]*num[1];
	p1[1] = abc[1]*sqrty_xz_y[1];
	temp = abc[5] - abc[0];
	p3[1] = temp*temp*(abc[1] + 2.0*abc[5])*atansqrt[1];
	out_quo[1] = (p1[1] - p2[0] + p3[1])*ONE_6_HI;
	ROUND_DOWN;
	temp = (p1[0] - p2[1] + p3[0])*ONE_6_LO;
	if( temp < 0.0 )
		temp = 0.0;
	out_quo[0] = temp;
}


void i_quoin_pars( double abc[6], double xyz[6],
	double out_quo[2], double out_quo_pars[4] )
{
	int i;
	double temp, num[2], den[2], p1[2], p2[2], p3[2];
	double y_x[2], z_y[2], z_x[2], sqrty_xz_y[2];
	double sqrtz_ydy_x[2], quot[2], atansqrt[2];
	
	/* all zero if b_hi <= a_lo or c_hi <= b_lo */
	if( abc[3] <= abc[0] || abc[5] <= abc[2] ) {	
		for( i=0; i<2; i++ ) {
			out_quo[i] = 0.0;
			out_quo_pars[i] = 0.0;
		}
		for( i=2; i<4; i++ ) {
			out_quo_pars[i] = 0.0;
		}
		return;
	}
	ROUND_DOWN;
	temp = xyz[2] - xyz[1];
	if( temp < 0.0 )
		temp = 0.0;
	y_x[0] = temp;
	temp = xyz[4] - xyz[3];
	if( temp < 0.0 )
		temp = 0.0;
	z_y[0] = temp;
	temp = xyz[4] - xyz[1];
	if( temp < 0.0 )
		temp = 0.0;
	z_x[0] = temp;
	ROUND_UP;
	y_x[1] = xyz[3] - xyz[0];
	z_y[1] = xyz[5] - xyz[2];
	z_x[1] = xyz[5] - xyz[0];
	sqrty_xz_y[1] = sqrt( y_x[1]*z_y[1] );
	sqrtz_ydy_x[1] = sqrt( z_y[1]/y_x[0] );
	ROUND_DOWN;
	sqrty_xz_y[0] = sqrt( y_x[0]*z_y[0] );
	sqrtz_ydy_x[0] = sqrt( z_y[0]/y_x[1] );
	i_atan( sqrtz_ydy_x, atansqrt );
	if( atansqrt[0] < 0.0 )
		atansqrt[0] = 0.0;
	
	/* Begin with quo */
	ROUND_DOWN;
	den[0] = (abc[0] + abc[2])*(abc[2] + abc[4]);
	ROUND_UP;
	den[1] = (abc[1] + abc[3])*(abc[3] + abc[5]);
	quot[1] = sqrty_xz_y[1]/den[0];
	ROUND_DOWN;
	quot[0] = sqrty_xz_y[0]/den[1];
	i_atan( quot, num );
	if( num[0] < 0.0 )
		num[0] = 0.0;
	ROUND_DOWN;
	p2[0] = 4.0*abc[4]*xyz[4]*num[0];
	p1[0] = abc[0]*sqrty_xz_y[0];
	temp = abc[4] - abc[1];
	if( temp < 0.0 )
		temp = 0.0;
	p3[0] = temp*temp*(abc[0] + 2.0*abc[4])*atansqrt[0];
	ROUND_UP;
	p2[1] = 4.0*abc[5]*xyz[5]*num[1];
	p1[1] = abc[1]*sqrty_xz_y[1];
	temp = abc[5] - abc[0];
	p3[1] = temp*temp*(abc[1] + 2.0*abc[5])*atansqrt[1];
	out_quo[1] = (p1[1] - p2[0] + p3[1])*ONE_6_HI;
	ROUND_DOWN;
	temp = (p1[0] - p2[1] + p3[0])*ONE_6_LO;
	if( temp < 0.0 )
		temp = 0.0;
	out_quo[0] = temp;
	
	/* Now do quo_a */
	den[0] = 3.0*xyz[0];
	ROUND_UP;
	den[1] = 3.0*xyz[1];
	num[1] = xyz[3] + 2.0*xyz[5] - den[0];
	ROUND_DOWN;
	temp = xyz[2] + 2.0*xyz[4] - den[1];
	if( temp < 0.0 )
		temp = 0.0;
	num[0] = temp;
	p1[0] = sqrtz_ydy_x[0]*num[0];
	p2[0] = 3.0*z_x[0]*atansqrt[0];
	ROUND_UP;
	p1[1] = sqrtz_ydy_x[1]*num[1];
	p2[1] = 3.0*z_x[1]*atansqrt[1];
	temp = p1[1] - p2[0];
	if( temp > 0.0 )
		out_quo_pars[1] = temp*ONE_6_HI;
	else
		out_quo_pars[1] = temp*ONE_6_LO;
	ROUND_DOWN;
	temp = p1[0] - p2[1];
	if( temp > 0.0 )
		out_quo_pars[0] = temp*ONE_6_HI;
	else
		out_quo_pars[0] = temp*ONE_6_LO;
	
	/* Now do quo_b */
	num[0] = abc[0]*z_y[0]*sqrtz_ydy_x[0];
	den[0] = 3.0*abc[2];
	ROUND_UP;
	num[1] = abc[1]*z_y[1]*sqrtz_ydy_x[1];
	den[1] = 3.0*abc[3];
	p1[1] = num[1]/den[0];
	out_quo_pars[2] = -p1[1];
	ROUND_DOWN;
	p1[0] = num[0]/den[1];
	out_quo_pars[3] = -p1[0];
}


/* beta = (2*t + h)*(t - h)^2 */
/* beta is decreasing in h, increasing in t */
void i_vor_beta( double hval[2], double out[2] )
{
	double h, t, temp;
	
	if( hval[0] < t0_val[1] ) {
		/* First compute min. */
		h = hval[1];
		t = t0_val[0];
		ROUND_DOWN;
		temp = t - h;
		if( temp < 0.0 )
			temp = 0.0;
		out[0] = temp*temp*(2.0*t + h);
		/* Now compute max. */
		h = hval[0];
		t = t0_val[1];
		ROUND_UP;
		temp = t - h;
		out[1] = temp*temp*(2.0*t + h);
	} else {
		out[0] = 0.0;
		out[1] = 0.0;
	}
}


/* beta_h = -3*(t + h)*(t - h) */
/* beta_h is increasing in h, decreasing in t */
void i_vor_beta_h( double hval[2], double out[2] )
{
	double h, t, temp, fred;
	
	if( hval[0] < t0_val[1] ) {
		/* First compute max. */
		h = hval[1];
		t = t0_val[0];
		ROUND_DOWN;
		temp = t - h;
		if( temp < 0.0 )
			temp = 0.0;
		fred = 3.0*temp*(t + h);
		out[1] = -fred;
		/* Now compute min. */
		h = hval[0];
		t = t0_val[1];
		ROUND_UP;
		fred = 3.0*(t - h)*(t + h);
		out[0] = -fred;
	} else {
		out[0] = 0.0;
		out[1] = 0.0;
	}
}


/* beta_hh = 6*h */
/* beta_hh is increasing in h, and discontinuous at h=t0 */
void i_vor_beta_hh( double hval[2], double out[2] )
{	
	if( hval[0] < t0_val[1] ) {
		/* First compute min. */
		if( hval[1] < t0_val[0] ) {
		ROUND_DOWN;
		out[0] = 6.0*hval[0];
		} else {
			out[0] = 0.0;
		}
		/* Now compute max. */
		ROUND_UP;
		out[1] = 6.0*hval[1];
	} else {
		out[0] = 0.0;
		out[1] = 0.0;
	}
}


void i_vor_beta_x( double hval[2], double oohval[2],
	double out[2] )
{
	double betah[2], temp[2];
	
	i_vor_beta_h( hval, betah );
	i_mult( oohval, betah, temp );
	out[0] = 0.125*temp[0];
	out[1] = 0.125*temp[1];
}


void i_vor_beta_xx( double hval[2], double oohval[2],
	double beta_x[2], double out[2] )
{
	double betahh[2], temp[2], temp2[2];
	
	i_vor_beta_hh( hval, betahh );
	ROUND_DOWN;
	temp[0] = 0.125*betahh[0] - beta_x[1];
	temp2[0] = 0.125*oohval[0]*oohval[0];
	ROUND_UP;
	temp[1] = 0.125*betahh[1] - beta_x[0];
	temp2[1] = 0.125*oohval[1]*oohval[1];
	i_mult( temp2, temp, out );
}


void i_vor0( double y[12], double x[12], 
	double sqrtdelta[2], double out[2] )
{
	int facelist[3][3] = {{0, 2, 10},{0, 4, 8},{2, 4, 6}};
	int i, j, k, n, index;
	double facex[6], abc[6], xyz[6];
	double quo[2], quosum[2], solmult[2];
	double h[2], wedgesum[2], beta_val[2], dih_val[2];
	
	i_solid( y, sqrtdelta, quo );
	i_mult( quo, sol_coeff, solmult );
/*	
	ROUND_NEAR;
	printf("sol_coeff = [%.18f, %.18f]\n", 
		sol_coeff[0], sol_coeff[1]);
	printf("solmult = [%.18f, %.18f]\n", solmult[0], solmult[1]);
*/
	
	/* First do quoins */
	abc[4] = t0_val[0];
	abc[5] = t0_val[1];
	ROUND_UP;
	xyz[5] = t0_val[1]*t0_val[1];
	ROUND_DOWN;
	xyz[4] = t0_val[0]*t0_val[0];
	
	quosum[0] = 0.0;
	quosum[1] = 0.0;
	for( i=0; i<3; i++ ) {
		for( j=0; j<3; j++ ) {
			index = facelist[i][j];
			k = 2*j;
			facex[k] = x[index];
			facex[k+1] = x[index+1];
		}
		i_crad3x2( facex, xyz + 2 );
/*		
		ROUND_NEAR;
		printf("crad2(%d %d %d) = [%.18f, %.18f]\n",
			facelist[i][0]/2+1,facelist[i][1]/2+1,
			facelist[i][2]/2+1,xyz[2],xyz[3]);
*/			
		i_sqrt( xyz + 2, abc + 2 );
		/* Check to see if b < c */
		if( abc[2] < t0_val[1] ) {
			for( k=0; k<2; k++ ) {
				index = facelist[i][k];
				n = index + 1;
				abc[0] = 0.5*y[index];
				abc[1] = 0.5*y[n];
				/* Check to see if a < b */
				if( abc[0] < abc[3] ) {
					xyz[0] = 0.25*x[index];
					xyz[1] = 0.25*x[n];
					i_quoin( abc, xyz, quo );
/*					
					ROUND_NEAR;
					printf("quoin = [%.18f, %.18f]\n", quo[0], quo[1]);
*/					
					ROUND_DOWN;
					quosum[0] += quo[0];
					ROUND_UP;
					quosum[1] += quo[1];
				}
			}
		}
	}
/*	
	ROUND_NEAR;
	printf("quosum = [%.18f, %.18f]\n", quosum[0], quosum[1]);
*/	
	
	/* Now do wedge terms */
	wedgesum[0] = 0.0;
	wedgesum[1] = 0.0;
	for( i=0; i<3; i++ ) {
		j = 2*i;
		h[0] = 0.5*y[j];
		h[1] = 0.5*y[j+1];
		/* Check to see if h_i < t_0 */
		if( h[0] < t0_val[1] ) {
			i_vor_beta( h, beta_val );
			i_index_dih( i, x, dih_val );
			if( dih_val[0] < 0.0 )
				dih_val[0] = 0.0;
			ROUND_DOWN;
			wedgesum[0] += beta_val[0]*dih_val[0];
			ROUND_UP;
			wedgesum[1] += beta_val[1]*dih_val[1];
		}
	}
/*	
	ROUND_NEAR;
	printf("wedgesum = [%.18f, %.18f]\n", wedgesum[0]/6, wedgesum[1]/6);
*/	
	ROUND_DOWN;
	h[0] = ONE_6_LO*wedgesum[0] - quosum[1];
	ROUND_UP;
	h[1] = ONE_6_HI*wedgesum[1] - quosum[0];
	quo[0] = DOCT_LO;
	quo[1] = DOCT_HI;
	i_mult( quo, h, quosum );
	ROUND_DOWN;
	out[0] = solmult[0] + 4.0*quosum[0];
	ROUND_UP;
	out[1] = solmult[1] + 4.0*quosum[1];
}


void i_vor0_partials( double y[12], double x[12], 
	double out_partials[12] )
{
	int facelist[3][3] = {{0, 2, 10},{0, 4, 8},{2, 4, 6}};
	int i, j, k, n, index;
	double face[6], facex[6], abc[6], xyz[6];
	double quo_pars[4], quosum_pars[12];
	double solmult_pars[12], partials[12];
	double h[2], wedgesum_pars[12], beta_val[2], dih_val[2];
	double delta[2], sqrtdelta[2], delta_part[12];
	double crad_pars[6], quob_sum[2], quo[2];
	double etaparquob[6], beta_h_val[2], dih_pars[12];
	double temp[2];
		
	i_bigdelta_partials_best( x, delta, delta_part );
	i_sqrt( delta, sqrtdelta );
	s_solid_partials( y, delta, sqrtdelta, delta_part, 
		partials );
	for( i=0; i<12; i+=2 ) {
		i_mult( sol_coeff, partials + i, solmult_pars + i );
	}

/*
	ROUND_NEAR;
	printf("solmult_pars = \n");
	for( i=0; i<12; i+=2 ) {
		printf("[%.18f, %.18f]\n", solmult_pars[i], solmult_pars[i+1]);
	}
*/

	/* First do quoins */
	abc[4] = t0_val[0];
	abc[5] = t0_val[1];
	ROUND_UP;
	xyz[5] = t0_val[1]*t0_val[1];
	ROUND_DOWN;
	xyz[4] = t0_val[0]*t0_val[0];
	
	for( i=0; i<12; i++ )
		quosum_pars[i] = 0.0;
	
	for( i=0; i<3; i++ ) {
		for( j=0; j<3; j++ ) {
			index = facelist[i][j];
			k = 2*j;
			facex[k] = x[index];
			face[k]  = y[index];
			facex[k+1] = x[index+1];
			face[k+1]  = y[index+1];
		}
		i_crad3x2( facex, xyz + 2 );
		crad3len_pars( face, facex, crad_pars );
/*
		ROUND_NEAR;
		printf("crad3len_pars = \n");
		for( j=0; j<6; j+=2 ) {
			printf("[%.18f, %.18f]\n", crad_pars[j], crad_pars[j+1]);
		}
*/
		i_sqrt( xyz + 2, abc + 2 );
		quob_sum[0] = 0.0;
		quob_sum[1] = 0.0;
		/* Check to see if b < c */
		if( abc[2] < t0_val[1] ) {
			for( k=0; k<2; k++ ) {
				index = facelist[i][k];
				n = index + 1;
				abc[0] = 0.5*y[index];
				abc[1] = 0.5*y[n];
/*				
				ROUND_NEAR;
				printf("abc = \n");
				for( j=0; j<6; j+=2 ) {
					printf("[%.18f, %.18f]\n", abc[j], abc[j+1]);
				}
*/				
				/* Check to see if a < b */
				if( abc[0] < abc[3] ) {
					xyz[0] = 0.25*x[index];
					xyz[1] = 0.25*x[n];
					i_quoin_pars( abc, xyz, quo, quo_pars );
/*					
					ROUND_NEAR;
					printf("quo_pars = \n");
					for( j=0; j<4; j+=2 ) {
						printf("[%.18f, %.18f]\n", quo_pars[j], quo_pars[j+1]);
					}
*/					
					ROUND_DOWN;
					quob_sum[0] += quo_pars[2];
					ROUND_UP;
					quob_sum[1] += quo_pars[3];
					quosum_pars[index] += 0.5*quo_pars[0];
					quosum_pars[n    ] += 0.5*quo_pars[1];
				}
			}
		}
		for( k=0; k<6; k+=2 ) {
			i_mult( quob_sum, crad_pars + k, etaparquob + k );
		}
		ROUND_DOWN;
		for( j=0; j<3; j++ ) {
			index = facelist[i][j];
			k = 2*j;
			quosum_pars[index] += etaparquob[k];
		}
		ROUND_UP;
		for( j=0; j<3; j++ ) {
			index = facelist[i][j] + 1;
			k = 2*j + 1;
			quosum_pars[index] += etaparquob[k];
		}
/*		
		if( i==0 ) {
		ROUND_NEAR;
			printf("orig quosum_pars = \n");
			for( k=0; k<12; k+=2 ) {
				printf("[%.18f, %.18f]\n", quosum_pars[k], quosum_pars[k+1]);
			}
		}
*/
	}
/*
	ROUND_NEAR;
	printf("quosum_pars = \n");
	for( k=0; k<12; k+=2 ) {
		printf("[%.18f, %.18f]\n", quosum_pars[k], quosum_pars[k+1]);
	}
*/
	
	/* Now do wedge terms */
	for( i=0; i<12; i++ ) {
		wedgesum_pars[i] = 0.0;
	}
	for( i=0; i<3; i++ ) {
		j = 2*i;
		h[0] = 0.5*y[j];
		h[1] = 0.5*y[j+1];
		/* Check to see if h_i < t_0 */
		if( h[0] < t0_val[1] ) {
			i_vor_beta( h, beta_val );
			i_vor_beta_h( h, beta_h_val );
			i_index_dih( i, x, dih_val );
			if( dih_val[0] < 0.0 )
				dih_val[0] = 0.0;
			i_index_dih_partials( i, y, x, dih_pars );
/*			
			ROUND_NEAR;
			printf("dih[%d]_partials = \n", i);
			for( k=0; k<12; k+=2 ) {
				printf("[%.18f, %.18f]\n", dih_pars[k], 
					dih_pars[k+1]);
			}
			printf("beta_val = [%.18f, %.18f]\n", beta_val[0],
				beta_val[1]);
			printf("beta_h_val = [%.18f, %.18f]\n", beta_h_val[0],
				beta_h_val[1]);
			printf("dih_val = [%.18f, %.18f]\n", dih_val[0],
				dih_val[1]);
*/			
			for( k=0; k<12; k+=2 ) {
				i_mult( beta_val, dih_pars + k, partials + k );
			}
			i_mult( beta_h_val, dih_val, temp );
			ROUND_DOWN;
			for( k=0; k<12; k+=2 ) {
				wedgesum_pars[k] += partials[k];
			}
			wedgesum_pars[j] += 0.5*temp[0];
			ROUND_UP;
			for( k=1; k<12; k+=2 ) {
				wedgesum_pars[k] += partials[k];
			}
			wedgesum_pars[j+1] += 0.5*temp[1];
/*
			ROUND_NEAR;
			printf("temp = [%.18f, %.18f]\n", temp[0],
				temp[1]);
*/
		}
	}
/*	
	ROUND_NEAR;
	printf("wedgesum_pars = \n");
	for( i=0; i<12; i+=2 ) {
		printf("[%.18f, %.18f]\n", wedgesum_pars[i], wedgesum_pars[i+1]);
	}
*/	
	temp[0] = ONE_6_LO;
	temp[1] = ONE_6_HI;
	
	for( i=0; i<12; i+=2 ) {
		i_mult( temp, wedgesum_pars + i, delta_part + i );
	}
	
	ROUND_DOWN;
	for( i=0; i<12; i+=2 ) {
		partials[i] = delta_part[i] - quosum_pars[i+1];
	}
	ROUND_UP;
	for( i=1; i<12; i+=2 ) {
		partials[i] = delta_part[i] - quosum_pars[i-1];
	}
	
	temp[0] = DOCT_LO;
	temp[1] = DOCT_HI;
	for( i=0; i<12; i+=2 ) {
		i_mult( temp, partials + i, delta_part + i );
	}

	ROUND_DOWN;
	for( i=0; i<12; i+=2 ) {
		out_partials[i] = solmult_pars[i] + 4.0*delta_part[i];
	}
	ROUND_UP;
	for( i=1; i<12; i+=2 ) {
		out_partials[i] = solmult_pars[i] + 4.0*delta_part[i];
	}
}


void clear_vor0_partials( double y[12], double x[12], 
	double out_partials[12] )
{
	int facelist[3][3] = {{0, 2, 10},{0, 4, 8},{2, 4, 6}};
	int i, j, k, n, index;
	double face[6], facex[6], abc[6], xyz[6];
	double quo_pars[4], quosum_pars[12];
	double solmult_pars[12], partials[12];
	double h[2], wedgesum_pars[12], beta_val[2], dih_val[2];
	double delta32[2], delta_part[12];
	double crad_pars[6], quob_sum[2], quo[2];
	double etaparquob[6], beta_h_val[2], dih_pars[12];
	double clear_quosum_pars[12];
	double temp[2];
	
	i_bigdelta_best( x, face );
	if( face[0] < 0.0 )
		face[0] = 0.0;
	i_sqrt( face, facex );
	ROUND_DOWN;
	delta32[0] = face[0]*facex[0];
	ROUND_UP;
	delta32[1] = face[1]*facex[1];
	
	clear_sol_i( y, x, partials );
	for( i=0; i<12; i+=2 ) {
		i_mult( sol_coeff, partials + i, solmult_pars + i );
	}
	
/*
	ROUND_NEAR;
	printf("solmult_pars = \n");
	for( i=0; i<12; i+=2 ) {
		printf("[%.18f, %.18f]\n", solmult_pars[i], solmult_pars[i+1]);
	}
*/

	/* First do quoins */
	abc[4] = t0_val[0];
	abc[5] = t0_val[1];
	ROUND_UP;
	xyz[5] = t0_val[1]*t0_val[1];
	ROUND_DOWN;
	xyz[4] = t0_val[0]*t0_val[0];
	
	for( i=0; i<12; i++ )
		quosum_pars[i] = 0.0;
	
	for( i=0; i<3; i++ ) {
		for( j=0; j<3; j++ ) {
			index = facelist[i][j];
			k = 2*j;
			facex[k] = x[index];
			face[k]  = y[index];
			facex[k+1] = x[index+1];
			face[k+1]  = y[index+1];
		}
#if ACUTEFACES
		s_crad3x2( facex, xyz + 2 );
		s_crad3len_pars( face, facex, crad_pars );
#else
		i_crad3x2( facex, xyz + 2 );
		crad3len_pars( face, facex, crad_pars );
#endif
		i_sqrt( xyz + 2, abc + 2 );
		quob_sum[0] = 0.0;
		quob_sum[1] = 0.0;
		/* Check to see if b < c */
		if( abc[2] < t0_val[1] ) {
			for( k=0; k<2; k++ ) {
				index = facelist[i][k];
				n = index + 1;
				abc[0] = 0.5*y[index];
				abc[1] = 0.5*y[n];
				/* Check to see if a < b */
				if( abc[0] < abc[3] ) {
					xyz[0] = 0.25*x[index];
					xyz[1] = 0.25*x[n];
/*
					ROUND_NEAR;
					printf("face = \n");
					printf("[%.18f, %.18f]\n", face[0], face[1]);
					printf("[%.18f, %.18f]\n", face[2], face[3]);
					printf("[%.18f, %.18f]\n", face[4], face[5]);
					printf("facex = \n");
					printf("[%.18f, %.18f]\n", facex[0], facex[1]);
					printf("[%.18f, %.18f]\n", facex[2], facex[3]);
					printf("[%.18f, %.18f]\n", facex[4], facex[5]);
					printf("abc = \n");
					printf("[%.18f, %.18f]\n", abc[0], abc[1]);
					printf("[%.18f, %.18f]\n", abc[2], abc[3]);
					printf("[%.18f, %.18f]\n", abc[4], abc[5]);
*/
					i_quoin_pars( abc, xyz, quo, quo_pars );
					ROUND_DOWN;
					quob_sum[0] += quo_pars[2];
					ROUND_UP;
					quob_sum[1] += quo_pars[3];
					quosum_pars[index] += 0.5*quo_pars[0];
					quosum_pars[n    ] += 0.5*quo_pars[1];
				}
			}
		}
		if( quob_sum[0] != 0.0 || quob_sum[1] != 0.0 ) {
			for( k=0; k<6; k+=2 ) {
				i_mult( quob_sum, crad_pars + k, etaparquob + k );
			}
		} else {
			for( k=0; k<6; k++ ) {
				etaparquob[k] = 0.0;
			}
		}
		ROUND_DOWN;
		for( j=0; j<3; j++ ) {
			index = facelist[i][j];
			k = 2*j;
			quosum_pars[index] += etaparquob[k];
		}
		ROUND_UP;
		for( j=0; j<3; j++ ) {
			index = facelist[i][j] + 1;
			k = 2*j + 1;
			quosum_pars[index] += etaparquob[k];
		}
	}
	/* clear quosum_pars */
	for( i=0; i<12; i+=2 ) {
		i_mult( delta32, quosum_pars + i, clear_quosum_pars + i );
	}
/*
	ROUND_NEAR;
	printf("quosum_pars = \n");
	for( k=0; k<12; k+=2 ) {
		printf("[%.18f, %.18f]\n", quosum_pars[k], quosum_pars[k+1]);
	}
*/
	
	/* Now do wedge terms */
	for( i=0; i<12; i++ ) {
		wedgesum_pars[i] = 0.0;
	}
	for( i=0; i<3; i++ ) {
		j = 2*i;
		h[0] = 0.5*y[j];
		h[1] = 0.5*y[j+1];
		/* Check to see if h_i < t_0 */
		if( h[0] < t0_val[1] ) {
			i_vor_beta( h, beta_val );
			i_vor_beta_h( h, beta_h_val );
			i_index_dih( i, x, temp );
			if( temp[0] < 0.0 )
				temp[0] = 0.0;
			i_mult( temp, delta32, dih_val );	/* clear dih_val */
			clear_index_dih_partials( i, y, x, dih_pars );
			for( k=0; k<12; k+=2 ) {
				i_mult( beta_val, dih_pars + k, partials + k );
			}
			i_mult( beta_h_val, dih_val, temp );
			ROUND_DOWN;
			for( k=0; k<12; k+=2 ) {
				wedgesum_pars[k] += partials[k];
			}
			wedgesum_pars[j] += 0.5*temp[0];
			ROUND_UP;
			for( k=1; k<12; k+=2 ) {
				wedgesum_pars[k] += partials[k];
			}
			wedgesum_pars[j+1] += 0.5*temp[1];
		}
	}
/*	
	for( i=0; i<12; i+=2 ) {
		i_div( wedgesum_pars + i, delta32, delta_part + i );
	}
	ROUND_NEAR;
	printf("new wedgesum_pars = \n");
	for( i=0; i<12; i+=2 ) {
		printf("[%.18f, %.18f]\n", delta_part[i], delta_part[i+1]);
	}
*/
	
	temp[0] = ONE_6_LO;
	temp[1] = ONE_6_HI;
	
	for( i=0; i<12; i+=2 ) {
		i_mult( temp, wedgesum_pars + i, delta_part + i );
	}
	
	ROUND_DOWN;
	for( i=0; i<12; i+=2 ) {
		partials[i] = delta_part[i] - clear_quosum_pars[i+1];
	}
	ROUND_UP;
	for( i=1; i<12; i+=2 ) {
		partials[i] = delta_part[i] - clear_quosum_pars[i-1];
	}
	
	temp[0] = DOCT_LO;
	temp[1] = DOCT_HI;
	for( i=0; i<12; i+=2 ) {
		i_mult( temp, partials + i, delta_part + i );
	}

	ROUND_DOWN;
	for( i=0; i<12; i+=2 ) {
		out_partials[i] = solmult_pars[i] + 4.0*delta_part[i];
	}
	ROUND_UP;
	for( i=1; i<12; i+=2 ) {
		out_partials[i] = solmult_pars[i] + 4.0*delta_part[i];
	}
}


void clear_vor0_xpars( double y[12], double x[12], 
	double oo2y[12], double out_partials[12] )
{
	int facelist[3][3] = {{0, 2, 10},{0, 4, 8},{2, 4, 6}};
	int i, j, k, n, index;
	double face[6], facex[6], abc[6], xyz[6];
	double quo_pars[4], quosum_pars[12];
	double solmult_pars[12], partials[12];
	double h[2], wedgesum_pars[12], beta_val[2], dih_val[2];
	double delta32[2], delta_part[12];
	double crad_pars[6], quob_sum[2], quo[2];
	double etaparquob[6], beta_h_val[2], dih_pars[12];
	double clear_quosum_pars[12], xpartials[12];
	double temp[2], ooh[2];
	
	i_bigdelta_best( x, face );
	if( face[0] < 0.0 )
		face[0] = 0.0;
	i_sqrt( face, facex );
	ROUND_DOWN;
	delta32[0] = face[0]*facex[0];
	ROUND_UP;
	delta32[1] = face[1]*facex[1];
	
	clear_sol_i( y, x, partials );
	sp_ytoxpars( oo2y, partials, xpartials );
	for( i=0; i<12; i+=2 ) {
		i_mult( sol_coeff, xpartials + i, solmult_pars + i );
	}
	
/*
	ROUND_NEAR;
	printf("solmult_pars = \n");
	for( i=0; i<12; i+=2 ) {
		printf("[%.18f, %.18f]\n", solmult_pars[i], solmult_pars[i+1]);
	}
*/

	/* First do quoins */
	abc[4] = t0_val[0];
	abc[5] = t0_val[1];
	ROUND_UP;
	xyz[5] = t0_val[1]*t0_val[1];
	ROUND_DOWN;
	xyz[4] = t0_val[0]*t0_val[0];
	
	for( i=0; i<12; i++ )
		quosum_pars[i] = 0.0;
	
	for( i=0; i<3; i++ ) {
		for( j=0; j<3; j++ ) {
			index = facelist[i][j];
			k = 2*j;
			facex[k] = x[index];
			face[k]  = y[index];
			facex[k+1] = x[index+1];
			face[k+1]  = y[index+1];
		}
#if ACUTEFACES
		s_crad3x2( facex, xyz + 2 );
		s_crad3len_pars( face, facex, crad_pars );
#else
		i_crad3x2( facex, xyz + 2 );
		crad3len_pars( face, facex, crad_pars );
#endif
		i_sqrt( xyz + 2, abc + 2 );
		quob_sum[0] = 0.0;
		quob_sum[1] = 0.0;
		/* Check to see if b < c */
		if( abc[2] < t0_val[1] ) {
			for( k=0; k<2; k++ ) {
				index = facelist[i][k];
				n = index + 1;
				abc[0] = 0.5*y[index];
				abc[1] = 0.5*y[n];
				/* Check to see if a < b */
				if( abc[0] < abc[3] ) {
					xyz[0] = 0.25*x[index];
					xyz[1] = 0.25*x[n];
/*
					ROUND_NEAR;
					printf("face = \n");
					printf("[%.18f, %.18f]\n", face[0], face[1]);
					printf("[%.18f, %.18f]\n", face[2], face[3]);
					printf("[%.18f, %.18f]\n", face[4], face[5]);
					printf("facex = \n");
					printf("[%.18f, %.18f]\n", facex[0], facex[1]);
					printf("[%.18f, %.18f]\n", facex[2], facex[3]);
					printf("[%.18f, %.18f]\n", facex[4], facex[5]);
					printf("abc = \n");
					printf("[%.18f, %.18f]\n", abc[0], abc[1]);
					printf("[%.18f, %.18f]\n", abc[2], abc[3]);
					printf("[%.18f, %.18f]\n", abc[4], abc[5]);
*/
					i_quoin_pars( abc, xyz, quo, quo_pars );
					ROUND_DOWN;
					quob_sum[0] += quo_pars[2];
					ROUND_UP;
					quob_sum[1] += quo_pars[3];
					quosum_pars[index] += 0.5*quo_pars[0];
					quosum_pars[n    ] += 0.5*quo_pars[1];
				}
			}
		}
		if( quob_sum[0] != 0.0 || quob_sum[1] != 0.0 ) {
			for( k=0; k<6; k+=2 ) {
				i_mult( quob_sum, crad_pars + k, etaparquob + k );
			}
		} else {
			for( k=0; k<6; k++ ) {
				etaparquob[k] = 0.0;
			}
		}
		ROUND_DOWN;
		for( j=0; j<3; j++ ) {
			index = facelist[i][j];
			k = 2*j;
			quosum_pars[index] += etaparquob[k];
		}
		ROUND_UP;
		for( j=0; j<3; j++ ) {
			index = facelist[i][j] + 1;
			k = 2*j + 1;
			quosum_pars[index] += etaparquob[k];
		}
	}
	/* clear quosum_pars */
	sp_ytoxpars( oo2y, quosum_pars, xpartials );
	for( i=0; i<12; i+=2 ) {
		i_mult( delta32, xpartials + i, clear_quosum_pars + i );
	}
/*
	ROUND_NEAR;
	printf("quosum_pars = \n");
	for( k=0; k<12; k+=2 ) {
		printf("[%.18f, %.18f]\n", quosum_pars[k], quosum_pars[k+1]);
	}
*/
	
	/* Now do wedge terms */
	for( i=0; i<12; i++ ) {
		wedgesum_pars[i] = 0.0;
	}
	for( i=0; i<3; i++ ) {
		j = 2*i;
		h[0] = 0.5*y[j];
		h[1] = 0.5*y[j+1];
		ooh[0] = 4.0*oo2y[j];
		ooh[1] = 4.0*oo2y[j+1];
		/* Check to see if h_i < t_0 */
		if( h[0] < t0_val[1] ) {
			i_vor_beta( h, beta_val );
			i_vor_beta_x( h, ooh, beta_h_val );
			i_index_dih( i, x, temp );
			if( temp[0] < 0.0 )
				temp[0] = 0.0;
			i_mult( temp, delta32, dih_val );	/* clear dih_val */
			clear_index_dih_xi( i, y, x, dih_pars );
			for( k=0; k<12; k+=2 ) {
				i_mult( beta_val, dih_pars + k, partials + k );
			}
			i_mult( beta_h_val, dih_val, temp );
			ROUND_DOWN;
			for( k=0; k<12; k+=2 ) {
				wedgesum_pars[k] += partials[k];
			}
			wedgesum_pars[j] += temp[0];
			ROUND_UP;
			for( k=1; k<12; k+=2 ) {
				wedgesum_pars[k] += partials[k];
			}
			wedgesum_pars[j+1] += temp[1];
		}
	}
/*	
	for( i=0; i<12; i+=2 ) {
		i_div( wedgesum_pars + i, delta32, delta_part + i );
	}
	ROUND_NEAR;
	printf("new wedgesum_pars = \n");
	for( i=0; i<12; i+=2 ) {
		printf("[%.18f, %.18f]\n", delta_part[i], delta_part[i+1]);
	}
*/
	
	temp[0] = ONE_6_LO;
	temp[1] = ONE_6_HI;
	
	for( i=0; i<12; i+=2 ) {
		i_mult( temp, wedgesum_pars + i, delta_part + i );
	}
	
	ROUND_DOWN;
	for( i=0; i<12; i+=2 ) {
		partials[i] = delta_part[i] - clear_quosum_pars[i+1];
	}
	ROUND_UP;
	for( i=1; i<12; i+=2 ) {
		partials[i] = delta_part[i] - clear_quosum_pars[i-1];
	}
	
	temp[0] = DOCT_LO;
	temp[1] = DOCT_HI;
	for( i=0; i<12; i+=2 ) {
		i_mult( temp, partials + i, delta_part + i );
	}

	ROUND_DOWN;
	for( i=0; i<12; i+=2 ) {
		out_partials[i] = solmult_pars[i] + 4.0*delta_part[i];
	}
	ROUND_UP;
	for( i=1; i<12; i+=2 ) {
		out_partials[i] = solmult_pars[i] + 4.0*delta_part[i];
	}
}


void i_index_dih( int index, double xp[12], double out[2] )
{
	int i;
	double x[12];
	
	switch( index ) {
		case 1:	/* Edge 2 */
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
			i_dih_old( x, out );
			break;
		case 2:	/* Edge 3 */
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
			i_dih_old( x, out );
			break;
		default:
			i_dih_old( xp, out );
			break;
	}
/*
	ROUND_NEAR;
	printf("dih[%d] = [%.18f, %.18f]\n", index, out[0], out[1]);
*/
}


void i_index_dih_partials( int index, double yp[12],
	double xp[12], double out_partials[12] )
{
	int i, j, k, n;
	int map[3][6] = 
		{{0, 2, 4, 6, 8, 10},
		{2, 0, 4, 8, 6, 10},
		{4, 2, 0, 10, 8, 6}};
	double x[12], y[12], partials[12];
	
	switch( index ) {
		case 0:
			i_dih_partials( xp, yp, out_partials );
			break;
		default:	/* Edges 2 and 3 */
			for( j=0; j<6; j++ ) {
				k = 2*j;
				n = map[index][j];
				for( i=0; i<2; i++ ) {
					x[k + i] = xp[n + i];
					y[k + i] = yp[n + i];
				}
			}
			i_dih_partials( x, y, partials );
			/* Now switch back. */
			for( j=0; j<6; j++ ) {
				k = 2*j;
				n = map[index][j];
				for( i=0; i<2; i++ ) {
					out_partials[k + i] = partials[n + i];
				}
			}
			break;
	}
}


void clear_index_dih_partials( int index, double yp[12],
	double xp[12], double out_partials[12] )
{
	int i, j, k, n;
	int map[3][6] = 
		{{0, 2, 4, 6, 8, 10},
		{2, 0, 4, 8, 6, 10},
		{4, 2, 0, 10, 8, 6}};
	double x[12], y[12], partials[12], xpartials[12];
	
	switch( index ) {
		case 0:
			clear_dih_i( yp, xp, xpartials );
			break;
		default:	/* Edges 2 and 3 */
			for( j=0; j<6; j++ ) {
				k = 2*j;
				n = map[index][j];
				for( i=0; i<2; i++ ) {
					x[k + i] = xp[n + i];
					y[k + i] = yp[n + i];
				}
			}
			clear_dih_i( y, x, partials );
			/* Now switch back. */
			for( j=0; j<6; j++ ) {
				k = 2*j;
				n = map[index][j];
				for( i=0; i<2; i++ ) {
					xpartials[k + i] = partials[n + i];
				}
			}
			break;
	}
	/* Switch to y partials */
	for( i=0; i<12; i+=2 ) {
		x[0] = 2.0*yp[i];
		x[1] = 2.0*yp[i+1];
		i_mult( x, xpartials + i, out_partials + i );
	}
}


void clear_index_dih_xi( int index, double yp[12],
	double xp[12], double out_partials[12] )
{
	int i, j, k, n;
	int map[3][6] = 
		{{0, 2, 4, 6, 8, 10},
		{2, 0, 4, 8, 6, 10},
		{4, 2, 0, 10, 8, 6}};
	double x[12], y[12], partials[12];
	
	switch( index ) {
		case 0:
			clear_dih_i( yp, xp, out_partials );
			break;
		default:	/* Edges 2 and 3 */
			for( j=0; j<6; j++ ) {
				k = 2*j;
				n = map[index][j];
				for( i=0; i<2; i++ ) {
					x[k + i] = xp[n + i];
					y[k + i] = yp[n + i];
				}
			}
			clear_dih_i( y, x, partials );
			/* Now switch back. */
			for( j=0; j<6; j++ ) {
				k = 2*j;
				n = map[index][j];
				for( i=0; i<2; i++ ) {
					out_partials[k + i] = partials[n + i];
				}
			}
			break;
	}
}


/* Formula: 2*y1*y2*y6( x1*x1*(x1 + x2 + x6) + 
					(x2-x6)*(x2-x6)*(3*(x2+x6)-5*x1) )/tomsu^(5/2) */
void crad3len_y1y1( double y[6], double x[6], double out[2] )
{
	double uval52[2], temp[2], temp2[2];
	double num1[2], num2[2], x2_x62[2], term[2];
	double prod[2], pterm[2];
	
	i_tomsu( x, temp2 );
	i_sqrt( temp2, temp );
	ROUND_DOWN;
	uval52[0] = temp2[0]*temp2[0]*temp[0];
	num1[0] = 2.0*y[0]*y[2]*y[4];
	num2[0] = x[2] - x[5];
	term[0] = 3.0*(x[2] + x[4]) + (-5.0)*x[1];
	pterm[0] = x[0]*x[0]*(x[0] + x[2] + x[4]);
	ROUND_UP;
	uval52[1] = temp2[1]*temp2[1]*temp[1];
	num1[1] = 2.0*y[1]*y[3]*y[5];
	num2[1] = x[3] - x[4];
	term[1] = 3.0*(x[3] + x[5]) + (-5.0)*x[0];
	pterm[1] = x[1]*x[1]*(x[1] + x[3] + x[5]);
	
	i_mult( num2, num2, x2_x62 );
	i_mult( x2_x62, term, prod );
	ROUND_DOWN;
	temp[0] = pterm[0] + prod[0];
	ROUND_UP;
	temp[1] = pterm[1] + prod[1];
	
	i_mult( num1, temp, temp2 );
	i_div( temp2, uval52, out );
}


void i_vor0_y1y1( double y[12], double x[12], 
	double out[2] )
{
	double sol_11[2], solmult[2];
	double quosum_y1y1[2], wedgesum_y1y1[2];
	double temp[2], temp2[2], sum[2];
		
	i_sol_y1y1( y, x, sol_11 );
	i_mult( sol_coeff, sol_11, solmult );
	
	/* First do quoins */
	i_quoinsum_y1y1( y, x, quosum_y1y1 );
	
	/* Now do wedge terms */
	i_betadihsum_y1y1( y, x, wedgesum_y1y1 );

	temp[0] = ONE_6_LO;
	temp[1] = ONE_6_HI;
	
	i_mult( temp, wedgesum_y1y1, temp2 );
	
	ROUND_DOWN;
	sum[0] = temp2[0] - quosum_y1y1[1];
	ROUND_UP;
	sum[1] = temp2[1] - quosum_y1y1[0];
	
	temp[0] = DOCT_LO;
	temp[1] = DOCT_HI;
	i_mult( temp, sum, temp2 );

	ROUND_DOWN;
	out[0] = solmult[0] + 4.0*temp2[0];
	ROUND_UP;
	out[1] = solmult[1] + 4.0*temp2[1];
	
/*
	ROUND_NEAR;
	printf("vor0_y1y1   = [%.18f, %.18f]\n", out[0], out[1]);
	printf("sol_11      = [%.18f, %.18f]\n", sol_11[0], sol_11[1]);
	printf("quosum_11   = [%.18f, %.18f]\n", quosum_y1y1[0], quosum_y1y1[1]);
	printf("wedgesum_11 = [%.18f, %.18f]\n", wedgesum_y1y1[0], wedgesum_y1y1[1]);
*/
}


void clear_vor0_y1y1( double y[12], double x[12], 
	double out[2] )
{
	double sol_11[2], solmult[2];
	double quosum_y1y1[2], wedgesum_y1y1[2];
	double temp[2], temp2[2], sum[2];
	double clear_i[12], clear_ij[6][12];
	double delta32[2];
		
	i_bigdelta_best( x, temp );
	if( temp[0] < 0.0 )
		temp[0] = 0.0;
	i_sqrt( temp, temp2 );
	ROUND_DOWN;
	delta32[0] = temp[0]*temp2[0];
	ROUND_UP;
	delta32[1] = temp[1]*temp2[1];

	clear_sol_ij( y, x, clear_i, clear_ij );
	sol_11[0] = clear_ij[0][0];
	sol_11[1] = clear_ij[0][1];
	i_mult( sol_coeff, sol_11, solmult );
	
	/* First do quoins */
	i_quoinsum_y1y1( y, x, temp );
	i_mult( temp, delta32, quosum_y1y1 );
	
	/* Now do wedge terms */
	clear_betadihsum_y1y1( y, x, wedgesum_y1y1 );

	temp[0] = ONE_6_LO;
	temp[1] = ONE_6_HI;
	
	i_mult( temp, wedgesum_y1y1, temp2 );
	
	ROUND_DOWN;
	sum[0] = temp2[0] - quosum_y1y1[1];
	ROUND_UP;
	sum[1] = temp2[1] - quosum_y1y1[0];
	
	temp[0] = DOCT_LO;
	temp[1] = DOCT_HI;
	i_mult( temp, sum, temp2 );

	ROUND_DOWN;
	out[0] = solmult[0] + 4.0*temp2[0];
	ROUND_UP;
	out[1] = solmult[1] + 4.0*temp2[1];
	
/*
	ROUND_NEAR;
	printf("clear_vor0_y1y1 = [%.18f, %.18f]\n", out[0], out[1]);
	printf("sol_11      = [%.18f, %.18f]\n", sol_11[0], sol_11[1]);
	printf("quosum_11   = [%.18f, %.18f]\n", quosum_y1y1[0], quosum_y1y1[1]);
	printf("wedgesum_11 = [%.18f, %.18f]\n", wedgesum_y1y1[0], wedgesum_y1y1[1]);
*/
}


void clear_vor0_x1x1( double y[12], double x[12], 
	double oo2y[12], double out[2] )
{
	double sol_11[2], solmult[2];
	double quosum_y1y1[2], quosum_x1x1[2], wedgesum_x1x1[2];
	double temp[2], temp2[2], sum[2];
	double clear_i[12], clear_ij[6][12];
	double delta32[2];
	
	i_bigdelta_best( x, temp );
	if( temp[0] < 0.0 )
		temp[0] = 0.0;
	i_sqrt( temp, temp2 );
	ROUND_DOWN;
	delta32[0] = temp[0]*temp2[0];
	ROUND_UP;
	delta32[1] = temp[1]*temp2[1];

	clear_sol_ij( y, x, clear_i, clear_ij );
/*
	ROUND_NEAR;
	printf("sol_y1 = [%.18f, %.18f]\n", clear_i[0],
		clear_i[1] );
	printf("sol_y1y1 = [%.18f, %.18f]\n", clear_ij[0][0],
		clear_ij[0][1]);
*/
	sol_11[0] = clear_ij[0][0];
	sol_11[1] = clear_ij[0][1];
	
	/* Now convert to x */
	i_mult( clear_i, oo2y, temp );
	ROUND_DOWN;
	temp2[0] = sol_11[0] - 2.0*temp[1];
	sum[0] = oo2y[0]*oo2y[0];
	ROUND_UP;
	temp2[1] = sol_11[1] - 2.0*temp[0];
	sum[1] = oo2y[1]*oo2y[1];
	i_mult( sum, temp2, sol_11 );
	i_mult( sol_coeff, sol_11, solmult );
/*	
	ROUND_NEAR;
	printf("sol_x1x1 = [%.18f, %.18f]\n", sol_11[0],
		sol_11[1]);
*/	
	/* First do quoins */
	i_quoinsum_y1( y, x, clear_i );
	i_quoinsum_y1y1( y, x, quosum_y1y1 );
	/* Now convert to x */
	i_mult( clear_i, oo2y, temp );
	ROUND_DOWN;
	temp2[0] = quosum_y1y1[0] - 2.0*temp[1];
	sum[0] = oo2y[0]*oo2y[0];
	ROUND_UP;
	temp2[1] = quosum_y1y1[1] - 2.0*temp[0];
	sum[1] = oo2y[1]*oo2y[1];
	i_mult( sum, temp2, temp );
	
	i_mult( temp, delta32, quosum_x1x1 );
	
	/* Now do wedge terms */
	clear_betadihsum_x1x1( y, x, oo2y, wedgesum_x1x1 );

	temp[0] = ONE_6_LO;
	temp[1] = ONE_6_HI;
	
	i_mult( temp, wedgesum_x1x1, temp2 );
	
	ROUND_DOWN;
	sum[0] = temp2[0] - quosum_x1x1[1];
	ROUND_UP;
	sum[1] = temp2[1] - quosum_x1x1[0];
	
	temp[0] = DOCT_LO;
	temp[1] = DOCT_HI;
	i_mult( temp, sum, temp2 );

	ROUND_DOWN;
	out[0] = solmult[0] + 4.0*temp2[0];
	ROUND_UP;
	out[1] = solmult[1] + 4.0*temp2[1];
	
/*
	ROUND_NEAR;
	printf("clear_vor0_x1x1 = [%.18f, %.18f]\n", out[0], out[1]);
	printf("sol_x1x1      = [%.18f, %.18f]\n", sol_11[0], sol_11[1]);
	printf("quosum_x1x1   = [%.18f, %.18f]\n", quosum_x1x1[0], quosum_x1x1[1]);
	printf("wedgesum_x1x1 = [%.18f, %.18f]\n", wedgesum_x1x1[0], wedgesum_x1x1[1]);
*/
}


void i_sol_y1y1( double y[12], double x[12], double out[2] )
{
	double a[2], c[2], d[2], a_i[12], c_i[12], d_xi[12];
	double a_ij[6][12], c_ij[6][12], d_ij[6][12];
	double delta[2], delta_part[12], delta_ij[6][12];
	double sol_ij[6][12], d_i[12], d_xij[6][12];
	
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
	sp_sol_ij( c, c_i, c_ij, sol_ij );
	
	out[0] = sol_ij[0][0];
	out[1] = sol_ij[0][1];
}


void i_betadihsum_y1y1( double y[12], double x[12], 
	double out[2] )
{
	int i, j, k, n, index;
	int map[3][6] = 
		{{0, 2, 4, 6, 8, 10},
		{2, 0, 4, 8, 6, 10},
		{4, 2, 0, 10, 8, 6}};
	double dih_xij[6][12], dih_xi[12], temp[2], temp2[2];
	double dih_yii[6], dih_val[2];
	double h[2], wedgesum[2], beta_val[2];
	double beta_h_val[2], beta_hh_val[2];
	double yp[12], xp[12];
	
	for( index=0; index<3; index++ ) {
		for( j=0; j<6; j++ ) {
			k = 2*j;
			n = map[index][j];
			for( i=0; i<2; i++ ) {
				xp[k + i] = x[n + i];
				yp[k + i] = y[n + i];
			}
		}
		i = 2*index;
		j = i + 1;
		sp_dihsecparbds( xp, dih_xij );
		i_dih_xpars( xp, dih_xi );
/*		
		ROUND_NEAR;
		printf("dih_xii = [%.18f, %.18f]\n",
			dih_xij[index][i], dih_xij[index][i+1]);
		printf("dih_xi  = [%.18f, %.18f]\n",
			dih_xi[i], dih_xi[i+1]);
*/		
		i_mult( x, dih_xij[index] + i, temp );
		ROUND_DOWN;
		dih_yii[i] = 2.0*dih_xi[i] + 4.0*temp[0];
		ROUND_UP;
		dih_yii[j] = 2.0*dih_xi[j] + 4.0*temp[1];
	}
	
	i_dih_xpars( x, dih_xi );
	i_dih_best( x, dih_xi, dih_val );
/*	
	ROUND_NEAR;
	printf("dih_yii = \n");
	for( i=0; i<6; i+=2 ) {
		printf("[%.18f, %.18f]\n", dih_yii[i], dih_yii[i+1]);
	}
*/

	/* Now do wedge terms */
	wedgesum[0] = 0.0;
	wedgesum[1] = 0.0;
	for( i=0; i<3; i++ ) {
		j = 2*i;
		h[0] = 0.5*y[j];
		h[1] = 0.5*y[j+1];
		/* Check to see if h_i < t_0 */
		if( h[0] < t0_val[1] ) {
			i_vor_beta( h, beta_val );
			i_mult( beta_val, dih_yii + j, temp );
			ROUND_DOWN;
			wedgesum[0] += temp[0];
			ROUND_UP;
			wedgesum[1] += temp[1];
		}
	}
/*	
	ROUND_NEAR;
	printf("wedgesum = [%.18f, %.18f]\n", 
		wedgesum[0], wedgesum[1]);
*/	
	h[0] = 0.5*y[0];
	h[1] = 0.5*y[1];
	/* Check to see if h_1 < t_0 */
	if( h[0] < t0_val[1] ) {
		i_vor_beta_h( h, beta_h_val );
		i_mult( beta_h_val, dih_xi, temp );
		i_mult( y, temp, temp2 );
		i_vor_beta_hh( h, beta_hh_val );
		i_mult( beta_hh_val, dih_val, temp );
		ROUND_DOWN;
		wedgesum[0] += 2.0*temp2[0] + 0.25*temp[0];
		ROUND_UP;
		wedgesum[1] += 2.0*temp2[1] + 0.25*temp[1];
	}
/*	
	ROUND_NEAR;
	printf("beta_h  = [%.18f, %.18f]\n", 
		beta_h_val[0], beta_h_val[1]);
	printf("beta_hh = [%.18f, %.18f]\n", 
		beta_hh_val[0], beta_hh_val[1]);
*/	
	out[0] = wedgesum[0];
	out[1] = wedgesum[1];
}


void clear_betadihsum_y1y1( double y[12], double x[12], 
	double out[2] )
{
	int i, j, k, n, index;
	int map[3][6] = 
		{{0, 2, 4, 6, 8, 10},
		{2, 0, 4, 8, 6, 10},
		{4, 2, 0, 10, 8, 6}};
	double dih_xij[6][12], dih_xi[12], temp[2], temp2[2];
	double dih_yii[6], dih_val[2];
	double h[2], wedgesum[2], beta_val[2];
	double beta_h_val[2], beta_hh_val[2];
	double yp[12], xp[12], delta32[2];
	
	i_bigdelta_best( x, temp );
	if( temp[0] < 0.0 )
		temp[0] = 0.0;
	i_sqrt( temp, temp2 );
	ROUND_DOWN;
	delta32[0] = temp[0]*temp2[0];
	ROUND_UP;
	delta32[1] = temp[1]*temp2[1];

	for( index=0; index<3; index++ ) {
		for( j=0; j<6; j++ ) {
			k = 2*j;
			n = map[index][j];
			for( i=0; i<2; i++ ) {
				xp[k + i] = x[n + i];
				yp[k + i] = y[n + i];
			}
		}
		i = 2*index;
		j = i + 1;
		clear_dih_ij( yp, xp, dih_xi, dih_xij );
/*		
		ROUND_NEAR;
		printf("dih_xii = [%.18f, %.18f]\n",
			dih_xij[index][i], dih_xij[index][i+1]);
		printf("dih_xi  = [%.18f, %.18f]\n",
			dih_xi[i], dih_xi[i+1]);
*/		
		i_mult( x, dih_xij[index] + i, temp );
		ROUND_DOWN;
		dih_yii[i] = 2.0*dih_xi[i] + 4.0*temp[0];
		ROUND_UP;
		dih_yii[j] = 2.0*dih_xi[j] + 4.0*temp[1];
	}
	
	clear_dih_i( y, x, dih_xi );
	i_dih_best( x, dih_xi, temp );
	i_mult( delta32, temp, dih_val );
/*	
	ROUND_NEAR;
	printf("dih_yii = \n");
	for( i=0; i<6; i+=2 ) {
		printf("[%.18f, %.18f]\n", dih_yii[i], dih_yii[i+1]);
	}
*/

	/* Now do wedge terms */
	wedgesum[0] = 0.0;
	wedgesum[1] = 0.0;
	for( i=0; i<3; i++ ) {
		j = 2*i;
		h[0] = 0.5*y[j];
		h[1] = 0.5*y[j+1];
		/* Check to see if h_i < t_0 */
		if( h[0] < t0_val[1] ) {
			i_vor_beta( h, beta_val );
			i_mult( beta_val, dih_yii + j, temp );
			ROUND_DOWN;
			wedgesum[0] += temp[0];
			ROUND_UP;
			wedgesum[1] += temp[1];
		}
	}
/*	
	ROUND_NEAR;
	printf("wedgesum = [%.18f, %.18f]\n", 
		wedgesum[0], wedgesum[1]);
*/	
	h[0] = 0.5*y[0];
	h[1] = 0.5*y[1];
	/* Check to see if h_1 < t_0 */
	if( h[0] < t0_val[1] ) {
		i_vor_beta_h( h, beta_h_val );
		i_mult( beta_h_val, dih_xi, temp );
		i_mult( y, temp, temp2 );
		i_vor_beta_hh( h, beta_hh_val );
		i_mult( beta_hh_val, dih_val, temp );
		ROUND_DOWN;
		wedgesum[0] += 2.0*temp2[0] + 0.25*temp[0];
		ROUND_UP;
		wedgesum[1] += 2.0*temp2[1] + 0.25*temp[1];
	}
/*	
	ROUND_NEAR;
	printf("beta_h  = [%.18f, %.18f]\n", 
		beta_h_val[0], beta_h_val[1]);
	printf("beta_hh = [%.18f, %.18f]\n", 
		beta_hh_val[0], beta_hh_val[1]);
*/	
	out[0] = wedgesum[0];
	out[1] = wedgesum[1];
}


void clear_betadihsum_x1x1( double y[12], double x[12], 
	double oo2y[12], double out[2] )
{
	int i, j, k, n, index;
	int map[3][6] = 
		{{0, 2, 4, 6, 8, 10},
		{2, 0, 4, 8, 6, 10},
		{4, 2, 0, 10, 8, 6}};
	double dih_xij[6][12], dih_xi[12], temp[2], temp2[2];
	double dih_xii[6], dih_val[2];
	double h[2], wedgesum[2], beta_val[2];
	double beta_h_val[2], beta_hh_val[2];
	double yp[12], xp[12], delta32[2];
	double ooh[2];
	
	i_bigdelta_best( x, temp );
	if( temp[0] < 0.0 )
		temp[0] = 0.0;
	i_sqrt( temp, temp2 );
	ROUND_DOWN;
	delta32[0] = temp[0]*temp2[0];
	ROUND_UP;
	delta32[1] = temp[1]*temp2[1];

	for( index=0; index<3; index++ ) {
		for( j=0; j<6; j++ ) {
			k = 2*j;
			n = map[index][j];
			for( i=0; i<2; i++ ) {
				xp[k + i] = x[n + i];
				yp[k + i] = y[n + i];
			}
		}
		i = 2*index;
		j = i + 1;
		clear_dih_ij( yp, xp, dih_xi, dih_xij );
/*		
		ROUND_NEAR;
		printf("dih_xii = [%.18f, %.18f]\n",
			dih_xij[index][i], dih_xij[index][i+1]);
		printf("dih_xi  = [%.18f, %.18f]\n",
			dih_xi[i], dih_xi[i+1]);
*/	
		dih_xii[i] = dih_xij[index][i];
		dih_xii[j] = dih_xij[index][j];
	}
	
	clear_dih_i( y, x, dih_xi );
	i_dih_best( x, dih_xi, temp );
	i_mult( delta32, temp, dih_val );
/*	
	ROUND_NEAR;
	printf("dih_xii = \n");
	for( i=0; i<6; i+=2 ) {
		printf("[%.18f, %.18f]\n", dih_xii[i], dih_xii[i+1]);
	}
	printf("dih_xi = \n");
	for( i=0; i<12; i+=2 ) {
		printf("[%.18f, %.18f]\n", dih_xi[i], dih_xi[i+1]);
	}
*/

	/* Now do wedge terms */
	wedgesum[0] = 0.0;
	wedgesum[1] = 0.0;
	for( i=0; i<3; i++ ) {
		j = 2*i;
		h[0] = 0.5*y[j];
		h[1] = 0.5*y[j+1];
		/* Check to see if h_i < t_0 */
		if( h[0] < t0_val[1] ) {
			i_vor_beta( h, beta_val );
			i_mult( beta_val, dih_xii + j, temp );
			ROUND_DOWN;
			wedgesum[0] += temp[0];
			ROUND_UP;
			wedgesum[1] += temp[1];
/*			
			ROUND_NEAR;
			printf("beta = [%.18f, %.18f]\n", 
				beta_val[0], beta_val[1]);
*/			
		}
	}
/*	
	ROUND_NEAR;
	printf("wedgesum = [%.18f, %.18f]\n", 
		wedgesum[0], wedgesum[1]);
*/	
	h[0] = 0.5*y[0];
	h[1] = 0.5*y[1];
	ooh[0] = 4.0*oo2y[0];
	ooh[1] = 4.0*oo2y[1];
	/* Check to see if h_1 < t_0 */
	if( h[0] < t0_val[1] ) {
		i_vor_beta_x( h, ooh, beta_h_val );
		i_mult( beta_h_val, dih_xi, temp2 );
		i_vor_beta_xx( h, ooh, beta_h_val, beta_hh_val );
		i_mult( beta_hh_val, dih_val, temp );
		ROUND_DOWN;
		wedgesum[0] += 2.0*temp2[0] + temp[0];
		ROUND_UP;
		wedgesum[1] += 2.0*temp2[1] + temp[1];
	}
/*	
	ROUND_NEAR;
	printf("beta_h  = [%.18f, %.18f]\n", 
		beta_h_val[0], beta_h_val[1]);
	printf("beta_hh = [%.18f, %.18f]\n", 
		beta_hh_val[0], beta_hh_val[1]);
*/	
	out[0] = wedgesum[0];
	out[1] = wedgesum[1];
}


void i_quoinsum_y1( double y[12], double x[12], 
	double out[2] )
{
	int facelist[3][3] = {{0, 2, 10},{0, 4, 8},{2, 4, 6}};
	int i, j, k, n, index;
	double face[6], facex[6], abc[6], xyz[6];
	double quo_pars[4], quosum_pars[12];
	double crad_pars[6], quob_sum[2], quo[2];
	double etaparquob[6];

	/* First do quoins */
	abc[4] = t0_val[0];
	abc[5] = t0_val[1];
	ROUND_UP;
	xyz[5] = t0_val[1]*t0_val[1];
	ROUND_DOWN;
	xyz[4] = t0_val[0]*t0_val[0];
	
	for( i=0; i<12; i++ )
		quosum_pars[i] = 0.0;
	
	/* Skip face (2 3 4) */
	for( i=0; i<2; i++ ) {
		for( j=0; j<3; j++ ) {
			index = facelist[i][j];
			k = 2*j;
			facex[k] = x[index];
			face[k]  = y[index];
			facex[k+1] = x[index+1];
			face[k+1]  = y[index+1];
		}
		i_crad3x2( facex, xyz + 2 );
		crad3len_pars( face, facex, crad_pars );
/*
		ROUND_NEAR;
		printf("crad3len_pars = \n");
		for( j=0; j<6; j+=2 ) {
			printf("[%.18f, %.18f]\n", crad_pars[j], crad_pars[j+1]);
		}
*/
		i_sqrt( xyz + 2, abc + 2 );
		quob_sum[0] = 0.0;
		quob_sum[1] = 0.0;
		/* Check to see if b < c */
		if( abc[2] < t0_val[1] ) {
			for( k=0; k<2; k++ ) {
				index = facelist[i][k];
				n = index + 1;
				abc[0] = 0.5*y[index];
				abc[1] = 0.5*y[n];
/*				
				ROUND_NEAR;
				printf("abc = \n");
				for( j=0; j<6; j+=2 ) {
					printf("[%.18f, %.18f]\n", abc[j], abc[j+1]);
				}
*/				
				/* Check to see if a < b */
				if( abc[0] < abc[3] ) {
					xyz[0] = 0.25*x[index];
					xyz[1] = 0.25*x[n];
					i_quoin_pars( abc, xyz, quo, quo_pars );
/*					
					ROUND_NEAR;
					printf("quo_pars = \n");
					for( j=0; j<4; j+=2 ) {
						printf("[%.18f, %.18f]\n", quo_pars[j], quo_pars[j+1]);
					}
*/					
					ROUND_DOWN;
					quob_sum[0] += quo_pars[2];
					ROUND_UP;
					quob_sum[1] += quo_pars[3];
					quosum_pars[index] += 0.5*quo_pars[0];
					quosum_pars[n    ] += 0.5*quo_pars[1];
				}
			}
		}
		for( k=0; k<6; k+=2 ) {
			i_mult( quob_sum, crad_pars + k, etaparquob + k );
		}
		ROUND_DOWN;
		for( j=0; j<3; j++ ) {
			index = facelist[i][j];
			k = 2*j;
			quosum_pars[index] += etaparquob[k];
		}
		ROUND_UP;
		for( j=0; j<3; j++ ) {
			index = facelist[i][j] + 1;
			k = 2*j + 1;
			quosum_pars[index] += etaparquob[k];
		}
	}
	/* This could be more efficient */
	out[0] = quosum_pars[0];
	out[1] = quosum_pars[1];
}


void i_quoinsum_y1y1( double y[12], double x[12], 
	double out[2] )
{
	int facelist[3][3] = {{0, 2, 10},{0, 4, 8},{2, 4, 6}};
	int i, j, k, n, index;
	double face[6], facex[6], abc[6], xyz[6];
	double quo_pars[4], quosum_y1y1[2];
	double crad_pars[6], quob_sum[2], quobb_sum[2], quo[2];
	double etasecparquob[2];
	double temp[2], crad_y1y1[2];
	double quo_secpars[6], inside[2], quoab_sum[2];

	abc[4] = t0_val[0];
	abc[5] = t0_val[1];
	ROUND_UP;
	xyz[5] = t0_val[1]*t0_val[1];
	ROUND_DOWN;
	xyz[4] = t0_val[0]*t0_val[0];
	
	quosum_y1y1[0] = 0.0;
	quosum_y1y1[1] = 0.0;
	
	/* Skip face (2 3 4) */
	for( i=0; i<2; i++ ) {
		for( j=0; j<3; j++ ) {
			index = facelist[i][j];
			k = 2*j;
			facex[k] = x[index];
			face[k]  = y[index];
			facex[k+1] = x[index+1];
			face[k+1]  = y[index+1];
		}
#if ACUTEFACES
		s_crad3x2( facex, xyz + 2 );
		s_crad3len_pars( face, facex, crad_pars );
#else
		i_crad3x2( facex, xyz + 2 );
		crad3len_pars( face, facex, crad_pars );
#endif
		crad3len_y1y1( face, facex, crad_y1y1 );
		i_sqrt( xyz + 2, abc + 2 );
		quob_sum[0] = 0.0;
		quob_sum[1] = 0.0;
		quoab_sum[0] = 0.0;
		quoab_sum[1] = 0.0;
		quobb_sum[0] = 0.0;
		quobb_sum[1] = 0.0;
/*
		ROUND_NEAR;
		printf("face = \n");
		printf("[%.18f, %.18f]\n", face[0], face[1]);
		printf("[%.18f, %.18f]\n", face[2], face[3]);
		printf("[%.18f, %.18f]\n", face[4], face[5]);
		printf("facex = \n");
		printf("[%.18f, %.18f]\n", facex[0], facex[1]);
		printf("[%.18f, %.18f]\n", facex[2], facex[3]);
		printf("[%.18f, %.18f]\n", facex[4], facex[5]);
		printf("abc = \n");
		printf("[%.18f, %.18f]\n", abc[2], abc[3]);
		printf("[%.18f, %.18f]\n", abc[4], abc[5]);
*/
		/* Check to see if b < c */
		if( abc[2] < t0_val[1] ) {
			for( k=0; k<2; k++ ) {
				index = facelist[i][k];
				n = index + 1;
				abc[0] = 0.5*y[index];
				abc[1] = 0.5*y[n];
				/* Check to see if a < b */
				if( abc[0] < abc[3] ) {
					xyz[0] = 0.25*x[index];
					xyz[1] = 0.25*x[n];
					i_quoin_all( abc, xyz, quo, quo_pars, quo_secpars );
/*					
					ROUND_NEAR;
					printf("abc = \n");
					printf("[%.18f, %.18f]\n", abc[0], abc[1]);
					printf("[%.18f, %.18f]\n", abc[2], abc[3]);
					printf("[%.18f, %.18f]\n", abc[4], abc[5]);
					printf("xyz = \n");
					printf("[%.18f, %.18f]\n", xyz[0], xyz[1]);
					printf("[%.18f, %.18f]\n", xyz[2], xyz[3]);
					printf("[%.18f, %.18f]\n", xyz[4], xyz[5]);
					printf("quo = [%.18f, %.18f]\n", quo[0], quo[1]);
					printf("quo_pars = \n");
					printf("[%.18f, %.18f]\n", quo_pars[0], quo_pars[1]);
					printf("[%.18f, %.18f]\n", quo_pars[2], quo_pars[3]);
					printf("quo_secpars = \n");
					printf("[%.18f, %.18f]\n", quo_secpars[0], quo_secpars[1]);
					printf("[%.18f, %.18f]\n", quo_secpars[2], quo_secpars[3]);
					printf("[%.18f, %.18f]\n", quo_secpars[4], quo_secpars[5]);
					printf("crad_pars = \n");
					printf("[%.18f, %.18f]\n", crad_pars[0], crad_pars[1]);
					printf("[%.18f, %.18f]\n", crad_pars[2], crad_pars[3]);
					printf("[%.18f, %.18f]\n", crad_pars[4], crad_pars[5]);
					printf("crad_y1y1 = [%.18f, %.18f]\n", crad_y1y1[0], crad_y1y1[1]);
*/
					if( k==0 ) {
						quoab_sum[0] = quo_secpars[2];
						quoab_sum[1] = quo_secpars[3];
					}
					ROUND_DOWN;
					quob_sum[0]  += quo_pars[2];
					quobb_sum[0] += quo_secpars[4];
					if( k==0 ) {
						quosum_y1y1[0] += 0.25*quo_secpars[0];
					}
					ROUND_UP;
					quob_sum[1]  += quo_pars[3];
					quobb_sum[1] += quo_secpars[5];
					if( k==0 ) {
						quosum_y1y1[1] += 0.25*quo_secpars[1];
					}
				}
			}
		}
		i_mult( quob_sum, crad_y1y1, etasecparquob );
		i_mult( quobb_sum, crad_pars, inside );
		ROUND_DOWN;
		quosum_y1y1[0] += etasecparquob[0];
		inside[0] += quoab_sum[0];
		ROUND_UP;
		quosum_y1y1[1] += etasecparquob[1];
		inside[1] += quoab_sum[1];
		i_mult( inside, crad_pars, temp );
		ROUND_DOWN;
		quosum_y1y1[0] += temp[0];
		ROUND_UP;
		quosum_y1y1[1] += temp[1];
	}
	out[0] = quosum_y1y1[0];
	out[1] = quosum_y1y1[1];
}


/* Clear denominators to avoid singularities near
x = (4 4 4 16 8 8), which is quite degenerate. */
void superclear_vor0_y1( double y[12], double x[12], 
	double out[2] )
{
	int i, ip;
	double beta_val[2], phi0pbeta[6], p1[2], p2[2], p3[2];
	double p4[2], h[2], coeff[6];
	double delta[2], delta_i[12], delta_1j[12];
	double u126[2], u135[2], facex[6], sqrtdelta[2];
	double sign_dihpars[12];
	
	/* (1 2 6) */
	facex[0] = x[0];
	facex[1] = x[1];
	facex[2] = x[2];
	facex[3] = x[3];
	facex[4] = x[10];
	facex[5] = x[11];

	super_u( facex, u126 );
	
	/* (1 3 5) */
	facex[2] = x[4];
	facex[3] = x[5];
	facex[4] = x[8];
	facex[5] = x[9];

	super_u( facex, u135 );
	
	/* load up super_coeff + 2.0*beta */
	for( i=0; i<6; i+=2 ) {
		ip = i + 1;
		h[0] = 0.5*y[i];
		h[1] = 0.5*y[ip];
		i_vor_beta( h, beta_val );
/*
		ROUND_NEAR;
		printf("h = [%.18f, %.18f]\n", h[0], h[1]);
		printf("beta_val = [%.18f, %.18f]\n", beta_val[0],
			beta_val[1]);
*/
		ROUND_DOWN;
		phi0pbeta[i]  = super_coeff[0] + 2.0*beta_val[0];
		ROUND_UP;
		phi0pbeta[ip] = super_coeff[1] + 2.0*beta_val[1];
	}

	i_bigdelta_partials_best( x, delta, delta_i );
	if( delta[0] < 0.0 )
		delta[0] = 0.0;
	i_sqrt( delta, sqrtdelta );
	i_delta1_j( x, delta_1j );
	
	/* first coeff: delta delta_4 + x1 delta_1*delta_4 
			- 2*x1*delta*delta_14 */
	i_mult( delta, delta_i + 6, p1 );
	i_mult( delta_i, delta_i + 6, p2 );
	i_mult( delta, delta_1j + 6, p3 );
	ROUND_DOWN;
	p4[0] = p2[0] - 2.0*p3[1];
	ROUND_UP;
	p4[1] = p2[1] - 2.0*p3[0];
	i_mult( x, p4, p2 );
	ROUND_DOWN;
	coeff[0] = p1[0] + p2[0];
	ROUND_UP;
	coeff[1] = p1[1] + p2[1];
	
	/* second coeff: y1*y2*delta_3*u135 */
	ROUND_DOWN;
	p1[0] = y[0]*y[2];
	ROUND_UP;
	p1[1] = y[1]*y[3];
	i_mult( p1, delta_i + 4, p2 );
	i_mult( p2, u135, coeff + 2 );
	
	/* third  coeff: y1*y3*delta_2*u126 */
	ROUND_DOWN;
	p1[0] = y[0]*y[4];
	ROUND_UP;
	p1[1] = y[1]*y[5];
	i_mult( p1, delta_i + 2, p2 );
	i_mult( p2, u126, coeff + 4 );
	
	/* compute dih term */
	i_mult( sqrtdelta, u135, p1 );
	i_mult( p1, u126, p2 );
	i_sign_dihpars( x, delta, delta_i, sign_dihpars );
	i_dih_best( x, sign_dihpars, p3 );
	i_mult( p2, p3, p1 );
	h[0] = 0.5*y[0];
	h[1] = 0.5*y[1];
	i_vor_beta_h( h, p2 );
	i_mult( p1, p2, h );
	
	/* now compute products */
	i_mult( phi0pbeta, coeff, p1 );
	i_mult( phi0pbeta + 2, coeff + 2, p2 );
	i_mult( phi0pbeta + 4, coeff + 4, p3 );
	
	ROUND_DOWN;
	out[0] = 0.5*h[0] + p1[0] - p2[1] - p3[1];
	ROUND_UP;
	out[1] = 0.5*h[1] + p1[1] - p2[0] - p3[0];
/*	
	ROUND_NEAR;
	printf("super_coeff = [%.18f, %.18f]\n",
		super_coeff[0], super_coeff[1]);
	printf("coeff = \n");
	for( i=0; i<6; i+=2 )
		printf("[%.18f, %.18f]\n", coeff[i], coeff[i+1]);
	printf("phi0pbeta = \n");
	for( i=0; i<6; i+=2 )
		printf("[%.18f, %.18f]\n", phi0pbeta[i], phi0pbeta[i+1]);
*/
}


void super_u( double fx[6], double out[2] )
{
	double pterms[2], nterms[2];

	ROUND_DOWN;
	pterms[0] = 2.0*(fx[0]*(fx[4] + fx[2]) + fx[2]*fx[4]);
	nterms[0] = fx[0]*fx[0] + fx[2]*fx[2] + fx[4]*fx[4];
	ROUND_UP;
	pterms[1] = 2.0*(fx[1]*(fx[5] + fx[3]) + fx[3]*fx[5]);
	nterms[1] = fx[1]*fx[1] + fx[3]*fx[3] + fx[5]*fx[5];
	out[1] = pterms[1] - nterms[0];
	ROUND_DOWN;
	out[0] = pterms[0] - nterms[1];
}


void superclear_vor0_x1x1( double y[12], double x[12],
	double out[2] )
{
	int i, ip;
	double delta[2], sqrtdelta[2], delta32[2];
	double delta_i[12], delta_1j[12];
	double p1[2], p2[2], p3[2], p4[2], p5[2];
	double u126[2], u135[2], u1261[2], u1351[2];
	double u1262[2], u1352[2], facex[6];
	double part1[2], part2[2], part3[2];
	double coeff[2], seccoeff[6];
	double h[2], alpha[6], beta_val[2];
	double beta_x[2], beta_xx[2];
	double ooh[2], dih_val[2];
	double sign_dihxpars[12];
#if DEBUG
	double prod[2];
#endif
	
	/* load up alphas = super_coeff + 2.0*beta */
	for( i=0; i<6; i+=2 ) {
		ip = i + 1;
		h[0] = 0.5*y[i];
		h[1] = 0.5*y[ip];
		i_vor_beta( h, beta_val );
		ROUND_DOWN;
		alpha[i]  = 0.5*super_coeff[0] + beta_val[0];
		ROUND_UP;
		alpha[ip] = 0.5*super_coeff[1] + beta_val[1];
	}

	/* 1 2 6 */
	for( i=0; i<2; i++ ) {
		facex[i    ] = x[     i];
		facex[2 + i] = x[2  + i];
		facex[4 + i] = x[10 + i];
	}
	i_tomsu( facex, u126 );
	ROUND_DOWN;
	u1262[0] = u126[0]*u126[0];
	u1261[0] = 2.0*(facex[2] + facex[4] - facex[1]);
	ROUND_UP;
	u1262[1] = u126[1]*u126[1];
	u1261[1] = 2.0*(facex[3] + facex[5] - facex[0]);
	/* 1 3 5 */
	for( i=0; i<2; i++ ) {
		facex[2 + i] = x[4 + i];
		facex[4 + i] = x[8 + i];
	}
	i_tomsu( facex, u135 );
	ROUND_DOWN;
	u1352[0] = u135[0]*u135[0];
	u1351[0] = 2.0*(facex[2] + facex[4] - facex[1]);
	ROUND_UP;
	u1352[1] = u135[1]*u135[1];
	u1351[1] = 2.0*(facex[3] + facex[5] - facex[0]);
	
	i_bigdelta_partials_best( x, delta, delta_i );
	if( delta[0] < 0.0 ) {
		delta[0] = 0.0;
	}
	ROUND_DOWN;
	sqrtdelta[0] = sqrt( delta[0] );
	delta32[0] = sqrtdelta[0]*delta[0];
	ROUND_UP;
	sqrtdelta[1] = sqrt( delta[1] );
	delta32[1] = sqrtdelta[1]*delta[1];
	
	i_delta1_j( x, delta_1j );
	
#if DEBUG
	i_mult( u1262, u1352, p1 );
	i_mult( delta32, p1, p2 );
	p3[0] = 2.0*DOCT_LO/3.0;
	p3[1] = 2.0*DOCT_HI/3.0;
	i_div( p3, p2, prod );
#endif
	
	/* do dih_x1x1 */
	
		/* do part1 */
	i_mult( x + 6, delta_i + 6, p3 );
	i_mult( delta_i, delta_1j + 6, p5 );
	ROUND_DOWN;
	p1[0] = 4.0*delta[0] - p5[1] - 2.0*p3[1];
	ROUND_UP;
	p1[1] = 4.0*delta[1] - p5[0] - 2.0*p3[0];
	i_mult( x, p1, p2 );
	i_mult( delta_i, delta_i + 6, p3 );
	i_mult( delta, delta_1j + 6, p4 );
	ROUND_DOWN;
	part1[0] = p2[0] + 2.0*p3[0] - p4[1];
	ROUND_UP;
	part1[1] = p2[1] + 2.0*p3[1] - p4[0];
	
		/* do part2 */
	i_mult( delta, delta_1j + 6, p1 );
	i_mult( delta_i, delta_i + 6, p2 );
	ROUND_DOWN;
	p3[0] = 2.0*p1[0] - p2[1];
	ROUND_UP;
	p3[1] = 2.0*p1[1] - p2[0];
	i_mult( x, p3, p1 );
	i_mult( delta, delta_i + 6, p2 );
	ROUND_DOWN;
	part2[0] = p1[0] - p2[1];
	ROUND_UP;
	part2[1] = p1[1] - p2[0];
	
		/* do dih_x1 while we are at it */
	p3[0] = - part2[1];
	p3[1] = - part2[0];
	p1[0] = u126[0]*u135[0]*delta[0]/y[1];
	ROUND_DOWN;
	p1[1] = u126[1]*u135[1]*delta[1]/y[0];
	i_mult( p1, p3, coeff );
	
		/* do part3 */
	i_mult( x, delta_i, p3 );
	ROUND_DOWN;
	p1[0] = x[0]*delta[0]*(u1261[0]*u135[0] + u126[0]*u1351[0]);
	p2[0] = u126[0]*u135[0];
	p4[0] = delta[0] + p3[0];
	ROUND_UP;
	p1[1] = x[1]*delta[1]*(u1261[1]*u135[1] + u126[1]*u1351[1]);
	p2[1] = u126[1]*u135[1];
	p4[1] = delta[1] + p3[1];
	i_mult( p2, p4, p3 );
	ROUND_DOWN;
	part3[0] = p1[0] + 0.5*p3[0];
	ROUND_UP;
	part3[1] = p1[1] + 0.5*p3[1];

#if DEBUG
	printf("part1 = [%.18f, %.18f]\n", part1[0], part1[1]);
	printf("part2 = [%.18f, %.18f]\n", part2[0], part2[1]);
	printf("part3 = [%.18f, %.18f]\n", part3[0], part3[1]);
#endif

	ROUND_DOWN;
	p5[0] = x[0]*delta[0]*u126[0]*u135[0];
	p4[0] = y[0]*x[0];
	ROUND_UP;
	p5[1] = x[1]*delta[1]*u126[1]*u135[1];
	p4[1] = y[1]*x[1];
	i_mult( p5, part1, p1 );
	i_mult( part2, part3, p2 );
	I_ADD( p1, p2, p3 );
#if DEBUG
	printf("poly1 = [%.18f, %.18f]\n", p3[0], p3[1]);
#endif
	i_div( p3, p4, seccoeff );
	
	/* do dih2_x1x1 */
	i_mult( delta_1j + 4, u126, p1 );
	i_mult( delta_i + 4, u1261, p2 );
	I_SUB( p2, p1, p3 );
	i_mult( p3, delta, p1 );
	i_mult( delta_i, delta_i + 4, p2 );
	i_mult( p2, u126, p3 );
	ROUND_DOWN;
	p2[0] = p1[0] + 0.5*p3[0];
	p4[0] = y[2]*u1352[0];
	ROUND_UP;
	p2[1] = p1[1] + 0.5*p3[1];
	p4[1] = y[3]*u1352[1];
	i_mult( p2, p4, seccoeff + 2 );
#if DEBUG
	printf("poly2 = [%.18f, %.18f]\n", p2[0], p2[1]);
#endif
	
	/* do dih3_x1x1 */
	i_mult( delta_1j + 2, u135, p1 );
	i_mult( delta_i + 2, u1351, p2 );
	I_SUB( p2, p1, p3 );
	i_mult( p3, delta, p1 );
	i_mult( delta_i, delta_i + 2, p2 );
	i_mult( p2, u135, p3 );
	ROUND_DOWN;
	p2[0] = p1[0] + 0.5*p3[0];
	p4[0] = y[4]*u1262[0];
	ROUND_UP;
	p2[1] = p1[1] + 0.5*p3[1];
	p4[1] = y[5]*u1262[1];
	i_mult( p2, p4, seccoeff + 4 );
#if DEBUG
	printf("poly3 = [%.18f, %.18f]\n", p2[0], p2[1]);
#endif

#if DEBUG
	printf("coeff = [%.18f, %.18f]\n", coeff[0], coeff[1] );
	printf("seccoeff = \n");
	for( i=0; i<6; i+=2 )
		printf("[%.18f, %.18f]\n", seccoeff[i], seccoeff[i+1]);
#endif
		
	h[0] = 0.5*y[0];
	h[1] = 0.5*y[1];
	ROUND_DOWN;
	ooh[0] = 1.0/h[1];
	ROUND_UP;
	ooh[1] = 1.0/h[0];
	
	i_vor_beta_x( h, ooh, beta_x );
	i_vor_beta_xx( h, ooh, beta_x, beta_xx );
	i_sign_dihpars( x, delta, delta_i, sign_dihxpars );
	i_dih_best( x, sign_dihxpars, dih_val );
	if( dih_val[0] < 0.0 )
		dih_val[0] = 0.0;
	
	ROUND_DOWN;
	p2[0] = delta32[0]*u1262[0]*u1352[0]*dih_val[0];
	ROUND_UP;
	p2[1] = delta32[1]*u1262[1]*u1352[1]*dih_val[1];
	
	i_mult( p2, beta_xx, p4 );
	i_mult( beta_x, coeff, p5 );
#if DEBUG
	printf("dih_term =    [%.18f, %.18f]\n", p4[0], p4[1]);
	printf("dih_x1_term = [%.18f, %.18f]\n", p5[0], p5[1]);
#endif
	i_mult( alpha    , seccoeff    , p1 );
	i_mult( alpha + 2, seccoeff + 2, p2 );
	i_mult( alpha + 4, seccoeff + 4, p3 );
	ROUND_DOWN;
	out[0] = p1[0] + p2[0] + p3[0] + p4[0] + 2.0*p5[0];
	ROUND_UP;
	out[1] = p1[1] + p2[1] + p3[1] + p4[1] + 2.0*p5[1];
	
#if DEBUG
	ROUND_NEAR;
	printf("alphas = \n");
	for( i=0; i<6; i+=2 ) {
		printf("[%.18f, %.18f]\n", alpha[i], alpha[i+1]);
	}
	printf("beta_x =  [%.18f, %.18f]\n", beta_x[0], beta_x[1]);
	printf("beta_xx = [%.18f, %.18f]\n", beta_xx[0], beta_xx[1]);
	printf("dih_x1x1_term  =    [%.18f, %.18f]\n", 
		p1[0], p1[1]);
	printf("dih2_x1x1_term =    [%.18f, %.18f]\n", 
		p2[0], p2[1]);
	printf("dih3_x1x1_term =    [%.18f, %.18f]\n", 
		p3[0], p3[1]);

	i_mult( prod, out, p1 );
	printf("vor0_x1x1 = [%.18f, %.18f]\n", p1[0], p1[1]);
#endif
}


void i_delta1_j( double x[12], double out[12] )
{
	/*	{-2*x[4], 
	x[4] + x[5] - x[6], 
	x[4] - x[5] + x[6], 
	-2*x[1] + x[2] + x[3] - 2*x[4] + x[5] + x[6], 
	x[2] - x[3] + x[4], 
	-x[2] + x[3] + x[4]} */
		out[0] = -2.0*x[7];
		out[1] = -2.0*x[6];
		ROUND_DOWN;
		out[2] = x[6] + x[8] - x[11];
		out[4] = x[6] - x[9] + x[10];
		out[6] = -2.0*x[1] + x[2] + x[4] - 2.0*x[7] + 
								x[8] + x[10];
		out[8] = x[2] - x[5] + x[6];
		out[10] = -x[3] + x[4] + x[6];
		ROUND_UP;
		out[3] = x[7] + x[9] - x[10];
		out[5] = x[7] - x[8] + x[11];
		out[7] = -2.0*x[0] + x[3] + x[5] - 2.0*x[6] + 
								x[9] + x[11];
		out[9] = x[3] - x[4] + x[7];
		out[11] = -x[2] + x[5] + x[7];
}


/*
newcrossdiag2[{x1_,x2_,x3_,x4_,x5_,x6_,xp1_,xp5_,
      xp6_}]:=(x1 - x6 - xp1 + xp6)^2/(4*x2) - 
  (x1*(x2 - x3 + x4) + (x3 - x4)*(x6 + xp1 - xp6) - 
      x2*(2*x5 - x6 + xp1 - 2*xp5 + xp6))^2/
   (4*x2*(x2^2 + (x3 - x4)^2 - 2*x2*(x3 + x4))) + 
  (Sqrt[x1 - (x1 + x2 - x6)^2/(4*x2) + 
       (-x2^2 + x1*(x2 - x3 + x4) + (x3 - x4)*x6 + 
           x2*(x3 + x4 - 2*x5 + x6))^2/
        (4*x2*(x2^2 + (x3 - x4)^2 - 2*x2*(x3 + x4)))] + 
     Sqrt[xp1 - (x2 + xp1 - xp6)^2/(4*x2) + 
       (x2^2 + (x3 - x4)*(xp1 - xp6) - 
           x2*(x3 + x4 + xp1 - 2*xp5 + xp6))^2/
        (4*x2*(x2^2 + (x3 - x4)^2 - 2*x2*(x3 + x4)))])^2
*/
/* den = 4*x2*(x2^2 + (x3 - x4)^2 - 2*x2*(x3 + x4)) */

/* cross_diag takes as arguments
	( x1, x2, x3, x4,  x5,  x6) and
	(xp1, x2, x3, x4, xp5, xp6). */
void cross_diag2( double x[12], double xp[12], 
	double out[2] )
{
	double p1[2], p2[2], p3[2], p4[2];
	double p5[2], p6[2], p7[2], p8[2];
	double den[2];

	/* first compute den */
	ROUND_DOWN;
	p1[0] = x[4] - x[7];
	p2[0] = x[2]*(x[4] + x[6]);
	ROUND_UP;
	p1[1] = x[5] - x[6];
	p2[1] = x[3]*(x[5] + x[7]);
	i_mult( p1, p1, p3 );
	ROUND_DOWN;
	p1[0] = x[2]*x[2] + p3[0] - 2.0*p2[1];
	ROUND_UP;
	p1[1] = x[3]*x[3] + p3[1] - 2.0*p2[0];
	p2[0] = 4.0*x[2];
	p2[1] = 4.0*x[3];
	i_mult( p2, p1, den );
/*
	ROUND_NEAR;
	printf("den = [%.18f, %.18f]\n", den[0], den[1]);
*/	
	aux_cross( x, den, p1 );
	aux_cross( xp, den, p2 );
	ROUND_UP;
	p3[1] = p1[1] + p2[1];
	p8[1] = p3[1]*p3[1];
	ROUND_DOWN;
	p3[0] = p1[0] + p2[0];
	p8[0] = p3[0]*p3[0];
	
	p1[0] = x[2] - x[5] + x[6];
	p2[0] = x[4] - x[7];
	p3[0] = x[10] + xp[0] - xp[11];
	p4[0] = 2.0*x[8] - x[11] + xp[0] - 2.0*xp[9] + xp[10];
	ROUND_UP;
	p1[1] = x[3] - x[4] + x[7];
	p2[1] = x[5] - x[6];
	p3[1] = x[11] + xp[1] - xp[10];
	p4[1] = 2.0*x[9] - x[10] + xp[1] - 2.0*xp[8] + xp[11];
	i_mult( x, p1, p5 );
	i_mult( p2, p3, p6 );
	i_mult( x + 2, p4, p7 );
	ROUND_DOWN;
	p1[0] = p5[0] + p6[0] - p7[1];
	ROUND_UP;
	p1[1] = p5[1] + p6[1] - p7[0];
	i_mult( p1, p1, p2 );
	i_div( p2, den, p7 );
	
	ROUND_DOWN;
	p1[0] = x[0] - x[11] - xp[1] + xp[10];
	ROUND_UP;
	p1[1] = x[1] - x[10] - xp[0] + xp[11];
	i_mult( p1, p1, p2 );
	p3[0] = 4.0*x[2];
	p3[1] = 4.0*x[3];
	i_div( p2, p3, p6 );
	
	ROUND_DOWN;
	out[0] = p6[0] - p7[1] + p8[0];
	ROUND_UP;
	out[1] = p6[1] - p7[0] + p8[1];
/*
	ROUND_NEAR;
	printf("p6 = [%.18f, %.18f]\n", p6[0], p6[1]);
	printf("p7 = [%.18f, %.18f]\n", p7[0], p7[1]);
	printf("p8 = [%.18f, %.18f]\n", p8[0], p8[1]);
*/
}


/* Sqrt[x1 - (x1 + x2 - x6)^2/(4*x2) + 
       (-x2^2 + x1*(x2 - x3 + x4) + (x3 - x4)*x6 + 
           x2*(x3 + x4 - 2*x5 + x6))^2/den)]  */
void aux_cross( double x[12], double den[2], double out[2] )
{
	double p1[2], p2[2], p3[2], p4[2];
	double p5[2], p6[2], p7[2], p8[2];
	double val;
	
	ROUND_DOWN;
	/* x1 + x2 - x6 */
	p1[0] = x[0] + x[2] - x[11];
	/* x2 - x3 + x4 */
	p2[0] = x[2] - x[5] + x[6];
	/* x3 - x4 */
	p3[0] = x[4] - x[7];
	/* x3 + x4 - 2*x5 + x6 */
	p4[0] = x[4] + x[6] - 2.0*x[9] + x[10];
	/* x2^2 */
	p5[0] = x[2]*x[2];
	ROUND_UP;
	/* x1 + x2 - x6 */
	p1[1] = x[1] + x[3] - x[10];
	/* x2 - x3 + x4 */
	p2[1] = x[3] - x[4] + x[7];
	/* x3 - x4 */
	p3[1] = x[5] - x[6];
	/* x3 + x4 - 2*x5 + x6 */
	p4[1] = x[5] + x[7] - 2.0*x[8] + x[11];
	/* x2^2 */
	p5[1] = x[3]*x[3];
	
	i_mult( x     , p2, p6 );
	i_mult( x + 10, p3, p7 );
	i_mult( x +  2, p4, p8 );
	
	ROUND_DOWN;
	p2[0] = -p5[1] + p6[0] + p7[0] + p8[0];
	ROUND_UP;
	p2[1] = -p5[0] + p6[1] + p7[1] + p8[1];
	
	i_mult( p2, p2, p3 );
	i_div( p3, den, p2 );
	
	i_mult( p1, p1, p5 );
	p1[0] = 4.0*x[2];
	p1[1] = 4.0*x[3];
	i_div( p5, p1, p3 );
	
	ROUND_DOWN;
	val = x[0] - p3[1] + p2[0];
	if( val < 0.0 )
		val = 0.0;
	out[0] = sqrt( val );
	ROUND_UP;
	val = x[1] - p3[0] + p2[1];
	if( val < 0.0 )
		val = 0.0;
	out[1] = sqrt( val );
/* When the upper value is negative, the cross-diagonal
doesn't exist.  However, I want to give a conservative
estimate on its size, so I can (appropriately) discard
the cell.  So the cross-diagonal will be somewhat bogus,
but only when there is no cross-diagonal. */
/*
	ROUND_NEAR;
	printf("aux = [%.18f, %.18f]\n", out[0], out[1]);
*/
}


