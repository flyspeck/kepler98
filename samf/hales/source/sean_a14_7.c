/* polynomials.c, by Samuel Ferguson, (c) 1997. */

#include "system_headers.h"
#include "sphere.h"
#include "interval.h"
#include "i_sphere.h"
#include "i_bounds.h"
#include "macros.h"


#define MAXDEPTH	32		/* Max recursion depth (was 32) */

#define SUBDIV		2			/* (4) 	number of subdivisions (1-d) */
#define SUBFRAC	0.5		/* (0.25) 	1/SUBDIV (make this exact . . .) */

#define MAXLEN		TWO51_HI

#define PARTIALS	0
#define TRASH		0

/* External variables */
/* External prototypes */

void i_init( void );

/* Global variables */

double cellcount;
double verifycount;
double partialcount;
double sccutcount;
double gmacount;
double vorcount;
double dualcount;
double celltotal;
int maxdepth;
int depthbd;
int failed;
double slop;
double t0_val[2];
double t0_val2[2];
double oot0_val[2];
double phi0_val[2];
double zetapt_val[2];
int boundary;
int constraints;
int y6constraint;
int dov1fun;
double original_cell[6];

struct Cell_Info {
	int depth; 
	int cell_type;
	};

double pcoeff[5][5];
double ncoeff[5][5];

/* Prototypes */

void verifypoly( void );
void verify_cell( double tet[6], struct Cell_Info info );
void i_poly( double xint[4], double out[2] );
void set_coeffs( void );
void get_trash( void );
double test_poly( double xin[2] );
void phi_fun( double h[2], double t[2], double out[2] );
void a_fun( double h[2], double out[2] );
void b_fun( double y[2], double out[2] );
void b0_fun( double y[2], double out[2] );
void v0_fun( double y[12], double x[12], double out[2] );
void v1_fun( double y[12], double x[12], double out[2] );
void quopar_fun( double y[6], double x[6], double out[2] );

/* Code */

void main( void )
{
	ROUND_NEAR;

	printf("Welcome to sphere packing routines, relation testing division.\n");
	verifypoly();
	/* get_trash(); */
}


void verifypoly( void )
{
	int i, n;
	double tet[12], temp[2], temp2[2];
	double diff, dt;
	time_t tstart, tstop;
	clock_t cstart, cstop;
	char *charptr;
	struct Cell_Info cell_info;
	double two189[2], three2[2];
	double sean_val[2];

/*
	dt = 		251682.0;
	diff =	100000.0;
	ROUND_DOWN;
	sean_val[0] = dt/diff;
	ROUND_UP;
	dt = 		251682.0;
	diff =	100000.0;
	sean_val[1] = dt/diff;
*/	
	dt = 2.51682;
	i_recognize( dt, sean_val );
	
	ROUND_NEAR;
	if( sean_val[1] <= sean_val[0] )
		printf("sean_val = [%.20f, %.20f]\n", sean_val[0], sean_val[1]);
	
	printf("Welcome to Sean A_14.7.\n");
	printf("Enter depth bound:  ");
	scanf("%d", &depthbd);
	printf("Enter slop:  ");
	scanf("%lf", &slop);
	printf("Depth bound = %d\n", depthbd );
	printf("slop = %g\n", slop);
	
	i_init();	/* initialize interval constants */
	
					/* initialize local constants */
	t0_val[0] = 0.5*sean_val[0];
	t0_val[1] = 0.5*sean_val[1];
	ROUND_DOWN;
	t0_val2[0] = t0_val[0]*t0_val[0];
	ROUND_UP;
	t0_val2[1] = t0_val[1]*t0_val[1];
	ROUND_DOWN;
	oot0_val[0] = 1.0/t0_val[1];
	ROUND_UP;
	oot0_val[1] = 1.0/t0_val[0];
	phi_fun( t0_val, t0_val, phi0_val );
	/* compute zetapt_val here */
	dt = 2.0;
	ROUND_DOWN;
	temp[0] = sqrt(dt)/5.0;
	ROUND_UP;
	temp[1] = sqrt(dt)/5.0;
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
	
	/* put constant on right-hand side */
	ROUND_DOWN;
	dt = 834.0;
	diff = dt/1000.0;
	slop = slop + diff/DOCT_HI;
	
	ROUND_UP;
	dt = 2189.0;
	two189[1] = dt/1000.0;
	dt = 32.0;
	three2[1] = dt/10.0;
	ROUND_DOWN;
	dt = 32.0;
	three2[0] = dt/10.0;
	dt = 2189.0;
	two189[0] = dt/1000.0;
	
#if TRASH
	get_trash();
	exit(1);
#endif
	
	time(&tstart);
	cstart = clock();
	
	charptr = ctime( &tstart );
	printf("Starting time:  \t");
	puts( charptr );
	
	celltotal = 0.0;
	for( n=1; n<2; n++ ) {
		printf("Starting case %d:\n", n);
		cellcount = 0.0;
		verifycount = 0.0;
		partialcount = 0.0;
		sccutcount = 0.0;
		vorcount = 0.0;
		gmacount = 0.0;
		dualcount = 0.0;
		maxdepth = 0;
		failed = 0;
		cell_info.depth = 0;
		cell_info.cell_type = -1;
		boundary = 0;
		constraints = 0;
		y6constraint = 0;
		dov1fun = 0;
		
		tet[0] = 2.0;
		tet[1] = sean_val[1];
		
		tet[2] = 2.0;
		tet[3] = sean_val[1];
		
		tet[4] = 2.0;
		tet[5] = sean_val[1];
		
		
		switch( n ) {
			case 1:
				break;
			default:
				printf("fell off end: this can't happen\n");
				break;
		}
		for( i=0; i<6; i+=2 )
			printf("%0.18f\t%0.18f\n", tet[i], tet[i+1]);
		for( i=0; i<6; i++ )
			original_cell[i] = tet[i];
		
		set_coeffs();
		verify_cell( tet, cell_info );
		printf("cellcount    = %16.0f\t\t", cellcount);
		printf("maxdepth = %d\n", maxdepth);
		printf("deltacut     = %16.0f\n", sccutcount);
		printf("impossible   = %16.0f\t", gmacount);
		printf("domaincut    = %16.0f\n", vorcount);
		printf("facecut      = %16.0f\t\t", dualcount);
		printf("verifycount  = %16.0f\n", verifycount);
		printf("partialcount = %16.0f\t\t", partialcount);
		diff = cellcount - verifycount - partialcount - sccutcount;
		printf("diff         = %16.0f\n", diff);
		/*	printf("frac = %g\n", ((double) i)/cellcount);	*/
		if( !failed )
			printf("Verification succeeded.\n\n");
		else
			printf("MAXDEPTH exceeded.\n\n");
		celltotal += cellcount;
	} /* end case loop */
		
	ROUND_NEAR;
		
	time(&tstop);
	cstop = clock();
	printf("Done.\n");
	diff = ( (double) (cstop - cstart) )/CLOCKS_PER_SEC;
	charptr = ctime( &tstop );
	printf("Ending time:    \t");
	puts( charptr );
#if DA_SYSTEM == 2
	dt = ( double ) tstop - tstart;
#else
	dt = difftime( tstop, tstart );
#endif
	printf("Elapsed time: %g seconds\n", dt);
	printf("diff = %g\t\t", diff);
	printf("%g cells per second\n\n\n", celltotal/dt);
}


void verify_cell( double y[6], struct Cell_Info info )
{
	int i, j, n0, n1, n2, bounds[3];
	double val;
	double x[6], delta[2];
	double yp[6], t[6], diff[3];
	
	/* If part of the verification has already failed, bail. */
	if( failed )
		return;
		
	cellcount += 1.0;
		
	/* adjust domain appropriately */
	ROUND_UP;
	/*
	if( constraints ) {
		if( y[6] > y[3] + y[5] ) {
			vorcount += 1.0;
			return;
		}
	}
	*/
	
	for( i=1; i<6; i+=2 ) {
		x[i] = y[i]*y[i];
	}
	ROUND_DOWN;
	for( i=0; i<6; i+=2 ) {
		x[i] = y[i]*y[i];
	}
		
	/* Score cell properly */
	
	quopar_fun( y, x, delta );

	val = delta[1];
	if( val < slop ) {
		verifycount += 1.0;
		return;
	}
	
	if( delta[0] > slop ) {
		/* We are in big trouble. */
		failed = 1;
		ROUND_NEAR;
		printf("Cell failed:\n");
		for( i=0; i<6; i+=2 )
			printf("%0.18f\t%0.18f\n", y[i], y[i+1]);
		printf("val = [%0.18g\t%0.18g]\n", delta[0], delta[1]);
		return;
	}
	/* if verification fails, try to bump off, or at least
		reduce dimension */

	for( i=0; i<6; i++ )
		yp[i] = y[i];	/* copy y */

	if( info.depth < depthbd ) {
		/* If all else fails, subdivide.  */
		info.depth++;
		if( info.depth > maxdepth )
			maxdepth = info.depth;
		ROUND_NEAR;
		j = 0;
		for( i=0; i<3; i++ ) {	/* scaled differentials */
			diff[i] = 0.5*(yp[j+1] - yp[j]);
			j += 2;
			if( diff[i] > 0.0 )
				bounds[i] = 2;
			else
				bounds[i] = 1;
		}
		for( n0=0; n0<bounds[0]; n0++ ) {
			val = yp[0];
			t[0] = val + n0*diff[0];
			t[1] = val + (n0+1)*diff[0];
			for( n1=0; n1<bounds[1]; n1++ ) {
				val = yp[2];
				t[2] = val + n1*diff[1];
				t[3] = val + (n1+1)*diff[1];
				for( n2=0; n2<bounds[2]; n2++ ) {
					val = yp[4];
					t[4] = val + n2*diff[2];
					t[5] = val + (n2+1)*diff[2];
					
					verify_cell( t, info );
					ROUND_NEAR;
				}
			}
		}
	} else {	/* exceeded MAXDEPTH */
		failed = 1;
		ROUND_NEAR;
		printf("Cell failed:\n");
		printf("val = [%0.18g\t%0.18g]\n", delta[0], delta[1]);
		for( i=0; i<6; i+=2 )
			printf("%0.18f\t%0.18f\n", y[i], y[i+1]);
	}
}


/*
bug =
1048576 - 393216*a + 53248*a^2 - 3072*a^3 + 64*a^4 - 
  393216*b + 167936*a*b - 22016*a^2*b + 832*a^3*b + 
  8*a^4*b + 53248*b^2 - 22016*a*b^2 + 3008*a^2*b^2 - 
  120*a^3*b^2 - a^4*b^2 - 3072*b^3 + 832*a*b^3 - 
  120*a^2*b^3 + 6*a^3*b^3 + 64*b^4 + 8*a*b^4 - a^2*b^4

dalist = CoefficientList[bug,{a,b}]
{{1048576, -393216, 53248, -3072, 64}, 
  {-393216, 167936, -22016, 832, 8}, 
  {53248, -22016, 3008, -120, -1}, 
  {-3072, 832, -120, 6, 0}, 
  {64, 8, -1, 0, 0}}
*/

void i_poly( double xint[4], double out[2] )
{
	double x, y, xsum, ysum, nterms[2], pterms[2];
	int i, j;
	
	x = xint[0];
	y = xint[2];
	xsum = 0;
	ROUND_DOWN;
	for( i=4; i>-1; i-- ) {
		ysum = 0;
		for( j=4; j>-1; j-- ) {
			ysum = pcoeff[i][j] + y*ysum;
		}
		xsum = ysum + x*xsum;
	}
	pterms[0] = xsum;
	
	xsum = 0;
	for( i=4; i>-1; i-- ) {
		ysum = 0;
		for( j=4; j>-1; j-- ) {
			ysum = ncoeff[i][j] + y*ysum;
		}
		xsum = ysum + x*xsum;
	}
	nterms[0] = xsum;
	
	x = xint[1];
	y = xint[3];
	xsum = 0;
	ROUND_UP;
	for( i=4; i>-1; i-- ) {
		ysum = 0;
		for( j=4; j>-1; j-- ) {
			ysum = pcoeff[i][j] + y*ysum;
		}
		xsum = ysum + x*xsum;
	}
	pterms[1] = xsum;
	
	xsum = 0;
	for( i=4; i>-1; i-- ) {
		ysum = 0;
		for( j=4; j>-1; j-- ) {
			ysum = ncoeff[i][j] + y*ysum;
		}
		xsum = ysum + x*xsum;
	}
	nterms[1] = xsum;
	
	out[1] = pterms[1] - nterms[0];
	ROUND_DOWN;
	out[0] = pterms[0] - nterms[1];
}


void set_coeffs( void )
{
	double coeffs[5][5] = 
		{{1048576, -393216, 53248, -3072, 64}, 
		{-393216, 167936, -22016, 832, 8}, 
		{53248, -22016, 3008, -120, -1}, 
		{-3072, 832, -120, 6, 0}, 
		{64, 8, -1, 0, 0}};
	int i, j;
	
	for( i=0; i<5; i++ ) {
		for( j=0; j<5; j++ ) {
			if( coeffs[i][j] >= 0 ) {
				pcoeff[i][j] = coeffs[i][j];
				ncoeff[i][j] = 0;
			} else {
				ncoeff[i][j] = -coeffs[i][j];
				pcoeff[i][j] = 0;
			}
		}
	}
}


void get_trash( void )
{
	double sy[6], sx[6], y[12], x[12], eps, val[2];
	int i, j;
	int interval;
	

	set_coeffs();
	ROUND_NEAR;
	j = 0;
	/*
	printf("pcoeff = \n");
	for( i=0; i<5; i++ ) {
		for( j=0; j<5; j++ ) {
			printf("\t%10.0f", pcoeff[i][j]);
		}
		printf("\n");
	}
	printf("ncoeff = \n");
	for( i=0; i<5; i++ ) {
		for( j=0; j<5; j++ ) {
			printf("\t%10.0f", ncoeff[i][j]);
		}
		printf("\n");
	}
	*/
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
		
#if OLD
		v0_fun( y, x, val );
		ROUND_NEAR;
		printf("val = [%.18g, %.18g]\n", val[0], val[1]);
		printf("\n\n");

		v1_fun( y, x, val );
		ROUND_NEAR;
		printf("val = [%.18g, %.18g]\n", val[0], val[1]);
		printf("\n\n");
#endif
		quopar_fun( y, x, val );
		ROUND_NEAR;
		printf("val = [%.18g, %.18g]\n", val[0], val[1]);
		printf("\n\n");

#if OLD
		s_crad3x2( x, val );
		ROUND_NEAR;
		printf("val = [%.18g, %.18g]\n", val[0], val[1]);
		printf("\n\n");

		s_crad3x( x, val );
		ROUND_NEAR;
		printf("val = [%.18g, %.18g]\n", val[0], val[1]);
		printf("\n\n");

		i_crad3len( y, val );
		ROUND_NEAR;
		printf("val = [%.18g, %.18g]\n", val[0], val[1]);
		printf("\n\n");
#endif
	}
}


double test_poly( double xin[2] )
{
	double x, y, xsum, ysum;
	double coeffs[5][5] = 
		{{1048576, -393216, 53248, -3072, 64}, 
		{-393216, 167936, -22016, 832, 8}, 
		{53248, -22016, 3008, -120, -1}, 
		{-3072, 832, -120, 6, 0}, 
		{64, 8, -1, 0, 0}};
	int i, j;
	
	x = xin[0];
	y = xin[1];
	xsum = 0;
	for( i=4; i>-1; i-- ) {
		ysum = 0;
		for( j=4; j>-1; j-- ) {
			ysum = coeffs[i][j] + y*ysum;
			/* printf("ysum = %g\n", ysum); */
		}
		xsum = ysum + x*xsum;
		/* printf("\txsum = %g\n", xsum); */
	}
	return( xsum );
}

/*
2.123 1.987
*/


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


void quopar_fun( double y[6], double x[6], double out[2] )
{
	double b2_val[2], t1, t2, temp[2], temp2[2], temp3[2];
	double den[2], den2[2];
	
	s_crad3x2( x, b2_val );
		
	ROUND_DOWN;
	t1 = -x[1] + x[2] + x[4];
	t2 =  x[0] - x[3] + x[4];
	temp[0] = t1*t2;
	den[0]  = sqrt( 4.0*b2_val[0] - x[1] );
	den2[0] = sqrt( 4.0*b2_val[0] - x[3] );
	ROUND_UP;
	t1 = -x[0] + x[3] + x[5];
	t2 =  x[1] - x[2] + x[5];
	temp[1] = t1*t2;
	den[1]  = sqrt( 4.0*b2_val[1] - x[0] );
	den2[1] = sqrt( 4.0*b2_val[1] - x[2] );
	
	temp2[1] = y[1]/den[0] + y[3]/den2[0];
	ROUND_DOWN;
	temp2[0] = y[0]/den[1] + y[2]/den2[1];
	
	t1 = t0_val2[0] - b2_val[1];
	if( t1 < 0.0 ) {
		t2 = 0.0;
	} else {
		t2 = sqrt( t1 )*t1;
	}

	temp3[0] = TWO_3_LO*temp[0]*temp2[0]*t2;
	
	ROUND_UP;
	t1 = t0_val2[1] - b2_val[0];
	if( t1 < 0.0 ) {
		t2 = 0.0;
	} else {
		t2 = sqrt( t1 )*t1;
	}
	
	temp3[1] = TWO_3_HI*temp[1]*temp2[1]*t2;
	out[1] = temp3[1]/x[4];
	
	ROUND_DOWN;
	out[0] = temp3[0]/x[5];
}


/*
2.123 2.29102 2.341248 2.011234 2.023309 2.03982134
2.021 2.034 2.098 2.021 2.034 2.098
*/
