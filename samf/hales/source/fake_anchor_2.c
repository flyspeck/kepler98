/* polynomials.c, by Samuel Ferguson, (c) 1997. */

#include "system_headers.h"
#include "sphere.h"
#include "interval.h"
#include "i_sphere.h"
#include "i_bounds.h"
#include "i_appendix.h"
#include "macros.h"


#define MAXDEPTH	32		/* Max recursion depth (was 32) */

#define SUBDIV		2			/* (4) 	number of subdivisions (1-d) */
#define SUBFRAC	0.5		/* (0.25) 	1/SUBDIV (make this exact . . .) */

#define MAXLEN		TWO51_HI

#define PARTIALS	1

#define TRASH		0
#define DEBUG		0

/* External variables */
extern int rog_dne;
extern int rog_inside;
extern double fake_anc_const[8];

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
int boundary;
double original_cell[6];

struct Cell_Info {
	int depth; 
	int cell_type;
	};

/* Prototypes */

void verifypoly( void );
void verify_cell( double tet[6], struct Cell_Info info );
void verify_cell2( double tet[6], struct Cell_Info info );

/* Code */

void main( void )
{
	ROUND_NEAR;

	printf("Welcome to sphere packing routines, relation testing division.\n");
	verifypoly();
	/* appendix_get_trash(); */
}


void verifypoly( void )
{
	int i, n;
	double tet[12];
	double diff, dt;
	time_t tstart, tstop;
	clock_t cstart, cstop;
	char *charptr;
	struct Cell_Info cell_info;
	
	printf("Welcome to Fake anchor 2\n");
	printf("Enter depth bound:  ");
	scanf("%d", &depthbd);
	printf("Enter slop:  ");
	scanf("%lf", &slop);
	printf("Depth bound = %d\n", depthbd );
	printf("slop = %g\n", slop);
	
	i_init();	/* initialize interval constants */
	appendix_init();
			
	/*
	ROUND_DOWN;
	dt = -0.003522;
	slop += dt;
	*/

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
		
		/* 1 */
		dt = 257.0;
		ROUND_UP;
		tet[0] = TWO51_LO;
		tet[1] = dt/100.0;
		/* 2 */
		tet[2] = 2.0;
		tet[3] = 2.25;
		/* 3 */
		tet[4] = 2.0;
		tet[5] = 2.25;
		/* 4 */
		dt = 277.0;
		ROUND_DOWN;
		tet[6] = dt/100.0;
		tet[7] = tet[6];
		/* 5 */
		tet[8] = 2.0;
		tet[9] = TWO51_HI;
		/* 6 */
		tet[10] = tet[2];
		tet[11] = tet[3];

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
		
		verify_cell2( tet, cell_info );
		printf("cellcount    = %16.0f\t\t", cellcount);
		printf("maxdepth = %d\n", maxdepth);
		printf("partialtries = %16.0f\n", sccutcount);
		printf("impossible   = %16.0f\t\t", gmacount);
		printf("domaincut    = %16.0f\n", vorcount);
		printf("acutecut     = %16.0f\t\t", dualcount);
		printf("verifycount  = %16.0f\n", verifycount);
		printf("partialcount = %16.0f\t\t", partialcount);
		/* diff = cellcount - verifycount - partialcount - sccutcount; */
		diff = cellcount - verifycount - partialcount;
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
	double x[6], yp[6], t[6], data[2], diff[3];
	double anc_val[2], fake_val[2], rel_part[2];
	double anc_y6[2], fake_y6[2];
	
	/* If part of the verification has already failed, bail. */
	if( failed )
		return;
		
	cellcount += 1.0;
		
	/* adjust domain appropriately */
	ROUND_UP;
	for( i=1; i<6; i+=2 ) {
		x[i] = y[i]*y[i];
	}
	ROUND_DOWN;
	for( i=0; i<6; i+=2 ) {
		x[i] = y[i]*y[i];
	}
	
	
	/* Score cell properly */
	
	/* Use y6 partial to improve bounds on anc_fun */
	anc_y6par( y, x, anc_y6 );
	max_anc_best( y, x, anc_y6, anc_val );
	
	fake_anchor2( y, fake_val );
	I_SUB( anc_val, fake_val, data );
	val = data[1];
	if( val < slop ) {
		verifycount += 1.0;
		return;
	}
	if( data[0] > slop ) {
		/* We are in big trouble. */
		failed = 1;
		printf("Cell failed:\n");
		for( i=0; i<6; i+=2 )
			printf("%0.18f\t%0.18f\n", y[i], y[i+1]);
		printf("val      = [%0.18g\t%0.18g]\n", data[0], data[1]);
		printf("anc_val  = [%0.18g\t%0.18g]\n", anc_val[0], anc_val[1]);
		printf("fake_val = [%0.18g\t%0.18g]\n", fake_val[0], fake_val[1]);
		return;
	}
	/* if verification fails, try to bump off, or at least
		reduce dimension */

	for( i=0; i<6; i++ )
		yp[i] = y[i];	/* copy y */

#if PARTIALS
	ROUND_NEAR;
	val = yp[5] - yp[4];
				/* rewrote anc_y6par so we don't need to check for discontinuity */
	if( info.depth > 10 && val > 0.0 ) {
		sccutcount += 1.0;
		fake_anchor2_y6( y, fake_y6 );
		I_SUB( anc_y6, fake_y6, rel_part );

		/* do y6 edge */
		i = 4;
		n0 = 5;
		if( rel_part[1] < 0.0 ) {
			if( yp[i] != original_cell[i] ) {
				partialcount += 1.0;
				return;	/* bumped off */
			}
			else
				yp[n0] = original_cell[i];	/* reduced dimension */
		}
		else {
			if( rel_part[0] > 0.0 ) {
				if( yp[n0] != original_cell[n0] ) {
					partialcount += 1.0;
					return;	/* bumped off */
				}
				else
					yp[i] = original_cell[n0];	/* reduced dimension */				
			}
		}
	}
#endif

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
		printf("val = [%0.18g\t%0.18g]\n", data[0], data[1]);
		printf("anc_val  = [%0.18g\t%0.18g]\n", anc_val[0], anc_val[1]);
		printf("fake_val = [%0.18g\t%0.18g]\n", fake_val[0], fake_val[1]);
		printf("anc_y6  = [%0.18g\t%0.18g]\n", anc_y6[0], anc_y6[1]);
		printf("fake_y6 = [%0.18g\t%0.18g]\n", fake_y6[0], fake_y6[1]);
		printf("y = \n");
		for( i=0; i<6; i+=2 )
			printf("%0.18f\t%0.18f\n", y[i], y[i+1]);
#if PARTIALS
		printf("rel_part = [%0.18g\t%0.18g]\n", rel_part[0], rel_part[1]);
		printf("yp = \n");
		for( i=0; i<6; i+=2 )
			printf("%0.18f\t%0.18f\n", yp[i], yp[i+1]);
#endif
	}
}


void verify_cell2( double y[6], struct Cell_Info info )
{
	int i, j, n0, n1, n2, bounds[3];
	double val;
	double x[6], yp[6], t[6], data[2], diff[3];
	double anc_val[2], fake_val[2], rel_part[6];
	double anc_2and6[4], fake_y6[2];
	
	/* If part of the verification has already failed, bail. */
	if( failed )
		return;
		
	cellcount += 1.0;
		
	/* adjust domain appropriately */
	ROUND_UP;
	for( i=1; i<6; i+=2 ) {
		x[i] = y[i]*y[i];
	}
	ROUND_DOWN;
	for( i=0; i<6; i+=2 ) {
		x[i] = y[i]*y[i];
	}
	
	/* Score cell properly */
	
	/* Use y2 and y6 partials to improve bounds on anc_fun */
	anc_2and6pars( y, x, anc_2and6 );
	max_anc_best( y, x, anc_2and6, anc_val );
	/* anc_fun( y, x, anc_val ); */
	
	fake_anchor2( y, fake_val );
	I_SUB( anc_val, fake_val, data );
	val = data[1];
	if( val < slop ) {
		verifycount += 1.0;
		return;
	}
	if( data[0] > slop ) {
		/* We are in big trouble. */
		failed = 1;
		printf("Cell failed:\n");
		for( i=0; i<6; i+=2 )
			printf("%0.18f\t%0.18f\n", y[i], y[i+1]);
		printf("val      = [%0.18g\t%0.18g]\n", data[0], data[1]);
		printf("anc_val  = [%0.18g\t%0.18g]\n", anc_val[0], anc_val[1]);
		printf("fake_val = [%0.18g\t%0.18g]\n", fake_val[0], fake_val[1]);
		return;
	}
	/* if verification fails, try to bump off, or at least
		reduce dimension */

	for( i=0; i<6; i++ )
		yp[i] = y[i];	/* copy y */

#if PARTIALS
	ROUND_NEAR;
	for( i=1; i<3; i++ ) {
		j = 2*i;
		diff[i] = yp[j+1] - yp[j];
	}
				/* rewrote anc_y6par so we don't need to check for discontinuity */
	if( info.depth > 11 ) {
	
		fake_anchor2_y6( y, fake_y6 );
		ROUND_DOWN;
		for( j=0; j<4; j+=2 ) {
			rel_part[j+2] = anc_2and6[j] - fake_y6[1];
		}
		ROUND_UP;
		for( j=1; j<4; j+=2 ) {
			rel_part[j+2] = anc_2and6[j] - fake_y6[0];
		}
		
		for( i=1; i<3; i++ ) {
			if( diff[i] > 1.0e-14 ) {
				sccutcount += 1.0;
				j = 2*i;
				n0 = j + 1;
				if( rel_part[n0] < 0.0 ) {
					if( yp[j] != original_cell[j] ) {
						partialcount += 1.0;
						return;	/* bumped off */
					}
					else
						yp[n0] = original_cell[j];	/* reduced dimension */
				}
				else {
					if( rel_part[j] > 0.0 ) {
						if( yp[n0] != original_cell[n0] ) {
							partialcount += 1.0;
							return;	/* bumped off */
						}
						else
							yp[j] = original_cell[n0];	/* reduced dimension */				
					}
				}
			}
		}
	}
#endif

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
			else {
				bounds[i] = 1;
				diff[i] = yp[j+1] - yp[j];
			}
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
					
					verify_cell2( t, info );
					ROUND_NEAR;
				}
			}
		}
	} else {	/* exceeded MAXDEPTH */
		failed = 1;
		ROUND_NEAR;
		printf("Cell failed:\n");
		printf("val = [%0.18g\t%0.18g]\n", data[0], data[1]);
		printf("anc_val  = [%0.18g\t%0.18g]\n", anc_val[0], anc_val[1]);
		printf("fake_val = [%0.18g\t%0.18g]\n", fake_val[0], fake_val[1]);
		printf("anc_y2  = [%0.18g\t%0.18g]\n", anc_2and6[0], anc_2and6[1]);
		printf("anc_y6  = [%0.18g\t%0.18g]\n", anc_2and6[2], anc_2and6[3]);
		printf("fake_y6 = [%0.18g\t%0.18g]\n", fake_y6[0], fake_y6[1]);
		printf("y = \n");
		for( i=0; i<6; i+=2 )
			printf("%0.18f\t%0.18f\n", y[i], y[i+1]);
#if PARTIALS
		printf("rel_part = [%0.18g\t%0.18g]\n", rel_part[2], rel_part[3]);
		printf("rel_part = [%0.18g\t%0.18g]\n", rel_part[4], rel_part[5]);
		printf("yp = \n");
		for( i=0; i<6; i+=2 )
			printf("%0.18f\t%0.18f\n", yp[i], yp[i+1]);
#endif
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

/*
From: Tom Hales <hales@math.lsa.umich.edu>
Date: Tue, 5 May 1998 17:52:22 -0400 (EDT)
Message-Id: <199805052152.RAA07083@diefledermaus.math.lsa.umich.edu>
Content-Type: text
Content-Length: 1543

-0.2695279326151798 + 0.1833230667013778*x - 0.02783887375001181*x^2


This bound works for case 1997... and case 2794...
but not for case 35722..

static double vorAnchor2(double y1,double y2,double y6)
    {
    double x = y2+y6;
     return (0.0002  
    -0.2695279326151798 + 0.1833230667013778*y1 - 0.02783887375001181*y1*y1
    -0.0152253
    +0.7557773828234548 - 0.3239044460786886*x + 0.0346916048615081*x*x);
    }

New Case = 19971603
numargs = 6
0.000434104367181605694
{2.50999999999999979,2,2.50999999999999979,3.20000000000000018,2.509999999999999
79,2}
 
New Case = 357223197
numargs = 6
-0.00246739112938332283
{2.50999999999999979,2,2.25,3.20000000000000018,2.25,2}
 
New Case = 279477953
numargs = 6
-6.37487851840590203e-05
{2.50999999999999979,2.25,2,3.20000000000000018,2,2.50999999999999979}
 
 smallest case = -0.00246739112938332283



This version works for 35722... (with constraints y2,y6<2.25
static double vorAnchor2(double y1,double y2,double y6)
    {
    double x = y2+y6;
     return (0.00011
    -0.2695279326151798 + 0.1833230667013778*y1 - 0.02783887375001181*y1*y1
    -0.0152253
    +0.4384851103526791 - 0.167484482694134*x + 0.01541738104479281*x*x);
    }

Numerical on fake anchor:
New Case = 108
numargs = 6
1.44577695544122198e-05
{2.50999999999999979,2,2.05988320025666116,2,2,2}
 
 smallest case = 1.44577695544122198e-05


Numerical on fake kappa:

New Case = 357223197
numargs = 6
6.34156375577399345e-05
{2.50999999999999979,2.13515069394793722,2,3.20000000000000018,2,2.1351506810856
6908}
*/

