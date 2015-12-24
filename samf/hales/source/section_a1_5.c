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

#define EPS			1.0e-14

/* External variables */
extern double t0_val[2];

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
double original_cell[12];
double one6_val[2];

struct Cell_Info {
	int depth; 
	int cell_type;
	};

/* Prototypes */

void verifypoly( void );
void verify_cell( double tet[12], struct Cell_Info info );

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
	double tet[12];
	double diff, dt;
	time_t tstart, tstop;
	clock_t cstart, cstop;
	char *charptr;
	struct Cell_Info cell_info;
	
	printf("Welcome to Section A_1.5.\n");
	printf("Enter depth bound:  ");
	scanf("%d", &depthbd);
	printf("Enter slop:  ");
	scanf("%lf", &slop);
	printf("Depth bound = %d\n", depthbd );
	printf("slop = %g\n", slop);
	
	i_init();	/* initialize interval constants */
	appendix_init();
	
	dt = 16.0;
	diff = 10.0;
	ROUND_DOWN;
	one6_val[0] = dt/diff;
	dt = 16.0;
	diff = 10.0;
	ROUND_UP;
	one6_val[1] = dt/diff;
	
	ROUND_NEAR;
	if( one6_val[1] <= one6_val[0] )
		printf("1.6   = [%.20f, %.20f]\n", one6_val[0], one6_val[1]);
	
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
		dt = 22.0;
		ROUND_UP;
		tet[0] = 2.0;
		tet[1] = dt/10.0;
		/* 2 */
		tet[2] = 2.0;
		tet[3] = TWO51_HI;
		/* 3 */
		tet[4] = TWO51_LO;
		tet[5] = TWO51_HI;
		/* 4 */
		dt = 32.0;
		ROUND_DOWN;
		tet[6] = dt/10.0;
		ROUND_UP;
		dt = 32.0;
		tet[7] = dt/10.0;
		/* 5 */
		ROUND_DOWN;
		dt = 32.0;
		tet[8] = dt/10.0;
		ROUND_UP;
		dt = 32.0;
		tet[9] = dt/10.0;
		/* 6 */
		tet[10] = 2.0;
		tet[11] = 2.0;

		switch( n ) {
			case 1:
				break;
			default:
				printf("fell off end: this can't happen\n");
				break;
		}
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", tet[i], tet[i+1]);
		for( i=0; i<12; i++ )
			original_cell[i] = tet[i];
		
		verify_cell( tet, cell_info );
		printf("cellcount    = %16.0f\t\t", cellcount);
		printf("maxdepth = %d\n", maxdepth);
		printf("deltacut     = %16.0f\n", sccutcount);
		printf("impossible   = %16.0f\t", gmacount);
		printf("domaincut    = %16.0f\n", vorcount);
		printf("acutecut     = %16.0f\t\t", dualcount);
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


void verify_cell( double y[12], struct Cell_Info info )
{
	int i, j, n0, n1, n2, n3, n4, n5, bounds[6];
	double val;
	double x[12], yp[12], t[12], data[2], diff[6];
	double cos2dih3_val[2], cos2betapsi_val[2];
	double cos2psi_val[2];
	double face[6], facex[6], sign[2], eta2_val[2];
	
	/* If part of the verification has already failed, bail. */
	if( failed )
		return;
		
	cellcount += 1.0;
		
	/* adjust domain appropriately */
	ROUND_UP;
	for( i=1; i<12; i+=2 ) {
		x[i] = y[i]*y[i];
	}
	ROUND_DOWN;
	for( i=0; i<12; i+=2 ) {
		x[i] = y[i]*y[i];
	}
	
	for( i=0; i<2; i++ ) {
		face[i] = y[i];
		facex[i] = x[i];
		face[i + 2] = t0_val[i];
		face[i + 4] = one6_val[i];
	}
	for( i=2; i<6; i+=2 ) {
		facex[i] = face[i]*face[i];
	}
	ROUND_UP;
	for( i=3; i<6; i+=2 ) {
		facex[i] = face[i]*face[i];
	}
	cos_arc_fun( face, facex, data );
	ROUND_DOWN;
	cos2psi_val[0] = data[0]*data[0];
	ROUND_UP;
	cos2psi_val[1] = data[1]*data[1];

	/* 1 2 6 */
	for( i=0; i<2; i++ ) {
		facex[i] = x[i];
		facex[i + 2] = x[i + 2];
		facex[i + 4] = x[i + 10];
	}
	i_crad3x2( facex, eta2_val );
	
	for( i=0; i<2; i++ ) {
		facex[i] = 0.25*x[i + 2];
		facex[i + 2] = eta2_val[i];
	}
	i_div( x, cos2psi_val, data );
	facex[4] = 0.25*data[0];
	facex[5] = 0.25*data[1];
	
	/* Score cell properly */
	cos2_dih2( x, sign, cos2dih3_val );
	cos2_rog_dih( facex, cos2betapsi_val );
	I_SUB( cos2dih3_val, cos2betapsi_val, data );
	val = data[1];
	if( val < slop ) {
		verifycount += 1.0;
		return;
	}
	if( data[0] > slop ) {
		/* We are in big trouble. */
		failed = 1;
		ROUND_NEAR;
		printf("Cell failed:\n");
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", y[i], y[i+1]);
		printf("val = [%0.18g\t%0.18g]\n", data[0], data[1]);
		printf("cos2_dih2 = [%0.18g\t%0.18g]\n", cos2dih3_val[0], cos2dih3_val[1]);
		printf("cos2_betapsi_val = [%0.18g\t%0.18g]\n", cos2betapsi_val[0], cos2betapsi_val[1]);
		printf("cos2psi_val = [%0.18g\t%0.18g]\n", cos2psi_val[0], cos2psi_val[1]);
		return;
	}
	/* if verification fails, try to bump off, or at least
		reduce dimension */

	for( i=0; i<12; i++ )
		yp[i] = y[i];	/* copy y */
	
	if( info.depth < depthbd ) {
		/* If all else fails, subdivide.  */
		info.depth++;
		if( info.depth > maxdepth )
			maxdepth = info.depth;
		ROUND_NEAR;
		j = 0;
		for( i=0; i<6; i++ ) {	/* scaled differentials */
			diff[i] = 0.5*(yp[j+1] - yp[j]);
			j += 2;
			if( diff[i] > EPS )
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
					for( n3=0; n3<bounds[3]; n3++ ) {
						val = yp[6];
						t[6] = val + n3*diff[3];
						t[7] = val + (n3+1)*diff[3];
						for( n4=0; n4<bounds[4]; n4++ ) {
							val = yp[8];
							t[8] = val + n4*diff[4];
							t[9] = val + (n4+1)*diff[4];
							for( n5=0; n5<bounds[5]; n5++ ) {
								val = yp[10];
								t[10] = val + n5*diff[5];
								t[11] = val + (n5+1)*diff[5];
					
								verify_cell( t, info );
								ROUND_NEAR;
							}
						}
					}
				}
			}
		}
	} else {	/* exceeded MAXDEPTH */
		failed = 1;
		ROUND_NEAR;
		printf("Cell failed:\n");
		printf("cos2_dih2 = [%0.18g\t%0.18g]\n", cos2dih3_val[0], cos2dih3_val[1]);
		printf("cos2_betapsi_val = [%0.18g\t%0.18g]\n", cos2betapsi_val[0], cos2betapsi_val[1]);
		printf("cos2psi_val = [%0.18g\t%0.18g]\n", cos2psi_val[0], cos2psi_val[1]);
		printf("facex = \n");
		for( i=0; i<6; i+=2 ) {
			printf("[%0.18g\t%0.18g]\n", facex[i],facex[i+1]);
		}
		printf("val = [%0.18g\t%0.18g]\n", data[0], data[1]);
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", y[i], y[i+1]);
	}
}


/*
2.123 2.29102 2.341248 2.011234 2.023309 2.03982134
2.0123 2.1765 2.0987 2.1567 2.0324 2.2187
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
