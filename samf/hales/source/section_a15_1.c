/* octahedra.c, by Samuel Ferguson, (c) 1998. */

#include "system_headers.h"
#include "interval.h"
#include "i_sphere.h"
#include "i_bounds.h"
#include "i_taylor.h"
#include "i_appendix.h"
#include "i_voronoi.h"
#include "second_partials.h"
#include "macros.h"


#define MAXDEPTH	32			/* Max recursion depth (was 10) */

#define SUBDIV		2			/* (4) 	number of subdivisions (1-d) */
#define SUBFRAC	0.5		/* (0.25) 	1/SUBDIV (make this exact . . .) */

#define COUNTLIMIT	1.0e5		/* 4.0e7 */

#define SEAN					0
#define SMARTSUBDIV		1
#define XPARS					1
#define SUPERCLEAR		1


/* External variables */
extern double t0_val[2];
extern double oot0_val[2];
extern double phi0_val[2];
extern double zetapt_val[2];
extern double sol_coeff[2];
extern double super_coeff[2];

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
int failed;
double slop;
double original_cell[12];
double countlimit;

struct Cell_Info {
	int depth; 
	int diff_skip[6];
	};

/* Prototypes */

void verifyocta( void );
void verify_cell( double tet[12], struct Cell_Info info );
void verify_cell2( double tet[12], struct Cell_Info info );
void verify_cell3( double tet[12], struct Cell_Info info );
void verify_cell4( double tet[12], struct Cell_Info info );
void testquoin( void );
void testvor0( void );
void sean_init( double sean_val[2] );

/* Code */

void main( void )
{
	ROUND_NEAR;

	printf("Welcome to sphere packing routines, relation testing division.\n");
	verifyocta();
	/*
	verifyocta();
	testquoin();
	testvor0();
	*/
}


void verifyocta( void )
{
	int i, m, n, casestart, casebound;
	double tet[12], tetdiff[6], max;
	double diff, dt;
	time_t tstart, tstop;
	clock_t cstart, cstop;
	char *charptr;
	struct Cell_Info cell_info;
	double sean_val[2];
	int subcase, onefailed;

#if SEAN
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
	casestart = 2;
#else
	sean_val[0] = TWO51_LO;
	sean_val[1] = TWO51_HI;
	casestart = 1;
#endif
	casebound = 2;
	
	onefailed = 0;
	
	ROUND_NEAR;
	if( sean_val[1] <= sean_val[0] )
		printf("sean_val = [%.20f, %.20f]\n", sean_val[0], sean_val[1]);

	printf("Welcome to Section A15.  \n");
	printf("Enter cell count limit:  ");
	scanf("%lf", &countlimit);

	printf("Enter slop:  ");
	scanf("%lf", &slop);
	printf("\t\tSlop = %g\n", slop );
	
	i_init();	/* initialize interval constants */
	sean_init( sean_val );
	
	ROUND_NEAR;

	time(&tstart);
	cstart = clock();
	
	charptr = ctime( &tstart );
	printf("Starting time:  \t");
	puts( charptr );
	
	celltotal = 0.0;
	
/*	for( n=1; n < 7; n++ ) { */
	for( n=1; n < 7; n++ ) { 
/* 
		n = 6; 
*/
		
		printf("Starting case %d:\n", n);
		
		for( subcase=casestart; subcase<=casebound; subcase++ ) {
			
			if( subcase==1 ) {
				printf("\tChecking  vor0.\n");
				sean_init( sean_val );
			} else {
				printf("\tChecking -tau0.\n");
				sean_init( sean_val );
				ROUND_DOWN;
				sol_coeff[0] = phi0_val[0] - zetapt_val[1];
				ROUND_UP;
				sol_coeff[1] = phi0_val[1] - zetapt_val[0];
				/* super_coeff = 3 phi0/doct */
				tet[0] = DOCT_LO;
				tet[1] = DOCT_HI;
				i_div( sol_coeff, tet, tet + 2 );
				tet[0] = 3.0;
				tet[1] = 3.0;
				i_mult( tet, tet + 2, super_coeff );
			}
			ROUND_NEAR;
			
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
			
			/* y1 */
			tet[0] = 2.0;
			tet[1] = sean_val[1];
			/* y2 */
			tet[2] = 2.0;
			tet[3] = sean_val[1];
			/* y3 */
			tet[4] = 2.0;
			tet[5] = sean_val[1];
			/* y5 */
			tet[8] = 2.0;
			tet[9] = 2.0;
			/* y6 */
			tet[10] = sean_val[0];
			tet[11] = sean_val[1];
			
			switch( n ) {
				case 1:
					/* y6 */
					tet[10] = 2.0;
					tet[11] = 2.0;
					break;
				case 2:
					break;
				case 3:
					/* y5 */
					tet[8] = sean_val[0];
					tet[9] = sean_val[1];
					break;
				case 4:
					/* y6 */
					tet[10] = TWOSQRT2_LO;
					tet[11] = TWOSQRT2_HI;
					break;
				case 5:
					/* y5 */
					tet[8] = sean_val[0];
					tet[9] = sean_val[1];
					/* y6 */
					tet[10] = TWOSQRT2_LO;
					tet[11] = TWOSQRT2_HI;
					break;
				case 6:
					/* y5 */
					tet[8] = TWOSQRT2_LO;
					tet[9] = TWOSQRT2_HI;
					/* y6 */
					tet[10] = TWOSQRT2_LO;
					tet[11] = TWOSQRT2_HI;
					break;
				case 7:
					/* y5 */
					tet[8] = TWOSQRT2_LO;
					tet[9] = TWOSQRT2_HI;
					/* y6 */
					tet[10] = TWOSQRT2_LO;
					tet[11] = TWOSQRT2_HI;
					break;
				default:
					printf("Somehow hit default case.\n");
					break;
			}
			
			/* y4 */
			/* diff = 2.0; */
			diff = TWOSQRT2_LO;
			if( n < 5 ) {
				tet[6] = diff;
				ROUND_UP;
				tet[7] = tet[9] + tet[11];
			} else {
				tet[6] = diff;
				ROUND_UP;
				tet[7] = tet[3] + tet[5];
			}
/*
			if( n==6 ) {
				tet[6] = 2.0;
				tet[7] = 3.9375;
			}
			if( n==7 ) {
				tet[6] = 3.9375;
				ROUND_UP;
				tet[7] = tet[3] + tet[5];
			}
*/
			/* Compute cell normalization information */
			max = 0.0;
			for( i=0; i<6; i++ ) {
				diff = tet[2*i+1] - tet[2*i];
				tetdiff[i] = diff;
				if( diff > max ) {
					max = diff;
				}
			}
			for( i=0; i<6; i++ ) {
				diff = tetdiff[i];
				if( diff > ATANERR ) {
					m = log(max/diff)/log(2.0) + 0.5;
					cell_info.diff_skip[i] = m;
				} else {
					cell_info.diff_skip[i] = 0;
				}
			}

			for( i=0; i<12; i++ )
				original_cell[i] = tet[i];
			
			ROUND_NEAR;
			for( i=0; i<12; i+=2 )
				printf("%0.18f\t%0.18f\n", tet[i], tet[i+1]);
			printf("\n");
/*			
			printf("sol_coeff = [%0.18g\t%0.18g]\n\n", sol_coeff[0], sol_coeff[1]);
*/		
/*
			if( n < 6 ) {
				verify_cell( tet, cell_info );
			} else if( n==6 ) {
				verify_cell2( tet, cell_info );
			} else {
				verify_cell3( tet, cell_info );
			}
*/
			if( n==1 || n==2 || n==4 ) {
				verify_cell3( tet, cell_info );
			} else {
				verify_cell4( tet, cell_info );
			}
			ROUND_NEAR;
			
			printf("cellcount    = %16.0f\t\t", cellcount);
			printf("maxdepth = %d\n", maxdepth);
			printf("discards2    = %16.0f\n", gmacount);
			printf("largecut     = %16.0f\t\t", dualcount);
			printf("discards     = %16.0f\n", vorcount);
			printf("smallcut     = %16.0f\t\t", sccutcount);
			printf("verifycount  = %16.0f\n", verifycount);
			printf("partialcount = %16.0f\t\t", partialcount);
			diff = cellcount - verifycount - gmacount - vorcount;
			printf("diff         = %16.0f\n", diff);
			if( !failed )
				printf("Verification succeeded.\n\n");
			else {
				printf("Verification FAILED.\n\n");
				onefailed++;
			}
			celltotal += cellcount;
		} /* end subcase loop */
	} /* end case loop */
		
	ROUND_NEAR;
	
	if( !onefailed )
		printf("All verifications succeeded.\n\n");
	else {
		if( onefailed == 1 )
			printf("One verification FAILED.\n\n");
		else
			printf("(%d) verifications FAILED.\n\n", onefailed);
	}

		
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
	printf("Elapsed time: %g seconds,\t\t", dt);
	printf("diff = %g\n", diff);
	printf("%.0f cells considered,\t\t", celltotal);
	printf("%g cells per second\n\n\n", celltotal/dt);
}


void verify_cell( double y[12], struct Cell_Info info )
{
	int i, j, bounds[6], n0, n1, n2, n3, n4, n5;
	double diff[6], t[12], val;
	double x[12], yp[12];
	double data[2], delta[2];
	
	
	/* If part of the verification has already failed, bail. */
	if( failed )
		return;
	
	cellcount += 1.0;
	
	/* Compute square of edge lengths */
	
	ROUND_UP;
	for( i=1; i<12; i+=2 ) {
		x[i] = y[i]*y[i];
	}
	ROUND_DOWN;
	for( i=0; i<12; i+=2 )
		x[i] = y[i]*y[i];
	if( x[9] > 8.0 )
		x[9] = 8.0;
	if( x[11] > 8.0 )
		x[11] = 8.0;

	/* First try to verify over the current cell.  */
	i_bigdelta_best( x, delta );
	
	if( delta[1] < 0.0 ) {	/* discard */
		vorcount += 1.0;
		return;
	}
		
	/* Score cell properly */
	/* i_vor0_y1y1( y, x, data ); */

	clear_vor0_y1y1( y, x, data );
	
	/* Attempt verification */	
	if( data[0] > slop ) {
		verifycount += 1.0;
		return;
	}
	if( data[1] < 0.0 && delta[0] > 0.0 ) {
		/* We are in big trouble. */
		failed = 1;
		printf("Cell failed:\n");
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", y[i], y[i+1]);
		printf("data = [%0.18g\t%0.18g]\n", data[0], data[1]);
		printf("delta = [%0.18g\t%0.18g]\n", delta[0], delta[1]);
		return;
	}
	
	/* if verification fails, try to bump off, or at least
		reduce dimension */

	for( i=0; i<12; i++ )
		yp[i] = y[i];	/* copy y */
	
	if( info.depth < MAXDEPTH && cellcount < countlimit ) {
		/* If all else fails, subdivide.  */
		info.depth++;
		if( info.depth > maxdepth )
			maxdepth = info.depth;
		ROUND_NEAR;
		j = 0;
		for( i=0; i<6; i++ ) {	/* scaled differentials */
			diff[i] = SUBFRAC*(yp[j+1] - yp[j]);
			if( diff[i] > 1.0e-14 )
				bounds[i] = SUBDIV;
			else {
				bounds[i] = 1;
				diff[i] = yp[j+1] - yp[j];
			}
			j += 2;
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
		printf("data = [%0.18g\t%0.18g]\n", data[0], data[1]);
		printf("delta = [%0.18g\t%0.18g]\n", delta[0], delta[1]);
		printf("y = \n");
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", y[i], y[i+1]);
	}
}


void verify_cell2( double y[12], struct Cell_Info info )
{
	int i, j, bounds[6], n0, n1, n2, n3, n4, n5;
	double diff[6], t[12], val;
	double x[12], yp[12];
	double data[2], data2[2], delta[2];
	
	/* If part of the verification has already failed, bail. */
	if( failed )
		return;
	
	cellcount += 1.0;
	
	/* Compute square of edge lengths */
	
	ROUND_UP;
	for( i=1; i<12; i+=2 ) {
		x[i] = y[i]*y[i];
	}
	ROUND_DOWN;
	for( i=0; i<12; i+=2 )
		x[i] = y[i]*y[i];
	if( x[9] > 8.0 )
		x[9] = 8.0;
	if( x[11] > 8.0 )
		x[11] = 8.0;

	/* First try to verify over the current cell.  */
	i_bigdelta_best( x, delta );
	
	if( delta[1] < 0.0 ) {	/* discard */
		vorcount += 1.0;
		return;
	}
		
	/* Score cell properly */
/*
	i_vor0_y1y1( y, x, data );
	i_vor0_partials( y, x, partials );
	clear_vor0_partials( y, x, partials );
*/
	superclear_vor0_y1( y, x, data2 );

	if( data2[0] > 0.0 ) {	/* discard */
		gmacount += 1.0;
		return;
	}
	
	if( data2[1] < 0.0 || info.depth > 6 ) {
		clear_vor0_y1y1( y, x, data );
		
		/* Attempt verification */	
		if( data[0] > slop ) {
			verifycount += 1.0;
			return;
		}
		if( data[1] < 0.0 && data2[1] < 0.0 && delta[0] > 0.0 ) {
			/* We are in big trouble. */
			failed = 1;
			printf("Cell failed:\n");
			for( i=0; i<12; i+=2 )
				printf("%0.18f\t%0.18f\n", y[i], y[i+1]);
			printf("data  = [%0.18g\t%0.18g]\n", data[0], data[1]);
			printf("data2 = [%0.18g\t%0.18g]\n", data2[0], data2[1]);
			printf("delta = [%0.18g\t%0.18g]\n", delta[0], delta[1]);
			return;
		}
	}
	
	/* if verification fails, try to bump off, or at least
		reduce dimension */

	for( i=0; i<12; i++ )
		yp[i] = y[i];	/* copy y */
	
	if( info.depth < MAXDEPTH && cellcount < countlimit ) {
		/* If all else fails, subdivide.  */
		info.depth++;
		if( info.depth > maxdepth )
			maxdepth = info.depth;
		ROUND_NEAR;
		j = 0;
		for( i=0; i<6; i++ ) {	/* scaled differentials */
			diff[i] = SUBFRAC*(yp[j+1] - yp[j]);
			if( diff[i] > 1.0e-14 )
				bounds[i] = SUBDIV;
			else {
				bounds[i] = 1;
				diff[i] = yp[j+1] - yp[j];
			}
			j += 2;
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
								verify_cell2( t, info );
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
		printf("data  = [%0.18g\t%0.18g]\n", data[0], data[1]);
		printf("data2 = [%0.18g\t%0.18g]\n", data2[0], data2[1]);
		printf("delta = [%0.18g\t%0.18g]\n", delta[0], delta[1]);
		printf("y = \n");
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", y[i], y[i+1]);
	}
}


void verify_cell3( double y[12], struct Cell_Info info )
{
	int i, j, bounds[6], n0, n1, n2, n3, n4, n5;
	double diff[6], t[12], val;
	double x[12], yp[12];
	double data[2], data2[2], delta[2], partials[12];
#if XPARS
	double oo2y[12];
#endif
	
	/* If part of the verification has already failed, bail. */
	if( failed )
		return;
	
	cellcount += 1.0;
	
	/* Compute square of edge lengths */
	
	ROUND_UP;
	for( i=1; i<12; i+=2 ) {
		x[i] = y[i]*y[i];
#if XPARS
		oo2y[i] = 0.5/y[i-1];
#endif
	}
	ROUND_DOWN;
	for( i=0; i<12; i+=2 ) {
		x[i] = y[i]*y[i];
#if XPARS
		oo2y[i] = 0.5/y[i+1];
#endif
	}
	if( x[9] > 8.0 )
		x[9] = 8.0;
	if( x[11] > 8.0 )
		x[11] = 8.0;

	/* First try to verify over the current cell.  */
	i_bigdelta_best( x, delta );
	
	if( delta[1] < 0.0 ) {	/* discard */
		vorcount += 1.0;
		return;
	}
		
	/* Score cell properly */
/*
	i_vor0_y1y1( y, x, data );
	i_vor0_partials( y, x, partials );
	clear_vor0_partials( y, x, partials );
*/
	clear_vor0_partials( y, x, partials );
	
	data2[0] = partials[0];
	data2[1] = partials[1];

	if( data2[0] > 0.0 || data2[1] < 0.0 ) {	/* discard */
		gmacount += 1.0;
		return;
	}

#if XPARS
	clear_vor0_x1x1( y, x, oo2y, data );
#else
	clear_vor0_y1y1( y, x, data );
#endif
	
	/* Attempt verification */	
	if( data[0] > slop ) {
		verifycount += 1.0;
		return;
	}
	if( data[1] < 0.0 && 
			(data2[1] < 0.0 || data2[0] > 0.0) && 
			delta[0] > 0.0 ) {
		/* We are in big trouble. */
		failed = 1;
		printf("Cell failed:\n");
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", y[i], y[i+1]);
		printf("data  = [%0.18g\t%0.18g]\n", data[0], data[1]);
		printf("data2 = [%0.18g\t%0.18g]\n", data2[0], data2[1]);
		printf("delta = [%0.18g\t%0.18g]\n", delta[0], delta[1]);
		return;
	}
	
	/* if verification fails, try to bump off, or at least
		reduce dimension */

	for( i=0; i<12; i++ )
		yp[i] = y[i];	/* copy y */
	
	if( info.depth < MAXDEPTH && cellcount < countlimit ) {
		/* If all else fails, subdivide.  */
		info.depth++;
		if( info.depth > maxdepth )
			maxdepth = info.depth;
		ROUND_NEAR;
		j = 0;
		
#if SMARTSUBDIV
		for( i=0; i<6; i++ ) {	/* differentials */
			diff[i] = yp[j+1] - yp[j];
			j += 2;
		}
		for( i=0; i<6; i++ ) {
			if( info.diff_skip[i] > 0 || diff[i] < 1.0e-14 ) { 
				/* skip */
				bounds[i] = 1;
				if( diff[i] > 1.0e-14 )
					info.diff_skip[i]--;
			} else {
			diff[i] *= SUBFRAC;
			bounds[i] = SUBDIV;
			}
		}
#else
		for( i=0; i<6; i++ ) {	/* scaled differentials */
			diff[i] = SUBFRAC*(yp[j+1] - yp[j]);
			if( diff[i] > 1.0e-14 )
				bounds[i] = SUBDIV;
			else {
				bounds[i] = 1;
				diff[i] = yp[j+1] - yp[j];
			}
			j += 2;
		}
#endif

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
								verify_cell3( t, info );
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
		printf("data  = [%0.18g\t%0.18g]\n", data[0], data[1]);
		printf("data2 = [%0.18g\t%0.18g]\n", data2[0], data2[1]);
		printf("delta = [%0.18g\t%0.18g]\n", delta[0], delta[1]);
		printf("y = \n");
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", y[i], y[i+1]);
	}
}


void verify_cell4( double y[12], struct Cell_Info info )
{
	int i, j, bounds[6], n0, n1, n2, n3, n4, n5;
	double diff[6], t[12], val;
	double x[12], yp[12];
	double data[2], data2[2], delta[2];
#if XPARS && !SUPERCLEAR
	double oo2y[12];
#endif
	
	/* If part of the verification has already failed, bail. */
	if( failed )
		return;
	
	cellcount += 1.0;
	
	/* Compute square of edge lengths */
	
	ROUND_UP;
	for( i=1; i<12; i+=2 ) {
		x[i] = y[i]*y[i];
#if XPARS && !SUPERCLEAR
		oo2y[i] = 0.5/y[i-1];
#endif
	}
	ROUND_DOWN;
	for( i=0; i<12; i+=2 ) {
		x[i] = y[i]*y[i];
#if XPARS && !SUPERCLEAR
		oo2y[i] = 0.5/y[i+1];
#endif
	}
	if( x[9] > 8.0 )
		x[9] = 8.0;
	if( x[11] > 8.0 )
		x[11] = 8.0;

	/* First try to verify over the current cell.  */
	i_bigdelta_best( x, delta );
	
	if( delta[1] < 0.0 ) {	/* discard */
		vorcount += 1.0;
		return;
	}
		
	/* Score cell properly */
/*
	i_vor0_y1y1( y, x, data );
	i_vor0_partials( y, x, partials );
	clear_vor0_partials( y, x, partials );
*/
	superclear_vor0_y1( y, x, data2 );

	if( data2[0] > 0.0 || data2[1] < 0.0 ) {	/* discard */
		gmacount += 1.0;
		return;
	}
	
#if XPARS

#if SUPERCLEAR
	superclear_vor0_x1x1( y, x, data );
#else
	clear_vor0_x1x1( y, x, oo2y, data );
#endif

#else
	clear_vor0_y1y1( y, x, data );
#endif
	
	/* Attempt verification */	
	if( data[0] > slop ) {
		verifycount += 1.0;
		return;
	}
	if( data[1] < 0.0 && 
			(data2[1] < 0.0 || data2[0] > 0.0) && 
			delta[0] > 0.0 ) {
		/* We are in big trouble. */
		failed = 1;
		printf("Cell failed:\n");
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", y[i], y[i+1]);
		printf("data  = [%0.18g\t%0.18g]\n", data[0], data[1]);
		printf("data2 = [%0.18g\t%0.18g]\n", data2[0], data2[1]);
		printf("delta = [%0.18g\t%0.18g]\n", delta[0], delta[1]);
		return;
	}
	
	/* if verification fails, try to bump off, or at least
		reduce dimension */

	for( i=0; i<12; i++ )
		yp[i] = y[i];	/* copy y */
	
	if( info.depth < MAXDEPTH && cellcount < countlimit ) {
		/* If all else fails, subdivide.  */
		info.depth++;
		if( info.depth > maxdepth )
			maxdepth = info.depth;
		ROUND_NEAR;
		j = 0;
		
#if SMARTSUBDIV
		for( i=0; i<6; i++ ) {	/* differentials */
			diff[i] = yp[j+1] - yp[j];
			j += 2;
		}
		for( i=0; i<6; i++ ) {
			if( info.diff_skip[i] > 0 || diff[i] < 1.0e-14 ) { 
				/* skip */
				bounds[i] = 1;
				if( diff[i] > 1.0e-14 )
					info.diff_skip[i]--;
			} else {
			diff[i] *= SUBFRAC;
			bounds[i] = SUBDIV;
			}
		}
#else
		for( i=0; i<6; i++ ) {	/* scaled differentials */
			diff[i] = SUBFRAC*(yp[j+1] - yp[j]);
			if( diff[i] > 1.0e-14 )
				bounds[i] = SUBDIV;
			else {
				bounds[i] = 1;
				diff[i] = yp[j+1] - yp[j];
			}
			j += 2;
		}
#endif

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
								verify_cell4( t, info );
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
		printf("data  = [%0.18g\t%0.18g]\n", data[0], data[1]);
		printf("data2 = [%0.18g\t%0.18g]\n", data2[0], data2[1]);
		printf("delta = [%0.18g\t%0.18g]\n", delta[0], delta[1]);
		printf("y = \n");
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", y[i], y[i+1]);
	}
}


void testquoin( void )
{
	int i;
	double abc[6], xyz[6];
	double quo[2], quo_pars[4], quo_secpars[6];
	
	while( 1 ) {
		ROUND_NEAR;
		printf("Enter interval values for a, b, c:\n");
		for( i=0; i<6; i++ )
			scanf("%lf", abc + i);
		printf("Got abc = \n");
		for( i=0; i<6; i+=2 )
			printf("[%.18f, %.18f]\n", abc[i], abc[i+1]);
		ROUND_DOWN;
		for( i=0; i<6; i+=2 )
			xyz[i] = abc[i]*abc[i];
		ROUND_UP;
		for( i=1; i<6; i+=2 )
			xyz[i] = abc[i]*abc[i];
		i_quoin_all( abc, xyz, quo, quo_pars, quo_secpars );
		ROUND_NEAR;
		printf("quoin =    [%.18f, %.18f]\n", quo[0], quo[1]);
		printf("quoin_a =  [%.18f, %.18f]\n", quo_pars[0], quo_pars[1]);
		printf("quoin_b =  [%.18f, %.18f]\n", quo_pars[2], quo_pars[3]);
		printf("quoin_aa = [%.18f, %.18f]\n", quo_secpars[0], quo_secpars[1]);
		printf("quoin_ab = [%.18f, %.18f]\n", quo_secpars[2], quo_secpars[3]);
		printf("quoin_bb = [%.18f, %.18f]\n", quo_secpars[4], quo_secpars[5]);
	}
}


/*
1.10234 1.32567 1.4142
1.01234 1.21432 1.255

1.10234 1.10234
1.32567 1.32567
1.4142  1.4142

1.01234 1.01234 
1.21432 1.21432 
1.255 1.255

*/


void testvor0( void )
{
	int i, j;
	double y[12], x[12], delta[2], sqrtdelta[2];
	double vor0_val[2], delta_part[12];
	double part[12], sph_pars[12];
	double dih_xij[6][12], sol_ij[6][12];
	double apart[12], delta32[2], clear_ij[6][12];
	double clear_i[12], oo2y[12];
	double temp[2], temp2[2], sum[2];
	double facex[6], u1262u1352[2];
	
	apart[0] = 0.0;
	part[0] = 0.0;
	temp[0] = 0.0;
	temp2[0] = 0.0;
	sum[0] = 0.0;
	
	i_init();
	appendix_init();
	apart[0] = TWO51_LO;
	apart[1] = TWO51_HI;
	sean_init( apart );
/*	
	ROUND_DOWN;
	t0_val[0] = sqrt(2.0);
	ROUND_UP;
	t0_val[1] = sqrt(2.0);
	phi_fun( t0_val, t0_val, phi0_val );
	sol_coeff[0] = phi0_val[0];
	sol_coeff[1] = phi0_val[1];
*/
/*
	ROUND_DOWN;
	sol_coeff[0] = phi0_val[0] - zetapt_val[1];
	ROUND_UP;
	sol_coeff[1] = phi0_val[1] - zetapt_val[0];
*/
	while( 1 ) {
		ROUND_NEAR;
		printf("Enter interval values for y:\n");
		for( i=0; i<12; i++ )
			scanf("%lf", y + i);
		printf("Got y = \n");
		for( i=0; i<12; i+=2 )
			printf("[%.18f, %.18f]\n", y[i], y[i+1]);
		ROUND_DOWN;
		for( i=0; i<12; i+=2 )
			x[i] = y[i]*y[i];
		ROUND_UP;
		for( i=1; i<12; i+=2 )
			x[i] = y[i]*y[i];
		ROUND_DOWN;
		for( i=0; i<12; i+=2 )
			oo2y[i] = 0.5/y[i+1];
		ROUND_UP;
		for( i=1; i<12; i+=2 )
			oo2y[i] = 0.5/y[i-1];
		
		/* 1 2 6 */
		for( i=0; i<2; i++ ) {
			facex[i    ] = x[     i];
			facex[2 + i] = x[2  + i];
			facex[4 + i] = x[10 + i];
		}
		super_u( facex, temp );
		i_mult( temp, temp, temp2 );
		/* 1 3 5 */
		for( i=0; i<2; i++ ) {
			facex[2 + i] = x[4 + i];
			facex[4 + i] = x[8 + i];
		}
		super_u( facex, temp );
		i_mult( temp, temp, sum );
		i_mult( temp2, sum, u1262u1352 );

		i_bigdelta_partials_best( x, delta, delta_part );
		ROUND_NEAR;
		printf("delta = [%.18f, %.18f]\n", delta[0], delta[1]);
		i_sqrt( delta, sqrtdelta );
		i_mult( delta, sqrtdelta, delta32 );
		i_vor0( y, x, sqrtdelta, vor0_val );
/*
		printf("delta32 = [%.18f, %.18f]\n", delta32[0], delta32[1]);
		ROUND_NEAR;
		printf("vor0 =      [%.18f, %.18f]\n", vor0_val[0], vor0_val[1]);
*/
/*
		i_vor_trunc( y, x, sqrtdelta, vor0_val );
		ROUND_NEAR;
		printf("vor_trunc = [%.18f, %.18f]\n", vor0_val[0], vor0_val[1]);
*/
/*
		i_crad( y, sqrtdelta );
		ROUND_NEAR;
		printf("crad  =     [%.18f, %.18f]\n", sqrtdelta[0], sqrtdelta[1]);
		i_squarelen( sqrtdelta, delta );
		ROUND_NEAR;
		printf("crad2 =     [%.18f, %.18f]\n", delta[0], delta[1]);
		ROUND_NEAR;
		for( i=1; i<4; i++ ) {
			j = 2*(i-1);
			printf("x[%d] =     [%.18f, %.18f]\n", i, 0.25*x[j], 0.25*x[j+1]);
		}
*/
/*
		crad3len_pars( y, x, part );
		ROUND_NEAR;
		printf("crad3len_pars = \n");
		for( j=0; j<6; j+=2 ) {
			printf("[%.18f, %.18f]\n", part[j], part[j+1]);
		}
		s_crad3len_pars( y, x, part );
		ROUND_NEAR;
		printf("s_crad3len_pars = \n");
		for( i=0; i<6; i+=2 ) {
			printf("[%.18f, %.18f]\n", part[i], part[i+1]);
		}
		crad3len2_xpars( x, part );
		ROUND_NEAR;
		printf("crad3len2_xpars = \n");
		for( j=0; j<6; j+=2 ) {
			printf("[%.18f, %.18f]\n", part[j], part[j+1]);
		}
		s_crad3len2_xpars( x, part );
		ROUND_NEAR;
		printf("s_crad3len2_xpars = \n");
		for( i=0; i<6; i+=2 ) {
			printf("[%.18f, %.18f]\n", part[i], part[i+1]);
		}
*/
		j = 0;
		sph_pars[0] = 0.0;
/*
		i_vor0_partials( y, x, part );
		ROUND_NEAR;
		printf("i_vor0_partials = \n");
		for( j=0; j<12; j+=2 ) {
			printf("[%.18f, %.18f]\n", part[j], part[j+1]);
		}
*/
/*
		s_solid_partials( y, delta, sqrtdelta, delta_part,
			sph_pars );
		clear_sol_i( y, x, apart );
		for( j=0; j<12; j+=2 ) {
			i_div( apart + j, delta32, part + j );
		}
		ROUND_NEAR;
		printf("old sol_i = \n");
		for( j=0; j<12; j+=2 ) {
			printf("[%.18f, %.18f]\n", sph_pars[j], sph_pars[j+1]);
		}
		printf("sol_i = \n");
		for( j=0; j<12; j+=2 ) {
			printf("[%.18f, %.18f]\n", part[j], part[j+1]);
		}
		for( j=0; j<12; j+=2 ) {
			i_sub( part + j, sph_pars + j, apart + j );
		}
		printf("diff_i = \n");
		for( j=0; j<12; j+=2 ) {
			printf("[%.18f, %.18f]\n", apart[j], apart[j+1]);
		}
*/
/*
		i_dih_xpars( x, sph_pars );
		clear_dih_i( y, x, apart );
		for( j=0; j<12; j+=2 ) {
			i_div( apart + j, delta32, part + j );
		}
		ROUND_NEAR;
		printf("old dih_i = \n");
		for( j=0; j<12; j+=2 ) {
			printf("[%.18f, %.18f]\n", sph_pars[j], sph_pars[j+1]);
		}
		printf("dih_i = \n");
		for( j=0; j<12; j+=2 ) {
			printf("[%.18f, %.18f]\n", part[j], part[j+1]);
		}
		for( j=0; j<12; j+=2 ) {
			i_sub( part + j, sph_pars + j, apart + j );
		}
		printf("diff_i = \n");
		for( j=0; j<12; j+=2 ) {
			printf("[%.18f, %.18f]\n", apart[j], apart[j+1]);
		}
*/
		i_vor0_partials( y, x, sph_pars );
		i_vor0_y1y1( y, x, sph_pars );
		clear_vor0_partials( y, x, apart );
		for( j=0; j<12; j+=2 ) {
			i_div( apart + j, delta32, part + j );
		}
/*
		ROUND_NEAR;
		printf("old vor0_i = \n");
		for( j=0; j<12; j+=2 ) {
			printf("[%.18f, %.18f]\n", sph_pars[j], sph_pars[j+1]);
		}
		printf("vor0_i = \n");
		for( j=0; j<12; j+=2 ) {
			printf("[%.18f, %.18f]\n", part[j], part[j+1]);
		}
		for( j=0; j<12; j+=2 ) {
			i_sub( part + j, sph_pars + j, apart + j );
		}
		printf("diff_i = \n");
		for( j=0; j<12; j+=2 ) {
			printf("[%.18f, %.18f]\n", apart[j], apart[j+1]);
		}
*/
		clear_vor0_xpars( y, x, oo2y, sph_pars );
		clear_vor0_partials( y, x, apart );
		sp_xtoypars( y, sph_pars, part );
/*
		ROUND_NEAR;
		printf("clear_vor0_i = \n");
		for( j=0; j<12; j+=2 ) {
			printf("[%.18f, %.18f]\n", apart[j], apart[j+1]);
		}
		printf("clear_vor0_xi = \n");
		for( j=0; j<12; j+=2 ) {
			printf("[%.18f, %.18f]\n", sph_pars[j], sph_pars[j+1]);
		}
		for( j=0; j<12; j+=2 ) {
			i_sub( part + j, apart + j, sph_pars + j );
		}
		printf("diff_i = \n");
		for( j=0; j<12; j+=2 ) {
			printf("[%.18f, %.18f]\n", sph_pars[j], sph_pars[j+1]);
		}
*/
/*
		s_solid_partials( y, delta, sqrtdelta, delta_part,
			sph_pars );
		trunc_pars( y, x, sph_pars, part );
		ROUND_NEAR;
		printf("trunc_pars = \n");
		for( j=0; j<12; j+=2 ) {
			printf("[%.18f, %.18f]\n", part[j], part[j+1]);
		}
*/
/*
		i_dih( x, part );
		ROUND_NEAR;
		printf("i_dih = [%.18f, %.18f]\n", part[0], part[1]);
		i_dih_xpars( x, part );
		ROUND_NEAR;
		printf("i_dih_xpars = \n");
		for( j=0; j<12; j+=2 ) {
			printf("[%.18f, %.18f]\n", part[j], part[j+1]);
		}
*/
/*
		i_dih_xpars( x, sph_pars );
		clear_dih_i( y, x, apart );
		for( j=0; j<12; j+=2 ) {
			i_div( apart + j, delta32, part + j );
		}
		ROUND_NEAR;
		printf("old dih_i = \n");
		for( j=0; j<12; j+=2 ) {
			printf("[%.18f, %.18f]\n", sph_pars[j], sph_pars[j+1]);
		}
		printf("dih_i = \n");
		for( j=0; j<12; j+=2 ) {
			printf("[%.18f, %.18f]\n", part[j], part[j+1]);
		}
*/
/*
		i_dih_xpars( x, part );
		ROUND_NEAR;
		printf("i_dih_xpars = \n");
		for( j=0; j<12; j+=2 ) {
			printf("[%.18f, %.18f]\n", part[j], part[j+1]);
		}
		old_dih_xpars( x, part );
		ROUND_NEAR;
		printf("old_dih_xpars = \n");
		for( j=0; j<12; j+=2 ) {
			printf("[%.18f, %.18f]\n", part[j], part[j+1]);
		}
*/
/*
		sp_dihsecparbds( x, dih_xij );
		sp_solsecyparbds( y, x, dih_xij );
		ROUND_NEAR;
		printf("old sol_ij = \n");
		for( i=0; i<6; i++ ) {
			for( j=2*i; j<12; j+=2 ) {
				printf("[%.18f, %.18f]\n", dih_xij[i][j],
					dih_xij[i][j+1] );
			}
			printf("\n");
		}
		clear_sol_ij( y, x, clear_i, clear_ij );
		for( i=0; i<6; i++ ) {
			for( j=2*i; j<12; j+=2 ) {
				i_div( clear_ij[i] + j, delta32, sol_ij[i] + j );
			}
		}
		ROUND_NEAR;
		printf("sol_ij = \n");
		for( i=0; i<6; i++ ) {
			for( j=2*i; j<12; j+=2 ) {
				printf("[%.18f, %.18f]\n", sol_ij[i][j],
					sol_ij[i][j+1] );
			}
			printf("\n");
		}
		for( i=0; i<6; i++ ) {
			for( j=2*i; j<12; j+=2 ) {
				i_sub( sol_ij[i] + j, dih_xij[i] + j, clear_ij[i] + j );
			}
		}
		ROUND_NEAR;
		printf("diff_ij = \n");
		for( i=0; i<6; i++ ) {
			for( j=2*i; j<12; j+=2 ) {
				printf("[%.18f, %.18f]\n", clear_ij[i][j],
					clear_ij[i][j+1] );
			}
			printf("\n");
		}
*/
		sp_dihsecparbds( x, dih_xij );
		clear_dih_ij( y, x, clear_i, clear_ij );
		for( i=0; i<6; i++ ) {
			for( j=2*i; j<12; j+=2 ) {
				i_div( clear_ij[i] + j, delta32, sol_ij[i] + j );
			}
		}
/*
		ROUND_NEAR;
		printf("old dih_ij = \n");
		for( i=0; i<6; i++ ) {
			for( j=2*i; j<12; j+=2 ) {
				printf("[%.18f, %.18f]\n", dih_xij[i][j],
					dih_xij[i][j+1] );
			}
			printf("\n");
		}
		ROUND_NEAR;
		printf("dih_ij = \n");
		for( i=0; i<6; i++ ) {
			for( j=2*i; j<12; j+=2 ) {
				printf("[%.18f, %.18f]\n", sol_ij[i][j],
					sol_ij[i][j+1] );
			}
			printf("\n");
		}
		for( i=0; i<6; i++ ) {
			for( j=2*i; j<12; j+=2 ) {
				i_sub( sol_ij[i] + j, dih_xij[i] + j, clear_ij[i] + j );
			}
		}
		ROUND_NEAR;
		printf("diff_ij = \n");
		for( i=0; i<6; i++ ) {
			for( j=2*i; j<12; j+=2 ) {
				printf("[%.18f, %.18f]\n", clear_ij[i][j],
					clear_ij[i][j+1] );
			}
			printf("\n");
		}
*/
/*
		ROUND_NEAR;
		printf("i_dih_x1x1 = [%.18f, %.18f]\n", dih_xij[0][0], dih_xij[0][1]);
*/
/*
		i_betadihsum_y1y1( y, x, part );
		ROUND_NEAR;
		printf("old i_betadihsum_y1y1 = [%.18f, %.18f]\n", part[0], part[1]);
		clear_betadihsum_y1y1( y, x, clear_i );
		i_div( clear_i, delta32, part );
		ROUND_NEAR;
		printf("i_betadihsum_y1y1     = [%.18f, %.18f]\n", part[0], part[1]);
*/
/*
		i_vor0_y1y1( y, x, part );
		ROUND_NEAR;
		printf("old i_vor0_y1y1 = [%.18f, %.18f]\n", part[0], part[1]);
		clear_vor0_y1y1( y, x, clear_i );
		i_div( clear_i, delta32, part );
		ROUND_NEAR;
		printf("i_vor0_y1y1     = [%.18f, %.18f]\n", part[0], part[1]);
*/
		clear_vor0_partials( y, x, apart );
		i_vor0_partials( y, x, part );
/*
		ROUND_NEAR;
		printf("clear_vor0_y1 = [%.18f, %.18f]\n", apart[0], apart[1]);
		printf("i_vor0_y1     = [%.18f, %.18f]\n", part[0], part[1]);
		clear_vor0_y1y1( y, x, part );
		i_vor0_y1y1( y, x, apart );
		ROUND_NEAR;
		printf("clear_vor0_y1y1 = [%.18f, %.18f]\n", part[0], part[1]);
		printf("i_vor0_y1y1     = [%.18f, %.18f]\n", apart[0], apart[1]);
*/
/*
		i_delta4( x, apart );
		sp_dih_a( x, part );
		ROUND_NEAR;
		printf("i_delta4 = [%.18f, %.18f]\n", apart[0], apart[1]);
		printf("sp_dih_a = [%.18f, %.18f]\n", part[0], part[1]);
*/
/*
		i_sol_y1y1( y, x, part );
		ROUND_NEAR;
		printf("i_sol_y1y1 = [%.18f, %.18f]\n", part[0], part[1]);
		i_betadihsum_y1y1( y, x, part );
		ROUND_NEAR;
		printf("i_betadihsum_y1y1 = [%.18f, %.18f]\n", part[0], part[1]);
		i_quoinsum_y1y1( y, x, part );
		ROUND_NEAR;
		printf("i_quoinsum_y1y1 = [%.18f, %.18f]\n", part[0], part[1]);
		i_vor0_y1y1( y, x, part );
		ROUND_NEAR;
		printf("i_vor0_y1y1 = [%.18f, %.18f]\n", part[0], part[1]);

		ROUND_DOWN;
		sol_coeff[0] = phi0_val[0] - zetapt_val[1];
		ROUND_UP;
		sol_coeff[1] = phi0_val[1] - zetapt_val[0];

		i_vor0_y1y1( y, x, part );
		ROUND_NEAR;
		printf("-i_tau0_y1y1 = [%.18f, %.18f]\n", part[0], part[1]);

		sol_coeff[0] = phi0_val[0];
		sol_coeff[1] = phi0_val[1];
*/
		clear_vor0_partials( y, x, sph_pars );
		clear_vor0_x1x1( y, x, oo2y, part );
		clear_vor0_y1y1( y, x, apart );
/*
		clear_vor0_partials( y, x, sph_pars );
		clear_vor0_x1x1( y, x, oo2y, part );
		clear_vor0_y1y1( y, x, apart );
		i_mult( sph_pars, oo2y, temp );
		ROUND_DOWN;
		temp2[0] = apart[0] - 2.0*temp[1];
		sum[0] = oo2y[0]*oo2y[0];
		ROUND_UP;
		temp2[1] = apart[1] - 2.0*temp[0];
		sum[1] = oo2y[1]*oo2y[1];
		i_mult( sum, temp2, temp );
		i_sub( temp, part, temp2 );
		ROUND_NEAR;
		printf("clear_vor0_y1y1 = [%.18f, %.18f]\n", apart[0], apart[1]);
		printf("clear_vor0_x1x1 = [%.18f, %.18f]\n", part[0], part[1]);
		printf("diff = [%.18g, %.18g]\n", temp2[0], temp2[1]);
*/
		superclear_vor0_y1( y, x, part );
		ROUND_NEAR;
		printf("superclear_vor0_y1 = [%.18f, %.18f]\n", part[0], part[1]);
		clear_vor0_x1x1( y, x, oo2y, part );
		clear_vor0_y1y1( y, x, apart );
		ROUND_NEAR;
		printf("clear_vor0_y1y1 = [%.18f, %.18f]\n", apart[0], apart[1]);
		printf("clear_vor0_x1x1 = [%.18f, %.18f]\n", part[0], part[1]);
		superclear_vor0_x1x1( y, x, part );
		clear_vor0_x1x1( y, x, oo2y, apart );
		i_div( part, u1262u1352, sph_pars );
		temp2[0] = 2.0/3.0*DOCT_LO;
		temp2[1] = 2.0/3.0*DOCT_HI;
		i_mult( temp2, sph_pars, temp );
		ROUND_NEAR;
		printf("superclear_vor0_x1x1 = [%.18f, %.18f]\n", part[0], part[1]);
		printf("clear_vor0_x1x1 = [%.18f, %.18f]\n", apart[0], apart[1]);
		printf("superclear_mod  = [%.18f, %.18f]\n", temp[0], temp[1]);
		i_vor0_partials( y, x, part );
		i_vor0_y1y1( y, x, apart );
		ROUND_NEAR;
		printf("vor0_y1 = [%.18f, %.18f]\n", part[0], part[1]);
		printf("vor0_y1y1 = [%.18f, %.18f]\n", apart[0], apart[1]);
		i_mult( y, apart, temp );
		i_sub( temp, part, apart );
		i_div( apart, x, part );
		i_div( part, y, temp );
		part[0] = 0.25*temp[0];
		part[1] = 0.25*temp[1];
		ROUND_NEAR;
		printf("vor0_x1x1 = [%.18f, %.18f] (from y partials)\n", 
			part[0], part[1]);
	}
}


void sean_init( double sean_val[2] )
{
	int i;
	double dt, temp[2], temp2[2], tet[12];
			
	/* initialize local constants */
	t0_val[0] = 0.5*sean_val[0];
	t0_val[1] = 0.5*sean_val[1];
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
	
	/* sol_coeff */
	ROUND_DOWN;
	sol_coeff[0] = phi0_val[0];
	ROUND_UP;
	sol_coeff[1] = phi0_val[1];
	
	/* compute super_coeff */
	/* super_coeff = 3 phi0/doct */
	temp[0] = DOCT_LO;
	temp[1] = DOCT_HI;
	i_div( sol_coeff, temp, temp2 );
	temp[0] = 3.0;
	temp[1] = 3.0;
	i_mult( temp, temp2, super_coeff );
}


/*
2.029102 		2.029102 
2.0341248 	2.0341248 
2.011234 		2.011234 
2.023309 		2.023309 
2.0123 			2.0123 
2.83982134 	2.83982134 

2.029102 		2.029102 
2.0341248 	2.0341248 
2.011234 		2.011234 
2.023309 		2.023309 
2.0123 			2.0123 
2.73982134 	2.73982134 

2.41				2.41
2.2291			2.2291
2.23412			2.23412
2.54				2.54
2.42331			2.42331
2.2123			2.2123

2.041				2.041
2.0291			2.0291
2.03412			2.03412
2.154				2.154
2.02331			2.02331
2.0123			2.0123

2.93982134 	2.93982134 
2.029102 		2.029102 
2.0341248 	2.0341248 
2.011234 		2.011234 
2.023309 		2.023309 
2.0123 			2.0123 

2.93982134 	2.9398213401
2.029102 		2.0291020001
2.0341248 	2.0341248001
2.011234 		2.0112340001
2.023309 		2.0233090001
2.0123 			2.0123000001

2.000062255859375249    2.000124511718750053
2.000186767578124858    2.000249023437500107
2.0    2.000062255859375249
3.463623046875    3.4638671875
2.0    2.0
2.0    2.0

2.0    2.0
2.0    2.0
2.0    2.0
3.975239257812499538    3.9752392578125
2.509999999999999787    2.510000000000000231
2.828427124746189847    2.828427124746190291

2.0 2.0
2.0 2.0
2.0 2.0
3.46410153 3.46410154
2.0 2.0
2.0 2.0

2.0    2.000000000118743682
2.0    2.000000000118743682
2.0    2.000000000118743682
2.0    2.000000000703148650
2.828427124746189847    2.828427124746190291
2.828427124746189847    2.828427124746190291

2.0    2.0
2.0    2.0
2.0    2.0
2.0    2.0
2.828427124746189847    2.828427124746190291
2.828427124746189847    2.828427124746190291

2.015937499999999716    2.031874999999999876
2.079687499999999911    2.095625000000000071
2.494062500000000071    2.510000000000000231
2.0    2.06250
2.0    2.0
2.0    2.0

2.0    2.000000000118743682
2.0    2.000000000118743682
2.0    2.000000000118743682
3.999999999543652596    4.000000000246800802
2.828427124746189847    2.828427124746190291
2.828427124746189847    2.828427124746190291

2.0    2.000000000000463629
2.034091385337014835    2.034091385337478464
2.254999999999536264    2.254999999999999893
3.9375    3.937500000000984546
2.828427124746189847    2.828427124746190291
2.828427124746189847    2.828427124746190291

2.041				2.041
2.0291			2.0291
2.03412			2.03412
2.154				2.154
2.52331			2.52331
2.5123			2.5123

2.041				2.0410001
2.0291			2.0291001
2.03412			2.0341201
2.154				2.1540001
2.52331			2.5233101
2.5123			2.5123001

2.041				2.041001
2.0291			2.029101
2.03412			2.034121
2.154				2.154001
2.52331			2.523311
2.5123			2.512301

Date: Thu, 21 May 1998 11:07:07 -0400 (EDT)
From: Tom Hales <hales@math.lsa.umich.edu>
To: Samuel Ferguson <samf@math.lsa.umich.edu>
cc: "Thomas C. Hales" <hales@math.lsa.umich.edu>
Subject: vor(sqrt2)
Message-ID: <Pine.SOL.3.96k.rh.980521102401.557A-100000@diefledermaus.math.lsa.u
mich.edu>
MIME-Version: 1.0
Content-Type: TEXT/PLAIN; charset=US-ASCII
Content-Length: 97


{2.41, 2.2291, 2.23412, 2.54, 2.42331, 2.2123}

vor truncated at sqrt2 = -0.2588174731801107

*/
