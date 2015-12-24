/* octahedra.c, by Samuel Ferguson, (c) 1998. */

#include "system_headers.h"
#include "interval.h"
#include "i_sphere.h"
#include "i_bounds.h"
#include "i_taylor.h"
#include "macros.h"


#define MAXDEPTH	32			/* Max recursion depth (was 10) */

#define SUBDIV		2			/* (4) 	number of subdivisions (1-d) */
#define SUBFRAC	0.5		/* (0.25) 	1/SUBDIV (make this exact . . .) */

#define PARTIALS	1

#define COUNTLIMIT	2.0e7		/* 4.0e7 */

#define INTERACTIVE	0
#define BADCHECK		0


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
int failed;
double slop;
int taylor_rel[7];
double taylor_relconst[7];
double y_coefficients[6];
int	faceinfo[2];
double original_cell[12];

struct Cell_Info {
	int depth; 
	};

/* Prototypes */

void verifyocta( void );
void verify_cell_original( double tet[12], struct Cell_Info info );
void verify_cell( double tet[12], struct Cell_Info info );
void verify_bdry126( double y[12], struct Cell_Info info );
void verify_bdry135( double y[12], struct Cell_Info info );
void verify_bdry( double y[12], struct Cell_Info info );
void set_relations( int casenumber );
int tomtomynum( int num );
int mynumtotom( int mynum );
void set_domain( int casenum );
void fix_domain( double tet[12], double fixtet[12] );
void octa_score( double y[12], double x[12],
	double out[2], double out_partials[12] );
void testvalues( void );
void constraint_functions( double x[12], double out[4] );
int find_lambdas( double tet[12], double out[2] );
void intersect( double xyz1[3], double xyz2[3], 
	double out[2] );

/* Code */

void main( void )
{
	ROUND_NEAR;

	printf("Welcome to sphere packing routines, relation testing division.\n");
	verifyocta();
	/*
	testvalues();
	verifyocta();
	*/
}


void verifyocta( void )
{
	int i, n, casenumber, mynumber, casebound, subcase;
	double tet[12];
	double diff, dt;
	time_t tstart, tstop;
	clock_t cstart, cstop;
	char *charptr;
	struct Cell_Info cell_info;
	

	printf("Welcome to Octahedra.  \n");
#if INTERACTIVE
	mynumber = -1;
	while( mynumber == -1 ) {
		printf("Enter case number:  ");
		scanf("%d", &casenumber );
		printf("\nCase = %d\n", casenumber );
		mynumber = tomtomynum( casenumber );
		if( mynumber == -1 )
			printf("Invalid case.  Try again.\n");
	}
#endif
	printf("Enter slop:  ");
	scanf("%lf", &slop);
	printf("\t\tSlop = %g\n", slop );
	
	i_init();	/* initialize interval constants */
	
	ROUND_NEAR;

	time(&tstart);
	cstart = clock();
	
	charptr = ctime( &tstart );
	printf("Starting time:  \t");
	puts( charptr );
	
	celltotal = 0.0;
	
#if INTERACTIVE
	casebound = 2;
#else
	casebound = 19;
#endif
	for( n=1; n < casebound; n++ ) { 

#if !INTERACTIVE
		casenumber = mynumtotom( n );
#endif

		/* printf("Starting case %d:\n", n); */
		printf("Starting case %d:\n", casenumber);
				
		mynumber = tomtomynum( casenumber );
		set_relations( mynumber );
		/* set_domain( casenumber ); */
		
		for( i=0; i<12; i++ )
			tet[i] = original_cell[i];
		
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", tet[i], tet[i+1]);
		printf("\n");
		
		for( subcase=1; subcase < 5; subcase ++ ) {
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
			
			for( i=0; i<12; i++ )
				tet[i] = original_cell[i];
		
			printf("Subcase %d:\n", subcase);
			switch( subcase ) {
				case 1:
					printf("\tNo hard boundary constraints.\n");
					verify_cell( tet, cell_info );
					break;
				case 2:
					printf("\tBoundary:  eta2(1 2 6) = 2.\n");
					verify_bdry126( tet, cell_info );
					break;
				case 3:
					printf("\tBoundary:  eta2(1 3 5) = 2.\n");
					verify_bdry135( tet, cell_info );
					break;
				case 4:
					printf("\tBoundary:  eta2(1 2 6) = eta2(1 3 5) = 2.\n");
					verify_bdry( tet, cell_info );
					break;
				default:
					printf("Somehow hit default case.\n");
					break;
			}
			
			printf("cellcount    = %16.0f\t\t", cellcount);
			printf("maxdepth = %d\n", maxdepth);
			printf("partialpushes= %16.0f\n", gmacount);
			printf("largecut     = %16.0f\t\t", dualcount);
			printf("discards     = %16.0f\n", vorcount);
			printf("smallcut     = %16.0f\t\t", sccutcount);
			printf("verifycount  = %16.0f\n", verifycount);
			printf("partialcount = %16.0f\t\t", partialcount);
			diff = cellcount - verifycount - partialcount - dualcount - sccutcount;
			printf("diff         = %16.0f\n", diff);
			if( !failed )
				printf("Verification succeeded.\n\n");
			else
				printf("Verification FAILED.\n\n");
			celltotal += cellcount;
		}	/* end subcase loop */
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


void verify_cell_original( double y[12], struct Cell_Info info )
{
	int i, j, bounds[6], n0, n1, n2, n3, n4, n5;
	int tight126, tight135, okaylist[6];
	double diff[6], t[12], val;
	double x[12], yp[12], facex[6];
	double data[2], rel_part[12];
	double crad1262[2], crad1352[2];
	
	
	/* If part of the verification has already failed, bail. */
	if( failed )
		return;
	
	cellcount += 1.0;
	
	/* Compute square of edge lengths */
	
	ROUND_UP;
	for( i=1; i<12; i+=2 )
		x[i] = y[i]*y[i];
	ROUND_DOWN;
	for( i=0; i<12; i+=2 )
		x[i] = y[i]*y[i];
	if( x[1] > 8.0 )
		x[1] = 8.0;

	/* First try to verify over the current cell.  */

	/* Check face constraints. */
	tight126 = 1;
	tight135 = 1;
	
	/* face:  (1 2 6) -> (0 2 10) */
	facex[0] = x[0];
	facex[1] = x[1];
	facex[2] = x[2];
	facex[3] = x[3];
	facex[4] = x[10];
	facex[5] = x[11];
	s_crad3x2( facex, crad1262 );
	
	if( faceinfo[0] < 1 ) {	/* should be small */
		if( crad1262[0] > 2.0 ) {
			dualcount += 1.0;
			return;
		}
		if( crad1262[1] < 2.0 ) {
			tight126 = 0;
		}
	} else {	/* should be large */
		if( crad1262[1] < 2.0 ) {
			sccutcount += 1.0;
			return;
		}
		if( crad1262[0] > 2.0 ) {
			tight126 = 0;
		}
	}
	
	/* face:  (1 3 5) -> (0 4 8) */
	facex[2] = x[4];
	facex[3] = x[5];
	facex[4] = x[8];
	facex[5] = x[9];
	s_crad3x2( facex, crad1352 );

	if( faceinfo[1] < 1 ) {	/* should be small */
		if( crad1352[0] > 2.0 ) {
			dualcount += 1.0;
			return;
		}
		if( crad1352[1] < 2.0 ) {
			tight135 = 0;
		}
	} else {	/* should be large */
		if( crad1352[1] < 2.0 ) {
			sccutcount += 1.0;
			return;
		}
		if( crad1352[0] > 2.0 ) {
			tight135 = 0;
		}
	}
	
	/* Score cell properly */
	octa_score( y, x, data, rel_part );
		
	
	/* Attempt verification */	
	if( data[1] < slop ) {
		verifycount += 1.0;
		return;
	}
#if BADCHECK
	if( data[0] > slop ) {
		/* We are in big trouble. */
		failed = 1;
		printf("Cell failed:\n");
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", y[i], y[i+1]);
		printf("data = [%0.18g\t%0.18g]\n", data[0], data[1]);
		return;
	}
#endif
	
	/* if verification fails, try to bump off, or at least
		reduce dimension */

	for( i=0; i<12; i++ )
		yp[i] = y[i];	/* copy y */
	
#if PARTIALS
	
	/* Score cell properly (find min score) */
	
	/* t contains the relation partials */
	ROUND_NEAR;
	j = 0;
	for( i=0; i<6; i++ ) {
		diff[i] = yp[j+1] - yp[j];	/* compute differentials */
		okaylist[i] = 1;
		j += 2;
	}
	if( tight126 ) {	/* eliminate 1 2 6 */
		okaylist[0] = 0;
		okaylist[1] = 0;
		okaylist[5] = 0;
	}
	if( tight135 ) {	/* eliminate 1 3 5 */
		okaylist[0] = 0;
		okaylist[2] = 0;
		okaylist[4] = 0;
	}
	for( j=0; j<6; j++ ) {
		if( okaylist[j] && diff[j] > 0.0 ) {	/* ignore tight spots */
			i = 2*j;
			n0 = i+1;
			if( rel_part[n0] < 0.0 ) {
				if( yp[i] != original_cell[i] ) {
					partialcount += 1.0;
					return;	/* bumped off */
				}
				else {
					yp[n0] = original_cell[i];	/* reduced dimension */
					gmacount += 1.0;
				}
			} else if( rel_part[i] > 0.0 ) {
				if( yp[n0] != original_cell[n0] ) {
					partialcount += 1.0;
					return;	/* bumped off */
				}
				else {
					yp[i] = original_cell[n0];	/* reduced dimension */
					gmacount += 1.0;
				}
			}
		}
	}
#endif
	if( info.depth < MAXDEPTH && cellcount < COUNTLIMIT) {
		/* If all else fails, subdivide.  */
		info.depth++;
		if( info.depth > maxdepth )
			maxdepth = info.depth;
		ROUND_NEAR;
		j = 0;
		for( i=0; i<6; i++ ) {	/* scaled differentials */
			diff[i] = SUBFRAC*(yp[j+1] - yp[j]);
			j += 2;
			if( diff[i] > 0.0 )
				bounds[i] = SUBDIV;
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
								verify_cell_original( t, info );
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
		printf("relation = [%g, %g]\n", 
			data[0], data[1]);
		printf("y = \n");
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", y[i], y[i+1]);
		printf("faceinfo = (%d %d)\n", faceinfo[0], faceinfo[1]);
		printf("crad1262 = [%0.18f\t%0.18f]\n", crad1262[0], crad1262[1]);
		printf("crad1352 = [%0.18f\t%0.18f]\n", crad1352[0], crad1352[1]);
#if PARTIALS
		printf("y partials = \n");
		for( i=0; i<12; i+=2 )
			printf("[%0.18f\t%0.18f]\n", rel_part[i], rel_part[i+1]);
		printf("yp = \n");
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", yp[i], yp[i+1]);
#endif
	}
}


void verify_cell( double y[12], struct Cell_Info info )
{
	int i, j, bounds[6], n0, n1, n2, n3, n4, n5;
	double diff[6], t[12], val;
	double x[12], yp[12], facex[6];
	double data[2], rel_part[12];
	double crad1262[2], crad1352[2];
	
	
	/* If part of the verification has already failed, bail. */
	if( failed )
		return;
	
	cellcount += 1.0;
	
	/* Compute square of edge lengths */
	
	ROUND_UP;
	for( i=1; i<12; i+=2 )
		x[i] = y[i]*y[i];
	ROUND_DOWN;
	for( i=0; i<12; i+=2 )
		x[i] = y[i]*y[i];
	if( x[1] > 8.0 )
		x[1] = 8.0;

	/* First try to verify over the current cell.  */

	/* Check face constraints. */
	
	/* face:  (1 2 6) -> (0 2 10) */
	facex[0] = x[0];
	facex[1] = x[1];
	facex[2] = x[2];
	facex[3] = x[3];
	facex[4] = x[10];
	facex[5] = x[11];
	s_crad3x2( facex, crad1262 );
	
	if( faceinfo[0] < 1 ) {	/* should be small */
		if( crad1262[0] > 2.0 ) {
			dualcount += 1.0;
			return;
		}
	} else {	/* should be large */
		if( crad1262[1] < 2.0 ) {
			sccutcount += 1.0;
			return;
		}
	}
	
	/* face:  (1 3 5) -> (0 4 8) */
	facex[2] = x[4];
	facex[3] = x[5];
	facex[4] = x[8];
	facex[5] = x[9];
	s_crad3x2( facex, crad1352 );

	if( faceinfo[1] < 1 ) {	/* should be small */
		if( crad1352[0] > 2.0 ) {
			dualcount += 1.0;
			return;
		}
	} else {	/* should be large */
		if( crad1352[1] < 2.0 ) {
			sccutcount += 1.0;
			return;
		}
	}
	
	/* Score cell properly */
	octa_score( y, x, data, rel_part );
		
	
	/* Attempt verification */	
	if( data[1] < slop ) {
		verifycount += 1.0;
		return;
	}
#if BADCHECK
	if( data[0] > slop ) {
		/* We are in big trouble. */
		failed = 1;
		printf("Cell failed:\n");
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", y[i], y[i+1]);
		printf("data = [%0.18g\t%0.18g]\n", data[0], data[1]);
		return;
	}
#endif
	
	/* if verification fails, try to bump off, or at least
		reduce dimension */

	for( i=0; i<12; i++ )
		yp[i] = y[i];	/* copy y */
	
#if PARTIALS
	
	/* Score cell properly (find min score) */
	
	/* t contains the relation partials */
	ROUND_NEAR;
	j = 0;
	for( i=0; i<6; i++ ) {
		diff[i] = yp[j+1] - yp[j];	/* compute differentials */
		j += 2;
	}
	for( j=0; j<6; j++ ) {
		if( diff[j] > 0.0 ) {	/* ignore tight spots */
			i = 2*j;
			n0 = i+1;
			if( rel_part[n0] < 0.0 ) {
				if( yp[i] != original_cell[i] ) {
					partialcount += 1.0;
					return;	/* bumped off */
				}
				else {
					yp[n0] = original_cell[i];	/* reduced dimension */
					gmacount += 1.0;
				}
			} else if( rel_part[i] > 0.0 ) {
				if( yp[n0] != original_cell[n0] ) {
					partialcount += 1.0;
					return;	/* bumped off */
				}
				else {
					yp[i] = original_cell[n0];	/* reduced dimension */
					gmacount += 1.0;
				}
			}
		}
	}
#endif
	if( info.depth < MAXDEPTH && cellcount < COUNTLIMIT) {
		/* If all else fails, subdivide.  */
		info.depth++;
		if( info.depth > maxdepth )
			maxdepth = info.depth;
		ROUND_NEAR;
		j = 0;
		for( i=0; i<6; i++ ) {	/* scaled differentials */
			diff[i] = SUBFRAC*(yp[j+1] - yp[j]);
			j += 2;
			if( diff[i] > 0.0 )
				bounds[i] = SUBDIV;
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
		printf("relation = [%g, %g]\n", 
			data[0], data[1]);
		printf("y = \n");
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", y[i], y[i+1]);
		printf("faceinfo = (%d %d)\n", faceinfo[0], faceinfo[1]);
		printf("crad1262 = [%0.18f\t%0.18f]\n", crad1262[0], crad1262[1]);
		printf("crad1352 = [%0.18f\t%0.18f]\n", crad1352[0], crad1352[1]);
#if PARTIALS
		printf("y partials = \n");
		for( i=0; i<12; i+=2 )
			printf("[%0.18f\t%0.18f]\n", rel_part[i], rel_part[i+1]);
		printf("yp = \n");
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", yp[i], yp[i+1]);
#endif
	}
}


void verify_bdry126( double y[12], struct Cell_Info info )
{
	int i, j, bounds[6], n0, n1, n2, n3, n4, n5;
	int okaylist[6];
	double diff[6], t[12], val;
	double x[12], yp[12], facex[6];
	double data[2], rel_part[12];
	double crad1352[2];
	
	
	/* If part of the verification has already failed, bail. */
	if( failed )
		return;
	
	cellcount += 1.0;
	
	/* Compute square of edge lengths */
	
	ROUND_UP;
	for( i=1; i<12; i+=2 )
		x[i] = y[i]*y[i];
	ROUND_DOWN;
	for( i=0; i<12; i+=2 )
		x[i] = y[i]*y[i];
	if( x[1] > 8.0 )
		x[1] = 8.0;

	/* First try to verify over the current cell.  */

	/* face:  (1 2 6) -> (0 2 10) */
	facex[0] = x[0];
	facex[1] = x[1];
	facex[2] = x[2];
	facex[3] = x[3];
	/* (1 2 6) is tight */
	
	/* Solve for x6 */
	face_edge( facex, facex + 2, facex + 4 );
	
	x[10] = facex[4];
	x[11] = facex[5];
	ROUND_DOWN;
	y[10] = sqrt( x[10] );
	ROUND_UP;
	y[11] = sqrt( x[11] );
	
	/* Check constraints on y6 */
	if( y[11] < original_cell[10] ) {
		vorcount += 1.0;
		return;
	}
	if( y[10] > original_cell[11] ) {
		vorcount += 1.0;
		return;
	}
	/* Adjust, if necessary */
	if( y[10] < original_cell[10] ) {
		y[10] = original_cell[10];
		ROUND_DOWN;
		x[10] = y[10]*y[10];
	}
	if( y[11] > original_cell[11] ) {
		y[11] = original_cell[11];
		ROUND_UP;
		x[11] = y[11]*y[11];
	}
	
	/* face:  (1 3 5) -> (0 4 8) */
	facex[2] = x[4];
	facex[3] = x[5];
	facex[4] = x[8];
	facex[5] = x[9];
	s_crad3x2( facex, crad1352 );

	if( faceinfo[1] < 1 ) {	/* should be small */
		if( crad1352[0] > 2.0 ) {
			dualcount += 1.0;
			return;
		}
	} else {	/* should be large */
		if( crad1352[1] < 2.0 ) {
			sccutcount += 1.0;
			return;
		}
	}
	
	/* Score cell properly */
	octa_score( y, x, data, rel_part );
		
	
	/* Attempt verification */	
	if( data[1] < slop ) {
		verifycount += 1.0;
		return;
	}
#if BADCHECK
	if( data[0] > slop && info.depth > 10 ) {
		/* We are in big trouble. */
		failed = 1;
		printf("Cell failed:\n");
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", y[i], y[i+1]);
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", x[i], x[i+1]);
		printf("data = [%0.18g\t%0.18g]\n", data[0], data[1]);
		return;
	}
#endif
	
	/* if verification fails, try to bump off, or at least
		reduce dimension */

	for( i=0; i<12; i++ )
		yp[i] = y[i];	/* copy y */
	
#if PARTIALS
	
	/* Score cell properly (find min score) */
	
	/* t contains the relation partials */
	ROUND_NEAR;
	j = 0;
	for( i=0; i<6; i++ ) {
		diff[i] = yp[j+1] - yp[j];	/* compute differentials */
		okaylist[i] = 1;
		j += 2;
	}
	/* eliminate 1 2 6 */
	okaylist[0] = 0;
	okaylist[1] = 0;
	okaylist[5] = 0;
	/* eliminate 1 3 5 */
	
	for( j=0; j<6; j++ ) {
		if( okaylist[j] && diff[j] > 0.0 ) {	/* ignore tight spots */
			i = 2*j;
			n0 = i+1;
			if( rel_part[n0] < 0.0 ) {
				if( yp[i] != original_cell[i] ) {
					partialcount += 1.0;
					return;	/* bumped off */
				}
				else {
					yp[n0] = original_cell[i];	/* reduced dimension */
					gmacount += 1.0;
				}
			} else if( rel_part[i] > 0.0 ) {
				if( yp[n0] != original_cell[n0] ) {
					partialcount += 1.0;
					return;	/* bumped off */
				}
				else {
					yp[i] = original_cell[n0];	/* reduced dimension */
					gmacount += 1.0;
				}
			}
		}
	}
#endif
	if( info.depth < MAXDEPTH && cellcount < COUNTLIMIT) {
		/* If all else fails, subdivide.  */
		info.depth++;
		if( info.depth > maxdepth )
			maxdepth = info.depth;
		ROUND_NEAR;
		j = 0;
		for( i=0; i<6; i++ ) {	/* scaled differentials */
			diff[i] = SUBFRAC*(yp[j+1] - yp[j]);
			j += 2;
			if( diff[i] > 0.0 )
				bounds[i] = SUBDIV;
			else
				bounds[i] = 1;
		}
		bounds[5] = 1;
		diff[5] = yp[11] - yp[10];
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
								verify_bdry126( t, info );
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
		printf("relation = [%g, %g]\n", 
			data[0], data[1]);
		printf("y = \n");
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", y[i], y[i+1]);
		printf("faceinfo = (%d %d)\n", faceinfo[0], faceinfo[1]);
		/* printf("crad1262 = [%0.18f\t%0.18f]\n", crad1262[0], crad1262[1]); */
		printf("crad1352 = [%0.18f\t%0.18f]\n", crad1352[0], crad1352[1]);
#if PARTIALS
		printf("y partials = \n");
		for( i=0; i<12; i+=2 )
			printf("[%0.18f\t%0.18f]\n", rel_part[i], rel_part[i+1]);
		printf("yp = \n");
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", yp[i], yp[i+1]);
#endif
	}
}


void verify_bdry135( double y[12], struct Cell_Info info )
{
	int i, j, bounds[6], n0, n1, n2, n3, n4, n5;
	int okaylist[6];
	double diff[6], t[12], val;
	double x[12], yp[12], facex[6];
	double data[2], rel_part[12];
	double crad1262[2];
	
	
	/* If part of the verification has already failed, bail. */
	if( failed )
		return;
	
	cellcount += 1.0;
	
	if( ISNAN( y[10] ) || ISNAN( y[11] ) ) {
		failed = 1;
		printf("Depth = %d\n", info.depth);
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", y[i], y[i+1]);
		return;
	}
	
	/* Compute square of edge lengths */
	
	ROUND_UP;
	for( i=1; i<12; i+=2 )
		x[i] = y[i]*y[i];
	ROUND_DOWN;
	for( i=0; i<12; i+=2 )
		x[i] = y[i]*y[i];
	if( x[1] > 8.0 )
		x[1] = 8.0;

	/* First try to verify over the current cell.  */

	/* face:  (1 2 6) -> (0 2 10) */
	facex[0] = x[0];
	facex[1] = x[1];
	facex[2] = x[2];
	facex[3] = x[3];
	facex[4] = x[10];
	facex[5] = x[11];
	s_crad3x2( facex, crad1262 );
	
	if( faceinfo[0] < 1 ) {	/* should be small */
		if( crad1262[0] > 2.0 ) {
			dualcount += 1.0;
			return;
		}
	} else {	/* should be large */
		if( crad1262[1] < 2.0 ) {
			sccutcount += 1.0;
			return;
		}
	}
	
	/* face:  (1 3 5) -> (0 4 8) */
	facex[2] = x[4];
	facex[3] = x[5];
	/* (1 3 5) is tight */
	
	/* Solve for x5 */
	face_edge( facex, facex + 2, facex + 4 );
	
	x[8] = facex[4];
	x[9] = facex[5];
	ROUND_DOWN;
	y[8] = sqrt( x[8] );
	ROUND_UP;
	y[9] = sqrt( x[9] );
	
	/* Check constraints on y5 */
	if( y[9] < original_cell[8] ) {
		vorcount += 1.0;
		return;
	}
	if( y[8] > original_cell[9] ) {
		vorcount += 1.0;
		return;
	}
	/* Adjust, if necessary */
	if( y[8] < original_cell[8] ) {
		y[8] = original_cell[8];
		ROUND_DOWN;
		x[8] = y[8]*y[8];
	}
	if( y[9] > original_cell[9] ) {
		y[9] = original_cell[9];
		ROUND_UP;
		x[9] = y[9]*y[9];
	}
	
	/* Score cell properly */
	octa_score( y, x, data, rel_part );
		
	
	/* Attempt verification */	
	if( data[1] < slop ) {
		verifycount += 1.0;
		return;
	}
#if BADCHECK
	if( data[0] > slop && info.depth > 10 ) {
		/* We are in big trouble. */
		failed = 1;
		printf("Cell failed:\n");
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", y[i], y[i+1]);
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", x[i], x[i+1]);
		printf("data = [%0.18g\t%0.18g]\n", data[0], data[1]);
		return;
	}
#endif

	/* if verification fails, try to bump off, or at least
		reduce dimension */

	for( i=0; i<12; i++ )
		yp[i] = y[i];	/* copy y */
	
#if PARTIALS
	
	/* Score cell properly (find min score) */
	
	/* t contains the relation partials */
	ROUND_NEAR;
	j = 0;
	for( i=0; i<6; i++ ) {
		diff[i] = yp[j+1] - yp[j];	/* compute differentials */
		okaylist[i] = 1;
		j += 2;
	}
	/* eliminate 1 2 6 */
	/* eliminate 1 3 5 */
	okaylist[0] = 0;
	okaylist[2] = 0;
	okaylist[4] = 0;
	
	for( j=0; j<6; j++ ) {
		if( okaylist[j] && diff[j] > 0.0 ) {	/* ignore tight spots */
			i = 2*j;
			n0 = i+1;
			if( rel_part[n0] < 0.0 ) {
				if( yp[i] != original_cell[i] ) {
					partialcount += 1.0;
					return;	/* bumped off */
				}
				else {
					yp[n0] = original_cell[i];	/* reduced dimension */
					gmacount += 1.0;
				}
			} else if( rel_part[i] > 0.0 ) {
				if( yp[n0] != original_cell[n0] ) {
					partialcount += 1.0;
					return;	/* bumped off */
				}
				else {
					yp[i] = original_cell[n0];	/* reduced dimension */
					gmacount += 1.0;
				}
			}
		}
	}
#endif
	if( info.depth < MAXDEPTH && cellcount < COUNTLIMIT) {
		/* If all else fails, subdivide.  */
		info.depth++;
		if( info.depth > maxdepth )
			maxdepth = info.depth;
		ROUND_NEAR;
		j = 0;
		for( i=0; i<6; i++ ) {	/* scaled differentials */
			diff[i] = SUBFRAC*(yp[j+1] - yp[j]);
			j += 2;
			if( diff[i] > 0.0 )
				bounds[i] = SUBDIV;
			else
				bounds[i] = 1;
		}
		bounds[4] = 1;
		diff[4] = yp[9] - yp[8];
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
								verify_bdry135( t, info );
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
		printf("relation = [%g, %g]\n", 
			data[0], data[1]);
		printf("y = \n");
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", y[i], y[i+1]);
		printf("faceinfo = (%d %d)\n", faceinfo[0], faceinfo[1]);
		printf("crad1262 = [%0.18f\t%0.18f]\n", crad1262[0], crad1262[1]);
		/* printf("crad1352 = [%0.18f\t%0.18f]\n", crad1352[0], crad1352[1]); */
#if PARTIALS
		printf("y partials = \n");
		for( i=0; i<12; i+=2 )
			printf("[%0.18f\t%0.18f]\n", rel_part[i], rel_part[i+1]);
		printf("yp = \n");
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", yp[i], yp[i+1]);
#endif
	}
}


void verify_bdry( double y[12], struct Cell_Info info )
{
	int i, j, bounds[6], n0, n1, n2, n3, n4, n5;
	double diff[6], t[12], val;
	double x[12], yp[12], facex[6];
	double data[2], rel_part[12];
	
	
	/* If part of the verification has already failed, bail. */
	if( failed )
		return;
	
	cellcount += 1.0;
	
	/* Compute square of edge lengths */
	
	ROUND_UP;
	for( i=1; i<12; i+=2 )
		x[i] = y[i]*y[i];
	ROUND_DOWN;
	for( i=0; i<12; i+=2 )
		x[i] = y[i]*y[i];
	if( x[1] > 8.0 )
		x[1] = 8.0;

	/* First try to verify over the current cell.  */

	/* face:  (1 2 6) -> (0 2 10) */
	facex[0] = x[0];
	facex[1] = x[1];
	facex[2] = x[2];
	facex[3] = x[3];
	/* (1 2 6) is tight */
	
	/* Solve for x6 */
	face_edge( facex, facex + 2, facex + 4 );
	
	x[10] = facex[4];
	x[11] = facex[5];
	ROUND_DOWN;
	y[10] = sqrt( x[10] );
	ROUND_UP;
	y[11] = sqrt( x[11] );
	
	/* Check constraints on y6 */
	if( y[11] < original_cell[10] ) {
		vorcount += 1.0;
		return;
	}
	if( y[10] > original_cell[11] ) {
		vorcount += 1.0;
		return;
	}
	/* Adjust, if necessary */
	if( y[10] < original_cell[10] ) {
		y[10] = original_cell[10];
		ROUND_DOWN;
		x[10] = y[10]*y[10];
	}
	if( y[11] > original_cell[11] ) {
		y[11] = original_cell[11];
		ROUND_UP;
		x[11] = y[11]*y[11];
	}
	
	/* face:  (1 3 5) -> (0 4 8) */
	facex[2] = x[4];
	facex[3] = x[5];
	/* (1 3 5) is tight */
	
	/* Solve for x5 */
	face_edge( facex, facex + 2, facex + 4 );
	
	x[8] = facex[4];
	x[9] = facex[5];
	ROUND_DOWN;
	y[8] = sqrt( x[8] );
	ROUND_UP;
	y[9] = sqrt( x[9] );
	
	/* Check constraints on y5 */
	if( y[9] < original_cell[8] ) {
		vorcount += 1.0;
		return;
	}
	if( y[8] > original_cell[9] ) {
		vorcount += 1.0;
		return;
	}
	/* Adjust, if necessary */
	if( y[8] < original_cell[8] ) {
		y[8] = original_cell[8];
		ROUND_DOWN;
		x[8] = y[8]*y[8];
	}
	if( y[9] > original_cell[9] ) {
		y[9] = original_cell[9];
		ROUND_UP;
		x[9] = y[9]*y[9];
	}
	
	/* Score cell properly */
	octa_score( y, x, data, rel_part );
		
	
	/* Attempt verification */	
	if( data[1] < slop ) {
		verifycount += 1.0;
		return;
	}
#if BADCHECK
	if( data[0] > slop && info.depth > 10 ) {
		/* We are in big trouble. */
		failed = 1;
		printf("Cell failed:\n");
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", y[i], y[i+1]);
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", x[i], x[i+1]);
		printf("data = [%0.18g\t%0.18g]\n", data[0], data[1]);
		return;
	}
#endif
	
	/* if verification fails, try to bump off, or at least
		reduce dimension */

	for( i=0; i<12; i++ )
		yp[i] = y[i];	/* copy y */
	
#if PARTIALS
	
	/* Score cell properly (find min score) */
	
	/* t contains the relation partials */
	ROUND_NEAR;
	i = 3;
	j = 6;
	diff[i] = yp[j+1] - yp[j];	/* compute differentials */

	/* eliminate 1 2 6 */
	/* eliminate 1 3 5 */
	
	/* only push at y4 */
	for( j=3; j<4; j++ ) {
		if( diff[j] > 0.0 ) {	/* ignore tight spots */
			i = 2*j;
			n0 = i+1;
			if( rel_part[n0] < 0.0 ) {
				if( yp[i] != original_cell[i] ) {
					partialcount += 1.0;
					return;	/* bumped off */
				}
				else {
					yp[n0] = original_cell[i];	/* reduced dimension */
					gmacount += 1.0;
				}
			} else if( rel_part[i] > 0.0 ) {
				if( yp[n0] != original_cell[n0] ) {
					partialcount += 1.0;
					return;	/* bumped off */
				}
				else {
					yp[i] = original_cell[n0];	/* reduced dimension */
					gmacount += 1.0;
				}
			}
		}
	}
#endif
	if( info.depth < MAXDEPTH && cellcount < COUNTLIMIT) {
		/* If all else fails, subdivide.  */
		info.depth++;
		if( info.depth > maxdepth )
			maxdepth = info.depth;
		ROUND_NEAR;
		j = 0;
		for( i=0; i<4; i++ ) {	/* scaled differentials */
			diff[i] = SUBFRAC*(yp[j+1] - yp[j]);
			j += 2;
			if( diff[i] > 0.0 )
				bounds[i] = SUBDIV;
			else
				bounds[i] = 1;
		}
		bounds[4] = 1;
		diff[4] = yp[9] - yp[8];
		bounds[5] = 1;
		diff[5] = yp[11] - yp[10];
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
								verify_bdry( t, info );
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
		printf("relation = [%g, %g]\n", 
			data[0], data[1]);
		printf("y = \n");
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", y[i], y[i+1]);
		printf("faceinfo = (%d %d)\n", faceinfo[0], faceinfo[1]);
/*
		printf("crad1262 = [%0.18f\t%0.18f]\n", crad1262[0], crad1262[1]);
		printf("crad1352 = [%0.18f\t%0.18f]\n", crad1352[0], crad1352[1]);
*/
#if PARTIALS
		printf("y partials = \n");
		for( i=0; i<12; i+=2 )
			printf("[%0.18f\t%0.18f]\n", rel_part[i], rel_part[i+1]);
		printf("yp = \n");
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", yp[i], yp[i+1]);
#endif
	}
}


void set_relations( int mynumber )
{
	int i, n;

	/* Order:  y1 y2 y3 1 y5 y6 dih1 dih2 dih3 */
	double rel_coeffs[18][9] = 
		{{-0.171086999999999989,0.282483000000000004,
	    0.0423758999999999996,-0.25889409228076996,
	    -0.0269357000000000024,-0.0813089000000000083,
	    0.0761185999999999918,0.198866999999999993,
	    0},
	    
	    {654.232435199999962,322.540000000000048,
	    0.449699000000000026,1998.90161565376978,
	    0.444429000000000051,-2000.,
	    0.35668299999999995,0,
	    0},
	    
	    {-2.78326999999999991,0.225892000000000026,
	    -0.147886999999999986,7.25363567174538204,
	    -0.206484000000000023,0.187246000000000005,
	    0.0767790999999999979,0.198866999999999993,
	    0.198866999999999993},
	    
	    {-2.96159999999999978,0.0112946000000000013,
	    0.311179000000000005,7.12463517774935528,
	    0.181643000000000007,-0.112644999999999995,
	    0.171553000000000022,0.198866999999999993,
	    0},
	    
	    {-623.135749999999966,-29.7859000000000015,
	    0.37367600000000003,1999.34328346878302,
	    0.39482800000000001,-89.7991000000000028,
	    0.296738999999999997,0,
	    0.198866999999999993},
	    
	    {2000.,304.735000000000022,
	    -0.171057999999999976,-2000.65671653121679,
	    -0.186587000000000013,-2000.,
	    0.296738999999999997,0,
	    0},
	    
	    {2000.,304.72699999999997,
	    -0.18345900000000002,-2000.65671653121679,
	    -0.167931000000000008,-2000.,
	    0.296738999999999997,0,
	    0},
	    
	    {527.638000000000051,122.604000000000001,
	    -0.228975000000000017,-2000.86671653121683,
	    -0.163413000000000003,122.596000000000016,
	    0.296738999999999997,0,
	    0},
	    
	    {-6.13844471300000016,0.184409999999999989,
	    -0.165877000000000007,16.9690956717453822,
	    -0.211226000000000002,0.138639000000000001,
	    0.0767790999999999979,0.198866999999999993,
	    0.198866999999999993},
	    
	    {-6.92706975000000024,0.190011999999999972,
	    -0.193186999999999997,19.303663429310891,
	    -0.25031500000000002,0.138883,
	    0.083739100000000004,0.198866999999999993,
	    0.198866999999999993},
	    
	    {-2.80038000000000018,0.226456999999999997,
	    -0.147743999999999982,7.30017317197146464,
	    -0.208187999999999995,0.187876999999999974,
	    0.0784740999999999999,0.198866999999999993,
	    0.198866999999999993},
	    
	    {4.23632000000000008,0.0999798999999999971,
	    0.700260000000000015,-15.5730468787901425,
	    0.749277000000000015,-0.00331063000000000062,
	    0.178919999999999994,0,
	    0.198866999999999993},
	    
	    {-8.19969999999999998,-0.165254000000000011,
	    0.261357000000000017,22.4461693810909723,
	    0.168102999999999997,-0.120040999999999997,
	    0.170066999999999985,0.198866999999999993,
	    0},
	    
	    {9.79630022499999952,0.192513999999999985,
	    0.471269000000000026,-30.7530748222506433,
	    0.37612899999999998,0.248019999999999996,
	    0.171553000000000022,0.198866999999999993,
	    0},
	    
	    {4.75617,0.741197000000000016,
	    -0.118935999999999996,-16.2236668280285334,
	    -0.0420990999999999981,0.736018999999999934,
	    0.0784740999999999999,0,
	    0},
	    
	    {-0.661767000000000038,-0.121748999999999996,
	    0.534628999999999976,-0.095233773206075849,
	    0.427205000000000012,-0.0958335000000000114,
	    0.174341999999999996,0.198866999999999993,
	    0},
	    
	    {-0.372773999999999983,0.0475567000000000028,
	    0.332564999999999999,0.529846182719230007,
	    -0.186431000000000004,0.0085941200000000002,
	    0.0761185999999999918,0,
	    0},
	    
	    {-1.71994000000000006,0.0332921999999999984,
	    0.0846380000000000087,4.40768817197146489,
	    0.0501668000000000002,-0.191493999999999999,
	    0.0784740999999999999,0.198866999999999993,
	    0.198866999999999993}};
    
/* Order:  octavor, gma, 126 135 */
  int rel_entries[18][4] = 
    {{1,0,0,1},{0,1,0,0},{1,0,0,1},{1,0,1,0},
    {1,0,1,0},{1,0,0,1},{1,0,0,1},{1,0,0,1},
    {1,0,0,1},{1,0,0,1},{1,0,0,1},{1,0,1,0},
    {1,0,1,0},{1,0,1,0},{1,0,0,1},{1,0,1,0},
    {1,0,1,0},{1,0,1,1}};
	
	n = mynumber - 1;
	
	/*
	printf("mynumber in set_relations = %d\n", mynumber);
	printf("n in set_relations = %d\n", n);
	*/
	
/* Order of relation constants:
	sol, gma, vor, octavor, dih1, dih2, dih3	
	 0		1		 2			3		 		4		  5		  6	*/

/* 	rel describes which values we want computed
independently.
		relconst gives the relation constants.  */

		for( i=0; i<7; i++ ) {
			taylor_rel[i] = 0;
			taylor_relconst[i] = 0.0;
		}
		
		for( i=0; i<6; i++ ) {
			y_coefficients[i] = rel_coeffs[n][i];
		}
		
		for( i=4; i<7; i++ ) {
			taylor_relconst[i] = rel_coeffs[n][i+2];
		}
	
		taylor_relconst[3] = rel_entries[n][0];
		taylor_relconst[1] = rel_entries[n][1];
		faceinfo[0] = rel_entries[n][2];
		faceinfo[1] = rel_entries[n][3];
		
		n = mynumtotom( mynumber );
		set_domain( n );
		/*
		ROUND_NEAR;
		printf("rel_coeffs = \n");
		for( i=0; i<18; i++ ) {
			for( n=0; n<9; n++ ) {
				if( n%3 == 0 )
					printf("\n");
				printf("%.18f\t", rel_coeffs[i][n]);
			}
			printf("\n");
		}
		printf("rel_entries = \n");
		for( i=0; i<18; i++ ) {
			for( n=0; n<4; n++ ) {
				printf("%d\t", rel_entries[i][n]);
			}
			printf("\n");
		}
		*/
		/*
		ROUND_NEAR;
		printf("taylor_relconst = \n");
		for( i=0; i<7; i++ )
			printf("%.18f\n", taylor_relconst[i]);
		printf("y coefficients = \n");
		for( i=0; i<6; i++ )
			printf("%.18f\n", y_coefficients[i]);
		*/
}


int tomtomynum( int num )
{
	int i;
	
	for( i=1; i<20; i++ ) {
		if( mynumtotom( i ) == num ) {
			return( i );
		}
	}
	return( -1 );
}


int mynumtotom( int mynum )
{
	int list[18] = 
		{355, 617, 688, 690, 703, 705,
		 707, 709, 761, 765, 837, 859,
		 882, 899, 969, 972, 1000, 1008};
	
	if( mynum > 0 && mynum < 19 )
		return( list[ mynum - 1 ] );
	else {
		/* printf("mynmtotom:  domain error\n"); */
		return( -1 );
	}
}


void set_domain( int casenum )
{
	int i;
	double xmin[6], xmax[6], tet[12];

	switch( casenum ) {
		case 355:
			xmin[0]=2.51;
			xmin[1]=2;
			xmin[2]=2;
			xmin[3]=2;
			xmin[4]=2;
			xmin[5]=2;

			xmax[0]=2.828427125;
			xmax[1]=2.255;
			xmax[2]=2.51;
			xmax[3]=2.255;
			xmax[4]=2.51;
			xmax[5]=2.255;
			break;
		case 617:
			xmin[0]=2.669213562;
			xmin[1]=2.255;
			xmin[2]=2;
			xmin[3]=2.255;
			xmin[4]=2;
			xmin[5]=2.255;

			xmax[0]=2.828427125;
			xmax[1]=2.3825;
			xmax[2]=2.1275;
			xmax[3]=2.3825;
			xmax[4]=2.1275;
			xmax[5]=2.3825;
			break;
		case 688:
			xmin[0]=2.669213562;
			xmin[1]=2;
			xmin[2]=2;
			xmin[3]=2;
			xmin[4]=2;
			xmin[5]=2;

			xmax[0]=2.828427125;
			xmax[1]=2.1275;
			xmax[2]=2.1275;
			xmax[3]=2.1275;
			xmax[4]=2.1275;
			xmax[5]=2.1275;
			break;
		case 690:
			xmin[0]=2.669213562;
			xmin[1]=2;
			xmin[2]=2;
			xmin[3]=2;
			xmin[4]=2;
			xmin[5]=2;

			xmax[0]=2.828427125;
			xmax[1]=2.1275;
			xmax[2]=2.1275;
			xmax[3]=2.1275;
			xmax[4]=2.1275;
			xmax[5]=2.1275;
			break;
		case 703:
			xmin[0]=2.669213562;
			xmin[1]=2;
			xmin[2]=2;
			xmin[3]=2;
			xmin[4]=2;
			xmin[5]=2;

			xmax[0]=2.828427125;
			xmax[1]=2.1275;
			xmax[2]=2.1275;
			xmax[3]=2.255;
			xmax[4]=2.1275;
			xmax[5]=2.1275;
			break;
		case 705:
			xmin[0]=2.669213562;
			xmin[1]=2.1275;
			xmin[2]=2;
			xmin[3]=2.1275;
			xmin[4]=2.1275;
			xmin[5]=2.1275;

			xmax[0]=2.828427125;
			xmax[1]=2.255;
			xmax[2]=2.1275;
			xmax[3]=2.255;
			xmax[4]=2.255;
			xmax[5]=2.255;
			break;
		case 707:
			xmin[0]=2.669213562;
			xmin[1]=2.1275;
			xmin[2]=2.1275;
			xmin[3]=2.1275;
			xmin[4]=2;
			xmin[5]=2.1275;

			xmax[0]=2.828427125;
			xmax[1]=2.255;
			xmax[2]=2.255;
			xmax[3]=2.255;
			xmax[4]=2.1275;
			xmax[5]=2.255;
			break;
		case 709:
			xmin[0]=2.669213562;
			xmin[1]=2.1275;
			xmin[2]=2.1275;
			xmin[3]=2.1275;
			xmin[4]=2.1275;
			xmin[5]=2.1275;

			xmax[0]=2.828427125;
			xmax[1]=2.255;
			xmax[2]=2.255;
			xmax[3]=2.255;
			xmax[4]=2.255;
			xmax[5]=2.255;
			break;
		case 761:
			xmin[0]=2.748820344;
			xmin[1]=2;
			xmin[2]=2;
			xmin[3]=2;
			xmin[4]=2;
			xmin[5]=2;

			xmax[0]=2.828427125;
			xmax[1]=2.06375;
			xmax[2]=2.06375;
			xmax[3]=2.06375;
			xmax[4]=2.06375;
			xmax[5]=2.06375;
			break;
		case 765:
			xmin[0]=2.748820344;
			xmin[1]=2;
			xmin[2]=2;
			xmin[3]=2;
			xmin[4]=2;
			xmin[5]=2;

			xmax[0]=2.828427125;
			xmax[1]=2.06375;
			xmax[2]=2.06375;
			xmax[3]=2.06375;
			xmax[4]=2.06375;
			xmax[5]=2.06375;
			break;
		case 837:
			xmin[0]=2.669213562;
			xmin[1]=2;
			xmin[2]=2;
			xmin[3]=2;
			xmin[4]=2;
			xmin[5]=2;

			xmax[0]=2.828427125;
			xmax[1]=2.1275;
			xmax[2]=2.1275;
			xmax[3]=2.1275;
			xmax[4]=2.1275;
			xmax[5]=2.1275;
			break;
		case 859:
			xmin[0]=2.51;
			xmin[1]=2;
			xmin[2]=2;
			xmin[3]=2;
			xmin[4]=2;
			xmin[5]=2;

			xmax[0]=2.828427125;
			xmax[1]=2.51;
			xmax[2]=2.255;
			xmax[3]=2.51;
			xmax[4]=2.255;
			xmax[5]=2.51;
			break;
		case 882:
			xmin[0]=2.748820344;
			xmin[1]=2;
			xmin[2]=2;
			xmin[3]=2;
			xmin[4]=2;
			xmin[5]=2;

			xmax[0]=2.828427125;
			xmax[1]=2.06375;
			xmax[2]=2.06375;
			xmax[3]=2.06375;
			xmax[4]=2.06375;
			xmax[5]=2.06375;
			break;
		case 899:
			xmin[0]=2.748820344;
			xmin[1]=2.06375;
			xmin[2]=2;
			xmin[3]=2;
			xmin[4]=2;
			xmin[5]=2;

			xmax[0]=2.828427125;
			xmax[1]=2.1275;
			xmax[2]=2.1275;
			xmax[3]=2.1275;
			xmax[4]=2.1275;
			xmax[5]=2.06375;
			break;
		case 969:
			xmin[0]=2.669213562;
			xmin[1]=2;
			xmin[2]=2;
			xmin[3]=2;
			xmin[4]=2;
			xmin[5]=2;

			xmax[0]=2.828427125;
			xmax[1]=2.255;
			xmax[2]=2.255;
			xmax[3]=2.255;
			xmax[4]=2.255;
			xmax[5]=2.255;
			break;
		case 972:
			xmin[0]=2.51;
			xmin[1]=2;
			xmin[2]=2;
			xmin[3]=2;
			xmin[4]=2;
			xmin[5]=2;

			xmax[0]=2.828427125;
			xmax[1]=2.255;
			xmax[2]=2.255;
			xmax[3]=2.255;
			xmax[4]=2.255;
			xmax[5]=2.255;
			break;
		case 1000:
			xmin[0]=2.51;
			xmin[1]=2;
			xmin[2]=2;
			xmin[3]=2;
			xmin[4]=2;
			xmin[5]=2;

			xmax[0]=2.828427125;
			xmax[1]=2.51;
			xmax[2]=2.51;
			xmax[3]=2.51;
			xmax[4]=2.51;
			xmax[5]=2.51;
			break;
		case 1008:
			xmin[0]=2.51;
			xmin[1]=2;
			xmin[2]=2;
			xmin[3]=2;
			xmin[4]=2;
			xmin[5]=2;

			xmax[0]=2.828427125;
			xmax[1]=2.255;
			xmax[2]=2.255;
			xmax[3]=2.255;
			xmax[4]=2.255;
			xmax[5]=2.255;
			break;
		default:
			printf("Unrecognized case: %d\n", casenum );
			return;
			break;
	}
	for( i=0; i<6; i++ ) {
		tet[2*i  ] = xmin[i];
		tet[2*i+1] = xmax[i];
	}
	/* Now correct the domain. */
	fix_domain( tet, original_cell );
}


void fix_domain( double tet[12], double fixtet[12] )
{
	int i, j, k, n, twotothen;
	double width[2], start;
	
	for( i=0; i<6; i++ ) {
		j = 2*i;
		ROUND_NEAR;
		if( i == 0 ) {
			width[0] = TWOSQRT2_LO - TWO51_HI;
			width[1] = TWOSQRT2_HI - TWO51_LO;
			start = TWO51_LO;
		} else {
			width[0] = TWO51_LO - 2.0;
			width[1] = TWO51_HI - 2.0;
			start = 2.0;
		}
		n = log( width[1]/(tet[j+1] - tet[j]) )/log(2.0) + 0.5;
		twotothen = (int) pow( 2.0, (double) n ) + 0.5;
		k = twotothen*(tet[j] - start)/width[1] + 0.5;
		ROUND_DOWN;
		fixtet[j] = start + ((double) k)*width[0]/twotothen;
		ROUND_UP;
		fixtet[j+1] = start + ((double) (k+1))*width[1]/twotothen;
	}
	/*
	ROUND_NEAR;
	printf("\ntet = \n");
	for( i=0; i<12; i+=2 )
		printf("%0.18f\t%0.18f\n", tet[i], tet[i+1]);
	printf("fixtet = \n");
	for( i=0; i<12; i+=2 )
		printf("%0.18f\t%0.18f\n", fixtet[i], fixtet[i+1]);
	printf("\n");
	*/
}


void octa_score( double y[12], double x[12],
	double out[2], double out_partials[12] )
{
	int i, j, jp;
	double t[12], compvals[14], rel_val[2];
	double yc, sum;
	
	t_composite( taylor_rel, taylor_relconst, y, x, 
		rel_val, t, compvals );
	
	ROUND_DOWN;
	sum = 0.0;
	for( i=0; i<6; i++ ) {
		j = 2*i;
		jp = j + 1;
		yc = y_coefficients[i];
		if( i != 3 ) {
			if( yc > 0.0 ) {
				sum += yc*y[j];
			} else {
				sum += yc*y[jp];
			}
		} else {
			sum += yc;
		}
	}
	out[0] = sum + rel_val[0];
	
	/*
	ROUND_NEAR;
	printf("ycontrib[0] = %.18f\n", sum);
	*/
	
	ROUND_UP;
	sum = 0.0;
	for( i=0; i<6; i++ ) {
		j = 2*i;
		jp = j + 1;
		yc = y_coefficients[i];
		if( i != 3 ) {
			if( yc > 0.0 ) {
				sum += yc*y[jp];
			} else {
				sum += yc*y[j];
			}
		} else {
			sum += yc;
		}
	}
	out[1] = sum + rel_val[1];
	
	/*
	ROUND_NEAR;
	printf("ycontrib[1] = %.18f\n", sum);
	printf("rel_val = [%.18f, %.18f]\n", rel_val[0], rel_val[1]);
	*/
	
	for( i=0; i<12; i+=2 ) {
		i_mult( y + i, t + i, compvals + i );
	}
	ROUND_DOWN;
	for( i=0; i<6; i++ ) {
		j = 2*i;
		if( i != 3 ) {
			out_partials[j] = 2.0*compvals[j] + y_coefficients[i];
		} else {
			out_partials[j] = 2.0*compvals[j];
		}
	}
	ROUND_UP;
	for( i=0; i<6; i++ ) {
		j = 2*i + 1;
		if( i != 3 ) {
			out_partials[j] = 2.0*compvals[j] + y_coefficients[i];
		} else {
			out_partials[j] = 2.0*compvals[j];
		}
	}
	
	/*
	ROUND_NEAR;
	printf("x partials (taylor)= \n");
	for( i=0; i<12; i+=2 )
		printf("[%0.18f\t%0.18f]\n", t[i], t[i+1]);
	printf("y partials (taylor)= \n");
	for( i=0; i<12; i+=2 )
		printf("[%0.18f\t%0.18f]\n", 2.0*compvals[i], 2.0*compvals[i+1]);
	*/
}


void testvalues( void )
{
	double testvals[18][6] = 
		{{2.82842461936835798,2.00266023434319074,2,2,2.00266023434303397,2},
		{2.70718338079173648,2.28729473854791321,2.12749999999999995,2.38249999999999984,
		2.12749999999999995,2.25499999999999989},
		{2.80384150402050469,2.12749999999999995,2.12749999999999995,2,2.12749999999999995,
		2.12749999999999995},
		{2.80384150402050514,2.12749999999999995,2.12749999999999995,2,2.12749999999999995,
		2.12749999999999995},
		{2.80384150402050603,2.12749999999999995,2.12749999999999995,2.25499999999999989,2,
		2.12749999999999995},
		{2.80384150402050514,2.12749999999999995,2,2.25499999999999989,2.2457434658651847,
		2.12749999999999995},
		{2.80384150402050514,2.12749999999999995,2.24574346586518425,2.25499999999999989,2,
		2.12749999999999995},
		{2.79631627722876708,2.14502135906075786,2.12749999999999995,2.25499999999999989,
		2.16230864522645394,2.14497844858591025},
		{2.82248877867752235,2,2.0637500000000002,2,2.0637500000000002,2.0637500000000002},
		{2.82248877867751613,2,2.0637500000000002,2,2.0637500000000002,2},
		{2.82842712474618985,2.00000000265083733,2.00000000044859316,2,2,2.00000002685834177},
		{2.82484517519379263,2.09808958271173918,2.0980895827117374,2,2,2},
		{2.82842712474618985,2,2.00000002921293918,2.01590414065797141,2.00000001363876256,2},
		{2.82805565102432821,2.0637500000000002,2.03215055649402077,2,2,2},
		{2.82687329199290893,2.04326649009163974,2,2,2.06518623120335132,2.02242106962795276},
		{2.82842711544692849,2.0001624324625471,2.00016208288528485,2.06035467348270096,2,2},
		{2.82842712474618985,2.00000000000000044,2.00000003798905635,2.02067129011180269,2,2},
		{2.82842712474619118,2,2,2,2,2}};
	double outvals[18] =
		{2.37994371599353609e-05,
		0.0559879938979149783,
		0.000540382517622034486,
		0.00963559576548286834,
		0.095305472624010304,
		0.0855911380603315841,
		0.0893324524985614765,
		0.0187314229706921492,
		0.000589604697884164169,
		0.00676353214323453453,
		0.000464848815941484373,
		0.0100683169525948996,
		0.000620748953891574196,
		0.0072860905003857418,
		0.00723663381365875855,
		0.014480362231012317,
		0.000367756480027485542,
		0.000602049142830632855};
	int casenumber, mynumber, i, n, found, start, stop;
	double y[12], x[12], val[2], partials[12];
	double eps, facex[6], crad1262[2], crad1352[2];
	double eps2, temp;
	double lambdas[2], origslop;
	
	eps2 = 1.0e-10;
	printf("Welcome to Octahedra.  \n");
	printf("Enter eps:  ");
	scanf("%lf", &eps);
	printf("Enter slop:  ");
	scanf("%lf", &slop);
	printf("\t\tSlop = %g\n", slop );
	
	origslop = slop;
	
	/* stop = 2; */
	stop = 19;
	for( start = 1; start < stop; start++ ) {
		/*
		mynumber = -1;
		while( mynumber == -1 ) {
			printf("Enter case number:  ");
			scanf("%d", &casenumber );
			printf("\nCase = %d\n", casenumber );
			mynumber = tomtomynum( casenumber );
			if( mynumber == -1 )
				printf("Invalid case.  Try again.\n");
		}
		*/
		casenumber = mynumtotom( start );
		mynumber = start;
		printf("\nCase = %d\n", casenumber );
		
		i_init();	/* initialize interval constants */
		
		ROUND_NEAR;
		
		n = mynumber - 1;
		maxdepth = n;
		slop = origslop + 0.5*outvals[n];
		printf("slop = %g\n", slop );
		
		/* printf("(mynumber, n) = (%d, %d)\n", mynumber, n); */

		set_relations( mynumber );
		/* set_domain( casenumber ); */
		for( i=0; i<6; i++ ) {
			temp = testvals[n][i] - eps/2;
			if( temp < 2.0 )
				temp = 2.0;
			y[2*i  ] = temp;
			y[2*i+1] = testvals[n][i] + eps/2;
		}
		ROUND_DOWN;
		for( i=0; i<12; i+=2 ) {
			x[i] = y[i]*y[i];
		}
		ROUND_UP;
		for( i=1; i<12; i+=2 ) {
			x[i] = y[i]*y[i];
		}
		octa_score( y, x, val, partials );
		
		/* face:  (1 2 6) -> (0 2 10) */
		facex[0] = x[0];
		facex[1] = x[1];
		facex[2] = x[2];
		facex[3] = x[3];
		facex[4] = x[10];
		facex[5] = x[11];
		s_crad3x2( facex, crad1262 );
		
		ROUND_NEAR;
		printf("faceinfo = (%d %d)\n", faceinfo[0], faceinfo[1]);
		if( faceinfo[0] < 1 ) {	/* should be small */
			if( crad1262[0] > 2.0 + eps2 ) {
				printf("crad126 is large, should be small.\n");
				printf("crad1262 = [%0.18f\t%0.18f]\n", crad1262[0], crad1262[1]);
			}
			if( crad1262[1] < 2.0 - eps2 ) {
				printf("crad126 is loose.\n");
				printf("crad1262 = [%0.18f\t%0.18f]\n", crad1262[0], crad1262[1]);
			}
		} else {	/* should be large */
			if( crad1262[1] < 2.0 - eps2 ) {
				printf("crad126 is small, should be large.\n");
				printf("crad1262 = [%0.18f\t%0.18f]\n", crad1262[0], crad1262[1]);
			}
			if( crad1262[0] > 2.0 + eps2 ) {
				printf("crad126 is loose.\n");
				printf("crad1262 = [%0.18f\t%0.18f]\n", crad1262[0], crad1262[1]);
			}
		}
		
		/* face:  (1 3 5) -> (0 4 8) */
		facex[2] = x[4];
		facex[3] = x[5];
		facex[4] = x[8];
		facex[5] = x[9];
		s_crad3x2( facex, crad1352 );

		ROUND_NEAR;
		if( faceinfo[1] < 1 ) {	/* should be small */
			if( crad1352[0] > 2.0 + eps2 ) {
				printf("crad135 is large, should be small.\n");
				printf("crad1352 = [%0.18f\t%0.18f]\n", crad1352[0], crad1352[1]);
			}
			if( crad1352[1] < 2.0 - eps2 ) {
				printf("crad135 is loose.\n");
				printf("crad1352 = [%0.18f\t%0.18f]\n", crad1352[0], crad1352[1]);
			}
		} else {	/* should be large */
			if( crad1352[1] < 2.0 - eps2 ) {
				printf("crad135 is small, should be large.\n");
				printf("crad1352 = [%0.18f\t%0.18f]\n", crad1352[0], crad1352[1]);
			}
			if( crad1352[0] > 2.0 + eps2 ) {
				printf("crad135 is loose.\n");
				printf("crad1352 = [%0.18f\t%0.18f]\n", crad1352[0], crad1352[1]);
			}
		}
		
		ROUND_NEAR;
		/*
		printf("taylor_relconst = \n");
		for( i=0; i<7; i++ )
			printf("%.18f\n", taylor_relconst[i]);
		printf("y coefficients = \n");
		for( i=0; i<6; i++ )
			printf("%.18f\n", y_coefficients[i]);
		printf("faceinfo = (%d %d)\n", faceinfo[0], faceinfo[1]);
		printf("original_cell = \n");
		for( i=0; i<12; i+=2 )
			printf("[%0.18f\t%0.18f]\n", original_cell[i], original_cell[i+1]);
		*/
		/*
		printf("tight cell = \n");
		for( i=0; i<12; i+=2 )
			printf("[%0.18f\t%0.18f]\n", y[i], y[i+1]);
		printf("val    = [%0.18f\t%0.18f]\n", val[0], val[1]);
		printf("tomval =  %0.18f\n", -outvals[n]);
		printf("y partials = \n");
		for( i=0; i<12; i+=2 )
			printf("[%0.18f\t%0.18f]\n", partials[i], partials[i+1]);
		printf("\n");
		*/
		
		found = find_lambdas( y, lambdas );
		ROUND_NEAR;
		if( found ) {
			printf("lambdas = (%.18g, %.18g)\n", lambdas[0], lambdas[1]);
		} else {
			printf("Feasible set was empty.\n");
		}
	}
}


void constraint_functions( double x[12], double out[4] )
{
	double facex[6], crad1262[2], crad1352[2];
	
		/* face:  (1 2 6) -> (0 2 10) */
		facex[0] = x[0];
		facex[1] = x[1];
		facex[2] = x[2];
		facex[3] = x[3];
		facex[4] = x[10];
		facex[5] = x[11];
		s_crad3x2( facex, crad1262 );
		
		/* face:  (1 3 5) -> (0 4 8) */
		facex[2] = x[4];
		facex[3] = x[5];
		facex[4] = x[8];
		facex[5] = x[9];
		s_crad3x2( facex, crad1352 );
	if( faceinfo[0] < 1 ) {	/* should be small */
		ROUND_DOWN;
		out[0] = 2.0 - crad1262[1];
		ROUND_UP;
		out[1] = 2.0 - crad1262[0];
	} else {	/* should be large */
		ROUND_DOWN;
		out[0] = crad1262[0] - 2.0;
		ROUND_UP;
		out[1] = crad1262[1] - 2.0;
	}
	if( faceinfo[1] < 1 ) {	/* should be small */
		ROUND_DOWN;
		out[2] = 2.0 - crad1352[1];
		ROUND_UP;
		out[3] = 2.0 - crad1352[0];
	} else {	/* should be large */
		ROUND_DOWN;
		out[2] = crad1352[0] - 2.0;
		ROUND_UP;
		out[3] = crad1352[1] - 2.0;
	}
}


int find_lambdas( double tet[12], double out[2] )
{
	int i, j, k, index, passed, num;
	int n0, n1, n2, n3, n4, n5;
	double points[64][3], x, y, z, pt1[3], pt2[3];
	double corner[2], *pot_corners, *corners;
	double test, yp[12], xp[12], chunk[6], partials[12];
	
	pot_corners = (double *) malloc( sizeof( double )*32*64*2 );
	corners = (double *) malloc( sizeof( double )*32*64*2 );
	ROUND_NEAR;
	
	/* Generate all of the points (64) */
	index = 0;
	for( n0=0; n0<2; n0++ ) {
		y = tet[n0];
		x = y*y;
		yp[0] = y;
		yp[1] = y;
		xp[0] = x;
		xp[1] = x;
		for( n1=2; n1<4; n1++ ) {
			y = tet[n1];
			x = y*y;
			yp[2] = y;
			yp[3] = y;
			xp[2] = x;
			xp[3] = x;
			for( n2=4; n2<6; n2++ ) {
				y = tet[n2];
				x = y*y;
				yp[4] = y;
				yp[5] = y;
				xp[4] = x;
				xp[5] = x;
				for( n3=6; n3<8; n3++ ) {
					y = tet[n3];
					x = y*y;
					yp[6] = y;
					yp[7] = y;
					xp[6] = x;
					xp[7] = x;
					for( n4=8; n4<10; n4++ ) {
						y = tet[n4];
						x = y*y;
						yp[8] = y;
						yp[9] = y;
						xp[8] = x;
						xp[9] = x;
						for( n5=10; n5<12; n5++ ) {
							y = tet[n5];
							x = y*y;
							yp[10] = y;
							yp[11] = y;
							xp[10] = x;
							xp[11] = x;
							octa_score( yp, xp, chunk, partials );
							constraint_functions( xp, chunk + 2 );
							ROUND_NEAR;
							points[index][0] = chunk[0];
							points[index][1] = chunk[2];
							points[index][2] = chunk[4];
							index++;
						}
					}
				}
			}
		}
	}
	
	/*
	for( i=0; i<64; i++ ) {
		printf("points[%d] = (%g, %g, %g)\n", 
			i, points[i][0], points[i][1], points[i][2]);
	}
	*/
	
	/* Generate all of the intersections (2016) */
	index = 0;
	for( i=0; i<63; i++ ) {
		for( k=0; k<3; k++ ) {
			pt1[k] = points[i][k];
		}
		for( j=i+1; j<64; j++ ) {
			for( k=0; k<3; k++ ) {
				pt2[k] = points[j][k];
			}
			intersect( pt1, pt2, corner );
			if( corner[0] > 0.0 && corner[1] > 0.0 ) {
				for( k=0; k<2; k++ ) {
					pot_corners[2*index+k] = corner[k];
				}
				index++;
			}
		}
	}
	num = index;

/*
	for( i=0; i<num; i++ ) {
		printf("pot_corners[%d] = (%g, %g)\n", 
			i, pot_corners[2*i], pot_corners[2*i+1]);
	}
*/
	
	/* Test all of the intersections, to find the corners */
	index = 0;
	for( i=0; i<num; i++ ) {
		for( k=0; k<2; k++ ) {
			corner[k] = pot_corners[2*i+k];
		}		/* Test to see if lambdas are feasible */
		/*
		if( i < 10 )
			printf("corner[%d] = (%g, %g)\n", i, corner[0], corner[1]);
		*/
		passed = 1;
		j = 0;
		while( passed && j < 64 ) {
			test = points[j][0] + points[j][1]*corner[0] +
				points[j][2]*corner[1];
			/*
			if( i < 10 )
				printf("test[%d] = %g\n", j, test);
			*/
			if( test > slop ) {
				passed = 0;
			}
			j++;
		}
		if( passed ) {
			for( k=0; k<2; k++ ) {
				corners[2*index+k] = corner[k];
			}
			index++;
		}
	}
	
	/* Got (index) corners. */
	if( index < 1 ) {		/* Empty feasible set */
		free( pot_corners );
		free( corners );
		return( 0 );
	} else {	/* Compute center. */
		printf("Got %d corners.\n", index);
		
		/*
		for( i=0; i<index; i++ ) {
			printf("corners[%d] = (%g, %g)\n", 
				i, corners[2*i], corners[2*i+1]);
		}
		*/
		
		x = 0.0;
		y = 0.0;
		z = 1.0/index;
		for( i=0; i<index; i++ ) {
			x += corners[2*i];
			y += corners[2*i+1];
		}
		/*
		printf("(x, y, z) = (%g, %g, %g)\n", x, y, z);
		*/
		out[0] = x*z;
		out[1] = y*z;
		free( pot_corners );
		free( corners );
		
		corner[0] = out[0];
		corner[1] = out[1];
		passed = 1;
		j = 0;
		x = 1.0e10;
		while( passed && j < 64 ) {
			test = points[j][0] + points[j][1]*corner[0] +
				points[j][2]*corner[1];
			if( test < x )
				x = test;
			/*
			if( i < 10 )
				printf("test[%d] = %g\n", j, test);
			*/
			if( test > slop ) {
				passed = 0;
			}
			j++;
		}
		if( passed ) {
			printf("Center passed.\n");
		} else {
			printf("Center failed.\n");
		}
		printf("Min test = %g\n", x);
		
		return( 1 );
	}
}


void intersect( double xyz1[3], double xyz2[3], 
	double out[2] )
{
	double det, oodet, u, v;
	
	det = xyz2[1]*xyz1[2] - xyz1[1]*xyz2[2];
	if( fabs( det ) < 1.0e-8 ) {
		out[0] = 0.0;
		out[1] = 0.0;
		return;
	}
	oodet = 1.0/det;
	
	u = (xyz1[0]*xyz2[2] - xyz2[0]*xyz1[2])*oodet;
	v = (xyz2[0]*xyz1[1] - xyz1[0]*xyz2[1])*oodet;
	if( fabs( u ) < 1.0e-12 )
		u = 0.0;
	if( fabs( v ) < 1.0e-12 )
		v = 0.0;
	out[0] = u;
	out[1] = v;
}


