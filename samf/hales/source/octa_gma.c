/* octa_gma.c, by Samuel Ferguson, (c) 1997. */

#include "system_headers.h"
#include "sphere.h"
#include "interval.h"
#include "i_sphere.h"
#include "i_bounds.h"
#include "i_taylor.h"
#include "macros.h"


#define MAXDEPTH	32			/* Max recursion depth (was 10) */
#define PARTIALS	1				/* use partials, or not */

#define SUBDIV		2			/* (4) 	number of subdivisions (1-d) */
#define SUBFRAC	0.5		/* (0.25) 	1/SUBDIV (make this exact . . . */

#define MAXLEN		TWO51_HI
#define SCORECUT 	-0.05758859149521		/* -1.04 pt */
#define CUTBD1		2.716
#define CUTBD2		2.81
#define Y4BOUND		TWO51_HI
#define PEELBD		2.2
#define ALPHA		0.14		/* 0.56/4 */ /* dihedral correction */

/* Even newer bounds:
line1[x_]:= 0.4922197796533495 - 0.3621*x
tetshift[x_]:= 0.253109 - 0.3621*x (eps = 0.07364)
*/

#define INTER1		0.4922197796533495
#define SLOPE1		0.3621

/* External variables */

extern double i_pi_const[2];
extern double i_pi_2_const[2];
extern double i_doct_const[2];
extern double i_two_pi_5_const[2];

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
double twosqrt2;

struct Cell_Info {
	int depth; 
	int cell_type;
	};

/* Prototypes */

void verifyoctagma( void );
void verify_cell( double tet[12], struct Cell_Info info );

/* Code */

void main( void )
{
	ROUND_NEAR;

	printf("Welcome to sphere packing routines, relation testing division.\n");
	verifyoctagma();
}


void verifyoctagma( void )
{
	int i, n;
	double tet[12];
	double diff, dt;
	time_t tstart, tstop;
	clock_t cstart, cstop;
	char *charptr;
	struct Cell_Info cell_info;
	
	printf("Welcome to OctaGma.  ");
#if	PARTIALS
	printf("Using partials to reduce complexity.\n");
#else
	printf("Not using partials.\n");
#endif
	printf("Subdivisions:	%d\n\n", SUBDIV);
	printf("Enter slop:  ");
	scanf("%lf", &slop);
	printf("\t\tSlop = %g\n", slop );
	
	i_init();	/* initialize interval constants */
	
	ROUND_DOWN;
	slop += 0.25*INTER1 + ALPHA*0.5*i_pi_const[0];
	ROUND_NEAR;

	time(&tstart);
	cstart = clock();
	
	charptr = ctime( &tstart );
	printf("Starting time:  \t");
	puts( charptr );
	
	celltotal = 0.0;
	for( n=1; n<6; n++ ) {
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

		for( i=0; i<12; i++ )
			tet[i] = 2.0;
		tet[0] = CUTBD1;
		tet[1] = TWOSQRT2_HI;
		
		switch( n ) {
			case 1:
				tet[1] = CUTBD1;		/* tight */
				tet[7] = Y4BOUND;
				tet[9] = PEELBD;
				tet[11] = PEELBD;
				break;
			case 2:
				tet[1] = CUTBD1;		/* tight */
				tet[3] = PEELBD;
				tet[5] = PEELBD;
				tet[9] = PEELBD;
				tet[11] = PEELBD;
				break;
			case 3:
				tet[3] = PEELBD;
				tet[7] = Y4BOUND;
				tet[9] = PEELBD;
				break;
			case 4:
				tet[3] = PEELBD;
				tet[5] = PEELBD;
				tet[7] = Y4BOUND;
				break;
			case 5:
				tet[3] = PEELBD;
				tet[5] = PEELBD;
				tet[11] = Y4BOUND;
				break;

			default:
				printf("fell off end: this can't happen\n");
				break;
		}
		verify_cell( tet, cell_info );
		printf("cellcount    = %16.0f\t\t", cellcount);
		printf("maxdepth = %d\n", maxdepth);
		printf("sccutcount   = %16.0f\n", sccutcount);
		printf("gmacount     = %16.0f\t\t", gmacount);
		printf("vorcount     = %16.0f\n", vorcount);
		printf("dualcount    = %16.0f\t\t", dualcount);
		printf("verifycount  = %16.0f\n", verifycount);
		printf("partialcount = %16.0f\t\t", partialcount);
		diff = cellcount - verifycount - partialcount - sccutcount;
		printf("diff         = %16.0f\n", diff);
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
	int i, j, bounds[6], n0, n1, n2, n3, n4, n5;
	int check, rel[5];
	double diff[6], t[12], val;
	double sc;
	double x[12], yp[12];
	double relconst[5], rel_val[2];
	double outvals[10];
	
	
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
		
	/* Determine cell type */
	check = info.cell_type;
	if( check == -1 ) {  /* indeterminate */
		check = octa_scoring( x );
		info.cell_type = check;
	}
	
	if( check == 1 ) {
		gmacount +=1.0;
	}
	else if( check == 0 ) {
		vorcount += 1.0;
		return;
	}
	else
		dualcount += 1.0;
	
	/* Score cell properly */
	/* Order:  sol, gma, vor, octavor, dih */
	/* 					0			1		2				3			4  */
	for( i=0; i<5; i++ ) {
		rel[i] = 0;
		relconst[i] = 0.0;
	}
	rel[1] = 1;
	relconst[0] = SLOPE1;
	relconst[1] = 1.0;
	relconst[4] = ALPHA;
	
	/* val = sc + SLOPE1*sph + ALPHA*dih */
	t_composite( rel, relconst, y, x, rel_val, t, 
		outvals );
	
	val = rel_val[1];
	sc = outvals[3];
	
	
	/* Discard cells with low score */
	if( sc < SCORECUT ) {
		sccutcount++;
		return;
	}
	
	/* Attempt verification */	
	if( val < slop ) {
		verifycount += 1.0;
		return;
	}
	
	/* if verification fails, try to bump off, or at least
		reduce dimension */

	for( i=0; i<12; i++ )
		yp[i] = y[i];	/* copy y */
	
#if PARTIALS
	
	/* Score cell properly (find min score) */
	if( check == 1 && outvals[2] > SCORECUT ) {
		/* t contains the relation partials */
		ROUND_NEAR;
		j = 0;
		for( i=0; i<6; i++ ) {
			diff[i] = yp[j+1] - yp[j];	/* compute differentials */
			j += 2;
		}
		/* Long edge is the first one, don't do that here. */
		/* Do edges 2 and 3 */
		j = 1;
		for( i=2; i<6; i+=2 ) {
			n0 = i + 1;
			if( diff[j] > 0.0 ) {	/* ignore tight spots */
				if( t[n0] < 0.0 ) {
					if( yp[i] != 2.0 ) {
						partialcount += 1.0;
						return;	/* bumped off */
					}
					else
						yp[n0] = 2.0;	/* reduced dimension */
				}
				else {
					if( t[i] > 0.0 ) {
						if( yp[n0] != PEELBD ) {
							partialcount += 1.0;
							return;	/* bumped off */
						}
						else
							yp[i] = PEELBD;	/* reduced dimension */				
					}
				}
			}
			j++;
		}
		/* Do edges 5 and 6 */
		j = 4;
		for( i=8; i<12; i+=2 ) {
			n0 = i + 1;
			if( diff[j] > 0.0 ) {	/* ignore tight spots */
				if( t[n0] < 0.0 ) {
					if( yp[i] != 2.0 ) {
						partialcount += 1.0;
						return;	/* bumped off */
					}
					else
						yp[n0] = 2.0;	/* reduced dimension */
				}
				else {
					if( t[i] > 0.0 ) {
						if( yp[n0] != PEELBD ) {
							partialcount += 1.0;
							return;	/* bumped off */
						}
						else
							yp[i] = PEELBD;	/* reduced dimension */				
					}
				}
			}
			j++;
		}
		/* Do edge 4 */
		j = 3;
		i = 6;
		n0 = 7;
		if( diff[j] > 0.0 ) {	/* ignore tight spots */
			if( t[n0] < 0.0 ) {
				if( yp[i] != 2.0 ) {
					partialcount += 1.0;
					return;	/* bumped off */
				}
				else
					yp[n0] = 2.0;	/* reduced dimension */
			}
			else {
				if( t[i] > 0.0 ) {
					if( yp[n0] != Y4BOUND ) {
						partialcount += 1.0;
						return;	/* bumped off */
					}
					else
						yp[i] = Y4BOUND;	/* reduced dimension */				
				}
			}
		}
		/* Now consider the long edge.  */
		/*
		j = 0;
		i = 0;
		n0 = 1;
		*/
		if( diff[0] > 0.0 ) {	/* ignore tight spots */
			if( t[1] < 0.0 ) {
				if( yp[0] != CUTBD1 ) {
					partialcount += 1.0;
					return;	/* bumped off */
				}
				else
					yp[1] = CUTBD1;	/* reduced dimension */
			}
			else {
				if( t[0] > 0.0 ) {
					if( yp[1] != TWOSQRT2_HI ) {
						partialcount += 1.0;
						return;	/* bumped off */
					}
					else
						yp[0] = TWOSQRT2_HI; /* reduced dimension */				
				}
			}
		}
	}
#endif
	if( info.depth < MAXDEPTH ) {
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
		printf("val = %g\n", val - slop);
		printf("sc = [%g, %g]\n", 
			outvals[2], outvals[3]);
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", y[i], y[i+1]);
	}
}

