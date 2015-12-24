/* octahedra.c, by Samuel Ferguson, (c) 1998. */

#include "system_headers.h"
#include "interval.h"
#include "i_sphere.h"
#include "i_bounds.h"
#include "macros.h"


#define MAXDEPTH	32			/* Max recursion depth (was 10) */

#define SUBDIV		2			/* (4) 	number of subdivisions (1-d) */
#define SUBFRAC	0.5		/* (0.25) 	1/SUBDIV (make this exact . . .) */

#define COUNTLIMIT	1.0e5		/* 4.0e7 */

#define SEAN					0
#define SMARTSUBDIV		1
#define PARTIALS			1
#define PARTIALSTART	1


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
double original_cell[12];
double countlimit;
double dih_bounds[2];
double sean_val[2];
int do_max;

struct Cell_Info {
	int depth; 
	int diff_skip[6];
	};

/* Prototypes */

void verifyocta( void );
void verify_cell( double tet[12], struct Cell_Info info );
void set_bounds( int casenum, double out_tet[12],
	int out_caselist[2] );
void testfuns( void );

/* Code */

void main( void )
{
	ROUND_NEAR;

	printf("Welcome to sphere packing routines, relation testing division.\n");
	verifyocta();
	/*
	verifyocta();
	testfuns();
	*/
}


void verifyocta( void )
{
	int i, m, n, casestart, casebound, caselist[2];
	double tet[12], starttet[12], tetdiff[6], max;
	double diff, dt, three02[2];
	double slop;
	time_t tstart, tstop;
	clock_t cstart, cstop;
	char *charptr;
	struct Cell_Info cell_info;
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
#else
	sean_val[0] = TWO51_LO;
	sean_val[1] = TWO51_HI;
#endif

	dt   = 302.0;
	diff = 100.0;
	ROUND_DOWN;
	three02[0] = dt/diff;
	ROUND_UP;
	dt   = 302.0;
	diff = 100.0;
	three02[1] = dt/diff;
	
	onefailed = 0;
	
	ROUND_NEAR;
	if( sean_val[1] <= sean_val[0] )
		printf("sean_val = [%.20f, %.20f]\n", sean_val[0], sean_val[1]);
	if( three02[1] <= three02[0] )
		printf("3.02 = [%.20f, %.20f]\n", three02[0], three02[1]);

	printf("Welcome to Dihedral Bounds.  \n");
	printf("Enter cell count limit:  ");
	scanf("%lf", &countlimit);

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
	
	for( n=0; n < 14; n++ ) { 
		
		printf("Starting case %d:\n", n);
		set_bounds( n, starttet, caselist );
		ROUND_DOWN;
		dih_bounds[0] = dih_bounds[0] - slop;
		ROUND_UP;
		dih_bounds[1] = dih_bounds[1] + slop;
		ROUND_NEAR;
		for( i=0; i<12; i++ ) {
			tet[i] = starttet[i];
		}
		if( caselist[0] == 1 ) {
			casestart = 1;
		} else {
			casestart = 2;
		}
		if( caselist[1] == 1 ) {
			casebound = 2;
		} else {
			casebound = 1;
		}
		
		for( subcase=casestart; subcase<=casebound; subcase++ ) {
			
			for( i=0; i<12; i++ ) {
				tet[i] = starttet[i];
			}
			if( subcase==1 ) {	/* min */
				do_max = 0;
				printf("\tVerifying dih > %.18f\n", dih_bounds[0]);
				/* dih is increasing in y4 */
				tet[7] = starttet[6];
			} else {						/* max */
				do_max = 1;
				printf("\tVerifying dih < %.18f\n", dih_bounds[1]);
				tet[6] = starttet[7];
				if( n==2 ) {
					tet[11] = three02[1];
				}
				if( n==6 ) {
					tet[9]  = three02[1];
					tet[11] = three02[1];
				}
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

			verify_cell( tet, cell_info );
			
			ROUND_NEAR;
			
			printf("cellcount    = %16.0f\t\t", cellcount);
			printf("maxdepth = %d\n", maxdepth);
			printf("partialpushes= %16.0f\n", gmacount);
			printf("largecut     = %16.0f\t\t", dualcount);
			printf("discards     = %16.0f\n", vorcount);
			printf("smallcut     = %16.0f\t\t", sccutcount);
			printf("verifycount  = %16.0f\n", verifycount);
			printf("partialcount = %16.0f\t\t", partialcount);
			diff = cellcount - verifycount - partialcount - vorcount;
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
	double data[2], delta[2], delta_part[12];
	double partials[12], rel_part[12];
	
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
	for( i=0; i<12; i+=2 ) {
		x[i] = y[i]*y[i];
	}

	/* First try to verify over the current cell.  */
	i_bigdelta_partials_best( x, delta, delta_part );
	
	if( delta[1] < 0.0 ) {	/* discard */
		vorcount += 1.0;
		return;
	}
		
	/* Score cell properly */
	i_sign_dihpars( x, delta, delta_part, partials );
	i_dih_best( x, partials, data );
	
	/* Attempt verification */
	if( do_max ) {
		if( data[1] < dih_bounds[1] ) {
			verifycount += 1.0;
			return;
		}
		if( data[0] > dih_bounds[1] && delta[0] > 0.0 ) {
			/* We are in big trouble. */
			failed = 1;
			printf("Cell failed:\n");
			for( i=0; i<12; i+=2 )
				printf("%0.18f\t%0.18f\n", y[i], y[i+1]);
			printf("data  = [%0.18g\t%0.18g]\n", data[0], data[1]);
			printf("delta = [%0.18g\t%0.18g]\n", delta[0], delta[1]);
			return;
		}
	} else {
		if( data[0] > dih_bounds[0] ) {
			verifycount += 1.0;
			return;
		}
		if( data[1] < dih_bounds[0] && delta[0] > 0.0 ) {
			/* We are in big trouble. */
			failed = 1;
			printf("Cell failed:\n");
			for( i=0; i<12; i+=2 )
				printf("%0.18f\t%0.18f\n", y[i], y[i+1]);
			printf("data  = [%0.18g\t%0.18g]\n", data[0], data[1]);
			printf("delta = [%0.18g\t%0.18g]\n", delta[0], delta[1]);
			return;
		}
	}
	
	/* if verification fails, try to bump off, or at least
		reduce dimension */

	for( i=0; i<12; i++ )
		yp[i] = y[i];	/* copy y */
	
#if PARTIALS
	if( delta[0] > 0.0 )	{
		if( do_max ) {
			for( i=0; i<12; i++ ) {
				rel_part[i] = partials[i];
			}
		} else {
			for( i=0; i<12; i+=2 ) {
				j = i + 1;
				rel_part[i] = -partials[j];
			}
			for( i=1; i<12; i+=2 ) {
				j = i - 1;
				rel_part[i] = -partials[j];
			}
		}
		
		ROUND_NEAR;
		j = 0;
		for( i=0; i<6; i++ ) {	/* differentials */
			diff[i] = yp[j+1] - yp[j];
			j += 2;
		}
		
		for( j=0; j<6; j++ ) {
			if( j != 3 ) {
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
		}
	}
#endif

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
		printf("data  = [%0.18g\t%0.18g]\n", data[0], data[1]);
		printf("delta = [%0.18g\t%0.18g]\n", delta[0], delta[1]);
		printf("y = \n");
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", y[i], y[i+1]);
		if( delta[0] > 0.0 )	{
			printf("partials = \n");
			for( i=0; i<12; i+=2 )
				printf("%0.18f\t%0.18f\n", partials[i], partials[i+1]);
			printf("yp = \n");
			for( i=0; i<12; i+=2 )
				printf("%0.18f\t%0.18f\n", yp[i], yp[i+1]);
		}
	}
}


void set_bounds( int casenum, double tet[12],
	int out_caselist[2] )
{
	int cases[14][3] = 
		{{0, 0, 1},
		 {0, 0, 3},
		 {0, 1, 0},
		 {0, 1, 1},
		 {0, 1, 2},
		 {0, 1, 3},
		 {1, 1, 0},
		 {1, 1, 1},
		 {1, 1, 2},
		 {1, 1, 3},
		 {4, 0, 4},
		 {4, 0, 5},
		 {4, 1, 4},
		 {4, 1, 5}};
	int casebd[14][2] = 
		{{1, 1},
		 {1, 0},
		 {1, 1},
		 {1, 1},
		 {1, 0},
		 {1, 0},
		 {1, 1},
		 {1, 1},
		 {1, 0},
		 {1, 0},
		 {1, 1},
		 {1, 0},
		 {0, 1},
		 {1, 0}};
	double dih_list[14][2] =
		{{1.153,	2.28},
		 {1.32,		8.0},
		 {0.633,	1.624},
		 {1.033,	1.929},
		 {1.033,	8.0},
		 {1.259,	8.0},
		 {0.817,	1.507},
		 {1.07,		1.761},
		 {1.07,		8.0},
		 {1.23,		8.0},
		 {0.633,	1.624},
		 {1.033,	8.0},
		 {0.0,		1.381},
		 {0.777,	8.0}};
	int i, j, k;
	
	if( casenum > 13 ) {
		printf("Bad casenum.\n");
		return;
	}
	
	for( i=0; i<6; i+=2 ) {
		tet[i] = 2.0;
		tet[i+1] = sean_val[1];
	}
	
	for( i=0; i<3; i++ ) {
		j = cases[casenum][i];
		k = 2*i + 8;
		if( i==2 ) {
			k = 6;
		}
		switch( j ) {
			case 0:
				tet[k] = 2.0;
				tet[k+1] = sean_val[1];
				break;
			case 1:
				tet[k] = sean_val[0];
				tet[k+1] = TWOSQRT2_HI;
				break;
			case 2:
				tet[k] = sean_val[0];
				ROUND_UP;
				tet[k+1] = tet[3] + tet[5];
				break;
			case 3:
				tet[k] = TWOSQRT2_LO;
				ROUND_UP;
				tet[k+1] = tet[3] + tet[5];
				break;
			case 4:
				tet[k] = 2.0;
				tet[k+1] = sean_val[1];
				tet[4] = sean_val[0];
				tet[5] = TWOSQRT2_HI;
				break;
			case 5:
				tet[k] = sean_val[0];
				ROUND_UP;
				tet[k+1] = tet[3] +tet[5];
				tet[4] = sean_val[0];
				tet[5] = TWOSQRT2_HI;
				break;
			default:
				printf("This can't happen.\n");
				break;
		}
	}
	dih_bounds[0] = dih_list[casenum][0];
	dih_bounds[1] = dih_list[casenum][1];
	
	out_caselist[0] = casebd[casenum][0];
	out_caselist[1] = casebd[casenum][1];
	/* Handle extra special cases elsewhere. */
}


void testfuns( void )
{
	int i, j;
	double y[12], x[12], temp[2];
	double delta[2], sqrtdelta[2], delta32[2], delta_part[12];
	double dih_part[12], partials[12];
		
	i_init();

	j = 0;
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
		
		i_bigdelta_partials_best( x, delta, delta_part );
		ROUND_NEAR;
		printf("delta = [%.18f, %.18f]\n", delta[0], delta[1]);
		i_sqrt( delta, sqrtdelta );
		i_mult( delta, sqrtdelta, delta32 );

		i_sign_dihpars( x, delta, delta_part, partials );
		ROUND_NEAR;
		printf("sign_dihpars = \n");
		for( i=0; i<12; i+=2 ) {
			printf("[%.18f, %.18f]\n", partials[i], partials[i+1]);
		}
		i_dih( x, temp );
		ROUND_NEAR;
		printf("i_dih      = [%.18f, %.18f]\n", temp[0], temp[1]);
		i_dih_alt( x, delta, temp );
		ROUND_NEAR;
		printf("i_dih_alt  = [%.18f, %.18f]\n", temp[0], temp[1]);
		i_dih_old( x, temp );
		ROUND_NEAR;
		printf("i_dih_old  = [%.18f, %.18f]\n", temp[0], temp[1]);
		i_dih_xpars( x, dih_part );
		i_dih_best( x, dih_part, temp );
		ROUND_NEAR;
		printf("i_dih_best = [%.18f, %.18f]\n", temp[0], temp[1]);
	}
}

/*
2.0    2.06375
2.0    2.06375
2.44625    2.51
2.828427124746189847    2.828427124746189847
2.0    2.06375
2.748820343559642332    2.828427124746190291

2.0		2.0
2.0		2.0
2.45	2.45
2.828427124746189847    2.828427124746189847
2.0		2.0
2.75	2.75

2.041				2.041
2.0291			2.0291
2.03412			2.03412
2.154				2.154
2.52331			2.52331
2.5123			2.5123

*/

