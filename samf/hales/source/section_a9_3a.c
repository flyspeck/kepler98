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

#define PARTIALS		0
#define REALPARTIALS	1
#define USEDELTA		1

#define TRASH		0
#define DEBUG		0

/* External variables */
/* extern double fake_anc_const[8]; */
extern int rog_dne;

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
int boundary;
double slop;
double original_cell[12];

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
	
	printf("Welcome to Section A_9.3a.\n");
	printf("Enter depth bound:  ");
	scanf("%d", &depthbd);
	printf("Enter slop:  ");
	scanf("%lf", &slop);
	printf("Depth bound = %d\n", depthbd );
	printf("slop = %g\n", slop);
	
	i_init();	/* initialize interval constants */
	appendix_init();
	
					/* initialize local constants */
	ROUND_DOWN;
	/* dt = -0.0290001; */
	dt = -0.022740001;
	slop += dt;

#if TRASH
	appendix_get_trash();
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
		ROUND_DOWN;
		tet[0] = dt/100.0;
		tet[1] = TWOSQRT2_HI;
		/* 2 */
		dt = 245.0;
		ROUND_DOWN;
		tet[2] = 2.0;
		tet[3] = TWO51_HI;
		/* 3 */
		tet[4] = 2.0;
		tet[5] = TWO51_HI;
		/* 4 */
		dt = 32.0;
		ROUND_DOWN;
		tet[6] = dt/10.0;
		tet[7] = tet[6];
		/* 5 */
		tet[8] = 2.0;
		tet[9] = TWO51_HI;
		/* 6 */
		tet[10] = 2.0;
		tet[11] = TWO51_HI;

		switch( n ) {
			case 0:
				dt = 0.001;
				/* 1 */
				ROUND_UP;
				tet[0] = TWO51_LO;
				tet[1] = tet[0] + dt;
				/* 2 */
				dt = 0.001;
				ROUND_DOWN;
				tet[2] = TWO51_LO - dt;
				tet[3] = TWO51_HI;
				/* 3 */
				tet[4] = 2.0;
				tet[5] = tet[4] + dt;
				/* 4 */
				dt = 32.0;
				ROUND_DOWN;
				tet[6] = dt/10.0;
				tet[7] = tet[6];
				/* 5 */
				dt = 0.001;
				tet[8] = 2.0;
				tet[9] = tet[8] + dt;
				/* 6 */
				ROUND_DOWN;
				tet[10] = TWO51_LO - dt;
				tet[11] = TWO51_HI;
				break;
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
		printf("partialpushes= %16.0f\t", gmacount);
		printf("domaincut    = %16.0f\n", vorcount);
		printf("smallcount   = %16.0f\t\t", dualcount);
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
#if USEDELTA
	double delta[2], delta_part[12], xp[12];
#endif
#if REALPARTIALS
	double rel_part[12], h[2], crown_val[2];
	double y6_par[2], face[6], temp[2];
#endif
#if PARTIALS
	double dih26_part[4], rel_part[12];
#endif
	double x[12], yp[12], t[12], data[2], diff[6];
	double dih_val[2], dih_part[12];
	
	/* If part of the verification has already failed, bail. */
	if( failed )
		return;
		
	cellcount += 1.0;
		
	/* adjust domain appropriately */
	ROUND_DOWN;
	for( i=0; i<12; i+=2 ) {
		x[i] = y[i]*y[i];
	}
	ROUND_UP;
	for( i=1; i<12; i+=2 ) {
		x[i] = y[i]*y[i];
	}

#if USEDELTA
	/* Generate best possible bounds for bigdelta */
	i_delta_partials(  x, delta_part );

	/* First find delta_max */
	for( i=0; i<12; i++ )
		xp[i] = x[i];		/* copy x, then modify copy */
	j = 0;
	for( i=0; i<6; i++ ) {
		data[0] = delta_part[j];
		data[1] = delta_part[j+1];
		if( data[0] > 0 )	/* positive sign */
			xp[j] = xp[j+1];
		if( data[1] < 0 )	/* negative sign */
			xp[j+1] = xp[j];
		j += 2;
	}
	delta[1] = rough_max_delta( xp );

	val = delta[1];
	if( val < 0.0 ) {
		sccutcount += 1.0;
		return;
	}
	
	/* Next find delta_min */
	for( i=0; i<12; i++ )
		xp[i] = x[i];		/* copy x, then modify copy */
	j = 0;
	for( i=0; i<6; i++ ) {
		data[0] = delta_part[j];
		data[1] = delta_part[j+1];
		/*	printf("partials data: %f\t%f\n", data[0], data[1]);	*/
		if( data[0] > 0 )	/* positive sign */
			xp[j+1] = xp[j];
		if( data[1] < 0 )	/* negative sign */
			xp[j] = xp[j+1];
		j += 2;
	}
	delta[0] = rough_min_delta( xp );
#endif

	i_dih_partials( x, y, dih_part );
	i_dih_best( x, dih_part, dih_val );
	
	if( dih_val[0] < 1.6499999 ) { /* see Section A_8 */
		dih_val[0] =  1.6499999;
	}
	
	/* Score cell properly */
	/* kappa_fun( y, x, data ); */
	fake_kappa_fun( y, dih_val, data );
	val = data[1];
	if( val < slop ) {
		verifycount += 1.0;
		return;
	}
		/* Check to see if we are in big trouble. */
#if USEDELTA
	if( data[0] > slop && delta[0] > 0.0 ) {
#else
	if( data[0] > slop  ) {
#endif
		failed = 1;
		printf("Cell failed:\n");
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", y[i], y[i+1]);
		printf("val = [%0.18g\t%0.18g]\n", data[0], data[1]);
		return;
	}
	/* if verification fails, try to bump off, or at least
		reduce dimension */

	for( i=0; i<12; i++ )
		yp[i] = y[i];	/* copy y */

#if PARTIALS

	ROUND_NEAR;
	j = 0;
	for( i=0; i<12; i+=2 ) {
		diff[j] = yp[i+1] - yp[i];
		j++;
	}
#if USEDELTA
	if( rog_dne && delta[0] > 0.0 )	{
#else
	if( rog_dne )	{
#endif
		/* Rogers associated with (y1, y2, y6) or (y1, y3, y5) 
		does not exist so we can use easy partials */
		if( rog_dne == 1 || rog_dne == 3 ) {
			/* anchor (y1, y2, y6) does not exist */
			
			/* (1 2 6) -> (0 2 10) */
			dih_sign26( x, dih26_part );
			rel_part[2] = - dih26_part[1];
			rel_part[3] = - dih26_part[0];
			rel_part[10] = - dih26_part[3];
			rel_part[11] = - dih26_part[2];
			
			/* do x2 edge first */
			j = 1;
			i = 2;
			n0 = 3;
			if( diff[j] > 0.0 ) {	/* ignore tight spots */
				if( rel_part[n0] < 0.0 ) {
					if( yp[i] != original_cell[i] ) {
						partialcount += 1.0;
						return;	/* bumped off */
					}
					else {
						yp[n0] = original_cell[i];	/* reduced dimension */
						gmacount += 1.0;
					}
				}
				else {
					if( rel_part[i] > 0.0 ) {
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
			/* now do x6 edge */
			j = 5;
			i = 10;
			n0 = 11;
			if( diff[j] > 0.0 ) {	/* ignore tight spots */
				if( rel_part[n0] < 0.0 ) {
					if( yp[i] != original_cell[i] ) {
						partialcount += 1.0;
						return;	/* bumped off */
					}
					else {
						yp[n0] = original_cell[i];	/* reduced dimension */
						gmacount += 1.0;
					}
				}
				else {
					if( rel_part[i] > 0.0 ) {
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
		if (rog_dne == 2 || rog_dne == 3) {
			/* anchor (y1, y3, y5) does not exist */
			
			/* (1 3 5) -> (0 4 8) */
			dih_sign35( x, dih26_part );
			rel_part[4] = - dih26_part[1];
			rel_part[5] = - dih26_part[0];
			rel_part[8] = - dih26_part[3];
			rel_part[9] = - dih26_part[2];

			/* do x3 edge first */
			j = 2;
			i = 4;
			n0 = 5;
			if( diff[j] > 0.0 ) {	/* ignore tight spots */
				if( rel_part[n0] < 0.0 ) {
					if( yp[i] != original_cell[i] ) {
						partialcount += 1.0;
						return;	/* bumped off */
					}
					else {
						yp[n0] = original_cell[i];	/* reduced dimension */
						gmacount += 1.0;
					}
				}
				else {
					if( rel_part[i] > 0.0 ) {
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
			/* now do x5 edge */
			j = 4;
			i = 8;
			n0 = 9;
			if( diff[j] > 0.0 ) {	/* ignore tight spots */
				if( rel_part[n0] < 0.0 ) {
					if( yp[i] != original_cell[i] ) {
						partialcount += 1.0;
						return;	/* bumped off */
					}
					else {
						yp[n0] = original_cell[i];	/* reduced dimension */
						gmacount += 1.0;
					}
				}
				else {
					if( rel_part[i] > 0.0 ) {
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

#if REALPARTIALS
	ROUND_NEAR;
	j = 0;
	for( i=0; i<6; i++ ) {	/* differentials */
		diff[i] = yp[j+1] - yp[j];
		j += 2;
	}
#if USEDELTA
	if( info.depth > 1 && delta[0] > 0.0 )	{
#else
	if( info.depth > 1 ) {
#endif
		h[0] = 0.5*y[0];
		h[1] = 0.5*y[1];
		crown_fun( h, crown_val );
		/* Use partials for y5, y6, for simplicity */
		
		/* y5 */
		/* face:  (1 3 5) -> (0 4 8) */
		face[0] = y[0];
		face[1] = y[1];
		face[2] = y[4];
		face[3] = y[5];
		face[4] = y[8];
		face[5] = y[9];
/*
		facex[0] = x[0];
		facex[1] = x[1];
		facex[2] = x[4];
		facex[3] = x[5];
		facex[4] = x[8];
		facex[5] = x[9];
		anc_y6par( face, facex, y6_par );
*/
		fake_anchor_y6( face, y6_par );

		i = 8;
		i_mult( crown_val, dih_part + i, temp );
		i_add( temp, y6_par, rel_part + i );

		i = 4;
		i_mult( crown_val, dih_part + i, temp );
		i_add( temp, y6_par, rel_part + i );
		
		/* y6 */
		/* face:  (1 2 6) -> (0 2 10) */
		face[2] = y[2];
		face[3] = y[3];
		face[4] = y[10];
		face[5] = y[11];
/*
		facex[2] = x[2];
		facex[3] = x[3];
		facex[4] = x[10];
		facex[5] = x[11];
		anc_y6par( face, facex, y6_par );
*/
		fake_anchor_y6( face, y6_par );

		i = 10;
		i_mult( crown_val, dih_part + i, temp );
		i_add( temp, y6_par, rel_part + i );
					
		i = 2;
		i_mult( crown_val, dih_part + i, temp );
		i_add( temp, y6_par, rel_part + i );
		
		for( j=1; j<6; j++ ) {
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
					for( n3=0; n3<bounds[3]; n3++ ) {
						t[6] = yp[6];
						t[7] = yp[7];
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
		printf("val    = [%0.18g\t%0.18g]\n", data[0], data[1]);
#if USEDELTA
		printf("delta  = [%0.18g\t%0.18g]\n", delta[0], delta[1]);
#endif
		i_dih( x, t );
		ROUND_NEAR;
		printf("i_dih  = [%0.18g\t%0.18g]\n", t[0], t[1]);
#if REALPARTIALS
		i = 8;
		printf("rel_part  = [%0.18g\t%0.18g]\n", rel_part[i], rel_part[i+1]);
		i = 10;
		printf("rel_part  = [%0.18g\t%0.18g]\n", rel_part[i], rel_part[i+1]);
#endif		
#if PARTIALS
		printf("rog_dne = %d\n", rog_dne );
		if( rog_dne == 1 || rog_dne == 3 ) {
			i = 2;
			printf("rel_part_%d  = [%0.18g\t%0.18g]\n", i/2+1, rel_part[i], rel_part[i+1]);
			i = 10;
			printf("rel_part_%d  = [%0.18g\t%0.18g]\n", i/2+1, rel_part[i], rel_part[i+1]);
		}
		if( rog_dne == 2 || rog_dne == 3 ) {
			i = 4;
			printf("rel_part_%d  = [%0.18g\t%0.18g]\n", i/2+1, rel_part[i], rel_part[i+1]);
			i = 8;
			printf("rel_part_%d  = [%0.18g\t%0.18g]\n", i/2+1, rel_part[i], rel_part[i+1]);
		}
#endif		
		printf("y = \n");
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", y[i], y[i+1]);
#if PARTIALS
		printf("yp = \n");
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", yp[i], yp[i+1]);
#endif
#if REALPARTIALS
		printf("yp = \n");
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", yp[i], yp[i+1]);
#endif
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

2.51 2.5158125
2.1115625 2.1275
2.0 2.0159375
2.769999999999999574 2.828427124746190291
2.031874999999999876 2.047812500000000036
2.239062500000000178 2.255

*/
