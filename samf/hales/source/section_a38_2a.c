/* octahedra.c, by Samuel Ferguson, (c) 1998. */

#include "system_headers.h"
#include "interval.h"
#include "i_sphere.h"
#include "i_bounds.h"
#include "i_taylor.h"
#include "i_appendix.h"
#include "i_voronoi.h"
#include "macros.h"


#define MAXDEPTH	32			/* Max recursion depth (was 10) */

#define SUBDIV		2			/* (4) 	number of subdivisions (1-d) */
#define SUBFRAC	0.5		/* (0.25) 	1/SUBDIV (make this exact . . .) */

#define SMARTSUBDIV		1
#define PARTIALS			0
#define PARTIALSTART	3


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
double original_cell[12];
double countlimit;
double sean_val[2];
double slop;
double origslop;
double y_coefficients[16];

struct Cell_Info {
	int depth; 
	int diff_skip[6];
	};

/* Prototypes */

void verifyocta( void );
void verify_cell( double tet[12], struct Cell_Info info );
void testfuns( void );
void sean_init( void );
void set_relation( int num );
void set_domain( int num, double tet[12] );
void octa_score( double y[12], double x[12],
	double out[2], double out_partials[12] );
void test_recognize( void );

/* Code */

void main( void )
{
	ROUND_NEAR;

	printf("Welcome to sphere packing routines, relation testing division.\n");
	verifyocta();
	/*
	verifyocta();
	testfuns();
	test_recognize();
	*/
}


void verifyocta( void )
{
	int i, m, n;
	double tet[12], tetdiff[6], max;
	double diff, dt;
	time_t tstart, tstop;
	clock_t cstart, cstop;
	char *charptr;
	struct Cell_Info cell_info;
	int onefailed, subcase;

	sean_val[0] = TWO51_LO;
	sean_val[1] = TWO51_HI;
	
	onefailed = 0;
	
	printf("Welcome to Section A.3.8.2a.  \n");
	/* dihedral boundary case */
	printf("Enter cell count limit:  ");
	scanf("%lf", &countlimit);

	printf("Enter slop:  ");
	scanf("%lf", &slop);
	printf("\t\tSlop = %g\n", slop );
	origslop = slop;
	
	i_init();	/* initialize interval constants */
	sean_init();
		
	ROUND_NEAR;

	time(&tstart);
	cstart = clock();
	
	charptr = ctime( &tstart );
	printf("Starting time:  \t");
	puts( charptr );
	
	celltotal = 0.0;
	
	for( n=1; n < 6; n++ ) { 	

		printf("Case %d:\n", n);
					
		set_relation( n );
		/*
		ROUND_NEAR;
		printf("y_coefficients = \n");
		for( i=0; i<16; i+=2 ) {
			printf("[%.18f, %.18f]\n", 
				y_coefficients[i], y_coefficients[i+1]);
		}
		*/
		ROUND_DOWN;
		slop = origslop + y_coefficients[12];
		
		ROUND_NEAR;
		printf("rhs = %.18f\n", y_coefficients[12]);

		/* turn vor0 into -tau0 */
		ROUND_DOWN;
		sol_coeff[0] = phi0_val[0] - zetapt_val[1];
		ROUND_UP;
		sol_coeff[1] = phi0_val[1] - zetapt_val[0];							

		for( subcase=1; subcase<8; subcase++ ) {
			
			set_domain( n, tet );
			/*
			ROUND_NEAR;
			printf("domain = \n");
			for( i=0; i<12; i+=2 ) {
				printf("[%.18f, %.18f]\n", 
					tet[i], tet[i+1]);
			}
			*/
			/*
			ROUND_NEAR;
			printf("slop = %.18f\n", slop);
			*/
			
			switch( subcase ) {
				case 1:
					tet[1] = tet[0];
					tet[3] = tet[2];
					tet[5] = tet[4];
					break;
				case 2:
					tet[3] = tet[2];
					tet[9] = tet[8];
					break;
				case 3:
					tet[5] = tet[4];
					tet[11] = tet[10];
					break;
				case 4:
					tet[7] = tet[6];
					tet[11] = tet[10];
					break;
				case 5:
					tet[9] = tet[8];
					tet[11] = tet[10];
					break;
				case 6:
					tet[1] = tet[0];
					tet[7] = tet[6];
					break;
				case 7:
					tet[7] = tet[6];
					tet[9] = tet[8];
					break;
				default:
					break;
			}
			
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
			diff = cellcount - verifycount - partialcount - 
				dualcount - sccutcount;
			printf("diff         = %16.0f\n", diff);
			if( !failed )
				printf("Verification succeeded.\n\n");
			else {
				printf("Verification FAILED.\n\n");
				onefailed++;
			}
			celltotal += cellcount;
		}	/* end subcase loop */
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
	double x[12], ynew[12];
	double data[2], delta[2], sqrtdelta[2], vor0_val[2];
	double delta_part[12], partials[12];
	double dih_val[2];
#if PARTIALS
	double rel_part[12];
#endif
	
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
		
	/* Score cell properly */
	i_bigdelta_partials_best( x, delta, delta_part );

	i_sign_dihpars( x, delta, delta_part, partials );
	i_dih_best( x, partials, dih_val );
	
	if( dih_val[1] < y_coefficients[0] ) {
		sccutcount += 1.0;	/* discard */
		return;
	}
	if( dih_val[0] > y_coefficients[1] ) {
		dualcount += 1.0;	/* discard */
		return;
	}
	
	i_sqrt( delta, sqrtdelta );
	i_vor0( y, x, sqrtdelta, vor0_val );
	
	if( dih_val[0] < 0.0 )
		dih_val[0] = 0.0;
	
	ROUND_DOWN;
	data[0] = vor0_val[0] + y_coefficients[14]*dih_val[0];
	ROUND_UP;
	data[1] = vor0_val[1] + y_coefficients[15]*dih_val[1];
	
	if( delta[1] < 0.0 ) {
		vorcount += 1.0;	/* discard */
		return;
	}
	
	/* Attempt verification */
	if( data[1] < slop ) {
		verifycount += 1.0;
		return;
	}
	
	if( data[0] > slop && dih_val[0] > y_coefficients[1] ) {
		/* We are in big trouble. */
		failed = 1;
		printf("Cell failed:\n");
		printf("y = \n");
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", y[i], y[i+1]);
		printf("data        = [%0.18g\t%0.18g]\n", 
			data[0], data[1]);
		printf("delta       = [%0.18g\t%0.18g]\n", 
			delta[0], delta[1]);
		printf("vor0_val    = [%0.18g\t%0.18g]\n", 
			vor0_val[0], vor0_val[1]);
		printf("dih_val     = [%0.18g\t%0.18g]\n", 
			dih_val[0], dih_val[1]);
		return;
	}
	
	/* if verification fails, try to bump off, or at least
		reduce dimension */

	for( i=0; i<12; i++ )
		ynew[i] = y[i];	/* copy y */
	
#if PARTIALS
	/* do delta==0 separately */
	if( dih_val[0] > y_coefficients[1] )	{
	
		i_vor0_partials( y, x, delta_part );
		i_dih_partials( x, y, partials );
		for( i=0; i<12; i+=2 ) {
			i_mult( y_coefficients + 14, partials + i, t + i );
		}
		ROUND_DOWN;
		for( i=0; i<12; i+=2 ) {
			rel_part[i] = delta_part[i] + t[i];
		}
		ROUND_UP;
		for( i=1; i<12; i+=2 ) {
			rel_part[i] = delta_part[i] + t[i];
		}

		ROUND_NEAR;
		j = 0;
		for( i=0; i<6; i++ ) {	/* differentials */
			diff[i] = ynew[j+1] - ynew[j];
			j += 2;
		}
		
		for( j=0; j<6; j++ ) {
			if( diff[j] > 0.0 ) {	/* ignore tight spots */
				i = 2*j;
				n0 = i+1;
				if( rel_part[n0] < 0.0 ) {
					if( ynew[i] != original_cell[i] ) {
						partialcount += 1.0;
						return;	/* bumped off */
					}
					else {
						ynew[n0] = original_cell[i];	/* reduced dimension */
						gmacount += 1.0;
					}
				} else if( rel_part[i] > 0.0 ) {
					if( ynew[n0] != original_cell[n0] ) {
						partialcount += 1.0;
						return;	/* bumped off */
					}
					else {
						ynew[i] = original_cell[n0];	/* reduced dimension */
						gmacount += 1.0;
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
			diff[i] = ynew[j+1] - ynew[j];
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
			diff[i] = SUBFRAC*(ynew[j+1] - ynew[j]);
			if( diff[i] > 1.0e-14 )
				bounds[i] = SUBDIV;
			else {
				bounds[i] = 1;
				diff[i] = ynew[j+1] - ynew[j];
			}
			j += 2;
		}
#endif

		for( n0=0; n0<bounds[0]; n0++ ) {
			val = ynew[0];
			t[0] = val + n0*diff[0];
			t[1] = val + (n0+1)*diff[0];
			for( n1=0; n1<bounds[1]; n1++ ) {
				val = ynew[2];
				t[2] = val + n1*diff[1];
				t[3] = val + (n1+1)*diff[1];
				for( n2=0; n2<bounds[2]; n2++ ) {
					val = ynew[4];
					t[4] = val + n2*diff[2];
					t[5] = val + (n2+1)*diff[2];
					for( n3=0; n3<bounds[3]; n3++ ) {
						val = ynew[6];
						t[6] = val + n3*diff[3];
						t[7] = val + (n3+1)*diff[3];
						for( n4=0; n4<bounds[4]; n4++ ) {
							val = ynew[8];
							t[8] = val + n4*diff[4];
							t[9] = val + (n4+1)*diff[4];
							for( n5=0; n5<bounds[5]; n5++ ) {
								val = ynew[10];
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
		printf("data        = [%0.18g\t%0.18g]\n", 
			data[0], data[1]);
		printf("delta       = [%0.18g\t%0.18g]\n", 
			delta[0], delta[1]);
		printf("vor0_val    = [%0.18g\t%0.18g]\n", 
			vor0_val[0], vor0_val[1]);
		printf("dih_val     = [%0.18g\t%0.18g]\n", 
			dih_val[0], dih_val[1]);
		printf("y = \n");
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", y[i], y[i+1]);
#if PARTIALS
	if( dih_val[0] > y_coefficients[1] )	{
			printf("partials = \n");
			for( i=0; i<12; i+=2 ) {
				printf("[%0.18f\t%0.18f]\n", 
					rel_part[i], rel_part[i+1]);
			}
			printf("ynew = \n");
			for( i=0; i<12; i+=2 )
				printf("%0.18f\t%0.18f\n", ynew[i], ynew[i+1]);
		}
#endif
	}
}


void testfuns( void )
{
	int i, j;
	double y[12], yp[12], x[12], xp[12], temp[2];
	double temp2[2];
	double delta[2], sqrtdelta[2], delta_part[12];
	double partials[12];
		
	i_init();
	
	sean_val[0] = TWO51_LO;
	sean_val[1] = TWO51_HI;
	ROUND_NEAR;
	if( sean_val[1] <= sean_val[0] ) {
		printf("sean_val = [%.20f, %.20f]\n", 
			sean_val[0], sean_val[1]);
	}
	
	sean_init();

	j = 0;
	while( 1 ) {
		ROUND_NEAR;
		printf("Enter interval values for y:\n");
		for( i=0; i<12; i++ )
			scanf("%lf", y + i);
/*		
		printf("Enter interval values for yp:\n");
		for( i=0; i<12; i++ )
			scanf("%lf", yp + i);
*/		
		printf("Got y = \n");
		for( i=0; i<12; i+=2 )
			printf("[%.18f, %.18f]\n", y[i], y[i+1]);
/*
		printf("Got yp = \n");
		for( i=0; i<12; i+=2 )
			printf("[%.18f, %.18f]\n", yp[i], yp[i+1]);
*/
		ROUND_DOWN;
		for( i=0; i<12; i+=2 )
			x[i] = y[i]*y[i];
		ROUND_UP;
		for( i=1; i<12; i+=2 )
			x[i] = y[i]*y[i];
		
		for( i=0; i<12; i++ )
			yp[i] = y[i];
		
		ROUND_DOWN;
		for( i=0; i<12; i+=2 )
			xp[i] = yp[i]*yp[i];
		ROUND_UP;
		for( i=1; i<12; i+=2 )
			xp[i] = yp[i]*yp[i];

		cross_diag2( x, xp, temp2 );
		/*
		ROUND_NEAR;
		printf("cross_diag2 = [%.18f, %.18f]\n", 
			temp2[0], temp2[1]);
		i_sqrt( temp2, temp );
		ROUND_NEAR;
		printf("cross_diag  = [%.18f, %.18f]\n", 
			temp[0], temp[1]);
		*/

		i_bigdelta_partials_best( x, delta, delta_part );
		i_sqrt( delta, sqrtdelta );
		i_vor0( y, x, sqrtdelta, temp );
		
		i_delta4( x, temp );
		printf("delta = [%.18f, %.18f]\n", 
			delta[0], delta[1]);
		printf("delta4 = [%.18f, %.18f]\n", 
			temp[0], temp[1]);
		
		i_dih( x, temp );
		ROUND_NEAR;
		printf("i_dih = [%.18f, %.18f]\n", 
			temp[0], temp[1]);

		i_dih_alt( x, delta, temp );
		ROUND_NEAR;
		printf("i_dih_alt = [%.18f, %.18f]\n", 
			temp[0], temp[1]);

		i_sign_dihpars( x, delta, delta_part, partials );
		i_dih_best( x, partials, temp );
		ROUND_NEAR;
		printf("i_dih_best = [%.18f, %.18f]\n", 
			temp[0], temp[1]);


		/*
		printf("vor0_val  = [%.18f, %.18f]\n", 
			temp[0], temp[1]);

		i_tomsrho_xpars( x, delta_part );
		ROUND_NEAR;
		printf("rho_pars = \n");
		for( i=0; i<12; i+=2 ) {
			printf("[%.18f, %.18f]\n", 
				delta_part[i], delta_part[i+1]);
		}
		i_crad2( x, temp );
		ROUND_NEAR;
		printf("crad2      = [%.18f, %.18f]\n", 
			temp[0], temp[1]);
		*/
		
		set_relation( 1 );
		
		i_crad2_best( x, temp );
		ROUND_NEAR;
		printf("crad2_best = [%.18f, %.18f]\n", 
			temp[0], temp[1]);
		octa_score( y, x, temp, delta_part );
		ROUND_NEAR;
		printf("rel_val = [%.18f, %.18f]\n", temp[0], temp[1]);
		printf("rel_pars = \n");
		for( i=0; i<12; i+=2 ) {
			printf("[%.18f, %.18f]\n", 
				delta_part[i], delta_part[i+1]);
		}

	}
}


void sean_init( void )
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


void set_relation( int num )
{
	int i, j;
	double bigmat[5][8] = 
		{{1.2, 0.0, 0.0, 0.0, 0.0, 0.0, 
				0.2391, 0.2529},
		 {1.2, 0.0, 0.0, 0.0, 0.0, 0.0, 
				0.1376, 0.2529},
		 {1.2, 0.0, 0.0, 0.0, 0.0, 0.0, 
				0.266, 0.2529},
		 {1.2, 0.0, 0.0, 0.0, 0.0, 0.0, 
				0.12, 0.2529},
		 {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
				0.1453, 0.2529}};
	
	if( num > 5 ) {
		printf("num > 5.\n");
		return;
	}
	if( num < 1 ) {
		printf("num < 1.\n");
		return;
	}
	
	i = num - 1;
	for( j=0; j<8; j++ ) {
		i_recognize( bigmat[i][j], y_coefficients + 2*j );
	}
}


void set_domain( int num, double tet[12] )
{
	int i;
	double temp[2];
	
	if( num > 5 ) {
		printf("num > 5.\n");
		return;
	}
	if( num < 1 ) {
		printf("num < 1.\n");
		return;
	}
	
	for( i=0; i<12; i+=2 ) {
		tet[i] = 2.0;
		tet[i+1] = sean_val[1];
	}
	
	i_recognize( 2.168, temp );
	
	tet[3] = temp[1];
	tet[5] = temp[1];
	
	switch( num ) {
		case 1:
			tet[10] = sean_val[0];
			ROUND_UP;
			tet[11] = tet[1] + tet[3];
			break;
		case 2:
			tet[8] = sean_val[0];
			tet[9] = sean_val[1];
			tet[10] = sean_val[0];
			ROUND_UP;
			tet[11] = tet[1] + tet[3];
			break;
		case 3:
			tet[6] = sean_val[0];
			tet[7] = TWOSQRT2_HI;
			tet[10] = sean_val[0];
			ROUND_UP;
			tet[11] = tet[1] + tet[3];
			break;
		case 4:
			tet[6] = sean_val[0];
			tet[7] = TWOSQRT2_HI;
			tet[8] = sean_val[0];
			tet[9] = sean_val[1];
			tet[10] = sean_val[0];
			ROUND_UP;
			tet[11] = tet[1] + tet[3];
			break;
		case 5:
			tet[8] = sean_val[0];
			i_recognize( 3.488, temp );
			tet[9] = temp[1];
			tet[10] = sean_val[0];
			tet[11] = sean_val[1];
			break;
	}
	
}


void octa_score( double y[12], double x[12],
	double out[2], double out_partials[12] )
{
	int j, jp;
	double t[12], compvals[14], delta_val[2], dih_val[2];
	
	t[0] = y[0];
	i_bigdelta_best( x, delta_val );
	i_dih_alt( x, delta_val, dih_val );
	i_mult( y_coefficients + 14, dih_val, out );

	i_dih_xpars( x, compvals );
	if( y_coefficients[14] < 0.0 ) {
		for( j=0; j<12; j+=2 ) {
			jp = j + 1;
			out_partials[j] = -compvals[jp];
			out_partials[jp] = -compvals[j];
		}
	} else {
		for( j=0; j<12; j++ ) {
			out_partials[j] = compvals[j];
		}
	}
	
}


void test_recognize( void )
{
	double x, bound[2];
	
	while( 1 ) {
		printf("Enter a number:  ");
		scanf("%lf", &x);
		printf("Got %.18f\n", x);
		i_recognize( x, bound );
		printf("bound = [%.18f, %.18f]\n", bound[0], bound[1]);
		if( bound[1] < bound[0] ) {
			printf("Interval is backwards.\n");
		}
		printf("\n\n");
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

2.0		2.0
2.0		2.0
2.0		2.0
2.828427124746189847    2.828427124746190291
2.0		2.0
2.0		2.0

2.041				2.041
2.0291			2.0291
2.03412			2.03412
2.6154			2.6154
2.52331			2.52331
2.5123			2.5123

2.0132			2.0132
2.0291			2.0291
2.03412			2.03412
2.6154			2.6154
2.42119			2.42119
2.1023			2.1023

2.51   2.51
2.067323847522493541    2.067323847641237222
2.255    2.255
2.51    2.51
2.255    2.255
2.255    2.255

2 2.01
2 2.01
2 2.01
2 2.01
2 2.01
2 2.01

2 2.005
2 2.005
2 2.005
2 2.005
2 2.005
2 2.005

2 2.001
2 2.001
2 2.001
2 2.001
2 2.001
2 2.001

2 2 
2.1275 2.1275 
2 2 
4.235625 4.235625 
2.255 2.255 
2.51 2.51

2 2
2 2
2 2
2.2265332 2.2265333
3.332299871 3.332299872
3.99996 4

Cell failed:
data        = [-1.25901535020628152     1.33528567645794372]
delta       = [-3.73618149751564488e-07 1.39499661599984393e-06]
vor0_val    = [-1.25901535020628152     1.03171540208955981]
dih_val     = [0        1.20035695677494569]
y = 
2.0    2.0
2.0    2.0
2.0    2.0
2.226533219516277029    2.226533223316073062
3.322998711809516337    3.322998715847731432
3.999999996334314467    4.000000000372528675

*/


