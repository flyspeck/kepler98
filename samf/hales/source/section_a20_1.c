/* octahedra.c, by Samuel Ferguson, (c) 1998. */

#include "system_headers.h"
#include "interval.h"
#include "i_sphere.h"
#include "i_bounds.h"
#include "i_appendix.h"
#include "i_voronoi.h"
#include "macros.h"


#define MAXDEPTH	32			/* Max recursion depth (was 10) */

#define SUBDIV		2			/* (4) 	number of subdivisions (1-d) */
#define SUBFRAC	0.5		/* (0.25) 	1/SUBDIV (make this exact . . .) */

#define COUNTLIMIT	1.0e5		/* 4.0e7 */

#define SEAN					0
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
double sean_val2[2];
double three2[2];
double three22[2];
double relation_const[2];
double slop;
int do_max;

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
void test_recognize( void );
void set_rhs( int index, int k0, int k1, int k2,
	double out[2] );

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
	double tet[14], tetdiff[6], max;
	double diff, dt, origslop;
	time_t tstart, tstop;
	clock_t cstart, cstop;
	char *charptr;
	struct Cell_Info cell_info;
	int onefailed, totalfailed;
	int a2ind, a4ind, b3ind, b4ind;
	int a2ran[2], a4ran[2], b3ran[2], b4ran[2];
	double k0, k1, k2;
	double intercellcount;
	int intermaxdepth;

#if SEAN
	dt = 2.51682;
	i_recognize( dt, sean_val );
#else
	sean_val[0] = TWO51_LO;
	sean_val[1] = TWO51_HI;
#endif
	i_mult( sean_val, sean_val, sean_val2 );

	dt = 3.2;
	i_recognize( dt, three2 );
	i_mult( three2, three2, three22 );
	
	onefailed = 0;
	
	ROUND_NEAR;
	if( sean_val[1] <= sean_val[0] )
		printf("sean_val = [%.20f, %.20f]\n", sean_val[0], sean_val[1]);
	if( three2[1] <= three2[0] )
		printf("3.2 = [%.20f, %.20f]\n", three2[0], three2[1]);

	printf("Welcome to Section A_20.1.  \n");
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
	totalfailed = 0;
	
	for( n=1; n < 7; n++ ) { 
		
		printf("Starting case %d:\n", n);

		if( n % 2 ) {
			printf("\tChecking  vor0.\n");
			sean_init();
		} else {
			printf("\tChecking -tau0.\n");
			sean_init();
			ROUND_DOWN;
			sol_coeff[0] = phi0_val[0] - zetapt_val[1];
			ROUND_UP;
			sol_coeff[1] = phi0_val[1] - zetapt_val[0];
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
		onefailed = 0;
		cell_info.depth = 0;
		intercellcount = 0.0;
		intermaxdepth = 0;
		
		b3ran[0] = 0;
		b3ran[1] = 2;
		b4ran[0] = 0;
		b4ran[1] = 1;
		a2ran[0] = 0;
		a2ran[1] = 1;
		a4ran[0] = 0;
		a4ran[1] = 1;
		
		switch( n ) {
			case 1:
				break;
			case 2:
				break;
			case 3:
				a4ran[1] = 0;
				b4ran[0] = 1;
				b3ran[0] = 1;
				break;
			case 4:
				a4ran[1] = 0;
				b4ran[0] = 1;
				b3ran[0] = 1;
				break;
			case 5:
				a4ran[1] = 0;
				b4ran[0] = 1;
				b3ran[1] = 0;
				break;
			case 6:
				a4ran[1] = 0;
				b4ran[0] = 1;
				b3ran[1] = 0;
				break;
			default:
				printf("Hit default case.\n");
				break;
		}
		for( i=0; i<14; i++ )
			tet[i] = 2.0;
		tet[6] = TWOSQRT2_LO;
		tet[7] = three2[1];
		
		/* Run through all of the cases now. */
		for( a2ind=a2ran[0]; a2ind<=a2ran[1]; a2ind++ ) {
			switch( a2ind ) {
				case 0:
					tet[12] = 2.0;
					tet[13] = 2.0;
					break;
				case 1:
					tet[12] = sean_val[0];
					tet[13] = sean_val[1];
					break;
				default:
					printf("Hit default in a2ind switch.\n");
					break;
			}
			for( a4ind=a4ran[0]; a4ind<=a4ran[1]; a4ind++ ) {
				switch( a4ind ) {
					case 0:
						tet[0] = 2.0;
						tet[1] = 2.0;
						break;
					case 1:
						tet[0] = sean_val[0];
						tet[1] = sean_val[1];
						break;
					default:
						printf("Hit default in a4ind switch.\n");
						break;
				}
				for( b3ind=b3ran[0]; b3ind<=b3ran[1]; b3ind++ ) {
					switch( b3ind ) {
						case 0:
							tet[8] = 2.0;
							tet[9] = 2.0;
							break;
						case 1:
							tet[8] = sean_val[0];
							tet[9] = sean_val[1];
							break;
						case 2:
							tet[8] = TWOSQRT2_LO;
							tet[9] = TWOSQRT2_HI;
							break;
						default:
							printf("Hit default in b3ind switch.\n");
							break;
					}
					i = b3ind - 1;
					if( i < b4ran[0] )
						i = b4ran[0];
					for( b4ind=i; b4ind<=b4ran[1]; b4ind++ ) {
						switch( b4ind ) {
							case 0:
								tet[10] = sean_val[0];
								tet[11] = sean_val[1];
								break;
							case 1:
								tet[10] = TWOSQRT2_LO;
								tet[11] = TWOSQRT2_HI;
								break;
							default:
								printf("Hit default in b4ind switch.\n");
								break;
						}
						
						/* Now compute k0, k1, k2 */
						k0 = 0;
						k1 = 0;
						k2 = 0;
						
						if( b3ind==0 )
							k0++;
						if( b3ind==1 )
							k1++;
						if( b3ind==2 )
							k2++;
						
						if( b4ind==0 )
							k1++;
						if( b4ind==1 )
							k2++;
						
						/* Compute right-hand-side constants */
						set_rhs( n, k0, k1, k2, tetdiff );

						ROUND_NEAR;
						printf("\nrhs = [%.18f, %.18f]\n",
							tetdiff[0], tetdiff[1]);

						ROUND_DOWN;
						slop = origslop + tetdiff[0];
						
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
						for( i=0; i<14; i+=2 )
							printf("%0.18f\t%0.18f\n", tet[i], tet[i+1]);

						maxdepth = 0;
						failed = 0;
						cell_info.depth = 0;
						cellcount = 0.0;

						if( !( b3ind==0 && b4ind==0 ) ) {
							verify_cell( tet, cell_info );
							ROUND_NEAR;
							if( failed ) {
								onefailed++;
							}
						} else {
							printf("\t(Case excluded.)\n");
						}
						if( maxdepth > intermaxdepth )
							intermaxdepth = maxdepth;
						intercellcount += cellcount;
					}
				}
			}
		}
		ROUND_NEAR;
		printf("\n");
		printf("cellcount    = %16.0f\t\t", intercellcount);
		printf("maxdepth = %d\n", intermaxdepth);
		printf("partialpushes= %16.0f\n", gmacount);
		printf("largecut     = %16.0f\t\t", dualcount);
		printf("discards     = %16.0f\n", vorcount);
		printf("smallcut     = %16.0f\t\t", sccutcount);
		printf("verifycount  = %16.0f\n", verifycount);
		printf("partialcount = %16.0f\t\t", partialcount);
		diff = intercellcount - verifycount - partialcount - vorcount;
		printf("diff         = %16.0f\n", diff);
		if( !onefailed )
			printf("Verifications succeeded.\n\n");
		else {
			printf("(%d) verifications FAILED.\n\n", onefailed);
			totalfailed += onefailed;
		}
		celltotal += intercellcount;
	} /* end case loop */
		
	ROUND_NEAR;
	
	if( !totalfailed )
		printf("All verifications succeeded.\n\n");
	else {
		if( totalfailed == 1 )
			printf("One verification FAILED.\n\n");
		else
			printf("(%d) verifications FAILED.\n\n", totalfailed);
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


void verify_cell( double y[14], struct Cell_Info info )
{
	int i, j, bounds[6], n0, n1, n2, n3, n4, n5;
	double diff[6], t[14], val;
	double x[12], xp[12], yp[12], ynew[14];
	double data[2], data2[2], vor0_val[2], vor0p_val[2];
	double delta[2], sqrtdelta[2], delta_part[12];
#if PARTIALS
	double partials[12], rel_part[12];
#endif
	
	/* If part of the verification has already failed, bail. */
	if( failed )
		return;
	
	cellcount += 1.0;
	
	/* Compute square of edge lengths */
	
	for( i=0; i<2; i++ ) {
		yp[i] = y[i + 12];
	}	

	ROUND_UP;
	for( i=1; i<12; i+=2 ) {
		x[i] = y[i]*y[i];
	}
	xp[1] = yp[1]*yp[1];
	ROUND_DOWN;
	for( i=0; i<12; i+=2 ) {
		x[i] = y[i]*y[i];
	}
	xp[0] = yp[0]*yp[0];
	
	for( i=2; i<8; i++ ) {
		yp[i] = y[i];
		xp[i] = x[i];
	}
	for( i=8; i<12; i++ ) {
		yp[i] = 2.0;
		xp[i] = 4.0;
	}
	

	/* First try to verify over the current cell.  */
	cross_diag2( x, xp, data2 );

	/*
	ROUND_NEAR;
	printf("cross_diag2 = [%0.18g\t%0.18g]\n", 
		data2[0], data2[1]);
	*/
	
	/* Check length of cross-diagonal. */
	if( data2[1] < three22[0] ) {	/* discard */
		vorcount += 1.0;
		return;
	}
		
	/* Score cell properly */
	i_bigdelta_best( xp, delta );
	i_sqrt( delta, sqrtdelta );
	i_vor0( yp, xp, sqrtdelta, vor0p_val );
	
	i_bigdelta_partials_best( x, delta, delta_part );
	i_sqrt( delta, sqrtdelta );
	i_vor0( y, x, sqrtdelta, vor0_val );
		
	ROUND_DOWN;
	data[0] = vor0_val[0] + vor0p_val[0];
	ROUND_UP;
	data[1] = vor0_val[1] + vor0p_val[1];
	
	/*
	ROUND_NEAR;
	printf("vor0_val    = [%0.18g\t%0.18g]\n", 
		vor0_val[0], vor0_val[1]);
	printf("vor0p_val   = [%0.18g\t%0.18g]\n", 
		vor0p_val[0], vor0p_val[1]);
	*/
	/*
	ROUND_NEAR;
	printf("Depth = %d\n", info.depth);
	printf("d1 = [%.18f, %.18f]\n", y[6], y[7]);
	printf("data        = [%0.18g\t%0.18g]\n", 
		data[0], data[1]);
	*/
	
	/* Attempt verification */
	if( data[1] < slop ) {
		verifycount += 1.0;
		return;
	}
	if( data[0] > slop && data2[0] > three22[1] ) {
		/* We are in big trouble. */
		failed = 1;
		printf("Cell failed:\n");
		printf("y = \n");
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", y[i], y[i+1]);
		printf("yp = \n");
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", yp[i], yp[i+1]);
		printf("data        = [%0.18g\t%0.18g]\n", 
			data[0], data[1]);
		printf("cross_diag2 = [%0.18g\t%0.18g]\n", 
			data2[0], data2[1]);
		printf("vor0_val    = [%0.18g\t%0.18g]\n", 
			vor0_val[0], vor0_val[1]);
		printf("vor0p_val   = [%0.18g\t%0.18g]\n", 
			vor0p_val[0], vor0p_val[1]);
		printf("delta       = [%0.18g\t%0.18g]\n", 
			delta[0], delta[1]);
		return;
	}
	
	/* if verification fails, try to bump off, or at least
		reduce dimension */

	for( i=0; i<14; i++ )
		ynew[i] = y[i];	/* copy y */
	
#if PARTIALS
	if( data2[0] > three22[1] && info.depth > PARTIALSTART )	{
		
		/* don't include yp partials for y1, y5, y6 */
		i_vor0_partials( y, x, t );
		i_vor0_partials( yp, xp, partials );
		for( i=0; i<2; i++ ) {
			rel_part[i] = t[i];
		}
		ROUND_DOWN;
		for( i=2; i<8; i+=2 ) {
			rel_part[i] = t[i] + partials[i];
		}
		ROUND_UP;
		for( i=3; i<8; i+=2 ) {
			rel_part[i] = t[i] + partials[i];
		}
		for( i=8; i<12; i++ ) {
			rel_part[i] = t[i];
		}
		
		ROUND_NEAR;
		j = 0;
		for( i=0; i<6; i++ ) {	/* differentials */
			diff[i] = ynew[j+1] - ynew[j];
			j += 2;
		}
		
		for( j=0; j<6; j++ ) {
			if( diff[j] > 1.0e-14 ) {	/* ignore tight spots */
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
		
		t[12] = ynew[12];
		t[13] = ynew[13];
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
		printf("cross_diag2 = [%0.18g\t%0.18g]\n", 
			data2[0], data2[1]);
		printf("vor0_val    = [%0.18g\t%0.18g]\n", 
			vor0_val[0], vor0_val[1]);
		printf("vor0p_val   = [%0.18g\t%0.18g]\n", 
			vor0p_val[0], vor0p_val[1]);
		printf("delta       = [%0.18g\t%0.18g]\n", 
			delta[0], delta[1]);
		printf("y = \n");
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", y[i], y[i+1]);
#if PARTIALS
	if( data2[0] > three22[1] && info.depth > PARTIALSTART )	{
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
		ROUND_NEAR;
		printf("cross_diag2 = [%.18f, %.18f]\n", 
			temp2[0], temp2[1]);
		i_sqrt( temp2, temp );
		ROUND_NEAR;
		printf("cross_diag  = [%.18f, %.18f]\n", 
			temp[0], temp[1]);

		i_bigdelta_partials_best( x, delta, delta_part );
		i_sqrt( delta, sqrtdelta );
		i_vor0( y, x, sqrtdelta, temp );
		printf("vor0_val  = [%.18f, %.18f]\n", 
			temp[0], temp[1]);
	}
}


void sean_init( void )
{
	int i;
	double temp[2], temp2[2], tet[12];
			
	/* initialize local constants */
	t0_val[0] = 0.5*sean_val[0];
	t0_val[1] = 0.5*sean_val[1];
	ROUND_DOWN;
	oot0_val[0] = 1.0/t0_val[1];
	ROUND_UP;
	oot0_val[1] = 1.0/t0_val[0];
	phi_fun( t0_val, t0_val, phi0_val );
	/* compute zetapt_val here */
	tet[0] = 2.0;
	tet[1] = 2.0;
	ROUND_DOWN;
	temp[0] = sqrt(tet[0])/5.0;
	ROUND_UP;
	temp[1] = sqrt(tet[1])/5.0;
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
	double relations[4] = 
		{-0.075, -0.235, -0.137, -0.3109};
	
	if( num > 4 ) {
		printf("num > 4.\n");
		return;
	}
	i_recognize( relations[num], relation_const );
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


void set_rhs( int index, int k0, int k1, int k2,
	double out[2] )
{
	int i;
	double vals[5] = 
		{0.009, 0.008, 0.04683, 0.0066, 0.0461};
	double fatvals[10];
	double p1[2], p2[2], p3[2], p4[2], p5[2];
	
	for( i=0; i<5; i++ ) {
		i_recognize( vals[i], fatvals + 2*i );
	}
	
	switch( index ) {
		case 1:
			tomZfun( 4, k1 + k2, p1 );
			p2[0] = (double) -k2;
			p2[1] = (double) -k2;
			i_mult( p2, fatvals, p3 );
			p2[0] = (double) -(k0 + 2*k2);
			p2[1] = (double) -(k0 + 2*k2);
			i_mult( p2, fatvals + 2, p4 );
			ROUND_DOWN;
			out[0] = p1[0] + p3[0] + p4[0]/3.0;
			ROUND_UP;
			out[1] = p1[1] + p3[1] + p4[1]/3.0;
			break;
		case 2:
			tomDfun( 4, k1 + k2, p1 );
			p2[0] = (double) -(k0 + 2*k2 - 3);
			p2[1] = (double) -(k0 + 2*k2 - 3);
			i_mult( p2, fatvals + 2, p3 );
			p4[0] = (double) -k2;
			p4[1] = (double) -k2;
			i_mult( p4, fatvals + 6, p5 );
			ROUND_DOWN;
			out[0] = -p1[1] - fatvals[5] + p3[0]/3.0 + p5[0];
			ROUND_UP;
			out[1] = -p1[0] - fatvals[4] + p3[1]/3.0 + p5[1];
			break;
		case 3:
			ROUND_NEAR;
			p1[0] = 0.0461 + 0.009 + 2.0*0.008;
			i_recognize( p1[0], p2 );
			tomZfun( 4, 2, p1 );
			i_sub( p1, p2, out );
			break;
		case 4:
			ROUND_NEAR;
			p1[0] = 0.04683 + 0.008 + 2.0*0.0066;
			i_recognize( p1[0], p2 );
			tomDfun( 5, 1, p1 );
			i_add( p1, p2, p3 );
			out[0] = -p3[1];
			out[1] = -p3[0];
			break;
		case 5:
			ROUND_NEAR;
			p1[0] = 0.0461 + 0.008;
			i_recognize( p1[0], p2 );
			sconst_fun( 5, p1 );
			i_sub( p1, p2, out );
			break;
		case 6:
			tconst_fun( 5, p1 );
			i_add( p1, fatvals + 2, p2 );
			out[0] = -p2[1];
			out[1] = -p2[0];
			break;
		default:
			printf("Hit default in set_rhs.\n");
			break;
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

*/


