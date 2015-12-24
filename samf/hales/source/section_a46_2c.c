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
#define PARTIALS			1
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
double relation_const[4];
double slop;
int do_max;

struct Cell_Info {
	int depth; 
	};

/* Prototypes */

void verifyocta( void );
void verify_cell( double tet[12], struct Cell_Info info );
void testfuns( void );
void sean_init( void );
void set_relation( int num );
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
	int i, n, casestart, casebound;
	double tet[12], tetdiff[6];
	double diff, dt, three2[2], origslop;
	time_t tstart, tstop;
	clock_t cstart, cstop;
	char *charptr;
	struct Cell_Info cell_info;
	int subcase, onefailed;
	int subsub;

	sean_val[0] = TWO51_LO;
	sean_val[1] = TWO51_HI;
	i_mult( sean_val, sean_val, sean_val2 );

	dt = 3.2;
	i_recognize( dt, three2 );
	
	onefailed = 0;
	
	printf("Welcome to Section A.4.6.2.c.  \n");
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
	
	for( n=1; n < 5; n++ ) { 
	
		sean_init();
		
		switch( n ) {
			case 1:
				i_recognize( -0.19, tetdiff );
				break;
			case 2:
				i_recognize( -0.281, tetdiff );
				break;
			case 3:
				i_recognize( -0.11, tetdiff );
				break;
			case 4:
				i_recognize( -0.205, tetdiff );
				break;
			default:
				break;
		}
			
		if( n%2==0 ) {
			/* turn vor0 into -tau0 */
			ROUND_DOWN;
			sol_coeff[0] = phi0_val[0] - zetapt_val[1];
			ROUND_UP;
			sol_coeff[1] = phi0_val[1] - zetapt_val[0];
		}
		ROUND_DOWN;
		slop = tetdiff[0] + origslop;
			
		printf("Starting case %d:\n", n);

		ROUND_NEAR;
		printf("\t\trhs = %.18f\n", tetdiff[0]);
		
		casestart = 1;
		casebound = 1;
		/*
		casebound = 2;
		*/
		
		for( subcase=casestart; subcase<=casebound; subcase++ ) {
			
			for( i=0; i<12; i+=2 ) {
				tet[i  ] = 2.0;
				tet[i+1] = sean_val[1];
			}
			
			tet[6] = 2.0;
			tet[7] = 2.0;
			if( n < 3 ) {
				tet[8] = TWOSQRT2_LO;
				tet[9] = TWOSQRT2_HI;
			} else {
				tet[8] = sean_val[0];
				tet[9] = sean_val[1];
			}
			tet[10] = three2[0];
			tet[11] = three2[1];
			
			for( subsub=1; subsub<2; subsub++ ) {				
				
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
			}	/* end subsub */
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
	double x[12], ynew[12], rel_part[12];
	double data[2];
	double delta[2];
	double sqrtdelta[2];
	
	/* If part of the verification has already failed, bail. */
	if( failed )
		return;
	
	cellcount += 1.0;
			
	/* Score cell properly */
	
	ROUND_UP;
	for( i=1; i<12; i+=2 ) {
		x[i] = y[i]*y[i];
	}
	ROUND_DOWN;
	for( i=0; i<12; i+=2 ) {
		x[i] = y[i]*y[i];
	}
	
	i_bigdelta_best( x, delta );
	
	if( delta[1] < 0.0 ) {	/* discard */
		vorcount += 1.0;
		return;
	}
	
	i_sqrt( delta, sqrtdelta );
	i_vor0( y, x, sqrtdelta, data );
	
	
	/* Attempt verification */
	if( data[1] < slop ) {
		verifycount += 1.0;
		return;
	}
	
	if( data[0] > slop ) {
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
		return;
	}
	
	/* if verification fails, try to bump off, or at least
		reduce dimension */

	for( i=0; i<12; i++ )
		ynew[i] = y[i];	/* copy y */
	
#if PARTIALS
	if( info.depth > PARTIALSTART && delta[0] > 0.0 )	{
		i_vor0_partials( y, x, rel_part );
		
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
		printf("y = \n");
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", y[i], y[i+1]);
		printf("data        = [%0.18g\t%0.18g]\n", 
			data[0], data[1]);
		printf("delta       = [%0.18g\t%0.18g]\n", 
			delta[0], delta[1]);
#if PARTIALS
	if( info.depth > PARTIALSTART && delta[0] > 0.0 )	{
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
/*
	double relations[6][4] = 
		{{-358.0,  228501.0,
			 100.0,	 100000.0},
		 {-2715.0, 167382.0,
		   1000.0, 100000.0},
		 {-1517.0,  8285.0,
		 	 1000.0, 10000.0},
		 {-858.0,   390925.0,
		  1000.0,	 1000000.0},
		 {-349.0,   12012.0,
		  1000.0,	 100000.0},
		 {-177.0,		  501.0,
		  1000.0,   10000.0}};
*/
	double relations[6][2] = 
		{{-3.58, -2.28501},
		 {-2.715, -1.67382},
		 {-1.517, -0.8285},
		 {-0.858, -0.390925},
		 {-0.358 + 0.009, -0.12012},
		 {-0.186 + 0.009, -0.0501}};
	
	if( num > 5 ) {
		printf("num > 5.\n");
		return;
	}
	i_recognize( relations[num][0], relation_const );
	i_recognize( relations[num][1], relation_const + 2 );
/*
	ROUND_DOWN;
	relation_const[0] = relations[num][0]/relations[num][2];
	relation_const[2] = (-relations[num][1])/relations[num][3];
	ROUND_UP;
	relation_const[1] = relations[num][0]/relations[num][2];
	relation_const[3] = (-relations[num][1])/relations[num][3];
*/
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

*/


