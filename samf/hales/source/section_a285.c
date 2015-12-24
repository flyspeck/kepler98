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
double slop;
double origslop;
double y_coefficients[6];
int taylor_rel[7];
int do_gma, do_tau, do_bdry, do_quarter;
double taylor_relconst[14];
double one41_2[2];

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
void set_tetrahedron( int num, double out[12] );
void get_tet( int num, double penta[32], double out[12] );
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
	int i, m, n, casestart, casebound;
	double tet[12], tetdiff[6], max;
	double diff, dt;
	time_t tstart, tstop;
	clock_t cstart, cstop;
	char *charptr;
	struct Cell_Info cell_info;
	int subcase, onefailed;
	int subsub, subsubbound;

	sean_val[0] = TWO51_LO;
	sean_val[1] = TWO51_HI;
	
	onefailed = 0;
	
	printf("Welcome to Section A.2.8.5.  \n");
	printf("Enter cell count limit:  ");
	scanf("%lf", &countlimit);

	printf("Enter slop:  ");
	scanf("%lf", &slop);
	printf("\t\tSlop = %g\n", slop );
	origslop = slop;
	
	i_init();	/* initialize interval constants */
	sean_init();
	
	dt = 1.41*1.41;
	i_recognize( dt, one41_2 );
	
	ROUND_NEAR;

	time(&tstart);
	cstart = clock();
	
	charptr = ctime( &tstart );
	printf("Starting time:  \t");
	puts( charptr );
	
	celltotal = 0.0;
	
	for( n=1; n < 41; n++ ) { 	

		printf("Case %d:\n", n);
		
		if( n%5==0 ) {
			do_quarter = 1;
			do_gma = 0;
		} else {
			do_quarter = 0;
			do_gma = 1;
		}
		
		do_tau = 1;
		
		casestart = 1;
		if( do_tau ) {
			casebound = 1;
			/* casebound = 2; */
		} else {
			casebound = 1;
		}
		
		for( subcase=casestart; subcase<=casebound; subcase++ ) {
			
			if( do_gma ) {
					printf("\tGma scoring:\n");
			} else {
					printf("\tVor scoring:\n");
			}
			
			set_relation( n );
			set_tetrahedron( n, tet );
			
			ROUND_NEAR;
			
			if( do_tau && !do_quarter ) {
				subsubbound = 2;
			} else {
				subsubbound = 1;
			}
			
			for( subsub=0; subsub<subsubbound; subsub++ ) {
			
				/* Use dimension-reduction?  Do we need it? */
				if( subsub==1 && !do_quarter ) {
					printf("\t\tBoundary:\n");
					do_bdry = 1;
				} else {
					do_bdry = 0;
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
	double x[12], ynew[12];
	double data[2], crad2_val[2];
	double rel_part[12];
	
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
	if( x[7] > 8.0 )
		x[7] = 8.0;

	/* First try to verify over the current cell.  */
	if( !do_quarter ) {
		i_crad2_best( x, crad2_val );
		if( do_bdry ) {
			if( crad2_val[0] > one41_2[1] ) {
				dualcount += 1.0;
				return;
			}
			if( crad2_val[1] < one41_2[0] ) {
				sccutcount += 1.0;
				return;
			}
		} else {
			if( do_gma ) {
				if( crad2_val[0] > one41_2[1] ) {
					dualcount += 1.0;
					return;
				}
				if( crad2_val[1] < one41_2[0] ) {
				}
			} else {
				if( crad2_val[1] < one41_2[0] ) {
					sccutcount += 1.0;
					return;
				}
				if( crad2_val[0] > one41_2[1] ) {
				}
			}
		}
	}
		
	/* Score cell properly */
	octa_score( y, x, data, rel_part );
	
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
		if( !do_quarter ) {
			printf("crad2 = [%0.18g\t%0.18g]\n", 
				crad2_val[0], crad2_val[1]);
		}
		return;
	}
	
	/* if verification fails, try to bump off, or at least
		reduce dimension */

	for( i=0; i<12; i++ )
		ynew[i] = y[i];	/* copy y */
	
#if PARTIALS
	if( !do_bdry )	{
	
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
		if( !do_quarter ) {
			printf("crad2 = [%0.18g\t%0.18g]\n", 
				crad2_val[0], crad2_val[1]);
		}
		printf("y = \n");
		for( i=0; i<12; i+=2 )
			printf("%0.18f\t%0.18f\n", y[i], y[i+1]);
#if PARTIALS
	if( !do_bdry )	{
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
		
		do_tau = 1;
		do_gma = 0;
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

	/* turn vor0 into -tau0 */
	ROUND_DOWN;
	sol_coeff[0] = phi0_val[0] - zetapt_val[1];
	ROUND_UP;
	sol_coeff[1] = phi0_val[1] - zetapt_val[0];
}


void set_relation( int num )
{
/*
cutval = 2.17667782820508515

2.8.5 (0,0,0) (0,0)

a23lo = 2.0
a23hi = cutval
a33lo = 2.0
a33hi = cutval
a12lo = 2.0
a12hi = 2.255
a43lo = 2.0
a43hi = 2.255

pracrange = {{2,2},{a12lo,a12hi},{2,cutval},{2,2},{2,2.2549},{2,2.51},
		{a23lo,a23hi},{2,2},{2,2.2549},
		{a33lo,a33hi},{2,2},{2,2.2549},
		{a43lo,a43hi},{2,2},{2,2.51},
		{2.7,2Sqrt[2]}}

{aval[1,2],aval[2,2],aval[3,2],aval[4,2],aval[5,2],aval[1,3],aval[2,3],
  aval[3,3],aval[4,3],aval[5,3],aval[1,5],aval[2,5],aval[3,5],aval[4,5],
  aval[5,5],aval[1,6],aval[2,6],aval[3,6],aval[4,6],aval[5,6],cval[1],cval[2],
  cval[3],cval[4],cval[5],bval}
  
List(0.0005824342923155256,0.01299678288138017,
   -0.01299678288138017,-0.01299678288138017,
   -0.0005824342923155256,-0.01299678288138017,
   0.01299678288138017,0.01299678288138017,
   0.0005824342923155256,-0.0005824342923155256,
   -0.0007597211794387881,0.02118533283902324,
   -0.02118533283902324,0.0007597211794387881,
   -0.0007597211794387881,0.0007597211794387881,
   0.0007597211794387881,-0.02118533283902324,
   0.02118533283902324,-0.0007597211794387881,
   0.5839166747660211,0.4632107380254502,
   0.6438293089439817,0.5400265667290888,
   0.867751414408088,0.4541887975443426)

room = 0.01468397099075227

-1.5 ptval - 0.06585- 0.022652 = -0.17156246850269623


2.8.5 (0,0,0) (0,1)

a23lo = 2.0
a23hi = cutval
a33lo = 2.0
a33hi = cutval
a12lo = 2.0
a12hi = 2.255
a43lo = 2.255
a43hi = 2.51

List(-0.04824643252537585,0.02329205714383042,
   -0.02329205714383042,0.1674983635677471,
   -0.1569377150264126,-0.02329205714383042,
   0.02329205714383042,-0.1674983635677471,
   0.1569377150264126,0.04824643252537585,
   0.01244415149606736,0.03235492928548255,
   -0.02186974787616857,-0.03702043276823618,
   -0.01244415149606736,0.01244415149606736,
   -0.01244415149606736,-0.03235492928548255,
   0.02186974787616857,0.03702043276823618,
   0.6014950905071683,0.3752049329988726,
   0.998224912899468,-0.2067879935185711,1.,
   0.4128444122478613)

room = 0.0005193058394593164


2.8.5 (0,0,0) (1,1)

a23lo = 2.0
a23hi = cutval
a33lo = 2.0
a33hi = cutval
a12lo = 2.255
a12hi = 2.51
a43lo = 2.255
a43hi = 2.51

List(0.02874944926886735,-0.1355829956515062,
   0.1070716351971938,0.1355829956515062,
   -0.02874944926886735,0.1355829956515062,
   -0.1070716351971938,-0.1355829956515062,
   0.02874944926886735,-0.02874944926886735,
   0.0870749172165697,-0.1232494555948707,
   0.1232494555948707,-0.1573359633463219,
   0.1573359633463219,-0.1573359633463219,
   -0.0870749172165697,0.1232494555948707,
   -0.1232494555948707,0.1573359633463219,
   -0.07664383380772021,1.,-0.3419331087911444,
   0.3440049118151675,-0.2412038690284744,
   0.07639731363595851)

room = 0.006528630627937043


2.8.5 (0,0,1)

a23lo = 2.0
a23hi = cutval
a33lo = cutval
a33hi = cutvalhi

List(0.01693860366506827,0.06650513085303666,
   0.1011072844100986,0.00933419519520839,
   -0.01693860366506827,-0.06650513085303666,
   -0.1011072844100986,-0.00933419519520839,
   0.01693860366506827,-0.01693860366506827,
   0.06148802693356714,0.07923785656245052,
   -0.1114108141186571,-0.0736489600671163,
   -0.06148802693356714,0.06148802693356714,
   -0.06148802693356714,-0.07923785656245052,
   0.1114108141186571,0.0736489600671163,
   0.1477578272601813,0.3282815284748696,
   0.5057271709696893,0.1703858231145003,
   0.5276230087005028,0.2393067362487741)

room = 0.000920864181930337

2.8.5 (0,1,0)

a23lo = cutval
a23hi = cutvalhi
a33lo = 2.0
a33hi = cutval

List(0.0169386036650676,0.13632725791818,
   -0.0093341951952074,0.06650513085303511,
   -0.0169386036650676,-0.13632725791818,
   0.0093341951952074,-0.06650513085303511,
   0.0169386036650676,-0.0169386036650676,
   0.06148802693356603,0.1109119698021565,
   -0.07923785656245163,0.06148802693356603,
   -0.06148802693356603,0.06148802693356603,
   -0.06148802693356603,-0.1109119698021565,
   0.07923785656245163,-0.06148802693356603,
   0.2874020813904749,-0.0917156630258051,
   0.839954312862954,-0.1537623554097296,
   0.7978969827018596,0.2393067362487742)

room = 0.000920864181930353

2.8.5 (0,1,1)

a23lo = cutval
a23hi = cutvalhi
a33lo =cutval
a33hi = cutvalhi

List(0.01693860366506894,-0.06650513085303555,
   0.02467148784967854,-0.00933419519520684,
   -0.02186989561511065,0.06650513085303555,
   -0.02467148784967854,0.00933419519520684,
   0.02186989561511065,-0.01693860366506894,
   0.06148802693356736,-0.110911969802157,
   0.110911969802157,-0.07740634150018744,
   0.0872689254002688,-0.0872689254002688,
   -0.06148802693356736,0.110911969802157,
   -0.110911969802157,0.07740634150018744,
   0.1792512085157029,0.83783899847087,
   -0.1934770173979916,0.6595412298021359,
   0.2324569250667658,0.2393067362487723)

room = 0.00808806136947969

2.8.5 (1,0,1)

List(0.01020244499071631,0.,-0.1060310277065948,0.,
   -0.01020244499071631,0.,0.1060310277065948,0.,
   0.01020244499071631,-0.01020244499071631,
   0.085976802039418,0.1019881007004708,
   -0.1031076114857755,-0.101471144262149,
   0.101471144262149,-0.101471144262149,
   -0.085976802039418,-0.1019881007004708,
   0.1031076114857755,0.101471144262149,
   0.2800313166479164,0.02536286944858234,
   0.89170100196956,0.2457696977552035,
   0.06357786944753662,0.2082364359685608)

room = 0.005298434373798893

2.8.5 (1,1,1)

List(0.02849679962536466,0.01836850892404773,
   0.01836850892404773,-0.01836850892404773,
   -0.03548719417712775,-0.01836850892404773,
   -0.01836850892404773,0.01836850892404773,
   0.03548719417712775,-0.02849679962536466,
   0.0711301335230865,0.102771474844137,
   0.102771474844137,-0.06337485423229527,
   0.07735564333582112,-0.07735564333582112,
   -0.0711301335230865,-0.102771474844137,
   -0.102771474844137,0.06337485423229527,
   0.3374611201122192,0.2957601381019636,
   0.2790781685197414,0.6433219695360872,
   0.3844509732745327,0.2693788350757662)

room = 0.01519055248572797
*/
	int i, j, k;
	double all_variables[8][26] = 
		{{0.0005824342923155256,0.01299678288138017,
   -0.01299678288138017,-0.01299678288138017,
   -0.0005824342923155256,-0.01299678288138017,
   0.01299678288138017,0.01299678288138017,
   0.0005824342923155256,-0.0005824342923155256,
   -0.0007597211794387881,0.02118533283902324,
   -0.02118533283902324,0.0007597211794387881,
   -0.0007597211794387881,0.0007597211794387881,
   0.0007597211794387881,-0.02118533283902324,
   0.02118533283902324,-0.0007597211794387881,
   0.5839166747660211,0.4632107380254502,
   0.6438293089439817,0.5400265667290888,
   0.867751414408088,0.4541887975443426},
   
   {-0.04824643252537585,0.02329205714383042,
   -0.02329205714383042,0.1674983635677471,
   -0.1569377150264126,-0.02329205714383042,
   0.02329205714383042,-0.1674983635677471,
   0.1569377150264126,0.04824643252537585,
   0.01244415149606736,0.03235492928548255,
   -0.02186974787616857,-0.03702043276823618,
   -0.01244415149606736,0.01244415149606736,
   -0.01244415149606736,-0.03235492928548255,
   0.02186974787616857,0.03702043276823618,
   0.6014950905071683,0.3752049329988726,
   0.998224912899468,-0.2067879935185711,1.,
   0.4128444122478613},
   
   {0.02874944926886735,-0.1355829956515062,
   0.1070716351971938,0.1355829956515062,
   -0.02874944926886735,0.1355829956515062,
   -0.1070716351971938,-0.1355829956515062,
   0.02874944926886735,-0.02874944926886735,
   0.0870749172165697,-0.1232494555948707,
   0.1232494555948707,-0.1573359633463219,
   0.1573359633463219,-0.1573359633463219,
   -0.0870749172165697,0.1232494555948707,
   -0.1232494555948707,0.1573359633463219,
   -0.07664383380772021,1.,-0.3419331087911444,
   0.3440049118151675,-0.2412038690284744,
   0.07639731363595851},
   
   {0.01693860366506827,0.06650513085303666,
   0.1011072844100986,0.00933419519520839,
   -0.01693860366506827,-0.06650513085303666,
   -0.1011072844100986,-0.00933419519520839,
   0.01693860366506827,-0.01693860366506827,
   0.06148802693356714,0.07923785656245052,
   -0.1114108141186571,-0.0736489600671163,
   -0.06148802693356714,0.06148802693356714,
   -0.06148802693356714,-0.07923785656245052,
   0.1114108141186571,0.0736489600671163,
   0.1477578272601813,0.3282815284748696,
   0.5057271709696893,0.1703858231145003,
   0.5276230087005028,0.2393067362487741},
   
   {0.0169386036650676,0.13632725791818,
   -0.0093341951952074,0.06650513085303511,
   -0.0169386036650676,-0.13632725791818,
   0.0093341951952074,-0.06650513085303511,
   0.0169386036650676,-0.0169386036650676,
   0.06148802693356603,0.1109119698021565,
   -0.07923785656245163,0.06148802693356603,
   -0.06148802693356603,0.06148802693356603,
   -0.06148802693356603,-0.1109119698021565,
   0.07923785656245163,-0.06148802693356603,
   0.2874020813904749,-0.0917156630258051,
   0.839954312862954,-0.1537623554097296,
   0.7978969827018596,0.2393067362487742},
   
   {0.01693860366506894,-0.06650513085303555,
   0.02467148784967854,-0.00933419519520684,
   -0.02186989561511065,0.06650513085303555,
   -0.02467148784967854,0.00933419519520684,
   0.02186989561511065,-0.01693860366506894,
   0.06148802693356736,-0.110911969802157,
   0.110911969802157,-0.07740634150018744,
   0.0872689254002688,-0.0872689254002688,
   -0.06148802693356736,0.110911969802157,
   -0.110911969802157,0.07740634150018744,
   0.1792512085157029,0.83783899847087,
   -0.1934770173979916,0.6595412298021359,
   0.2324569250667658,0.2393067362487723},
   
   {0.01020244499071631,0.,-0.1060310277065948,0.,
   -0.01020244499071631,0.,0.1060310277065948,0.,
   0.01020244499071631,-0.01020244499071631,
   0.085976802039418,0.1019881007004708,
   -0.1031076114857755,-0.101471144262149,
   0.101471144262149,-0.101471144262149,
   -0.085976802039418,-0.1019881007004708,
   0.1031076114857755,0.101471144262149,
   0.2800313166479164,0.02536286944858234,
   0.89170100196956,0.2457696977552035,
   0.06357786944753662,0.2082364359685608},
   
   {0.02849679962536466,0.01836850892404773,
   0.01836850892404773,-0.01836850892404773,
   -0.03548719417712775,-0.01836850892404773,
   -0.01836850892404773,0.01836850892404773,
   0.03548719417712775,-0.02849679962536466,
   0.0711301335230865,0.102771474844137,
   0.102771474844137,-0.06337485423229527,
   0.07735564333582112,-0.07735564333582112,
   -0.0711301335230865,-0.102771474844137,
   -0.102771474844137,0.06337485423229527,
   0.3374611201122192,0.2957601381019636,
   0.2790781685197414,0.6433219695360872,
   0.3844509732745327,0.2693788350757662}};
	double all_room[8] = 
		{0.01468397099075227,0.0005193058394593164,
		0.006528630627937043,0.000920864181930337,
		0.000920864181930353,0.00808806136947969,
		0.005298434373798893,0.01519055248572797};
	double bval, room, avals[6];

	sean_init();
	
	if( num > 40 ) {
		printf("num > 40.\n");
		return;
	}
	if( num < 1 ) {
		printf("num < 1.\n");
		return;
	}
	
	k = (num - 1)/5;

	bval = all_variables[k][25];
	room = all_room[k];
	
	i = (num - 1)%5;
	/* y1 coeff */
	avals[0] = 0.0;
	/* y2 coeff */
	avals[1] = all_variables[k][0*5+i];
	/* y3 coeff */
	avals[2] = all_variables[k][1*5+i];
	/* y4 coeff */
	avals[3] = 0.0;
	/* y5 coeff */
	avals[4] = all_variables[k][2*5+i];
	/* y6 coeff */
	avals[5] = all_variables[k][3*5+i];

	for( j=0; j<6; j++ ) {
		y_coefficients[j] = avals[j];
	}

	ROUND_DOWN;
	slop = origslop - all_variables[k][i+20] + room;
	
	ROUND_NEAR;
	printf("rhs = %.18f\n", room - all_variables[k][i+20]);
	
	/* Order of relation constants:
		sol, gma, vor, octavor, dih1, dih2, dih3	
		 0		1		 2			3		 		4		  5		  6	*/
	for( j=0; j<7; j++ ) {
		taylor_rel[j] = 0;
		taylor_relconst[2*j] = 0.0;
		taylor_relconst[2*j+1] = 0.0;
	}
	taylor_relconst[8] = -bval;
	taylor_relconst[9] = -bval;
	if( do_tau ) {
		if( do_gma ) {
			taylor_relconst[2] = 1.0;
			taylor_relconst[3] = 1.0;
		} else {
			taylor_relconst[4] = 1.0;
			taylor_relconst[5] = 1.0;
		}
		/* This gives us -tau */
		taylor_relconst[0] = -zetapt_val[1];
		taylor_relconst[1] = -zetapt_val[0];
	}
}


void set_tetrahedron( int num, double out[12] )
{
/*
cutval = 2.17667782820508515

2.8.5 (0,0,0) (0,0)

a23lo = 2.0
a23hi = cutval
a33lo = 2.0
a33hi = cutval
a12lo = 2.0
a12hi = 2.255
a43lo = 2.0
a43hi = 2.255

pracrange = {{2,2},{a12lo,a12hi},{2,cutval},
		{2,2},{2,2.2549},{2,2.51},
		{a23lo,a23hi},{2,2},{2,2.2549},
		{a33lo,a33hi},{2,2},{2,2.2549},
		{a43lo,a43hi},{2,2},{2,2.51},
		{2.7,2Sqrt[2]}}

2.8.5 (0,0,0) (0,1)

a23lo = 2.0
a23hi = cutval
a33lo = 2.0
a33hi = cutval
a12lo = 2.0
a12hi = 2.255
a43lo = 2.255
a43hi = 2.51

2.8.5 (0,0,0) (0,1)

a23lo = 2.0
a23hi = cutval
a33lo = 2.0
a33hi = cutval
a12lo = 2.255
a12hi = 2.51
a43lo = 2.255
a43hi = 2.51

2.8.5 (0,0,1)

a23lo = 2.0
a23hi = cutval
a33lo = cutval
a33hi = cutvalhi

2.8.5 (0,1,0)

a23lo = cutval
a23hi = cutvalhi
a33lo = 2.0
a33hi = cutval

2.8.5 (0,1,1)

a23lo = cutval
a23hi = cutvalhi
a33lo =cutval
a33hi = cutvalhi

2.8.5 (1,0,1)

2.8.5 (1,1,1)
*/
	int i, j, k;
	double full_pent[32], temp[2], cutval, cutval2;

	if( num > 40 ) {
		printf("num > 40.\n");
		return;
	}
	if( num < 1 ) {
		printf("num < 1.\n");
		return;
	}

	cutval = 2.17667782820508515;
	cutval2 = 2.255;
	
	k = (num - 1)/5;
	
	for( i=0; i<32; i++ ) {
		full_pent[i] = 0.0;
	}
	
	/* y1 */
	full_pent[0] = 2.0;
	full_pent[1] = 2.245;
	
	/* y4's */
	full_pent[6] = 2.0;
	full_pent[7] = TWO51_HI;
	for( i=7; i<15; i+=3 ) {
		full_pent[2*i] = 2.0;
		full_pent[2*i+1] = TWO51_HI;
	}
	i_recognize( 2.7, temp );
	full_pent[30] = temp[0];
	full_pent[31] = TWOSQRT2_HI;
	
	/* y2 and y3 */
	for( i=1; i<3; i++ ) {
		full_pent[2*i] = 2.0;
		full_pent[2*i+1] = TWO51_HI;
	}
	for( i=6; i<15; i+=3 ) {
		full_pent[2*i] = 2.0;
		full_pent[2*i+1] = TWO51_HI;
	}
	
	/* y5 and y6 */
	for( i=4; i<6; i++ ) {
		full_pent[2*i] = 2.0;
		full_pent[2*i+1] = TWO51_HI;
	}
	for( i=8; i<15; i+=3 ) {
		full_pent[2*i] = 2.0;
		full_pent[2*i+1] = TWO51_HI;
	}
	full_pent[11] = 2.4016;
	full_pent[17] = 2.4016;
	full_pent[23] = 2.4016;
	
	if( k > 2 ) {
		if( k < 6 ) {
			full_pent[5] = cutval;
		} else {
			full_pent[4] = cutval;
		}
		if( k==4 ) {
			full_pent[19] = cutval;
		} else {
			full_pent[18] = cutval;
		}
		if( k==3 || k==6 ) {
			full_pent[13] = cutval;
		} else {
			full_pent[12] = cutval;
		}
	}
	if( k < 3 ) {
		full_pent[5] = cutval;
		full_pent[13] = cutval;
		full_pent[19] = cutval;
		if( k==2 ) {
			full_pent[2] = cutval2;
		} else {
			full_pent[3] = cutval2;
		}
		if( k!=0 ) {
			full_pent[24] = cutval2;
		} else {
			full_pent[25] = cutval2;
		}
	}
	
	i = (num - 1)%5;
	get_tet( i, full_pent, out );
	
	j = 0;
	/*
	ROUND_NEAR;
	for( i=0; i<10; i+=2 ) {
		for( j=0; j<6; j+=2 ) {
			printf("[%.4f, %.4f]\t", 
				full_pent[3*i+j], full_pent[3*i+j+1]);
		}
		printf("\n");
	}
	printf("[%.4f, %.4f]\n\n", 
		full_pent[30], full_pent[31]);
	*/
}


void get_tet( int num, double penta[32], double out[12] )
{
	int i, j;
	double tets[5][12];
	
	for( i=0; i<5; i++ ) {
		for( j=0; j<12; j++ ) {
			tets[i][j] = -1.0;
		}
	}
	
	for( i=0; i<5; i++ ) {
		tets[i][0] = penta[0];
		tets[i][1] = penta[1];
	}
	for( i=0; i<6; i++ ) {
		tets[0][2*i] = penta[2*i];
		tets[0][2*i+1] = penta[2*i+1];
	}
	for( j=1; j<4; j++ ) {
		for( i=2; i<5; i++ ) {
			tets[j][2*i] = penta[2*(6 + 3*(j-1) + i-2)];
			tets[j][2*i+1] = penta[2*(6 + 3*(j-1) + i-2)+1];
		}
	}
	tets[4][6] = penta[30];
	tets[4][7] = penta[31];
	for( j=1; j<5; j++ ) {
		tets[j][2*1] = tets[j-1][2*2];
		tets[j][2*1+1] = tets[j-1][2*2+1];
		tets[j][2*5] = tets[j-1][2*4];
		tets[j][2*5+1] = tets[j-1][2*4+1];
	}
	tets[4][4] = tets[0][2];
	tets[4][5] = tets[0][3];
	tets[4][8] = tets[0][10];
	tets[4][9] = tets[0][11];
	
	for( i=0; i<6; i++ ) {
		out[2*i] = tets[num][2*i];
		out[2*i+1] = tets[num][2*i+1];
	}
}


void octa_score( double y[12], double x[12],
	double out[2], double out_partials[12] )
{
	int i, j, jp;
	double t[12], compvals[14], rel_val[2], partials[12];
	double yc, sum;
	
	if( do_gma ) {
		fat_composite( taylor_rel, taylor_relconst, y, x, 
			rel_val, t, compvals );
	} else {
		i_bigdelta_best( x, t );
		i_sqrt( t, t + 2 );
		i_dih_alt( x, t, compvals );
		i_mult( taylor_relconst + 8, compvals, compvals + 2 );
		i_vor0( y, x, t + 2, compvals + 4 );
		i_add( compvals + 2, compvals + 4, rel_val );
		
		i_vor0_partials( y, x, compvals );
		s_dih_partials( x, y, partials );
		ROUND_DOWN;
		for( i=0; i<12; i+=2 ) {
			t[i] = taylor_relconst[8]*partials[i+1] + compvals[i];
		}
		ROUND_UP;
		for( i=1; i<12; i+=2 ) {
			t[i] = taylor_relconst[8]*partials[i-1] + compvals[i];
		}
	}

	/*
	ROUND_NEAR;
	printf("taylor partials = \n");
	for( i=0; i<12; i+=2 )
		printf("[%0.18f\t%0.18f]\n", t[i], t[i+1]);
	*/
	
	ROUND_DOWN;
	sum = 0.0;
	for( i=1; i<6; i++ ) {
		j = 2*i;
		jp = j + 1;
		yc = y_coefficients[i];
		if( i != 3 ) {
			if( yc > 0.0 ) {
				sum += yc*y[j];
			} else {
				sum += yc*y[jp];
			}
		}
	}
	out[0] = sum + rel_val[0];
	
	/*
	ROUND_NEAR;
	printf("ycontrib[0] = %.18f\n", sum);
	*/
	
	ROUND_UP;
	sum = 0.0;
	for( i=1; i<6; i++ ) {
		j = 2*i;
		jp = j + 1;
		yc = y_coefficients[i];
		if( i != 3 ) {
			if( yc > 0.0 ) {
				sum += yc*y[jp];
			} else {
				sum += yc*y[j];
			}
		}
	}
	out[1] = sum + rel_val[1];
		
	/*
	ROUND_NEAR;
	printf("ycontrib[1] = %.18f\n", sum);
	printf("rel_val = [%.18f, %.18f]\n", rel_val[0], rel_val[1]);
	*/
	
	if( do_gma ) {
		for( i=0; i<12; i+=2 ) {
			i_mult( y + i, t + i, compvals + i );
		}
		ROUND_DOWN;
		for( i=0; i<6; i++ ) {
			j = 2*i;
			if( i != 0 && i != 3 ) {
				out_partials[j] = 2.0*compvals[j] + y_coefficients[i];
			} else {
				out_partials[j] = 2.0*compvals[j];
			}
		}
		ROUND_UP;
		for( i=0; i<6; i++ ) {
			j = 2*i + 1;
			if( i != 0 && i != 3 ) {
				out_partials[j] = 2.0*compvals[j] + y_coefficients[i];
			} else {
				out_partials[j] = 2.0*compvals[j];
			}
		}
	} else {
		ROUND_DOWN;
		for( i=0; i<6; i++ ) {
			j = 2*i;
			if( i != 0 && i != 3 ) {
				out_partials[j] = t[j] + y_coefficients[i];
			} else {
				out_partials[j] = t[j];
			}
		}
		ROUND_UP;
		for( i=0; i<6; i++ ) {
			j = 2*i + 1;
			if( i != 0 && i != 3 ) {
				out_partials[j] = t[j] + y_coefficients[i];
			} else {
				out_partials[j] = t[j];
			}
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
*/


/*
cutval = 2.17667782820508515

2.8.5 (0,0,0) (0,0)

a23lo = 2.0
a23hi = cutval
a33lo = 2.0
a33hi = cutval
a12lo = 2.0
a12hi = 2.255
a43lo = 2.0
a43hi = 2.255

pracrange = {{2,2},{a12lo,a12hi},{2,cutval},{2,2},{2,2.2549},{2,2.51},
		{a23lo,a23hi},{2,2},{2,2.2549},
		{a33lo,a33hi},{2,2},{2,2.2549},
		{a43lo,a43hi},{2,2},{2,2.51},
		{2.7,2Sqrt[2]}}

{aval[1,2],aval[2,2],aval[3,2],aval[4,2],aval[5,2],aval[1,3],aval[2,3],
  aval[3,3],aval[4,3],aval[5,3],aval[1,5],aval[2,5],aval[3,5],aval[4,5],
  aval[5,5],aval[1,6],aval[2,6],aval[3,6],aval[4,6],aval[5,6],cval[1],cval[2],
  cval[3],cval[4],cval[5],bval}
  
List(-0.00823296394650475,0.004414106930963113,
   0.004414106930963113,-0.004414106930963113,
   0.00823296394650475,-0.004414106930963113,
   -0.004414106930963113,0.004414106930963113,
   -0.00823296394650475,0.00823296394650475,
   -0.03796881260612205,-0.03796881260612205,
   0.03796881260612205,-0.00898102748994389,
   0.00898102748994389,-0.00898102748994389,
   0.03796881260612205,0.03796881260612205,
   -0.03796881260612205,0.00898102748994389,
   0.7207090677395251,0.6015152457924683,
   0.43198356764416,0.7207090677395251,0.849421416626457,
   0.4886556269189748)

room = 0.0164924083423007576

-1.5 ptval - 0.06585- 0.022652 = -0.17156246850269623


2.8.5 (0,0,0) (0,1)

a23lo = 2.0
a23hi = cutval
a33lo = 2.0
a33hi = cutval
a12lo = 2.0
a12hi = 2.255
a43lo = 2.255
a43hi = 2.51

List(0.005566207946332779,0.01784898729958728,
   0.01784898729958728,0.1677799590658981,
   -0.1000787252571603,-0.01784898729958728,
   -0.01784898729958728,-0.1677799590658981,
   0.1000787252571603,-0.005566207946332779,
   0.05176749878922826,-0.05176749878922826,
   -0.05223270272861679,-0.03069118952105154,
   -0.00626664046699632,0.00626664046699632,
   -0.05176749878922826,0.05176749878922826,
   0.05223270272861679,0.03069118952105154,
   0.4435990582068949,0.7421717731697626,
   0.835894129424219,-0.1142629861583313,1.,
   0.4347030214601293)

room = 0.000903973742999059126


2.8.5 (0,0,0) (0,1)

a23lo = 2.0
a23hi = cutval
a33lo = 2.0
a33hi = cutval
a12lo = 2.255
a12hi = 2.51
a43lo = 2.255
a43hi = 2.51

List(0.02874944926886158,0.1070716351971938,
   -0.1070716351971938,0.1355829956515035,
   -0.02874944926886158,-0.1070716351971938,
   0.1070716351971938,-0.1355829956515035,
   0.02874944926886158,-0.02874944926886158,
   -0.1434062246271567,-0.1434062246271567,
   0.1434062246271567,0.122671774674842,
   0.1573359633463213,-0.1573359633463213,
   0.1434062246271567,0.1434062246271567,
   -0.1434062246271567,-0.122671774674842,
   0.869627711577136,-0.3342445481090603,
   0.005726355868492127,-0.1756970261626013,
   0.3188116070138338,0.07639731363595436)

room = 0.00652863062793700432


2.8.5 (0,0,1)

a23lo = 2.0
a23hi = cutval
a33lo = cutval
a33hi = cutvalhi

List(0.01693860366506827,0.06650513085303599,
   0.1363272579181817,-0.00933419519520795,
   -0.05115092026814959,-0.06650513085303599,
   -0.1363272579181817,0.00933419519520795,
   0.05115092026814959,-0.01693860366506827,
   -0.07923785656245241,-0.0935496987309429,
   0.1109119698021592,-0.07364896006711751,
   0.0838680336254256,-0.0838680336254256,
   0.07923785656245241,0.0935496987309429,
   -0.1109119698021592,0.07364896006711751,
   0.7199217153702022,0.4628448190857721,
   -0.4017911143548266,0.5934644176299197,
   0.3053355207886716,0.2393067362487727)

room = 0.000920864181930308767

2.8.5 (0,1,0)

a23lo = cutval
a23hi = cutvalhi
a33lo = 2.0
a33hi = cutval

List(0.01693860366506894,0.1363272579181793,
   -0.00933419519520917,-0.1011072844100964,
   -0.01693860366506894,-0.1363272579181793,
   0.00933419519520917,0.1011072844100964,
   0.01693860366506894,-0.01693860366506894,
   0.07923785656245119,-0.0959163210189348,
   0.0923342764958428,0.06148802693356736,
   -0.06148802693356736,0.06148802693356736,
   -0.07923785656245119,0.0959163210189348,
   -0.0923342764958428,-0.06148802693356736,
   0.2519024221326904,0.3574405778741409,
   -0.2520713654220794,0.5246067412331177,
   0.7978969827018572,0.2393067362487726)

room = 0.000920864181930308767

2.8.5 (0,1,1)

a23lo = cutval
a23hi = cutvalhi
a33lo =cutval
a33hi = cutvalhi

List(0.0169386036650696,0.1011072844100999,
   0.00933419519520617,-0.00933419519520617,
   -0.0180389236477455,-0.1011072844100999,
   -0.00933419519520617,0.00933419519520617,
   0.0180389236477455,-0.0169386036650696,
   -0.0923342764958428,0.0959163210189364,
   -0.0959163210189364,-0.0861686054175912,
   0.0872689254002688,-0.0872689254002688,
   0.0923342764958428,-0.0959163210189364,
   0.0959163210189364,0.0861686054175912,
   0.822120645900794,-0.2520713654220926,
   0.6672204907520682,0.2710711199294792,
   0.2072704532972271,0.2393067362487717)

room = 0.00808806136947968212

2.8.5 (1,0,1)

List(0.004996749481564588,0.,-0.1060310277065953,0.,
   -0.004996749481564588,0.,0.1060310277065953,0.,
   0.004996749481564588,-0.004996749481564588,
   -0.1031076114857766,-0.1019881007004714,
   -0.1031076114857766,0.07931614955685595,
   -0.07931614955685595,0.07931614955685595,
   0.1031076114857766,0.1019881007004714,
   0.1031076114857766,-0.07931614955685595,
   0.3070369470785845,0.05514644520007938,
   0.4837485991676782,-0.1053934988645204,
   0.7659042626869791,0.2082364359685616)

room = 0.00529843437379895476

2.8.5 (1,1,1)

List(0.03876164308154517,0.01836850892404906,
   -0.01836850892404906,-0.01836850892404906,
   -0.04538054589400819,-0.01836850892404906,
   0.01836850892404906,0.01836850892404906,
   0.04538054589400819,-0.03876164308154517,
   -0.102771474844139,0.102771474844139,
   -0.0889560665102218,-0.0822073757100353,
   0.0888262785224989,-0.0888262785224989,
   0.102771474844139,-0.102771474844139,
   0.0889560665102218,0.0822073757100353,
   0.6876759203076707,-0.1320077308568168,
   0.7424979034527878,0.2777452263490909,
   0.3641610502918215,0.2693788350757678)

room = 0.0151905524857279488
*/
