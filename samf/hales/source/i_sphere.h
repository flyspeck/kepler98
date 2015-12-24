/*  This file contains header information for the routines 
in i_sphere.c */

#define PI 		3.14159265358979323846
#define DOCT	0.7209029495174650928
#define POINT	0.055373645668463869730
#define SQRT2	1.414213562
#define SEED	10L
#define BUF		1000
#define PERT	10000.0


void i_bvol( double y[12], double sqrtdelta[2], 
	double out[2] );
void i_afunc( double y[12], double out[2] );
void i_solid( double y[12], double sqrtdelta[2],
	double out[2] );
void i_tvol( double y[12], double out[2] );
void i_gma( double y[12], double sqrtdelta[2], 
	double out[2] );
void i_dih( double x[12], double out[2] );
void i_dih_alt( double x[12], double delta[2], double out[2] );
void i_dih_old( double x[12], double out[2] );
void i_dih_best( double x[12], double dih_part[12], 
	double out[2] );
void i_dbvol( double y[12], double out[2] );
void i_dbvol_alt( double y[12], double sqrtdelta[2],
	double out[2] );
void i_bigdelta( double x[12], double out[2] );
void i_bigdelta_best( double x[12], double out[2] );
void i_bigdelta_partials_best( double x[12], double delta[2],
	double delta_part[12] );
int i_istet( double y[12] );
void i_tomsrho( double x[12], double out[2] );
void i_crad( double y[12], double out[2] );
void i_crad2( double x[12], double out[2] );
void i_crad2_best( double x[12], double out[2] );
void i_tomsu( double x[6], double out[2] );
void i_tomsv( double x[12], double out[2] );
void i_auxPfun( double x[6], double out[2] );
void i_tomsP( double x[12], double out[2] );
void i_tomsPwod( double x[12], double out[2] );
void i_tomschi( double x[12], double out[2] );
double max_tomschi( double x[12] );
double s_max_tomschi( double x[12] );
double s_min_tomschi( double x[12] );
int s_crad2_6( double x[12], double out[2] );
void s_crad2_1( double x[12], double out[2], int finfo[2] );
void i_rogersvol2( double xyz[6], double out[2] );
void i_raw_rogersvol2( double xyz[6], double out[2] );
void i_rogers_density( double abc[6], double out[2] );
void old_i_rogers_density( double abc[6], double out[2] );
void i_voronoivol( double x[12], double sqrtdelta[2], 
	double out[2] );
void s_density_vor_6( double y[12], double x[12], 
	double out[2] );
void s_density_vor_1( double y[12], double x[12], 
	double out[2] );
void i_vor( double y[12], double x[12], double sqrtdelta[2], 
	double out[2] );
void i_vor_alt( double x[12],
	double sqrtdelta[2], double sol[2], double out[2] );
void best_i_vor_1( double y[12], double x[12], 
	double sqrtdelta[2], double out[2] );
void best_i_vor_6( double y[12], double x[12], 
	double sqrtdelta[2], double out[2] );
void i_crad3len( double y[6], double out[2] );
void i_crad3x2( double x[6], double out[2] );
int i_isqrtet( double y[12] );
int scoring_system( double x[12] );
int old_scoring_system( double y[12], double x[12] );
int octa_scoring( double x[12] );
int old_octa_scoring( double y[12], double x[12] );
void i_score( double y[12], double x[12],
	double sqrtdelta[2], double out[2] );
void s_octa_vor( double y[12], double x[12], double out[2] );
void best_octa_vor( double y[12], double x[12], 
	double sqrtdelta[2], double out[2] );
void i_octavor( double y[12], double x[12], 
	double sqrtdelta[2], double sol[2], double out[2] );
void i_rog_sph( double ab[4], double out[2] );
void i_wedge_vol( double alpha[2], double a[2], double out[2] );
void obtuse_wedge_vol( double alpha[2], double a[2], 
	double out[2] );
void i_wedge_sph( double alpha[2], double a[2], double out[2] );
void obtuse_wedge_sph( double alpha[2], double a[2], 
	double out[2] );
void i_dih_rog( double xy[4], double out[2] );
void old_dih_rog( double xy[4], double out[2] );
void i_vor_trunc( double y[12], double x[12], 
	double sqrtdelta[2], double out[2] );
void alt_vor_trunc( double y[12], double x[12], 
	double sph_in[2], double out[2] );
void obtuse_vor_trunc( double y[12], double x[12], 
	double sph_in[2], double out[2] );
void test_rog_wed( double y[12], double x[12], 
	double vol[2], double sol[2] );
void test_fullwed( double y[12], double x[12], 
	double vol[2], double sol[2] );



/* end of prototypes */
