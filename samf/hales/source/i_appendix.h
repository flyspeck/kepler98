/* i_appendix.h, by Samuel Ferguson, (c) 1998. 			*/
/* Headers for auxiliary routines for doing 			*/
/* calculations in the appendix of Sphere Packings IV.	*/


/* External variables */

/* Global variables */

/* Prototypes */

void appendix_init( void );
void appendix_get_trash( void );
void phi_fun( double h[2], double t[2], double out[2] );
void a_fun( double h[2], double out[2] );
void b_fun( double y[2], double out[2] );
void b0_fun( double y[2], double out[2] );
void v0_fun( double y[12], double x[12], double out[2] );
void v1_fun( double y[12], double x[12], double out[2] );
void rog_vor( double abc[6], double xyz[6], double out[2] );
void rog_sol( double abc[6], double out[2] );
void rog_dih( double xyz[6], double out[2] );
void eta0_fun( double h[2], double out[2] );
void crown_fun( double h[2], double out[2] );
void anc_den( double abc[6], double xyz[6], double out[2] );
void anc_fun( double y[6], double x[6], double out[2] );
void old_anc_fun( double y[6], double x[6], double out[2] );
void older_anc_fun( double y[6], double x[6], double out[2] );
void max_anc_best( double y[6], double x[6], 
	double anc_2and6[2], double out[2] );
void kappa_fun( double y[12], double x[12], 
	double dih_val[2], double out[2] );
void dih_sign26( double x[12], double out[4] );
void dih_sign35( double x[12], double out[4] );
void set_fake_anc_const( void );
void fake_anc_fun( double y[6], double out[2] );
void rog_vol_ab( double abc[6], double xyz[6], double part[4] );
void rog_sol_ab( double abc[6], double xyz[6], double part[4] );
void rog_dih_ab( double abc[6], double xyz[6], double part[4] );
void rog_vor_ab( double abc[6], double xyz[6], double part[4] );
void anc_y6par( double y[6], double x[6], double out[2] );
void anc_2and6pars( double y[6], double x[6], double out[4] );
void fake_kappa_fun( double y[12], double dih_val[2], 
	double out[2] );
void fake_kappa_fun2( double y[12], double dih_val[2], 
	double out[2] );
void fake_anchor( double y[6], double out[2] );
void fake_anchor2( double y[6], double out[2] );
void fake_anchor_y6( double y[6], double out[2] );
void fake_anchor2_y6( double y[6], double out[2] );
void cos_arc_fun( double y[6], double x[6], double out[2] );
void cos2_beta_psi( double cos2psi[2], double cos2theta[2],
	double out[2] );
void cos2_dih3( double xp[12], double sign[2], double out[2] );
void cos2_dih2( double xp[12], double sign[2], double out[2] );
void cos2_rog_dih( double xyz[6], double out[2] );
void tconst_fun( int i, double out[2] );
void sconst_fun( int index, double out[2] );
void tomDfun( int n, int k, double out[2] );
void tomZfun( int n, int k, double out[2] );

