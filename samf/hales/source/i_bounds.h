/* i_bounds.h contains prototypes for i_bounds.c  */
/* (c) 1997, Samuel Ferguson  */




/* Prototypes */

double rough_min_delta( double x[12] );
double rough_max_delta( double x[12] );
double min_a( double y[12] );
double max_a( double y[12] );
double s_min_u( double x[6] );
double s_max_u( double x[6] );
double s_min_delta4( double x[12] );
double s_max_delta4( double x[12] );
void i_delta4( double x[12], double out[2] );
double rough_min_tvol( double y[12] );
double rough_max_tvol( double y[12] );
double rough_max_solid( double y[12] );
double max_solid( double y[12], double sqrtdelta[2] );
double rough_min_bvol( double y[12] );
double rough_max_bvol( double y[12] );
double s_min_bvol( double y[12], double sqrtdelta[2] );
double s_max_bvol( double y[12], double sqrtdelta[2] );
double s_max_dih( double x[12] );
double s_min_dih( double x[12] );
void s_dih( double x[12], double out[2] );
double max_acos( double xv );
double min_acos( double xv );
int max_qr_crad_test( double y[12] );
int min_qr_crad_test( double y[12] );
double max_qr_crad( double y[12] );
double min_qr_crad( double y[12] );
double rough_min_gma( double y[12] );
double rough_max_gma( double y[12] );
double s_min_gma( double y[12], double sqrtdelta[2] );
double s_max_gma( double y[12], double sqrtdelta[2] );
void delta_partial( int n, double x[12], double out[2] );
void i_delta_partials( double x[12], double part[12] );
void s_delta_partials( double x[12], double part[12] );
int s_delta_partial_sign( int n, double x[12] );
void a_partials( double y[12], double part[12] );
void s_solid_partials( double y[12], double delta[2], 
  double sqrtdelta[2], double delta_part[12], double sol_part[12] );
void s_solid_xpars( double y[12], double ooy[12],
	double delta[2], double sqrtdelta[2], 
	double delta_part[12], double sol_part[12] );
void s_bvol_partials( double y[12], double delta[2], 
	double sqrtdelta[2], double delta_part[12], 
	double sol_part[12], double bvol_part[12] );
void s_bvol_xpars( double y[12], double ooy[12],
	double delta[2], double sqrtdelta[2], 
	double delta_part[12], double sol_xpart[12], 
	double bvol_part[12] );
void s_uu_partials( double x[12], double u_2[2],
	double u_3[2], double uu_part[12] );
void i_uu_partials( double x[12], double u2[2],
	double u3[2], double uu_part[12] );
void delta4_partials( double x[12], double part[12] );
void s_dih_partials( double x[12], double y[12], 
	double dih_part[12] );
void i_dih_partials( double x[12], double y[12], 
	double dih_part[12] );
void i_dih_xpars( double x[12], double dih_part[12] );
void old_dih_xpars( double x[12], double dih_part[12] );
void s_dih_xpars( double x[12], double dih_part[12] );
void s_gma_partials( double y[12], double sqrtdelta[2], 
	double delta_part[12], double bvol_part[12], 
	double part[12] );
void s_gma_xpars( double oosqrtdelta[2], double delta_part[12],
	double bvol_xpart[12], double part[12] );
void auxp_partials( double x[6], double part[6] );
void tomsu_partials( double x[6], double part[6] );
void tomsv_partials( double x[12], double part[12] );
void i_auxp( double x[6], double out[2] );
void tomsP_partials( double x[12], double delta[2],
	double delta32[2], double deltapart[12], 
	double part[12] );
void vorvol_partials( double y[12], double x[12], 
	double delta[2], double sqrtdelta[2],
	double deltapart[12], double part[12] );
void vorvol_xpars( double x[12], double delta[2], 
	double sqrtdelta[2], double delta_part[12], 
	double part[12] );
void vor_partials( double y[12], double x[12], 
	double delta[2], double sqrtdelta[2], 
	double delta_part[12], double sol_part[12], 
	double part[12] );
void vor_xpars( double x[12], double delta[2], 
	double sqrtdelta[2], double delta_part[12], 
	double sol_xpart[12], double part[12] );
void octa_vor_partials( double y[12], double x[12], 
	double delta[2], double sqrtdelta[2], 
	double delta_part[12], double sol_part[12], 
	double part[12] );
void octa_vor_xpars( double y[12], double ooy[12],
	double x[12], double delta[2], double sqrtdelta[2], 
	double delta_part[12], double sol_xpart[12], 
	double part[12] );
void s_crad3x( double x[6], double out[2] );
void s_crad3x2( double x[6], double out[2] );
void face_edge( double x2[2], double x6[2], double out[2] );
void face_x6_x2( double x1[2], double x2[2], double out[2] );
void crad3len_pars( double y[6], double x[6],
	double part[6] );
void s_crad3len_pars( double y[6], double x[6],
	double part[6] );
void crad3len2_xpars( double x[6], double part[6] );
void s_crad3len2_xpars( double x[6], double part[6] );
void alpha_ab( double xy[4], double part[4] );
void wedgevol_ab( double ab[4], double xy[4], double al[2], 
	double al_ab[4], double part[4] );
void wedgesol_ab( double ab[4], double al[2],
	double al_ab[4], double part[4] );
void rogvol_ab( double ab[4], double xy[4], double part[4] );
void rogsol_ab( double ab[4], double xy[4], double part[4] );
void old_rogsol_ab( double ab[4], double xy[4], double part[4] );
void wedge_pars( double y[12], double x[12],
	double vol_pars[12], double sol_pars[12] );
void trunc_pars( double y[12], double x[12], 
	double sph_pars[12], double part[12] );
void approx_trunc_pars( double y[12], double part[12] );
void approx_rog_wed_pars( double y[12], double vpart[12],
	double spart[12] );
void approx_fullwed_pars( double y[12], double vpart[12],
	double spart[12] );
void test_ab_pars( double ab[4] );
void i_sign_dihpars( double x[12], double delta[2],
	double delta_i[12], double out[12] );
void i_tomsrho_xpars( double x[12], double out[12] );
void i_sign_crad2_xpars( double x[12], double delta[2], 
	double delta_partials[12], double out[12] );
