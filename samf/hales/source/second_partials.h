/*  This file contains header information for the routines 
in second_partials.c */


void sp_afunction( double x[12], double apars[12],
	double out[2] );
void sp_newafun( double x[12], double apars[12], 
	double out[2] );
void sp_bfunction( double x[12], double bpars[12],
	double out[2] );
void sp_cfunction( double x[12], double cpars[12],
	double out[2] );
void sp_dfunction( double x[12], double delta_part[12],
	double out_delta[2], double out_sqrtdelta[2] );
void sp_ufunction( double aval[2], double bval[2], 
	double out[2] );
void sp_ubv_i( double uval[2], double vval[2],
	double uval_i[12], double vval_i[12], double out[12] );
void sp_u_i( double aval[2], double bval[2],
	double aval_i[12], double bval_i[12], double out[12] );
void sp_u_ij( double aval[2], double bval[2],
	double aval_i[12], double bval_i[12],
	double aval_ij[6][12], double bval_ij[6][12], 
	double out[6][12] );
void sp_ubv_ij( double uval[2], double vval[2],
	double uval_i[12], double vval_i[12],
	double uval_ij[6][12], double vval_ij[6][12], 
	double out[6][12] );
void sp_apars( double x[12], double out[12] );
void sp_newapars( double x[12], double out[12] );
void sp_bpars( double x[12], double b_ij[6][12],
	double out[12] );
void sp_cpars( double x[12], double out[12] );
void sp_dpars( double delta_pars[12],
	double sqrt_delta[2], double out[12] );
void sp_deltapars( double rx[12], double delta_ij[6][12],
	double part[12] );
void sp_asecondpars( double bigout[6][12] );
void sp_newasecpars( double bigout[6][12] );
void sp_bsecondpars( double x[12], double bigout[6][12] );
void sp_csecondpars( double bigout[6][12] );
void sp_dsecondpars( double delta[2], double sqrtdelta[2], 
	double delta_pars[12], double delta_ij[6][12],
	double out[6][12] );
void sp_deltasecondpars( double x[12], double bigout[6][12] );
void sp_secparbounds( double x[12], double out[6][12] );
void sp_vorsecparbds( double x[12], double out[6][12] );
void sp_findmaxmin( double ubv_ij[6][12], double out[2] );
void sp_solasecpars( double y[12], double a_ij[6][12] );
void sp_sold_i( double y[12], double d_xi[12], 
	double out[12] );
void sp_sold_ij( double y[12], double d_xi[12], 
	double d_xij[6][12], double out[6][12] );
void sp_sol_ij( double c[2], double c_i[12],
	double c_ij[6][12], double sol_ij[6][12] );
void sp_solsecparbds( double y[12], double x[12],
	double out[6][12] );
void sp_solsecyparbds( double y[12], double x[12],
	double out[6][12] );
void sp_dih_a( double x[12], double out[2] );
void sp_dih_c( double x[12], double cpars[12],
	double out[2] );
void sp_dih_apars( double x[12], double out[12] );
void sp_dih_cpars( double x[12], double out[12] );
void sp_dihasecpars( double a_ij[6][12] );
void sp_dihcsecpars( double c_ij[6][12] );
void sp_dihsecparbds( double x[12], double out[6][12] );
void sp_ytoxpars( double oo2y[12], double ypars[12],
	double xpars[12] );
void sp_ytoxsecpars( double oo2y[12], double xpars[12],
	double ysecpars[6][12], double xsecpars[6][12] );
void sp_xtoypars( double y[12], double xpars[12],
	double ypars[12] );
void sp_xtoysecpars( double y[12], double xpars[12],
	double xsecpars[6][12], double ysecpars[6][12] );
void sp_gmavolsecparbds( double x[12], double out[6][12] );
void dih_bfun( double x[12], double delta[2], double out[2] );
void dih_bfun_i( double x[12], double delta[2],
	double delta_i[12], double out[12] );
void dih_bfun_ij( double x[12], double delta_i[12],
	double delta_ij[6][12], double out[6][12] );
void dih_afun( double delta_i[12], double out[2] );
void dih_afun_i( double delta_ij[6][12], double out[12] );
void dih_afun_ij( double out[6][12] );
void dih_cfun( double a[2], double b[2], double out[2] );
void dih_cfun_i( double a[2], double b[2], double a_i[12],
	double b_i[12], double out[12] );
void sol_bfun( double delta[2], double out[2] );
void sol_bfun_i( double y[12], double delta_i[12], 
	double out[12] );
void sol_bfun_ij( double y[12], double delta_i[12],
	double delta_ij[6][12], double out[6][12] );
void sol_afun( double y[12], double out[2] );
void sol_afun_i( double y[12], double x[12], double part[12] );
void sol_afun_ij( double y[12], double out[6][12] );
void sol_cfun( double a[2], double b[2], double out[2] );
void sol_cfun_i( double a[2], double b[2], double a_i[12],
	double b_i[12], double out[12] );
void clear_pars( double a[2], double b[2], double c[2],
	double a_i[12], double b_i[12], double out[12] );
void clear_secpars( double a[2], double b[2], double c[2],
	double a_i[12], double b_i[12], double c_i[12],
	 double a_ij[6][12], double b_ij[6][12], 
	 double out[6][12] );
void clear_dih_i( double y[12], double x[12], 
	double clear_i[12] );
void clear_dih_ij( double y[12], double x[12], 
	double clear_i[12], double clear_ij[6][12] );
void clear_sol_i( double y[12], double x[12], 
	double clear_i[12] );
void clear_sol_ij( double y[12], double x[12], 
	double clear_i[12], double clear_ij[6][12] );


/* end of prototypes */
