/* i_voronoi.h contains prototypes for i_voronoi.c  */
/* (c) 1998, Samuel Ferguson  */




/* Prototypes */
void i_quoin_all( double abc[6], double xyz[6],
	double out_quo[2], double out_quo_pars[4], 
	double out_quo_secpars[6] );
void i_quoin( double abc[6], double xyz[6],
	double out_quo[2] );
void i_quoin_pars( double abc[6], double xyz[6],
	double out_quo[2], double out_quo_pars[4] );
void i_vor_beta( double hval[2], double out[2] );
void i_vor_beta_h( double hval[2], double out[2] );
void i_vor_beta_hh( double hval[2], double out[2] );
void i_vor_beta_x( double hval[2], double oohval[2],
	double out[2] );
void i_vor_beta_xx( double hval[2], double oohval[2],
	double beta_x[2], double out[2] );
void i_vor0( double y[12], double x[12], 
	double sqrtdelta[2], double out[2] );
void i_vor0_partials( double y[12], double x[12], 
	double out_partials[12] );
void clear_vor0_partials( double y[12], double x[12], 
	double out_partials[12] );
void clear_vor0_xpars( double y[12], double x[12], 
	double oo2y[12], double out_partials[12] );
void i_index_dih( int index, double xp[12], double out[2] );
void i_index_dih_partials( int index, double yp[12],
	double xp[12], double out_partials[12] );
void clear_index_dih_partials( int index, double yp[12],
	double xp[12], double out_partials[12] );
void clear_index_dih_xi( int index, double yp[12],
	double xp[12], double out_partials[12] );
void crad3len_y1y1( double y[6], double x[6], double out[2] );
void i_vor0_y1y1( double y[12], double x[12], 
	double out[2] );
void clear_vor0_y1y1( double y[12], double x[12], 
	double out[2] );
void clear_vor0_x1x1( double y[12], double x[12], 
	double oo2y[12], double out[2] );
void i_sol_y1y1( double y[12], double x[12], double out[2] );
void i_betadihsum_y1y1( double y[12], double x[12], 
	double out[2] );
void clear_betadihsum_y1y1( double y[12], double x[12], 
	double out[2] );
void clear_betadihsum_x1x1( double y[12], double x[12], 
	double oo2y[12], double out[2] );
void i_quoinsum_y1( double y[12], double x[12], 
	double out[2] );
void i_quoinsum_y1y1( double y[12], double x[12], 
	double out[2] );
void superclear_vor0_y1( double y[12], double x[12], 
	double out[2] );
void super_u( double fx[6], double out[2] );
void superclear_vor0_x1x1( double y[12], double x[12],
	double out[2] );
void i_delta1_j( double x[12], double out[12] );
void cross_diag2( double x[12], double xp[12], 
	double out[2] );
void aux_cross( double x[12], double den[2], double out[2] );

