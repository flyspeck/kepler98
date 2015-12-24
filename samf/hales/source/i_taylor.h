/* i_taylor.h  This file contains headers for i_taylor.c */


/* Prototypes */

void t_dih( double t[12], double h[12], double h2, 
	double out[2] );
void t_dih_base( double xt[12], double out[2], 
	double out_i[12] );
void t_solid( double yt[12], double ooyt[12], 
	double delta[2],double sqrtdelta[2], 
	double delta_part[12],double h[12], 
	double h2, double out[2] );
void t_solid_base( double yt[12], double ooyt[12], 
	double delta[2], double sqrtdelta[2], double delta_part[12],
	double out[2], double out_i[12] );
void t_vorvol( double xt[12], double delta[2],
	double sqrtdelta[2], double delta_part[12], 
	double h[12], double h2, double out[2] );
void t_vor( double xt[12], double yt[12], double ooyt[12], 
	double delta[2],double sqrtdelta[2], double delta_part[12], 
	double h[12], double h2, double out[2] );
void t_vor_base( double xt[12], double yt[12],
	double delta[2],double sqrtdelta[2], double delta_part[12], 
	double solid_i[12], double out[2], double out_i[12] );
void t_gma( double yt[12], double ooyt[12],
	double delta[2], double sqrtdelta[2], 
	double delta_part[12], double h[12], double h2, 
	double out[2] );
void t_gma_base( double yt[12], double ooyt[12],
	double delta[2], double sqrtdelta[2], 
	double delta_part[12], double solid_i[12],
	double out[2], double out_i[12] );
void t_composite( int rel[7], double relconst[7], 
	double y[12], double x[12], double out[2], 
	double out_i[12], double outvals[14] );
void fat_composite( int rel[7], double relconst[14], 
	double y[12], double x[12], double out[2], 
	double out_i[12], double outvals[14] );
void t_print_const( void );
