/* interval.h	(c) 1997, Samuel Ferguson.
Headers for routines in interval.c */



/* Prototypes */
void i_init( void );
void i_add( double x[2], double y[2], double out[2] );
void i_sub( double x[2], double y[2], double out[2] );
void i_negate( double x[2], double out[2] );
void i_smult( double x, double y[2], double out[2] );
void i_mult( double x[2], double y[2], double out[2] );
void i_div( double x[2], double y[2], double out[2] );
void i_sqrt( double x[2], double out[2] );
void i_atan( double x[2], double out[2] );
double max_atanquot( double x );
double min_atanquot( double x );
void i_acos( double x[2], double out[2] );
void i_squarelen( double len[2], double out[2] );
void set_infinity( double out[2] );
void i_old_recognize( double num, double out[2] );
void i_recognize( double num, double out[2] );

/* Interval constants */

/* Pi = 
[3.141592653589793116, 
 3.14159265358979356]	*/
#define PI_LO	3.14159265358979311600e+0
#define PI_HI	3.14159265358979356009e+0

/* Pi/2 = 
[1.570796326794896558, 
 1.57079632679489678]	*/
#define PI_2_LO	1.57079632679489655800e+0
#define PI_2_HI	1.57079632679489678004e+0

/* delta_oct = 
[0.7209029495174650304, 
 0.7209029495174651414]	*/
#define DOCT_LO	7.20902949517465030382e-1
#define DOCT_HI	7.20902949517465141405e-1

/* 2 Pi/5 = 
[1.256637061435917246, 
 1.256637061435917468]	*/
#define TWOPI_5_LO	1.25663706143591724640e+0
#define TWOPI_5_HI	1.25663706143591746844e+0

/* 2 atan( 0.25 sqrt(2) ) = 
[0.6796738189082438542, 
 0.6796738189082439652] */
#define SPHCUT_LO	6.79673818908243854153e-1
#define SPHCUT_HI	6.79673818908243965176e-1

/* Sqrt[2] = 
[1.414213562373094923, 
 1.414213562373095145]	*/
#define SQRT2_LO	1.41421356237309492343e+0
#define SQRT2_HI	1.41421356237309514547e+0

/* 2 Sqrt[2] = 
[2.828427124746189847, 
 2.828427124746190291]	*/
#define TWOSQRT2_LO	2.82842712474618984686e+0
#define TWOSQRT2_HI	2.82842712474619029095e+0

/* 1/3 = 
[0.3333333333333333148, 
 0.3333333333333333703]	*/
#define ONE_3_LO	3.33333333333333314830e-1
#define ONE_3_HI	3.33333333333333370341e-1

/* 2/3 = 
[0.6666666666666666297, 
 0.6666666666666667407] */
#define TWO_3_LO	6.66666666666666629659e-1
#define TWO_3_HI	6.66666666666666740682e-1

/* 1/6 = 
[0.1666666666666666574, 
 0.1666666666666666852]	*/
#define ONE_6_LO	1.66666666666666657415e-1
#define ONE_6_HI	1.66666666666666685170e-1

/* 1/12 = 
[0.08333333333333332871, 
 0.08333333333333334259]	*/
#define ONE_12_LO	8.33333333333333287074e-2
#define ONE_12_HI	8.33333333333333425852e-2

/* 1/36 = 
[0.02777777777777777624, 
 0.02777777777777777971]	*/
#define ONE_36_LO	2.77777777777777762358e-2
#define ONE_36_HI	2.77777777777777797052e-2

/* 1/48 =  */
#define ONE_48_LO	2.08333333333333321769e-2
#define ONE_48_HI	2.08333333333333356463e-2

/* 1/sqrt(2) = 
[0.707106781186547, 
 0.707106781186548]	*/
#define ONE_SQRT2_LO	7.07106781186547461715e-1
#define ONE_SQRT2_HI	7.07106781186547683760e-1

/* 1/(2*sqrt(2)) =  */
#define ONE_2SQRT2_LO	3.53553390593273730858e-1
#define ONE_2SQRT2_HI	3.53553390593273841880e-1

/* 2.51 =   
2.50999999999999978684e+0
 2.51000000000000023093e+0 */
#define TWO51_LO			2.50999999999999978684
#define TWO51_HI			2.51000000000000023093

/* 4/2.51 = 
 1.59362549800796804433e+0
 1.59362549800796826638e+0 */
#define FB251_LO			1.59362549800796804433
#define FB251_HI			1.59362549800796826638

/* 2^-49 = 
 1.77635683940025046468e-15 */
#define ATANERR				1.77635683940025046468e-15

 
