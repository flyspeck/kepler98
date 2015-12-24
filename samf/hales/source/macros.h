/* macros.h, a collection of macros useful for interval computations.
(c) 1997, Samuel Ferguson.  */

#if	DA_SYSTEM == 1	/* PPC, PowerMac */

#define	ROUND_UP		fesetround(FE_UPWARD)
#define ROUND_DOWN		fesetround(FE_DOWNWARD)
#define ROUND_NEAR		fesetround(FE_TONEAREST)
#define ROUND_ZERO		fesetround(FE_TOWARDZERO)
#define ISNAN( x )		isnan( x )

#elif DA_SYSTEM == 2	/* Sparc */

/*
#define	ROUND_UP		ieee_flags("set","direction","positive",0)
#define ROUND_DOWN		ieee_flags("set","direction","negative",0)
#define ROUND_NEAR		ieee_flags("set","direction","nearest",0)
#define ROUND_ZERO		ieee_flags("set","direction","tozero",0)
*/
#define ROUND_UP		fpsetround( FP_RP )
#define ROUND_DOWN		fpsetround( FP_RM )
#define ROUND_NEAR		fpsetround( FP_RN )
#define ROUND_ZERO		fpsetround( FP_RZ )
#define ISNAN( x )		isnan( x )

#elif DA_SYSTEM == 3	/* SGI */

#define ROUND_UP		fpsetround( FP_RP )
#define ROUND_DOWN		fpsetround( FP_RM )
#define ROUND_NEAR		fpsetround( FP_RN )
#define ROUND_ZERO		fpsetround( FP_RZ )
#define ISNAN( x )		isnand( x )

#else					/* default */

#define ROUND_UP		fpsetround( FP_RP )
#define ROUND_DOWN		fpsetround( FP_RM )
#define ROUND_NEAR		fpsetround( FP_RN )
#define ROUND_ZERO		fpsetround( FP_RZ )
#define ISNAN( x )		isnand( x )

#endif


#define MIN( a, b )		( (a) < (b) ? a : b )
#define MAX( a, b )		( (a) > (b) ? a : b )

#define I_ADD( x, y, out )		ROUND_DOWN;	\
								out[0] = (x)[0] + (y)[0];		\
								ROUND_UP;	\
								out[1] = (x)[1] + (y)[1]

#define I_SUB( x, y, out )		ROUND_DOWN;	\
								out[0] = (x)[0] - (y)[1];		\
								ROUND_UP;	\
								out[1] = (x)[1] - (y)[0]

#define I_NEGATE( x, out )		out[0] = -(x)[1];		\
								out[1] = -(x)[0]

#define I_SMULT( x, y, out )	if( (x) >= 0.0 ) {	\
									ROUND_DOWN;	\
									out[0] = (x)*(y)[0];			\
									ROUND_UP;	\
									out[1] = (x)*(y)[1];			\
									}							\
								else {							\
									ROUND_DOWN;	\
									out[0] = (x)*(y)[1];			\
									ROUND_UP;	\
									out[1] = (x)*(y)[0];			\
									}

#define I_MULT( x, y, out )		if( ((x)[0] >= 0.0) && (y)[0] >= 0.0 ) {	\
									ROUND_DOWN;	\
									out[0] = (x)[0]*(y)[0];			\
									ROUND_UP;	\
									out[1] = (x)[1]*(y)[1];			\
									}							\
								else							\
									i_mult( x, y, out )

#define I_SQRT( x, out )		ROUND_DOWN;	\
								if( (x)[0] < 0.0 )			\
									out[0] = 0.0;			\
								else						\
									out[0] = sqrt( (x)[0] );\
								ROUND_UP;	\
								out[1] = sqrt( (x)[1] )

#define I_ATAN( x, out )		ROUND_DOWN;	\
								out[0] = atan( (x)[0] ) - ATANERR;	\
								ROUND_UP;	\
								out[1] = atan( (x)[1] ) + ATANERR

#define I_SQUARELEN( x, out )	ROUND_DOWN;	\
								out[0] = (x)[0]*(x)[0];		\
								ROUND_UP;	\
								out[1] = (x)[1]*(x)[1]



