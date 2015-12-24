/* system_headers.h, a header file which includes the appropriate header information,
depending on the target system. */

/* Currently, this file only supports PPC machines (PowerMacs)
and SGI's (blackbox).  */

#include <stdio.h>

#define DA_SYSTEM	2
						/* 1:  PowerPC	*/
						/* 2:  Sparc	*/
						/* 3:	SGI		*/
						
#if DA_SYSTEM == 1	/* PPC code */

#include <fp.h>
#include <fenv.h>
#include <stdlib.h>
#include <time.h>

#define GOT_FMADD	1	/* have fmadd instruction */

#elif DA_SYSTEM == 2	/* Sparc (Sun) code */

#include <math.h>
/* #include <sys/ieeefp.h> */
#include <ieeefp.h>
#include <stdlib.h>
#include <time.h>

/* #define CLOCKS_PER_SEC	1.0e6 */	/* microseconds */
#define GOT_FMADD	0	/* no fmadd instruction */

#elif DA_SYSTEM == 3	/* generating SGI code */

#include <math.h>
#include <ieeefp.h>
#include <stdlib.h>
#include <time.h>

#define GOT_FMADD	1	/* have fmadd instruction */

#else					/* default */

#include <math.h>
#include <ieeefp.h>
#include <stdlib.h>
#include <time.h>

#define GOT_FMADD	0	/* no fmadd instruction */

#endif
