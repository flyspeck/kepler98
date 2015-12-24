/*  This file contains header information for the routines in sphere.c */

#define PI 		3.14159265358979323846
#define DOCT	0.7209029495174650928
#define POINT	0.055373645668463869730
#define SQRT2	1.414213562
#define SEED	10L
#define BUF		1000
#define PERT	10000.0


double aindf( double ls[7], int n1, int n2, int n3, int n4, int n5, int n6 );
double bvol( double lens[7] );
double afunc( double lens[7] );
double solid( double lens[7] );
double pf( int j, int i, int k, double len[5][5] );
double tvol( double lens[7] );
double tetvol( double lenlist[7] );
double norm( double v[4] );
void vtolen( double vects[4][4], double lens[7] );
void lentov( double lens[7], double vects[4][4] );
double det3( double a[4][4] );
void solve3( double a[4][4], double b[4], double c[4] );
double volvects( double vlist[4][4] );
void cross_product( double avector[4], double bvector[4], double result[4] );
double gma( double edges[7] );
double gmapt( double edges[7] );
double old_dih( double edges[7] );
double dih( double edges[7] );
int ccenter( double vects[4][4], double cent[4] );
int circumcent( double lens[7], double cent[4] );
double circumrad( double lens[7] );
double dbvol( double lens[7] );
double bigdelta( double lens[7] );
double tetvolume( double yvect[7] );
int istet( double yvect[7] );
double tomsrho( double xvect[7] );
double circumradius( double yvect[7] );
double tomsu( double xvect[4] );
double tomsv( double xvect[7] );
double auxPfun( double xvect[4] );
double tomsP( double xvect[7] );
double tomschi( double xvect[7] );
double voronoivol( double xvect[7] );
double vor( double edges[7] );
double vorpt( double edges[7] );
double crad3len( double lens[4] );
int isqrtet( double edges[7] );
int issmall( double edges[7] );
int oppval( int val );
int areedgefacessmall( double edges[7], int longedge );
int areothersshort( double edges[7], int longedge );
double scoresmall( double edges[7], int isgroup );
double genscore( double edges[7], int isgroup );
double score( double tet[7] );
double gscore( double tet[7] );
double scorept( double tet[7] );
double gscorept( double tet[7] );
double crossdiag( double len[10] );
int iscrossok( double len[10] );
double anglesum( double octa[13] );
double seciter( double x0, double x1, double roct[13] );
double octadiag( double roct[13] );
int isoctaok( double octa[13] );
double octascore( double octa[13] );
double octasph( double octa[13] );
void octapair( double octa[13], double xy[2] );
int isquadok( double len[10] );
void sph( double vect[4], double rvect[4] );
void penta( double **vects );
void hcp( double **vects );
double quaddiff( double quad[10] );
/* double mrand( void ); */
/* end of prototypes */
