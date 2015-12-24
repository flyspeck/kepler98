/* sphere.c, by Samuel Ferguson, (c) 1997. */
/* Routines for doing stuff with tetrahedrons and sphere packings */
/* Note: I'm doing Mathematica type indexing, not C indexing, normally.  Be warned. */

/* This presentation of the routines for sphere packing computations is taken from the	*/
/* presentation in sphere.m, the collection of Mathematica routines.  The intent of   	*/
/* this duplication is, to put it simply, speed.  The composite scoring function, 			*/
/* score(), is quite complex, and could be prohibitively slow for doing extensive				*/
/* computations.																																				*/


#include "system_headers.h"
/*	#include "dim.h"	*/
#include "sphere.h"

#define SQUARE( x )		( x )*( x )

/*
#define DEBUG			1
*/

/* Useful definitions */
#define TWOPI				6.283185307179586476925286766559
#define TWOSQRT2		2.8284271247461900976033774484193961571


/* These are constants which appear in the modified definition for small simplices. */
#define CIRCUMCUT		1.39			/* originally 2.78/2 */
#define EDGELENCUT	2.12			/* originally 2.12   */
#define LONGCUT			2.8284272			/* something less than 2 Sqrt[2] */


double aindf( double ls[7], int n1, int n2, int n3, int n4, int n5, int n6 )
{
  double ans, temp[7];
  
  temp[1] = ls[n1];
  temp[2] = ls[n2];
  temp[3] = ls[n3];
  temp[4] = ls[n4];
  temp[5] = ls[n5];
  temp[6] = ls[n6];
  
  ans = afunc( temp );
  return( ans );
}


double bvol( double ls[7] )
{
  double a[4], tp, sum, ans;
  int i;
  
  tp = 6*tvol( ls );
  
  a[0] = afunc( ls );
  a[1] = aindf( ls, 1, 5, 6, 4, 2, 3 );
  a[2] = aindf( ls, 2, 4, 6, 5, 1, 3 );
  a[3] = aindf( ls, 3, 4, 5, 6, 1, 2 );
  
  sum = 0.0;
  for( i=0; i<4; i++ )
    sum += atan( tp/a[i] );
  ans = 2*sum/3.0;
  return( ans );
}


double afunc( double y[7] )
{
  double p1, p2, p3, ans;
  
  p1 = y[2]*y[2] + y[3]*y[3] - y[4]*y[4];
  p2 = y[1]*y[1] + y[3]*y[3] - y[5]*y[5];
  p3 = y[1]*y[1] + y[2]*y[2] - y[6]*y[6];
  ans = y[1]*y[2]*y[3] + (y[1]*p1 + y[2]*p2 + y[3]*p3)/2.0;
  return( ans );
}


double solid( double lens[7] )
{
  double ans;
  ans = 2.0*atan( 6.0*tvol( lens )/afunc( lens ) );
  return( ans );
}


double pf( int j, int i, int k, double len[5][5] ) 
{
  double a, b, c;

  a = len[j][i];
  b = len[i][k];
  c = len[j][k];

return( (a*a + b*b - c*c)/2.0 );
}


double tvol( double lenlist[7] )
{
  double vol, det, v[4][4];
  
  lentov( lenlist, v );
  det = v[1][1]*v[2][2]*v[3][3];
  vol = det/6.0;
  return(vol);
}


double tetvol( double lenlist[7] )
{
  double len[5][5], p213, p214, p314, vol, t1, t2, t3, t4, dets, slist[7];
  int i,j;

  for( i=1; i<7; i++ )
    slist[i] = lenlist[i];
  slist[4] = lenlist[6];
  slist[6] = lenlist[4];
  for( i=1; i<4; i++ )
    for( j=i+1; j<5; j++ )
      len[j][i] = len[i][j] = slist[j+(i*(7-i))/2 - 4];

				/* format the list of lengths properly */

  p213 = pf(2,1,3, len);
  p214 = pf(2,1,4, len);
  p314 = pf(3,1,4, len);

  t1 = len[1][2]*len[1][3]*len[1][4];
  t2 = len[1][3]*p214;
  t3 = len[1][4]*p213;
  t4 = len[1][2]*p314;

  dets = t1*t1 + 2.0*p213*p214*p314 - t2*t2 - t3*t3 - t4*t4;

  vol = sqrt(dets)/6.0;

  return(vol);
}


double norm( double v[4] )
{
  int i;
  double sum;

  sum = 0.0;
  for( i=1; i<4; i++ )
    sum += v[i]*v[i];

  return( sqrt(sum) );
}


void vtolen( double vs[4][4], double len[7] )
{
  double v1[4], v2[4], v3[4], v12[4], v13[4], v23[4];
  int i;
  
  for( i=1; i<4; i++ )
    v1[i] = vs[1][i];
  for( i=1; i<4; i++ )
    v2[i] = vs[2][i];
  for( i=1; i<4; i++ )
    v3[i] = vs[3][i];
  for( i=1; i<4; i++ )
    v12[i] = vs[1][i] - vs[2][i];
  for( i=1; i<4; i++ )
    v13[i] = vs[1][i] - vs[3][i];
  for( i=1; i<4; i++ )
    v23[i] = vs[2][i] - vs[3][i];
  len[1] = norm( v1 );
  len[2] = norm( v2 );
  len[3] = norm( v3 );
  len[4] = norm( v23 );
  len[5] = norm( v13 );
  len[6] = norm( v12 );
}


void lentov( double len[7], double vs[4][4] )
{
  double x1, x2, x3, y2, y3, z3;
  
  x1 = len[1];
  x2 = (len[1]*len[1] + len[2]*len[2] - len[6]*len[6])/(2*len[1]);
  x3 = (len[1]*len[1] + len[3]*len[3] - len[5]*len[5])/(2*len[1]);
  y2 = sqrt( len[2]*len[2] - x2*x2 );
  y3 = -x2*x3/y2 + (len[2]*len[2] + len[3]*len[3] - len[4]*len[4])/(2*y2);
  z3 = sqrt( len[3]*len[3] - x3*x3 - y3*y3 );
  vs[1][1] = x1;
  vs[1][2] = 0.0;
  vs[1][3] = 0.0;
  vs[2][1] = x2;
  vs[2][2] = y2;
  vs[2][3] = 0.0;
  vs[3][1] = x3;
  vs[3][2] = y3;
  vs[3][3] = z3;
}


double det3( double a[4][4] )
{
  double d, p1, p2, p3, m1, m2, m3, t1, t2;

  p1 = a[1][2]*a[2][3]*a[3][1];
  p2 = a[1][3]*a[2][1]*a[3][2];
  p3 = a[1][1]*a[2][2]*a[3][3];

  m1 = a[1][3]*a[2][2]*a[3][1];
  m2 = a[1][1]*a[2][3]*a[3][2];
  m3 = a[1][2]*a[2][1]*a[3][3];

  t1 = p1 + p2 + p3;
  t2 = m1 + m2 + m3;

  d = t1 - t2;

  return(d);
}


void solve3( double a[4][4], double b[4], double c[4] )
{
					/* gives c=Inverse[a] . b 	*/
  double d, ia[4][4], sum;
  int i, k;

  d = det3(a);

  ia[1][1] = -(a[2][3]*a[3][2]) + a[2][2]*a[3][3];
  ia[1][2] =  a[1][3]*a[3][2] - a[1][2]*a[3][3];
  ia[1][3] = -(a[1][3]*a[2][2]) + a[1][2]*a[2][3];

  ia[2][1] = a[2][3]*a[3][1] - a[2][1]*a[3][3];
  ia[2][2] =  -(a[1][3]*a[3][1]) + a[1][1]*a[3][3];
  ia[2][3] = a[1][3]*a[2][1] - a[1][1]*a[2][3];

  ia[3][1] = -(a[2][2]*a[3][1]) + a[2][1]*a[3][2];
  ia[3][2] =  a[1][2]*a[3][1] - a[1][1]*a[3][2];
  ia[3][3] = -(a[1][2]*a[2][1]) + a[1][1]*a[2][2];

  for( i=1; i<4; i++ ) {
    sum = 0.0;
    for( k=1; k<4; k++ )
      sum += ia[i][k]*b[k];
    c[i] = sum/d;
  };

}


double volvects( double vlist[4][4] )
{
  return( fabs( det3(vlist) )/6.0 );
}


double gma( double edges[7] )
{
  double fred;

  fred = -DOCT*tvol(edges) + bvol(edges);
  return(fred);
}


double gmapt( double edges[7] )
{
  double fred;
  fred = gma( edges )/POINT;
  return(fred);
}


void cross_product( double avector[4], double bvector[4], double result[4] )
{
	result[1] = avector[2]*bvector[3] - avector[3]*bvector[2];
	result[2] = avector[3]*bvector[1] - avector[1]*bvector[3];
	result[3] = avector[1]*bvector[2] - avector[2]*bvector[1];
}


double old_dih( double lens[7] )
{
	double vects[4][4], v[7], vp[7], tp, bot, ans, sum;
	int i;
	
	lentov( lens, vects );
	cross_product( vects[1], vects[2], v );
	cross_product( vects[1], vects[3], vp );
	sum = 0.0;
	for( i=1; i<4; i++ )
		sum += v[i]*vp[i];
	tp = sum;
	bot = norm( v )*norm( vp );
	ans = acos( tp/bot );
	return( ans );
}


double dih( double lens[7] )
{
	double dx4, u1, u2;
	double x1, x2, x3, x4, x5, x6;
	double xv[4];
	
	x1 = SQUARE(lens[1]);
	x2 = SQUARE(lens[2]);
	x3 = SQUARE(lens[3]);
	x4 = SQUARE(lens[4]);
	x5 = SQUARE(lens[5]);
	x6 = SQUARE(lens[6]);
	
	dx4 = x1*(-x1 + x2 + x3 - 2*x4 + x5 + x6) +
		(x5 - x3)*(x2 - x6);
	xv[1] = x1;
	xv[2] = x2;
	xv[3] = x6;
	u1 = tomsu( xv );
	xv[2] = x3;
	xv[3] = x5;
	u2 = tomsu( xv );
	x1 = sqrt(u1*u2);
	x2 = dx4/x1;
	return( acos( x2 ) );
}

 /*
dihedral[lens__]:=
Module[{vects,v,vp, tp, bot, ans},
	vects = lentov[lens];
	v = crossprod[ vects[[1]], vects[[2]] ];
	vp = crossprod[ vects[[1]], vects[[3]] ];
	tp = v . vp;
	bot = norm[v] norm[vp];
	ans = ArcCos[tp/bot];
	Return[ans]
	]
*/


/*
double density( double **config )
{
  return( DOCT/(1.0 - score(config)*POINT*3.0/(16.0*PI)) );
}
*/

int ccenter( double vects[4][4], double cent[4] )
{
  double vects2[4], sum;
  int i, j;

  for( i=1; i<4; i++ ) {
    sum = 0.0;
    for( j=1; j<4; j++ )
      sum += vects[i][j]*vects[i][j];
    vects2[i] = sum/2.0;
  }; 
  if( fabs( det3(vects) ) > 0.01 ) {
    solve3(vects, vects2, cent);
    return( 1 );
    } 
  else
    return( 0 );
}


int circumcent( double lens[7], double cent[4] )
{
	double vect[4][4];
	int val;
	
	lentov( lens, vect );
	val = ccenter( vect, cent );
	return( val );
}

/*
circumrad[lens__]:=
	norm[circumcent[lens]]
*/

double circumrad( double lens[7] )
{
	double cent[4], rad;
	
	circumcent( lens, cent );
	rad = norm( cent );
	return( rad );
}


/*
dbvol[ls__]:=
Module[{vects, tp, l, a, i},
	vects = lentov[ls];
	tp = Abs[Det[vects]];
	a = afunc[ ls ];
	v = ArcTan[ tp/a ];
	ans = 2 v/3;
	Return[ ans ];
	]
*/

double dbvol( double lens[7] )
{
	double vects[4][4], tp, a, v, ans;
	
	lentov( lens, vects );
	tp = det3( vects );	/* Don't need fabs() because of the way lentov is set up */
	a = afunc( lens );
	v = atan2( tp, a );
	ans = 2.0*v/3.0;
	return( ans );
}


/*
bigdelta[{x1_,x2_,x3_,x4_,x5_,x6_}] :=
-(x1^2*x4) + x1*x2*x4 + x1*x3*x4 - 
	x2*x3*x4 - x1*x4^2 + 
	x1*x2*x5 - x2^2*x5 - x1*x3*x5 + 
	x2*x3*x5 + x1*x4*x5 + 
	x2*x4*x5 - x2*x5^2 - x1*x2*x6 + 
	x1*x3*x6 + x2*x3*x6 - 
	x3^2*x6 + x1*x4*x6 + x3*x4*x6 + 
	x2*x5*x6 + x3*x5*x6 - 
	x4*x5*x6 - x3*x6^2
*/

double bigdelta( double lens[7] )
{
	double x1, x2, x3, x4, x5, x6, val;
	
	x1 = lens[1];
	x2 = lens[2];
	x3 = lens[3];
	x4 = lens[4];
	x5 = lens[5];
	x6 = lens[6];
	
	val = - x1*x1*x4 + x1*x2*x4 + x1*x3*x4 - 
	x2*x3*x4 - x1*x4*x4 + 
	x1*x2*x5 - x2*x2*x5 - x1*x3*x5 + 
	x2*x3*x5 + x1*x4*x5 + 
	x2*x4*x5 - x2*x5*x5 - x1*x2*x6 + 
	x1*x3*x6 + x2*x3*x6 - 
	x3*x3*x6 + x1*x4*x6 + x3*x4*x6 + 
	x2*x5*x6 + x3*x5*x6 - 
	x4*x5*x6 - x3*x6*x6;

	return( val );
}


/*
tetvolume[{y1_,y2_,y3_,y4_,y5_,y6_}]:=
Block[{yvector, xvector, volume},
	yvector = {y1, y2, y3, y4, y5, y6};
	xvector = yvector^2;
	volume = Sqrt[bigdelta[xvector]]/12;
	Return[volume]
]
*/

double tetvolume( double yvect[7] )
{
	double xvect[7], volume;
	int i;
	
	for( i=1; i<7; i++ )
		xvect[i] = yvect[i]*yvect[i];
	volume = sqrt( bigdelta( xvect ) )/12.0;
	return( volume );
}


int istet( double yvect[7] )
{
	double xvect[7], volume;
	int i;
	
	for( i=1; i<7; i++ )
		xvect[i] = yvect[i]*yvect[i];
	volume = bigdelta( xvect );
	if( volume > 0.0 )
		return( 1 );
	else
		return( 0 );
}


/*
tomsrho[{x1_,x2_,x3_,x4_,x5_,x6_}]:=
	-(x1^2*x4^2) + 2*x1*x2*x4*x5 - x2^2*x5^2 + 
	2*x1*x3*x4*x6 + 2*x2*x3*x5*x6 - x3^2*x6^2
*/

double tomsrho( double xvect[7] )
{
	double x1, x2, x3, x4, x5, x6, val;
	
	x1 = xvect[1];
	x2 = xvect[2];
	x3 = xvect[3];
	x4 = xvect[4];
	x5 = xvect[5];
	x6 = xvect[6];
	
	val = 	- x1*x1*x4*x4 + 2*x1*x2*x4*x5 - x2*x2*x5*x5 + 
	2*x1*x3*x4*x6 + 2*x2*x3*x5*x6 - x3*x3*x6*x6;
	
	return( val );
}


/*
circumradius[{y1_,y2_,y3_,y4_,y5_,y6_}]:=
Block[{yvector, xvector, top, bot},
	yvector = {y1, y2, y3, y4, y5, y6};
	xvector = yvector^2;
	top = tomsrho[xvector];
	bot = bigdelta[xvector];
	Return[Sqrt[top/bot]/2]
]
*/

double circumradius( double yvect[7] )
{
	double xvect[7], top, bot;
	int i;
	
	for( i=1; i<7; i++ )
		xvect[i] = yvect[i]*yvect[i];
	top = tomsrho( xvect );
	bot = bigdelta( xvect );
	return( sqrt( top/bot )/2.0 );
}


/*
tomsu[{x1_,x2_,x3_}]:=
	-x1^2 - x2^2 - x3^2 + 2 x1 x3 + 2 x1 x2 + 2 x2 x3
*/

double tomsu( double xvect[4] )
{
	double x1, x2, x3, val;
	
	x1 = xvect[1];
	x2 = xvect[2];
	x3 = xvect[3];
	val = -x1*x1 - x2*x2 - x3*x3 + 2*x1*x3 + 2*x1*x2 + 2*x2*x3;
	return( val );
}


/*
tomsv[{x1_,x2_,x3_,x4_,x5_,x6_}]:=
	x1 x4(-x1 + x2 + x6) + x2 x5(x1 - x2 + x6) +
	x3 x6(x1 + x2 - x6) - 2 x1 x2 x6
*/

double tomsv( double xvect[7] )
{
	double x1, x2, x3, x4, x5, x6, val;
	
	x1 = xvect[1];
	x2 = xvect[2];
	x3 = xvect[3];
	x4 = xvect[4];
	x5 = xvect[5];
	x6 = xvect[6];
	
	val = 	x1*x4*(-x1 + x2 + x6) + x2*x5*(x1 - x2 + x6) +
	x3*x6*(x1 + x2 - x6) - 2*x1*x2*x6;

	return( val );
}


/*
auxPfun[{x1_,x2_,x3_}]:=
Block[{p1, p2},
	p1 = -x1^2 + 2 x1 x2 - x2^2 + x1 x3 + x2 x3;
	p2 = 48 tomsu[{x1, x2, x3}];
	Return[p1/p2]
]
*/

double auxPfun( double xvect[4] )
{
	double p1, p2, x1, x2, x3;
	
	x1 = xvect[1];
	x2 = xvect[2];
	x3 = xvect[3];
	
	p1  = -x1*x1 + 2*x1*x2 - x2*x2 + x1*x3 + x2*x3;
	p2 = 48*tomsu( xvect );
	return( p1/p2 );
}


/*
tomsP[{x1_,x2_,x3_,x4_,x5_,x6_}]:=
Block[{p1, p2, p3, demall, xvect},
	xvect = {x1, x2, x3, x4, x5, x6};
	p1 = auxPfun[{x1,x2,x6}];
	p2 = tomsv[xvect];
	p3 = Sqrt[bigdelta[xvect]];
	demall = p1 p2/p3;
	Return[demall]
]
*/

double tomsP( double xvect[7] )
{
	double p1, p2, p3, demall, xv[4];
	
	xv[1] = xvect[1];
	xv[2] = xvect[2];
	xv[3] = xvect[6];
	
	p1 = auxPfun( xv );
	p2 = tomsv( xvect );
	p3 = sqrt( bigdelta( xvect ) );
	demall = p1*p2/p3;
	return( demall );
}

/*
tomschi[{x1_,x2_,x3_,x4_,x5_,x6_}]:=
	x1 x4 (x5 + x6) + x2 x5 (x4 + x6) + 
		x3 x6 (x4 + x5) - 2 x4 x5 x6 - 
		x1 x4^2 - x2 x5^2 - x3 x6^2
*/

double tomschi( double x[7] )
{
	double p1, p2;
	
	p1 = x[1]*x[4]*(x[5] + x[6]) + x[2]*x[5]*(x[4] + x[6]) +
		x[3]*x[6]*(x[4] + x[5]);
	p2 = 2.0*x[4]*x[5]*x[6] + x[1]*x[4]*x[4] + 
		x[2]*x[5]*x[5] + x[3]*x[6]*x[6];
	return( p1 - p2 );
}

/*
voronoivol[{x1_,x2_,x3_,x4_,x5_,x6_}]:=
Block[{p1, p2, p3, xv, yv},
	xv = {x1, x2, x3, x4, x5, x6};
	yv = xv^2;
	p1 = tomsP[yv];
	p2 = tomsP[{yv[[2]], yv[[3]], yv[[1]], 
		yv[[5]], yv[[6]], yv[[4]]}];
	p3 = tomsP[{yv[[3]], yv[[1]], yv[[2]], 
		yv[[6]], yv[[4]], yv[[5]]}];
	Return[p1 + p2 + p3]
]
*/

double voronoivol( double xvect[7] )
{
	double yvect[7], p1, p2, p3, yv2[7], yv3[7];
	int i;
	
	for( i=1; i<7; i++ )
		yvect[i] = xvect[i]*xvect[i];
	p1 = tomsP( yvect );
	
	yv2[1] = yvect[2];
	yv2[2] = yvect[3];
	yv2[3] = yvect[1];
	yv2[4] = yvect[5];
	yv2[5] = yvect[6];
	yv2[6] = yvect[4];
	p2 = tomsP( yv2 );
	
	yv3[1] = yvect[3];
	yv3[2] = yvect[1];
	yv3[3] = yvect[2];
	yv3[4] = yvect[6];
	yv3[5] = yvect[4];
	yv3[6] = yvect[5];
	p3 = tomsP( yv3 );
	
	return( p1 + p2 + p3 );
}


/*
vor[edges__]:= 
	4( dbvol[edges] - doct voronoivol[edges] )
*/

double vor( double edges[7] )
{
	double val;
	
	val = dbvol( edges ) - DOCT*voronoivol( edges );
	return( 4*val );
}


/*
vorpt[edges__]:=
	vor[edges]/pt
*/

double vorpt( double edges[7] )
{
	return( vor( edges )/POINT );
}


/*
crad3len[{l1_,l2_,l3_}]:=
Block[{x,y,a,b,r},
	x = (l1^2 + l2^2 - l3^2)/(2 l1);
	y = Sqrt[l2^2 - x^2];
	a = l1/2;
	b = (l2^2 - l1 x)/(2 y);
	r = Sqrt[a^2 + b^2];
	Return[r]
]
*/

double crad3len( double lens[4] )
{
	double l1, l2, l3, x, y, a, b, r, l12, l22, l32;
	
	l1 = lens[1];
	l2 = lens[2];
	l3 = lens[3];
	l12 = l1*l1;
	l22 = l2*l2;
	l32 = l3*l3;
	x = (l12 + l22 - l32)/(2*l1);
	y = sqrt( l22 - x*x );
	a = l1/2.0;
	b = (l22 - l1*x)/(2*y);
	r = sqrt( a*a + b*b );
	return( r );
}


/*
isqrtet[edges__]:=
  Block[ {ok, i},
    ok = True;
    i = 1;
    While[ (ok == True && i <= 6),
      If[ Not[( 2 <= edges[[i]] <= 2.51 )], 
        ok = False ];
      i = i+1
    ];
    Return[ ok ]
  ]
*/

int isqrtet( double edges[7] )
{
	int ok, i;
	
	ok = 1;
	i = 1;
	while( (ok == 1 && i < 7) ) {
		if( !( edges[i] >= 2.0 && edges[i] <= 2.51 ) )
			ok = 0;
		i += 1;
		}
	return( ok );
}


/*
issmall[edges__]:=
Module[{c1, c2, c3, c4, clist, mx},
	c1 = crad3len[{edges[[1]],edges[[2]],edges[[6]]}];
	c2 = crad3len[{edges[[1]],edges[[3]],edges[[5]]}];
	c3 = crad3len[{edges[[2]],edges[[3]],edges[[4]]}];
	c4 = crad3len[{edges[[4]],edges[[5]],edges[[6]]}];
	clist = N[{c1, c2, c3, c4}];
	mx = Max[clist];
	If[mx < N[Sqrt[2]],
		Return[True],
		Return[False]
		]
	]
*/

int issmall( double edges[7] )
{
	double cval, old, list[4], max, temp;
	int i;
	
	max = 0.0;
	for( i=1; i<7; i++ ) {
		temp = edges[i];
		if( temp > max )
			max = temp;
		}
	
	if( max < LONGCUT ) {		/* We simply modify the definition of small--that seems easiest. */
		old = 0.0;
		list[1] = edges[1];
		list[2] = edges[2];
		list[3] = edges[6];
		cval = crad3len( list );
		if( cval > old )
			old = cval;
		
		list[1] = edges[1];
		list[2] = edges[3];
		list[3] = edges[5];
		cval = crad3len( list );
		if( cval > old )
			old = cval;
		
		list[1] = edges[2];
		list[2] = edges[3];
		list[3] = edges[4];
		cval = crad3len( list );
		if( cval > old )
			old = cval;
		
		list[1] = edges[4];
		list[2] = edges[5];
		list[3] = edges[6];
		cval = crad3len( list );
		if( cval > old )
			old = cval;
		
		if( old < SQRT2 )
			return( 1 );
		else
			return( 0 );
	} else
		return( 0 );	/* If the long edge is too long, we define the tet to be large. */
}


int oppval( int val )
{
	if( val > 3 )
		return( val - 3 );
	else
		return( val + 3 );
}

/*
edgefaces[long_]:=
Module[{faces},
	faces = {{{1,2,6},{1,3,5}},
		{{1,2,6},{2,3,4}},
		{{1,3,5},{2,3,4}},
		{{2,3,4},{4,5,6}},
		{{1,3,5},{4,5,6}},
		{{1,2,6},{4,5,6}}};
	Return[faces[[long]]]
	]
	

areedgefacessmall[edges__,long_]:=
Module[{faces, f1, f2, len1, len2, i, c1, c2, val, ans},
	faces = edgefaces[long];
	f1 = faces[[1]];
	f2 = faces[[2]];
	len1 = Table[edges[[ f1[[i]] ]], {i,1,3}];
	len2 = Table[edges[[ f2[[i]] ]], {i,1,3}];
	c1 = N[crad3len[len1]];
	c2 = N[crad3len[len2]];
	val = 2.78/2;
	If[ c1 < val && c2 < val, ans = True, ans = False ];
	Return[ans]
	]
*/


int areedgefacessmall( double edges[7], int longedge )
{
	int faces[6][2][3] = 
		{{{1,2,6},{1,3,5}},
		{{1,2,6},{2,3,4}},
		{{1,3,5},{2,3,4}},
		{{2,3,4},{4,5,6}},
		{{1,3,5},{4,5,6}},
		{{1,2,6},{4,5,6}}};
	double lens[4], crad[2], val;
	int i, j;
	
	val = CIRCUMCUT;
	for( j=0; j<2; j++ ) {
		for( i=1; i<4; i++ )
			lens[i] = edges[faces[longedge-1][j][i-1]];
		crad[j] = crad3len(lens);
		}
	if( crad[0] < val && crad[1] < val )
		return( 1 );
	else
		return( 0 );
}


/*
areothersshort[edges__, long_]:=
Module[{i, opp, ret},
	If[long < 4,
		opp = long + 3,
		opp = long - 3
		];
	ret = True;
	Do[If[edges[[i]]>2.12 && i!=long && i!=opp,
		ret = False],{i,1,6}];
	Return[ret]
	]
*/


int areothersshort( double edges[7], int longedge )
{
	int i, opp, ret;
	
	opp = oppval( longedge );
	ret = 1;
	for( i=1; i<7; i++ )
		if( edges[i] > EDGELENCUT && i != longedge && i != opp )
			ret = 0;
	return( ret );
}


/*
scoresmall[edges__, isgroup_]:=
Module[{sc, i, ind, ct, long, opp, long2, opp2, 
	foo, foo2, foo3, foo4, foo5, rad},
	ind = Table[0,{i,1,6}];
	ct = 0;
	Do[ If[edges[[i]] > 2.51,
			ct = ct + 1;
			ind[[ct]] = i;
			],
		{i,1,6}];
	If[ ct == 1,
		long = ind[[1]];
		If[long < 4,
			opp = long + 3,
			opp = long - 3
			];
		foo = areothersshort[edges, long];
		foo2 = areedgefacessmall[edges, long];
		foo3 = foo || foo2;
		foo4 = edges[[opp]] < 2.06;
		If[ isgroup, foo5 = foo3, foo5 = foo4 && foo3 ];
		If[foo5,
			sc = gma[edges];
			Return[sc], (*else*)
			sc = vor[edges];
			Return[sc]
			],
		(*else*)
		If[ ct == 2,
			long = ind[[1]];
			If[long < 4,
				opp = long + 3,
				opp = long - 3
				];
			long2 = ind[[2]];
			If[long2 < 4,
				opp2 = long2 + 3,
				opp2 = long2 - 3
				];
			If[ opp != long2 && edges[[long]] < 2.56
				&& edges[[long2]] < 2.56,
				sc = gma[edges];
				Return[sc],
				sc = vor[edges];
				Return[sc]
				];
			];
		];
	sc = vor[edges];
	Return[sc]
]
*/


double scoresmall( double edges[7], int isgroup )
{
	double sc;
	int ct, i, ind[7], big1, opp1, big2, opp2;
	int foo, foo2, foo3, foo4, foo5;
	
	ct = 0;
	for( i=1; i<7; i++ )
		if( edges[i] > 2.51 ) {
			ct++;
			ind[ct] = i;
			}
	if( ct == 1 ) {
		big1 = ind[1];
		opp1 = oppval( big1 );
		/* This is the part which is modified under the new scheme */
		foo = areothersshort( edges, big1 );
		foo2 = areedgefacessmall( edges, big1 );
		foo3 = (foo || foo2);
		foo4 = (edges[ opp1 ] < 2.06);
		if( isgroup )
			foo5 = foo3;
		else
			foo5 = (foo4 && foo3);
			
		if( foo5 ) {
			sc = gma( edges );
			return( sc );
			}
		else {
			sc = vor( edges );
			return( sc );
			}
		}
	else if( ct == 2 ) {
		big1 = ind[1];
		opp1 = oppval( big1 );
		big2 = ind[2];
		opp2 = oppval( big2 );
		if( opp1 != big2 && (edges[big1] < 2.56) && (edges[big2] < 2.56) ) {
			sc = gma( edges );
			return( sc );
			}
		else {
			sc = vor( edges );
			return( sc );
			}
		}
	sc = vor( edges );
	return( sc );
}


/*
genscore[tet__, isgroup_]:=
Module[{sc},
	If[ isqrtet[tet],
		If[ N[circumrad[tet]] <= 1.41,
			sc = gma[tet];
			Return[sc], (* Else *)
			sc = vor[tet];
			Return[sc]
			],
		If[ issmall[tet],
			Return[scoresmall[tet, isgroup]]
			]
		];
	sc = vor[tet];
	Return[sc]
	]
*/


double genscore( double tet[7], int isgroup )
{
	double sc;
	
	if( istet(tet) ) {
		if( isqrtet( tet ) ) {
			if( circumradius( tet ) <= 1.41 ) {
				sc = gma( tet );
				return( sc );
				}
			else {
				sc = vor( tet );
				return( sc );
				}
			}
		else {
			if( issmall( tet ) ) {
				sc = scoresmall( tet, isgroup );
				return( sc );
				}
			else {
				sc = vor( tet );
				return( sc );
				}
			}
	} else {
		return( -1.0 );
		}
}


double score( double tet[7] )
{
	return( genscore( tet, 0 ) );
}


double gscore( double tet[7] )
{
	return( genscore( tet, 1 ) );
}


double scorept( double tet[7] )
{
	return( score( tet )/POINT );
}


double gscorept( double tet[7] )
{
	return( gscore( tet )/POINT );
}


/*
crossdiag[{e1_,e2_,e3_,e4_,e5_,e6_,ep3_,ep4_,ep5_}]:=
Module[{i, x, y, z, diag, vect},
	x[1] = e1;
	x[2] = (e1^2 + e2^2 - e6^2)/(2 e1);
	x[3] = (e1^2 + e3^2 - e5^2)/(2 e1);
	y[2] = Sqrt[e2^2 - x[2]^2];
	y[3] = ((e2^2 + e3^2 - e4^2)/2 - x[2] x[3])/y[2];
	z[3] = Sqrt[e3^2 - x[3]^2 - y[3]^2];
	(* replace x[3],y[3],z[3] by x[4],y[4],z[4], 
		e3, e4, e5 by ep3, ep4, ep5, and change
		the sign of z. *)
	x[4] = (e1^2 + ep3^2 - ep5^2)/(2 e1);
	y[4] = ((e2^2 + ep3^2 - ep4^2)/2 - x[2] x[4])/y[2];
	z[4] = - Sqrt[ep3^2 - x[4]^2 - y[4]^2];
	vect = {x[3] - x[4], y[3] - y[4], z[3] - z[4]};
	diag = Sqrt[vect . vect];
	Return[diag]
	]
*/


double crossdiag( double len[10] )
{
  double x1, x2, x3, y2, y3, z3, x4, y4, z4, diag, vect[3];
  int i;
  
  x1 = len[1];
  x2 = (len[1]*len[1] + len[2]*len[2] - len[6]*len[6])/(2*len[1]);
  x3 = (len[1]*len[1] + len[3]*len[3] - len[5]*len[5])/(2*len[1]);
  y2 = sqrt( len[2]*len[2] - x2*x2 );
  y3 = -x2*x3/y2 + (len[2]*len[2] + len[3]*len[3] - len[4]*len[4])/(2*y2);
  z3 = sqrt( len[3]*len[3] - x3*x3 - y3*y3 );
  
  x4 = (len[1]*len[1] + len[7]*len[7] - len[9]*len[9])/(2*len[1]);
  y4 = -x2*x4/y2 + (len[2]*len[2] + len[7]*len[7] - len[8]*len[8])/(2*y2);
  z4 = -sqrt( len[7]*len[7] - x4*x4 - y4*y4 );
  diag = 0.0;
  vect[0] = x3 - x4;
  vect[1] = y3 - y4;
  vect[2] = z3 - z4;
  for( i=0; i<3; i++ )
  	diag += vect[i]*vect[i];
  return( sqrt( diag ) );
}


/*
crossokay[{e1_,e2_,e3_,e4_,e5_,e6_,ep3_,ep4_,ep5_}]:=
Module[{diag, face, crad2, val},
	diag = crossdiag[{e1,e2,e3,e4,e5,e6,ep3,ep4,ep5}];
	face = {e3, ep3, diag};
	crad2 = crad3len[face]^2;
	If[ crad2 > 2.0,
		val = True,
		val = False];
	Return[val]
	]
*/


int iscrossok( double len[10] )
{
	double diag, face[4], crad, crad2;
	int val;
	
	diag = crossdiag( len );
	if( diag < 2.51 )
		return( 0 );
	face[1] = len[3];
	face[2] = len[7];
	face[3] = diag;
	crad = crad3len( face );
	crad2 = crad*crad;
	if( crad2 > 2.0 )
		val = 1;
	else
		val = 0;
	return( val );
}


/*
anglesum[{e1_, e2_, ep2_, e3_, ep3_, e4_, ep4_, epp4_,
	eppp4_, e5_, ep5_, e6_, ep6_}]:=
Module[{tet1, tet2, tet3, tet4, sum},
	tet1 = {e1, e2, e3, e4, e5, e6};
	tet2 = {e1, e2, ep3, ep4, ep5, e6};
	tet3 = {e1, ep2, ep3, epp4, ep5, ep6};
	tet4 = {e1, ep2, e3, eppp4, e5, e6};
	sum = ndih[tet1] + ndih[tet2] + ndih[tet3] +
		ndih[tet4];
	Return[sum]
	]
*/


double anglesum( double octa[13] )
{
	double e1, e2, e3, e4, e5, e6, ep2, ep3, ep4, ep5, ep6, epp4, eppp4;
	double tet1[7], tet2[7], tet3[7], tet4[7];
	double sum;
	
	/* This is kinda ugly, but what the hell.  */
	e1 = 		octa[0];
	e2 =		octa[1];
	ep2 = 	octa[2];
	e3 = 		octa[3];
	ep3 = 	octa[4];
	e4 = 		octa[5];
	ep4 = 	octa[6];
	epp4 = 	octa[7];
	eppp4 = octa[8];
	e5 = 		octa[9];
	ep5 = 	octa[10];
	e6 = 		octa[11];
	ep6 = 	octa[12];
	
	/* tet1 = {e1, e2, e3, e4, e5, e6}; */
	tet1[1] = e1;
	tet1[2] = e2;
	tet1[3] = e3;
	tet1[4] = e4;
	tet1[5] = e5;
	tet1[6] = e6;
	
	/* tet2 = {e1, e2, ep3, ep4, ep5, e6}; */
	tet2[1] = e1;
	tet2[2] = e2;
	tet2[3] = ep3;
	tet2[4] = ep4;
	tet2[5] = ep5;
	tet2[6] = e6;
	
	/* tet3 = {e1, ep2, ep3, epp4, ep5, ep6}; */
	tet3[1] = e1;
	tet3[2] = ep2;
	tet3[3] = ep3;
	tet3[4] = epp4;
	tet3[5] = ep5;
	tet3[6] = ep6;
	
	/* tet4 = {e1, ep2, e3, eppp4, e5, ep6}; */
	tet4[1] = e1;
	tet4[2] = ep2;
	tet4[3] = e3;
	tet4[4] = eppp4;
	tet4[5] = e5;
	tet4[6] = ep6;
	
	sum = dih( tet1 ) + dih( tet2 ) + dih( tet3 ) + dih( tet4 );
	return( sum );
}


/*
seciter[x0_,x1_,roct__]:=
Module[{oct0, oct1, y0, y1, m, x2},
	oct0 = Join[{x0}, roct];
	oct1 = Join[{x1}, roct];
	y0 = anglesum[oct0];
	y1 = anglesum[oct1];
	m = (x1 - x0)/(y1 - y0);
	x2 = x0  - (y0 - twopi) m;
	Return[x2]
	]
*/


double seciter( double x0, double x1, double roct[13] )
{
	double y0, y1, m, x2, octa[13];
	int i;
	
	for( i=1; i<13; i++ )
		octa[i] = roct[i];
	octa[0] = x0;
	y0 = anglesum( octa );
	octa[0] = x1;
	y1 = anglesum( octa );
	m = (x1 - x0)/(y1 - y0);
	x2 = x0 - (y0 - TWOPI)*m;
	return( x2 );
}


/*
secdiag[{e2_, ep2_, e3_, ep3_, e4_, ep4_, epp4_,
	eppp4_, e5_, ep5_, e6_, ep6_}]:=
Module[{roct, eps, eps2, x0, x1, x2, done, sol, num},
	roct = {e2, ep2, e3, ep3, e4, ep4, epp4,
		eppp4, e5, ep5, e6, ep6};
	eps = 1.0 10^(-4);
	eps2 = 1.0 10^(-8);
	x0 = twosqrt2;
	x1 = twosqrt2 + eps;
	num = 0;
	done = False;
	While[ !done,
		x2 = seciter[x0, x1, roct];
		x0 = x1;
		x1 = x2;
		diff = Abs[x0 - x1];
		num = num + 1;
		Print["x1 = ", x1];
		If[ diff < eps2 || num > 20, done = True]
		];
	sol = Abs[ anglesum[Join[{x1}, roct]] - twopi ];
	If[ sol < eps2,
		Return[x1],
		Return[-1]]
	]
*/


double octadiag( double roct[13] )
{
	double octa[13], eps, eps2, x0, x1, x2, diff, sol, val;
	int i, num, notdone;
	
	eps = 1.0e-4;
	eps2 = 1.0e-8;
	num = 0;
	notdone = 1;
	x0 = TWOSQRT2;
	x1 = TWOSQRT2 + eps;
	
	for( i=1; i<13; i++ )
		octa[i] = roct[i];
	
	while( notdone ) {
		x2 = seciter( x0, x1, roct );
		x0 = x1;
		x1 = x2;
		diff = fabs( x0 - x1 );
		num += 1;
		if( diff < eps2 || num > 20 )
			notdone = 0;
		octa[0] = x1;
		sol = fabs( anglesum( octa ) - TWOPI );
		if( sol < eps2 )
			val = x1;
		else
			val = -100.0;
		}
		return( val );
}


/*
octacheck[{e1_, e2_, ep2_, e3_, ep3_, e4_, ep4_, 
	epp4_, eppp4_, e5_, ep5_, e6_, ep6_}]:=
Module[{c1, c2, val},
	c1 = niscrossok[{e5, e3, e4, e2, e6, e1, 
		eppp4, ep2, ep6}];
	c2 = niscrossok[{ep6, ep2, eppp4, e3, e5, e1, 
		epp4, ep3, ep5}];
	If[ c1 + c2 == 2,
		val = True,
		val = False];
	Return[val]
	]
*/


int isoctaok( double octa[13] )
{
	double len[10];
	double e1, e2, e3, e4, e5, e6, ep2, ep3, ep4, ep5, ep6, epp4, eppp4;
	
	/* This is kinda ugly, but what the hell.  */
	e1 = 		octa[0];
	e2 =		octa[1];
	ep2 = 	octa[2];
	e3 = 		octa[3];
	ep3 = 	octa[4];
	e4 = 		octa[5];
	ep4 = 	octa[6];
	epp4 = 	octa[7];
	eppp4 = octa[8];
	e5 = 		octa[9];
	ep5 = 	octa[10];
	e6 = 		octa[11];
	ep6 = 	octa[12];

	if( e1 < 2.51 )
		return( 0 );
		
	/* {e2, e6, e4, e5, e3, e1, ep4, ep5, ep3} */	
	len[1] = e2;
	len[2] = e6;
	len[3] = e4;
	len[4] = e5;
	len[5] = e3;
	len[6] = e1;
	len[7] = ep4;
	len[8] = ep5;
	len[9] = ep3;
	
	if( !isquadok( len ) )
		return( 0 );
	
	/* {e3, e5, e4, e6, e2, e1, eppp4, ep6, ep2} */
	len[1] = e3;
	len[2] = e5;
	len[3] = e4;
	len[4] = e6;
	len[5] = e2;
	len[6] = e1;
	len[7] = eppp4;
	len[8] = ep6;
	len[9] = ep2;
	
	if( !isquadok( len ) )
		return( 0 );
	
/* {ep3, ep5, ep4, e6, e2, e1, epp4, ep6, ep2} */
	len[1] = ep3;
	len[2] = ep5;
	len[3] = ep4;
	len[4] = e6;
	len[5] = e2;
	len[6] = e1;
	len[7] = epp4;
	len[8] = ep6;
	len[9] = ep2;
	
	if( !isquadok( len ) )
		return( 0 );
	
/* {ep2, ep6, epp4, ep5, ep3, e1, eppp4, e5, e3} */		
	len[1] = ep2;
	len[2] = ep6;
	len[3] = epp4;
	len[4] = ep5;
	len[5] = ep3;
	len[6] = e1;
	len[7] = eppp4;
	len[8] = e5;
	len[9] = e3;
	
	if( !isquadok( len ) )
		return( 0 );
	
	return( 1 );
}


/*
octascore[{e1_, e2_, ep2_, e3_, ep3_, e4_, ep4_, 
	epp4_, eppp4_, e5_, ep5_, e6_, ep6_}]:=
Module[{tet1, tet2, tet3, tet4, dsum},
	tet1 = {e1, e2, e3, e4, e5, e6};
	tet2 = {e1, e2, ep3, ep4, ep5, e6};
	tet3 = {e1, ep2, ep3, epp4, ep5, ep6};
	tet4 = {e1, ep2, e3, eppp4, e5, ep6};
	If[ issmall[tet1] && issmall[tet2] &&
		issmall[tet3] && issmall[tet4],
		dsum = gscore[tet1] + gscore[tet2] + 
			gscore[tet3] + gscore[tet4],
		dsum = score[tet1] + score[tet2] + 
			score[tet3] + score[tet4]
		];
	Return[dsum]
	]
*/


double octascore( double octa[13] )
{
	double e1, e2, e3, e4, e5, e6, ep2, ep3, ep4, ep5, ep6, epp4, eppp4;
	double tet1[7], tet2[7], tet3[7], tet4[7];
	double sum;
	
	/* This is kinda ugly, but what the hell.  */
	e1 = 		octa[0];
	e2 =		octa[1];
	ep2 = 	octa[2];
	e3 = 		octa[3];
	ep3 = 	octa[4];
	e4 = 		octa[5];
	ep4 = 	octa[6];
	epp4 = 	octa[7];
	eppp4 = octa[8];
	e5 = 		octa[9];
	ep5 = 	octa[10];
	e6 = 		octa[11];
	ep6 = 	octa[12];
	
	/* tet1 = {e1, e2, e3, e4, e5, e6}; */
	tet1[1] = e1;
	tet1[2] = e2;
	tet1[3] = e3;
	tet1[4] = e4;
	tet1[5] = e5;
	tet1[6] = e6;
	
	/* tet2 = {e1, e2, ep3, ep4, ep5, e6}; */
	tet2[1] = e1;
	tet2[2] = e2;
	tet2[3] = ep3;
	tet2[4] = ep4;
	tet2[5] = ep5;
	tet2[6] = e6;
	
	/* tet3 = {e1, ep2, ep3, epp4, ep5, ep6}; */
	tet3[1] = e1;
	tet3[2] = ep2;
	tet3[3] = ep3;
	tet3[4] = epp4;
	tet3[5] = ep5;
	tet3[6] = ep6;
	
	/* tet4 = {e1, ep2, e3, eppp4, e5, ep6}; */
	tet4[1] = e1;
	tet4[2] = ep2;
	tet4[3] = e3;
	tet4[4] = eppp4;
	tet4[5] = e5;
	tet4[6] = ep6;
	
	if( issmall( tet1 ) && issmall( tet2 ) && issmall( tet3 ) && issmall( tet4 ) )
		sum = gscore(tet1) + gscore(tet2) + gscore(tet3) + gscore(tet4);
	else
		sum = score(tet1) + score(tet2) + score(tet3) + score(tet4);
		
	return( sum );
}


/*
octasph[{e1_, e2_, ep2_, e3_, ep3_, e4_, ep4_, 
	epp4_, eppp4_, e5_, ep5_, e6_, ep6_}]:=
Module[{tet1, tet2, tet3, tet4, dsum},
	tet1 = {e1, e2, e3, e4, e5, e6};
	tet2 = {e1, e2, ep3, ep4, ep5, e6};
	tet3 = {e1, ep2, ep3, epp4, ep5, ep6};
	tet4 = {e1, ep2, e3, eppp4, e5, ep6};
	dsum = dbvol[tet1] + dbvol[tet2] + 
		dbvol[tet3] + dbvol[tet4];
	Return[3 dsum]
	]
*/


double octasph( double octa[13] )
{
	double e1, e2, e3, e4, e5, e6, ep2, ep3, ep4, ep5, ep6, epp4, eppp4;
	double tet1[7], tet2[7], tet3[7], tet4[7];
	double sum;
	
	/* This is kinda ugly, but what the hell.  */
	e1 = 		octa[0];
	e2 =		octa[1];
	ep2 = 	octa[2];
	e3 = 		octa[3];
	ep3 = 	octa[4];
	e4 = 		octa[5];
	ep4 = 	octa[6];
	epp4 = 	octa[7];
	eppp4 = octa[8];
	e5 = 		octa[9];
	ep5 = 	octa[10];
	e6 = 		octa[11];
	ep6 = 	octa[12];
	
	/* tet1 = {e1, e2, e3, e4, e5, e6}; */
	tet1[1] = e1;
	tet1[2] = e2;
	tet1[3] = e3;
	tet1[4] = e4;
	tet1[5] = e5;
	tet1[6] = e6;
	
	/* tet2 = {e1, e2, ep3, ep4, ep5, e6}; */
	tet2[1] = e1;
	tet2[2] = e2;
	tet2[3] = ep3;
	tet2[4] = ep4;
	tet2[5] = ep5;
	tet2[6] = e6;
	
	/* tet3 = {e1, ep2, ep3, epp4, ep5, ep6}; */
	tet3[1] = e1;
	tet3[2] = ep2;
	tet3[3] = ep3;
	tet3[4] = epp4;
	tet3[5] = ep5;
	tet3[6] = ep6;
	
	/* tet4 = {e1, ep2, e3, eppp4, e5, ep6}; */
	tet4[1] = e1;
	tet4[2] = ep2;
	tet4[3] = e3;
	tet4[4] = eppp4;
	tet4[5] = e5;
	tet4[6] = ep6;
	
	sum = dbvol(tet1) + dbvol(tet2) + dbvol(tet3) + dbvol(tet4);
		
	return( 3*sum );
}


void octapair( double octa[13], double xy[2] )
{
	double e1, e2, e3, e4, e5, e6, ep2, ep3, ep4, ep5, ep6, epp4, eppp4;
	double tet1[7], tet2[7], tet3[7], tet4[7];
	double sum;
	
	/* This is kinda ugly, but what the hell.  */
	e1 = 		octa[0];
	e2 =		octa[1];
	ep2 = 	octa[2];
	e3 = 		octa[3];
	ep3 = 	octa[4];
	e4 = 		octa[5];
	ep4 = 	octa[6];
	epp4 = 	octa[7];
	eppp4 = octa[8];
	e5 = 		octa[9];
	ep5 = 	octa[10];
	e6 = 		octa[11];
	ep6 = 	octa[12];
	
	/* tet1 = {e1, e2, e3, e4, e5, e6}; */
	tet1[1] = e1;
	tet1[2] = e2;
	tet1[3] = e3;
	tet1[4] = e4;
	tet1[5] = e5;
	tet1[6] = e6;
	
	/* tet2 = {e1, e2, ep3, ep4, ep5, e6}; */
	tet2[1] = e1;
	tet2[2] = e2;
	tet2[3] = ep3;
	tet2[4] = ep4;
	tet2[5] = ep5;
	tet2[6] = e6;
	
	/* tet3 = {e1, ep2, ep3, epp4, ep5, ep6}; */
	tet3[1] = e1;
	tet3[2] = ep2;
	tet3[3] = ep3;
	tet3[4] = epp4;
	tet3[5] = ep5;
	tet3[6] = ep6;
	
	/* tet4 = {e1, ep2, e3, eppp4, e5, ep6}; */
	tet4[1] = e1;
	tet4[2] = ep2;
	tet4[3] = e3;
	tet4[4] = eppp4;
	tet4[5] = e5;
	tet4[6] = ep6;

	if( issmall( tet1 ) && issmall( tet2 ) && issmall( tet3 ) && issmall( tet4 ) )
		sum = gscore(tet1) + gscore(tet2) + gscore(tet3) + gscore(tet4);
	else
		sum = score(tet1) + score(tet2) + score(tet3) + score(tet4);
	xy[1] = sum;	
	
	sum = dbvol(tet1) + dbvol(tet2) + dbvol(tet3) + dbvol(tet4);
	xy[0] = 3*sum;
}


/*  This routine checks to see if the given decomposition of a quadrilateral
actually squares with the Delaunay decomposition.  It checks to see if the
opposite vertex is closer to the circumcenter than the rest of the vertices
or not. */

int isquadok( double len[10] )
{
  double x1, x2, x3, y2, y3, z3, x4, y4, z4;
  double vect3p[4], vect[4][4], cent[4], dist1, dist2, sum, diff;
  int i, val;
  
  x1 = len[1];
  x2 = (len[1]*len[1] + len[2]*len[2] - len[6]*len[6])/(2*len[1]);
  x3 = (len[1]*len[1] + len[3]*len[3] - len[5]*len[5])/(2*len[1]);
  y2 = sqrt( len[2]*len[2] - x2*x2 );
  y3 = -x2*x3/y2 + (len[2]*len[2] + len[3]*len[3] - len[4]*len[4])/(2*y2);
  z3 = sqrt( len[3]*len[3] - x3*x3 - y3*y3 );
  
  x4 = (len[1]*len[1] + len[7]*len[7] - len[9]*len[9])/(2*len[1]);
  y4 = -x2*x4/y2 + (len[2]*len[2] + len[7]*len[7] - len[8]*len[8])/(2*y2);
  z4 = -sqrt( len[7]*len[7] - x4*x4 - y4*y4 );

	vect[1][1] = x1;
	vect[1][2] = 0.0;
	vect[1][3] = 0.0;
	
	vect[2][1] = x2;
	vect[2][2] = y2;
	vect[2][3] = 0.0;
	
	vect[3][1] = x3;
	vect[3][2] = y3;
	vect[3][3] = z3;
	
	vect3p[1] = x4;
	vect3p[2] = y4;
	vect3p[3] = z4;
	
	val = ccenter( vect, cent );
	if( val == 0 )
		return( 0 );	/* Got an error in ccenter() */
	
	sum = 0.0;
	for( i=1; i<4; i++ ) {
		diff = vect[1][i] - cent[i];
		sum += diff*diff;
	}
	dist1 = sum;
	
	sum = 0.0;
	for( i=1; i<4; i++ ) {
		diff = vect3p[i] - cent[i];
		sum += diff*diff;
	}
	dist2 = sum;
	
#ifdef DEBUG
	printf("dist1 = %lf, dist2 = %lf\n", dist1, dist2);
#endif
	
	if( dist2 > dist1 )
		return( 1 );
	else
		return( 0 );
}



/* The following routines, penta and hcp, produce the local packings associated with
the pentahedral prism and hexagonal close packings. */

void sph( double vect[4], double rvect[4] )
{
  double rho, theta, phi, x, y, z, r;

  rho = vect[1];
  theta = vect[2];
  phi = vect[3];
  r = rho*sin(phi);
  x = r*cos(theta);
  y = r*sin(theta);
  z = rho*cos(phi);
  rvect[1] = x;
  rvect[2] = y;
  rvect[3] = z;
}


void penta( double **vects )
{
  double th, tab[13][4], temp[4], rtemp[4];
  int j, k;

  th = 2*asin( 1.0/sqrt(3.0) );
  tab[1][1] = 2.0;
  tab[1][2] = 0.0;
  tab[1][3] = 0.0;
  tab[12][1] = 2.0;
  tab[12][2] = 0.0;
  tab[12][3] = PI;
  for( k=0; k<5; k++ ) {
    tab[k+2][1] = 2.0;
    tab[k+2][2] = k*th;
    tab[k+2][3] = PI/3.0;
  }
  for( k=0; k<5; k++ ) {
    tab[k+7][1] = 2.0;
    tab[k+7][2] = k*th;
    tab[k+7][3] = 2*PI/3.0;
  }
  for( k=1; k<=12; k++ ) {
    for( j=1; j<4; j++ )
      temp[j] = tab[k][j];
    sph(temp, rtemp);
    for( j=1; j<4; j++ )
      vects[k][j] = rtemp[j];
  }
  vects[0][1] = 12.0;
}


void hcp( double **vects )
{
  double th, phi, tab[13][4], temp[4], rtemp[4];
  int j, k;

  th = 2.0*PI/3.0;
  phi = asin( 1.0/sqrt(3.0) );
  for( k=0; k<3; k++ ) {
    tab[k+1][1] = 2.0;
    tab[k+1][2] = PI/6.0 + k*th;
    tab[k+1][3] = phi;
  }
  for( k=0; k<3; k++ ) {
    tab[k+10][1] = 2.0;
    tab[k+10][2] = PI/2.0 + k*th;
    tab[k+10][3] = PI - phi;
  }
  for( k=0; k<6; k++ ) {
    tab[k+4][1] = 2.0;
    tab[k+4][2] = k*PI/3.0;
    tab[k+4][3] = PI/2.0;
  }
  for( k=1; k<=12; k++ ) {
    for( j=1; j<4; j++ )
      temp[j] = tab[k][j];
    sph(temp, rtemp);
    for( j=1; j<4; j++ )
      vects[k][j] = rtemp[j];
  }
  vects[0][1] = 12.0;
}

/*
double mrand( void )
{
  double phil;

  phil = 2.0*drand48() - 1.0;
  return(phil);
}
*/

/*
compfun[{e1_,e2_,e3_,e4_,e5_,e6_,ep3_,ep4_,ep5_}]:=
Module[{auxc1, auxc2, auxc3, auxc4, auxc5, auxc6,
	x1, x2, x3, x4, x5, x6, y3, y4, y5, tot, auxc32},
	x1 = e1^2;
	x2 = e2^2;
	x3 = e3^2;
	x4 = e4^2;
	x5 = e5^2;
	x6 = e6^2;
	y3 = ep3^2;
	y4 = ep4^2;
	y5 = ep5^2;
	auxc1 = -x1^2 + 2*x1*x2 - x2^2 + 2*x1*x6 + 
		2*x2*x6 - x6^2;
	auxc2 = x2 - (x1 + x2 - x6)^2/(4*x1);
	auxc3 = (x2 + y3 - y4)/2 - 
		((x1 + x2 - x6)*(x1 + y3 - y5))/(4*x1);
	auxc4 = -(x1 + y3 - y5)^2/(4*x1);
	auxc5 = -(x1^2*x4) + x1*x2*x4 + x1*x2*x5 - x2^2*x5 - 
		2*x1*x2*x6 + x1*x3*x6 + x2*x3*x6 + x1*x4*x6 + 
		x2*x5*x6 - x3*x6^2;
	auxc6 = x1^2*x4 - x1*x2*x4 - x1*x3*x4 + x2*x3*x4 + 
		x1*x4^2 - x1*x2*x5 + x2^2*x5 + x1*x3*x5 - 
		x2*x3*x5 - x1*x4*x5 - x2*x4*x5 + 
  		x2*x5^2 + x1*x2*x6 - x1*x3*x6 - x2*x3*x6 + 
  		x3^2*x6 - x1*x4*x6 - x3*x4*x6 - x2*x5*x6 - 
  		x3*x5*x6 + x4*x5*x6 + x3*x6^2;
  	auxc32 = auxc3*auxc3;
  	tot = -(auxc32/auxc2) + auxc4 + 
		(auxc3*(auxc3/Sqrt[auxc2] + 
       	(x1 - x2 - x6)/Sqrt[auxc1/x1]))/Sqrt[auxc2] + 
       	y3 - 
  		(auxc5*Sqrt[(-(auxc6/auxc1))]*
  		Sqrt[(-(auxc32/auxc2) + 
  		auxc4 + y3)])/auxc6 + 
  		((-x1 + y3 - y5)*(x1 + y3 - y5))/(4*x1);
  	Return[tot]
  	]
*/

#ifdef WANTOLD
double quaddiff( double quad[10] )
{
	double auxc1, auxc2, auxc3, auxc4, auxc5, auxc6;
	double x1, x2, x3, x4, x5, x6, y3, y4, y5, tot, auxc32;
	/*	double e1, e2, e3, e4, e5, e6, ep3, ep4, ep5;	*/
	double x12, x22, x32, x42, x52, x62;
	double temp;
	
	/*
	e1 = quad[1];
	e2 = quad[2];
	e3 = quad[3];
	e4 = quad[4];
	e5 = quad[5];
	e6 = quad[6];
	ep3 = quad[7];
	ep4 = quad[8];
	ep5 = quad[9];
	*/
	x1 = SQUARE(quad[1]);
	x2 = SQUARE(quad[2]);
	x3 = SQUARE(quad[3]);
	x4 = SQUARE(quad[4]);
	x5 = SQUARE(quad[5]);
	x6 = SQUARE(quad[6]);
	y3 = SQUARE(quad[7]);
	y4 = SQUARE(quad[8]);
	y5 = SQUARE(quad[9]);
	/*
	x1 = SQUARE(e1);
	x2 = SQUARE(e2);
	x3 = SQUARE(e3);
	x4 = SQUARE(e4);
	x5 = SQUARE(e5);
	x6 = SQUARE(e6);
	y3 = SQUARE(ep3);
	y4 = SQUARE(ep4);
	y5 = SQUARE(ep5);
	*/
	x12 = SQUARE( x1 );
	x22 = SQUARE( x2 );
	x32 = SQUARE( x3 );
	x42 = SQUARE( x4 );
	x52 = SQUARE( x5 );
	x62 = SQUARE( x6 );
	auxc1 = -x12 + 2*x1*x2 - x22 + 2*x1*x6 + 
		2*x2*x6 - x62;
	temp = (x1 + x2 - x6);
	auxc2 = x2 - temp*temp/(4*x1);
	auxc3 = (x2 + y3 - y4)/2 - 
		((x1 + x2 - x6)*(x1 + y3 - y5))/(4*x1);
	temp = (x1 + y3 - y5);
	auxc4 = -temp*temp/(4*x1);
	auxc5 = -(x12*x4) + x1*x2*x4 + x1*x2*x5 - x22*x5 - 
			2*x1*x2*x6 + x1*x3*x6 + x2*x3*x6 + x1*x4*x6 + 
			x2*x5*x6 - x3*x62;
	auxc6 = x12*x4 - x1*x2*x4 - x1*x3*x4 + x2*x3*x4 + 
			x1*x42 - x1*x2*x5 + x22*x5 + x1*x3*x5 - 
			x2*x3*x5 - x1*x4*x5 - x2*x4*x5 + 
  		x2*x52 + x1*x2*x6 - x1*x3*x6 - x2*x3*x6 + 
  		x32*x6 - x1*x4*x6 - x3*x4*x6 - x2*x5*x6 - 
  		x3*x5*x6 + x4*x5*x6 + x3*x62;
  auxc32 = SQUARE(auxc3);
  tot = -(auxc32/auxc2) + auxc4 + 
				(auxc3*(auxc3/sqrt(auxc2) + 
       	(x1 - x2 - x6)/sqrt(auxc1/x1)))/sqrt(auxc2) + y3 - 
  			(auxc5*sqrt((-(auxc6/auxc1)))*
  			sqrt((-(auxc32/auxc2) + 
  			auxc4 + y3)))/auxc6 + 
  			((-x1 + y3 - y5)*(x1 + y3 - y5))/(4*x1);
  return(tot);
}
#endif

