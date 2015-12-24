//  copyright (c) 1997, Thomas C. Hales, all rights reserved.

// This and seconds.cc give essential bounds on functions used
// in the Kepler conjecture.  This file contains the functions
// of three variables such as eta and quoin, and seconds.cc contains
// the others such as solid,dih,vor_analytic, etc.

// The output has been saved in several files out.quoin.* out.eta2.
// This will be fed into the inequality verifier to produce good
// Taylor approximations to our functions.

// This program only needs to be run once, which is fortunate, because
// these verifications took a few days each on UofM computers supernova
// and blackbox.

// 2nd derivative verifications of 3d guys:
#include <iomanip.h>
#include <iostream.h>
#include <math.h>
#include <stdlib.h>
#include "error.h"
#include "interval.h"
#include "io.h"


interval eta2xx(double x[3],double z[3])
	{
	down();
	double umin = (-z[0])*z[0] + 2.0*x[0]*x[1] +(-z[1])*z[1] +2.0*x[0]*x[2] 
		+2.0*x[1]*x[2] + (-z[2])*z[2];
	if (umin<0.0) umin=0.0;
	double umin3 = umin*umin*umin;
	cout << umin3 << endl;
	double nmin =  x[0]*x[0]*x[0]  
		+ 2.0*x[1]*x[1]*x[1] + 6.0*x[0]*x[1]*x[2] 
		+ 2.0*x[2]*x[2]*x[2]
		+3.0*(((-z[0])*z[1])*z[1]) +3.0*(((-z[0])*z[2])*z[2])
		+ 2.0*(((-z[2])*z[1])*z[1]) + 2.0*((-z[2]*z[1])*z[2]);
	cout << nmin << endl;
	double x12 = x[1]*x[2]*2.0;
	cout << x12 << endl;
 
	up();
	double umax = (-x[0])*x[0] + 2.0*z[0]*z[1] +(-x[1])*x[1] +2.0*z[0]*z[2] 
		+2.0*z[1]*z[2] + (-x[2])*x[2];
	if (umax<0.0) umax=0.0;
	double umax3 = umax*umax*umax;
	double nmax =  z[0]*z[0]*z[0] + 3.0*(((-x[0])*x[1])*x[1]) + 
		2.0*z[1]*z[1]*z[1] + 6.0*z[0]*z[1]*z[2] 
		+2.0*(((-x[1])*x[1])*x[2])
		 + 3.0*(((-x[0])*x[2])*x[2]) + 2.0*(((-x[1])*x[2])*x[2])
		+ 2.0*z[2]*z[2]*z[2];
	double z12 = z[1]*z[2]*2.0;
 
	interval n = interval(nmin,nmax);
	interval u3 = interval(umin3,umax3);
	return interval(x12,z12)*n/u3;
	}

// second derivative of eta^2, wrt x[0],x[1].
interval eta2xy(double x[3],double z[3])
	{
	up();
	double z0sq = z[0]*z[0], z1sq = z[1]*z[1], z2sq=z[2]*z[2];
 
	down();
        double umin = -z1sq-z2sq-z0sq 
		+2.0*(x[0]*x[1] +x[0]*x[2] +x[1]*x[2]);
        if (umin<0.0) umin=0.0;
        double umin3 = umin*umin*umin;
	double x0sq = x[0]*x[0], x1sq = x[1]*x[1], x2sq=x[2]*x[2];
        double nmin =    
		+ 6.0*x0sq*x1sq
		+ 2.0*x[0]*x0sq*x[2] 
		+ 2.0*x[1]*x1sq*x[2] + 10.0*x[0]*x[1]*x2sq
		+ x2sq*x2sq+
		(-z0sq)*z0sq +2.0*((-z[0])*z0sq)*z[1] 
		 +2.0*((-z[0])*z2sq)*z[2]+  2.0*((-z[1])*z2sq)*z[2] +
		 2.0*((-z[0])*z1sq)*z[1] + (-z1sq)*z1sq +
		 6.0*((-z0sq)*z[1])*z[2] + 6.0*((-z[0])*z1sq)*z[2];
 
 
        up();
        double umax = -x0sq-x1sq-x2sq 
		+ 2.0*(z[0]*z[1] +z[0]*z[2] +z[1]*z[2] );
        if (umax<0.0) umax=0.0;
        double umax3 = umax*umax*umax;
        double nmax =  
		+ 6.0*z0sq*z1sq
		+ 2.0*z[0]*z0sq*z[2] 
		+ 2.0*z[1]*z1sq*z[2] + 10.0*z[0]*z[1]*z2sq
		+ z2sq*z2sq +
		(-x0sq)*x0sq +2.0*((-x[0])*x0sq)*x[1] 
		 +2.0*((-x[0])*x2sq)*x[2]  + 2.0*((-x[1])*x2sq)*x[2] +
		 2.0*((-x[0])*x1sq)*x[1] + (-x1sq)*x1sq +
		 6.0*((-x0sq)*x[1])*x[2] + 6.0*((-x[0])*x1sq)*x[2];
 
	return interval(nmin,nmax)*interval(x[2],z[2])/interval(umin3,umax3);
	}
 

void runeta2()
	{
	int i,j,k,q;
	double x[3],z[3];
	interval tempxx = zero;
	interval tempxy = zero;
	double eminxx=10000.0, emaxxx = -100000.0;
	double eminxy=10000.0, emaxxy = -100000.0;
	interval w0 = interval("2.3001")/interval("100.0"); 
	interval w1 = interval("1.6999")/interval("100.0"); 
	interval v[3]; 
	interval x0[3]; 
	char* tag;
	
	int top=100;
 
	for (q=0;q<3;q++)
	{
	x0[0]=four; x0[1]=four; x0[2]=four;
	v[0]=v[1]=v[2]=w0;
	switch(q) {
		case 0 :  // quasi-regular;
			tag = "QUASI-REG : ";
			break;
		case 1 : x0[0]=interval("6.3001"); v[0]=w1;
			tag = "upright X-long";
			break;
		case 2 : x0[2]=interval("6.3001");  v[2]=w1;
			tag = "upright Z-long";
			break;
		  }
	for (i=0;i<top;i++) 
	for (j=0;j<top;j++)
	for (k=0;k<top;k++)
		{ 
		down();
		x[0]= inf(x0[0]) + inf(v[0])*double(i);
		x[1]= inf(x0[1]) + inf(v[1])*double(j);
		x[2]= inf(x0[2]) + inf(v[2])*double(k);
		up();
		z[0]= sup(x0[0]) + sup(v[0])*double(i+1);
		z[1]= sup(x0[1]) + sup(v[1])*double(j+1);
		z[2]= sup(x0[2]) + sup(v[2])*double(k+1);
		tempxx = eta2xx(x,z);
		if (emaxxx< sup(tempxx)) emaxxx = sup(tempxx);
		if (eminxx> inf(tempxx)) eminxx = inf(tempxx);
		tempxy = eta2xy(x,z);
		if (emaxxy< sup(tempxy)) emaxxy = sup(tempxy);
		if (eminxy> inf(tempxy)) eminxy = inf(tempxy);
		}
	cout.precision(20);
	cout << tag << endl;
	cout << "XX: " << eminxx << " " << emaxxx << endl;
	cout << "XY: " << eminxy << " " << emaxxy << endl;
	}
	cout << "*" << endl;
	}

void half_array_multiply(interval v,interval DDu[3][3],interval vDDu[3][3])
	{
	int i,j;
	for (i=0;i<3;i++) for (j=i;j<3;j++)
		vDDu[i][j] = v*DDu[i][j];
	}

void product3(interval& u,interval Du[3],interval DDu[3][3],
             interval& v,interval Dv[3],interval DDv[3][3],
             interval& uv,interval Duv[3],interval DDuv[3][3])
        {
        int i,j;
        uv = u*v;
        for (i=0;i<3;i++) Duv[i] = u*Dv[i] + v*Du[i];
        interval DuDv[3][3];
        for (i=0;i<3;i++) for (j=0;j<3;j++) DuDv[i][j]=Du[i]*Dv[j];
        interval uDDv[3][3],vDDu[3][3];
        half_array_multiply(v,DDu,vDDu);
        half_array_multiply(u,DDv,uDDv);
        for (i=0;i<3;i++) for (j=i;j<3;j++)
                DDuv[i][j]= uDDv[i][j]+DuDv[i][j]+DuDv[j][i]+vDDu[i][j];
        for (i=0;i<3;i++) for (j=0;j<i;j++)
                DDuv[i][j]= DDuv[j][i];
        }

void product2(interval& u,interval Du[2],interval DDu[2][2],
             interval& v,interval Dv[2],interval DDv[2][2],
             interval& uv,interval Duv[2],interval DDuv[2][2])
        {
        int i,j;
        uv = u*v;
        for (i=0;i<2;i++) Duv[i] = u*Dv[i] + v*Du[i];
        for (i=0;i<2;i++) for (j=i;j<2;j++)
                DDuv[i][j]= u*DDv[i][j]+Du[i]*Dv[j]+Du[j]*Dv[i]+v*DDu[i][j];
	DDuv[1][0]= DDuv[0][1];
        }

void Dsqrt3(interval&u,interval Du[3],interval DDu[3][3],
           interval& sqrt_u,interval Dsqrt_u[3],interval DDsqrt_u[3][3])
        {
        int i,j;
        sqrt_u = sqrt(u);
	if (inf(u)<=0.0) error_msg("derivative at zero encountered");
        interval t = interval(
		(down(), 0.5/sup(sqrt_u)),
		(up(), 0.5/inf(sqrt_u)))
;
 
        for (i=0;i<3;i++) Dsqrt_u[i]= Du[i]*t;
 
        interval v = interval(
		(down(), -0.5/inf(u)),
		(up(), -0.5/sup(u)));
 
        interval vDu[3];
        for (i=0;i<3;i++) vDu[i] = Du[i]*v;
 
        for (i=0;i<3;i++) for (j=i;j<3;j++)
                DDsqrt_u[i][j] =(Du[i]*vDu[j] + DDu[i][j])*t;
        for (i=0;i<3;i++) for (j=0;j<i;j++)
                DDsqrt_u[i][j] = DDsqrt_u[j][i];
        }

void Dfivehalves(interval& u,interval Du[2],interval DDu[2][2],
	interval& fh,interval Dfh[2],interval DDfh[2][2])
	{
	int i,j;
	interval sqrt_u = sqrt(u);
	static interval fiveh = five/two;
	interval x = fiveh*u*sqrt_u;
	fh = sqrt_u*u*u;			// Def' is u^(5/2) & derivs.
	for (i=0;i<2;i++) Dfh[i]= x*Du[i];
	static interval fifq = three*five/four;
	interval y = fifq*sqrt_u;
	for (i=0;i<2;i++) for (j=i;j<2;j++)
		DDfh[i][j]= y*Du[i]*Du[j] + x*DDu[i][j];
	DDfh[1][0]= DDfh[0][1];
	}

void atan(interval&u,interval Du[3],interval DDu[3][3],
	interval& atanu,interval Datanu[3],interval DDatanu[3][3])
	{
	int i,j;
	atanu = atan(u);
	interval t = one/(one+u*u);
	interval t1 = -two*t*t*u;
	interval v[3] = {t1*Du[0],t1*Du[1],t1*Du[2]};
	for (i=0;i<3;i++) Datanu[i]= Du[i]*t;
	for (i=0;i<3;i++) for (j=i;j<3;j++)
		DDatanu[i][j]= DDu[i][j]*t + t1*Du[i]*Du[j];
	for (i=0;i<3;i++) for (j=0;j<i;j++)
		DDatanu[i][j]= DDatanu[j][i];
	}

// modified arctan is defined as matan(x) = atan(sqrt(x))-sqrt(x)+sqrt(x)^3/3.
void matan(interval&u,interval Du[2],interval DDu[2][2],
	interval& ma,interval Dma[2],interval DDma[2][2])
	{
	int i,j;
	interval t = sqrt(u);
	interval ru1=one/(two*(one+u));
	static interval third= one/three;
	// make use of monotonicity of matan in computation of its interval:
	interval atant = atan(t);
	double mau;
	up(); mau = sup(atant) + (-sup(t)) + sup(t)*sup(u)*sup(third);
	//                         ^^ sup instead of inf by monotonicity.
	double man;
	down(); man = inf(atant) + (-inf(t)) + inf(t)*inf(u)*inf(third);
	ma = interval(man,mau); 
	interval x = t*u*ru1;
	for (i=0;i<2;i++) Dma[i]=Du[i]*x;
	interval y = t*(three+u)*ru1*ru1;
	interval Df[2]= {Du[0]*y,Du[1]*y};
	for (i=0;i<2;i++) for (j=i;j<2;j++)
		DDma[i][j]= DDu[i][j]*x + Df[i]*Du[j];
	DDma[1][0]= DDma[0][1];
	}


void quotient3(interval& a,interval Da[3],interval DDa[3][3],
                interval& b,interval Db[3],interval DDb[3][3],
                interval& v,interval Dv[3],interval DDv[3][3])
        {
        int i,j;
        v = a/b;
        interval b2 = one/(b*b);
        interval t = -two*b2/b;
        interval r[3];
        for (i=0;i<3;i++) r[i] = (Da[i]*b- a*Db[i]);
        for (i=0;i<3;i++) Dv[i] = r[i]*b2;
        for (i=0;i<3;i++) r[i] = r[i]*t;
        interval DaDb[3][3];
        for (i=0;i<3;i++) for (j=0;j<3;j++) DaDb[i][j]=Da[i]*Db[j];
 
        interval bDDa[3][3]; half_array_multiply(b,DDa,bDDa);
        interval aDDb[3][3]; half_array_multiply(a,DDb,aDDb);
 
        for (i=0;i<3;i++) for (j=i;j<3;j++)
                DDv[i][j]= (bDDa[i][j]+DaDb[i][j]-DaDb[j][i]-aDDb[i][j]);
        half_array_multiply(b2,DDv,DDv);
        for (i=0;i<3;i++) for (j=i;j<3;j++)
                DDv[i][j] = DDv[i][j] + r[i]*Db[j];
        for (i=0;i<3;i++) for (j=0;j<i;j++)
                DDv[i][j]= DDv[j][i];
        }

void quotient2(interval& a,interval Da[2],interval DDa[2][2],
                interval& b,interval Db[2],interval DDb[2][2],
                interval& v,interval Dv[2],interval DDv[2][2])
        {
        int i,j;
        v = a/b;
        interval b2 = one/(b*b);
        interval t = -two*b2/b;
        interval r[2];
        for (i=0;i<2;i++) r[i] = (Da[i]*b- a*Db[i]);
        for (i=0;i<2;i++) Dv[i] = r[i]*b2;
        for (i=0;i<2;i++) r[i] = r[i]*t;
 
        for (i=0;i<2;i++) for (j=i;j<2;j++)
                DDv[i][j]= (b*DDa[i][j]+Da[i]*Db[j]-Da[j]*Db[i]-a*DDb[i][j]);
        for (i=0;i<2;i++) for (j=i;j<2;j++)
                DDv[i][j] = b2*DDv[i][j] + r[i]*Db[j];
	DDv[1][0]= DDv[0][1];
        }


void output3(interval d,interval Dd[3],interval DDd[3][3])
        {
        int i,j;
        cout << " ----- \n" << band(d) << endl;
        for (i=0;i<3;i++) cout << i << " " << band(Dd[i]) << endl;
        for (i=0;i<3;i++) for (j=0;j<3;j++) 
		cout << i << " " << j << " " << band(DDd[i][j]) << endl;
	cout << flush;
        }

void output2(interval d,interval Dd[2],interval DDd[2][2])
        {
        int i,j;
        cout << " ----- \n" << band(d) << endl;
        for (i=0;i<2;i++) cout << i << " " << band(Dd[i]) << endl;
        for (i=0;i<2;i++) for (j=0;j<2;j++) 
		cout << i << " " << j << " " << band(DDd[i][j]) << endl;
	cout << flush;
        }

void compute_rogersterm(interval x,interval y,interval z,
	interval& f,interval Df[2],interval DDf[2][2])
	{// the nasty polynomial that shows up: Rogers5d
	double xxmin[2][2],xxmax[2][2];
	double  xn=inf(x), yn=inf(y), zn=inf(z), 
		xu=sup(x), yu=sup(y), zu=sup(z);
	up();
	double xu2=xu*xu, yu2=yu*yu, zu2=zu*zu;
	double xu3=xu*xu2, yu3=yu*yu2, zu3=zu*zu2;
	down();
	double xn2=xn*xn, yn2=yn*yn, zn2=zn*zn;
	double xn3=xn2*xn, yn3=yn2*yn, zn3=zn2*zn;
 
	double fmin,fmax;
	down();
	fmin = //=Rogers5d
		2.0*xn2*xn*yn2 + 3.0*(xu*(yu2*(-yu2))) + 8.0*xn2*xn*yn*zn 
		+ 12.0*(xu*(yu2*(yu*(-zu)))) + 8.0*xn2*xn*zn2 
		+ 12.0*(xu*(yu2*(-zu2))) + 
		6.0*xn2*zn3 + 12.0*(xu*(yu*(-zu3))) + 8.0*yn2*zn3 
		+ 3.0*(xu*(zu2*(-zu2))) + 8.0*yn*zn2*zn2 + 
		   2.0*zn2*zn2*zn;
 
	up();
	fmax = //=Rogers5d
		2.0*xu2*xu*yu2 + 3.0*(xn*(yn2*(-yn2))) + 8.0*xu2*xu*yu*zu 
		+ 12.0*(xn*(yn2*(yn*(-zn)))) + 8.0*xu2*xu*zu2 
		+ 12.0*(xn*(yn2*(-zn2))) + 
		6.0*xu2*zu3 + 12.0*(xn*(yn*(-zn3))) + 8.0*yu2*zu3 
		+ 3.0*(xn*(zn2*(-zn2))) + 8.0*yu*zu2*zu2 + 
		   2.0*zu2*zu2*zu;
	f = interval(fmin,fmax);
 
	// Warning!! We never use the z-derivatives, so we set them to zero.
	down(); 
	xxmin[0][0] = 12.0*(xn*yn2 + 4.0*xn*zn*(yn+zn) + zn3);
	xxmin[0][1]=xxmin[1][0]= 12.0*(xn2*yn -yu3 + 
		2.0*xn2*zn + 3.0*(yu2*(-zu)) + 2.0*(yu*(-zu2)) - zu3);
	xxmin[1][1]=4.0*(xn3+9.0*(xu*(-yu2))+18.0*(xu*(yu*(-zu))) 
		+ 6.0*(xu*(-zu2))
		+4.0*zn3);
	up();
	xxmax[0][0] = 12.0*(xu*yu2 + 4.0*xu*zu*(yu+zu) + zu3);
	xxmax[0][1]=xxmax[1][0]= 12.0*(xu2*yu -yn3 + 
		2.0*xu2*zu + 3.0*(yn2*(-zn)) + 2.0*(yn*(-zn2)) - zn3);
	xxmax[1][1]=4.0*(xu3+9.0*(xn*(-yn2))+18.0*(xn*(yn*(-zn))) 
		+ 6.0*(xn*(-zn2))
		+4.0*zu3);
	int i,j; for (i=0;i<2;i++) for (j=0;j<2;j++)
		DDf[i][j]=interval(xxmin[i][j],xxmax[i][j]);
 
	// now for the first derivatives:
	double xmin[2],xmax[2];
 
	down();
	xmin[0]= 3.0*(2.0*xn2*yn2 +yu2*(-yu2)+8.0*xn2*yn*zn+4.0*(yu3*(-zu))
		+8.0*xn2*zn2 + 4.0*(yu2*(-zu2))+4.0*xn*zn3 + 4.0*(yu*(-zu3))
		+zu2*(-zu2));
	xmin[1]= 4.0*(xn3*yn + 3.0*(xu*(-yu3))+2.0*xn3*zn +9.0*(xu*(yu2*(-zu)))
		+6.0*(xu*(yu*(-zu2))) + 3.0*(xu*(-zu3))+4.0*yn*zn3+2.0*zn2*zn2);
 
	up();
	xmax[0]= 3.0*(2.0*xu2*yu2 +yn2*(-yn2)+8.0*xu2*yu*zu +4.0*(yn3*(-zn))
		+8.0*xu2*zu2 + 4.0*(yn2*(-zn2))+4.0*xu*zu3 + 4.0*(yn*(-zn3))
		+(zn2*(-zn2)));
	xmax[1]= 4.0*(xu3*yu + 3.0*(xn*(-yn3))+2.0*xu3*zu +9.0*(xn*(yn2*(-zn)))
		+6.0*(xn*(yn*(-zn2))) + 3.0*(xn*(-zn3))+4.0*yu*zu3+2.0*zu2*zu2);
	for (i=0;i<2;i++) Df[i]=interval(xmin[i],xmax[i]);
	}


// x < y < z : rogers simplex variables.
// Warning!!, we don't compute z derivatives, they are not needed (z constant).
void compute_quoin(interval x,interval y,interval z,interval& f,
		interval Df[2],interval DDf[2][2]) 
	// method two : tangent expansion.
	{
	int i,j;
	interval x2 = x*x;
	interval y2 = y*y;
	interval z2 = z*z;
	interval DDf1[2][2];
	DDf1[0][0] = -six*x;
	DDf1[0][1] = DDf1[1][0]= zero;
	DDf1[1][1] = zero;
	interval Df1[2];
	Df1[0]= three*(-x*x+z*z);
	Df1[1]= zero;
	interval f1 = (-x+z)*(x*x+x*z-two*z*z);
 
	interval unum = z2-y2;
	interval Dunum[2] = {zero,-two*y};
	interval DDunum[2][2]= {{zero,zero},
			  	{zero,-two}};
	interval uden = y2-x2;
	interval Duden[2]= {-two*x,two*y};
	interval DDuden[2][2] = {{-two,zero},
				{zero,two}};
 
	interval u2, Du2[2], DDu2[2][2];
	interval atu,Datu[2],DDatu[2][2];
	quotient2(unum,Dunum,DDunum,uden,Duden,DDuden,u2,Du2,DDu2);
		
	matan(u2,Du2,DDu2,atu,Datu,DDatu);
 
	interval g1,Dg1[2],DDg1[2][2];
	product2(atu,Datu,DDatu,f1,Df1,DDf1,g1,Dg1,DDg1);
 
 
	interval d5,Dd5[2],DDd5[2][2];
	compute_rogersterm(x,y,z,d5,Dd5,DDd5);
	interval fh,Dfh[2],DDfh[2][2];
	Dfivehalves(u2,Du2,DDu2,fh,Dfh,DDfh);
	//interval lnum,Dlnum[2],DDlnum[2][2]; // same as uden //
	interval lden;
	interval yz = y+z; interval yz2=yz*yz; interval yz3=yz*yz2;
	lden = yz2*yz2*three;
	interval ldenf = four*three*yz3;
	interval Dlden[2] = {zero,ldenf};
	interval ldenf2 = interval("36.0")*yz2;
	interval DDlden[2][2]= {{zero,zero},
			{zero,ldenf2}};
	interval quo,Dquo[2],DDquo[2][2];
	quotient2(uden,Duden,DDuden,lden,Dlden,DDlden,quo,Dquo,DDquo);
	interval xx,Dxx[2],DDxx[2][2];
	product2(quo,Dquo,DDquo,fh,Dfh,DDfh,xx,Dxx,DDxx);
	interval g2,Dg2[2],DDg2[2][2];
	product2(xx,Dxx,DDxx,d5,Dd5,DDd5,g2,Dg2,DDg2);
 
	// Term g3.
	interval fz3 = four*z*z2;
	interval vnum = (-x+y)*(-y+z);
	interval Dvnum[2]= {y-z,-two*y+z+x};
	interval DDvnum[2][2]= {{zero,one},{one,-two}};
 
	interval vden = (x+y)*(y+z);
	interval Dvden[2]= {y+z,two*y+z+x};
	interval DDvden[2][2] = {{zero,one},{one,two}};
	interval v,Dv[2],DDv[2][2];
	quotient2(vnum,Dvnum,DDvnum,vden,Dvden,DDvden,v,Dv,DDv);
	interval mv,Dmv[2],DDmv[2][2];
	matan(v,Dv,DDv,mv,Dmv,DDmv);
	interval g3,Dg3[2],DDg3[2][2];
	g3 = mv*fz3;
	for (i=0;i<2;i++) Dg3[i]=fz3*Dmv[i];
	for (i=0;i<2;i++) for (j=0;j<2;j++) DDg3[i][j]=fz3*DDmv[i][j];
 
	// Combine terms. m16 because we consistently left out a factor of -1/6.
	static interval m16 = -one/interval("6.0");
	f = m16*(g1+g2+g3);
	for (i=0;i<2;i++) Df[i] = m16*(Dg1[i]+Dg2[i]+Dg3[i]);
	for (i=0;i<2;i++) for (j=i;j<2;j++) DDf[i][j]=
			m16*(DDg1[i][j]+DDg2[i][j]+DDg3[i][j]);
	for (i=0;i<2;i++) for (j=0;j<i;j++)
			DDf[i][j]=DDf[j][i];
	}

void compute_chainquoin(interval trunc,
	interval x1,interval x2,interval x6,
		interval DDf[3][3]) 
	{
	interval a2 = x1/four;
	interval a = sqrt(a2);
	interval Da2[3] = {one/four,zero,zero};
	interval DDa2[3][3]= 
		{{zero,zero,zero},{zero,zero,zero},{zero,zero,zero}};
	interval Da[3],DDa[3][3]; Dsqrt3(a2,Da2,DDa2,a,Da,DDa);
 
	interval u126 = -x1*x1+two*x1*x2-x2*x2+two*x1*x6-x6*x6+
		two*x2*x6;
	interval Du126[3] = {two*(-x1+x2+x6),two*(x1-x2+x6),two*(x1+x2-x6)};
	interval DDu126[3][3]={{-two,two,two},{two,-two,two},{two,two,-two}};
	interval num = x1*x2*x6;
	interval Dnum[3] = {x2*x6,x1*x6,x1*x2};
	interval DDnum[3][3] = {{zero,x6,x2},{x6,zero,x1},{x2,x1,zero}};
	interval etas,Detas[3],DDetas[3][3]; 
	quotient3(num,Dnum,DDnum,u126,Du126,DDu126,etas,Detas,DDetas);
	interval eta,Deta[3],DDeta[3][3];
	Dsqrt3(etas,Detas,DDetas,eta,Deta,DDeta);
	interval g,Dg[2],DDg[2][2];
 
	// output from out.quoin.1255;
	// CHANGEABLE:
	/*
	g = interval("-1.0e-8","0.00219");
	Dg[0]= interval("-1.0e-6","0.00797");
	Dg[1]= interval("-0.0601","1.0e-8");
	DDg[0][0]= interval("-1.0e-4","0.0673");
	DDg[0][1]=DDg[1][0]= interval("-0.252","0.000288");
	DDg[1][1]= interval("-0.000291","1.1596");
	*/
 
 
	// output from out.quoin.sqrt;
	DDg[0][0]= interval("-1.078","7.024");
	DDg[0][1]=DDg[1][0]= interval("-9.048","0.317");
	DDg[1][1]= interval("-7.5e-06","12.73");
	Dg[0]= interval("-4.42e-10","0.12");
	Dg[1]= interval("-0.283","1.13e-10");
	g = interval("-1.996e-10","0.0231");
 
	int i,j;
	for (i=0;i<3;i++) for (j=i;j<3;j++)
		DDf[i][j]= DDg[0][0]*Da[i]*Da[j]+
			DDg[0][1]*Da[i]*Deta[j]+
			DDg[1][0]*Deta[i]*Da[j]+
			DDg[1][1]*Deta[i]*Deta[j]+
			Dg[0]*DDa[i][j]+Dg[1]*DDeta[i][j];
	for (i=0;i<3;i++) for (j=0;j<i;j++) DDf[i][j]=DDf[j][i];
	}

void test_quoin()
	{
	interval f,Df[3],DDf[3][3];
	interval x= "4.04", y = "4.17", z = "4.255";
	interval eps = interval(0.0, 0.02);
	interval trunc = "1.255";
	x = x+eps; y = y+eps; z = z+eps;
	compute_chainquoin(trunc,x,y,z,DDf);
	output3(f,Df,DDf);
	}

void test_eta2()
	{
	double x[3]={4.04,4.15,4.225};
	double z[3]={4.04,4.15,4.225};
	cout << eta2xx(x,z) << endl;
	cout << eta2xy(x,z) << endl;
	}



void run_quoin_2D()
	{
	cout << "*" << flush;
	interval a,b,c;
	interval q,Dq[2],f,fD[2],DD[2][2],fDD[2][2];
	double emin=10000.0, emax = -100000.0;
	interval biz = interval(emin,emax);
 
	//CHANGEABLE: CASE SQRT(2).
	int top=10000;  
	c=interval("1.41421356237309504880"); 
	// // for 2sq  let a go out to 1.255 (as y1 is short);
	interval wa = interval("0.255")/interval(double(top),double(top)); 
	char* tag = "QUASI-REG : sqrt(2)-trunc : rogers variables (2 dim'l)";
 
	// CHANGEABLE: CASE 1.255.
	//int top=2000;  
	//c=interval("1.255"); 
	// for 2.51 : rad(2(1.20846),2,2)>1.255);
	//interval wa = interval("0.20846")/interval(double(top),double(top)); 
	//char* tag = "QUASI-REG : 1.255-trunc : rogers variables (2 dim'l)";
	// END OF CHANGES:
 
	cout << tag << endl << flush;	
	int count = 0;
	double an,au,bn,bu,bmin,temp;
	int i,j,jmax,r,s; 
	for (i=0;i<3;i++) for (j=0;j<3;j++) fDD[i][j]=biz;
	for (i=0;i<3;i++) fD[i]=biz;
	f = biz;
 
	// MAIN LOOP
	for (i=0;i<top;i++) 
		{
		down();
		an = 1.0 + inf(wa)*double(i);
		up();
		au = 1.0 + sup(wa)*double(i+1);
		a = interval(an,au); 
		temp = sqrt(4.0 + an*(-an));
		down();
		bmin = 2.0/temp;     // radf(2a,2,2);
		up(); jmax = int(ceil((sup(c)+(-bmin))/inf(wa)));
		for (j=0;j<jmax;j++)
			{ 
			down();
			bn= bmin + inf(wa)*double(j); 
			if (bn>sup(c)) bn = sup(c);
			up();
			bu= bmin + sup(wa)*double(j+1);
			if (bu>sup(c)) bu = sup(c);
			b=interval(bn,bu);
			compute_quoin(a,b,c,q,Dq,DD);
			for (r=0;r<2;r++) for (s=0;s<2;s++)
				fDD[r][s]=combine(fDD[r][s],DD[r][s]);
			for (r=0;r<2;r++)
				fD[r]=combine(fD[r],Dq[r]);
			f=combine(f,q);
 
			if (0 == (count++ % 1000000))
				{
				cout << "*" << i << " " 
				  << band(fDD[0][0]) << endl << flush;
				}
			}
		}
	cout.precision(20);
	cout << tag << endl;
	for (i=0;i<2;i++) for (j=0;j<2;j++)
		cout << "DD[" << i << "," << j << "]: " 
			<< band(fDD[i][j]) << endl;
	for (i=0;i<2;i++)
		cout << "Df[" << i << "]: "
			<< band(fD[i]) << endl;
	cout << "f: " << band(f) << endl;
	cout << "*" << endl;
	cout << count << endl;
	}


void run_quoin()
	{
	cout << "*" << flush;
	interval DD[3][3],fDD[3][3];
	interval x1,x2,x3;
	double x1n,x1u,x2n,x2u,x3n,x3u;
	int i,j,k,r,s;
	double emin=10000.0, emax = -100000.0;
	interval biz = interval(emin,emax);
	int top=200;  // CHANGEABLE use 200 for 1.255 truncation;
	interval wa = interval("2.3001")/interval(double(top),double(top)); 
 
	// CHANGEABLE 1.255:
	//interval trunc = "1.255";
	//interval wb = wa;
	//char* tag = "QUASI-REG : 1.255-trunc : simplex variables (3 dim'l)";
	//int top2=top;
 
	// CHANGEABLE sqrt:
	interval trunc = "1.4142135623730950488";
	interval wb = four/interval(double(2*top),double(2*top)); 
	char* tag = "FLAT&QRTET : sq : simplex variables (3 dim'l)";
	int top2=2*top;
	//
 
	int count = 0;
	for (i=0;i<3;i++) for (j=0;j<3;j++) fDD[i][j]=biz;
 
	// MAIN LOOP
	for (i=0;i<top;i++) 
	for (j=0;j<top;j++)
	for (k=0;k<top2;k++) 
		{ 
		down();
		x1n = 4.0 + inf(wa)*double(i);
		x2n = 4.0 + inf(wa)*double(j);
		x3n = 4.0 + inf(wb)*double(k);
		up();
		x1u = 4.0 + sup(wa)*double(i+1);
		x2u = 4.0 + sup(wa)*double(j+1);
		x3u = 4.0 + sup(wb)*double(k+1);
		x1 = interval(x1n,x1u); 
		x2 = interval(x2n,x2u); 
		x3 = interval(x3n,x3u); 
		compute_chainquoin(trunc,x1,x2,x3,DD);
		for (r=0;r<3;r++) for (s=0;s<3;s++)
		fDD[r][s]=combine(fDD[r][s],DD[r][s]);
		if (0 == (count++ % 1000000))  // 10^6
			{
			cout << "*" << band(fDD[0][0]) << endl << flush;
			}
		}
	cout.precision(20);
	cout << tag << endl;
	for (i=0;i<3;i++) for (j=0;j<3;j++)
		cout << "DD[" << i << "," << j << "]: " 
			<< band(fDD[i][j]) << endl;
	cout << "*" << endl;
	cout << count << endl;
	}





main()
	{
	//test_quoin();
	run_quoin();
	//runeta2();
	}
