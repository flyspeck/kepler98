#include <iomanip.h>
#include <iostream.h>
#include <fstream.h>
#include <math.h>
#include <stdlib.h>
#include "numerical.h"
#include "gradient.h"


const static double pi = 3.141592653589793;
const static double zeta = 1.813941816806566;
const static double pt = 0.0553736456684637;
const static double tildemax = 0.6670811193203888;

double vorAV(double y[6])
	{
	return 0.5*(vor_analytic(y[0],y[1],y[2],y[3],y[4],y[5])+
			vor_analytic(y[0],y[5],y[4],y[3],y[2],y[1]));
	}

double vorVcAV(double y[6])
	{
	return 0.5*(vorVc(y[0],y[1],y[2],y[3],y[4],y[5])+
			vorVc(y[0],y[5],y[4],y[3],y[2],y[1]));
	}

double crownV(double h) { return -crown(h); }

double tauAnchor1(double y1,double y2,double y6,double eta1)
	{
	double eta = radf(y1,y2,y6);
	if (eta>eta1) return 0.0;
	double h=y1/2.0;
	return -dihR(h,eta,eta1)*crown(h)/(2*pi) +
		tauR(h,eta,eta1)-solR(h,eta,eta1)*tildemax;
	}

// CHAFF: (used in 1.18 truncation argument in Part IV, April 97)
double tauK(double y1,double y2,double y3,double y4,double y5,double y6)
	{
	double h = (y1>2.36 ? 1.18 : y1/2.0);
	double pr = pretilde(h,1.18);
	double dd = dihedraly(y1,y2,y3,y4,y5,y6);
	double r1 = radf(y1,y2,y6);
	double r2 = radf(y1,y3,y5);
	double d1 = (r1>1.18 ? 0.0 : dihR(h,r1,1.18));
	double d2 = (r2>1.18 ? 0.0 : dihR(h,r2,1.18));
	double t1 = (r1>1.18 ? 0.0 : tauR(h,r1,1.18)+tauR(y2/2.0,r1,1.18));
	double t2 = (r2>1.18 ? 0.0 : tauR(h,r2,1.18)+tauR(y3/2.0,r2,1.18));
	return pr*(dd-d1-d2) + t1+t2;
	}




double tauAnchor1(double y1,double y2,double y6)
	{
	return tauAnchor1(y1,y2,y6,radf(y1,2.0,2.51));
	}

double tauAnchor2(double y1,double y2,double y6,double eta1)
	{
	double eta0 = radf(y1,y2,y6);
	if (eta0>eta1) return 0.0;
	double neuf = tauR(y2/2.0,eta0,eta1);
	double old = 
		dihR(y2/2.0,eta0,eta1)*(1-y2/2.51)*(tilde(y2/2.0)-tildemax) +
		solR(y2/2.0,eta0,eta1)*tildemax;
	return neuf-old;
	}

double tauAnchor2(double y1,double y2,double y6)
	{
	return tauAnchor2(y1,y2,y6,radf(y1,2.0,2.51));
	}

double tauAnchor(double y1,double y2,double y6)
	{
	return tauAnchor1(y1,y2,y6)+tauAnchor2(y1,y2,y6);
	}

double vorAnchor(double y1,double y2,double y6)
	{
	return -tauAnchor(y1,y2,y6);
	}

double tauU(double y1,double y2,double y3,double y4,double y5,double y6)
	{
	return tau(y1,y2,y3,y4,y5,y6)+ // should be just tau
		(2.0*pi-dihedraly(y1,y2,y3,y4,y5,y6))*crown(y1/2.0)/(2.0*pi)
		+tauAnchor(y1,y2,y6)+tauAnchor(y1,y3,y5);
	}

double gammaU(double y1,double y2,double y3,double y4,double y5,double y6)
        {
        return solid(y1,y2,y3,y4,y5,y6)* 0.1004445714270561
                        - tauU(y1,y2,y3,y4,y5,y6);
        }



double octatau(double y1,double y2,double y3,double y4,double y5,
	double y6)
	{
	return solid(y1,y2,y3,y4,y5,y6)* 0.1004445714270561
                        - octavor(y1,y2,y3,y4,y5,y6);
	}

double tauU2(double y1,double y2,double y3,double y4,double y5,double y6)
	{
	return octatau(y1,y2,y3,y4,y5,y6)+ // should be tau -BUGBUG
		(pi-dihedraly(y1,y2,y3,y4,y5,y6))*crown(y1/2.0)/(2.0*pi)
		+tauAnchor(y1,2.0,2.0); // BUGBUG should be (y1,y2,y6);
	}


double gammaU2(double y1,double y2,double y3,double y4,double y5,double y6)
	{
	return solid(y1,y2,y3,y4,y5,y6)* 0.1004445714270561
			- tauU2(y1,y2,y3,y4,y5,y6);
	}


double vorU2(double y1,double y2,double y3,double y4,double y5,double y6)
	{
	return vor_analytic(y1,y2,y3,y4,y5,y6)-gamma(y1,y2,y3,y4,y5,y6)
		+solid(y1,y2,y3,y4,y5,y6)* 0.1004445714270561
			- tauU2(y1,y2,y3,y4,y5,y6);
	}


double eta2side(double z[6])
	{
	return eta2(z[1]*z[1],z[2]*z[2],z[3]*z[3]);
	}

double triangle(double y1,double y2,double y6)
    {
    return acos( (y1*y1+y2*y2-y6*y6)/(2*y1*y2) );
    }
 
double psi (double y) { return triangle(y,1.255,1.6); }
 
double beta (double psi0, double y0,double y1,double y6)
    {
    double p = cos(psi0);
    double t = cos(triangle(y0,y1,y6));
    return acos( sqrt((p*p-t*t)/(1.0-t*t)) );
    }
 

