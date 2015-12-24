//  copyright (c) 1997, Thomas C. Hales, all rights reserved.


#include <iomanip.h>
#include <iostream.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "numerical.h"
#include "morefn.h"

/******************* DELAUNAY STUFF ***********************/

double mabs(double a) { return (a>0? a : -a); }
double max(double a,double b) { return (a>b? a : b); }
double min(double a,double b) { return (a<b? a : b); }

double myrand()
	{
	static int w =0;
	if (w==0) { srand(time(0)); w = 1; }
	double u = double(rand());
	double v =  double(/*stdlib.h*/RAND_MAX); 
	return u/v;
	}

double safesqrt(double t)
	{
	if (t<0.0) return 0.0; else return sqrt(t);
	}

double pos(double t)
	{
	return (t>0.0 ? t : 0.0);
	}

double safeacos(double t)
	{
	if (t<-1.0) return 3.14159265358979323846264338328;
	if (t>1.0) return 0.0;
	return acos(t);
	}

double delta(double x1,double x2,double x3,double x4,double x5,double x6)
        {
        return (-x1+x2+x3-x4+x5+x6)*(x1*x4) +
                (x1-x2+x3+x4-x5+x6)*(x2*x5)+
                (x1+x2-x3+x4+x5-x6)*(x3*x6) -
                (x2*x3*x4 +x1*x3*x5+x1*x2*x6+x4*x5*x6);
        }


double U(double a,double b,double c)
       {
        return -a*a-b*b-c*c+2.0*(a*b+b*c+c*a);
       }


double AA(double y1,double y2,double y3,double y4,double y5,double y6)
        {
        return y1*y2*y3 + (y1*(y2*y2+y3*y3-y4*y4)+
                          y2*(y1*y1+y3*y3-y5*y5)+
                          y3*(y1*y1+y2*y2-y6*y6))/2.0;
        }

double presolid(      double y1,double y2,double y3,
                        double y4,double y5,double y6)
        {
        double a = AA(y1,y2,y3,y4,y5,y6);
        return delta(y1*y1,y2*y2,y3*y3,y4*y4,y5*y5,y6*y6)/(a*a);
        }

double solid(double y1,double y2,double y3,double y4,double y5,double y6)
        {
        return 2.0*atan(safesqrt(presolid(y1,y2,y3,y4,y5,y6))/2.0);
        }

double solR(double y1,double y2,double y3)
	{
	if (y1>=y2) return 0.0;
	if (y2>=y3) return 0.0;
	return solid(y1,y2,y3,safesqrt(y3*y3-y2*y2),safesqrt(y3*y3-y1*y1),
		safesqrt(y2*y2-y1*y1));
	}

const double doct = 0.72090294951746509284124483502;
double gamma(double y1,double y2,double y3,double y4,
                        double y5, double y6)
        {
        double x1=y1*y1, x2 = y2*y2, x3=y3*y3, x4=y4*y4, x5=y5*y5, x6=y6*y6;
        double t = safesqrt(delta(x1,x2,x3,x4,x5,x6))/2.0;
        double a1,a2,a3,a4;
        a1 = AA(y1,y2,y3,y4,y5,y6);
        a2 = AA(y1,y5,y6,y4,y2,y3);
        a3 = AA(y4,y5,y3,y1,y2,y6);
        a4 = AA(y4,y2,y6,y1,y5,y3);
        return -doct*t/6.0 + (2.0/3.0)*(atan(t/a1)+atan(t/a2)+
                                          atan(t/a3)+atan(t/a4));
        };

double dihedral(double x1,double x2,double x3,double x4,
                double x5, double x6)
        {
        double d4 = -(x2*x3) - x1*x4 + x2*x5 + x3*x6 - x5*x6 +
                        x1*(-x1 + x2 + x3 - x4 + x5 + x6);
        return safeacos(d4/safesqrt(U(x1,x2,x6)*U(x1,x3,x5)));
        }

double dihedraly(double y1,double y2,double y3,double y4,double y5,
        double y6)
        {
        return dihedral(y1*y1,y2*y2,y3*y3,y4*y4,y5*y5,y6*y6);
        }

double dihR(double y1,double y2,double y3)
	{
	static const double pi2 = 1.57079632679489661923132169164;
	if (y1>=y2) return pi2;
	if (y2>=y3) return 0.0;
	return dihedral(y1*y1,y2*y2,y3*y3,y3*y3-y2*y2,y3*y3-y1*y1,
			y2*y2-y1*y1);
	}

double eta2(double x1,double x2,double x3)
        {
        return x1*x2*x3/U(x1,x2,x3);
        }

double radf(double y1,double y2,double y3)
        {
        return safesqrt(eta2(y1*y1,y2*y2,y3*y3));
        }

double chi(   double x1,double x2,double x3,
                double x4,double x5,double x6)
        {
        return x3*x6*(x4+x5-x6) + x2*x5*(x4-x5+x6)+x1*x4*(-x4+x5+x6)
               -2.0*x4*x5*x6;
        }

double rad(double x1,double x2,double x3,double x4,
                double x5,double x6)
        {
        double c = chi(x1,x2,x3,x4,x5,x6);
        return safesqrt(eta2(x4,x5,x6)+
                        c*c/(4.0*U(x4,x5,x6)*delta(x1,x2,x3,x4,x5,x6)));
        }

double rady(double y1,double y2,double y3,double y4,double y5,double y6)
	{
	return rad(y1*y1,y2*y2,y3*y3,y4*y4,y5*y5,y6*y6);
	}

double arc(double y1,double y2,double y6)
    {
    return acos( (y1*y1+y2*y2-y6*y6)/(2.*y1*y2));
    }

double circum2(double x1,double x2,double x3,double x4,
                double x5,double x6)
        {
        double c = chi(x1,x2,x3,x4,x5,x6);
        return (eta2(x4,x5,x6)+
                        c*c/(4.0*U(x4,x5,x6)*delta(x1,x2,x3,x4,x5,x6)));
        };

double vor_analytic(double y1,double y2,double y3,double y4,
	double y5,double y6)
	{
	double x1=y1*y1, x2=y2*y2, x3=y3*y3, 
		x4=y4*y4, x5=y5*y5, x6=y6*y6;
	double u126 = U(x1,x2,x6), u135 = U(x1,x3,x5);
	double u234 = U(x2,x3,x4);
	double vol = (x1*(x2 + x6 - x1) + x2*(x1 + x6 - x2))*
              chi(x4, x5, x3, x1, x2, x6)/u126 +   
           (x2*(x3 + x4 - x2) + x3*(-x3 + x4 + x2))*
              chi(x6, x5, x1, x3, x2, x4)/u234 + 
           (x1*(-x1 + x3 + x5) + x3*(x1 - x3 + x5))*
              chi(x4, x6, x2, x1, x3, x5)/u135;
	vol = vol/(48.0*safesqrt(delta(x1,x2,x3,x4,x5,x6)));
	return -4.0*doct*vol + 4.0*solid(y1,y2,y3,y4,y5,y6)/3.0;
	}

double octavor(double y1,double y2,double y3,double y4,
        double y5,double y6)
        {
        return (vor_analytic(y1,y2,y3,y4,y5,y6) +
                vor_analytic(y1,y5,y6,y4,y2,y3))/2.0;
        }


double octavorVc(double y1,double y2,double y3,double y4,
        double y5,double y6)
        {
        return (vorVc(y1,y2,y3,y4,y5,y6) +
                vorVc(y1,y5,y6,y4,y2,y3))/2.0;
        }

double vorR(double y1,double y2,double y3)
	{
	double sol = solR(y1,y2,y3);
	double vol = y1*safesqrt(y2*y2-y1*y1)*safesqrt(y3*y3-y2*y2)/6.0;
	return -4.0*doct*vol + 4.0*sol/3.0;
	}

double pretilde(double h,double t)
        {
        return -1.232888761906277 + 0.4806019663449*t*(t+h)*h;
        }

double tilde(double h)
        {
        return pretilde(h,1.255);
	}

double crown(double h,double t)
        {
        return 2.0*3.141592653589793*(1.0-h/t)*(pretilde(h,t)-tilde(1.255));
        }

double kappa(double y1,double y2,double y3,double y4,double y5,
            double y6)
    {
	double pi = 3.141592653589793;
    return (crownV(y1/2.0)*dihedraly(y1,y2,y3,y4,y5,y6)/(2.0*pi))
        + vorAnchor(y1,y2,y6)+ vorAnchor(y1,y3,y5);
    }


double crown(double h)
        {
        return crown(h,radf(2.0*h,2.0,2.51));
        }

double phi(double h,double t)
    {
    static const double doc = 0.72090294951746509284124;
    return 2.0*(2.0 - doc*h*t*(t+h))/3.0;
    }

double quoin(double a,double b,double c)
    {
    if ((a>b)||(b>c)) return 0;
    double u = sqrt((c*c - b*b)/(b*b - a*a));
    return -(a*a + a*c - 2*c*c)*(c - a)*atan(u)/6.0 +
        (a*(b*b - a*a)*u)/6.0 - (2.0/3.0)*c*c*c*atan((u*(b - a))/(b + c));
    }


// function A from SPIII:
double A_of_vorVc(double h)
    {
    static const double t0=1.255;
    return (1.0-h/t0)*(phi(h,t0)-phi(t0,t0));
    }

// function B from SPIV: 
double B_of_vorVc(double y)
    {
    static const double t0=1.255;
    return phi(t0,t0) + A_of_vorVc(y/2.0);
    }


double sigtilde(double h,double t)
	{
	if (h>t) h=t;
	return 4.0/3.0 - 2.0*doct*h*t*(h+t)/3.0;
	}

double skelV(double y1,double y2,double y3,double y4,double y5,double y6)
	{
	double sq =  1.41421356237309504880;
	double eta126 = radf(y1,y2,y6);
	double eta135 = radf(y1,y3,y5);
	return vorR(y1/2.0,eta126,sq)+ vorR(y2/2.0,eta126,sq)+
		vorR(y1/2.0,eta135,sq)+ vorR(y3/2.0,eta135,sq)+
		(dihedraly(y1,y2,y3,y4,y5,y6)-dihR(y1/2.0,eta126,sq)-
			dihR(y1/2.0,eta135,sq))*(1-y1/(2.0*sq))*
				sigtilde(y1/2.0,sq);
	}

	
double sigma_star(double h1,double h2,double h3,
		double a1,double a2,double a3,double trunc)
	{
	static double pi = 3.14159265358979323;
	double sigmamax = sigtilde(trunc,trunc);
	h1 = min(h1,trunc); h2 = min(h2,trunc); h3 = min(h3,trunc);
	return a1*pos(1.0-h1/trunc)*(sigtilde(h1,trunc)- sigmamax) +
	       a2*pos(1.0-h2/trunc)*(sigtilde(h2,trunc)- sigmamax) +
	       a3*pos(1.0-h3/trunc)*(sigtilde(h3,trunc)- sigmamax) +
		(a1+a2+a3-pi)*sigmamax;
	}

double one_sigma_star(double h1,double a1,double trunc)
	// one term from sigma_star
	{
	static double pi = 3.14159265358979323;
	double sigmamax = sigtilde(trunc,trunc);
	h1 = min(h1,trunc); 
	return a1*pos(1.0-h1/trunc)*(sigtilde(h1,trunc)- sigmamax);
	}

double dihVc(double x1,double x2,double x3,double x4,double x5,double x6)
	{
	double y1 = safesqrt(x1);
	double h1 = y1/2.0; 
	double dih1 = dihedral(x1,x2,x3,x4,x5,x6);
	return one_sigma_star(h1,dih1,1.255);
	}

double dihSqc(double x1,double x2,double x3,double x4,double x5,double x6)
	{
	const static double sq = 1.414213562373095048801;
	double y1 = safesqrt(x1);
	double h1 = y1/2.0; 
	double dih1 = dihedral(x1,x2,x3,x4,x5,x6);
	return one_sigma_star(h1,dih1,sq);
	}

double sfacelift(double y1,double y2,double y6,double t)
	{
	double eta = radf(y1,y2,y6);
	if (eta>=t) return 0.0;
	double a = y1/2.0; 
	double sigmamax = sigtilde(t,t);
	double d = dihR(a,eta,t);
	return vorR(a,eta,t) 
		+ sigmamax*(d*(1.0-a/t)- solR(a,eta,t))
		- sigtilde(a,t)*(1.0-a/t)*d;
	}

double sigmafacelift(double y1,double y2,double y6,double t)
	{
	return sfacelift(y1,y2,y6,t) + sfacelift(y2,y1,y6,t);
	}

double vorVc(double y1,double y2,double y3,double y4,double y5,
	double y6,double trunc) // y1 is the long edge.
	{
	double h1 = y1/2.0, h2 = y2/2.0, h3 = y3/2.0;
	double dih1 = dihedraly(y1,y2,y3,y4,y5,y6);
	double dih2 = dihedraly(y2,y3,y1,y5,y6,y4);
	double dih3 = dihedraly(y3,y1,y2,y6,y4,y5);
	return sigma_star(h1,h2,h3,dih1,dih2,dih3,trunc) +
		sigmafacelift(y1,y2,y6,trunc) +
		sigmafacelift(y2,y3,y4,trunc) +
		sigmafacelift(y1,y3,y5,trunc);
	}
double vorVc(double y1,double y2,double y3,double y4,double y5,double y6) 
	{ return vorVc(y1,y2,y3,y4,y5,y6,1.255); }
double vorSqc(double y1,double y2,double y3,double y4,double y5,double y6)
	{ return vorVc(y1,y2,y3,y4,y5,y6,1.414213562373095); }


double tauVc(double y1,double y2,double y3,double y4,double y5,double y6,
	double t)
	{
	return solid(y1,y2,y3,y4,y5,y6)*0.1004445714270561 
		- vorVc(y1,y2,y3,y4,y5,y6,t);
	}
double tauVc(double y1,double y2,double y3,double y4,double y5,double y6)
	{ return tauVc(y1,y2,y3,y4,y5,y6,1.255); }
double tauSqc(double y1,double y2,double y3,double y4,double y5,double y6)
	{ return tauVc(y1,y2,y3,y4,y5,y6,1.414213562373095); }

double lng = 1.6;
double sol_cor(double y1,double y2,double y6)
	// the negative area of the snippet cut by perp to F1 from a knob.
	{
	double trunc=1.255;
	double cospsi = (y1*y1+trunc*trunc-lng*lng)/(2.0*y1*trunc);
	if (cospsi>1.0) return 0.0;
	double eta = radf(y1,y2,y6);
	double costheta = y1/(2.0*eta);
	if (costheta>1.0) costheta=1.0;
	double a = trunc*cospsi;
	double b = a/costheta;
	if (b>=trunc) return 0.0;
	return -dihR(a,b,trunc)*(1.0-cospsi) + solR(a,b,trunc);
	}

double corsigma(double y1,double y2,double y3,double y4,double y5,double y6)
	{
	double trunc = 1.255;
	double sigmamax = sigtilde(trunc,trunc);
	double h = y1/2.0;
	double cospsi = (y1*y1+trunc*trunc-lng*lng)/(2.0*y1*trunc);
	if (cospsi>1.0) cospsi=1.0;
	double costhet = h/trunc;
	if (costhet>1.0) costhet=1.0;
	return dihedraly(y1,y2,y3,y4,y5,y6)*((1.0-cospsi)*sigmamax
		+(1.0-costhet)*(sigtilde(h,trunc)-sigmamax)) +
		sigmamax*(sol_cor(y1,y2,y6)+sol_cor(y1,y3,y5))+
		sfacelift(y1,y2,y6,trunc)+sfacelift(y1,y3,y5,trunc);
	}

double corsolid(double y1,double y2,double y3,double y4,double y5,double y6)
	{
	double trunc = 1.255;
	double cospsi = (y1*y1+trunc*trunc-lng*lng)/(2.0*y1*trunc);
	if (cospsi>1.0) return 0.0;
	return dihedraly(y1,y2,y3,y4,y5,y6)*(1.0-cospsi) +
		sol_cor(y1,y2,y6)+ sol_cor(y1,y3,y5);
	}

double cortau(double y1,double y2,double y3,double y4,double y5,double y6)
	{
	return corsolid(y1,y2,y3,y4,y5,y6)*0.1004445714270561 -
		corsigma(y1,y2,y3,y4,y5,y6);
	}

double tau(double y1,double y2,double y3,double y4,double y5,double y6)
	{
	return solid(y1,y2,y3,y4,y5,y6)*0.1004445714270561 
		- gamma(y1,y2,y3,y4,y5,y6);
	}

double tauR(double y1,double y2,double y3)
	{
	return solR(y1,y2,y3)*0.1004445714270561 
		- vorR(y1,y2,y3);
	}

double tau_analytic(double y1,double y2,double y3,double y4,double y5,double y6)
	{
	return solid(y1,y2,y3,y4,y5,y6)*0.1004445714270561 
		- vor_analytic(y1,y2,y3,y4,y5,y6);
	}

// Various scoring on upright quarters.

double vorNu(double x1,double x2,double x3,double x4,double x5,
        double x6)
    {
    return octavor(x1,x2,x3,x4,x5,x6) + 0.5*(
        vorVc(x1,x2,x3,x4,x5,x6)-vorVc(x1,x6,x5,x4,x3,x2));
    }

double gammaNu(double x1,double x2,double x3,double x4,double x5,
        double x6)
    {
    return gamma(x1,x2,x3,x4,x5,x6) + 0.5*(
        vorVc(x1,x2,x3,x4,x5,x6)-vorVc(x1,x6,x5,x4,x3,x2));
    }

double tauVnu(double x1,double x2,double x3,double x4,double x5,
        double x6)
    {
    return octatau(x1,x2,x3,x4,x5,x6) - 0.5*(
        vorVc(x1,x2,x3,x4,x5,x6)-vorVc(x1,x6,x5,x4,x3,x2));
    }

double tauGnu(double x1,double x2,double x3,double x4,double x5,
        double x6)
    {
    return tau(x1,x2,x3,x4,x5,x6) - 0.5*(
        vorVc(x1,x2,x3,x4,x5,x6)-vorVc(x1,x6,x5,x4,x3,x2));
    }

//---- stuff for quad clusters.
double vorVc9(double y[9])
	{
	return vorVc(y[0],y[1],y[2],y[3],y[4],y[5])+
		vorVc(y[6],y[1],y[2],y[3],y[7],y[8]);
	}

double solid9(double y[9])
	{
	return solid(y[0],y[1],y[2],y[3],y[4],y[5])+
		solid(y[6],y[1],y[2],y[3],y[7],y[8]);
	}

double dih9_0(double y[9])
	{
	return dihedraly(y[0],y[1],y[2],y[3],y[4],y[5]);
	}

double dih9_2(double y[9])
	{
	return dihedraly(y[2],y[0],y[1],y[5],y[3],y[4])+
		dihedraly(y[2],y[6],y[1],y[8],y[3],y[7]);
	}

double dih9_1(double y[9])
	{
	return dihedraly(y[1],y[0],y[2],y[4],y[3],y[5])+
		dihedraly(y[1],y[6],y[2],y[7],y[3],y[8]);
	}

double dih9_6(double y[9])
	{
	return dihedraly(y[6],y[1],y[2],y[3],y[7],y[8]);
	}

double crossdiag(double y[9])
	{
	double    x3=y[3]*y[3], 
		x4=y[4]*y[4], x5=y[5]*y[5], 
		x7=y[7]*y[7], x8=y[8]*y[8];
	double dd = dihedraly(y[3],y[1],y[5],y[0],y[4],y[2])+
			dihedraly(y[3],y[1],y[8],y[6],y[7],y[2]);
	static double pi = 3.14159265358979323;
	double c0;
	double c1 = -x3*x3-x5*x8-x4*x7+x3*(x5+x8+x4+x7)+x5*x7+x4*x8;
	if (dd>pi)
		{
		c0 = cos(2.0*pi-dd)*sqrt(U(x3,x5,x4)*(U(x3,x8,x7)));
		return sqrt((c0-c1)/(-2.0*x3));
		}
	else
		{
		c0 = cos(dd)*sqrt(U(x3,x5,x4)*(U(x3,x8,x7)));
		return  sqrt((c0-c1)/(-2.0*x3));
		}
	}
