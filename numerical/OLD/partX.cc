#include <iomanip.h>
#include <stdlib.h>
#include "numerical.h"
#include "constants.h"
#include "morefn.h"
#include <math.h>
#include "quoinfn.h"

class iter {
    double* xmin;
	double* xmax;
	double* x;
    int numiter,numargs,nconstr;
    void (*func)(int numargs,int whichFn,double* x,double* ret,void*);
    void (*constraintfunc)
        (int numargs,int which,double* x,double* ret,void*);
    iter(int);
    ~iter();
    };
static double kappa(double y1,double y2,double y3,double y4,double y5,
            double y6)
    {
    return (crownV(y1/2.0)*dihedraly(y1,y2,y3,y4,y5,y6)/(2.0*global::pi))
        + vorAnchor(y1,y2,y6)+ vorAnchor(y1,y3,y5);
    }
double U(double a,double b,double c); // in numerical.cc

double SeanVol(double y1,double y2,double y3,double y4,
	double y5,double y6)
	// dodecrad truncated volume of VOronoi
	{
	double dr = global::dodecrad/2.;
	double dt = 0.7209029495174648;
	return
		(-1./dt)*(vorVc(y1,y2,y3,y4,y5,y6,dr)/4.
			-solid(y1,y2,y3,y4,y5,y6)/3.);
	}

// cross diag in terms of edges . Equivalent to version in numerical.cc

double crossdiagNew(double y1,double y2,double y3,double y4,double
		y5,double y6,double y7,double y8,double y9)
	{
	double x1=y1*y1, x2=y2*y2, x3=y3*y3, x4=y4*y4, x5=y5*y5, x6=y6*y6,
			x7=y7*y7, x8=y8*y8, x9=y9*y9;
	double s = dihedraly(y4,y2,y6,y1,y5,y3)+dihedraly(y4,y2,y9,y7,y8,y3);
	if (s>global::pi) { s = 2*global::pi - s; }
	double uc = cos(s)*safesqrt(U(x4,x8,x9)*U(x4,x5,x6));
	double c= -x4*x4+x4*(x6+x9+x8+x5) - x6*x9-x5*x8+x5*x9+x6*x8;
	double a= -2.0*x4;
	return safesqrt((uc-c)/a);
	}

static double quoin(double a,double b,double c)
  	{
	double u = sqrt((c*c - b*b)/(b*b - a*a));
    if ((a>b)||(b>c)) return 0;
    return -(a*a + a*c - 2*c*c)*(c - a)*atan(u)/6.0 + 
		(a*(b*b - a*a)*u)/6.0 - (2.0/3.0)*c*c*c*atan((u*(b - a))/(b + c));
	}
static double phi(double h,double t)
    {
    static const double doc = 0.72090294951746509284124;
    return 2.0*(2.0 - doc*h*t*(t+h))/3.0;
    }

static double quoinH(double y1,double y2,double y6)
	{
	return quoin(y1/2.0,radf(y1,y2,y6),1.255);
	}

static double u135M(double y1,double y3,double y5)
	{
	static const double doc = 0.72090294951746509284124;
	double x1 = y1*y1, x3 = y3*y3, x5=y5*y5;
	double a = y1/2.0, ap= y3/2.0, b = radf(y1,y3,y5), c = 1.255;
	double c2b2 = c*c-b*b;
	if (c2b2<0) return 0;
	double root = -c2b2*sqrt(c2b2);
	double Dquo135 = a*root/(3.0*b*sqrt(b*b-a*a));
	double Dquo315 = ap*root/(3.0*b*sqrt(b*b-ap*ap));
	double u135 = U(x1,x3,x5);
	return (Dquo135+Dquo315)*x1*x3*(x5+x3-x1)*(x5+x1-x3)*(-4.0*doc)/
				(2.0*b*u135);
	}

static void nofunc(int numargs,int whichFn,double* x,double* ret,void*)
    {
    cout << "nofunc should not be called" << endl << flush;
    }

static double xix = 0.2066818092722632;
static double taudih(double x1,double x2,double x3,double x4,double x5,
		double x6)
	{
	return tau(x1,x2,x3,x4,x5,x6)
			+xix*dihedraly(x1,x2,x3,x4,x5,x6);
	}

static double tauXg(double x1,double x2,double x3,double x4,double x5,
		double x6)
	{
	return tau(x1,x2,x3,x4,x5,x6) - 0.5*(
		vorVc(x1,x2,x3,x4,x5,x6)-vorVc(x1,x6,x5,x4,x3,x2)); 
	}

static double sigXg(double x1,double x2,double x3,double x4,double x5,
		double x6)
	{
	return gamma(x1,x2,x3,x4,x5,x6) + 0.5*(
		vorVc(x1,x2,x3,x4,x5,x6)-vorVc(x1,x6,x5,x4,x3,x2)); 
	}

static double tauXo(double x1,double x2,double x3,double x4,double x5,
		double x6)
	{
	return octatau(x1,x2,x3,x4,x5,x6) - 0.5*(
		vorVc(x1,x2,x3,x4,x5,x6)-vorVc(x1,x6,x5,x4,x3,x2)); 
	}

static double sigXo(double x1,double x2,double x3,double x4,double x5,
		double x6)
	{
	return octavor(x1,x2,x3,x4,x5,x6) + 0.5*(
		vorVc(x1,x2,x3,x4,x5,x6)-vorVc(x1,x6,x5,x4,x3,x2)); 
	}



double fa(double y1,double y2,double y3,double y4,double y5,double y6,
    double fA[2],int type[4])
    {
    double g = (type[0]||type[1] ? vor_analytic(y1,y2,y3,y4,y5,y6) :
        gamma(y1,y2,y3,y4,y5,y6) );
    return g
    + fA[0]*(dihedraly(y2,y1,y3,y5,y4,y6)+ dihedraly(y3,y1,y2,y6,y4,y5))
    + fA[1]/4.0;
    }
 
double fb(double y1,double y2,double y3,double y4,double y5,double y6,
    double fA[2],int type[4])
    {
    double g = (type[1]||type[2] ? vor_analytic(y1,y2,y3,y4,y5,y6) :
        gamma(y1,y2,y3,y4,y5,y6) );
    return g
    + fA[0]*(dihedraly(y2,y1,y3,y5,y4,y6))+
    + fA[1]/4.0;
    }
 
double fc(double y1,double y2,double y3,double y4,double y5,double y6,
    double fA[2],int type[4])
    {
    double g = (type[2]||type[3] ? vor_analytic(y1,y2,y3,y4,y5,y6) :
        gamma(y1,y2,y3,y4,y5,y6) );
    return g
    +fA[1]/4.0;
    }
 
double fd(double y1,double y2,double y3,double y4,double y5,double y6,
    double fA[2],int type[4])
    {
    double g = (type[3]||type[0] ? vor_analytic(y1,y2,y3,y4,y5,y6) :
        gamma(y1,y2,y3,y4,y5,y6) );
    return g
    + fA[0]*(dihedraly(y3,y1,y2,y6,y4,y5))
    + fA[1]/4.0;
    }

double dpi(double y1,double y2,double y3,double y4,double y5,double y6)
    {
    return dihedraly(y1,y2,y3,y4,y5,y6)-global::pi2;
    }


static double fx(double y1,double y2,double y3,double y4,double y5,double y6)
	{
    return gamma(y1,y2,y3,y4,y5,y6)
        - 3.0508*(dihedraly(y2,y1,y3,y5,y4,y6)+ dihedraly(y3,y1,y2,y6,y4,y5))
        + 9.494/4.0;
    }

static double fA(double y[6])
	{
	static const double X[7]=
	{-1.466301, -0.759678,  0.074232,  0.578203,  1.768768,  0.910089, -1.954001};
	static const double cA = X[0];
	static const double cy[6] = {X[1],X[2],X[3],0,X[4],X[5]};
	static const double d = X[6];
	double total  =cA;
	int i;
	for (i=0;i<6;i++) total += cy[i]*y[i];
	total += fx(y[0],y[1],y[2],y[3],y[4],y[5]);
	total += d*(dihedraly(y[0],y[1],y[2],y[3],y[4],y[5])- global::pi2);
	return total;
	}

static double V(double y1,double y2,double y3,double y4,double y5,double y6)
	{
	return octavorVc(y1,y2,y3,y4,y5,y6)-gamma(y1,y2,y3,y4,y5,y6);
	}

static double VV(double y1,double y2,double y3,double y4,double y5,double y6)
	{
	return octavorVc(y1,y2,y3,y4,y5,y6)-octavor(y1,y2,y3,y4,y5,y6);
	}

static double overlap(double y1,double y2,double y5)
	{
	return delta(y1*y1,y2*y2,1.255*1.255,1.6*1.6,1.6*1.6,y5*y5);
	}

static double areadiff(double y1,double y2,double y5)
	{
	if (overlap(y1,y2,y5)<=0.0) return 0;
	double a1 = dihedraly(y1,y2,1.255,1.6,1.6,y5);
	double a2 = dihedraly(y2,y1,1.255,1.6,1.6,y5);
	double cospsi1 = (y1*y1+1.255*1.255-1.6*1.6)/(2.51*y1);
	double cospsi2 = (y2*y2+1.255*1.255-1.6*1.6)/(2.51*y2);
	double s =  solid(y1,y2,1.255,1.6,1.6,y5);
	return a1*(1.0-cospsi1)+a2*(1.0-cospsi2) - s;
	}

static double A(double h)
    {
    static const double t0=1.255;
    return (1.0-h/t0)*(phi(h,t0)-phi(t0,t0));
    }
 
static double B(double y)
    {
    static const double t0=1.255;
    return phi(t0,t0) + A(y/2.0);
    }
 
static double Vee0(double y1,double y2,double y3,double y4,double y5,double y6)
    {
    double x1=y1*y1, x2=y2*y2,x3=y3*y3, x4=y4*y4, x5=y5*y5, x6=y6*y6;
    double delta4 = -(x2*x3) - x1*x4 + x2*x5 + x3*x6 - x5*x6 + x1*(-x1 + x2 + x3 - x4 + x5 + x6);
    double delta6 = -(x1*x2) + x1*x4 + x2*x5 - x4*x5 + x3*(x1 + x2 - x3 + x4 + x5 - x6) - x3*x6;
    double u135 = U(x1,x3,x5);
    return -B(y1)*y1*delta6 + B(y2)*y2*u135 - B(y3)*y3*delta4;
    }


static double Vee1(double y1,double y2,double y3,double y4,double y5,double y6)
    {
    static const double zp = global::zetapt;
    double x1=y1*y1, x2=y2*y2,x3=y3*y3, x4=y4*y4, x5=y5*y5, x6=y6*y6;
    double delta4 = -(x2*x3) - x1*x4 + x2*x5 + x3*x6 - x5*x6 + x1*(-x1 + x2 + x3
 - x4 + x5 + x6);
    double delta6 = -(x1*x2) + x1*x4 + x2*x5 - x4*x5 + x3*(x1 + x2 - x3 + x4 + x5 - x6) - x3*x6;
    double u135 = U(x1,x3,x5);
    return -(B(y1)-zp)*y1*delta6 + (B(y2)-zp)*y2*u135 - (B(y3)-zp)*y3*delta4;
    }

static double arc(double y1,double y2,double y6)
	{
	return acos( (y1*y1+y2*y2-y6*y6)/(2.*y1*y2));
	}
 
static double tempAngle(double y1,double y2,double y3)
	{
	double h = y2/2.;
	double ar = arc(y1,1.255,1.6);
	double eta = radf(y1,y2,y3);
	double c = h/cos(ar);
	return dihR(h,eta,c);
	}

int INEQ_NUMBER=0;
static void generic(int numargs,int whichFn,double* x, double* ret,void*)
	{
	switch (INEQ_NUMBER) {
		// start of first page of inequalities for Section 2, SPIV.
		case 248:
			*ret= (SeanVol(x[0],x[1],x[2],x[3],x[4],x[5])+
				SeanVol(x[6],x[1],x[2],x[3],x[7],x[8]))
			-0.68*
			(solid(x[0],x[1],x[2],x[3],x[4],x[5])+
				solid(x[6],x[1],x[2],x[3],x[7],x[8]))
					//-0.10957*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]);
			-0.10957*
			(dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])+
				dihedraly(x[6],x[1],x[2],x[3],x[7],x[8]));
			break;
		case 247:
			*ret= -gamma(x[0],x[1],x[2],x[3],x[4],x[5])
				+0.252*dihedraly(x[1],x[0],x[2],x[4],x[3],x[5]); break;
		case 246:
			*ret = -vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 245:
			*ret = -dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 244:
			*ret = dihedraly(x[1],x[0],x[2],x[4],x[3],x[5])-
				beta(acos(x[0]/2.51),x[1],x[0],x[5]); break;
		case 2442:
			*ret = dihedraly(x[1],x[0],x[2],x[4],x[3],x[5])-
				tempAngle(x[0],x[1],x[5]); break;
		case 2441:
			*ret = dihedraly(x[1],x[0],x[2],x[4],x[3],x[5])-
				beta(arc(x[0],1.255,1.6),x[1],x[0],x[5]); break;
		case 243:
			*ret= dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 242:
			*ret=-vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
				-vorVc(x[6],x[1],x[2],x[3],x[7],x[8]);
			break;
		case 240:
			*ret=tauVc(x[0],x[1],x[2],x[3],x[4],x[5])+
				tauVc(x[6],x[1],x[2],x[3],x[7],x[8]);
			break;
		case 241:
			*ret=-vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
				-vorVc(x[6],x[1],x[2],x[3],x[7],x[8]);
			break;
		case 238:
			*ret=D1Vee(x[0],x[1],x[2],x[3],x[4],x[5]);
			break;
		case 239:
			*ret=D1Vee1(x[0],x[1],x[2],x[3],x[4],x[5]);
			break;
		case 91:
			*ret=-gamma(x[0],x[1],x[2],x[3],x[4],x[5])
				-0.145*x[0]-0.081*x[1]-0.081*x[2]
				-0.133*x[4]-0.133*x[5]+1.17401; break;
		case 92:
			*ret=-vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
				-0.145*x[0]-0.081*x[1]-0.081*x[2]
				-0.133*x[4]-0.133*x[5]+1.17401; break;
		case 93:
			*ret=-vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
				-0.145*x[0]-0.081*x[1]-0.081*x[2]
				-0.133*x[4]-0.133*x[5]+1.17401; break;
		case 94: case 95:
			*ret=-vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
				-0.145*x[0]-0.081*x[1]-0.081*x[2]
				-0.133*x[4]-0.133*x[5]+1.17401; break;
		case 96:
			*ret=-gamma(x[0],x[1],x[2],x[3],x[4],x[5])
				+0.00897
				+1.0e-4 -0.153*(x[3]+x[4]+x[5]-4.0-global::sqrt8);
			break;
		case 97:
			*ret=-vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
				+0.00897
				+1.0e-4 -0.153*(x[3]+x[4]+x[5]-4.0-global::sqrt8);
			break;
		case 98:
			*ret=-vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
				+0.00897
				+1.0e-4 -0.153*(x[3]+x[4]+x[5]-4.0-global::sqrt8);
			break;
		case 99: case 100:
			*ret=-vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
				+0.00897
				+1.0e-4 -0.153*(x[3]+x[4]+x[5]-4.0-global::sqrt8);
			break;
		case 8102: case 8103:
			*ret= -vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
				-0.0571 - 0.0624*(x[4]+x[5]-2.0*2.51) -0.115*(x[3]-2.);
				break;
		case 8104:
			*ret= -vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5])
				-0.0571 - 0.0624*(x[4]+x[5]-2.0*2.51) -0.115*(x[3]-2.);
				break;
		case 8101:
			*ret= -gamma(x[0],x[1],x[2],x[3],x[4],x[5])
				- 0.5*(vorVc(x[0],x[1],x[2],x[3],x[4],x[5])-
						vorVc(x[0],x[5],x[4],x[3],x[2],x[1]))
				-0.0538*(x[1]+x[2]+x[4]+x[5])
				-0.083*x[3]+0.59834; break;

		case 90: *ret=vorVc(x[0],x[1],x[2],x[3],x[4],x[5])-
					vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 89 : *ret=dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 891: *ret=tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-
					0.2529*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]); break;

		case 1968: *ret=- dih(x[0],x[1],x[2],x[3],x[4],x[5])
				-dih(x[0],x[6],x[7],x[8],x[9],x[10]);
				break;

		case 892: *ret=dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
				//-1.4412
				-1.6966
				-0.291*(x[0]-2.0)
				+0.586*(x[1]-2.0)
				+0.393*(x[2]-2.0)
				//-0.79*(x[3]-2.51)
				+0.397*(x[4]-2.51)
				+0.321*(x[5]-2.0);
				break;

		case 893: *ret=dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
				-1.072
				-0.115*(x[0]-2.0)
				+0.452*(x[1]-2.0)
				+0.452*(x[2]-2.0)
				-0.613*(x[3]-2.)
				+0.15*(x[4]-2.51)
				+0.15*(x[5]-2.51);
				break;

		case 894: *ret=dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
				//-1.394
				-0.115*(x[0]-2.0)
				+0.452*(x[1]-2.0)
				+0.452*(x[2]-2.0)
				//-0.618*(x[3]-2.51)
				+0.15*(x[4]-2.51)
				+0.15*(x[5]-2.51);
				break;

		case 895: *ret=dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
				//-1.8718
				-0.47*(x[0]-2.51)
				+0.522*(x[1]-2.0)
				+0.522*(x[2]-2.0)
				//-0.812*(x[3]-2.51)
				+0.522*(x[4]-2.0)
				+0.522*(x[5]-2.0);
				break;
		case 896: *ret=dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
				-2.5481
				+0.4*x[1]-0.15*x[0]+0.09*x[2]
				+0.631*x[4]-0.57*x[3]+0.234*x[5];
				break;
		case 897: *ret=dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
				-0.863
				-0.289*(x[0]-2.0)
				+1.36*(x[1]-2.51)
				+0.148*(x[2]-2.0)
				-0.688*(x[3]-2.)
				+1.36*(x[4]-2.51)
				+0.148*(x[5]-2.);
				break;
		case 898: *ret=dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
				-1.214
				-0.289*(x[0]-2.0)
				+0.723*(x[1]-2.51)
				+0.148*(x[2]-2.0)
				+0.723*(x[4]-2.51)
				+0.148*(x[5]-2.);
				break;
		case 899: *ret=solid(x[0],x[1],x[2],x[3],x[4],x[5])
				-0.8176
				+0.492*(x[0]-2.)
				+0.492*(x[1]+x[2]-4.)
				-0.43*(x[3]-2.0)
				-0.038*(x[4]+x[5]-2.*2.51);
				break;

		case 830: case 831 :*ret=-vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
				-0.0571
				-0.058*(x[0]-2.0)
				-0.105*(x[1]+x[2]-4.0)
				-0.115*(x[3]-2.0)
				-0.062*(x[4]+x[5]-2.*2.51);
				break;
		case 832: *ret=-vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5])
				-0.0571
				-0.058*(x[0]-2.0)
				-0.105*(x[1]+x[2]-4.0)
				-0.115*(x[3]-2.0)
				-0.062*(x[4]+x[5]-2.*2.51);
				break;

		case 88 : *ret= dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			+0.457*(x[1]+x[2]+x[4]+x[5]-8)-0.34*(x[0]-2)
			-1.9089; break;
		case -18 : *ret=cortau(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case -15 : *ret= 0.0268+vorVc(x[0],x[1],x[2],x[3],x[4],x[5])-
			gamma(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case -16 : case -17 :
			*ret= 0.0268+vorVc(x[0],x[1],x[2],x[3],x[4],x[5])-
			vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case -14 :
			*ret=dihedraly(x[1],x[0],x[2],x[4],x[3],x[5])+
				dihedraly(x[1],x[6],x[2],x[7],x[3],x[8])+
			dihedraly(x[6],x[1],x[2],x[3],x[7],x[8]);
			break;

		case -13: *ret = -x[7]; break;
		case -11: *ret = -dihedraly(x[0],x[1],1.255,1.6,1.6,2)+dihedraly(x[0],x[1],2.51,
					3.2,global::sqrt8,2); break;
		case -9 : *ret = -x[5]; break;
		case -7 : *ret = -u135M(x[0],x[1],x[2]); break;
		case -5 : *ret = solid(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case -3 : *ret = -corsigma(x[0],x[1],x[2],x[3],x[4],x[5])/global::pt; break;
		case -4 : *ret = tauVc(x[0],x[1],x[2],x[3],x[4],x[5])+
			tauVc(x[6],x[1],x[2],x[3],x[7],x[8]); break;
		case -1 : *ret = x[1]*x[1]+x[2]*x[2];
			break;
		case -2 : *ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]);
			break;

		case 1200 :
			{
			double t = crossdiag(x);
			*ret = cortau(x[0],x[1],x[2],x[3],x[4],x[5])+
				cortau(x[1],x[0],x[6],t,x[8],x[5])
				-areadiff(x[0],x[1],x[3])*2.0*(global::zetapt - phi(1.255,1.255));
			/*
			//*ret = -corsigma(x[0],x[1],x[2],x[3],x[4],x[5])
			//	-corsigma(x[1],x[0],x[6],t,x[8],x[5]);
			//*ret = tauVc(x[0],x[1],x[2],x[3],x[4],x[5])+
				//tauVc(x[6],x[1],x[2],x[3],x[7],x[8]);
			*ret = -vorVc(x[0],x[1],x[2],x[3],x[4],x[5])+
				-vorVc(x[6],x[1],x[2],x[3],x[7],x[8]);
			*/
			}
			break;


		case 900 : *ret = tau(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 901 : 
		case 902 : *ret = tau_analytic(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 903 : 
		case 904 : *ret = tauVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;

		case 800 : *ret = 0.0
			+ tau(x[0],x[1],x[2],x[3],x[4],x[5]) 
			+ tau(x[0],x[1],x[8],x[6],x[7],x[5])
			+ tau(x[0],x[8],x[10],x[9],x[11],x[7])
				+ tau(x[0],x[10],x[14],x[12],x[13],x[11])
				+ 0*tau(x[0],x[14],x[2],x[15],x[4],x[13])
			;
			break;

		case 801 : *ret = 0.0
			+ tau(x[0],x[1],x[2],x[3],x[4],x[5]) 
			+ tau(x[0],x[1],x[8],x[6],x[7],x[5])
			+ tau(x[0],x[8],x[10],x[9],x[11],x[7])
				+ tau(x[0],x[10],x[14],x[12],x[13],x[11])
				+ tau(x[0],x[14],x[2],x[15],x[4],x[13])
			;
			break;

		case 802 : case 803 : *ret = 0.0
			+ tau(x[0],x[1],x[2],x[3],x[4],x[5]) 
			+ tau(x[0],x[1],x[8],x[6],x[7],x[5])
			+ tau(x[0],x[8],x[10],x[9],x[11],x[7])
				+ tau(x[0],x[10],x[14],x[12],x[13],x[11])
				+ tau_analytic(x[0],x[14],x[2],x[15],x[4],x[13])
			;
			break;

		case 804 : case 805 : case 806 : *ret = 0.0
			+ tau(x[0],x[1],x[2],x[3],x[4],x[5]) 
			+ tau(x[0],x[1],x[8],x[6],x[7],x[5])
			+ tau(x[0],x[8],x[10],x[9],x[11],x[7])
				+ tau(x[0],x[10],x[14],x[12],x[13],x[11])
				+ tauVc(x[0],x[14],x[2],x[15],x[4],x[13])
			;
			break;

		case 6 : *ret = vorVc(x[0],x[1],x[2],x[3],x[4],x[5]) -
					gamma(x[0],x[1],x[2],x[3],x[4],x[5]);
			break;
		case 7 : *ret = vorVc(x[0],x[1],x[2],x[3],x[4],x[5])-
				vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5]);
			break;

		// INEQUALITIES FOR FLAT QUARTERS INSIDE A QUAD CLUSTER...10--18
		case 10 : *ret = 0.0123 
				- 0.15*x[0] + 0.35*x[1] - 0.15*x[2] - 0.17*x[3] + 0.7022*x[4]
				-dihedraly(x[1],x[0],x[2],x[4],x[3],x[5]);
				break;
		case 11 : *ret = -2.63363 + 0.631*x[0] - 0.13*x[1] + 0.31*x[2] + 
				0.413*x[3] - 0.58*x[4] + 0.025*x[5]
				+dihedraly(x[1],x[0],x[2],x[4],x[3],x[5]);
				break;
		case 12 : *ret = -0.3482 + 0.714*x[0] - 0.221*x[1] - 0.221*x[2] + 
				0.92*x[3] - 0.221*x[4] - 0.221*x[5]
				-dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]);
				break;
		case 13 : *ret = -2.37095 - 0.315*x[0] + 0.3972*x[1] + 0.3972*x[2] - 
				0.715*x[3] + 0.3972*x[4] + 0.3972*x[5]
				+dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]);
				break;
		case 14 : *ret = -0.437235 - 0.187*x[0] - 0.187*x[1] 
				- 0.187*x[2] + 0.1185*x[3] + 0.479*x[4] + 0.479*x[5]
				-solid(x[0],x[1],x[2],x[3],x[4],x[5]);
				break;
		case 15 : *ret = -2.244 +0.488*x[0] + 0.488*x[1] + 0.488*x[2] - 
				0.334*x[4] - 0.334*x[5]
				+solid(x[0],x[1],x[2],x[3],x[4],x[5]);
				break;
		case 16 : *ret = 1.17401 - 0.159*x[0] - 0.081*x[1] - 0.081*x[2] - 
				0.133*x[4] - 0.133*x[5]
				-gamma(x[0],x[1],x[2],x[3],x[4],x[5]);
				break;
		case 17 : 
		case 18 : *ret = 1.17401 - 0.159*x[0] - 0.081*x[1] - 0.081*x[2] - 
				0.133*x[4] - 0.133*x[5]
				-vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5]);
				break;
		// END OF FLAT QUARTERS 

		case 49 : /*experimental flat*/ *ret = 
				1.0e-6 -0.197*(x[3]+x[4]+x[5]-4-global::sqrt8)
				-gamma(x[0],x[1],x[2],x[3],x[4],x[5]);
				
				break;

		case 50 : 
		case 51 : /*experimental flat*/ *ret = 
				1.0e-6 -0.197*(x[3]+x[4]+x[5]-4-global::sqrt8)
				-vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5]);
				break;

		case 53 : /* experimental qrtet*/ *ret = 
				gamma(2.0,2.0,2.0,2.51,x[3]+x[4]+x[5]-4.51,2.0) 
					-gamma(x[0],x[1],x[2],x[3],x[4],x[5]);
				break;

		// UPRIGHT QUARTERS IN QUADS
		case 20 : *ret = -1.82419 - 0.636*x[0] + 0.462*x[1] + 0.462*x[2] - 
				0.82*x[3] + 0.462*x[4] + 0.462*x[5]
				+dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]);
				break;
		case 21 : *ret = -0.75281 + 0.55*x[0] - 0.214*x[1] - 0.214*x[2] + 
				1.24*x[3] - 0.214*x[4] - 0.214*x[5]
				-dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]);
				break;
		case 22 : *ret = -2.5481 + 0.4*x[0] - 0.15*x[1] + 0.09*x[2] + 
				0.631*x[3] - 0.57*x[4] + 0.23*x[5]
				+dihedraly(x[1],x[0],x[2],x[4],x[3],x[5]);
				break;
		case 23 : *ret = 0.3429 - 0.454*x[0] + 0.34*x[1] + 0.154*x[2] - 
				0.346*x[3] + 0.805*x[4]
				-dihedraly(x[1],x[0],x[2],x[4],x[3],x[5]);
				break;
		case 24 : *ret = -0.2618 + 0.065*x[1] + 0.065*x[2] + 
				0.061*x[3] - 0.115*x[4] - 0.115*x[5]
				+ solid(x[0],x[1],x[2],x[3],x[4],x[5]);
				break;
		case 25 : *ret = -0.2514 - 0.293*x[0] - 0.03*x[1] - 0.03*x[2] + 
				0.12*x[3] + 0.325*x[4] + 0.325*x[5]
				- solid(x[0],x[1],x[2],x[3],x[4],x[5]);
				break;

		case 26 : 
		case 27 : *ret = 0.59834 - 0.054*x[1] - 0.054*x[2] - 0.083*x[3] - 
				0.054*x[4] - 0.054*x[5]
				- octavor(x[0],x[1],x[2],x[3],x[4],x[5]);
				break;
		case 28 : *ret = 0.59834 - 0.054*x[1] - 0.054*x[2] - 0.083*x[3] - 
				0.054*x[4] - 0.054*x[5]
				- gamma(x[0],x[1],x[2],x[3],x[4],x[5]);
				break;
		case 261 : 
		case 271 : 
				*ret = - 0.133*(x[1]+x[2]+x[4]+x[5]-8) 
				- 0.135*(x[3]-2.0)
				+ 0.07*(x[0]-2.51)
				- octavor(x[0],x[1],x[2],x[3],x[4],x[5]);
				break;
		case 281 : 
				*ret = 0 - 0.133*(x[1]+x[2]+x[4]+x[5]-8) 
				- 0.135*(x[3]-2.0)
				+ 0.07*(x[0]-2.51)
				- gamma(x[0],x[1],x[2],x[3],x[4],x[5]);
				break;

		// END OF UPRIGHT QUARTERS.

		case 69 : *ret = /*experimental*/
				-0.419351*solid(x[0],x[1],x[2],x[3],x[4],x[5])
				+0.079431*dihedraly(x[1],x[0],x[2],x[4],x[3],x[5])
				+0.06904
				-0.0846*(x[0]-2.8)
				-gamma(x[0],x[1],x[2],x[3],x[4],x[5]);
				break;

		case 70 : case 71 :
				*ret = /*experimental*/
				-0.419351*solid(x[0],x[1],x[2],x[3],x[4],x[5])
				+0.06904
				-0.0846*(x[0]-2.8)
				+0.079431*dihedraly(x[1],x[0],x[2],x[4],x[3],x[5])
				-octavor(x[0],x[1],x[2],x[3],x[4],x[5]);
				break;

		// QUAD VC 
		case 30 : *ret = -4.885 - 0.372*x[0] + 0.465*x[1] + 0.465*x[2] + 
				0.465*x[4] + 0.465*x[5]
				+dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]);
				break;
		case 31: *ret = 0.9978 -0.06*x[1] - 0.06*x[2] - 0.185*x[4] - 0.185*x[5]
				-vorVc(x[0],x[1],x[2],x[3],x[4],x[5]);
				break;
		case 32: *ret = -0.06333 + quoinH(x[0],x[1],x[2])+
				0.00758*x[0]+0.0115*x[1]+0.0115*x[2];
				break;
		// END QUAD VC 

		// QRTET edge inequalities.
		case 200 :
				*ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
					- 1.23095 - 0.283*(x[4]+x[5]-4.0);
				break;

		case 201 :
				*ret = solid(x[0],x[1],x[2],x[3],x[4],x[5])
				-0.551285 - 0.221*(x[3]+x[4]+x[5]-6) 
				+0.377076*(x[0]+x[1]+x[2]-6);
				break;
		case 202 :
				*ret = -solid(x[0],x[1],x[2],x[3],x[4],x[5])
				+0.55778 + 0.221*(x[3]+x[4]+x[5]-6) 
				;
				break;
		case 203 :
				*ret = -dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
				+1.23116
				+0.731*(x[3]-2)+0.498*(x[0]-2)-0.212*(x[4]+x[5]-4)
				;
				break;
		case 204 : 
				*ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
				-1.23095
				+0.34*(x[1]+x[2]-4)
				+0.27*(x[4]+x[5]-4)
				-0.689*(x[3]-2)
				;
				break;
		case 205 :
				*ret = -gamma(x[0],x[1],x[2],x[3],x[4],x[5])
				+0.0553737
				-0.109*(x[0]+x[1]+x[2]-6)
				-0.14135*(x[3]+x[4]+x[5]-6)
				;
				break;
		case 206 :
				*ret = -gamma(x[0],x[1],x[2],x[3],x[4],x[5])
				-0.419351*solid(x[0],x[1],x[2],x[3],x[4],x[5])
				+0.28665
				-0.2*(x[0]+x[1]+x[2]-6)
				-0.048*(x[3]+x[4]+x[5]-6)
				;
				break;
		case 207 :
				*ret = tau(x[0],x[1],x[2],x[3],x[4],x[5])
				+1.0e-6 
				-0.0845696*(x[0]+x[1]+x[2]-6)
				-0.163*(x[3]+x[4]+x[5]-6)
				;
				break;

		case 208 :
				*ret = solid(x[0],x[1],x[2],x[3],x[4],x[5])
				-0.60657
				-0.1781*(x[3]+x[4]+x[5]-6.25)
				+0.378*(x[0]+x[1]+x[2]-6)
				;
				break;

		case 2091 :
				*ret = -solid(x[0],x[1],x[2],x[3],x[4],x[5])
				+0.61298
				+0.28*(x[3]+x[4]+x[5]-6.25)
				-0.254*(x[0]+x[1]+x[2]-6)
				;
				break;

		case 209 :
				*ret = -solid(x[0],x[1],x[2],x[3],x[4],x[5])
				+0.61298
				+0.3405*(x[3]+x[4]+x[5]-6.25)
				-0.171*(x[0]+x[1]+x[2]-6)
				;
				break;

		case 210 :
				*ret = -gamma(x[0],x[1],x[2],x[3],x[4],x[5])
				+0.02004
				-0.0781*(x[3]+x[4]+x[5]-6.25)
				-0.1208*(x[0]+x[1]+x[2]-6)
				;
				break;

		case 2101 :
				*ret = -vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5])
				+0.02004
				-0.0781*(x[3]+x[4]+x[5]-6.25)
				-0.167*(x[0]+x[1]+x[2]-6)
				;
				break;

		case 211 :
				*ret = -gamma(x[0],x[1],x[2],x[3],x[4],x[5])
				-0.419351*solid(x[0],x[1],x[2],x[3],x[4],x[5])
				+0.27441
				+0.0106*(x[3]+x[4]+x[5]-6.25)
				-0.2*(x[0]+x[1]+x[2]-6)
				;
				break;
		case 212 :
				*ret = -vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5])
				-0.419351*solid(x[0],x[1],x[2],x[3],x[4],x[5])
				+0.27441
				+0.0106*(x[3]+x[4]+x[5]-6.25)
				-0.2*(x[0]+x[1]+x[2]-6)
				;
				break;

		// (4,1) stuff:
		case 40 : *ret = tau(x[0],x[1],x[2],x[3],x[4],x[5])-
					tau(x[0],x[1],x[2],2.0,x[4],x[5]);
			break;
		case 41 : *ret = tau(x[0],x[1],x[2],x[3],x[4],x[5])-
					tau(2.0,x[1],x[2],x[3],x[4],x[5]);
			break;
		case 42 : case 43 :
			*ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-
					dihedraly(2.0,x[1],x[2],x[3],x[4],x[5]);
			break;

		case 44 :
			*ret = tau(x[0],x[1],x[2],x[3],x[4],x[5])
			-2.7209127*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			//-(tau(x[0],x[1]+x[5]-2.0,x[2]+x[4]-2.0,x[3],2.0,2.0)
			//-2.7209127*dihedraly(x[0],x[1]+x[5]-2.0,x[2]+x[4]-2.0,x[3],2.0,2.0));
			-(tau(x[0],x[1],x[2],x[3],2.0,2.0)
			-2.7209127*dihedraly(x[0],x[1],x[2],x[3],2.0,2.0));
			break;
		case 45 :
			*ret = taudih(x[0],x[1],x[2],x[3],x[4],x[5])-
				taudih(x[0],x[1],x[2],x[3],x[4],2.0);
			break;
		case 46 :
			*ret = 
				+xix*dihedraly(2.0,x[1],x[2],global::sqrt8,x[3],x[0])+
				taudih(2.0,x[2],x[4],2,2,x[3])
				-(xix*dihedraly(2.0,x[1],x[2],global::sqrt8,2.0,x[0])+
				taudih(2.0,x[2],x[4],2,2,2.0));
			break;
		case 47 :
			*ret = 
				+xix*dihedraly(2.0,x[0],x[1],global::sqrt8,2,2)
				+taudih(2,x[1],x[2],2,2,2)
				+taudih(2,x[2],x[3],2,2,2)
				+taudih(2,x[3],x[4],2,2,2)
				+taudih(2,x[4],x[0],2,2,2);
		break;
		

		case 101: *ret = tauXg(x[0],x[1],x[2],x[3],x[4],x[5])
				-1.01104 
				+ 0.78701*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]);
				break;

		case 102: *ret = tauXo(x[0],x[1],x[2],x[3],x[4],x[5])
				-1.01104 
				+ 0.78701*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]);
				break;

		case 103: *ret = -sigXg(x[0],x[1],x[2],x[3],x[4],x[5])
				-0.9871
				+ 0.80449*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]);
				break;

		case 104: *ret = -sigXo(x[0],x[1],x[2],x[3],x[4],x[5])
				-0.9871
				+ 0.80449*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]);
				break;
	
		case 11103: *ret = -sigXg(x[0],x[1],x[2],x[3],x[4],x[5])
				-0.3429 +0.24573*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]);
				break;
		case 11104: *ret = -sigXo(x[0],x[1],x[2],x[3],x[4],x[5])
				-0.3429 +0.24573*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]);
				break;
		case 11105: 
			*ret = -1.08*global::pt - vor_analytic(x[0],x[1],x[2],x[3],x[4],
			x[5]); break;
		case 11106: case 11107: case 11108:
			*ret = -1.08*global::pt - vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); 
			break;

		case 300: *ret = -fA(x); break;

		case 400: *ret = -x[3]; break;

		case 500: *ret = V(x[0],x[1],x[2],x[3],x[4],x[5]) 
				-kappa(x[0],x[2],x[7],x[6],x[8],x[4])+
				V(x[0],x[7],x[11],x[9],x[10],x[8])+
				V(x[0],x[11],x[1],x[12],x[5],x[10]); break;

		case 501: *ret = V(x[0],x[1],x[2],x[3],x[4],x[5]) 
				-kappa(x[0],x[2],x[7],x[6],x[8],x[4])+
				VV(x[0],x[7],x[11],x[9],x[10],x[8])+
				V(x[0],x[11],x[1],x[12],x[5],x[10]); break;

		// Appendix on contraction (IV).

		// embedded pentagon stuff:
		case 700: *ret = tau(x[0],x[1],x[2],x[3],x[4],x[5])-
				0.2529*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]);
				break;
		case 701: *ret = tau_analytic(x[0],x[1],x[2],x[3],x[4],x[5])-
				0.2529*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]);
				break;
		case 702: case 703 :
				*ret = tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-
				0.2529*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]);
				break;
		case 704: 
				*ret = tau(x[0],x[1],x[2],x[3],x[4],x[5])-
				0.2529*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]);
				break;
		case 705: 
				*ret = tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-
				0.2529*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]);
				break;
		case 706: case 707 :
				*ret = tau_analytic(x[0],x[1],x[2],x[3],x[4],x[5])-
				0.2529*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]);
				break;
		case 708: 
				*ret = tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-
				0.2529*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]);
				break;



		case 507: *ret = -dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]);
				break;

		case 503: *ret= -delta(x[0]*x[0],x[1]*x[1],x[2]*x[2],x[3]*x[3],x[4]*x[4],x[5]*x[5]);
				break;
		case 504: *ret= -U(x[0]*x[0],x[2]*x[2],x[4]*x[4]);
				break;

		default : cout << "generic default" << endl << flush;
			*ret=0;
		}
	}

iter::~iter() 
	{ delete[] xmin; 
	  delete[] xmax; 
		delete[] x; }


static void ConstraintPage1(int numargs,int whichFn,double* x,double* ret,void*)
    {
	*ret = 0;
	switch (INEQ_NUMBER) {

		case 246:
			*ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])- 2.588; break;
		case 245:
			*ret = x[1]+x[2]+x[4]+x[5]-8.57; break;
		case 248:
			*ret =x[3]-crossdiag(x); break;

		
		case 240: case 241:
			*ret=-crossdiag(x) + x[3]; break;
		case 243:
			*ret=x[1]+x[2]-4.6; break;
		case 242:
			switch(whichFn)
			{
			case 1:
			*ret= dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.743; break;
			case 2:
			*ret= 3.5-crossdiag(x); break;
			}
			break;
		case 238: 
			switch(whichFn)
			{
			case 1:
			*ret= 0.000001-delta(x[0]*x[0],x[1]*x[1],x[2]*x[2],x[3]*x[3],x[4]*x[4],x[5]*x[5]);
			break;
			case 2: *ret=D2Vee(x[0],x[1],x[2],x[3],x[4],x[5]); break;
			}
			break;

		case 239:
			switch(whichFn)
			{
			case 1:
			*ret= 0.000001-delta(x[0]*x[0],x[1]*x[1],x[2]*x[2],x[3]*x[3],x[4]*x[4],x[5]*x[5]);
			break;
			case 2: *ret=D2Vee1(x[0],x[1],x[2],x[3],x[4],x[5]); break;
			}
			break;

		case 91: case 96 :
			switch(whichFn) {
			case 1: *ret=radf(x[3],x[4],x[5])-global::sqrt2; break;
			case 2: *ret=radf(x[1],x[2],x[3])-global::sqrt2; break;
			}
			break;
		case 11105 : case 11106 :
			*ret=-2.2+dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 11107:
			switch(whichFn) {
			case 1:
			*ret=global::sqrt2-radf(x[3],x[4],x[5]); break;
			case 2:
			*ret=-2.2+dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]); break;
			}
			break;
		case 11108:
			switch(whichFn) {
			case 1:
			*ret=global::sqrt2-radf(x[3],x[2],x[1]); break;
			case 2:
			*ret=-2.2+dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]); break;
			}
			break;
		
		case 94: case 99: 
			*ret=global::sqrt2-radf(x[3],x[4],x[5]); break;
		case 95: case 100: 
			*ret=global::sqrt2-radf(x[3],x[2],x[1]); break;
		case 9101:
			*ret=global::sqrt2-radf(x[0],x[2],x[4]); break;

		case 90: switch(whichFn)
			{
			case 1 : *ret=global::sqrt2-radf(x[3],x[2],x[1]); break;
			}
		break;
			
		case -14:
			{
			double t = crossdiag(x);
			switch(whichFn) {
			case 1 : *ret= -delta(x[0]*x[0],x[1]*x[1],x[2]*x[2],x[3]*x[3],x[4]*x[4],x[5]*x[5]);
				break;
			case 2 : *ret= -delta(x[1]*x[1],x[0]*x[0],x[6]*x[6],t*t,x[8]*x[8],x[5]*x[5]);
				break;
			case 3 : *ret= -delta(x[6]*x[6],x[1]*x[1],x[2]*x[2],x[3]*x[3],x[7]*x[7],x[8]*x[8]);
				break;
			case 4 : *ret= 3.2/*x[3]*/-t; 
				break;
			case 5 : *ret= 0*(-vorVc(x[6],x[1],x[2],x[3],x[7],x[8])); break;

			case 6 :
			{
			double dA = safesqrt(delta(x[0]*x[0],x[1]*x[1],x[2]*x[2],x[3]*x[3],
					x[4]*x[4],x[5]*x[5]));
			double dB = safesqrt(delta(x[1]*x[1],x[6]*x[6],x[2]*x[2],x[7]*x[7],
					x[3]*x[3],x[8]*x[8]));
			*ret = 
				-Vee0(x[1],x[0],x[2],x[4],x[3],x[5])*dB
				-Vee0(x[1],x[6],x[2],x[7],x[3],x[8])*dA;
			}
			break;
			}
			}
			break;

		case -13: 
			{
			double t = crossdiag(x);
			switch(whichFn) {
			case 1 : *ret= -overlap(x[1],x[2],x[3]); break;
			case 2 : *ret= -overlap(x[0],x[6],t); break;
			case 3 : *ret= -delta(x[0]*x[0],x[1]*x[1],x[2]*x[2],x[3]*x[3],x[4]*x[4],x[5]*x[5]);
				break;
			case 4 : *ret= -delta(x[1]*x[1],x[0]*x[0],x[6]*x[6],t*t,x[8]*x[8],x[5]*x[5]);
				break;
			case 5 : *ret= -vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
			case 6 : *ret= -vorVc(x[1],x[0],x[6],t,x[8],x[5]); break;
			}
			}
			break;
		case -9 : *ret = 2.46/2.0 - dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case -7 : *ret = -1.255 + radf(x[0],x[1],x[2]); break;
		case -5 : 
			*ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-2.89;
				break;
			
		case 7 : switch(whichFn) {
			case 2 : *ret = -radf(x[3],x[4],x[5])+global::sqrt2; break;
			case 1 : *ret = -radf(x[1],x[2],x[3])+global::sqrt2; break;
			}
			break;

		// QUAD STUFF
		case 49 : case -15 :
		case 16 : switch(whichFn) {
			case 1 : *ret = +radf(x[3],x[4],x[5])-global::sqrt2; break;
			case 2 : *ret = +radf(x[1],x[2],x[3])-global::sqrt2; break;
			}
			break;
		case 52 : *ret = x[3]+x[4]+x[5]- 6.51; break;
		case 53 : switch(whichFn) {
			case 1 : *ret = 6.51-x[3]-x[4]-x[5]; break;
			case 2 : *ret = x[3]+x[4]+x[5]-7.02; break;
				}
			break;
		case 26 : case 261 :
		case 70 : *ret = -radf(x[0],x[1],x[2])+global::sqrt2; break;
		case 27 : case 271 :
		case 71 : *ret = -radf(x[0],x[2],x[4])+global::sqrt2; break;
		case 28 : case 281 :
		case 69 : switch(whichFn) {
			case 1 : *ret = +radf(x[0],x[1],x[2])-global::sqrt2; break;
			case 2 : *ret = +radf(x[0],x[2],x[5])-global::sqrt2; break;
			}
			break;
 
		case 17 : case 50 : case -16 :
			*ret = -radf(x[3],x[4],x[5])+global::sqrt2; break;
		case 18 : case 51 : case -17 :
			*ret = -radf(x[3],x[2],x[1])+global::sqrt2; break;

		case 31 : *ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-2.12; break;
		// END QUAD STUFF

		// QRTET edge inequalities. 
		case 200 : *ret = x[0]+x[1]+x[2]-6.13; break;
		case 201 : case 202 : case 203 : 
				*ret = x[3]+x[4]+x[5]-6.25; break;
		case 204 :  case 205 : case 206 : case 207 :
			switch(whichFn) {
			case 1 : *ret = x[3]+x[4]+x[5]-6.25; break;
			case 2 : *ret = 0*(x[0]+x[1]+x[2]-6.13); break;
			}
			break;
		case 208 :  case 209 : case 210 : 
			switch(whichFn) {
			case 1 : *ret = 6.25-(x[3]+x[4]+x[5]); break;
			case 2 : *ret = 0*(x[0]+x[1]+x[2]-6.13); break;
			}
			break ;
		case 2091 :
			switch(whichFn) {
			case 1 : *ret = 6.22-(x[3]+x[4]+x[5]); break;
			case 2 : *ret = 0*(x[0]+x[1]+x[2]-6.13); break;
			}
			break ;

		case 211 : 
			switch(whichFn) {
			case 1 : *ret = 6.25-(x[3]+x[4]+x[5]); break;
			case 2 : *ret = rady(x[0],x[1],x[2],x[3],x[4],x[5])-1.41; break;
			case 3 : *ret = 0*(x[0]+x[1]+x[2]-6.13); break;
			}
			break;
		case 212 : case 2101 :
			switch(whichFn) {
			case 1 : *ret = 6.25-(x[3]+x[4]+x[5]); break;
			case 2 : *ret = -rady(x[0],x[1],x[2],x[3],x[4],x[5])+1.41; break;
			case 3 : *ret = x[0]+x[1]+x[2]-6.13; break;
			}
			break;

		case 102 : 
		case 104 : case 11104 : *ret = -radf(x[0],x[1],x[2])+global::sqrt2; break;
		case 300 :
		case 8103: case 831:
			*ret=global::sqrt2-radf(x[3],x[4],x[5]); break;
		case 8104: case 832:
			*ret=radf(x[3],x[4],x[5])-global::sqrt2; break;
		case 8101 : case 103 : case 11103: switch(whichFn) {
			case 1 : *ret = +radf(x[0],x[1],x[2])-global::sqrt2; break;
			case 2 : *ret = +radf(x[0],x[2],x[5])-global::sqrt2; break;
			}
			break;
		case 400 :
			*ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-2.12;
			break;

		case 500 :
			switch (whichFn) {
			case 1 : *ret = +radf(x[0],x[1],x[5])-global::sqrt2; break;
			case 2 : *ret = +radf(x[0],x[2],x[4])-global::sqrt2; break;
			case 3 : *ret = +radf(x[0],x[7],x[8])-global::sqrt2; break;
			case 4 : *ret = +radf(x[0],x[10],x[11])-global::sqrt2; break;
			case 5: *ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]) +
				dihedraly(x[0],x[2],x[7],x[6],x[8],x[4])+
				dihedraly(x[0],x[7],x[11],x[9],x[10],x[8])+
				dihedraly(x[0],x[11],x[1],x[12],x[5],x[10]) - global::pi*2.0; break;
			}
			break;

		case 501 :
			switch (whichFn) {
			case 1 : *ret = +radf(x[0],x[1],x[5])-global::sqrt2; break;
			case 2 : *ret = +radf(x[0],x[2],x[4])-global::sqrt2; break;
			case 3 : *ret = -radf(x[0],x[7],x[8])+global::sqrt2; break;
			case 4 : *ret = +radf(x[0],x[10],x[11])-global::sqrt2; break;
			case 5: *ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]) +
				dihedraly(x[0],x[2],x[7],x[6],x[8],x[4])+
				dihedraly(x[0],x[7],x[11],x[9],x[10],x[8])+
				dihedraly(x[0],x[11],x[1],x[12],x[5],x[10]) - global::pi*2.0; break;
			}
			break;


		case 700:
			switch(whichFn) {
				case 1 : *ret=1.51-dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]); break;
				case 2 : *ret=circum2(x[0],x[1],x[2],x[3],x[4],x[5])-1.41*1.41;
					break;
				}
			break;
		case 701:
			switch(whichFn) {
				case 1 : *ret=1.51-dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]); break;
				case 2 : *ret=-circum2(x[0],x[1],x[2],x[3],x[4],x[5])+1.41*1.41;
					break;
				}
			break;
		case 702:
			switch(whichFn) {
				case 1 : *ret=1.26-dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]); break;
				case 2 : *ret=-1.63+dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]); break;
				}
			break;
		case 703:
			switch(whichFn) {
				case 1 : *ret=1.14-dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]); break;
				case 2 : *ret=-1.51+dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]); break;
				}
			break;

		case 704:
			switch(whichFn) {
				case 1 : *ret=1.26-dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]); break;
				case 2 : *ret=-1.63+dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]); break;
				case 3 : *ret=-global::sqrt2+radf(x[0],x[1],x[5]); break;
				case 4 : *ret=-global::sqrt2+radf(x[3],x[4],x[5]); break;
				}
			break;

		case 705:
			switch(whichFn) {
				case 1 : *ret=1.26-dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]); break;
				case 2 : *ret=-1.63+dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]); break;
				}
			break;

		case 706: case 708 :
			switch(whichFn) {
				case 1 : *ret=1.26-dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]); break;
				case 2 : *ret=-1.63+dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]); break;
				case 3 : *ret=global::sqrt2-radf(x[3],x[4],x[5]); break;
				}
			break;

		case 707:
			switch(whichFn) {
				case 1 : *ret=1.26-dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]); break;
				case 2 : *ret=-1.63+dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]); break;
				case 3 : *ret=global::sqrt2-radf(x[0],x[1],x[5]); break;
				}
			break;

		//case 503 : *ret= 0.00758*x[0]+0.0115*(x[2]+x[4])-0.06333; break;
		//case 504 : *ret= 0.00758*x[0]+0.0115*(x[2]+x[4])-0.06333; break;
		case 503 : case 504 : *ret= 
				-1.255 + radf(x[0],x[1],x[4]); break;

		case 44 : *ret = -0.0830604685026962//-1.5pt 
			+tau(x[0],x[1],x[2],x[3],x[4],x[5])+
			tau(2.0,x[1],2,2,2,x[5])+tau(2,2,x[2],2,x[4],2); break;
		case 45 : *ret = -0.0830604685026962//-1.5pt 
			+tau(x[0],x[1],x[2],x[3],x[4],x[5])+
			tau(2.0,x[1],2,2,2,x[5]); break;
		case 46 : *ret = -0.0830604685026962//-1.5pt
			+tau(2,2,x[1],2,x[0],2)
			+tau(2,x[2],x[4],2,2,x[3])
			+tau(2,2,x[4],2,2,2);
			break;
		case 47 : *ret =  -0.0830604685026962//-1.5pt
			+tau(2,x[1],x[2],2,2,2)
			+tau(2,x[2],x[3],2,2,2)
			+tau(2,x[3],x[4],2,2,2)
			+tau(2,x[4],x[0],2,2,2);
			break;

		case 2 : *ret = -radf(x[2],x[13],x[15])+global::sqrt2; break;

		case 1200 : switch(whichFn) {
			case 1 :
			*ret= -delta(x[0]*x[0],x[1]*x[1],x[2]*x[2],
									x[3]*x[3],x[4]*x[4],x[5]*x[5]);
				break;
			case 2 : 
			{
			double t = crossdiag(x);
			*ret= -delta(x[1]*x[1],x[0]*x[0],x[6]*x[6],
									t*t,x[8]*x[8],x[5]*x[5]);
			}
				break;
			case 3 : *ret= 3.2-crossdiag(x); break;

		//	case 4 : *ret= -x[3]*x[3]+x[1]*x[1]+x[2]*x[2]; break;

//			case 4 : *ret=-0.0498362811016177 + tauVc(x[0],x[1],x[2],x[3],x[4],x[5]);
			//case 4 : *ret=0.0110747 -vorVc(x[0],x[1],x[2],x[3],x[4],x[5]);
						// = 0.9 pt. // 0.2 pt.
				break;
			}
		break;


		case 999 :
			*ret = -3.14159 + dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])+ dihedraly(x[0],x[7],x[2],x[6],x[4],x[8])+ dihedraly(x[0],x[7],x[11],x[9],x[10],x[8]);
				break;

		case 1100 : case 1101 : case 1102 : case 1103 : case 1104 : case 1105 :
		case 1106 : case 1107 : case 1108 :

		case 1000: case 1001 : case 1002 : case 1003 : case 1004 : case 1005 :
		case 1006: case 1007 : case 1008 :
			*ret= -delta(x[0]*x[0],x[1]*x[1],x[2]*x[2],
									x[3]*x[3],x[4]*x[4],x[5]*x[5]);
			break;

		case 1010 : case 1011 : case 1012 : case 1013 : case 1014 : case 1015 :
		case 1016 : case 1017 : case 1018 : case 1019 : case 1020 : case 1021:
		case 1022 : case 1023 : case 1024 : case 1025 : case 1026 : case 1027:

		case 1028 : case 1029 : case 1030 : case 1031 : case 1032 : case 1033 :
		case 1034 : case 1035 : case 1036 : case 1037 : case 1038 : case 1039 :
		case 1040 : case 1041 : case 1042 : case 1043 : case 1044 : case 1045 :

		switch (whichFn)
			{
			case 1 :
			*ret= -delta(x[0]*x[0],x[1]*x[1],x[2]*x[2],
									x[3]*x[3],x[4]*x[4],x[5]*x[5]);
					break;
			case 2 :
			*ret= -delta(x[0]*x[0],x[7]*x[7],x[2]*x[2],
									x[6]*x[6],x[4]*x[4],x[8]*x[8]);
					break;
			}
			break;

		case 900 : 
		case 903 :
		case 904 : switch(whichFn)
			{
			case 1 : *ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.32;
				break;
			}
			break;

		case 901 :  switch(whichFn)
			{
			case 1 : *ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.32;
				break;
			case 2 : *ret = -radf(x[3],x[4],x[5])+global::sqrt2; break;
			}
			break;

		case 902 :  switch(whichFn)
			{
			case 1 : *ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.32;
				break;
			case 2 : *ret = -radf(x[3],x[2],x[1])+global::sqrt2; break;
			}
			break;


		case 800 : case 801 : case 804 : case 806 : switch(whichFn) {
			case 1 : *ret = 
			dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]) 
			+ dihedraly(x[0],x[1],x[8],x[6],x[7],x[5])
			+ dihedraly(x[0],x[8],x[10],x[9],x[11],x[7])
				+ dihedraly(x[0],x[10],x[14],x[12],x[13],x[11])
				+ dihedraly(x[0],x[14],x[2],x[15],x[4],x[13])
			 - global::pi*2.0;
				break;
			}
			break;

		case 802 : case 805 : switch(whichFn) {
			case 1 : *ret = 
			dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]) 
			+ dihedraly(x[0],x[1],x[8],x[6],x[7],x[5])
			+ dihedraly(x[0],x[8],x[10],x[9],x[11],x[7])
				+ dihedraly(x[0],x[10],x[14],x[12],x[13],x[11])
				+ dihedraly(x[0],x[14],x[2],x[15],x[4],x[13])
			 - global::pi*2.0;
				break;
			case 2 : *ret = -radf(x[4],x[13],x[15])+global::sqrt2; break;
			}
			break;

		case 803 : switch(whichFn) {
			case 1 : *ret = 
			dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]) 
			+ dihedraly(x[0],x[1],x[8],x[6],x[7],x[5])
			+ dihedraly(x[0],x[8],x[10],x[9],x[11],x[7])
				+ dihedraly(x[0],x[10],x[14],x[12],x[13],x[11])
				+ dihedraly(x[0],x[14],x[2],x[15],x[4],x[13])
			 - global::pi*2.0;
				break;
			case 2 : *ret = -radf(x[2],x[14],x[15])+global::sqrt2; break;
			}
			break;

		

		case -1 : *ret = -radf(x[0],x[1],x[2])+radf(x[0],2,2.51); break;

		case -4 : switch(whichFn) {
			case 1 : *ret= 
		 -delta(x[0]*x[0],x[1]*x[1],x[2]*x[2],x[3]*x[3],x[4]*x[4],x[5]*x[5]);
		 break;
			case 2 : *ret= global::sqrt8-crossdiag(x); break;
			case 7 : *ret= global::sqrt8-crossdiagNew(x[0],x[1],x[2],x[3],x[4],
				x[5],x[6],x[7],x[8]); break;
			}
			break;

		case 1 : switch(whichFn) {
			case 1 :*ret = -dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]) -
			dihedraly(x[0],x[1],x[8],x[6],x[7],x[5])-
			dihedraly(x[0],x[8],x[10],x[9],x[11],x[7])-
			dihedraly(x[0],x[2],x[10],x[12],x[11],x[4]) + global::pi*2.0;
			break;
			case 2 :
				global::pi*2.0 
					-dihedraly(global::sqrt8,x[0],x[2], x[4],2.0,2.0)
					-dihedraly(global::sqrt8,x[0],x[10],x[11],2.0,2.0)
					-dihedraly(global::sqrt8,x[0],x[10],x[2],2.0,2.0);
			break;
			}
			break;

		default : cout << "unexpected case in constraint" << endl;
		}
    }

iter::iter(int ineqSwitch) {
	numiter = 20; numargs = 6; nconstr=0;
	if (ineqSwitch==500) numargs=13;
	if (ineqSwitch==501) numargs=13;
	if ((ineqSwitch>=1009)&&(ineqSwitch<1082)) numargs = 9;
	switch(ineqSwitch)
		{
		case 800 : case 801 : case 802 : case 803 :
		case 804 : case 805 : case 806 : numargs=16; break;

		case 999 : numargs = 12; break;

		case -4 :
		case 1200 : numargs = 9; break;
		case -7 : numargs = 3; break;
		case -13: numargs =9; break;
		case -14: numargs =9; break;
		case 240: case 241 : case 242 : case 248:
				numargs= 9; break;
		}
	// temp: 
	if (ineqSwitch==5) numargs =20;
	xmin = new double[numargs];
	xmax = new double[numargs];
	x = new double[numargs];
	constraintfunc = nofunc;
	func = generic;
	int i;
	for (i=0;i<numargs;i++) { xmin[i]=x[i]=2.0; xmax[i]=2.51; }
	INEQ_NUMBER = ineqSwitch;
	switch (ineqSwitch)
		{
		case 248:
			xmin[3]=global::dodecrad; xmax[3]=3.2;
			xmax[0]=xmax[1]=xmax[2]=xmax[4]=xmax[5]=xmax[6]=xmax[7]
			=xmax[8]=global::dodecrad;
			nconstr=1;
			break;

		case 247:
			xmin[0]=2.51; xmax[0]=global::sqrt8;
			break;
		case 246:
			xmin[0]=2.51; xmax[0]=global::sqrt8;
			xmin[3]=2.51; xmax[0]=global::sqrt8;
			nconstr=1;
			break;
		case 245:
			xmin[3]=xmax[3]=global::sqrt8;
			xmin[4]=2.; xmax[4]=global::sqrt8;
			xmin[5]=2.51; xmax[5]=global::sqrt8;
			//nconstr=1;
			break;
		case 2442:
			xmin[3]=xmax[3]=3.2;
			xmin[4]=xmax[4]=3.2;
			xmin[5]=xmax[5]=2.;
			xmax[0]=2.2;
			xmin[2]=xmax[2]=2.51;
			break;
		case 2441:
			xmin[3]=xmax[3]=3.2;
			xmin[4]=xmax[4]=3.2;
			xmin[5]=xmax[5]=2.;
			xmin[0]=2.2;
			break;
		case 244:
			xmin[5]=2.51;
			xmin[4]=xmax[4]=2.;
			xmin[3]=xmax[3]=3.2;
			break;
		case 243:
			xmin[3]=xmax[3]=3.2;
			xmin[0]=2.51; xmax[0]=global::sqrt8;
			nconstr=1;
			break;
		case 240: case 241:
			xmin[0]=xmax[0]=2.;
			xmin[1]=xmax[1]=2.;
			xmin[3]=global::sqrt8; xmax[3]=3.2;
			xmin[7]=global::sqrt8; xmax[7]=3.2;
			nconstr=1;
			break;	
		case 242:
			xmin[0]=2.51; xmax[0]=global::sqrt8;
			xmin[1]=2.4; xmin[5]=2.4;
			xmin[3]=3.2; xmax[3]=3.3;
			xmax[7]=2.0;
			xmin[8]=xmax[8]=global::sqrt8;
			nconstr=2;
			break;
			
		case 238: case 239:
			xmin[4]=xmax[4]=global::sqrt8;
			xmin[5]=xmax[5]=global::sqrt8;
			xmin[3]=3.7;  xmax[3]=xmax[1]+xmax[2];
			nconstr=2;
			break;
		case 91: case 96:
			xmin[3]=2.51; xmax[3]=global::sqrt8;
			nconstr=2;
			break;
		case 92: case 97:
			xmin[3]=2.6; xmax[3]=global::sqrt8;
			xmin[0]=2.2;
			break;
		case 93: case 98:
			xmin[3]=2.7; xmax[3]=global::sqrt8;
			xmax[0]=2.2;
			break;
		case 94: case 95: case 99: case 100:
			xmin[3]=2.51; xmax[3]=global::sqrt8;
			nconstr=1;
			break;
		case 8101:
			xmin[0]=2.51; xmax[0]=global::sqrt8;
			nconstr=2;
			break;
		case 8102:
			xmin[4]=2.51; xmax[4]=global::sqrt8;
			xmin[5]=2.77; xmax[5]=global::sqrt8;
			break;
		case 8103:
			xmin[4]=2.51; xmax[4]=global::sqrt8;
			xmin[5]=2.51; xmax[5]=global::sqrt8;
			nconstr=1;
			break;
		case 8104:
			xmin[4]=2.51; xmax[4]=2.77;
			xmin[5]=2.51; xmax[5]=2.77;
			nconstr=1;
			break;
		case 90:
			xmin[3]=2.51; xmax[0]=global::sqrt8;
			nconstr=1;
			break;
		case 892:
			xmin[3]=global::sqrt8; xmax[3]=global::sqrt8;
			xmin[4]=2.51; xmax[4]=global::sqrt8;
			break;
		case 893:
			xmin[3]=2.; xmax[3]=2.51;
			xmin[4]=2.51; xmax[4]=global::sqrt8;
			xmin[5]=2.51; xmax[5]=global::sqrt8;
			break;
		case 894:
			xmin[3]=global::sqrt8; xmax[3]=global::sqrt8;
			xmin[4]=2.51; xmax[4]=global::sqrt8;
			xmin[5]=2.51; xmax[5]=global::sqrt8;
			break;
		case 895:
			xmin[0]=2.51; xmax[0]=global::sqrt8;
			xmin[3]=global::sqrt8; xmax[3]=global::sqrt8;
			break;
		case 896:
			xmin[1]=2.51; xmax[1]=global::sqrt8;
			xmin[3]=xmax[3]=2.51;
			break;
		case 897:
			xmin[1]=2.51; xmax[1]=global::sqrt8;
			xmin[4]=2.51; xmax[4]=global::sqrt8;
			break;
		case 898:
			xmin[1]=2.51; xmax[1]=global::sqrt8;
			xmin[4]=2.51; xmax[4]=global::sqrt8;
			xmin[3]=xmax[3]=2.51;
			break;
		case 899:
			xmin[4]=2.51; xmax[4]=global::sqrt8;
			xmin[5]=2.51; xmax[5]=global::sqrt8;
			break;
		case 891:
			xmin[4]=2.51; xmax[4]=3.2;
			xmin[5]=2.51; xmax[5]=3.2;
			break;
		case 89: 
			//xmin[0]=2.51; xmax[0]=global::sqrt8;
			xmin[3]=global::sqrt8; xmax[3]=global::sqrt8;
			xmin[4]=2.51; xmax[4]=global::sqrt8;
			xmin[5]=2.51; xmax[5]=global::sqrt8;
			break;
		case 830:
			xmin[4]=2.77; xmax[4]=global::sqrt8;
			xmin[5]=2.51; xmax[5]=global::sqrt8;
			break;
		case 831:
			xmin[4]=2.51; xmax[4]=global::sqrt8;
			xmin[5]=2.51; xmax[5]=global::sqrt8;
			nconstr=1;
			break;
		case 832:
			xmin[4]=2.51; xmax[4]=2.77;
			xmin[5]=2.51; xmax[5]=2.77;
			nconstr=1;
			break;
		case 88 :
			xmin[3]=xmax[3]=global::sqrt8;
			break;
		case -18 :
			xmin[0]=2.1;
			xmin[5]=xmax[5]=2.51;
			xmin[3]=3.2; xmax[3]=3.2;
			break;
		case -15 : case -16 : case -17 :
			xmin[3]=2.51; xmax[3]=global::sqrt8;
			nconstr=1;
			break;
		case -14 :
			xmax[5]=xmax[7]=xmax[8]=2.0;
			xmin[4]=global::sqrt8; xmax[4]=xmax[0]+xmax[2];
			xmin[3]=3.2;
			xmin[5]=xmax[5]=2;
			xmin[8]=xmax[8]=2.51;
			xmin[7]=xmax[7]=2.;
			xmin[3]=global::sqrt8; xmax[3]=xmax[7]+xmax[8];
			if (xmax[3]>xmax[1]+xmax[2]) xmax[3]=xmax[1]+xmax[2];
			nconstr=6;
			break;
		case -13 :
			xmin[3]=global::sqrt8; xmax[3]=3.2;
			xmax[7]=2.51;
			nconstr=6;
			break;
		case -11 :
			break;
		case -9 :
			xmin[0]=2.51; xmax[0]=global::sqrt8;
			xmin[5]=2.51; xmax[5]=3.49;
			xmax[3]=2.0;
			nconstr=1; constraintfunc=ConstraintPage1;
			break;
		case -7 :
			xmax[2]=2.417;
			xmin[2]=2.189;
			nconstr=1; constraintfunc=ConstraintPage1;
			break;
		case -5 :
			xmin[0]=2.51; xmax[0]=global::sqrt8;
			//xmin[3]=global::sqrt8;
			//xmax[3]=3.2;
			//nconstr=1; constraintfunc=ConstraintPage1;
			break;
		case -3 :
			xmin[3]=xmax[3]=global::sqrt8;
			xmin[4]=xmax[4]=2;
			xmin[5]=xmax[5]=2;
			break;
		case -2 :
			xmin[3]=xmax[3]=global::sqrt8; 
			break;
		case -1 :
			xmin[0]=2.51; xmax[0]=global::sqrt8;
			xmax[3]=xmax[4]=xmax[5]=2;
			nconstr=1; constraintfunc=ConstraintPage1;
			break;
		case -4 :
			xmin[4]=xmax[4]=2.0;
			xmin[5]=xmax[5]=2.0;
			xmin[7]=xmax[7]=2.0;
			xmin[8]=global::sqrt8; xmax[8]=3.2;
			xmin[3]=global::sqrt8; xmax[3]=3.2;
			nconstr=2; constraintfunc=ConstraintPage1;
			break;
			

		case 1200 :
			xmin[4]=xmax[4]=2.0;
			xmin[5]=xmax[5]=2.0;
			xmin[3]=global::sqrt8; xmax[3]=3.2;
			xmin[7]=xmax[7]=3.2;

			xmin[8]=xmax[8]=2.51;
			nconstr=3; constraintfunc=ConstraintPage1;
			break;

		case 1100 :
			xmin[4]=xmax[4]=2.0;
			xmin[5]=xmax[5]=2.0;
			xmin[3]=global::sqrt8; xmax[3]=xmax[4]+xmax[5];
			nconstr=1; constraintfunc=ConstraintPage1;
			break;

		case 1101 :
			xmin[4]=xmax[4]=2.0;
			xmin[5]=xmax[5]=2.51;
			xmin[3]=3.2; xmax[3]=xmax[4]+xmax[5];
			nconstr=1; constraintfunc=ConstraintPage1;
			break;

		case 1102 :
			xmin[4]=xmax[4]=2.0;
			xmin[5]=xmax[5]=global::sqrt8;
			xmin[3]=3.2; xmax[3]=xmax[4]+xmax[5];
			nconstr=1; constraintfunc=ConstraintPage1;
			break;

		case 1103 :
			xmin[4]=xmax[4]=2.51;
			xmin[5]=xmax[5]=2.51;
			xmin[3]=3.2; xmax[3]=xmax[4]+xmax[5];
			nconstr=1; constraintfunc=ConstraintPage1;
			break;

		case 1104 :
			xmin[4]=xmax[4]=2.51;
			xmin[5]=xmax[5]=global::sqrt8;
			xmin[3]=3.2; xmax[3]=xmax[4]+xmax[5];
			nconstr=1; constraintfunc=ConstraintPage1;
			break;

		case 1105 :
			xmin[4]=xmax[4]=global::sqrt8;
			xmin[5]=xmax[5]=global::sqrt8;
			xmin[3]=3.2; xmax[3]=xmax[4]+xmax[5];
			nconstr=1; constraintfunc=ConstraintPage1;
			break;

		case 1000 :
			xmin[3]=xmax[3]=2.0;
			xmin[4]=global::sqrt8; xmax[4]= 3.2;
			xmin[5]=xmax[5]=2.0;
			nconstr=1; constraintfunc=ConstraintPage1;
			break;

		case 999 :
			xmin[0]=2.51; xmax[0]=global::sqrt8;
			nconstr=1; constraintfunc=ConstraintPage1;
			break;

		case 1001 :
			xmin[3]=xmax[3]=2.0;
			xmin[5]=xmax[5]=2.51;
			xmin[4]=3.2; xmax[4]=2.0+2.51;
			nconstr=1; constraintfunc=ConstraintPage1;
			break;

		case 1002 :
			xmin[3]=xmax[3]=2.0;
			xmin[5]=xmax[5]=global::sqrt8;
			xmin[4]=3.2; xmax[4]=2.0+global::sqrt8;
			nconstr=1; constraintfunc=ConstraintPage1;
			break;

		case 1003 :
			xmin[3]=xmax[3]=2.51;
			xmin[5]=xmax[5]=2.51;
			xmin[4]=3.2; xmax[4]=xmax[3]+xmax[5];
			nconstr=1; constraintfunc=ConstraintPage1;
			break;

		case 1004 :
			xmin[3]=xmax[3]=2.51;
			xmin[5]=xmax[5]=global::sqrt8;
			xmin[4]=3.2; xmax[4]=xmax[3]+xmax[5];
			nconstr=1; constraintfunc=ConstraintPage1;
			break;

		case 1005 :
			xmin[3]=xmax[3]=global::sqrt8;
			xmin[5]=xmax[5]=global::sqrt8;
			xmin[4]=3.2; xmax[4]=xmax[3]+xmax[5];
			nconstr=1; constraintfunc=ConstraintPage1;
			break;

		case 1006 :
			xmin[3]=3.2; xmax[3]=5.02;
			xmin[4]=global::sqrt8; xmax[4]=3.2;
			xmin[5]=xmax[5]=2;
			nconstr=1; constraintfunc=ConstraintPage1;
			break;

		case 1007 :
			xmin[3]=3.2; xmax[3]=5.02;
			xmin[4]=global::sqrt8; xmax[4]=3.2;
			xmin[5]=xmax[5]=2.51;
			nconstr=1; constraintfunc=ConstraintPage1;
			break;

		case 1008 :
			xmin[3]=3.2; xmax[3]=5.02;
			xmin[4]=global::sqrt8; xmax[4]=3.2;
			xmin[5]=xmax[5]=global::sqrt8;
			nconstr=1; constraintfunc=ConstraintPage1;
			break;

		// 2 2 2 
		case 1010 : case 1010+18 : case 1010+36 : case 1010+54 :
			xmin[3]=xmax[3]=2;
			xmin[5]=xmax[5]=2;
			xmin[6]=xmax[6]=2;
			xmin[4]=global::sqrt8; xmax[4]=4;
			xmin[8]=global::sqrt8; xmax[8]=xmax[5]+xmax[3]+xmax[6];
			nconstr=2; constraintfunc=ConstraintPage1;
			break;

		// 2 2 2.51
		case 1011 : case 1011+18 : case 1011+36 : case 1011+54 :
			xmin[5]=xmax[5]=2;
			xmin[3]=xmax[3]=2;
			xmin[6]=xmax[6]=2.51;
			xmin[4]=global::sqrt8; xmax[4]=4;
			xmin[8]=global::sqrt8; xmax[8]=xmax[5]+xmax[3]+xmax[6];
			nconstr=2; constraintfunc=ConstraintPage1;
			if (ineqSwitch>1045) nconstr += 1;
			break;

		// 2 2 sqrt(8)
		case 1012 : case 1012+18 : case 1012+36 : case 1012+54 :
			xmin[5]=xmax[5]=2;
			xmin[3]=xmax[3]=2;
			xmin[6]=xmax[6]=global::sqrt8;
			xmin[4]=global::sqrt8; xmax[4]=4;
			xmin[8]=global::sqrt8; xmax[8]=xmax[5]+xmax[3]+xmax[6];
			nconstr=2; constraintfunc=ConstraintPage1;
			if (ineqSwitch>1045) nconstr += 1;
			break;

		// 2 2.51 2
		case 1013 : case 1013+18 : case 1013+36 : case 1013+54 :
			xmin[5]=xmax[5]=2;
			xmin[3]=xmax[3]=2.51;
			xmin[6]=xmax[6]=2;
			xmin[4]=3.2;  xmax[4]=xmax[5]+xmax[3];
			xmin[8]=global::sqrt8; xmax[8]=xmax[5]+xmax[3]+xmax[6];
			nconstr=2; constraintfunc=ConstraintPage1;
			if (ineqSwitch>1045) nconstr += 1;
			break;

		// 2 2.51 2.51
		case 1014 : case 1014+18 : case 1014+36 : case 1014+54 :
			xmin[5]=xmax[5]=2;
			xmin[3]=xmax[3]=2.51;
			xmin[6]=xmax[6]=2.51;
			xmin[4]=3.2;  xmax[4]=xmax[5]+xmax[3];
			xmin[8]=global::sqrt8; xmax[8]=xmax[5]+xmax[3]+xmax[6];
			nconstr=2; constraintfunc=ConstraintPage1;
			if (ineqSwitch>1045) nconstr += 1;
			break;

		// 2 2.51 sqrt(8)
		case 1015 : case 1015+18 : case 1015+36 : case 1015+54 :
			xmin[5]=xmax[5]=2;
			xmin[3]=xmax[3]=2.51;
			xmin[6]=xmax[6]=global::sqrt8;
			xmin[4]=3.2;  xmax[4]=xmax[5]+xmax[3];
			xmin[8]=global::sqrt8; xmax[8]=xmax[5]+xmax[3]+xmax[6];
			nconstr=2; constraintfunc=ConstraintPage1;
			if (ineqSwitch>1045) nconstr += 1;
			break;

		// 2 sqrt(8) 2
		case 1016 : case 1016+18 : case 1016+36 : case 1016+54 :
			xmin[5]=xmax[5]=2;
			xmin[3]=xmax[3]=global::sqrt8;
			xmin[6]=xmax[6]=2;
			xmin[4]=3.2;  xmax[4]=xmax[5]+xmax[3];
			xmin[8]=global::sqrt8; xmax[8]=xmax[5]+xmax[3]+xmax[6];
			nconstr=2; constraintfunc=ConstraintPage1;
			if (ineqSwitch>1045) nconstr += 1;
			break;

		// 2 sqrt(8) 2.51
		case 1017 : case 1017+18 : case 1017+36 : case 1017+54 :
			xmin[5]=xmax[5]=2;
			xmin[3]=xmax[3]=global::sqrt8;
			xmin[6]=xmax[6]=2.51;
			xmin[4]=3.2;  xmax[4]=xmax[5]+xmax[3];
			xmin[8]=global::sqrt8; xmax[8]=xmax[5]+xmax[3]+xmax[6];
			nconstr=2; constraintfunc=ConstraintPage1;
			if (ineqSwitch>1045) nconstr += 1;
			break;

		// 2 sqrt(8) sqrt(8)
		case 1018 : case 1018+18 : case 1018+36 : case 1018+54 :
			xmin[5]=xmax[5]=2;
			xmin[3]=xmax[3]=global::sqrt8;
			xmin[6]=xmax[6]=global::sqrt8;
			xmin[4]=3.2;  xmax[4]=xmax[5]+xmax[3];
			xmin[8]=global::sqrt8; xmax[8]=xmax[5]+xmax[3]+xmax[6];
			nconstr=2; constraintfunc=ConstraintPage1;
			if (ineqSwitch>1045) nconstr += 1;
			break;

		// 2.51 2 2.51 
		case 1019 : case 1019+18 : case 1019+36 : case 1019+54 :
			xmin[5]=xmax[5]=2.51;
			xmin[3]=xmax[3]=2;
			xmin[6]=xmax[6]=2.51;
			xmin[4]=3.2;  xmax[4]=xmax[5]+xmax[3];
			xmin[8]=global::sqrt8; xmax[8]=xmax[5]+xmax[3]+xmax[6];
			nconstr=2; constraintfunc=ConstraintPage1;
			if (ineqSwitch>1045) nconstr += 1;
			break;

		// 2.51 2 sqrt(8)
		case 1020 : case 1020+18 : case 1020+36 : case 1020+54 :
			xmin[5]=xmax[5]=2.51;
			xmin[3]=xmax[3]=2;
			xmin[6]=xmax[6]=global::sqrt8;
			xmin[4]=3.2;  xmax[4]=xmax[5]+xmax[3];
			xmin[8]=global::sqrt8; xmax[8]=xmax[5]+xmax[3]+xmax[6];
			nconstr=2; constraintfunc=ConstraintPage1;
			if (ineqSwitch>1045) nconstr += 1;
			break;

		// 2.51 2.51 2.51
		case 1021 : case 1021+18 : case 1021+36 : case 1021+54 :
			xmin[5]=xmax[5]=2.51;
			xmin[3]=xmax[3]=2.51;
			xmin[6]=xmax[6]=2.51;
			xmin[4]=3.2;  xmax[4]=xmax[5]+xmax[3];
			xmin[8]=global::sqrt8; xmax[8]=xmax[5]+xmax[3]+xmax[6];
			nconstr=2; constraintfunc=ConstraintPage1;
			if (ineqSwitch>1045) nconstr += 1;
			break;

		// 2.51 2.51 sqrt(8)
		case 1022 : case 1022+18 : case 1022+36 : case 1022+54 :
			xmin[5]=xmax[5]=2.51;
			xmin[3]=xmax[3]=2.51;
			xmin[6]=xmax[6]=global::sqrt8;
			xmin[4]=3.2;  xmax[4]=xmax[5]+xmax[3];
			xmin[8]=global::sqrt8; xmax[8]=xmax[5]+xmax[3]+xmax[6];
			nconstr=2; constraintfunc=ConstraintPage1;
			if (ineqSwitch>1045) nconstr += 1;
			break;

		// 2.51 sqrt(8) 2.51 
		case 1023 : case 1023+18 : case 1023+36 : case 1023+54 :
			xmin[5]=xmax[5]=2.51;
			xmin[3]=xmax[3]=global::sqrt8;
			xmin[6]=xmax[6]=2.51;
			xmin[4]=3.2;  xmax[4]=xmax[5]+xmax[3];
			xmin[8]=global::sqrt8; xmax[8]=xmax[5]+xmax[3]+xmax[6];
			nconstr=2; constraintfunc=ConstraintPage1;
			if (ineqSwitch>1045) nconstr += 1;
			break;

		// 2.51 sqrt(8) sqrt(8)
		case 1024 : case 1024+18 : case 1024+36 : case 1024+54 :
			xmin[5]=xmax[5]=2.51;
			xmin[3]=xmax[3]=global::sqrt8;
			xmin[6]=xmax[6]=global::sqrt8;
			xmin[4]=3.2;  xmax[4]=xmax[5]+xmax[3];
			xmin[8]=global::sqrt8; xmax[8]=xmax[5]+xmax[3]+xmax[6];
			nconstr=2; constraintfunc=ConstraintPage1;
			if (ineqSwitch>1045) nconstr += 1;
			break;

		// sqrt(8) 2 sqrt(8)
		case 1025 : case 1025+18 : case 1025+36 : case 1025+54 :
			xmin[5]=xmax[5]=global::sqrt8;
			xmin[3]=xmax[3]=2;
			xmin[6]=xmax[6]=global::sqrt8;
			xmin[4]=3.2;  xmax[4]=xmax[5]+xmax[3];
			xmin[8]=global::sqrt8; xmax[8]=xmax[5]+xmax[3]+xmax[6];
			nconstr=2; constraintfunc=ConstraintPage1;
			if (ineqSwitch>1045) nconstr += 1;
			break;

		// sqrt(8) 2.51 sqrt(8)
		case 1026 : case 1026+18 : case 1026+36 : case 1026+54 :
			xmin[5]=xmax[5]=global::sqrt8;
			xmin[3]=xmax[3]=2.51;
			xmin[6]=xmax[6]=global::sqrt8;
			xmin[4]=3.2;  xmax[4]=xmax[5]+xmax[3];
			xmin[8]=global::sqrt8; xmax[8]=xmax[5]+xmax[3]+xmax[6];
			nconstr=2; constraintfunc=ConstraintPage1;
			if (ineqSwitch>1045) nconstr += 1;
			break;

		// sqrt(8) sqrt(8) sqrt(8)
		case 1027 : case 1027+18 : case 1027+36 : case 1027+54 :
			xmin[5]=xmax[5]=global::sqrt8;
			xmin[3]=xmax[3]=global::sqrt8;
			xmin[6]=xmax[6]=global::sqrt8;
			xmin[4]=3.2;  xmax[4]=xmax[5]+xmax[3];
			xmin[8]=global::sqrt8; xmax[8]=xmax[5]+xmax[3]+xmax[6];
			nconstr=2; constraintfunc=ConstraintPage1;
			if (ineqSwitch>1045) nconstr += 1;
			break;


		case 900 : case 903 :
			xmin[3]=2.51; xmax[3]=global::sqrt8;
			nconstr=1; constraintfunc=ConstraintPage1;
			break;

		case 904 : 
			xmin[3]=2.6; xmax[3]=global::sqrt8;
			nconstr=1; constraintfunc=ConstraintPage1;
			break;

		case 901 : case 902 :
			xmin[3]=2.51; xmax[3]=global::sqrt8;
			nconstr=2; constraintfunc=ConstraintPage1;
			break;

		case 800 : 
			xmin[15]=global::sqrt8; xmax[15]=global::sqrt8;
			nconstr=1; constraintfunc=ConstraintPage1;
			break;

		case 801 : 
			xmin[15]=2.51; xmax[15]=global::sqrt8;
			nconstr=1; constraintfunc=ConstraintPage1;
			break;

		case 806 : 
			xmin[15]=2.6; xmax[15]=global::sqrt8;
			nconstr=1; constraintfunc=ConstraintPage1;
			break;

		case 802 : case 803 : case 805 :
			xmin[15]=2.51; xmax[15]=global::sqrt8;
			nconstr=2; constraintfunc=ConstraintPage1;
			break;

		case 804 : 
			xmin[15]=2.72; xmax[15]=global::sqrt8;
			nconstr=1; constraintfunc=ConstraintPage1;
			break;

		case 6 :
			xmin[3]=2.51; xmax[3]=global::sqrt8;
			break;
		case 7 :
			xmin[3]=2.51; xmax[3]=global::sqrt8;
			nconstr = 1; constraintfunc=ConstraintPage1;
			break;


		// START OF FLAT CASES.
		case 10 : case 11 : case 12 : case 13 : case 14 : case 15 :
			xmin[3]=2.51; xmax[3]=global::sqrt8;
			break;
		case 16 : case 49 :
			xmin[3]=2.51; xmax[3]=global::sqrt8;
			nconstr = 2; constraintfunc=ConstraintPage1;
			break;
		case 17 : case 18 :
		case 50 : case 51 :
			xmin[3]=2.51; xmax[3]=global::sqrt8;
			nconstr = 1; constraintfunc=ConstraintPage1;
			break;
		// END OF FLAT CASES.

		case 52 : 
			nconstr=1; constraintfunc=ConstraintPage1;
			break;
		case 53 : 
			nconstr=2; constraintfunc=ConstraintPage1;
			break;

		// UPRIGHT QUARTERS IN QUADS
		case 20 : case 21 : case 22 : case 23 : case 24 : case 25 :
			xmin[0]=2.51; xmax[0]=global::sqrt8; break;
		case 26 : case 261 :
			xmin[0]=2.51; xmax[0]=global::sqrt8;
			xmax[1]= 2.13; xmax[2]= 2.13;
			nconstr=1; constraintfunc=ConstraintPage1;
			break;
		case 27 : case 70 : case 71 : case 271 :
			xmin[0]=2.51; xmax[0]=global::sqrt8;
			xmax[1]= 2.13; xmax[2]= 2.13;
			nconstr=1; constraintfunc=ConstraintPage1;
			break;
		case 28 : case 69 : case 281 :
			xmin[0]=2.51; xmax[0]=global::sqrt8;
			xmax[1]= 2.13; xmax[2]= 2.13;
			nconstr=2; constraintfunc=ConstraintPage1;
			break;
		// END UPRIGHT QUARTERS IN QUADS

		// QUAD VC
		case 30 :
			xmin[3]=global::sqrt8;  xmax[3]=3.0;
			break;
		case 31 :
			xmin[3]=global::sqrt8;  xmax[3]=3.2;
			nconstr=1; constraintfunc=ConstraintPage1;
			break;
		case 32 :
			xmin[3]=xmax[3]=xmin[4]=xmax[4]=xmin[5]=xmax[5]=2;
			break;
		// END QUAD VC

		// QRTET edge inequalities.
		case 200 : 
			//xmax[0]=xmax[1]=xmax[2]=2.0;
			nconstr=1; constraintfunc=ConstraintPage1;
			break;
		case 201 : case 202 : case 203 : 
			nconstr=1; constraintfunc=ConstraintPage1;
			break;
		case 204 : case 205 : case 206 : case 207 :
		case 208 : case 209 : case 210 : 
			nconstr=2; constraintfunc=ConstraintPage1;
			break;
		case 2091 :
			nconstr=2; constraintfunc=ConstraintPage1;
			xmax[3]=xmax[4]=xmax[5]=2.2;
			break;
		case 211 : case 212 : case 2101 :
			nconstr=3; constraintfunc=ConstraintPage1;
			break;

		case 40 :
			xmin[3]=2.001; xmax[0]=2.2;
			break;
		case 41 : case 42 :
			xmin[0]=2.001; xmax[0]=2.2; xmin[3]=xmax[3]=2.0; break;
		case 43 :
			xmin[0]=2.001; xmax[0]=2.2; xmin[3]=xmax[3]=global::sqrt8; break;
		case 44 :
			xmin[0]=xmax[0]=2.0; xmin[3]=xmax[3]=2.0; 
			xmin[4]=2.001;
			nconstr = 1; constraintfunc=ConstraintPage1; break;
		case 45 :
			xmin[0]=xmax[0]=2.0; xmin[3]=xmax[3]=2.0; 
			xmin[5]=2.00001;
			nconstr = 1; constraintfunc=ConstraintPage1; break;
		case 46 :
			xmin[5]=xmax[5]=2.0;
			xmin[3]=2.0001;
			nconstr = 1; constraintfunc=ConstraintPage1; break;
		case 47 :
			xmin[5]=xmax[5]=2.0;
			break;//nconstr = 1; constraintfunc=ConstraintPage1; break;


		case 101: case 103 : case 11103:
			xmin[0]=2.75; xmax[0]=global::sqrt8; // modified.
			nconstr=2; constraintfunc=ConstraintPage1;
			break;

		case 102: case 104 : case 11104:
			xmin[0]=2.75; xmax[0]=global::sqrt8; // modified.
			nconstr=1; constraintfunc=ConstraintPage1;
			break;

		case 11105 :
			xmin[0]=2.51; xmax[0]=2.75;
			xmin[3]=2.51; xmax[3]=global::sqrt8;
			nconstr=1;
			break;
		case 11106:
			xmin[0]=2.51; xmax[0]=2.75;
			xmin[3]=2.51; xmax[3]=global::sqrt8;
			xmin[3]=2.77;
			nconstr=1;
			break;
		case 11107: case 11108:
			xmin[0]=2.51; xmax[0]=2.75;
			xmin[3]=2.51; xmax[3]=global::sqrt8;
			nconstr = 2;
			xmin[3]=2.51;
			break;
			

		case 300 : 
			xmin[0]=2.51; xmax[0]=global::sqrt8;
			xmax[0]=2.65; 
			//{ for (int q = 1;q<6;q++) xmax[q]=2.2; }
			nconstr=1; constraintfunc=ConstraintPage1;
			break;

		case 400 :
			xmin[3]=global::sqrt8; xmax[3]=4.5;
			nconstr=1; constraintfunc=ConstraintPage1;
			break;

		case 500 :
			xmin[6]=3.2; xmax[6]=3.2;
			xmin[0]=2.51; xmax[0]=global::sqrt8;
			nconstr=5; constraintfunc=ConstraintPage1;
			break;
		case 501 :
			xmin[6]=3.2; xmax[6]=3.2;
			xmin[0]=2.51; xmax[0]=global::sqrt8;
			nconstr=5; constraintfunc=ConstraintPage1;
			break;

		case 700 : 
			xmin[0]=2.3;
			nconstr=2; constraintfunc=ConstraintPage1;
			break;

		case 701 : 
			xmin[0]=2.3;
			nconstr=1; constraintfunc=ConstraintPage1;
			break;

		case 702 : 
			xmin[0]=2.3;
			xmin[5]=global::sqrt8; xmax[5]=3.02;
			nconstr=2; constraintfunc=ConstraintPage1;
			break;
		case 703 : 
			xmin[0]=2.3;
			xmin[5]=2.51; xmax[5]=3.02;
			xmin[4]=2.51; xmax[4]=global::sqrt8;
			nconstr=2; constraintfunc=ConstraintPage1;
			break;
		case 704 : 
			xmin[0]=2.3;
			xmin[5]=2.51; xmax[5]=global::sqrt8;
			nconstr=4; constraintfunc=ConstraintPage1;
			break;
		case 705 : 
			xmin[0]=2.3;
			xmin[5]=2.6; xmax[5]=global::sqrt8;
			nconstr=2; constraintfunc=ConstraintPage1;
			break;
		case 706 : case 707 : case 708 :
			xmin[0]=2.3;
			xmin[5]=2.51; xmax[5]=global::sqrt8;
			nconstr=3; constraintfunc=ConstraintPage1;
			break;

		case 503 : case 504 :
			xmin[3]=global::sqrt8; xmax[3]=3.9;
			xmin[4]=2; xmax[4]=2.417;   //2.189;
			xmin[5]=2; xmax[5]=3.2;
			nconstr=1; constraintfunc=ConstraintPage1;
			break;

		default : cout << "error " << ineqSwitch << ": not installed " << endl;
		}

	if (nconstr>0) constraintfunc=ConstraintPage1;
	}

void /*ineq.cc*/minimize2(int);

void page1()
	{
	/*
		minimize2(10);
		minimize2(11);
		minimize2(12);
		minimize2(13);
		minimize2(14);
		minimize2(15);
		minimize2(16);
		minimize2(17);
		minimize2(18);
		minimize2(47);

		minimize2(20);
		minimize2(21);
		minimize2(22);
		minimize2(23);
		minimize2(24);
		minimize2(25);
		minimize2(26);
		minimize2(27);
		minimize2(28);
		minimize2(30);
		minimize2(31);
		minimize2(32);
		minimize2(49);
		minimize2(50);
		minimize2(49);
		minimize2(50);
		minimize2(51);
		minimize2(261);
		minimize2(271);
		minimize2(281);

// contraction lemmas, appendix to IV.
		minimize2(600);
		minimize2(601);
		minimize2(602);
		minimize2(603);
		minimize2(604);
		minimize2(605);

		minimize2(700);
		minimize2(701);
		minimize2(702);
		minimize2(703);
		minimize2(704);
		minimize2(705);
		minimize2(706);
		minimize2(707);
		minimize2(800);
		minimize2(801);
		minimize2(802);
		minimize2(803);
		minimize2(804);
		minimize2(805);
		minimize2(806);


	// (K) checking that a pentagon with dih<1.32 squanders>5.66 :
		minimize2(900);
		minimize2(901);
		minimize2(902);
		minimize2(903);
		minimize2(904);

	//
		minimize2(604);
		minimize2(605);
		minimize2(606);
		minimize2(607);
		minimize2(608);

	// 1/5/98: III.Appendix. check that second derivative is positive.
	// Some of these give negative results.  So they won't be
	// used.
		minimize2(1000);
		minimize2(1001);
		minimize2(1002);
		minimize2(1003);
		minimize2(1004);
		minimize2(1005);
		minimize2(1006);
		minimize2(1007);
		minimize2(1008);

	// 1/6/98. Since some of 1/5/98 gives negative results we
	// have the following coupled systems for III.Appendix, second
	// derivative is positive (at critical points).

		minimize2(1010);
		minimize2(1011);
		minimize2(1012);
		minimize2(1013);
		minimize2(1014);
		minimize2(1015);
		minimize2(1016);
		minimize2(1017);
		minimize2(1018);
		minimize2(1019);

		minimize2(1020);
		minimize2(1021);
		minimize2(1022);
		minimize2(1023);
		minimize2(1024);
		minimize2(1025);
		minimize2(1026);
		minimize2(1027);
		minimize2(1028);
		minimize2(1029);

		minimize2(1030);
		minimize2(1031);
		minimize2(1032);
		minimize2(1033);
		minimize2(1034);
		minimize2(1035);
		minimize2(1036);
		minimize2(1037);
		minimize2(1038);
		minimize2(1039);

		minimize2(1040);
		minimize2(1041);
		minimize2(1042);
		minimize2(1043);
		minimize2(1044);
		minimize2(1045);
		minimize2(1046);
		minimize2(1047);
		minimize2(1048);
		minimize2(1049);

		minimize2(1050);
		minimize2(1051);
		minimize2(1052);
		minimize2(1053);
		minimize2(1054);
		minimize2(1055);
		minimize2(1056);
		minimize2(1057);
		minimize2(1058);
		minimize2(1059);

		minimize2(1060);
		minimize2(1061);
		minimize2(1062);
		minimize2(1063);
		minimize2(1064);
		minimize2(1065);
		minimize2(1066);
		minimize2(1067);
		minimize2(1068);
		minimize2(1069);

		minimize2(1070);
		minimize2(1071);
		minimize2(1072);
		minimize2(1073);
		minimize2(1074);
		minimize2(1075);
		minimize2(1076);
		minimize2(1077);
		minimize2(1078);
		minimize2(1079);

		minimize2(1080);
		minimize2(1081);


		minimize2(600);
		minimize2(601);
		minimize2(6011);
		minimize2(604);
		minimize2(605);

		minimize2(700);
		minimize2(701);
		minimize2(702);
		minimize2(703);
		minimize2(704);
		minimize2(705);
		minimize2(706);
		minimize2(707);
		minimize2(708);

		minimize2(91);
		minimize2(92);
		minimize2(93);
		minimize2(94);
		minimize2(95);
		minimize2(96);
		minimize2(97);
		minimize2(98);
		minimize2(99);
		minimize2(8102);
		minimize2(8103);
		minimize2(8104);
		minimize2(91);
		minimize2(92);
		minimize2(93);
		minimize2(94);
		minimize2(95);
		minimize2(96);
		minimize2(97);
		minimize2(98);
		minimize2(99);
		minimize2(100);

		minimize2(11103);
		minimize2(11104);
		minimize2(11105);
		minimize2(11106);
		minimize2(11107);
		minimize2(11108);
	*/
		minimize2(201);
		minimize2(202);
		minimize2(203);
		minimize2(204);
		minimize2(205);
		minimize2(206);
		minimize2(207);
		minimize2(208);
		minimize2(209);
		minimize2(210);
		minimize2(2101);
		minimize2(211);
		minimize2(212);
	/*
	*/
	}
