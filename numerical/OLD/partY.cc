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

double vol_analytic(double y1,double y2,double y3,double y4,
	double y5,double y6)
	// dodecrad truncated volume of VOronoi
	{
	double dt = 0.7209029495174648;
	return
		(-1./dt)*(vor_analytic(y1,y2,y3,y4,y5,y6)/4.
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

static double octatau(double x1,double x2,double x3,double x4,double x5,
		double x6)
	{
	return solid(x1,x2,x3,x4,x5,x6)*global::zeta*global::pt - 
		0.5*(vor_analytic(x1,x2,x3,x4,x5,x6)+vor_analytic(x1,x6,x5,x4,x3,x2));
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
		case 1 :
			*ret = vol_analytic(x[0],x[1],x[2],x[3],x[4],x[5])
				+0.06229*(x[0]+x[1]+x[2]-6)
				-0.235702;
				break;
		case 2 :
			*ret = SeanVol(x[0],x[1],x[2],x[3],x[4],x[5])
				+0.243*(x[3]-global::dodecrad)
				+0.164*(x[0]-2.)
				+0.164*(x[1]+x[2]-4.)
				-0.198*(x[4]+x[5]-4.)
				-0.295;
		case 3 :
			*ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
				-0.503*(x[0]-2.)
				+0.347*(x[1]+x[2]-4)
				+0.347*(x[4]+x[5]-4)
				-1.145*(x[3]-global::dodecrad)
				-1.62655;
			break;
		case 4 :
			*ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
				-2.339
				-0.59*(x[0]-2)
				+0.65*(x[1]+x[2]+x[4]+x[5]-8.);
				break;
		case 5 :
			*ret = -x[3];
				break;


		// start of first page of inequalities for Section 2, SPIV.
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
		case 2 : case 3: case 4:
			*ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-2.5;
			break;
		break;
		case 5 :
			switch(whichFn)
			{
			case 1: *ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.715;
			break;
			case 2: *ret = x[0]+x[1]+x[2]-6.28; 
			break;
			case 3: *ret = x[4]+x[5]-4.94;
			break;
			}
			break;

		default : cout << "unexpected case in constraint" << endl;
		}
    }

iter::iter(int ineqSwitch) {
	numiter = 20; numargs = 6; nconstr=0;
	switch(ineqSwitch)
		{
			
		}
	// temp: 
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
		case 1 :
			xmax[0]=xmax[1]=xmax[2]=xmax[3]=xmax[4]=xmax[5]=
			global::dodecrad;
			break;
		case 2 :
			xmax[0]=xmax[1]=xmax[2]=xmax[4]=xmax[5]=
			global::dodecrad;
			xmin[3]=global::dodecrad; xmax[3]=3.5;
			nconstr=1;
			break;
		case 4 : 
			xmax[0]=xmax[1]=xmax[2]=xmax[4]=xmax[5]=
			global::dodecrad;
			xmin[3]=xmax[3]=3.2;
			nconstr=1;
			break;
		case 5 :
			xmin[3]=global::sqrt8; xmax[3]=3.5;
			xmin[5]=2.51; xmax[5]=global::sqrt8;
			nconstr=3;
			break;
			
			

		default : cout << "error " << ineqSwitch << ": not installed " << endl;
		}

	if (nconstr>0) constraintfunc=ConstraintPage1;
	}

void /*ineq.cc*/minimize2(int);

void page1()
	{
	minimize2(5);
	}
