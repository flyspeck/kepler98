#include <iomanip.h>
#include <stdlib.h>
#include "numerical.h"
#include "constants.h"
#include "morefn.h"
#include <math.h>

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

// cross diag in terms of edges . Equivalent to version in numerical.cc

static double quoin(double a,double b,double c)
  	{
	double u = sqrt((c*c - b*b)/(b*b - a*a));
    if ((a>b)||(b>c)) return 0;
    return -(a*a + a*c - 2*c*c)*(c - a)*atan(u)/6.0 + 
		(a*(b*b - a*a)*u)/6.0 - (2.0/3.0)*c*c*c*atan((u*(b - a))/(b + c));
	}

static double quoinH(double y1,double y2,double y6)
	{
	return quoin(y1/2.0,radf(y1,y2,y6),1.255);
	}

//////////
// u135M = - 4 doct u135 D[quoi135+quoin315,x5]
//
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


static double phi(double h,double t)
	{
	static const double doc = 0.72090294951746509284124;
	return 2.0*(2.0 - doc*h*t*(t+h))/3.0; 
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
	double delta4 = -(x2*x3) - x1*x4 + x2*x5 + x3*x6 - x5*x6 + x1*(-x1 + x2 + x3 - x4 + x5 + x6);
	double delta6 = -(x1*x2) + x1*x4 + x2*x5 - x4*x5 + x3*(x1 + x2 - x3 + x4 + x5 - x6) - x3*x6;
	double u135 = U(x1,x3,x5);
	return -(B(y1)-zp)*y1*delta6 + (B(y2)-zp)*y2*u135 - (B(y3)-zp)*y3*delta4;
	}

static double cortauRoundLong(double y1,double lon)
	{
	double h = y1/2.0;
	double t0=1.255;
	double cospsi = (y1*y1+t0*t0-lon*lon)/(2.0*y1*t0);
	double costheta= h/t0;
	double tildemax = global::zetapt - phi(t0,t0);
	return ((1.0-cospsi)*tildemax +
			(1.0-costheta)*(phi(t0,t0)-phi(h,t0)))*2.0*global::pi;
	}

static double cortauRound(double y1)
	{
	return cortauRoundLong(y1,1.6);
	}

static double crossTerm(double y1,double y2,double lambda)
	{
	double lon=3.2;
	double zp = global::zetapt;
	double t0=1.255;
	double alpha1=dihedraly(y1,t0,y2,lambda,lon,lambda);
	double alpha2=dihedraly(y2,t0,y1,lambda,lon,lambda);
	double sol=solid(y2,t0,y1,lambda,lon,lambda);
	double phi0=phi(t0,t0);
	double phi1=phi(y1/2.0,t0);
	double phi2=phi(y2/2.0,t0);
	double pi=global::pi;
	return 2.0*(sol*(zp-phi0)+
				alpha1*(1-y1/(2*t0))*(-phi1+phi0)+
				alpha2*(1-y2/(2*t0))*(-phi2+phi0))+
		(pi-2.0*alpha2)*cortauRoundLong(y2,lambda)/(2.0*pi)+
		(pi-2.0*alpha1)*cortauRoundLong(y1,lambda)/(2.0*pi);
	}

int INEQ_NUMBER=0;
static void generic(int numargs,int whichFn,double* x, double* ret,void*)
	{
	switch (INEQ_NUMBER) {
		// start of first page of inequalities for Section 2, SPIV.

		// IV contraction:
		// Section 1.
		case 705592875 :
			*ret = tau(x[0],x[1],x[2],x[3],x[4],x[5])-0.61*global::pt;
			break;
		case 705592876 : case 705592877 :
			*ret = tau_analytic(x[0],x[1],x[2],x[3],x[4],x[5]) - 0.61*global::pt;
			break;
		case 600340212 :
			*ret = delta(x[0]*x[0],x[1]*x[1],x[2]*x[2],x[3]*x[3],x[4]*x[4],x[5]*x[5]);
			break;
		// Appendix on contraction (IV).
		case 602317600: 
		  *ret = -Vee0(x[0],x[1],x[2],x[3],x[4],x[5]);
				break;
		case 602317601: 
		  *ret = -Vee1(x[0],x[1],x[2],x[3],x[4],x[5]);
				break;
		case 602317602: case 602317606 : 
		  *ret = -Vee0(x[0],x[1],x[2],x[3],x[4],x[5])-0.82*sqrt(421.0);
				break;
		case 602317603: case 602317607 : 
		  *ret = -Vee1(x[0],x[1],x[2],x[3],x[4],x[5])-0.82*sqrt(421.0);
				break;
		case 602317604: case 602317608 : 
		  *ret = -Vee0(x[0],x[1],x[2],x[3],x[4],x[5])-0.5*sqrt(421.0);
				break;
		case 602317605: case 602317609 : 
		  *ret = -Vee1(x[0],x[1],x[2],x[3],x[4],x[5])-0.5*sqrt(421.0);
				break;
		case 602317610: *ret= 421.0-delta(x[0]*x[0],x[1]*x[1],x[2]*x[2],x[3]*x[3],x[4]*x[4],x[5]*x[5]);
				break;
		case 602317611: *ret= 0.82 -u135M(x[0],x[2],x[4]);
				break;
		case 602317612: *ret= 0.5 -u135M(x[0],x[2],x[4]);
				break;
		case 308471379: *ret= cortauRound(x[0]);
				break;
		case 290054199: 
			{
			double t0 = 1.255;
			double y=x[0];
			double lambda=1.945;
			double psi = acos((t0*t0+y*y-lambda*lambda)/(2.0*t0*y));
			*ret=dihedraly(2.51,2.51,x[0],3.2,3.2,3.621)-beta(psi,x[0],2.51,3.2);
			}
			break;
		case 301657403 :
			*ret= 0.75*global::pi - dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]); 
			break;
		case 802813268 :
			*ret = crossTerm(x[0],x[1],1.945)-(14.8+1.21)*global::pt; 
			break;
		case 531903371 : case 881912007 :
			*ret= dihedraly(x[1],x[0],1.255,1.6,1.6,x[5])-
				dihedraly(x[1],x[0],x[2],x[4],x[3],x[5]);
			break;
		case 74877405 : case 502078137 :
			*ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-
				dihedraly(x[0],x[1],1.255,1.6,1.6,x[5])-
				dihedraly(x[0],x[2],1.255,1.6,1.6,x[4]);
			break;

		case 972512649 :
			*ret = cortau(x[0],x[1],x[2],x[3],x[4],x[5])-1.5488*global::pt; break;
		case 669958327 :
			*ret = cortau(x[0],x[1],x[2],x[3],x[4],x[5])-1.724*global::pt; break;
		case 7355789 :
			*ret = cortau(x[0],x[1],x[2],x[3],x[4],x[5])-1.787*global::pt; break;
		case 250333482 : case 336976176 : case 257461153 :
			*ret = cortau(x[0],x[1],x[2],x[3],x[4],x[5])-2.44*global::pt; break;

		case 702115356 :
			{
			double t = crossdiag(x);
			*ret =cortau(x[0],x[1],x[2],x[3],x[4],x[5])+
				cortau(x[1],x[0],x[6],t,x[8],x[5])-3.1*global::pt;
			}
			break;
		case 335682356 :
			{
			double t = crossdiag(x);
			*ret =cortau(x[0],x[1],x[2],x[3],x[4],x[5])+
				cortau(x[1],x[0],x[6],t,x[8],x[5])-(1.54+1.72)*global::pt;
			}
			break;
		case 858804156 :
			{
			double t = crossdiag(x);
			*ret =cortau(x[0],x[1],x[2],x[3],x[4],x[5])+
				cortau(x[1],x[0],x[6],t,x[8],x[5])-(1.54+1.787)*global::pt;
			}
			break;
		case 398156231 :
			{
			double t = crossdiag(x);
			*ret =cortau(x[0],x[1],x[2],x[3],x[4],x[5])+
				cortau(x[1],x[0],x[6],t,x[8],x[5])-(2.335)*global::pt;
			}
			break;
		case 631313278 : case 968194929 :
			{
			double t = crossdiag(x);
			*ret =cortau(x[0],x[1],x[2],x[3],x[4],x[5])+
				cortau(x[1],x[0],x[6],t,x[8],x[5])-(3.32)*global::pt;
			}
			break;
		case 452966150 :
			*ret = cortau(x[0],x[1],x[2],x[3],x[4],x[5])-0.3*global::pt; break;
			break;

		case 825350859 :
			*ret =tauVc(x[0],x[1],x[2],x[3],x[4],x[5])+
				tauVc(x[1],x[2],x[6],x[7],x[8],x[3])-(4.25)*global::pt;
			break;
		case 44912956 :
			*ret =tauVc(x[0],x[1],x[2],x[3],x[4],x[5])+
				tauVc(x[1],x[2],x[6],x[7],x[8],x[3])-(5.6)*global::pt;
			break;
		case 726107690 :
			*ret =tauVc(x[0],x[1],x[2],x[3],x[4],x[5])+
				tauVc(x[1],x[2],x[6],x[7],x[8],x[3])-(6.72)*global::pt;
			break;
		case 938593293 :
			*ret =tauVc(x[0],x[1],x[2],x[3],x[4],x[5])+
				tauVc(x[1],x[2],x[6],x[7],x[8],x[3])-(8.21)*global::pt;
			break;
		case 118374693 :
			*ret =tauVc(x[0],x[1],x[2],x[3],x[4],x[5])+
				tauVc(x[1],x[2],x[6],x[7],x[8],x[3])-(9.0)*global::pt;
			break;
		case 352630087 :
			*ret =tauVc(x[0],x[1],x[2],x[3],x[4],x[5])+
				tauVc(x[1],x[2],x[6],x[7],x[8],x[3])-(10.3)*global::pt;
			break;
		case 678021661 :
			*ret = tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-3.7*global::pt;
			break;
		case 11542431 :
			*ret = tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-6.58*global::pt;
			break;
		case 735258244 :
			{
			double psi = acos(x[2]/2.51);
			*ret =dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-beta(psi,x[0],x[2],x[4]); 
			}
			break;

		case 600951623 : *ret=0; break;
		case 333730624 : *ret=tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-1.07*global::pt;
			break;
		case 219059227 : *ret=tau_analytic(x[0],x[1],x[2],x[3],x[4],x[5])-2.518*global::pt;
			break;
		case 219059228 : *ret=tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-2.518*global::pt;
			break;
		case 219059229 : *ret=tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-2.518*global::pt;
			break;
		case 376614286 : *ret=tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-3.17*global::pt;
			break;
		case 688359402 : *ret=tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-4.25*global::pt;
			break;
		case 226629472 : *ret=tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-4.26*global::pt;
			break;
		case 757798626 : *ret=tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-5.64*global::pt;
			break;
		case 710070524 : *ret=tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-6.81*global::pt;
			break;
		case 939838329 : *ret=tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-8.54*global::pt;
			break;
		case 73379960 : *ret=-corsigma(x[0],x[1],x[2],x[3],x[4],x[5])-0.656*global::pt;
			break;
		case 630331617 : *ret=-corsigma(x[0],x[1],x[2],x[3],x[4],x[5])-0.919*global::pt;
			break;
		case 510614908 : *ret=-corsigma(x[0],x[1],x[2],x[3],x[4],x[5])-1.02*global::pt;
			break;
		case 979368646 : *ret=-corsigma(x[0],x[1],x[2],x[3],x[4],x[5])-1.579*global::pt;
			break;
		case 301227566 : *ret=-corsigma(x[0],x[1],x[2],x[3],x[4],x[5])-1.818*global::pt;
			break;
		case 82624302 : *ret=-corsigma(x[0],x[1],x[2],x[3],x[4],x[5])-2.20*global::pt;
			break;

		case 614416493 :
			{
			double t = crossdiag(x);
			*ret =-corsigma(x[0],x[1],x[2],x[3],x[4],x[5])+
				-corsigma(x[1],x[0],x[6],t,x[8],x[5])-(2.0*0.656)*global::pt;
			}
			break;
		case 792901760 :
			{
			double t = crossdiag(x);
			*ret =-corsigma(x[0],x[1],x[2],x[3],x[4],x[5])+
				-corsigma(x[1],x[0],x[6],t,x[8],x[5])-(0.656+0.919)*global::pt;
			}
			break;
		case 954734867 :
			{
			double t = crossdiag(x);
			*ret =-corsigma(x[0],x[1],x[2],x[3],x[4],x[5])+
				-corsigma(x[1],x[0],x[6],t,x[8],x[5])-(0.656+1.02)*global::pt;
			}
			break;
		case 861020178 :
			{
			double t = crossdiag(x);
			*ret =-corsigma(x[0],x[1],x[2],x[3],x[4],x[5])+
				-corsigma(x[1],x[0],x[6],t,x[8],x[5])-(0.74)*global::pt;
			}
			break;
		case 379125722 :
			*ret = 0.3*global::pt-corsigma(x[0],x[1],x[2],x[3],x[4],x[5]);
			break;

		case 713466234 :
			*ret =-vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
				-vorVc(x[1],x[2],x[6],x[7],x[8],x[3])-(1.36)*global::pt;
			break;

		case 582437770 :
			*ret =-vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
				-vorVc(x[1],x[2],x[6],x[7],x[8],x[3])-(3.34)*global::pt;
			break;

		case 92640954 :
			*ret = 0.163*global::pt - vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 162618659 :
			*ret = -1.0319*global::pt - vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 162618660 :
			*ret = -1.0319*global::pt - vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 162618661 :
			*ret = -1.0319*global::pt - vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 609725228 :
			*ret = -1.52*global::pt - vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 726130636 :
			*ret = -2.64*global::pt - vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 642934483 :
			*ret = -2.37*global::pt - vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 143431740 :
			*ret = -3.51*global::pt - vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 366126947 :
			*ret = tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-(1.189/2.0)*global::pt; break;
		case 697089973 :
			*ret = -vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 519972222 :
			*ret = 2.46/2.0 - dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]); break;

		case 906280948:
			*ret= -3.58 + 2.28501*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
				-vorVc(x[0],x[1],x[2],x[3],x[4],x[5])-vorVc(x[6],x[1],x[2],x[3],x[7],x[8]); 
			break;
		case 931179769:
			*ret= -2.715 + 1.67382*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
				-vorVc(x[0],x[1],x[2],x[3],x[4],x[5])-vorVc(x[6],x[1],x[2],x[3],x[7],x[8]); 
			break;
		case 3452357:
			*ret= -1.517+ 0.8285*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
				-vorVc(x[0],x[1],x[2],x[3],x[4],x[5])-vorVc(x[6],x[1],x[2],x[3],x[7],x[8]); 
			break;
		case 323155964:
			*ret= -0.858+ 0.390925*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
				-vorVc(x[0],x[1],x[2],x[3],x[4],x[5])-vorVc(x[6],x[1],x[2],x[3],x[7],x[8]); 
			break;
		case 690309687:
			*ret= -0.358+ 0.12012*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]) +0.009
				-vorVc(x[0],x[1],x[2],x[3],x[4],x[5])-vorVc(x[6],x[1],x[2],x[3],x[7],x[8]); 
			break;
		case 894530411:
			*ret= -0.186 + 0.0501*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]) +0.009
				-vorVc(x[0],x[1],x[2],x[3],x[4],x[5])-vorVc(x[6],x[1],x[2],x[3],x[7],x[8]); 
			break;
		case 14074184:
			*ret= -3.48 + 2.1747*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
				-1.189*global::pt+
				tauVc(x[0],x[1],x[2],x[3],x[4],x[5])+tauVc(x[6],x[1],x[2],x[3],x[7],x[8]); 
			break;
		case 87853018:
			*ret= -3.06+ 1.87427*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
				-1.189*global::pt+
				tauVc(x[0],x[1],x[2],x[3],x[4],x[5])+tauVc(x[6],x[1],x[2],x[3],x[7],x[8]); 
			break;
		case 717179286:
			*ret= -1.58+ 0.83046*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
				-1.189*global::pt+
				tauVc(x[0],x[1],x[2],x[3],x[4],x[5])+tauVc(x[6],x[1],x[2],x[3],x[7],x[8]); 
			break;
		case 249519771:
			*ret= -1.06 + 0.48263*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
				-1.189*global::pt+
				tauVc(x[0],x[1],x[2],x[3],x[4],x[5])+tauVc(x[6],x[1],x[2],x[3],x[7],x[8]); 
			break;
		case 706069213:
			*ret= -0.83+ 0.34833*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
				-1.189*global::pt+
				tauVc(x[0],x[1],x[2],x[3],x[4],x[5])+tauVc(x[6],x[1],x[2],x[3],x[7],x[8]); 
			break;
		case 451756659:
			*ret= -0.5+ 0.1694*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
				-1.189*global::pt+
				tauVc(x[0],x[1],x[2],x[3],x[4],x[5])+tauVc(x[6],x[1],x[2],x[3],x[7],x[8]); 
			break;
		case 609025895: 
			*ret= -0.29+ 0.0822*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]) + 0.0014
				-1.189*global::pt+
				tauVc(x[0],x[1],x[2],x[3],x[4],x[5])+tauVc(x[6],x[1],x[2],x[3],x[7],x[8]); 
			break;



		default : cout << INEQ_NUMBER << " generic default : function not installed " << endl << flush;
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
	
		case 705592875 : switch(whichFn) {
			case 1 : *ret= -global::sqrt2 + radf(x[0],x[1],x[5]); break;
			case 2 : *ret= -global::sqrt2 + radf(x[0],x[2],x[4]); break;
			}
			break;
		case 705592876 : 
			*ret= global::sqrt2-radf(x[0],x[1],x[5]); break;
		case 705592877 :
			*ret= global::sqrt2-radf(x[0],x[2],x[4]); break;
		

		case 602317600 : case 602317601: switch(whichFn) {
			case 1 : *ret= -delta(x[0]*x[0],x[1]*x[1],x[2]*x[2],x[3]*x[3],x[4]*x[4],x[5]*x[5]); break;
			case 2 : *ret = x[4]-x[5]; break;
			}
			break;
		case 602317602: case 602317603: case 602317604 : case 602317605 : 
		case 602317606: case 602317607: case 602317608 : case 602317609 : 
			*ret= -delta(x[0]*x[0],x[1]*x[1],x[2]*x[2],x[3]*x[3],x[4]*x[4],x[5]*x[5]);
			break;
		case 602317610 : case 602317611 : case 602317612 : 
			*ret= -1.255 + radf(x[0],x[2],x[4]); 
			break;
		case 531903371 : case 74877405 :
			*ret= -vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 881912007 : case 502078137 :
			*ret= tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-1.189*global::pt; break;

		case 825350859:
		case 44912956:
			{
			double t = crossdiag(x);
			switch(whichFn) {
			case 1 : *ret=-delta(x[5]*x[5],x[8]*x[8],x[3]*x[3],x[7]*x[7],x[4]*x[4],t*t); break;
			case 2 : *ret=x[3]-x[4]-x[5]; break;
			case 3 : *ret=x[3]-x[1]-x[2]; break;
			case 4 : *ret=global::sqrt8-t; break;
			}
			}
			break;

		case 726107690:
		case 938593293:
		case 118374693: 
		case 352630087:
			{
			double t = crossdiag(x);
			switch(whichFn) {
			case 1 : *ret=-delta(x[5]*x[5],x[8]*x[8],x[3]*x[3],x[7]*x[7],x[4]*x[4],t*t); break;
			case 2 : *ret=x[3]-x[4]-x[5]; break;
			case 3 : *ret=x[3]-x[1]-x[2]; break;
			case 4 : *ret=3.2-t; break;
			}
			}
			break;

		case 219059227 : case 162618659 :
			*ret = radf(x[3],x[4],x[5])-global::sqrt2; break;
		case 219059229 : case 162618661 :
			*ret = -radf(x[3],x[4],x[5])+global::sqrt2; break;

		case 713466234:
			{
			double t = crossdiag(x);
			switch(whichFn) {
			case 1 : *ret=-delta(x[5]*x[5],x[8]*x[8],x[3]*x[3],x[7]*x[7],x[4]*x[4],t*t); break;
			case 2 : *ret=x[3]-x[4]-x[5]; break;
			case 3 : *ret=x[3]-x[1]-x[2]; break;
			case 4 : *ret=global::sqrt8-t; break;
			}
			}
			break;
		case 582437770:
			{
			double t = crossdiag(x);
			switch(whichFn) {
			case 1 : *ret=-delta(x[5]*x[5],x[8]*x[8],x[3]*x[3],x[7]*x[7],x[4]*x[4],t*t); break;
			case 2 : *ret=x[3]-x[4]-x[5]; break;
			case 3 : *ret=x[3]-x[1]-x[2]; break;
			case 4 : *ret=3.2-t; break;
			}
			}
			break;

		case 906280948:
		case 931179769:
		case 3452357:
		case 323155964:
		case 690309687:
		case 894530411:
		case 14074184:
		case 87853018:
		case 717179286:
		case 249519771:
		case 706069213:
		case 451756659:
		case 609025895: 
			{
			double t = crossdiag(x);
			switch(whichFn) {
			case 1 : *ret=-delta(x[5]*x[5],x[8]*x[8],x[3]*x[3],x[7]*x[7],x[4]*x[4],t*t); break;
			case 2 : *ret=x[3]-x[4]-x[5]; break;
			case 3 : *ret=x[3]-x[1]-x[2]; break;
			case 4 : *ret=2.51-t; break;
			case 5 : *ret=dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-2.46; break;
			}
			}
			break;

		default : cout << INEQ_NUMBER << " constraints not installed" << endl;
		}
    }

iter::iter(int ineqSwitch) {
	numiter = 20; numargs = 6; nconstr=0;
	switch(ineqSwitch)
		{
		case 702115356:
		case 335682356:
		case 858804156:
		case 398156231:
		case 631313278:
		case 968194929: numargs=9;
		case 825350859:
		case 44912956:
		case 726107690:
		case 938593293:
		case 118374693: 
		case 352630087: numargs=9;
		case 614416493 :
		case 792901760 :
		case 954734867 :
		case 861020178 : numargs=9;
		case 713466234 :
		case 582437770 : numargs=9;

		case 906280948:
		case 931179769:
		case 3452357:
		case 323155964:
		case 690309687:
		case 894530411:
		case 14074184:
		case 87853018:
		case 717179286:
		case 249519771:
		case 706069213:
		case 451756659:
		case 609025895: numargs=9;
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
		case 705592875 :
		case 705592876 :
		case 705592877 :
			xmin[0]=2.51; xmax[0]=global::sqrt8;
			nconstr=1;
			if (ineqSwitch==705592875) nconstr=2;
			break;
		case 600340212 :
			xmin[3]=2.51; xmax[3]=3.2;
			xmin[4]=xmax[4]=2.0;
			xmin[5]=xmax[5]=2.51;
			break;
		case 602317600 : case 602317601 :
			xmin[3]=3.2; xmax[3]=3.5; //xmax[1]+xmax[2];
			xmin[4]=3.2; xmax[4]=3.4;   
			xmin[5]=2.51; xmax[5]=2.51;
			nconstr=1; 
			break;
		case 602317602 : case 602317603 :
			xmin[3]=global::sqrt8; xmax[3]=3.2;
			xmax[4]=2.189;
			nconstr=1; 
			break;
		case 602317604 : case 602317605 :
			xmin[3]=global::sqrt8; xmax[3]=3.2;
			xmin[4]=2.189; xmax[4]=global::sqrt8;
			nconstr=1; 
			break;
		case 602317606 : case 602317607 :
			xmin[3]=3.2; xmax[3]=xmax[1]+xmax[2];
			xmax[4]=xmax[5]=3.2;
			xmax[4]=2.189;
			nconstr=1; 
			break;
		case 602317608 : case 602317609 :
			xmin[3]=3.2; xmax[3]=xmax[1]+xmax[2];
			xmax[4]=xmax[5]=3.2;
			xmin[4]=2.189; 
			nconstr=1; 
			break;
		case 602317610 : case 602317611 : 
			xmin[3]=global::sqrt8; xmax[3]=xmax[1]+xmax[2];
			xmin[4]=2; xmax[4]=2.417;  
			xmin[5]=2; xmax[5]=3.2;
			nconstr=1; 
			break;
		case 602317612 : 
			xmin[3]=global::sqrt8; xmax[3]=xmax[1]+xmax[2];
			xmin[4]=2.189; xmax[4]=2.417;  
			xmin[5]=2; xmax[5]=3.2;
			nconstr=1; 
			break;
		case 308471379 :
			break;
		case 290054199 :
			break;
		case 301657403 :
			xmin[3]=xmax[3]=3.2;
			xmin[4]=xmin[5]=3.2;
			xmax[4]=xmax[5]=3.7;
			break;
		case 802813268 :
			break;
		case 531903371 : case 881912007 :
		case 74877405 : case 502078137 :
			xmin[3]=global::sqrt8; xmax[3]=3.2;
			nconstr=1;
			break;

		case 972512649 :
			xmin[3]=xmax[3]=3.2;
			xmin[4]=xmax[4]=2.0;
			xmin[5]=xmax[5]=2.0;
			break;

		case 669958327 :
			xmin[3]=xmax[3]=3.2;
			xmin[4]=xmax[4]=2.51;
			xmin[5]=xmax[5]=2.0;
			break;

		case 7355789 :
			xmin[3]=xmax[3]=3.2;
			xmin[4]=xmax[4]=global::sqrt8;
			xmin[5]=xmax[5]=2.0;
			break;

		case 250333482 :
			xmin[3]=xmax[3]=3.2;
			xmin[4]=xmax[4]=2.51;
			xmin[5]=xmax[5]=2.51;
			break;

		case 336976176 :
			xmin[3]=xmax[3]=3.2;
			xmin[4]=xmax[4]=global::sqrt8;
			xmin[5]=xmax[5]=2.51;
			break;

		case 257461153 :
			xmin[3]=xmax[3]=3.2;
			xmin[4]=xmax[4]=global::sqrt8;
			xmin[5]=xmax[5]=global::sqrt8;
			break;

		case 702115356 :
			xmin[7]=xmax[7]=3.2;
			xmin[4]=xmax[4]=2.0;
			xmin[5]=xmax[5]=2.0;
			xmin[8]=xmax[8]=2.0;
			xmin[3]=global::sqrt8; xmax[3]=3.2;
			break;
		case 335682356 :
			xmin[7]=xmax[7]=3.2;
			xmin[4]=xmax[4]=2.0;
			xmin[5]=xmax[5]=2.0;
			xmin[8]=xmax[8]=2.51;
			xmin[3]=global::sqrt8; xmax[3]=3.2;
			break;
		case 858804156 :
			xmin[7]=xmax[7]=3.2;
			xmin[4]=xmax[4]=2.0;
			xmin[5]=xmax[5]=2.0;
			xmin[8]=xmax[8]=global::sqrt8;
			xmin[3]=global::sqrt8; xmax[3]=3.2;
			break;
		case 398156231 :
			xmin[7]=xmax[7]=global::sqrt8;
			xmin[4]=xmax[4]=2.0;
			xmin[5]=xmax[5]=2.0;
			xmin[8]=xmax[8]=2.0;
			xmin[3]=global::sqrt8; xmax[3]=3.2;
			break;
		case 631313278 :
			xmin[7]=xmax[7]=global::sqrt8;
			xmin[4]=xmax[4]=2.0;
			xmin[5]=xmax[5]=2.0;
			xmin[8]=xmax[8]=2.51;
			xmin[3]=global::sqrt8; xmax[3]=3.2;
			break;
		case 968194929 :
			xmin[7]=xmax[7]=global::sqrt8;
			xmin[4]=xmax[4]=2.0;
			xmin[5]=xmax[5]=2.0;
			xmin[8]=xmax[8]=global::sqrt8;
			xmin[3]=global::sqrt8; xmax[3]=3.2;
			break;
		case 452966150 :
			xmin[3]=global::sqrt8; xmax[3]=3.2;
			xmin[4]=xmax[4]=2.0;
			xmin[5]=xmax[5]=2.0;
			break;

		case 825350859 :
			xmin[4]=xmax[4]=2.0;
			xmin[5]=xmax[5]=2.0;
			xmin[7]=xmax[7]=2.51;
			xmin[8]=xmax[8]=2.0;
			xmin[3]=global::sqrt8; xmax[3]=3.17; 
			nconstr=4;
			break;
		case 44912956 :
			xmin[4]=xmax[4]=2.0;
			xmin[5]=xmax[5]=2.0;
			xmin[7]=xmax[7]=global::sqrt8;
			xmin[8]=xmax[8]=2.0;
			xmin[3]=global::sqrt8; xmax[3]=3.37;
			nconstr=4;
			break;
		case 726107690 :
			xmin[4]=xmax[4]=2.0;
			xmin[5]=xmax[5]=2.0;
			xmin[7]=xmax[7]=2.51;
			xmin[8]=xmax[8]=2.51;
			xmin[3]=global::sqrt8; xmax[3]=xmax[4]+xmax[5];
			nconstr=4;
			break;
		case 938593293 :
			xmin[4]=xmax[4]=2.0;
			xmin[5]=xmax[5]=2.0;
			xmin[7]=xmax[7]=2.51;
			xmin[8]=xmax[8]=global::sqrt8;
			xmin[3]=global::sqrt8; xmax[3]=xmax[4]+xmax[5];
			nconstr=4;
			break;
		case 118374693 :
			xmin[4]=xmax[4]=2.0;
			xmin[5]=xmax[5]=2.0;
			xmin[7]=xmax[7]=global::sqrt8;
			xmin[8]=xmax[8]=global::sqrt8;
			xmin[3]=global::sqrt8; xmax[3]=xmax[4]+xmax[5];
			nconstr=4;
			break;
		case 352630087 :
			xmin[4]=xmax[4]=2.0;
			xmin[5]=xmax[5]=2.51;
			xmin[7]=xmax[7]=2.51;
			xmin[8]=xmax[8]=2.51;
			xmin[3]=3.2; xmax[3]=3.53;
			nconstr=4;
			break;
		case 678021661 :
			xmin[3]=xmax[3]=2.0;
			xmin[4]=xmax[4]=2.51;
			xmin[5]=xmax[5]=3.2;
			break;
		case 11542431 :
			xmin[3]=xmax[3]=2.51;
			xmin[4]=xmax[4]=2.51;
			xmin[5]=xmax[5]=3.2;
			break;
		case 735258244 :
			xmin[3]=xmax[3]=2.0;
			xmin[4]=xmax[4]=2.51;
			xmin[5]=xmax[5]=3.2;
			break;

		case 600951623 : 
			cout << "DONE IN part4sec2.cc ";
			break;

		case 333730624 :
			xmin[3]=xmax[3]=2.0;
			xmin[4]=xmax[4]=2.0;
			xmin[5]=xmax[5]=global::sqrt8;
			break;
		case 219059227 :
			xmin[3]=2.51; xmax[3]=2.77;
			xmin[4]=2.51; xmax[4]=2.77;
			nconstr=1;
			break;
		case 219059228 :
			xmin[3]=2.77; xmax[3]=global::sqrt8;
			xmin[4]=2.51; xmax[4]=global::sqrt8;
			break;
		case 219059229 :
			xmin[3]=2.51; xmax[3]=2.77;
			xmin[4]=2.51; xmax[4]=2.77;
			nconstr=1;
			break;

		case 376614286 :
			xmin[4]=2.51; xmax[4]=global::sqrt8;
			xmin[5]=global::sqrt8; xmax[5]=3.2;	
			break;
		case 688359402 :
			xmin[4]=global::sqrt8; xmax[4]=3.2;
			xmin[5]=global::sqrt8; xmax[5]=3.2;
			break;
		case 226629472 :
			xmin[3]=2.51; xmax[3]=global::sqrt8;
			xmin[4]=2.51; xmax[4]=global::sqrt8;
			xmin[5]=2.51; xmax[5]=global::sqrt8;
			break;
		case 757798626 :
			xmin[3]=2.51; xmax[3]=global::sqrt8;
			xmin[4]=2.51; xmax[4]=global::sqrt8;
			xmin[5]=global::sqrt8; xmax[5]=3.2;
			break;
		case 710070524 :
			xmin[3]=2.51; xmax[3]=global::sqrt8;
			xmin[4]=global::sqrt8; xmax[4]=3.2;
			xmin[5]=global::sqrt8; xmax[5]=3.2;
			break;
		case 939838329 :
			xmin[3]=global::sqrt8; xmax[3]=3.2;
			xmin[4]=global::sqrt8; xmax[4]=3.2;
			xmin[5]=global::sqrt8; xmax[5]=3.2;
			break;

		case 73379960 :
			xmin[3]=xmax[3]=3.2;
			xmin[4]=xmax[4]=2.0;
			xmin[5]=xmax[5]=2.0;
			break;
		case 630331617 :
			xmin[3]=xmax[3]=3.2;
			xmin[4]=xmax[4]=2.51;
			xmin[5]=xmax[5]=2.0;
			break;
		case 510614908 :
			xmin[3]=xmax[3]=3.2;
			xmin[4]=xmax[4]=global::sqrt8;
			xmin[5]=xmax[5]=2.0;
			break;
		case 979368646 :
			xmin[3]=xmax[3]=3.2;
			xmin[4]=xmax[4]=2.51;
			xmin[5]=xmax[5]=2.51;
			break;
		case 301227566 :
			xmin[3]=xmax[3]=3.2;
			xmin[4]=xmax[4]=2.51;
			xmin[5]=xmax[5]=global::sqrt8;
			break;
		case 82624302 :
			xmin[3]=xmax[3]=3.2;
			xmin[4]=xmax[4]=global::sqrt8;
			xmin[5]=xmax[5]=global::sqrt8;
			break;

		case 614416493 :
			xmin[7]=xmax[7]=3.2;
			xmin[4]=xmax[4]=2.0;
			xmin[5]=xmax[5]=2.0;
			xmin[8]=xmax[8]=2.0;
			xmin[3]=global::sqrt8; xmax[3]=3.2;
			break;
		case 792901760 :
			xmin[7]=xmax[7]=3.2;
			xmin[4]=xmax[4]=2.0;
			xmin[5]=xmax[5]=2.0;
			xmin[8]=xmax[8]=2.51;
			xmin[3]=global::sqrt8; xmax[3]=3.2;
			break;
		case 954734867 :
			xmin[7]=xmax[7]=3.2;
			xmin[4]=xmax[4]=2.0;
			xmin[5]=xmax[5]=2.0;
			xmin[8]=xmax[8]=global::sqrt8;
			xmin[3]=global::sqrt8; xmax[3]=3.2;
			break;
		case 861020178 :
			xmin[7]=xmax[7]=global::sqrt8;
			xmin[4]=xmax[4]=2.0;
			xmin[5]=xmax[5]=2.0;
			xmin[8]=xmax[8]=2.0;
			xmin[3]=global::sqrt8; xmax[3]=3.2;
			break;
		case 379125722 :
			xmin[3]=xmax[3]=global::sqrt8;
			xmin[4]=xmax[4]=2.0;
			xmin[5]=xmax[5]=2.0;
			break;
		case 713466234 :
			xmin[4]=xmax[4]=2.0;
			xmin[5]=xmax[5]=2.0;
			xmin[8]=xmax[8]=2.0;
			xmin[7]=2.51; xmax[7]=3.2;
			xmin[3]=global::sqrt8;
			xmax[3]=3.25;
			nconstr=4;
			break;
		
		case 582437770 :
			xmin[4]=xmax[4]=2.0;
			xmin[5]=xmax[5]=2.0;
			xmin[8]=xmax[8]=2.51;
			xmin[7]=xmax[8]=2.51;
			xmin[3]=global::sqrt8;
			xmax[3]=3.14;
			nconstr=4;
			break;

		case 92640954 :
			xmin[3]=2.0; xmax[3]=2.0;
			xmin[4]=2.0; xmax[4]=2.0;
			xmin[5]=global::sqrt8; xmax[5]=global::sqrt8;
			break;
		case 162618659 :
			xmin[3]=2.51; xmax[3]=2.77;
			xmin[4]=2.51; xmax[4]=2.77;
			nconstr=1;
			break;
		case 162618660 :
			xmin[3]=2.77; xmax[3]=global::sqrt8;
			xmin[4]=2.51; xmax[4]=global::sqrt8;
			break;
		case 162618661 :
			xmin[3]=2.51; xmax[3]=2.77;
			xmin[4]=2.51; xmax[4]=2.77;
			nconstr=1;
			break;

		case 609725228 :
			xmin[4]=2.51; xmax[4]=global::sqrt8;
			xmin[5]=global::sqrt8; xmax[5]=3.2;
			break;

		case 726130636 :
			xmin[4]=global::sqrt8; xmax[4]=3.2;
			xmin[5]=global::sqrt8; xmax[5]=3.2;
			break;
		case 642934483 :
			xmin[3]=2.51; xmax[3]=global::sqrt8;
			xmin[4]=2.51; xmax[4]=global::sqrt8;
			xmin[5]=2.51; xmax[5]=global::sqrt8;
			break;
		case 143431740 :
			xmin[3]=global::sqrt8; xmax[3]=3.2;
			xmin[4]=2.51; xmax[4]=3.2;
			xmin[5]=2.51; xmax[5]=3.2;
			break;

		case 366126947 : case 697089973 :
			xmin[0]=2.51; xmax[0]=global::sqrt8;
			xmin[5]=2.51; xmax[5]=2.74;
			break;
		case 519972222 :
			xmin[0]=2.51; xmax[0]=global::sqrt8;
			xmin[5]=2.74; xmax[5]=3.4;
			xmin[3]=xmax[3]=2.0;
			break;

		case 906280948:
		case 931179769:
		case 3452357:
		case 323155964:
		case 690309687:
		case 894530411:
		case 14074184:
		case 87853018:
		case 717179286:
		case 249519771:
		case 706069213:
		case 451756659:
		case 609025895:

			xmin[0]=2.51; xmax[0]=global::sqrt8;	
			xmin[3]=global::sqrt8; xmax[3]=3.4; //xmax[4]+xmax[5];
			nconstr=5;
			break;
			

		default : cout << "error " << ineqSwitch << ": bounds not installed " << endl;
		}
		if (nconstr>0) constraintfunc=ConstraintPage1;
	}

void /*ineq.cc*/minimize2(int);

void page1()
	{
	// test code goes here....
/*
		minimize2(735258244); // A1 
		minimize2(705592875); // A13
		minimize2(705592876); // A13
		minimize2(705592877); // A13
		minimize2(600340212); // (not found in SPIV)

		minimize2(602317600);
		minimize2(602317600);
		minimize2(602317601);
		minimize2(602317602);
		minimize2(602317603);
		minimize2(602317604);
		minimize2(602317605);
		minimize2(602317606);
		minimize2(602317607);
		minimize2(602317608);
		minimize2(602317609);
		minimize2(602317610);
		minimize2(602317611);
		minimize2(602317612);
		minimize2(308471379);
		minimize2(290054199);

		//minimize2(301657403); // doesn't check.
		minimize2(802813268);
		//minimize2(531903371); // doesn't check
		//minimize2(881912007); // doesn't check
		//minimize2(74877405);  // doesn't check
		//minimize2(502078137); // doesn't check

		minimize2(972512649);
		minimize2(669958327);
		minimize2(7355789);
		minimize2(250333482);
		minimize2(336976176);
		minimize2(257461153);

		minimize2(702115356);
		minimize2(335682356);
		minimize2(858804156);
		minimize2(398156231);
		minimize2(631313278);
		minimize2(968194929);
		minimize2(452966150);

		minimize2(825350859);
		minimize2(44912956);
		minimize2(726107690);
		minimize2(938593293);
		minimize2(118374693);

		minimize2(352630087);
		minimize2(678021661);
		minimize2(11542431);  // modified! original doesn't work 
		minimize2(735258244);

		minimize2(600951623);
		minimize2(333730624);
		minimize2(219059227);
		minimize2(219059228);
		minimize2(219059229);
		minimize2(376614286);
		minimize2(688359402); // constant adjusted
		minimize2(226629472);
		minimize2(757798626);// constant adjusted
		minimize2(710070524);// constant adjusted
		minimize2(939838329);// constant adjusted

		minimize2(73379960);
		minimize2(630331617);
		minimize2(510614908);
		minimize2(979368646);
		minimize2(301227566);
		minimize2(82624302);

		minimize2(614416493);
		minimize2(792901760);
		minimize2(954734867);
		minimize2(861020178);
		minimize2(379125722);

		minimize2(713466234);
		minimize2(582437770);
		minimize2(92640954);
		minimize2(162618659);
		minimize2(162618660);
		minimize2(162618661);
		minimize2(609725228);
		minimize2(726130636);
		minimize2(642934483);
		minimize2(143431740);

		minimize2(366126947);
		minimize2(697089973);
		minimize2(519972222);

		minimize2(906280948);
		minimize2(931179769);
		minimize2(3452357);
		minimize2(323155964);
		minimize2(690309687);
		minimize2(894530411);
		minimize2(14074184);
		minimize2(87853018);
		minimize2(717179286);
		minimize2(249519771);
		minimize2(706069213);
		minimize2(451756659);
		minimize2(609025895);

		minimize2(602317600);
		minimize2(602317601);
	*/



	// Section A23:
	minimize2(4591018);		//A23
	minimize2(193728878);		//A23
	minimize2(2724096);		//A23
	minimize2(213514168);		//A23
	minimize2(750768322);		//A23
	minimize2(371464244);		//A23
	minimize2(657011065);		//A23
	minimize2(953023504);		//A23
	minimize2(887276655);		//A23
	minimize2(246315515);		//A23
	minimize2(784421604);		//A23
	minimize2(258632246);		//A23
	minimize2(404164527);		//A23
	minimize2(163088471);		//A23
	}