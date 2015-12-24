/*  A function that may be used in Sphere Packings IV.
	Second derivative wrt to x1 of VorVc,tauVc.
				*/
#include <math.h>
#include <iomanip.h>
#include "numerical.h"


double DaQuoin(double a,double b,double c) 
	{
	if ((b<=a)||(c<=b)) return 0.;
	double a2=a*a,b2=b*b,c2=c*c;
	double m = sqrt((b2-c2)/(a2-b2));
	return  (m*(-3.*a2 + b2 + 2.*c2) + 3.*(a2 - c2)*atan(m))/6.;
	}

double DaDaQuoin(double a,double b,double c)
	{
	if ((b<=a)||(c<=b)) return 0.;
	double a2=a*a,b2=b*b,c2=c*c;
	double m = sqrt((b2-c2)/(a2-b2));
	return
	  -(a*(m*(3.*a2 - 4.*b2 + c2) - 3.*(a2 - b2)*atan(m)))/
	   (3.*(a2 - b2));
	}


double DbQuoin(double a,double b,double c)
	{
	if ((b<=a)||(c<=b)) return 0.;
	double a2=a*a,b2=b*b,c2=c*c;
	double m = sqrt((b2-c2)/(a2-b2));
	return (a*(b2 - c2)*m)/(3.*b);
	}

double DaDbQuoin(double a,double b,double c)
	{
	if ((b<=a)||(c<=b)) return 0.;
	double a2=a*a,b2=b*b,c2=c*c;
	double m = sqrt((b2-c2)/(a2-b2));
	return  (b*(b2 - c2)*m)/(3.*(-a2 + b2));
	}

double DbDbQuoin(double a,double b,double c)
	{
	if ((b<=a)||(c<=b)) return 0.;
	double a2=a*a,b2=b*b,c2=c*c;
	double m = sqrt((b2-c2)/(a2-b2));
	return    (a*m*(2.*a2*b2 - b2*b2 + a2*c2 - 2.*b2*c2))/
   (3.*b2*(a2 - b2));
	}


void DDeta(double x1,double x3,double x5,
	double& eta,double& Deta,double& DDeta) // D[eta,{x1,2}];
	{
	double u = -x1*x1 + 2.*x1*x3 - x3*x3 + 2.*x1*x5 + 2.*x3*x5 - x5*x5;
	double eta2 = x1*x3*x5/u;
	double Deta2 = -x3*(-x1+x3-x5)*(x1+x3-x5)*x5/(u*u);
	double DDeta2 =   (2.*x3*x5*(x1*x1*x1 - 3.*x1*x3*x3 + 2.*x3*x3*x3 + 
		6.*x1*x3*x5 - 2.*x3*x3*x5 - 3.*x1*x5*x5 - 2.*x3*x5*x5 + 2.*x5*x5*x5))/
			(u*u*u);
	eta = sqrt(eta2);
	Deta = Deta2/(2.*eta);
	DDeta = (DDeta2 - 2.*Deta*Deta)/(2.*eta);
	}

double D2Quoin(double x1,double x3,double x5)
	{
	double e,De,DDe;
	DDeta(x1,x3,x5,e,De,DDe);
	double t0=1.255;
	double a = sqrt(x1)/2.;
	double Da = 1./(4.*sqrt(x1));
	double DDa= -1./(8.*sqrt(x1)*x1);
	
	return DaDaQuoin(a,e,t0)*Da*Da + DaQuoin(a,e,t0)*DDa
			+ DbDbQuoin(a,e,t0)*De*De+ DbQuoin(a,e,t0)*DDe
			+2.*DaDbQuoin(a,e,t0)*Da*De;
	}

double D1Quoin(double x1,double x3,double x5)
	{
	double e,De,DDe;
	DDeta(x1,x3,x5,e,De,DDe);
	double t0=1.255;
	double a = sqrt(x1)/2.;
	double Da = 1./(4.*sqrt(x1));
	return DaQuoin(a,e,t0)*Da + DbQuoin(a,e,t0)*De;
	}

double D2Quoin315(double x1,double x3,double x5)
	{
	double e,De,DDe;
	DDeta(x1,x3,x5,e,De,DDe);
	double t0=1.255;
	double a = sqrt(x3)/2.;
	return DbDbQuoin(a,e,t0)*De*De+ DbQuoin(a,e,t0)*DDe;
	}


double D1Quoin315(double x1,double x3,double x5)
	{
	double e,De,DDe;
	DDeta(x1,x3,x5,e,De,DDe);
	double t0=1.255;
	double a = sqrt(x3)/2.;
	return DbQuoin(a,e,t0)*De;
	}


void Ddih(double x1,double x2,double x3,double x4,double x5,double x6,
	double& sqrtdDd,double& sqrtd32DDd)
	{
	double delta =   -(x2*x3*x4) - x1*x3*x5 - x1*x2*x6 - x4*x5*x6 + 
   x3*(x1 + x2 - x3 + x4 + x5 - x6)*x6 + 
   x2*x5*(x1 - x2 + x3 + x4 - x5 + x6) + x1*x4*(-x1 + x2 + x3 - x4 + x5 + x6);
	double d1 =  -(x1*x4) + x2*x5 - x3*x5 - x2*x6 + x3*x6 + x4*(-x1 + x2 + x3 - x4 + x5 + x6);
	double d4 = -(x2*x3) - x1*x4 + x2*x5 + x3*x6 - x5*x6 + x1*(-x1 + x2 + x3 - x4 + x5 + x6);
	double d11= -2.*x4;
	double d14= -2.*x1 + x2 + x3 - 2.*x4 + x5 + x6;
	double d44= -2.*x1;
	double y1 = sqrt(x1);
	double u135= -x1*x1 + 2*x1*x3 - x3*x3 + 2*x1*x5 + 2*x3*x5 - x5*x5;
	double u126= -x1*x1 + 2*x1*x2 - x2*x2 + 2*x1*x6 + 2*x2*x6 - x6*x6;
	double sqrtd= (delta>0? sqrt(delta): 0);
	double ND = -2.*x1*delta*d14+ x1*d1*d4+delta*d4;
	double den= u126*u135*y1;
	// Sqrt[delta] D[dih,x1]:
	sqrtdDd = ND/den;

	// Now Sqrt[Delta]^3 *D[dih,{x1,2}]:	
	double num1 = -2.*(delta*d14+x1*d1*d14+x1*delta*(-2.))+
			(d1*d4+x1*d11*d4+x1*d1*d14)+(d1*d4+delta*d14);
	double term1 = num1*delta/den;
	double sqrtdDden= 2.*(-x1+x2+x6)*u135*y1*delta +
				2.*(-x1+x3+x5)*u126*y1*delta +
				u135*u126*(0.5/y1)*delta +
				u135*u126*y1*(0.5)*d1;
	sqrtd32DDd = term1 -ND*sqrtdDden/(den*den);
	}


double BB(double x) // x = y^2. // B[Sqrt[x]]
	{
	double y = sqrt(x);
	return 1.333333333333333 - 1.135440168063744*y 
				+ 0.06007524579312208*x*y;
	}

double dB(double x) // x = y^2   // D[B[Sqrt[x]],x].
	{
	double y = sqrt(x);
	 return 
		(-0.5677200840318724 + 0.0901128686896831*x)/y;
	}

double ddB(double x) // x = y^2 // D[B[Sqrt[x]],{x,2}].
	{
	double y = sqrt(x);
	 return 
		 (0.2838600420159362 + 0.04505643434484155*x)/(x*y);
	}

double sqrt32DDdih2(double x1,double x2,double x3,double x4,double x5,double x6)
	// Sqrt[delta]^3 D[dih2,{x1,2}];
	{
	double delta =   -(x2*x3*x4) - x1*x3*x5 - x1*x2*x6 - x4*x5*x6 + 
   x3*(x1 + x2 - x3 + x4 + x5 - x6)*x6 + 
   x2*x5*(x1 - x2 + x3 + x4 - x5 + x6) + x1*x4*(-x1 + x2 + x3 - x4 + x5 + x6);
	double d1 =  -(x1*x4) + x2*x5 - x3*x5 - x2*x6 + x3*x6 + x4*(-x1 + x2 + x3 - x4 + x5 + x6);
	double d13 =  x4 - x5 + x6;
	double d3=x1*x4 - x2*x4 - x1*x5 + x2*x5 - x3*x6 + (x1 + x2 - x3 + x4 + x5 - x6)*x6;
	double u126= -x1*x1 + 2*x1*x2 - x2*x2 + 2*x1*x6 + 2*x2*x6 - x6*x6;
	double Du126 = 2.*(-x1+x2+x6);
	double y2 = sqrt(x2);
	return -y2*d13*delta/u126 + y2*d3*(Du126*delta+u126*d1/2.)/(u126*u126);
	}

double sqrtDdih2(double x1,double x2,double x3,double x4,double x5,double x6)
	// Sqrt[delta] D[dih2,x1];
	{
	double d3=
	x1*x4 - x2*x4 - x1*x5 + x2*x5 - x3*x6 + (x1 + x2 - x3 + x4 + x5 - x6)*x6;
	double u126= -x1*x1 + 2*x1*x2 - x2*x2 + 2*x1*x6 + 2*x2*x6 - x6*x6;
	double y2 = sqrt(x2);
	return -y2*d3/u126;
	}

double sqrt32DDdih3(double x1,double x2,double x3,double x4,double x5,double x6)
	// Sqrt[delta]^3 D[dih3,{x1,2}];
	{
	double delta =   -(x2*x3*x4) - x1*x3*x5 - x1*x2*x6 - x4*x5*x6 + 
   x3*(x1 + x2 - x3 + x4 + x5 - x6)*x6 + 
   x2*x5*(x1 - x2 + x3 + x4 - x5 + x6) + x1*x4*(-x1 + x2 + x3 - x4 + x5 + x6);
	double d1 =  -(x1*x4) + x2*x5 - x3*x5 - x2*x6 + x3*x6 + x4*(-x1 + x2 + x3 - x4 + x5 + x6);
	double d12= x4 + x5 - x6;
	double d2= x1*x4 - x3*x4 - x2*x5 - x1*x6 + x3*x6 + x5*(x1 - x2 + x3 + x4 - x5 + x6);
	double u135= -x1*x1 + 2*x1*x3 - x3*x3 + 2*x1*x5 + 2*x3*x5 - x5*x5;
	double Du135 = 2.*(-x1+x3+x5);
	double y3 = sqrt(x3);
	return -y3*d12*delta/u135 + y3*d2*(Du135*delta+u135*d1/2.)/(u135*u135);
	}

double sqrtDdih3(double x1,double x2,double x3,double x4,double x5,double x6)
	// Sqrt[delta] D[dih3,x1];
	{
	double d2= 
	x1*x4 - x3*x4 - x2*x5 - x1*x6 + x3*x6 + x5*(x1 - x2 + x3 + x4 - x5 + x6);
	double u135= -x1*x1 + 2*x1*x3 - x3*x3 + 2*x1*x5 + 2*x3*x5 - x5*x5;
	double y3 = sqrt(x3);
	return -y3*d2/u135;
	}

double D2Vee0(double y1,double y2,double y3,double y4,double y5,double y6)
	{
	double x1=y1*y1,x2=y2*y2,x3=y3*y3,x4=y4*y4,x5=y5*y5,x6=y6*y6;
	double delta =   -(x2*x3*x4) - x1*x3*x5 - x1*x2*x6 - x4*x5*x6 + 
   x3*(x1 + x2 - x3 + x4 + x5 - x6)*x6 + 
   x2*x5*(x1 - x2 + x3 + x4 - x5 + x6) + x1*x4*(-x1 + x2 + x3 - x4 + x5 + x6);
	double sdelta= (delta>0? sqrt(delta) : delta);
	double term1 = sdelta*delta*ddB(x1)*dihedraly(y1,y2,y3,y4,y5,y6);
	double Ddih1,DDdih1;
	Ddih(x1,x2,x3,x4,x5,x6,Ddih1,DDdih1);
	double term2 = 2.*delta*dB(x1)*Ddih1;
	double term345= BB(x1)*DDdih1+BB(x2)*sqrt32DDdih2(x1,x2,x3,x4,x5,x6)+
			BB(x3)*sqrt32DDdih3(x1,x2,x3,x4,x5,x6);
	double quo = (D2Quoin(x1,x3,x5)+D2Quoin315(x1,x3,x5)
                        +D2Quoin(x1,x2,x6)+D2Quoin315(x1,x2,x6));
	double dt=0.7209029495174648;
	double termQUOIN =
		-4.*dt*delta*sdelta*(D2Quoin(x1,x3,x5)+D2Quoin315(x1,x3,x5)
						+D2Quoin(x1,x2,x6)+D2Quoin315(x1,x2,x6));
	return term1 + term2 + term345 + termQUOIN;
	}

double D2Vee1(double y1,double y2,double y3,double y4,double y5,double y6)
	{
	double zp = 0.1004445714270568;
	double x1=y1*y1,x2=y2*y2,x3=y3*y3,x4=y4*y4,x5=y5*y5,x6=y6*y6;
	double delta =   -(x2*x3*x4) - x1*x3*x5 - x1*x2*x6 - x4*x5*x6 + 
   x3*(x1 + x2 - x3 + x4 + x5 - x6)*x6 + 
   x2*x5*(x1 - x2 + x3 + x4 - x5 + x6) + x1*x4*(-x1 + x2 + x3 - x4 + x5 + x6);
	double sdelta= (delta>0? sqrt(delta) : delta);
	double term1 = sdelta*delta*ddB(x1)*dihedraly(y1,y2,y3,y4,y5,y6);
	double Ddih1,DDdih1;
	Ddih(x1,x2,x3,x4,x5,x6,Ddih1,DDdih1);
	double term2 = 2.*delta*dB(x1)*Ddih1;
	double term345= (BB(x1)-zp)*DDdih1+(BB(x2)-zp)*sqrt32DDdih2(x1,x2,x3,x4,x5,x6)+
			(BB(x3)-zp)*sqrt32DDdih3(x1,x2,x3,x4,x5,x6);
	double quo = (D2Quoin(x1,x3,x5)+D2Quoin315(x1,x3,x5)
                        +D2Quoin(x1,x2,x6)+D2Quoin315(x1,x2,x6));
	double termQUOIN =
		-4.*doct*delta*sdelta*(D2Quoin(x1,x3,x5)+D2Quoin315(x1,x3,x5)
						+D2Quoin(x1,x2,x6)+D2Quoin315(x1,x2,x6));
	return term1 + term2 + term345 + termQUOIN;
	}

double D1Vee0(double y1,double y2,double y3,double y4,double y5,double y6)
	{
	double x1=y1*y1,x2=y2*y2,x3=y3*y3,x4=y4*y4,x5=y5*y5,x6=y6*y6;
	double delta =   -(x2*x3*x4) - x1*x3*x5 - x1*x2*x6 - x4*x5*x6 + 
   x3*(x1 + x2 - x3 + x4 + x5 - x6)*x6 + 
   x2*x5*(x1 - x2 + x3 + x4 - x5 + x6) + x1*x4*(-x1 + x2 + x3 - x4 + x5 + x6);
	double sdelta= (delta>0? sqrt(delta) : delta);
	double term1 = sdelta*dB(x1)*dihedraly(y1,y2,y3,y4,y5,y6);
	double Ddih1,DDdih1;
	Ddih(x1,x2,x3,x4,x5,x6,Ddih1,DDdih1);
	double term345= BB(x1)*Ddih1+BB(x2)*sqrtDdih2(x1,x2,x3,x4,x5,x6)+
			BB(x3)*sqrtDdih3(x1,x2,x3,x4,x5,x6);
	double quo = (D1Quoin(x1,x3,x5)+D1Quoin315(x1,x3,x5)
                        +D1Quoin(x1,x2,x6)+D1Quoin315(x1,x2,x6));
	double dt=0.7209029495174648;
	double termQUOIN =
		-4.*dt*sdelta*(D1Quoin(x1,x3,x5)+D1Quoin315(x1,x3,x5)
						+D1Quoin(x1,x2,x6)+D1Quoin315(x1,x2,x6));
	return term1 + term345 + termQUOIN;
	}

double D1Vee1(double y1,double y2,double y3,double y4,double y5,double y6)
	{
	double zp = 0.1004445714270568;
	double x1=y1*y1,x2=y2*y2,x3=y3*y3,x4=y4*y4,x5=y5*y5,x6=y6*y6;
	double delta =   -(x2*x3*x4) - x1*x3*x5 - x1*x2*x6 - x4*x5*x6 + 
   x3*(x1 + x2 - x3 + x4 + x5 - x6)*x6 + 
   x2*x5*(x1 - x2 + x3 + x4 - x5 + x6) + x1*x4*(-x1 + x2 + x3 - x4 + x5 + x6);
	double sdelta= (delta>0? sqrt(delta) : delta);
	double term1 = sdelta*dB(x1)*dihedraly(y1,y2,y3,y4,y5,y6);
	double Ddih1,DDdih1;
	Ddih(x1,x2,x3,x4,x5,x6,Ddih1,DDdih1);
	double term345= (BB(x1)-zp)*Ddih1+(BB(x2)-zp)*sqrtDdih2(x1,x2,x3,x4,x5,x6)+
			(BB(x3)-zp)*sqrtDdih3(x1,x2,x3,x4,x5,x6);
	double quo = (D1Quoin(x1,x3,x5)+D1Quoin315(x1,x3,x5)
                        +D1Quoin(x1,x2,x6)+D1Quoin315(x1,x2,x6));
	double dt=0.7209029495174648;
	double termQUOIN =
		-4.*dt*sdelta*(D1Quoin(x1,x3,x5)+D1Quoin315(x1,x3,x5)
						+D1Quoin(x1,x2,x6)+D1Quoin315(x1,x2,x6));
	return term1 + term345 + termQUOIN;
	}
