#include <iomanip.h>
#include <iostream.h>
#include <fstream.h>
#include <math.h>
#include <stdlib.h>
#include "numerical.h"
#include "gradient.h"
#include "morefn.h"
#include "constants.h"
#include "float.h"

// from CFSQP:
double constrained_min(double xmin[], double xmax[], double x[],
		int numargs,int nconstr,
        void (*func)(int numargs,int whichFn,double* x,double* ret,void*),
        void (*constraintfunc)
			(int numargs,int which,double* x,double* ret,void*));

double iter_constrained_min(double xmin[], double xmax[], double x[],
		int numiter,int numargs,int nconstr,
        void (*func)(int numargs,int whichFn,double* x,double* ret,void*),
        void (*constraintfunc)
			(int numargs,int which,double* x,double* ret,void*))
	{
	double currentMin=DBL_MAX;
	int i;
	double* y;
	y = new double[numargs];
	for (i=0;i<numargs;i++) x[i]=xmin[i];
	for (i=0;i<numiter;i++) 
			{
			double t = constrained_min(xmin,xmax,y,numargs,nconstr,
				func,constraintfunc);
			t = constrained_min(xmin,xmax,y,numargs,nconstr,
				func,constraintfunc);
			if (t<currentMin) {
				currentMin=t;
				for (int j=0;j<numargs;j++) { x[j]=y[j]; }
				}
			}
	delete [] y;
	return currentMin;
	}

class iter 
{
public:
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


double vorvol(double y[6])
	{
	return (vor_analytic(y[0],y[1],y[2],y[3],y[4],y[5])
		-(4.0/3.0)*solid(y[0],y[1],y[2],y[3],y[4],y[5]))/(-4.0*doct);
	}

double vorIP(double y[6]) { return
	vorVc(y[0],y[1],y[2],y[3],y[4],y[5],1.385);
	}

double tauIP(double y[6]) { return
	-vorIP(y)+global::zetapt*solid(y[0],y[1],y[2],y[3],y[4],y[5]);
	}

double mustbepositive(double x)
	{
	return (x>0.0 ? 0.0 : -5000.0*x*x);
	}


double betahack(double y0,double y1, double y6)
	{
	double p = y0/2.77;
	double t = cos(triangle(y0,y1,y6));
	return acos ( sqrt((p*p-t*t)/(1.0-t*t)));
	}
double beta0(double y0,double y1,double y5)
	{
	if (delta(y0*y0,y1*y1,1.255*1.255,1.6*1.6,1.6*1.6,y5*y5)>0.002)
	return dihedraly(y0,y1,1.255,1.6,1.6,y5);
	return 0;
	}
double rad(double y[6])
	{
	return sqrt(circum2(y[0]*y[0],y[1]*y[1],y[2]*y[2],y[3]*y[3],
			y[4]*y[4],y[5]*y[5]));
	}

double fnTemp1(double y[4])   // put your favorite inequality here.
	{
	double u = 0.0,v=0.0;
	double t = 10.0;
	t=delta(y[0]*y[0],y[1]*y[1],y[2]*y[2],y[3]*y[3],y[4]*y[4],y[5]*y[5]);
	double f126= radf(y[0],y[1],y[5])-global::sqrt2;
	double y1=y[0],y2=y[1],y3=y[2],y4=y[3],y5=y[4],y6=y[5];
	if ((t > 0.05) && (t > 0.05) )
		{
		u = -tau_analytic(y[0],y[1],y[2],y[3],y[4],y[5])
			+0.0658 + 0.1*(y1-2)+0.061*(y2+y3-4) + 0.166*(y5+y6-4);
		}
	return u 
		+mustbepositive(t)
		+mustbepositive(-global::sqrt2+radf(y[3],y[4],y[5]))
		//+mustbepositive(-global::sqrt2+radf(y[1],y[2],y[3]))
		;
	}


void PARTIV(double y[])
	{
	int numargs=6;
	int numiters= 20; // was 20;
	double ymin[15], ymax[15];
	int i;
	for (i=0;i<15;i++) { ymin[i]=2.0; ymax[i]=2.51; }
	// ADJUST HERE: Temp
	ymin[3]=2.51; ymax[3]=global::sqrt8;
	double t = unconstrained_iterated_max(numiters,ymin,
			ymax,y,numargs,*fnTemp1);
	cout << t << endl << "{";
	for (i=0;i<numargs;i++) cout << y[i] << (i+1<numargs ? ",": "}");
	cout << endl;
	}

double lowest = DBL_MAX;
double minimize2(int j)
	{
	cout << "\nNew Case = " << j << "\n";
	iter X(j);
	int i;
	double t = iter_constrained_min(X.xmin,X.xmax,X.x,
			X.numiter,X.numargs,X.nconstr,
			X.func,
			X.constraintfunc);
	cout.precision(18);
	cout << "numargs = " << X.numargs << endl;
	cout << t << endl << "{";
	for (i=0;i<X.numargs;i++) cout << X.x[i] << (i+1<X.numargs ? ",": "}") ;
	cout << endl;
	if (t<lowest) lowest=t;
	if (X.numargs==9) { cout << "cd = " << crossdiag(X.x) << endl; }
	return t;
	}

void /*part4sec2.cc*/page1();
double crossdiagNew(double y1,double y2,double y3,double y4,double
        y5,double y6,double y7,double y8,double y9);

int main() 
	{
	/* test:
	double y[9]={2.0,2.05,2.07,2.09, 2.10,2.11,2.12,2.14,2.18};
	cout << crossdiagNew(y[0],y[1],y[2],y[3],y[4],y[5],y[6],y[7],y[8]) << endl;
	cout << crossdiag(y) << endl;
	*/
	
	page1();
	cout << "\n smallest case = " << lowest << endl << endl;
	cout << "---------\n\n";
	return 1;
	}
