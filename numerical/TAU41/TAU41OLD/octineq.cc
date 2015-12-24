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
#include "octeqns.h"

// from CFSQP:
double constrained_min(double xmin[], double xmax[], double x[],
		int numargs,int nconstr,
        void (*func)(int numargs,int whichFn,double* x,double* ret,void*),
        void (*constraintfunc)
			(int numargs,int which,double* x,double* ret,void*));

static double iter_constrained_min(double xmin[], double xmax[], double x[],
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


double lowest = DBL_MAX;
int PRINTING = 0;
double minimize2(int j)
	{
	if (PRINTING) cout << "\nNew Case = " << j << "\n";
	iter X(j);
	int i;
	double t = iter_constrained_min(X.xmin,X.xmax,X.x,
			X.numiter,X.numargs,X.nconstr,
			X.func,
			X.constraintfunc);
	if (PRINTING)
		{
		cout.precision(18);
		cout << "numargs = " << X.numargs << endl;
		cout << t << endl << "{";
		for (i=0;i<X.numargs;i++) cout << X.x[i] << (i+1<X.numargs ? ",": "}") ;
		cout << endl;
		}
	if (t<lowest) lowest=t;
	return t;
	}

iter::~iter()
    { delete[] xmin;
      delete[] xmax;
        delete[] x; }
 

static void nofunc(int numargs,int whichFn,double* x,double* ret,void*)
    {
    cout << "nofunc should not be called" << endl << flush;
    }
 

// DETAILS OF IMPLEMENTATION.

double fA[2];
int type[4]={0,0,0,0};
double C[29];
double umin[13],umax[13];
int INEQ_NUMBER=0;

static void ConstraintPage1(int numargs,int whichFn,double* x,double* ret,void*)
    {
	static double eps = 1.0e-10;
    *ret = 0;
    switch (INEQ_NUMBER % 4) {
        case 0 : switch(whichFn) {
            case 1 : *ret = -eps+radf(x[0],x[1],x[5])-global::sqrt2; break;
            case 2 : *ret = -eps+radf(x[0],x[2],x[4])-global::sqrt2; break;
            }
            break;
        case 1 : switch(whichFn) {
            case 1 : *ret = -eps-radf(x[0],x[1],x[5])+global::sqrt2; break;
            case 2 : *ret = -eps+radf(x[0],x[2],x[4])-global::sqrt2; break;
            }
            break;
        case 2 : switch(whichFn) {
            case 1 : *ret = -eps+radf(x[0],x[1],x[5])-global::sqrt2; break;
            case 2 : *ret = -eps-radf(x[0],x[2],x[4])+global::sqrt2; break;
            }
            break;
        case 3 : switch(whichFn) {
            case 1 : *ret = -eps-radf(x[0],x[1],x[5])+global::sqrt2; break;
            case 2 : *ret = -eps-radf(x[0],x[2],x[4])+global::sqrt2; break;
            }
            break;
        default : cout << "unexpected case in constraint" << endl;
        }
    }

static void generic(int numargs,int whichFn,double* x, double* ret,void*)
    {
    switch (INEQ_NUMBER) {
        // start of first page of inequalities for Section 2, SPIV.
 
        case 0: case 1 : case 2 : case 3 :
			*ret = -fa(x[0],x[1],x[2],x[3],x[4],x[5],fA,type)
			 -C[0]-C[1]*x[0]-C[2]*x[1]-C[3]*x[2]-C[4]*x[4]-C[5]*x[5]
			 -C[6]*dpi(x[0],x[1],x[2],x[3],x[4],x[5])+C[28];
			 break;
        case 4: case 5 : case 6 : case 7 :
			*ret = -fb(x[0],x[1],x[2],x[3],x[4],x[5],fA,type)
			 -C[7+0]-C[7+1]*x[0]-C[7+2]*x[1]-C[7+3]*x[2]-C[7+4]*x[4]-C[7+5]*x[5]
			 -C[7+6]*dpi(x[0],x[1],x[2],x[3],x[4],x[5])+C[28];
			break;
        case 8: case 9 : case 10 : case 11 :
			*ret = -fc(x[0],x[1],x[2],x[3],x[4],x[5],fA,type) 
			 -C[14+0]-C[14+1]*x[0]-C[14+2]*x[1]-C[14+3]*x[2]
			 -C[14+4]*x[4]-C[14+5]*x[5]
			 -C[14+6]*dpi(x[0],x[1],x[2],x[3],x[4],x[5])+C[28];
			break;
        case 12: case 13 : case 14 : case 15 :
			*ret = -fd(x[0],x[1],x[2],x[3],x[4],x[5],fA,type) 
			 -C[21+0]-C[21+1]*x[0]-C[21+2]*x[1]-C[21+3]*x[2]
				-C[21+4]*x[4]-C[21+5]*x[5]
			 -C[21+6]*dpi(x[0],x[1],x[2],x[3],x[4],x[5])+C[28];
			break;
 
        default : cout << "generic default" << endl << flush;
            *ret=0;
        }
    }

iter::iter(int ineqSwitch) {
    numiter = 20; numargs = 6; nconstr=2;
    // temp:
    xmin = new double[numargs];
    xmax = new double[numargs];
    x = new double[numargs];
    constraintfunc = ConstraintPage1;
    func = generic;
    int i;
    for (i=0;i<numargs;i++) { xmin[i]=x[i]=2.0; xmax[i]=2.51; }
    INEQ_NUMBER = ineqSwitch;
    switch (ineqSwitch/4)
        {
        case 0 : 
			{
			int inv[6]={0,1,2,3,4,5};
			for (int i=0;i<6;i++)
				{
				xmin[i]=umin[inv[i]]; xmax[i]=umax[inv[i]]; 
				}
			}
            break;
        case 1 : 
			{
			int inv[6]={0,2,7,6,8,4};
			for (int i=0;i<6;i++)
				{
				xmin[i]=umin[inv[i]]; xmax[i]=umax[inv[i]]; 
				}
			}
            break;
        case 2 : 
			{
			int inv[6]={0,7,11,9,10,8};
			for (int i=0;i<6;i++)
				{
				xmin[i]=umin[inv[i]]; xmax[i]=umax[inv[i]]; 
				}
			}
            break;
        case 3 : 
			{
			int inv[6]={0,11,1,12,5,10};
			for (int i=0;i<6;i++)
				{
				xmin[i]=umin[inv[i]]; xmax[i]=umax[inv[i]]; 
				}
			}
            break;
		}
		//if (0 == (ineqSwitch  % 4)) nconstr = 2;
		//if (3 == (ineqSwitch  % 4)) nconstr = 2;
    }

int smallR(double y1,double y2,double y3)
	{
	return (radf(y1,y2,y3)<1.4142135623730950488);
	}

int bigR(double y1,double y2,double y3)
	{
	return (radf(y1,y2,y3)>1.4142135623730950488);
	}

int hasPoints(double ymin[13],double ymax[13],int tx[4])
	{
	/* return 1 if all four faces are compatible with types */
	if (tx[0]&&(smallR(ymax[0],ymax[1],ymax[5]))) return 0;
	if (tx[1]&&(smallR(ymax[0],ymax[2],ymax[4]))) return 0;
	if (tx[2]&&(smallR(ymax[0],ymax[7],ymax[8]))) return 0;
	if (tx[3]&&(smallR(ymax[0],ymax[10],ymax[11]))) return 0;
	if (tx[0]&&(bigR(ymin[0],ymin[1],ymin[5]))) return 0;
	if (tx[1]&&(bigR(ymin[0],ymin[2],ymin[4]))) return 0;
	if (tx[2]&&(bigR(ymin[0],ymin[7],ymin[8]))) return 0;
	if (tx[3]&&(bigR(ymin[0],ymin[10],ymin[11]))) return 0;
	return 1;
	}


void separateChecked(double ymin[13],double ymax[13],double CC[29],
	double f[2],int tx[4])
    {
	int i;
	if (!hasPoints(ymin,ymax,tx)) 
		{
		for (i=0;i<29;i++) CC[i]=0;
		CC[28]= -5000;
		return;
		}
	lowest = DBL_MAX;
	// PRINTING=1;
	for (i=0;i<29;i++) C[i]=CC[i];
	for (i=0;i<2;i++) fA[i]=f[i];
	for (i=0;i<4;i++) type[i]=tx[i];
	for (i=0;i<13;i++) { umin[i]=ymin[i]; umax[i]=ymax[i]; }
	double total = 0;
	separate(ymin,ymax,C,fA,type);
    double t;
	/* the offset 0,4,8,12 indicates whether we are at face a,b,c, or d.
	   the type+2*type' indicates which circumradius constraints hold */
	t = minimize2(type[0]+2*type[1]);
	total += t; C[0] += t;
	t = minimize2(4+type[1]+2*type[2]);
	total += t; C[7] += t;
	t = minimize2(8+type[2]+2*type[3]);
	total += t; C[14] += t;
	t = minimize2(12+type[3]+2*type[0]);
	total += t; C[21] += t;
	total = total/4.0;
	C[0] -= total; C[7] -= total; C[14] -= total; C[21] -= total;
	C[28] -= total;
	if (PRINTING)
		{
		cout << "ymin= {";
		for (i=0;i<13;i++) cout << ymin[i] << ","; cout << endl;
		cout << "ymax= {";
		for (i=0;i<13;i++) cout << ymax[i] << ","; cout << endl;
		 cout << "\n smallest case = " << lowest << endl << endl;
		 cout << "---------\n\n";
		}
	for (i=0;i<29;i++) CC[i]=C[i];
    }

static void test()
	{
	double f[2] = {-3.0508,9.494};
	int t[4]={0,0,0,0};
	double CC[29];
	double ymin[13]={2.51,2,2,2,2,2,2,2,2,2,2,2,2};
	double ymax[13]={2.51,2.2,2.2,2.2,2.2,2.2,2.2,2.2,2.2,2.2,2.2,2.2,2.2};
	separateChecked(ymin,ymax,CC,f,t);
	cout << " slack = " << C[28] << endl;
	}
