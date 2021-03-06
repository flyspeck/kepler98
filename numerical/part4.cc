#include <iomanip.h>
#include "numerical.h"

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

static void func0(int numargs,int whichFn,double* x,double* ret,void*)
    {
    if (ret) *ret= -gamma(x[0],x[1],x[2],x[3],x[4],x[5]);
    else { cout << "ret expected" << endl << flush; }
    if (whichFn!=1) cout << "whichFn=" << whichFn  << endl << flush;
    }

static void nofunc(int numargs,int whichFn,double* x,double* ret,void*)
    {
    cout << "nofunc should not be called" << endl << flush;
    }


static void temp(int numargs,int whichFn,double* x,double* ret,void*)
    {
	*ret = tau(x[0],x[1],x[2],x[3],x[4],x[5]) 
		- (x[0]+x[1]-4)*0.08456 - (x[4]-2)*0.15973 - (x[5]-2)*0.129;
    }

iter::~iter() 
	{ delete[] xmin; 
	  delete[] xmax; 
		delete[] x; }

iter::iter(int ineqSwitch) {
	numiter = 20; numargs = 6; nconstr=0;
	xmin = new double[numargs];
	xmax = new double[numargs];
	x = new double[numargs];
	int i;
	for (i=0;i<numargs;i++) { xmin[i]=x[i]=2.0; xmax[i]=2.51; }
	switch (ineqSwitch) {
		case 0 : func = func0; constraintfunc=nofunc; break;
		case 1 : func = temp; constraintfunc=nofunc; 
				xmax[0]=2.12; xmin[5]=2.25; xmax[2]=2.1;
				break;
		default : func = func0; constraintfunc=nofunc; break;
		}
	}


