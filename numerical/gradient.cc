#include <iomanip.h>
#include <iostream.h>
#include <fstream.h>
#include <math.h>
#include <stdlib.h>
#include "numerical.h"
#include "gradient.h"

//#define VERBOSE

static const int MAXARG = 20;

double gradient(int ith,
        double x[],int numargs,double (*func)(double x[])  )
        {
        double eps = 1.0e-2;  //BUGBUG 1.0e-3.  using 1.0e-2 for 2nd partials.
	if (numargs>MAXARG) { cout << "array overflow"; return 0.0; }
        double y[MAXARG]; 
	for (int i=0;i<numargs;i++) y[i]=x[i]; // make a copy:
        y[ith] += eps;
        return ((*func)(y) - (*func)(x))/eps;
        }

double second_partial(int i,int j,int numargs,
            double (*f)(double []),double x[])
        {
        double eps = 1.0e-4;
        double y[MAXARG];
	if (numargs>MAXARG) { cout << "array overflow"; return 0.0; }
	for (int k=0;k<numargs;k++) y[k]=x[k];
        double f1 = (*f)(x);
        y[i-1] += eps; double f2 = (*f)(y);
        y[j-1] += eps; double f4 = (*f)(y);
        y[i-1] -= eps; double f3 = (*f)(y);
        return (f4+f1-f2-f3)/(eps*eps);
        }


double unconstrained_max(double xmin[],double xmax[],
		double x[],int numargs,double (*func)(double x[]) )
	{
        double eps = 5.0e-2;
        int i,j,icount=0;

	for (i=0;i<numargs;i++) x[i] = xmin[i] + (xmax[i]-xmin[i])*myrand(); 
	if (numargs>MAXARG) { cout << "array overflow"; return 0.0; }
        double gra[MAXARG], z[MAXARG];
        double gx = (*func)(x), gz;
	int imax = 10000;
        for (i=0;((icount<imax)||(i<500));i++)
                {
                for (j=0;j<numargs;j++)
                        {
                        gra[j]=eps*gradient(j,x,numargs,*func);
                        if (((x[j]>=xmax[j])&&(gra[j]>0.0))
                            ||((x[j]<=xmin[j])&&(gra[j]<0.0))) gra[j]=0.0;
                        z[j]= x[j]+ gra[j];
                        if (z[j] > xmax[j]) z[j]=xmax[j];
                        if (z[j] < xmin[j]) z[j]=xmin[j];
                        }
                gz = (*func)(z);
                if (gz > gx)
                        {
			icount ++;
                        gx = gz;
                        for (j=0;j<numargs;j++) x[j]=z[j];
                        }
                else { eps = 2.0*eps/3.0; if (eps < 1.0e-4) icount = imax; }
                }
#ifdef VERBOSE
	for (j=0;j<numargs;j++) cout << x[j] << " "; cout << endl;
	cout << " max = " << gx << endl << flush;
#endif
        return gx;
        }

valuator* globalHACK;

static double HACKfunction(double x[])
	{
	return globalHACK->evalf(x);
	}

double unconstrained_iterated_max(int num_iter,double xmin[],double xmax[],
		double x[],int numargs,valuator& u)
	{
	globalHACK=&u;
	return unconstrained_iterated_max(num_iter,xmin,xmax,x,numargs,HACKfunction);
	}

double unconstrained_iterated_max(int num_iter,double xmin[],double xmax[],
		double x[],int numargs,double (*func)(double x[]) )
	{
	if (numargs>MAXARG) { cout << "array overflow"; return 0.0; }
	double z[MAXARG];
        double gx = unconstrained_max(xmin,xmax,x,numargs,*func),gz;
        for (int i = 0; i<num_iter; i++)
                { 
		gz = unconstrained_max(xmin,xmax,z,numargs,*func);
                if (gz > gx)
                  {
                  gx = gz;
		  for (int j=0;j<numargs;j++) x[j]=z[j];
                  }
                }
        return gx;
        }
