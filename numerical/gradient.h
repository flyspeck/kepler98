#ifndef grad_c
#define grad_c

#include <iomanip.h>
#include <iostream.h>
#include <fstream.h>
#include <math.h>
#include <stdlib.h>
#include "numerical.h"

class valuator {
public:
		virtual double evalf(double y[])=0;
		};

double gradient(int ith,
        double x[],int numargs,double (*func)(double x[])  );

double second_partial(int i,int j,int numargs,double (*f)(double[]),double x[]);

double unconstrained_max(double xmin[],double xmax[],
                double x[],int numargs,double (*func)(double x[]) );

double unconstrained_iterated_max(int num_iter,double xmin[],double xmax[],
                double x[],int numargs,double (*func)(double x[]) );

double unconstrained_iterated_max(int num_iter,double xmin[],double xmax[],
                double x[],int numargs,valuator& );

#endif

/* SAMPLE:

double gradfunction(double x[6])
        {
        return gamma(x[0],x[1],x[2],x[3],x[4],x[5]);
        }

int main()
        {
        double y[6];
        double ymin[6] = {2.0,2.0,2.0,2.0,2.0,2.0};
        double ymax[6] = {2.51,2.51,2.51,2.51,2.51,2.51};
        cout << unconstrained_iterated_max(30,ymin,ymax,y,6,gradfunction);
        cout << endl;
        for (int i=0;i<6;i++) cout << y[i] << " ";
        cout << endl;
        }

*/
