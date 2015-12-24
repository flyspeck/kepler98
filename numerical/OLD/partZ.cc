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


static void nofunc(int numargs,int whichFn,double* x,double* ret,void*)
    {
    cout << "nofunc should not be called" << endl << flush;
    }

int INEQ_NUMBER=0;
static void generic(int numargs,int whichFn,double* x, double* ret,void*)
	{
	switch (INEQ_NUMBER) {
		// start of first page of inequalities for Section 2, SPIV.
		case 809220761 : 
			*ret = tau(x[0],x[1],x[2],x[3],x[4],x[5])-3.07*global::pt; 
			break;
		case 809220762 : case 809220763 :
			*ret = tau_analytic(x[0],x[1],x[2],x[3],x[4],x[5])-3.07*global::pt; 
			break;
		case 809220764 : 
			*ret = tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-3.07*global::pt; 
			break;
		case 809220765 : 
			*ret = -0.035+tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-3.07*global::pt;
			 break;

		case 316093823 :
			*ret= -1.4*global::pt+tau(x[0],x[1],x[2],x[3],x[4],x[5])+
					tau(x[0],x[2],x[7],x[6],x[8],x[4])+
					0*tau(x[0],x[7],x[11],x[9],x[10],x[8])+
					tau(x[0],x[11],x[13],x[12],x[14],x[10])+
					0*tau(x[0],x[1],x[13],x[15],x[14],x[5]);
				break;

		case 170976188 :
			*ret= -1.4*global::pt+tau(x[0],x[1],x[2],x[3],x[4],x[5])+
					tau(x[0],x[2],x[7],x[6],x[8],x[4])+
					tau(x[0],x[7],x[11],x[9],x[10],x[8])+
					0*tau(x[0],x[11],x[13],x[12],x[14],x[10])+
					0*tau(x[0],x[1],x[13],x[15],x[14],x[5]);
				break;

		case 671740346 :
			*ret= -1.5*global::pt+0*tau(x[0],x[1],x[2],x[3],x[4],x[5])+
					tau(x[0],x[2],x[7],x[6],x[8],x[4])+
					tau(x[0],x[7],x[11],x[9],x[10],x[8])+
					tau(x[0],x[11],x[13],x[12],x[14],x[10])+
					tau(x[0],x[1],x[13],x[15],x[14],x[5]);
				break;
		case 22430965 :
			*ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.32; break;

		case 317122931 :
			*ret= -1.4*global::pt+
					tau(x[0],x[1],x[2],x[3],x[4],x[5])+
					tau(x[0],x[2],x[7],x[6],x[8],x[4])+
					tau(x[0],x[7],x[11],x[9],x[10],x[8])+
					tau(x[0],x[11],x[13],x[12],x[14],x[10])+
					0*tau(x[0],x[1],x[13],x[15],x[14],x[5]);
				break;

		case 317122932 : case 317122933 :
			*ret= -1.4*global::pt+
					tau(x[0],x[1],x[2],x[3],x[4],x[5])+
					tau(x[0],x[2],x[7],x[6],x[8],x[4])+
					tau(x[0],x[7],x[11],x[9],x[10],x[8])+
					tau_analytic(x[0],x[11],x[13],x[12],x[14],x[10])+
					0*tau(x[0],x[1],x[13],x[15],x[14],x[5]);
				break;

		case 317122934 : 
			*ret= -1.4*global::pt+
					tau(x[0],x[1],x[2],x[3],x[4],x[5])+
					tau(x[0],x[2],x[7],x[6],x[8],x[4])+
					tau(x[0],x[7],x[11],x[9],x[10],x[8])+
					tauVc(x[0],x[11],x[13],x[12],x[14],x[10])+
					0*tau(x[0],x[1],x[13],x[15],x[14],x[5]);
				break;

		case 317122935 : 
			*ret= -1.4*global::pt+
					tau(x[0],x[1],x[2],x[3],x[4],x[5])+
					tau(x[0],x[2],x[7],x[6],x[8],x[4])+
					tau(x[0],x[7],x[11],x[9],x[10],x[8])+
					tauVc(x[0],x[11],x[13],x[12],x[14],x[10])-2.0*(0.00935)+
					0*tau(x[0],x[1],x[13],x[15],x[14],x[5]);
				break;


		case 359077171 :
			*ret= -1.4*global::pt+
					tau(x[0],x[1],x[2],x[3],x[4],x[5])+
					tau(x[0],x[2],x[7],x[6],x[8],x[4])+
					tau(x[0],x[7],x[11],x[9],x[10],x[8])+
					tau(x[0],x[11],x[13],x[12],x[14],x[10])+
					0*tau(x[0],x[1],x[13],x[15],x[14],x[5]);
				break;

		case 359077172 : case 359077173 :
			*ret= -1.4*global::pt+
					tau(x[0],x[1],x[2],x[3],x[4],x[5])+
					tau(x[0],x[2],x[7],x[6],x[8],x[4])+
					tau_analytic(x[0],x[7],x[11],x[9],x[10],x[8])+
					tau(x[0],x[11],x[13],x[12],x[14],x[10])+
					0*tau(x[0],x[1],x[13],x[15],x[14],x[5]);
				break;

		case 359077174 : 
			*ret= -1.4*global::pt+
					tau(x[0],x[1],x[2],x[3],x[4],x[5])+
					tau(x[0],x[2],x[7],x[6],x[8],x[4])+
					tauVc(x[0],x[7],x[11],x[9],x[10],x[8])+
					tau(x[0],x[11],x[13],x[12],x[14],x[10])+
					0*tau(x[0],x[1],x[13],x[15],x[14],x[5]);
				break;

		case 359077175 : 
			*ret= -1.4*global::pt+
					tau(x[0],x[1],x[2],x[3],x[4],x[5])+
					tau(x[0],x[2],x[7],x[6],x[8],x[4])+
					tauVc(x[0],x[7],x[11],x[9],x[10],x[8])+
					tau(x[0],x[11],x[13],x[12],x[14],x[10])-2.0*(0.00935)+
					0*tau(x[0],x[1],x[13],x[15],x[14],x[5]);
				break;

		case 103362341 :
			*ret= -1.5*global::pt+
					tau(x[0],x[1],x[2],x[3],x[4],x[5])+
					tau(x[0],x[2],x[7],x[6],x[8],x[4])+
					tau(x[0],x[7],x[11],x[9],x[10],x[8])+
					tau(x[0],x[11],x[13],x[12],x[14],x[10])+
					tau(x[0],x[1],x[13],x[15],x[14],x[5]);
				break;

		case 103362342 : case 103362343 :
			*ret= -1.5*global::pt+
					tau_analytic(x[0],x[1],x[2],x[3],x[4],x[5])+
					tau(x[0],x[2],x[7],x[6],x[8],x[4])+
					tau(x[0],x[7],x[11],x[9],x[10],x[8])+
					tau(x[0],x[11],x[13],x[12],x[14],x[10])+
					tau(x[0],x[1],x[13],x[15],x[14],x[5]);
				break;

		case 103362344 : 
			*ret= -1.5*global::pt+
					tauVc(x[0],x[1],x[2],x[3],x[4],x[5])+
					tau(x[0],x[2],x[7],x[6],x[8],x[4])+
					tau(x[0],x[7],x[11],x[9],x[10],x[8])+
					tau(x[0],x[11],x[13],x[12],x[14],x[10])+
					tau(x[0],x[1],x[13],x[15],x[14],x[5]);
				break;

		case 103362345 : 
			*ret= -1.5*global::pt+
					tauVc(x[0],x[1],x[2],x[3],x[4],x[5])+
					tau(x[0],x[2],x[7],x[6],x[8],x[4])+
					tau(x[0],x[7],x[11],x[9],x[10],x[8])+
					tau(x[0],x[11],x[13],x[12],x[14],x[10])-2.0*(0.00935)+
					tau(x[0],x[1],x[13],x[15],x[14],x[5]);
				break;



		default : cout << INEQ_NUMBER << " generic default : function not installed " << endl << flush;
			*ret=0;
		}
	}

iter::~iter() 
	{ delete[] xmin; 
	  delete[] xmax; 
		delete[] x; }

double dih15(double x[15])
	{
	return dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])+
			dihedraly(x[0],x[2],x[7],x[6],x[8],x[4])+
			dihedraly(x[0],x[7],x[11],x[9],x[10],x[8])+
			dihedraly(x[0],x[11],x[13],x[12],x[14],x[10])+
			dihedraly(x[0],x[1],x[13],x[15],x[14],x[5]);
	}



static void ConstraintPage1(int numargs,int whichFn,double* x,double* ret,void*)
    {
	*ret = 0;
	switch (INEQ_NUMBER) {
		case 809220761 : switch(whichFn) {
			case 1 : *ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.32;
				break;
			case 2 : *ret = radf(x[3],x[4],x[5])-global::sqrt2; break;
			case 3 : *ret = radf(x[1],x[2],x[3])-global::sqrt2; break;
			}
			break;
		case 809220762 :
		case 809220765 :
			switch(whichFn) {
			case 1 : *ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.32;
				break;
			case 2 : *ret = -radf(x[3],x[4],x[5])+global::sqrt2; break;
			}
			break;
		case 809220763 :
			switch(whichFn) {
			case 1 : *ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.32;
				break;
			case 2 : *ret = -radf(x[1],x[2],x[3])+global::sqrt2; break;
			}
			break;
		case 809220764 :
			*ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.32; 
			break;
		case 316093823 : case 170976188 : case 671740346 :
			*ret = -global::pi*2.0+dih15(x); break;

		case 317122931 : switch(whichFn) {
			case 1 : *ret = dih15(x)-2.0*global::pi;
				break;
			case 2 : *ret = -dihedraly(x[0],x[13],x[1],x[15],x[5],x[14])+1.32;
			case 3 : *ret = radf(x[10],x[12],x[14])-global::sqrt2; break;
			case 4 : *ret = radf(x[11],x[12],x[13])-global::sqrt2; break;
			}
			break;

		case 317122932 : case 317122935 : switch(whichFn) {
			case 1 : *ret = dih15(x)-2.0*global::pi;
				break;
			case 2 : *ret = -dihedraly(x[0],x[13],x[1],x[15],x[5],x[14])+1.32;
			case 3 : *ret = -radf(x[10],x[12],x[14])+global::sqrt2; break;
			}
			break;

		case 317122933 : switch(whichFn) {
			case 1 : *ret = dih15(x)-2.0*global::pi;
				break;
			case 2 : *ret = -dihedraly(x[0],x[13],x[1],x[15],x[5],x[14])+1.32;
			case 3 : *ret = -radf(x[11],x[12],x[13])+global::sqrt2; break;
			}
			break;

		case 317122934 : switch(whichFn) {
			case 1 : *ret = dih15(x)-2.0*global::pi;
				break;
			case 2 : *ret = -dihedraly(x[0],x[13],x[1],x[15],x[5],x[14])+1.32;
			}
			break;

		case 359077171 : switch(whichFn) {
			case 1 : *ret = dih15(x)-2.0*global::pi;
				break;
			case 2 : *ret = -dihedraly(x[0],x[13],x[1],x[15],x[5],x[14])+1.32;
			case 3 : *ret = radf(x[8],x[9],x[10])-global::sqrt2; break;
			case 4 : *ret = radf(x[7],x[9],x[11])-global::sqrt2; break;
			}
			break;

		case 359077172 : case 359077175 : switch(whichFn) {
			case 1 : *ret = dih15(x)-2.0*global::pi;
				break;
			case 2 : *ret = -dihedraly(x[0],x[13],x[1],x[15],x[5],x[14])+1.32;
			case 3 : *ret = -radf(x[8],x[9],x[10])+global::sqrt2; break;
			}
			break;

		case 359077173 : switch(whichFn) {
			case 1 : *ret = dih15(x)-2.0*global::pi;
				break;
			case 2 : *ret = -dihedraly(x[0],x[13],x[1],x[15],x[5],x[14])+1.32;
			case 3 : *ret = -radf(x[7],x[9],x[11])+global::sqrt2; break;
			}
			break;

		case 359077174 : switch(whichFn) {
			case 1 : *ret = dih15(x)-2.0*global::pi;
				break;
			case 2 : *ret = -dihedraly(x[0],x[13],x[1],x[15],x[5],x[14])+1.32;
			}
			break;

		case 103362341 : switch(whichFn) {
			case 1 : *ret = dih15(x)-2.0*global::pi;
				break;
			case 2 : *ret = radf(x[3],x[4],x[5])-global::sqrt2; break;
			case 3 : *ret = radf(x[1],x[2],x[3])-global::sqrt2; break;
			}
			break;

		case 103362342 : case 103362345 : switch(whichFn) {
			case 1 : *ret = dih15(x)-2.0*global::pi;
				break;
			case 2 : *ret = -radf(x[3],x[4],x[5])+global::sqrt2; break;
			}
			break;

		case 103362343 : switch(whichFn) {
			case 1 : *ret = dih15(x)-2.0*global::pi;
				break;
			case 2 : *ret = -radf(x[1],x[2],x[3])+global::sqrt2; break;
			}
			break;

		case 103362344 : switch(whichFn) {
			case 1 : *ret = dih15(x)-2.0*global::pi;
				break;
			}
			break;

		default : cout << INEQ_NUMBER << " constraints not installed" << endl;
		}
    }

iter::iter(int ineqSwitch) {
	numiter = 20; numargs = 6; nconstr=0;
	switch(ineqSwitch)
		{
		case 316093823 :
		case 170976188 :
		case 671740346 : numargs=16; break;

		case 317122931 :
		case 317122932 :
		case 317122933 :
		case 317122934 :
		case 317122935 : numargs=16; break;

		case 359077171 :
		case 359077172 :
		case 359077173 :
		case 359077174 :
		case 359077175 : numargs=16; break;

		case 103362341 :
		case 103362342 :
		case 103362343 :
		case 103362344 :
		case 103362345 : numargs=16; break;
		};

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
		case 809220761 :
			xmin[3]=2.51; xmax[3]=global::sqrt8;
			nconstr=3;
			break;
		case 809220762 : case 809220763 : case 809220765 :
			xmin[3]=2.51; xmax[3]=global::sqrt8;
			nconstr=2;
			break;
		case 809220764 : 
			xmin[3]=2.6; xmax[3]=global::sqrt8;
			nconstr= 1;
			break;

		case 316093823 :
			xmin[9]=xmax[9]=global::sqrt8;
			xmin[15]=xmax[15]=global::sqrt8;
			nconstr=1;
			break;
		case 170976188 :
			xmin[12]=xmax[12]=global::sqrt8;
			xmin[15]=xmax[15]=global::sqrt8;
			nconstr=1;
			break;
		case 22430965 :
			xmin[3]=xmax[3]=global::sqrt8;
			break;
		case 671740346 :
			xmin[3]=xmax[3]=global::sqrt8;
			nconstr=1;
			break;

		case 317122931 :
			xmin[12]=2.51; xmax[12]=global::sqrt8;
			xmin[15]=2.51; xmax[15]=3.2;
			nconstr=4;
			break;
		case 317122932 : case 317122933 : case 317122935 :
			xmin[12]=2.51; xmax[12]=global::sqrt8;
			xmin[15]=2.51; xmax[15]=3.2;
			nconstr=3;
			break;
		case 317122934 : 
			xmin[12]=2.6; xmax[12]=global::sqrt8;
			xmin[15]=2.51; xmax[15]=3.2;
			nconstr=2;
			break;

		case 359077171 :
			xmin[9]=2.51; xmax[9]=global::sqrt8;
			xmin[15]=2.51; xmax[15]=3.2;
			nconstr=4;
			break;
		case 359077172 : case 359077173 : case 359077175 :
			xmin[9]=2.51; xmax[9]=global::sqrt8;
			xmin[15]=2.51; xmax[15]=3.2;
			nconstr=3;
			break;
		case 359077174 : 
			xmin[9]=2.6; xmax[9]=global::sqrt8;
			xmin[15]=2.51; xmax[15]=3.2;
			nconstr=2;
			break;

		case 103362341 :
			xmin[3]=2.51; xmax[3]=global::sqrt8;
			nconstr=3;
			break;
		case 103362342 : case 103362343 : case 103362345 :
			xmin[3]=2.51; xmax[3]=global::sqrt8;
			nconstr=2;
			break;
		case 103362344 : 
			xmin[3]=2.6; xmax[3]=global::sqrt8;
			nconstr=1;
			break;


		


		default : cout << "error " << ineqSwitch << ": bounds not installed " << endl;
		}
		if (nconstr>0) constraintfunc=ConstraintPage1;
	}

void /*ineq.cc*/minimize2(int);

void page1()
	{
	// test code goes here....

	
	minimize2(809220761);
	minimize2(809220762);
	minimize2(809220763);
	minimize2(809220764);
	minimize2(809220765);

	return;
	minimize2(316093823);
	minimize2(170976188);
	minimize2(22430965);
	minimize2(671740346 );

	minimize2(317122931);
	minimize2(317122932);
	minimize2(317122933);
	minimize2(317122934);
	minimize2(317122935);

	minimize2(359077171);
	minimize2(359077172);
	minimize2(359077173);
	minimize2(359077174);
	minimize2(359077175);
	
	minimize2(103362341);
	minimize2(103362342);
	minimize2(103362343);
	minimize2(103362344);
	minimize2(103362345);

	}
