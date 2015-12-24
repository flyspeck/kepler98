#include <iomanip.h>
#include <stdlib.h>
#include "numerical.h"
#include "numerical/constants.h"
#include "numerical/morefn.h"
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

	case 1 : case 2 : 

	  *ret = -x[3];break;

	case 3 : case 4 : case 5:

	  *ret = truncVol(x[0],x[1],x[2],x[3],x[4],x[5]);break;




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
	
	case 1 : case 3 : case 5 :

	  *ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.705;break;


	case 5 : break;


	default : cout << "unexpected case in constraint" << endl;
		
	  


	}
    }

iter::iter(int ineqSwitch) {
	numiter = 20; numargs = 6; nconstr=1;
	// temp: 
	
	
	xmin = new double[numargs];
	xmax = new double[numargs];
	x = new double[numargs];
	constraintfunc = ConstraintPage1;
	func = generic;
	int i;
	

	for (i=0;i<numargs;i++) { xmin[i]=x[i]=2.0; xmax[i]=2.51682; }
	INEQ_NUMBER = ineqSwitch;
	switch (ineqSwitch)
		{
		
		case 1 :

		xmax[0]=2.207;xmax[1]=2.199;xmax[2]=2.222;
		xmin[3]=2.51682;xmax[3]=3.4;break;
		
		case 2 :

		xmax[0]=2.161;xmax[1]=2.207;xmax[2]=2.280;
		xmin[3]=2.51682;xmax[3]=3.4;break;

		case 3 : xmax[0]=2.207;xmax[1]=2.51682;xmax[2]=2.219;
                xmin[3]=2.51682;xmax[3]=3.373;break;

		case 4 : xmax[2]=2.219;xmin[3]=2.51682;xmax[3]=3.373;break;

		case 5 : xmax[0]=2.161;xmax[1]=2.207;xmax[2]=2.280;
		  xmin[3]=2.51682;xmax[3]=3.6;break;


		default : cout << "error " << ineqSwitch << ": not installed " << endl;
		}
	}

void /*ineq.cc*/minimize2(int);

void page1()
{ 
  

  minimize2(1);

}








