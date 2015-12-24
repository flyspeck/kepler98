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
	
	case 1 : case 3 :

	  *ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.705;break;

	case 2 : 

	  *ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.681;break;


	case 4 : 

	  *ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-2.652;break;

	case 5 : case 7 :

	  *ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.681;break;

	case 6 : case 8 :

	  *ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-2.074;break;

	case 9 : case 14 : case 15 : break;


	case 10 : case 12 :

	  *ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.664;break;

	case 11 : case 13 :

	  *ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.917;break;
	  
	case 16 : case 17 : switch(whichFn)
	  {
	  case 1 : *ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.618;break;
	    
	  case 2 : *ret = 2.51681-crossdiag(x);break;
	  }
	break;



	case 18 : switch(whichFn)
	  {
	  case 2 : *ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-2.074;break;
	  case 1 : *ret = 2.51681-crossdiag(x);break;
	  }
	break;

	default : cout << "unexpected case in constraint" << endl;
		
	  


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
	

	for (i=0;i<numargs;i++) { xmin[i]=x[i]=2.0; xmax[i]=2.51682; }
	INEQ_NUMBER = ineqSwitch;
	switch (ineqSwitch)
		{
		
		case 1 :

		xmax[0]=2.329;xmax[1]=2.207;xmax[2]=2.51682;
		xmin[3]=2.51682;xmax[3]=3.4;break;
		
		case 2 :

		xmax[0]=2.161;xmax[1]=2.207;xmax[2]=2.280;
		xmin[3]=2.51682;xmax[3]=3.4;break;

		case 3 : xmax[0]=2.207;xmax[1]=2.51682;xmax[2]=2.219;
                xmin[3]=2.51682;xmax[3]=3.373;break;

		case 4 : xmax[2]=2.219;xmin[3]=2.51682;xmax[3]=3.373;break;

		case 5 : xmax[0]=2.161;xmax[1]=2.207;xmax[2]=2.280;
		  xmin[3]=2.51682;xmax[3]=3.6;break;

		case 6 : xmax[0]=2.33;xmax[2]=2.207;xmin[3]=2.51682;
		         xmax[3]=3.9;break;

		case 7 : xmax[0]=2.161;xmax[1]=2.207;xmax[2]=2.280;
                  xmin[3]=2.51682;xmax[3]=3.189;break;

		case 8 : xmax[0]=2.33;xmax[2]=2.207;xmin[3]=2.51682;
		         xmax[3]=3.827;break;

		case 9 : xmax[0]=2.207;xmax[1]=2.28;xmin[4]=2.51682;
                         xmax[4]=3.6;xmin[5]=2.51682;xmax[5]=3.189;break;



		case 10 : xmax[0]=2.149;xmax[1]=2.164;xmax[2]=2.222;
                  xmin[3]=2.51682;xmax[3]=3.6;break;

                case 11 : xmax[0]=2.31;xmax[2]=2.164;xmin[3]=2.51682;
		  xmax[3]=3.6;break;

                case 12 : xmax[0]=2.149;xmax[1]=2.164;xmax[2]=2.222;
                  xmin[3]=2.51682;xmax[3]=3.112;break;

                case 13 : xmax[0]=2.31;xmax[2]=2.164;xmin[3]=2.51682;
                         xmax[3]=3.6;break;

                case 14 : xmax[0]=2.164;xmax[1]=2.222;xmin[4]=2.51682;
                         xmax[4]=3.6;xmin[5]=2.51682;xmax[5]=3.112;break;


		case 15 : xmax[0]=2.207;xmax[1]=2.28;xmin[4]=3.6;
		  xmax[4]=3.827;xmin[5]=2.51682;xmax[5]=3.189;break;

		case 16 : case 17 : xmax[0]=2.161;xmax[1]=2.207;xmax[2]=2.28;
		  xmin[3]=2.51682;xmax[3]=3.189;xmin[7]=3.6;xmax[7]=3.827;
		  break;

		case 18 : xmax[0]=2.33;xmax[1]=2.207;xmin[3]=2.51681;
		  xmax[3]=3.827;xmax[6]=2.28;xmin[8]=2.51682;
		  xmax[8]=3.189;break;



		default : cout << "error " << ineqSwitch << ": not installed " << endl;
		}
	}

void /*ineq.cc*/minimize2(int);

void page1()
{ 
  

  minimize2(1);
  minimize2(2); 
  //minimize2(14);

  /*

  minimize2(12);
  minimize2(9);
  minimize2(10);
  minimize2(11);
  minimize2(12);

  
  minimize2(1);
  minimize2(2);
  minimize2(3);
  minimize2(4);
  */
}








