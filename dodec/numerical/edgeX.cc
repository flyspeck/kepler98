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


	  //Solid edge
	
	case 1 : *ret = solid(x[0],x[1],x[2],x[3],x[4],x[5])+ 
		   .245*(x[0]+x[1]+x[2]-6)-.063*(x[3]+x[4]+x[5]-6)
		   - 0.551285;break;

	

	case 2 : *ret = solid(x[0],x[1],x[2],x[3],x[4],x[5])+ 
		   0.3798*(x[0]+x[1]+x[2]-6)-.198*(x[3]+x[4]+x[5]-6)
		   - 0.551285;break;

	case 3 : *ret = -1*solid(x[0],x[1],x[2],x[3],x[4],x[5]) 
		   -0.151*(x[0]+x[1]+x[2]-6)+0.323*(x[3]+x[4]+x[5]-6)
		   + 0.551286;break;

	

		   // Dih edge

        
	case 4 : *ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
		   -0.237*(x[0]-2) - 0.708 *(x[3]-2)+
		   .372*(x[1]+x[2]+x[4]+x[5]-8) - 1.23095;break;		   
	case 5 : *ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
		   -0.237*(x[0]-2) - 0.688 * (x[3]-2)+
		   .363*(x[1]+x[2]+x[4]+x[5]-8) - 1.23095;break;		   
	case 6 : *ret = -1*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
		   + 0.505* (x[0]-2) + .766 * (x[3]-2)
		   -0.152*(x[1]+x[2]+x[4]+x[5]-8) + 1.23096;break;		   
		   // Vol and Squander edge


	case 7 : *ret = squander(x[0],x[1],x[2],x[3],x[4],x[5])
		   - 0.0392* (x[0]+x[1]+x[2]-6) 
		   -0.0101 * (x[3]+x[4]+x[5]-6);break;  
		   
	case 8 : *ret = squander2(x[0],x[1],x[2],x[3],x[4],x[5])
		   - 0.0392* (x[0]+x[1]+x[2]-6) 
		   -.0101* (x[3]+x[4]+x[5]-6);break; 
		   
	
	case 9 : *ret = volume(x[0],x[1],x[2],x[3],x[4],x[5])+
		   .107 * (x[0]+x[1]+x[2]-6)- 0.116 *(x[3]+x[4]+x[5]-6)
		   -0.235702;break;

	case 10 : *ret = volume(x[0],x[1],x[2],x[3],x[4],x[5])+
		   0.0623 * (x[0]+x[1]+x[2]-6) -0.0722 *(x[3]+x[4]+x[5]-6)
		   -0.235702;break;	

	
		   // Quad Vol edge


        case 11 : *ret = QuadVarSquander(0,0,x[0],x[1],x[2],x[3],x[4],x[5]
                    ,x[6],x[7],x[8])
		    +0.166 * (x[0]+x[1]+x[2]+x[6]-8)
		    -.143 *(x[4]+x[5]+x[7]+x[8]-8) -0.590491;break;   
	
	
	case 12 : *ret = truncVol(x[0],x[1],x[2],x[3],x[4],x[5])		                   +0.166 * (x[0]+x[1]/2+x[2]/2 - 4)
		    -.143*(x[4]+x[5]-4) -0.590491/2;break;   	    


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
	

	
	case 11 : *ret = 3 - crossdiag(x) ; break;
               
	
	
	default : cout << "unexpected case in constraint" << endl;
		
	  


	}
    }

iter::iter(int ineqSwitch) {
	numiter = 20; numargs = 6; nconstr=0;
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
		
		case 11 : xmin[3]=3; xmax[3]=3.56;break;
		 
		case 12 : xmin[3]= 2.51682; xmax[3]=3;break;   

		
		  

		default : cout << "error " << ineqSwitch << ": not installed " << endl;
		}
	}

void /*ineq.cc*/minimize2(int);

void page1()
{ 

  minimize2(12);

  /*
  minimize2(1);
  minimize2(2);
  minimize2(3);
  minimize2(4);
  minimize2(5);
  minimize2(6);
  minimize2(7);
  minimize2(8);
  minimize2(9);
  minimize2(10);
  minimize2(11);
  minimize2(12);
  */
}








