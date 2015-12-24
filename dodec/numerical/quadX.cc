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

	case 1 : *ret = QuadVarSquander
	(0.42775, 0.15098,x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;

	case 2 : *ret = QuadVarSquander2
       (0.42775, 0.15098, x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;

	case 3 : *ret = QuadVarSquander
       (0.42775, 0.09098, x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;
			    
        case 4 : *ret = QuadVarSquander2
       (0.42775, 0.09098, x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;
			 
      	case 5 : *ret = QuadVarSquander
       (0.42775, 0, x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;
			    
        case 6 : *ret = QuadVarSquander2
       (0.42775, 0, x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;

	case 7 : *ret = QuadVarSquander
      (0.42775, -0.18519, x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]);break;
			    
        case 8 : *ret = QuadVarSquander2
      (0.42775, -0.18519, x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;
		
        case 9 : *ret = QuadVarSquander
      (0.42775, -0.20622, x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]);break;
			    
        case 10 : *ret = QuadVarSquander2
      (0.42775, -0.20622, x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;

	case 11 : *ret = QuadVarSquander
       (0.55792, 0.30124, x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]);break;
		
        case 12 : *ret = QuadVarSquander2
       (0.55792, 0.30124, x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;
			 
      	case 13 : *ret = QuadVarSquander
       (0.55792, 0.02921, x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]);break;
			    
        case 14 : *ret = QuadVarSquander2
       (0.55792, 0.02921, x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;

	case 15 : *ret = QuadVarSquander
      (0.55792, 0, x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;
			    
        case 16 : *ret = QuadVarSquander2
      (0.55792, 0, x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;
		
	case 17 : *ret = QuadVarSquander
      (0.55792, -0.05947, x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]);break;
			    
        case 18 : *ret = QuadVarSquander2
      (0.55792, -0.05947, x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;
			 
      	case 19 : *ret = QuadVarSquander
      (0.55792, -0.39938, x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]);break;
			    
        case 20 : *ret = QuadVarSquander2
      (0.55792, -0.39938, x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;

	case 21 : *ret = QuadVarSquander
      (0.55792, -2.5021, x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;
			    
        case 22 : *ret = QuadVarSquander2
      (0.55792, -2.5021, x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;

	case 23 : *ret = QuadVarSquander
       (0.68, 0.44194, x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;
			    
        case 24 : *ret = QuadVarSquander2
       (0.68, 0.44194, x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;

	case 25 : *ret = QuadVarSquander
       (0.68, 0.10957, x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;
			    
        case 26 : *ret = QuadVarSquander2
       (0.68, 0.10957, x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;
			 
      	case 27 : *ret = QuadVarSquander
       (0.68, 0, x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;	
		    
        case 28 : *ret = QuadVarSquander2
       (0.68, 0, x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;

	case 29 : *ret = QuadVarSquander
      (0.68, -0.86096, x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;
			    
        case 30 : *ret = QuadVarSquander2
      (0.68,-0.86096 , x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;
		
        case 31 : *ret = QuadVarSquander
      (0.68, -2.44439, x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;
			    
        case 32 : *ret = QuadVarSquander2
      (0.68, -2.44439, x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;

	case 33 : *ret = QuadVarSquander
       (0.3,0.12596 ,x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;
			    
	case 34 : *ret = QuadVarSquander2
       (0.3,0.12596 ,x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;

	case 35 : *ret = QuadVarSquander
       (0.3,0.02576 ,x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;

	case 36 : *ret = QuadVarSquander2
       (0.3,0.02576 ,x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;

	case 37 : *ret = QuadVarSquander
       (0.3,0 ,x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;

	case 38 : *ret = QuadVarSquander2
       (0.3,0 ,x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;

	case 39 : *ret = QuadVarSquander
       (0.3,-0.037 ,x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;

	case 40 : *ret = QuadVarSquander2
       (0.3,-0.037 ,x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;

	case 41 : *ret = QuadVarSquander
       (0.3,-0.22476 ,x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;

	case 42 : *ret = QuadVarSquander2
       (0.3,-0.22476 ,x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;

	case 43 : *ret = QuadVarSquander
       (0.3,-2.31852 ,x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;

	case 44 : *ret = QuadVarSquander2
       (0.3,-2.31852 ,x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;

	case 45 : *ret = QuadVarSquander
       (0,0.23227 ,x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;

	case 46 : *ret = QuadVarSquander2
       (0,0.23227 ,x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;

	case 47 : *ret = QuadVarSquander
       (0,-0.07448 ,x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;

	case 48 : *ret = QuadVarSquander2
       (0,-0.07448 ,x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;

	case 49 : *ret = QuadVarSquander
       (0,-0.22019 ,x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;

	case 50 : *ret = QuadVarSquander2
       (0,-0.22019 ,x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;

	case 51 : *ret = QuadVarSquander
       (0,-0.80927 ,x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;
				    // checking sigma ineq
	case 52 : *ret = QuadVarSquander2
       (0, -0.80927, x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;

	case 53 : *ret = QuadVarSquander
       (0,-5.8438 ,x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;

	case 54 : *ret = QuadVarSquander2
       (0,-5.8438 ,x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); break;



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
	
	case 1 : case 2 : case 3 : case 4 : case 5 : case 6 : case 7 : 
	case 8 : case 9 : case 10 : case 11 : case 12 : case 13 : case 14 : 
	case 15 : case 16 : case 17 : case 18 : case 19 : case 20 : case 21 : 
	case 22 : case 23 : case 24 : case 25 : case 26 : case 27 : case 28 : 
	case 29 : case 30 : case 31 : case 32 : case 33 : case 34 : case 35 : 
	case 36 : case 37 : case 38 : case 39 : case 40 : case 41 : case 42 : 
	case 43 : case 44 : case 45 : case 46 : case 47 : case 48 : case 49 : 
	case 50 : case 51 : case 52 : case 53 : case 54 : 

	  *ret = 2.51682 - crossdiag(x); break;

	 

	default : cout << "unexpected case in constraint" << endl;
		
	  


	}
    }

iter::iter(int ineqSwitch) {
	numiter = 20; numargs = 9; nconstr=1;
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
		
		case 1 : case 2 : case 3 : case 4 : case 5 : case 6 : case 7 :
		case 8 : case 9 : case 10 : case 11 : case 12 : case 13 :
		case 14 : case 15 : case 16 : case 17 : case 18 : case 19 : 
		case 20 : case 21 : case 22 : case 23 : case 24 : case 25 : 
		case 26 : case 27 : case 28 : case 29 : case 30 : case 31 : 
		case 32 : case 33 : case 34 : case 35 : case 36 : case 37 : 
		case 38 : case 39 : case 40 : case 41 : case 42 : case 43 : 
		case 44 : case 45 : case 46 : case 47 : case 48 : case 49 : 
		case 50 : case 51 : case 52 : case 53 : case 54 : 

		  xmin[3]=2.51682; xmax[3]=3.55933; break;
		 


      default : cout << "error " << ineqSwitch << ": not installed " << endl;
		
		}
	}

void /*ineq.cc*/minimize2(int);

void page1()
{ 
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
  minimize2(13);
  minimize2(14);
  minimize2(15);
  minimize2(16);
  minimize2(17);
  minimize2(18);
  minimize2(19);
  minimize2(20);
  minimize2(21);
  minimize2(22);
  minimize2(23);
  minimize2(24);
  minimize2(25);
  minimize2(26);
  minimize2(27);
  minimize2(28);
  minimize2(29);
  minimize2(30);
  minimize2(31);
  minimize2(32);
  minimize2(33);
  minimize2(34);
  minimize2(35);
  minimize2(36);
  minimize2(37);
  minimize2(38);
  minimize2(39);
  minimize2(40);
  minimize2(41);
  minimize2(42);
  minimize2(43);
  minimize2(44);
  minimize2(45);
  minimize2(46);
  minimize2(47);
  minimize2(48);
  minimize2(49);
  minimize2(50);
  minimize2(51);
  minimize2(52);
  minimize2(53);
  minimize2(54);   
}








