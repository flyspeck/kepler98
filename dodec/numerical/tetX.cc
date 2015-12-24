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

		
	case 1 : *ret = VarSquander(0.68, -1.88718,
				    x[0],x[1],x[2],x[3],x[4],x[5]);break;

	case 2 : *ret = VarSquander2(0.68, -1.88718,
				     x[0],x[1],x[2],x[3],x[4],x[5]);break;
	


	case 3 : *ret =  VarSquander(0.68, -0.90746, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;		    
	case 4 : *ret =  VarSquander2(0.68, -0.90746, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;



	case 5 : *ret =  VarSquander(0.68, -0.46654, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;

	case 6 : *ret =  VarSquander2(0.68, -0.46654, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;
				    

     
	case 7 : *ret =  VarSquander(0.55889, 0, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;

	case 8 : *ret =  VarSquander2(0.55889, 0, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;


	
	case 9 : *ret =  VarSquander(0.63214, 0, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;
	
	case 10 : *ret =  VarSquander2(0.63214, 0, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;		    

	
	case 11 : *ret =  VarSquander(0.73256, 0, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;		
	case 12 : *ret =  VarSquander2(0.73256, 0, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;	    

	
	case 13 : *ret =  VarSquander(0.89346, 0, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;		
	case 14 : *ret =  VarSquander2(0.89346, 0, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;		    
	

	case 15 : *ret =  VarSquander(0.3, 0.5734, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;		
	case 16 : *ret =  VarSquander2(0.3, 0.5734, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;
	


	case 17 : *ret =  VarSquander(0.3, 0.03668, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			
	case 18 : *ret =  VarSquander2(0.3, 0.03668, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;		    
	

	case 19 : *ret =  VarSquander(0.3, -0.04165, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			   
	case 20 : *ret =  VarSquander(0.3, -0.04165, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;
	


	case 21 : *ret =  VarSquander(0.3, -0.1234, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			  
	case 22 : *ret =  VarSquander2(0.3, -0.1234, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;
	


	case 23 : *ret =  VarSquander(0.42755, 0.11509, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			
	case 24 : *ret =  VarSquander2(0.42755, 0.11509, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;


    
	case 25 : *ret =  VarSquander(0.42755, 0.04078, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;		
	case 26 : *ret =  VarSquander2(0.42755, 0.04078, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			    
	    

	case 27 : *ret =  VarSquander(0.42755,-0.11031, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;		
	case 28 : *ret =  VarSquander(0.42755,-0.11031, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;

			
    
	case 29 : *ret =  VarSquander(0.42755,-0.13091, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			
	case 30 : *ret =  VarSquander2(0.42755,-0.13091, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			    

    
	case 31 : *ret =  VarSquander(0.55792,0.21394, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			
	case 32 : *ret =  VarSquander2(0.55792,0.21394, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			    


	case 33 : *ret =  VarSquander(0.55792,0.0068, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;		
	case 34 : *ret =  VarSquander2(0.55792,0.0068, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			    


	case 35 : *ret =  VarSquander(0.55792,-0.0184, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;		
	case 36 : *ret =  VarSquander2(0.55792,-0.0184, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			    
	

	case 37 : *ret =  VarSquander(0.55792,-0.24335, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;		
	case 38 : *ret =  VarSquander2(0.55792,-0.24335, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			    
	

	case 39 : *ret =  VarSquander(0.68,0.30651, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			    
	case 40 : *ret =  VarSquander2(0.68,0.30651, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;



	case 41 : *ret =  VarSquander(0.68,0.06965, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;		
	case 42 : *ret =  VarSquander2(0.68,0.06965, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			    			    

	    
	case 43 : *ret =  VarSquander(0.68,-0.0172, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			    
	case 44 : *ret =  VarSquander2(0.68,-0.0172, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;



	case 45 : *ret =  VarSquander(0.68,-0.41812, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			    
	case 46 : *ret =  VarSquander2(0.68,-0.41812, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;



	case 47 : *ret =  VarSquander(0.64934,0, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			    
	case 48 : *ret =  VarSquander2(0.64934,0, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;



	case 49 : *ret =  VarSquander(0.6196,0, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			 				    			    
	case 50 : *ret =  VarSquander2(0.6196,0, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;



	case 51 : *ret =  VarSquander(0.58402,0, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			
	case 52 : *ret =  VarSquander2(0.58402,0, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			    

 			    
	case 53 : *ret =  VarSquander(0.25181,0, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			 			    
	case 54 : *ret =  VarSquander2(0.25181,0, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;


	case 55 : *ret =  VarSquander(0.00909,0, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			 			    
	case 56 : *ret =  VarSquander2(0.00909,0, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;


	case 57 : *ret =  VarSquander(-0.93877,0, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			     
	case 58 : *ret =  VarSquander2(-0.93877,0, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;



	case 59 : *ret =  VarSquander(-0.93877,0.20211, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			     
        case 60 : *ret =  VarSquander2(-0.93877,0.20211, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;



	case 61 : *ret =  VarSquander(-0.93877,-0.63517, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			     
	case 62 : *ret =  VarSquander2(-0.93877,-0.63517, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;



	case 63 : *ret =  VarSquander(-1.93877,0, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			     
	case 64 : *ret =  VarSquander2(-1.93877,0, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;



	case 65 : *ret =  VarSquander(-1.93877,0.20211, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			     
        case 66 : *ret =  VarSquander2(-1.93877,0.20211, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;



	case 67 : *ret =  VarSquander(-1.93877,-0.63517, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			     
        case 68 : *ret =  VarSquander2(-1.93877,-0.63517, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;



	case 69 : *ret =  VarSquander(0.42775,0, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			     
        case 70 : *ret =  VarSquander2(0.42775,0, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;



	case 71 : *ret =  VarSquander(0.55792,0, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			     
        case 72 : *ret =  VarSquander2(0.55792,0, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;



	case 73 : *ret =  VarSquander(0,0.07853, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			     
        case 74 : *ret =  VarSquander2(0,0.07853, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;



	case 75 : *ret =  VarSquander(0,0.00339, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			     
        case 76 : *ret =  VarSquander2(0,0.00339, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;

	

	case 77 : *ret =  VarSquander(0,-0.18199, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			     
	case 78 : *ret =  VarSquander2(0,-0.18199, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;



	case 79 : *ret =  VarSquander(0.42755,0.2, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			     
        case 80 : *ret =  VarSquander2(0.42755,0.2, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;



	case 81 : *ret =  VarSquander(0.3,0.36373, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			     
        case 82 : *ret =  VarSquander2(0.3,0.36373, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;



	case 83 : *ret =  VarSquander(0.3,-0.20583, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			     
        case 84 : *ret =  VarSquander2(0.3,-0.20583, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;



	case 85 : *ret =  VarSquander(0.3,-0.40035, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			     
        case 86 : *ret =  VarSquander2(0.3,-0.40035, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;



	case 87 : *ret =  VarSquander(0.3,-0.83259, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			     
        case 88 : *ret =  VarSquander2(0.3,-0.83259, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;



	case 89 : *ret =  VarSquander(0.42755,0.51838, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			     
        case 90 : *ret =  VarSquander2(0.42755,0.51838, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;



	case 91 : *ret =  VarSquander(0.42755,-0.29344, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			     
        case 92 : *ret =  VarSquander2(0.42755,-0.29344, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;



	case 93 : *ret =  VarSquander(0.42755,-0.57056, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			     
        case 94 : *ret =  VarSquander2(0.42755,-0.57056, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;



	case 95 : *ret =  VarSquander(0.42755,-1.18656, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			     
        case 96 : *ret =  VarSquander2(0.42755,-1.18656, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;



	case 97 : *ret =  VarSquander(0.55792,0.67644, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			     
        case 98 : *ret =  VarSquander2(0.55792,0.67644, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;



	case 99 : *ret =  VarSquander(0.55792,-0.38278, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			     
        case 100 : *ret =  VarSquander2(0.55792,-0.38278, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;



	case 101 : *ret =  VarSquander(0.55792,-0.74454, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			     
        case 102 : *ret =  VarSquander2(0.55792,-0.74454, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;



	case 103 : *ret =  VarSquander(0.55792,-1.54837, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;			     
        case 104 : *ret =  VarSquander2(0.55792,-1.54837, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;



	case 105 : *ret =  VarSquander(0.68,0.82445, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;

	case 106 : *ret =  VarSquander2(0.68,0.82445, 
				    x[0],x[1],x[2],x[3],x[4],x[5]); break;  
	

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
	
	case 1 : case 3 : case 5 : case 7 : case 9 : case 11 : case 13 : 
	case 15 : case 17 : case 19 : case 21 : case 23 : case 25 : case 27 : 
	case 29 : case 31 : case 33 : case 35 : case 37 : case 39 : case 41 : 
	case 43 : case 45 : case 47 : case 49 : case 51 : case 53 : case 55 : 
	case 57 : case 59 : case 61 : case 63 : case 65 : case 67 : case 69 : 
	case 71 : case 73 : case 75 : case 77 : case 79 : case 81 : case 83 : 
	case 85 : case 87 : case 89 : case 91 : case 93 : case 95 : case 97 : 
	case 99 : case 101 : case 103 : case 105 : 

	  *ret = rady(x[0],x[1],x[2],x[3],x[4],x[5])-T; break;
	  

case 2 : case 4 : case 6 : case 8 : case 10 : case 12 : case 14 : case 16 : 
	case 18 : case 20 : case 22 : case 24 : case 26 : case 28 : case 30 : 
	case 32 : case 34 : case 36 : case 38 : case 40 : case 42 : case 44 : 
	case 46 : case 48 : case 50 : case 52 : case 54 : case 56 : case 58 : 
	case 60 : case 62 : case 64 : case 66 : case 68 : case 70 : case 72 : 
	case 74 : case 76 : case 78 : case 80 : case 82 : case 84 : case 86 :
	case 88 : case 90 : case 92 : case 94 : case 96 : case 98 : case 100 :
	case 102 : case 104 : case 106 : 

	  *ret = T-rady(x[0],x[1],x[2],x[3],x[4],x[5]); break;
	

	


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
		
		

		  //default : cout << "error " << ineqSwitch << ": not installed " << endl;
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
  minimize2(55);
  minimize2(56);
  minimize2(57);
  minimize2(58);
  minimize2(59);
  minimize2(60);
  minimize2(61);
  minimize2(62);
  minimize2(63);
  minimize2(64);
  minimize2(65);
  minimize2(66);
  minimize2(67);
  minimize2(68);
  minimize2(69);
  minimize2(70);
  minimize2(71);
  minimize2(72);
  minimize2(73);
  minimize2(74);
  minimize2(75);
  minimize2(76);
  minimize2(77);
  minimize2(78);
  minimize2(79);
  minimize2(80);
  minimize2(81);
  minimize2(82);
  minimize2(83);
  minimize2(84);
  minimize2(85);
  minimize2(86);
  minimize2(87);
  minimize2(88);
  minimize2(89);
  minimize2(90);
  minimize2(91);
  minimize2(92);
  minimize2(93);
  minimize2(94);
  minimize2(95);
  minimize2(96);
  minimize2(97);
  minimize2(98);
  minimize2(99);
  minimize2(100);
  minimize2(101);
  minimize2(102);
  minimize2(103);
  minimize2(104);
  minimize2(105);
  minimize2(106);

}








