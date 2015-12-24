#include <iomanip.h>
#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>
#include "numerical.h"
#include "constants.h"
#include "morefn.h"
#include "quoinfn.h"
#include <math.h>
#include "temp355.cc" 


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

iter::~iter() 
	{ delete[] xmin; 
	  delete[] xmax; 
		delete[] x; }

static double vorAnchor3(double y1,double y2,double y6) // was 2
    {
	double x = y2+y6;
     return (0.0002  
	-0.2695279326151798 + 0.1833230667013778*y1 - 0.02783887375001181*y1*y1
	-0.0152253
	+0.7557773828234548 - 0.3239044460786886*x + 0.0346916048615081*x*x);
	//(- 0.027*(y2-2.)-0.0264*(y6-2.));
    }

static double vorAnchor2(double y1,double y2,double y6) 
    {
	double x = y2+y6;
     return (0.00011
	-0.2695279326151798 + 0.1833230667013778*y1 - 0.02783887375001181*y1*y1
	-0.0152253
	+0.4384851103526791 - 0.167484482694134*x + 0.01541738104479281*x*x);
    }

static double dih2R(double y1,double y2,double y6) // used in A1...
	{
	return dihR(y2/2.,radf(y1,y2,y6),y1/(2.*cos(arc(y1,1.255,1.6))));
	}

static void nofunc(int numargs,int whichFn,double* x,double* ret,void*)
    {
    cout << "nofunc should not be called" << endl << flush;
    }

double dips(double x[16])
    {
    return
        dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
            +dihedraly(x[0],x[2],x[7],x[6],x[3],x[4])
            +dihedraly(x[0],x[7],x[11],x[9],x[10],x[8])
            +dihedraly(x[0],x[11],x[13],x[12],x[14],x[10])
            +dihedraly(x[0],x[13],x[1],x[15],x[5],x[14])
            -2.0*global::pi;
    }

double gammaAX(double y1,double y2,double y3,double y4,double y5,double y6)
	{
	return gamma(y1,y2,y3,y4,y5,y6)+0.419351*solid(y1,y2,y3,y4,y5,y6);
	}

double dihConstraint(double x[6],double dihmax)
	{
	double x1=x[0]*x[0], x2=x[1]*x[1], x3=x[2]*x[2],
		x4=x[3]*x[3], x5=x[4]*x[4], x6=x[5]*x[5];
	double d4 = -(x2*x3) - x1*x4 + x2*x5 + x3*x6 - x5*x6 + 
			x1*(-x1 + x2 + x3 - x4 + x5 + x6);
	double t = tan(dihmax-global::pi/2.0);
	return (d4*d4- delta(x1,x2,x3,x4,x5,x6)*
					4.*x1*t*t);
	}

static double deltay(double y1,double y2,double y3,double y4,double y5,
    double y6)
    {
    return delta(y1*y1,y2*y2,y3*y3,y4*y4,y5*y5,y6*y6);
    }


double XMIN[6],XMAX[6],
	   TAUC,XC,y2C,y3C,y5C,y6C,xDC;

int INEQ_NUMBER=0;
static void generic(int numargs,int whichFn,double* x, double* ret,void*)
	{
	switch (INEQ_NUMBER) {

	case 0: 
		*ret = -(TAUC*tau(x[0],x[1],x[2],x[3],x[4],x[5])+
				XC+y2C*x[1]+y3C*x[2]+y5C*x[4]+y6C*x[5]+
				xDC*(dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-2*global::pi/5.)
				);
		break;

	//Z-ineq
		// start of first page of inequalities for Section 2, SPIV.
		default : cout << "generic default "<< INEQ_NUMBER << endl << flush;
			*ret=0;
		}
	}


static void ConstraintPage1(int numargs,int whichFn,double* x,double* ret,void*)
    {
	*ret = 0;
	switch (INEQ_NUMBER) {



		//Z-con
		default : cout << "unexpected case in constraint " << INEQ_NUMBER<< endl;
		}
    }

iter::iter(int ineqSwitch) {
	numiter = 80; numargs = 6; nconstr=0; // numiter was 20;
	switch(ineqSwitch)
		{
		case 348940660+1 : numargs=16;
		}
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
		case 0: 
		
		for (i=0;i<6;i++) 
			{ 
			xmin[i]=XMIN[i];
			xmax[i]=XMAX[i];
			}
		if (ineqSwitch==1) { xmin[3]=xmax[3]=global::sqrt8; }
		break;

			
		//Z-vars
		default : cout << "error " << ineqSwitch << ": not installed " << endl;
		}

	if (nconstr>0) constraintfunc=ConstraintPage1;
	}

double /*ineq.cc*/minimize2(int);

void page1()
	{
	ifstream inFile;
	inFile.open("/tmp/subdiv_cfsqp.dat");
	XMIN[0]=XMIN[3]=0.; XMAX[0]=XMAX[3]=1;
	int j;
	double amin = 10.0;
	double slack;
	for (j=0;j<5;j++)
		{
		double x[6];
		// order is essential here , see Subdivide.m:varSetTauA, etc. 
		inFile >> XMIN[0] >> XMIN[1] >> XMIN[2] >> 
				XMIN[3] >> XMIN[4] >> XMIN[5];
		inFile >> XMAX[0] >> XMAX[1] >> XMAX[2] >> 
				XMAX[3] >> XMAX[4] >> XMAX[5];
		inFile >> TAUC >> slack >> XC >> y2C >> y3C >> y5C >> y6C >> xDC;
		if (j<4) { XMAX[3]=XMAX[0]=2.51; }
		if (j==4){ XMAX[0]=2.51; }
		x[j]=minimize2(0);
		if (x[j]<amin) amin = x[j];
		cout << "slack = " << slack << endl;
		}
	ofstream oFile;
	oFile.open("/tmp/subdiv_cfsqp.m");
	oFile << "cfsqpVal = " << amin << endl;
	}
