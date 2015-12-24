#include <iomanip.h>
#include <stdlib.h>
#include "numerical.h"
#include "constants.h"
#include "morefn.h"
#include <math.h>
#include "quoinfn.h"


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


double vol_analytic(double y1,double y2,double y3,double y4,
	double y5,double y6)
	// dodecrad truncated volume of VOronoi
	{
	double dt = 0.7209029495174648;
	return
		(-1./dt)*(vor_analytic(y1,y2,y3,y4,y5,y6)/4.
			-solid(y1,y2,y3,y4,y5,y6)/3.);
	}

// cross diag in terms of edges . Equivalent to version in numerical.cc

double crossdiagNew(double y1,double y2,double y3,double y4,double
		y5,double y6,double y7,double y8,double y9)
	{
	double x1=y1*y1, x2=y2*y2, x3=y3*y3, x4=y4*y4, x5=y5*y5, x6=y6*y6,
			x7=y7*y7, x8=y8*y8, x9=y9*y9;
	double s = dihedraly(y4,y2,y6,y1,y5,y3)+dihedraly(y4,y2,y9,y7,y8,y3);
	if (s>global::pi) { s = 2*global::pi - s; }
	double uc = cos(s)*safesqrt(U(x4,x8,x9)*U(x4,x5,x6));
	double c= -x4*x4+x4*(x6+x9+x8+x5) - x6*x9-x5*x8+x5*x9+x6*x8;
	double a= -2.0*x4;
	return safesqrt((uc-c)/a);
	}

static double quoinH(double y1,double y2,double y6)
	{
	return quoin(y1/2.0,radf(y1,y2,y6),1.255);
	}

static void nofunc(int numargs,int whichFn,double* x,double* ret,void*)
    {
    cout << "nofunc should not be called" << endl << flush;
    }

static double deltay(double y1,double y2,double y3,double y4,double y5,
    double y6)
    {
    return delta(y1*y1,y2*y2,y3*y3,y4*y4,y5*y5,y6*y6);
    }
 

int INEQ_NUMBER=0;
static void generic(int numargs,int whichFn,double* x, double* ret,void*)
	{
	switch (INEQ_NUMBER) {

	/* SAMPLES
	// Section VI.A.4.4.2  (Kepler)
	case 867359387+0:
		*ret= 0.114 
			-gamma(x[0],x[1],x[2],x[3],x[4],x[5])
			-gamma(x[0],x[2],x[7],x[6],x[8],x[4])
			-gamma(x[0],x[7],x[11],x[9],x[10],x[8])
			-gamma(x[0],x[11],x[13],x[12],x[14],x[10])
			-gamma(x[0],x[13],x[1],x[15],x[5],x[14]);
			break;
	*/

	case 222907973: // 1/30/98 test for second generation proof.
		*ret = -(gamma(x[0],x[1],x[2],x[3],x[4],x[5])
			+gamma(x[0],x[2],x[6],x[7],x[8],x[4])
			+gamma(x[0],x[6],x[9],x[10],x[11],x[8]))
		 +(octavorVc(x[0],x[1],x[2],x[3],x[4],x[5])
			+octavorVc(x[0],x[2],x[6],x[7],x[8],x[4])
			+octavorVc(x[0],x[6],x[9],x[10],x[11],x[8]));
		break;

	case 201741996: // 1/30/98 test for second generation proof.
		*ret = -(gamma(x[0],x[1],x[2],x[3],x[4],x[5])
			+gamma(x[0],x[2],x[6],x[7],x[8],x[4])
			+gamma(x[0],x[6],x[9],x[10],x[11],x[8]))
		 +(octavorVc(x[0],x[1],x[2],x[3],x[4],x[5])
			+octavorVc(x[0],x[2],x[6],x[7],x[8],x[4])
			+octavorVc(x[0],x[6],x[9],x[10],x[11],x[8]));
		break;

	case 779007829: // 1/30/98 test for second generation proof.
		*ret = -(gamma(x[0],x[1],x[2],x[3],x[4],x[5])
			+gamma(x[0],x[2],x[6],x[7],x[8],x[4])
			+gamma(x[0],x[6],x[9],x[10],x[11],x[8])
			+gamma(x[0],x[9],x[12],x[13],x[14],x[11]))
			
		 +(octavorVc(x[0],x[1],x[2],x[3],x[4],x[5])
			+octavorVc(x[0],x[2],x[6],x[7],x[8],x[4])
			+octavor(x[0],x[9],x[12],x[13],x[14],x[11])
			+octavorVc(x[0],x[6],x[9],x[10],x[11],x[8]));
		break;



	//Z-ineq
		// start of first page of inequalities for Section 2, SPIV.
		default : cout << "generic default "<<INEQ_NUMBER << endl << flush;
			*ret=0;
		}
	}

iter::~iter() 
	{ delete[] xmin; 
	  delete[] xmax; 
		delete[] x; }


static void ConstraintPageK2
	(int numargs,int whichFn,double* x,double* ret,void*)
    {
	*ret = 0;
	switch (INEQ_NUMBER) {

	/* SAMPLES:
	// Section VI.A.2.8 (Kepler)
	case 314974315+3:
		*ret=dips(x); break;
	case 314974315+1:
		switch(whichFn) {
		case 1 : *ret=dips(x); break;
		case 2 : *ret=global::sqrt2-radf(x[1],x[13],x[15]); break;
		}
		break;
	case 314974315+4:
		switch(whichFn) {
		case 1 : *ret=dips(x); break;
		case 2 : *ret=global::sqrt2-radf(x[5],x[14],x[15]); break;
		}
		break;
	*/

	case 222907973:
		 *ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			+dihedraly(x[0],x[2],x[6],x[7],x[8],x[4])
			+dihedraly(x[0],x[6],x[9],x[10],x[11],x[8])- global::pi;
		break;

	case 201741996:
		 *ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			+dihedraly(x[0],x[2],x[6],x[7],x[8],x[4])
			+dihedraly(x[0],x[6],x[9],x[10],x[11],x[8])
			+dihedraly(x[0],x[1],x[9],x[12],x[11],x[5])
			- global::pi*2.0;
		break;

	case 779007829:
		 *ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			+dihedraly(x[0],x[2],x[6],x[7],x[8],x[4])
			+dihedraly(x[0],x[6],x[9],x[10],x[11],x[8])
			+dihedraly(x[0],x[9],x[12],x[13],x[14],x[11])
			+dihedraly(x[0],x[1],x[12],x[15],x[14],x[5])
			- global::pi*2.0;
		break;

	//Z-con
	default : cout << "unexpected case in constraint " << INEQ_NUMBER<< endl;
		}
    }

iter::iter(int ineqSwitch) {
	numiter = 20; numargs = 6; nconstr=0;
	switch(ineqSwitch)
		{

	/* SAMPLES:
	// Section VI.A.4.12 (Kepler)
	case 836331201+1:
			numargs=9; break;
	*/
	case 222907973: numargs=12; break;
	case 201741996: numargs=13; break;
	case 779007829: numargs=16; break;

	// Z-NUMARGS
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

	/* SAMPLES:
	// Section VI.A.2.5. (Kepler)
	case 352079526+1:
		xmin[3]=2.51; xmax[3]=global::sqrt8;
		nconstr=2;
		break;
	case 352079526+3:
		xmin[3]=2.7; xmax[3]=global::sqrt8;
		nconstr=1;
		break;
	*/
	case 222907973:
		xmin[0]=2.51; xmax[0]=2.67;
		nconstr=1;
		break;

	case 201741996:
		xmin[0]=2.51; xmax[0]=2.67;
		xmin[12]=xmax[12]=3.2;
		nconstr=1;
		break;

	case 779007829:
		xmin[0]=2.51; xmax[0]=2.67;
		xmin[15]=xmax[15]=3.2;
		nconstr=1;
		break;

		//Z-vars
		default : cout << "error " << ineqSwitch << ": not installed " << endl;
		}

	if (nconstr>0) constraintfunc=ConstraintPageK2;
	}

void /*ineq.cc*/minimize2(int);

void page1 ()
	{
	/*
	cout << "------ Section VI.A.2.5 of Kepler. --------" << endl;
	minimize2(293389419+1);


	// BREAK
	*/
	//minimize2(222907973);
	//minimize2(201741996);
	minimize2(779007829);
	return;
	}
