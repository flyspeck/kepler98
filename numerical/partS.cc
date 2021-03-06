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


double SeanVol(double y1,double y2,double y3,double y4,
	double y5,double y6)
	// dodecrad truncated volume of VOronoi
	{
	double dr = global::dodecrad/2.;
	double dt = 0.7209029495174648;
	return
		(-1./dt)*(vorVc(y1,y2,y3,y4,y5,y6,dr)/4.
			-solid(y1,y2,y3,y4,y5,y6)/3.);
	}

double SeanSquander(double y1,double y2,double y3,double y4,
	double y5,double y6)
	{
	return SeanVol(y1,y2,y3,y4,y5,y6)
		- 0.42755*solid(y1,y2,y3,y4,y5,y6);
	}

double dips(double x[16])
	{
	return 
		dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			+dihedraly(x[0],x[2],x[7],x[6],x[8],x[4])
			+dihedraly(x[0],x[7],x[11],x[9],x[10],x[8])
			+dihedraly(x[0],x[11],x[13],x[12],x[14],x[10])
			+dihedraly(x[0],x[13],x[1],x[15],x[5],x[14])
			-2.0*global::pi;
	}

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

	// Section A.2.5.

	case 12920:
		*ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]);
		break;
	case 12921:
		*ret = SeanSquander(x[0],x[1],x[2],x[3],x[4],x[5]);
		break;
	case 12638:
		*ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]);
		break;
	case 12786: case 13035:
		*ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]);
		break;
	case 402251:
		*ret = x[0];
		break;
	case 735559098:
		{
		int NumSphere=24;
		int i;
		double t =0;
		for (i=0;i<4*(NumSphere+1);i++) t+= x[i]*x[i];
		*ret = t;
		}
		break;
		




	//Z-ineq
		// start of first page of inequalities for Section 2, SPIV.
		default : cout << "generic default"<<INEQ_NUMBER << endl << flush;
			*ret=0;
		}
	}

iter::~iter() 
	{ delete[] xmin; 
	  delete[] xmax; 
		delete[] x; }

double ovol(double * x)
	{
	double q1[3]= {2,0,0};
	double q2[3]= {0,2,0};
	double q3[3]= {0,0,2};
	return 33.51032163829112/(q1[0]*(q2[1]+x[2])*(q3[2]+x[5]));
	}

double dist(double* r1,double* r2)
	{
	double t=0;

	for (int i=0;i<3;i++) t += (r1[i]-r2[i])*(r1[i]-r2[i]);
	return sqrt(t);
	}

double radd(double* r1,double* r2,double* r3,double* r4)
	{
	return rady(dist(r1,r2),dist(r1,r3),dist(r1,r4),
				dist(r3,r4),dist(r2,r4),dist(r2,r3));
	}

double radix(int which,double* x) 
	//	(* circumradius in phelan-weaire, which=1..46 *)
	{
	double q0[3]= {0,0,0};
	double q1[3]= {2,0,0};
	double q2[3]= {0+x[1],2+x[2],0};
	double q3[3]= {0+x[3],0+x[4],2+x[5]};
	double r0[3]= {1+x[6],1+x[7],1+x[8]};
	double p1[3]= {1+x[9],0.5+x[10],0+x[11]};
	double p2[3]= {1+x[12],1.5+x[13],0+x[14]};
	double p3[3]= {0+x[15],1+x[16],0.5+x[17]};
	double p4[3]= {0+x[18],1+x[19],1.5+x[20]};
	double p5[3]= {1.5+x[21],0+x[22],1+x[23]};
	double p6[3]= {0.5+x[24],0+x[25],1+x[26]};
	double r1[3],r2[3],r3[3],r4[3];
	int i;
	switch (which) {

		case 1:
			for (i=0;i<3;i++)
				{
				r1[i]= p1[i];
				r2[i]= p2[i];
				r3[i]= p3[i]+q1[i];
				r4[i]= p4[i]+q1[i]-q3[i];
				}
			break;

		case 2:
			for (i=0;i<3;i++)
				{
				r1[i]= p1[i];
				r2[i]= p2[i];
				r3[i]= p3[i];
				r4[i]= p4[i]-q3[i];
				}
			break;

		case 3:
			for (i=0;i<3;i++)
				{
				r1[i]= p5[i];
				r2[i]= p6[i];
				r3[i]= p1[i];
				r4[i]= p2[i]-q2[i];
				}
			break;

		case 4:
			for (i=0;i<3;i++)
				{
				r1[i]= p5[i];
				r2[i]= p6[i];
				r3[i]= p1[i]+q3[i];
				r4[i]= p2[i]+q3[i]-q2[i];
				}
			break;

		case 5:
			for (i=0;i<3;i++)
				{
				r1[i]= p3[i];
				r2[i]= p4[i];
				r3[i]= p6[i];
				r4[i]= p5[i]-q1[i];
				}
			break;

		case 6:
			for (i=0;i<3;i++)
				{
				r1[i]= p3[i];
				r2[i]= p4[i];
				r3[i]= p6[i]+q2[i];
				r4[i]= p5[i]+q2[i]-q1[i];
				}
			break;

		case 7:
			for (i=0;i<3;i++)
				{
				r1[i]= r0[i];
				r2[i]= p5[i]+q2[i];
				r3[i]= p6[i]+q2[i];
				r4[i]= p2[i];
				}
			break;

		case 8:
			for (i=0;i<3;i++)
				{
				r1[i]= r0[i];
				r2[i]= p5[i]+q2[i];
				r3[i]= p2[i];
				r4[i]= p3[i]+q1[i];
				}
			break;

		case 9:
			for (i=0;i<3;i++)
				{
				r1[i]= r0[i];
				r2[i]= p5[i]+q2[i];
				r3[i]= p3[i]+q1[i];
				r4[i]= p4[i]+q1[i];
				}
			break;

		case 10:
			for (i=0;i<3;i++)
				{
				r1[i]= r0[i];
				r2[i]= p5[i]+q2[i];
				r3[i]= p4[i]+q1[i];
				r4[i]= p2[i]+q3[i];
				}
			break;

		case 11:
			for (i=0;i<3;i++)
				{
				r1[i]= r0[i];
				r2[i]= p5[i]+q2[i];
				r3[i]= p6[i]+q2[i];
				r4[i]= p2[i]+q3[i];
				}
			break;

		case 12:
			for (i=0;i<3;i++)
				{
				r1[i]= r0[i];
				r2[i]= p6[i]+q2[i];
				r3[i]= p2[i]+q3[i];
				r4[i]= p4[i];
				}
			break;

		case 13:
			for (i=0;i<3;i++)
				{
				r1[i]= r0[i];
				r2[i]= p6[i]+q2[i];
				r3[i]= p4[i];
				r4[i]= p3[i];
				}
			break;

		case 14:
			for (i=0;i<3;i++)
				{
				r1[i]= r0[i];
				r2[i]= p6[i]+q2[i];
				r3[i]= p3[i];
				r4[i]= p2[i];
				}
			break;

		case 15:
			for (i=0;i<3;i++)
				{
				r1[i]= r0[i];
				r2[i]= p5[i];
				r3[i]= p6[i];
				r4[i]= p1[i];
				}
			break;

		case 16:
			for (i=0;i<3;i++)
				{
				r1[i]= r0[i];
				r2[i]= p5[i];
				r3[i]= p1[i];
				r4[i]= p3[i]+q1[i];
				}
			break;

		case 17:
			for (i=0;i<3;i++)
				{
				r1[i]= r0[i];
				r2[i]= p5[i];
				r3[i]= p3[i]+q1[i];
				r4[i]= p4[i]+q1[i];
				}
			break;

		case 18:
			for (i=0;i<3;i++)
				{
				r1[i]= r0[i];
				r2[i]= p5[i];
				r3[i]= p4[i]+q1[i];
				r4[i]= p1[i]+q3[i];
				}
			break;

		case 19:
			for (i=0;i<3;i++)
				{
				r1[i]= r0[i];
				r2[i]= p5[i];
				r3[i]= p6[i];
				r4[i]= p1[i]+q3[i];
				}
			break;

		case 20:
			for (i=0;i<3;i++)
				{
				r1[i]= r0[i];
				r2[i]= p6[i];
				r3[i]= p1[i]+q3[i];
				r4[i]= p4[i];
				}
			break;

		case 21:
			for (i=0;i<3;i++)
				{
				r1[i]= r0[i];
				r2[i]= p6[i];
				r3[i]= p4[i];
				r4[i]= p3[i];
				}
			break;

		case 22:
			for (i=0;i<3;i++)
				{
				r1[i]= r0[i];
				r2[i]= p6[i];
				r3[i]= p3[i];
				r4[i]= p1[i];
				}
			break;

		case 23:
			for (i=0;i<3;i++)
				{
				r1[i]= r0[i];
				r2[i]= p3[i];
				r3[i]= p1[i];
				r4[i]= p2[i];
				}
			break;

		case 24:
			for (i=0;i<3;i++)
				{
				r1[i]= r0[i];
				r2[i]= p1[i];
				r3[i]= p2[i];
				r4[i]= p3[i]+q1[i];
				}
			break;

		case 25:
			for (i=0;i<3;i++)
				{
				r1[i]= r0[i];
				r2[i]= p1[i]+q3[i];
				r3[i]= p2[i]+q3[i];
				r4[i]= p4[i];
				}
			break;

		case 26:
			for (i=0;i<3;i++)
				{
				r1[i]= r0[i];
				r2[i]= p1[i]+q3[i];
				r3[i]= p2[i]+q3[i];
				r4[i]= p4[i]+q1[i];
				}
			break;

		case 27:  
			for (i=0;i<3;i++)
				{
				r1[i]= q0[i];
				r2[i]= p4[i]-q3[i];
				r3[i]= p3[i];
				r4[i]= p1[i];
				}
			break;

		case 28:
			for (i=0;i<3;i++)
				{
				r1[i]= q0[i];
				r2[i]= p3[i];
				r3[i]= p1[i];
				r4[i]= p6[i];
				}
			break;

		case 29:
			for (i=0;i<3;i++)
				{
				r1[i]= q0[i];
				r2[i]= p3[i];
				r3[i]= p6[i];
				r4[i]= p5[i]-q1[i];
				}
			break;

		case 30:
			for (i=0;i<3;i++)
				{
				r1[i]= q0[i];
				r2[i]= p3[i];
				r3[i]= p5[i]-q1[i];
				r4[i]= p1[i]-q1[i];
				}
			break;

		case 31:
			for (i=0;i<3;i++)
				{
				r1[i]= q0[i];
				r2[i]= p3[i];
				r3[i]= p4[i]-q3[i];
				r4[i]= p1[i]-q1[i];
				}
			break;

		case 32:
			for (i=0;i<3;i++)
				{
				r1[i]= q0[i];
				r2[i]= p4[i]-q3[i];
				r3[i]= p1[i]-q1[i];
				r4[i]= p5[i]-q1[i]-q3[i];
				}
			break;

		case 33:
			for (i=0;i<3;i++)
				{
				r1[i]= q0[i];
				r2[i]= p4[i]-q3[i];
				r3[i]= p5[i]-q1[i]-q3[i];
				r4[i]= p6[i]-q3[i];
				}
			break;

		case 34:
			for (i=0;i<3;i++)
				{
				r1[i]= q0[i];
				r2[i]= p4[i]-q3[i];
				r3[i]= p6[i]-q3[i];
				r4[i]= p1[i];
				}
			break;

		case 35:
			for (i=0;i<3;i++)
				{
				r1[i]= q0[i];
				r2[i]= p4[i]-q3[i]-q2[i];
				r3[i]= p3[i]-q2[i];
				r4[i]= p2[i]-q2[i];
				}
			break;

		case 36:
			for (i=0;i<3;i++)
				{
				r1[i]= q0[i];
				r2[i]= p3[i]-q2[i];
				r3[i]= p2[i]-q2[i];
				r4[i]= p6[i];
				}
			break;

		case 37:
			for (i=0;i<3;i++)
				{
				r1[i]= q0[i];
				r2[i]= p3[i]-q2[i];
				r3[i]= p6[i];
				r4[i]= p5[i]-q1[i];
				}
			break;

		case 38:
			for (i=0;i<3;i++)
				{
				r1[i]= q0[i];
				r2[i]= p3[i]-q2[i];
				r3[i]= p5[i]-q1[i];
				r4[i]= p2[i]-q2[i]-q1[i];
				}
			break;

		case 39:
			for (i=0;i<3;i++)
				{
				r1[i]= q0[i];
				r2[i]= p3[i]-q2[i];
				r3[i]= p4[i]-q3[i]-q2[i];
				r4[i]= p2[i]-q2[i]-q1[i];
				}
			break;

		case 40:
			for (i=0;i<3;i++)
				{
				r1[i]= q0[i];
				r2[i]= p4[i]-q3[i]-q2[i];
				r3[i]= p2[i]-q2[i]-q1[i];
				r4[i]= p5[i]-q3[i]-q1[i];
				}
			break;

		case 41:
			for (i=0;i<3;i++)
				{
				r1[i]= q0[i];
				r2[i]= p4[i]-q3[i]-q2[i];
				r3[i]= p5[i]-q3[i]-q1[i];
				r4[i]= p6[i]-q3[i];
				}
			break;

		case 42:
			for (i=0;i<3;i++)
				{
				r1[i]= q0[i];
				r2[i]= p4[i]-q3[i]-q2[i];
				r3[i]= p6[i]-q3[i];
				r4[i]= p2[i]-q2[i];
				}
			break;

		case 43:
			for (i=0;i<3;i++)
				{
				r1[i]= q0[i];
				r2[i]= p1[i];
				r3[i]= p2[i]-q2[i];
				r4[i]= p6[i]-q3[i];
				}
			break;

		case 44:
			for (i=0;i<3;i++)
				{
				r1[i]= q0[i];
				r2[i]= p1[i];
				r3[i]= p2[i]-q2[i];
				r4[i]= p6[i];
				}
			break;

		case 45:
			for (i=0;i<3;i++)
				{
				r1[i]= q0[i];
				r2[i]= p1[i]-q1[i];
				r3[i]= p2[i]-q1[i]-q2[i];
				r4[i]= p5[i]-q3[i]-q1[i];
				}
			break;

		case 46:
		default:
			for (i=0;i<3;i++)
				{
				r1[i]= q0[i];
				r2[i]= p1[i]-q1[i];
				r3[i]= p2[i]-q1[i]-q2[i];
				r4[i]= p5[i]-q1[i];
				}
			break;
		}


	return radd(r1,r2,r3,r4);
	};

double negsquares(int i,int j,double* x)
	{
	int k;
	double t =0;
	for (k=0;k<4;k++) 
		t += (x[4*i+k]-x[4*j+k])*(x[4*i+k]-x[4*j+k]);
	return 4-t;
	}

static void ConstraintPageK2
	(int numargs,int whichFn,double* x,double* ret,void*)
    {
	*ret = 0;
	switch (INEQ_NUMBER) {
	case 12921:
		*ret=dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.53; 
		break;
	case 402251:
		{
		double r = radix(whichFn,x);
		double v = ovol(x);
		*ret= -(x[0]- v*r*r*r);
		}
		break;
	case 735559098:
		{
		int NumSphere=24;
		int i=0;
		int j = whichFn;
		int start = NumSphere;
		while ((j>start)&&(start>0)) { j -= start; i++; start--; }
		*ret = negsquares(i,i+j,x);
		}
		break;

	//Z-con
	default : cout << "unexpected case in constraint " << INEQ_NUMBER<< endl;
		}
    }

iter::iter(int ineqSwitch) {
	numiter = 10; numargs = 6; nconstr=0;
	int NumSphere=24;
	switch(ineqSwitch)
		{
		case 402251:  numargs = 27; break;
		case 735559098: numargs = 4*(NumSphere+1); break;

	// Z-NUMARGS
		}

	// temp: 
	xmin = new double[numargs];
	xmax = new double[numargs];
	x = new double[numargs];
	constraintfunc = nofunc;
	func = generic;
	int i;
	for (i=0;i<numargs;i++) { xmin[i]=x[i]=2.0; xmax[i]=2.51681714472963793; }
	INEQ_NUMBER = ineqSwitch;
	switch (ineqSwitch)
		{

	case 12920: 
		xmax[0]= 2.072;
		xmax[1]= 2.11;
		xmax[2]= 2.11;
		xmin[3]= 2.82;
		xmax[3]= 2.82;
		xmax[4]= 2.37;
		xmax[5]= 2.37;
		break;

	case 12921: 
		xmax[0]= 2.072;
		xmax[1]= 2.11;
		xmax[2]= 2.11;
		xmin[3]= 2.51;
		xmax[3]= 2.82;
		xmax[4]= 2.37;
		xmax[5]= 2.37;
		nconstr = 1;
		break;

	case 12638:
		xmax[0]=2.03;
		xmax[1]=2.03;
		xmax[2]=2.03;
		xmin[3]=2.51;
		xmax[4]=2.13;
		xmax[5]=2.13;
		break;

	case 12786:
		xmax[0]=2.03;
		xmax[1]=2.01;
		xmax[2]=2.03;
		xmin[3]=2.51;
		xmax[4]=2.181;
		xmax[5]=2.21;
		break;
	case 13035:
		xmax[0]=xmax[1]=2.01;
		xmax[2]=2.03;
		xmin[3]=2.51;
		xmax[4]=2.2;
		xmax[5]=2.01;	
		break;
	case 402251:
		{
		int i;
		for (i=0;i<27;i++) { xmin[i]= -0.5; xmax[i]=0.5; }
		}
		xmin[0]=-10; xmax[0]= 100.0;
		nconstr=46; 
		break;
	case 735559098:
		{
		int NumSphere=24;
		int i;
		for (i=0;i<4*(NumSphere+1);i++) { xmin[i]= -10; xmax[i]=10; }
		for (i=0;i<4;i++) { xmin[i]=xmax[i]=0; }
		nconstr=(NumSphere*(NumSphere+1))/2;
		}
		break;



	case 480930831:
		xmin[3]=2.51; xmax[3]=global::sqrt8;
		xmin[4]=2.77; xmax[4]=global::sqrt8;
		break;
	case 463544803:
		xmin[3]=2.7; xmax[3]=global::sqrt8;
		break;
	case 594246986:
	case 381970727:
		xmin[3]=2.51; xmax[3]=global::sqrt8;
		xmax[4]=xmax[5]=2.3;
		nconstr=1;
		break;
	case 951798877:
	case 923397705:
		xmin[3]=2.51; xmax[3]=global::sqrt8;
		nconstr=1;
	case 495568072:
		xmin[3]=2.51; xmax[3]=2.7;
		nconstr=1;
		xmax[0]=xmax[1]=xmax[2]=2.0;
		break;
	case 378020227:
		xmin[4]=2.51; xmax[4]=2.77; 
		xmin[5]=2.51; xmax[5]=2.77; break;
	case 256893386:
		xmin[4]=2.51; xmax[4]=2.77;
		xmin[5]=2.77; xmax[5]=global::sqrt8;
		break;
	case 749955642:
		xmax[0]=xmax[1]=xmax[2]=2.;
		nconstr=1;
		xmin[4]=2.51; xmax[4]=2.77; 
		xmin[5]=2.51; xmax[5]=2.77; break;
	case 653849975:
		nconstr=1;
		xmin[4]=2.51; xmax[4]=2.77; 
		xmin[5]=2.51; xmax[5]=2.77; break;

		//Z-vars
		default : cout << "error " << ineqSwitch << ": not installed " << endl;
		}

	if (nconstr>0) constraintfunc=ConstraintPageK2;
	}

void /*ineq.cc*/minimize2(int);

void testRadix()
	{
	int i;
	double x[27]={1.25444752074805277,-0.25,0.25,0,
    -0.113397449804600231,
-0.11610429144945987,0,-0.181698671918555571,
0.0669478250568131816,0,-0.0533634667495881121,
0.0605299985482918926,0,-0.196636493819324437,
0.189469966213560326,0,-0.136204267709799598,
0.0811489434338724669,0,-0.227193132524489111,
0.0527467222257731461,-0.00694156613224826413,-0.0566986908117247598,
-0.0580521534422987623,0.00694153855748141056,-0.0566986907019925435,
-0.0580521587327554114};

	// for (i=0;i<27;i++) x[i]=0.0;

	if (1==1) {
	for (i=1;i<47;i++) cout << radix(i,x) << endl;
	for (i=1;i<47;i++)
		{
		double r=0;
		void* a;
		INEQ_NUMBER=402251;
		ConstraintPageK2(27,i,x,&r,a);
		cout << r << endl;
		}
	cout << endl << ovol(x) << endl;
	}

	
	}

void page1 ()
	{
	/*
	cout << "------ Section VI.A.2.5 of Kepler. --------" << endl;
	minimize2(269048407+0 );
	minimize2(269048407+1 );
	minimize2(553285469+0);
	minimize2(553285469+1);
	minimize2(293389419+0);
	minimize2(293389419+1);
	*/

	/*
	// Appendices, listed by archive number
	minimize2(12920);
	minimize2(12920+1);
	minimize2(12638);
	minimize2(12786);
	minimize2(13035);
	minimize2(402251);
	//testRadix();
	*/


	minimize2(735559098);


	// BREAK
	return;
	}
