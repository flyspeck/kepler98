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

	case 269048407+0 :
		*ret = vorVc(x[0],x[1],x[2],x[3],x[4],x[5]) +
			0.01*(global::pi/2. - dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]))
			- gammaNu(x[0],x[1],x[2],x[3],x[4],x[5]); break;
	case 269048407+1 :
		*ret = vorVc(x[0],x[1],x[2],x[3],x[4],x[5]) +
			0.01*(global::pi/2. - dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]))
			- vorNu(x[0],x[1],x[2],x[3],x[4],x[5]); break;
	case 553285469+0:
		*ret = vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
			- gammaNu(x[0],x[1],x[2],x[3],x[4],x[5]); break;
	case 553285469+1:
		*ret = vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
			- vorNu(x[0],x[1],x[2],x[3],x[4],x[5]); break;
	case 293389419+0:
		*ret = vorVc(x[0],x[1],x[2],x[3],x[4],x[5])+0.0268
            - gamma(x[0],x[1],x[2],x[3],x[4],x[5]); break;
	case 293389419+1:
		*ret = vorVc(x[0],x[1],x[2],x[3],x[4],x[5])+0.0268
            - vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5]); break;
	case 695069283+0:
		*ret = vorVc(x[0],x[1],x[2],x[3],x[4],x[5])+0.02
            - gamma(x[0],x[1],x[2],x[3],x[4],x[5]); break;
	case 695069283+1:
		*ret = vorVc(x[0],x[1],x[2],x[3],x[4],x[5])+0.02
            - vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5]); break;
	case 814398901 :
		*ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.32; break;
	case 352079526+0:
		*ret = tau(x[0],x[1],x[2],x[3],x[4],x[5]) -3.07*global::pt; break;
	case 352079526+1:
		*ret = tau_analytic(x[0],x[1],x[2],x[3],x[4],x[5]) -3.07*global::pt; 
		break;
	case 352079526+2:
	case 352079526+3:
	case 352079526+4:
		*ret = tauVc(x[0],x[1],x[2],x[3],x[4],x[5]) -3.07*global::pt; break;
	case 179025673:
		*ret = tauVc(x[0],x[1],x[2],x[3],x[4],x[5]) -3.07*global::pt
			- 0.003521 - 2.*0.00935; break;




	// Section VI.A.2.7 
	case 551665569:
		*ret= -1.4*global::pt
			+tau(x[0],x[1],x[2],x[3],x[4],x[5])
			+tau(x[0],x[2],x[7],x[6],x[8],x[4])
			+0*tau(x[0],x[7],x[11],x[9],x[10],x[8])
			+tau(x[0],x[11],x[13],x[12],x[14],x[10])
			+0*tau(x[0],x[13],x[1],x[15],x[5],x[14]);
			break;
	case 824762926:
		*ret= -1.4*global::pt+tau(x[0],x[1],x[2],x[3],x[4],x[5])
			+tau(x[0],x[2],x[7],x[6],x[8],x[4])
			+tau(x[0],x[7],x[11],x[9],x[10],x[8])
			+0*tau(x[0],x[11],x[13],x[12],x[14],x[10])
			+0*tau(x[0],x[13],x[1],x[15],x[5],x[14]);
			break;
	// 0,1,2,3,4, gamma,vor,vor0,vor0,vor0 as in VI.2.5 listing.
	case 675785884+0:
		*ret= -1.4*global::pt-0.06585
			+tau(x[0],x[1],x[2],x[3],x[4],x[5])
			+tau(x[0],x[2],x[7],x[6],x[8],x[4])
			+tau(x[0],x[7],x[11],x[9],x[10],x[8])
			+tau(x[0],x[11],x[13],x[12],x[14],x[10])
			+0*tau(x[0],x[13],x[1],x[15],x[5],x[14]);
			break;
	case 675785884+1:
		*ret= -1.4*global::pt-0.06585
			+tau(x[0],x[1],x[2],x[3],x[4],x[5])
			+tau(x[0],x[2],x[7],x[6],x[8],x[4])
			+tau_analytic(x[0],x[7],x[11],x[9],x[10],x[8])
			+tau(x[0],x[11],x[13],x[12],x[14],x[10])
			+0*tau(x[0],x[13],x[1],x[15],x[5],x[14]);
			break;
	case 675785884+2:
		*ret= -1.4*global::pt-0.06585
			+tau(x[0],x[1],x[2],x[3],x[4],x[5])
			+tau(x[0],x[2],x[7],x[6],x[8],x[4])
			+tauVc(x[0],x[7],x[11],x[9],x[10],x[8])-0.01561-2*0.003521
			+tau(x[0],x[11],x[13],x[12],x[14],x[10])
			+0*tau(x[0],x[13],x[1],x[15],x[5],x[14]);
			break;
	case 675785884+3:
		*ret= -1.4*global::pt-0.06585
			+tau(x[0],x[1],x[2],x[3],x[4],x[5])
			+tau(x[0],x[2],x[7],x[6],x[8],x[4])
			+tauVc(x[0],x[7],x[11],x[9],x[10],x[8])-0.01561-2*0.003521
			+tau(x[0],x[11],x[13],x[12],x[14],x[10])
			+0*tau(x[0],x[13],x[1],x[15],x[5],x[14]);
			break;
	case 675785884+4:
		*ret= -1.4*global::pt-0.06585
			+tau(x[0],x[1],x[2],x[3],x[4],x[5])
			+tau(x[0],x[2],x[7],x[6],x[8],x[4])
			+tauVc(x[0],x[7],x[11],x[9],x[10],x[8])
			+tau(x[0],x[11],x[13],x[12],x[14],x[10])
			+0*tau(x[0],x[13],x[1],x[15],x[5],x[14]);
			break;

	case 193592217+0:
		*ret= -1.4*global::pt-0.06585
			+tau(x[0],x[1],x[2],x[3],x[4],x[5])
			+tau(x[0],x[2],x[7],x[6],x[8],x[4])
			+tau(x[0],x[7],x[11],x[9],x[10],x[8])
			+tau(x[0],x[11],x[13],x[12],x[14],x[10])
			+0*tau(x[0],x[13],x[1],x[15],x[5],x[14]);
			break;

	case 193592217+1:
		*ret= -1.4*global::pt-0.06585
			+tau(x[0],x[1],x[2],x[3],x[4],x[5])
			+tau(x[0],x[2],x[7],x[6],x[8],x[4])
			+tau(x[0],x[7],x[11],x[9],x[10],x[8])
			+tau_analytic(x[0],x[11],x[13],x[12],x[14],x[10])
			+0*tau(x[0],x[13],x[1],x[15],x[5],x[14]);
			break;
	case 193592217+2:
		*ret= -1.4*global::pt-0.06585
			+tau(x[0],x[1],x[2],x[3],x[4],x[5])
			+tau(x[0],x[2],x[7],x[6],x[8],x[4])
			+tau(x[0],x[7],x[11],x[9],x[10],x[8])
			+tauVc(x[0],x[11],x[13],x[12],x[14],x[10])-2.*0.003521-0.01561
			+0*tau(x[0],x[13],x[1],x[15],x[5],x[14]);
			break;
	case 193592217+3:
		*ret= -1.4*global::pt-0.06585
			+tau(x[0],x[1],x[2],x[3],x[4],x[5])
			+tau(x[0],x[2],x[7],x[6],x[8],x[4])
			+tau(x[0],x[7],x[11],x[9],x[10],x[8])
			+tauVc(x[0],x[11],x[13],x[12],x[14],x[10])-2.*0.003521-0.01561
			+0*tau(x[0],x[13],x[1],x[15],x[5],x[14]);
			break;
	case 193592217+4:
		*ret= -1.4*global::pt-0.06585
			+tau(x[0],x[1],x[2],x[3],x[4],x[5])
			+tau(x[0],x[2],x[7],x[6],x[8],x[4])
			+tau(x[0],x[7],x[11],x[9],x[10],x[8])
			+tauVc(x[0],x[11],x[13],x[12],x[14],x[10])
			+0*tau(x[0],x[13],x[1],x[15],x[5],x[14]);
			break;

	// Section VI.A.2.8 
	case 325738864:
		*ret= -1.5*global::pt+tau(x[0],x[1],x[2],x[3],x[4],x[5])
			+tau(x[0],x[2],x[7],x[6],x[8],x[4])
			+tau(x[0],x[7],x[11],x[9],x[10],x[8])
			+tau(x[0],x[11],x[13],x[12],x[14],x[10])
			+0*tau(x[0],x[13],x[1],x[15],x[5],x[14]);
			break;
	case 314974315+0:
		*ret= -1.5*global::pt-0.06585
			+tau(x[0],x[1],x[2],x[3],x[4],x[5])
			+tau(x[0],x[2],x[7],x[6],x[8],x[4])
			+tau(x[0],x[7],x[11],x[9],x[10],x[8])
			+tau(x[0],x[11],x[13],x[12],x[14],x[10])
			+tau(x[0],x[13],x[1],x[15],x[5],x[14]);
			break;
	case 314974315+1:
		*ret= -1.5*global::pt-0.06585
			+tau(x[0],x[1],x[2],x[3],x[4],x[5])
			+tau(x[0],x[2],x[7],x[6],x[8],x[4])
			+tau(x[0],x[7],x[11],x[9],x[10],x[8])
			+tau(x[0],x[11],x[13],x[12],x[14],x[10])
			+tau_analytic(x[0],x[13],x[1],x[15],x[5],x[14]);
			break;
	case 314974315+2:
		*ret= -1.5*global::pt-0.06585
			+tau(x[0],x[1],x[2],x[3],x[4],x[5])
			+tau(x[0],x[2],x[7],x[6],x[8],x[4])
			+tau(x[0],x[7],x[11],x[9],x[10],x[8])
			+tau(x[0],x[11],x[13],x[12],x[14],x[10])
			+tauVc(x[0],x[13],x[1],x[15],x[5],x[14])
			-2.*0.003521-0.01561;
			break;
	case 314974315+3:
		*ret= -1.5*global::pt-0.06585
			+tau(x[0],x[1],x[2],x[3],x[4],x[5])
			+tau(x[0],x[2],x[7],x[6],x[8],x[4])
			+tau(x[0],x[7],x[11],x[9],x[10],x[8])
			+tau(x[0],x[11],x[13],x[12],x[14],x[10])
			+tauVc(x[0],x[13],x[1],x[15],x[5],x[14])
			-2.*0.003521-0.01561;
			break;
	case 314974315+4:
		*ret= -1.5*global::pt-0.06585
			+tau(x[0],x[1],x[2],x[3],x[4],x[5])
			+tau(x[0],x[2],x[7],x[6],x[8],x[4])
			+tau(x[0],x[7],x[11],x[9],x[10],x[8])
			+tau(x[0],x[11],x[13],x[12],x[14],x[10])
			+tauVc(x[0],x[13],x[1],x[15],x[5],x[14]);
			break;

	// Section VI.A.3.1 of Kepler
	case 572068135+0 :
		*ret=tau(x[0],x[1],x[2],x[3],x[4],x[5])
			-0.2529*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			+0.3442;
			break;
	case 572068135+1 :
		*ret=tau_analytic(x[0],x[1],x[2],x[3],x[4],x[5])
			-0.2529*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			+0.3442;
			break;
	case 723700608 :
		*ret=tauVc(x[0],x[1],x[2],x[3],x[4],x[5])
			-0.2529*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			+0.1787;
			break;
	case 560470084+0:
		*ret=tau(x[0],x[1],x[2],x[3],x[4],x[5])
			-0.2529*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			+0.2137 ;
			break;
	case 560470084+1:
		*ret=tau_analytic(x[0],x[1],x[2],x[3],x[4],x[5])
			-0.2529*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			+0.2137 ;
			break;
	case 560470084+2:
		*ret=tauVc(x[0],x[1],x[2],x[3],x[4],x[5])
			-0.2529*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			+0.2137 ;
			break;
	case 560470084+3:
		*ret=tauVc(x[0],x[1],x[2],x[3],x[4],x[5])
			-0.2529*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			+0.2137 ;
			break;
	case 560470084+4:
		*ret=tauVc(x[0],x[1],x[2],x[3],x[4],x[5])
			-0.2529*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			+0.2137 ;
			break;
	case 535502975 :
		*ret=tauVc(x[0],x[1],x[2],x[3],x[4],x[5])
			-0.2529*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			+0.1371 ;
			break;

	// Section VI.A.3.8 (Kepler)
	case 821707685:
		*ret = 1.63-dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]);
		break;
	case 115383627:
		*ret = 1.51-dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]);
        break;
	case 576221766:
		*ret = 1.93-dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]);
        break;
	case 122081309:
		*ret = 1.77-dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]);
        break;
	case 644534985:
		*ret = tauVc(x[0],x[1],x[2],x[3],x[4],x[5])
			-0.2529*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			+0.2391;
			break;
	case 467530297:
		*ret = tauVc(x[0],x[1],x[2],x[3],x[4],x[5])
			-0.2529*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			+0.1376;
        break;
	case 603910880:
		*ret = tauVc(x[0],x[1],x[2],x[3],x[4],x[5])
			-0.2529*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			+0.266;
			break;
	case 135427691:
		*ret = tauVc(x[0],x[1],x[2],x[3],x[4],x[5])
			-0.2529*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			+0.12;
			break;
	case 60314528:
		*ret = 1.16- dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]);
			break;
	case 312132053:
		*ret = tauVc(x[0],x[1],x[2],x[3],x[4],x[5])
			-0.2529*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			+0.1453;
			break;

	// Section VI.A.3.8. Case 2-b (Kepler)
	case 751442360:
		*ret= dihedraly(x[1],x[0],x[2],x[4],x[3],x[5])
			- 0.74; break;
	case 893059266:
		*ret = tauVc(x[0],x[1],x[2],x[3],x[4],x[5])
			-0.2529*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			+0.2391;
			break;
	case 690646028:
		*ret= -dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			+global::pi/2. 
			-0.5*(2.402-x[3]);
			break;

	// Section VI.A.3.9.  (Kepler)
	case 161665083:
		*ret= dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			- 1.78; break;

	// Section VI.A.4.4.  (Kepler)
	case 867513567+1:
		*ret= -dihedraly(x[1],x[0],x[2],x[4],x[3],x[5])
			+0.35*x[1]-0.15*x[0]-0.15*x[2]+0.7022*x[4]-0.17*x[3]
			+0.0123;
			break;
	case 867513567+2:
		*ret= dihedraly(x[1],x[0],x[2],x[4],x[3],x[5])
			-0.13*x[1]+0.631*x[0]+0.31*x[2]-0.58*x[4]
			+0.413*x[3]+0.025*x[5]-2.63363;
			break;
	case 867513567+3:
		*ret= -dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			+0.714*x[0]-0.221*x[1]-0.221*x[2]+0.92*x[3]
			-0.221*x[4]-0.221*x[5]-0.3482;
			break;
	case 867513567+4:
		*ret= dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			-0.315*x[0]+0.3972*x[1]+0.3972*x[2]
			-0.715*x[3]+0.3972*x[4]+0.3972*x[5]
			-2.37095;
			break;
	case 867513567+5:
		*ret= -solid(x[0],x[1],x[2],x[3],x[4],x[5])
			-0.187*x[0]-0.187*x[1]-0.187*x[2]
			+0.1185*x[3]+0.479*x[4]+0.479*x[5]-0.437235;
			break;
	case 867513567+6:
		*ret= solid(x[0],x[1],x[2],x[3],x[4],x[5])
			+0.488*x[0]+0.488*x[1]+0.488*x[2]
			-0.334*x[4]-0.334*x[5]-2.244;
			break;
	case 867513567+70+0:
	case 867513567+70+1:
	case 867513567+70+2:
	case 867513567+70+3:
	case 867513567+70+4:
		{
		double t = vorVc(x[0],x[1],x[2],x[3],x[4],x[5]);
		switch (INEQ_NUMBER-( 867513567+70))
			{
			case 0 : t = gamma(x[0],x[1],x[2],x[3],x[4],x[5]); break;
			case 1 : t = vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5]); break;
			}
		*ret = -t
			-0.145*x[0]-0.081*x[1]-0.081*x[2]-0.133*x[4]-0.133*x[5]
			+1.17401;
		}
		break;
	case 867513567+80+0:
	case 867513567+80+1:
	case 867513567+80+2:
	case 867513567+80+3:
	case 867513567+80+4:
		{
		double t = vorVc(x[0],x[1],x[2],x[3],x[4],x[5]);
		switch (INEQ_NUMBER-( 867513567+80))
			{
			case 0 : t = gamma(x[0],x[1],x[2],x[3],x[4],x[5]); break;
			case 1 : t = vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5]); break;
			}
		*ret = -t
			-0.12*x[0]-0.081*x[1]-0.081*x[2]-0.113*x[4]-0.113*x[5]
			+0.029*x[3] + 0.94903;
		}
		break;
	case 867513567+90+0:
	case 867513567+90+1:
	case 867513567+90+2:
	case 867513567+90+3:
	case 867513567+90+4:
		{
		double t = vorVc(x[0],x[1],x[2],x[3],x[4],x[5]);
		switch (INEQ_NUMBER-( 867513567+90))
			{
			case 0 : t = gamma(x[0],x[1],x[2],x[3],x[4],x[5]); break;
			case 1 : t = vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5]); break;
			}
		*ret = -t
			+1.05382-0.153*(x[3]+x[4]+x[5]);
		}
		break;
	case 867513567+100+0:
	case 867513567+100+1:
	case 867513567+100+2:
	case 867513567+100+3:
	case 867513567+100+4:
		{
		double t = vorVc(x[0],x[1],x[2],x[3],x[4],x[5]);
		switch (INEQ_NUMBER-( 867513567+100))
			{
			case 0 : t = gamma(x[0],x[1],x[2],x[3],x[4],x[5]); break;
			case 1 : t = vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5]); break;
			}
		*ret = 0.0114-t;
		}
		break;

	case 867513567+110+0:
	case 867513567+110+1:
	case 867513567+110+2:
	case 867513567+110+3:
	case 867513567+110+4:
		{
		double t = tauVc(x[0],x[1],x[2],x[3],x[4],x[5]);
		switch (INEQ_NUMBER-( 867513567+110))
			{
			case 0 : t = tau(x[0],x[1],x[2],x[3],x[4],x[5]); break;
			case 1 : t = tau_analytic(x[0],x[1],x[2],x[3],x[4],x[5]); break;
			}
		*ret = t-1.019*global::pt;
		}
		break;

	case 867513567+120+0:
	case 867513567+120+1:
	case 867513567+120+2:
	case 867513567+120+3:
	case 867513567+120+4:
		{
		double t = vorVc(x[0],x[1],x[2],x[3],x[4],x[5]);
		switch (INEQ_NUMBER-( 867513567+120))
			{
			case 0 : t = gamma(x[0],x[1],x[2],x[3],x[4],x[5]); break;
			case 1 : t = vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5]); break;
			}
		*ret = -t+1.449 -0.419351*solid(x[0],x[1],x[2],x[3],x[4],x[5])
			-0.19*x[0]-0.19*x[1]-0.19*x[2];
		}
		break;

	case 867513567+130+0:
	case 867513567+130+1:
	case 867513567+130+2:
	case 867513567+130+3:
	case 867513567+130+4:
		{
		double t = vorVc(x[0],x[1],x[2],x[3],x[4],x[5]);
		switch (INEQ_NUMBER-( 867513567+130))
			{
			case 0 : t = gamma(x[0],x[1],x[2],x[3],x[4],x[5]); break;
			case 1 : t = vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5]); break;
			}
		*ret = -t-0.01465-0.419351*solid(x[0],x[1],x[2],x[3],x[4],x[5])
			+0.0436*(x[4]+x[5])
			+0.079431*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]);
		}
		break;

	// Section VI.A.4.4.2  (Kepler)
	case 867359387+0:
		*ret= 0.114 
			-gamma(x[0],x[1],x[2],x[3],x[4],x[5])
			-gamma(x[0],x[2],x[7],x[6],x[8],x[4])
			-gamma(x[0],x[7],x[11],x[9],x[10],x[8])
			-gamma(x[0],x[11],x[13],x[12],x[14],x[10])
			-gamma(x[0],x[13],x[1],x[15],x[5],x[14]);
			break;
	case 867359387+1:
		*ret= 0.114 
			-gamma(x[0],x[1],x[2],x[3],x[4],x[5])
			-gamma(x[0],x[2],x[7],x[6],x[8],x[4])
			-gamma(x[0],x[7],x[11],x[9],x[10],x[8])
			-gamma(x[0],x[11],x[13],x[12],x[14],x[10])
			-vor_analytic(x[0],x[13],x[1],x[15],x[5],x[14]);
			break;
	case 867359387+2:
		*ret= 0.0875 
			-gamma(x[0],x[1],x[2],x[3],x[4],x[5])
			-gamma(x[0],x[2],x[7],x[6],x[8],x[4])
			-gamma(x[0],x[7],x[11],x[9],x[10],x[8])
			-gamma(x[0],x[11],x[13],x[12],x[14],x[10])
			-vorVc(x[0],x[13],x[1],x[15],x[5],x[14]);
			break;
	case 867359387+3:
		*ret= 0.0875 
			-gamma(x[0],x[1],x[2],x[3],x[4],x[5])
			-gamma(x[0],x[2],x[7],x[6],x[8],x[4])
			-gamma(x[0],x[7],x[11],x[9],x[10],x[8])
			-gamma(x[0],x[11],x[13],x[12],x[14],x[10])
			-vorVc(x[0],x[13],x[1],x[15],x[5],x[14]);
			break;
	case 867359387+4:
		*ret= 0.0875 
			-gamma(x[0],x[1],x[2],x[3],x[4],x[5])
			-gamma(x[0],x[2],x[7],x[6],x[8],x[4])
			-gamma(x[0],x[7],x[11],x[9],x[10],x[8])
			-gamma(x[0],x[11],x[13],x[12],x[14],x[10])
			-vorVc(x[0],x[13],x[1],x[15],x[5],x[14]);
			break;

	// Section VI.A.4.5.1  (Kepler)
	case 498839271+1:
		*ret= dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			-0.636*x[0]+0.462*x[1]+0.462*x[2]
			-0.82*x[3]+0.462*x[4]+0.462*x[5]-1.82419;
			break;
	case 498839271+2:
		*ret= -dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			+0.55*x[0]-0.214*x[1]-0.214*x[2]+1.24*x[3]
			-0.214*x[4]-0.214*x[5]-0.75281;
			break;
	case 498839271+3:
		*ret= dihedraly(x[1],x[0],x[2],x[4],x[3],x[5])
			+0.4*x[0]-0.15*x[1]+0.09*x[2]
			+0.631*x[3]-0.57*x[4]+0.23*x[5]-2.5481;
			break;
	case 498839271+4:
		*ret= -dihedraly(x[1],x[0],x[2],x[4],x[3],x[5])
			-0.454*x[0]+0.34*x[1]+0.154*x[2]
			-0.346*x[3]+0.805*x[4]+0.3429;
			break;
	case 498839271+5:
		*ret= solid(x[1],x[0],x[2],x[4],x[3],x[5])
			+0.065*x[1]+0.065*x[2]+0.061*x[3]-0.115*x[4]
			-0.115*x[5]-0.2618;
			break;
	case 498839271+6:
		*ret= -solid(x[1],x[0],x[2],x[4],x[3],x[5])
			-0.293*x[0]-0.03*x[1]-0.03*x[2]+0.12*x[3]
			+0.325*x[4]+0.325*x[5]-0.2514;
			break;
	case 498839271+7:
		*ret= -gammaNu(x[1],x[0],x[2],x[4],x[3],x[5])
			-0.0538*x[1]-0.0538*x[2]-0.083*x[3]
			-0.0538*x[4]-0.0538*x[5]+0.5995;
			break;
	case 498839271+8:
		*ret= -vorNu(x[1],x[0],x[2],x[4],x[3],x[5])
			-0.0538*x[1]-0.0538*x[2]-0.083*x[3]
			-0.0538*x[4]-0.0538*x[5]+0.5995;
			break;
	case 498839271+9:
		*ret= -gammaNu(x[1],x[0],x[2],x[4],x[3],x[5]);
		break;
	case 498839271+10:
		*ret= -vorNu(x[1],x[0],x[2],x[4],x[3],x[5]);
		break;
	case 498839271+11:
		*ret= tauGnu(x[1],x[0],x[2],x[4],x[3],x[5])
			-0.5945*global::pt;
		break;
	case 498839271+12:
		*ret= tauVnu(x[1],x[0],x[2],x[4],x[3],x[5])
			-0.5945*global::pt;
		break;

	// Section VI.A.4.5.4  (Kepler)
	case 319046543+1 :
		*ret= dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
		-0.49*x[0]+0.44*x[1]+0.44*x[2]-0.82*x[3]+0.44*x[4]+0.44*x[5]
		-2.0421;
		break;
	case 319046543+2 :
		*ret= -dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
		+0.495*x[0]-0.214*x[1]-0.214*x[2]
		+1.05*x[3]-0.214*x[4]-0.214*x[5]-0.2282;
		break;
	case 319046543+3:
		*ret= dihedraly(x[1],x[0],x[2],x[4],x[3],x[5])
		+0.38*x[0]-0.15*x[1]+0.09*x[2]+0.54*x[3]-0.57*x[4]+0.24*x[5]
		-2.3398;
		break;
	case 319046543+4:
		*ret= -dihedraly(x[1],x[0],x[2],x[4],x[3],x[5])
		-0.375*x[0]+0.33*x[1]+0.11*x[2]-0.36*x[3]+0.72*x[4]
		+0.034*x[5] + 0.36135;
		break;
	case 319046543+5:
		*ret= solid(x[0],x[1],x[2],x[3],x[4],x[5])
		+0.42*x[0]+0.165*x[1]+0.165*x[2]-0.06*x[3]-0.135*x[4]-0.135*x[5]
		-1.479;
		break;
	case 319046543+6:
		*ret= -solid(x[0],x[1],x[2],x[3],x[4],x[5])
		-0.265*x[0]-0.06*x[1]-0.06*x[2]+0.124*x[3]+0.296*x[4]+0.296*x[5]
		-0.0997;
		break;
	case 319046543+7:
		*ret= -gammaNu(x[0],x[1],x[2],x[3],x[4],x[5])
		+0.112*x[0]-0.142*x[1]-0.142*x[2]-0.16*x[3]-0.074*x[4]-0.074*x[5]
		+0.9029;
		break;
	case 319046543+8:
		*ret= -vorNu(x[0],x[1],x[2],x[3],x[4],x[5])
		+0.112*x[0]-0.142*x[1]-0.142*x[2]-0.16*x[3]-0.074*x[4]-0.074*x[5]
		+0.9029;
		break;
	case 319046543+9:
		*ret= -gammaNu(x[0],x[1],x[2],x[3],x[4],x[5])
		+0.11-0.07611*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]);
		break;
	case 319046543+10:
		*ret= -vorNu(x[0],x[1],x[2],x[3],x[4],x[5])
		+0.11-0.07611*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]);
		break;
	case 319046543+11:
		*ret= tauGnu(x[0],x[1],x[2],x[3],x[4],x[5])
		-0.07106*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
		+0.06429;
		break;
	case 319046543+12:
		*ret= tauVnu(x[0],x[1],x[2],x[3],x[4],x[5])
		-0.07106*dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
		+0.06429;
		break;
	case 319046543+13:
		*ret= tauGnu(x[0],x[1],x[2],x[3],x[4],x[5])
		-0.0414;
		break;
	case 319046543+14:
		*ret= tauVnu(x[0],x[1],x[2],x[3],x[4],x[5])
		-0.0414;
		break;
	case 319046543+15:
		*ret= -gammaNu(x[0],x[1],x[2],x[3],x[4],x[5])
			-0.015*x[0] -0.16*(x[1]+x[2]+x[3])- 0.0738*(x[4]+x[5])
			+1.29285; break;
	case 319046543+16 : case 319046543+17 :
		*ret= -dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
		+0.495*x[0]-0.214*x[1]-0.214*x[2]
		+1.05*x[3]-0.214*x[4]-0.214*x[5]-0.23545;
		break;
	case 319046543+18:
		*ret= vorVc(x[0],x[1],x[2],x[3],x[4],x[5])-0.03122
		-gammaNu(x[0],x[1],x[2],x[3],x[4],x[5]);
		break;
	case 319046543+19:
		*ret= vorVc(x[0],x[1],x[2],x[3],x[4],x[5])-0.03122
		-vorNu(x[0],x[1],x[2],x[3],x[4],x[5]);
		break;
	case 533270809 :
		*ret= vorVc(x[0],x[1],x[2],x[3],x[4],x[5])+0.007805
			-gammaNu(x[0],x[1],x[2],x[3],x[4],x[5]);
		break;

	// Section VI.A.4.5.5  (Kepler)
	case 365179082:
		*ret= -vorVc(x[0],x[1],x[2],x[3],x[4],x[5])-0.05;
		break;
	case 365179082+1:
		*ret= -vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5])-0.119;
		break;
	case 365179082+2:
		*ret= -vorVc(x[0],x[1],x[2],x[3],x[4],x[5])-0.119;
		break;
	case 368244553:
		*ret= -0.043/2. - vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
	case 820900672:
		*ret=-0.043-vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
		-vorVc(2.,x[1],x[2],x[3],2.,2.);
		break;
	case 961078136:
		*ret=-0.043-vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
		-vorVc(2.51,x[1],x[2],x[3],2.,2.);
		break;

	case 424186517+1:
		*ret = -0.033-vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
	case 424186517+2:
		*ret = -0.058-vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
	case 424186517+3:
		*ret = -0.073-vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;

	// Section VI.A.4.6.1.a  (Kepler)
	case 725257062:
		*ret= -0.212-0.0461+0.137
		-vorVc(x[0],x[1],x[2],x[3],x[4],x[5]);
		break;
	case 977272202:
		*ret= tauVc(x[0],x[1],x[2],x[3],x[4],x[5])
		-0.54525+0.31;
		break;
	// Section VI.A.4.6.1.b  (Kepler)
	case 583626763:
		*ret= -0.212-0.0461
		-vorVc(x[0],2,2,global::sqrt8,x[4],x[5])
		-vorVc(x[2],x[0],2,x[4],2,2)
		-vorVc(x[1],x[0],2,x[5],2,2);
		break;
	case 390951718:
		*ret= -0.54525
		+tauVc(x[0],2,2,global::sqrt8,x[4],x[5])
		+tauVc(x[2],x[0],2,x[4],2,2)
		+tauVc(x[1],x[0],2,x[5],2,2);
		break;

	// Section VI.A.4.6.1.c  (Kepler)
	case 621852152+1:
		*ret = -vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
		-vorVc(2,x[0],x[1],x[5],2,2)
		-vorVc(2,x[0],x[2],x[4],2,2)
		-0.212-0.0461+0.137;
		break;
	case 621852152+2:
		*ret = -vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
		-vorVc(2.51,x[0],x[1],x[5],2,2)
		-vorVc(2,x[0],x[2],x[4],2,2)
		-0.212-0.0461+0.137;
		break;
	case 621852152+3:
		*ret = -vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
		-vorVc(2.51,x[0],x[1],x[5],2,2)
		-vorVc(2.51,x[0],x[2],x[4],2,2)
		-0.212-0.0461+0.137;
		break;
	case 207203174+1:
		*ret = tauVc(x[0],x[1],x[2],x[3],x[4],x[5])
		+tauVc(2,x[0],x[1],x[5],2,2)
		+tauVc(2,x[0],x[2],x[4],2,2)
		-0.54525+0.31;
		break;
	case 207203174+2:
		*ret = tauVc(x[0],x[1],x[2],x[3],x[4],x[5])
		+tauVc(2.51,x[0],x[1],x[5],2,2)
		+tauVc(2,x[0],x[2],x[4],2,2)
		-0.54525+0.31;
		break;
	case 207203174+3:
		*ret = tauVc(x[0],x[1],x[2],x[3],x[4],x[5])
		+tauVc(2.51,x[0],x[1],x[5],2,2)
		+tauVc(2.51,x[0],x[2],x[4],2,2)
		-0.54525+0.31;
		break;

	// Section VI.A.4.6.1.d  (Kepler)
	case 368258024+1:
		*ret= -vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
		-vorVc(x[6],x[1],x[2],x[3],x[7],x[8])
		-vorVc(2,x[0],x[1],x[5],2,2)
		-vorVc(2,x[6],x[1],x[8],2,2)
		-0.212;
		break;
	case 368258024+2:
		*ret= -vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
		-vorVc(x[6],x[1],x[2],x[3],x[7],x[8])
		-vorVc(2.51,x[0],x[1],x[5],2,2)
		-vorVc(2,x[6],x[1],x[8],2,2)
		-0.212;
		break;
	case 368258024+3:
		*ret= -vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
		-vorVc(x[6],x[1],x[2],x[3],x[7],x[8])
		-vorVc(2.51,x[0],x[1],x[5],2,2)
		-vorVc(2.51,x[6],x[1],x[8],2,2)
		-0.212;
		break;
	case 564618342+1:
		*ret= +tauVc(x[0],x[1],x[2],x[3],x[4],x[5])
		+tauVc(x[6],x[1],x[2],x[3],x[7],x[8])
		+tauVc(2,x[0],x[1],x[5],2,2)
		+tauVc(2,x[6],x[1],x[8],2,2)
		-0.54525;
		break;
	case 564618342+2:
		*ret= tauVc(x[0],x[1],x[2],x[3],x[4],x[5])
		+tauVc(x[6],x[1],x[2],x[3],x[7],x[8])
		+tauVc(2.51,x[0],x[1],x[5],2,2)
		+tauVc(2,x[6],x[1],x[8],2,2)
		-0.54525;
		break;
	case 564618342+3:
		*ret= tauVc(x[0],x[1],x[2],x[3],x[4],x[5])
		+tauVc(x[6],x[1],x[2],x[3],x[7],x[8])
		+tauVc(2.51,x[0],x[1],x[5],2,2)
		+tauVc(2.51,x[6],x[1],x[8],2,2)
		-0.54525;
		break;

	// Section VI.A.4.6.1.e  (Kepler)
	case 498774382+1:
		*ret= -vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
		-vorVc(x[6],x[1],x[2],x[3],x[7],x[8])
		-vorVc(2,x[0],x[1],x[5],2,2)
		-vorVc(2,x[6],x[2],x[7],2,2)
		-0.212;
		break;
	case 498774382+2:
		*ret= -vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
		-vorVc(x[6],x[1],x[2],x[3],x[7],x[8])
		-vorVc(2.51,x[0],x[1],x[5],2,2)
		-vorVc(2,x[6],x[2],x[7],2,2)
		-0.212;
		break;
	case 498774382+3:
		*ret= -vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
		-vorVc(x[6],x[1],x[2],x[3],x[7],x[8])
		-vorVc(2.51,x[0],x[1],x[5],2,2)
		-vorVc(2.51,x[6],x[2],x[7],2,2)
		-0.212;
		break;
	case 544865225+1:
		*ret= +tauVc(x[0],x[1],x[2],x[3],x[4],x[5])
		+tauVc(x[6],x[1],x[2],x[3],x[7],x[8])
		+tauVc(2,x[0],x[1],x[5],2,2)
		+tauVc(2,x[6],x[2],x[7],2,2)
		-0.54525;
		break;
	case 544865225+2:
		*ret= +tauVc(x[0],x[1],x[2],x[3],x[4],x[5])
		+tauVc(x[6],x[1],x[2],x[3],x[7],x[8])
		+tauVc(2.51,x[0],x[1],x[5],2,2)
		+tauVc(2,x[6],x[2],x[7],2,2)
		-0.54525;
		break;
	case 544865225+3:
		*ret= +tauVc(x[0],x[1],x[2],x[3],x[4],x[5])
		+tauVc(x[6],x[1],x[2],x[3],x[7],x[8])
		+tauVc(2.51,x[0],x[1],x[5],2,2)
		+tauVc(2.51,x[6],x[2],x[7],2,2)
		-0.54525;
		break;

	// Section VI.A.4.6.2.a  (Kepler)
	case 234734606:
		*ret= -vorVc(x[0],x[1],x[2],x[3],x[4],x[5])-0.221-2*0.009; 
		break;
	case 791682321:
		*ret= tauVc(x[0],x[1],x[2],x[3],x[4],x[5])-0.486+2*0.05925; 
		break;

	// Section VI.A.4.6.2.b  (Kepler)
	case 995351614:
		*ret= -vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
		-vorVc(x[6],x[1],x[2],x[3],x[7],x[8])
		-0.221-0.009;
		break;
	case 321843503:
		*ret= +tauVc(x[0],x[1],x[2],x[3],x[4],x[5])
		+tauVc(x[6],x[1],x[2],x[3],x[7],x[8])
		-0.486+0.05925;
		break;

	// Section VI.A.4.6.2.c  (Kepler)
	case 818985272:
		*ret = -0.19-vorVc(x[0],x[1],x[2],x[3],x[4],x[5]);
		break;
	case 906443824:
		*ret = -0.281+tauVc(x[0],x[1],x[2],x[3],x[4],x[5]);
		break;
	case 547486831:
		*ret = -0.11-vorVc(x[0],x[1],x[2],x[3],x[4],x[5]);
		break;
	case 683897354:
		*ret = -0.205+tauVc(x[0],x[1],x[2],x[3],x[4],x[5]);
		break;

	// Section VI.A.4.6.2.d  (Kepler)
	case 109046923:
		*ret= -0.221-0.0461
		-vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
        -vorVc(x[6],x[1],x[2],x[3],x[7],x[8]);
		break;
	case 642590101:
		*ret= -0.486
		+tauVc(x[0],x[1],x[2],x[3],x[4],x[5])
        +tauVc(x[6],x[1],x[2],x[3],x[7],x[8]);
		break;

	// Section VI.A.4.6.2.e  (Kepler)
	case 160800042:
		*ret = -0.221
		-vorVc(x[0],2,2,x[3],2,2)
		-vorVc(x[1],2,2,x[2],2,2)
		-vorVc(2,2,2,2.51,x[2],x[3]);
		break;
	case 690272881:
		*ret = -0.486
		+tauVc(x[0],2,2,x[3],2,2)
		+tauVc(x[1],2,2,x[2],2,2)
		+tauVc(2,2,2,2.51,x[2],x[3]);
		break;

	// Section VI.A.4.6.2.f  (Kepler)
	case 713930036:
		*ret= -vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
		-vorVc(2,x[0],x[2],x[4],2,2)
		-vorVc(2,x[0],x[1],x[5],2,2)
		-0.221;
		break;
	case 724922588:
		*ret= +tauVc(x[0],x[1],x[2],x[3],x[4],x[5])
		+tauVc(2,x[0],x[2],x[4],2,2)
		+tauVc(2,x[0],x[1],x[5],2,2)
		-0.486;
		break;

	// Section VI.A.4.6.2.g  (Kepler)
	case 821730621+1:
		*ret= -vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
		-vorVc(x[6],x[1],x[2],x[3],x[7],x[8])
		-vorVc(2,x[1],x[6],x[8],2,2)
		-0.221;
		break;
	case 821730621+2:
		*ret= -vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
		-vorVc(x[6],x[1],x[2],x[3],x[7],x[8])
		-vorVc(2.51,x[1],x[6],x[8],2,2)
		-0.221;
		break;
	case 890642961+1:
		*ret= +tauVc(x[0],x[1],x[2],x[3],x[4],x[5])
		+tauVc(x[6],x[1],x[2],x[3],x[7],x[8])
		+tauVc(2,x[1],x[6],x[8],2,2)
		-0.486;
		break;
	case 890642961+2:
		*ret= +tauVc(x[0],x[1],x[2],x[3],x[4],x[5])
		+tauVc(x[6],x[1],x[2],x[3],x[7],x[8])
		+tauVc(2.51,x[1],x[6],x[8],2,2)
		-0.486;
		break;

	// Section VI.A.4.6.3  (Kepler)
	case 341667126:
		*ret= -0.168-0.009
		-vorVc(x[0],x[1],x[2],x[3],x[4],x[5]);
		break;
	case 535906363:
		*ret= -0.352+0.05925
		+tauVc(x[0],x[1],x[2],x[3],x[4],x[5]);
		break;
	case 302085207:
		*ret= -0.168
		-vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
		-vorVc(x[6],x[1],x[2],x[3],x[7],x[8]);
		break;
	case 411491283:
		*ret= -0.352
		+tauVc(x[0],x[1],x[2],x[3],x[4],x[5])
		+tauVc(x[6],x[1],x[2],x[3],x[7],x[8]);
		break;

	// Section VI.A.4.6.8  (Kepler)
	case 516537931:
		*ret= -0.146-vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
	case 130008809+1:
		*ret= -0.31+tauVc(x[0],x[1],x[2],x[3],x[4],x[5])
			+tauVc(x[0],2,x[2],2,x[4],2); break;
	case 130008809+2:
		*ret= -0.31+tauVc(x[0],x[1],x[2],x[3],x[4],x[5])
			+tauVc(x[0],2.51,x[2],2,x[4],2); break;

	// Section VI.A.4.6.10  (Kepler)
	case 286122364:
		*ret= -0.176+tauVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
	case 531861442:
		*ret= -0.084-vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;

	// Section VI.A.4.7.1  (Kepler)
	case 131574415:
		*ret= 1.01 - 0.1*x[0]-0.05*x[1]-0.05*x[2]-0.15*x[4]-0.15*x[5]
			-vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
	case 929773933:
		*ret= 1.1227-0.1*x[0]-0.1*x[1]-0.03*x[2]-0.17*x[4]-0.16*x[5]
			-vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
	case 223261160:
		*ret= 1.0159-0.1*x[0]-0.08*x[1]-0.08*x[2]+0.04*x[3]-0.15*x[4]
			-0.15*x[5]
			-vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
	case 135018647:
		*ret= 1.01054 - 0.1*x[0]-0.06*x[1]-0.06*x[2]-0.04*x[3]
				-0.12*(x[4]+x[5])
			-vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
	case 559676877:
		*ret=-vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
			-vorVc(x[6],x[1],x[2],x[3],x[7],x[8])
		-0.419351*solid(x[0],x[1],x[2],x[3],x[4],x[5])
			-0.419351*solid(x[6],x[1],x[2],x[3],x[7],x[8])
		+0.4542+0.0238*(x[4]+x[5]+x[8]); break;

	// Section VI.A.4.7.2  (Kepler)
	case 587781327:
		*ret = 
			-vorVc(2,2,2,x[0],2,2)
			-vorVc(2,2,2,x[1],2,2)
			-vorVc(2,2,2,2,x[1],x[0]) -0.128;
		break;
	case 807067544:
		*ret = 
			tauVc(2,2,2,x[0],2,2)
			+tauVc(2,2,2,x[1],2,2)
			+tauVc(2,2,2,2,x[1],x[0]) -0.36925;
		break;

	// Section VI.A.4.8  (Kepler)
	case 853728973+1: *ret=
		dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.153; break;
	case 853728973+2: *ret=
		dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.32; break;
	case 853728973+3: *ret=
		dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-0.633; break;
	case 853728973+4: *ret=
		dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.033; break;
	case 853728973+5: *ret=
		dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.033; break;
	case 853728973+6: *ret=
		dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.259; break;
	case 853728973+7: *ret=
		dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-0.817; break;
	case 853728973+8: *ret=
		dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.07; break;
	case 853728973+9: *ret=
		dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.07; break;
	case 853728973+10: *ret=
		dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.23; break;
	case 853728973+11: *ret=
		-dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])+2.28; break;
	case 853728973+12: *ret=
		-dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])+1.624; break;
	case 853728973+13: *ret=
		-dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])+1.929; break;
	case 853728973+14: *ret=
		-dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])+1.507; break;
	case 853728973+15: *ret=
		-dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])+1.761; break;
	case 853728973+16: *ret=
		dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-0.633; break;
	case 853728973+17: *ret=
		dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.033; break;
	case 853728973+18: *ret=
		dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-0.777; break;
	case 853728973+19: *ret=
		-dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])+1.624; break;
	case 853728973+20: *ret=
		-dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])+1.381; break;

	// Section VI.A.4.9  (Kepler)
	case 529738375+1:
		*ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			-0.372*x[0]+0.465*x[1]+0.465*x[2]+0.465*x[4]
			+0.465*x[5]-4.885;
			break;
	case 529738375+2:
		*ret = -2.47277+dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			 -(
			0.291*x[0]-0.393*x[1]-0.586*x[2]+0.79*x[3]-0.321*x[4]
			-0.397*x[5]); break;
	case 529738375+3:
		*ret = -4.45567+dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			-(
			0.291*x[0]-0.393*x[1]-0.586*x[2]-0.321*x[4]-0.397*x[5]); break;
	case 529738375+4:
		*ret = -4.71107+dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			-(
			0.291*x[0]-0.393*x[1]-0.586*x[2]-0.321*x[4]-0.397*x[5]); break;
	case 529738375+5:
		*ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			-0.214*x[0]+0.4*x[1]+0.58*x[2]+0.155*x[4]+0.395*x[5]
			-4.52345;
			break;
	case 529738375+6:
		*ret = +solid(x[0],x[1],x[2],x[3],x[4],x[5])
		-2.71884
		-(
		-0.492*x[0]-0.492*x[1]-0.492*x[2]+0.43*x[3]+0.038*x[4]
		+0.038*x[5]);
		break;
	case 529738375+7:
		*ret = -vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5])
		-0.058*x[0]-0.105*x[1]-0.105*x[2]-0.115*x[3]
		-0.062*x[4]-0.062*x[5]+1.02014; break;
	case 529738375+8:
	case 529738375+9:
		*ret = -vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
		-0.058*x[0]-0.105*x[1]-0.105*x[2]-0.115*x[3]
		-0.062*x[4]-0.062*x[5]+1.02014; break;
	case 529738375+10:
		*ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
		-2.177
		-(
		0.115*x[0]-0.452*x[1]-0.452*x[2]+0.613*x[3]-0.15*x[4]
		-0.15*x[5]); break;
	case 529738375+11:
		*ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
		-2.17382
		-(
		0.115*x[0]-0.452*x[1]-0.452*x[2]+0.618*x[3]-0.15*x[4]
		-0.15*x[5]); break;
	case 529738375+12:
		*ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
		-3.725
		-(
		0.115*x[0]-0.452*x[1]-0.452*x[2]-0.15*x[4]
		-0.15*x[5]); break;
	case 529738375+13:
		*ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
		-3.927
		-(
		0.115*x[0]-0.452*x[1]-0.452*x[2]-0.15*x[4]
		-0.15*x[5]); break;
	case 529738375+14:
		*ret= 0.3085-0.419351*solid(x[0],x[1],x[2],x[3],x[4],x[5])
			-vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5]); break;
	case 529738375+15:
	case 529738375+16:
		*ret= 0.3085-0.419351*solid(x[0],x[1],x[2],x[3],x[4],x[5])
			-vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
	case 529738375+17:
		*ret= -0.121-vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;

	// Section VI.A.4.9.2  (Kepler)
	case 456320257+1:
		*ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
		-2.82998
		-(
		0.47*x[0]-0.522*x[1]-0.522*x[2]+0.812*x[3]-0.522*x[4]
		-0.522*x[5]);
		break;
	case 456320257+2:
		*ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
		-4.8681
		-(
		0.47*x[0]-0.522*x[1]-0.522*x[2]-0.522*x[4]-0.522*x[5]);
		break;
	case 456320257+3:
		*ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
		-5.1623
		-(
		0.47*x[0]-0.522*x[1]-0.522*x[2]-0.522*x[4]-0.522*x[5]);
		break;

	// Section VI.A.4.9.3  (Kepler)
	case 664959245+1:
		*ret = -3.9788+dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
		-(
		-0.4*x[2]+0.15*x[0]-0.09*x[1]-0.631*x[5]-0.23*x[4]);break;
	case 664959245+2:
		*ret = -6.3282+dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
		-(
		0.289*x[0]-0.148*x[1]-1.36*x[2]+0.688*x[3]-0.148*x[4]-1.36*x[5]);
		break;
	case 664959245+3:
		*ret = -4.85746+dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
		-(
		0.289*x[0]-0.148*x[1]-0.723*x[2]-0.148*x[4]-0.723*x[5]);
		break;

	// Section VI.A.4.10 (Kepler)
	case 615073260:
		*ret= x[1]+x[4]+x[7]+x[10]+x[13]
			 +x[2]+x[5]+x[8]+x[11]+x[14]-20.42;
		break;
	case 844430737:
		*ret= x[1]+x[4]+x[7]+x[10]+x[13]
			 +x[2]+x[5]+x[8]+x[11]+x[14]-20.76;
		break;

	// Section VI.A.4.12  (Kepler)
	case 704795925+1:
		*ret= -0.055-gammaNu(x[0],x[1],x[2],x[3],x[4],x[5]); break;
	case 704795925+2:
		*ret= -0.055-vorNu(x[0],x[1],x[2],x[3],x[4],x[5]); break;
	case 704795925+3:
		*ret= -0.092+tauGnu(x[0],x[1],x[2],x[3],x[4],x[5]); break;
	case 704795925+4:
		*ret= -0.092+tauVnu(x[0],x[1],x[2],x[3],x[4],x[5]); break;
	case 332919646+1:
		*ret= -0.039 - gamma(x[0],x[1],x[2],x[3],x[4],x[5]); break;
	case 332919646+2:
		*ret= -0.039 - vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5]); break;
	case 332919646+3:
	case 332919646+4:
	case 332919646+5:
		*ret= -0.039 - vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
	case 332919646+6:
		*ret= -0.094 + tau(x[0],x[1],x[2],x[3],x[4],x[5]); break;
	case 332919646+7:
		*ret= -0.094 + tau_analytic(x[0],x[1],x[2],x[3],x[4],x[5]); break;
	case 332919646+8:
	case 332919646+9:
	case 332919646+10:
		*ret= -0.094 + tauVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;

	case 335795137+1:
		*ret = -0.197-vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5]); break;
	case 335795137+2:
		*ret = -0.239+tau_analytic(x[0],x[1],x[2],x[3],x[4],x[5]); break;
	case 967376139:
		*ret = -vorVc(x[0],x[1],x[2],x[3],x[4],x[5])-
			vorVc(2,x[1],x[2],x[3],2,2) -0.136; break;
	case 666869244:
		*ret = -vorVc(x[0],x[1],x[2],x[3],x[4],x[5])-
			vorVc(2.51,x[1],x[2],x[3],2,2) -0.136; break;
	case 268066802:
		*ret = tauVc(x[0],x[1],x[2],x[3],x[4],x[5])+
			tauVc(2,x[1],x[2],x[3],2,2) -0.224; break;
	case 508108214:
		*ret = tauVc(x[0],x[1],x[2],x[3],x[4],x[5])+
			tauVc(2.51,x[1],x[2],x[3],2,2) -0.224; break;
	case 322505397:
		*ret = -vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
			 -0.125; break;
	case 736616321:
		*ret = -vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
			 +0.011 ;  break;
	case 689417023:
		*ret = tauVc(x[0],x[1],x[2],x[3],x[4],x[5])
			 -0.17 ;  break;
	case 748466752:
		*ret = tauVc(x[0],x[1],x[2],x[3],x[4],x[5])
			 -0.054 ;  break;
	case 369386367+1:
		*ret = -vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
			 -0.126 ;  break;
	case 369386367+2:
		*ret = tauVc(x[0],x[1],x[2],x[3],x[4],x[5])
			 -0.16 ;  break;
	case 724943459+1:
	case 724943459+3:
		*ret = -vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
			 -0.114 ;  break;
	case 724943459+2:
	case 724943459+4:
		*ret = tauVc(x[0],x[1],x[2],x[3],x[4],x[5])
			 -0.186 ;  break;
	case 605071818+1:
		*ret = -vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
			 -0.089 ;  break;
	case 605071818+2:
		*ret = tauVc(x[0],x[1],x[2],x[3],x[4],x[5])
			 -0.154 ;  break;
	case 642806938+1:
		*ret = -vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
			 -0.089 ;  break;
	case 642806938+2:
		*ret = tauVc(x[0],x[1],x[2],x[3],x[4],x[5])
			 -0.154 ;  break;
	case 836331201+1:
		*ret = -vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
			-vorVc(x[6],x[1],x[2],x[3],x[7],x[8])
			 -0.149 ;  break;
	case 836331201+2:
		*ret = tauVc(x[0],x[1],x[2],x[3],x[4],x[5])
			+tauVc(x[6],x[1],x[2],x[3],x[7],x[8])
			 -0.281 ;  break;
	case 327474205+1:
		*ret = -vorVc(x[0],x[1],x[2],x[3],x[4],x[5])
			-vorVc(x[6],x[1],x[2],x[3],x[7],x[8])
			 -0.128 ;  break;
	case 327474205+2:
		*ret = tauVc(x[0],x[1],x[2],x[3],x[4],x[5])
			+tauVc(x[6],x[1],x[2],x[3],x[7],x[8])
			 -0.26625 ;  break;


	// Section VI.A.7  (Kepler)
	case 104506452:
		*ret=octavorVc(x[0],x[1],x[2],x[3],x[4],x[5])
			-0.008-octavor(x[0],x[1],x[2],x[3],x[4],x[5]); break;
	case 601083647:
		*ret= 3-x[3]; break;
	case 543730647:
		*ret= 0.3138- 0.157*x[4] - gamma(x[0],x[1],x[2],x[3],x[4],x[5]);
		break;
	case 163030624:
		*ret= -0.06-gamma(x[0],x[1],x[2],x[3],x[4],x[5]); break;
	case 181462710:
		{
		double t = gamma(x[0],x[1],x[2],x[3],x[4],x[5]);
		*ret = -t+1.4 - 0.1*x[0]-0.15*(x[1]+x[2]+x[4]+x[5]);
		}
		break;

	// Section VI.Appendix 2 (Kepler)
	case 480930831:
		*ret = -0.077 - vorVc(x[0],x[1],x[2],x[3],x[4],x[5]); break;
	case 463544803:
		*ret =vorVc(x[0],x[1],x[2],x[3],x[4],x[5])-
			vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5]); break;
	case 594246986:
		*ret = -gamma(x[0],x[1],x[2],x[3],x[4],x[5])-
			0.145*x[0]-0.08*(x[1]+x[2]) - 0.133*(x[4]+x[5])+
			1.146; break;
	case 381970727:
		*ret = -gamma(x[0],x[1],x[2],x[3],x[4],x[5])-
			0.145*x[0]-0.081*(x[1]+x[2]) - 0.16*(x[4]+x[5])+
			1.255; break;
	case 951798877:
		*ret = -gamma(x[0],x[1],x[2],x[3],x[4],x[5])-
			0.03*x[0]-0.03*(x[1]+x[2]) - 0.094*(x[4]+x[5])+
			0.5361; break;
	case 923397705:
		*ret = -gamma(x[0],x[1],x[2],x[3],x[4],x[5])-
			0.03*x[0]-0.03*(x[1]+x[2]) - 0.16*(x[4]+x[5])+
			0.82; break;
	case 495568072:
		*ret = 1.69*x[3]+x[4]+x[5]-9.0659; break;
	case 378020227:
		*ret = -vor_analytic(x[0],x[1],x[2],x[3],x[4],x[5])-
			0.058*x[0]-0.08*(x[1]+x[2]) - 0.16*x[3]-0.21*(x[4]+x[5])+
			1.7531; break;
	case 256893386:
		*ret = -vorVc(x[0],x[1],x[2],x[3],x[4],x[5])-
			0.058*x[0]-0.1*(x[1]+x[2]) - 0.165*x[3]
			-0.115*x[5] - 0.12*x[4] + 1.38875; break;
			 break;
	case 749955642:
		*ret = x[3]+x[4]+x[5]-7.206; break;
	case 653849975:
		*ret = -vorVc(x[0],x[1],x[2],x[3],x[4],x[5])-
			0.058*x[0]-0.05*(x[1]+x[2]) - 0.16*x[3]-0.13*(x[4]+x[5])+
			1.24547; break;
	



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


static void ConstraintPageK2
	(int numargs,int whichFn,double* x,double* ret,void*)
    {
	*ret = 0;
	switch (INEQ_NUMBER) {

	// Section VI.A.2.5. (Kepler)
	case 269048407+0 :  // upright gamma
	case 553285469+0:
		switch (whichFn) {
		case 1: *ret=radf(x[0],x[1],x[5])-global::sqrt2; break;
		case 2: *ret=radf(x[0],x[2],x[4])-global::sqrt2; break;
		}
		break;
	case 269048407+1 :
	case 553285469+1:
		*ret = global::sqrt2-radf(x[0],x[1],x[5]); break;
	case 293389419+0: // flat gamma
	case 695069283+0:
	case 352079526+0:
		switch (whichFn) {
		case 1: *ret=radf(x[1],x[2],x[3])-global::sqrt2; break;
		case 2: *ret=radf(x[3],x[4],x[5])-global::sqrt2; break;
		case 3: *ret=dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.32; break;
		}
		break;
	case 293389419+1:
	case 695069283+1:
	case 352079526+1:
		switch(whichFn) {
		case 1: *ret = global::sqrt2-radf(x[1],x[2],x[3]); break;
		case 2: *ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.32; break;
		}
		break;
	case 352079526+2:
	case 352079526+3:
		*ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.32; break;
	case 352079526+4:
	case 179025673:
		switch(whichFn) {
		case 1: *ret = global::sqrt2-radf(x[3],x[4],x[5]); break;
		case 2: *ret = dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.32; break;
		}
		break;

	// Section VI.A.2.7 (Kepler)
	// 0,1,2,3,4, gamma,vor,vor0,vor0,vor0 as in VI.2.5 listing.
	case 551665569 : *ret=dips(x); break;
	case 824762926 : *ret= dips(x); break;

	case 675785884+0:
	case 675785884+2:
	case 675785884+3:
		switch(whichFn) {
		case 1 : *ret=dips(x); break;
		case 2 : *ret= 1.32-dihedraly(x[0],x[1],x[13],x[15],x[14],x[5]); break;
		}
		break;
	case 675785884+1:
		switch(whichFn) {
		case 1 : *ret=dips(x); break;
		case 2 : *ret= 1.32-dihedraly(x[0],x[1],x[13],x[15],x[14],x[5]); break;
		case 3 : *ret=global::sqrt2-radf(x[7],x[9],x[11]); break;
		}
		break;
	case 675785884+4:
		switch(whichFn) {
		case 1 : *ret=dips(x); break;
		case 2 : *ret= 1.32-dihedraly(x[0],x[1],x[13],x[15],x[14],x[5]); break;
		case 3 : *ret=global::sqrt2-radf(x[8],x[9],x[10]); break;
		}
		break;

	case 193592217+0:
	case 193592217+2:
	case 193592217+3:
		switch(whichFn) {
		case 1 : *ret=dips(x); break;
		case 2 : *ret= 1.32-dihedraly(x[0],x[1],x[13],x[15],x[14],x[5]); break;
		}
		break;
	case 193592217+1:
		switch(whichFn) {
		case 1 : *ret=dips(x); break;
		case 2 : *ret= 1.32-dihedraly(x[0],x[1],x[13],x[15],x[14],x[5]); break;
		case 3 : *ret=global::sqrt2-radf(x[11],x[12],x[13]); break;
		}
		break;
	case 193592217+4:
		switch(whichFn) {
		case 1 : *ret=dips(x); break;
		case 2 : *ret= 1.32-dihedraly(x[0],x[1],x[13],x[15],x[14],x[5]); break;
		case 3 : *ret=global::sqrt2-radf(x[10],x[12],x[14]); break;
		}
		break;

	// Section VI.A.2.8 (Kepler)
	case 325738864: 
	case 314974315+0:
	case 314974315+2:
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

	// Section VI.A.3.1 of Kepler
	case 572068135+0 :
		*ret=1.51-dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]); break;

	case 560470084+4:
		switch(whichFn) {
		case 1 : *ret=global::sqrt2-radf(x[3],x[4],x[5]);  break;
		case 2 : *ret=1.26-dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		case 3 : *ret=-1.63+dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]); break;
		}
		break;

	// Section VI.A.3.8 (Kepler)
	case 644534985:
	case 467530297:
	case 603910880:
	case 135427691:
		*ret = 1.2 - dihedraly(x[0],x[1],x[2],x[3],x[4],x[5]); break;

	// Section VI.A.3.8. Case 2-b (Kepler)
	case 893059266:
		*ret = -deltay(x[4],2,2,global::sqrt8,2.51,x[5]); break;

	// Section VI.A.3.9.  (Kepler)
	case 161665083:
		*ret= x[1]+x[2]-4.6; break;

	// Section VI.A.4.4.  (Kepler)
	case 867513567+70+0:
	case 867513567+80+0:
	case 867513567+90+0:
	case 867513567+100+0:
	case 867513567+110+0:
	case 867513567+120+0:
	case 867513567+130+0:
		switch(whichFn) {
		case 1 : *ret= -global::sqrt2+radf(x[1],x[2],x[3]); break;
		case 2 : *ret= -global::sqrt2+radf(x[3],x[4],x[5]); break;
		}
		break;
	case 867513567+70+1:
	case 867513567+80+1:
	case 867513567+90+1:
	case 867513567+100+1:
	case 867513567+110+1:
	case 867513567+120+1:
	case 867513567+130+1:
		*ret= global::sqrt2-radf(x[1],x[2],x[3]); break;
	case 867513567+70+4:
	case 867513567+80+4:
	case 867513567+90+4:
	case 867513567+100+4:
	case 867513567+110+4:
	case 867513567+120+4:
	case 867513567+130+4:
		*ret= global::sqrt2-radf(x[3],x[4],x[5]); break;

	// Section VI.A.4.4.2  (Kepler)
	case 867359387+1:
		switch(whichFn) {
		case 1 : *ret= dips(x);  break;
		case 2 : *ret= global::sqrt2-radf(x[1],x[13],x[15]); break;
		}
		break;
	case 867359387+0:
	case 867359387+2:
	case 867359387+3:
		*ret=dips(x); break;
	case 867359387+4:
		switch(whichFn) {
		case 1 : *ret= dips(x);  break;
		case 2 : *ret= global::sqrt2-radf(x[5],x[14],x[15]); break;
		}
		break;

	// Section VI.A.4.5.1  (Kepler)
	case 498839271+7: 
	case 498839271+9:
	case 498839271+11:
		switch(whichFn) {
		case 1 : *ret= -global::sqrt2+radf(x[0],x[1],x[5]); break;
		case 2 : *ret= -global::sqrt2+radf(x[0],x[2],x[4]); break;
		}
		break;
	case 498839271+8:
	case 498839271+10:
	case 498839271+12:
		*ret= global::sqrt2-radf(x[0],x[1],x[5]); break;


	// Section VI.A.4.5.4  (Kepler)
	case 319046543+7:
	case 319046543+9:
	case 319046543+11:
	case 319046543+13:
	case 319046543+15:
		switch(whichFn) {
		case 1 : *ret= -global::sqrt2+radf(x[0],x[1],x[5]); break;
		case 2 : *ret= -global::sqrt2+radf(x[0],x[2],x[4]); break;
		}
		break;
	case 319046543+8:
	case 319046543+10:
	case 319046543+12:
	case 319046543+14:
	case 319046543+19:
		*ret= global::sqrt2-radf(x[0],x[1],x[5]); break;

	case 365179082+1:
	case 365179082+2:
		*ret= global::sqrt2-radf(x[0],x[1],x[5]); break;

	case 820900672:
		{
		double v[9]={x[0],x[1],x[2],x[3],x[4],x[5],2.,2.,2.};
		*ret = -crossdiag(v) + 2.51;
		}
		break;

	case 961078136:
		{
		double v[9]={x[0],x[1],x[2],x[3],x[4],x[5],2.51,2.,2.};
		*ret = -crossdiag(v) + 2.51;
		}
		break;
	case 424186517+1:
		*ret= dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.8; break;
	case 424186517+2:
		*ret= dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-2.5; break;
	case 424186517+3:
		*ret= global::sqrt2-radf(x[0],x[1],x[5]); break;

	// Section VI.A.4.6.1.a  (Kepler)
	case 725257062:
	case 977272202: break;
	// Section VI.A.4.6.1.b  (Kepler)
	case 583626763:
	case 390951718:
		switch(whichFn) {
		case 1 : *ret= -deltay(x[2],x[0],2,x[4],2,2); break;
		case 2 : *ret= -deltay(x[1],x[0],2,x[5],2,2); break;
		}
		break;

	// Section VI.A.4.6.1.c  (Kepler)
	case 621852152+1:
	case 207203174+1:
		switch(whichFn) {
		case 1 : *ret=-deltay(2,x[0],x[1],x[5],2,2); break;
		case 2 : *ret=-deltay(2,x[0],x[2],x[4],2,2); break;
		}
		break;
	case 621852152+2:
	case 207203174+2:
		switch(whichFn) {
		case 1 : *ret=-deltay(2.51,x[0],x[1],x[5],2,2); break;
		case 2 : *ret=-deltay(2,x[0],x[2],x[4],2,2); break;
		}
		break;
	case 621852152+3:
	case 207203174+3:
		switch(whichFn) {
		case 1 : *ret=-deltay(2.51,x[0],x[1],x[5],2,2); break;
		case 2 : *ret=-deltay(2.51,x[0],x[2],x[4],2,2); break;
		}
		break;

	// Section VI.A.4.6.1.d  (Kepler)
	case 368258024+1:
	case 564618342+1:
		{
		double cd = crossdiag(x);
		switch(whichFn) {
		case 1 : *ret=-deltay(2,x[0],x[1],x[5],2,2); break;
		case 2 : *ret=-deltay(2,x[6],x[1],x[8],2,2); break;
		case 3 : *ret=-deltay(x[2],x[0],x[6],cd,2,2); break;
		case 4 : *ret=global::sqrt8-cd; break;
		}
		}
		break;
	case 368258024+2:
	case 564618342+2:
		{
		double cd = crossdiag(x);
		switch(whichFn) {
		case 1 : *ret=-deltay(2.51,x[0],x[1],x[5],2,2); break;
		case 2 : *ret=-deltay(2,x[6],x[1],x[8],2,2); break;
		case 3 : *ret=-deltay(x[2],x[0],x[6],cd,2,2); break;
		case 4 : *ret=global::sqrt8-cd; break;
		}
		}
		break;
	case 368258024+3:
	case 564618342+3:
		{
		double cd = crossdiag(x);
		switch(whichFn) {
		case 1 : *ret=-deltay(2.51,x[0],x[1],x[5],2,2); break;
		case 2 : *ret=-deltay(2.51,x[6],x[1],x[8],2,2); break;
		case 3 : *ret=-deltay(x[2],x[0],x[6],cd,2,2); break;
		case 4 : *ret=global::sqrt8-cd; break;
		}
		}
		break;

	// Section VI.A.4.6.1.e  (Kepler)
	case 498774382+1:
	case 544865225+1:
		{
		double cd = crossdiag(x);
		switch(whichFn) {
		case 1 : *ret=-deltay(2,x[0],x[1],x[5],2,2); break;
		case 2 : *ret=-deltay(2,x[6],x[2],x[7],2,2); break;
		case 3 : *ret=x[3]-cd; break;
		}
		}
		break;
	case 498774382+2:
	case 544865225+2:
		{
		double cd = crossdiag(x);
		switch(whichFn) {
		case 1 : *ret=-deltay(2,x[0],x[1],x[5],2,2); break;
		case 2 : *ret=-deltay(2,x[6],x[2],x[7],2,2); break;
		case 3 : *ret=x[3]-cd; break;
		}
		}
		break;
	case 498774382+3:
	case 544865225+3:
		{
		double cd = crossdiag(x);
		switch(whichFn) {
		case 1 : *ret=-deltay(2,x[0],x[1],x[5],2,2); break;
		case 2 : *ret=-deltay(2,x[6],x[2],x[7],2,2); break;
		case 3 : *ret=x[3]-cd; break;
		}
		}
		break;

	// Section VI.A.4.6.2.d  (Kepler)
	case 109046923:
	case 642590101:
		*ret=x[3]-crossdiag(x); break;

	// Section VI.A.4.6.2.e  (Kepler)
	case 160800042:
	case 690272881:
		switch(whichFn) {
		case 1 : *ret=-deltay(x[0],2,2,x[3],2,2); break;
		case 2 : *ret=-deltay(x[1],2,2,x[2],2,2); break;
		}
		break;

	// Section VI.A.4.6.2.f  (Kepler)
	case 713930036:
	case 724922588:
		switch(whichFn) {
		case 1 : *ret=-deltay(2,x[0],x[2],x[4],2,2); break;
		case 2 : *ret=-deltay(2,x[0],x[1],x[5],2,2); break;
		case 3 : *ret= dihedraly(2,2,2,3.2,2,2)
			- dihedraly(2,2,x[2],2,x[4],2)
			- dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])
			- dihedraly(2,2,x[1],2,x[5],2); break;
		}
		break;

	// Section VI.A.4.6.2.g  (Kepler)
	case 821730621+1:
	case 890642961+1:
		switch(whichFn) {
		case 1 : *ret=-deltay(2,x[1],x[6],x[8],2,2); break;
		case 2 : *ret=3.2-crossdiag(x);  break;
		}
		break;
	case 821730621+2:
	case 890642961+2:
		switch(whichFn) {
		case 1 : *ret=-deltay(2.51,x[1],x[6],x[8],2,2); break;
		case 2 : *ret=3.2-crossdiag(x);  break;
		}
		break;

	// Section VI.A.4.6.3  (Kepler)
	case 302085207:
	case 411491283:
		*ret=3.2-crossdiag(x); 
		break;

	// Section VI.A.4.7  (Kepler)
	case 131574415:
		*ret= dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.9; break;
	case 929773933:
		switch(whichFn)
		{
		case 1:*ret= dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-2; break;
		case 2:*ret= x[1]+x[2]-4.67; break;
		}
		break;
	case 559676877:
		*ret= global::sqrt8-crossdiag(x); break;

	// Section VI.A.4.9 (Kepler)
	case 529738375+9:
	case 529738375+15:
		*ret = global::sqrt2-radf(x[3],x[4],x[5]); break;

	// Section VI.A.4.6.2.b  (Kepler)
	case 995351614:
	case 321843503:
		*ret = global::sqrt8-crossdiag(x); break;

	// Section VI.A.4.10 (Kepler)
	case 615073260:
	case 844430737:
		*ret=dips(x); break;

	// Section VI.A.4.12  (Kepler)
	case 704795925+2:
	case 704795925+4:
		*ret= global::sqrt2-radf(x[0],x[1],x[5]); break;
	case 332919646+2:
	case 332919646+7:
		*ret= global::sqrt2-radf(x[3],x[2],x[1]); break;
	case 332919646+5:
	case 332919646+10:
		*ret= global::sqrt2-radf(x[3],x[4],x[5]); break;
	case 836331201+1:
	case 836331201+2:
	case 327474205+1:
	case 327474205+2:
		*ret = -crossdiag(x) +global::sqrt8;
		break;

	// Section VI.A.7 (Kepler)
	case 104506452:
		*ret = global::sqrt2-radf(x[0],x[2],x[4]); break;
	case 601083647:
		switch(whichFn) {
		case 1: *ret= dihedraly(x[0],x[1],x[2],x[3],x[4],x[5])-1.678; break;
		case 2: *ret=x[1]+x[2]+x[4]+x[5]-8.77; break;
		}
		break;
	case 181462710:
		switch(whichFn) {
		case 1 : *ret= -global::sqrt2+radf(x[1],x[2],x[3]); break;
		case 2 : *ret= -global::sqrt2+radf(x[3],x[4],x[5]); break;
		}
		break;

	// Section VI.Appendix 2 (Kepler)
	case 594246986:
	case 381970727:
		*ret= -global::sqrt2+radf(x[3],x[4],x[5]); break;
	
	case 951798877:
		*ret= 4.3-x[4]-x[5]; break;
	case 923397705:
		*ret = -4.3+x[4]+x[5]; break;
	case 495568072:
		*ret = global::sqrt2-radf(x[3],x[4],x[5]); break;
	case 749955642:
	case 653849975:
		*ret = global::sqrt2-radf(x[3],x[4],x[5]); break;

	//Z-con
	default : cout << "unexpected case in constraint " << INEQ_NUMBER<< endl;
		}
    }

iter::iter(int ineqSwitch) {
	numiter = 20; numargs = 6; nconstr=0;
	switch(ineqSwitch)
		{

	// Section VI.A.2.7 of Kepler.
	case 551665569 :
	case 824762926 :
	case 675785884+0:
	case 675785884+1:
	case 675785884+2:
	case 675785884+3:
	case 675785884+4:
	case 193592217+0:
	case 193592217+1:
	case 193592217+2:
	case 193592217+3:
	case 193592217+4:
			numargs=16; break;

	// Section VI.A.2.8 (Kepler)
	case 325738864: 
	case 314974315+0:
	case 314974315+2:
	case 314974315+3:
	case 314974315+1:
	case 314974315+4:
			numargs=16; break;

	// Section VI.A.4.4.2  (Kepler)
	case 867359387+0:
	case 867359387+1:
	case 867359387+2:
	case 867359387+3:
	case 867359387+4:
			numargs=16; break;
	
	// Section VI.A.4.6.1.d (Kepler)
	case 368258024+1:
	case 564618342+1:
	case 368258024+2:
	case 564618342+2:
	case 368258024+3:
	case 564618342+3:
			numargs=9; break;

	// Section VI.A.4.6.1.e  (Kepler)
	case 498774382+1:
	case 544865225+1:
	case 498774382+2:
	case 544865225+2:
	case 498774382+3:
	case 544865225+3:
			numargs=9; break;

	// Section VI.A.4.6.2.b (Kepler)
	case 995351614:
	case 321843503:
			numargs=9; break;

	// Section VI.A.4.6.2.d  (Kepler)
	case 109046923:
	case 642590101:
			numargs=9; break;

	// Section VI.A.4.6.2.g  (Kepler)
	case 821730621+1:
	case 821730621+2:
	case 890642961+1:
	case 890642961+2:
			numargs=9; break;

	// Section VI.A.4.6.3  (Kepler)
	case 302085207:
	case 411491283:
			numargs=9; break;

	// Section VI.A.4.7.1  (Kepler)
	case 559676877:
			numargs=9; break;

	// Section VI.A.4.10 (Kepler)
	case 615073260:
	case 844430737:
			numargs=16; break;

	// Section VI.A.4.12 (Kepler)
	case 836331201+1:
	case 836331201+2:
	case 327474205+1:
	case 327474205+2:
			numargs=9; break;

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

	// Section VI.A.2.5. (Kepler)
	case 269048407+0 :
		xmin[0]=2.696; xmax[0]=global::sqrt8;
		nconstr=2;
		break;
	case 269048407+1 :
		xmin[0]=2.696; xmax[0]=global::sqrt8;
		nconstr=1;
		break;
	case 553285469+0:
		xmin[0]=2.6; xmax[0]=2.696;
		xmin[3]=2.1;
		nconstr=2;
		break;
	case 553285469+1:
		xmin[0]=2.6; xmax[0]=2.696;
		xmin[3]=2.1;
		nconstr=2;
		break;
	case 293389419+0:
		xmin[3]=2.51; xmax[3]=global::sqrt8;
		nconstr=2;
		break;
	case 293389419+1:
		xmin[3]=2.51; xmax[3]=global::sqrt8;
		nconstr=1;
		break;
	case 695069283+0:
		xmin[3]=2.51; xmax[3]=global::sqrt8;
		xmax[0]=2.17;
		nconstr=2;
		break;
	case 695069283+1:
		xmin[3]=2.51; xmax[3]=global::sqrt8;
		xmax[0]=2.17;
		nconstr=1;
		break;
	case 814398901 :
		xmin[3]=xmax[3]=global::sqrt8;
		break;
	case 352079526+0:
		xmin[3]=2.51; xmax[3]=global::sqrt8;
		nconstr=3;
		break;
	case 352079526+1:
		xmin[3]=2.51; xmax[3]=global::sqrt8;
		nconstr=2;
		break;
	case 352079526+2:
		xmin[3]=2.6; xmax[3]=global::sqrt8;
		xmin[0]=2.2; 
		nconstr=1;
		break;
	case 352079526+3:
		xmin[3]=2.7; xmax[3]=global::sqrt8;
		nconstr=1;
		break;
	case 352079526+4:
		xmin[3]=2.51; xmax[3]=global::sqrt8;
		nconstr=2;
		break;
	case 179025673:
		xmin[3]=2.51; xmax[3]=global::sqrt8;
		nconstr=2;
		break;

	// Section VI.A.2.7 (Kepler)
	case 551665569:
		xmin[9]=xmax[9]=global::sqrt8;
		xmin[15]=xmax[15]=global::sqrt8;
		nconstr=1;
		break;
	case 824762926:
		xmin[12]=xmax[12]=global::sqrt8;
		xmin[15]=xmax[15]=global::sqrt8;
		nconstr=1;
		break;
	case 675785884+0:
		xmin[9]=2.51; xmax[9]=global::sqrt8;
		xmin[15]=xmax[15]=global::sqrt8;
		nconstr=2;
		break;
	case 675785884+1:
	case 675785884+4:
		xmin[9]=2.51; xmax[9]=global::sqrt8;
		xmin[15]=xmax[15]=global::sqrt8;
		nconstr=3;
		break;
	case 675785884+2:
		xmin[9]=2.6; xmax[9]=global::sqrt8;
		xmin[15]=xmax[15]=global::sqrt8;
		xmin[0]=2.2;
		nconstr=2;
		break;
	case 675785884+3:
		xmin[9]=2.7; xmax[9]=global::sqrt8;
		xmin[15]=xmax[15]=global::sqrt8;
		nconstr=2;
		break;
	case 193592217+0:
		xmin[12]=2.51; xmax[12]=global::sqrt8;
		xmin[15]=xmax[15]=global::sqrt8;
		nconstr=2;
		break;
	case 193592217+1:
	case 193592217+4:
		xmin[12]=2.51; xmax[12]=global::sqrt8;
		xmin[15]=xmax[15]=global::sqrt8;
		nconstr=3;
		break;
	case 193592217+2:
		xmin[12]=2.6; xmax[12]=global::sqrt8;
		xmin[15]=xmax[15]=global::sqrt8;
		xmin[0]=2.2;
		nconstr=2;
		break;
	case 193592217+3:
		xmin[12]=2.7; xmax[12]=global::sqrt8;
		xmin[15]=xmax[15]=global::sqrt8;
		nconstr=2;
		break;

	// Section VI.A.2.8 (Kepler)
	case 325738864: 
		xmin[15]=xmax[15]=global::sqrt8;
		nconstr=1;
		break;
	case 314974315+0:
		xmin[15]=2.51; xmax[15]=global::sqrt8;
		nconstr=1;
		break;
	case 314974315+2:
		xmin[15]=2.6; xmax[15]=global::sqrt8;
		xmin[0]=2.2;
		nconstr=1;
		break;
	case 314974315+3:
		xmin[15]=2.7; xmax[15]=global::sqrt8;
		nconstr=1;
		break;
	case 314974315+1:
	case 314974315+4:
		xmin[15]=2.51; xmax[15]=global::sqrt8;
		nconstr=2; 
		break;


	// Section VI.A.3.1 of Kepler
	case 572068135+0 :
		xmin[0]=2.3;
		nconstr=1;
		break;
	case 572068135+1 :
		xmin[0]=2.3;
		break;
	case 723700608 :
		xmin[0]=2.3;
		xmin[5]=global::sqrt8; xmax[5]=3.02;
		break;
	case 560470084+0:
		xmin[0]=2.3;
		xmin[5]=2.51; xmax[5]=global::sqrt8;
		break;
	case 560470084+1:
		xmin[0]=2.3;
		xmin[5]=2.51; xmax[5]=global::sqrt8;
		break;
	case 560470084+2:
		xmin[0]=2.3;
		xmin[2]=2.2;
		xmin[5]=2.6; xmax[5]=global::sqrt8;
		break;
	case 560470084+3:
		xmin[0]=2.3;
		xmin[5]=2.7; xmax[5]=global::sqrt8;
		break;
	case 560470084+4:
		xmin[0]=2.3;
		xmin[5]=2.51; xmax[5]=global::sqrt8;
		nconstr=1;
		break;
	case 535502975 :
		xmin[0]=2.3;
		xmin[5]=2.51; xmax[5]=3.02;
		xmin[4]=2.51; xmax[4]=3.02;
		break;


	// Section VI.A.3.8 (Kepler)
	case 821707685:
		xmin[5]=2.51; xmax[5]=3.2;
		xmax[1]=xmax[2]=2.168;
		break;
	case 115383627:
		xmin[5]=2.51; xmax[5]=3.2;
		xmin[4]=2.51; xmax[4]=3.2;
		xmax[1]=xmax[2]=2.168;
		break;
	case 576221766:
		xmin[5]=2.51; xmax[5]=3.2;
		xmin[3]=xmax[3]=2.51;
		xmax[1]=xmax[2]=2.168;
		break;
	case 122081309:
		xmin[5]=2.51; xmax[5]=3.2;
		xmin[4]=2.51; xmax[4]=3.2;
		xmin[3]=xmax[3]=2.51;
		xmax[1]=xmax[2]=2.168;
		break;
	case 644534985:
		xmin[5]=2.51; xmax[5]=3.2;
		xmax[1]=xmax[2]=2.168;
		nconstr=1;
		break;
	case 467530297:
		xmin[5]=2.51; xmax[5]=3.2;
		xmin[4]=2.51; xmax[4]=3.2;
		xmax[1]=xmax[2]=2.168;
		nconstr=1;
	case 603910880:
		xmin[5]=2.51; xmax[5]=3.2;
		xmin[3]=2.51; xmax[3]=global::sqrt8;
		xmax[1]=xmax[2]=2.168;
		nconstr=1;	
		break;
	case 135427691:
		xmin[5]=2.51; xmax[5]=3.2;
		xmin[4]=2.51; xmax[4]=3.2;
		xmin[3]=2.51; xmax[3]=global::sqrt8;
		xmax[1]=xmax[2]=2.168;
		nconstr=1;
		break;
	case 60314528:
		xmin[5]=2.51; xmax[5]=3.2;
		xmin[4]=2.51; xmax[4]=3.2;
		xmax[3]=2; 
		xmax[1]=xmax[2]=2.168;
		break;
	case 312132053:
		xmin[5]=2.51; xmax[5]=3.488;
		xmin[4]=2.51; xmax[4]=2.51;
		xmax[1]=xmax[2]=2.168;
		break;


	// Section VI.A.3.8. Case 2-b (Kepler)
	case 751442360: 
		xmin[0]=2.51; xmax[0]=2.696;
		xmax[1]=xmax[2]=2.168;
		break;
	case 893059266:
		xmin[4]=2.51; xmax[4]=3.488;
		nconstr=1;
		break;
	case 690646028:
		xmin[4]=2.51; xmax[4]=3.5;
		xmax[1]=xmax[2]=2.168;
		break;

	// Section VI.A.3.9.  (Kepler)
	case 161665083:
		xmin[0]=2.51; xmax[0]=global::sqrt8;
		xmin[3]=xmax[3]=3.2;
		nconstr=1;
		break;

	// Section VI.A.4.4.  (Kepler)
	case 867513567+1:
	case 867513567+2:
	case 867513567+3:
	case 867513567+4:
	case 867513567+5:
	case 867513567+6: 
		xmin[3]=2.51; xmax[3]=global::sqrt8;
		break;
	case 867513567+70+0:
	case 867513567+80+0:
	case 867513567+90+0:
	case 867513567+100+0:
	case 867513567+110+0:
	case 867513567+120+0:
	case 867513567+130+0:
		xmin[3]=2.51; xmax[3]=global::sqrt8;
		nconstr=2;
		break;
	case 867513567+70+1:
	case 867513567+80+1:
	case 867513567+90+1:
	case 867513567+100+1:
	case 867513567+110+1:
	case 867513567+120+1:
	case 867513567+130+1:
		xmin[3]=2.51; xmax[3]=global::sqrt8;
		nconstr=1;
		break;
	case 867513567+70+2:
	case 867513567+80+2:
	case 867513567+90+2:
	case 867513567+100+2:
	case 867513567+110+2:
	case 867513567+120+2:
	case 867513567+130+2:
		xmin[3]=2.6; xmax[3]=global::sqrt8;
		xmin[0]=2.2;
		break;
	case 867513567+70+3:
	case 867513567+80+3:
	case 867513567+90+3:
	case 867513567+100+3:
	case 867513567+110+3:
	case 867513567+120+3:
	case 867513567+130+3:
		xmin[3]=2.7; xmax[3]=global::sqrt8;
		break;
	case 867513567+70+4:
	case 867513567+80+4:
	case 867513567+90+4:
	case 867513567+100+4:
	case 867513567+110+4:
	case 867513567+120+4:
	case 867513567+130+4:
		xmin[3]=2.51; xmax[3]=global::sqrt8;
		nconstr=1;
		break;

	// Section VI.A.4.4.2  (Kepler)
	case 867359387+0:
		xmin[15]=2.51; xmax[15]=global::sqrt8;
		nconstr=1;
		break;
	case 867359387+1:
	case 867359387+4:
		xmin[15]=2.51; xmax[15]=global::sqrt8;
		nconstr=2;
		break;
	case 867359387+2:
		xmin[15]=2.6; xmax[15]=global::sqrt8;
		xmin[0]=2.2;
		nconstr=1;
		break;
	case 867359387+3:
		xmin[15]=2.7; xmax[15]=global::sqrt8;
		nconstr=1;
		break;

	// Section VI.A.4.5.1  (Kepler)
	case 498839271+1:
	case 498839271+2:
	case 498839271+3:
	case 498839271+4:
	case 498839271+5:
	case 498839271+6:
		xmin[0]=2.51; xmax[0]=global::sqrt8;
		break;
	case 498839271+7: 
	case 498839271+9:
	case 498839271+11:
		xmin[0]=2.51; xmax[0]=global::sqrt8;
		nconstr=2;
		break;
	case 498839271+8:
	case 498839271+10:
	case 498839271+12:
		xmin[0]=2.51; xmax[0]=global::sqrt8;
		nconstr=1;
		break;


	// Section VI.A.4.5.4  (Kepler)
	case 319046543+1 :
	case 319046543+2 :
	case 319046543+3:
	case 319046543+4:
	case 319046543+5:
	case 319046543+6:
		xmin[0]=2.51; xmax[0]=2.696;
		break;
	case 319046543+7:
	case 319046543+9:
	case 319046543+11:
	case 319046543+13:
	case 319046543+15:
		xmin[0]=2.51; xmax[0]=2.696;
		nconstr=2;
		break;
	case 319046543+16:
		xmin[0]=2.51; xmax[0]=2.696;
		xmax[3]=2.475;
		break;
	case 319046543+17:
		xmin[0]=2.51; xmax[0]=2.68;
		xmin[3]=2.475;
		break;
	case 319046543+18:
		xmin[0]=2.68; xmax[0]=2.696;
		xmin[3]=2.475;
		break;
	case 319046543+19:
		xmin[0]=2.68; xmax[0]=2.696;
		xmin[3]=2.475;
		nconstr=1;
		break;
	case 319046543+8:
	case 319046543+10:
	case 319046543+12:
	case 319046543+14:
		xmin[0]=2.51; xmax[0]=2.696;
		nconstr=1;
		break;
	case 533270809:
		xmin[0]=2.68; xmax[0]=2.696;
		break;
	

	// Section VI.A.4.5.5  (Kepler)
	case 365179082:
		xmin[0]=2.51; xmax[0]=2.696;
		xmin[3]=2.51; xmax[3]=global::sqrt8;
		break;
	case 365179082+1:
		xmin[0]=2.51; xmax[0]=2.696;
		xmin[3]=2.51; xmax[3]=2.77;
		nconstr=1;
		break;
	case 365179082+2:
		xmin[0]=2.51; xmax[0]=2.696;
		xmin[3]=2.77; xmax[3]=global::sqrt8; 
		nconstr=1;
		break;
	case 368244553:
		xmin[0]=2.51; xmax[0]=2.696;
		xmin[5]=xmax[5]=2.51;
		break;
	case 820900672:
	case 961078136:
		xmin[0]=2.51; xmax[0]=2.696;
		xmin[3]=global::sqrt8; xmax[3]=3.2;
		nconstr=1;
		break;
	case 424186517+1:
		xmax[0]=2.12;
		xmin[3]=global::sqrt8; xmax[3]=3.2;
		nconstr=1;
		break;
	case 424186517+2:
	case 424186517+3:
		xmin[0]=2.51; xmax[0]=2.696;
		xmin[3]=global::sqrt8; xmax[3]=3.2;
		nconstr=1;
		break;
		

	// Section VI.A.4.6.1.a  (Kepler)
	case 725257062:
	case 977272202: 
		xmin[3]=xmax[3]=3.2;
		xmin[4]=global::sqrt8; xmax[4]=3.2;
		break;
	// Section VI.A.4.6.1.b  (Kepler)
	case 583626763:
	case 390951718:
		xmin[4]=3.2; xmax[4]=3.6;
		xmin[5]=3.2; xmax[5]=3.6;
		xmin[3]=xmax[3]=2.;
		nconstr=2;
		break;

	// Section VI.A.4.6.1.c  (Kepler)
	case 621852152+1:
	case 207203174+1:
	case 621852152+2:
	case 207203174+2:
	case 621852152+3:
	case 207203174+3:
		xmin[3]=xmin[4]=xmin[5]=global::sqrt8;
		xmax[3]=xmax[4]=xmax[5]=3.2;
		break;

	// Section VI.A.4.6.1.d  (Kepler)
	case 368258024+1:
	case 564618342+1:
	case 368258024+2:
	case 564618342+2:
	case 368258024+3:
	case 564618342+3:
		xmin[5]=global::sqrt8; xmax[5]=3.2;
		xmin[8]=global::sqrt8; xmax[8]=3.2;
		xmin[3]=3.2; xmax[3]=3.8;
		nconstr=4;
		break;

	// Section VI.A.4.6.1.e  (Kepler)
	case 498774382+1:
	case 544865225+1:
	case 498774382+2:
	case 544865225+2:
	case 498774382+3:
	case 544865225+3:
		xmin[3]=3.2; xmax[3]=3.78;
		xmin[5]=global::sqrt8; xmax[5]=3.2;
		xmin[7]=global::sqrt8; xmax[7]=3.2;
		nconstr=3;
		break;

	// Section VI.A.4.6.2.a  (Kepler)
	case 234734606:
	case 791682321:
		xmin[3]=xmax[3]=global::sqrt8;
		xmin[4]=global::sqrt8; xmax[4]=3.2;
		xmin[5]=2.51; xmax[5]=global::sqrt8;
		break;

	// Section VI.A.4.6.2.b  (Kepler)
	case 995351614:
	case 321843503:
		xmin[3]=3.2; xmax[3]=3.7;
		xmin[5]=xmax[5]=2.51;
		xmin[8]=global::sqrt8; xmax[8]=3.2;
		nconstr=1;
		break;

	// Section VI.A.4.6.2.c  (Kepler)
	case 818985272:
	case 906443824:
		xmax[3]=2.;
		xmin[4]=xmax[4]=global::sqrt8;
		xmin[5]=xmax[5]=3.2;
		break;
	case 547486831:
	case 683897354:
		xmax[3]=2.;
		xmin[4]=xmax[4]=2.51;
		xmin[5]=xmax[5]=3.2;
		break;

	// Section VI.A.4.6.2.d  (Kepler)
	case 109046923:
	case 642590101:
		xmin[5]=xmax[5]=2.;
		xmin[4]=xmax[4]=2.51;
		xmin[8]=xmax[8]=global::sqrt8;
		xmin[7]=xmax[7]=2.;
		xmin[3]=3.2; xmax[3]=3.7;
		nconstr=1;
		break;

	// Section VI.A.4.6.2.e  (Kepler)
	case 160800042:
	case 690272881:
		xmin[2]=xmin[3]=global::sqrt8;
		xmax[2]=xmax[3]=3.2;
		xmax[4]=xmax[5]=2.;
		nconstr=2;
		break;

	// Section VI.A.4.6.2.f  (Kepler)
	case 713930036:
	case 724922588:
		xmax[0]=2.;
		xmin[3]=xmax[3]=2.51;
		xmin[4]=3.2; xmax[4]=3.4;
		xmin[5]=3.2; xmax[5]=3.4;
		nconstr=3;
		break;

	// Section VI.A.4.6.2.g  (Kepler)
	case 821730621+1:
	case 821730621+2:
	case 890642961+1:
	case 890642961+2:
		xmin[5]=xmax[5]=2.51;
		xmax[7]=2.;
		xmax[4]=2.;
		xmin[8]=global::sqrt8; xmax[8]=3.2;
		xmin[3]=3.2; xmax[3]=3.7;
		nconstr=2;
		break;
		

	// Section VI.A.4.6.3  (Kepler)
	case 341667126:
	case 535906363:
		xmin[3]=xmax[3]=global::sqrt8;
		xmin[4]=2.51; xmax[4]=global::sqrt8;
		xmin[5]=2.51; xmax[5]=global::sqrt8;
		break;
	case 302085207:
	case 411491283:
		xmin[8]=xmax[8]=2.51;
		xmax[4]=2;
		xmax[5]=2;
		xmin[7]=xmax[7]=2.51;
		xmin[3]=global::sqrt8; xmax[3]=3.2;
		nconstr=1;
		break;

	// Section VI.A.4.6.8  (Kepler)
	case 516537931:
	case 130008809+1:
	case 130008809+2:
		xmin[4]=global::sqrt8; xmax[4]=3.2;
		xmin[5]=global::sqrt8; xmax[5]=3.2;
		break;

	// Section VI.A.4.6.10  (Kepler)
	case 286122364:
	case 531861442:
		xmin[4]=2.51; xmax[4]=global::sqrt8;
		xmin[5]=global::sqrt8; xmax[5]=3.2;
		break;

	// Section VI.A.4.7.1  (Kepler)
	case 131574415:
		xmin[3]=global::sqrt8; xmax[3]= 3.4;
		nconstr=1;
		xmax[0]=2.2;
		break;
	case 929773933:
		xmin[3]=global::sqrt8; xmax[3]= 3.7; // was 3.4;
		xmin[4]=2.51; xmax[4]=global::sqrt8;
		nconstr=2;
		break;
	case 223261160:
		xmin[3]=global::sqrt8; xmax[3]= 3.;
		xmax[0]=2.08;
		break;
	case 135018647:
		xmin[3]=global::sqrt8; xmax[3]= 3.;
		xmin[4]=2.51; xmax[4]=global::sqrt8;
		break;
	case 559676877:
		xmin[3]=global::sqrt8; xmax[3]=3.2;
		xmin[7]=2.51; xmax[7]=global::sqrt8;
		nconstr=1;
		break;

	// Section VI.A.4.7.2  (Kepler)
	case 587781327:
	case 807067544:
		xmax[2]=xmax[3]=xmax[4]=xmax[5]=2.;
		xmin[0]=global::sqrt8; xmax[0]=3.2;
		xmin[1]=global::sqrt8; xmax[1]=3.23;
		break;

	// Section VI.A.4.8  (Kepler)
	case 853728973+1:
	case 853728973+11:
		xmin[3]=2.51; xmax[3]=global::sqrt8;
		break;
	case 853728973+2:
		xmin[3]= xmax[3]=global::sqrt8;
		break;
	case 853728973+3:
	case 853728973+12:
		xmin[5]=2.51; xmax[5]=global::sqrt8;
		break;
	case 853728973+4:
	case 853728973+13:
		xmin[5]=2.51; xmax[5]=global::sqrt8;
		xmin[3]=2.51; xmax[3]=global::sqrt8;
		break;
	case 853728973+5:
		xmin[5]=2.51; xmax[5]=global::sqrt8;
		xmin[3]=2.51; xmax[3]=2.51;
		break;
	case 853728973+6:
		xmin[5]=2.51; xmax[5]=global::sqrt8;
		xmin[3]= xmax[3]=global::sqrt8;
		break;
	case 853728973+7:
	case 853728973+14:
		xmin[4]=2.51; xmax[4]=global::sqrt8;
		xmin[5]=2.51; xmax[5]=global::sqrt8;
		break;
	case 853728973+8:
	case 853728973+15:
		xmin[4]=2.51; xmax[4]=global::sqrt8;
		xmin[5]=2.51; xmax[5]=global::sqrt8;
		xmin[3]=2.51; xmax[3]=global::sqrt8;
		break;
	case 853728973+9:
		xmin[4]=2.51; xmax[4]=global::sqrt8;
		xmin[5]=2.51; xmax[5]=global::sqrt8;
		xmin[3]=2.51; xmax[3]=2.51;
		break;
	case 853728973+10:
		xmin[4]=2.51; xmax[4]=global::sqrt8;
		xmin[5]=2.51; xmax[5]=global::sqrt8;
		xmin[3]= xmax[3]=global::sqrt8;
		break;
	case 853728973+16:
	case 853728973+19:
		xmin[2]=2.51; xmax[2]=global::sqrt8;
		break;
	case 853728973+17:
		xmin[2]=2.51; xmax[2]=global::sqrt8;
		xmin[3]= xmax[3]=2.51;
		break;
	case 853728973+18:
		xmin[2]=2.51; xmax[2]=global::sqrt8;
		xmin[5]=2.51; xmax[5]=global::sqrt8;
		xmin[3]= xmax[3]=2.51;
		break;
	case 853728973+20:
		xmin[2]=2.51; xmax[2]=global::sqrt8;
		xmin[5]=2.51; xmax[5]=global::sqrt8;
		break;


	// Section VI.A.4.9  (Kepler)
	case 529738375+1:
		xmin[3]=xmax[3]=global::sqrt8;
		break;
	case 529738375+2:
		xmin[3]=2.51; xmax[3]=global::sqrt8;
		xmin[5]=2.51; xmax[5]=global::sqrt8;
		break;
	case 529738375+3:
		xmin[3]=2.51; xmax[3]=2.51;
		xmin[5]=2.51; xmax[5]=global::sqrt8;
		break;
	case 529738375+4:
	case 529738375+5:
		xmin[3]= xmax[3]=global::sqrt8;
		xmin[5]=2.51; xmax[5]=global::sqrt8;
		break;
	case 529738375+6:
	case 529738375+10:
		xmin[4]=2.51; xmax[4]=global::sqrt8;
		xmin[5]=2.51; xmax[5]=global::sqrt8;
		break;
	case 529738375+7:
		xmin[4]=2.51; xmax[4]=2.77;
		xmin[5]=2.51; xmax[5]=2.77;
		break;
	case 529738375+8:
		xmin[4]=2.77; xmax[4]=global::sqrt8;
		xmin[5]=2.51; xmax[5]=global::sqrt8;
		break;
	case 529738375+9:
		xmin[4]=2.51; xmax[4]=2.77;
		xmin[5]=2.51; xmax[5]=2.77;
		nconstr=1;
		break;
	case 529738375+11:
	case 529738375+12:
		xmin[3]=2.51; xmax[3]=global::sqrt8;
		xmin[4]=2.51; xmax[4]=global::sqrt8;
		xmin[5]=2.51; xmax[5]=global::sqrt8;
		break;
	case 529738375+13:
		xmin[3]= xmax[3]=global::sqrt8;
		xmin[4]=2.51; xmax[4]=global::sqrt8;
		xmin[5]=2.51; xmax[5]=global::sqrt8;
		break;
	case 529738375+14:
		xmin[4]=2.51; xmax[4]=2.77;
		xmin[5]=2.51; xmax[5]=2.77;
		break;
	case 529738375+15:
		xmin[4]=2.51; xmax[4]=2.77;
		xmin[5]=2.51; xmax[5]=2.77;
		nconstr=1;
		break;
	case 529738375+16:
		xmin[4]=2.77; xmax[4]=global::sqrt8;
		xmin[5]=2.51; xmax[5]=2.77;
		break;
	case 529738375+17:
		xmin[3]=2.51; xmax[3]=global::sqrt8;
		xmin[4]=2.51; xmax[4]=global::sqrt8;
		xmin[5]=2.51; xmax[5]=global::sqrt8;
		break;

	// Section VI.A.4.9.2  (Kepler)
	case 456320257+1:
	case 456320257+2:
		xmin[0]=2.51; xmax[0]=global::sqrt8;
		xmin[3]=2.51; xmax[3]=global::sqrt8; 
		break;
	case 456320257+3:
		xmin[0]=2.51; xmax[0]=global::sqrt8;
		xmin[3]= xmax[3]=global::sqrt8; 
		break;

	// Section VI.A.4.9.3  (Kepler)
	case 664959245+1:
		xmin[2]=2.51; xmax[2]=global::sqrt8;
		xmin[3]=xmax[3]=2.51;
		break;
	case 664959245+2:
		xmin[2]=2.51; xmax[2]=global::sqrt8;
		xmin[5]=2.51; xmax[5]=global::sqrt8;
		break;
	case 664959245+3:
		xmin[2]=2.51; xmax[2]=global::sqrt8;
		xmin[5]=2.51; xmax[5]=global::sqrt8;
		xmin[3]=xmax[3]=2.51;
		break;

	// Section VI.A.4.10 (Kepler)
	case 615073260:
		xmin[15]=xmax[15]=2.51;
		nconstr=1;
		break;
	case 844430737:
		xmin[15]=xmax[15]=global::sqrt8;
		nconstr=1;
		break;

	// Section VI.A.4.12  (Kepler)
	case 704795925+1:
	case 704795925+3:
		xmin[0]=2.696; xmax[0]=global::sqrt8;
		xmin[1]=2.45; xmax[1]=2.51;
		xmin[5]=2.45; xmax[5]=2.51;
		break;
	case 704795925+2:
	case 704795925+4:
		xmin[0]=2.696; xmax[0]=global::sqrt8;
		xmin[1]=2.45; xmax[1]=2.51;
		xmin[5]=2.45; xmax[5]=2.51;
		nconstr=1;
		break;
	case 332919646+1:
	case 332919646+6:
		xmin[1]=2.45; xmax[1]=2.51;
		xmin[3]=2.51; xmax[3]=global::sqrt8;
		break;
	case 332919646+2:
	case 332919646+7:
		xmin[1]=2.45; xmax[1]=2.51;
		xmin[3]=2.51; xmax[3]=global::sqrt8;
		nconstr=1;
		break;
	case 332919646+3:
	case 332919646+8:
		xmin[0]=2.2;
		xmin[1]=2.45; xmax[1]=2.51;
		xmin[3]=2.6; xmax[3]=global::sqrt8;
		break;
	case 332919646+4:
	case 332919646+9:
		xmin[1]=2.45; xmax[1]=2.51;
		xmin[3]=2.7; xmax[3]=global::sqrt8;
		break;
	case 332919646+5:
	case 332919646+10:
		xmin[1]=2.45; xmax[1]=2.51;
		xmin[3]=2.51; xmax[3]=global::sqrt8;
		nconstr=1;
		break;
	case 335795137+1:
	case 335795137+2:
		xmin[0]=2.696; xmax[0]=global::sqrt8;
		xmin[1]=2.45; xmax[1]=2.51;
		xmin[5]=2.45; xmax[5]=2.51;
		xmin[3]=2.51;  xmax[3]=global::sqrt8;
		break;

	case 967376139:
	case 666869244:
	case 268066802:
	case 508108214:
		xmin[0]=2.696; xmax[0]=global::sqrt8;
		xmin[1]=2.45; xmax[1]=2.51;
		xmin[5]=2.45; xmax[5]=2.51;
		xmin[3]=global::sqrt8;   xmax[3]=3.2;
		break;
	case 322505397:
	case 689417023:
		xmin[0]=2.696; xmax[0]=global::sqrt8;
		xmin[1]=2.45; xmax[1]=2.51;
		xmin[5]=2.45; xmax[5]=2.51;
		xmin[4]=2.51;
		break;
	case 736616321:
	case 748466752:
		xmin[0]=2.696;xmax[0]=global::sqrt8;
		xmin[5]=2.51;
		break;
	case 369386367+1:
	case 369386367+2:
		xmin[0]=2.696; xmax[0]=global::sqrt8;
		xmin[1]=2.45; xmax[1]=2.51;
		xmin[5]=2.45; xmax[5]=2.51;
		xmin[3]=global::sqrt8;   xmax[3]=3.2;
		break;
	case 724943459+1:
	case 724943459+2:
		xmin[0]=2.45; xmax[0]=2.51;
		xmin[4]=2.51; xmax[4]=global::sqrt8;
		xmin[5]=global::sqrt8; xmax[5]=3.2;
		break;
	case 724943459+3:
	case 724943459+4:
		xmin[0]=2.45; xmax[0]=2.51;
		xmin[3]=2.51; xmax[3]=global::sqrt8;
		xmin[5]=global::sqrt8;  xmax[5]=3.2;
		break;
	case 605071818+1:
	case 605071818+2:
		xmin[0]=2.45; xmax[0]=2.51;
		xmin[4]=2.51; xmax[4]=global::sqrt8;
		xmin[5]=2.51; xmax[5]=global::sqrt8;
		break;
	case 642806938+1:
	case 642806938+2:
		xmin[0]=2.45; xmax[0]=2.51;
		xmin[3]=2.51; xmax[3]=global::sqrt8;
		xmin[5]=2.51; xmax[5]=global::sqrt8;
		break;
	case 836331201+1:
	case 836331201+2:
		xmin[1]=2.45;
		xmin[3]=global::sqrt8; xmax[3]=3.5;
		xmin[5]=2.51; xmax[5]=3.2;
		nconstr=1;
		break;
	case 327474205+1:
	case 327474205+2:
		xmin[1]=2.45;
		xmin[3]=global::sqrt8;  xmax[3]=3.5;  
		xmin[5]=global::sqrt8; xmax[5]=3.2;
		nconstr=1;
		break;

	// Section VI.A.7  (Kepler)
	case 104506452:
		xmin[0]=2.51; xmax[0]=2.696;
		nconstr=1;
		break;
	case 601083647:
		xmin[3]=2.51; xmax[3]=3.3;
		nconstr=2;
		break;
	case 543730647:
		xmin[3]=2.51; xmax[3]=2.6;
		xmax[4]=2.138;
		break;
	case 163030624:
		xmin[3]=2.51; xmax[3]=global::sqrt8;
		xmin[4]=2.22; xmax[4]=2.238;
		xmin[1]=2.121;
		break;
	case 181462710:
		xmin[3]=2.51; xmax[3]=global::sqrt8;
		xmax[0]=xmax[1]=xmax[2]=2.2;
		xmax[4]=xmax[5]=2.35;
		nconstr=2;
		break;

	// Section VI.Appendix 2 (Kepler)


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
	minimize2(695069283+0);
	minimize2(695069283+1);
	/*
	minimize2(814398901 );

	minimize2(352079526+0);
	minimize2(352079526+1);
	minimize2(352079526+2);
	minimize2(352079526+3);
	minimize2(352079526+4);
	minimize2(179025673);

	cout << "\n\n------ Section VI.A.2.7 of Kepler. ----- " << endl;

	// Section VI.A.2.7 
	minimize2(551665569);
	minimize2(824762926);
	minimize2(675785884+0);
	minimize2(675785884+1);
	minimize2(675785884+2);
	minimize2(675785884+3);
	minimize2(675785884+4);
	minimize2(193592217+0);
	minimize2(193592217+1);
	minimize2(193592217+2);
	minimize2(193592217+3);
	minimize2(193592217+4);

	cout << "\n\n------ Section VI.A.2.8 of Kepler ----------" << endl;
	minimize2(325738864); 
	minimize2(314974315+0);
	minimize2(314974315+1);
	minimize2(314974315+2);
	minimize2(314974315+3);
	minimize2(314974315+4);
	cout << "\n\n------ Section VI.A.3.1 of Kepler ----------" << endl;
	// Section VI.A.3.1 of Kepler
	minimize2(572068135+0 );
	minimize2(572068135+1 );
	minimize2(723700608 );
	minimize2(560470084+0);
	minimize2(560470084+1);
	minimize2(560470084+2);
	minimize2(560470084+3);
	minimize2(560470084+4);
	minimize2(535502975 );

	cout << "\n\n------ Section VI.A.3.8 of Kepler ----------" << endl;
	// Section VI.A.3.8 (Kepler)
	minimize2(821707685);
	minimize2(115383627);
	minimize2(576221766);
	minimize2(122081309);
	minimize2(644534985);
	minimize2(467530297);
	minimize2(603910880);
	minimize2(135427691);
	minimize2(60314528);
	minimize2(312132053);

	cout << "\n\n------ Section VI.A.3.8 Case 2-b of Kepler --" << endl;
	// Section VI.A.3.8. Case 2-b (Kepler)
	minimize2(751442360); 
	minimize2(893059266);
	minimize2(690646028);


	cout << "\n\n------ Section VI.A.3.9 of Kepler --" << endl;
	// Section VI.A.3.9.  (Kepler)
	minimize2(161665083);

	cout << "\n\n------ Section VI.A.4.4 of Kepler --" << endl;
	// Section VI.A.4.4.  (Kepler)
	minimize2(867513567+1);
	minimize2(867513567+2);
	minimize2(867513567+3);
	minimize2(867513567+4);
	minimize2(867513567+5);
	minimize2(867513567+6); 
	minimize2(867513567+70+0);
	minimize2(867513567+70+1);
	minimize2(867513567+70+2);
	minimize2(867513567+70+3);
	minimize2(867513567+70+4);
	minimize2(867513567+80+0);
	minimize2(867513567+80+1);
	minimize2(867513567+80+2);
	minimize2(867513567+80+3);
	minimize2(867513567+80+4);
	minimize2(867513567+90+0);
	minimize2(867513567+90+1);
	minimize2(867513567+90+2);
	minimize2(867513567+90+3);
	minimize2(867513567+90+4);
	minimize2(867513567+100+0);
	minimize2(867513567+100+1);
	minimize2(867513567+100+2);
	minimize2(867513567+100+3);
	minimize2(867513567+100+4);
	minimize2(867513567+110+0);
	minimize2(867513567+110+1);
	minimize2(867513567+110+2);
	minimize2(867513567+110+3);
	minimize2(867513567+110+4);
	minimize2(867513567+120+0);
	minimize2(867513567+120+1);
	minimize2(867513567+120+2);
	minimize2(867513567+120+3);
	minimize2(867513567+120+4);
	minimize2(867513567+130+0);
	minimize2(867513567+130+1);
	minimize2(867513567+130+2);
	minimize2(867513567+130+3);
	minimize2(867513567+130+4);

	cout << "\n\n------ Section VI.A.4.4.2 of Kepler --" << endl;
	// Section VI.A.4.4.2  (Kepler)
	minimize2(867359387+0);
	minimize2(867359387+1);
	minimize2(867359387+2);
	minimize2(867359387+3);
	minimize2(867359387+4);

	cout << "\n\n------ Section VI.A.4.5.1 of Kepler --" << endl;
	// Section VI.A.4.5.1  (Kepler)
	minimize2(498839271+1);
	minimize2(498839271+2);
	minimize2(498839271+3);
	minimize2(498839271+4);
	minimize2(498839271+5);
	minimize2(498839271+6);
	minimize2(498839271+7); 
	minimize2(498839271+8);
	minimize2(498839271+9);
	minimize2(498839271+10);
	minimize2(498839271+11);
	minimize2(498839271+12);

	cout << "\n\n------ Section VI.A.4.5.4 of Kepler --" << endl;
	// Section VI.A.4.5.4  (Kepler)
	minimize2(319046543+1 );
	minimize2(319046543+2 );
	minimize2(319046543+3);
	minimize2(319046543+4);
	minimize2(319046543+5);
	minimize2(319046543+6);
	minimize2(319046543+7);
	minimize2(319046543+8);
	minimize2(319046543+9);
	minimize2(319046543+10);
	minimize2(319046543+11);
	minimize2(319046543+12);
	minimize2(319046543+13);
	minimize2(319046543+14);
	minimize2(319046543+15);
	minimize2(319046543+16);
	minimize2(319046543+17);
	minimize2(319046543+18);
	minimize2(319046543+19);
	minimize2(533270809);
	

	cout << "\n\n------ Section VI.A.4.5.5 of Kepler --" << endl;
	// Section VI.A.4.5.5  (Kepler)
	minimize2(365179082);
	minimize2(365179082+1);
	minimize2(365179082+2);
	minimize2(368244553);
	minimize2(820900672);
	minimize2(961078136);
	minimize2(424186517+1);
	minimize2(424186517+2);
	minimize2(424186517+3);
	

	cout << "\n\n------ Section VI.A.4.6.1.a of Kepler --" << endl;
	// Section VI.A.4.6.1.a  (Kepler)
	minimize2(725257062);
	minimize2(977272202); 
	cout << "\n\n------ Section VI.A.4.6.1.b of Kepler --" << endl;
	// Section VI.A.4.6.1.b  (Kepler)
	minimize2(583626763);
	minimize2(390951718);


	cout << "\n\n------ Section VI.A.4.6.1.c of Kepler --" << endl;
	// Section VI.A.4.6.1.c  (Kepler)
	minimize2(621852152+1);
	minimize2(207203174+1);
	minimize2(621852152+2);
	minimize2(207203174+2);
	minimize2(621852152+3);
	minimize2(207203174+3);


	cout << "\n\n------ Section VI.A.4.6.1.d of Kepler --" << endl;
	// Section VI.A.4.6.1.d  (Kepler)
	minimize2(368258024+1);
	minimize2(564618342+1);
	minimize2(368258024+2);
	minimize2(564618342+2);
	minimize2(368258024+3);
	minimize2(564618342+3);

	cout << "\n\n------ Section VI.A.4.6.1.e of Kepler --" << endl;
	// Section VI.A.4.6.1.e  (Kepler)
	minimize2(498774382+1);
	minimize2(544865225+1);
	minimize2(498774382+2);
	minimize2(544865225+2);
	minimize2(498774382+3);
	minimize2(544865225+3);

	cout << "\n\n------ Section VI.A.4.6.2.a of Kepler --" << endl;
	// Section VI.A.4.6.2.a  (Kepler)
	minimize2(234734606);
	minimize2(791682321);

	cout << "\n\n------ Section VI.A.4.6.2.b of Kepler --" << endl;
	// Section VI.A.4.6.2.b  (Kepler)
	minimize2(995351614);
	minimize2(321843503);

	cout << "\n\n------ Section VI.A.4.6.2.c of Kepler --" << endl;
	// Section VI.A.4.6.2.c  (Kepler)
	minimize2(818985272);
	minimize2(906443824);
	minimize2(547486831);
	minimize2(683897354);

	cout << "\n\n------ Section VI.A.4.6.2.d of Kepler --" << endl;
	// Section VI.A.4.6.2.d  (Kepler)
	minimize2(109046923);
	minimize2(642590101);

	cout << "\n\n------ Section VI.A.4.6.2.e of Kepler --" << endl;
	// Section VI.A.4.6.2.e  (Kepler)
	minimize2(160800042);
	minimize2(690272881);

	cout << "\n\n------ Section VI.A.4.6.2.f of Kepler --" << endl;
	// Section VI.A.4.6.2.f  (Kepler)
	minimize2(713930036);
	minimize2(724922588);

	cout << "\n\n------ Section VI.A.4.6.2.g of Kepler --" << endl;
	// Section VI.A.4.6.2.g  (Kepler)
	minimize2(821730621+1);
	minimize2(821730621+2);
	minimize2(890642961+1);
	minimize2(890642961+2);

	cout << "\n\n------ Section VI.A.4.6.3 of Kepler --" << endl;
	// Section VI.A.4.6.3  (Kepler)
	minimize2(341667126);
	minimize2(535906363);
	minimize2(302085207);
	minimize2(411491283);

	cout << "\n\n------ Section VI.A.4.6.8 of Kepler --" << endl;
	// Section VI.A.4.6.8  (Kepler)
	minimize2(516537931);
	minimize2(130008809+1);
	minimize2(130008809+2);


	cout << "\n\n------ Section VI.A.4.6.10 of Kepler --" << endl;
	// Section VI.A.4.6.10  (Kepler)
	minimize2(286122364);
	minimize2(531861442);


	// Section VI.A.4.7.1  (Kepler)
	minimize2(131574415);
	minimize2(929773933);
	minimize2(223261160);
	minimize2(135018647);
	minimize2(559676877);


	cout << "\n\n------ Section VI.A.4.7.2 of Kepler --" << endl;
	// Section VI.A.4.7.2  (Kepler)
	minimize2(587781327);
	minimize2(807067544);
	cout << "\n\n------ Section VI.A.4.8 of Kepler --" << endl;
	// Section VI.A.4.8  (Kepler)
	minimize2(853728973+1);
	minimize2(853728973+2);
	minimize2(853728973+3);
	minimize2(853728973+4);
	minimize2(853728973+5);
	minimize2(853728973+6);
	minimize2(853728973+7);
	minimize2(853728973+8);
	minimize2(853728973+9);
	minimize2(853728973+10);
	minimize2(853728973+11);
	minimize2(853728973+12);
	minimize2(853728973+13);
	minimize2(853728973+14);
	minimize2(853728973+15);
	minimize2(853728973+16);
	minimize2(853728973+17);
	minimize2(853728973+18);
	minimize2(853728973+19);
	minimize2(853728973+20);

	cout << "\n\n------ Section VI.A.4.9.1 of Kepler --" << endl;
	// Section VI.A.4.9  (Kepler)
	minimize2(529738375+1);
	minimize2(529738375+2);
	minimize2(529738375+3);
	minimize2(529738375+4);
	minimize2(529738375+5);
	minimize2(529738375+6);
	minimize2(529738375+7);
	minimize2(529738375+8);
	minimize2(529738375+9);
	minimize2(529738375+10);
	minimize2(529738375+11);
	minimize2(529738375+12);
	minimize2(529738375+13);
	minimize2(529738375+14);
	minimize2(529738375+15);
	minimize2(529738375+16);
	minimize2(529738375+17);

	cout << "\n\n------ Section VI.A.4.9.2 of Kepler --" << endl;
	// Section VI.A.4.9.2  (Kepler)
	minimize2(456320257+1);
	minimize2(456320257+2);
	minimize2(456320257+3);

	cout << "\n\n------ Section VI.A.4.9.3 of Kepler --" << endl;
	// Section VI.A.4.9.3  (Kepler)
	minimize2(664959245+1);
	minimize2(664959245+2);
	minimize2(664959245+3);

	cout << "\n\n------ Section VI.A.4.10 of Kepler --" << endl;
	// Section VI.A.4.10 (Kepler)
	minimize2(615073260);
	minimize2(844430737);

	// Section VI.A.4.12  (Kepler)
	minimize2(704795925+1);
	minimize2(704795925+2);
	minimize2(704795925+3);
	minimize2(704795925+4);
	minimize2(332919646+1);
	minimize2(332919646+2);
	minimize2(332919646+3);
	minimize2(332919646+4);
	minimize2(332919646+5);
	minimize2(332919646+6);
	minimize2(332919646+7);
	minimize2(332919646+8);
	minimize2(332919646+9);
	minimize2(332919646+10);
	minimize2(335795137+1);
	minimize2(335795137+2);
	minimize2(967376139);
	minimize2(666869244);
	minimize2(268066802);
	minimize2(508108214);
	minimize2(322505397);
	minimize2(736616321);
	minimize2(689417023);
	minimize2(748466752);
	minimize2(369386367+1);
	minimize2(369386367+2);
	minimize2(724943459+1);
	minimize2(724943459+2);
	minimize2(724943459+3);
	minimize2(724943459+4);
	minimize2(605071818+1);
	minimize2(605071818+2);
	minimize2(642806938+1);
	minimize2(642806938+2);
	minimize2(836331201+1);
	minimize2(836331201+2);
	minimize2(327474205+1);
	minimize2(327474205+2);

	// Section VI.A.7  (Kepler)
	minimize2(104506452);
	minimize2(601083647);
	minimize2(543730647);
	minimize2(163030624);
	minimize2(181462710);


	// Section VI.Appendix 2 (Kepler)
	minimize2(480930831);
	minimize2(463544803);
	minimize2(594246986);
	minimize2(381970727);
	minimize2(951798877);
	minimize2(923397705);
	minimize2(495568072);
	minimize2(378020227);
	minimize2(256893386);
	minimize2(749955642);
	minimize2(653849975);
	
	// BREAK
	*/
	return;
	}