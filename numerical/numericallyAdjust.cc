#include <iomanip.h>
#include <iostream.h>
#include <fstream.h>
#include <math.h>
#include <stdlib.h>
#include "numerical.h"
#include "gradient.h"
#include "morefn.h"
#include "constants.h"

// copied from ineq.cc on 8/21/1997, to give a version of
// ineq.cc that works with the callable CPLEX library

static double gamma141(double y[6]);
static double vor141(double y[6]);
class trace {
 
public:
    virtual double getY(int index)=0;         // 0..5
    virtual double getDih(int index)=0;       // 0..2
    virtual double getS()=0;
    virtual void setY(int index,double y)=0;
    virtual void setDih(int index,double d)=0;
    virtual void setS(double s)=0;
};

class option {
public:
	enum fnChoice { gamma,vor,sean };
private:
	static fnChoice r;
public:
	static fnChoice getChoice() { return r; }
	static void setChoice(fnChoice x) { r = x; }
	static double sigmafn(double,double,double,double,double,double);
};

static double vorvol(double y[6]);
option::fnChoice option::r = option::gamma;

double option::sigmafn(double y1,double y2,double y3,double y4,double y5,
	double y6)
	{
	double y[6] = {y1,y2,y3,y4,y5,y6};
	fnChoice opt = option::getChoice();
	switch(opt)
	{
	case option::gamma : return gamma141(y); 
	case option::vor : return vor141(y);
	case option::sean : return -vorvol(y);
	default : cout << "error, cases not found" << endl;
	}
	return 0;
	}

class tData : public valuator {
public: 
	double y[6];
	double d[3];
	double s;
	tData() {}
	tData(trace& t)
		{
		int i;
		for (i=0;i<6;i++) y[i]=t.getY(i);
		for (i=0;i<3;i++) d[i]=t.getDih(i);
		s = t.getS();
		}

	double evalf(double z[6])
		{
		return -y[0]*z[0]-y[1]*z[1]-y[2]*z[2]-y[3]*z[3]-y[4]*z[4]-y[5]*z[5]-
			d[0]*dihedraly(z[0],z[1],z[2],z[3],z[4],z[5])
			-d[1]*dihedraly(z[1],z[0],z[2],z[4],z[3],z[5])
			-d[2]*dihedraly(z[2],z[1],z[0],z[5],z[4],z[3])
			- s + option::sigmafn(z[0],z[1],z[2],z[3],z[4],z[5]);
		}
	};



static double vorvol(double y[6])
	{
	return (vor_analytic(y[0],y[1],y[2],y[3],y[4],y[5])
		-(4.0/3.0)*solid(y[0],y[1],y[2],y[3],y[4],y[5]))/(-4.0*doct);
	}


static double vorIP(double y[6]) { return
	vorVc(y[0],y[1],y[2],y[3],y[4],y[5],1.385);
	}

static double mustbepositive(double x)
	{
	return (x>0.0 ? 0.0 : -5000.0*x*x);
	}





/********************************** MEMO *******************************/
/*
static double triangle(double y1,double y2,double y6)
	{
	return acos( (y1*y1+y2*y2-y6*y6)/(2*y1*y2) );
	}
static double psi (double y) { return triangle(y,1.255,1.6); }
static double beta (double y0,double y1,double y6)
	{
	double p = cos(psi(y1));
	double t = cos(triangle(y0,y1,y6));
	return acos( sqrt((p*p-t*t)/(1.0-t*t)) );
	}

static double betahack(double y0,double y1, double y6)
	{
	double p = y0/2.77;
	double t = cos(triangle(y0,y1,y6));
	return acos ( sqrt((p*p-t*t)/(1.0-t*t)));
	}
*/
static double beta0(double y0,double y1,double y5)
	{
	if (delta(y0*y0,y1*y1,1.255*1.255,1.6*1.6,1.6*1.6,y5*y5)>0.002)
	return dihedraly(y0,y1,1.255,1.6,1.6,y5);
	return 0;
	}
static double rad(double y[6])
	{
	return sqrt(circum2(y[0]*y[0],y[1]*y[1],y[2]*y[2],y[3]*y[3],
			y[4]*y[4],y[5]*y[5]));
	}

static double gamma141(double y[6])
	{
	return mustbepositive(-rad(y)+1.41) + gamma(y[0],y[1],y[2],y[3],y[4],y[5]);
	}

static double vor141(double y[6])
	{
	return mustbepositive(+rad(y)-1.41) 
		+ vor_analytic(y[0],y[1],y[2],y[3],y[4],y[5]);
	}

static void PARTIV(tData& u,double y[],double ymin[6],double ymax[6])
	{
	int i;
	int numargs=6;
	int numiters= 20; // was 20;
	double t = unconstrained_iterated_max(numiters,ymin,
			ymax,y,numargs,u);
	cout << "adjustment factor (max) = " << t << ", {";
	for (i=0;i<numargs;i++) cout << y[i] << (i+1<numargs ? ",": "} ");
	for (i=0;i<6;i++) cout << ymin[i] << " ";
	cout << ";";
	for (i=0;i<6;i++) cout << ymax[i] << " ";
	cout << endl;
	u.s -= t;
	}

void numericallyAdjust(trace& t,double ymin[6],double ymax[6])
	{
	tData u(t);
	double y[6];
	option::setChoice(option::sean);
	cout << "WARNING: using McLaughlin's parameters.\n\n\n" << flush;
	PARTIV(u,y,ymin,ymax);
	tData v(t);
	//BUG:
	t.setS(u.s);
	/*
	option::setChoice(option::vor);
	PARTIV(v,y,ymin,ymax);
	if (u.s<v.s) t.setS(u.s); else t.setS(v.s);
	*/
	}


/*
class traceD : public trace {
	private:
	double y[6];
	double d[3];
	double s;
	public:
	double getY(int i) { if ((i<0)||(i>5)) return 0; return y[i]; }
	double getDih(int i) { if ((i<0)||(i>2)) return 0; return d[i]; }
	double getS() { return s; }
	void setY(int i,double y0) { if ((i<0)||(i>5)) return; y[i]=y0; }
	void setDih(int i,double d0) { if ((i<0)||(i>2)) return; d[i]=d0; }
	void setS(double s0) { s = s0; }
	traceD(tData t) 
		{
		int i;
		for (i=0;i<6;i++) y[i]= t.y[i];
		for (i=0;i<3;i++) d[i]= t.d[i];
		s = t.s;
		}
	};
*/
