#include <iomanip.h>
#include <stdlib.h>
#include "numerical.h"
#include "constants.h"
#include "morefn.h"
#include <math.h>



double fa(double y1,double y2,double y3,double y4,double y5,double y6)
	{
	return gamma(y1,y2,y3,y4,y5,y6) 
		- 3.0508*(dihedraly(y2,y1,y3,y5,y4,y6)+ dihedraly(y3,y1,y2,y6,y4,y5)) 
		+ 9.494/4.0;
	}

double fb(double y1,double y2,double y3,double y4,double y5,double y6)
	{
	return gamma(y1,y2,y3,y4,y5,y6) 
		- 3.0508*(dihedraly(y2,y1,y3,y5,y4,y6))+
		9.494/4.0;
	}

double fc(double y1,double y2,double y3,double y4,double y5,double y6)
	{
	return gamma(y1,y2,y3,y4,y5,y6) 
		 + 9.494/4.0;
	}

double fd(double y1,double y2,double y3,double y4,double y5,double y6)
	{
	return gamma(y1,y2,y3,y4,y5,y6) 
		- 3.0508*(dihedraly(y3,y1,y2,y6,y4,y5)) 
		+ 9.494/4.0;
	}

double dpi(double y1,double y2,double y3,double y4,double y5,double y6)
	{
	return dihedraly(y1,y2,y3,y4,y5,y6)-global::pi2;
	}

class sign {
	public:
	double x;
	sign(double u) { x = u; }
	sign() {}
	};
ostream &operator<< (ostream& stream,sign u)
	{
	if (u.x>=0.0) stream << "+";
	stream << u.x;
	return stream;
	}


void outA(double y[6])
	{
	cout << " cA - t ";
	cout << sign(y[0]) << " La1 ";
	cout << sign(y[1]) << " La2 ";
	cout << sign(y[2]) << " La3 ";
	cout << sign(y[4]) << " La5 ";
	cout << sign(y[5]) << " La6 ";
	cout << sign(dpi(y[0],y[1],y[2],y[3],y[4],y[5])) << " Lad ";
	cout << " < " << -fa(y[0],y[1],y[2],y[3],y[4],y[5]);
	cout << endl;
	}

void outB(double y[6])
	{
	cout << " cB - t ";
	cout << sign(y[0]) << " Lb1 ";
	cout << sign(y[1]) << " Lb2 ";
	cout << sign(y[2]) << " Lb3 ";
	cout << sign(y[4]) << " Lb5 ";
	cout << sign(y[5]) << " Lb6 ";
	cout << sign(dpi(y[0],y[1],y[2],y[3],y[4],y[5])) << " Lbd ";
	cout << " < " << -fb(y[0],y[1],y[2],y[3],y[4],y[5]);
	cout << endl;
	}

void outC(double y[6])
	{
	cout << " cC - t ";
	cout << sign(y[0]) << " Lc1 ";
	cout << sign(y[1]) << " Lc2 ";
	cout << sign(y[2]) << " Lc3 ";
	cout << sign(y[4]) << " Lc5 ";
	cout << sign(y[5]) << " Lc6 ";
	cout << sign(dpi(y[0],y[1],y[2],y[3],y[4],y[5])) << " Lcd ";
	cout << " < " << -fc(y[0],y[1],y[2],y[3],y[4],y[5]);
	cout << endl;
	}

void outD(double y[6])
	{
	cout << " cD - t ";
	cout << sign(y[0]) << " Ld1 ";
	cout << sign(y[1]) << " Ld2 ";
	cout << sign(y[2]) << " Ld3 ";
	cout << sign(y[4]) << " Ld5 ";
	cout << sign(y[5]) << " Ld6 ";
	cout << sign(dpi(y[0],y[1],y[2],y[3],y[4],y[5])) << " Ldd ";
	cout << " < " << -fd(y[0],y[1],y[2],y[3],y[4],y[5]);
	cout << endl;
	}

void out(double y[6])
	{ outA(y); outB(y); outC(y); outD(y); }

void constraint()
	{
	const int size=25;
	char c[size][100]= {
		"cA+cB+cC+cD>0",
		"Lad+Lbd+Lcd+Ldd=0",
		"La1+Lb1+Lc1+Ld1>0",
		"La3+Lb2>0",
		"Lb3+Lc2>0",
		"Lc3+Ld2>0",
		"Ld3+La2>0",

		"La5+Lb6>0",
		"Lb5+Lc6>0",
		"Lc5+Ld6>0",
		"Ld5+La6>0"};

	for (int i=0;i<size;i++)
		cout << c[i] << endl;
	}

void bound()
	{
	int i;
	const int size = 29;
 	char c[size][5]={
		"cA","La1","La2","La3","La5","La6","Lad",
		"cB","Lb1","Lb2","Lb3","Lb5","Lb6","Lbd",
		"cC","Lc1","Lc2","Lc3","Lc5","Lc6","Lcd",
		"cD","Ld1","Ld2","Ld3","Ld5","Ld6","Ldd",
		"t"};
	cout << "BOUNDS\n\n";
	for (i=0;i<size;i++)
		cout << " -50 < " << c[i] << " < 50 " << endl;
	cout << "\n\nEND\n" << flush;
	}


int compression( double y[6])
	{
	return
		((radf(y[0],y[1],y[5])< 1.41422)&&
		   (radf(y[0],y[2],y[4])< 1.41422));
	}

int vortype( double y[6])
	{
	return
		!((radf(y[0],y[1],y[5])< 1.4142)&&
		   (radf(y[0],y[2],y[4])< 1.4142));
	}

int main()
	{
	cout << "MINIMIZE\n t\n\nST\n\n";

	double y[6]={2,2,2,2.51,2,2};
	double u[6]={2.25,2.25,2.25,2.65,2.25,2.25};
	for (int i0=0;i0<3;i0++)
	for (int i1=0;i1<3;i1++)
	for (int i2=0;i2<3;i2++)
	for (int i3=0;i3<3;i3++)
	for (int i4=0;i4<3;i4++)
	for (int i5=0;i5<3;i5++)
		{
		y[0]= 2.51 + (u[0]-2.51)*double(i0)/2.0;
		y[1]= 2.0 + (u[1]-2.0)*double(i1)/2.0;
		y[2]= 2.0 + (u[2]-2.0)*double(i2)/2.0;
		y[3]= 2.0 + (u[3]-2.0)*double(i3)/2.0;
		y[4]= 2.0 + (u[4]-2.0)*double(i4)/2.0;
		y[5]= 2.0 + (u[5]-2.0)*double(i5)/2.0;
		if (compression(y))
			out(y);
		}
	constraint();
	bound();
	}
