#include <iomanip.h>
#include "constants.h"
#include "octeqns.h"

double fTest[2];
int tType[4];
static int test(double ymin[13],double ymax[13])
    {
    double CC[29];
    separateChecked(ymin,ymax,CC,fTest,tType); 
	int i;
	if (CC[28]>=0) return 0;

	// report answer ...
	cout << "\nf= {" << fTest[0] << "," << fTest[1] << "}\n";
	cout << "T= {";
	for (i=0;i<4;i++) cout << tType[i] << (i==3 ? "}\n" : ",");
	cout << "Dm = {"; 
	for (i=0;i<13;i++) cout << ymin[i] << (i==12 ? "}\n" : ",");
	cout << "DM = {"; 
	for (i=0;i<13;i++) cout << ymax[i] << (i==12 ? "}\n" : ",");
	cout << "C= {";
	for (i=0;i<29;i++) cout << CC[i] << (i==28 ? "}\n" : ", ") 
		<< ((i+1)%7 ? "" : "\n   " ); 
	cout << flush;
	return 1; 
    }

void getABC(double ymin[13],double ymax[13],int r[3])
	{
	double fat0=-1;
	double fat1=-1;
	double fat2=-1;
	double temp =-1;
	int i;
	int pos0=-1,pos1=-1,pos2=-1;
	for (i=0;i<13;i++) 
		if (ymax[i]-ymin[i] > fat0)
			{
			fat0=ymax[i]-ymin[i]; pos0=i;
			}
	for (i=0;i<13;i++)
		if ((ymax[i]-ymin[i] > fat1)&&(i!=pos0))
			{
			fat1=ymax[i]-ymin[i]; pos1=i;
			}
	for (i=0;i<13;i++)
		if ((ymax[i]-ymin[i] > fat2)&&(i!=pos0)&&(i!=pos1))
			{
			fat2=ymax[i]-ymin[i]; pos2=i;
			}
	r[0]=pos0; r[1]=pos1; r[2]=pos2;
	return;
	}

void setABC(double ymin[13],double ymax[13],double zmin[13],double zmax[13],
		int r[3],int c)
	{
	int i;
	double w[3];
	for (i=0;i<13;i++) { zmin[i]=ymin[i]; zmax[i]=ymax[i]; }
	for (i=0;i<3;i++) 
		{
		w[i]= (ymax[r[i]]-ymin[r[i]])/2.0;
		zmax[r[i]] -= w[i];
		}
	int k=1;
	for (i=0;i<3;i++) 
		{
		if (1 == ((c/k) % 2)) { zmin[r[i]] += w[i]; zmax[r[i]] += w[i]; }
		k = 2*k;
		}
	}

void subdivide(double ymin[13],double ymax[13])
    {
	int i;
	double zmin[13], zmax[13];
	int r[3];
	int c;
	getABC(ymin,ymax,r);
	cout << "\n\\\\ NewDivision: ";
	for (i=0;i<3;i++) cout << r[i] << " "; cout << endl;
	for (c=0;c<8;c++) 
		{
		setABC(ymin,ymax,zmin,zmax,r,c);
		//cout << " "; for (i=0;i<13;i++) cout << zmin[i] << " "; cout << endl;
		//cout << " ";for (i=0;i<13;i++) cout << zmax[i] << " "; cout << endl;
		if (!test(zmin,zmax))
			{
			subdivide(zmin,zmax);
			}
		}
	cout << "\\\\ EndDivision: ";
	for (i=0;i<3;i++) cout << r[i] << " "; cout << endl;
    }

main()
	{
    double f[4][2] = {{-3.0508,9.494},
					 {-0.27605,1.0472},
					 {0.198867,-0.7624},
					 {0.844,-3.5926}};
    double ymin[13];
    double ymax[13];
    int i;
    for (i=0;i<13;i++) { ymin[i]=2; ymax[i]=2.51; }
    ymin[0]=2.51; ymax[0]=global::sqrt8;
	for (i=2;i<4;i++) // BUGBUG was i=0;...
		{
		cout << "\n\n\\********* STARTING NEW INEQUALITY *************\n";
		cout << "\\parameters = {" << f[i][0] << "," << f[i][1] << "}\n";
		fTest[0]= f[i][0]; fTest[1]=f[i][1];
		for (int i0=0;i0<2;i0++)
		for (int i1=0;i1<2;i1++)
		for (int i2=0;i2<2;i2++)
		for (int i3=0;i3<2;i3++)
			{
			if ((i==2)&&(i0+i1+i2==0)) continue; // BUGBUG. pick up..
			tType[0]=i0; tType[1]=i1; tType[2]=i2; tType[3]=i3;
			subdivide(ymin,ymax);
			}
		}
	}

/*



*/
