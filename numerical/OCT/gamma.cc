#include <iomanip.h>
#include "constants.h"
#include "octeqns.h"
#include <math.h>

/* This is a fully functional version, that gave data for
three of the four inequalities in 4.2.  I'm about to redo
it with three purposes in mind.  First to make use of the
symmetries to reduce the number of cases.  Second to tweak
things so that all four inequalities go through. Third,
I would like the output to produce C++ code that can be
compiled directly, rather than the output data is currently
does in gamma.out.

The data for the three cases can be found in gamma.out.

Each inequality in this version involves 2^4 cases depending
on whether each face is small or not.  Not all these cases
are necessary.
*/

static int printIfSeparate
	(const double ymin[13],const double ymax[13],
	const double ineqCoeff[2], const int etaIndex[4])
    {
	static int count=0;
	int i;
    double CC[29];
	double coeff[2]={ineqCoeff[0],ineqCoeff[1]};
	double yminX[13],ymaxX[13];
	for (i=0;i<13;i++) { yminX[i]=ymin[i]; ymaxX[i]=ymax[i]; }
	int index[4]={etaIndex[0],etaIndex[1],etaIndex[2],etaIndex[3]};
    separateChecked(yminX,ymaxX,CC,coeff,index); 
	if (CC[28]>=0) return 0;

	// report answer ...
	cout << "\nf"<<count<<"[2]= {" << ineqCoeff[0] << "," << ineqCoeff[1] << "};\n";
	cout << "T"<<count<<"[4]= {";
	for (i=0;i<4;i++) cout << etaIndex[i] << (i==3 ? "};\n" : ",");
	cout << "Dm"<<count<<"[13] = {"; 
	for (i=0;i<13;i++) cout << ymin[i] << (i==12 ? "};\n" : ",");
	cout << "DM"<<count<<"[13] = {"; 
	for (i=0;i<13;i++) cout << ymax[i] << (i==12 ? "};\n" : ",");
	cout << "C"<<count<<"[29]= {";
	for (i=0;i<29;i++) cout << CC[i] << (i==28 ? "};\n" : ", ") 
		<< ((i+1)%7 ? "" : "\n   " ); 
	cout << flush;
	count++;
	return 1; 
    }

void pickSubIndex(const double ymin[13],const double ymax[13],
	int subIndex[3])
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
	subIndex[0]=pos0; subIndex[1]=pos1; subIndex[2]=pos2;
	return;
	}

void setz(const double ymin[13],const double ymax[13],
	double zmin[13],double zmax[13],const int subIndex[3],
	int binaryCode/*<8*/)
	{
	int i;
	double w[3];
	for (i=0;i<13;i++) { zmin[i]=ymin[i]; zmax[i]=ymax[i]; }
	for (i=0;i<3;i++) 
		{
		w[i]= (ymax[subIndex[i]]-ymin[subIndex[i]])/2.0;
		zmax[subIndex[i]] -= w[i];
		}
	int k=1;
	for (i=0;i<3;i++) 
		{
		if (1 == ((binaryCode/k) % 2)) { zmin[subIndex[i]] += w[i]; zmax[subIndex[i]] += w[i]; }
		k = 2*k;
		}
	}

void subdivide(const double ymin[13],const double ymax[13],
	const double ineqCoeff[2],const int etaIndex[4])
    {
	int i;
	double zmin[13], zmax[13];
	int subIndex[3];
	int binaryCode;
	pickSubIndex(ymin,ymax,subIndex);
	cout << "\n// NewDivision: ";
	for (i=0;i<3;i++) cout << subIndex[i] << " "; cout << endl;
	for (binaryCode=0;binaryCode<8;binaryCode++) 
		{
		setz(ymin,ymax,zmin,zmax,subIndex,binaryCode);
		if (!printIfSeparate(zmin,zmax,ineqCoeff,etaIndex))
			{
			subdivide(zmin,zmax,ineqCoeff,etaIndex);
			}
		}
	cout << "//// EndDivision: ";
	for (i=0;i<3;i++) cout << subIndex[i] << " "; cout << endl;
    }

oldmain()
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
	for (i=0;i<4;i++) // BUGBUG was i=0;i<4;i++
		{
		cout << "\n\n//********* STARTING NEW INEQUALITY *************\n";
		cout << "// parameters = {" << f[i][0] << "," << f[i][1] << "}\n";
		for (int i0=0;i0<2;i0++)
		for (int i1=i0;i1<2;i1++) // use symmetry.
		for (int i2=0;i2<2;i2++)
		for (int i3=0;i3<2;i3++)
			{
			int etaIndex[4] = {i0,i1,i2,i3};
			double ineqCoeff[2] = {f[i][0],f[i][1]};
			subdivide(ymin,ymax,ineqCoeff,etaIndex);
			}
		}
	return 1;
	}

////////////////////////////////////////////////////////////
// Everything that follows is here because there were 18 cases that
// didn't not go through easily on the interval verifications.
// So I have come back to subdivide them further.
/* This is defined in part3oData.cc. The file part3oData.cc is a C++-file
obtained by converting the output data of oldmain() to C++-format.
 setFTC is defined in part3oData.cc.  It reads one record.  */
int setFTC(double f[2],int T[4],double Dm[13],double DM[13],double C[29],
    int index);

// copied from nint::part3o.cc.
static int Casenum(int index)
    {
    if (index<96) { return  0; }
    else if (index <311) { return 1; }
    else if (index <1037) { return 2; }
    else return 3;
    }

// copied.
static int checkIntegrity(double f[2],int casenum)
    {
    static const double fvalue[4][2]=
        {{-3.0508,9.494},{-0.27605,1.0472},
        {0.198867,-0.7624},{0.844,-3.5926}};
    if (fabs(f[0]-fvalue[casenum][0])> 1.0e-14) return 0;
    if (fabs(f[1]-fvalue[casenum][1])> 1.0e-14) return 0;
    return 1;
    }

// copied.
static int nearest(double x)
    {
    return (int)(x< 0 ? ceil(x-0.5) : floor (x+0.5));
    }
 


// adapted from part3o.cc:setInterval
 
    //////////
    // There has been rounding in the file part3oData.cc.
    // We need to reconstruct the intervals from the
    // rounded data.  
    // ym and yM are equal to cleaned up versions of Dm and DM.
    //
static void setEndpoint(const double Dm[13],const double DM[13],
    double ym[13],double yM[13])
    {
    static const double t51 = 2.51;
    static const double sqrt8 = 2.8284271247461900976033774;
    int i;
    double x0[13]={t51,2,2,2,2, 2,2,2,2,2, 2,2,2};
    double z0[13];
    z0[0]=sqrt8; for (i=1;i<13;i++) z0[i]=t51;
 
    int segment[13];
    for (i=0;i<13;i++)  { segment[i] = nearest((z0[i]-x0[i])/(DM[i]-Dm[i])); }
 
    int segN[13];
    for (i=0;i<13;i++)
        { segN[i]= nearest (((Dm[i]-x0[i])*segment[i])/(z0[i]-x0[i])); }
    int segU[13];
    for (i=0;i<13;i++)
        { segU[i]= nearest (((DM[i]-x0[i])*segment[i])/(z0[i]-x0[i])); }
 
    for (i=0;i<13;i++)
        {
        double s=segment[i];
        double n=segN[i];
        double u=segU[i];
        ym[i]= x0[i]+ (n/s)*(z0[i]-x0[i]);
        yM[i]= x0[i]+ (u/s)*(z0[i]-x0[i]);
        }
    }

// adapted from part3o.cc:testSetInterval.
static void testSetEndpoints(int index)
    {
    double f[2];int T[4];double Dm[13];double DM[13];double C[29];
    setFTC(f,T,Dm,DM,C,index);
    double ym[13], yM[13];
    setEndpoint(Dm,DM,ym,yM);
    int i;
	cout.precision(16);
    for (i=0;i<13;i++)
     cout <<Dm[i]<< " " << ym[i] << endl; cout << endl;
    for (i=0;i<13;i++)
     cout <<DM[i]<< " " << yM[i] << endl; cout << endl;
    }
 



main()
	{
	// The following cases are not going through rigorously with intervals.
	// Let's rerun them.
	int Retry[18]={355,617,688,690,703,705,
				   707,709,761,765,837,859,
				   882,899,969,972,1000,1008};
	double ineqCoeff[2];
	int etaIndex[4];
	double Dm[13],DM[13],ym[13],yM[13];
	double C[29];
	cout.precision(16);
	for (int i=0;i<18;i++)
		{
		setFTC(ineqCoeff,etaIndex,Dm,DM,C,Retry[i]);
		cout << "\n\n//********* STARTING NEW INEQUALITY(" << Retry[i] 
				<< ") *************\n";
		cout << "// parameters = {" << ineqCoeff[0] 
				<< "," << ineqCoeff[1] << "}\n";
		setEndpoint(Dm,DM,ym,yM); // beautify endpoints.
		subdivide(ym,yM,ineqCoeff,etaIndex);
		}
	}

