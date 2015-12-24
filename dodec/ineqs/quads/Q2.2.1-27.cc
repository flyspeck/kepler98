#include "error.h"
#include <iostream.h>
#include<fstream.h>
#include <iomanip.h>
#include "interval.h"
#include "lineInterval.h"
#include "secondDerive.h"
#include "taylorInterval.h"
#include "recurse.h"

// stuff for truncation....

void selfTest()
    {
    interMath::selfTest();
    linearization::selfTest();
    secondDerive::selfTest();
    taylorFunction::selfTest();
    }

	////////
	// This runs the asymmetric verification on a quad cluster
	// with the short edge drawn.
	//
int openQ(const taylorFunction& FA,const taylorFunction& FB,
	cellOption opt)
	{
	domain xA,xB,zA,zB;
	interval x[6]={"4","4","4","6.3343685400050472861798","4","4"};
	xB = xA = domain::lowerD(x);
	interval v("6.3343685400050472861798");
	interval z[6]={v,v,v,"12.668737080010094572359",v,v};
	zB = zA = domain::upperD(z);
	if (!FA.getReducibleState() ||
		!FB.getReducibleState())
		{
		error::message(
		 "Program is unlikely to terminate without reducibility");
		}
	const taylorFunction* IA[1]={&FA};
	const taylorFunction* IB[1]={&FB};
	prove::recursiveVerifierQ(0,xA,xB,zA,zB,IA,IB,1,opt);
	return 1;
	}

	////////
	// The verification of a "standard" quad cluster inequality
	// using truncatedVoronoiVolume.
	// truncatedVoronoiVolume > rhsC + rhsDIH dih + rhsSOL sol.
	//
int standardQ(interval rhsSOL,interval rhsDIH,interval rhsC)
	{
	static const interval zero("0");
	cout << "running standardQ(sol,dih,const)= " << rhsSOL << " "
		 << rhsDIH << " " << rhsC << endl << flush;
	/*asymmetric*/{
	taylorFunction FA= taylorSimplex::truncatedVoronoiVolume*("-1")
		+taylorSimplex::unit*(rhsC)
		+taylorSimplex::dih*(rhsDIH)
		+taylorSimplex::sol*(rhsSOL);
	taylorFunction FB= taylorSimplex::truncatedVoronoiVolume*("-1")
		+taylorSimplex::sol*(rhsSOL);
	FA.setReducibleState(1);
	FB.setReducibleState(1);
	cellOption opt;
	// an examination of the symmetries of the 18 cases gives...
	int skippable[7] = {7,8,9,10,12,14,16};
	opt.setSkipCases(skippable,7);
	openQ(FA,FB,opt);
	error::printTime("asymmetric phase complete ");
	}
 
	/*split*/{
	taylorFunction FA=
	         taylorSimplex::truncatedVoronoiVolume*("-1")
		+taylorSimplex::unit*(rhsC)
		+taylorSimplex::dih2*(rhsDIH)
		+taylorSimplex::sol*(rhsSOL);
	taylorFunction FB= taylorSimplex::truncatedVoronoiVolume*("-1")
		+taylorSimplex::dih2*(rhsDIH)
		+taylorSimplex::sol*(rhsSOL);
	FA.setReducibleState(1);
	FB.setReducibleState(1);
	cellOption opt;
	// an examination of the split symmetries of the 18 cases gives...
	int skippable[6] = {3,5,7,9,15,16};
	opt.setSkipCases(skippable,6);
	openQ(FA,FB,opt);
	error::printTime("split phase complete ");
	}
	return 1;
	}

typedef char str[31];

int main()
{
  ofstream outFile;
  outFile.open("quadChex27.dat");
  if(!outFile) {cout<<"Can't open output file."<<endl; return 1;}
 
  int i;
  const int NUM_QUAD_EQS=27;
  interval quads[NUM_QUAD_EQS][3]={
                                {"0.42775","0.15098","-0.3670"},
				{"0.42775","0.09098","-0.1737"},
				{"0.42775","0.00000","0.0310"},
				{"0.42775","-0.18519","0.3183"},
				{"0.42775","-0.20622","0.3438"},
				{"0.55792","0.30124","-1.0173"},
				{"0.55792","0.02921","-0.2101"},
				{"0.55792","0.00000","-0.1393"},
				{"0.55792","-0.05947","-0.0470"},
		      /* 10 */  {"0.55792","-0.39938","0.4305"},
				{"0.55792","-2.50210","2.8976"},
				{"0.68000","0.44194","-1.6264"},
				{"0.68000","0.10957","-0.6753"},
				{"0.68000","0.00000","-0.4029"},
				{"0.68000","-0.86096","0.8262"},
				{"0.68000","-2.44439","2.7002"},
				{"0.30000","0.12596","-0.1279"},
				{"0.30000","0.02576","0.1320"},
				{"0.30000","-0.00000","0.1945"},
		      /* 20 */  {"0.30000","-0.03700","0.2480"},
				{"0.30000","-0.22476","0.5111"},
				{"0.30000","-2.31852","2.9625"},
				{"0.00000","0.23227","-0.1042"},
				{"0.00000","-0.07448","0.5591"},
				{"0.00000","-0.22019","0.7627"},
				{"0.00000","-0.80927","1.5048"},
				{"0.00000","-5.84380","7.3468"}
                               };

  
  for(i=26;i<27;i++)
    {
      if( standardQ(quads[i][0],quads[i][1],quads[i][2]) ) outFile<<"quad"<<i+1<<" complete!"<<endl;
      else outFile<<"quad"<<i+1<<" FAILED!!!"<<endl;
    }

  return 0;
}


