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


void standardQ()
	{
	static const interval zero("0");
	cout << "running standardQ(QuadMuLB)= "<<endl;
	/*asymmetric*/{
	taylorFunction FA= taylorSimplex::truncatedVoronoiVolume*"-1"
	  +taylorSimplex::sol*"0.42755"
	  +taylorSimplex::unit*"0.031350";
		
	taylorFunction FB= taylorSimplex::truncatedVoronoiVolume*"-1"+
	  taylorSimplex::sol*"0.42755";
		
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
	         taylorSimplex::truncatedVoronoiVolume*"-1"+
	  taylorSimplex::sol*"0.42755"
		+taylorSimplex::unit*"0.031350";
		
	taylorFunction FB= taylorSimplex::truncatedVoronoiVolume*"-1"+
	  taylorSimplex::sol*"0.42755";
		
	FA.setReducibleState(1);
	FB.setReducibleState(1);
	cellOption opt;
	// an examination of the split symmetries of the 18 cases gives...
	int skippable[6] = {3,5,7,9,15,16};
	opt.setSkipCases(skippable,6);
	openQ(FA,FB,opt);
	error::printTime("split phase complete ");
	}
	return ;
	}

int main()
{
  standardQ();

  return 0;
}

