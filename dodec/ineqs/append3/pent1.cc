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
int openQ(const taylorFunction& FA1,const taylorFunction& FA2,
	  const taylorFunction& FB1,const taylorFunction& FB2,
	cellOption opt)
	{
	
	  interval v("6.3343685400050472861798");
	  domain xA,xB,zA,zB;
	interval x1[6]={"4","4","4","6.3343685400050472861798","4","4"};
	interval x2[6]={"4","4","4","6.3343685400050472861798","4",v};
	xA = domain::lowerD(x1);xB = domain::lowerD(x2);
	
	interval z1[6]={"5.4289","4.870849",v,"14.645929",v,v};
	interval z2[6]={"5.1984","4.870849",v,"14.645929",v,"10.169721"};

	zA = domain::upperD(z1);zB = domain::upperD(z2);
	
	/*if (!FA.getReducibleState() ||
		!FB.getReducibleState())
		{
		error::message(
		 "Program is unlikely to terminate without reducibility");
		}*/
	const taylorFunction* IA[2]={&FA1,&FA2};
	const taylorFunction* IB[2]={&FB1,&FB2};
	prove::recursiveVerifierQ(0,xA,xB,zA,zB,IA,IB,2,opt);
	return 1;
	}


void standardQ()
	{
	static const interval zero("0");
	cout << "running standardQ(pent1)= "<<endl;
	/*asymmetric*/{
	taylorFunction FA1= taylorSimplex::truncatedVoronoiVolume*"-1"
	 	  +taylorSimplex::unit*"0.55";
		
	taylorFunction FA2= taylorSimplex::dih*"-1"
	 	  +taylorSimplex::unit*"1.917";

	taylorFunction FB1= taylorSimplex::truncatedVoronoiVolume*"-1";
	 
	taylorFunction FB2= taylorSimplex::unit*"0";
		
	FA1.setReducibleState(1);FA2.setReducibleState(1);
	FB1.setReducibleState(1);FB2.setReducibleState(1);
	

	cellOption opt;


	openQ(FA1,FA2,FB1,FB2,opt);
	error::printTime("asymmetric phase complete ");
	}
 
	/*split*/{
	
	 taylorFunction FA1= taylorSimplex::truncatedVoronoiVolume*"-1"
	 	  +taylorSimplex::unit*"0.55";
		
	taylorFunction FA2= taylorSimplex::dih*"-1"
	 	  +taylorSimplex::unit*"1.917";

	taylorFunction FB1= taylorSimplex::truncatedVoronoiVolume*"-1";
	 
	taylorFunction FB2= taylorSimplex::dih;

		
	FA1.setReducibleState(1);FA2.setReducibleState(1);
	FB1.setReducibleState(1);FB2.setReducibleState(1);

	cellOption opt;


	openQ(FA1,FA2,FB1,FB2,opt);
	error::printTime("split phase complete ");
	}
	return ;
	}

int main()
{
  standardQ();

  return 0;
}

