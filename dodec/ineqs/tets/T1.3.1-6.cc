#include<iostream.h>
#include<iomanip.h>
#include<fstream.h>
#include "error.h"
#include "interval.h"
#include "lineInterval.h"
#include "secondDerive.h"
#include "taylorInterval.h"
#include "recurse.h"

typedef char str[31];

int proveTet1(str,str,str,str,str,str,str,str,str);
int proveTet2(str,str,str,str,str,str,str,str,str);
int proveTet3(str,str,str,str,str,str,str,str,str);
int proveTet4(str,str,str,str,str,str,str,str,str);
int proveTet5(str,str,str,str,str,str,str,str,str);

void selfTest()
    {
    interMath::selfTest();
    linearization::selfTest();
    secondDerive::selfTest();
    taylorFunction::selfTest();
    }

static void setqr(domain& x,domain& z)
	{
	interval tx[6]={"4","4","4","4","4","4"};
	interval tz[6]={"6.334368540005047286","6.334368540005047286","6.334368540005047286","6.334368540005047286","6.334368540005047286","6.334368540005047286"};
	x = domain::lowerD(tx);
	z = domain::upperD(tz);
	}

int main()
{
 
  int i;
  const int NUM_TET_EQS=1;
  str tets[NUM_TET_EQS][9]={
  
    /* 1.3.1. */    {"0","-1","-0.245","-0.245","-0.245","0.063",
		     "0.063","0.063","1.6432"},
    /* 1.3.2. */    {"0","-1","-0.3798","-0.3798","-0.3798","0.198",
		     "0.198","0.198","1.642"},
    /* 1.3.3. */    {"0","1","0.151","0.151","0.151","-0.323","-0.323",
		     "-0.323","0.4807"},
    /* 1.3.4. */    {"1","0.42755","0.0392","0.0392","0.0392","0.0101",
		     "0.0101","0.0101","-0.2958"},
    /* 1.3.5. */    {"1","0","-0.107","-0.107","-0.107","0.116","0.116",
		     "0.116","0.1817"},  
    /* 1.3.6. */    {"1","0","-0.0623","-0.0623","-0.0623","0.0722",
		     "0.0722","0.0722","0.1763"}

                            };

  for(i=0;i<NUM_TET_EQS;i++)
    {
      
      if( proveTet1(tets[i][0],tets[i][1],tets[i][2],tets[i][3],
		    tets[i][4],tets[i][5],tets[i][6],tets[i][7],
		    tets[i][8]) )
	outFile<<"Tet Equation # "<<i+1<<" part 1 PASSED!!!"<<endl;
      else outFile<<"Tet Equation # "<<i+1<<" part 1 FAILED!!!"<<endl;
     
      if( proveTet2(tets[i][0],tets[i][1],tets[i][2],tets[i][3],
		    tets[i][4],tets[i][5],tets[i][6],tets[i][7],
		    tets[i][8]) )
	outFile<<"Tet Equation # "<<i+1<<" part 2 PASSED!!!"<<endl;
      else outFile<<"Tet Equation # "<<i+1<<" part 2 FAILED!!!"<<endl;

      if( proveTet3(tets[i][0],tets[i][1],tets[i][2],tets[i][3],
		    tets[i][4],tets[i][5],tets[i][6],tets[i][7],
		    tets[i][8]))
	outFile<<"Tet Equation # "<<i+1<<" part 3 PASSED!!!"<<endl;
      else outFile<<"Tet Equation # "<<i+1<<" part 3 FAILED!!!"<<endl;

      if( proveTet4(tets[i][0],tets[i][1],tets[i][2],tets[i][3],
		    tets[i][4],tets[i][5],tets[i][6],tets[i][7],
		    tets[i][8]) )
	outFile<<"Tet Equation # "<<i+1<<" part 4 PASSED!!!"<<endl;
      else outFile<<"Tet Equation # "<<i+1<<" part 4 FAILED!!!"<<endl;

      if( proveTet5(tets[i][0],tets[i][1],tets[i][2],tets[i][3],
		    tets[i][4],tets[i][5],tets[i][6],tets[i][7],
		    tets[i][8]) )
	outFile<<"Tet Equation # "<<i+1<<" part 5 PASSED!!!"<<endl;
      else outFile<<"Tet Equation # "<<i+1<<" part 5 FAILED!!!"<<endl;

    }

  return 0;
}
//{"0","-1","-0.245","-0.245","-0.245","0.063","0.063","0.063","1.6432"},
int proveTet1(str sigma,str sol,str y1,str y2,str y3,str y4,str y5,
	      str y6,str num)
{
  cellOption opt;
  domain x,z,x0,z0;
  setqr(x0,z0);
  setqr(x,z);
  taylorFunction truncSq=taylorSimplex::unit*"1.583592135001261822";
  taylorFunction two=taylorSimplex::unit*"2.0";
 
  taylorFunction F=taylorSimplex::voronoiVolume*"-1"*sigma+
                   taylorSimplex::sol*sol+taylorSimplex::y1*y1+
                   taylorSimplex::y2*y2+taylorSimplex::y3*y3+
                   taylorSimplex::y4*y4+
                   taylorSimplex::y5*y5+taylorSimplex::y6*y6+
                   taylorSimplex::unit*num;



  
  taylorFunction G=taylorSimplex::eta2_126*"-1"+two;
  taylorFunction H=taylorSimplex::eta2_135*"-1"+two;
  taylorFunction I=taylorSimplex::eta2_234*"-1"+two;
  taylorFunction J=taylorSimplex::eta2_456*"-1"+two;
 
  const taylorFunction* K[5]={&F,&G,&H,&I,&J};
  
  
  
  G.setReducibleState(1);
  H.setReducibleState(1);
  I.setReducibleState(1);
  J.setReducibleState(1);

  return prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt);
}

int proveTet2(str sigma,str sol,str y1,str y2,str y3,str y4,str y5,
	      str y6,str num)
{
  cellOption opt;
  domain x,z,x0,z0;
  setqr(x0,z0);
  setqr(x,z);
  taylorFunction truncSq=taylorSimplex::unit*"1.583592135001261822";
  taylorFunction two=taylorSimplex::unit*"2.0";

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"*sigma+
                   taylorSimplex::sol*sol+taylorSimplex::y1*y1+
                   taylorSimplex::y2*y2+taylorSimplex::y3*y3+
                   taylorSimplex::y4*y4+
                   taylorSimplex::y5*y5+taylorSimplex::y6*y6+
                   taylorSimplex::unit*num;
   
  taylorFunction G=taylorSimplex::eta2_126+two*"-1";  
  const taylorFunction* K[2]={&F,&G};  
  
  return prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt);
}

int proveTet3(str sigma,str sol,str y1,str y2,str y3,str y4,str y5,
	      str y6,str num)
{
  cellOption opt;
  domain x,z,x0,z0;
  setqr(x0,z0);
  setqr(x,z);
  taylorFunction truncSq=taylorSimplex::unit*"1.583592135001261822";
  taylorFunction two=taylorSimplex::unit*"2.0";

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"*sigma+
                   taylorSimplex::sol*sol+taylorSimplex::y1*y1+
                   taylorSimplex::y2*y2+taylorSimplex::y3*y3+
                   taylorSimplex::y4*y4+
                   taylorSimplex::y5*y5+taylorSimplex::y6*y6+
                   taylorSimplex::unit*num;
  
  
  taylorFunction H=taylorSimplex::eta2_135+two*"-1";
  
  const taylorFunction* K[2]={&F,&H};
  
  
  
  return prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt);
}

int proveTet4(str sigma,str sol,str y1,str y2,str y3,str y4,str y5,
	      str y6,str num)
{
  cellOption opt;
  domain x,z,x0,z0;
  setqr(x0,z0);
  setqr(x,z);
  taylorFunction truncSq=taylorSimplex::unit*"1.583592135001261822";
  taylorFunction two=taylorSimplex::unit*"2.0";

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"*sigma+
                   taylorSimplex::sol*sol+taylorSimplex::y1*y1+
                   taylorSimplex::y2*y2+taylorSimplex::y3*y3+
                   taylorSimplex::y4*y4+
                   taylorSimplex::y5*y5+taylorSimplex::y6*y6+
                   taylorSimplex::unit*num;
  
  
  taylorFunction I=taylorSimplex::eta2_234+two*"-1";
  
  const taylorFunction* K[2]={&F,&I};
  
  
  
  return prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt);
}

int proveTet5(str sigma,str sol,str y1,str y2,str y3,str y4,str y5,
	      str y6,str num)
{
  cellOption opt;
  domain x,z,x0,z0;
  setqr(x0,z0);
  setqr(x,z);
  taylorFunction truncSq=taylorSimplex::unit*"1.583592135001261822";
  taylorFunction two=taylorSimplex::unit*"2.0";

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"*sigma+
                   taylorSimplex::sol*sol+taylorSimplex::y1*y1+
                   taylorSimplex::y2*y2+taylorSimplex::y3*y3+
                   taylorSimplex::y4*y4+
                   taylorSimplex::y5*y5+taylorSimplex::y6*y6+
                   taylorSimplex::unit*num;
  
  
  taylorFunction J=taylorSimplex::eta2_456+two*"-1";
  const taylorFunction* K[2]={&F,&J};
  
  
  F.setReducibleState(1);
  return prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt);
}







