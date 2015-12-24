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

int proveTet1(str,str,str);
int proveTet2(str,str,str);
int proveTet3(str,str,str);
int proveTet4(str,str,str);
int proveTet5(str,str,str);

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
	interval tz[6]={"6.33436854000504728617","6.33436854000504728617","6.33436854000504728617","6.33436854000504728617","6.33436854000504728617","6.33436854000504728617"};
	x = domain::lowerD(tx);
	z = domain::upperD(tz);
	}

int main()
{
 
  int i;
  const int NUM_TET_EQS=53;
  str tets[NUM_TET_EQS][3]={
                             {"0.68","-1.88718","1.545510"},
                             {"0.68","-0.90746","0.706725"},
    		             {"0.68","-0.46654","0.329233"},
                             {"0.55889","0","-0.073648"},
                             {"0.63214","0","-0.13034"},
    	                     {"0.73256","0","-0.23591"},
                             {"0.89346","0","-0.40505"},
                             {"0.3","0.5734","-0.978221"},
                             {"0.3","0.03668","0.024767"},
                  /* 10 */   {"0.3","-0.04165","0.121199"},
			     {"0.3","-0.1234","0.209279"},
			     {"0.42755","0.11509","-0.171859"},
			     {"0.42755","0.04078","-0.050713"},
			     {"0.42755","-0.11031","0.135633"},
			     {"0.42755","-0.13091","0.157363"},
			     {"0.55792","0.21394","-0.417998"},
			     {"0.55792","0.0068","-0.081902"},
			     {"0.55792","-0.0184","-0.051224"},
			     {"0.55792","-0.24335","0.193993"},
		  /* 20 */   {"0.68","0.30651","-0.648496"},
			     {"0.68","0.06965","-0.27800"},
			     {"0.68","-0.0172","-0.15662"},
			     {"0.68","-0.41812","0.287778"},
			     {"0.64934","0","-0.14843"},
			     {" 0.6196","0","-0.11800"},
			     {"0.58402","0","-0.090290"},
			     {"0.25181","0","0.096509"},
			     {"0.00909","0","0.199559"},
		             {"-0.93877","0","0.537892"},
		   /* 30 */  {"-0.93877","0.20211","0.27313"},
			     {"-0.93877","-0.63517","1.20578"},
			     {"-1.93877","0","0.854804"},
			     {"-1.93877","0.20211","0.621886"},
			     {"-1.93877","-0.63517","1.57648"},
			     {"0.42775","0","-0.000111"},
			     {"0.55792","0","-0.073037"},
			     {"0","0.07853","0.08865"},
			     {"0","0.00339","0.198693"},
		             {"0","-0.18199","0.396670"},
	          /* 40 */   {"0.42755","0.20000","-0.332061"},
			     {"0.3","0.36373","-0.582630"},
			     {"0.3","-0.20583","0.279851"},
			     {"0.3","-0.40035","0.446389"},
			     {"0.3","-0.83259","0.816450"},
			     {"0.42755","0.51838","-0.932759"},
			     {"0.42755","-0.29344","0.296513"},
			     {"0.42755","-0.57056","0.533768"},
			     {"0.42755","-1.18656","1.06115"},
			     {"0.55792","0.67644","-1.29062"},
		   /* 50 */  {"0.55792","-0.38278","0.313365"},
			     {"0.55792","-0.74454","0.623085"},
			     {"0.55792","-1.54837","1.31128"},
			     {"0.68","0.82445","-1.62571"}
                            };

  for(i=0;i<NUM_TET_EQS;i++)
    {
      
      if( proveTet1(tets[i][0],tets[i][1],tets[i][2]) )
	outFile<<"Tet Equation # "<<i+1<<" part 1 PASSED!!!"<<endl;
      else outFile<<"Tet Equation # "<<i+1<<" part 1 FAILED!!!"<<endl;
     
      if( proveTet2(tets[i][0],tets[i][1],tets[i][2]) )
	outFile<<"Tet Equation # "<<i+1<<" part 2 PASSED!!!"<<endl;
      else outFile<<"Tet Equation # "<<i+1<<" part 2 FAILED!!!"<<endl;

      if( proveTet3(tets[i][0],tets[i][1],tets[i][2]) )
	outFile<<"Tet Equation # "<<i+1<<" part 3 PASSED!!!"<<endl;
      else outFile<<"Tet Equation # "<<i+1<<" part 3 FAILED!!!"<<endl;

      if( proveTet4(tets[i][0],tets[i][1],tets[i][2]) )
	outFile<<"Tet Equation # "<<i+1<<" part 4 PASSED!!!"<<endl;
      else outFile<<"Tet Equation # "<<i+1<<" part 4 FAILED!!!"<<endl;

      if( proveTet5(tets[i][0],tets[i][1],tets[i][2]) )
	outFile<<"Tet Equation # "<<i+1<<" part 5 PASSED!!!"<<endl;
      else outFile<<"Tet Equation # "<<i+1<<" part 5 FAILED!!!"<<endl;

    }

  return 0;
}

int proveTet1(str sol,str dih,str num)
{
  cellOption opt;
  domain x,z,x0,z0;
  setqr(x0,z0);
  setqr(x,z);
  taylorFunction truncSq=taylorSimplex::unit*"1.583592135001261822";
  taylorFunction two=taylorSimplex::unit*"2.0";
 
  taylorFunction F=taylorSimplex::voronoiVolume*"-1"+taylorSimplex::sol*sol+taylorSimplex::dih*dih+taylorSimplex::unit*num;
  
  taylorFunction G=taylorSimplex::eta2_126*"-1"+two;
  taylorFunction H=taylorSimplex::eta2_135*"-1"+two;
  taylorFunction I=taylorSimplex::eta2_234*"-1"+two;
  taylorFunction J=taylorSimplex::eta2_456*"-1"+two;
 
  const taylorFunction* K[5]={&F,&G,&H,&I,&J};
  
  
  F.setReducibleState(1);
  G.setReducibleState(1);
  H.setReducibleState(1);
  I.setReducibleState(1);
  J.setReducibleState(1);

  return prove::recursiveVerifier(0,x,z,x0,z0,K,5,opt);
}

int proveTet2(str sol,str dih,str num)
{
  cellOption opt;
  domain x,z,x0,z0;
  setqr(x0,z0);
  setqr(x,z);
  taylorFunction truncSq=taylorSimplex::unit*"1.583592135001261822";
  taylorFunction two=taylorSimplex::unit*"2.0";

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+taylorSimplex::sol*sol+taylorSimplex::dih*dih+taylorSimplex::unit*num;
  
  taylorFunction G=taylorSimplex::eta2_126+two*"-1";
  
  const taylorFunction* K[2]={&F,&G};
  
  
  F.setReducibleState(1);
  return prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt);
}

int proveTet3(str sol,str dih,str num)
{
  cellOption opt;
  domain x,z,x0,z0;
  setqr(x0,z0);
  setqr(x,z);
  taylorFunction truncSq=taylorSimplex::unit*"1.583592135001261822";
  taylorFunction two=taylorSimplex::unit*"2.0";

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+taylorSimplex::sol*sol+taylorSimplex::dih*dih+taylorSimplex::unit*num;
  
  
  taylorFunction H=taylorSimplex::eta2_135+two*"-1";
  
  const taylorFunction* K[2]={&F,&H};
  
  
  F.setReducibleState(1);
  return prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt);
}

int proveTet4(str sol,str dih,str num)
{
  cellOption opt;
  domain x,z,x0,z0;
  setqr(x0,z0);
  setqr(x,z);
  taylorFunction truncSq=taylorSimplex::unit*"1.583592135001261822";
  taylorFunction two=taylorSimplex::unit*"2.0";

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+taylorSimplex::sol*sol+taylorSimplex::dih*dih+taylorSimplex::unit*num;
  
  
  taylorFunction I=taylorSimplex::eta2_234+two*"-1";
  
  const taylorFunction* K[2]={&F,&I};
  
  
  F.setReducibleState(1);
  return prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt);
}

int proveTet5(str sol,str dih,str num)
{
  cellOption opt;
  domain x,z,x0,z0;
  setqr(x0,z0);
  setqr(x,z);
  taylorFunction truncSq=taylorSimplex::unit*"1.583592135001261822";
  taylorFunction two=taylorSimplex::unit*"2.0";

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+taylorSimplex::sol*sol+taylorSimplex::dih*dih+taylorSimplex::unit*num;
  
  
  taylorFunction J=taylorSimplex::eta2_456+two*"-1";
  const taylorFunction* K[2]={&F,&J};
  
  
  F.setReducibleState(1);
  return prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt);
}







