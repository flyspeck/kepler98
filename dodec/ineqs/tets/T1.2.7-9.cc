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

int proveTet1(str,str,str,str,str,str,str,str);


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
  str tets[NUM_TET_EQS][8]={

    /* 1.3.7. */  {"-1","0.237","-0.372","-0.372","0.708","-0.372",
		   "-0.372","2.3169"},
    /* 1.3.8. */  {"-1","0.237","-0.363","-0.363","0.688","-0.363",
		  "-0.363","2.2849"},
    /* 1.3.9. */  {"1","-0.505","0.152","0.152","-0.766","0.152",
		   "0.152","0.09503"},
                            };

  for(i=0;i<NUM_TET_EQS;i++)
    {
      
      if( proveTet1(tets[i][0],tets[i][1],tets[i][2],tets[i][3],
		    tets[i][4],tets[i][5],tets[i][6],tets[i][7]))
	outFile<<"Tet Equation # "<<i+1<<" part 1 PASSED!!!"<<endl;
      else outFile<<"Tet Equation # "<<i+1<<" part 1 FAILED!!!"<<endl;
     
      

    }

  return 0;
}

int proveTet1(str dih,str y1,str y2,str y3,str y4,str y5,
	      str y6,str num)
{
  cellOption opt;
  domain x,z,x0,z0;
  setqr(x0,z0);
  setqr(x,z);
  taylorFunction truncSq=taylorSimplex::unit*"1.583592135001261822";
  taylorFunction two=taylorSimplex::unit*"2.0";
 
  taylorFunction F=taylorSimplex::dih*dih+
                   taylorSimplex::y1*y1+
                   taylorSimplex::y2*y2+taylorSimplex::y3*y3+
                   taylorSimplex::y4*y4+
                   taylorSimplex::y5*y5+taylorSimplex::y6*y6+
                   taylorSimplex::unit*num;



  
  
 
  const taylorFunction* K[1]={&F};
  
  
  


  return prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt);
}
