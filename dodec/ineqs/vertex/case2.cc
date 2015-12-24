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

void proveTet1();
void proveTet2();
void proveTet3();
void proveTet4();
void proveTet5();

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
 
  proveTet1();
  proveTet2();
  proveTet3();
  proveTet4();
  proveTet5();
  return 0;
}


 
void proveTet1()
{
  cellOption opt;
  domain x,z,x0,z0;
  

  interval tx[6]={"4","4","4","4","4","4"};

  interval tz[6]={"4.2025","4.41",
		  "4.41","6.33436854000504728617",
		  "5.2441","5.2441"};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
  
  taylorFunction truncSq=taylorSimplex::unit*"1.583592135001261822";
  taylorFunction two=taylorSimplex::unit*"2.0";
 
  taylorFunction F=taylorSimplex::voronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::y2*"-0.037846"+
                   taylorSimplex::y3*"-0.05782"+ 
                   taylorSimplex::y5*"0"+
                   taylorSimplex::y6*"-0.014907"+
                   taylorSimplex::dih*"-0.110014"+
                   taylorSimplex::unit*"0.354003"; 
  


   
  const taylorFunction* K[1]={&F};  
  
  
  

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt)) 
    cout<<"xst1 complete!"<<endl;
  else
    cout<<"xst1 failed!!"<<endl;

}

void proveTet2()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","4","4","4"};

  interval tz[6]={"4.4025","4.41",
		  "4.41","6.33436854000504728617",
		  "5.2441","5.2441"};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);




  taylorFunction truncSq=taylorSimplex::unit*"1.583592135001261822";
  taylorFunction two=taylorSimplex::unit*"2.0";

  taylorFunction F=taylorSimplex::voronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::y2*"0.05782"+
                   taylorSimplex::y3*"0.05782"+ 
                   taylorSimplex::y5*"0"+
                   taylorSimplex::y6*"0"+
                   taylorSimplex::dih*"-0.110014"+
                   taylorSimplex::unit*"-0.098422";
 
  const taylorFunction* K[1]={&F};
  
  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt))
    cout<<"xst2 complete!"<<endl;
  else
    cout<<"xst 2 failed!!"<<endl;
}

void proveTet3()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","4","4","4"};

  interval tz[6]={"4.4025","4.41",
		  "4.41","6.33436854000504728617",
		  "5.2441","5.2441"};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction truncSq=taylorSimplex::unit*"1.583592135001261822";
  taylorFunction two=taylorSimplex::unit*"2.0";

  taylorFunction F=taylorSimplex::voronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::y2*"-0.05782"+
                   taylorSimplex::y3*"-0.05782"+ 
                   taylorSimplex::y5*"0"+
                   taylorSimplex::y6*"0"+
                   taylorSimplex::dih*"-0.110014"+
                   taylorSimplex::unit*"0.364137";
  
  
    
  const taylorFunction* K[1]={&F};
  
  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt))
    cout<<"xst3 complete!"<<endl;
  else
    cout<<"xst3 failed!!"<<endl;
}

void proveTet4()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","4","4","4"};

  interval tz[6]={"4.4025","4.41",
		  "4.41","6.33436854000504728617",
		  "5.2441","5.2441"};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);




  taylorFunction truncSq=taylorSimplex::unit*"1.583592135001261822";
  taylorFunction two=taylorSimplex::unit*"2.0";

  taylorFunction F=taylorSimplex::voronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::y2*"0.05782"+
                   taylorSimplex::y3*"0.05782"+ 
                   taylorSimplex::y5*"0"+
                   taylorSimplex::y6*"0"+
                   taylorSimplex::dih*"-0.110014"+
                   taylorSimplex::unit*"-0.098422";
  
  const taylorFunction* K[1]={&F};
  
  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt))
    cout<<"xst4 complete!"<<endl;
  else
    cout<<"xst4 failed!!"<<endl;
}

void proveTet5()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","8","4","4"};

  interval tz[6]={"4.2025","4.41",
		  "4.41","10.24",
		  "5.2411","5.2411"};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
  


  taylorFunction truncSq=taylorSimplex::unit*"1.583592135001261822";
  taylorFunction two=taylorSimplex::unit*"2.0";

  taylorFunction F=taylorSimplex::voronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::y2*"-0.05782"+
                   taylorSimplex::y3*"0.037846"+ 
                   taylorSimplex::y5*"0.014907"+
                   taylorSimplex::y6*"0"+
                   taylorSimplex::dih*"-0.110014"+
                   taylorSimplex::unit*"0.18605";

  const taylorFunction* K[1]={&F};
  
  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt))
    cout<<"xst5 complete!"<<endl;
  else
    cout<<"xst5 failed!!"<<endl;
}







