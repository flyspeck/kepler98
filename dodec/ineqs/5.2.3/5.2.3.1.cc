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
	  interval tx[6]={"4","4","4","4.74063529","4","4"}; //2.1773^2
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
  setqr(x0,z0);
  setqr(x,z);
  taylorFunction truncSq=taylorSimplex::unit*"1.583592135001261822";
  taylorFunction two=taylorSimplex::unit*"2.0";
 
  taylorFunction F=taylorSimplex::voronoiVolume*"-1"
                  +taylorSimplex::sol*"0.42755"
                  +taylorSimplex::unit*"0.0047";
  
  taylorFunction G= two*"-1" + taylorSimplex::eta2_126;
  taylorFunction H= two*"-1" + taylorSimplex::eta2_135;
  taylorFunction I= two*"-1" + taylorSimplex::eta2_234;
  taylorFunction J= two*"-1" + taylorSimplex::eta2_456;

  const taylorFunction* K[5]={&F,&G,&H,&I,&J};  
  
  F.setReducibleState(1);
  

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,5,opt)) 
    cout<<"5.4.1 Part 1 complete!"<<endl;
  else
    cout<<"5.4.1 Part 1 failed!!"<<endl;

}

void proveTet2()
{
  cellOption opt;
  domain x,z,x0,z0;
  setqr(x0,z0);
  setqr(x,z);
  taylorFunction truncSq=taylorSimplex::unit*"1.583592135001261822";
  taylorFunction two=taylorSimplex::unit*"2.0";

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"
                  +taylorSimplex::sol*"0.42755"
                  +taylorSimplex::unit*"0.0047";
  
  taylorFunction G=taylorSimplex::eta2_126+two*"-1";
 
  const taylorFunction* K[2]={&F,&G};
  
  
  F.setReducibleState(1);
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"5.4.1 Part 2 complete!"<<endl;
  else
    cout<<"5.4.1 Part 2 failed!!"<<endl;
}

void proveTet3()
{
  cellOption opt;
  domain x,z,x0,z0;
  setqr(x0,z0);
  setqr(x,z);
  taylorFunction truncSq=taylorSimplex::unit*"1.583592135001261822";
  taylorFunction two=taylorSimplex::unit*"2.0";

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"
                  +taylorSimplex::sol*"0.42755"
                  +taylorSimplex::unit*"0.0047";  
  
  taylorFunction H=taylorSimplex::eta2_135+two*"-1";
  
  const taylorFunction* K[2]={&F,&H};
  
  
  F.setReducibleState(1);
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<" 5.4.1 Part 3 complete!"<<endl;
  else
    cout<<" 5.4.1 Part 3 failed!!"<<endl;
}

void proveTet4()
{
  cellOption opt;
  domain x,z,x0,z0;
  setqr(x0,z0);
  setqr(x,z);
  taylorFunction truncSq=taylorSimplex::unit*"1.583592135001261822";
  taylorFunction two=taylorSimplex::unit*"2.0";

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"
                  +taylorSimplex::sol*"0.42755"
                  +taylorSimplex::unit*"0.0047";
  
  
  taylorFunction I=taylorSimplex::eta2_234+two*"-1";
  
  const taylorFunction* K[2]={&F,&I};
  
  
  F.setReducibleState(1);
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<" 5.4.1 Part 4 complete!"<<endl;
  else
    cout<<" 5.4.1 Part 4 failed!!"<<endl;
}

void proveTet5()
{
  cellOption opt;
  domain x,z,x0,z0;
  setqr(x0,z0);
  setqr(x,z);
  taylorFunction truncSq=taylorSimplex::unit*"1.583592135001261822";
  taylorFunction two=taylorSimplex::unit*"2.0";

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"
                  +taylorSimplex::sol*"0.42755"
                  +taylorSimplex::unit*"0.0047";
  
  
  taylorFunction J=taylorSimplex::eta2_456+two*"-1";
  const taylorFunction* K[2]={&F,&J};
  
  
  F.setReducibleState(1);
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<" 5.4.1 Part 5 complete!"<<endl;
  else
    cout<<" 5.4.1 Part 5 failed!!"<<endl;
}







