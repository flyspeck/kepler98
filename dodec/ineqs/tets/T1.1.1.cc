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

  int i;
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
 
  taylorFunction F=taylorSimplex::voronoiVolume*"-1"+
                   taylorSimplex::unit*"0.202804";
  
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

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,5,opt)) 
    cout<<"Part 1, TetVolLB complete!"<<endl;
  else
    cout<<"Part 1, TetVolLB failed!!"<<endl;

}

void proveTet2()
{
  cellOption opt;
  domain x,z,x0,z0;
  setqr(x0,z0);
  setqr(x,z);
  taylorFunction truncSq=taylorSimplex::unit*"1.583592135001261822";
  taylorFunction two=taylorSimplex::unit*"2.0";

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.202804";
  
  taylorFunction G=taylorSimplex::eta2_126+two*"-1";
 
  const taylorFunction* K[2]={&F,&G};
  
  
  F.setReducibleState(1);
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"Part 2, TetVolLB complete!"<<endl;
  else
    cout<<"Part 2, TetVolLB failed!!"<<endl;
}

void proveTet3()
{
  cellOption opt;
  domain x,z,x0,z0;
  setqr(x0,z0);
  setqr(x,z);
  taylorFunction truncSq=taylorSimplex::unit*"1.583592135001261822";
  taylorFunction two=taylorSimplex::unit*"2.0";

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.202804";
  
  
  taylorFunction H=taylorSimplex::eta2_135+two*"-1";
  
  const taylorFunction* K[2]={&F,&H};
  
  
  F.setReducibleState(1);
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"Part 3, TetVolLB complete!"<<endl;
  else
    cout<<"Part 3, TetVolLB failed!!"<<endl;
}

void proveTet4()
{
  cellOption opt;
  domain x,z,x0,z0;
  setqr(x0,z0);
  setqr(x,z);
  taylorFunction truncSq=taylorSimplex::unit*"1.583592135001261822";
  taylorFunction two=taylorSimplex::unit*"2.0";

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.202804";
  
  
  taylorFunction I=taylorSimplex::eta2_234+two*"-1";
  
  const taylorFunction* K[2]={&F,&I};
  
  
  F.setReducibleState(1);
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"Part 4, TetVolLB complete!"<<endl;
  else
    cout<<"Part 4, TetVolLB failed!!"<<endl;
}

void proveTet5()
{
  cellOption opt;
  domain x,z,x0,z0;
  setqr(x0,z0);
  setqr(x,z);
  taylorFunction truncSq=taylorSimplex::unit*"1.583592135001261822";
  taylorFunction two=taylorSimplex::unit*"2.0";

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.202804";
  
  
  taylorFunction J=taylorSimplex::eta2_456+two*"-1";
  const taylorFunction* K[2]={&F,&J};
  
  
  F.setReducibleState(1);
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"Part 5, TetVolLB complete!"<<endl;
  else
    cout<<"Part 5, TetVolLB failed!!"<<endl;
}







