#include<iostream.h>
#include<iomanip.h>
#include<fstream.h>
#include "error.h"
#include "interval.h"
#include "lineInterval.h"
#include "secondDerive.h"
#include "taylorInterval.h"
#include "recurse.h"



int proveTet1();
int proveTet2();


void selfTest()
    {
    interMath::selfTest();
    linearization::selfTest();
    secondDerive::selfTest();
    taylorFunction::selfTest();
    }

static void setqr1(domain& x,domain& z)
	{
	interval tx[6]={"4","4","4","6.33436854000504728617","4","4"};
	interval tz[6]={"6.33436854000504728617","6.33436854000504728617","6.33436854000504728617","6.33436854000504728617","6.33436854000504728617","6.33436854000504728617"};
	x = domain::lowerD(tx);
	z = domain::upperD(tz);
	}


static void setqr2(domain& x,domain& z)
	{
	interval tx[6]={"4","4","4","6.33436854000504728617","4","6.33436854000504728617"};
	interval tz[6]={"6.33436854000504728617","6.33436854000504728617","6.33436854000504728617","6.33436854000504728617","6.33436854000504728617","8"};
	x = domain::lowerD(tx);
	z = domain::upperD(tz);
	}






int main()
{
  if(proveTet1()) cout<<"DihLB Part 1 : Passed!"<<endl;
  else cout<<"DihLB Part 1 : Failed!"<<endl;

  if(proveTet2()) cout<<"DihUB Part 1 : Passed!"<<endl;
  else cout<<"DihUB Part 1 : Failed!"<<endl;

  return 0;
}


int proveTet1()
{
  cellOption opt;
  domain x,z,x0,z0;
  setqr1(x0,z0);
  setqr1(x,z);
  taylorFunction truncSq=taylorSimplex::unit*"1.583592135001261822";
  taylorFunction two=taylorSimplex::unit*"2.0";
 
  taylorFunction F=taylorSimplex::dih*"-1"+taylorSimplex::unit*"1.15242";
  
     
  const taylorFunction* K[1]={&F};
  
  
  F.setReducibleState(1);
 
  return prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt);
}


int proveTet2()
{
  cellOption opt;
  domain x,z,x0,z0;
  setqr2(x0,z0);
  setqr2(x,z);
  taylorFunction truncSq=taylorSimplex::unit*"1.583592135001261822";
  taylorFunction two=taylorSimplex::unit*"2.0";
 
  taylorFunction F=taylorSimplex::dih+taylorSimplex::unit*"-1.629435";
  // 1.629435 = bound/2
     
  const taylorFunction* K[1]={&F};
  
  
  F.setReducibleState(1);
 
  return prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt);
}







