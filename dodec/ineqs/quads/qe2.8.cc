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
     interval tx[6]={"4","4","4","6.33436854000504728617","4","4"};
	interval tz[6]={"6.33436854000504728617","6.33436854000504728617","6.33436854000504728617","7.84","6.33436854000504728617","6.33436854000504728617"};
	x = domain::lowerD(tx);
	z = domain::upperD(tz);
	}

int main()
{
 
  proveTet1();
  
 
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
 
  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"
                  +taylorSimplex::y1*"-0.166"
                  +taylorSimplex::y2*"-0.083" 
                  +taylorSimplex::y3*"-0.083"
                  +taylorSimplex::y5*"0.143"
                  +taylorSimplex::y6*"0.143"
                  +taylorSimplex::unit*"0.3872455";
  
  
  const taylorFunction* K[1]={&F};
  

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt)) 
    cout<<"quad edge Part 1 complete!"<<endl;
  else
    cout<<"quad edge Part 1 failed!!"<<endl;

}


