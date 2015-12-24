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
 
  taylorFunction F=taylorSimplex::sol*"-1"+
                   taylorSimplex::unit*"0.315696";
  
   
  const taylorFunction* K[1]={&F};  
  
  F.setReducibleState(1);
  

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt)) 
    cout<<"SolLB complete!"<<endl;
  else
    cout<<"SolLB failed!!"<<endl;

}


 
void proveTet2()
{
  cellOption opt;
  domain x,z,x0,z0;
  setqr(x0,z0);
  setqr(x,z);
  taylorFunction truncSq=taylorSimplex::unit*"1.583592135001261822";
  taylorFunction two=taylorSimplex::unit*"2.0";
 
  taylorFunction F=taylorSimplex::sol+
                   taylorSimplex::unit*"-1.051232";
  
   
  const taylorFunction* K[1]={&F};  
  
  F.setReducibleState(1);
  

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt)) 
    cout<<"SolUB complete!"<<endl;
  else
    cout<<"SolUB failed!!"<<endl;

}



 
void proveTet3()
{
  cellOption opt;
  domain x,z,x0,z0;
  setqr(x0,z0);
  setqr(x,z);
  taylorFunction truncSq=taylorSimplex::unit*"1.583592135001261822";
  taylorFunction two=taylorSimplex::unit*"2.0";
 
  taylorFunction F=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"0.856147";
  
   
  const taylorFunction* K[1]={&F};  
  
  F.setReducibleState(1);
  

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt)) 
    cout<<"DihLB complete!"<<endl;
  else
    cout<<"DihLB failed!!"<<endl;

}



 
void proveTet4()
{
  cellOption opt;
  domain x,z,x0,z0;
  setqr(x0,z0);
  setqr(x,z);
  taylorFunction truncSq=taylorSimplex::unit*"1.583592135001261822";
  taylorFunction two=taylorSimplex::unit*"2.0";
 
  taylorFunction F=taylorSimplex::dih+
                   taylorSimplex::unit*"-1.886730";
  
   
  const taylorFunction* K[1]={&F};  
  
  F.setReducibleState(1);
  

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt)) 
    cout<<"DihUB complete!"<<endl;
  else
    cout<<"DihUB failed!!"<<endl;

}
