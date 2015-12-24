

// This file includes all simplex verifications needed for the proof of 
// Theorem 1 of Appendix 2 of the proof of the dodecahedral conjecture. *)


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



void long1();
void long2();

const interval v="6.33436854000504728617";  // = (2 Sqrt[3]Tan[Pi/5])^2

void selfTest()
    {
    interMath::selfTest();
    linearization::selfTest();
    secondDerive::selfTest();
    taylorFunction::selfTest();
    }


int main()
{
  long1();
  long2();
 
  return 0;
}

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  // 
  // dih<1.57, y1<2.040, y2<2.220, y3<2.222 -> y4<3.1
  //
  /////////////////////////////////////////////////////////////////

void long1()
{
  cellOption opt;
  domain x,z,x0,z0;  

  interval tx[6]={"4","4","4",v,"4","4"};
  interval tz[6]={"4.1616","4.9282","4.937284","15",v,v};

  //2.04^2=4.1616, 2.22^2=4.9282, 2.222^2=4.937284

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
      
  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-3.1";
   
  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.57";

  opt.setDihMax(1.57);

  const taylorFunction* K[2]={&F,&G};  

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt)) 
    cout<<"long1 complete!"<<endl;
  else
    cout<<"long1 failed!!"<<endl;

}

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  // 
  // dih<1.35, y1<2.005, y2<2.005, y3<2.005 -> y4<2.8
  //
  /////////////////////////////////////////////////////////////////


void long2()
{
  cellOption opt;
  domain x,z,x0,z0;

  interval tx[6]={"4","4","4","12.96","4","4"};
  interval tz[6]={"4.020025","4.020025","4.020025","15",v,v};
  
    //2.005^2=4.020025

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
  
  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-2.8";
 
  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.535";

  opt.setDihMax(1.535);

  const taylorFunction* K[2]={&F,&G};  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"long2 complete!"<<endl;
  else
    cout<<"long2 failed!!"<<endl;
}

