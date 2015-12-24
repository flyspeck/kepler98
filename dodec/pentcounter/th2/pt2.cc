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

void diag1();
void vol1();
void vol2();

void selfTest()
    {
    interMath::selfTest();
    linearization::selfTest();
    secondDerive::selfTest();
    taylorFunction::selfTest();
    }

const interval v="6.33436854000504728617"; // 2T^2

int main()
{
  diag1();
  vol1();  
  vol2();

  return 0;
}

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  //  dih < 1.705, y1 < 2.207, y3 < 2.222, y2 < 2.199 
  //
  //  -> y4 < 3.188
  //
  /////////////////////////////////////////////////////////////////
 
void diag1()
{
  cellOption opt;
  domain x,z,x0,z0;

  interval tx[6]={"4","4","4",v,"4","4"};
  interval tz[6]={"4.870849","4.835601","4.937284","12",v,v};
  
  // 2.207^2=4.870849, 2.222^2=4.937284, 2.199^2=4.835601, 
  // 12 > 3.188^2=10.163344
  
  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
  
  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-3.188";

  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.705";

  opt.setDihMax(1.705);

  const taylorFunction* K[2]={&F,&G};  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"diag1 complete!"<<endl;
  else
    cout<<"diag1 failed!!"<<endl;
}

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  //  dih < 1.705, y1 < 2.207, y3 < 2.222, y2 < 2.199, 
  //
  //  y4 < 3.188 -> vol > .277
  //
  /////////////////////////////////////////////////////////////////
 
void vol1()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4",v,"4","4"};
  interval tz[6]={"4.870849","4.835601","4.937284","10.163344",v,v};

  // 2.207^2=4.870849, 2.222^2=4.937284, 2.199^2=4.835601, 
  // 3.188^2=10.163344  

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.277";

  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.705";

  opt.setDihMax(1.705);
    
  const taylorFunction* K[2]={&F,&G};  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"vol1 complete!"<<endl;
  else
    cout<<"vol1 failed!!"<<endl;
}

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  //  y1 < 2.222, y2,y3 < 2.199, 
  //
  //  y5,y6 < 3.188 -> vol > .362
  //
  /////////////////////////////////////////////////////////////////
 

void vol2()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","4",v,v};
  interval tz[6]={"4.937284","4.835601","4.835601",v,
		  "10.163344","10.163344"};

  // 2.222^2=4.937284, 2.199^2=4.835601, 
  // 3.188^2=10.163344  

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.362";
      
  const taylorFunction* K[1]={&F};
    
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt))
    cout<<"vol2 complete!"<<endl;
  else
    cout<<"vol2 failed!!"<<endl;
}

