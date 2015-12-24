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
void diag2();
void vol1();
void vol2();
void vol3();

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
  diag2();  
  vol1();  
  vol2();
  vol3();
  
  return 0;
}

//*****************************************************************
 
  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  //  dih < 1.319, y1 < 2.243, y2 < 2.269 -> y4 < 2.757
  //
  /////////////////////////////////////////////////////////////////
  
void diag1()
{
  cellOption opt;
  domain x,z,x0,z0;  

  interval tx[6]={"4","4","4",v,"4","4"};
  interval tz[6]={"5.031049","5.148361",v,"12",v,v};

  // 2.243^2=5.031049, 2.269^2=5.148361, 2.169^2=4.704561, 
  // 12 > 2.757^2 = 7.601049

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);  
  
  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-2.757";

  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.319";

  opt.setDihMax(1.319);

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
  //  dih < 1.713, y1 < 2.169, y2 < 2.269, y3 < 2.360 -> y4 < 3.308
  //
  /////////////////////////////////////////////////////////////////
   
void diag2()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4",v,"4","4"};
  interval tz[6]={"4.704561","5.148361","5.5696","12",v,v};

  // 2.36^2=5.5696, 2.269^2=5.148361, 2.169^2=4.704561, 
  // 12 > 3.308^2 = 10.942864

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-3.308";  

  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.713";

  opt.setDihMax(1.713);
    
  const taylorFunction* K[2]={&F,&G};
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"diag2 complete!"<<endl;
  else
    cout<<"diag2 failed!!"<<endl;
}

//*****************************************************************
 
  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  //  dih < 1.319, y1 < 2.243, y2 < 2.269, y4 < 2.757 -> vol > .287
  //
  /////////////////////////////////////////////////////////////////
  
void vol1()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4",v,"4","4"};
  interval tz[6]={"5.031049","5.148361",v,"7.601049",v,v};

  // 2.243^2=5.031049, 2.269^2=5.148361, 2.169^2=4.704561, 
  // 2.757^2 = 7.601049

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.287";

  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.319";

  opt.setDihMax(1.319);
    
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
  //  dih < 1.713, y1 < 2.169, y2 < 2.269, y3 < 2.360,y4 < 3.308
  //
  //  -> vol > .264
  //
  /////////////////////////////////////////////////////////////////
   
void vol2()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4",v,"4","4"};
  interval tz[6]={"4.704561","5.148361","5.5696","10.942864",v,v};

  // 2.36^2=5.5696, 2.269^2=5.148361, 2.169^2=4.704561, 
  // 3.308^2 = 10.942864

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.247";
  
  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.713";
    
  opt.setDihMax(1.713);

  const taylorFunction* K[2]={&F,&G};  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"vol2 complete!"<<endl;
  else
    cout<<"vol2 failed!!"<<endl;
}

//*****************************************************************
 
  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  //  y1 < 2.269, y2 < 2.360, y5 < 3.308, y6 < 2.757 -> vol > .318
  //
  /////////////////////////////////////////////////////////////////

void vol3()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","4",v,v};
  interval tz[6]={"5.148361","5.5696",v,v,"10.942864","7.601049"};

  // 2.36^2=5.5696, 2.269^2=5.148361, 2.169^2=4.704561, 
  // 3.308^2 = 10.942864, 2.757^2 = 7.601049

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.300";
      
  const taylorFunction* K[1]={&F};  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt))
    cout<<"vol3 complete!"<<endl;
  else
    cout<<"vol3 failed!!"<<endl;
}


