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
void diag3();

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
  diag3();  
 
  return 0;
}

//*****************************************************************
 
  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  //  dih < 1.637, y1 < 2.130, y2 < 2.329, y3 < 2.179 -> y4 < 3.139
  //
  /////////////////////////////////////////////////////////////////
  
void diag1()
{
  cellOption opt;
  domain x,z,x0,z0;  

  interval tx[6]={"4","4","4",v,"4","4"};
  interval tz[6]={"4.5369","5.424241","4.748041","12",v,v};

  // 2.13^2=4.5369, 2.329^2=5.424241, 2.179^2=4.748041

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);  
  
  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-3.139";

  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.637";

  opt.setDihMax(1.637);

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
  //  dih < 1.625, y1 < 2.189, y2 < 2.373, y3 < 2.179 -> y4 < 3.145
  //
  /////////////////////////////////////////////////////////////////
   
void diag2()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4",v,"4","4"};
  interval tz[6]={"4.791721","5.631129","4.748041","12",v,v};

  // 2.189^2=4.791721, 2.373^2=5.631129, 2.179^2=4.748041

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-3.145";  

  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.625";

  opt.setDihMax(1.625);
    
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
  //  dih < 1.637, y1 < 2.130, y2 < 2.329, y3 < 2.179, y4 < 3.139
  //  
  //  -> vol > .271
  //
  /////////////////////////////////////////////////////////////////
  
void vol1()
{
  cellOption opt;
  domain x,z,x0,z0;

  interval tx[6]={"4","4","4",v,"4","4"};
  interval tz[6]={"4.5369","5.424241","4.748041","9.853321",v,v};

  // 2.13^2=4.5369, 2.329^2=5.424241, 2.179^2=4.748041,
  // 3.139^2=9.853321

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.271";

  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.637";

  opt.setDihMax(1.637);
    
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
  //  dih < 1.625, y1 < 2.189, y2 < 2.373, y3 < 2.179,y4 < 3.145
  //
  //  -> vol > .267
  //
  /////////////////////////////////////////////////////////////////
   
void vol2()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4",v,"4","4"};
  interval tz[6]={"4.791721","5.631129","4.748041","9.891025",v,v};

  // 2.189^2=4.791721, 2.373^2=5.631129, 2.179^2=4.748041,
  // 3.145^2=9.891025

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.267";
  
  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.625";
    
  opt.setDihMax(1.625);

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
  //  y1 < 2.179, y2 < 2.329, y3<2.373, y5 < 3.145, y6 < 3.139
  //
  // -> vol > .346
  //
  /////////////////////////////////////////////////////////////////

void vol3()
{
  cellOption opt;
  domain x,z,x0,z0;

  interval tx[6]={"4","4","4","4",v,v};
  interval tz[6]={"4.748041","5.424241","5.631129",v,"9.891025","9.853321"};

  // 2.373^2=5.631129, 2.179^2=4.748041,3.139^2=9.853321
  // 3.145^2=9.891025, 2.329^2=5.424241

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.346";
      
  const taylorFunction* K[1]={&F};  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt))
    cout<<"vol3 complete!"<<endl;
  else
    cout<<"vol3 failed!!"<<endl;
}

//*****************************************************************
 
  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  //  dih < 1.545, y1 < 2.078, y2 < 2.167, y3 > 2.076 -> y4 < 2.87
  //
  /////////////////////////////////////////////////////////////////
  
void diag3()
{
  cellOption opt;
  domain x,z,x0,z0;  

  interval tx[6]={"4","4","4","5.1076","4","4"};
  interval tz[6]={"4.318084","4.695889","4.309776","12",v,v};

  // 2.078^2=4.318084, 2.167^2=4.695889, 2.076^2=4.309776

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);  
  
  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-2.87";

  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.545";

  opt.setDihMax(1.545);

  const taylorFunction* K[2]={&F,&G};
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt)) 
    cout<<"diag3 complete!"<<endl;
  else
    cout<<"diag3 failed!!"<<endl;
}


