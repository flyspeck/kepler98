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
void diag3();
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
  diag3();
  vol1();  
  vol2();
  vol3();
  
  return 0;
}

//*****************************************************************
 
  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  //  dih < 1.676, y1 < 2.176, y2 < 2.343, y3 < 2.205 -> y4 < 3.126
  //
  /////////////////////////////////////////////////////////////////
  
void diag1()
{
  cellOption opt;
  domain x,z,x0,z0;  

  interval tx[6]={"4","4","4",v,"4","4"};
  interval tz[6]={"4.734976","5.489649","4.862025","12",v,v};

  // 2.176^2=4.734976, 2.343^2=5.489649, 2.205^2=4.862025

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);  
  
  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-3.126";

  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.676";

  taylorFunction H=taylorSimplex::y5*"-1"+
                   taylorSimplex::y6*"-1"+
                   taylorSimplex::unit*"4.965";

  taylorFunction I=taylorSimplex::y2*"-1"+
                   taylorSimplex::y3*"-1"+
                   taylorSimplex::unit*"4.386";

  opt.setDihMax(1.676);

  const taylorFunction* K[4]={&F,&G,&H,&I};
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,4,opt)) 
    cout<<"diag1 complete!"<<endl;
  else
    cout<<"diag1 failed!!"<<endl;
}

//*****************************************************************
 
  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  //  dih < 1.685, y1 < 2.165, y2 < 2.221, y3 < 2.202 -> y4 < 3.161
  //
  /////////////////////////////////////////////////////////////////
   
void diag2()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4",v,"4","4"};
  interval tz[6]={"4.687225","4.932841","4.848804","12",v,v};

  // 2.165^2=4.687225, 2.221^2=4.932841, 2.202^2=4.848804
  
  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-3.161";  

  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.685";

  opt.setDihMax(1.685);
    
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
  //  dih < 1.828, y1 < 2.343, y2 < 2.221, y3 < 2.23 -> y4 < 3.370
  //
  /////////////////////////////////////////////////////////////////
   
void diag3()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4",v,"4","4"};
  interval tz[6]={"5.489649","4.932841","4.9729","12",v,v};

  // 2.343^2=5.489649, 2.221^2=4.932841, 2.23^2=4.9729
  
  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-3.37";  

  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.828";

  opt.setDihMax(1.828);
    
  const taylorFunction* K[2]={&F,&G};
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"diag3 complete!"<<endl;
  else
    cout<<"diag3 failed!!"<<endl;
}



//*****************************************************************
 
  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  //  dih < 1.685, y1 < 2.165, y2 < 2.221, y3 < 2.202 
  //
  //  y4 < 3.161 -> vol > .277
  //
  /////////////////////////////////////////////////////////////////
  
void vol1()
{
  cellOption opt;
  domain x,z,x0,z0;

  interval tx[6]={"4","4","4",v,"4","4"};
  interval tz[6]={"4.687225","4.932841","4.848804","9.991921",v,v};

  // 2.165^2=4.687225, 2.221^2=4.932841, 2.202^2=4.848804, 
  // 3.161^2=9.991921

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.277";

  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.685";

  opt.setDihMax(1.85);
    
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
  //  dih < 1.828, y1 < 2.343, y2 < 2.221, y3 < 2.23,y4 < 3.37
  //
  //  -> vol > .267
  //
  /////////////////////////////////////////////////////////////////
   
void vol2()
{
  cellOption opt;
  domain x,z,x0,z0;

  interval tx[6]={"4","4","4",v,"4","4"};
  interval tz[6]={"5.489649","4.932841","4.9729","11.3569",v,v};

  // 2.343^2=5.489649, 2.221^2=4.932841, 2.23^2=4.9729, 3.37^2=11.3569
  
  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.267";
  
  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.828";
    
  opt.setDihMax(1.828);

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
  //  y1 < 2.221, y2 < 2.23, y3 < 2.202, y5 < 3.161, y6 < 3.37
  //
  //  -> vol > .36
  //
  /////////////////////////////////////////////////////////////////

void vol3()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","4",v,v};
  interval tz[6]={"4.932841","4.9729","4.848804",v,"9.991921","11.3569"};

  // 2.165^2=4.687225, 2.221^2=4.932841, 2.202^2=4.848804, 
  // 3.161^2=9.991921, 2.23^2=4.9729, 3.37^2=11.3669

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.36";
      
  const taylorFunction* K[1]={&F};  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt))
    cout<<"vol3 complete!"<<endl;
  else
    cout<<"vol3 failed!!"<<endl;
}


