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


const interval v="6.33436854000504728617";


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
  

  interval tx[6]={"5.29","4","4","4",v,"4"};

  interval tz[6]={v,v,v,v,"9.3636",v};


  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
  
  
  
  taylorFunction F=taylorSimplex::dih+
                   taylorSimplex::unit*"-1.630";
   



  const taylorFunction* K[1]={&F};
  
  
  

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt)) 
    cout<<"exc1 complete!"<<endl;
  else
    cout<<"exc1 failed!!"<<endl;

}


void proveTet2()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"5.29","4","4","4",v,v};

  interval tz[6]={v,v,v,v,"8","9.3636"};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  


  taylorFunction F=taylorSimplex::dih+
                   taylorSimplex::unit*"-1.512";
  

  const taylorFunction* K[1]={&F};
  

  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt))
    cout<<"exc2 complete!"<<endl;
  else
    cout<<"exc2 failed!!"<<endl;
}


void proveTet3()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"5.29","4","4","4","4","4"};

  interval tz[6]={v,v,v,v,v,v};




  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::voronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::dih*"0.02365"+
                   taylorSimplex::unit*"-0.0126";

  
    
  const taylorFunction* K[1]={&F};
  
  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt))
    cout<<"exc3 complete!"<<endl;
  else
    cout<<"exc3 failed!!"<<endl;
}

void proveTet4()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"5.29","4","4","4","4",v};

  interval tz[6]={v,v,v,v,v,"9.3636"};


  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::dih*"0.02365"+
                   taylorSimplex::unit*"0.0026";
  
 
    
  const taylorFunction* K[1]={&F};
  
  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt))
    cout<<"exc4 complete!"<<endl;
  else
    cout<<"exc4 failed!!"<<endl;
}




void proveTet5()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"5.29","4","4","4",v,v};

  interval tz[6]={v,v,v,v,"8","9.3636"};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::dih*"0.02365"+
                   taylorSimplex::unit*"0.0274";
  
  
 
    
  const taylorFunction* K[1]={&F};
  
  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt))
    cout<<"exc5 complete!"<<endl;
  else
    cout<<"exc5 failed!!"<<endl;
}


