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
	interval tx[6]={"4","4","4","12.96","4","4"};
	interval tz[6]={"5.166529","5.035536","6.33436854000504728617","12.96","6.33436854000504728617","6.33436854000504728617"};
	x = domain::lowerD(tx);
	z = domain::upperD(tz);
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
  

  interval tx[6]={"4","4","4",v,"4","4"};

  interval tz[6]={"4.704561","5.128361","5.5696","12",
		  v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
  
  
  
  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-2.683";
   


  //   taylorFunction I=taylorSimplex::dih*"-1"+
  //    taylorSimplex::unit*"1.701";

  opt.setDihMax(1.319);

  const taylorFunction* K[1]={&F};
  
  
  

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt)) 
    cout<<"pent1 complete!"<<endl;
  else
    cout<<"pent1 failed!!"<<endl;

}


void proveTet2()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4",v,"4","4"};

  interval tz[6]={"5.031049","5.148361",
		  v,"12",
		  v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  


  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-3.390";
  

    opt.setDihMax(1.713);



    
  const taylorFunction* K[1]={&F};
  

  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt))
    cout<<"pent2 complete!"<<endl;
  else
    cout<<"pent2 failed!!"<<endl;
}


void proveTet3()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4",v,"4","4"};

  interval tz[6]={"4.4.704561","5.128361","5.5696","7.198489",v,v};




  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.20";

  
  opt.setDihMax(1.319);


    
  const taylorFunction* K[1]={&F};
  
  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt))
    cout<<"pent3 complete!"<<endl;
  else
    cout<<"pent3 failed!!"<<endl;
}

void proveTet4()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4",v,"4","4"};

  interval tz[6]={"5.031049","5.148361",v,"11.4921",v,v};


  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  opt.setDihMax(1.713);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.247";
  
 
    
  const taylorFunction* K[1]={&F};
  
  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt))
    cout<<"pent4 complete!"<<endl;
  else
    cout<<"pent4 failed!!"<<endl;
}




void proveTet5()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","4",v,v};

  interval tz[6]={"5.148361","5.5696",
		    v,v,"11.4921","7.198489"};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.300";
  
 
    
  const taylorFunction* K[1]={&F};
  
  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt))
    cout<<"pent5 complete!"<<endl;
  else
    cout<<"pent5 failed!!"<<endl;
}


