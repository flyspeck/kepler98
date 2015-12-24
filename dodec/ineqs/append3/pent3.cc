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
void proveTet6();

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
  proveTet6();
  return 0;
}


 
void proveTet1()
{
  cellOption opt;
  domain x,z,x0,z0;
  

  interval tx[6]={"4","4","4",v,"4","4"};

  interval tz[6]={"4.734976","5.489649","4.862025","12",
		  v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
  
  
  
  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-3.126";
   
  taylorFunction G=taylorSimplex::y2*"-1"+
    taylorSimplex::y2*"-1"+
    taylorSimplex::unit*"4.386";

  taylorFunction H=taylorSimplex::y6*"-1"+
    taylorSimplex::y5*"-1"+
    taylorSimplex::unit*"4.965";


  //   taylorFunction I=taylorSimplex::dih*"-1"+
  //    taylorSimplex::unit*"1.701";

  opt.setDihMax(1.676);

  const taylorFunction* K[3]={&F,&G,&H};  
  
  
  

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,3,opt)) 
    cout<<"pent1 complete!"<<endl;
  else
    cout<<"pent1 failed!!"<<endl;

}

void proveTet2()
{
  cellOption opt;
  domain x,z,x0,z0;

  interval tx[6]={"4","4","4",v,"4","4"};

  interval tz[6]={"4.687225","4.932841","4.848804","12",
		  v,v};
  
  opt.setDihMax(1.685);

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  
  taylorFunction F=taylorSimplex::y4+
    taylorSimplex::unit*"-3.161";

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

  interval tz[6]={"5.489649","4.932841",
		  "4.9729","12",
		  v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  


  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-3.370";
  

    opt.setDihMax(1.828);



    
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

  interval tz[6]={"4.687225","4.932841","4.848804","9.991921",v,v};


  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.277";

  
  opt.setDihMax(1.685);


    
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
  
  interval tx[6]={"4","4","4",v,"4","4"};

    interval tz[6]={"5.489649","4.932841",
		    "4.9729","11.3569",v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  opt.setDihMax(1.828);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.260";
  
 
    
  const taylorFunction* K[1]={&F};
  
  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt))
    cout<<"pent5 complete!"<<endl;
  else
    cout<<"pent5 failed!!"<<endl;
}




void proveTet6()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","4",v,v};

  interval tz[6]={"4.932841","4.9729",
		    "4.848804",v,"9.991921","11.3569"};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.360";
  
 
    
  const taylorFunction* K[1]={&F};
  
  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt))
    cout<<"pent6 complete!"<<endl;
  else
    cout<<"pent6 failed!!"<<endl;
}


