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


  //  proveTet1();
  //proveTet2();  
  proveTet3();  
  //proveTet4();
  //proveTet5();
  //proveTet6();
  
  return 0;
}


 
void proveTet1()
{
  cellOption opt;
  domain x,z,x0,z0;
  

  interval tx[6]={"4","4","4",v,"4","4"};

  interval tz[6]={"4.5369","5.424241","5.4.748041","12",
		  v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
  
  
  
  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-3.139";
   


  //   taylorFunction I=taylorSimplex::dih*"-1"+
  //    taylorSimplex::unit*"1.701";

  opt.setDihMax(1.637);

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

  interval tz[6]={"4.791721","4.748041",
		  "5.631129","12",
		  v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  


  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-3.145";
  

    opt.setDihMax(1.625);



    
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

    interval tz[6]={"4.5369","5.424241","4.748041","12",
		  v,v};

  




  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.271";

  
  opt.setDihMax(1.637);


    
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

    interval tz[6]={"4.791721","4.748041",
		  "5.631129","12",
		  v,v};




  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  opt.setDihMax(1.713);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.267";
  
 
    
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

  interval tz[6]={"4.748041","5.424241",
		    "5.631129",v,"9.891025","9.853321"};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.346";
  
 
    
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
  
  interval tx[6]={"4","4","4",v,"4","4"};

  interval tz[6]={"4.318084","4.695889",
		  "4.309776","12",
		  v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  


  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-2.87";
  

    opt.setDihMax(1.545);



    
  const taylorFunction* K[1]={&F};
  

  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt))
    cout<<"pent6 complete!"<<endl;
  else
    cout<<"pent6 failed!!"<<endl;
}



