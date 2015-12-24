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
void proveTet7();
void proveTet8();
void proveTet9();
void proveTet10();


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
  proveTet7();  
  proveTet8();  
  proveTet9();
  proveTet10();

  return 0;
}


 
void proveTet1()
{
  cellOption opt;
  domain x,z,x0,z0;
  

  interval tx[6]={"4","4","4",v,"4","4"};

  interval tz[6]={"4.813636","5.354596","5.1076","12",
		  v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
  
  
  
  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-2.697";
   


  taylorFunction G=taylorSimplex::dih*"-1"+
    taylorSimplex::unit*"1.341";

  opt.setDihMax(1.341);

  const taylorFunction* K[2]={&F,&G};
  
  
  

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt)) 
    cout<<"pent1 complete!"<<endl;
  else
    cout<<"pent1 failed!!"<<endl;

}


void proveTet2()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4",v,"4","4"};

  interval tz[6]={"4.635409","5.354596",
		  "4.726276","12",
		  v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  


  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-3.365";

  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.684";
  

    opt.setDihMax(1.684);



    
  const taylorFunction* K[2]={&F,&G};
  

  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"pent2 complete!"<<endl;
  else
    cout<<"pent2 failed!!"<<endl;
}


void proveTet3()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4",v,"4","4"};

  interval tz[6]={"4.813636","5.354596","5.1076","7.273809",v,v};




  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.328";

  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.341";

  
  opt.setDihMax(1.341);


    
  const taylorFunction* K[2]={&F,&G};
  
  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"pent3 complete!"<<endl;
  else
    cout<<"pent3 failed!!"<<endl;
}

void proveTet4()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4",v,"4","4"};

  interval tz[6]={"4.635409","5.354596",
		  "4.726276","11.323225",v,v};


  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  opt.setDihMax(1.684);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.272";
  
   taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.684";
    
  const taylorFunction* K[2]={&F,&G};
  
  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"pent4 complete!"<<endl;
  else
    cout<<"pent4 failed!!"<<endl;
}




void proveTet5()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","4",v,v};

  interval tz[6]={"5.354596","4.726276",
		    "5.1076",v,"7.273809","11.323225"};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.340";
  
 
    
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
  

  interval tx[6]={"4","4","5.1076",v,"4","4"};

  interval tz[6]={"4.813636","5.354596",v,"12",
		  v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
  
  
  
  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-2.807";
   


  taylorFunction G=taylorSimplex::dih*"-1"+
    taylorSimplex::unit*"1.341";

  opt.setDihMax(1.341);

  const taylorFunction* K[2]={&F,&G};
  
  
  

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt)) 
    cout<<"pent6 complete!"<<endl;
  else
    cout<<"pent6 failed!!"<<endl;

}


void proveTet7()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4",v,"4","4"};

  interval tz[6]={"4.635409","5.354596",
		  "4.726276","12",
		  v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  


  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-3.365";

  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.684";
  

    opt.setDihMax(1.684);



    
  const taylorFunction* K[2]={&F,&G};
  

  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"pent7 complete!"<<endl;
  else
    cout<<"pent7 failed!!"<<endl;
}


void proveTet8()
{
  cellOption opt;
  domain x,z,x0,z0;

  interval tx[6]={"4","4","5.1076",v,"4","4"};

  interval tz[6]={"4.813636","5.354596",v,"7.879249",
		  v,v};



  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
taylorSimplex::sol*"0.42755"+
taylorSimplex::unit*"0.0453";

  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.341";

  
  opt.setDihMax(1.341);


    
  const taylorFunction* K[2]={&F,&G};
  
  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"pent8 complete!"<<endl;
  else
    cout<<"pent8 failed!!"<<endl;
}

void proveTet9()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4",v,"4","4"};

  interval tz[6]={"4.635409","5.354596",
		  "4.726276","11.323225",
		  v,v};




  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  opt.setDihMax(1.684);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.0156";
  
   taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.684";
    
  const taylorFunction* K[2]={&F,&G};
  
  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"pent9 complete!"<<endl;
  else
    cout<<"pent9 failed!!"<<endl;
}




void proveTet10()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","5.1076","4",v,v};

  interval tz[6]={"5.354596","4.726276",
		    v,v,"7.879249","11.323225"};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.0480";
  
 
    
  const taylorFunction* K[1]={&F};
  
  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt))
    cout<<"pent10 complete!"<<endl;
  else
    cout<<"pent10 failed!!"<<endl;
}


