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
	interval tx[6]={"4","4","4","4","4","4"};
	interval tz[6]={"6.33436854000504728617","6.33436854000504728617","6.33436854000504728617","6.33436854000504728617","6.33436854000504728617","6.33436854000504728617"};
	x = domain::lowerD(tx);
	z = domain::upperD(tz);
	}

int main()
{
 
  //  proveTet1();
  //proveTet2();
  //proveTet3();
  //proveTet4();
  proveTet5();
  proveTet6();
  //proveTet7();
  // proveTet8();
  // proveTet9();
  //proveTet10();


  return 0;
}


 
void proveTet1()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval v="6.33436854000504728617";// = (2 Sqrt[3]Tan[Pi/5])^2

  interval tx[6]={"5.0932365124","4","4","7.84","4","4"};

  interval tz[6]={v,v,v,"10.1124",v,v};

  //5.0932365124=2.25841^2 


  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
  
  
  taylorFunction two=taylorSimplex::unit*"2.0";
 
  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::y1*"-0.166"+
                   taylorSimplex::y2*"-0.174893"+
                   taylorSimplex::y3*"-0.174893"+ 
                   taylorSimplex::y4*"0.306136"+
                   taylorSimplex::y5*"0.143"+
                   taylorSimplex::y6*"0.143"+
                   taylorSimplex::dih*"-0.263488"+
                   taylorSimplex::unit*"-0.038250"; 

  //  taylorFunction G=taylorSimplex::dih*"-1"+taylorSimplex::unit*"1.701";

   
  const taylorFunction* K[1]={&F};  
  
    opt.setDihMax(1.701);
  

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt)) 
    cout<<"qe1 complete!"<<endl;
  else
    cout<<"qe1 failed!!"<<endl;

}




void proveTet2()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval v="6.33436854000504728617";// = (2 Sqrt[3]Tan[Pi/5])^2

  interval tx[6]={"4","4","4","7.84","4","4"};

  interval tz[6]={v,v,v,"10.1124",v,v};

  //5.0932365124=2.25841^2, 10.1124=3.179^2, 7.84=2.8^2


  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
  
  
  taylorFunction two=taylorSimplex::unit*"2.0";
 
  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::y1*"-0.166"+
                   taylorSimplex::y2*"0.008893"+
                   taylorSimplex::y3*"0.008893"+ 
                   taylorSimplex::y4*"-0.306136"+
                   taylorSimplex::y5*"0.143"+
                   taylorSimplex::y6*"0.143"+
                   taylorSimplex::dih*"0"+
                   taylorSimplex::unit*"0.80825"; 

    opt.setDihMax(2.652);

   
  const taylorFunction* K[1]={&F};  
  
  
  

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt)) 
    cout<<"qe2 complete!"<<endl;
  else
    cout<<"qe2 failed!!"<<endl;

}



void proveTet3()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval v="6.33436854000504728617";// = (2 Sqrt[3]Tan[Pi/5])^2

  interval tx[6]={"4","4","4","7.84","4","4"};

  interval tz[6]={"5.0932365124",v,v,"8.7025",v,v};

  //5.0932365124=2.25841^2, 10.1124=3.179^2, 7.84=2.8^2, 8.7025=2.95^2


  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
  
  
  taylorFunction two=taylorSimplex::unit*"2.0";
 
  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::y1*"-0.166"+
                   taylorSimplex::y2*"-0.120535"+
                   taylorSimplex::y3*"-0.120535"+ 
                   taylorSimplex::y4*"0.107880"+
                   taylorSimplex::y5*"0.143"+
                   taylorSimplex::y6*"0.143"+
                   taylorSimplex::dih*"-0.104704"+
                   taylorSimplex::unit*"0.251657"; 


    opt.setDihMax(1.701);

   
  const taylorFunction* K[1]={&F};  
  
  
  

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt)) 
    cout<<"qe3 complete!"<<endl;
  else
    cout<<"qe3 failed!!"<<endl;

}


void proveTet4()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval v="6.33436854000504728617"; // = (2 Sqrt[3]Tan[Pi/5])^2

  interval tx[6]={"4","4","4","7.84","4","4"};

  interval tz[6]={v,v,v,"8.7025",v,v};

  //5.0932365124=2.25841^2, 10.1124=3.179^2, 7.84=2.8^2, 8.7025=2.95^2


  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
  
  
  taylorFunction two=taylorSimplex::unit*"2.0";
 
  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::y1*"-0.166"+
                   taylorSimplex::y2*"-0.045465"+
                   taylorSimplex::y3*"-0.045465"+ 
                   taylorSimplex::y4*"-0.107880"+
                   taylorSimplex::y5*"0.143"+
                   taylorSimplex::y6*"0.143"+
                   taylorSimplex::dih*"0"+
                   taylorSimplex::unit*"0.518343"; 

   
  const taylorFunction* K[1]={&F};  
  
      opt.setDihMax(2.652);
  

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt)) 
    cout<<"qe4 complete!"<<endl;
  else
    cout<<"qe4 failed!!"<<endl;

}


void proveTet5()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval v="6.33436854000504728617"; // = (2 Sqrt[3]Tan[Pi/5])^2

  interval tx[6]={"4","4","4","8.7025","4","4"};

  interval tz[6]={"5.0932365124",v,v,"10.1124",v,v};

  //5.0932365124=2.25841^2, 10.1124=3.179^2, 7.84=2.8^2, 8.7025=2.95^2


  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
  
  
  taylorFunction two=taylorSimplex::unit*"2.0";
 
  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::y1*"-0.166"+
                   taylorSimplex::y2*"-0.114552"+
                   taylorSimplex::y3*"-0.114552"+ 
                   taylorSimplex::y4*"0.115382"+
                   taylorSimplex::y5*"0.143"+
                   taylorSimplex::y6*"0.143"+
                   taylorSimplex::dih*"-0.153420"+
                   taylorSimplex::unit*"0.22"; 

   
  const taylorFunction* K[1]={&F};  
  
    opt.setDihMax(1.701);
  

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt)) 
    cout<<"qe5 complete!"<<endl;
  else
    cout<<"qe5 failed!!"<<endl;

}


void proveTet6()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval v="6.33436854000504728617"; // = (2 Sqrt[3]Tan[Pi/5])^2

  interval tx[6]={"4","4","4","8.7025","4","4"};

  interval tz[6]={"5.0932365124",v,v,"10.1124",v,v};

  //5.0932365124=2.25841^2, 10.1124=3.179^2, 7.84=2.8^2, 8.7025=2.95^2


  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
  
  
  taylorFunction two=taylorSimplex::unit*"2.0";
 
  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::y1*"-0.166"+
                   taylorSimplex::y2*"-0.051448"+
                   taylorSimplex::y3*"-0.051448"+ 
                   taylorSimplex::y4*"-0.115382"+
                   taylorSimplex::y5*"0.143"+
                   taylorSimplex::y6*"0.143"+
                   taylorSimplex::dih*"0"+
                   taylorSimplex::unit*"0.55"; 

   
  const taylorFunction* K[1]={&F};  
  
  
      opt.setDihMax(2.652);

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt)) 
    cout<<"qe6 complete!"<<endl;
  else
    cout<<"qe6 failed!!"<<endl;

}


void proveTet7()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval v="6.33436854000504728617"; // = (2 Sqrt[3]Tan[Pi/5])^2

  interval tx[6]={"4","4","4","8.7025","4","4"};

  interval tz[6]={"5.0932365124","5.0932365124",v,"10.1124",v,v};

  //5.0932365124=2.25841^2, 10.1124=3.179^2, 7.84=2.8^2, 8.7025=2.95^2


  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
  
  
  taylorFunction two=taylorSimplex::unit*"2.0";
 
  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::y1*"-0.166"+
                   taylorSimplex::y2*"-0.279805"+
                   taylorSimplex::y3*"-0.340136"+ 
                   taylorSimplex::y4*"0.422343"+
                   taylorSimplex::y5*"0.143"+
                   taylorSimplex::y6*"0.143"+
                   taylorSimplex::dih*"-0.380615"+
                   taylorSimplex::unit*"0.147006"; 

   
  const taylorFunction* K[1]={&F};  
  
    opt.setDihMax(1.701);
  

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt)) 
    cout<<"qe7 complete!"<<endl;
  else
    cout<<"qe7 failed!!"<<endl;

}


void proveTet8()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval v="6.33436854000504728617"; // = (2 Sqrt[3]Tan[Pi/5])^2

  interval tx[6]={"5.0932365124","4","4","8.7025","4","4"};

  interval tz[6]={v,"5.0932365124",v,"10.1124",v,v};

  //5.0932365124=2.25841^2, 10.1124=3.179^2, 7.84=2.8^2, 8.7025=2.95^2


  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
  
  
  taylorFunction two=taylorSimplex::unit*"2.0";
 
  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::y1*"-0.166"+
                   taylorSimplex::y2*"0.113805"+
                   taylorSimplex::y3*"0.174136"+ 
                   taylorSimplex::y4*"-0.422343"+
                   taylorSimplex::y5*"0.143"+
                   taylorSimplex::y6*"0.143"+
                   taylorSimplex::dih*"0"+
                   taylorSimplex::unit*"0.622994"; 

   
  const taylorFunction* K[1]={&F};  
  
      opt.setDihMax(2.652);
  

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt)) 
    cout<<"qe8 complete!"<<endl;
  else
    cout<<"qe8 failed!!"<<endl;

}


void proveTet9()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval v="6.33436854000504728617"; // = (2 Sqrt[3]Tan[Pi/5])^2

  interval tx[6]={"4","5.0932365124","4","8.7025","4","4"};

  interval tz[6]={"5.0932365124",v,v,"10.1124",v,v};

  //5.0932365124=2.25841^2, 10.1124=3.179^2, 7.84=2.8^2, 8.7025=2.95^2


  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
  
  
  taylorFunction two=taylorSimplex::unit*"2.0";
 
  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::y1*"-0.166"+
                   taylorSimplex::y2*"-0.182394"+
                   taylorSimplex::y3*"-0.147301"+ 
                   taylorSimplex::y4*"0.250100"+
                   taylorSimplex::y5*"0.143"+
                   taylorSimplex::y6*"0.143"+
                   taylorSimplex::dih*"-0.206213"+
                   taylorSimplex::unit*"0.001504"; 

   
  const taylorFunction* K[1]={&F};  
  
    opt.setDihMax(1.701);
  

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt)) 
    cout<<"qe9 complete!"<<endl;
  else
    cout<<"qe9 failed!!"<<endl;

}


void proveTet10()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval v="6.33436854000504728617"; // = (2 Sqrt[3]Tan[Pi/5])^2

  interval tx[6]={"5.0932365124","5.0932365124","4","8.7025","4","4"};

  interval tz[6]={v,v,v,"10.1124",v,v};

  //5.0932365124=2.25841^2, 10.1124=3.179^2, 7.84=2.8^2, 8.7025=2.95^2


  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
  
  
  taylorFunction two=taylorSimplex::unit*"2.0";
 
  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::y1*"-0.166"+
                   taylorSimplex::y2*"0.016394"+
                   taylorSimplex::y3*"0.018699"+ 
                   taylorSimplex::y4*"-0.250100"+
                   taylorSimplex::y5*"0.143"+
                   taylorSimplex::y6*"0.143"+
                   taylorSimplex::dih*"-0.206213"+
                   taylorSimplex::unit*"0.768496"; 

   
  const taylorFunction* K[1]={&F};  
  
      opt.setDihMax(2.652);
  

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt)) 
    cout<<"qe10 complete!"<<endl;
  else
    cout<<"qe10 failed!!"<<endl;

}
