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
void proveTet6a();
void proveTet7();
void proveTet8();
void proveTet9();
void proveTet10();
void proveTet11();
void proveTet12();



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
  proveTet6a();
  proveTet7();
  proveTet8();
  proveTet9();
  proveTet10();  
  proveTet11();
  proveTet12();

  return 0;
}


 
void proveTet1()
{
  cellOption opt;
  domain x,z,x0,z0;
  

  interval tx[6]={"4","4","4",v,"4","4"};

  interval tz[6]={"4.923961","5.4289","4.870849","12",
		  v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
  
  
  
  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-3.179";
   
  taylorFunction G=taylorSimplex::y2*"-1"+
        taylorSimplex::y3*"-1"+             
    taylorSimplex::unit*"4.399";

  taylorFunction H=taylorSimplex::y2*"-1"+
    taylorSimplex::y3*"-1"+
    taylorSimplex::y1*"-1"+
    taylorSimplex::unit*"6.406";

     taylorFunction I=taylorSimplex::dih*"-1"+
      taylorSimplex::unit*"1.701";

  opt.setDihMax(1.701);

  const taylorFunction* K[4]={&F,&G,&H,&I};  
  
  
  

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,4,opt)) 
    cout<<"pent1 complete!"<<endl;
  else
    cout<<"pent1 failed!!"<<endl;

}

void proveTet2()
{
  cellOption opt;
  domain x,z,x0,z0;

  interval tx[6]={"4","4","4",v,"4","4"};

  interval tz[6]={"4.669921","4.870849","5.1984","12",
		  v,v};
  
  opt.setDihMax(1.681);

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  
  taylorFunction F=taylorSimplex::y4+
    taylorSimplex::unit*"-3.101";

  taylorFunction G=taylorSimplex::y2*"-1"+
    taylorSimplex::y3*"-1"+
    taylorSimplex::unit*"4.307";

  taylorFunction H=taylorSimplex::y2*"-1"+
    taylorSimplex::y3*"-1"+
    taylorSimplex::y1*"-1"+
    taylorSimplex::unit*"6.313";



 
  const taylorFunction* K[3]={&F,&G,&H};
  
  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,3,opt))
    cout<<"pent2 complete!"<<endl;
  else
    cout<<"pent2 failed!!"<<endl;
}

void proveTet3()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4",v,"4","4"};

  interval tz[6]={"5.424241","4.870849",
		  v,"14",
		  v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  


  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-3.611";
  

    opt.setDihMax(1.917);


   taylorFunction G=taylorSimplex::y2*"-1"+
     taylorSimplex::y3*"-1";
     taylorSimplex::unit*"4.657";

   taylorFunction H=taylorSimplex::y2*"-1"+
     taylorSimplex::y3*"-1"+
     taylorSimplex::y1*"-1"+
     taylorSimplex::unit*"6.757";



    
  const taylorFunction* K[3]={&F,&G,&H};
  

  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,3,opt))
    cout<<"pent3 complete!"<<endl;
  else
    cout<<"pent3 failed!!"<<endl;
}


void proveTet4()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4",v,"4","4"};

  interval tz[6]={"4.669921","4.870849",
		  "5.1984","9.616201",
		  v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.279";

  taylorFunction G=taylorSimplex::y2*"-1"+
    taylorSimplex::y3*"-1"+             
    taylorSimplex::unit*"4.307";

  taylorFunction H=taylorSimplex::y2*"-1"+
    taylorSimplex::y3*"-1"+
    taylorSimplex::y1*"-1"+
    taylorSimplex::unit*"6.313";

  
  opt.setDihMax(1.681);


    
  const taylorFunction* K[3]={&F,&G,&H};
  
  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,3,opt))
    cout<<"pent4 complete!"<<endl;
  else
    cout<<"pent4 failed!!"<<endl;
}

void proveTet5()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4",v,"4","4"};

  interval tz[6]={"4.519876","4.553956",
		  "4.704561","15",
		  v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-3.024";
  
  opt.setDihMax(1.634);
  
  //taylorFunction G=taylorSimplex::dih*"-1"+
  //               taylorSimplex::unit*"1.719";
    
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

  interval tz[6]={"5.121169","4.553956",
		  "5.1076","14",
		  v,v};


  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-3.399";
  
  opt.setDihMax(1.875);

  taylorFunction G=taylorSimplex::dih*"-1"+
                 taylorSimplex::unit*"1.875";
    
  const taylorFunction* K[2]={&F,&G};
  
  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"pent6 complete!"<<endl;
  else
    cout<<"pent6 failed!!"<<endl;
}


void proveTet6a()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","5.1076",v,"4","4"};

  interval tz[6]={"5.121169","4.553956",
		  v,"14",
		  v,v};


  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-3.554";
  
  opt.setDihMax(1.875);

  taylorFunction G=taylorSimplex::dih*"-1"+
                 taylorSimplex::unit*"1.875";
    
  const taylorFunction* K[2]={&F,&G};
  
  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"pent6a complete!"<<endl;
  else
    cout<<"pent6a failed!!"<<endl;
}


void proveTet7()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4",v,"4","4"};

  interval tz[6]={"4.519876","4.553956","4.704561","9.144576",v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.283";
  
  opt.setDihMax(1.634);



  const taylorFunction* K[1]={&F};
  

  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt))
    cout<<"pent7 complete!"<<endl;
  else
    cout<<"pent7 failed!!"<<endl;
} 


void proveTet8()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4",v,"4","4"};

  interval tz[6]={"5.121169","4.553956","5.1076","12.630916",v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.277";
  
  opt.setDihMax(1.875);


  const taylorFunction* K[1]={&F};
  

  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt))
    cout<<"pent8 complete!"<<endl;
  else
    cout<<"pent8 failed!!"<<endl;
} 

void proveTet9()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","4",v,v};

  interval tz[6]={"4.553956","4.704561","5.1076",v,"11.553201","9.144576"};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.335";
  

  const taylorFunction* K[1]={&F};
  

   
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt))
    cout<<"pent9 complete!"<<endl;
  else
    cout<<"pent9 failed!!"<<endl;
} 


void proveTet10()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4",v,"4","4"};

  interval tz[6]={"4.519876","4.553956","4.704561","9.144576",v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.0156";
  
    opt.setDihMax(1.634);

  //taylorFunction G=taylorSimplex::dih*"-1"+
  //               taylorSimplex::unit*"1.719";
    

  const taylorFunction* K[1]={&F};
  

    
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt))
    cout<<"pent10 complete!"<<endl;
  else
    cout<<"pent10 failed!!"<<endl;
} 


void proveTet11()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","5.1076",v,"4","4"};

  interval tz[6]={"5.121169","4.553956",v,"12.630916",v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.0304";
  
    opt.setDihMax(1.875);


    //  taylorFunction G=taylorSimplex::dih*"-1"+
    //           taylorSimplex::unit*"2.419";
    


  const taylorFunction* K[1]={&F};
  

  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt))
    cout<<"pent11 complete!"<<endl;
  else
    cout<<"pent11 failed!!"<<endl;
} 

void proveTet12()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","5.1076","4",v,v};

  interval tz[6]={"4.553956","4.704561",v,v,"12.630916","9.144576"};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.0406";
  


  const taylorFunction* K[1]={&F};
  

  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt))
    cout<<"pent12 complete!"<<endl;
  else
    cout<<"pent12 failed!!"<<endl;
} 

