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
void proveTet11();
void proveTet12();
void proveTet13();
void proveTet14();
void proveTet15();
void proveTet16();
void proveTet17();
void proveTet18();
void proveTet19();
void proveTet20();
void proveTet21();
void proveTet22();
void proveTet23();
void proveTet24();
void proveTet25();
void proveTet26();
void proveTet27();
void proveTet28();
void proveTet29();



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
  proveTet6();
  proveTet7();
  proveTet8();
  proveTet9();
  proveTet10();  
  proveTet11();
  proveTet12();
  proveTet13();
  proveTet14();
  proveTet15();
  proveTet16();
  proveTet17();
  proveTet18();
  proveTet19();
  proveTet20();
  proveTet21();
  proveTet22();
  proveTet23();
  proveTet24();
  proveTet25();
  proveTet26();
  proveTet27();
  proveTet28();
  proveTet29();
  

  return 0;
}


 
void proveTet1()
{
  cellOption opt;
  domain x,z,x0,z0;
  

  interval tx[6]={"4","4","4","12.96","4","4"};

  interval tz[6]={"5.166529","5.035536",
		  v,"12.96",
		  v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
  
  opt.setDihMax(1.86);
  
  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-3.6";
   
  //taylorFunction G=taylorSimplex::dih*"-1"+
  //               taylorSimplex::unit*"1.86";

  const taylorFunction* K[1]={&F};  
  
  
  

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt)) 
    cout<<"quad1 complete!"<<endl;
  else
    cout<<"quad1 failed!!"<<endl;

}

void proveTet2()
{
  cellOption opt;
  domain x,z,x0,z0;

  interval tx[6]={"4","4","4","12.96","4","4"};

  interval tz[6]={"5.035536","5.16653",
		  v,"15",
		  v,v};
  
  opt.setDihMax(1.86);

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  
  taylorFunction F=taylorSimplex::y4+
    taylorSimplex::unit*"-3.6";

 
  const taylorFunction* K[1]={&F};
  
  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt))
    cout<<"quad2 complete!"<<endl;
  else
    cout<<"quad2 failed!!"<<endl;
}

void proveTet3()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","9.61","4","4"};

  interval tz[6]={"5.16653","5.035536",
		  v,"10.24",
		  v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  


  taylorFunction F=taylorSimplex::y2*"-1"+
                   taylorSimplex::y3*"-1"+
                   taylorSimplex::unit*"4.21";
  

    opt.setDihMax(1.719);


   taylorFunction G=taylorSimplex::dih*"-1"+
            taylorSimplex::unit*"1.719";


    
  const taylorFunction* K[2]={&F,&G};
  

  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"quad3 complete!"<<endl;
  else
    cout<<"quad3 failed!!"<<endl;
}


void proveTet4()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","10.24","4","4"};

  interval tz[6]={"5.16653","5.035536",
		  v,"10.89",
		  v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::y2*"-1"+
                   taylorSimplex::y3*"-1"+
                   taylorSimplex::unit*"4.37";
  
  opt.setDihMax(1.719);

  //taylorFunction G=taylorSimplex::dih*"-1"+
  //               taylorSimplex::unit*"1.719";
    
  const taylorFunction* K[1]={&F};
  
  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt))
    cout<<"quad4 complete!"<<endl;
  else
    cout<<"quad4 failed!!"<<endl;
}

void proveTet5()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","10.89","4","4"};

  interval tz[6]={"5.16653","5.035536",
		  v,"11.56",
		  v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::y2*"-1"+
                   taylorSimplex::y3*"-1"+
                   taylorSimplex::unit*"4.53";
  
  opt.setDihMax(1.719);
  
  //taylorFunction G=taylorSimplex::dih*"-1"+
  //               taylorSimplex::unit*"1.719";
    
  const taylorFunction* K[1]={&F};
  
  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt))
    cout<<"quad5 complete!"<<endl;
  else
    cout<<"quad5 failed!!"<<endl;
}

void proveTet6()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","11.56","4","4"};

  interval tz[6]={"5.16653","5.035536",
		  v,"12.96",
		  v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::y2*"-1"+
                   taylorSimplex::y3*"-1"+
                   taylorSimplex::unit*"4.8";
  
  opt.setDihMax(1.719);

  //taylorFunction G=taylorSimplex::dih*"-1"+
  //               taylorSimplex::unit*"1.719";
    
  const taylorFunction* K[1]={&F};
  
  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt))
    cout<<"quad6 complete!"<<endl;
  else
    cout<<"quad6 failed!!"<<endl;
}



void proveTet7()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","9.61","4","4"};

  interval tz[6]={"5.16653","5.035536",v,"10.24",v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.057";
  
  opt.setDihMax(1.719);


  //taylorFunction G=taylorSimplex::dih*"-1"+
  //               taylorSimplex::unit*"1.719";
    
  taylorFunction G=taylorSimplex::y2+
                   taylorSimplex::y3+
                   taylorSimplex::unit*"-4.21";

  const taylorFunction* K[2]={&F,&G};
  

  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"quad7 complete!"<<endl;
  else
    cout<<"quad7 failed!!"<<endl;
} 


void proveTet8()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","10.24","4","4"};

  interval tz[6]={"5.16653","5.035536",v,"10.89",v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.069";
  
  opt.setDihMax(1.719);


  //taylorFunction G=taylorSimplex::dih*"-1"+
  //               taylorSimplex::unit*"1.719";
    
  taylorFunction G=taylorSimplex::y2+
                   taylorSimplex::y3+
                   taylorSimplex::unit*"-4.37";

  const taylorFunction* K[2]={&F,&G};
  

  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"quad8 complete!"<<endl;
  else
    cout<<"quad8 failed!!"<<endl;
} 

void proveTet9()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","10.89","4","4"};

  interval tz[6]={"5.16653","5.035536",v,"11.56",v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.083";
  
    opt.setDihMax(1.719);

  //taylorFunction G=taylorSimplex::dih*"-1"+
  //               taylorSimplex::unit*"1.719";
    
  taylorFunction G=taylorSimplex::y2+
                   taylorSimplex::y3+
                   taylorSimplex::unit*"-4.53";

  const taylorFunction* K[2]={&F,&G};
  

   
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"quad9 complete!"<<endl;
  else
    cout<<"quad9 failed!!"<<endl;
} 


void proveTet10()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","11.56","4","4"};

  interval tz[6]={"5.16653","5.035536",v,"12.96",v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.11";
  
    opt.setDihMax(1.719);

  //taylorFunction G=taylorSimplex::dih*"-1"+
  //               taylorSimplex::unit*"1.719";
    
  taylorFunction G=taylorSimplex::y2+
                   taylorSimplex::y3+
                   taylorSimplex::unit*"-4.8";

  const taylorFunction* K[2]={&F,&G};
  

    
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"quad10 complete!"<<endl;
  else
    cout<<"quad10 failed!!"<<endl;
} 


void proveTet11()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","9.61","4","4"};

  interval tz[6]={v,"5.035536",v,"10.24",v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.0196";
  
    opt.setDihMax(2.419);


  taylorFunction G=taylorSimplex::dih*"-1"+
                 taylorSimplex::unit*"2.419";
    
  taylorFunction H=taylorSimplex::y2+
                   taylorSimplex::y3+
                   taylorSimplex::unit*"-4.21";

  const taylorFunction* K[3]={&F,&G,&H};
  

  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,3,opt))
    cout<<"quad11 complete!"<<endl;
  else
    cout<<"quad11 failed!!"<<endl;
} 

void proveTet12()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","10.24","4","4"};

  interval tz[6]={v,"5.035536",v,"10.89",v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.016";
  
    opt.setDihMax(2.419);


  //taylorFunction G=taylorSimplex::dih*"-1"+
  //               taylorSimplex::unit*"2.419";
    
  taylorFunction G=taylorSimplex::y2+
                   taylorSimplex::y3+
                   taylorSimplex::unit*"-4.37";

  const taylorFunction* K[2]={&F,&G};
  

  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"quad12 complete!"<<endl;
  else
    cout<<"quad12 failed!!"<<endl;
} 


void proveTet13()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","10.89","4","4"};

  interval tz[6]={v,"5.035536",v,"11.56",v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.009";
  
    opt.setDihMax(2.419);


  //taylorFunction G=taylorSimplex::dih*"-1"+
  //               taylorSimplex::unit*"2.419";
    
  taylorFunction G=taylorSimplex::y2+
                   taylorSimplex::y3+
                   taylorSimplex::unit*"-4.53";

  const taylorFunction* K[2]={&F,&G};
  

  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"quad13 complete!"<<endl;
  else
    cout<<"quad13 failed!!"<<endl;
} 


void proveTet14()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","10.24","4","4"};

  interval tz[6]={"5.166529","5.035536",v,"10.89",v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.12";
  
    opt.setDihMax(1.532);


  //taylorFunction G=taylorSimplex::dih*"-1"+
  //               taylorSimplex::unit*"1.532";
    
  taylorFunction G=taylorSimplex::y2+
                   taylorSimplex::y3+
                   taylorSimplex::unit*"-4.37";

  const taylorFunction* K[2]={&F,&G};
  

 
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"quad14 complete!"<<endl;
  else
    cout<<"quad14 failed!!"<<endl;
} 


void proveTet15()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","9.61","4","4"};

  interval tz[6]={"5.166529","5.035536",v,"10.24",v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.076";
  
    opt.setDihMax(1.618);

  //taylorFunction G=taylorSimplex::dih*"-1"+
  //               taylorSimplex::unit*"1.614";
    
  taylorFunction G=taylorSimplex::y2+
                   taylorSimplex::y3+
                   taylorSimplex::unit*"-4.21";

  const taylorFunction* K[2]={&F,&G};
  

   
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"quad15 complete!"<<endl;
  else
    cout<<"quad15 failed!!"<<endl;
} 

//***

void proveTet16()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","9.61","4","4"};

  interval tz[6]={"5.16653","5.035536",
		  v,"10.24",
		  v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::y2*"-1"+
                   taylorSimplex::y3*"-1"+
                   taylorSimplex::unit*"4.134";
  
    opt.setDihMax(1.76);


  taylorFunction G=taylorSimplex::dih*"-1"+
                 taylorSimplex::unit*"1.76";
    
  const taylorFunction* K[2]={&F,&G};
  
 
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"quad16 complete!"<<endl;
  else
    cout<<"quad16 failed!!"<<endl;
}


void proveTet17()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","10.24","4","4"};

  interval tz[6]={"5.16653","5.035536",
		  v,"10.89",
		  v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::y2*"-1"+
                   taylorSimplex::y3*"-1"+
                   taylorSimplex::unit*"4.29";
  
  opt.setDihMax(1.76);


  taylorFunction G=taylorSimplex::dih*"-1"+
                 taylorSimplex::unit*"1.76";
    
  const taylorFunction* K[2]={&F,&G};
  
  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"quad17 complete!"<<endl;
  else
    cout<<"quad17 failed!!"<<endl;
}

void proveTet18()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","10.89","4","4"};

  interval tz[6]={"5.16653","5.035536",
		  v,"11.56",
		  v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
  
  opt.setDihMax(1.76);

  taylorFunction F=taylorSimplex::y2*"-1"+
                   taylorSimplex::y3*"-1"+
                   taylorSimplex::unit*"4.44";
  
  taylorFunction G=taylorSimplex::dih*"-1"+
                 taylorSimplex::unit*"1.76";
    
  const taylorFunction* K[2]={&F,&G};
  
  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"quad18 complete!"<<endl;
  else
    cout<<"quad18 failed!!"<<endl;
}

void proveTet19()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","11.56","4","4"};

  interval tz[6]={"5.16653","5.035536",
		  v,"12.96",
		  v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::y2*"-1"+
                   taylorSimplex::y3*"-1"+
                   taylorSimplex::unit*"4.64";
  
    opt.setDihMax(1.76);

  //taylorFunction G=taylorSimplex::dih*"-1"+
  //               taylorSimplex::unit*"1.76";
    
  const taylorFunction* K[1]={&F};
  
  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt))
    cout<<"quad19 complete!"<<endl;
  else
    cout<<"quad19 failed!!"<<endl;
}

void proveTet20()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","9.61","4","4"};

  interval tz[6]={"5.16653","5.035536",v,"10.24",v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.049";
  
    opt.setDihMax(1.76);


  //taylorFunction G=taylorSimplex::dih*"-1"+
  //               taylorSimplex::unit*"1.76";
    
  taylorFunction G=taylorSimplex::y2+
                   taylorSimplex::y3+
                   taylorSimplex::unit*"-4.134";

  const taylorFunction* K[2]={&F,&G};
  

 
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"quad20 complete!"<<endl;
  else
    cout<<"quad20 failed!!"<<endl;
} 


void proveTet21()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","10.24","4","4"};

  interval tz[6]={"5.16653","5.035536",v,"10.89",v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.0617";
  
    opt.setDihMax(1.76);


  //taylorFunction G=taylorSimplex::dih*"-1"+
  //               taylorSimplex::unit*"1.76";
    
  taylorFunction G=taylorSimplex::y2+
                   taylorSimplex::y3+
                   taylorSimplex::unit*"-4.29";

  const taylorFunction* K[2]={&F,&G};
  

  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"quad21 complete!"<<endl;
  else
    cout<<"quad21 failed!!"<<endl;
} 

void proveTet22()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","10.89","4","4"};

  interval tz[6]={"5.16653","5.035536",v,"11.56",v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.074";
  
    opt.setDihMax(1.76);

  //taylorFunction G=taylorSimplex::dih*"-1"+
  //               taylorSimplex::unit*"1.76";
    
  taylorFunction G=taylorSimplex::y2+
                   taylorSimplex::y3+
                   taylorSimplex::unit*"-4.44";

  const taylorFunction* K[2]={&F,&G};
  

    
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"quad22 complete!"<<endl;
  else
    cout<<"quad22 failed!!"<<endl;
} 


void proveTet23()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","11.56","4","4"};

  interval tz[6]={"5.16653","5.035536",v,"12.96",v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.093";
  
  opt.setDihMax(1.76);

  //taylorFunction G=taylorSimplex::dih*"-1"+
  //               taylorSimplex::unit*"1.76";
    
  taylorFunction G=taylorSimplex::y2+
                   taylorSimplex::y3+
                   taylorSimplex::unit*"-4.64";

  const taylorFunction* K[2]={&F,&G};
  

  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"quad23 complete!"<<endl;
  else
    cout<<"quad23 failed!!"<<endl;
} 


void proveTet24()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","9.61","4","4"};

  interval tz[6]={v,"5.035536",v,"10.24",v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.024";
  
  opt.setDihMax(2.11);

  taylorFunction G=taylorSimplex::dih*"-1"+
                 taylorSimplex::unit*"2.11";
    
  taylorFunction H=taylorSimplex::y2+
                   taylorSimplex::y3+
                   taylorSimplex::unit*"-4.134";

  const taylorFunction* K[3]={&F,&G,&H};
  

 
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,3,opt))
    cout<<"quad24 complete!"<<endl;
  else
    cout<<"quad24 failed!!"<<endl;
} 

void proveTet25()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","10.24","4","4"};

  interval tz[6]={v,"5.035536",v,"10.89",v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.02";
  

  opt.setDihMax(2.11);

  taylorFunction G=taylorSimplex::dih*"-1"+
                 taylorSimplex::unit*"2.11";
    
  taylorFunction H=taylorSimplex::y2+
                   taylorSimplex::y3+
                   taylorSimplex::unit*"-4.29";

  const taylorFunction* K[3]={&F,&G,&H};
  

 
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,3,opt))
    cout<<"quad25 complete!"<<endl;
  else
    cout<<"quad25 failed!!"<<endl;
} 


void proveTet26()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","10.89","4","4"};

  interval tz[6]={v,"5.035536",v,"11.56",v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.02";
  

  opt.setDihMax(2.11);

 taylorFunction G=taylorSimplex::dih*"-1"+
                 taylorSimplex::unit*"2.11";
    
  taylorFunction H=taylorSimplex::y2+
                   taylorSimplex::y3+
                   taylorSimplex::unit*"-4.44";

  const taylorFunction* K[3]={&F,&G,&H};
  

    
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,3,opt))
    cout<<"quad26 complete!"<<endl;
  else
    cout<<"quad26 failed!!"<<endl;
} 

void proveTet27()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","9.61","4","4"};

  interval tz[6]={"5.166529","5.035536",v,"10.24",v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.074";
  
  opt.setDihMax(1.625);

  //taylorFunction G=taylorSimplex::dih*"-1"+
  //               taylorSimplex::unit*"1.625";
    
  taylorFunction G=taylorSimplex::y2+
                   taylorSimplex::y3+
                   taylorSimplex::unit*"-4.134";

  const taylorFunction* K[2]={&F,&G};
  

   
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"quad27 complete!"<<endl;
  else
    cout<<"quad27 failed!!"<<endl;
} 

void proveTet28()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","11.56","4","4"};

  interval tz[6]={"5.166529","5.035536",v,"12.96",v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"-0.02";
  
  opt.setDihMax(2.419);

  //taylorFunction G=taylorSimplex::dih*"-1"+
  //               taylorSimplex::unit*"1.625";
    
  taylorFunction G=taylorSimplex::y2+
                   taylorSimplex::y3+
                   taylorSimplex::unit*"-4.8";

  const taylorFunction* K[2]={&F,&G};
  

   
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"quad28 complete!"<<endl;
  else
    cout<<"quad28 failed!!"<<endl;
} 


void proveTet29()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","11.56","4","4"};

  interval tz[6]={"5.166529","5.035536",v,"12.96",v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"-0.01";
  
  opt.setDihMax(2.11);

  //taylorFunction G=taylorSimplex::dih*"-1"+
  //               taylorSimplex::unit*"1.625";
    
  taylorFunction G=taylorSimplex::y2+
                   taylorSimplex::y3+
                   taylorSimplex::unit*"-4.64";

  const taylorFunction* K[2]={&F,&G};
  

   
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"quad29 complete!"<<endl;
  else
    cout<<"quad29 failed!!"<<endl;
} 

