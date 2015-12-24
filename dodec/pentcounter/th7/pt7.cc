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
void vol1();
void vol2();
void diag2();
void mu1();
void mu2();
void mu3();
void diag3();
void diag4();
void vol3();
void vol4();
void vol5();
void mu4();

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
  vol1();  
  vol2();
  diag2();  
  mu1();
  mu2();
  mu3();
  diag3();  
  diag4();  
  vol3();
  vol4();
  vol5();
  mu4();

  return 0;
}

//*****************************************************************
 
  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  //  dih < 1.735, y1 < 2.199, y2 < 2.26, y3 < 2.206 -> y4 < 3.254
  //
  /////////////////////////////////////////////////////////////////
  
void diag1()
{
  cellOption opt;
  domain x,z,x0,z0;  

  interval tx[6]={"4","4","4",v,"4","4"};
  interval tz[6]={"4.835601","5.1076","4.866436","12",v,v};

  // 2.199^2=4.835601, 2.26^2=5.1076, 2.206^2=4.866436

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);  
  
  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-3.254";

  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.735";

  opt.setDihMax(1.735);

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
  //  dih < 1.735, y1 < 2.199, y2 < 2.26, y3 < 2.206,y4 < 3.254
  //
  //  -> vol > .275
  //
  /////////////////////////////////////////////////////////////////
  
void vol1()
{
  cellOption opt;
  domain x,z,x0,z0;

  interval tx[6]={"4","4","4",v,"4","4"};
  interval tz[6]={"4.835601","5.1076","4.866436","10.588516",v,v};

  // 2.199^2=4.835601, 2.26^2=5.1076, 2.206^2=4.866436, 
  // 3.254^2=10.588516
  
  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.275";

  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.735";

  opt.setDihMax(1.735);
    
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
  //  y1 < 2.206, y2,y3 < 2.26, y5,y6 < 3.254
  //
  //  -> vol > .357
  //
  /////////////////////////////////////////////////////////////////
   
void vol2()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","4",v,v};
  interval tz[6]={"4.866436","5.1076","5.1076",v,"10.588516","10.588516"};

  // 2.206^2=4.866436, 2.26^2=5.1076, 3.254^2=10.588516, 
 
  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.357";
  
  const taylorFunction* K[1]={&F};  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt))
    cout<<"vol2 complete!"<<endl;
  else
    cout<<"vol2 failed!!"<<endl;
}

//*****************************************************************
 
  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  //  dih < 1.735, y1 < 2.199, y2 > 2.26, y3 < 2.206 -> y4 < 3.396
  //
  /////////////////////////////////////////////////////////////////
  
void diag2()
{
  cellOption opt;
  domain x,z,x0,z0;  

  interval tx[6]={"4","5.1076","4",v,"4","4"};
  interval tz[6]={"4.835601",v,"4.866436","12",v,v};

  // 2.199^2=4.835601, 2.26^2=5.1076, 2.206^2=4.866436

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);  
  
  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-3.396";

  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.735";

  opt.setDihMax(1.735);

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
  //  dih < 1.735, y1 < 2.199, y2 < 2.26, y3 < 2.206,y4 < 3.254
  //
  //  -> mu > -.0027
  //
  /////////////////////////////////////////////////////////////////
  
void mu1()
{
  cellOption opt;
  domain x,z,x0,z0;

  interval tx[6]={"4","4","4",v,"4","4"};
  interval tz[6]={"4.835601","5.1076","4.866436","10.588516",v,v};

  // 2.199^2=4.835601, 2.26^2=5.1076, 2.206^2=4.866436, 
  // 3.254^2=10.588516
  
  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"-0.0027";

  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.735";

  opt.setDihMax(1.735);
    
  const taylorFunction* K[2]={&F,&G};  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"mu1 complete!"<<endl;
  else
    cout<<"mu1 failed!!"<<endl;
}

//*****************************************************************
 
  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  //  dih < 1.735, y1 < 2.199, y2 > 2.26, y3 < 2.206,y4 < 3.396
  //
  //  -> mu > .0313
  //
  /////////////////////////////////////////////////////////////////
  
void mu2()
{
  cellOption opt;
  domain x,z,x0,z0;

  interval tx[6]={"4","5.1076","4",v,"4","4"};
  interval tz[6]={"4.835601",v,"4.866436","11.532816",v,v};

  // 2.199^2=4.835601, 2.26^2=5.1076, 2.206^2=4.866436, 
  // 3.396^2=11.532816
  
  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.0313";

  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.735";

  opt.setDihMax(1.735);
    
  const taylorFunction* K[2]={&F,&G};  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"mu2 complete!"<<endl;
  else
    cout<<"mu2 failed!!"<<endl;
}



//*****************************************************************
 
  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  //  y1 < 2.206, y2 < 2.26, y3 > 2.26, y5 < 3.396, y6 < 3.254
  //
  //  -> mu > .0627
  //
  /////////////////////////////////////////////////////////////////
   
void mu3()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","5.1076","4",v,v};
  interval tz[6]={"4.866436","5.1076",v,v,"11.532816","10.588516"};

  // 2.206^2=4.866436, 2.26^2=5.1076, 3.254^2=10.588516,
  // 3.396^2=10.588516 
 
  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.0627";
  
  const taylorFunction* K[1]={&F};  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt))
    cout<<"mu3 complete!"<<endl;
  else
    cout<<"mu3 failed!!"<<endl;
}

//*****************************************************************
 
  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  //  dih < 1.574, y1 < 2.199, y2 < 2.26, y3 < 2.206 -> y4 < 3.021
  //
  /////////////////////////////////////////////////////////////////
  
void diag3()
{
  cellOption opt;
  domain x,z,x0,z0;  

  interval tx[6]={"4","4","4",v,"4","4"};
  interval tz[6]={"4.835601","5.1076","4.866436","12",v,v};

  // 2.199^2=4.835601, 2.26^2=5.1076, 2.206^2=4.866436

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);  
  
  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-3.021";

  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.574";

  opt.setDihMax(1.574);

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
  //  dih < 1.574, y1 < 2.199, y2 > 2.26, y3 < 2.206 -> y4 < 3.157
  //
  /////////////////////////////////////////////////////////////////
  
void diag4()
{
  cellOption opt;
  domain x,z,x0,z0;  

  interval tx[6]={"4","5.1076","4",v,"4","4"};
  interval tz[6]={"4.835601",v,"4.866436","12",v,v};

  // 2.199^2=4.835601, 2.26^2=5.1076, 2.206^2=4.866436

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);  
  
  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-3.157";

  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.574";

  opt.setDihMax(1.574);

  const taylorFunction* K[2]={&F,&G};
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt)) 
    cout<<"diag4 complete!"<<endl;
  else
    cout<<"diag4 failed!!"<<endl;
}


//*****************************************************************
 
  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  //  dih < 1.574, y1 < 2.199, y2 < 2.26, y3 < 2.206,y4 < 3.021
  //
  //  -> vol > .275
  //
  /////////////////////////////////////////////////////////////////
  
void vol3()
{
  cellOption opt;
  domain x,z,x0,z0;

  interval tx[6]={"4","4","4",v,"4","4"};
  interval tz[6]={"4.835601","5.1076","4.866436","9.126441",v,v};

  // 2.199^2=4.835601, 2.26^2=5.1076, 2.206^2=4.866436, 
  // 3.021^2=9.126441
  
  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.275";

  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.574";

  opt.setDihMax(1.574);
    
  const taylorFunction* K[2]={&F,&G};  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"vol3 complete!"<<endl;
  else
    cout<<"vol3 failed!!"<<endl;
}


//*****************************************************************
 
  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  //  dih < 1.574, y1 < 2.199, y2 > 2.26, y3 < 2.206,y4 < 3.157
  //
  //  -> vol > .249
  //
  /////////////////////////////////////////////////////////////////
  
void vol4()
{
  cellOption opt;
  domain x,z,x0,z0;

  interval tx[6]={"4","4","4",v,"4","4"};
  interval tz[6]={"4.835601","5.1076","4.866436","9.966649",v,v};

  // 2.199^2=4.835601, 2.26^2=5.1076, 2.206^2=4.866436, 
  // 3.157^2=9.966649
  
  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.249";

  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.574";

  opt.setDihMax(1.574);
    
  const taylorFunction* K[2]={&F,&G};  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"vol4 complete!"<<endl;
  else
    cout<<"vol4 failed!!"<<endl;
}

//*****************************************************************
 
  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  //  y1 < 2.206, y2 < 2.26, y3 > 2.26, y5 < 3.157,
  //
  //  y6 < 3.021 -> vol > .321
  //
  /////////////////////////////////////////////////////////////////
   
void vol5()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","5.1076","4",v,v};
  interval tz[6]={"4.866436","5.1076",v,v,"9.966649","9.126441"};

  // 2.206^2=4.866436, 2.26^2=5.1076, 3.157^2=9.966649,  
  // 3.021^2=9.126441
 
  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.321";
  
  const taylorFunction* K[1]={&F};  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt))
    cout<<"vol5 complete!"<<endl;
  else
    cout<<"vol5 failed!!"<<endl;
}

//*****************************************************************
 
  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  //  y1 < 2.206, y2,y3 > 2.26, y5,y6 < 3.396
  //
  //  -> mu > .0336
  //
  /////////////////////////////////////////////////////////////////
   
void mu4()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","5.1076","5.1076","4",v,v};
  interval tz[6]={"4.866436",v,v,v,"10.588516","10.588516"};

  // 2.206^2=4.866436, 2.26^2=5.1076,3.396^2=10.588516 
 
  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.0336";
  
  const taylorFunction* K[1]={&F};  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt))
    cout<<"mu4 complete!"<<endl;
  else
    cout<<"mu4 failed!!"<<endl;
}

