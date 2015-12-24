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
void vol1();
void vol2();
void vol3();
void diag3();
void mu1();
void mu2();
void mu3();

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
  vol1();
  vol2();
  vol3();
  diag3();
  mu1();
  mu2();
  mu3();

  return 0;
}

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  //  dih < 1.342, y1 < 2.194, y2 < 2.314, y3 < 2.26  -> y4 < 2.699
  //
  /////////////////////////////////////////////////////////////////
  
void diag1()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4",v,"4","4"};
  interval tz[6]={"4.813636","5.354596","5.1076","12",v,v};

  // 2.194^2=4.813636, 2.314^2=5.354596, 12 > 2.699^2 = 7.284601

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
  
  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-2.809";

  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.342";

  opt.setDihMax(1.342);

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
  //  dih < 1.684, y1 < 2.153, y2 < 2.314, y3 < 2.174 -> y4 < 3.196
  //
  /////////////////////////////////////////////////////////////////
 
void diag2()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4",v,"4","4"};

  interval tz[6]={"4.635409","5.354596","4.726276","12",v,v};

  // 2.153^2=4.635409, 2.314^2=5.354596, 2.174^2=4.726276, 
  // 12 > 3.196^2 = 10.214416

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-3.196";

  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.684";  

  opt.setDihMax(1.684);
    
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
  //  dih < 1.342, y1 < 2.194, y2 < 2.314, y3 < 2.26, y4 < 2.699
  //
  //  -> vol > .328
  //
  /////////////////////////////////////////////////////////////////
  
void vol1()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4",v,"4","4"};
  interval tz[6]={"4.813636","5.354596","5.1076","7.284601",v,v};

  // 2.194^2=4.813636, 2.314^2=5.354596, 2.26^2=5.1076,
  // 2.699^2 = 7.284601

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.328";

  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.342";
  
  opt.setDihMax(1.342);
    
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
  //  dih < 1.684, y1 < 2.153, y2 < 2.314, y3 < 2.174, y4 < 3.196
  //
  //  -> vol > .272
  //
  /////////////////////////////////////////////////////////////////
 
void vol2()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4",v,"4","4"};
  interval tz[6]={"4.635409","5.354596","4.726276","10.214416",v,v};

  // 2.153^2=4.635409, 2.314^2=5.354596, 2.174^2=4.726276, 
  // 3.196^2 = 10.214416

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.272";
  
  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.684";

  opt.setDihMax(1.684);
    
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
  //  y1 < 2.314, y2 < 2.174, y3 < 2.26, y5 < 2.699^2, y6 < 3.196 
  //
  //  -> vol > .35
  //
  /////////////////////////////////////////////////////////////////
 
void vol3()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","4",v,v};
  interval tz[6]={"5.354596","4.726276","5.1076",v,"7.284601","10.214416"};

  // 2.314^2=5.354596, 2.174^2=4.726276, 2.26^2=5.1076,
  // 3.196^2 = 10.214416,  2.699^2 = 7.284601

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.35"; 
    
  const taylorFunction* K[1]={&F};  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt))
    cout<<"vol3 complete!"<<endl;
  else
    cout<<"vol3 failed!!"<<endl;
}

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  //  dih < 1.342, y1 < 2.194, y2 < 2.314, y3 > 2.26 -> y4 < 2.809
  //
  /////////////////////////////////////////////////////////////////
  
void diag3()
{
  cellOption opt;
  domain x,z,x0,z0;  

  interval tx[6]={"4","4","5.1076",v,"4","4"};
  interval tz[6]={"4.813636","5.354596",v,"12",v,v};

  // 2.194^2=4.813636, 2.314^2=5.354596, 2.26^2=5.1076,
  // 3.196^2 = 10.214416, 12 > 2.809^2=7.890481
  
  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);  
  
  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-2.809";

  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.342";

  opt.setDihMax(1.342);

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
  //  dih < 1.342, y1 < 2.194, y2 < 2.314, y3 > 2.26, y4 < 2.809 
  //
  //  -> mu > .0451
  //
  /////////////////////////////////////////////////////////////////
  
void mu1()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","5.1076",v,"4","4"};
  interval tz[6]={"4.813636","5.354596",v,"7.890481",v,v};

  // 2.194^2=4.813636, 2.314^2=5.354596, 2.26^2=5.1076,
  // 2.809^2 = 7.890481

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.0451";

  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.342";
  
  opt.setDihMax(1.342);
    
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
  //  dih < 1.684, y1 < 2.153, y2 < 2.314, y3 < 2.174, y4 < 3.196
  //
  //  -> mu > .0156
  //
  /////////////////////////////////////////////////////////////////
 
void mu2()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4",v,"4","4"};
  interval tz[6]={"4.635409","5.354596","4.726276","10.214416",v,v};

  // 2.153^2=4.635409, 2.314^2=5.354596, 2.174^2=4.726276, 
  // 3.196^2 = 10.214416

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
 
  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.0156";
  
   taylorFunction G=taylorSimplex::dih*"-1"+
                    taylorSimplex::unit*"1.684";
    
   opt.setDihMax(1.684);

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
  //  y1 < 2.314, y2 < 2.174, y3 > 2.26, y5 < 2.809, y6 < 3.196 
  //
  //  -> mu > .0627
  //
  /////////////////////////////////////////////////////////////////
 
void mu3()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","5.1076","4",v,v};
  interval tz[6]={"5.354596","4.726276",v,v,"7.890481","10.214416"};

  // 2.314^2=5.354596, 2.174^2=4.726276, 2.26^2=5.1076,
  // 3.196^2 = 10.214416, 2.809^2=7.890481
  
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


