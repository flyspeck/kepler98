


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




void pentStar6();
void pentStar7();             
void pentStar8();
void pentStar9();             
void pentStar10();
void pentStar11();


const interval v="6.33436854000504728617";  // = (2 Sqrt[3]Tan[Pi/5])^2
const interval a="5.1004157281"; // =2.25841^2


void selfTest()
    {
    interMath::selfTest();
    linearization::selfTest();
    secondDerive::selfTest();
    taylorFunction::selfTest();
    }


int main()
{

  pentStar6();
  pentStar7();
  pentStar8();
  pentStar9();
  pentStar10();
  pentStar11();
  
  return 0;
}

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  // if y1<2.25841, 3.179>y4>2.95, 2.25841>y7>2 than  
  //    
  // -vol(S) - .166 y1 -.114552 (y2+y3)+.115382 y4 + .143 (y5+y6)
  // -.15342(dih-1.701)+.21 < 0
  //
  // .15342*1.701+.21=0.47096742
  //
  /////////////////////////////////////////////////////////////////



void pentStar6()
{
  cellOption opt;
  domain x,z,x0,z0;

  interval tx[6]={"4","4","4","8.7025","4","4"};
  interval tz[6]={a,v,v,"10.106041",v,v};

  // 10.106041=3.179^2, 8.7025=2.95^2

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
 
  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::y1*"-0.166"+
                   taylorSimplex::y2*"-0.114552"+
                   taylorSimplex::y3*"-0.114552"+ 
                   taylorSimplex::y4*"0.115382"+
                   taylorSimplex::y5*"0.143"+
                   taylorSimplex::y6*"0.143"+
                   taylorSimplex::dih*"-0.153420"+
                   taylorSimplex::unit*"0.47096742"; 
   
  const taylorFunction* K[1]={&F};    

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt)) 
    cout<<"pentStar6 complete!"<<endl;
  else
    cout<<"pentStar6 failed!!"<<endl;

}

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  // if y1<2.25841, y2 <2.33, y3<2.207, 3.179>y4>2.95, 2.25841>y7>2 than  
  //    
  // -vol(S) - .166 y7 -.051448 (y2+y3)-.115382 y4 + .143 (y8+y9)
  // +.56 < 0
  //
  /////////////////////////////////////////////////////////////////

void pentStar7()
{
  cellOption opt;
  domain x,z,x0,z0;

  interval tx[6]={"4","4","4","8.7025","4","4"};
  interval tz[6]={a,"5.4289","4.870849","10.106041",v,v};

  // 10.106041=3.179^2, 7.84=2.8^2, 8.7025=2.95^2, 2.33^2=5.4289
  // 2.207^2=4.870849

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
 
  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::y1*"-0.166"+
                   taylorSimplex::y2*"-0.051448"+
                   taylorSimplex::y3*"-0.051448"+ 
                   taylorSimplex::y4*"-0.115382"+
                   taylorSimplex::y5*"0.143"+
                   taylorSimplex::y6*"0.143"+
                   taylorSimplex::unit*"0.56"; 
   
  const taylorFunction* K[1]={&F};    

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt)) 
    cout<<"pentStar7 complete!"<<endl;
  else
    cout<<"pentStar7 failed!!"<<endl;

}

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  // if y1<2.25841, 3.179>y4>2.95, 2T>y7>2.25841, 2.25841>y2>2 than  
  //    
  // -vol(S) - .166 y1 -.179514 y2 - .257750 y3 +.169516 y4 + .143 (y5+y6)
  // -.273372 (dih-1.701)+0.491150 < 0
  //
  // .273372*1.701+.49115=.956155772
  //
  /////////////////////////////////////////////////////////////////

void pentStar8()
{
  cellOption opt;
  domain x,z,x0,z0;

  interval tx[6]={"4","4","4","8.7025","4","4"};
  interval tz[6]={a,a,v,"10.106041",v,v};

  //10.106041=3.179^2, 8.7025=2.95^2

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
 
  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::y1*"-0.166"+
                   taylorSimplex::y2*"-0.179514"+
                   taylorSimplex::y3*"-0.257750"+ 
                   taylorSimplex::y4*"0.169516"+
                   taylorSimplex::y5*"0.143"+
                   taylorSimplex::y6*"0.143"+
                   taylorSimplex::dih*"-0.273372"+
                   taylorSimplex::unit*"0.956155772"; 
   
  const taylorFunction* K[1]={&F};    

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt)) 
    cout<<"pentStar8 complete!"<<endl;
  else
    cout<<"pentStar8 failed!!"<<endl;

}

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  // if y1<2.25841, 3.179>y4>2.95, 2T>y7>2.25841, 2.25841>y2>2 than 
  //    
  // -vol(S) - .166 y7 +.013514 y2 + .09175 y3 -.169516 y4 + .143 (y8+y9)
  // +.278850 < 0
  //
  /////////////////////////////////////////////////////////////////

void pentStar9()
{
  cellOption opt;
  domain x,z,x0,z0;

  interval tx[6]={a,"4","4","8.7025","4","4"};

  interval tz[6]={v,a,"4.870849","10.106041",v,v};

  // 10.106041=3.179^2, 8.7025=2.95^2, 2.33^2=5.4289
  // 2.207^2=4.870849



  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
 
  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::y1*"-0.166"+
                   taylorSimplex::y2*"0.013514"+
                   taylorSimplex::y3*"0.091750"+ 
                   taylorSimplex::y4*"-0.169516"+
                   taylorSimplex::y5*"0.143"+
                   taylorSimplex::y6*"0.143"+
                   taylorSimplex::unit*"0.278850"; 

  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"2.652";

  opt.setDihMax(2.652);
 
  const taylorFunction* K[2]={&F,&G};    

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt)) 
    cout<<"pentStar9 complete!"<<endl;
  else
    cout<<"pentStar9 failed!!"<<endl;

}

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  // if y1<2.25841, 3.179>y4>2.95, 2T>y7>2.25841, 2T>y2>2.25841 than  
  //    
  // -vol(S) - .166 y1 -.195794 y2 - .147301 y3 +.227016 y4 + .143 (y5+y6)
  // -.206212 (dih-1.701)+.107554 < 0
  //
  // .206213*1.701+.107554=.458322313
  //
  /////////////////////////////////////////////////////////////////

void pentStar10()
{
  cellOption opt;
  domain x,z,x0,z0;

  interval tx[6]={"4",a,"4","8.7025","4","4"};
  interval tz[6]={a,v,v,"10.106041",v,v};

  // 10.106041=3.179^2,  8.7025=2.95^2

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
 
  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::y1*"-0.166"+
                   taylorSimplex::y2*"-0.195794"+
                   taylorSimplex::y3*"-0.147301"+ 
                   taylorSimplex::y4*"0.227016"+
                   taylorSimplex::y5*"0.143"+
                   taylorSimplex::y6*"0.143"+
                   taylorSimplex::dih*"-0.206212"+
                   taylorSimplex::unit*".458322313"; 
 
  
  const taylorFunction* K[1]={&F};    

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt)) 
    cout<<"pentStar10 complete!"<<endl;
  else
    cout<<"pentStar10 failed!!"<<endl;

}

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  // 
  // if y1<2.25841, 3.179>y4>2.95, 2T>y7>2.25841, 2T>y2>2.25841 than  
  // 
  // -vol(S) - .166 y7 +.029794 y2 -.018699 y3 -.227016 y4 + .143 (y8+y9)
  // +.662446 < 0
  //
  /////////////////////////////////////////////////////////////////

void pentStar11()
{
  cellOption opt;
  domain x,z,x0,z0;

  interval tx[6]={a,a,"4","8.7025","4","4"};
  interval tz[6]={v,v,v,"10.106041",v,v};

  //10.106041=3.179^2, 8.7025=2.95^2

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
 
  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::y1*"-0.166"+
                   taylorSimplex::y2*"0.029794"+
                   taylorSimplex::y3*"-0.018699"+ 
                   taylorSimplex::y4*"-0.227016"+
                   taylorSimplex::y5*"0.143"+
                   taylorSimplex::y6*"0.143"+
                   taylorSimplex::unit*"0.662446"; 
   
  const taylorFunction* K[1]={&F};    

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt)) 
    cout<<"pentStar11 complete!"<<endl;
  else
    cout<<"pentStar11 failed!!"<<endl;

}

