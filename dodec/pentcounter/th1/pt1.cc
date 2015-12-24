//  This file contains the proof of all tetrahedral inequalities 
//  needed for Pentagon Theorem 1 of the proof of the dodecahedral
//  conjecture.


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
void diag3();
void diag4();
void vol1();
void vol2();
void vol3();
void diag5();
void vol4();
void mu1();
void mu2();
void mu3();
void diag6();
void mu4();

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
  diag1();
  diag2();  
  diag3();  
  diag4();
  vol1();
  vol2();
  vol3();
  diag5();
  vol4();
  mu1();
  mu2();
  mu3();
  diag6();
  mu4();

  return 0;
}

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  // dih<1.701, y1<2.219, y2<2.330,y3<2.207, y1+y2+y3<6.406
  //
  // y2+y3<4.399 -> y4<3.179
  //
  //  
  /////////////////////////////////////////////////////////////////

 
void diag1()
{
  cellOption opt;
  domain x,z,x0,z0;
  interval tx[6]={"4","4","4",v,"4","4"};
  interval tz[6]={"4.923961","5.4289","4.870849","12",v,v};

  // 2.219^2=4.923961, 2.33^2=5.4289, 2.207^2=4.870849

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
    cout<<"diag1 complete!"<<endl;
  else
    cout<<"diag1 failed!!"<<endl;

}

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  // dih<1.681, y1<2.161, y2<2.207,y3<2.280, y1+y2+y3<6.313
  //
  // y2+y3<4.307 -> y4<3.101
  //
  //  
  /////////////////////////////////////////////////////////////////


void diag2()
{
  cellOption opt;
  domain x,z,x0,z0;

  interval tx[6]={"4","4","4",v,"4","4"};
  interval tz[6]={"4.669921","4.870849","5.1984","12",v,v};
  
  // 2.161^2=4.669921,2.207^2=4.870849, 2.28^2=5.1984

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

  taylorFunction I=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.681";

  opt.setDihMax(1.681);
 
  const taylorFunction* K[4]={&F,&G,&H,&I};  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,4,opt))
    cout<<"diag2 complete!"<<endl;
  else
    cout<<"diag2 failed!!"<<endl;
}

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  // dih<1.917, y1<2.329, y2<2.207,y3<2T, y1+y2+y3<6.757
  //
  // y2+y3<4.657 -> y4<3.611
  //
  //  
  /////////////////////////////////////////////////////////////////


void diag3()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4",v,"4","4"};
  interval tz[6]={"5.424241","4.870849",v,"14",v,v};

  // 2.329^2=4.669921,2.207^2=4.870849, 

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-3.611";

  taylorFunction G=taylorSimplex::y2*"-1"+
                   taylorSimplex::y3*"-1";
                   taylorSimplex::unit*"4.657";

  taylorFunction H=taylorSimplex::y2*"-1"+
                   taylorSimplex::y3*"-1"+
                   taylorSimplex::y1*"-1"+
                   taylorSimplex::unit*"6.757";

  taylorFunction I=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.917";  

  opt.setDihMax(1.917);
    
  const taylorFunction* K[4]={&F,&G,&H,&I};
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,4,opt))
    cout<<"diag3 complete!"<<endl;
  else
    cout<<"diag3 failed!!"<<endl;
}

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  // dih<1.917, y1<2.329, y2<2.207,y3<2.15, y1+y2+y3<6.757
  //
  // y2+y3<4.657 -> y4<3.3
  //
  //  
  /////////////////////////////////////////////////////////////////

void diag4()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4",v,"4","4"};
  interval tz[6]={"4.669921","4.870849","4.6225","11",v,v};

  // 2.329^2=4.669921,2.207^2=4.870849,2.15^2= 4.6225

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-3.427";

  taylorFunction G=taylorSimplex::y2*"-1"+
                   taylorSimplex::y3*"-1"+             
                   taylorSimplex::unit*"4.657";

  taylorFunction H=taylorSimplex::y2*"-1"+
                   taylorSimplex::y3*"-1"+
                   taylorSimplex::y1*"-1"+
                   taylorSimplex::unit*"6.757";

  taylorFunction I=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.917";  

  opt.setDihMax(1.917);
    
  const taylorFunction* K[4]={&F,&G,&H,&I};  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,4,opt))
    cout<<"diag4 complete!"<<endl;
  else
    cout<<"diag4 failed!!"<<endl;
}


//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  // dih<1.681, y1<2.161, y2<2.207,y3<2.28, y1+y2+y3<6.313
  //
  // y2+y3<4.307,y4<3.101 -> vol(S)>0.273
  //
  //  
  /////////////////////////////////////////////////////////////////


void vol1()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4",v,"4","4"};

  interval tz[6]={"4.669921","4.870849","5.1984","9.616201",v,v};
  
  // 2.161^2=4.669921,2.207^2=4.870849, 2.28^2=5.1984, 3.101^2=9.616201
  
  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.273";

  taylorFunction G=taylorSimplex::y2*"-1"+
                   taylorSimplex::y3*"-1"+             
                   taylorSimplex::unit*"4.307";

  taylorFunction H=taylorSimplex::y2*"-1"+
                   taylorSimplex::y3*"-1"+
                   taylorSimplex::y1*"-1"+
                   taylorSimplex::unit*"6.313";
  
  taylorFunction I=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.681";

  opt.setDihMax(1.681);
    
  const taylorFunction* K[4]={&F,&G,&H,&I};  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,4,opt))
    cout<<"vol1 complete!"<<endl;
  else
    cout<<"vol1 failed!!"<<endl;
}


//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  // dih<1.917, y1<2.329, y2<2.207,y3<2.15, y1+y2+y3<6.757
  //
  // y2+y3<4.657,y4<3.427 -> vol(S)>0.270
  //
  //  
  /////////////////////////////////////////////////////////////////

void vol2()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4",v,"4","4"};
  interval tz[6]={"5.424241","4.870849","4.6225","11.744329",v,v};

  // 2.329^2=5.424241,2.207^2=4.870849,2.15^2=4.6225,3.427^2=11.744329 

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.270"; 

  taylorFunction G=taylorSimplex::y2*"-1"+
                   taylorSimplex::y3*"-1"+             
                   taylorSimplex::unit*"4.657";

  taylorFunction H=taylorSimplex::y2*"-1"+
                   taylorSimplex::y3*"-1"+
                   taylorSimplex::y1*"-1"+
                   taylorSimplex::unit*"6.757";
  
  taylorFunction I=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.917";
    
  opt.setDihMax(1.917);
      
  const taylorFunction* K[4]={&F,&G,&H,&I};  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,4,opt))
    cout<<"vol2 complete!"<<endl;
  else
    cout<<"vol2 failed!!"<<endl;
}


//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  // y1<2.207, y2<2.28,y3<2.15, 
  //
  // 2T<y5<3.427, 2T<y6<3.101 -> vol(S)>0.328
  //
  //  
  /////////////////////////////////////////////////////////////////


void vol3()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","4",v,v};
  interval tz[6]={"4.870849","5.1984","4.6225",v,"11.744329","9.616201"};

  // 2.207^2=4.870849,2.28^2=5.1984, 2.15^2=4.6225,3.427^2=11.744329  
  // 3.101^2=9.616201

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.328";
  
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
  // dih<1.817, y1<2.329, y2<2.207,y3<2.15,-> y4<3.303
  //
  //  
  /////////////////////////////////////////////////////////////////

void diag5()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4",v,"4","4"};
  interval tz[6]={"4.669921","4.870849","4.6225","11",v,v};

  // 2.329^2=4.669921,2.207^2=4.870849,2.15^2= 4.6225

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-3.427";

  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.817";  

  opt.setDihMax(1.817);
    
  const taylorFunction* K[2]={&F,&G};  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"diag5 complete!"<<endl;
  else
    cout<<"diag5 failed!!"<<endl;
}


//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  // y1<2.207, y2<2.28,y3<2.15, 
  //
  // 2T<y5<3.303, 2T<y6<3.101 -> vol(S)>0.352
  //
  //  
  /////////////////////////////////////////////////////////////////

void vol4()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","4",v,v};
  interval tz[6]={"4.870849","5.1984","4.6225",v,"10.909809","9.616201"};

  // 2.207^2=4.870849,2.28^2=5.1984, 2.15^2=4.6225,3.303^2=10.909809  
  // 3.101^2=9.616201

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::unit*"0.352";
  
  const taylorFunction* K[1]={&F};  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt))
    cout<<"vol4 complete!"<<endl;
  else
    cout<<"vol4 failed!!"<<endl;
}

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  // dih<1.681, y1<2.161, y2<2.207,y3<2.28, y1+y2+y3<6.313
  //
  // y2+y3<4.307,y4<3.101 -> mu(S)>0.156
  //
  //  
  /////////////////////////////////////////////////////////////////


void mu1()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4",v,"4","4"};

  interval tz[6]={"4.669921","4.870849","5.1984","9.616201",v,v};
  
  // 2.161^2=4.669921,2.207^2=4.870849, 2.28^2=5.1984, 3.101^2=9.616201
  
  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.0156";

  taylorFunction G=taylorSimplex::y2*"-1"+
                   taylorSimplex::y3*"-1"+             
                   taylorSimplex::unit*"4.307";

  taylorFunction H=taylorSimplex::y2*"-1"+
                   taylorSimplex::y3*"-1"+
                   taylorSimplex::y1*"-1"+
                   taylorSimplex::unit*"6.313";
  
  taylorFunction I=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.681";

  opt.setDihMax(1.681);
    
  const taylorFunction* K[4]={&F,&G,&H,&I};  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,4,opt))
    cout<<"mu1 complete!"<<endl;
  else
    cout<<"mu1 failed!!"<<endl;
}


//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  // dih<1.917, y1<2.329, y2<2.207,y3>2.15, y1+y2+y3<6.757
  //
  // y2+y3<4.657,y4<3.611 -> mu(S)>0.0263
  //
  //  
  /////////////////////////////////////////////////////////////////

void mu2()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4.6225",v,"4","4"};
  interval tz[6]={"5.424241","4.870849",v,"13.039321",v,v};

  // 2.329^2=5.424241,2.207^2=4.870849,2.15^2=4.6225,3.611^2=13.039321 

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+  
                   taylorSimplex::unit*"0.0263"; 

  taylorFunction G=taylorSimplex::y2*"-1"+
                   taylorSimplex::y3*"-1"+             
                   taylorSimplex::unit*"4.657";

  taylorFunction H=taylorSimplex::y2*"-1"+
                   taylorSimplex::y3*"-1"+
                   taylorSimplex::y1*"-1"+
                   taylorSimplex::unit*"6.757";
  
  taylorFunction I=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.917";
    
  opt.setDihMax(1.917);
      
  const taylorFunction* K[4]={&F,&G,&H,&I};  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,4,opt))
    cout<<"mu2 complete!"<<endl;
  else
    cout<<"mu2 failed!!"<<endl;
}


//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  // y1<2.207, y2<2.28,y3>2.15, 
  //
  // 2T<y5<3.611, 2T<y6<3.101 -> mu(S)>0.0546
  //
  //  
  /////////////////////////////////////////////////////////////////


void mu3()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4.6225","4",v,v};
  interval tz[6]={"4.870849","5.1984",v,v,"13.039321","9.616201"};

  // 2.207^2=4.870849,2.28^2=5.1984, 2.15^2=4.6225,3.611^2=13.039321  
  // 3.101^2=9.616201

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.0546";
  
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
  // dih<1.57, y1<2.329, y2<2.207,y3<2.15,-> y4<3.2
  //
  //  
  /////////////////////////////////////////////////////////////////

void diag6()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4",v,"4","4"};
  interval tz[6]={"4.669921","4.870849","4.6225","11",v,v};

  // 2.329^2=4.669921,2.207^2=4.870849,2.15^2= 4.6225

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-3.2";

  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.57";  

  opt.setDihMax(1.57);
    
  const taylorFunction* K[2]={&F,&G};  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"diag6 complete!"<<endl;
  else
    cout<<"diag6 failed!!"<<endl;
}

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  // y1<2.207, y2<2.28,y3>2.15, 
  //
  // 2T<y5<3.2, 2T<y6<3.101 -> mu(S)>0.057
  //
  //  
  /////////////////////////////////////////////////////////////////


void mu4()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4.6225","4",v,v};
  interval tz[6]={"4.870849","5.1984",v,v,"10.24","9.616201"};

  // 2.207^2=4.870849,2.28^2=5.1984, 2.15^2=4.6225,3.2^2=10.24
  // 3.101^2=9.616201

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.057";
  
  const taylorFunction* K[1]={&F};  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt))
    cout<<"mu6 complete!"<<endl;
  else
    cout<<"mu6 failed!!"<<endl;
}
