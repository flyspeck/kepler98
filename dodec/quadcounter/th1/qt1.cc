

// This file includes all simplex verifications needed for the proof of 
// Theorem 1 of Appendix 2 of the proof of the dodecahedral conjecture. *)


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



void star1();
void star2();
void star2();
void star3();
void star4();
void star5();
void star6();
void star7();             
void star8();
void star9();             
void star10();
void star11();
void long1();
void long2();
void table1();
void table2();
void table3();
void table4();
void table5();
void table6();
void table7();
void table8();
void table9();
void table10();
void table11();
void table12();
void spec1();
void spec2();
void tableA1();
void tableA2();
void tableA3();
void tableA4();
void tableA5();
void tableA6();
void tableA7();
void tableA8();
void tableA9();
void tableA10();
void tableA11();
void tableA12();
void specA1();


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
 

  star1();
  star2();
  star3();
  star4();
  star5();
  star6();
  star7();
  star8();
  star9();
  star10();
  star11();
  long1();
  long2();
  table1();
  table2();
  table3();
  table4();
  table5();
  table6();
  table7();
  table8();
  table9();
  table10();
  table11();
  table12();
  spec1();
  spec2();
  tableA1();
  tableA2();
  tableA3();
  tableA4();
  tableA5();
  tableA6();
  tableA7();
  tableA8();
  tableA9();
  tableA10();
  tableA11();
  tableA12();
  specA1();

  return 0;
}



//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  // -vol(S)-.166(y1-2)-0.083(y2+y3-4)+.143(y5+y6-4)+.590491/2 < 0
  //
  //  2 * .166 + 4 * .083 - 4 * .143 +.2952455 = .387246
  //
  /////////////////////////////////////////////////////////////////

void star1()
{
  cellOption opt;
  domain x,z,x0,z0;
    
  interval tx[6]={"4","4","4",v,"4","4"};
  interval tz[6]={v,v,v,"7.84",v,v};

  // 7.84=2.8^2

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);  
   
  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::y1*"-0.166"+
                   taylorSimplex::y2*"-0.083"+
                   taylorSimplex::y3*"-0.083"+
                   taylorSimplex::y5*"0.143"+
                   taylorSimplex::y6*"0.143"+
                   taylorSimplex::unit*"0.387245"; 
   
  const taylorFunction* K[1]={&F};  

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt)) 
    cout<<"star1 complete!"<<endl;
  else
    cout<<"star1 failed!!"<<endl;

}

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  // if y1>2.25841, 3.1>y4>2.8, than  
  //    
  // -vol(S) - .166 y1 - 0.174893 (y2+y3) + .306136 y4 + .143 (y5+y6)
  // -.263488 (dih-1.79)-.038250 < 0
  //
  //  1.719*.263488-.038250=.414686
  //
  /////////////////////////////////////////////////////////////////

void star2()
{
  cellOption opt;
  domain x,z,x0,z0;
    
  interval tx[6]={a,"4","4","7.84","4","4"};
  interval tz[6]={v,v,v,"9.61",v,v};

  // 7.84=2.8^2, 9.61=3.1^2, 

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);  
   
  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::y1*"-0.166"+
                   taylorSimplex::y2*"-0.174893"+
                   taylorSimplex::y3*"-0.174893"+
                   taylorSimplex::y4*"0.306136"+
                   taylorSimplex::y5*"0.143"+
                   taylorSimplex::y6*"0.143"+
                   taylorSimplex::dih*"-.263488"+
                   taylorSimplex::unit*"0.414686"; 
   
  const taylorFunction* K[1]={&F};  

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt)) 
    cout<<"star2 complete!"<<endl;
  else
    cout<<"star2 failed!!"<<endl;

}

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  // if y1>2.25841, 3.1>y4>2.8, than  
  //    
  // -vol(S) - .166 y7 +.008893 (y2+y3) -.036136 y4 + .143 (y8+y9)
  // +.80825 < 0
  //
  /////////////////////////////////////////////////////////////////

void star3()
{
  cellOption opt;
  domain x,z,x0,z0;

  interval tx[6]={a,"4","4","7.84","4","4"};
  interval tz[6]={v,v,v,"9.61",v,v};

  //9.61=3.1^2, 7.84=2.8^2

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
 
  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::y1*"-0.166"+
                   taylorSimplex::y2*"0.008893"+
                   taylorSimplex::y3*"0.008893"+ 
                   taylorSimplex::y4*"-0.306136"+
                   taylorSimplex::y5*"0.143"+
                   taylorSimplex::y6*"0.143"+
                   taylorSimplex::unit*"0.80825"; 
   
  const taylorFunction* K[1]={&F};  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt)) 
    cout<<"star3 complete!"<<endl;
  else
    cout<<"star3 failed!!"<<endl;

}

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  // if y1<2.25841, 2.95>y4>2.8, than  
  //    
  // -vol(S) - .166 y1 -.120535 (y2+y3)+.10788 y4 + .143 (y5+y6)
  // -.104704(dih-1.719)+.251657 < 0
  //
  // .104704*1.719+.251657=.431643
  //
  /////////////////////////////////////////////////////////////////

void star4()
{
  cellOption opt;
  domain x,z,x0,z0;

  interval tx[6]={"4","4","4","7.84","4","4"};
  interval tz[6]={a,v,v,"8.7025",v,v};

  // 9.61=3.1^2, 7.84=2.8^2, 8.7025=2.95^2

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
 
  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::y1*"-0.166"+
                   taylorSimplex::y2*"-0.120535"+
                   taylorSimplex::y3*"-0.120535"+ 
                   taylorSimplex::y4*"0.107880"+
                   taylorSimplex::y5*"0.143"+
                   taylorSimplex::y6*"0.143"+
                   taylorSimplex::dih*"-0.104704"+
                   taylorSimplex::unit*"0.431643"; 
   
  const taylorFunction* K[1]={&F};    

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt)) 
    cout<<"star4 complete!"<<endl;
  else
    cout<<"star4 failed!!"<<endl;

}

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  // if y1<2.25841, 2.95>y4>2.8, than  
  //    
  // -vol(S) - .166 y7 -.045465 (y2+y3)-.107880 y4 + .143 (y8+y9)
  // +.518343 < 0
  //
  /////////////////////////////////////////////////////////////////

void star5()
{
  cellOption opt;
  domain x,z,x0,z0;

  interval tx[6]={"4","4","4","7.84","4","4"};
  interval tz[6]={a,v,v,"8.7025",v,v};

  // 9.61=3.1^2, 7.84=2.8^2, 8.7025=2.95^2

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
    
  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::y1*"-0.166"+
                   taylorSimplex::y2*"-0.045465"+
                   taylorSimplex::y3*"-0.045465"+ 
                   taylorSimplex::y4*"-0.107880"+
                   taylorSimplex::y5*"0.143"+
                   taylorSimplex::y6*"0.143"+
                   taylorSimplex::unit*"0.518343"; 
   
  const taylorFunction* K[1]={&F};    

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt)) 
    cout<<"star5 complete!"<<endl;
  else
    cout<<"star5 failed!!"<<endl;

}

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  // if y1<2.25841, 3.1>y4>2.95, 2.25841>y7>2 than  
  //    
  // -vol(S) - .166 y1 -.114552 (y2+y3)+.115382 y4 + .143 (y5+y6)
  // -.15342(dih-1.719)+.193572 < 0
  //
  // .15342*1.719+.193572=.457301
  //
  /////////////////////////////////////////////////////////////////



void star6()
{
  cellOption opt;
  domain x,z,x0,z0;

  interval tx[6]={"4","4","4","8.7025","4","4"};
  interval tz[6]={a,v,v,"9.61",v,v};

  // 9.61=3.1^2, 8.7025=2.95^2

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
                   taylorSimplex::unit*"0.457301"; 
   
  const taylorFunction* K[1]={&F};    

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt)) 
    cout<<"star6 complete!"<<endl;
  else
    cout<<"star6 failed!!"<<endl;

}

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  // if y1<2.25841, 3.1>y4>2.95, 2.25841>y7>2 than  
  //    
  // -vol(S) - .166 y7 -.051448 (y2+y3)-.115382 y4 + .143 (y8+y9)
  // +.576428 < 0
  //
  /////////////////////////////////////////////////////////////////

void star7()
{
  cellOption opt;
  domain x,z,x0,z0;

  interval tx[6]={"4","4","4","8.7025","4","4"};
  interval tz[6]={a,v,v,"9.61",v,v};

  // 9.61=3.1^2, 7.84=2.8^2, 8.7025=2.95^2

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
                   taylorSimplex::unit*"0.576428"; 
   
  const taylorFunction* K[1]={&F};    

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt)) 
    cout<<"star7 complete!"<<endl;
  else
    cout<<"star7 failed!!"<<endl;

}

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  // if y1<2.25841, 3.1>y4>2.95, 2T>y7>2.25841, 2.25841>y2>2 than  
  //    
  // -vol(S) - .166 y1 -.279805 y2 - .340136 y3 +.422343 y4 + .143 (y5+y6)
  // -.380615 (dih-1.719)+.147006 < 0
  //
  // .380615*1.719+.147006=.801283
  //
  /////////////////////////////////////////////////////////////////

void star8()
{
  cellOption opt;
  domain x,z,x0,z0;

  interval tx[6]={"4","4","4","8.7025","4","4"};
  interval tz[6]={a,a,v,"9.61",v,v};

  //9.61=3.1^2, 7.84=2.8^2, 8.7025=2.95^2

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
 
  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::y1*"-0.166"+
                   taylorSimplex::y2*"-0.279805"+
                   taylorSimplex::y3*"-0.340136"+ 
                   taylorSimplex::y4*"0.422343"+
                   taylorSimplex::y5*"0.143"+
                   taylorSimplex::y6*"0.143"+
                   taylorSimplex::dih*"-0.380615"+
                   taylorSimplex::unit*"0.801283"; 
   
  const taylorFunction* K[1]={&F};    

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt)) 
    cout<<"star8 complete!"<<endl;
  else
    cout<<"star8 failed!!"<<endl;

}

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  // if y1<2.25841, 3.1>y4>2.95, 2T>y7>2.25841, 2.25841>y2>2 than 
  //    
  // -vol(S) - .166 y7 +.113805 y2 + .174136 y3 -.422343 y4 + .143 (y8+y9)
  // +.622994 < 0
  //
  /////////////////////////////////////////////////////////////////

void star9()
{
  cellOption opt;
  domain x,z,x0,z0;

  interval tx[6]={a,"4","4","8.7025","4","4"};

  interval tz[6]={v,a,v,"9.61",v,v};

  // 9.61=3.1^2, 8.7025=2.95^2

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
 
  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::y1*"-0.166"+
                   taylorSimplex::y2*"0.113805"+
                   taylorSimplex::y3*"0.174136"+ 
                   taylorSimplex::y4*"-0.422343"+
                   taylorSimplex::y5*"0.143"+
                   taylorSimplex::y6*"0.143"+
                   taylorSimplex::unit*"0.622994"; 
   
  const taylorFunction* K[1]={&F};    

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt)) 
    cout<<"star9 complete!"<<endl;
  else
    cout<<"star9 failed!!"<<endl;

}

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  //
  // if y1<2.25841, 3.1>y4>2.95, 2T>y7>2.25841, 2T>y2>2.25841 than  
  //    
  // -vol(S) - .166 y1 -.195794 y2 - .147301 y3 +.227016 y4 + .143 (y5+y6)
  // -.206212 (dih-1.719)+.107554 < 0
  //
  // .206212*1.719+.107554=0.462032428
  //
  /////////////////////////////////////////////////////////////////

void star10()
{
  cellOption opt;
  domain x,z,x0,z0;

  interval tx[6]={"4",a,"4","8.7025","4","4"};
  interval tz[6]={a,v,v,"9.61",v,v};

  // 9.61=3.1^2,  8.7025=2.95^2

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
                   taylorSimplex::unit*"0.462032428"; 
   
  const taylorFunction* K[1]={&F};    

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt)) 
    cout<<"star10 complete!"<<endl;
  else
    cout<<"star10 failed!!"<<endl;

}

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  // 
  // if 2T>y1>2.25841, 3.1>y4>2.95, 2T>y7>2.25841, 2T>y2>2.25841 than  
  // 
  // -vol(S) - .166 y7 +.029794 y2 -0.018699 y3 -.227016 y4 + .143 (y8+y9)
  // +.662446 < 0
  //
  /////////////////////////////////////////////////////////////////

void star11()
{
  cellOption opt;
  domain x,z,x0,z0;

  interval tx[6]={a,a,"4","8.7025","4","4"};
  interval tz[6]={v,v,v,"9.61",v,v};

  //9.61=3.1^2, 8.7025=2.95^2

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
    cout<<"star11 complete!"<<endl;
  else
    cout<<"star11 failed!!"<<endl;

}

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  // 
  // dih<1.86, y1,y2<2.273,2.244 resp than y4<3.6
  //
  /////////////////////////////////////////////////////////////////
  
  
void long1()
{
  cellOption opt;
  domain x,z,x0,z0;  

  interval tx[6]={"4","4","4",v,"4","4"};
  interval tz[6]={"5.166529","5.035536",v,"15",v,v};

  //2.273^2=5.166529, 2.244^2=5.035536, 15 > 3.6^2

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
      
  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-3.6";
   
  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.86";

  opt.setDihMax(1.86);

  const taylorFunction* K[2]={&F,&G};  

  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt)) 
    cout<<"long1 complete!"<<endl;
  else
    cout<<"long1 failed!!"<<endl;

}

void long2()
{
  cellOption opt;
  domain x,z,x0,z0;

  interval tx[6]={"4","4","4","12.96","4","4"};
  interval tz[6]={"5.035536","5.16653",v,"15",v,v};
  
    //2.273^2=5.166529, 2.244^2=5.035536, 15 > 3.6^2

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
  
  taylorFunction F=taylorSimplex::y4+
                   taylorSimplex::unit*"-3.6";
 
  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.86";

  opt.setDihMax(1.86);

  const taylorFunction* K[2]={&F,&G};  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"long2 complete!"<<endl;
  else
    cout<<"long2 failed!!"<<endl;
}


//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  // 
  // 3.1<y4<3.2, dih<1.719 -> y2+y3>4.2
  //
  /////////////////////////////////////////////////////////////////

void table1()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","9.61","4","4"};
  interval tz[6]={"5.16653","5.035536",v,"10.24",v,v};
  
  // 3.1^2=9.61, 3.2^2=10.24, 2.273^2=5.16653, 2.244^2=5.035536, 9.61=3.1^2

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::y2*"-1"+
                   taylorSimplex::y3*"-1"+
                   taylorSimplex::unit*"4.21";
  
  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.719";

  opt.setDihMax(1.719);
    
  const taylorFunction* K[2]={&F,&G};
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"table1 complete!"<<endl;
  else
    cout<<"table1 failed!!"<<endl;
}

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  // 
  // 3.2<y4<3.3, dih<1.719 -> y2+y3>4.37
  //
  /////////////////////////////////////////////////////////////////

void table2()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","10.24","4","4"};
  interval tz[6]={"5.16653","5.035536",v,"10.89",v,v};

  // 3.2^2=10.24, 3.3^2=10.89, 2.273^2=5.16653, 2.244^2=5.035536

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::y2*"-1"+
                   taylorSimplex::y3*"-1"+
                   taylorSimplex::unit*"4.37";
    
  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.719";
    
  opt.setDihMax(1.719);

  const taylorFunction* K[2]={&F,&G};  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,1,opt))
    cout<<"table2 complete!"<<endl;
  else
    cout<<"table2 failed!!"<<endl;
}

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  // 
  // 3.3<y4<3.4, dih<1.719 -> y2+y3>4.53
  //
  /////////////////////////////////////////////////////////////////


void table3()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","10.89","4","4"};
  interval tz[6]={"5.16653","5.035536",v,"11.56",v,v};

  // 3.4^2=11.56, 3.3^2=10.89, 2.273^2=5.16653, 2.244^2=5.035536

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);


  taylorFunction F=taylorSimplex::y2*"-1"+
                   taylorSimplex::y3*"-1"+
                   taylorSimplex::unit*"4.53";
  
  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.719";

  opt.setDihMax(1.719);
      
  const taylorFunction* K[2]={&F,&G};
  
  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"table3 complete!"<<endl;
  else
    cout<<"table3 failed!!"<<endl;
}

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  // 
  // 3.4<y4<3.6, dih<1.719 -> y2+y3>4.8
  //
  /////////////////////////////////////////////////////////////////

void table4()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","11.56","4","4"};
  interval tz[6]={"5.16653","5.035536",v,"12.96",v,v};

  // 3.4^2=11.56, 3.6^2=12.96, 2.273^2=5.16653, 2.244^2=5.035536

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::y2*"-1"+
                   taylorSimplex::y3*"-1"+
                   taylorSimplex::unit*"4.8";
    
  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.719";
    
  opt.setDihMax(1.719);

  const taylorFunction* K[2]={&F,&G};  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"table4 complete!"<<endl;
  else
    cout<<"table4 failed!!"<<endl;
}

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  // 
  // 3.1<y4<3.2, y2+y3>4.21, dih<1.719 -> mu(A)>.057
  //
  /////////////////////////////////////////////////////////////////


void table5()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","9.61","4","4"};
  interval tz[6]={"5.16653","5.035536",v,"10.24",v,v};

  // 3.2^2=10.24, 3.1^2=9.61, 2.273^2=5.16653, 2.244^2=5.035536

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.057";
        
  taylorFunction G=taylorSimplex::y2+
                   taylorSimplex::y3+
                   taylorSimplex::unit*"-4.21";

  taylorFunction H=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.719";

  opt.setDihMax(1.719);

  const taylorFunction* K[3]={&F,&G,&H};  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,3,opt))
    cout<<"table5 complete!"<<endl;
  else
    cout<<"table5 failed!!"<<endl;
} 

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  // 
  // 3.2<y4<3.3, y2+y3>4.37, dih<1.719 -> mu(A)>.069
  //
  /////////////////////////////////////////////////////////////////

void table6()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","10.24","4","4"};
  interval tz[6]={"5.16653","5.035536",v,"10.89",v,v};

  // 3.2^2=10.24, 3.3^2=10.89, 2.273^2=5.16653, 2.244^2=5.035536

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.069";
    
  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.719";
    
  taylorFunction H=taylorSimplex::y2+
                   taylorSimplex::y3+
                   taylorSimplex::unit*"-4.37";

  opt.setDihMax(1.719);

  const taylorFunction* K[3]={&F,&G,&H};
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,3,opt))
    cout<<"table6 complete!"<<endl;
  else
    cout<<"table6 failed!!"<<endl;
} 

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  // 
  // 3.3<y4<3.4, y2+y3>4.53, dih<1.719 -> mu(A)>.083
  //
  /////////////////////////////////////////////////////////////////

void table7()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","10.89","4","4"};
  interval tz[6]={"5.16653","5.035536",v,"11.56",v,v};

  // 3.3^2=10.89, 3.4^2=11.56, 2.273^2=5.16653, 2.244^2=5.035536

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.083";
  
  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.719";
    
  taylorFunction H=taylorSimplex::y2+
                   taylorSimplex::y3+
                   taylorSimplex::unit*"-4.53";
  
  opt.setDihMax(1.719);


  const taylorFunction* K[3]={&F,&G,&H};
   
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,3,opt))
    cout<<"table7 complete!"<<endl;
  else
    cout<<"table7 failed!!"<<endl;
} 

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  // 
  // 3.4<y4<3.6, y2+y3>4.8, dih<1.719 -> mu(A)>.11
  //
  /////////////////////////////////////////////////////////////////

void table8()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","11.56","4","4"};
  interval tz[6]={"5.16653","5.035536",v,"12.96",v,v};

  // 3.6^2=12.96, 3.4^2=11.56, 2.273^2=5.16653, 2.244^2=5.035536

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.11";
  
  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.719";
    
  taylorFunction H=taylorSimplex::y2+
                   taylorSimplex::y3+
                   taylorSimplex::unit*"-4.8";

  opt.setDihMax(1.719);

  const taylorFunction* K[3]={&F,&G,&H};
    
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,3,opt))
    cout<<"table8 complete!"<<endl;
  else
    cout<<"table8 failed!!"<<endl;
} 

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  // 
  // 3.1<y4<3.2, y2+y3>4.21, dih<2.419 -> mu(B)>.02
  //
  /////////////////////////////////////////////////////////////////


void table9()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","9.61","4","4"};
  interval tz[6]={v,"5.035536",v,"10.24",v,v};

  // 3.1^2=9.61, 3.2^2=10.24, 2.244^2=5.035536

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.0196";   // was .0196
    
  taylorFunction G=taylorSimplex::dih*"-1"+
                 taylorSimplex::unit*"2.419";
    
  taylorFunction H=taylorSimplex::y2+
                   taylorSimplex::y3+
                   taylorSimplex::unit*"-4.21";

  opt.setDihMax(2.419);

  const taylorFunction* K[3]={&F,&G,&H};  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,3,opt))
    cout<<"table9 complete!"<<endl;
  else
    cout<<"table9 failed!!"<<endl;
} 

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  // 
  // 3.2<y4<3.3, y2+y3>4.37, dih<2.419 -> mu(B)>.016
  //
  /////////////////////////////////////////////////////////////////



void table10()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","10.24","4","4"};
  interval tz[6]={v,"5.035536",v,"10.89",v,v};

  // 3.3^2=10.89, 3.2^2=10.24, 2.244^2=5.035536

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.016";

  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"2.419";
    
  taylorFunction H=taylorSimplex::y2+
                   taylorSimplex::y3+
                   taylorSimplex::unit*"-4.37";
  
  opt.setDihMax(2.419);

  const taylorFunction* K[3]={&F,&G,&H};  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,3,opt))
    cout<<"table10 complete!"<<endl;
  else
    cout<<"table10 failed!!"<<endl;
} 

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  // 
  // 3.3<y4<3.4, y2+y3>4.53, dih<2.419 -> mu(B)>.009
  //
  /////////////////////////////////////////////////////////////////


void table11()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","10.89","4","4"};
  interval tz[6]={v,"5.035536",v,"11.56",v,v};

  // 3.3^2=10.89, 3.4^2=11.56, 2.244^2=5.035536

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.009";
  
  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"2.419";
    
  taylorFunction H=taylorSimplex::y2+
                   taylorSimplex::y3+
                   taylorSimplex::unit*"-4.53";

  opt.setDihMax(2.419);

  const taylorFunction* K[3]={&F,&G,&H};  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,3,opt))
    cout<<"table11 complete!"<<endl;
  else
    cout<<"table11 failed!!"<<endl;
} 

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  // 
  // 3.4<y4<3.6, y2+y3>4.8, dih<2.419 -> mu(B)>-.02
  //
  /////////////////////////////////////////////////////////////////

void table12()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","11.56","4","4"};
  interval tz[6]={"5.166529","5.035536",v,"12.96",v,v};

  // 3.6^2=12.96, 3.4^2=11.56, 2.244^2=5.035536

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"-0.02";
    
  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.625";
    
  taylorFunction H=taylorSimplex::y2+
                   taylorSimplex::y3+
                   taylorSimplex::unit*"-4.8";

  const taylorFunction* K[3]={&F,&G,&H};
  
  opt.setDihMax(2.419);
   
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,3,opt))
    cout<<"table12 complete!"<<endl;
  else
    cout<<"table12 failed!!"<<endl;
} 



//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  // 
  // 3.2<y4<3.3, y2+y3>4.37, dih<1.532 -> mu(A)>.129
  //
  /////////////////////////////////////////////////////////////////

void spec1()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","10.24","4","4"};
  interval tz[6]={"5.166529","5.035536",v,"10.89",v,v};

  // 3.2^2=10.24, 3.3^2=10.89, 2.273^2=5.16653, 2.244^2=5.035536

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.12";
  
  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.532";
    
  taylorFunction H=taylorSimplex::y2+
                   taylorSimplex::y3+
                   taylorSimplex::unit*"-4.37";

  opt.setDihMax(1.532);

  const taylorFunction* K[3]={&F,&G,&H}; 
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,3,opt))
    cout<<"spec1 complete!"<<endl;
  else
    cout<<"spec1 failed!!"<<endl;
} 

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  // 
  // 3.1<y4<3.2, y2+y3>4.21, dih<1.614 -> mu(A)>.076
  //
  /////////////////////////////////////////////////////////////////

void spec2()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","9.61","4","4"};
  interval tz[6]={"5.166529","5.035536",v,"10.24",v,v};

  // 3.2^2=10.24, 3.1^2=9.61, 2.273^2=5.16653, 2.244^2=5.035536

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.076";

  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.614";
    
  taylorFunction H=taylorSimplex::y2+
                   taylorSimplex::y3+
                   taylorSimplex::unit*"-4.21";

  const taylorFunction* K[3]={&F,&G,&H};
   
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,3,opt))
    cout<<"spec2 complete!"<<endl;
  else
    cout<<"spec2 failed!!"<<endl;
} 


//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  // 
  // 3.1<y4<3.2, dih<1.76 -> y2+y3>4.134
  //
  /////////////////////////////////////////////////////////////////

void tableA1()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","9.61","4","4"};
  interval tz[6]={"5.16653","5.035536",v,"10.24",v,v};

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  // 3.2^2=10.24, 3.1^2=9.61, 2.273^2=5.16653, 2.244^2=5.035536

  taylorFunction F=taylorSimplex::y2*"-1"+
                   taylorSimplex::y3*"-1"+
                   taylorSimplex::unit*"4.134";
    
  taylorFunction G=taylorSimplex::dih*"-1"+
                 taylorSimplex::unit*"1.76";
    
  opt.setDihMax(1.76);

  const taylorFunction* K[2]={&F,&G}; 
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"tableA1 complete!"<<endl;
  else
    cout<<"tableA1 failed!!"<<endl;
}

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  // 
  // 3.2<y4<3.3, dih<1.76 -> y2+y3>4.29
  //
  /////////////////////////////////////////////////////////////////

void tableA2()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","10.24","4","4"};
  interval tz[6]={"5.16653","5.035536",v,"10.89",v,v};

  // 3.2^2=10.24, 3.3^2=10.89, 2.273^2=5.16653, 2.244^2=5.035536

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::y2*"-1"+
                   taylorSimplex::y3*"-1"+
                   taylorSimplex::unit*"4.29";
  
  taylorFunction G=taylorSimplex::dih*"-1"+
                 taylorSimplex::unit*"1.76";
    
  opt.setDihMax(1.76);

  const taylorFunction* K[2]={&F,&G};  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"tableA2 complete!"<<endl;
  else
    cout<<"tableA2 failed!!"<<endl;
}

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  // 
  // 3.3<y4<3.4, dih<1.76 -> y2+y3>4.44
  //
  /////////////////////////////////////////////////////////////////

void tableA3()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","10.89","4","4"};
  interval tz[6]={"5.16653","5.035536",v,"11.56",v,v};

  // 3.4^2=11.56, 3.3^2=10.89, 2.273^2=5.16653, 2.244^2=5.035536

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);
 
  taylorFunction F=taylorSimplex::y2*"-1"+
                   taylorSimplex::y3*"-1"+
                   taylorSimplex::unit*"4.44";
  
  taylorFunction G=taylorSimplex::dih*"-1"+
                 taylorSimplex::unit*"1.76";
 
  opt.setDihMax(1.76);
   
  const taylorFunction* K[2]={&F,&G};  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"tableA3 complete!"<<endl;
  else
    cout<<"tableA3 failed!!"<<endl;
}

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  // 
  // 3.4<y4<3.6, dih<1.76 -> y2+y3>4.64
  //
  /////////////////////////////////////////////////////////////////

void tableA4()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","11.56","4","4"};
  interval tz[6]={"5.16653","5.035536",v,"12.96",v,v};

  // 3.4^2=11.56, 3.6^2=12.96, 2.273^2=5.16653, 2.244^2=5.035536

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::y2*"-1"+
                   taylorSimplex::y3*"-1"+
                   taylorSimplex::unit*"4.64";

  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.76";
    
  opt.setDihMax(1.76);

  const taylorFunction* K[2]={&F,&G};  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,2,opt))
    cout<<"tableA4 complete!"<<endl;
  else
    cout<<"tableA4 failed!!"<<endl;
}

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  // 
  // 3.1<y4<3.2, dih<1.76, y2+y3>4.134 -> mu(C)>.049
  //
  /////////////////////////////////////////////////////////////////

void tableA5()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","9.61","4","4"};
  interval tz[6]={"5.16653","5.035536",v,"10.24",v,v};

  // 3.1^2=9.61, 3.2^2=10.24, 2.273^2=5.16653, 2.244^2=5.035536

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.049";

  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.76";
    
  taylorFunction H=taylorSimplex::y2+
                   taylorSimplex::y3+
                   taylorSimplex::unit*"-4.134";

  const taylorFunction* K[3]={&F,&G,&H}; 
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,3,opt))
    cout<<"tableA5 complete!"<<endl;
  else
    cout<<"tableA5 failed!!"<<endl;
} 

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  // 
  // 3.2<y4<3.3, dih<1.76, y2+y3>4.29 -> mu(C)>.0617
  //
  /////////////////////////////////////////////////////////////////

void tableA6()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","10.24","4","4"};
  interval tz[6]={"5.16653","5.035536",v,"10.89",v,v};

  // 3.2^2=10.24, 3.3^2=10.89, 2.273^2=5.16653, 2.244^2=5.035536

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.0617";
  
  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.76";
    
  taylorFunction H=taylorSimplex::y2+
                   taylorSimplex::y3+
                   taylorSimplex::unit*"-4.29";

  const taylorFunction* K[3]={&F,&G,&H};
  
  opt.setDihMax(1.76);  
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,3,opt))
    cout<<"tableA6 complete!"<<endl;
  else
    cout<<"tableA6 failed!!"<<endl;
} 

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  // 
  // 3.3<y4<3.4, dih<1.76, y2+y3>4.44 -> mu(C)>.074
  //
  /////////////////////////////////////////////////////////////////

void tableA7()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","10.89","4","4"};
  interval tz[6]={"5.16653","5.035536",v,"11.56",v,v};

  // 3.4^2=11.56, 3.3^2=10.89, 2.273^2=5.16653, 2.244^2=5.035536

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.074";
     
  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.76";
    
  taylorFunction H=taylorSimplex::y2+
                   taylorSimplex::y3+
                   taylorSimplex::unit*"-4.44";

  opt.setDihMax(1.76);

  const taylorFunction* K[3]={&F,&G,&H};
    
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,3,opt))
    cout<<"tableA7 complete!"<<endl;
  else
    cout<<"tableA7 failed!!"<<endl;
} 

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  // 
  // 3.4<y4<3.6, dih<1.76, y2+y3>4.64 -> mu(C)>.093
  //
  /////////////////////////////////////////////////////////////////

void tableA8()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","11.56","4","4"};
  interval tz[6]={"5.16653","5.035536",v,"12.96",v,v};

  // 3.4^2=11.56, 3.6^2=12.96, 2.273^2=5.16653, 2.244^2=5.035536

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.093";
  
  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.76";
    
  taylorFunction H=taylorSimplex::y2+
                   taylorSimplex::y3+
                   taylorSimplex::unit*"-4.64";

  opt.setDihMax(1.76);

  const taylorFunction* K[3]={&F,&G,&H};
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,3,opt))
    cout<<"tableA8 complete!"<<endl;
  else
    cout<<"tableA8 failed!!"<<endl;
} 

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  // 
  // 3.1<y4<3.2, dih<2.11, y2+y3>4.134 -> mu(D)>.024
  //
  /////////////////////////////////////////////////////////////////

void tableA9()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","9.61","4","4"};
  interval tz[6]={v,"5.035536",v,"10.24",v,v};

  // 3.1^2=9.61, 3.2^2=10.24, 2.244^2=5.035536

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.024";
  
  taylorFunction G=taylorSimplex::dih*"-1"+
                 taylorSimplex::unit*"2.11";
    
  taylorFunction H=taylorSimplex::y2+
                   taylorSimplex::y3+
                   taylorSimplex::unit*"-4.134";

  opt.setDihMax(2.11); 

  const taylorFunction* K[3]={&F,&G,&H}; 
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,3,opt))
    cout<<"tableA9 complete!"<<endl;
  else
    cout<<"tableA9 failed!!"<<endl;
} 

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  // 
  // 3.2<y4<3.3, dih<2.11, y2+y3>4.29 -> mu(D)>.020
  //
  /////////////////////////////////////////////////////////////////

void tableA10()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","10.24","4","4"};
  interval tz[6]={v,"5.035536",v,"10.89",v,v};

  // 3.2^2=10.24, 3.3^2=10.89, 2.244^2=5.035536

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.02";

  taylorFunction G=taylorSimplex::dih*"-1"+
                 taylorSimplex::unit*"2.11";
    
  taylorFunction H=taylorSimplex::y2+
                   taylorSimplex::y3+
                   taylorSimplex::unit*"-4.29";

  opt.setDihMax(2.11);

  const taylorFunction* K[3]={&F,&G,&H}; 
  
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,3,opt))
    cout<<"tableA10 complete!"<<endl;
  else
    cout<<"tableA10 failed!!"<<endl;
} 

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  // 
  // 3.3<y4<3.4, dih<2.11, y2+y3>4.44 -> mu(D)>.020
  //
  /////////////////////////////////////////////////////////////////

void tableA11()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","10.89","4","4"};
  interval tz[6]={v,"5.035536",v,"11.56",v,v};

  // 3.4^2=11.56, 3.3^2=10.89, 2.244^2=5.035536

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.02";
  
  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"2.11";
    
  taylorFunction H=taylorSimplex::y2+
                   taylorSimplex::y3+
                   taylorSimplex::unit*"-4.44";

  opt.setDihMax(2.11);

  const taylorFunction* K[3]={&F,&G,&H};
    
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,3,opt))
    cout<<"tableA11 complete!"<<endl;
  else
    cout<<"tableA11 failed!!"<<endl;
} 

//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  // 
  // 3.4<y4<3.6, dih<2.11, y2+y3>4.64 -> mu(D)>-.01
  //
  /////////////////////////////////////////////////////////////////

void tableA12()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","11.56","4","4"};
  interval tz[6]={v,"5.035536",v,"12.96",v,v};

  // 3.4^2=11.56, 3.6^2=12.96, 2.244^2=5.035536

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"-0.01";
  
  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"2.11";
    
  taylorFunction H=taylorSimplex::y2+
                   taylorSimplex::y3+
                   taylorSimplex::unit*"-4.64";

  opt.setDihMax(2.11);

  const taylorFunction* K[3]={&F,&G,&H};
   
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,3,opt))
    cout<<"tableA12 complete!"<<endl;
  else
    cout<<"tableA12 failed!!"<<endl;
} 


//*****************************************************************

  /////////////////////////////////////////////////////////////////
  // Want to show: 
  // 
  // 3.1<y4<3.2, dih<1.652, y2+y3>4.134 -> mu(C)>.074
  //
  /////////////////////////////////////////////////////////////////

void specA1()
{
  cellOption opt;
  domain x,z,x0,z0;
  
  interval tx[6]={"4","4","4","9.61","4","4"};
  interval tz[6]={"5.166529","5.035536",v,"10.24",v,v};

  // 3.1^2=9.61, 3.2^2=10.24, 2.244^2=5.035536,2.273^2=5.166529

  x = domain::lowerD(tx);
  z = domain::upperD(tz);
  x0 = domain::lowerD(tx);
  z0 = domain::upperD(tz);

  taylorFunction F=taylorSimplex::truncatedVoronoiVolume*"-1"+
                   taylorSimplex::sol*"0.42755"+
                   taylorSimplex::unit*"0.074";
  
  taylorFunction G=taylorSimplex::dih*"-1"+
                   taylorSimplex::unit*"1.652";
    
  taylorFunction H=taylorSimplex::y2+
                   taylorSimplex::y3+
                   taylorSimplex::unit*"-4.134";

  opt.setDihMax(1.652);

  const taylorFunction* K[3]={&F,&G,&H};
   
  if(prove::recursiveVerifier(0,x,z,x0,z0,K,3,opt))
    cout<<"specA1 complete!"<<endl;
  else
    cout<<"specA1 failed!!"<<endl;
} 

