#include<string.h>
#include<iostream.h>
#include<iomanip.h>
#include<fstream.h>
#include<stdlib.h>
#include<math.h>
#include<assert.h>
#include<time.h>

//***************************************************************************
//
//                                CONSTANTS
//
//***************************************************************************

const int NUM_FACES=40;        // Max faces in graph
const int MAX_EDGES=9;         // Max edges in a face
const int NUM_VERT=30;         // Max spheres within 2T of 0
const int MAX_DIH=12;          // Max faces around a vertex + other info
const int MAX_CON_EDGE=100;    // Max number of edges of the form y(1,2)
const double TET_LB=2;         // Min edge length
const double TET_UB=2.51682;   // Max edge length
const double M=0.42755;        // Slope


//***************************************************************************
//
//                               FUNCTIONS
//
//***************************************************************************


void GetData(int faceArr[NUM_FACES][MAX_EDGES],
             int vertArr[NUM_VERT][MAX_DIH],
             int& configFaceNum,
             int& configVertNum,
             ofstream& outFile,
             ifstream& inFile,
             int counter);




void PrintGraph(const int faceArr[NUM_FACES][MAX_EDGES],
                      int configFaceNum,
                      ofstream& outFile);


void GetEdgeList(int list[MAX_CON_EDGE][2],
                 int& edgeNum,
                 const int faceArr[NUM_FACES][MAX_EDGES],
                 int configFaceNum,
                 int configVertNum,
                 ofstream& outFile);

void ArrangeVert(const int faceArr[NUM_FACES][MAX_EDGES],
		 int vertArr[NUM_VERT][MAX_DIH],
		 const int configVertNum,
		 ofstream& outFile);

static char SignOf(double x)
  {
  if(x>=0) return '+';
  else return '-';
  }

void PrintProblem(int configFaceNum, ofstream& outFile, int num);

void PrintSpecial(const int faceArr[NUM_FACES][MAX_EDGES],
               int configFaceNum,
               ofstream& outFile);


void PrintEdgesEqual(const int list[MAX_CON_EDGE][2],
                     int edgeNum,
                     ofstream& outFile);

void PrintDihSum(const int faceArr[NUM_FACES][MAX_EDGES],
                 const int vertArr[NUM_VERT][MAX_DIH],
                 int configVertNum,
                 ofstream& outFile);

void PrintTets(const int faceArr[NUM_FACES][MAX_EDGES],
               int configFaceNum,
               ofstream& outFile);

void PrintQuads(const int faceArr[NUM_FACES][MAX_EDGES],
               int configFaceNum,
               ofstream& outFile);

void PrintSolid(const int faceArr[NUM_FACES][MAX_EDGES],
                int configFaceNum,
                ofstream& outFile);

void Print5Sets(const int faceArr[NUM_FACES][MAX_EDGES],
                const int vertArr[NUM_VERT][MAX_DIH],
                const int edgeList[MAX_CON_EDGE][2],
                int configFaceNum,
                int configVertNum,
                int configEdgeNum,
                ofstream& outFile);

void PrintDihEdge(const int faceArr[NUM_FACES][MAX_EDGES],
                  int configFaceNum,
                  ofstream& outFile);

void PrintVolEdge(const int faceArr[NUM_FACES][MAX_EDGES],
                  int configFaceNum,
                  ofstream& outFile);


void PrintQuadEdge(const int faceArr[NUM_FACES][MAX_EDGES],
                  int configFaceNum,
                  ofstream& outFile);

void PrintXST(const int faceArr[NUM_FACES][MAX_EDGES],
                  int configFaceNum,
                  ofstream& outFile);


void PrintSquanderFace(const int faceArr[NUM_FACES][MAX_EDGES],
                  int configFaceNum,
                  ofstream& outFile);

void PrintXST(const int vertArr[NUM_VERT][MAX_DIH],
                  int configVertNum,
                  ofstream& outFile);



void PrintBreakQuads(const int faceArr[NUM_FACES][MAX_EDGES],
                     int configFaceNum,
                     ofstream& outFile);

void PrintBoundSigma(const int faceArr[NUM_FACES][MAX_EDGES],
                     int configFaceNum,
                     ofstream& outFile);

void PrintBoundSolid(const int faceArr[NUM_FACES][MAX_EDGES],
                     int configFaceNum,
                     ofstream& outFile);

void PrintBoundDih(const int faceArr[NUM_FACES][MAX_EDGES],
                   int configFaceNum,
                   ofstream& outFile);

void PrintBoundEL(const int list[MAX_CON_EDGE][2],
                  int edgeNum,
                  int configVertNum,
                  ofstream& outFile);

void PrintBrokenStuff(const int faceArr[NUM_FACES][MAX_EDGES],
                      int configFaceNum,
                      ofstream& outFile);



//***************************************************************************
//
//                                 MAIN
//
//***************************************************************************


int main()
{
  ofstream cplexFile;                   
  int faces;                            
  int vertices;
  int edges;
  
  /* The arrays faceArr, vertArr, and edgeList are the main information 
     holders in this program.  The array faceArr has the same number of 
     rows as there are faces of the graph.  The first entry in a row is 
     the face's number.  So faceArr[3][0] will be 4, as c++ begins indexing
     at 0.  The second entry is the number of edges that particular face had.
     In particular, if the face numbered 4 is a pentagon, then faceArr[3][1]
     will be 5.  The remaining entries in a row are a cyclic ordering of the
     vertices of that face.  Thus, faceArr[3] will have length 7, the first 
     two I just described and the last 5 are the vertices of the pentagon. */


  int faceArr[NUM_FACES][MAX_EDGES]; 

  /* The array vertArr has the same number of rows as there are vertices
     of the graph.  As with faceArr, the first entry is the numbering of 
     the vertex.  So vertArr[n][0] = n+1.  The second entry is the number
     of faces around a particular vertex.  The next entries are the numbers
     of the faces around that vertex.  */


  int vertArr[NUM_VERT][MAX_DIH];
  

  /*  The array edgeList has the same number of rows as there are edges in
      the entire graph.  The first entry, as before, is the number of the 
      edge in my ordering.  So edgeList[n][0]=n+1.  The next two are the 
      numbers of the vertices which are the endpoints of the edge. */

  int edgeList[MAX_CON_EDGE][2];


  cplexFile.setf(ios::fixed, ios::floatfield);
  cplexFile.setf(ios::showpoint);
  cplexFile<<setprecision(6);   // Sets up output style

  ifstream inFile;  //the Java file
  
  char file[30]="simpleEx.dat";
  //cout<<"Enter file name"<<endl;
  //cin>>file;


  inFile.open(file);
  if(!inFile)
    {cout<<"Can't open input file!"<<endl;  //Make sure the file opens
    return 0;}

  int i,j;
  char title[20]="ceh.lp";
  char title2[20];
  char buffer[20];
  //cout<<"Enter a general name for the output files. (Ex. cplex4.lp)"<<endl;
  //cin>>title;
  int number=1297;
  //cout<<"Enter number of graphs"<<endl;
  //cin>>number;

  


  for(i=0;i<number;i++)
    {
      strcpy(title2,title);
      sprintf(buffer,"%d",i);
      strcat(title2,buffer);
      //cout<<endl<<"File is : "<<title2<<endl;
      cplexFile.open(title2);
      if(!cplexFile)
        {
          cout<<"**Can't open output file**"<<endl;
          return 1;
        }
      //That routine just opens the input and output files.

      
      GetData(faceArr, vertArr, faces, vertices, cplexFile, inFile,i);
      GetEdgeList(edgeList, edges, faceArr, faces, vertices, cplexFile);
      ArrangeVert(faceArr, vertArr, vertices, cplexFile);
      PrintProblem(faces, cplexFile,i);
      cplexFile<<endl<<endl<<"ST"<<endl<<endl;
      


      PrintEdgesEqual(edgeList, edges, cplexFile);
      PrintDihSum(faceArr, vertArr, vertices, cplexFile);
      PrintTets(faceArr, faces, cplexFile);
      PrintQuads(faceArr, faces, cplexFile);
      PrintSolid(faceArr, faces, cplexFile);
      PrintDihEdge(faceArr, faces, cplexFile);
      PrintVolEdge(faceArr, faces, cplexFile);
      //PrintQuadEdge(faceArr, faces, cplexFile);
      //PrintXST(vertArr,vertices,cplexFile);

      PrintSquanderFace(faceArr, faces, cplexFile);

 Print5Sets(faceArr, vertArr, edgeList, faces, vertices, edges, cplexFile);
      
      


      cplexFile<<endl<<endl<<"BOUNDS"<<endl<<endl;
      cplexFile<<"c=1"<<endl<<endl;
      cplexFile<<endl<<"pi = 3.141592653589793"<<endl<<endl;

      PrintBoundSigma(faceArr, faces, cplexFile);
      PrintBoundSolid(faceArr, faces, cplexFile);
      PrintBoundDih(faceArr, faces, cplexFile);
      PrintBoundEL(edgeList, edges, vertices, cplexFile);
      
      cplexFile<<endl<<endl<<"END"<<endl;
      PrintGraph(faceArr, faces, cplexFile);
      

      cplexFile.close();
      //cout<<endl<<endl<<i+1<< " complete"<<endl;
    }

  cout<<endl<<"All graphs complete."<<endl<<endl;

  return 0;
}


//***************************************************************************
//
//                         FUNCTION BODIES
//
//***************************************************************************
//
//                          SET UP ARRAYS
//
//***************************************************************************

void GetData(/* inout */ int faceArr[NUM_FACES][MAX_EDGES],
               /* inout */ int vertArr[NUM_VERT][MAX_DIH],
               /* inout */ int& configFaceNum,
               /* inout */ int& configVertNum,
               /* inout */ ofstream& outFile,
			   ifstream& inFile,
			   int counter)
{
  

  int i,j,k;
  char ch='a';
  char ch1='a';
  char str[3]="aa";

  for(j=0;j<NUM_FACES;j++)    //Init. faceArr to -1
    for(k=0;k<MAX_EDGES;k++)
      faceArr[j][k]=-1;

  for(j=0;j<NUM_VERT;j++)    //Init. vertArr to -1
    for(k=0;k<MAX_DIH;k++)
      vertArr[j][k]=-1;

  for(j=0;j < (NUM_VERT) ;j++)    //Init. first in vertArr to vert number
    vertArr[j][0]=j+1;

  /* cout<<"PrePrint for faceArr"<<endl<<endl; //print to check initialization

  for(i=0;i<NUM_FACES;i++)
    { for(j=0;j<MAX_EDGES;j++)
      {cout<<faceArr[i][j]<<" ";}
    cout<<endl;}

  cout<<endl<<"PrePrint for vertArr"<<endl<<endl; //print to check 


    for(i=0;i<NUM_VERT;i++)
      { for(j=0;j<MAX_DIH;j++)
	{cout<<vertArr[i][j]<<" ";}
      cout<<endl;}*/



  /*  TESTING
    
    inFile>>ch;
    while(ch != '"') inFile>>ch;   // skip to first parenthesis
    inFile>>ch;  // skip the 0
    assert(ch=='0');
    while(ch != '"'){ cout<<ch<<" ";inFile>>ch;}  TEST STUFF
    return; */

  int faceNum,vertNum;
  configVertNum=0;
  
  inFile>>ch;
  while(ch != '"') inFile>>ch;  //get to first 0
  inFile>>ch;                   // skip the 0
  assert(ch == '0');            // Check to be sure
  if(ch != '0') {cout<<"Error in "<<counter<<" th loop."<<endl; return;}
  inFile>>str[0];
  inFile.get(ch1);
  if(ch1 != ' ') str[1]=ch1;
  configFaceNum=atoi(str);      // Array str is used to store the number. 
  if(configFaceNum>NUM_FACES)
    {
      cout<<"** Max number of faces (40) exceeded in graph "
	  <<counter<<" !  **"<<endl; 
      return; 
    }

  str[0]='a';str[1]='a';
  
  for(i=0;i<configFaceNum;i++)
    {
      inFile>>str[0];
      faceArr[i][0]=i+1;        // Set the numbering
      faceArr[i][1]=atoi(str);
      
      assert(faceArr[i][1] < MAX_EDGES);      

      str[0]='a';
      for(j=2;j<faceArr[i][1]+2;j++)
	{
	  inFile>>str[0];
	  inFile.get(ch1);
	  if(ch1 != ' ') str[1]=ch1;
	  faceArr[i][j]=atoi(str);  // Set the value from Java
	  faceArr[i][j]++;          // Re-index
	  
	  assert(faceArr[i][j] < NUM_VERT);

	  str[0]='a';str[1]='a';    // Reset string
	  
	  // Next line scrolls through vertices, checking for bigger ones.
	  // Note:  It's already been reindexed from faceArr.

	  if(configVertNum<faceArr[i][j]) configVertNum=faceArr[i][j]; 

	  // This sticks the face next to its vertex in the vertArr.

	  
	  for(k=2;k<MAX_DIH+2;k++)
	    {
	      if( vertArr[(faceArr[i][j]-1)][k]==-1)
		{
		  vertArr[(faceArr[i][j]-1)][k]=faceArr[i][0];
		  break;
		}
	    }
	}
    }

  if(configVertNum>NUM_VERT) 
    {
      cout<<" ** Max vertices exceeded in "<<counter
	  <<" !! **"<<endl;
      return;
    }

  int numAround;

  // numAround counts the faces around the vertex.  This must be done 
  // after entering the faces for obvious reasons.

  for(i=0;i<configVertNum;i++)
    {
      numAround=0;
      for(j=2;j<MAX_DIH;j++)
        {
          while(vertArr[i][j] != -1 )
            {
              numAround++;
              j++;
            }
          vertArr[i][1]=numAround;
        }
    }

  
  //cout<<endl<<"Vertices="<<configVertNum<<endl;
  /*cout<<endl<<"Final faceArr"<<endl<<endl;

     for(i=0;i<NUM_FACES;i++)
       { for(j=0;j<MAX_EDGES;j++)
	 {cout<<faceArr[i][j]<<" ";}
       cout<<endl;}

     cout<<endl<<"final vertArr"<<endl<<endl;

    for(i=0;i<NUM_VERT;i++)
      { for(j=0;j<MAX_DIH;j++)
	{cout<<vertArr[i][j]<<" ";}
      cout<<endl;}*/
      
    inFile>>ch;
    while(ch != ',') inFile>>ch;  //Skip to the comma to be ready for next file

  return;


}

//***************************************************************************
//
//                        SET EDGE ARRAY
//
//***************************************************************************

void GetEdgeList(int list[MAX_CON_EDGE][2],
		 int& edgeNum,
		 const int faceArr[NUM_FACES][MAX_EDGES], 
		 int configFaceNum,
		 int configVertNum,
		 ofstream& outFile)



{
  int a,b,c,d;
  int index=0;
  int i;
  int j;
  int sum;
  int counter;
  int counter2;

  for(counter=0;counter<configFaceNum;counter++)
    {
      for(counter2=2;counter2<faceArr[counter][1]+1;counter2++)
	{
	  a=faceArr[counter][counter2];
	  b=faceArr[counter][counter2+1];
	  if(a<b)
	    {
	      c=a;
	      d=b;
	    }
	  else if(b<a)
	    {
	      c=b;
	      d=a;
	    }
	  else if(b==a)
	    cout<<"** Error in edge list, same numbers **"<<endl;

	  sum=0;
	  for(i=0;i<index;i++)
	    {
	      if(  c == list[i][0] && d == list[i][1] )
	        {sum++;}
	    }
	  if(sum==0)
	    {
	      list[index][0]=c;list[index][1]=d;
	      index++;
	    }
	}
       a=faceArr[counter][2];
       b=faceArr[counter][(faceArr[counter][1]+1)];
       if(a<b)
	 {
	   c=a;
	   d=b;
	 }
          else if(b<a)
            {
              c=b;
              d=a;
            }
          else if(b==a)
            cout<<"** Error in edge list, same numbers **"<<endl;

       sum=0;  
       for(i=0;i<index;i++)
	    {
	      if(  c == list[i][0] && d == list[i][1] ) 
	        sum++;
	    }


       if(sum==0)
	    {
	      list[index][0]=c;list[index][1]=d;	
	      index++;
	    }
    }

  /*cout<<endl<<"Edge List : "<<endl;
  for(i=0;i<index;i++)
    {
      for(j=0;j<2;j++)
	cout<<list[i][j]<<" ";
      if(list[i][j]==-1)
	{ 
	  cout<<endl<<endl<<"*** error!! edge is -1 in edgeArr"<<endl<<endl; 
	  return;
	}
      cout<<" ";
    } */
    

  edgeNum=index;
  //cout<<"Number of Edges = "<<edgeNum<<endl;

  return;
}

//***************************************************************************
//
//                  Sticks the type of vertex at the end of vertArr  
//
//***************************************************************************

void ArrangeVert(const int faceArr[NUM_FACES][MAX_EDGES],
                 int vertArr[NUM_VERT][MAX_DIH],
                 const int configVertNum,
                 ofstream& outFile)
{
  int i,j;
  int tet,quad,other;
  int sides;
  int faces;
  int F;
  
  for(i=0;i<configVertNum;i++)
    {
      tet=0;
      quad=0;
      other=0;
      faces=vertArr[i][1];
      for(j=2;j<faces+2;j++)
	{
	  F=vertArr[i][j];
	  sides=faceArr[F-1][1];
	  if(sides==3) tet++;
	  else if(sides==4) quad++;
	  else if(sides>4) other++;
	}
      vertArr[i][9]=tet;
      vertArr[i][10]=quad;
      vertArr[i][11]=other;
    }
 
  /*cout<<endl<<"VertArr with added 3 : "<<endl;
  for(i=0;i<configVertNum;i++)
   {
     for(j=0;j<MAX_DIH;j++)
       cout<<vertArr[i][j]<<" ";
     cout<<endl;
   }*/
   
  return;
}

//***************************************************************************
//
//                       PRINT OPTIMIZATION PROBLEM
//
//***************************************************************************

void PrintProblem(/* in */ int configFaceNum,
		 /* inout */ ofstream& outFile,
                 /* in */ int num)
{
  outFile<<"\\ LP Format CPLEX file generated by C++ on "<<time(0)<<endl
         <<endl<<"\\ Problem: Bound the volume of configuration"<<endl
         <<"\\ "<<num+1<< " arising in the Dodecahedral conjecture."         
         <<endl<<endl<<"MAXIMIZE"<<endl<<endl<<"score : ";     



  int counter;
 
  for(counter=0;counter<configFaceNum;counter++)
    outFile<<"+sigma{"<<counter+1<<"}";

  

  return;

}

//***************************************************************************
//
//         The deep observation that edge (2,3)=(3,2)
//
//***************************************************************************

void PrintEdgesEqual(/* in */ const int list[MAX_CON_EDGE][2],
		     /* in */ int edgeNum,
		     /* inout */ ofstream& outFile)
{
  int i;
  
  
  outFile<<endl<<"\\ edges equal "<<endl<<endl;

  for(i=0;i<edgeNum;i++)
    {
      outFile<<"y("<<list[i][0]<<","<<list[i][1]<<") - y("<<list[i][1]<<","
             <<list[i][0]<<") = 0"<<endl;
    }
  outFile<<endl<<endl;
  return;
}

//***************************************************************************
//
//         Sum of dihedral angles about a vertex = 2pi
//
//***************************************************************************

void PrintDihSum(/* in */ const int faceArr[NUM_FACES][MAX_EDGES],
/* in */ const int vertArr[NUM_VERT][MAX_DIH],
/* in */ int configVertNum,
/* inout */ ofstream& outFile)

{
  int counter;
  int counter2;

   outFile<<"\\ The dihedral angles around each vertex sum to 2pi."
          <<endl<<endl;

   for(counter=0;counter<configVertNum;counter++)
     {
       outFile<<"dihsum"<<vertArr[counter][0]<<" : -2 pi ";
       counter2=2;
       while(vertArr[counter][counter2] != -1 )
         {
           outFile<<"+dih{" << vertArr[counter][counter2] <<
	     "}"<<vertArr[counter][0];
           counter2++;
         }
       outFile<<" = 0"<<endl;
     }
   outFile<<endl;

 return;

}

//***************************************************************************
//
//                        Print Tet ineqs
//
//***************************************************************************

void PrintTets(/* in */ const int faceArr[NUM_FACES][MAX_EDGES], 
               /* in */ int configFaceNum,
               /* inout */ ofstream& outFile)

{
  const int NUM_TET_EQNS=54;   // Must be changed if you add any tet eqns
  
  double tets[NUM_TET_EQNS][3]={
                            {0.68,-1.88718,1.545510},
                            {0.68,-0.90746,0.706725},
    		            {0.68,-0.46654,0.329233},
                            {0.55889,0,-0.073648},
                            {0.63214,0,-0.13034},
         	            {0.73256,0,-0.23591},
                            {0.89346,0,-0.40505},
                            {0.3,0.5734,-0.978221},
                            {0.3,0.03668,0.024767},
	         /* 10 */   {0.3,-0.04165,0.121199},
			    {0.3,-0.1234,0.209279},
			    {0.42755,0.11509,-0.171859},
			    {0.42755,0.04078,-0.050713},
			    {0.42755,-0.11031,0.135633},
			    {0.42755,-0.13091,0.157363},
			    {0.55792,0.21394,-0.417998},
		            {0.55792,0.0068,-0.081902},
		            {0.55792,-0.0184,-0.051224},
		            {0.55792,-0.24335,0.193993},
		 /* 20 */   {0.68,0.30651,-0.648496},
		            {0.68,0.06965,-0.27800},
		            {0.68,-0.0172,-0.15662},
		            {0.68,-0.41812,0.287778},
		            {0.64934,0,-0.14843},
		            { 0.6196,0,-0.11800},
		            {0.58402,0,-0.090290},
		            {0.25181,0,0.096509},
		            {0.00909,0,0.199559},
		            {-0.93877,0,0.537892},
		 /* 30 */   {-0.93877,0.20211,0.27313},
		            {-0.93877,-0.63517,1.20578},
		            {-1.93877,0,0.854804},
		            {-1.93877,0.20211,0.621886},
		            {-1.93877,-0.63517,1.57648},
		            {0.42775,0,-0.000111},
		            {0.55792,0,-0.073037},
		            {0,0.07853,0.08865},
		            {0,0.00339,0.198693},
		            {0,-0.18199,0.396670},
	         /* 40 */   {0.42755,0.20000,-0.332061},
		            {0.3,0.36373,-0.582630},
		            {0.3,-0.20583,0.279851},
		            {0.3,-0.40035,0.446389},
		            {0.3,-0.83259,0.816450},
		            {0.42755,0.51838,-0.932759},
		            {0.42755,-0.29344,0.296513},
		            {0.42755,-0.57056,0.533768},
		            {0.42755,-1.18656,1.06115},
		            {0.42755,0,0},
			    {0.55792,0.67644,-1.29062},
		 /* 50 */   {0.55792,-0.38278,0.313365},
	        	    {0.55792,-0.74454,0.623085},
			    {0.55792,-1.54837,1.31128},
			    {0.68,0.82445,-1.62571}
                            };
                             


  int counter;
  int counter2;
  int counter3;

  for(counter=0; counter<NUM_TET_EQNS; counter++)
    {
      outFile<<"\\Tetrahedral Equation "<<counter+1<<endl<<endl;
      
      for(counter2=0; counter2<configFaceNum; counter2++)
	{
	  if(faceArr[counter2][1]!=3) continue;
	  for(counter3=2;counter3<5;counter3++)
	    {
	      outFile<<"TetE"<<counter+1<<"{"<<counter2+1<<"}"
                     <<faceArr[counter2][counter3]<<" : sigma{"
		     <<counter2+1<<"}"<<SignOf(tets[counter][0])
		     <<fabs(tets[counter][0])<<" solid{"<<faceArr[counter2][0]
		     <<"}"
		     <<SignOf(tets[counter][1])<<fabs(tets[counter][1])
		     <<" dih{"<<faceArr[counter2][0]<<"}"
		     <<faceArr[counter2][counter3]
		     <<SignOf(tets[counter][2])
		     <<fabs(tets[counter][2])<<" c <=0"<<endl;	      
	   }
      }
      outFile<<endl;
   }  

  return;

}

//***************************************************************************
//
//         Prints quad ineqs
//
//***************************************************************************

void PrintQuads(/* in */ const int faceArr[NUM_FACES][MAX_EDGES], 
                /* in */ int configFaceNum,
                /* inout */ ofstream& outFile)

{
  const int NUM_QUAD_EQNS=28;   // Must be changed if you add any quad eqns
  
  double quads[NUM_QUAD_EQNS][3]={
                                {0.42775,0,0.031350},     
                                {0.42775,0.15098,-0.3670},
				{0.42775,0.09098,-0.1737},
				{0.42775,0.00000,0.0310},
				{0.42775,-0.18519,0.3183},
				{0.42775,-0.20622,0.3438},
				{0.55792,0.30124,-1.0173},
				{0.55792,0.02921,-0.2101},
				{0.55792,0.00000,-0.1393},
				{0.55792,-0.05947,-0.0470},
		      /* 10 */  {0.55792,-0.39938,0.4305},
				{0.55792,-2.50210,2.8976},
				{0.68000,0.44194,-1.6264},
				{0.68000,0.10957,-0.6753},
				{0.68000,0.00000,-0.4029},
				{0.68000,-0.86096,0.8262},
				{0.68000,-2.44439,2.7002},
				{0.30000,0.12596,-0.1279},
				{0.30000,0.02576,0.1320},
				{0.30000,-0.00000,0.1945},
		      /* 20 */  {0.30000,-0.03700,0.2480},
				{0.30000,-0.22476,0.5111},
				{0.30000,-2.31852,2.9625},
				{0.00000,0.23227,-0.1042},
				{0.00000,-0.07448,0.5591},
				{0.00000,-0.22019,0.7627},
				{0.00000,-0.80927,1.5048},
				{0.00000,-5.84380,7.3468}
                               };


  int counter;
  int counter2;
  int counter3;

  for(counter=0; counter<NUM_QUAD_EQNS; counter++)
    {
      outFile<<"\\Quadrilateral Equation "<<counter+1<<endl<<endl;
      
      for(counter2=0; counter2<configFaceNum; counter2++)
	{
	  if(faceArr[counter2][1]!=4) continue;
	  for(counter3=2;counter3<6;counter3++)
	    {
	      outFile<<"QuadE"<<counter+1<<"{"<<faceArr[counter2][0]<<"}"
                     <<faceArr[counter2][counter3]<<" : sigma{"
		     <<faceArr[counter2][0]<<"}"<<SignOf(quads[counter][0])
		     <<fabs(quads[counter][0])<<" solid{"<<faceArr[counter2][0]
		     <<"}"<<SignOf(quads[counter][1])<<fabs(quads[counter][1])
		     <<" dih{"<<faceArr[counter2][0]<<"}"
		     <<faceArr[counter2][counter3]<<SignOf(quads[counter][2])
		     <<fabs(quads[counter][2])<<" c <= 0"<<endl;	      
	    }
	}
      outFile<<endl;
    }  

  return;

}

//***************************************************************************
//
// Solid angle is sum of dihedrals - (num edges - 2 pi)
//
//***************************************************************************

void PrintSolid(/* in */ const int faceArr[NUM_FACES][MAX_EDGES], 
                /* in */ int configFaceNum,
                /* inout */ ofstream& outFile)

{
  outFile<<"\\ Solid angles in terms of dihedral angles"<<endl<<endl;
  
  int counter;
  int counter2;

  for(counter=0;counter<configFaceNum;counter++)
    {
      outFile<<"DS"<<counter+1<<" : - solid{"<<counter+1<<"} ";
      for(counter2=2;counter2<faceArr[counter][1]+2;counter2++)
    	{outFile<<" + dih{" << counter+1 << "}" << faceArr[counter][counter2];}
      outFile<<" - "<<faceArr[counter][1]-2<<" pi = 0"<<endl;
    }

  outFile<<endl<<endl;
  return;

}

//***************************************************************************
//
//            $$$$ 5 Sets $$$$  Extra squander around a (5,0)
//
//***************************************************************************
void Print5Sets(/* in */ const int faceArr[NUM_FACES][MAX_EDGES], 
	        /* in */ const int vertArr[NUM_VERT][MAX_DIH],
        	/* in */ const int edgeList[MAX_CON_EDGE][2],
                /* in */ int configFaceNum,
		/* in */ int configVertNum,
	        /* in */ int configEdgeNum,
                /* inout */ ofstream& outFile)

{
  const int MAX_5_CONFIG=100;
  int counter;
  int counter2;
  int counter3;
  int counter4;
  int i;
  int j;
  int index=0;
  int sum;

  int list[MAX_5_CONFIG][2];   //[vertex,number of faces]
  int list2[MAX_5_CONFIG][3];  //[vertex,vertex, num faces]

  for(i=0;i<MAX_5_CONFIG;i++)
    for(j=0;j<2;j++)
      list[i][j]=-1;

  for(i=0;i<MAX_5_CONFIG;i++)
    for(j=0;j<3;j++)
      list2[i][j]=-1;
  

  for(counter=0;counter<configVertNum;counter++)
    {
      if(vertArr[counter][9]==5 && vertArr[counter][10]==0 && 
	 vertArr[counter][11]==0)
	{
	  list[index][0]=vertArr[counter][0];
	  list[index][1]=5;
	  index++;
	 
	}
    }	  
  
  /*cout<<endl<<"LIST 1"<<endl;
  for(i=0;i<15;i++)
    {
      for(j=0;j<2;j++)
	cout<<list[i][j]<<" ";
      cout<<" ";
    }*/
    
  for(counter=0;counter<index;counter++)
     {
       
       outFile<<"Vs"<<list[counter][0]<<" : ";
       for(counter2=2;counter2<7;counter2++)
	 {
	   outFile<<"+sigma{"<<vertArr[(list[counter][0])-1][counter2]
		  <<"} + 0.42755 solid{"
		  <<vertArr[(list[counter][0])-1][counter2]<<"}";
		  
	 }
       outFile<<" + 0.004 c <= 0"<<endl;

     }


   int index2=0;
   int a,b,x,y;
 
   //cout<<endl<<"Index="<<index<<endl;
   outFile<<"\\ 2 sets"<<endl<<endl;
   for(counter=0;counter<index-1;counter++)
     {
       for(counter2=counter+1;counter2<index;counter2++)
	 {
	   a=list[counter][0]; //cout<<a<<" ";
	   b=list[counter2][0]; // cout<<b<<" ";
	   if(a<b)
	     {
	       x=a; 
	       y=b;
	     }
	   else if(b<a)
	     {
	       x=b;
	       y=a;
	     }
	   else if(a==b) cout<<"** Error in 5set set2 edge check**"<<endl;
	   for(counter3=0;counter3<configEdgeNum;counter3++)
	     {
	       if( edgeList[counter3][0]==x && edgeList[counter3][1]==y) 
		 {
		   list2[index2][0]=x; 
		   list2[index2][1]=y;
		   list2[index2][2]=8;
		   index2++;
		 }
	     }
	 }
     }
   
   /*cout<<endl<<endl<<"List 2:"<<endl;
    
   for(counter=0;counter<index2;counter++)
      {
	for(counter2=0;counter2<3;counter2++)
	  {cout<<list2[counter][counter2]<<" ";}
	cout<<"   ";
      }
      */ 
   

  return;

}





//***************************************************************************

//***************************************************************************
	
void PrintBreakQuads(/* in */ const int faceArr[NUM_FACES][MAX_EDGES], 
                     /* in */ int configFaceNum,
                     /* inout */ ofstream& outFile)
{
  int counter;
  int counter2;
  int i;
  int j;
  
  outFile<<endl<<"\\ Division of quads into simplices"<<endl<<endl;

  for(counter=0;counter<configFaceNum;counter++)
    {
      if(faceArr[counter][1] != 4) continue;
      outFile<<"Fs"<<faceArr[counter][0]<<"a"<<" : solid{"<<faceArr[counter][0]
	     <<"} - solid{"<<faceArr[counter][0]<<"}"<<faceArr[counter][2]
	     <<" - solid{"<<faceArr[counter][0]<<"}"<<faceArr[counter][4]
	     <<" = 0"<<endl;
      outFile<<"Fs"<<faceArr[counter][0]<<"b"<<" : solid{"<<faceArr[counter][0]
	     <<"} - solid{"<<faceArr[counter][0]<<"}"<<faceArr[counter][3]
	     <<" - solid{"<<faceArr[counter][0]<<"}"<<faceArr[counter][5]
	     <<" = 0"<<endl;
      outFile<<"F"<<faceArr[counter][0]<<"a"<<" : sigma{"<<faceArr[counter][0]
	     <<"} - sigma{"<<faceArr[counter][0]<<"}"<<faceArr[counter][2]
	     <<" - sigma{"<<faceArr[counter][0]<<"}"<<faceArr[counter][4]
	     <<" = 0"<<endl;
      outFile<<"F"<<faceArr[counter][0]<<"b"<<" : sigma{"<<faceArr[counter][0]
	     <<"} - sigma{"<<faceArr[counter][0]<<"}"<<faceArr[counter][3]
	     <<" - sigma{"<<faceArr[counter][0]<<"}"<<faceArr[counter][5]
	     <<" = 0"<<endl;      
      
      outFile<<"F"<<faceArr[counter][0]<<",1"<<" : solid{"<<faceArr[counter][0]<<"}"
             <<faceArr[counter][2]<<"-"<<" dih{"<<faceArr[counter][0]<<"}"
             <<faceArr[counter][2]<<"-"<<" dih{"<<faceArr[counter][0]<<"}"
             <<faceArr[counter][3]<<".s("<<faceArr[counter][2]<<") - dih{"
	     <<faceArr[counter][0]<<"}"<<faceArr[counter][5]<<".s("
	     <<faceArr[counter][2]<<") + pi = 0"<<endl;
      outFile<<"F"<<faceArr[counter][0]<<",2"<<" : dih{"<<faceArr[counter][0]<<"}"
             <<faceArr[counter][2]<<"-"<<" dih{"<<faceArr[counter][0]<<"}"
             <<faceArr[counter][2]<<".s("<<faceArr[counter][3]<<") - dih{"
	     <<faceArr[counter][0]<<"}"<<faceArr[counter][2]<<".s("
	     <<faceArr[counter][5]<<") = 0"<<endl;
      outFile<<"F"<<faceArr[counter][0]<<",3"<<" : solid{"<<faceArr[counter][0]<<"}"
             <<faceArr[counter][3]<<"-"<<"dih{"<<faceArr[counter][0]<<"}"
             <<faceArr[counter][3]<<"-"<<"dih{"<<faceArr[counter][0]<<"}"
             <<faceArr[counter][2]<<".s("<<faceArr[counter][2]<<") - dih{"
	     <<faceArr[counter][0]<<"}"<<faceArr[counter][4]<<".s("
	     <<faceArr[counter][2]<<") + pi = 0"<<endl;
      outFile<<"F"<<faceArr[counter][0]<<",4"<<" : dih{"<<faceArr[counter][0]<<"}"
             <<faceArr[counter][3]<<"-"<<" dih{"<<faceArr[counter][0]<<"}"
             <<faceArr[counter][3]<<".s("<<faceArr[counter][2]<<") - dih{"
	     <<faceArr[counter][0]<<"}"<<faceArr[counter][3]<<".s("
	     <<faceArr[counter][4]<<") = 0"<<endl;
      outFile<<"F"<<faceArr[counter][0]<<",5"<<" : solid{"<<faceArr[counter][0]<<"}"
             <<faceArr[counter][4]<<"-"<<"dih{"<<faceArr[counter][0]<<"}"
             <<faceArr[counter][4]<<"-"<<"dih{"<<faceArr[counter][0]<<"}"
             <<faceArr[counter][5]<<".s("<<faceArr[counter][4]<<") - dih{"
	     <<faceArr[counter][0]<<"}"<<faceArr[counter][3]<<".s("
	     <<faceArr[counter][4]<<") + pi = 0"<<endl;
      outFile<<"F"<<faceArr[counter][0]<<",6"<<" : dih{"<<faceArr[counter][0]<<"}"
             <<faceArr[counter][4]<<"-"<<" dih{"<<faceArr[counter][0]<<"}"
             <<faceArr[counter][4]<<".s("<<faceArr[counter][5]<<") - dih{"
	     <<faceArr[counter][0]<<"}"<<faceArr[counter][4]<<".s("
	     <<faceArr[counter][3]<<") = 0"<<endl;
      outFile<<"F"<<faceArr[counter][0]<<",7"<<" : solid{"<<faceArr[counter][0]<<"}"
             <<faceArr[counter][5]<<"-"<<"dih{"<<faceArr[counter][0]<<"}"
             <<faceArr[counter][5]<<"-"<<"dih{"<<faceArr[counter][0]<<"}"
             <<faceArr[counter][2]<<".s("<<faceArr[counter][5]<<") - dih{"
	     <<faceArr[counter][0]<<"}"<<faceArr[counter][4]<<".s("
	     <<faceArr[counter][5]<<") + pi = 0"<<endl;
      outFile<<"F"<<faceArr[counter][0]<<",8"<<" : dih{"<<faceArr[counter][0]<<"}"
             <<faceArr[counter][5]<<"-"<<" dih{"<<faceArr[counter][0]<<"}"
             <<faceArr[counter][5]<<".s("<<faceArr[counter][2]<<") - dih{"
	     <<faceArr[counter][0]<<"}"<<faceArr[counter][5]<<".s("
	     <<faceArr[counter][4]<<") = 0"<<endl;
    }
  outFile<<endl<<endl;
  return;
}

//***************************************************************************

//***************************************************************************

void PrintBoundSigma(/* in */ const int faceArr[NUM_FACES][MAX_EDGES],
                     /* in */ int configFaceNum,
                     /* inout */ ofstream& outFile)

{  
  int counter;
  const char TET_LB[10]="-infinity";
  const double TET_UB=0;
  const char QUAD_LB[10]="-infinity";
  const double QUAD_UB=0;

  outFile<<endl<<"\\Give bounds on scores"<<endl<<endl;
  for(counter=0;counter<configFaceNum;counter++)
    {
      if(faceArr[counter][1]==3)
        {
          outFile<<TET_LB<<" <= sigma{"<<faceArr[counter][0]<<"} <= "
                 <<TET_UB<<endl;
        }
      else
        {
          outFile<<QUAD_LB<<" <= sigma{"<<faceArr[counter][0]<<"} <= "
                 <<QUAD_UB<<endl;
        }
    }

  outFile<<endl<<endl;
  return;
}

//***************************************************************************

void PrintBoundSolid(/* in */ const int faceArr[NUM_FACES][MAX_EDGES], 
                     /* in */ int configFaceNum,
                     /* inout */ ofstream& outFile)
{
  int counter;
  const double LB=0;
  const double UB=12.57;
  
  outFile<<"\\Give bounds on solid angle (4pi < 12.57)"<<endl<<endl;
  for(counter=0;counter<configFaceNum;counter++)
    outFile<<LB<<" <= solid{"<<faceArr[counter][0]<<"} <= "<<UB<<endl;
  outFile<<endl<<endl;
  return;
}

//***************************************************************************


//***************************************************************************
void PrintBoundDih(/* in */ const int faceArr[NUM_FACES][MAX_EDGES],
/* in */ int configFaceNum,
/* inout */ ofstream& outFile)

{
  const double TET_LB=0.856147;
  const double TET_UB=1.886730;
  const double QUAD_LB=1.15242;
  const double QUAD_UB=3.25887;  // 2pi
  int counter;
  int counter2;

  for(counter=0;counter<configFaceNum;counter++)
    {
      for(counter2=2;counter2<faceArr[counter][1]+2;counter2++)
        {
          if(faceArr[counter][1]==3)
	    outFile<<TET_LB<<" <= dih{"<<faceArr[counter][0]<<"}"
		   <<faceArr[counter][counter2]<<" <= "<<TET_UB<<endl;
	  else if(faceArr[counter][1]==4)
	    outFile<<QUAD_LB<<" <= dih{"<<faceArr[counter][0]<<"}"
		   <<faceArr[counter][counter2]<<" <= "<<QUAD_UB<<endl;
	  
	  else outFile<<QUAD_LB<<" <= dih{"<<faceArr[counter][0]<<"}"
		      <<faceArr[counter][counter2]<<endl;

	}
    }
  outFile<<endl<<endl;

  return;
}

//***************************************************************************


void PrintBoundEL(/* in */ const int list[MAX_CON_EDGE][2],
		  /* in */ int numEdges,
                  /* in */ int configVertNum,
                  /* inout */ ofstream& outFile)



{
  int counter;
  
  
  outFile<<"\\Edge length bounds"<<endl<<endl;
  
  for(counter=0;counter<configVertNum;counter++)    
    outFile<<TET_LB<<" <= y("<<counter+1<<") <= "<<TET_UB<<endl;

  for(counter=0;counter<numEdges;counter++)
    {
      outFile<< TET_LB << " <= y(" << list[counter][0] << "," << list[counter][1]
	     << ") <= " << TET_UB <<endl;  	
    }
  
  outFile<<endl<<endl;
  return;
} 

//***************************************************************************


void PrintBrokenStuff(/* in */ const int faceArr[NUM_FACES][MAX_EDGES], 
                      /* in */ int configFaceNum,
                      /* inout */ ofstream& outFile)
{
  int counter;
  int counter2;
  int i;
  int j;  

  for(counter=0;counter<configFaceNum;counter++)
    {
      if(faceArr[counter][1] != 4) continue;
      for(counter2=2;counter2<5+1;counter2++)
	{
	  outFile<<"solid{"<<faceArr[counter][0]<<"}"
		 <<faceArr[counter][counter2]<<" free"<<endl;
	}
      for(counter2=2;counter2<5+1;counter2++)
	{
	  outFile<<"sigma{"<<faceArr[counter][0]<<"}"
		 <<faceArr[counter][counter2]<<" free"<<endl;
	}
      outFile<<"dih{"<<faceArr[counter][0]<<"}"<<faceArr[counter][2]<<".s("
	     <<faceArr[counter][3]<<") free"<<endl;
      outFile<<"dih{"<<faceArr[counter][0]<<"}"<<faceArr[counter][2]<<".s("
	     <<faceArr[counter][5]<<") free"<<endl;
      outFile<<"dih{"<<faceArr[counter][0]<<"}"<<faceArr[counter][3]<<".s("
	     <<faceArr[counter][2]<<") free"<<endl;
      outFile<<"dih{"<<faceArr[counter][0]<<"}"<<faceArr[counter][3]<<".s("
	     <<faceArr[counter][4]<<") free"<<endl;
      outFile<<"dih{"<<faceArr[counter][0]<<"}"<<faceArr[counter][4]<<".s("
	     <<faceArr[counter][3]<<") free"<<endl;
      outFile<<"dih{"<<faceArr[counter][0]<<"}"<<faceArr[counter][4]<<".s("
	     <<faceArr[counter][5]<<") free"<<endl;
      outFile<<"dih{"<<faceArr[counter][0]<<"}"<<faceArr[counter][5]<<".s("
	     <<faceArr[counter][4]<<") free"<<endl;
      outFile<<"dih{"<<faceArr[counter][0]<<"}"<<faceArr[counter][5]<<".s("
	     <<faceArr[counter][2]<<") free"<<endl;
    }
  outFile<<endl<<endl;
  return;
}

//***************************************************************************





//***************************************************************************
//
//                           PRINT JAVA GRAPH 
//
//***************************************************************************


void PrintGraph(/* in */ const int faceArr[NUM_FACES][MAX_EDGES], 
                /* in */ int configFaceNum,
                /* inout */ ofstream& outFile)
{
  int i;
  int j;
  
  
  // Prints Java code for graph.  (Used for searches)

  outFile<<endl<<endl<<endl<<" # \" 0 "<<configFaceNum<<" ";
  for(i=0;i<configFaceNum;i++)
    {
      
      outFile<<faceArr[i][1]<<" ";
      for(j=2;j<faceArr[i][1]+2;j++)
      {
	outFile<<faceArr[i][j]-1<<" ";
      }
     
    }
  outFile<<"\","<<endl;
  

  // Prints Java code reindexed beginning at 1.  (Used for drawing.)
  
  outFile<<endl<<endl<<endl<<" # \" 0 "<<configFaceNum<<" ";
  for(i=0;i<configFaceNum;i++)
    {
      
      outFile<<faceArr[i][1]<<" ";
      for(j=2;j<faceArr[i][1]+2;j++)
      {
	outFile<<faceArr[i][j]<<" ";	
      }
     
    }
  outFile<<"\","<<endl;




  return;
}

//***************************************************************************
//
//                  Dihedral with edge length ineqs
//
//***************************************************************************

void PrintDihEdge(/* in */ const int faceArr[NUM_FACES][MAX_EDGES], 
                  /* in */ int configFaceNum,
                  /* inout */ ofstream& outFile)

{
  const int NUM_DE_INEQS=6;

  int counter;
  int counter2;
  int a;int b;int c;
  int i;

  double dihEdge[NUM_DE_INEQS][8]={
    {-1,0.237,-0.372,-0.372,0.708,-0.372,-0.372,2.3169},
    {-1,0.237,-0.363,-0.363,0.688,-0.363,-0.363,2.2849},
    {1,-0.505,0.152,0.152,-0.766,0.152,0.152,0.09503},
        
                                    };

  outFile<<"\\ Dihedral, edge length ineqs"<<endl<<endl;

  for(counter=0;counter<configFaceNum;counter++)
    {
      if(faceArr[counter][1]!=3) continue;
      a=faceArr[counter][2];b=faceArr[counter][3];c=faceArr[counter][4];
      for(counter2=0;counter2<NUM_DE_INEQS;counter2++)
	{
	  outFile<<"DihEdgeE"<<counter2+1<<"F"<<faceArr[counter][0]
		 <<"D"<<a<<" : "
	         <<dihEdge[counter2][0]<<" dih{"<<faceArr[counter][0]
		 <<"}"<<a<<" "
		 <<SignOf(dihEdge[counter2][1])<<" "
		 <<fabs(dihEdge[counter2][1])
	         <<" y("<<a<<") "<<SignOf(dihEdge[counter2][2])<<" "
		 <<fabs(dihEdge[counter2][2])<<" y("<<b<<") "
		 <<SignOf(dihEdge[counter2][3])<<" "
		 <<fabs(dihEdge[counter2][3])<<" y("<<c<<") "
		 <<SignOf(dihEdge[counter2][4])<<" "
		 <<fabs(dihEdge[counter2][4])<<" y("<<b<<","<<c<<") "
		 <<SignOf(dihEdge[counter2][5])<<" "
		 <<fabs(dihEdge[counter2][5])<<" y("<<a<<","<<b<<") "
		 <<SignOf(dihEdge[counter2][6])<<" "
		 <<fabs(dihEdge[counter2][6])<<" y("<<a<<","<<c<<")"
		 <<SignOf(dihEdge[counter2][7])
	         <<fabs(dihEdge[counter2][7])<<" c <= 0"<<endl;

	  outFile<<"DihEdgeE"<<counter2+1<<"F"<<faceArr[counter][0]
		 <<"D"<<b<<" : "
	         <<dihEdge[counter2][0]<<" dih{"<<faceArr[counter][0]
		 <<"}"<<b<<" "
		 <<SignOf(dihEdge[counter2][1])<<" "
		 <<fabs(dihEdge[counter2][1])
	         <<" y("<<b<<") "<<SignOf(dihEdge[counter2][2])<<" "
		 <<fabs(dihEdge[counter2][2])<<" y("<<a<<") "
		 <<SignOf(dihEdge[counter2][3])<<" "
		 <<fabs(dihEdge[counter2][3])<<" y("<<c<<") "
		 <<SignOf(dihEdge[counter2][4])<<" "
		 <<fabs(dihEdge[counter2][4])<<" y("<<a<<","<<c<<") "
		 <<SignOf(dihEdge[counter2][5])<<" "
		 <<fabs(dihEdge[counter2][5])<<" y("<<a<<","<<b<<") "
		 <<SignOf(dihEdge[counter2][6])<<" "
		 <<fabs(dihEdge[counter2][6])<<" y("<<b<<","<<c<<")"
		 <<SignOf(dihEdge[counter2][7])
	         <<fabs(dihEdge[counter2][7])<<" c <= 0"<<endl;

	  outFile<<"DihEdgeE"<<counter2+1<<"F"<<faceArr[counter][0]
		 <<"D"<<c<<" : "
	         <<dihEdge[counter2][0]<<" dih{"<<faceArr[counter][0]
		 <<"}"<<c<<" "
		 <<SignOf(dihEdge[counter2][1])<<" "
		 <<fabs(dihEdge[counter2][1])
	         <<" y("<<c<<") "<<SignOf(dihEdge[counter2][2])<<" "
		 <<fabs(dihEdge[counter2][2])<<" y("<<a<<") "
		 <<SignOf(dihEdge[counter2][3])<<" "
		 <<fabs(dihEdge[counter2][3])<<" y("<<b<<") "
		 <<SignOf(dihEdge[counter2][4])<<" "
		 <<fabs(dihEdge[counter2][4])<<" y("<<a<<","<<b<<") "
		 <<SignOf(dihEdge[counter2][5])<<" "
		 <<fabs(dihEdge[counter2][5])<<" y("<<c<<","<<b<<") "
		 <<SignOf(dihEdge[counter2][6])<<" "
		 <<fabs(dihEdge[counter2][6])<<" y("<<a<<","<<c<<")"
		 <<SignOf(dihEdge[counter2][7])
	         <<fabs(dihEdge[counter2][7])<<" c <= 0"<<endl;


	}
    }
  
  outFile<<endl<<endl;
  return;
}

//***************************************************************************

void PrintVolEdge(/* in */ const int faceArr[NUM_FACES][MAX_EDGES], 
                  /* in */ int configFaceNum,
                  /* inout */ ofstream& outFile)

{
  const int NUM_VE_INEQS=6;
  
  int i,j;
  
  double volEdge[NUM_VE_INEQS][9]={
    {0,-1,-0.245,-0.245,-0.245,0.063,0.063,0.063,1.6432},
    {0,-1,-0.3798,-0.3798,-0.3798,0.198,0.198,0.198,1.642},
    {0,1,0.151,0.151,0.151,-0.323,-0.323,-0.323,0.4807},
    {1,0.42755,0.0392,0.0392,0.0392,0.0101,0.0101,0.0101,-0.2958},
    {1,0,-0.107,-0.107,-0.107,0.116,0.116,0.116,0.1817},
    {1,0,-0.0623,-0.0623,-0.0623,0.0722,0.0722,0.0722,0.1763}

    
    

  };
  
    outFile<<"\\ Vol, Sol edge ineqs"<<endl<<endl;

  for(i=0;i<configFaceNum;i++)
    {
      for(j=0;j<NUM_VE_INEQS;j++)
        {  	
	  if(faceArr[i][1]!=3) continue;
	  outFile<<"VE"<<i+1<<" : "<<volEdge[j][0]<<" sigma{"
	  <<faceArr[i][0]<<"} "
	  <<SignOf(volEdge[j][1])<<fabs(volEdge[j][1])
	  <<" solid{"<<faceArr[i][0]<<"} "<<SignOf(volEdge[j][2])
	  <<fabs(volEdge[j][2])<<" y("<<faceArr[i][2]<<") "
	  <<SignOf(volEdge[j][3])
	  <<fabs(volEdge[j][3])<<" y("<<faceArr[i][3]<<") "
	  <<SignOf(volEdge[j][4])
	  <<fabs(volEdge[j][4])<<" y("<<faceArr[i][4]<<") "
	  <<SignOf(volEdge[j][5])
	  <<fabs(volEdge[j][5])<<" y("<<faceArr[i][2]<<","
	  <<faceArr[i][3]<<") "
	  <<SignOf(volEdge[j][6])
	  <<fabs(volEdge[j][6])<<" y("<<faceArr[i][2]<<","<<faceArr[i][4]<<") "
	  <<SignOf(volEdge[j][7])
	  <<fabs(volEdge[j][7])<<" y("<<faceArr[i][3]<<","<<faceArr[i][4]<<") "
	  <<SignOf(volEdge[j][8])<<fabs(volEdge[j][8])<<" c <= 0"<<endl;
	}
    }

    return;
}

//***************************************************************************

void PrintQuadEdge(const int faceArr[NUM_FACES][MAX_EDGES],
                  int configFaceNum,
                  ofstream& outFile)
{
  int i,j;

  for(i=0;i<configFaceNum;i++)
    {
      if(faceArr[i][1]!=4) continue;
      outFile<<"sigma{"<<faceArr[i][0]<<"} -0.166 y("<<faceArr[i][2]
	     <<") -0.166 y("<<faceArr[i][3]<<") -0.166 y("<<faceArr[i][4]
	     <<") -0.166 y("<<faceArr[i][5]
	     <<") + 0.143 y("<<faceArr[i][2]<<","<<faceArr[i][3]
	     <<") +0.143 y("<<faceArr[i][3]<<","<<faceArr[i][4]
	     <<") +0.143 y("<<faceArr[i][4]<<","<<faceArr[i][5]
	     <<") +0.143 y("<<faceArr[i][2]<<","<<faceArr[i][5]
	     <<") + 0.77 c < 0"<<endl;
      // Note: Original c=0.774491

  }

  return;
}











//***************************************************************************

void PrintSquanderFace(const int faceArr[NUM_FACES][MAX_EDGES],
                  int configFaceNum,
                  ofstream& outFile)
{
  
  int i,j;
  
  outFile<<"\\ Squander Face inequalities"<<endl<<endl;
  for(i=0;i<configFaceNum;i++)
    {
      if(faceArr[i][1]==5)
	{	  
	  outFile<<"sigma{"<<faceArr[i][0]<<"} + 0.42755 solid{"
		 <<faceArr[i][0]<<"} + 0.076 c < 0"<<endl;
	}

      if(faceArr[i][1]==6)
	{	  
	  outFile<<"sigma{"<<faceArr[i][0]<<"} + 0.42755 solid{"
		 <<faceArr[i][0]<<"} + 0.121 c < 0"<<endl;
	}

      if(faceArr[i][1]==7)
	{	  
	  outFile<<"sigma{"<<faceArr[i][0]<<"} + 0.42755 solid{"
		 <<faceArr[i][0]<<"} + 0.166 c < 0"<<endl;
	}
    }

  outFile<<endl<<endl;
  return;
}

//***************************************************************************

//Note:  This procedure is only temporary.  It is not correct.  If there is overlap, it is false.

void PrintXST(const int vertArr[NUM_VERT][MAX_DIH],
                  int configVertNum,
                  ofstream& outFile)
{
  int i,j;

  outFile<<endl<<endl<<"// XSTC stuff"<<endl<<endl;

  for(i=0;i<configVertNum;i++)
    {
      if(vertArr[i][1] != 5) continue;
      if(vertArr[i][9]!=4 || vertArr[i][10]!=0 || vertArr[i][11]!=1) 
	continue;
      for(j=2;j<7;j++)
	{
	  outFile<<"+ sigma{"<<vertArr[i][j]<<"} + 0.42755 solid{"
		 <<vertArr[i][j]<<"}";
	}
      outFile<<" + 0.092 c < 0"<<endl;
    }

  return;
}
