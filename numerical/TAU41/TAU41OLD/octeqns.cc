#include <iomanip.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "../numerical.h"
#include "../constants.h"
#include "../morefn.h"

extern "C" 
{
#include "cplex.h"
}

#define COLSPACE   30
#define ROWSPACE   9500
#define NZSPACE    285000

double fa(double y1,double y2,double y3,double y4,double y5,double y6,
	double fA[2],int type[4])
    {
	double g = (type[0]||type[1] ? octavor(y1,y2,y3,y4,y5,y6) :
		gamma(y1,y2,y3,y4,y5,y6) );
	return g 
	+ fA[0]*(dihedraly(y2,y1,y3,y5,y4,y6)+ dihedraly(y3,y1,y2,y6,y4,y5))
    + fA[1]/4.0;
    }

double fb(double y1,double y2,double y3,double y4,double y5,double y6,
	double fA[2],int type[4])
    {
	double g = (type[1]||type[2] ? octavor(y1,y2,y3,y4,y5,y6) :
		gamma(y1,y2,y3,y4,y5,y6) );
	return g 
    + fA[0]*(dihedraly(y2,y1,y3,y5,y4,y6))+
	+ fA[1]/4.0;
    }

double fc(double y1,double y2,double y3,double y4,double y5,double y6,
	double fA[2],int type[4])
    {
	double g = (type[2]||type[3] ? octavor(y1,y2,y3,y4,y5,y6) :
		gamma(y1,y2,y3,y4,y5,y6) );
	return g 
	+fA[1]/4.0;
    }

double fd(double y1,double y2,double y3,double y4,double y5,double y6,
	double fA[2],int type[4])
    {
	double g = (type[3]||type[0] ? octavor(y1,y2,y3,y4,y5,y6) :
		gamma(y1,y2,y3,y4,y5,y6) );
	return g 
	+ fA[0]*(dihedraly(y3,y1,y2,y6,y4,y5))
	+ fA[1]/4.0;
    }

double dpi(double y1,double y2,double y3,double y4,double y5,double y6)
    {
    return dihedraly(y1,y2,y3,y4,y5,y6)-global::pi2;
    }

static int indomain( double y[6],int type1,int type2)
    {
	int b1,b2;
	if (type1) b1 = (radf(y[0],y[1],y[5])>1.41421356);
	else b1 = (radf(y[0],y[1],y[5])<1.41421357);
	if (type2) b2 = (radf(y[0],y[2],y[4])>1.41421356);
	else b2 = (radf(y[0],y[2],y[4])<1.41421357);
	return (b1 && b2);
    }
 

static void addrow(int matcnt[],int index[29][10000],double value[29][10000],
	int* numrows_p,char sense[],double rhs[],double entries[30],char S)
	{
	sense[*numrows_p]= S;
	rhs[*numrows_p] = entries[29];
	int i;
	for (i=0;i<29;i++) if (entries[i]!=0)
		{
		value[i][matcnt[i]] = entries[i];
		index[i][matcnt[i]] = *numrows_p;
		matcnt[i]++;
		}
	(*numrows_p)++;
	if (*numrows_p > ROWSPACE)
		{
		cout << " row overflow; increase allocation " << flush;
		exit(0);
		}
	}

static void clear(double x[],int L)
	{
	for (int i=0;i<L;i++) x[i]=0;
	}

	//////////
	// Set the arrays matbeg, matind, matval
	//
static void synthesize(int matcnt[],int index[29][10000],double value[29][10000],
	int numrows,int matbeg[],int matind[],double matval[])
	{
	int i,j;
	matbeg[0]=0;
	for (i=1;i<29;i++) matbeg[i]=matbeg[i-1]+matcnt[i-1];
	for (i=0;i<29;i++) for (j=0;j<matcnt[i];j++)
		{
		//cout << "setting " << j+matbeg[i] << " " << value[i][j] << " " << index[i][j] << endl;
		matval[j+matbeg[i]] = value[i][j];
		matind[j+matbeg[i]] = index[i][j];
		}
	}

static void setentryA(double y[6],double entry[30],double fA[2],int type[4])
	{
	clear(entry,30);
	entry[28]= -1;
	entry[29]= -fa(y[0],y[1],y[2],y[3],y[4],y[5],fA,type);
	entry[0]=1; 
	entry[1]= y[0]; 
	entry[2]= y[1]; 
	entry[3]= y[2];
	entry[4]=y[4]; 
	entry[5]=y[5]; 
	entry[6]= dpi(y[0],y[1],y[2],y[3],y[4],y[5]);
	}

static void setentryB(double y[6],double entry[30],double fA[2],int type[4])
	{
	clear(entry,30);
	entry[28]= -1;
	entry[29]= -fb(y[0],y[1],y[2],y[3],y[4],y[5],fA,type);
	entry[7+0]=1; 
	entry[7+1]= y[0]; 
	entry[7+2]= y[1]; 
	entry[7+3]= y[2];
	entry[7+4]=y[4]; 
	entry[7+5]=y[5]; 
	entry[7+6]= dpi(y[0],y[1],y[2],y[3],y[4],y[5]);
	}

static void setentryC(double y[6],double entry[30],double fA[2],int type[4])
	{
	clear(entry,30);
	entry[28]= -1;
	entry[29]= -fc(y[0],y[1],y[2],y[3],y[4],y[5],fA,type);
	entry[14+0]=1; 
	entry[14+1]= y[0]; 
	entry[14+2]= y[1]; 
	entry[14+3]= y[2];
	entry[14+4]=y[4]; 
	entry[14+5]=y[5]; 
	entry[14+6]= dpi(y[0],y[1],y[2],y[3],y[4],y[5]);
	}

static void setentryD(double y[6],double entry[30],double fA[2],int type[4])
	{
	clear(entry,30);
	entry[28]= -1;
	entry[29]= -fd(y[0],y[1],y[2],y[3],y[4],y[5],fA,type);
	entry[21+0]=1; 
	entry[21+1]= y[0]; 
	entry[21+2]= y[1]; 
	entry[21+3]= y[2];
	entry[21+4]=y[4]; 
	entry[21+5]=y[5]; 
	entry[21+6]= dpi(y[0],y[1],y[2],y[3],y[4],y[5]);
	}

double edge(double y0,double z0) // eta(sqrt(return),sqrt(y),sqrt(z))=2.
	{
	double y = y0*y0, z = z0*z0;
	double x = (4*y + 4*z - y*z + sqrt((8.0 - y)*y*z*(8.0 - z)))/4.0;
	return sqrt(x);
	}

int radAdjust(int k[6],int inv[6],double ymin[6],double ymax[13],double y[6]) 
	{
	int n4=0; int pos=-1;
	for (int j=0;j<6;j++) if (k[j]>2) { n4++; pos = j; }
	if (n4>1) return 0;
	if (n4<1) return 1;
	switch(pos) 
		{
		case 0 : 
			if (k[0]==3) y[0] = edge(y[1],y[5]); 
			else y[0] = edge(y[2],y[4]); break;
		case 1 : y[1] = edge(y[0],y[5]); break;
		case 2 : y[2] = edge(y[0],y[4]); break;
		case 4 : y[4] = edge(y[0],y[2]); break; 
		case 5 : y[5] = edge(y[0],y[1]); break; 
		default : break;
		}
	double eps = 1.0e-8;
	if ((y[pos]>ymax[inv[pos]]+eps)||(y[pos]<ymin[inv[pos]]-eps))
		return 0;
	return 1;
	}


static void
setproblemdata (char *probname, int *numcols_p, int *numrows_p,
                int *objsen_p, double *obj, double *rhs, char *sense,
                int *matbeg, int *matcnt, int *matind, double *matval,
                double *lb, double *ub,double ymin[13],double ymax[13],
				double fA[2],int type[4])
{
	static const numVar = 29;
   strcpy (probname, "example");
   *numcols_p   = numVar;
   *numrows_p   = 0;
   *objsen_p = 1;   /* The problem is minimization */

	int i;
	for (i=0;i<numVar;i++) obj[i]=0; obj[numVar-1]=1;

	for (i=0;i<numVar;i++) { lb[i]= -2000; ub[i]= 2000; }

	/* matcnt[0]...matcnt[28] keeps the current number of nonzeros in each col
	   index[0][r]..index[28][s] keeps the row num of value[0][r]..value[28][s]
	   We'll set matcnt, matind, matval after matcnt, index, value are done.
	*/

	for (i=0;i<numVar;i++) matcnt[i]=0;
	double value[29][10000];
	int index[29][10000];
	double entry[30]; 
	double y[6];

	int i0,i1,i2,i3,i4,i5;
	for (i0=0;i0<5;i0++)
    for (i1=0;i1<4;i1++)
    for (i2=0;i2<4;i2++)
    for (i3=0;i3<3;i3++)
    for (i4=0;i4<4;i4++)
    for (i5=0;i5<4;i5++)
        {
        y[0]= ymin[0] + (ymax[0]-ymin[0])*double(i0)/2.0;
        y[1]= ymin[1] + (ymax[1]-ymin[1])*double(i1)/2.0;
        y[2]= ymin[2] + (ymax[2]-ymin[2])*double(i2)/2.0;
        y[3]= ymin[3] + (ymax[3]-ymin[3])*double(i3)/2.0;
        y[4]= ymin[4] + (ymax[4]-ymin[4])*double(i4)/2.0;
        y[5]= ymin[5] + (ymax[5]-ymin[5])*double(i5)/2.0;
		int inv[6]={0,1,2,3,4,5};
		int k[6]={i0,i1,i2,i3,i4,i5};
		if (!radAdjust(k,inv,ymin,ymax,y)) continue;
        if (indomain(y,type[0],type[1] ))
			{
            setentryA(y,entry,fA,type); 
			addrow(matcnt,index,value,numrows_p,sense,rhs,entry,'L');
			}
		}

	for (i0=0;i0<5;i0++)
    for (i1=0;i1<4;i1++)
    for (i2=0;i2<4;i2++)
    for (i3=0;i3<3;i3++)
    for (i4=0;i4<4;i4++)
    for (i5=0;i5<4;i5++)
        {
        y[0]= ymin[0] + (ymax[0]-ymin[0])*double(i0)/2.0;
        y[1]= ymin[2] + (ymax[2]-ymin[2])*double(i1)/2.0;
        y[2]= ymin[7] + (ymax[7]-ymin[7])*double(i2)/2.0;
        y[3]= ymin[6] + (ymax[6]-ymin[6])*double(i3)/2.0;
        y[4]= ymin[8] + (ymax[8]-ymin[8])*double(i4)/2.0;
        y[5]= ymin[4] + (ymax[4]-ymin[4])*double(i5)/2.0;
		int inv[6]={0,2,7,6,8,4};
		int k[6]={i0,i1,i2,i3,i4,i5};
		if (!radAdjust(k,inv,ymin,ymax,y)) continue;
        if (indomain(y,type[1],type[2]))
			{
            setentryB(y,entry,fA,type); 
			addrow(matcnt,index,value,numrows_p,sense,rhs,entry,'L');
			}
		}

	for (i0=0;i0<5;i0++)
    for (i1=0;i1<4;i1++)
    for (i2=0;i2<4;i2++)
    for (i3=0;i3<3;i3++)
    for (i4=0;i4<4;i4++)
    for (i5=0;i5<4;i5++)
        {
        y[0]= ymin[0] + (ymax[0]-ymin[0])*double(i0)/2.0;
        y[1]= ymin[7] + (ymax[7]-ymin[7])*double(i1)/2.0;
        y[2]=ymin[11] + (ymax[11]-ymin[11])*double(i2)/2.0;
        y[3]= ymin[9] + (ymax[9]-ymin[9])*double(i3)/2.0;
        y[4]=ymin[10] + (ymax[10]-ymin[10])*double(i4)/2.0;
        y[5]= ymin[8] + (ymax[8]-ymin[8])*double(i5)/2.0;
		int inv[6]={0,7,11,9,10,8};
		int k[6]={i0,i1,i2,i3,i4,i5};
		if (!radAdjust(k,inv,ymin,ymax,y)) continue;
        if (indomain(y,type[2],type[3]))
			{
            setentryC(y,entry,fA,type); 
			addrow(matcnt,index,value,numrows_p,sense,rhs,entry,'L');
			}
		}

	for (i0=0;i0<5;i0++)
    for (i1=0;i1<4;i1++)
    for (i2=0;i2<4;i2++)
    for (i3=0;i3<3;i3++)
    for (i4=0;i4<4;i4++)
    for (i5=0;i5<4;i5++)
        {
        y[0]= ymin[0] + (ymax[0]-ymin[0])*double(i0)/2.0;
        y[1]=ymin[11] + (ymax[11]-ymin[11])*double(i1)/2.0;
        y[2]= ymin[1] + (ymax[1]-ymin[1])*double(i2)/2.0;
        y[3]=ymin[12] + (ymax[12]-ymin[12])*double(i3)/2.0;
        y[4]= ymin[5] + (ymax[5]-ymin[5])*double(i4)/2.0;
        y[5]=ymin[10] + (ymax[10]-ymin[10])*double(i5)/2.0;
		int inv[6]={0,11,1,12,5,10};
		int k[6]={i0,i1,i2,i3,i4,i5};
		if (!radAdjust(k,inv,ymin,ymax,y)) continue;
        if (indomain(y,type[3],type[0]))
			{
            setentryD(y,entry,fA,type); 
			addrow(matcnt,index,value,numrows_p,sense,rhs,entry,'L');
			}
		}

	// canned constraints:
	clear(entry,30);
	entry[0]=entry[7]=entry[14]=entry[21]=1;
	addrow(matcnt,index,value,numrows_p,sense,rhs,entry,'G');
	clear(entry,30);
	entry[1]=entry[8]=entry[15]=entry[22]=1;
	addrow(matcnt,index,value,numrows_p,sense,rhs,entry,'G');

	clear(entry,30);
	entry[6]= 1; entry[13]= -1;
	addrow(matcnt,index,value,numrows_p,sense,rhs,entry,'E');

	clear(entry,30);
	entry[13]= 1; entry[20]= -1;
	addrow(matcnt,index,value,numrows_p,sense,rhs,entry,'E');

	clear(entry,30);
	entry[20]= 1; entry[27]= -1;
	addrow(matcnt,index,value,numrows_p,sense,rhs,entry,'E');

	for (i=0;i<4;i++)
		{
		clear(entry,30);
		entry[7*i +3] = entry[7*((i+1) % 4) + 2] = 1;
		addrow(matcnt,index,value,numrows_p,sense,rhs,entry,'G');
		}
	for (i=0;i<4;i++)
		{
		clear(entry,30);
		entry[7*i +4] = entry[7*((i+1) % 4) + 5] = 1;
		addrow(matcnt,index,value,numrows_p,sense,rhs,entry,'G');
		}
	// end of canned constraints.
	synthesize(matcnt,index,value,*numrows_p,matbeg,matind,matval);



	}
	
	


static void setproblemdata (char *probname, int *numcols_p, int *numrows_p, 
                   int *objsen_p, double *obj, double *rhs, char *sense, 
                   int *matbeg, int *matcnt, int *matind, double *matval, 
                   double *lb, double *ub,double ymin[13],double ymax[13],
					double fA[2],int type[4]);



static int terminate(CPXENVptr env,CPXLPptr lp)
	{
	int status;

   if ( lp != NULL ) {
      status = CPXfreeprob (env, &lp);
      if ( status ) {
         fprintf (stderr, "CPXfreeprob failed, error code %d.\n", status);
      }
   }

   if ( env != NULL ) {
      status = CPXcloseCPLEX (&env);
      if ( status ) {
      char  errmsg[1024];
         fprintf (stderr, "Could not close CPLEX environment.\n");
         CPXgeterrorstring (env, status, errmsg);
         fprintf (stderr, "%s", errmsg);
      }
   }
   return (status);
	}


int separate(double ymin[13],double ymax[13],double c[29],
	double fA[2],int type[4])
{

char     probname[16];  /* Problem name is max 16 characters */
int      numcols;
int      numrows;
int      objsen;
double   obj[COLSPACE];
double   rhs[ROWSPACE];
char     sense[ROWSPACE];
int      matbeg[COLSPACE];
int      matcnt[COLSPACE];
int      matind[NZSPACE];
double   matval[NZSPACE];
double   lb[COLSPACE];
double   ub[COLSPACE];

int      solstat;
double   objval;
double   x[COLSPACE];
double   pi[ROWSPACE];
double   slack[ROWSPACE];
double   dj[COLSPACE];


/* Initialize the CPLEX environment */
int           status;
CPXENVptr     env = CPXopenCPLEX (&status);
CPXLPptr      lp = NULL;
int           i, j;
int           cur_numrows, cur_numcols;



   if ( env == NULL ) {
   char  errmsg[1024];
      fprintf (stderr, "Could not open CPLEX environment.\n");
      CPXgeterrorstring (env, status, errmsg);
      fprintf (stderr, "%s", errmsg);
      return terminate(env,lp);
   }

   setproblemdata (probname, &numcols, &numrows, &objsen, obj, 
                   rhs, sense, matbeg, matcnt, matind, matval, lb, ub,
					ymin,ymax,fA,type);

   /* Load the problem. */

   lp = CPXloadlp (env, probname, numcols, numrows, objsen, obj, rhs, 
                   sense, matbeg, matcnt, matind, matval,
                   lb, ub, NULL, COLSPACE, ROWSPACE, NZSPACE);
   if ( lp == NULL ) {
      fprintf (stderr, "Failed to load LP.\n");
      return terminate(env,lp);
   }

   status = CPXoptimize (env, lp);
   if ( status ) {
      fprintf (stderr, "Failed to optimize LP.\n");
      return terminate(env,lp);
   }
   status = CPXsolution (env, lp, &solstat, &objval, x, pi, slack, dj);
   if ( status ) {
      fprintf (stderr, "Failed to obtain solution.\n");
      return terminate(env,lp);
   }

   for (j = 0; j < 29; j++) c[j]=x[j];

	// status = CPXlpwrite (env, lp, "lpex1.lp");

	return terminate(env,lp);
}  /* END separate */
   

/*
static void round(double fA[2],int type[4],double c[29],double ymin[13],double
	ymax[13])
	{

	int i,j;
	for (i=0;i<13;i++) w[i]= (umax[i]-umin[i])/2.0;

	for (i=0;i< 8;i++) // 8192.
	{
	int k=1;
	for (j=0;j<13;j++) { ymin[j]= ((i/k) % 2 ? w[j] : 0 )+umin[j] ; k=2*k; }
	for (j=0;j<13;j++) { ymax[j]= ymin[j]+w[j]; }
	k=1;
	cout << i << "; ";
	for (j=0;j<4;j++) cout << type[j] << " ";
	cout << "; ";
	for (j=0;j<2;j++) cout << fA[j] << " ";
	cout << "; ";
	for (j=0;j<13;j++) { cout << ((i/k) % 2 ? 1 : 0 ) << " "; k=2*k; }
	cout << "; ";
	 separate(ymin,ymax,c,fA,type);
	//cout << i << " "<< c[28] << endl;
	for (int j=0;j<29;j++) cout << c[j] << " "; 
	cout << "; \n\n" << flush;
	}
	}
*/

static void test()
	{
	double ymin[13]={2.6,2,2,2,2,2,2,2,2,2,2,2,2};
	double ymax[13]=
		{2.82842712474619,2.2,2.2,2.2,2.2,2.2,2.2,2.2,2.2,2.2,2.2,2.2,2.2};
	int type[4]={0,0,0,0};
	double fA[2]={-3.0508,9.494};
	double c[29];
	separate(ymin,ymax,c,fA,type);
	for (int i=0;i<29;i++) cout << c[i] << " ";
	}


