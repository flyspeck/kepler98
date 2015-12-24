extern "C"
{
#include "cplex.h"
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
}
#include <iomanip.h>
#include "cplexInit.h"


class globalInit {
public:
	static CPXENVptr env;
	globalInit() { }

	~globalInit() { CPXcloseCPLEX(&env); cerr << "closing..." << endl; }
	};



int status;
CPXENVptr globalInit::env = /* page 329 */ CPXopenCPLEX (&status);

globalInit CX;

static int readproblemdata(CPXENVptr env, char *filename, char **probname_p,
        int *numcols_p, int *numrows_p, int *objsen_p,
        double **obj_p, double **rhs_p, char **sense_p,
        int **matbeg_p, int **matcnt_p,
        int **matind_p, double **matval_p,
        double ** lb_p, double ** ub_p, double **rngval_p,
        char ***colname_p,char **colnamestore_p,
        char ***rowname_p,char **rownamestore_p,
        int *colspace_p, int *rowspace_p, int*nzspace_p,
        unsigned *colnamespace_p, unsigned *rownamespace_p)
 
/* page 333 */
    {
    int namelen;
    int status = 0;
	int numrims = 0;
	int *freerowind  =0;
	int *rimtype = 0;
	int *rimbeg = 0;
	int *rimcnt = 0;
	int *rimind = 0;
	double *rimval = 0;
	char *dataname = 0;
	char *objname = 0;
	char *rhsname = 0;
	char *rngname = 0;
	char *bndname = 0;
	char **rimname = 0;
	char *rimnamestore = 0;
	int rimspace = 0;
	int rimnzspace =0;
	unsigned rimnamespace = 0;
	int *cstat = 0;
	int *rstat = 0;

	namelen = strlen(filename);
	*probname_p = (char*)malloc(namelen+1);
	if (!(*probname_p)) {
		cerr << "failed to get memory for problem.\n";
		status = -2;
		return 0;
		}

	strncpy(*probname_p,filename, namelen);
	(*probname_p)[namelen]= '\0';

	*colspace_p=0;
	*rowspace_p=0;
	*nzspace_p =0;
	*colnamespace_p=0;
	*rownamespace_p=0;
	status=CPXlpread(CX.env,filename,numcols_p,numrows_p,
			objsen_p,obj_p,rhs_p,sense_p,matbeg_p,matcnt_p,
			matind_p,matval_p,lb_p,ub_p,&objname,&rhsname,
			colname_p,colnamestore_p,rowname_p,rownamestore_p,
			colspace_p,rowspace_p,nzspace_p,
			colnamespace_p,rownamespace_p);
	if (status) { cerr << "Failed to read file " << filename << endl; }
	cout << "objective name = " << objname << endl;
	}

void Dump(const CPXLPptr& lp)
	{
	CPXlpwrite(CX.env,lp,"/tmp/dumpfile.lp");
	}
	


void loadAndSolve(char* filename, CPXLPptr& lp,int extrarow)
	{
	int status;

	char *probname = 0;
    int numcols=0;
    int numrows=0;
    int objsen=0;
    double *obj=0;
    double *rhs = 0;
    char *sense = 0;
    int *matbeg=0;
    int *matcnt = 0;
    int *matind=0;
    double *matval=0;
    double *lb = 0;
    double *ub = 0;
    double *rngval = 0;
    char **colname = 0;
    char *colnamestore = 0;
    char **rowname = 0;
    char *rownamestore=0;
    int colspace =0;
    int rowspace =0;
    int nzspace = 0;
    unsigned colnamespace =0;
    unsigned rownamespace =0;

    /* page 329 */

    // env = CPXopenCPLEX (&status);
	lp = 0;

    if (!CX.env) { char errmsg[1024]; 
		cerr << "Error: cplex initialization failed.\n"; 
		CPXgeterrorstring(CX.env,status,errmsg);
		cerr << errmsg; 
		exit(1);
		return ;
		}

    status = CPXsetintparam(CX.env,CPX_PARAM_SCRIND,CPX_ON);
    if (status!=0) { cout << "screen " << status << endl ;  return ;}
 
    status = readproblemdata(CX.env,filename,&probname,&numcols,&numrows,
        &objsen,&obj,&rhs,&sense,&matbeg,&matcnt,&matind,&matval,
        &lb,&ub,&rngval,&colname,&colnamestore,&rowname,&rownamestore,
        &colspace,&rowspace,&nzspace,&colnamespace,&rownamespace);
    if (status) { cout << "problem not read."; exit(1); return ; }

    lp = CPXloadlpwnames(CX.env,probname,numcols,numrows,objsen,
        obj,rhs,sense,matbeg,matcnt,matind,matval,lb,ub,rngval,
        colname,colnamestore,rowname,rownamestore,
		colspace,rowspace,nzspace,colnamespace,rownamespace);
    if (!lp) { cout << "Problem not loaded.\n"; exit(1); return ; }

	status = CPXreallocprob(CX.env,lp,
		&obj,&rhs,&sense,&matbeg,&matcnt,&matind,&matval,&lb,&ub,&rngval,
		&colname,&colnamestore,&rowname,&rownamestore,
		NULL,
		colspace,rowspace+extrarow,nzspace+10*extrarow,
		colnamespace,rownamespace+10*extrarow);
	if (status) { cout << "failed to reallocate\n" << status; exit(1); return; }
			

	status = CPXoptimize(CX.env,lp);
	if (status) { cerr << "Failed to optimize" << endl; exit(1); return ; }
	}

void printSolution(CPXLPptr lp)
	{
	double objval=0.0;
	int status = CPXgetobjval(CX.env,lp,&objval);
	if (status) { cerr << "Failed to obtain objective value" << endl;
			return ; }
	cout << "Optimal value is " << objval;
	cout <<  " (or " << objval/0.0553736456684637 << " pt)" << endl;
	}

double getSolution(CPXLPptr lp)
	{
	double objval=0.0;
	int status = CPXgetobjval(CX.env,lp,&objval);
	if (status) { cerr << "Failed to obtain objective value" << endl;
			return 0; }
	return objval;
	}

void printInfo(CPXLPptr lp)
	{
	cout << "Problem name = ";
	char buf[50];
	int status = CPXgetprobname(CX.env,lp,buf);
	if (status) cout << "(unknown)" << endl;
	else cout << buf << endl;
	cout << "Number of constraints = " << CPXgetnumrows(CX.env,lp) << endl;
	cout << "Number of variables = " << CPXgetnumcols(CX.env,lp) << endl;
	}
	

void printVarNames(CPXLPptr lp)
	{
	// print out variable names:
	int status ;
	char* cur_colnamestore=0;
	char** cur_colname = 0;
	int surplus;
	int cur_numcols = CPXgetnumcols(CX.env,lp);
	int cur_colnamespace = CPXgetcolnamespace(CX.env,lp);
	cur_colname = (char **) malloc(sizeof(char*)*cur_numcols);
	cur_colnamestore = (char*) malloc(cur_colnamespace);
	status = CPXgetcolname(CX.env,lp,cur_colname,cur_colnamestore,
			cur_colnamespace,&surplus,0,cur_numcols-1);
	int i;
	for (i=0;i<cur_numcols;i++) cout << cur_colname[i] << " ";
	cout << "\n";
	}

void printVarName(CPXLPptr lp,int VarNumber)
	{
	// print out variable names:
	int status ;
	char* cur_colnamestore=0;
	char** cur_colname = 0;
	int surplus;
	int cur_numcols = CPXgetnumcols(CX.env,lp);
	int cur_colnamespace = CPXgetcolnamespace(CX.env,lp);
	cur_colname = (char **) malloc(sizeof(char*)*cur_numcols);
	cur_colnamestore = (char*) malloc(cur_colnamespace);
	status = CPXgetcolname(CX.env,lp,cur_colname,cur_colnamestore,
			cur_colnamespace,&surplus,0,cur_numcols-1);
	cout << cur_colname[VarNumber];
	}

static void printRowNames(CPXLPptr lp)
	{
	// print out row names:
	int status ;
	char* cur_namestore=0;
	char** cur_name = 0;
	int surplus;
	int cur_num = CPXgetnumrows(CX.env,lp);
	int cur_namespace = CPXgetrownamespace(CX.env,lp);
	cur_name = (char **) malloc(sizeof(char*)*cur_num);
	cur_namestore = (char *)malloc (cur_namespace);
	status = CPXgetrowname(CX.env,lp,cur_name,cur_namestore,
			cur_namespace,&surplus,0,cur_num-1);
	int i;
	for (i=0;i<cur_num;i++) cout << cur_name[i] << " ";
	cout << "\n";
	}

static void reportSpace(CPXLPptr lp);
int varNameIndex(CPXLPptr lp,char* varName)
	{
	// match variable names:
	int status ;
	char* cur_colnamestore=0;
	char** cur_colname = 0;
	int surplus;
	int cur_numcols = CPXgetnumcols(CX.env,lp);
	int cur_colnamespace = CPXgetcolnamespace(CX.env,lp);
	cur_colname = (char **) malloc(sizeof(char*)*cur_numcols);
	cur_colnamestore = (char*) malloc(cur_colnamespace);
	status = CPXgetcolname(CX.env,lp,cur_colname,cur_colnamestore,
			cur_colnamespace,&surplus,0,cur_numcols-1);
	int i;
	for (i=0;i<cur_numcols;i++) 
		if (!strcmp(cur_colname[i],varName)) 
			{
			return i;
			}
	cout << "not found: [" << varName << "]";
	return 0;
	}

void printVarValues(CPXLPptr lp,char* st)
	{
	// print out variable names:
	int status ;
	char* cur_colnamestore=0;
	char** cur_colname = 0;
	int surplus;
	int cur_numcols = CPXgetnumcols(CX.env,lp);
	int cur_colnamespace = CPXgetcolnamespace(CX.env,lp);
	cur_colname = (char **) malloc(sizeof(char*)*cur_numcols);
	cur_colnamestore = (char*) malloc(cur_colnamespace);
	status = CPXgetcolname(CX.env,lp,cur_colname,cur_colnamestore,
			cur_colnamespace,&surplus,0,cur_numcols-1);
	int i;
	int modifier=2;
	for (i=0;i<cur_numcols;i++) 
			{
			int j=0;
			double x[1];
			int ok=1;
			for (;st[j];j++)
				if (st[j]!=cur_colname[i][j]) ok=0;
			if (ok) 
					{
					cout << cur_colname[i] << " = ";
					status = CPXgetx(CX.env,lp,x,i,i);
					cout << x[0] << ";   ";
					if (0== (modifier++ % 3)) cout << endl;
					}
					
			}
	cout << "\n";
	}

void addSigmaConstraint(const CPXLPptr& lp,double target,faceData& d)
	{
	int fc = d.getFaceCount();
	double rhs[1] = {target};
	char sense[1]= {'G'};
	int rmatbeg[1]= {0};
	int* rmatind = (int*)malloc(sizeof(int)*fc);
	double* rmatval = (double*)malloc(sizeof(double)*fc);
	int i;
	for (i=0;i<fc;i++) rmatind[i]= sigmaVariableLpIndex(lp,i,d);
	for (i=0;i<fc;i++) rmatval[i]= 1.0;
	int status=CPXaddrows(CX.env,lp,0,1,fc,rhs,sense,rmatbeg,rmatind,
		rmatval,0,0);
	if (status) { cerr << "failed to add sigma" << endl; return; }
	}

void setYRange(const CPXLPptr& lp,int face,double ymin[6],double ymax[6],
	faceData& d)
	{
	for (int i=0;i<6;i++) { ymin[i]=2.0; ymax[i]=2.51; }
	ymax[0]=heightData::getHeight(lp,d.getVertexAt(face,0));
	ymax[1]=heightData::getHeight(lp,d.getVertexAt(face,1));
	ymax[2]=heightData::getHeight(lp,d.getVertexAt(face,2));
	}

double varValues(CPXLPptr lp,int variableIndex)
	{
	double x[1];
	int status = CPXgetx(CX.env,lp,x,variableIndex,variableIndex);
	if (status) cerr << "varValues failed " << endl;
	return x[0];
	}

CPXLPptr loadConvexHullFinder()
	{
	CPXLPptr lp;
	loadAndSolve("qrtetDual.lp",lp);
	return lp;
	}

CPXLPptr loadWeakHullFinder()
	{
	CPXLPptr lp;
	loadAndSolve("qrtet.lp",lp);
	return lp;
	}

int dihVariableLpIndex(const CPXLPptr& lp,int face,int pos,
	const faceData& d)
	{
	int vnum = d.getVertexAt(face,pos);
	// form dih{ + face + } + vnum;
	char facestr[30];
	char vstr[30];
	facestr[0]='d'; facestr[1]='i'; facestr[2]='h'; facestr[3]='{';
	sprintf(facestr+4,"%d",face+1);
	strcat(facestr,"}");
	sprintf(vstr,"%d",vnum+1);
	strcat(facestr,vstr);
	return varNameIndex(lp,facestr);
	}

static int yVariableLpIndexTop
	(const CPXLPptr& lp,int face,int pos,const faceData& d)
	{
	char facestr[30];
	facestr[0]='y'; facestr[1]='(';
	int v[3]= {d.getVertexAt(face,0),d.getVertexAt(face,1),d.getVertexAt(face,2)};
	int m =100;
	int i;
	for (i=0;i<3;i++) if ((v[i]<m)&&(i!=pos-3)) m = v[i];
	int M=-1;
	for (i=0;i<3;i++) if ((v[i]>M)&&(i!=pos-3)) M = v[i];
	sprintf(facestr+2,"%d",1+m);
	strcat(facestr,",");
	sprintf(facestr+strlen(facestr),"%d",1+M);
	strcat(facestr,")");
	return varNameIndex(lp,facestr);
	}

int yVariableLpIndex(const CPXLPptr& lp,int face,int pos,const faceData& d)
	{
	if (d.getVertexCount(face)!=3) return -1;
	if ((pos>5)||(pos<0)) return -1;
	if (pos>2) return yVariableLpIndexTop(lp,face,pos,d);
	char facestr[30];
	facestr[0]='y'; facestr[1]='(';
	sprintf(facestr+2,"%d",1+d.getVertexAt(face,pos));
	strcat(facestr,")");
	return varNameIndex(lp,facestr);
	}

int yVariableLpIndex(const CPXLPptr& lp,int vertexNumber)
	{
	char facestr[30];
	facestr[0]='y'; facestr[1]='(';
	sprintf(facestr+2,"%d",1+vertexNumber);
	strcat(facestr,")");
	return varNameIndex(lp,facestr);
	}

double* heightArray; 

static void clearObjective(CPXLPptr& lp)
	{
	int cnt = CPXgetnumcols(CX.env,lp);
	int* index = (int*)malloc(sizeof(int)*cnt);
	int i;
	for (i=0;i<cnt;i++) index[i]=i;
	double* values = (double*)malloc(sizeof(double)*cnt);
	for (i=0;i<cnt;i++) values[i]=0.0;
	CPXchgobj(CX.env,lp,cnt,index,values);
	}
	

void heightData::setHeight(const CPXLPptr& lp,faceData& d)
	{
	static int initialized=0;
	int vc = d.getVertexCount();
	if (!initialized)
		{
		heightArray = (double*)malloc(sizeof(double)*vc);
		initialized++;
		}
	double values[1]={1.0};
	int j;
	for (j=0;j<vc;j++)
		{
		clearObjective(lp);
		int index[1] = {yVariableLpIndex(lp,j)};
		//cout << "\nnew optimization: index is " << index[0] << endl;
		//cout << "index name is ";
		//printVarName(lp,index[0]);
		//cout  << endl;

		CPXchgobj(CX.env,lp,1,index,values);
		int status = CPXoptimize(CX.env,lp);
		if (status) { cerr << "Failed to optimize" << endl; exit(1); return ; }
		double objval;
		status = CPXgetobjval(CX.env,lp,&objval);
		if (status) { cerr << "Failed to obtain objective value" << endl;
				return ; }
		heightArray[j]=objval;
		//cout << "height at vertex " << j << " = " << objval << endl;
		}
	cout << "heights= ";
	for (j=0;j<vc;j++) cout << "(" << j+1<< ","<< heightArray[j] << ") ";
	cout << endl;
	}

double heightData::getHeight(const CPXLPptr& lp,int i)
	{
	return heightArray[i];
	}

int sigmaVariableLpIndex(const CPXLPptr& lp,int face,const faceData& d)
	{
	char facestr[30];
	facestr[0]='s'; facestr[1]='i'; facestr[2]='g'; facestr[3]='m';
	facestr[4]='a'; facestr[5]='{';
	sprintf(facestr+6,"%d",1+face);
	strcat(facestr,"}");
	return varNameIndex(lp,facestr);
	}

class traceD : public trace {
private:
    double y[6];
    double d[3];
    double s;
public:
    double getY(int i)const { if ((i<0)||(i>5)) return 0; return y[i]; }
    double getDih(int i)const { if ((i<0)||(i>2)) return 0; return d[i]; }
    double getS()const { return s; }
    void setY(int i,double y0) { if ((i<0)||(i>5)) return; y[i]=y0; }
    void setDih(int i,double d0) { if ((i<0)||(i>2)) return; d[i]=d0; }
    void setS(double s0) { s = s0; }
    };
 
traceD globalTrace;

trace& setTriangularFace(const CPXLPptr& lp,int face,const faceData& d)
	{
	int i;
	for (i=0;i<6;i++) 
		globalTrace.setY(i,varValues(lp,yVariableLpIndex(lp,face,i,d)));
	for (i=0;i<3;i++)
		globalTrace.setDih(i,varValues(lp,dihVariableLpIndex(lp,face,i,d)));
	globalTrace.setS(varValues(lp,sigmaVariableLpIndex(lp,face,d)));
	return globalTrace;
	}

traceD td;
trace& setHullTrace(const CPXLPptr& lp,trace& t,double ymin[],double ymax[])
	{
	int index[10]= {0,1,2,3,4,5,6,7,8,9};
	char* names[10]= {"a1","a2","a3","a4","a5","a6","d1","d2","d3","b"};
	double values[10];
	int i;
	for (i=0;i<6;i++) values[varNameIndex(lp,names[i])]=t.getY(i);
	for (i=0;i<3;i++) values[varNameIndex(lp,names[6+i])]=t.getDih(i);
	values[varNameIndex(lp,names[9])]= 1.0;
	CPXchgobj(CX.env,lp,10,index,values);
	int status = CPXoptimize(CX.env,lp);
	if (status) { cerr << "Failed to optimize" << endl; return td; }
	printSolution(lp);
	cout << "{";
	for (i=0;i<10;i++) 
		{
		cout << varValues(lp,varNameIndex(lp,names[i])) << ",";
		}
	for (i=0;i<6;i++) td.setY(i,varValues(lp,varNameIndex(lp,names[i])));
	for (i=0;i<3;i++) td.setDih(i,varValues(lp,varNameIndex(lp,names[i+6])));
	td.setS(varValues(lp,varNameIndex(lp,names[9])));
	numericallyAdjust(td,ymin,ymax);
	return td;
	}

static void reportSpace(CPXLPptr lp)
	{
	int colspace_p=0,rowspace_p=0,nzspace_p=0;
	unsigned colnamespace_p=0,rownamespace_p=0;	

	status = CPXgetspace(CX.env,lp,&colspace_p,&rowspace_p,&nzspace_p,
		&colnamespace_p,&rownamespace_p);
	cout << "colspace = " << colspace_p << endl;
	cout << "rowspace = " << rowspace_p << endl;
	cout << "nzspace = " << nzspace_p << endl;
	cout << "colnamespace = " << colnamespace_p << endl;
	cout << "rownamespace = " << rownamespace_p << endl << flush;
	if (status) { cerr << "failed to analyze space " << endl; return ; }
	}

void addNewIneq(const CPXLPptr& lp,int face,trace& t,faceData& d,int doOpt)
	{
	double rhs[1] = {t.getS()};
	char sense[1]= {'L'};
	int rmatbeg[1]= {0};
	int rmatind[10];
	double rmatval[10];
	char* rowname[1]={"XX"};
	int i;
	for (i=0;i<6;i++) rmatind[i]= yVariableLpIndex(lp,face,i,d);
	for (i=0;i<3;i++) rmatind[i+6]= dihVariableLpIndex(lp,face,i,d);
	rmatind[9]= sigmaVariableLpIndex(lp,face,d);
	//for (i=0;i<10;i++) cout << rmatind[i] << " ";
	for (i=0;i<6;i++) rmatval[i]= -t.getY(i);
	for (i=0;i<3;i++) rmatval[i+6]= -t.getDih(i);
	rmatval[9]=1.0;
	//for (i=0;i<10;i++) cout << rmatval[i] << " ";
	int status=0;
	status=CPXaddrows(CX.env,lp,0,1,10,rhs,sense,rmatbeg,rmatind,
		rmatval,NULL,NULL );
	if (status) { cerr << "failed to add " << endl; return; }
	if (doOpt)
		{
		status = CPXoptimize(CX.env,lp);
		if (status) { cerr << "Failed to optimize" << endl; return ; }
		cout << "## " ; printSolution(lp);
		}
	}
