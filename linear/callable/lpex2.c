/* lpex2.c based on example programs lpex2.c in cplex manual */
#define CPX_PROTOTYPE_ANSI
#include "cplex.h"
#include <ctype.h>
#include <stdlib.h>
#include <string.h>

static int readproblemdata(CPXENVptr env, char *filename, char **probname_p,
		int *numcols_p, int *numrows_p, int *objsen_p,
		double **obj_p, double **rhs_p, char **sense_p,
		int **matbeg_p, int **matcnt_p,
		int **matind_p, double **matval_p,
		double ** lb_p, double ** ub_p, double **rngval_p,
		char ***colname_p,char **colnamestore_p,
		char ***rowname_p,char **rownamestore_p,
		char *colspace_p, int *rowspace_p, int*nzspace_p,
		unsigned *colnamespace_p, unsigned *rownamespace_p);

static void
	freereaddata (char **probname_p, double **obj_p, double **rhs_p,
		char **sense_p, int **matbeg_p, int **matcnt_p,
		int **matind_p, double ** matval_p, double ** lb_p,
		double **ub_p, double **rngval_p,
		char ***colname_p, char **colnamestore_p,
		char ***rowname_p, char **rownamestore_p),
	free_and_null(char **ptr), usage(char *progname);


int main(int argc,char *argv[])
	{
	char *probname = NULL;
	int numcols=0;
	int numrows=0;
	int objsen=0;
	double *obj=NULL;
	double *rhs = NULL;
	char *sense = NULL;
	int *matbeg=NULL;
	int *matcnt = NULL;
	int *matind=NULL;
	double *matval=NULL;
	double *lb = NULL;
	double *ub = NULL;
	double *rngval = NULL;
	char **colname = NULL;
	char *colnamestore = NULL;
	char **rowname = NULL;
	char *rownamestore=NULL;
	int colspace =0;
	int rowspace =0;
	int nzspace = 0;
	unsigned colnamespace =0;
	unsigned rownamespace =0;

	int solstat;
	double objval;
	double *x = NULL;
	int *cstat = NULL;
	int *rstat = NULL;

	CPXENVptr env=NULL;
	CPXLPptr lp=NULL;
	int status;
	int i,j;
	int cur_numrows,cur_numcols;
	char *basismsg;

	if (( argc !=3) 
		||(strchr("odthb",argv[2][0])==NULL) ) 
		{ usage(argv[0]); goto TERMINATE;
		}

	/* page 329 */
	env = CPXopenCPLEX (&status);
	if (env==NULL) { char errmsg[1024]; fprintf(stderr,"cplex envir\n");
	CPXgeterrorstring(env,status,errmsg);
	fprintf(stderr,"%s",errmsg); goto TERMINATE;
	}

	status = CPXsetintparam(env,CPX_PARAM_SCRIND,CPX_ON);
	if (status!=0) { fprintf(stderr,"screen %d\n",status); goto TERMINATE;}

	status = readproblemdata(env,argv[1],&probname,&numcols,&numrows,
		&objsen,&obj,&rhs,&sense,&matbeg,&matcnt,&matind,&matval,
		&lb,&ub,&rngval,&colname,&colnamestore,&rowname,&rownamestore,
		&colspace,&rowspace,&nzspace,&colnamespace,&rownamespace);
	if (status) goto TERMINATE;

	lp = CPXloadlpwnames(env,probname,numcols,numrows,objsen,
		obj,rhs,sense,matbeg,matcnt,matind,matval,lb,ub,rngval,
		colspace,rowspace,nzspace,colnamespace,rownamespace);

	if (lp==NULL) { fprintf(stderr,"no load.\n"); goto TERMINATE; }
	
	/* page 330 */
	switch (argv[2][0]) {
		case 'o' : status=CPXoptimize(env,lp); break;
		case 't':case 'd': status=CPXdualopt(env,lp); break;
		case 'b': status=CPXbaropt(env,lp); break;
		case 'h': status=CPXhybbaropt(env,lp,'o'); break;
		default: status=-1; break; 
		}

	if (status) { fprintf(stderr,"no opt\n"); goto TERMINATE; }

	solstat = CPXgetstat(env,lp);
	status = CPXgetobjval(env,lp,&objval);

	if (status) { fprintf(stderr,"no obj\n"); goto TERMINATE; }
	printf("Solution status %d. Objective value %.10g\n",solstat,objval);

	cur_numcols=CPXgetnumcols(env,lp);
	cur_numrows=CPXgetnumrows(env,lp);

	cstat = (int *) malloc(cur_numcols*sizeof(int));
	rstat = (int *) malloc(cur_numrows*sizeof(int));
	x = (double *) malloc (cur_numcols*sizeof(double));

	if (cstat==NULL || rstat==NULL ) { fprintf(stderr,"no mem\n");
		goto TERMINATE; }

	CPXsetintparam(env,CPX_PARAM_SCRIND,CPX_OFF);
	status = CPXgetbase(env,lp,cstat,rstat);

	/* page 331 */
	CPXsetintparam(env,CPX_PARAM_SCRIND,CPX_ON);

	if ( status == CPXERR_NO_BASIS) {
		printf("No basis exists.\n"); free((char*)cstat);
		free ((char*) rstat); cstat=NULL;rstat=NULL; }
	else if (status) { fprintf(stderr,"basis %d\n",status);
		goto TERMINATE; }
	
	status = CPXgetx(env,lp,x,0,cur_numcols-1);
	if (status) { fprintf("stderr,"primal\n"); goto TERMINATE; }

	/* write out solution */

	TERMINATE:
	free_and_null((char **) &cstat);
	free_and_null((char **) &rstat);
	free_and_null((char **) &x);
	if (lp!=NULL) {
		status = CPXunloadprob(env,&lp);
		if (status) {
		fprintf("stderr,"unload\n",status); }
		}

	freereaddata(&probname,&obj,&rhs,&sense,&matbeg,
		&matcnt,&matind,&matval,&lb,&ub,
		&rngval,&colname,&colnamestore,
		&rowname,&rownamestore);

	if (env!=NULL) { status = CPXcloseCPLEX(&env);
	  if (status) {
		char errmsg[1024];
		fprintf(stderr,"Could not close\n");
		CPXgeterrorstring(env,status,errmsg);
		fprintf(stderr,"%s",errmsg);
		}
	}

} /* end main */




static int readproblemdata(CPXENVptr env, char *filename, char **probname_p,
		int *numcols_p, int *numrows_p, int *objsen_p,
		double **obj_p, double **rhs_p, char **sense_p,
		int **matbeg_p, int **matcnt_p,
		int **matind_p, double **matval_p,
		double ** lb_p, double ** ub_p, double **rngval_p,
		char ***colname_p,char **colnamestore_p,
		char ***rowname_p,char **rownamestore_p,
		char *colspace_p, int *rowspace_p, int*nzspace_p,
		unsigned *colnamespace_p, unsigned *rownamespace_p)

/* page 333 */
	{
	int namelen;
	char filetype[4];
	int i;
	int status = 0;
	}
