#include "cplex.h"

class faceData {

public:
	virtual int getVertexCount() const=0 ;			// 0..
	virtual int getFaceCount() const=0;
	virtual int getVertexCount(int face) const=0;	// 0..faceCount()-1.
	virtual int getVertexAt(int face,int pos) const=0; // 0..fC-1,0..2, return 0..vC-1.
};

class trace {

public:
	virtual double getY(int index)const=0; 		// 0..5
	virtual double getDih(int index)const=0;		// 0..2
	virtual double getS()const=0;
	virtual void setY(int index,double y)=0;
	virtual void setDih(int index,double d)=0;
	virtual void setS(double s)=0;
};

//class CPXLPptr;

void Dump(const CPXLPptr&);
void loadAndSolve(char* filename,  CPXLPptr& lp,int extrarow =0 );
faceData& loadFaceData(char* filename);

int yVariableLpIndex(const CPXLPptr& lp,int vertexNumber);
int yVariableLpIndex(const CPXLPptr&,int face, int pos, const faceData& d);
int dihVariableLpIndex(const CPXLPptr&,int face, int pos, const faceData& d);
int sigmaVariableLpIndex(const CPXLPptr&,int face, const faceData& d);

double varValue(const CPXLPptr&,int variableIndex);

class heightData {
public:
	static void setHeight(const CPXLPptr& lp,faceData& d);
	static double getHeight(const CPXLPptr& lp,int vertexNumber); // 0..vN-1.
	};

void addSigmaConstraint(const CPXLPptr&,double target,faceData& d);


void printSolution(CPXLPptr lp);
double getSolution(CPXLPptr lp);

trace& setTriangularFace(const CPXLPptr&,int face,const faceData& d);


CPXLPptr loadConvexHullFinder(); // absolute min
CPXLPptr loadWeakHullFinder();   // uses target and minimizes coefficients.

trace& setHullTrace(const CPXLPptr&,trace&,double ymin[],double ymax[]);
void numericallyAdjust(trace&,double ymin[],double ymax[]); // from ineq.cc variant.


void addNewIneq(const CPXLPptr& lp,int face,trace& t,faceData& d,int doOpt=1);

void setYRange(const CPXLPptr& lp,int face,double ymin[6],double ymax[6],
	faceData& d);

