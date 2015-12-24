//  copyright (c) 1997, Thomas C. Hales, all rights reserved.
#ifndef RECURSE_H
#define RECURSE_H

#include <iomanip.h>
#include <iostream.h>
extern "C"
{
#include <math.h>
#include <stdlib.h>
}
#include "error.h"
#include "interval.h"
#include "taylorInterval.h"

/* I'm in a verifying frenzy.  I'm not responding to outside noises.
	It's a heightened state of consciousness. S.McL.
*/


/*
CLASS
    cellOption
 
    A cellOption contains various parameter settings used to
	fine-tune the behavior of the recursive verification procedures in
	the class prove.
 
OVERVIEW TEXT
	The proof class is designed to prove many different inequalities
	in many different contexts.  It is necessary to customize 
	certain actions the recursive verification procedure.  This
	is done through various parameter settings.

AUTHOR
    Thomas C. Hales
*/

class cellOption {
public:
    enum X { silent, verbose };
    enum cellStatus { /*likely*/counterexample, cellPasses, inconclusive };
private:
	int iterationLimit;
	int iterationCount;
    double dihmax;
    int usingDihMax;
	int usingBigFace126;
    int chiShortCut;
	int startingIndex;
	int recursionDepth;
    X printingMode;
	double widthCutoff;
	int usingWidthCutoff;
	enum { Nskip = 30 };
	int skip_these[Nskip];
	int skipped;
public:
	//////////
	// 
	//
	void setRecursionDepth(int rd) { recursionDepth=rd; }

	//////////
	// This only has an effect if it is positive.
	//
	int getRecursionDepth() const { return recursionDepth; }

	//////////
	// quit after so many tries.
	void setIterationLimit(int lim) { iterationLimit=lim; }
	
	//////////
	// 
	int getIterationLimit() const { return iterationLimit; }

	//////////
	// 
	void resetIterationCount() { iterationCount=0; }

	//////////
	//
	int getIterationCount() const { return iterationCount; }

	//////////
	//
	void augmentIterationCount() { iterationCount++; }

	

	//////////
	// setSkipCases allows the user to enter an array that
	// specifies what input integers caseNumber will return
	// a nonzero value skip(caseNumber).
	//
	void setSkipCases(const int skiplist[],int len)
		{
		int i;
		if (len>Nskip) { error::message("skip array out of bounds"); len=Nskip; }
		for (i=0;i<len;i++) skip_these[i]=skiplist[i];
		skipped=len;
		cout << "\nWarning: skipping R-cases ";
		for (i=0;i<len;i++) cout << skiplist[i] << " ";
		cout << endl << flush;
		}
	//////////
	// skip is called by recursiveVerifierQ when there is 
	// dimension reduction.  Dimension reduction breaks a quad
	// inequality into many cases.  Before running each case,
	// recurseiveVerifierQ calls skip(casenumber) and omits the
	// case if a nonzero value is returned.  The nonzero return values
	// are set using setSkipCases.
	//
	int skip(int i) const
		{
		for (int j=0;j<skipped;j++)  if (i==skip_these[j]) return 1;
		return 0;
		}

	//////////
	// 
	void setWidthCutoff(double x) { widthCutoff=x; 
		usingWidthCutoff=1; }

	//////////
	//
	double getWidthCutoff() const { return widthCutoff; }

	//////////
	//
	int hasWidthCutoff() const { return usingWidthCutoff; }

	//////////
	//
    void setChiShortCut(int i) { chiShortCut=i; }

	//////////
	//
    int getChiShortCut() const { return chiShortCut; }

	//////////
	//
    void setDihMax(double theta) { dihmax=theta; usingDihMax=1;  }

	//////////
	//
    double getDihMax() { return dihmax; }

	//////////
	//
    int isUsingDihMax() { return usingDihMax; }

	//////////
	//
	void setBigFace126() { usingBigFace126=1; }

	//////////
	//
	int isUsingBigFace126() { return usingBigFace126; }

	//////////
	//
    void setPrintMode(X u) { printingMode = u; }

	//////////
	//
    X getPrintMode() { return printingMode; }

	//////////
	//
	void setStartingIndex(int i) { startingIndex=i; }

	//////////
	//
	int getStartingIndex() const { return startingIndex; }

	//////////
	//
    cellOption() { dihmax = 0.0; usingDihMax = 0; chiShortCut=0; 
				usingBigFace126=0;
                printingMode=verbose; usingWidthCutoff=0; 
				recursionDepth=0;
				startingIndex =0; for (int i=0;i<Nskip;i++) skip_these[i]=-1;
				iterationCount=0;
				iterationLimit=0;
				skipped=0; }
};




/*
CLASS
    prove

	Given a list of functions, prove that at least one is negative
	at every point of a given domain.
 
OVERVIEW TEXT
	This class gathers together a number of routines that
	take a taylorFunction (or more generally a list of taylorFunctions F)
	and proves that at least one of the functions F is negative
	on the given domain.

AUTHOR
    Thomas C. Hales
*/

class prove {
public:

	//////////
	// recursiveVerifier is the main inequality verification
	// procedure for simplices.  
	// It starts with a list of taylorFunctions I[],
	// and attempts to prove that at every point of the domain
	// x--z, at least one of the functions is negative.
	// Start with depth=0.  Each step of the recursion will
	// increase the depth by 1.
	// Start with x=x0, z=z0.  As the domain is subdivided,
	// x will increase and z will decrease.  x0, z0 remain fixed,
	// unless there is dimension reduction.
	//
static int recursiveVerifier(int depth,
	const domain& x,const domain& z,     /// current cell
	const domain& x0,const domain& z0,   // boundary
	const taylorFunction* I[],
	int count,cellOption& options);

	//////////
	// recursiveVerifierQ is the main inequality verification
	// procedure for quad clusters.
	// Each quad cluster is divided into two simplices A,B along
	// the shortest diagonal.  A list of taylorFunctions is
	// given IA, IB for the two simplices.  The recursiveVerifierQ
	// procedure attempts to show that at every point in the domain
	// xA--zA (on A), xB--zB (on B), there is an index for which
	// the sum of the values of IA[index] on A and IB[index] on B
	// is negative.  IA[index] or IB[index] is allowed to be positive
	// as long as the sum is negative.  Dimension reduction is
	// always used in recursiveVerifierQ.  Otherwise the dimensions
	// are too great to be handled by computer.  This means that
	// if unreducible taylorFunctions are used, the
	// results are unreliable.  
	//
	// The depth starts out at 0.
	//
static void recursiveVerifierQ(int depth, 
	const domain& xA,const domain& xB,
	const domain& zA,const domain& zB,
	const taylorFunction* IA[],const taylorFunction* IB[],int Nineq,
	cellOption& options);

	//////////
	// Prove that F<0 on the domain x--z.  This function calls
	// recursiveVerifier.
	//
static int generic(const domain& x,const domain& z,const taylorFunction& F);

	//////////
	// Prove that sigma+FlessScore <0 on the the domain x--z.
	// Two cases are considered according to whether the score
	// is compression or the analytic voronoi function.
	// qrtet calls recursiveVerifier twice.
	//
static int qrtet(const domain& x,const domain& z,const taylorFunction& FlessScore);

	//////////
	// Prove that sigma+FlessScore <0 on the the domain x--z, using
	// dimension reduction.
	// Two cases are considered according to whether the score
	// is compression or the analytic voronoi function.
	// It is assumed that FlessScore is reducible.  The procedure
	// keeps track of when the score sigma is reducible and acts
	// accordingly.
	// qrtetReduce calls recursiveVerifier twice.
	//
static int qrtetReduce
	(const domain& x,const domain& z,const taylorFunction& FlessScore);

	//////////
	//
static int flat(const domain& x,const domain& z,const taylorFunction& FlessScore);

	//////////
	//
static int flatVc(const domain& x,const domain& z,const taylorFunction& FlessScore);

	//////////
	//  This uses gamma or octavor depending on the context.
	//  This verification function is suitable for quad clusters.
	//
static int upright(const domain& x,const domain& z,const taylorFunction& FlessScore);

	//////////
	//
static int octa(const domain& x,const domain& z,const taylorFunction& FlessScore);

	//////////
	//
static int octaWithoutEta(const domain& x,const domain& z,
		const taylorFunction& FlessScore);

	//////////
	//
static int quadcluster(const domain& xA,const domain& xB,
		const domain& zA, const domain& zB, 
		const taylorFunction& fA, const taylorFunction& fB);
};


#endif
