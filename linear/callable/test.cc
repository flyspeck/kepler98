#include <iomanip.h>
#include "cplex.h"
#include <string.h>
#include "cplexInit.h"

void loadAndSolve(char* filename,  CPXLPptr& lp );
void printVarNames(CPXLPptr lp);
void printVarValues(CPXLPptr lp,char* st);
void printInfo(CPXLPptr lp);

int printMenu(CPXLPptr lp)
	{
	cout << "\n\n";
	cout << " 0= quit, ";
	cout << " 1= solution, ";
	cout << " 2= y,";
	cout << " 3= dih,";
	cout << " 4= sol,"; 
	cout << " 5= sig,";
	cout << " 6= vars,";
	cout << " 7= info,";
	cout << "\n\nCPLEX-option-> ";
	int choice=0;
	cin >> choice;
	switch (choice)
			{
			case 0: return 0;
			case 1: printSolution(lp); break;
			case 2: printVarValues(lp,"y"); break;
			case 3: printVarValues(lp,"dih"); break;
			case 4: printVarValues(lp,"sol"); break;
			case 5: printVarValues(lp,"sig"); break;
			case 6: printVarNames(lp); break;
			case 7: printInfo(lp); break;
			default: return 0;
			}
	return 1;
	}


main()
	{
	CPXLPptr lp,lpheight;
	char str[80];
	CPXLPptr qrtet = loadConvexHullFinder();
	int addRowSpace = 1500; 
	cout << "rowSpace = " << addRowSpace << endl;
	//CPXLPptr qrtetWeak = loadWeakHullFinder();
	cout << "\n\nCPLEX-file-> ";
	cin >> str;
	//strcpy(str,"/tmp/Z/cplex.lp70");
	loadAndSolve(str,lp,addRowSpace);
	loadAndSolve(str,lpheight,addRowSpace+1);
	//double target = 0.44298; /*8pt*/
	double target = -5.550; // sean's target.
	cout << "target = -5.550 ";
	//cin >> target;
	cout << endl;
	printSolution(lp);
	faceData& d = loadFaceData(str);
	cout << "There are " << d.getVertexCount() << " vertices " << endl;
	cout << "There are " << d.getFaceCount() << " faces" << endl;
	addSigmaConstraint(lpheight,target,d);
	int i;
	int iter;
	double ymin[6]={2,2,2,2,2,2};
	double ymax[6]={2,2,2,2,2,2};
	for (iter=0;iter<addRowSpace;iter++)
		{	
		int face = (iter/5) % d.getFaceCount();
		if (0 == (iter % 100)) 
				heightData::setHeight(lpheight,d);
		setYRange(lp,face,ymin,ymax,d);
		if (d.getVertexCount(face)!=3) continue;

		trace& t = setTriangularFace(lp,face,d);
		cout << endl << endl << flush;
		cout << "{";
		for (i=0;i<6;i++)
				cout << t.getY(i) << ",";
		cout << "}, dih={" << flush;
		for (i=0;i<3;i++)
				cout << t.getDih(i) << ",";
		cout << "} S=" << t.getS() << endl;
		trace& tCoeff = setHullTrace(qrtet,t,ymin,ymax);
		cout << "face = " << face << endl;
		cout << "iteration = " << iter << endl;
		addNewIneq(lp,face,tCoeff,d,1);	
		addNewIneq(lpheight,face,tCoeff,d,0);
		/*ContinuePrompt*/{
		double x = getSolution(lp);
		if ((x< target /*8pt*/)|| (-1 == (iter % 10)))
			{
			char strx[1024];
			cout << "\n\nCPLEX-continue? -> ";
			cin >> strx;
			if ((strx[0] == 'n')||(strx[0]=='N')) break;
			}
		}
		}
	// save
	cout << "\n\nCPLEX-save? -> ";
	char ctrx[1024];
	cin >> ctrx;
	if ((ctrx[0]=='y')||(ctrx[0]=='Y'))
		{
		Dump(lp);
		cout << "CPLEX lp file written to " << "/tmp/dumpfile.lp" << endl << flush;
		}
	
	
	//while (printMenu(lp)) {};
	cout << "\n\n";
	}
