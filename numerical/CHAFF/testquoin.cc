
#include <math.h>
#include "quoinfn.h"
#include <iomanip.h>
#include "numerical.h"

main () {
	double x1 = 4.034, x3= 4.213, x5=4.135;
	double x2= 4.4,  x4 = 4.31, x6=4.227;
	cout.precision(16);
	cout << D2Quoin(x1,x3,x5) << endl;
	double e,De,DDe;
	DDeta(x1,x3,x5,e,De,DDe);
	/*
	cout << "e="<<e << endl;
	cout << "De="<<De << endl;
	cout << "DDe="<<DDe << endl;
	cout << "Da=" << DaQuoin(1.02,1.15,1.255) << endl;
	cout << "DDa=" << DaDaQuoin(1.02,1.15,1.255) << endl;
	cout << "Db=" << DbQuoin(1.02,1.15,1.255) << endl;
	cout << "DaDb="<< DaDbQuoin(1.02,1.15,1.255) << endl;
	cout << "DbDb=" <<DbDbQuoin(1.02,1.15,1.255) << endl;
	double d,DD;
	Ddih(x1,x2,x3,x4,x5,x6,d,DD);
	cout << "d=" << d << endl;
	cout << "DD=" << DD << endl<<endl;
	cout << "dd2="<<sqrt32DDdih2(x1,x2,x3,x4,x5,x6) << endl;
	cout << "dd3="<<sqrt32DDdih3(x1,x2,x3,x4,x5,x6) << endl;
	cout << "B=" << BB(6.32) << endl;
	cout << "dB=" << dB(6.32) << endl;
	cout << "ddB=" << ddB(6.32) << endl;
	cout << "q315=" << D2Quoin315(x1,x3,x5) << endl;
	*/ //OK
	cout << "D1Vee1=" << D1Vee1(sqrt(x1),sqrt(x2),sqrt(x3),sqrt(x4),sqrt(x5),
					sqrt(x6))<<endl;
	}
