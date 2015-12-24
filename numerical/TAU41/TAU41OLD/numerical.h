//  copyright (c) 1997, Thomas C. Hales, all rights reserved.

#ifndef ineq_c
#define ineq_c

// These routines use ordinary floating point arithmetic without any
// reference to IEEE-754 standards, rounding modes, interval methods, etc.
// They cannot be used in the rigorous verification of inequalities.
// They can be used to propose and test inequalities that are to be
// proved later by rigorous means.

// simplex:
double delta(double x1,double x2,double x3,double x4,double x5,double x6);
double solid(double y1,double y2,double y3,double y4,double y5,double y6);
double dihedraly(double y1,double y2,double y3,double y4,
                        double y5, double y6);

// circumradius: 
double eta2(double x1,double x2,double x3);
double radf(double y1,double y2,double y3);
double chi(double x1,double x2,double x3,double x4,double x5,double x6);
double circum2(double x1,double x2,double x3,double x4,double x5,double x6);
double rady(double y1,double y2,double y3,double y4,double y5,double y6);

// scoring:
extern const double doct;
double gamma(double y1,double y2,double y3,double y4,double y5, double y6);
double vor_analytic(double y1,double y2,double y3,double y4,
        double y5,double y6);
double tau_analytic(double y1,double y2,double y3,double y4,
        double y5,double y6);
double octavor(double y1,double y2,double y3,double y4,
        double y5,double y6); // average vor_analytic on an upright simplex 
double octavorVc(double y1,double y2,double y3,double y4,
        double y5,double y6); // average VorVc on an upright
double vorVc(double y1,double y2,double y3,double y4,
        double y5,double y6);
double vorVc(double y1,double y2,double y3,double y4,double y5,
        double y6,double trunc);
double vorSqc(double y1,double y2,double y3,double y4,
        double y5,double y6);
double tauVc(double y1,double y2,double y3,double y4,
        double y5,double y6);
double tauSqc(double y1,double y2,double y3,double y4,
        double y5,double y6);
double tau(double y1,double y2,double y3,double y4,
        double y5,double y6);
double dihVc(double x1,double x2,double x3,double x4,double x5,double x6);
double dihSqc(double x1,double x2,double x3,double x4,double x5,double x6);
double skelV(double x1,double x2,double x3,double x4,double x5,double x6);
	
double corsigma(double y1,double y2,double y3,double y4,double y5,double y6);
double corsolid(double y1,double y2,double y3,double y4,double y5,double y6);
double cortau(double y1,double y2,double y3,double y4,double y5,double y6);
double sol_cor(double y1,double y2,double y3);



//
double tilde(double h);
double pretilde(double h,double t);
double crown(double h,double t);
double crown(double h);
double sfacelift(double y1,double y2,double y6,double t);


// rogers:
double dihR(double y1,double y2,double y3);
double solR(double y1,double y2,double y3);
double vorR(double y1,double y2,double y3);
double tauR(double y1,double y2,double y3);

// misc:
double myrand();
double mabs(double a);
double safesqrt(double x);

// full quad cluster:
// we use mod three numbering:
//  	indices on back simplex are standard order mod 3.
//
//y[0],y[1],y[2],y[3],y[4],y[5] = 1st simplex. y[3]=diag.
//y[6],y[1],y[2],y[3],y[7],y[8] = 2nd simplex.
//
//opposite pairs: (y[4],y[8]),(y[5],y[7]),(y[0],y[6]),
//			(y[1],y[2])
//diagonals: y[3], other is between (0,6), computed by crossdiag(y);
// corners are 0,1,2,6.
//
double vorVc9(double y[9]);
double solid9(double y[9]);
double dih9_0(double y[9]);
double dih9_1(double y[9]);
double dih9_2(double y[9]);
double dih9_6(double y[9]);
double crossdiag(double y[9]);


#endif
