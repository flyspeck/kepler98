
// The file was created to compute the second derivative wrt x1
// of vorVc and tauVc.  These are D2Vee0 and -D2Vee1, respectively,
// up to a power of delta.

/* second derivatives */
// D[vorVc,{x1,2}] Sqrt[delta]^3  // x1 deriv. but input yi.
double D2Vee0(double y1,double y2,double y3,double y4,double y5,double y6);
// D[-tauVc,{x1,2}] Sqrt[delta]^3  // x1 deriv. but input yi.
double D2Vee1(double y1,double y2,double y3,double y4,double y5,double y6);


/* and the first derivatives */
// D[vorVc,x1] Sqrt[delta] // x1 deriv. but input yi.
double D1Vee0(double y1,double y2,double y3,double y4,double y5,double y6);
// D[-tauVc,x1] Sqrt[delta] // x1 deriv. but input yi.
double D1Vee1(double y1,double y2,double y3,double y4,double y5,double y6);



// These are intermediate functions used to construct D2Vee0 and D2Vee1.

double DaQuoin(double a,double b,double c);
double DaDaQuoin(double a,double b,double c);
double DbQuoin(double a,double b,double c);
double DaDbQuoin(double a,double b,double c);
double DbDbQuoin(double a,double b,double c);
void DDeta(double x1,double x3,double x5,
    double& eta,double& Deta,double& DDeta);
double D2Quoin(double x1,double x3,double x5);

// set sqrt(delta)D[dih,x1] and sqrt(delta)^3 D[dih,{x1,2}]:
void Ddih(double x1,double x2,double x3,double x4,double x5,double x6,
    double& sqrtdDd,double& sqrtd32DDd);


double sqrt32DDdih3(double x1,double x2,double x3,double x4,double x5,double x6);
double sqrt32DDdih2(double x1,double x2,double x3,double x4,double x5,double x6);

double ddB(double x);
double dB(double x);  // D[B[Sqrt[x]],x]
double BB(double x); // x = y^2. // B[Sqrt[x]]
double D2Quoin315(double x1,double x3,double x5);

