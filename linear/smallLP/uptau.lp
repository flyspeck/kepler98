\\ upright quarters.

\\ d = upright quarter, d2= [2.51,2sq2], d3 = [2sqrt2,3.2], d4 > 2.46

MINIMIZE 

2 tauQuarter + 1 tau2+ 0 tau3 + 1 tau4 - 0 pen

ST

2 d + 1 d2 + 0 d3 + 1 d4 < 6.28318530717958
special = 0 \\ special=0 no special simplices, special=1 specials.
pen =  0.06585
tau4=0
 d4>2.46
\\d4 > 1.65

tauQuarter - 0.05537364566846414 tauQuarterpt = 0
tau2 -  0.05537364566846414 tau2pt = 0
tau3 -  0.05537364566846414 tau3pt = 0

\\ IV.A3... quarters with "nu" scoring.
-tauQuarter -4.16523 d < -4.42873
-tauQuarter -0.78701 d < -1.01104
-tauQuarter -0.77627 d < -0.99937
-tauQuarter -0.21916 d < -0.34877
-tauQuarter -0.05107 d < -0.11434
-tauQuarter +0.07106 d <  0.07749

\\ IV.A5... 2.51<y4<2sq, tau_x scoring.
-tau2 -2.1747 d2 < -3.3407
-tau2 -1.87427 d2 < -2.945
-tau2 -0.83046 d2 < -1.5035
-tau2 -0.48263 d2 < -1.0009
-tau2 -0.34833 d2 < -0.7787  \\ was -0.7829 until 3/22/98
-tau2 -0.1694 d2 < -0.4475
-tau2 -0.0822 d2 < -0.2568

\\ IV.A7... y4>2sq, tau0 scoring.
-tau3 -2.1747 d3 < -3.48
-tau3 -1.87427 d3 < -3.06
-tau3 -0.83046 d3 < -1.58
-tau3 -0.48263 d3 < -1.06
-tau3 -0.34833 d3 < -0.83
-tau3 -0.1694 d3 < -0.50
-tau3 -0.0822 d3 - 0.0014 special < -0.29
\\ patch: -tau3 -0.0822 d3  < -0.318 \\ for th (4,1) case when d3 > 2.46
BOUNDS

tauQuarter free
tau2 free
tau3 free
tauQuarterpt free
tau2pt free
tau3pt free
0< d < 6.7
0< d2 < 6.7
0< d3 < 6.7

END
