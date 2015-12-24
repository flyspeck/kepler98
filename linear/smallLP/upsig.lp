\\ anchored simplex : bounds on the score
\\ sig and d are for upright quarters.
\\ sig2 and d2 are for y4 in [2.51,2sq],
\\ sig3 and d3 are for y4 in [2sq,3.2].
\\ sig4 is a simplex with angle > 2.46

MAXIMIZE


     2 sig + 0 sig2 + 2 sig3 + 0 sig4
ST
sum: 2 d +  0 d2 + 2 d3  +  0 d4 < 6.283185307179586

special = 2 \\ special = 1 if there is a special, special =0 if no special
sig4 =0
d4 > 2.46
\\ d2>2.2
\\ d3>2.46

ptdef: pt = 0.05537364566846414

s1pt: sig -   0.05537364566846414 sigpt = 0
s2pt: sig2 -  0.05537364566846414 sig2pt = 0
s3pt: sig3 -  0.05537364566846414 sig3pt = 0

\\ IV.A2
sig1: sig -4.10113 d < -4.3223 
sig2: sig -0.80449 d < -0.9871
sig3: sig -0.70186 d < -0.8756
sig4: sig -0.24573 d < -0.3404 \\ changed from -0.3429 on 3/28/98.
\\ sig4old: sig -0.24573 d < -0.3429 
sig5: sig -0.00154 d < -0.0024
sig6: sig +0.07611 d < 0.1196

\\ IV.A4
sig21: sig2 -2.28501 d2< -3.421
sig22: sig2 -1.67382 d2< -2.616
sig23: sig2 -0.8285  d2< -1.4486
sig24: sig2 -0.390925 d2< -0.79
sig25: sig2 -0.12012 d2< -0.3088
sig26: sig2 -0.05007 d2 < -0.1558

\\ IV.A6
sig31: sig3 -2.28501 d3< -3.58
sig32: sig3 -1.67382 d3< -2.715
sig33: sig3 -0.8285  d3< -1.517
sig34: sig3 -0.390925 d3< -0.858
sig35: sig3 -0.12012 d3 - 0.009 special < -0.358
sig36: sig3 -0.05007 d3 - 0.009 special < -0.186

\\ IV.A22 gives special correction.

BOUNDS

sig free
sig2 free
sigpt free
sig2pt free
sig3 free
sig3pt free
0< d < 6.7
0< d2 < 6.7
0< d3 < 6.7

END
