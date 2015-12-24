
Angle[y1_, y2_, y6_] := ArcCos[(y1^2 + y2^2 - y6^2)/(2*y1*y2)]
APart[y__] := y[[{1, 2, 3, 4, 5, 6}]]
BPart[y__] := y[[1 + {0, 2, 7, 6, 8, 4}]]
CPart[y__] := y[[1 + {0, 7, 11, 9, 10, 8}]]
DPart[y__] := y[[1 + {0, 11, 1, 12, 5, 10}]]
Rogify[y1_, y2_, y6_] := {y1/2, eta[y1, y2, y6], TRUNC}
B[y_] := phi[TRUNC, TRUNC] + A[y/2]
A[h_] := (1 - h/TRUNC)*(phi[h, TRUNC] - phi[TRUNC, TRUNC])

dihanc[y1_,y2_,y6_]:= (* The dihedral angle of an anchor *)
	Module[{h1,h2,b,c}, (* Equation F.4.5 *)
        h1 = y1/2;
        h2 = y2/2;
        b = eta[y1,y2,y6];
        c = eta0[h1];
        If[N[b]>N[c],Return[0]];
        dihR[h1,b,c]];
kappaCycle[y1_,y2_,y6_]:= (* 1/2-kappa for an angle pi *)
	anc[y1,y2,y6] + (Pi/2 - dihanc[y1,y2,y6])crown[y1/2]/(2Pi);
sq:= Sqrt[2];

CorTauRoundLong[y1_, long_] := 
  Module[{Cospsi, h, Costhet}, h = y1/2; 
    Cospsi = (y1^2 + TRUNC^2 - long^2)/(2*y1*TRUNC); Costhet = h/TRUNC; 
    N[((1 - Cospsi)*tildemax + (1 - Costhet)*(tilde[h] - tildemax))*2*Pi]]

tildemax := zeta*pt - phi[TRUNC, TRUNC]
tilde[h_] := zeta*pt - phi[h, TRUNC]

overlap[y1_, y3_, y5_] := Delta[y1, y3, TRUNC, 1.6, 1.6, y5]

CrossDiag[y0_, y1_, y2_, y3_, y4_, y5_, y6_, y7_, y8_] := 
  Module[{x0, x1, x2, x3, x4, x5, x6, x7, x8, dd, c1, c0}, 
   {x0, x1, x2, x3, x4, x5, x6, x7, x8} = 
     {y0, y1, y2, y3, y4, y5, y6, y7, y8}^2; 
    dd = Dihedral[y3, y1, y5, y0, y4, y2] + Dihedral[y3, y1, y8, y6, y7, y2]; 
    c1 = -(x3*x3) - x5*x8 - x4*x7 + x3*(x5 + x8 + x4 + x7) + x5*x7 + x4*x8; 
    If[N[dd] > N[Pi], c0 = 
       Cos[2*Pi - dd]*Sqrt[U[x3, x5, x4]]*Sqrt[U[x3, x8, x7]]; 
      Return[Sqrt[(c0 - c1)/(-2*x3)]]]; 
    c0 = Cos[dd]*Sqrt[U[x3, x5, x4]]*Sqrt[U[x3, x8, x7]]; 
    Sqrt[(c0 - c1)/(-2*x3)]]

CrossTerm[y1_, y2_, lambda_] := 
  Module[{long, alpha1, alpha2, sol, t0, phi0, phi1, phi2}, 
   t0 = 1.254999999999999; long = 3.2; 
    alpha1 = Dihedral[y1, t0, y2, lambda, long, lambda]; 
    alpha2 = Dihedral[y2, t0, y1, lambda, long, lambda]; 
    sol = Solid[y2, t0, y1, lambda, long, lambda]; phi0 = phi[t0, t0]; 
    phi1 = phi[y1/2, t0]; phi2 = phi[y2/2, t0]; 
    2*(sol*(zeta*pt - phi0) + alpha1*(1 - y1/(2*t0))*(-phi1 + phi0) + 
        alpha2*(1 - y2/(2*t0))*(-phi2 + phi0)) + 
     ((Pi - 2*alpha2)*CorTauRoundLong[y2, lambda])/(2*Pi) + 
     ((Pi - 2*alpha1)*CorTauRoundLong[y1, lambda])/(2*Pi)]

Coss[x_, y_, z_] := (x^2 + y^2 - z^2)/(2*x*y)
Ac[x_, y_, z_] := ArcCos[Coss[x, y, z]]
arc[x1_, x2_, x3_] := ArcCos[(x1^2 + x2^2 - x3^2)/(2*x1*x2)];
dodec := 2*Sqrt[3]*Tan[Pi/5]
SeanVol[x__] := (VorVc[x, dodec]/4 - Solid[x]/3)/-doct
SeanVol[x__] := (NewTrunc[VorVc, x, dodec/2]/4 - Solid[x]/3)/-doct
Vol[x__] := (VorAnalytic[x]/4 - Solid[x]/3)/-doct
xiV = 0.003520999999999999
xiG = 0.01561
xiGG = 0.009350000000000001
Alambda[h_] := (A[h] - A[1])/(h - 1)
VorNu[x__] := (VorAnalytic[x] + VorAnalytic @@ Hat[x])/2 + VorVc[x]/2 - 
   VorVc @@ Hat[x]/2
Hat[x1_, x2_, x3_, x4_, x5_, x6_] := {x1, x6, x5, x4, x3, x2}
GammaNu[x__] := Gamma[x] + VorVc[x]/2 - VorVc @@ Hat[x]/2
Octavor[y1_, y2_, y3_, y4_, y5_, y6_] := 
  (VorAnalytic[y1, y2, y3, y4, y5, y6] + VorAnalytic[y1, y6, y5, y4, y3, y2])/
   2
OctavorVc[y1_, y2_, y3_, y4_, y5_, y6_] := 
  (VorVc[y1, y2, y3, y4, y5, y6] + VorVc[y1, y6, y5, y4, y3, y2])/2
xsub = {x[0] -> y1, x[1] -> y2, x[2] -> y3, x[3] -> y4, x[4] -> y5, 
   x[5] -> y6}
tauNu[x__] := Solid[x]*zeta*pt - GammaNu[x]
tauAnalyticNu[x__] := Solid[x]*zeta*pt - VorNu[x]
xikG = -0.01339
xik = -0.029
pimax = 0.06687999999999999
target := (4*Pi*zeta - 8)*pt

delta = -(x2*x3*x4) - x1*x3*x5 - x1*x2*x6 - x4*x5*x6 + 
   x3*(x1 + x2 - x3 + x4 + x5 - x6)*x6 + 
   x2*x5*(x1 - x2 + x3 + x4 - x5 + x6) + x1*x4*(-x1 + x2 + x3 - x4 + x5 + x6)
sqx = {Sqrt[x1], Sqrt[x2], Sqrt[x3], Sqrt[x4], Sqrt[x5], Sqrt[x6]}
MakeRandom := Random[Integer, {0, 10^9}]
MakeRandom := (SeedRandom; MakeRandom := Random[Integer, {0, 10^9}]; 
    MakeRandom)
DodecVol = (10*Sqrt[2]*(-5 + Sqrt[5])^2)/(5 + Sqrt[5])^(3/2)
Vc[i1_, i2_, i3_, i4_, i5_] := 
  CorTau[m[i1], m[i2], m[i3], 3.2, 2, 2] + 
   CorTau[m[i2], m[i3], m[i4], 3.2, 2, 2] + 
   CorTau[m[i3], m[i4], m[i5], 3.2, 2, 2] + 
   CorTau[m[i4], m[i5], m[i1], 3.2, 2, 2.509999999999998] + 
   CorTau[m[i5], m[i1], m[i2], 3.2, 2, 2.509999999999998]
Vc[i1_, i2_, i3_, i4_, i5_] := 
  CorTau[m[i1], m[i2], m[i3], 2, 3.2, 2] + 
   CorTau[m[i2], m[i3], m[i4], 2, 3.2, 2] + 
   CorTau[m[i3], m[i4], m[i5], 2, 3.2, 2] + 
   CorTau[m[i4], m[i5], m[i1], 2.509999999999999, 3.2, 2] + 
   CorTau[m[i5], m[i1], m[i2], 2, 3.2, 2.509999999999999]
Vc[i1_, i2_, i3_, i4_, i5_] := 
  CorTau[m[i2], m[i1], m[i3], 3.2, 2, 2] + 
   CorTau[m[i3], m[i2], m[i4], 3.2, 2, 2] + 
   CorTau[m[i4], m[i3], m[i5], 3.2, 2, 2] + 
   CorTau[m[i5], m[i4], m[i1], 3.2, 2.509999999999998, 2] + 
   CorTau[m[i1], m[i5], m[i2], 3.2, 2, 2.509999999999998]
fax[x__, fn_] := fn[x] + 0.4193509999999999*Solid[x]
(* T := Sqrt[3]*Tan[Pi/5] *)
muS[x__] := Vol[x] - 0.42755*Solid[x]
SeanM = 0.42755
SeanSquanderTarget := DodecVol - 4*Pi*SeanM
SeanSquander[x__] := SeanVol[x] - SeanM*Solid[x]
SeanT := Sqrt[3]*Tan[Pi/5]
dist[x_, y_] := Sqrt[(x - y) . (x - y)]
SeanDodecVol := DodecVol
NewGamma[x__] := Gamma[x] + (nr*ArcTanx @@ ({x} - {2, 2, 2, 2, 2, 2}))/nv
ArcTanx[x__] := {1, 1, 1, -1, -1, -1} . ArcTan /@ {x}
(* nr := 0.004 nv := 0.25 *)
NewTau[x__] := Solid[x]*zeta*pt - NewGamma[x]
NewTau[x__] := Solid[x]*zeta*pt - NewGamma[x]
nr := 0.005
nv := 0.1
NewGamma[x__] := Gamma[x] + nr*ArcTanx @@ (({x} - {2, 2, 2, 2, 2, 2})/nv)
tt[y1_, y2_, y3_, y4_, y5_, y6_] := 
  ((y1 - 2)*Dihedral[y1, y2, y3, y4, y5, y6])/(2*Pi) + 
   ((y2 - 2)*Dihedral[y2, y1, y3, y5, y4, y6])/(2*Pi) + 
   (Dihedral[y3, y2, y1, y6, y5, y4]*(y3 - 2))/(2*Pi) + 
   ttc*(Solid[y1, y2, y3, y4, y5, y6] - Solid[2, 2, 2, 2, 2, 2])
ttc := 0.5





Bisect[{v1_, v2_}, eps_] := 
  Module[{y1, y2, ym, x1, x2, xm, reps, f}, 
   {x1, x2} = Sort[N[{v1, v2}]]; f = Bisectf; xm := (x1 + x2)/2; y1 = f[x1]; 
    y2 = f[x2]; ym := f[xm]; If[y1*y2 > 0, Return[]]; 
    reps = Ceiling[Log[(x2 - x1)/eps]/Log[2]]; 
    Do[ym = f[xm]; If[ym*y1 < 0, x2 = xm, x1 = xm]; , {i, 1, reps}]; {x1, x2}]
 
Bisect[{v1_, v2_}, eps_, f_] := 
  Module[{y1, y2, ym, x1, x2, xm, reps}, 
   {x1, x2} = Sort[N[{v1, v2}]]; xm := (x1 + x2)/2; y1 = f[x1]; y2 = f[x2]; 
    ym := f[xm]; If[y1*y2 > 0, Return[]]; 
    reps = Ceiling[Log[(x2 - x1)/eps]/Log[2]]; 
    Do[ym = f[xm]; If[ym*y1 < 0, x2 = xm, x1 = xm]; , {i, 1, reps}]; {x1, x2}]



GammaY[y1_, y2_, y3_, y4_, y5_, y6_] := 
  Gamma[y1, y2, y3, y4, y5, y6] + Y[y1] + Y[y2] + Y[y3] - Y[y4] - Y[y5] - 
   Y[y6]
Y[y_] := 0.006*ArcTan[10*(y - 2)]

NewDihedral[y1_,y2_,y3_,y4_,y5_,y6_]:= Module[
	{d4,x1=y1*y1,x2=y2*y2,x3=y3*y3,x4=y4*y4,x5=y5*y5,x6=y6*y6},
	d4 = -x1^2 + x1*x2 + x1*x3 - x2*x3 - 2*x1*x4 + x1*x5 + x2*x5 + x1*x6 + x3*x6 - x5*x6;
	Pi/2 + ArcTan[-d4/Sqrt[4 x1 Delta[y1,y2,y3,y4,y5,y6]]]
	];
A[h_] := (1 - h/TRUNC)*(phi[h, TRUNC] - phi[TRUNC, TRUNC])
 
A[h_, t_] := (1 - h/t)*(phi[h, t] - phi[t, t])
phiV[h_, t_] := (h*t*(h + t))/6
AV[h_, t_] := (1 - h/t)*(phi[h, t] - phi[t, t])
AV[h_, t_] := (1 - h/t)*(phi[t, t] - phi[h, t])
AV[h_, t_] := (1 - h/t)*(phiV[h, t] - phiV[t, t])
red[f_, x_, i_] := 
  Module[{j, k}, k = Series[f, {x, 0, i}]; 
    k[[3]] . Table[x^j, {j, k[[4]], k[[4]] + Length[k[[3]]] - 1}]]
degree[f_, x_] := -Series[f /. x -> 1/x, {x, 0, 1}][[4]]
degree[f_, x_] := -Series[f /. x -> 1/x, {x, 0, 1}][[4]]
RED[f_, x_] := Complement[Table[red[f, x, i], {i, 0, degree[f, x]}], {0}]
RED[f_, x_] := Complement[Table[red[f, x, i], {i, 0, degree[f, x]}], {0}]
red[f_, x_, i_] := 
  Module[{j, k}, k = CoefficientList[f, x]; k = Drop[k, -i]; 
    k . Table[x^j, {j, 0, Length[k] - 1}]]
dtet := Sqrt[8]*ArcTan[Sqrt[2]/5]
taut[x__] := Volt[x] - Solid[x]/(3*dtet)
phit := phi[trunc, trunc] - 1/(3*dtet)
arc[a_, b_, c_] := ArcCos[(a^2 + b^2 - c^2)/(2*a*b)]
taut[x__] := Volt[x] - Solid[x]/(3*dtet)
eSolid[x__] := Solid[x] - Solid @@ base
eVolt[x__] := Volt[x] - Volt @@ base
sec[x_, a_, b_, f_] := 
  (f /. x -> a) + ((x - a)*((f /. x -> b) - (f /. x -> a)))/(b - a)
tan[x_, a_, f_] := (f /. x -> a) + (D[f, x] /. x -> a)*(x - a)
L[x_] := 0.035*(x - 2) - 0.06*(x - 2)*(x - edgeMax)
fcc[x__] := VolAn[x] - VolAn[2, 2, 2, 2, 2, 2] - 
   AA*(Solid[x] - Solid[2, 2, 2, 2, 2, 2]) + 
   ({x} - {2, 2, 2, 2, 2, 2}) . {-aa, -aa, -aa, bb, bb, bb}
quad[x__] := VolAn[x] - VolAn[2, 2, 2, Sqrt[8], 2, 2] - 
   AA*(Solid[x] - Solid[2, 2, 2, Sqrt[8], 2, 2]) + 
   ({x} - {2, 2, 2, Sqrt[8], 2, 2}) . 
    {-aa, -(aa/2), -(aa/2), -0.003, 0.06061999999999999, 0.06061999999999999}
As[y1_, y2_, y3_, y4_, y5_, y6_] := 
  Module[{d1, d2, di, t1, t2}, d1 = dihR[y1/2, eta[y1, y2, y5], eM/2]; 
    d2 = dihR[y1/2, eta[y1, y3, y6], eM/2]; 
    t1 = tauR[y1/2, eta[y1, y2, y5], eM/2]; 
    t2 = tauR[y1/2, eta[y1, y3, y6], eM/2]; 
    di = Dihedral[y1, y2, y3, y4, y5, y6]; 
    If[N[eta[y1, y2, y5]] > eM/2, d1 = 0; t1 = 0]; 
    If[N[eta[y1, y3, y5]] > eM/2, d2 = 0; t2 = 0]; di = di - d1 - d2; 
    t1 + t2 + di*(1 - Cos[arc[y1, eM/2, eM/2]])*phit - di*A[y1/2] - L[y1]]
Bs[y1_, y2_, y3_, y4_, y5_, y6_] := 
  Module[{d1, d2, di, t1, t2}, d1 = dihR[y1/2, eta[y1, y2, y5], eM/2]; 
    t1 = tauR[y1/2, eta[y1, y2, y5], eM/2]; 
    di = Dihedral[y1, y2, y3, y4, y5, y6]; 
    If[N[eta[y1, y2, y5]] > eM/2, d1 = 0; t1 = 0]; di = di - d1; 
    t1 + di*(1 - Cos[arc[y1, eM/2, cc]])*phit - di*A[y1/2] - L[y1]/2]
Bs[y1_, y2_, y3_, y4_, y5_, y6_] := 
  Module[{d1, d2, di, t1, t2}, d1 = dihR[y1/2, eta[y1, y2, y5], eM/2]; 
    t1 = tauR[y1/2, eta[y1, y2, y5], eM/2]; 
    di = Dihedral[y1, y2, y3, y4, y5, y6]; 
    If[N[eta[y1, y2, y5]] > eM/2, d1 = 0; t1 = 0]; di = di - d1; 
    t1 + di*(1 - Cos[arc[y1, eM/2, cc]])*phit - di*A[y1/2] - L[y1]/2]
jacn := 3
jacA := (jacn - 3)/2
J[n_, x_] := JacobiP[n, jacA, jacA, x]
aPoly[n_] := Sum[J[i, x]*a[i], {i, 0, n}]
eqn[p_, n_] := Table[(D[p, {x, i}] /. x -> 0) == 0, {i, 0, n}]
solveJ[p_, n_] := Solve[eqn[p, n], Table[a[i], {i, 0, n}]]
Ac = {0, 2/Sqrt[3]}
Bc = {1, -(1/Sqrt[3])}
