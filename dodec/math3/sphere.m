(* 
<html>
<br><br>
<body background="paperbackground.JPG" alink=blue vlink=red>
	Mathematica package for sphere packings
<PRE>
*)


SPhi[h_]:=1/6 T (T+h) h


dode=5.550291041
M=.42755
target=dode-4 Pi M//N

S={2,2,2,2,2,2}

L=1.94159

Squander[x1_,x2_,x3_,x4_,x5_,x6_]:=Volume[x1,x2,x3,x4,x5,x6]-M Solid[x1,x2,x3,x4,x5,x6]//N

Squander2[x1_,x2_,x3_,x4_,x5_,x6_]:=TruncVol[x1,x2,x3,x4,x5,x6]-M Solid[x1,x2,x3,x4,x5,x6]//N



VarSquander[y1_,y2_,x1_,x2_,x3_,x4_,x5_,x6_]:=Volume[x1,x2,x3,x4,x5,x6]-y1 Solid[x1,x2,x3,x4,x5,x6]- y2 Dihedral[x1,x2,x3,x4,x5,x6]//N

VarSquander2[y1_,y2_,x1_,x2_,x3_,x4_,x5_,x6_]:=TruncVol[x1,x2,x3,x4,x5,x6]-y1 Solid[x1,x2,x3,x4,x5,x6]- y2 Dihedral[x1,x2,x3,x4,x5,x6]//N



Keep[x_]:= PutAppend[Definition[x],"more.m"];
SetAttributes[Keep,HoldFirst];
 
Fix[x_] := Edit[Definition[x]]
SetAttributes[Fix,HoldFirst];
(Unprotect[In,Out]; Format[In]=SphereIn; 
	Format[Out]=SphereOut; Protect[In,Out];)

$SpherePrecision:= 16;
Ns[x_]:= N[x,$SpherePrecision];
zeta := 1/(2*ArcTan[2^(1/2)/5])//Ns;
pt := (-Pi/3 + 4*ArcTan[2^(1/2)/5])//Ns;
doct := ((Pi - 4*ArcTan[2^(1/2)/5])/(2*2^(1/2)))//Ns;

chi = -(x1*x4^2) + x1*x4*x5 + x2*x4*x5 - x2*x5^2 + x1*x4*x6 + x3*x4*x6 +
   x2*x5*x6 + x3*x5*x6 - 2*x4*x5*x6 - x3*x6^2

eta[x_,y_,z_]:= Module[{s= (x+y+z)/2}, x y z/4/Sqrt[s (s-x)(s-y)(s-z)]]//Ns;

volR[a_, b_, c_] := Sqrt[a^2*(b^2 - a^2)*(c^2 - b^2)]/6//Ns;
solR[x_, y_, z_] := 2*ArcTan[Sqrt[((z - y)*(y - x))/((z + y)*(x + y))]]//Ns;
dihR[x_, y_, z_] := ArcTan[Sqrt[(z^2 - y^2)/(y^2 - x^2)]]//Ns;
vorR[a_, b_, c_] := 4*(-(doct*volR[a, b, c]) + solR[a, b, c]/3)//Ns;
denR[a_, b_, c_] := solR[a, b, c]/(3*volR[a, b, c])
tauR[x__] := (-vorR[x] + solR[x]*zeta*pt)//Ns;
ttauR[x__] := (vorR[x] - 3.2*solR[x]*zeta*pt)//Ns;
  
T=N[Sqrt[3]Tan[Pi/5],20];
TRUNC= T;

Dihedral2[y1_,y2_,y3_,y4_,y5_,y6_] :=   (* Dihedral along edge 2 *)
  Dihedral[y2,y1,y3,y5,y4,y6];
 
 
Dihedral3[y1_,y2_,y3_,y4_,y5_,y6_] :=   (* Dihedral along edge 3 *)
	Dihedral[y3,y1,y2,y6,y4,y5];


Chi[x__] := Module[{x1, x2, x3, x4, x5, x6}, 
   {x1, x2, x3, x4, x5, x6} = {x}^2; 
    -(x1*x4^2) + x1*x4*x5 + x2*x4*x5 - x2*x5^2 + x1*x4*x6 + x3*x4*x6 + 
     x2*x5*x6 + x3*x5*x6 - 2*x4*x5*x6 - x3*x6^2];


tauAnalytic[x__] := Solid[x]*zeta*pt - VorAnalytic[x]
ttauAnalytic[x__]:= VorAnalytic[x] - 3.2 Solid[x] zeta pt;
tau[x__] := Solid[x]*zeta pt - Gamma[x];
ttau[x__] := Gamma[x] - 3.2*Solid[x]*zeta pt;


VorAnalytic[x__]:= Module[{y1,y2,y3,y4,y5,y6,x1,x2,x3,x4,x5,x6,
    del,u126,u135,u234,vol},
    {y1,y2,y3,y4,y5,y6} = {x};
    {x1,x2,x3,x4,x5,x6} = {x}^2;
    del = Sqrt[Delta[x]];
    u126 = U[x1,x2,x6];
    u135 = U[x1,x3,x5];
    u234 = U[x2,x3,x4];
    vol=1/(48 del) (
     (x1 (x2+x6-x1)+ x2 (x1+x6-x2)) Chi[y4,y5,y3,y1,y2,y6]/u126 +
     (x2 (x3+x4-x2) + x3 (-x3+x4+x2)) Chi[y6,y5,y1,y3,y2,y4]/u234 +
     (x1 (-x1+x3+x5) + x3 (x1-x3+x5)) Chi[y4,y6,y2,y1,y3,y5]/u135 );
     4(-doct vol + Solid[x]/3)
    ];


Volume[x__]:= Module[{y1,y2,y3,y4,y5,y6,x1,x2,x3,x4,x5,x6,
    del,u126,u135,u234,vol},
    {y1,y2,y3,y4,y5,y6} = {x};
    {x1,x2,x3,x4,x5,x6} = {x}^2;
    del = Sqrt[Delta[x]];
    u126 = U[x1,x2,x6];
    u135 = U[x1,x3,x5];
    u234 = U[x2,x3,x4];
    vol=1/(48 del) (
     (x1 (x2+x6-x1)+ x2 (x1+x6-x2)) Chi[y4,y5,y3,y1,y2,y6]/u126 +
     (x2 (x3+x4-x2) + x3 (-x3+x4+x2)) Chi[y6,y5,y1,y3,y2,y4]/u234 +
     (x1 (-x1+x3+x5) + x3 (x1-x3+x5)) Chi[y4,y6,y2,y1,y3,y5]/u135 );
     N[vol,20]
    ];

Delta[y1_,y2_,y3_,y4_,y5_,y6_]:= Module[{x1,x2,x3,x4,x5,x6},
	x1= y1*y1; x2=y2*y2; x3=y3*y3;
	x4= y4*y4; x5=y5*y5; x6=y6*y6;
	(x1*x4*(-x1+x2+x3-x4+x5+x6)+
			x2*x5*(x1-x2+x3+x4-x5+x6)
			+x3*x6*(x1+x2-x3+x4+x5-x6)
			-x2*x3*x4-x1*x3*x5-x1*x2*x6-x4*x5*x6)];


HPhi[h_,t_]:=2/3*(2-doct * h * t * (h+t))
P0=HPhi[T,T]
A[h_]:=(1-h/T)*(HPhi[h,T]-HPhi[T,T])





U[x1_,x2_,x6_]:= -x1*x1-x2*x2-x6*x6+2*x1*x6+2*x1*x2+2*x2*x6;


Rho[y1_,y2_,y3_,y4_,y5_,y6_]:= Module[{x1,x2,x3,x4,x5,x6},
x1= y1*y1; x2=y2*y2; x3=y3*y3;
x4= y4*y4; x5=y5*y5; x6=y6*y6;
(-x1*x1*x4*x4-x2*x2*x5*x5-x3*x3*x6*x6+
        2*x1*x2*x4*x5+2*x1*x3*x4*x6+2*x2*x3*x5*x6)];

Rad[y1_,y2_,y3_,y4_,y5_,y6_]:= 
	Sqrt[Rho[y1,y2,y3,y4,y5,y6]/Delta[y1,y2,y3,y4,y5,y6]]/2;


aSolid[y1_,y2_,y3_,y4_,y5_,y6_]:=
	y1*y2*y3 + y1*(y2*y2+y3*y3-y4*y4)/2 +
					y2*(y1*y1+y3*y3-y5*y5)/2 +
					y3*(y1*y1+y2*y2-y6*y6)/2;


Unprotect[Gamma];

Attributes[Gamma]= {};
Gamma::usage="Gamma[x] is the compression of the simplex with edge lengths x";
Clear[Gamma];
Gamma[y1_,y2_,y3_,y4_,y5_,y6_]:= Module[{u,t,a1,a2,a3,a4},
	u=Delta[y1,y2,y3,y4,y5,y6];
	If[u<=0,Return[0]];
	t = Sqrt[u]/2;
	a1 = aSolid[y1,y2,y3,y4,y5,y6];
	a2 = aSolid[y1,y5,y6,y4,y2,y3];
	a3 = aSolid[y2,y4,y6,y5,y1,y3];
	a4 = aSolid[y4,y5,y3,y1,y2,y6];
	(-doct*t/6 + (2/3)*(ArcTan[t/a1]+ArcTan[t/a2]+ArcTan[t/a3]+ArcTan[t/a4]))
];
Protect[Gamma];


Solid[y1_,y2_,y3_,y4_,y5_,y6_]:= Module[{u,t,a1},
	u = Delta[y1,y2,y3,y4,y5,y6];
	If [N[u]<=0, Return[0]];
	t = Sqrt[u]/2;
	a1 = aSolid[y1,y2,y3,y4,y5,y6];
	2*ArcTan[t/a1]
]//Ns;


Dihedral[y1_,y2_,y3_,y4_,y5_,y6_]:= Module[{ux,x1,x2,x3,x4,x5,x6,t},
	x1= y1*y1; x2=y2*y2; x3=y3*y3;
	x4= y4*y4; x5=y5*y5; x6=y6*y6;
	ux = U[x1,x2,x6]*U[x1,x3,x5];
	t=(-2*x1*x4+x1*(x2+x3+x5+x6)+x2*x5+x3*x6-x1*x1 -x2*x3-x5*x6)/Sqrt[ux];
	ArcCos[t]
];

Density[x_]:= doct/(1-3 x/(16 Pi));

Norm[x_]:= x.x;
Norm[x_,y_]:= Norm[x-y];
Distance[x_,y_]:= Sqrt[(x-y).(x-y)];

FarFrom[x_,pair_]:= If[Norm[x,pair[[1]]]<Norm[x,pair[[2]]],
                        pair[[2]],pair[[1]]];

ExtremePoint[{x_,y_,z_}]:= LinearSolve[{x,y,z},{x.x/2,y.y/2,z.z/2}];

Vertex[{X_, dx_}, {Y_, dy_}, {Z_, dz_}] :=
  Simplify[Block[{w, x, y, z}, w = {x, y, z};
     w /. Solve[{Norm[X, w] == dx^2,
             Norm[Y, w] == dy^2, Norm[Z, w] == dz^2}, w]]];

SimplexCoordinates[y1_,y2_,y3_,y4_,y5_,y6_]:=
	Module[{x1,x2,x3,x4,x5,x6,X,Y,Z},
		x1 = y1^2; x2 = y2^2; x3 = y3^2;
		x4 = y4^2; x5 = y5^2; x6 = y6^2;
		 X = {x1^(1/2), 0, 0}; 
    Y = {(x1 + x2 - x6)/(2*x1^(1/2)), 
      (-x1^2 + 2*x1*x2 - x2^2 + 2*x1*x6 + 2*x2*x6 - x6^2)^(1/2)/
       (2*x1^(1/2)), 0}; Z = 
     {(x1 + x3 - x5)/(2*x1^(1/2)), 
      -((x1^2 - x1*x2 - x1*x3 + x2*x3 + 2*x1*x4 - x1*x5 - x2*x5 - x1*x6 - 
           x3*x6 + x5*x6)/
         (2*x1^(1/2)*(-x1^2 + 2*x1*x2 - x2^2 + 2*x1*x6 + 2*x2*x6 - x6^2)^
            (1/2))), (-(x2*x3*x4) - x1*x3*x5 - x1*x2*x6 - x4*x5*x6 + 
          x3*(x1 + x2 - x3 + x4 + x5 - x6)*x6 + 
          x2*x5*(x1 - x2 + x3 + x4 - x5 + x6) + 
          x1*x4*(-x1 + x2 + x3 - x4 + x5 + x6))^(1/2)/
   (-x1^2 + 2*x1*x2 - x2^2 + 2*x1*x6 + 2*x2*x6 - x6^2)^(1/2)}; {X, Y, Z}]//Ns


Enclosed[x__, z1_, z2_, z3_] :=
  Module[{X, Y, Z}, If[Length[{x}] != 6, Return["invalid input"]];
    {X, Y, Z} = SimplexCoordinates[x];
    Sqrt[Norm[FarFrom[NullVector, Vertex[{X, z1}, {Y, z2}, {Z, z3}]]]]];

NullVector={0,0,0};


Quoin[a_, b_, c_] :=
  Module[{u = Sqrt[(c^2 - b^2)/(b^2 - a^2)]},
	If[N[a]>N[b],Return[0]];
	If[N[b]>N[c],Return[0]];
   -(a^2 + a*c - 2*c^2)*(c - a)*ArcTan[u]/6 + a*(b^2 - a^2)*u/6 -
    (2/3)*c^3*ArcTan[(u*(b - a))/(b + c)]]//Ns;

Qn[y1_,y2_,z_]:= -4 doct (Quoin[y1/2,eta[y1,y2,z],T] +
		Quoin[y2/2,eta[y1,y2,z],T])//Ns;

(* VorVc is the truncation vor(S,1.255) of the Voronoi cell at t0=1.255 *)
VorVc[y1_,y2_,y3_,y4_,y5_,y6_]:= Module[{h1,h2,h3,a1,a2,a3},
        h1 = y1/2; h2 = y2/2; h3 = y3/2;
        a1 = Dihedral[y1,y2,y3,y4,y5,y6];
        a2 = Dihedral[y2,y3,y1,y5,y6,y4];
        a3 = Dihedral[y3,y1,y2,y6,y4,y5];
        vorstar[h1,h2,h3,a1,a2,a3]+
           Qn[y1,y2,y6]+Qn[y2,y3,y4]+Qn[y1,y3,y5]
        ];

TruncVol[y1_,y2_,y3_,y4_,y5_,y6_]:=N[(VorVc[y1,y2,y3,y4,y5,y6]-
	(4/3) Solid[y1,y2,y3,y4,y5,y6])/(-4 doct),20]

phi[h_,t_]:= 2 (2 - doct h t (h+t))/3//Ns;
phi0 := phi[TRUNC,TRUNC]//Ns;

(* two function from formulation *)
crown[h_]:= Module[{e}, e = eta0[h];
		2Pi (1-h/e) (phi[h,e] - phi0)
		];
eta0[h_]:= eta[ 2h,2 TRUNC,2];
anc[y1_,y2_,y6_]:= Module[{h1,h2,b,c}, (* Equation F.4.5 *)
		h1 = y1/2;
		h2 = y2/2;
		b = eta[y1,y2,y6];
		c = eta0[h1];
		If[N[b]>N[c],Return[0]];
		-dihR[h1,b,c] crown[h1]/(2Pi) -solR[h1,b,c] phi0 + vorR[h1,b,c] -
		dihR[h2,b,c](1-h2/TRUNC)(phi[h2,TRUNC]-phi0)-
		solR[h2,b,c] phi0 + vorR[h2,b,c]
		];
(* from SP IV, section 2.5. *)
kappa[y1_,y2_,y3_,y4_,y5_,y6_]:= crown[y1/2] Dihedral[y1,y2,y3,y4,y5,y6]/(2Pi)+
	anc[y1,y2,y6]+anc[y1,y3,y5];



vorstar[h1_,h2_,h3_,a1_,a2_,a3_]:=  Module[{Mx},
        Mx[h_]:= If[N[h]<N[TRUNC],h,TRUNC];
        a1 (1-Mx[h1]/TRUNC)(phi[h1,TRUNC]-phi0) +
        a2 (1-Mx[h2]/TRUNC)(phi[h2,TRUNC]-phi0) +
        a3 (1-Mx[h3]/TRUNC)(phi[h3,TRUNC]-phi0) +
        (a1+a2+a3-Pi)phi0//Ns
        ];




tauVc[x__]:= Solid[x] zeta pt - VorVc[x];
VorSqc[x__] := NewTrunc[VorVc, x, Sqrt[2]];

NewTrunc[f_,x__,t_]:= Module[{a,tbak}, tbak = TRUNC; TRUNC=t;
    a = f[x]; TRUNC=tbak; a]

(* checked this far *)

phiTau[h_]:= zeta pt - phi[h,TRUNC]//Ns;
phiTauMax:= phiTau[TRUNC];

beta[psi_,y1_,y3_,y5_]:= (* function IV.2.8 *)
	Module[{cb2,ctheta2},
		ctheta2 = ((y1^2+y3^2-y5^2)/(2 y1 y3))^2;
		cb2 = (Cos[psi]^2-ctheta2)/(1-ctheta2);
		ArcCos[Sqrt[cb2]]
		];

Arc[a_,b_,c_]:=ArcCos[((a*a+b*b-c*c)/(2*a*b))];


Psi2[y1_]:=ArcCos[(y1/(2*T))];

beta2[y1_,y3_,y5_]:=Module[{a,b,ppb,pb},
  a=Psi2[(y1)]; b=Arc[y1,y3,y5];
  ppb=(Cos[(a)]*Cos[(a)]-Cos[(b)]*Cos[(b)])/(1-Cos[(b)]*Cos[(b)]);
  pb =Sqrt[(ppb)];Let $w_1,w_2$ be the vertices adjacent to $v$.
  ArcCos[(pb)] ];





(* Calculate value of f at a, and all the partials of f at a *)
HyperInterval[f_,a_]:= Module[{f0,var,fs,i,j,sub},
		f0 = f @@ a;
		var = Table[x[i],{i,1,Length[a]}];
		fs = f @@ var;
		sub = Table[x[j]->a[[j]],{j,1,Length[a]}];
		{f0,Table[D[fs,x[i]]/.sub,{i,1,Length[var]}],
		Table[D[fs,x[i],x[j]]/.sub,{i,1,Length[var]},{j,1,Length[var]}]}
		];











(*
</PRE>
	copyright (c) 1997, Thomas C. Hales, all rights reserved.
</HTML>
*)


Print[" -- Sphere Packings initialized -- "];










