(* 2001 package for sphere packings *)

Keep[x_]:= PutAppend[Definition[x],"more.m"];
SetAttributes[Keep,HoldFirst];

Fix[x_] := Edit[Definition[x]]
SetAttributes[Fix,HoldFirst];
(Unprotect[In,Out]; Format[In]=SphereIn;
    Format[Out]=SphereOut; Protect[In,Out];) 

$SpherePrecision:= 16;
Ns[x_]:= N[x,$SpherePrecision];   

Delta[y1_,y2_,y3_,y4_,y5_,y6_]:= Module[{x1,x2,x3,x4,x5,x6},
    x1= y1*y1; x2=y2*y2; x3=y3*y3;
    x4= y4*y4; x5=y5*y5; x6=y6*y6;
    (x1*x4*(-x1+x2+x3-x4+x5+x6)+
            x2*x5*(x1-x2+x3+x4-x5+x6)
            +x3*x6*(x1+x2-x3+x4+x5-x6)
            -x2*x3*x4-x1*x3*x5-x1*x2*x6-x4*x5*x6)];    

(* edgeMax := 2.558913201377738; (old value)*)
edgeMax:= Ns[25589/10^4];
trunc := edgeMax/2;

phi[h_,t_]:= h t (h+t)/6;

phi0 := phi[trunc,trunc];

A[h_,t_]:= (1-h/t) (phi[t,t]-phi[h,t]);
A[h_]:= A[h,trunc];

Solid[y1_,y2_,y3_,y4_,y5_,y6_]:= Module[{},
	Dihedral[y1,y2,y3,y4,y5,y6]+Dihedral[y2,y3,y1,y5,y6,y4]
	+Dihedral[y3,y1,y2,y6,y4,y5]-Pi];

Volt[y1_,y2_,y3_,y4_,y5_,y6_]:= Module[{},
	phi0 Solid[y1,y2,y3,y4,y5,y6] - A[y1/2] Dihedral[y1,y2,y3,y4,y5,y6]
	-A[y2/2] Dihedral[y2,y1,y3,y5,y4,y6]
	-A[y3/2] Dihedral[y3,y1,y2,y6,y4,y5]
	+Quoin[y1/2,eta[y1,y2,y6],trunc]+Quoin[y2/2,eta[y1,y2,y6],trunc]
	+Quoin[y2/2,eta[y2,y3,y4],trunc]+Quoin[y3/2,eta[y2,y3,y4],trunc]
	+Quoin[y3/2,eta[y3,y1,y5],trunc]+Quoin[y1/2,eta[y3,y1,y5],trunc]
	];

Volq[y1_,y2_,y3_,y4_,y5_,y6_,t_]:= Module[{},
	phi[t,t] Solid[y1,y2,y3,y4,y5,y6] 
	- A[y1/2,t] Dihedral[y1,y2,y3,y4,y5,y6]
	-A[y2/2,t] Dihedral[y2,y1,y3,y5,y4,y6]
	-A[y3/2,t] Dihedral[y3,y1,y2,y6,y4,y5]
	+Quoin[y1/2,eta[y1,y2,y6],t]+Quoin[y2/2,eta[y1,y2,y6],t]
	+Quoin[y2/2,eta[y2,y3,y4],t]+Quoin[y3/2,eta[y2,y3,y4],t]
	+Quoin[y3/2,eta[y3,y1,y5],t]+Quoin[y1/2,eta[y3,y1,y5],t]
	];

assert[x_]:= If[!x,Print["ASSERTION FAILED"]];

volR[a_,b_,c_]:=Sqrt[a^2*(b^2 - a^2)*(c^2 - b^2)]/6;

VolAn[y1_,y2_,y3_,y4_,y5_,y6_]:= Module[{r = Rad[y1,y2,y3,y4,y5,y6]},
	volR[y1/2,eta[y1,y2,y6],r]+volR[y2/2,eta[y1,y2,y6],r]
	+volR[y2/2,eta[y2,y3,y4],r]+volR[y3/2,eta[y2,y3,y4],r]
	+volR[y3/2,eta[y3,y1,y5],r]+volR[y1/2,eta[y3,y1,y5],r]
	];

Dihedral[y1_, y2_, y3_, y4_, y5_, y6_] :=
  Module[{d4, x1 = y1*y1, x2 = y2*y2, x3 = y3*y3, x4 = y4*y4, x5 = y5*y5,
    x6 = y6*y6}, d4 =
     -x1^2 + x1*x2 + x1*x3 - x2*x3 - 2*x1*x4 + x1*x5 + x2*x5 + x1*x6 +
      x3*x6 - x5*x6; Pi/2 + ArcTan[-(d4/
         Sqrt[4*x1*Delta[y1, y2, y3, y4, y5, y6]])]]  ;

Rho[y1_,y2_,y3_,y4_,y5_,y6_]:= Module[{x1,x2,x3,x4,x5,x6},
x1= y1*y1; x2=y2*y2; x3=y3*y3;
x4= y4*y4; x5=y5*y5; x6=y6*y6;
(-x1*x1*x4*x4-x2*x2*x5*x5-x3*x3*x6*x6+
        2*x1*x2*x4*x5+2*x1*x3*x4*x6+2*x2*x3*x5*x6)];

Rad[y1_,y2_,y3_,y4_,y5_,y6_]:=
    Sqrt[Rho[y1,y2,y3,y4,y5,y6]/Delta[y1,y2,y3,y4,y5,y6]]/2; 

eta[x_,y_,z_]:= Module[{s= (x+y+z)/2}, x y z/4/Sqrt[s (s-x)(s-y)(s-z)]]//Ns;

Quoin[a_, b_, c_] :=
  Module[{u = Sqrt[(c^2 - b^2)/(b^2 - a^2)]},
    If[N[a]>N[b],Return[0]];
    If[N[b]>N[c],Return[0]];
   -(a^2 + a*c - 2*c^2)*(c - a)*ArcTan[u]/6 + a*(b^2 - a^2)*u/6 -
    (2/3)*c^3*ArcTan[(u*(b - a))/(b + c)]]//Ns; 

Qy[y1_,y2_,y3_]:= Quoin[y1/2,eta[y1,y2,y3],trunc];





(* Rogers's density *)
dtet := Sqrt[8]*ArcTan[Sqrt[2]/5]

(* Squander of a simplex *)
taut[x__] := Volt[x] - Solid[x]/(3*dtet)

(* adjustment to phi to deal with squander rather than c-vol *)
phit := phi[trunc, trunc] - 1/(3*dtet)

(* law of cosines, c opposite *)
arc[a_, b_, c_] := ArcCos[(a^2 + b^2 - c^2)/(2*a*b)] 

(* chiMinus *)
chiMinus[y1_,y2_,y3_,y4_,y5_,y6_,lambda_]:= Module[{psi,eta126,eta135,y1p,dihS},
	psi = arc[y1,eM/2,lambda];
	eta126 = eta[y1,y2,y6];
	eta135 = eta[y1,y3,y5];
	y1p = y1/(2 Cos[psi]);
	dihS = Dihedral[y1,y2,y3,y4,y5,y6];
	dihS (1-Cos[psi]) phit
		- solp[y1,eta126,y1p,psi]phit - solp[y1,eta135,y1p,psi]phit
		- A[y1/2] dihS + Qy[y1,y2,y6] + Qy[y1,y3,y5]
		- mu[y5] L[y1]/2 - mu[y6] L[y1]/2
	];

mu[y_]:= If[y>=eM-10^-10,0,1];  (* fixed May 26, 2001 *)

solp[y1_,eta_,yp_,psi_]:= Module[{},
	dihR[y1/2,eta,yp] (1-Cos[psi]) - solR[y1/2,eta,yp]
	];

(* chiPlus at angle Pi *)
chiPi[y1_,y2_,y3_,y5_,y6_,lambda_]:= Module[{psi,eta126,eta135,y1p},
	psi = arc[y1,eM/2,lambda];
	eta126 = eta[y1,y2,y6];
	eta135 = eta[y1,y3,y5];
	y1p = y1/(2 Cos[psi]);
	Pi (1-Cos[psi])phit - Pi A[y1/2] - solp[y1,eta126,y1p,psi]phit
		-solp[y1,eta135,y1p,psi]phit
	];

