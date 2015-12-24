(* Functions for a proof of the honeycomb conjecture *)

PI[n_]:= n Tan[Pi/n];

eps[n_,T_]:= T Sqrt[PI[6]] + (n-6) 0.0505 - 2 Sqrt[PI[6]];

L[n_,x_]:= 2 Sqrt[(1-x) PI[n]];
Lplus = 2 Sqrt[Pi]; 
Lminus[x_]:= 2 Sqrt[(1-2x) Pi];
tau0 = 1/2;
Ld[x_]:= -x Sqrt[2Pi/tau0];  (* for x<0 *)
L1[n_,x_]:= L[n,x] HoneyEdge[1,Abs[x]/L[n,x]];

HoneyEdge[polyedge_, x_] := 
	(* called arc in paper,
		arclength of an arc cut by a chord of length polyedge and
		enclosing area x *)
  Module[{ChordArea, PolygonRadius, PolygonEdge, rho, theta, Bisectf}, 
   ChordArea = x; PolygonEdge = polyedge; 
    Bisectf = #1/Sin[#1]^2 - Cos[#1]/Sin[#1] - 
       (4*ChordArea)/PolygonEdge^2 & ; 
    theta := Bisect[{10^(-8), 3.14}, 0.0001, Bisectf][[1]]; 
    rho := PolygonEdge/(2*Sin[theta]); Return[2*theta*rho]; ]

ConstrainedNgonArea[N_,r_] := 
	(* Area of an N-gon  inscribed in a circle, sides t (N-1 times) and max(1,t)
		*)
	Module[{t0,t,perim,area,alpha,beta},
		t0 = Max[1,2 r Sin[Pi/N]];
		alpha = ArcSin[t0/(2r)];
		beta = (Pi-alpha)/(N-1);
		t = 2r Sin[beta];
		perim = (N-1)t + t0;
		area = r^2 ((N-1) Cos[beta]Sin[beta] + Cos[alpha]Sin[alpha]);
		{perim,area}
		];


XNplus[N_]:= 1 - Pi/PI[N];
XNminus[N_]:= (1-Pi/PI[N])/(1-2Pi/PI[N]);

Fp[t_]:= (t-Sin[t] Cos[t])/Sin[t]^2;

Fq[t_]:= t/Sin[t];

Fr[t_]:= (Sin[t]- t Cos[t])/(t Sin[t]^2);


