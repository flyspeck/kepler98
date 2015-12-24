As[y1_,y2_,y3_,y4_,y5_,y6_]:= Module[{d1,d2,di,t1,t2},
	d1 = dihR[y1/2,eta[y1,y2,y5],eM/2];
	d2 = dihR[y1/2,eta[y1,y3,y6],eM/2];
	t1 = tauR[y1/2,eta[y1,y2,y5],eM/2];
	t2 = tauR[y1/2,eta[y1,y3,y6],eM/2];
	di = Dihedral[y1,y2,y3,y4,y5,y6];
	If[N[eta[y1,y2,y5]]>eM/2,d1=0;t1=0];
	If[N[eta[y1,y3,y5]]>eM/2,d2=0;t2=0];
	di = di - d1-d2;
	t1+t2+di*(1-Cos[arc[y1,eM/2,eM/2]])*phit - di*A[y1/2] - L[y1]
	];

L[x_] := 0.035*(x - 2) - 0.06*(x - 2)*(x - edgeMax);
s8:= Sqrt[8];
eM:= edgeMax;
tauR[x__]:= volR[x] - solR[x]/3/dtet;
volR[a_, b_, c_] := Sqrt[a^2*(b^2 - a^2)*(c^2 - b^2)]/6//Ns;
solR[x_, y_, z_] := 2*ArcTan[Sqrt[((z - y)*(y - x))/((z + y)*(x + y))]]//Ns;
dihR[x_, y_, z_] := ArcTan[Sqrt[(z^2 - y^2)/(y^2 - x^2)]]//Ns;  

FarFrom[x_,pair_]:= If[Norm[x,pair[[1]]]< Norm[x,pair[[2]]],
                        pair[[2]],pair[[1]]];
NearTo[x_,pair_]:= If[Norm[x,pair[[1]]]< Norm[x,pair[[2]]],
                        pair[[1]],pair[[2]]];
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

Norm[x_]:= x.x;
Norm[x_,y_]:= Norm[x-y];
Distance[x_,y_]:= Sqrt[(x-y).(x-y)];

Enclosed[x__, z1_, z2_, z3_] :=
  Module[{X, Y, Z}, If[Length[{x}] != 6, Return["invalid input"]];
    {X, Y, Z} = SimplexCoordinates[x];
    Sqrt[Norm[FarFrom[NullVector, Vertex[{X, z1}, {Y, z2}, {Z, z3}]]]]];
 
NullVector={0,0,0};

dist[x_, y_] := Sqrt[(x - y) . (x - y)]

beta[psi_,y1_,y3_,y5_]:= (* function IV.2.8 *)
    Module[{cb2,ctheta2},
        ctheta2 = ((y1^2+y3^2-y5^2)/(2 y1 y3))^2;
        cb2 = (Cos[psi]^2-ctheta2)/(1-ctheta2);
        ArcCos[Sqrt[cb2]]
        ];

beta[psi_,theta_]:= (* function IV.2.8 *)
    Module[{cb2,ctheta},
		ctheta = Cos[theta];
        cb2 = (Cos[psi]^2-ctheta^2)/(1-ctheta^2);
        ArcCos[Sqrt[cb2]]
        ];

uT[x1_,x2_,x6_]:= -x1*x1-x2*x2-x6*x6+2*x1*x6+2*x1*x2+2*x2*x6;

(* height in radians of spherical triangle with base alpha, sides beta, gamma*)
ht[beta_,gamma_,alpha_]:= Module[{alpha1},
	alpha1 = ArcTan[ -(Cos[alpha] Cos[beta]-Cos[gamma])/(Cos[beta] Sin[alpha])];
	ArcCos[ Cos[beta]/Cos[alpha1] ]
	];
