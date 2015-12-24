
(* Routines for Cylindrical algebraic decompositions :
	Code by Thomas C. Hales
	March 2001
	copyright Thomas C. Hales
	licensed under GPL

	References are to the book Quantifier Elimination and Cylindrical
	Algebraic Decomposition, Caviness and Johnson (eds.)

*)

KeepCAD[x_] := (Definition[x] >>> "cad.m")

degree[f_, x_] := Length[CoefficientList[f, x]] - 1

(* LeadingCoefficient *)

ldcf[f_,x_]:= CoefficientList[f,x]//Last;

(* Reductum *)

red[f_, x_, i_] := 
  Module[{j, k}, k = CoefficientList[f, x]; k = Drop[k, -i]; 
    k . Table[x^j, {j, 0, Length[k] - 1}]]

RED[f_, x_] := Module[{i},
	Complement[Table[red[f, x, i], {i, 0, degree[f, x]}], {0}]];

REDLIST[A_,x_]:= Table[RED[A[[i]],x],{i,1,Length[A]}]//Flatten//Union;

(* PSC: Principal Subresultant Coefficient *)

psc[f_, g_, x_, j_] := 
  Module[{f1,f2,m,n,k,matrix,keep}, 
	f1 = CoefficientList[f, x]//Reverse; 
	f2 = CoefficientList[g, x]//Reverse; 
    m = Length[f1] - 1; 
	n = Length[f2] - 1; 
	f1 = Join[f1,Table[0,{k,1,n-1}]];
	f2 = Join[f2,Table[0,{k,1,m-1}]];
	matrix = Join[Table[RotateRight[f1,k],{k,0,n-j-1}],
	Table[RotateRight[f2,k],{k,0,m-j-1}]
				];
	If[Length[matrix]==0,Return[1]]; (* I assume --- *)
	matrix = Transpose[matrix];
	keep = matrix[[m+n-2j]];
	Join[Drop[matrix,-1-2j],{keep}]//Det
			]

PSC[f_,g_,x_]:= Complement[
	Table[psc[f,g,x,j],{j,0,Min[degree[f,x],degree[g,x]]}],{0}];

PSC[f_,x_]:= PSC[f,D[f,x],x];

(* QECAD p 140 Projection Phase *)

proj1[f_,x_]:= Union[{ldcf[f,x]},PSC[f,x]];

PROJ1[A_,x_]:= Map[proj1[#,x]&,REDLIST[A,x]]//Flatten;

PROJ2[A_,x_]:= Module[{i,j},
	Table[Proj2[A[[i]],A[[j]],x],{i,1,Length[A]},{j,i+1,Length[A]}]//
	Flatten
	];

Proj2[f_,g_,x_]:= Map[PSC[#,g,x]&,RED[f,x]];  (* see QECAD, p 168 *)

PROJ[A_,x_]:= Union[PROJ1[A,x],PROJ2[A,x]];

(* QECAD p 146: Base Phase *)

RealRoots[f_,x_]:= Module[{k},
	k=x/.NSolve[f==0,x];
	If[k==x,Return[{}]];
	Select[k,(Head[#]==Real)&]];

RealRootList[A_,x_]:= Map[RealRoots[#,x]&,A]//Flatten;

BaseList[A_,x_]:= Module[{k,k2,j},
	k = RealRootList[A,x]//Union//Sort;
	If[Length[k]==0,Return[{{x->0}}]];
	k = Join[{k[[1]]-20},k,{Last[k]+20}];
	k2 = Table[k[[j]]/2 + k[[j+1]]/2, {j,1,Length[k]-1}];
	k2 = Join[Take[k,{2,Length[k]-1}],k2]//Union;
	Table[{x -> k2[[j]]},{j,1,Length[k2]}]
	];

(* QECAD p 147 : Extension Phase *)

CylinderExtend[A_,x_,cellPoint_]:= 
	Flatten[Map[CylinderExtendPoly[A,x,#]&,cellPoint],1];

CylinderExtendPoly[A_,x_,sub_]:= Module[{i},
	k = BaseList[A/.sub,x];
	Table[Join[k[[i]],sub], {i,1,Length[k]}]
	];



