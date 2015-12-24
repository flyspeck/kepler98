(* Mathematica Code that calls CPLEX to deal
	with octahedral inequalities in SPhere Packings III,
	and vertexAdjustment (excessTcount in java parlance) in "Kepler"

	This program calls CPLEX.
	and cfsqp
*)

(* Call cT to generate the list.
	This program will have to be adapted for the cases (3,2),
	or (4,1) with a flat quarter.
*)

rS:= << "/afs/math.lsa.umich.edu/group/fac/hal\
es/pub/html/countdown/numerical/TAU41/Subdivide.m";

(* This sets the octahedral parameters *)
	(* Currently only one particular inequality is being used (the 1st of 18),
		and that with a slack of 0.002; *)

SetOct355:=(
	mode="oct355";
	CleanList=CleanList355;
	FF[x__] := {F[x], 2 - eta[{x}[[1]], {x}[[2]], {x}[[6]]]^2,
	   eta[{x}[[1]], {x}[[3]], {x}[[5]]]^2 - 2};
 
	F[x__] := Module[{y1, y2, y3, y4, y5, y6, Fsig, dih2, dih1, sig, slack},
	   {y1, y2, y3, y4, y5, y6} = {x}; 
		slack = 0.002;
		sig = Octavor[x];
		Fsig = sig - 0.1905999999999998 + dih2*0.1988669999999999;
		dih2 = Dihedral2[x]; 
		pi = Pi; 
		dih1 = Dihedral[x];
		-slack+Fsig + 0.051272725 + y1*-0.1710869999999999 + y2*0.282483 +
		 y3*0.04237589999999999 + y5*-0.0269357 + y6*-0.0813089 +
		 (dih1 - pi/2)*0.07611859999999999
		];
	z0 = Array[0&,6]; z1 = Array[1&,6];
	TextLineFn=TextLine;
	header = "MAX t\nST\n";
	xmin ={2.51, 2, 2, 2, 2, 2};
	xmax ={2.8284271247462, 2.255, 2.51, 2.255, 2.51, 2.255};
	tail = "x1<10\nx2<10\nBOUNDS\nt free\nend\n";
 
	outofrange[cmin_,cmax_]:= Module[{u},
		u={Drop[(FF @@ interp[cmin]),1], Drop[(FF @@ interp[cmax]),1]}
				//Transpose;
		If[Min[Map[Max,u]]<0, True, False]];
 
	TextLine[x_,j_]:= Module[{v = Chop[FF @@ interp[x]]},
		StringJoin[S[v[[2]]], " x1 ", SR[v[[3]]], " x2 + t < ", S[-v[[1]]]]];

	CleanList355[x_]:= Join[x[[2]],x[[3]],varSet355[x[[1]]]]//N;

	varSet355[x_]:= Module[{r,i,xx},
		If[x=={{1,0}},Return[{0,0,-1}]];
		r = {0,0,0};
		Do[
			xx=x[[i]];
			If[SameQ[xx[[1]],vx0],r[[3]]=Last[xx]];
			If[SameQ[xx[[1]],vx1],r[[1]]=Last[xx]];
			If[SameQ[xx[[1]],vx2],r[[2]]=Last[xx]];
			,{i,1,Length[x]}];
		r
		];
);



S[x_] := ToString[CForm[x]]
SR[x_] := StringJoin[If[N[x] >= 0, "+", ""], S[x]]


(*
xmin = domain lower bound
xmax = domain upper bound
ratmin = current lower = (1-ratmin) xmin + (ratmin) xmax;
*)

interp[rat_] := Array[(1-rat[[#]]) xmin[[#]]+(rat[[#]])xmax[[#]]&,Length[rat]];
MaxPos[x_]:= Position[x,Max[x]]//First;

Subdivide[ratmin_,ratmax_]:= Module[{pos,r,ww,siz},
	siz= Apply[Times,#]&;
	pos = MaxPos[interp[ratmax]-interp[ratmin]];
	r = (ratmin[[pos]]+ratmax[[pos]])/2;
	ww={{ratmin,ratmax},{ratmin,ratmax}};
	ww[[1,2,pos]]=r;
	ww[[2,1,pos]]=r;
	If[siz[ww[[1,2]]-ww[[1,1]]]+siz[ww[[2,2]]-ww[[2,1]]]!=siz[ratmax-ratmin],
		Print[ratmin," ",ratmax," ",ww];
		Print["siz ",siz[ww[[1,2]]-ww[[1,1]]]," ",
			siz[ww[[2,2]]-ww[[2,1]]]," ",siz[ratmax-ratmin]];
		Input["ERROR in Subdivide"];
		];
	ww
	];

(* file:
/tmp/subdiv.exec:
set logfile /tmp/subdiv.log
set preprocessing presolve no
set output results n /tmp/y.log
set output dialog n /tmp/x.log
read /tmp/subdiv.lp lp
opt
di so -
quit

/tmp/getAddValues  (add the following):

# If starts with x or t keep as a string
/^[xt]/{
s/^x\([^ ]*\) /{ vx\1, /g
s/e/*10^/g
s/^t/{ vx0,/g
s/$/},/g
p
}


*)
(* binary10 = Table[Mod[Floor[i/2^(j-1)],2],{i,0,1023},{j,1,10}]; *)
array[i_,cmin_,cmax_,r_]:= Module[{j,t},
	t =cmin;
	If[Length[r]>10,Print["ERROR!!"]];
	Do[If[binary10[[i+1,j]]>0,t[[r[[j]]]]=cmax[[r[[j]]]] ] ,
		{j,1,Length[r]}];
	t
	];

CorrectionTerms[cmin_,cmax_] := 
  Module[{stream, x, clopen, out, vv,r,type,rj,j,i}, 
	(* If[outofrange[cmin,cmax],(values = {{1,0}}); Return[1]]; *)
	clopen[x_] := Close[OpenWrite[x]]; 
    stream = OpenWrite["/tmp/subdiv"]; 
	streamABORT=stream;
    out[x_] := WriteString[stream, StringJoin[x, "\n"]]; 
	out[header];
	r = Select[Range[Length[cmax]],(cmax-cmin)[[#]]>10^-6&];

	type[1]={2,3,5,6};
	type[2]={3,8,9,5};
	type[3]={8,12,11,9};
	type[4]={12,14,15,11};
	type[5]={14,2,6,15};
	Do[
		rj = Intersection[r,type[j]];
		vv = Table[
			 TextLineFn[array[i,cmin,cmax,rj],j]<>"\n",
			{i,0,2^Length[rj] -1} ]//StringJoin;
		out[vv]; 
		,{j,1,5}];
	out[tail];
	Close[stream];
    clopen["/tmp/subdiv.log"]; clopen["/tmp/x.log"]; clopen["/tmp/y.log"]; 
    Run["cplex < /tmp/subdiv.exec"]; clopen["/tmp/subdiv.sum"]; 
	Run["sed -f /tmp/getAddValues /tmp/subdiv.log > /tmp/subdiv.sum"];
	Run["cat /tmp/getAddTail >> /tmp/subdiv.sum"];
    << "/tmp/subdiv.sum"; 
	{If[values[[1,1]]>epsilon,
		values[[1,1]]=
				{ values[[1,1]],
				 cfsqp[CleanList[{values,interp[cmin],interp[cmax]}]]
				}//Min
		]};
    Print[values]; 
	 values[[1,1]]
	];

epsilon = 10^-5;

FindCorrections[cmin_,cmax_]:= Module[{w},
	If[CorrectionTerms[cmin,cmax]>epsilon,
		AppendTo[CorrectionList,{values,cmin,cmax}],
		w = Subdivide[cmin,cmax];
		Print["Subdivide = ",N[cmin]," ",N[cmax-cmin]];
		Print["Length[CorrectionList]= ",Length[CorrectionList]];
		FindCorrections[w[[1,1]],w[[1,2]]];
		FindCorrections[w[[2,1]],w[[2,2]]]
		];
	r
	];

cT:= (CorrectionList={};FindCorrections[z0,z1];
	Array[CleanList[CorrectionList[[#]]]&,Length[CorrectionList]]);

size[x_]:= Apply[Times,x[[3]]-x[[2]]];
CheckList:= Apply[Plus,Map[size,CorrectionList]];

(* Start on the VertexAdjustment[4] inequalities *)

CleanListTau[x_]:=
	{
	Join[x[[2]][[{1,2,3, 4, 5,6}]],
		 x[[3]][[{1,2,3, 4, 5,6}]],varSetTauA[x[[1]]] ],
	Join[x[[2]][[{1,3,8, 7, 9,5}]],
		 x[[3]][[{1,3,8, 7, 9,5}]],varSetTauB[x[[1]]] ],
	Join[x[[2]][[{1,8,12,10,11,9}]],
		 x[[3]][[{1,8,12,10,11,9}]],varSetTauC[x[[1]]] ],
	Join[x[[2]][[{1,12,14,13,15,11}]],
		 x[[3]][[{1,12,14,13,15,11}]],varSetTauD[x[[1]]] ],
	Join[x[[2]][[{1,14,2, 16, 6,15}]],
		 x[[3]][[{1,14,2, 16, 6,15}]],varSetTauE[x[[1]]] ]
	}//N;
(*
*)

varSetTauA[x_]:= Module[{r,i,xx},
	r = Array[0&,8];
	r[[1]]= -1; (* tau coeff *)
	Do[
		xx=x[[i]];
		If[SameQ[xx[[1]],vx0],r[[2]]=Last[xx]];
		If[SameQ[xx[[1]],vxA],r[[3]]=Last[xx]];
		If[SameQ[xx[[1]],vxA2],r[[4]]=Last[xx]];
		If[SameQ[xx[[1]],vxA3],r[[5]]=Last[xx]];
		If[SameQ[xx[[1]],vxA5],r[[6]]=Last[xx]];
		If[SameQ[xx[[1]],vxA6],r[[7]]=Last[xx]];
		If[SameQ[xx[[1]],vxdih],r[[8]]=Last[xx]];
		,{i,1,Length[x]}];
	r
	];

varSetTauB[x_]:= Module[{r,i,xx},
	r = Array[0&,8];
	r[[1]]= -1; (* tau coeff *)
	Do[
		xx=x[[i]];
		If[SameQ[xx[[1]],vx0],r[[2]]=Last[xx]];
		If[SameQ[xx[[1]],vxB],r[[3]]=Last[xx]];
		If[SameQ[xx[[1]],vxB3],r[[4]]=Last[xx]];
		If[SameQ[xx[[1]],vxB8],r[[5]]=Last[xx]];
		If[SameQ[xx[[1]],vxB9],r[[6]]=Last[xx]];
		If[SameQ[xx[[1]],vxB5],r[[7]]=Last[xx]];
		If[SameQ[xx[[1]],vxdih],r[[8]]=Last[xx]];
		,{i,1,Length[x]}];
	r
	];

varSetTauC[x_]:= Module[{r,i,xx},
	r = Array[0&,8];
	r[[1]]= -1; (* tau coeff *)
	Do[
		xx=x[[i]];
		If[SameQ[xx[[1]],vx0],r[[2]]=Last[xx]];
		If[SameQ[xx[[1]],vxC],r[[3]]=Last[xx]];
		If[SameQ[xx[[1]],vxC8],r[[4]]=Last[xx]];
		If[SameQ[xx[[1]],vxC12],r[[5]]=Last[xx]];
		If[SameQ[xx[[1]],vxC11],r[[6]]=Last[xx]];
		If[SameQ[xx[[1]],vxC9],r[[7]]=Last[xx]];
		If[SameQ[xx[[1]],vxdih],r[[8]]=Last[xx]];
		,{i,1,Length[x]}];
	r
	];

varSetTauD[x_]:= Module[{r,i,xx},
	r = Array[0&,8];
	r[[1]]= -1; (* tau coeff *)
	Do[
		xx=x[[i]];
		If[SameQ[xx[[1]],vx0],r[[2]]=Last[xx]];
		If[SameQ[xx[[1]],vxD],r[[3]]=Last[xx]];
		If[SameQ[xx[[1]],vxD12],r[[4]]=Last[xx]];
		If[SameQ[xx[[1]],vxD14],r[[5]]=Last[xx]];
		If[SameQ[xx[[1]],vxD15],r[[6]]=Last[xx]];
		If[SameQ[xx[[1]],vxD11],r[[7]]=Last[xx]];
		If[SameQ[xx[[1]],vxdih],r[[8]]=Last[xx]];
		,{i,1,Length[x]}];
	r
	];

varSetTauE[x_]:= Module[{r,i,xx},
	r = Array[0&,8];
	r[[1]]= 0; (* tau coeff *)
	Do[
		xx=x[[i]];
		If[SameQ[xx[[1]],vx0],r[[2]]=Last[xx]];
		If[SameQ[xx[[1]],vxE],r[[3]]=Last[xx]];
		If[SameQ[xx[[1]],vxE14],r[[4]]=Last[xx]];
		If[SameQ[xx[[1]],vxE2],r[[5]]=Last[xx]];
		If[SameQ[xx[[1]],vxE6],r[[6]]=Last[xx]];
		If[SameQ[xx[[1]],vxE15],r[[7]]=Last[xx]];
		If[SameQ[xx[[1]],vxdih],r[[8]]=Last[xx]];
		,{i,1,Length[x]}];
	r
	];




DihedralA[x_]:= (Dihedral @@ x[[{1,2,3,4,5,6}]])-2.Pi/5;
DihedralB[x_]:= (Dihedral @@ x[[{1,3,8,7,9,5}]])-2.Pi/5;
DihedralC[x_]:= (Dihedral @@ x[[{1,8,12,10,11,9}]])-2.Pi/5;
DihedralD[x_]:= (Dihedral @@ x[[{1,12,14,13,15,11}]])-2.Pi/5;
DihedralE[x_]:= (Dihedral @@ x[[{1,2,14,16,15,6}]])-2.Pi/5;


tauA[x_]:= (tau @@ x[[{1,2,3,4,5,6}]]);
tauB[x_]:= (tau @@ x[[{1,3,8,7,9,5}]]);
tauC[x_]:= (tau @@ x[[{1,8,12,10,11,9}]]);
tauD[x_]:= (tau @@ x[[{1,12,14,13,15,11}]]);
tauE[x_]:= 0;

TextLineTau[y_,j_]:= Module[{x},
	x = interp[y];
	Switch[j,
	1,{"t + xA  ",SR[x[[2]]]," xA2 ",SR[x[[3]]]," xA3 ",
	 SR[x[[5]]], " xA5 ",SR[x[[6]]]," xA6 ",
	 SR[DihedralA[x]]," xdih < ",
	 SR[tauA[x]],"\n"},
	 2,{"t + xB  ",SR[x[[3]]]," xB3 ",SR[x[[8]]]," xB8 ",
	 SR[x[[5]]], " xB5 ",SR[x[[9]]]," xB9 ",
	 SR[DihedralB[x]]," xdih < ",
	 SR[tauB[x]],"\n"},
	 3,{"t + xC  ",SR[x[[8]]]," xC8 ",SR[x[[11]]]," xC11 ",
	 SR[x[[9]]], " xC9 ",SR[x[[12]]]," xC12 ",
	 SR[DihedralC[x]]," xdih < ",
	 SR[tauC[x]],"\n"},
	 4,{"t + xD  ",SR[x[[12]]]," xD12 ",SR[x[[15]]]," xD15 ",
	 SR[x[[11]]], " xD11 ",SR[x[[14]]]," xD14 ",
	 SR[DihedralD[x]]," xdih < ",
	 SR[tauD[x]],"\n"},
	 5,{"t + xE  ",SR[x[[2]]]," xE2 ",SR[x[[15]]]," xE15 ",
	 SR[x[[6]]], " xE6 ",SR[x[[14]]]," xE14 ",
	 SR[DihedralE[x]]," xdih < ",
	 SR[tauE[x]],"\n"}
	]//StringJoin
	];

SetTau:=
	(
	mode="tau5";
	CleanList=CleanListTau;
	TextLineFn=TextLineTau;
	xmin=Array[2&,16]; xmin[[16]]=2.sq;
	xmax=Array[2.51&,16]; 
	xmax[[5]]=xmax[[9]]=xmax[[11]]=2.26;
	xmax[[1]]=2;
	xmax[[16]]=2.sq;
	xmax[[4]]=xmax[[7]]=xmax[[10]]=xmax[[13]]=2;
	xmax[[16]]=2.sq;
	xmax[[3]]=xmax[[8]]=xmax[[12]]=2.45;
	z0 = Array[0&,16]; z1= Array[1&,16];
	header="MAX t\nST\n";
	tail = StringJoin[
		{
		"xA+xB+xC+xD+xE >0.0831\n", (* 0.0831 approx 1.5 pt, adjust as needed *)
		"xA2+xE2>0\n",
		"xA3+xB3>0\n",
		"xB8+xC8>0\n",
		"xC12+xD12>0\n",
		"xD14+xE14>0\n",
		"xA5+xB5>0\n",
		"xB9+xC9>0\n",
		"xC11+xD11>0\n",
		"xD15+xE15>0\n",
		"xE6+xA6>0\n",

		"BOUNDS\n",
		"-10< xA < 10\n","-10< xA2 < 10\n","-10< xA3 < 10\n",
				"-10< xA5 < 10\n","-10< xA6 < 10\n",
		"-10< xB < 10\n","-10< xB3 < 10\n","-10< xB9 < 10\n",
				"-10< xB5 < 10\n","-10< xB8 < 10\n",
		"-10< xC < 10\n","-10< xC8 < 10\n","-10< xC11 < 10\n",
				"-10< xC9 < 10\n","-10< xC12 < 10\n",
		"-10< xD < 10\n","-10< xD12 < 10\n","-10< xD15 < 10\n",
				"-10< xD11 < 10\n","-10< xD14 < 10\n",
		"-10< xE < 10\n","-10< xE2 < 10\n","-10< xE15 < 10\n",
				"-10< xE6 < 10\n","-10< xE14 < 10\n",
		"t free\n",
		"-10 < xdih < 10\n",
		"end\n"
		}];
	);


ConvertString[x_]:= StringReplace[S[x],
	{"{"->" ","("->" ",")"->" ","}"->"\n",","->" ","List"->"\n"}];

cfsqp[x_]:= Module[{stream,i},
	stream = OpenWrite["/tmp/subdiv_cfsqp.dat"];
	Do[WriteString[stream,ConvertString[x[[i]]]],{i,1,Length[x]}];
	streamABORT=stream;
	Close[stream];
	Run["/tmp/partU"];
	<< /tmp/subdiv_cfsqp.m;
	cfsqpVal
	];  
