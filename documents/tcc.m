(* This file contains the procedures needed for the bounds on
	convex polygons in Sphere Packings IV *)

(* Config = {{v1,v2,...,vr},{e1,e2,...,er}}  the a's and b's of SPIV.5.8.
	vi,ei integers 0..4.
	converted to length by FixedLength,
	3 represents the range (2sq,3.2).   *)

FirstPass:= SameQ[Head[FirstPassToken],Symbol];

(*                 0,  1,  2,   3,      4   *)
FixedLength[i_]:= {2,2.51,2sq,"invalid",3.2}[[i+1]];
IPart[v_,i_]:= v[[IMod[i,Length[v] ] ]];
IMod[i_,r_]:= 1+Mod[i-1,r];



(******************   CONSTANTS ***********************)

(* The sum of tauD+pitau.   tauD is the target
	for what we want squandered, penalty comes from the
	erasure of upright quarters, eps comes from flats that
	are score with vorVc that arise by deformation from
	specials.  *)

PenalizedD[config_]:= Module[{i,k0,v,e,n,k,k2,pitau},
	{v,e}=config;
	k2 = Select[e,#>1&]//Length;
	k = Select[e,#>0&]//Length;
	n = Length[e];
	k0 = n-k;
	pitau= 0.04683 + (k0+2 k2-3)/3 0.008 + 0.0066 k2;
	tauD[n,k]+ pitau
	];

PenalizedZ[config_]:= Module[{i,v,e,n,k,k2,k0,pisigma},
	{v,e}=config;
	k2 = Select[e,#>1&]//Length;
	k = Select[e,#>0&]//Length;
	n = Length[e];
	k0 = n-k;
	pisigma =  (k0+2k2)/3 0.008 + 0.009 k2; 
	If[k0+2k2 < 3, pisigma = 0.009 k2];
	sigmaZ[n,k]- pisigma
	];

(* The bounds desired from Sphere Packings IV *)

tauD[n_, k_] := tauT[n + k] - k*0.06585;
tauT[4] = 0.13170;
tauT[5] = 0.27113;
tauT[6] = 0.41056;
tauT[7] = 0.54999;
tauT[8] = 0.6045;
tauT[9] = 0.6978;
tauT[10] = 0.7891;
tauT[n_] := (4Pi zeta -8) pt/; n > 10

sigmaZ[3,1]= 0.00005;
sigmaZ[n_, k_] := sigmaS[n + k] - k*0.00005;

sigmaS[4] := 0;
sigmaS[5] := -0.05704;
sigmaS[6] := -0.11408;
sigmaS[7] := -0.17112;
sigmaS[8] := -0.22816;
sigmaS[9] := -0.19720;

sigmaS[n_] := 0 /; n > 9

(* These are the bounds on a special simplex that has been
	severed from the convex polygon.
	See SP.IV.A13 for these bounds.  *)
SpecialTau[0]= 0.014; (* <tauVc[2,2,2,/2sq,3.2/,2,2] *)
SpecialTau[1]= 0; (* Min of tauVc[2.51,2,2,/2sq,3.114466.../,2,2] *)
SpecialSigma[0]= 0.0461; (* Max of VorVc[2,2,2,/2sq,3.2/,2,2] *)
SpecialSigma[1]= 0; (* Max of VorVc[2.51,2,2,/2sq,3.114466.../,2,2] *)

(******** END CONSTANTS ************)




(************* SCORING TCCS ************************)

(* What a truncated corner cell squanders: *)
CorSigma[x__] := Module[{y1, y2, y3, y4, y5, y6, long, Cospsi, h, Costhet}, 
   {y1, y2, y3, y4, y5, y6} = Ns[{x}]; 
	long = Ns[16/10]; 
	h = y1/2; 
    If[Delta[y1, y2, y3, y4, y5, y6] <= 0, Return[0]]; 
    Cospsi = Cos[arc[y1,TRUNC,long]];
    Costhet = Min[h/TRUNC, 1]; 
    Dihedral[x]*((1 - Cospsi)*phi0 + (1 - Costhet)*(phi[h,TRUNC]-phi0))
	- TccSolPrime[y1, y2, y6]*phi0
	- TccSolPrime[y1, y3, y5]*phi0 
	- 4 doct (Quoin[y1/2,eta[y1,y2,y6],TRUNC]+Quoin[y1/2,eta[y1,y3,y5],TRUNC])
     ];

CorTauFixedLength[x__]:=  CorTauFixedLength[x] = CorTau @@ Map[FixedLength,{x}];

(* sol' of SP.IV.5.3. *)
TccSolPrime[y1_, y2_, y6_] := Module[{h,long, cospsi, e}, 
	long = Ns[16/10]; 
	h = y1/2;
	e = eta[y1, y2, y6]; 
	cospsi = Cos[arc[y1,TRUNC,long]];
	If[e >= h/cospsi, Return[0]]; 
    +(dihR[h, e, h/cospsi]*(1 - cospsi)) - solR[h,e,h/cospsi]
	];

(* What a truncated corner cell scores *)
(* sol(C0) of SP.IV.5.3.  *)
CorTau[x__] := CorSolid[x]*zeta*pt - CorSigma[x];
CorSigmaFixedLength[x__]:=  CorSigmaFixedLength[x] = CorSigma @@ Map[FixedLength,{x}];

CorSolid[x__] := Module[{y1, y2, y3, y4, y5, y6, long, Cospsi}, 
   {y1, y2, y3, y4, y5, y6} = Ns[{x}]; 
	long = Ns[16/10]; 
    If[Delta[y1, y2, y3, y4, y5, y6] <= 0, Return[0]]; 
    Cospsi = Cos[arc[y1,TRUNC,long]];
    Dihedral[x]*(1 - Cospsi) - TccSolPrime[y1, y2, y6] - TccSolPrime[y1, y3, y5]
	]


(* The sum of what the truncated corner cells in a convex polygon
	squander *)
TccSquanderFixedLength[config_]:= Module[{i,v,e},
	{v,e}=config;
	If[Max[v]>1 || Min[v]<0 || Max[e]>4 || Min[e]<0, Return["bad data"]];
	Apply[Plus,Table[CorTauFixedLength[
		v[[i]],IPart[v,i-1],IPart[v,i+1],4,IPart[e,i],IPart[e,i-1]
		],
		{i,1,Length[v]}
	   ]]
	];

(* The same as TccSquanderFixedLength, but becoming tauVc on simplices *)
TccOrTauVcFixedLength[config_]:= Module[{i,v,e},
	{v,e}=config;
	If[Length[v]>3,Return[TccSquanderFixedLength[config]]];
	tauVcFixedLength[v[[1]],v[[2]],v[[3]],e[[2]],e[[3]],e[[1]]]
	];

tauVcFixedLength[x__]:= tauVcFixedLength[x] = tauVc @@ Map[FixedLength,{x}];
VorVcFixedLength[x__]:= VorVcFixedLength[x] = VorVc @@ Map[FixedLength,{x}];

(* The sum of what a truncated corner cells in a convex
	polygon  score *)
TccScoreFixedLength[config_]:= Module[{i,v,e},
	{v,e}=config;
	If[Max[v]>1 || Min[v]<0 || Max[e]>4 || Min[e]<0, Return["bad data"]];
	Apply[Plus,Table[CorSigmaFixedLength[
		v[[i]],IPart[v,i-1],IPart[v,i+1],4,IPart[e,i],IPart[e,i-1]
		],
		{i,1,Length[v]}
	   ]]
	];

(* The same as TccScoreFixedLength, but becoming VorVc on simplices *)
TccOrVorVcFixedLength[config_]:= Module[{i,v,e},
	{v,e}=config;
	If[Length[v]>3,Return[TccScoreFixedLength[config]]];
	VorVcFixedLength[v[[1]],v[[2]],v[[3]],e[[2]],e[[3]],e[[1]]]
	];

(* Draw an edge of length 3.2 from corner 1 to 3..config, breaking the
	convex polygon in 2 (or more by recursion) pieces.
	Calculate the sum of what the two or more pieces squander *)
RecursiveTau[config_]:= Module[{i,c1,c2},
	If[Length[config[[1]]]<4,Return[TccOrTauVcFixedLength[config]]];
	{TccOrTauVcFixedLength[config],
	Table[
		{c1,c2}=Cut[config,i];
		TccOrTauVcFixedLength[c1]+RecursiveTau[c2],
		{i,3,Length[config[[1]]]-1}
		]}//Min
	];

(* Same as RecursiveTau, but for score rather than squander. *)
RecursiveSigma[config_]:= Module[{i,c1,c2},
	If[Length[config[[1]]]<4,Return[TccOrVorVcFixedLength[config]]];
	{TccOrVorVcFixedLength[config],
	Table[
		{c1,c2}=Cut[config,i];
		TccOrVorVcFixedLength[c1]+RecursiveSigma[c2],
		{i,3,Length[config[[1]]]-1}
		]}//Max
	];


(***************** END OF TCC SCORING ********************)

(**************** CHECKING *************************)

(* Break config into two pieces by a new edge running from
	corner 1 to corner i.  
	The length of the new edge is 3.2==FixedLength[4].
	This is used in RecursiveTau/Sigma *)
(* Example: Cut[{{a1,b1,c1,d1,e1},{a2,b2,c2,d2,e2}},3] ==
{{{a1, b1, c1}, {a2, b2, 4}}, {{c1, d1, e1, a1}, {c2, d2, e2, 4}}}  *)

Cut[config_,i_]:= Module[{c1,c2,v,e,j},
	{v,e}= config;
	c1 = {v[[Range[i] ]],e[[Range[i] ]]}; c1[[2,i]]=4;
	v = RotateLeft[v,i-1];
	e = RotateLeft[e,i-1];
	j = Length[v]-i+2;
	c2 = {v[[Range[j] ]],e[[Range[j] ]]}; c2[[2,j]]=4;
	{c1,c2}
	];

DoesntExist[config_]:= Module[{v,e,i},
	(* Delta[2.51,2,2,3.2,2,2]<0 *)
	{v,e}=config;
	Count[Table[{v[[i]],IPart[v,i-1],IPart[v,i+1],
				IPart[e,i-1],e[[i]]},{i,1,Length[v]}],
		{1,0,0,0,0}]>0];

CheckTcc[config_]:= Module[{eps=10.0^-6},
	(DoesntExist[config]) ||
	((TccSquanderFixedLength[config]>PenalizedD[config]+eps)&&
	 (TccScoreFixedLength[config]<PenalizedZ[config]-eps)) 
	];

(* This is the first major verification that is needed.
	It verifies that with pure truncated corner cell scoring
	(no specials), all the bounds are met. *)

CheckAllTcc:= Module[{veri,list},
	veri[list_]:= Array[CheckTcc[list[[#]]]&,Length[list]]//Union;
	Join[
		veri[QuadList],veri[PentList],veri[HexList],veri[HeptList]
		]
	];

(* Now we turn to the procedures needed for the bounds when
	there are special simplices.  
	Calculation of the vertices at which a special simplex is
	centered *)
(* Examples: CutSpecial[{{a1,b1,c1,d1,e1},{a2,b2,c2,d2,e2}},3]==
		{{{d1, e1, a1, b1}, {d2, e2, a2, 2}}, c1} 
   CutSpecial[{{a1,b1,c1,0,e1,0,g1},{a2,b2,c2,0,0,f2,g2}},3]
		== {{{b1, a1, g1, 0, e1, 0}, {a2, g2, f2, 0, 0, 2}}, c1} *)

CutSpecial[config_,at_Integer(* center of special*)]:= Module[{},
	{v,e}=config;
	v=RotateLeft[v,at-1];
	specialHt=v[[1]];
	v=Drop[v,1]//RotateRight;
	e=Drop[RotateLeft[e,at-1],1]//RotateRight;
	e[[1]]=2;
	(* try to put the newly exposed edge last,
		and a nonspecial at 2, The reason is that we
		cut from the first vertex to the ith (3...len-1),
		and we don't want a special interfering. *)
	v=RotateLeft[v]; e=RotateLeft[e];
	If[v[[1]]==0 && v[[3]]==0 && e[[1]]==0 && e[[2]]==0,
		v=Reverse[v]; e=Reverse[e]//RotateLeft]; 
	{{v,e},specialHt}
	];


(* called by CheckConfigWithSpecial *)
WithSpecial[config_,defig_,ht_]:= Module[{eps=10.0^-6},
	(DoesntExist[defig]) ||
	((RecursiveTau[defig]+SpecialTau[ht]>PenalizedD[config]+eps)&&
	 (RecursiveSigma[defig]+SpecialSigma[ht]<
						PenalizedZ[config]-eps)) 
	];


(* Sample call CheckConfigWithSpecial @@ SpecialHexList *)
CheckConfigWithSpecial[configwithSpecial_]:= Module[
	{eps=10.0^-6,config,ht,i,specialList,tab},
	{config,specialList}=configwithSpecial;
	tab=Table[
		{defig,ht}=CutSpecial[config,specialList[[i]]];
		WithSpecial[config,defig,ht],
		{i,1,Length[specialList]}];
	Union[tab]=={True}
	];

FullCheck:= Module[{r6},
	Print[" -- Generating Lists -- "];
	QuadList;
	Print[" -- QuadList -- "];
	PentList;
	Print[" -- PentList -- "];
	HexList;
	Print[" -- HexList -- "];
	HeptList;
	Print[" -- HeptList -- "];
	Print[" -- Cases with no specials -- ",CheckAllTcc];
	Print[" -- Special Pent -- ",Union[CheckConfigWithSpecial /@ SpecialPentList]];
	Print[" -- Special Hex -- ",
		r6 = (CheckConfigWithSpecial /@ SpecialHexList);
		Union[r6],
		Position[r6,False]//Flatten];
	Print[" -- Special Hept -- ",Union[CheckConfigWithSpecial /@ SpecialHeptList]];
	Print[" -- finished! -- "]
	];

(* To do the verifications it is necessary to call FullCheck;
	It takes a few minutes to run.
	They all return "True" except for CheckConfigWithSpecial[SpecialHexList[[i]]],
		i = 175,176.  There are some addition notes in the paper
		about how to deal with these two cases. 
*)

(****************** DOUBLE SPECIALS ***********************)

(* Select[SpecialHeptList,Length[#[[2]]>1&], gives several isomorphic cases.
	SpecialSigma[1]=SpecialTau[1]=0.
	reduced7 = {{0,0,1,1,0},{2,0,0,0,2}};
	RecursiveTau[reduced7]== 0.75...
	RecursiveSigma[reduced7]== -0.557...
    not even close!

  Hexagon, two specials on opposite vertices.

	Table[
		RecursiveTau[{{0,0,0,0},{i,2,j,2}}]
		-PenalizedD[{"filler",{i,0,0,j,0,0}}],{i,0,2},{j,0,2}]//Min
	== 0.0202339 > 0.

	Table[
		RecursiveSigma[{{0,0,0,0},{i,2,j,2}}]+SpecialSigma[0] 2
        -PenalizedZ[{"filler",{i,0,0,j,0,0}}],{i,0,2},{j,0,2}]//Max
	== -0.0684666 <0.  

  Hexagon, two specials at nonopposite nonadjacent vertices.
	
	Table[
		RecursiveTau[{{0,0,j,0},{2,i,k,2}}]
		-PenalizedD[{"filler",{0,0,0,0,i,k}}],{i,0,2},{j,0,1},{k,0,2}]//Min
	== 0.02078 > 0.

	Table[
		RecursiveSigma[{{0,0,j,0},{2,i,k,2}}]+SpecialSigma[0] 2
        -PenalizedZ[{"filler",{0,0,0,0,i,k}}],{i,0,2},{j,0,1},{k,0,2}]//Max
	== -0.10842 <0.  

	So all the double special cases pass easily!
		

*)

(********************** LISTS ************************)

(* The arc-length of the perimeter *)
Arc[config_]:= Module[{i,v,e},
	{v,e}=config;
	Apply[Plus,Table[arcFixedLength[v[[i]],IPart[v,i+1],e[[i]]],
			{i,1,Length[v]}]]//Ns
	];
HasShortArc[config_]:=(Arc[config]<2.Pi+10.0^-6);

arcFixedLength[x__]:= arcFixedLength[x] = arc @@ Map[FixedLength,{x}];

(* Quad/Pent/Hex/HeptLists are the lists describing vertex/edge
	lengths of the possible convex polygons.  The format
	is described at the top of the file.  config={v,e}, v
	giving vertex heights, and e giving edge lengths.  *)

If[FirstPass,
HeptList:= HeptList = Module[
	{i1,i2,i3,i4,i5,i6,i7,twoTable,c1,c2,c3,c4,c5,c6},
	(* 4 arc[2.51,2.51,2]+3 arc[2.51,2.51,2.51] > 2.Pi
	   5 arc[2.51,2.51,2]+  arc[2.51,2.51,2.51]+arc[2.51,2.51,2sq]>2.Pi,
	  so there is at most two variables on top >2, if two both are 2.51.
	*)
	twoTable = Table[{i1,i2,i3,i4,i5,i6,i7}-1,
		{i1,1,2},{i2,1,2},{i3,1,2},{i4,1,2},{i5,1,2},{i6,1,2},{i7,1,2}]~
		Flatten~6;
	c1 = Map[{#,{0,0,0,0,0,0,0}}&,twoTable];
	c2 = Map[{#,{1,0,0,0,0,0,0}}&,twoTable];
	c3 = Map[{#,{2,0,0,0,0,0,0}}&,twoTable];
	c4 = Map[{#,{1,1,0,0,0,0,0}}&,twoTable];
	c5 = Map[{#,{1,0,1,0,0,0,0}}&,twoTable];
	c6 = Map[{#,{1,0,0,1,0,0,0}}&,twoTable];
	Select[Join[c1,c2,c3,c4,c5,c6],HasShortArc]
	]];

If[FirstPass,
HexList:=HexList = Module[{vlist,elist,i1,i2,i3,i4,i5,i6},
	(* If there are two consecutive 2.51s in the vlist,
		then either they are all 2.51s, or we can begin
		{2,2.51,2.51,...}.  This leads to the enumeration that follows.
	*)
	vlist = Join[{{1,1,1,1,1,1}},
		Table[{0,1,1,i1,i2,i3},{i1,0,1},{i2,0,1},{i3,0,1}]~Flatten~2,
		{{0,0,0,0,0,0},
		 {1,0,0,0,0,0},
		 {1,0,1,0,0,0},
		 {1,0,0,1,0,0},
		 {1,0,1,0,1,0}}];
	elist = Table[{i1,i2,i3,i4,i5,i6},
		{i1,0,2},{i2,0,2},{i3,0,2},{i4,0,2},{i5,0,2},{i6,0,2}]~
		Flatten~5;
	Select[Outer[{#1,#2}&,vlist,elist,1]~Flatten~1,HasShortArc]
	]];

If[FirstPass,
PentList:=PentList = Module[{vlist,elist,i1,i2,i3,i4,i5},
	(* If there are two consecutive 2.51s in the vlist,
		then either they are all 2.51s, or we can begin
		{2,2.51,2.51,...}.  This leads to the enumeration that follows.
	*)
	vlist = Join[{{1,1,1,1,1}},
		Table[{0,1,1,i1,i2},{i1,0,1},{i2,0,1}]~Flatten~1,
		{{0,0,0,0,0},
		 {1,0,0,0,0},
		 {1,0,1,0,0}
		 }];
	elist = Table[{i1,i2,i3,i4,i5},
		{i1,0,2},{i2,0,2},{i3,0,2},{i4,0,2},{i5,0,2}]~
		Flatten~4;
	Select[Outer[{#1,#2}&,vlist,elist,1]~Flatten~1,HasShortArc]
	]];

If[FirstPass,
QuadList:=QuadList = Module[{vlist,elist,i1,i2,i3,i4},
	vlist = {
		{0,0,0,0},
		{1,0,0,0},
		{1,1,0,0},
		{1,0,1,0},
		{1,1,1,0},
		{1,1,1,1}};
	elist = Table[{i1,i2,i3,i4},
		{i1,0,2},{i2,0,2},{i3,0,2},{i4,0,2}]~
		Flatten~3;
	Select[Outer[{#1,#2}&,vlist,elist,1]~Flatten~1,HasShortArc]
	]];

If[FirstPass, Module[{},
	SpecialPosition[config_]:= Module[{v,e,i},
		{v,e}=config;
		Position[Table[{IPart[e,i-1],IPart[e,i],IPart[v,i-1],IPart[v,i+1]},
				{i,1,Length[v]}],{0,0,0,0}]//Flatten
		];
	Specialize[list_]:= Module[{i},
		Select[
		Table[{list[[i]],SpecialPosition[list[[i]] ]},{i,1,Length[list]}],
		Length[#[[2]]]>0&]
		];
	SpecialQuadList:=SpecialQuadList=Specialize[QuadList];
	SpecialPentList:=SpecialPentList=Specialize[PentList]; 
	SpecialHexList:=SpecialHexList=Specialize[HexList]; 
	SpecialHeptList:=SpecialHeptList=Specialize[HeptList];
	];
];

FirstPassToken=0; 
