(* Division into simplices 0 represents an enclosed vertex. *)
MathToCplexExcept=1;
rE:= << MathToCplexExcept.m;

(* This files is only for the particular list
	ConfigurationList=(*SHORT/shortlist.m:*)shortlist;
*)


(* typically x_List is a subregion AllPent[[i,j]],
	The order puts the smallest first, and the next smallest
	among the adjacent second *)

(* DihOrder[{3,5,7,2,8}] == {2, 7, 5, 3, 8}; *)

DihOrder[cycle_List]:= Module[{t},
		If[Length[cycle]<3,Return[cycle]];
		t=MoveFirst[cycle,Min[cycle]];
		If[t[[2]]>Last[t],t=Reverse[t]];
		t=MoveFirst[t,Min[t]]
		];

(* Typically x_List is a decomposition AllPent[[i]] *)
DihListOrder[x_List]:= Sort[Map[DihOrder,x]];

(* subCyclic[3,4] == {1 -> 4, 2 -> 1, 3 -> 2, 4 -> 3}; *)
subCyclic[k_,n_]:= Array[#-> IMod[#+k,n]&,n];

(* subDih[3,4] == {1 -> 2, 2 -> 1, 3 -> 4, 4 -> 3}; *)
subDih[k_,n_]:= Table[i->IMod[k-i,n],{i,1,n}];

(* A full list of decompositions of the pentagon that we use,
	up to symmetry.  The order and orientation of the faces
	do not matter in this listing.
	vertex 0 = enclosed. *)

(* Expand PentDecompositions according to symmetries *)
AllPent:= AllPent =  Module[{i},
	PentDecompositions = {
		{{1,2,3,4,5}},
		{{1,2,5},{2,3,4,5}},
		{{1,2,5},{5,2,4},{4,2,3}},
		{{1,2,5},{5,2,0},{0,2,3},{0,3,4}, {0,4,5}},
		{{1,2,0,5},{0,2,3},{0,3,4},{0,4,5}},
		{{0,1,2},{0,2,3},{0,3,4},{0,4,5},{0,5,1}}
	};
	Map[DihListOrder,
		Flatten[Join[Table[PentDecompositions/.subCyclic[i,5],{i,1,5}],
			   Table[PentDecompositions/.subDih[i,5],{i,1,5}]],
 		1]
	]//Union
	];

(* a list of decompositions of the hexagon, up to dihedral
	symmetry,
	vertex 0 = enclosed.
*)

(* Expand HexDecompositions using symmetry *)
AllHex:= AllHex = Module[{},

	HexDecompositions = {
	{{1,2,3,4,5,6}},
	{{1,2,3,4,5},{5,6,1}},
	{{1,5,6},{1,2,4,5},{2,3,4}},
	{{1,5,6},{1,3,4,5},{1,2,3}},
	{{1,5,6},{1,2,3},{3,4,5},{1,3,5}},
	{{1,2,3,6},{3,4,5,6}},
	{{1,2,6},{2,3,6},{3,4,5,6}},
	{{1,5,6},{1,2,5},{2,4,5},{2,3,4}},
	{{1,5,6},{1,2,3},{1,3,4},{1,4,5}},

	{{0,5,6,1},{0,5,4},{0,4,3,2},{0,2,1}},
	{{0,1,6},{0,6,5,4},{0,4,3,2},{0,2,1}},
	{{1,5,6},{0,1,5},{0,5,4},{0,4,3,2},{0,2,1}},
	{{0,1,6},{0,6,5,4},{0,4,2},{0,2,1},{2,3,4}},
	{{1,5,6},{0,1,5},{0,5,4},{0,4,2},{0,2,1},{2,3,4}},
	{{0,1,6},{0,6,4},{0,4,2},{0,2,1},{6,5,4},{4,3,2}},

	{{0,6,3},{0,3,2},{0,2,1},{0,1,6},{3,4,5,6}},
	{{0,1,6},{0,6,5,4,3},{0,3,2},{0,2,1}},
	{{0,1,6},{0,6,3},{0,3,2},{0,2,1},{3,4,6},{5,4,6}},
	{{0,1,6},{0,6,4,3},{0,3,2},{0,2,1},{4,5,6}}
	};

	Map[DihListOrder,
		Flatten[Join[Table[HexDecompositions/.subCyclic[i,6],{i,1,6}],
			   Table[HexDecompositions/.subDih[i,6],{i,1,6}]],
 		1]
	]//Union
	];

(* A full list of heptagon decompositions mod symmetry.  
	All enclosed vertices
	have been erased on heptagons and octagons, but not on
	pentagons and hexagons *)

HeptDecompositions = {
{{1,2,3,4,5,6,7}},
{{1,2,3,4,5,7},{5,6,7}},
{{1,2,7},{2,3,5,6,7},{3,4,5}},
{{1,2,3,4,6},{4,5,6},{6,7,1}},
{{1,2,3},{3,4,5},{5,6,7},{1,3,5,7}}
};

(* Expand HeptDecompositions by symmetry *)
AllHept:= AllHept = 
	Map[DihListOrder,
		Flatten[Join[Table[HeptDecompositions/.subCyclic[i,7],{i,1,7}],
			   Table[HeptDecompositions/.subDih[i,7],{i,1,7}]],
 		1]
	]//Union;

(* A full list of octagon decompositions up to symmetry.
	All enclosed vertices have been erased *)

OctDecompositions = {
{{1,2,3,4,5,6,7,8}},
{{1,2,3,4,5,7,8},{5,6,7}},
{{1,2,3,5,7,8},{5,6,7},{3,4,5}},

{{1,2,4,5,7,8},{2,3,4},{5,6,7}},
{{1,3,4,5,7,8},{1,2,3},{5,6,7}},
{{1,3,5,7,8},{1,2,3},{3,4,5},{5,6,7}},

{{2,3,5,7,8},{8,1,2},{3,4,5},{5,6,7}},
{{1,3,5,7},{1,2,3},{3,4,5},{5,6,7},{7,8,1}}
};


(* Expanded form of OctDecompositions *)
AllOct:= AllOct = 
	Map[DihListOrder,
		Flatten[Join[Table[OctDecompositions/.subCyclic[i,8],{i,1,8}],
			   Table[OctDecompositions/.subDih[i,8],{i,1,8}]],
 		1]
	]//Union;

(* consistency: we get True from 
	Map[CheckConsistency,AllPent],...AllHex,AllHept,AllOct.. *)


CheckConsistency[x_]:= Module[{e,m,one,edgeCounts},
	e = Edges[x];
	edgeCounts = Map[Count[e,#]&,e];
	If[Max[edgeCounts]>2,Return[False]];
	one = Select[e,Count[e,#]==1&]//Sort;
	m = Max[x];
	If[Length[one]!= m,Return[False]];
	If[one!= Sort[Array[Sort[{#,AugMod[#,m]}]&,m]],
			Return[False]];
	True
	];

(* The variables expressing the internal structure of a exceptional region 
	pen is the penalty attached to the region.
	cquo is the "combined quoin" along an edge, it represents two quoins.
	cquo is expressed in terms of the numbering wrt a face-position.
	cquoI is expressed in terms of the numbering wrt vertex number
	sigE,solE,tauE are score,solid,squander of a subregion. The vertex_List
		is a list of the corners of the subregion.
	dihE is the dihedral angle at i, between edges extending to j and k
	Adih is the product A[h] dihE occuring in the formula for truncated Vor
	
*)

ListToString[x_List]:= S[Mod[x.Array[17^#&,Length[x]],(*prime*)882377]];
	(*  17 > 15  ==Max[Array[(Initialize[#]; Max[GBLregions])&,180]] *)
CheckListToString[f_,p_]:=Module[{},
	c = Map[ListToString,p/.VertexSub[f]];
	Length[Union[c]]==Length[p]
	];
	



StringVar["pen",f_]:= "pen"<>S[f]; 
StringVar["cquo",f_,p1_,p2_]:= Module[{i,j},
	i=GBLregions[[f,p1]];
	j=GBLregions[[f,p2]];
	StringVar["cquoI",f,i,j]
	];
StringVar["cquoI",f_,i_,j_]:= "cquo{"<>S[f]<>"}("<>S[Min[i,j]]<>","<>S[Max[i,j]]<>")";
StringVar["sigE",f_,vertex_List]:= StringVar["sig",f]<>
			ListToString[vertex//DihOrder];
StringVar["solE",f_,vertex_List]:= StringVar["sol",f]<>
			ListToString[vertex//DihOrder];
StringVar["tauE",f_,vertex_List]:= StringVar["tau",f]<>
			ListToString[vertex//DihOrder];
StringVar["dihE",f_,i_,j_,k_]:= "dih{"<>S[f]<>"}"<>ListToString[{i,Min[j,k],Max[j,k]}];
StringVar["AdihE",f_,i_,j_,k_]:= "A"<>StringVar["dihE",f,i,j,k];

(* VertexSub converts AllPent,AllHex,etc. to vertex numbering *)
VertexSub[f_]:= Array[#->GBLregions[[f,#]]&,Length[GBLregions[[f]]]];




(************** Section VI.4.3. VARIABLE RELATIONS ********************)

ExceptSolDihText[f_,p_]:= (* f is the face number, p = AllPent[[i]] *)
	Module[{AreaText,pp},
	AreaText[pp_]:= (* pp = AllPent[[i,j]], 
					or AllHex or whatever.
					Produce the text for the area of subregions in p *)
		Module[{pv,lP},
			pv = (pp/.VertexSub[f])//DihOrder;
			lP = pv//Length;
			{
				StringVar["solE",f,pv],
				Array[{" -",StringVar["dihE",f,Period[pv,#+1],pv[[#]],
						Period[pv,#+2]]}&,lP],
				" +",S[lP-2]," pi=0\n"
			}
		];
	Map[AreaText,p]//StringJoin
	];

ExceptDih2pi[f_,p_]:= (* f face number, p=AllPent[[i]], angles sum to 0 around
						enclosed *)
	Module[{pv},
		pv = (Select[p,Count[#,0]>0&]/.VertexSub[f])//DihListOrder;
		If[Length[pv]==0,Return[""]];
		{
		Map[{" -",StringVar["dihE",f,0,#[[2]],Last[#]]}&,pv],
		" + 2 pi = 0\n"
		}//StringJoin
	];

ExceptDihRelationText[f_,p_]:=
	(* now relate the dihedral angles at each vertex. *)
	Module[{pv,i},
	Table[
		pv = Map[MoveFirst[#,i]&,
				Select[p,Count[#,i]>0&]
				]/.VertexSub[f];
		{
		Map[{" -",StringVar["dihE",f,#[[1]],#[[2]],Last[#]]}&,pv],
					" + ",StringVar["dih",f,i]," = 0\n"
		},
	{i,1,Lg[GBLregions[[f]]]}]
	//StringJoin
	];

Penalty5[m_]:= If[m==5,0.008,0];

Penalty6[parameter_]:= Module[{X=0.01561,Y=0.003521}, 
	Switch[parameter,
	{0,0,0,0},2 0.008,
	{1,0,1,0},0.008,
	{0,3,4,1},3X,
	{0,3,3,0},3X,
	{1,3,6,0},X+2Y,
	{1,3,4,1},X+2Y,
	_,0]
	];

Penalty7[m_]:= Module[{X=0.01561,Y=0.003521},Switch[m,
	7,6X,
	6,5X,
	5,3X+2Y,
	4,X+4Y,
	_,0]
	];

Penalty8[m_]:= Module[{X=0.01561,Y=0.003521},Switch[m,
	8,6X,
	7,6X,
	6,4X+2Y,
	5,2X+4Y,
	_,0]
	];

Penalty[f_,p_]:= Module[{m},
	m = Max[Map[Length,p]];
	Switch[Max[p],
		5,Penalty5[m],
		6,Penalty6[HexParameters[p]], 
		7,Penalty7[m],
		8,Penalty8[m],
		_,0]
	];


(* HexParameters Input = AllHex[[n]]; *)

HexParameters[p_]:= 
	{
	FlatQuarters[p]//Length,
	UprightQuarters[p]//Length,
	Select[p,Length[#]==3&]//Length,
	Select[p,Length[#]==4&]//Length
	};

FlatQuarters[p_]:= Module[{r=Max[p]},Select[p,
	  Min[#]>0
	&&Length[#]==3
	&&Inside[#+1-Min[#],{{1,2,3},{1,r-1,r},{1,2,r}}]&]];

UprightQuarters[p_]:= Select[p,
	  Min[#]==0
	&&Length[#]==3
	&&Min[
		Mod[#[[3]]-#[[2]],Max[p]],
		Mod[#[[2]]-#[[3]],Max[p]]
		 ]==1&];

ExceptScoreText[f_,p_]:= Module[{pv},
	pv = (p/.VertexSub[f])//DihListOrder;
	{
		" -",StringVar["sig",f],
		Map[{" +",StringVar["sigE",f,#]}&,pv],
		" +",StringVar["pen",f]," >0\n",
	Map[{StringVar["sigE",f,#]," +",StringVar["tauE",f,#],
		" -0.1004445714270568 ",StringVar["solE",f,#]," =0\n"}&,pv]
	} 
	//StringJoin
	];

(* Use the table of penalties for the face in question, and
	the crude bound 0.01561 N penalty on other faces *)

ExceptPenaltyText[f_,p_]:= StringVar["pen",f]<>" ="<>S[Penalty[f,p]]<>"\n";
ExceptSecondaryPenaltyText[f_,p_]:= 
	StringVar["pen",f]<> " = "<>S[0.01561 Length[GBLregions[[f]]]]<>" \n";

(* Use this... *)
ExceptRelationText[f_,p_]:= 
	StringJoin[
		"\n\n\\ Relations among exceptional variables\n",
		ExceptSolDihText[f,p],
		ExceptDih2pi[f,p],
		ExceptDihRelationText[f,p],
		ExceptScoreText[f,p]
	];

(****************** END VARIABLE RELATIONS **********************)

(****************** BOUNDS ******************************)

ExceptBounds[f_,p_]:= Module[{pv},
	pv = (p/.VertexSub[f])//DihListOrder;
	Do[AppendTo[freeVar,StringVar["sigE",f,pv[[i]]]],{i,1,Length[pv]}];
	Do[AppendTo[freeVar,StringVar["tauE",f,pv[[i]]]],{i,1,Length[pv]}];
	];

(******************* Section VI.4.4. FLAT QUARTERS *************************)


FlatQuarterInequalities = Map[(#<>"\n")&,
    {
	"\n\n\\ Section VI.4.4. Flat Quarters, hat-sig scoring \n",
    " y4 > 2.51",
    " y4 < 2.8284271247462",
    "- dih2 + 0.35 y2 - 0.15 y1 - 0.15 y3 + 0.7022 y5 - 0.17 y4 > -0.0123",
    "- dih3 + 0.35 y3 - 0.15 y1 - 0.15 y2 + 0.7022 y6 - 0.17 y4 > -0.0123",
    "  dih2 - 0.13 y2 + 0.631 y1 + "<>
               "0.31 y3 - 0.58 y5 + 0.413 y4 + 0.025 y6 > 2.63363 ",
    "  dih3 - 0.13 y3 + 0.631 y1 + "<>
               "0.31 y2 - 0.58 y6 + 0.413 y4 + 0.025 y5 > 2.63363 ",
    " -dih1 + 0.714 y1 - 0.221 y2 - 0.221 y3 + "<>
               "0.92 y4 - 0.221 y5 - 0.221 y6 > 0.3482",
    "  dih1 - 0.315 y1 + 0.3972 y2 + 0.3972 y3 - "<>
               "0.715 y4 +  0.3972 y5 + 0.3972 y6 > 2.37095",
    "- solid - 0.187 y1 - 0.187 y2 - "<>
               "0.187 y3 + 0.1185 y4 + 0.479 y5 + 0.479 y6 > 0.437235 ",
    "+ solid + 0.488 y1 + 0.488 y2 + "<>
               "0.488 y3 - 0.334 y5 - 0.334 y6 > 2.244 ",
	"- sigma - 0.145 y1 - 0.081 y2 - 0.081 y3 - "<>
				"0.133 y5 - 0.133 y6 > -1.17401",
	"- sigma - 0.12 y1 - 0.081 y2 - 0.081 y3 - "<>
				"0.113 y5 - 0.113 y6 + 0.029 y4 > -0.94903", 
    " sigma + 0.153 y4 + 0.153 y5 + 0.153 y6 < 1.05382", 
	" sigma + 0.419351 solid + 0.19 y1 + 0.19 y2 + 0.19 y3 < 1.449", (* 2 NEW*)
(* constant fixed on 7/8/98, error introduced on 6/29 *)
	" sigma + 0.419351 solid -0.079431 dih1 -0.0436 y5 -0.0436 y6 < -0.01465",
	" sigma < 0.0114 ",
	" tau - 1.019 pt > 0 ",
	"\n\n"
    }
	]//StringJoin;


OneFlat[f_,{i_,j_,k_}]:= Module[{},
		 (* j=central vertex *)
		"\\ {i,j,k} = {"<>S[i]<>","<>S[j]<>","<>S[k]<>"}\n"<>
      StringReplace[FlatQuarterInequalities,
        {"y2"->StringVar["y",i],
         "y1"->StringVar["y",j],
         "y3"->StringVar["y",k],
         "y5"->StringVar["y",j,k],
         "y4"->StringVar["y",i,k],
         "y6"->StringVar["y",i,j],
         "dih2"->StringVar["dihE",f,i,j,k],
         "dih1"->StringVar["dihE",f,j,i,k],
         "dih3"->StringVar["dihE",f,k,i,j],
         "solid"->StringVar["solE",f,{i,j,k}],
         "tau"->StringVar["tauE",f,{i,j,k}],
		 "slack"->"slack{"<>S[f]<>"}",
         "sigma"->StringVar["sigE",f,{i,j,k}]}
    ]
    ];

				
ExceptFlatText[f_,p_]:= Module[{max,flatv},
	(* build list of flat quarters *)
	max = Max[p];
	flatv= (FlatQuarters[p]/.
		{{1,max-1,max}->{max-1,max,1}, {1,2,max}->{max,1,2}})
		/.VertexSub[f];
	{
		"\n\n\\ Vertex Sum < 0.114, at Type (4,1)-vertices VI.4.4.2  \n",
		Map[VertexSumType41[f,#]&,flatv],
		Map[OneFlat[f,#]&,flatv]
	}
		//StringJoin
	];


(* fv = {fv[[1]],fv[[2]],fv[[3]]}, fv[[2]] = central of flat *)

(* VI.4.4.2 *)
VertexSumType41[f_,fv_]:= Module[{pos},
	If[VertexType[fv[[2]]]!={4,0,1},Return[""]];
	pos = PSelect[GBLregions,Inside[fv[[2]],#]&]//Sort;
	{
	Map[{" +",StringVar["sig",#]}&,Drop[pos,-1]],
	" +",StringVar["sigE",f,fv],
	" < 0.114\n",
	Map[{" +",StringVar["sig",#]}&,Drop[pos,-1]],
	" +",StringVar["sigE",f,fv],
	" -",StringReplace[StringVar["sigE",f,fv],"sigma"->"slk41"], 
	" < 0.0875\n" (* slack 4.4.3*)
	} //StringJoin
	];

(***************** Section VI.4.5. UPRIGHT NU-QUARTERS ********************)

ExceptUprightEquations:= 
	Join[{
	" \\ Section VI.4.5. Upright Nu-Quarters in Exceptional Regions  ",
    " y1 > 2.51 ",
    " y1 < 2.8284271247462 ",
    " y2 > 2 ",
    " y3 > 2 ",
    " y4 > 2 ",
    " y5 > 2 ",
    " y6 > 2 ",
    " y2 < 2.51 ",
    " y3 < 2.51 ",
    " y4 < 2.51 ",
    " y5 < 2.51 ",
    " y6 < 2.51 ",
    " dih1 - 0.636 y1 + 0.462 y2 + 0.462 y3 - 0.82 y4 + 0.462 y5 + "<>
        " 0.462 y6 > 1.82419 ", 
    " - dih1 + 0.55 y1 - 0.214 y2 - 0.214 y3 + 1.24 y4 - 0.214 y5 "<>
        " - 0.214 y6 > 0.75281 ",  
    " dih2 + 0.4 y1 - 0.15 y2 + 0.09 y3 + 0.631 y4 - 0.57 y5 + 0.23 y6 " <>
    "  > 2.5481", 
    " - dih2 - 0.454 y1 + 0.34 y2 + 0.154 y3 - 0.346 y4 + " <>
        "0.805 y5 > -0.3429", 
    " dih3 + 0.4 y1 - 0.15 y3 + 0.09 y2 + 0.631 y4 - 0.57 y6 + 0.23 y5 " <>
    "  > 2.5481", 
    " - dih3 - 0.454 y1 + 0.34 y3 + 0.154 y2 - 0.346 y4 + " <>
        "0.805 y6 > -0.3429", 
    " sol + 0.065 y2 + 0.065 y3 + 0.061 y4 - 0.115 y5 - "<>
        "0.115 y6 > 0.2618", 
    " - sol - 0.293 y1 - 0.03 y2 - 0.03 y3 + 0.12 y4 + " <>
        "0.325 y5 + 0.325 y6 > 0.2514", 
	" -sig - 0.0538 y2 - 0.0538 y3 -0.083 y4 - 0.0538 y5 - "<>
		"0.0538 y6 > -0.5995",
	" sig < 0 ",
	" tau - 0.5945 pt > 0 "  ,
	(* Section VI.4.5.3 *)
	(* Part IV.A2  inequalities, 9052168, 746202672, etc. *)
	" sig -4.10113 dih1< -4.3223 ",
	" sig -0.80449 dih1< -0.9871 ",
	" sig -0.70186 dih1< -0.8756 ",
	" sig -0.24573 dih1< -0.3404 ", 
	" sig -0.00154 dih1< -0.0024 ",
	" sig +0.07611 dih1<  0.1196 ",
	(* Part IV.A3 inequalities, Section VI.4.5.3. *)
	" tau +4.16523 dih1>  4.42873 ",
	" tau +0.78701 dih1>  1.01104 ",
	" tau +0.77627 dih1>  0.99937 ",
	" tau +0.21916 dih1>  0.34877 ",
	" tau +0.05107 dih1>  0.11434 ",
	" tau -0.07106 dih1> -0.07749 "
    },ExceptUprightEquationsSlack
	];

							(* Section VI.4.5.4 *)
ExceptUprightEquationsSlack= 
	StringReplace[{
	" \\ Modified (these hold if diag is less than  2.696) ",
    " y1 < 2.696 ",
	"dih1 - 0.49 y1 + 0.44 y2 + 0.44 y3 - 0.82 y4 + 0.44 y5 + 0.44 y6 > 2.0421",
	"-dih1 + 0.495 y1 - 0.214 y2 - 0.214 y3 + 1.05 y4 - 0.214 y5 - "<>
		" 0.214 y6 > 0.2282 ", (* was 0.23545, changed 7/17/98 *)
	"-dih1 + 0.495 y1 - 0.214 y2 - 0.214 y3 + 1.05 y4 - 0.214 y5 - "<>
		" 0.214 y6 + oldA454 > 0.23545 ",  (* was this *)
	" dih2 + 0.38 y1 - 0.15 y2 + 0.09 y3 + 0.54 y4 - 0.57 y5 + 0.24 y6 > 2.3398",
	"-dih2 - 0.375 y1 + 0.33 y2 + 0.11 y3 - 0.36 y4 + 0.72 y5 +"<>
		"0.034 y6 > -0.36135",
	" dih3 + 0.38 y1 - 0.15 y3 + 0.09 y2 + 0.54 y4 - 0.57 y6 + 0.24 y5 > 2.3398",
	"-dih3 - 0.375 y1 + 0.33 y3 + 0.11 y2 - 0.36 y4 + 0.72 y6 +"<>
		"0.034 y5 > -0.36135",
	" sol+ 0.42 y1 + 0.165 y2 + 0.165 y3 - 0.06 y4 - 0.135 y5 -"<>
		"0.135 y6 > 1.479",
	"-sol - 0.265 y1 - 0.06 y2 - 0.06 y3 + 0.124 y4 + 0.296 y5 +"<>
		"0.296 y6 > 0.0997",
	"-sig + 0.112 y1 - 0.142 y2 - 0.142 y3 - 0.16 y4 - "<>
		"0.074 y5 - 0.074 y6 > -0.9029 ",
	" sig +0.07611 dih1<  0.11 ",  
	"-sig -0.015 y1 - 0.16 y2 -0.16 y3 -0.16 y4 -0.0738 y5 - 0.0738 y6 "<>
		" + slacG > -1.29285 ",
	" tau -0.07106 dih1> -0.06429 ", 
	" tau > 0.0414 "  
    },  {">"-> " +slack2696 > ",
		 "<"-> " -slack2696 < "}];


OneUprightText[f_,i_,j_]:= (* quad f, vertices i,j in 1..numvertices *)
		(* 0 = vertex 1, i = vertex 2, j = vertex 3 *)
    Module[{list},
		list=Apply[StringJoin,Map[(#<>"\n")&,ExceptUprightEquations]];
		"\n\n\\ {0,i,j} = {"<>S[0]<>","<>S[i]<>","<>S[j]<>"}\n"<>
      StringReplace[list,
        {"y1"->StringVar["y",10f+0],
         "y2"->StringVar["y",i],
         "y3"->StringVar["y",j],
         "y4"->StringVar["y",i,j],
         "y5"->StringVar["y",10f+0,j],
         "y6"->StringVar["y",10f+0,i],
         "dih1"->StringVar["dihE",f,0,i,j],
         "dih2"->StringVar["dihE",f,i,0,j],
         "dih3"->StringVar["dihE",f,j,0,i],
         "sol"->StringVar["solE",f,{0,i,j}],
         "sig"->StringVar["sigE",f,{0,i,j}],
         "tau"->StringVar["tauE",f,{0,i,j}],
		 "slack"->"slack{"<>S[f]<>"}",
		 "slacG"->"slkG{"<>S[f]<>"}"<>ListToString[{0,i,j}//DihOrder]
		}
    ]
    ];

ExceptUprightText[f_,p_]:= Module[{upv},
	(* build list of upright quarters *)
	upv = UprightQuarters[p]/.VertexSub[f];
	{
	"\n\n\\ Section VI.4.5. Upright nu-Quarters\n",
	Map[OneUprightText[f,#[[2]],#[[3]]]&,upv],
	ExceptUprightTextMore[f,p]
	}
	]//StringJoin;

									(* Section VI.4.5.5, VI.4.5.6  *)

ExceptUprightTextMore[f_,p_]:= Module[{p1,p2},
	If[Min[p]>0,Return[""]]; 
	p2 = Select[p,Min[#]==0&&Length[#]==4&]/.VertexSub[f]; 
	p1 = Complement[
				Select[p,Min[#]==0&&Length[#]==3&],
				UprightQuarters[p]
		]/.VertexSub[f];
	{
	"\n\\ Section VI.4.5.5: -0.05 bound\n",
	Map[{StringVar["sigE",f,#]," - slack{",S[f],"}2696 < -0.05\n"}&,p1],
	Map[{StringVar["sigE",f,#]," - slack{",S[f],"}2696 - "<>
		"slkI{"<>S[f]<>"}"<>ListToString[#//DihOrder]<>" < -0.119\n"}&,p1],
	Map[{StringVar["sigE",f,#]," < 0\n"}&,p1],

	"\n\n\\ Section VI.4.5.6: -0.043 bound \n",
	(* according to VI.4.6.11, the quad not covered by VI.4.5.6 satisfies
		a stronger bound, so it is not necessary to separate it out. *)
	Map[{StringVar["sigE",f,#]," - slack{",S[f],"}2696 < -0.043\n"}&,p2],
	Map[{StringVar["sigE",f,#]," - slack{",S[f],"}2696 - "<>
		" slkH{"<>S[f]<>"}"<>ListToString[#//DihOrder]  <>
		" < -0.091\n"}&,p2],
	Map[{StringVar["sigE",f,#]," < 0\n"}&,p2]
	}//StringJoin
	];

	
(******************* Section VI.4.6. Hexagonal Regions *******************)

HexagonText[f_,p_]:= Module[{sv,tv,i,psub}, 
	sv[i_]:= "+"<>StringVar["sigE",f,i];
	tv[i_]:= "+"<>StringVar["tauE",f,i];
	psub = (p/.VertexSub[f])//Sort; (*shorter expressions come first in Sort*)
	pr = (Select[p,Min[#]==0&]/.VertexSub[f])//Sort;
	"\n\\ Hex Inequalities  --  HexagonText (VI.4.6) \n"<>
	Switch[HexParameters[p],
	{0,0,0,0}, 	{ sv[psub[[1]]],"  < -0.212\n",
				tv[psub[[1]]]," > 0.54525\n" },
	{1,0,1,0}, { sv[psub//Last]," < -0.221\n",
				tv[psub//Last],"  > 0.486\n"},
	{2,0,2,1}, sv[psub//Last]<>" < -0.168 \n"<>
				tv[psub//Last]<>" >0.352\n",
	{0,0,0,2}, sv[psub[[1]]]<>" < -0.075 \n"<>
				tv[psub[[1]]]<>" > 0.176\n"<>
			   sv[psub[[2]]]<>" < -0.075 \n"<>
				tv[psub[[2]]]<>" >0.176\n",
	{1,0,2,1}, 	sv[psub//Last]<>" < -0.075 \n"<>
				tv[psub//Last]<>" >0.176\n",
	{0,2,2,2}, StringVar["sig",f]<>" < - 0.297\n"<>
			   StringVar["tau",f]<>" > 0.504\n",
	{1,2,4,1}, {Map[sv,pr]," < -0.253\n",
			    Map[tv,pr],"> 0.4686\n"
				},
	{2,2,6,0}, {Map[sv,pr]," < -0.2\n",
                Map[tv,pr] ,">0.3992\n"
                },
	{0,3,4,1}, sv[psub//Last]<>" < -0.075\n"<>
				tv[psub//Last]<>" >0.176\n",
	{0,3,3,0}, {Map[sv,pr]," < -0.2187\n",
                Map[tv,pr]," >0.518\n",
				sv[psub//Last]," < -0.137\n",
				tv[psub//Last]," > 0.31\n"
                } ,
	{1,3,4,1}, {Map[sv,pr ]," < -0.1657\n",
                Map[tv,pr ],"> 0.384\n",
				sv[psub//Last]," < -0.084\n",
				tv[psub//Last]," > 0.176\n"
                },
	_,""
		]<>HexagonDiagLength[f,p]
		]//StringJoin;

HexagonDiagLength[f_,p_]:= Module[{v,i,e},
	e = Edges[p];
	v= Intersection[e,{{1,4},{2,5},{3,6}}];
	v=v/.VertexSub[f];
	Table[StringVar["y",v[[i,1]],v[[i,2]] ]<>" > 2.51\n"<>
		  StringVar["y",v[[i,1]],v[[i,2]] ]<>" < 2.82842712474619\n",
			{i,1,Length[v]}]//StringJoin
		];



(*************** Section VI.4.7.  Pentagonal Regions ****************)

PentagonText[f_,p_]:= Module[{s},
	If[Max[p]!=5,Return[""]];
	{
	 If[Min[p]>0&&Length[p]==2,
		{
		"\n\\ Inequality VI.4.7.1 \n",
		StringVar["sigE",f,Last[Sort[p]]/.VertexSub[f]]<>" < -0.075\n",
		StretchQuadEqn[f,Last[Sort[p]]],
		StringVar["tauE",f,Last[Sort[p]]/.VertexSub[f]]<>" > 0.176\n"
		},
		""
		],
	 If[Min[p]==0&&Length[p//Last]==4,
		{
		"\n\\ Inequality VI.4.3 (quad region with uprights, sig<0) \n",
		StringVar["sigE",f,Last[p]/.VertexSub[f]]<>" <0\n"
		},
		""
		],
	If[Length[p]==1,
		{ "\n\\ Inequality VI.4.7.2 \n",
		StringVar["sigE",f,Last[p]/.VertexSub[f]]<>" < -0.128\n",
		StringVar["tauE",f,Last[p]/.VertexSub[f]]<>" > 0.36925\n"},
		""
		]
	}//StringJoin];

(* VI.4.7.1 and VI.4.6.4, a quad with one edge stretched past 2.51,
	both diagonals ge 2sqrt2. *)

StretchQuadEqn[f_,cycle_]:= Module[{eqn,cycf,g}, 
	If[Length[cycle]!=4,Return[""]];
	skip = Complement[Range[5],cycle];
	If[Length[skip]!=1,Return[""]];
	skip = skip//First;
	g = RotateLeft[GBLregions[[f]],skip-1];
	cycf = cycle/.VertexSub[f];
	StretchOne[f,cycf,g]
	]

	(* f is the usual face number,
	   cycf is the subregion.
	   g = {_,g2,g3,g4,g5}, with g2,g3,g4,g5 the vertices of the
				quad subregion, with the diagonal along (g2,g5)
	*)
StretchOne[f_,cycf_,g_]:= Module[{eqn},
	eqn = {
		" sig +0.1 y1+0.15 y2 +0.08 y3 +0.15 y5 +0.15 y6 +0.1 y7",
		" +0.17 y8 +0.16 y9 - slkJ.K1 < 2.1327","\n",
		" sig +0.1 y2+0.15 y1 +0.08 y7 +0.15 y9 +0.15 y6 +0.1 y3",
		" +0.17 y8 +0.16 y5 - slkJ.K2 < 2.1327","\n",
		" sig + 0.419351 sol - 0.0238 y5 - 0.0238 y6 -0.0238 y9",
		" < 0.4542 ","\n"
		};
	StringReplace[eqn,
		{
		"sig" -> StringVar["sigE",f,cycf],
		"sol" -> StringVar["solE",f,cycf],
		"slkJ"-> StringReplace[StringVar["sigE",f,cycf],"sigma"->"sJ"],
		"K1"-> S[g[[3]]],
		"K2"-> S[g[[4]]],
		"y1"-> StringVar["y",g[[3]]],
		"y2"-> StringVar["y",g[[4]]],
		"y3"-> StringVar["y",g[[2]]],
		"y5"-> StringVar["y",g[[2]],g[[3]]],
		"y6"-> StringVar["y",g[[3]],g[[4]]],
		"y7"-> StringVar["y",g[[5]]],
		"y8"-> StringVar["y",g[[2]],g[[5]]],
		"y9"-> StringVar["y",g[[4]],g[[5]]]
		}
		]
	];

	


(************************** EDGE LENGTH y(i,j) *******************)

(* The flat and upright sections have edge length inequalites for
edges on quarters.  We need edge length constraints 2-2.51 for
edges not on a quarter. By construction these are anchors not on a
quarter. This only arises in a few hexagonal cases (context (4.2)). *)

				
AnchorLengthText[f_,p_]:= Module[{v},
	v = Complement[
			Map[{#[[2]],Last[#]}&,Select[p,Min[#]==0&]]//Flatten,
			UprightQuarters[p]//Flatten
			]/.VertexSub[f];
	{"\n\\  AnchorLength\n",
	Map[StringVar["y",10f+0,# ]<>">2\n"<>
			StringVar["y",10f+0,# ]<>"<2.51\n"&,v]
	}//StringJoin
	]



(******************** Section VI.4.8. Dihederal Bounds ***********)


(*
	edge length codes:
		0 = short 2-2.51.
		1 = long, 2.51-2sq.
		2 = undrawn short,  2.51++
		3 = undrawn long,   2sq++
*)

EdgeCode[p_,v_]:= Module[{m,u,edge1,edge2,pedge},
	(* v = {a,b} subset {0,..,m} *)
	m = Max[p];
	u = Sort[v];
	edge1 = Array[Sort[{#,IMod[#+1,m]}]&,m];
	edge2 = Array[Sort[{#,IMod[#+2,m]}]&,m];
	pedge = Edges[p]//Union;
	If[Inside[u,edge1],Return[0] ];  (* boundary edge *)
	If[Inside[u,pedge]&&Min[u]==0,Return[0]];  (* anchor *)
	If[Inside[u,pedge],Return[1]];  (* drawn edge between corners *)
	If[Min[u]==0,Return[2]];  (* nonanchor *)
	If[Max[p]<7, Return[3]];  (* undrawn between corners, >2sq *)
	If[Inside[u,edge2],Return[3]]; (* no flat *)
	2
	];


NPi=N[Pi];
DihBound0[{y5_,y6_,y4_}]  := Switch[{y5,y6,y4},
	(* first section 5&6 short. *)
(*  y5,y6,y4,	min, max *)
	{0,0,1},	{1.153,2.28},
	{0,0,3},	{1.32,2NPi},
 
	(* second section 5 short, 6 long *)
	{0,1,0},	{0.633,1.624},
	{0,1,1},	{1.033,1.929},
	{0,1,2},	{1.033,2NPi},
	{0,1,3},	{1.259,2NPi},
 
	(* 5 long, 6 short *)
	{1,0,0},	{0.633,1.624},
	{1,0,1},	{1.033,1.929},
	{1,0,2},	{1.033,2NPi},
	{1,0,3},	{1.259,2NPi},
	
	(* 5 long, 6 long *)
	{1,1,0},	{0.817,1.507},
	{1,1,1},	{1.07,1.761},
	{1,1,2},	{1.07,2NPi},
	{1,1,3},	{1.23,2NPi},
	_,			{0,2NPi}
	];

DihBound1[{y5_,y6_,y4_}]:= Switch[{y5,y6,y4},
	(* upright diagonal at vertex 1, 5 drawn, 6 drawn *)
(*  y5,y6,y4,	min, max *)
	{0,0,0},	{0.956,2.184},
	{0,0,1},	{1.23,NPi},
	{0,0,2},	{1.23,NPi},
	{0,0,3},	{1.416,NPi},
	_,			{0,2NPi}
	];

DihBound2[{y5_,y6_,y4_}]:= Switch[{y5,y6,y4},
	(* upright diagonal at vertex 2, 6 drawn *)
(*  y5,y6,y4,	min, max *)
	{0,0,0},	{0.633,1.624},
	{0,0,2},	{1.033,2NPi},
	{1,0,0},	{0,1.381},
	{1,0,2},	{0.777,2NPi},
	_,			{0,2NPi}
	];

DihBound3[{y5_,y6_,y4_}]:= Switch[{y5,y6,y4},
	(* upright diagonal at vertex 3, 5 drawn *)
(*  y5,y6,y4,	min, max *)
	{0,0,0},	{0.633,1.624},
	{0,0,2},	{1.033,2NPi},
	{0,1,0},	{0,1.381},
	{0,1,2},	{0.777,2NPi},
	_,			{0,2NPi}
	];

DihBoundText[f_,p_]:= Module[
	{triple,c,tlist,v1,v2,v3,y4,y5,y6,i1,i2,i3,dih,switch,bds},
	triple[c_]:= Array[{Period[c,#],Period[c,#+1],Period[c,#+2]}&,Length[c]];
	tlist = Flatten[Map[triple,p],1]; 
	{
	"\n\n\\ Section VI.4.8. DihBoundText\n",
	Map[
		({v2,v1,v3} = #;
		{y5,y6,y4} = 
			{EdgeCode[p,{v1,v3}],EdgeCode[p,{v1,v2}],EdgeCode[p,{v2,v3}]};
		switch = DihBound0;
		If[v1==0,switch=DihBound1]; 
		If[v2==0,switch=DihBound2]; 
		If[v3==0,switch=DihBound3]; 
		bds = switch[{y5,y6,y4}];
		{i1,i2,i3}={v1,v2,v3}/.VertexSub[f];
		dih = StringVar["dihE",f,i1,i2,i3];
		{ 
		  dih," - 0.00001 dihSlack < ",S[bds[[2]]],"\n",
		  dih," + 0.00001 dihSlack > ",S[bds[[1]]],"\n"
		}
		)&,tlist]
	}//StringJoin
	];

		

(**************** Section VI.4.9. Additional Inequalities *******************)

(* 
DihEqn[{1,2,3,4,5,6},u]==
 +1 y1 +2 y2 +3 y3 +4 y4 +5 y5 +6 y6 -dih -0.000001 Buffer < - u
*)

DihEqn[coeff_,c_]:= 
	{Array[{SR[coeff[[#]]]," y",S[#]}&,6], " -dih",
	" -0.000001 Slack < ", S[-c ],"\n"}//StringJoin;

(* extra inequalities for dih,sigma depending on shape *)

AdditionalList0[{y5_,y6_,y4_}] := Switch[{y5,y6,y4},
	(* first section 5&6 short. *)
(*  {y5,y6,y4},		addtext, *)
	{0,0,1},		"", (* flat *)
	{0,0,3},
		"dih-0.372 y1 +0.465 y2 +0.465 y3 + 0.465 y5 + 0.465 y6 >4.885\n",
 
	(* second section 5 short, 6 long *)
	{0,1,0}, 		"",	(* flat *)
	{0,1,1}, DihEqn[{0.291,-0.393,-0.586,0.79,-0.321,-0.397},2.47477],
	{0,1,2},
			DihEqn[{0.291,-0.393,-0.586,0.0,-0.321,-0.397},4.45567],
	{0,1,3},
			{
			DihEqn[{0.291,-0.393,-0.586,0.0,-0.321,-0.397},4.71107],
			"dih -0.214 y1 +0.4 y2 +0.58 y3 +0.155 y5 +0.395 y6 > 4.52345\n"
			},
 
	(* 5 long, 6 short *)
	{1,0,0}, 		"", (* flat *)
	{1,0,1}, 
			DihEqn[{0.291,-0.586,-0.393,0.79,-0.397,-0.321},2.47477],
	{1,0,2}, 
			DihEqn[{0.291,-0.586,-0.393,0.0,-0.397,-0.321},4.45567],
	{1,0,3},
			{
			DihEqn[{0.291,-0.586,-0.393,0.0,-0.397,-0.321},4.71107],
			"dih -0.214 y1 +0.4 y3 +0.58 y2 +0.155 y6 +0.395 y5 > 4.52345\n"
			},
	
	(* 5 long, 6 long *)
	{1,1,0},
		{
		" tau  > 0.13943 \n",
		" sig  < -0.05714 \n",
		" -sol -0.492 y1 -0.492 y2 -0.492 y3 +0.43 y4 +0.038 y5+0.038 y6 < ",
			 S[-2.71884],"\n",
		" -sig -0.058 y1 -0.105 y2 -0.105 y3 -0.115 y4 -0.062 y5 -0.062 y6>",
			S[-1.02014],"\n",
		" sig+ 0.419351 sol < 0.3085","\n",
		DihEqn[{0.115,-0.452,-0.452,0.613,-0.15,-0.15},2.177]
		},
	{1,1,1},
		{
		DihEqn[{0.115,-0.452,-0.452,0.618,-0.15,-0.15},2.17382],
		" sig < -0.121\n", (* changed 7/14/98 *)
		" tau > 0.21301\n"
		},
	{1,1,2},
		DihEqn[{0.115,-0.452,-0.452,0.,-0.15,-0.15},3.725],
	{1,1,3},
		DihEqn[{0.115,-0.452,-0.452,0.,-0.15,-0.15},3.927],
	_,			""
	];

AdditionalList1[{y5_,y6_,y4_}] := Switch[{y5,y6,y4},
	(* upright diagonal, 5 drawn, 6 drawn *)
(*  {y5,y6,y4},		addtext, *)
	{0,0,0},	"", (* upright quarter *)
	{0,0,1},
		"sig < 0\n"<> 
		DihEqn[{0.47,-0.522,-0.522,0.812,-0.522,-0.522},2.82998],
	{0,0,2},
		DihEqn[{0.47,-0.522,-0.522,0.,-0.522,-0.522},4.8681],
	{0,0,3},
		DihEqn[{0.47,-0.522,-0.522,0.,-0.522,-0.522},5.1623],
	_,			""
	];

AdditionalList2[{y5_,y6_,y4_}] := Switch[{y5,y6,y4},
	(* upright diagonal at vertex 2, 6 drawn *)
(*  {y5,y6,y4},		addtext, *)
	{0,0,0},	"", (* upright quarter *)
	{0,0,2},
		" -0.4 y2 +0.15 y1 -0.09 y3 -0.631 y5-0.23 y6-dih < "<>
			" -3.9788\n",
	{1,0,0},
		DihEqn[{0.289,-1.36,-0.148,0.688,-1.36,-0.148},6.3282],
	{1,0,2},
		DihEqn[{0.289,-0.723,-0.148,0.,-0.723,-0.148},4.85746],
	_,			""
	];

AdditionalList3[{y5_,y6_,y4_}] := Switch[{y5,y6,y4},
	(* upright diagonal at vertex 3, 5 drawn *)
(*  {y5,y6,y4},		addtext, *)
	{0,0,0},"", (* upright quarter *)
	{0,0,2},
		" -0.4 y3 +0.15 y1 -0.09 y2 -0.631 y6-0.23 y5-dih < "<>
			" -3.9788\n",
	{0,1,0},
		DihEqn[{0.289,-0.148,-1.36,0.688,-0.148,-1.36},6.3282],
	{0,1,2},
		DihEqn[{0.289,-0.148,-0.723,0.,-0.148,-0.723},4.85746],
	_,			""
	];

AdditionalText[f_,p_]:= Module[
	{triple,c,tlist,v1,v2,v3,y4,y5,y6,i1,i2,i3,j1,j2,j3,switch,text},
	triple[c_]:= Array[{Period[c,#],Period[c,#+1],Period[c,#+2]}&,Length[c]];
	tlist = Flatten[Map[triple,p],1]; 
	{
	"\n\n\\ Section VI.4.9. Additional Inequalities\n",
	Map[
		({v2,v1,v3} = #;
		{y5,y6,y4} = 
			{EdgeCode[p,{v1,v3}],EdgeCode[p,{v1,v2}],EdgeCode[p,{v2,v3}]};
		switch = AdditionalList0;
		If[v1==0,switch=AdditionalList1]; 
		If[v2==0,switch=AdditionalList2]; 
		If[v3==0,switch=AdditionalList3]; 
		text = switch[{y5,y6,y4}];
		{i1,i2,i3}={v1,v2,v3}/.Join[VertexSub[f],{0-> (10f)}];
		{j1,j2,j3}={v1,v2,v3}/.VertexSub[f];
		cycle = Select[p,Subset[{v1,v2,v3},#]&]/.VertexSub[f];
		If[Length[cycle]!=1, Print["cycle error",cycle];Error["cycleX"]];
		StringReplace[text,
			{"y1"->StringVar["y",i1],
			 "y2"->StringVar["y",i2],
			 "y3"->StringVar["y",i3],
			 "y4"->StringVar["y",i2,i3],
			 "y5"->StringVar["y",i1,i3],
			 "y6"->StringVar["y",i1,i2],
			 "dih"->StringVar["dihE",f,j1,j2,j3],
			 "Slack"->"dihSlack",
			 "sig"->StringVar["sigE",f,First[cycle]],
			 "sol"->StringVar["solE",f,First[cycle]],
			 "tau"->StringVar["tauE",f,First[cycle]]
			}]
		)&,tlist]
	}//StringJoin
	];


(******* Section VI.4.10. Miscellaneous Inequalities (VC EQUATIONS) ********)


(* Formula for VorVc as a function of Adih, quoins, etc.*)
VcDefText[f_,p_]:= Module[{x,fdoct,phi0,i,AdihVertex,v,quoinEdge,boundaryEdge},
	x = Complement[p,FlatQuarters[p]~Union~UprightQuarters[p]];  
	x = (x/.VertexSub[f])//DihListOrder;
	boundaryEdge= Module[{face},
		face = GBLregions[[f]];
		Array[Sort[{Period[face,#],Period[face,#+1]}]&,Length[face]]
		];
	fdoct = 2.883611798069859;
	phi0 = -0.5666365478933329;
	{"\n\n\\ Truncated vor expansion (VI.4.10) : \n",
	 Table[
		v = x[[i]];
		AdihVertex = 
		  Array[{Period[v,#+1],Period[v,#],Period[v,#+2]}&,Length[v]];
		AdihVertex = Select[AdihVertex,#[[1]]>0&];
		quoinEdge= Intersection[ boundaryEdge,
			Array[Sort[{v[[#]],Period[v,#+1]}]&,Length[v]]
			];
		{StringVar["sigE",f,v],
		 SR[-phi0], " ",StringVar["solE",f,v], 
		 Map[{"\n",SR[fdoct]," ",StringVar["cquoI",f,#[[1]],#[[2]]]}&,quoinEdge],
		"\n",
		 Map[{" - ",StringVar["AdihE",f,#[[1]],#[[2]],#[[3]]]}&,AdihVertex],
		" < 0 \n"}
 
		,{i,1,Length[x]}
		]}//StringJoin
	]; 


UpdateAdihText[basefile_,f_,p_]:= If[basefile=="","",
	(* else *)
	Map[UpdateSubregion[basefile,f,#]&,
		Complement[p,UprightQuarters[p]~Union~FlatQuarters[p]]
		]
	  ];

UpdateSubregion[basefile_,f_,c_]:= 
	UpdateSubregion[basefile,f,c,""];

UpdateSubregion[basefile_,f_,c_,add_]:= Module[{v,AdihVertex},
	v = (c/.VertexSub[f])//DihOrder;
	AdihVertex = 
	  Array[{Period[v,#+1],Period[v,#],Period[v,#+2]}&,Length[v]];
	AdihVertex = Select[AdihVertex,#[[1]]>0&];
	{
	"\n\n\\ UpdateSubegion ",S[v],"\n",
	Map[  LPmUpdateAdih[basefile,StringVar["y",#[[1]]],
		 StringVar["dihE",f,#[[1]],#[[2]],#[[3]]],
		 StringVar["AdihE",f,#[[1]],#[[2]],#[[3]]],
			add
				] &,AdihVertex
	   ]
	}//StringJoin
	];

(********** Section VI.4.10.1. vertexAdjustment (1.4 and 1.5 excesses) *******)

(* if there is a flat whose first edge is at a vertex with excess,
	then the excess can be formed with the flat and the surrounding
	regions, rather than drawing in the entire face.
	
	If the fourth edge has length greater than 2sq, then the excess
	can be formed using nothing from the face. *)

VertexAdjustmentText[f_,p_]:= 
	Module[{eVertices,centeredFlat,uncenteredFlat,str},
    eVertices = Module[{vertexRange,vertex401,vertex311,vertex302},
        vertexRange = GBLregions[[f]]//Sort;
        vertex401 = Select[vertexRange,(VertexType[#]=={4,0,1})&];
        vertex311 = Select[vertexRange,(VertexType[#]=={3,1,1})&];
        vertex302 = Select[vertexRange,(VertexType[#]=={3,0,2})&];
        Union[vertex401,vertex311,vertex302]];
    (*build flats centered on eVertices *)
	{centeredFlat,uncenteredFlat} = 
	   Module[{pr,x,max,triples,reorder,flatv,flatr,flats,centered,lenG3},
		x = Select[p,Count[#,0]==0&];
		x = Select[x,Length[#]==3&];
		max = Max[p];
		triples = Table[{i,IMod[i+1,max],IMod[i+2,max]},{i,1,max}]//DihListOrder;
		x = Intersection[triples,x];
		reorder[k_]:= (* place first edge as central variable *)
			Switch[k,
				{1,max-1,max},{max-1,max,1},
				{1,2,max},{max,1,2},_,k];
		x = Map[reorder,x];
		flatv= x/.VertexSub[f];
		centered= Intersection[Map[#[[2]]&,flatv],eVertices];
		flatr= Intersection[Map[{#[[1]],#[[3]]}&,flatv]//Flatten,
				eVertices];
		(* Add open vertices to uncentered *)
		pr = p/.VertexSub[f];
		lenG3 = Select[pr,Length[#]>3&]//Flatten//Union;
		flats=Select[eVertices,Inside[#,lenG3]&&Count[pr//Flatten,#]==1&];
		{centered,Union[flatr,flats]}
		];
	(* Print[p/.VertexSub[f]," ",centeredFlat," ",uncenteredFlat];   *)
 
 
    (*build uncentered string*)
	ustr = {
    "\n\n\\ Vertex Adjustment Text\n\n",
	Module[{i,vertex,faces},
        Table[
			vertex = uncenteredFlat[[i]];
			faces = Select[Range[Length[GBLregions]],
				Length[Intersection[GBLregions[[#]],{vertex}]]>0&];
			faces=Complement[faces,{f}];
		{
        Map[{" +",StringVar["tau",#]}&,faces],
			 " > ", S[vSquander1[vertex]+fSquander[faces]],"\n"
		}
        ,{i,1,Length[uncenteredFlat]}
        ]]
		};
 
    (*build centered string*)
    str = {
	"\n\n\\ Vertex Centered Adjustment Text\n\n",
	Module[{i,vertex,faces,C,mf},
		C =     0.06585 ; (* D(3,1) *)
        Table[
			vertex = centeredFlat[[i]];
			faces = Select[Range[Length[GBLregions]],
				Length[Intersection[GBLregions[[#]],{vertex}]]>0&];
			faces=Complement[faces,{f}];
			mf = MoveFirst[GBLregions[[f]],vertex];
			mf = mf[[{1,2,Length[mf]}]]//DihOrder;
		{
        Map[{" +",StringVar["tau",#]}&,faces],
		 " +", StringVar["tauE",f,mf],"\n",
        " > ", S[C+vSquander1[vertex]+fSquander[faces]], "\n"
		}
        ,{i,1,Length[centeredFlat]}
        ]]
		};
 
    {str,ustr}//StringJoin
    ];


(******** Section VI.4.10.2, VI.4.10.3 Edge Distortion Inequalities *****)


(* If there are four qrtets and something else at a vertex.
	The ten edges around the common edge have total length at least const.*)
EdgeDistortionOneText[f_,v_,const_]:= 
	Module[{s,vertices,vars},
	s=PSelect[GBLregions,Inside[v,#]&];
	s=Complement[s,{f}];
	If[Length[s]!=4
		|| Union[Map[Length,GBLregions[[s]] ]]!={3},
		Return[""]];
	vertices = Complement[GBLregions[[s]]//Flatten//Union,{v}];
	vars = Map[StringVar["y",#]&,vertices]
			~Join~
			Map[StringVar["y",v,#]&,vertices];
	If[Length[vars]!=10,Error["EdgeDistortionError"];Return[""]];
	{Map[{"+",#}&,vars]," > ",S[const],"\n"}//StringJoin
	];

EdgeDistortion251:= Module[{fc,j},
	fc=PSelect[GBLregions,Length[#]>3&];
	{
	"\n\n\\ Edge Distortion 2.51 (See Sphere Packings VI.4.10.2) \n",
	Map[
		Table[EdgeDistortionOneText[#,GBLregions[[# ,j]],20.42]
		,{j,1,Lg[GBLregions[[#]] ]}]&,fc
		]
	}//StringJoin
	]


EdgeDistortionSqrt2[f_,p_]:= Module[{i,pc,pcf,v,s,str},
	pc = p/.VertexSub[f];
	pcf = Complement[Flatten[Union[pc]],{0}];
	If[Length[GBLregions[[f]]]<5,Return[""]];
	str="\n\n\\ EdgeDistortionSqrt2 (See Sphere Packings VI.4.10.3)\n";
	Do[
		v = pcf[[i]];
		s=Select[pc,Inside[v,#]&];
		(* a nontriangular region has diag greater than 2sqrt2 *)
		If[Length[s]==1&&Length[s[[1]]]>3,
			str=str<>EdgeDistortionOneText[f,v,20.76]
			];
		(* a vertex with a cross edge has diag > 2sqrt2 *)
		If[Length[s]>1&&Min[s]>0,
			str=str<>EdgeDistortionOneText[f,v,20.76]
			];
		,{i,1,Length[pcf]}
		];
	str
	];

SnText[f_,p_]:= Module[{h,k,v},
	h = "\n\n \\ Sn Text\n";
	If[Max[p]<7,Return[""]];
	k = Length[FlatQuarters[p]];
	v = Last[p]/.VertexSub[f];
	Switch[Max[p],
		8,
			{
			h<>StringVar["sig",f]<>" < -0.22816\n",
			StringVar["sigE",f,v]," < -0.22816\n",
			StringVar["tauE",f,v]," > ",S[0.6045-k 0.06585]
			},
		7,	{
			h<>StringVar["sig",f]<>" < -0.17112\n",
			StringVar["sigE",f,v]," < -0.17112\n",
			StringVar["tauE",f,v]," > ",S[0.54999-k 0.06585]
			},
		_,""]//StringJoin
	];
		

(**************** PUT IT ALL TOGETHER (OUTPUT) *********************)

Setp[f_,ppos_]:= Module[{p},
	p = Switch[GBLregions[[f]]//Length,
		5,AllPent[[ppos]],
		6,AllHex[[ppos]],
		7,AllHept[[ppos]],
		8,AllOct[[ppos]],
		_,Print["Setp out of range!"]; {}
		];
	If[!CheckListToString[f,p],Error["Setp Error ",i," ",f," ",ppos]];
	p
	];

fileText[stem_,n_,i_,f_,ppos_]:= stem<>S[n]<>"."<>S[i]<>".F"<>S[f]<>".C"<>S[ppos];

WRITEOUTtwoPhase[i_,f_,ppos_]:= Module[{exceptBak,basefile},
	exceptBak=exceptStem;
	exceptStem = "/tmp/cplexE.lp";
	WRITEOUTexcept[i,f,ppos];
	basefile = exceptStem<>S[i]<>".F"<>S[f]<>".C"<>S[ppos];
	exceptStem= exceptBak;
	WRITEOUTexcept[basefile,i,f,ppos];
	];

Initialize[n_,config_]:= Module[{list,arrang},
		ConfigurationList = Switch[n,
			5,lpent,
			6,lhex,
			7,lhept,
			8,loct,
			_,{}
			];
        list = ConfigurationList[[config,3]];
        arrang = Map[prune,list ];
        GBLconfig = config;
        GBLregions  = StdRegions[arrang];
        ];

exceptStem ="SHORT/cplexE.lp";
WRITEOUTexcept[n_,i_,f_,ppos_]:= WRITEOUTexcept["",n,i,f,ppos];
(* basefile is used for LPbounds used in constructing the file *)
WRITEOUTexcept[basefile_String,n_,i_,f_,ppos_]:= Module[{stream,p,file},
	Initialize[n,i];
	freeVar={"X","sigsum"}; (* global *)
	file=exceptStem<>S[n]<>"."<>S[i]<>".F"<>S[f]<>".C"<>S[ppos];
	stream=OpenWrite[file];
	streamBAK=stream;
	WriteString[stream,"\n\\ basefile = "<>basefile<>"\n"];
	WRITEOUTstd[stream,basefile];
	p = Setp[f,ppos]; 
	WriteString[stream,AnchorLengthText[f,p]];

	(* Section VI.4.3 *)
	WriteString[stream,ExceptRelationText[f,p]];
	WriteString[stream,ExceptPenaltyText[f,p]];
	ExceptBounds[f,p];

	(* Section VI.4.4, VI.4.5 *)
	WriteString[stream,ExceptFlatText[f,p]];
	WriteString[stream,ExceptUprightText[f,p]];

	(* Section VI.4.6, VI.4.7 *)
	If[Max[p]==6, WriteString[stream,HexagonText[f,p]]]; 
	If[Max[p]==5, WriteString[stream,PentagonText[f,p]]];

	(* Section VI.4.8, VI.4.9 *)
	WriteString[stream,DihBoundText[f,p]];
	WriteString[stream,AdditionalText[f,p]];

	(* Section VI.4.10 *)
	WriteString[stream,QuoinText[f]];
	WriteString[stream,VcDefText[f,p]];
	WriteString[stream,UpdateAdihText[basefile,f,p]];
	WriteString[stream,SnText[f,p]];  (* 4.10.0 *)
	WriteString[stream,VertexAdjustmentText[f,p]]; (* 4.10.1 *)
	WriteString[stream,EdgeDistortion251]; (* 4.10.2 *)
	WriteString[stream,EdgeDistortionSqrt2[f,p]]; (* 4.10.3 *)

	(* wrapup *)
	WriteString[stream,bounds];
	WriteString[stream,"\n\nEND\n\n"];
	WriteString[stream,faceCode];
	Close[stream];
	If[Length[ErrorLog]>0,Print["Warning: there are errors"]];
	file
	];

(******************  BRANCH AND BOUND STUFF ***********************)





