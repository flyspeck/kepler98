rP:= << MathToCplexPent.m;

filePentText[stem_,n_, i_,fplist_] := 
	{
	stem,S[n],".",S[i],
	Map[{".F",S[#[[1]]],".C",S[#[[2]]]}&,fplist]
	}//StringJoin;



(* Based on WRITEOUTexcept *)
(* For example the pent-hex case, WRITEOUTmany[167,{{15,68},{14,1}}]; *)
WRITEOUTmany[n_,i_,fplist_]:= Module[{file,stream},
	Initialize[n,i];
    freeVar={"X","sigsum"}; (* global *)
	file = {
			exceptStem,S[n],".",S[i],
			Map[{".F",S[#[[1]]],".C",S[#[[2]]]}&,fplist]
			}//StringJoin;
	Print[file];
    stream=OpenWrite[file];
	streamBAK=stream;
	WRITEOUTstd[stream];
	Map[WRITEOUTonePart[stream,#[[1]],#[[2]]]&,fplist];
	WriteString[stream,bounds];
    WriteString[stream,"\n\nEND\n\n"];
    WriteString[stream,faceCode];
    Close[stream];
    If[Length[ErrorLog]>0,Print["Warning: there are errors"]];
	file
    ];

(* based on WRITEOUTexcept *)
WRITEOUTonePart[stream_,f_,ppos_]:= Module[{p},
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
    (* WriteString[stream,UpdateAdihText[basefile,f,p]]; *)
    WriteString[stream,SnText[f,p]];  (* 4.10.0 *)
    WriteString[stream,VertexAdjustmentText[f,p]]; (* 4.10.1 *)
    WriteString[stream,EdgeDistortion251]; (* 4.10.2 *)
    WriteString[stream,EdgeDistortionSqrt2[f,p]]; (* 4.10.3 *)
	];


LiftCandidates[p_]:=
		Map[{#[[2]],Last[#]}&,UprightQuarters[p]]//Flatten//Union;

(* Lift vertex GBLregions[[f,i]] to 2.45 and see what happens to pattern p *)
PentLiftAdjust[f_,p_,i_]:= Module[{v,guys,ups},
	If[Min[p]>0,Return[""]];
	If[!Inside[i,Flatten[p]],Return[""]];
	v = i/.VertexSub[f];
	guys = Select[p,Inside[i,#]&];
	ups = UprightQuarters[p]/.VertexSub[f];
	{
	StringVar["y",10f]," > 2.696\n",
	StringVar["y",v,10f] ," > 2.45\n",
	StringVar["y",v], " > 2.45\n",
	Map[ Switch[True,
	Length[#]==4,
		{
		StringVar["sigE",f,#]," < -0.136\n",
		StringVar["tauE",f,#]," > 0.224\n"
		},
	Min[#]>0,
		{
		StringVar["sigE",f,#]," < -0.039\n",
		StringVar["tauE",f,#]," > 0.094\n"
		},
	Inside[#,ups],
		{
		StringVar["sigE",f,#]," < -0.055\n",
		StringVar["tauE",f,#]," > 0.092\n"
		},
	True,	(* Else anchored simplex *)
		{
		StringVar["sigE",f,#]," < -0.197\n",
		StringVar["tauE",f,#]," > 0.239\n"
		}
		] (* end of Switch *)
		&,
		guys/.VertexSub[f]
		] (* end of Map *)
	}//StringJoin
	];

LiftDoublet[i_,p_,q_]:=Module[{bf,f},
	Initialize[5,i];
	f = GBLregions//Length;
    bf = "SHORT/PENT/LPmany/cplexE.lp5."<>S[i]<>".F"<>S[f-1]<>".C"<>
			S[p]<>".F"<>S[f]<>".C"<>S[q];
        {
        Map[LPdisplay[bf,PentLiftAdjust[f-1,AllPent[[p]],#]]&,
                LiftCandidates[AllPent[[p]]]
            ],
        Map[LPdisplay[bf,PentLiftAdjust[f,AllPent[[q]],#]]&,
                LiftCandidates[AllPent[[q]]]
            ] 
        }
    ];

(* These are the cases I'm interested in running on 6/28/98 *)
PentDoublets =
	{895, 947, 949, 951, 953, 957, 997, 1894, 1983, 2016, 2246};
Infeasible=0;
If[Length[DoubletValues]<5,<< SHORT/PENT/cplexMany.sum; DoubletValues=values];
pick[i_]:= Select[DoubletValues,#[[3]]==i&&#[[1]]>0.4429&];


DisplayLiftDoublet[i_]:= Module[{pairs,pairLift},
	Initialize[5,i];
	pairs = Map[{#[[5]],#[[7]]}&,pick[i]];
	pairLift[i]=Map[LiftDoublet[i,#[[1]],#[[2]]]&,pairs];
	Definition[pairLift] >>> SHORT/PENT/displayDoublet.m;
	Max[pairLift[i]]
	];

(* Map[DisplayLiftDoublet,PentDoublets]==
  {0.4558582787999999, 0.39383321374, 0.3219080457399999, 0.3416384307199999, 
   0.4007247711599999, 0.43097382691, 0.2996838499999999, 0.25913797899, 
   0.29771575471, 0.3308552866199999, 0.39002672923};
	PSelect[%,#>0.4429&] == {1},
	PentDoublets[[1]] == 895, 
	so in all others with have both slacks{}2696=0
	 *)


(***************************** SINGLETS ********************)

(* I'm particularly interested in these cases today 6/28/98 *)
PentSinglets= {663,1454,1590};

LiftSinglet[i_,p_]:=Module[{bf,f},
	Initialize[5,i];
	f = GBLregions//Length;
    bf = "SHORT/PENT/LP/cplexE.lp5."<>S[i]<>".F"<>S[f]<>".C"<>S[p];
	Map[LPdisplay[bf,PentLiftAdjust[f,AllPent[[p]],#]]&,
			LiftCandidates[AllPent[[p]]]
		]
    ];
Infeasible=0;
(* This is for 1454 and 1590, 663 will be based on cplexE.sum *)
If[Length[SingletValues]<5,<< SHORT/PENT/cplexQ.sum; SingletValues=values];
pickSinglet[i_]:= Select[SingletValues,#[[3]]==i&&#[[1]]>0.4429&];


DisplayLiftSinglet[i_]:= Module[{singles,singleLift},
	Initialize[5,i];
	singles = Map[#[[5]]&,pickSinglet[i]]//Union;
	singleLift[i]=Map[LiftSinglet[i,#]&,singles];
	Definition[singleLift] >>> SHORT/PENT/displaySinglet.m;
	Max[Flatten[singleLift[i]]]
	];

(* DisplayLiftSinglet[1454]; DisplayLiftSinglet[1590]; *)
(* Now for 663...
	<< SHORT/PENT/cplexE.sum; 
	set663= Select[values,#[[3]]==i&&#[[1]]>0.4429&];
	singles = Map[#[[5]]&,set663]//Union;
	Initialize[5,663];
	Clear[singleLift];
	singleLift[663]=Map[LiftSinglet[663,#]&,singles];
	Definition[singleLift] >>> SHORT/PENT/displaySinglet.m;
	Max[Flatten[singleLift[663]]] == 0.317562;
*)


(********************* HEXAGONAL 2.696 UPRIGHT DIAGONAL TESTS *************)


(* Lift vertex GBLregions[[f,i]] to 2.45 and see what happens to pattern p *)
HexLiftAdjust[f_,p_,i_]:= Module[{v,guys,ups,sigVar,tauVar,flats},
	If[Min[p]>0,Return[""]];
	If[!Inside[i,Flatten[p]],Return[""]];
	v = i/.VertexSub[f];
	guys = Select[p,Inside[i,#]&];
	ups = UprightQuarters[p]/.VertexSub[f];
	flats = FlatQuarters[p]/.VertexSub[f];
	{
	StringVar["y",10f]," > 2.696\n",
	StringVar["y",v,10f] ," > 2.45\n",
	StringVar["y",v], " > 2.45\n",
	Map[ 
	(sigVar = StringVar["sigE",f,#];
	tauVar = StringVar["tauE",f,#];
	Switch[True,
	Length[#]==5&&Min[#]==0,
		{
		sigVar," < -0.254\n",
		tauVar," > 0.42625\n"
		},
	Length[#]==4&&Min[#]>0,
		{
		sigVar," < -0.149\n",
		tauVar," > 0.281\n"
		},
	Length[#]==3&&Min[#]>0&&(!Inside[#,flats]),
		{
		sigVar," < -0.089\n",
		tauVar," > 0.154\n"
		},
	Inside[#,flats],
		{
		sigVar," < -0.039\n",
		tauVar," > 0.094\n"
		},
	Length[#]==4&&Min[#]==0&&Length[ups]==3,
		{
		sigVar," < -0.24\n",
		tauVar," > 0.346\n"
		},
	Length[#]==4&&Min[#]==0&&Length[ups]==2,
		{
		sigVar," < -0.136\n",
		tauVar," > 0.224\n"
		},
	Length[#]==3&&Min[#]==0&&(!Inside[#,ups]),
		{
		sigVar," < -0.197\n",
		tauVar," > 0.239\n"
		},
	Inside[#,ups],
		{
		sigVar," < -0.055\n",
		tauVar," > 0.092\n"
		},
	True,	(* Else anchored simplex *)
		{
		"\\ unexpected default. \n"
		}
		] (* end of Switch *)
		)&,
		guys/.VertexSub[f]
		] (* end of Map *)
	}//StringJoin
	];

LiftHexSinglet[i_,ppos_]:=Module[{bf,f},
	Initialize[6,i];
	f = GBLregions//Length;
    bf = "SHORT/HEX/LP/cplexE.lp6."<>S[i]<>".F"<>S[f]<>".C"<>S[ppos];
	Map[LPdisplay[bf,HexLiftAdjust[f,AllHex[[ppos]],#]]&,
			LiftCandidates[AllHex[[ppos]]]
		]
    ];
Infeasible=0;

(* see SHORT/HEX/index.html for derivation of upPair *)

upPair = {{59, 35}, {59, 39}, {59, 61}, {59, 73}, {59, 79}, {59, 84}, 
   {59, 98}, {59, 103}, {59, 110}, {70, 64}, {70, 94}, {104, 61}, {104, 79}, 
   {111, 88}, {111, 96}, {131, 33}, {131, 41}, {131, 68}, {131, 71}, 
   {131, 85}, {131, 88}, {131, 105}, {131, 112}, {131, 114}, {226, 66}, 
   {226, 68}, {226, 89}, {226, 94}, {248, 32}, {248, 44}, {248, 64}, 
   {248, 66}, {248, 94}, {248, 96}, {248, 102}, {248, 117}, {296, 64}, 
   {296, 66}, {296, 94}, {296, 102}, {296, 116}, {302, 61}, {310, 64}, 
   {310, 116}, {363, 68}, {363, 85}, {363, 112}, {385, 94}}

DisplayHexLiftSinglet[i_]:= Module[{singles,singleLift},
	singles = Select[upPair,#[[1]]==i&];
	singleLift[i]=Map[LiftHexSinglet[#[[1]],#[[2]]]&,singles];
	Definition[singleLift] >>> SHORT/HEX/displayHexSinglet.m;
	Max[Flatten[singleLift[i]]]
	];

