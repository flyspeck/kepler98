(* Generate equations for CPLEX from final.m *)
(* 
	$Id: MathToCplex.m,v 1.26 1998/07/31 05:58:57 hales Exp hales $
	$Log: MathToCplex.m,v $
	Revision 1.26  1998/07/31 05:58:57  hales
	*** empty log message ***

	Revision 1.25  1998/06/25 00:51:32  hales
	Another bug in VertexExcessText fixed.

	Revision 1.24  1998/06/23 23:53:02  hales
	*** empty log message ***

	Revision 1.23  1998/06/23 23:46:28  hales
	bug in VertexExcessText fixed.

	Revision 1.22  1998/06/23 14:01:49  hales
	*** empty log message ***

	Revision 1.21  1998/06/04 21:25:36  hales
	Constant in SquanderFace corrected.  Doesn't affect anything
	I've run with this version of MathToCplex.m

	Revision 1.20  1998/06/02 17:52:06  hales
	Version used in 6/2/98 version of Sphere Packings III.

	Revision 1.19  1998/05/27 02:45:24  hales
	Another major revision.  MathToCplexQuadInExcept.m and
	MathToCplexQuadCase.m have been eliminated.  The file
	MathToCplex.m now has everything (I hope) needed to
	do the verifications of SPIII.

	Revision 1.18  1998/05/25 22:47:39  hales
	Major revision of MathToCplex.m.   There was one small error
	in a constant in the previous version (1.486615 vs. 1.48665),
	but all the other revisions were made to simplify the code.
	It is much much cleaner now.  I checked a cplex output file
	to verify that all the inequalities are correct.  It looks good.

	Revision 1.17  1998/05/03 23:33:19  hales
	another small bug fixed.

	Revision 1.16  1998/04/22 16:58:46  hales
	*** empty log message ***

	Revision 1.15  1998/04/03 02:19:39  hales
	textheaderX added, AnyPosition modified.

	Revision 1.14  1998/02/24 05:43:57  hales
	dihSlack variable added. This relaxes allowing dihedral
	bounds to fail with a penalty.

	Revision 1.13  1998/02/23 03:50:01  hales
	Some stuff for exceptionals that was never seriously used
	was moved to another file MathToPenalty.old.m

	Revision 1.12  1998/02/12 15:19:58  hales
	MakeExec rewritten as Exec.  Doesn't affect any inequalities.

	Revision 1.11  1998/02/12 03:35:56  hales
	revision number added.

	Revision 1.10  1998/02/12 03:33:15  hales
	This version has a conservative but presumably correct
	list of inequalities for EXCEPTIONALS.  I have cut back enough
	that things should never have to be rerun with weaker inequalities.

	Revision 1.9  1998/02/06 23:01:14  hales
	several obsolete functions have been deleted.

	Revision 1.8  1997/12/03 00:33:36  hales
	This is the version that was used in the treatment of 18 cases
	that remained for Part III.

	Revision 1.7  1997/10/24 01:10:13  hales
	carriage returns added because cplex complained of long lines

	Revision 1.6  1997/10/23 12:51:53  hales
	"drastic corrections" factors changed to a good number
	n x 0.01561 x pt

	Revision 1.5  1997/10/23 00:43:25  hales
	Major changes have been made.
	Code for excesses around crowded vertices type (4,1) or (3,2)...
	This is not rigorous! and must be changed because there is a
	factor called the "drastic correction" measuring the effect of
	replacing everything by VorVc, and this factor is not what it needs
	to be.
	Code for the score of an exceptional cluster scored by VorVc has
	been added ("drastic correction" needs adjustment).

	Revision 1.4  1997/09/19 12:48:30  hales
	An error was found Sept 18, 1997 that affects everything with
	and exceptional cluster.  The error was in qedges as it relates to
	the edgeineqIsoPeri string.  It was applying it to exceptionals,
	when it should only apply to quads.  Fixed in this new version.

	Revision 1.3  1997/09/16 03:01:49  hales
	There was a bad constant in quad42.CC that has now been fixed
	(A counterexample was found in part3openSqc)
	The constant 0.7602884 was changed to 0.7624

	Revision 1.2  1997/09/14 19:00:27  hales
	Height cascade material was deleted.
	Two constants in group3 qrtet inequalities were changed,
		because counterexamples were found 3.15 and 3.16.
	group10 was renamed group5 for compatibility with the paper SPIII.

*)
rM:= << MathToCplex.m;
ErrorLog={};
Error[x_]:= AppendTo[ErrorLog,{datestring,"\n","config=",S[GBLconfig],"\n",
	"Length[ConfigurationList]=",S[Length[ConfigurationList]],"\n",x}//StringJoin];
	
(* AnyPosition[{a,b,c},b] -> 2 *)
AnyPosition[list_,x_]:= If[Count[list,x]==0,0,
	Min[Map[First,Position[list,x]]]];
PSelect[x_,f_]:= Select[Range[Length[x]],f[x[[#]]]&];
Lg[x_]:= Length[x];
Inside[x_,list_]:= (Count[list,x]>0);
Subset[sub_,C_]:= (Length[Complement[sub,C]]==0);

(* MoveFirst[{3,4,5},5] -> {5,3,4} *)
MoveFirst[list_,i_]:= RotateLeft[list,AnyPosition[list,i]-1];
Nonempty[x__]:=(Length[Intersection[x]]>0);
Deselect[x_,f_]:= Select[x,!f[#]&];
Flatten1[x_]:= Flatten[x,1];

(* Period[{a,b,c},7] -> a *)
Period[x_,i_]:= x[[ 1 + Mod[i-1,Lg[x]] ]];
IMod[i_,p_]:= 1 + Mod[i-1,p];
AugMod[i_,p_]:= IMod[i+1,p];
CR[i_,mod_]:= If[0==Mod[i,mod],"\n",""];

(* SR[33.44] -> +33.44 *)
SR[x_]:= If[N[x]>=0," +"<>S[x]," " <> S[x]];
S[x_]:= ToString[CForm[x]];


(* StdRegions[finalX[[107,3]] -> {{1, 3, 2}, {2, 3, 8}, ..., {1, 5, 4, 3}}; *)

StdRegions[arrang_]:= Module[{OneStdRegion,OneStdRegionAt},
 
	OneStdRegionAt[i_,j_,a_]:= Module[{b,s,u,t},
				b = Map[#[[2]]&,a];
				t = {i,u=b[[i,j]]};
				While[Last[t]!=i,s = b[[Last[t]]];
					q = AnyPosition[s,t[[-2]] ];
					AppendTo[t,Period[s,q-1] ];
					];
				t = Drop[t,-1];
				MoveFirst[t,Min[t]]
				];
 
	OneStdRegion[i_, a_] := Module[{j},
	  Table[OneStdRegionAt[i, j, a], {j, 1, Length[a[[i,2]]]}]];
 
	
	Array[OneStdRegion[#,arrang]&,Lg[arrang]]//Flatten1//Union
	];

prune[list_]:= Module[{r=list},
		If[Length[list[[1]]]>2 && list[[1,3]]==0,
		r[[1]]=list[[1,{1,2}]] ]; r];



(* Initialize[107] -> { GBLconfig -> 107, 
	GBLregions -> {{1, 3, 2}, {2, 3, 8},... },
*)

Initialize[config_]:= Module[{list,arrang},
		list = ConfigurationList[[config,3]];
		arrang = Map[prune,list ];  
		GBLconfig = config;
		GBLregions  = StdRegions[arrang];
		];

VertexType[i_]:= Module[{r},
	r=Map[Length,Select[GBLregions,Inside[i,#]&]];
	{Count[r,3],Count[r,4],Length[Select[r,#>4&]]}
	];

(* datestring -> "7/28/1997" *)
datestring:=
	Module[{d=Date[]},S[d[[2]]]<>"/"<>S[d[[3]]]<>"/"<>S[d[[1]]]<>
			" "<>S[d[[4]]]<>":"<>S[d[[5]]]<>":"<>S[d[[6]]] ];

(* sigsum[2] -> "score :  + sigma{1} + sigma{2}" *)
sigsum[n_]:= Module[{s="score : "},
	Array[(s = s<> " + "<>StringVar["sig",#]<> 
			If[0==Mod[#,5],"\n",""])&,n];
	s];

textheader:= {
	"\\ LP Format CPLEX file generated by Mathematica on ",
	datestring,
	"\n\n\n",
	"\\ Problem: Bound the score of configuration number ", S[GBLconfig],
	"/",S[Length[ConfigurationList]],
	"\n\\   arising in Part III/IV/V of the Kepler Conjecture.",
	"\n\n\\ $Revision: 1.26 $",
	"\n \n \n",
	(* if dihSlack!=0, get relaxation allowing dihedral bounds to break*)
	"MAXIMIZE \nX - dihSlack\n",
	"\n\nST \n",
	sigsum[Length[GBLregions]]," -sigsum= 0\n",
	"sigsumX: X-sigsum=0\n"
	}//StringJoin;

exec4429:= {
	"change delete constraints sigsumX\n",
	"add\n",
	"sig4429: sigsum>0.4429\n",
	"end\n"
	}//StringJoin;

(* twopi ->
	 \ The dihedral angles around each vertex sum to 2pi : 
         dihsum1 : - 2 pi  + dih{1}1 + dih{19}1 + dih{20}1 = 0
         dihsum2 : - 2 pi  + dih{1}2 + dih{2}2 + dih{3}2 + dih{4}2 + dih{19}2\
         >   = 0
		....		
		dihsum13 : - 2 pi  + dih{14}13 + dih{15}13 + dih{16}13 + dih{17}13 +\
          
         >   dih{18}13 = 0
		*)

twopi:= Module[{i,j,dihList,t},
	{"\n\n",
	"\\ The dihedral angles around each vertex sum to 2pi : ",
	Table[
	t= Position[GBLregions,i]; 
	dihList = Map[" + "<>StringVar["dih",#[[1]],#[[2]]]&,t];
	{
	"\n",
	"dihsum" <> S[i] <> " : -2 pi ",
	dihList,
	" =0\n"
	}
	,{i,1,Max[GBLregions]}]
	}//StringJoin
	];

(* tauText:
	- tau{1} +0.100444571427056 solid{1} - sigma{1}=0
	- tau{2} +0.100444571427056 solid{2} - sigma{2}=0
	- tau{3} +0.100444571427056 solid{3} - sigma{3}=0
	etc.
*)

tauText:= Array[
	{" - ",StringVar["tau",#], SR[0.100444571427056] , " ", 
	StringVar["sol",#]," - ", StringVar["sig",#], "=0\n"}&,
	 Length[GBLregions]]//StringJoin;

(* StringVar .... Generate Variable names *)

StringVar["y",i_]:= "y("<>S[i]<>")";
StringVar["y",i_,j_]:= "y("<>S[Min[i,j]]<>","<>S[Max[i,j]]<>")";
StringVar["yG",f_,i_]:= StringVar["y",Period[GBLregions[[f]],i]];
StringVar["yG",f_,i_,j_]:= 
	StringVar["y",Period[GBLregions[[f]],i],Period[GBLregions[[f]],j]];
StringVar["dih",f_,pos_]:= "dih{"<>S[f]<>"}"<>S[Period[GBLregions[[f]],pos]];
StringVar["Adih",f_,pos_]:= "A"<>StringVar["dih",f,pos];
StringVar["sol",f_Integer]:= "solid{"<>S[f]<>"}";
StringVar["sig",f_Integer]:= "sigma{"<>S[f]<>"}";
StringVar["sol",f_,pos_Integer]:= StringVar["sol",f]<>S[GBLregions[[f,pos]]];
StringVar["sig",f_,pos_Integer]:= StringVar["sig",f]<>S[GBLregions[[f,pos]]];
StringVar["dih",f_,pos_,pos2_]:= StringVar["dih",f,pos]<>".s("<>
			S[GBLregions[[f,pos2]]]<>")";
StringVar["tau",f_]:= "tau{"<>S[f]<>"}";
StringVar["ht",f_]:= "ht{"<>S[f]<>"}213";



(********************************  BOUNDS ********************************)


(* bounds ->
		 BOUNDS
 
 
          pi = 3.141592653589793
          pt = 0.05537364566846414
 
         \Give bounds on dihedral angles
         0.8538 < dih{1}1 < 1.874445
         0.8538 < dih{1}3 < 1.874
 
         \ Give bounds on edge lengths   *)

bounds:= Module[{r,i,j,pt,str,dihmin,dihmax,dihQmin,dihQmax,d1,d2,
	edgeNames},
	dihmin = S[0.8538];
	dihmax = S[1.874445];
	dihQmin= S[1.153];
	dihQmax= S[6.29];
	pt = S[0.0553736456684637];
	str="\n\nBOUNDS\n\n" <>
	"\n pi = 3.141592653589793" <>
	"\n pt = 0.05537364566846414" 
	<> "\n\n\\Give bounds on dihedral angles";
	Do[  If[Length[GBLregions[[i]]]<4,
			d1=dihmin; d2=dihmax, (* else *)
			d1=dihQmin; d2=dihQmax
		];
	   str = str <> "\n" <> d1 <> " < " <> StringVar["dih",i,j] <> " < " <> d2,
		{i,1,Lg[GBLregions]},{j,1,Lg[GBLregions[[i]]]}];
	str = str <> "\n\n\\ Give bounds on scores";
	Do[ If[Length[GBLregions[[i]]]<4,d1=pt,d1= "0"];
	  str = str<>"\n -infinity < "<>StringVar["sig",i] <> " < " <> d1,
		{i,1,Lg[GBLregions]}
		];
	str = str <> "\n\n\\ Give bounds on edge lengths";
	edgeNames = MakeEdgeNames;
	Do[str = str<> "\n 2 < "<> edgeNames[[i]]<> " < 2.51",
		{i,1,Lg[edgeNames]}
		];
	Do[str = str<> "\n "<>freeVar[[i]]<>" free",{i,1,Length[freeVar]}];
	str
	];

MakeEdgeNames:= Module[{i,f,t1,t2},
	t1 = Array[StringVar["y",#]&,Max[GBLregions]];
	t2 = Table[StringVar["yG",f,i,i+1],
			{f,1,Length[GBLregions]},
			{i,1,Length[GBLregions[[f]] ]}];
	Union[t1,t2]//Flatten
	];

(**************************************  PART III QRTET *******************)

(* textQR[{"label",4,5,6},{f,v}] -> 
		labelf.v  :  sigma{f} -4 solid{f} -5 dih{f}v < 6 *)

textQR[{eqnid_,solid_,dih_,const_},{f_,v_}]:=   
	(* sigma < {solid,const,dih} *) 
	If[eqnid=="","",eqnid <> S[f]<>"."<>S[v] <> " :  "] <> 
	StringVar["sig",f]<>
	If[solid!=0,SR[-solid]<>" "<> StringVar["sol",f],""]<>
	If[dih!=0,SR[-dih]<>" " <> StringVar["dih",f,v],""]<>
	" < " <> S[const];

(* trieqn["title","label",{3,4,5}] ->
		 \ title
         label{1}1 :  sigma{1} -3 solid{1} -4 dih{1}1 < 5
         label{1}3 :  sigma{1} -3 solid{1} -4 dih{1}3 < 5
		 ....
		 label{18}13 :  sigma{18} -3 solid{18} -4 dih{18}13 < 5  *)

trieqn[tag_,lab_,{sol_,dih_,con_}]:= Module[{r,f,j,str},
	str = "\n\n\\ "<> tag;
	Do[r = Lg[GBLregions[[f]]];
		If[r<4,
	      str=str<>"\n"<>
		textQR[{lab,sol,dih,con},{f,j}]],
		{f,1,Lg[GBLregions]},{j,1,Lg[GBLregions[[f]]]}
		];
	str
	];

group1:= (
	trieqn["Group 1. Equation 4.","Grp1E4",{-0.37642101,0,0.287389}] <>
	trieqn["Group 1. Equation 5.","Grp1E5",{ 0.446634,0,-0.190249}] <>
	trieqn["Group 1. Equation 6.","Grp1E6",{-0.419351,0,0.2856354+0.001}]
	  );

group3:= Module[{zp=0.1004445714270568,zp32=0.321422628566582,
		ax=-0.419351},
	trieqn["Group 3. Equation 1.","Grp3E1.",{0, 0.37898, -0.4111}]<>
	trieqn["Group 3. Equation 2.","Grp3E2.",{0, -0.142, 0.23021}]<>
	trieqn["Group 3. Equation 3.","Grp3E3.",{0, -0.3302, 0.5353}]<>
	trieqn["Group 3. Equation 4.","Grp3E4.",{zp, 0.3897, -0.4666}]<>
	trieqn["Group 3. Equation 5.","Grp3E5.",{zp, 0.2993, -0.3683}]<>
	"\n\n\\ Group 3. Equation 6 appears in the bounds section" <>
	trieqn["Group 3. Equation 7.","Grp3E7.",{zp, -0.1689, 0.208}]<>
	trieqn["Group 3. Equation 8.","Grp3E8.",{zp, -0.2529, 0.3442}]<>
	trieqn["Group 3. Equation 9.","Grp3E9.",{zp32, 0.4233, -0.5974}]<>
	trieqn["Group 3. Equation 10.","Grp3E10.",{zp32, 0.1083, -0.255}]<>
	trieqn["Group 3. Equation 11.","Grp3E11.",{zp32, -0.0953, -0.0045}]<>
	trieqn["Group 3. Equation 12.","Grp3E12.",{zp32, -0.1966, 0.1369}]<>
	trieqn["Group 3. Equation 13.","Grp3E13.",{ax, 0.796456, -0.5786316}]<>
	trieqn["Group 3. Equation 14.","Grp3E14.",{ax, 0.0610397, 0.211419}]<>
	trieqn["Group 3. Equation 15.","Grp3E15.",{ax, -0.0162028, 0.308526}]<>
	trieqn["Group 3. Equation 16.","Grp3E16.",{ax, -0.0499559, 0.35641}]<>
	trieqn["Group 3. Equation 17.","Grp3E17.",{ax, -0.64713719, 1.3225}]
	  ];


(******************************** PART III QUAD **************************)

(* text42[{"label",4,5},{"X","Y","Z"}] ->
			 labelXY :  sigmaX -4 dihXY -4 dihXZ < 5 *)

text42[{eqnid_,dih_,const_},{f_,v1_,v2_}]:=   
	(* sigma < {const,dih} *) 
	eqnid <> S[f] <> "." <> S[v1] <> " :  " <> 
	StringVar["sig",f]<>
	If[dih!=0,SR[-dih] <> StringVar["dih",f,v1],""]<>
	If[dih!=0,SR[-dih] <> StringVar["dih",f,v2],""]<>
	" < " <> S[const];

(* quad["title","label",1] ->
 
         \ title
         label{19}1 :  sigma{19} -4.56766 dih{19}1 < -5.7906
         label{19}2 :  sigma{19} -4.56766 dih{19}2 < -5.7906
         label{19}6 :  sigma{19} -4.56766 dih{19}6 < -5.7906
         label{19}5 :  sigma{19} -4.56766 dih{19}5 < -5.7906
         label{20}1 :  sigma{20} -4.56766 dih{20}1 < -5.7906
         label{20}5 :  sigma{20} -4.56766 dih{20}5 < -5.7906
         label{20}4 :  sigma{20} -4.56766 dih{20}4 < -5.7906
         label{20}3 :  sigma{20} -4.56766 dih{20}3 < -5.7906   *)
			

CCHASH = 2032330977; 
quad[tag_,lab_,q_]:= Module[{r,f,j,str,sol,dih,con,CC,zp,zp32},
	zp=0.1004445714270568;
	zp32=0.321422628566582;
	CC={ 
	{0,-5.7906,4.56766}, 
        {0,-2.0749,1.5094},
        {0,-0.8341,0.5301},
        {0,-0.6284,0.3878},  
        {0,0.4124,-0.1897},
        {0,1.5707,-0.5905},
        {-0.3,0.41717,0}, 
        {zp,-5.81446,4.49461},  
        {zp,-2.955,2.1406}, 
        {zp,-0.6438,0.316},
        {zp,-0.1317,0.0},
        {zp,0.3825,-0.2365},
        {zp,1.071,-0.4747},
        {zp32,-5.77942,4.25863}, 
        {zp32,-4.893,3.5294},
        {zp32,-0.4126,0.0},
        {zp32,0.33,-0.316}, 
        {-0.419351,-5.350181,4.611391},
        {-0.419351,-1.66174,1.582508}, 
        {-0.419351,0.0895,0.342747}, 
        {-0.419351,3.36909,-0.974137}  
        };

	{sol,con,dih}=CC[[q]];
	str = "\n\n\\ "<> tag;
	Do[r = Lg[GBLregions[[f]]];
		If[r==4,
	      str=str<>"\n"<>
		textQR[{lab,sol,dih,con}, {f,j} ]],
		{f,1,Lg[GBLregions]},{j,1,Lg[GBLregions[[f]]]}
		];
	str
	];



quad42[tag_,lab_,q_]:= Module[{r,i,j,str,dih,con,CC},
	CC={ (* 5/24/97 *)
		{-9.494,3.0508},
		{-1.0472,0.27605},
(*		{0.7624,-0.198867}, *)  (* removed 5/25/98 *)
		{3.5926,-0.844}
        };
	If[q>Length[CC],Return[""]];
	{con,dih}=CC[[q]];
	str = "\n\n\\ "<> tag;
	Do[r = Lg[GBLregions[[i]]];
		If[r==4,
	      str=str<>"\n"<>
		text42[{lab,dih,con}, {i,j,j+1} ]],
		{i,1,Lg[GBLregions]},{j,1,Lg[GBLregions[[i]]]}
		];
	str
	];

(* quadeqn42[1]
 
         \ Quad 4.2.1
         Qu421{19}1 :  sigma{19} -3.0508 dih{19}1 -3.0508 dih{19}2 < -9.494
         Qu421{19}2 :  sigma{19} -3.0508 dih{19}2 -3.0508 dih{19}6 < -9.494
         Qu421{19}6 :  sigma{19} -3.0508 dih{19}6 -3.0508 dih{19}5 < -9.494
         Qu421{19}5 :  sigma{19} -3.0508 dih{19}5 -3.0508 dih{19}1 < -9.494
         Qu421{20}1 :  sigma{20} -3.0508 dih{20}1 -3.0508 dih{20}5 < -9.494
         Qu421{20}5 :  sigma{20} -3.0508 dih{20}5 -3.0508 dih{20}4 < -9.494
         Qu421{20}4 :  sigma{20} -3.0508 dih{20}4 -3.0508 dih{20}3 < -9.494
         Qu421{20}3 :  sigma{20} -3.0508 dih{20}3 -3.0508 dih{20}1 < -9.494 *)



hyp41:= Module[{i,quadeqn,str=""},
	quadeqn[i_]:= quad["Quad 4.1."<>S[i],"Qu41"<>S[i],i];
	Do[str = str<> quadeqn[i],{i,1,21}];
	str
	];

hyp42:= Array[ quad42["Quad 4.2."<>S[#],"Qu42"<>S[#],#]&, 4]//StringJoin

solid:= Module[{r,i},
	{
	"\n\n\\ Solid angles in terms of dihedral angles\n",
	Table[r = Lg[GBLregions[[i]]];
		{
		"DS",S[i]," :  - ",StringVar["sol",i],
		Array[" + "<>StringVar["dih",i,#]&,r],
		" - "<>S[r-2]<>" pi = 0\n"
		},
	{i,1,Lg[GBLregions]}]
	}//StringJoin
	];


group55:= Module[{v,q,str,i,j,s,t},
	v = Select[Range[1,Max[GBLregions]],VertexType[#]=={5,0,0}&];
	If[Length[v]<1,Return];
	{"\n\n\\ The group 4.5 (=I.5.1.1) inequalities for type (5,0) vertices \n",
	Table[ i = v[[q]]; t= Map[First,Position[GBLregions,i]]; 
	{"\nGrp105Vertex" <> S[i] <> " : ",
	Map[{" + ",StringVar["sig",#],"+0.419351",StringVar["sol",#]}&,t],
	 " <  1.428177 \\ = 5(0.2856354)\n"
	}
	,{q,1,Length[v]}]
	}//StringJoin
	];

(* group54 -> \ The group 5 inequalities for type (4,0) vertices  *)

group54:= Module[{v,q,str,i,s,t},
	v = Select[Range[1,Max[GBLregions]],VertexType[#]=={4,0,0}&];
	If[Length[v]<1,Return[""]];
	{"\n\n\\ The group 5 inequalities for type (4,0) vertices \n",
	Table[ i = v[[q]]; t= Map[First,Position[GBLregions,i]]; 
		{ "\nGrp101Vertex" <> S[i] <> " :  -0.33 pt ",
		Map[{" + ",StringVar["sig",#]}&,t],
		" < 0\n"
		}
	,{q,1,Length[v]}]
	}//StringJoin
	];

group5:= (group55<>group54);

WRITEOUTstd[stream_]:=WRITEOUTstd[stream,""];
WRITEOUTstd[stream_,basefile_]:=
	(WriteString[stream,textheader];
        WriteString[stream,twopi];
        WriteString[stream,tauText];
        WriteString[stream,group1];
        WriteString[stream,group3];
        WriteString[stream,group5];
        WriteString[stream,hyp41];
        WriteString[stream,hyp42];
        WriteString[stream,solid];
        WriteString[stream,edges];
		WriteString[stream,QRText]; 
        WriteString[stream,qedges];
		WriteString[stream,QuadVarRelations];
		WriteString[stream,quadInstall[basefile]];
        WriteString[stream,ValenceFiveText];
        WriteString[stream,ExceptFaceText];
        WriteString[stream,VertexExcessText]; 
	);

WRITEOUT[i_]:= Module[{stream},
		Initialize[i];
		freeVar={"X","sigsum"};
		stream=OpenWrite["/tmp/cplex8.lp"<>S[i]];
		(* on first pass use WRITEOUTstd[stream]; *)
		WRITEOUTstd[stream,"TEMP/cplex8.lp"<>S[i]];
		WriteString[stream,bounds];
		WriteString[stream,"\n\nEND\n\n"];
		WriteString[stream,faceCode];
		Close[stream];
		];

faceCode:= Module[{str},
		str = Array[{Length[GBLregions[[#]]],GBLregions[[#]]}&,
				Length[GBLregions]];
		str = Flatten[PrependTo[str,Length[GBLregions]]]//ToString;
		str = "# "<>str;
		str = StringReplace[str,{","->"","{"->"","}"->""}];
		str
		];


(*********************** Valence Five STUFF ************************)

(* ValenceFiveVertices ->
   {{7, 11}, {7, 12}, {8, 9}, {8, 10}, {8, 13}, {9, 12}, {9, 13}, 
    {10, 11}, {10, 13}, {11, 12}, {11, 13}, {12, 13}, {7, 9, 12}, ...., 
    {10, 11, 12, 13}}  *)

(* Valence 5 stuff *)
(* This gives the various collections of at most four vertices that are all
	of type (5,0) and are all connected together *)
ValenceFiveVertices:=
	Module[{v5,rawlist,v}, 
	v5 = Select[Range[1,Max[GBLregions]],VertexType[#]=={5,0,0}&];
	rawlist= Module[{x,i,j},
		Table[j=v5[[i]]; 
		x=Intersection[Flatten[Select[GBLregions,Count[#,j]>0&]],v5]~
			Complement~{j};
		Map[Join[{j},#]&,subsets[Length[x]]/.Array[#->x[[#]]&,Length[x]]],
		{i,1,Length[v5]}] ~Flatten~ 1
		];
	v = Select[Union[rawlist,Map[{#}&,v5]],Length[#]<5&];
	v = Map[Sort,Select[v,Subset[#,v5]&]]//Union;
	If[Union[Flatten[v]]!=v5,Error["ValenceFiveVertices"]];
	v
	];

(* Valence5txt[{7,11}]
         +sigma{7}+sigma{8}+sigma{9}+sigma{10}+sigma{12}
         +sigma{13}+sigma{17}+sigma{18} -7.04 pt <0
         +sigma{7}- 0.1004445714270561 solid{7}+sigma{8}- 0.1004445714270561\
         >   solid{8}+sigma{9}- 0.1004445714270561 solid{9}+sigma{10}-\
         >   0.1004445714270561 solid{10}+sigma{12}- 0.1004445714270561\
          
         >   solid{12} +sigma{13}- 0.1004445714270561 solid{13}+sigma{17}-\
         >   0.1004445714270561 solid{17}+sigma{18}- 0.1004445714270561\
         >   solid{18} +1.1 pt < 0  *)

Valence5txt[v_]:= Module[{fac,s,i},
	fac=Select[Range[1,Lg[GBLregions]],Nonempty[GBLregions[[#]],v]&];
	s ="";
	Do[s = s<> "+"<>StringVar["sig",fac[[i]]]<>CR[i,5],{i,1,Length[fac]}];
	s = s<> SR[Length[v] 0.48- Length[fac]]<>" pt < 0\n";
	Do[s = s<> "+"<>StringVar["tau",fac[[i]]]<>CR[i,5],{i,1,Length[fac]}];
	s = s<> SR[-Length[v] 0.55]<>" pt > 0\n\n";
	s
	];

ValenceFiveText:= {
	"\n\n\\ Valence 5 inequalities: \n\n",
	Map[Valence5txt,ValenceFiveVertices]
	}//StringJoin

(********************* EDGE LENGTHS ***************************)

(* edge length stuff *)
edgeineq = 
	{"-solid + 0.199235 y4 + 0.199235 y5 + 0.199235 y6 -\
		0.377076 y1 - 0.377076 y2 - 0.377076 y3 < -1.618331\n",
	"solid - 0.320937 y4 - 0.320937 y5 - 0.320937 y6 +\
		0.152679 y1 + 0.152679 y2 + 0.152679 y3 < -0.458262\n",
	"sigma + 0.10857 y1 + 0.10857 y2 + 0.10857 y3 + \
		0.10857 y4 + 0.10857 y5 + 0.10857 y6 < 1.3582137\n",
	"sigma + 0.419351 solid + 0.2 y1 + 0.2 y2 + 0.2 y3 < 1.48665\n",
	"sigma - 0.1004445714270561 solid + 0.129119 y4 + 0.129119 y5 +\
		0.129119 y6 + 0.0845696 y1 + 0.0845696 y2 + \
		0.0845696 y3 < 1.2821326\n",
	"dih1 + 0.153598 y2 + 0.153598 y3 + 0.153598 y5 + 0.153598 y6 -\
		0.498 y1 - 0.76446 y4 < -0.065176\n",
	"dih2 + 0.153598 y1 + 0.153598 y3 + 0.153598 y4 + 0.153598 y6 -\
		0.498 y2 - 0.76446 y5 < -0.065176\n",
	"dih3 + 0.153598 y1 + 0.153598 y2 + 0.153598 y4 + 0.153598 y5 -\
		0.498 y3 - 0.76446 y6 < -0.065176\n",
	"-dih1 - 0.359894 y2 - 0.359894 y3 - 0.359894 y5 - 0.359894 y6 +\
		0.003 y1 + 0.685 y4 < -2.734102\n",
	"-dih2 - 0.359894 y1 - 0.359894 y3 - 0.359894 y4 - 0.359894 y6 +\
		0.003 y2 + 0.685 y5 < -2.734102\n",
	"-dih3 - 0.359894 y1 - 0.359894 y2 - 0.359894 y4 - 0.359894 y5 +\
		0.003 y3 + 0.685 y6 < -2.734102\n\n\n"}//StringJoin


EdgeInequality[f_]:= "\\ Group 2, face "<>S[f]<>"\n"<>
	StringReplace[edgeineq,
	{"solid"->StringVar["sol",f],
	 "dih1"->StringVar["dih",f,1],
	 "dih2"->StringVar["dih",f,2],
	 "dih3"->StringVar["dih",f,3],
	 "sigma"->StringVar["sig",f],
	 "y1"->StringVar["yG",f,1],
	 "y2"->StringVar["yG",f,2],
	 "y3"->StringVar["yG",f,3],
	 "y4"->StringVar["yG",f,2,3],
	 "y5"->StringVar["yG",f,1,3],
	 "y6"->StringVar["yG",f,1,2]
	}];


edges:= Module[{li},
	(* edgeNames={}; *)
	li = Select[Range[1,Length[GBLregions]],Length[GBLregions[[#]]]==3&];
	{
	"\\ Group 2 \n",
	Array[EdgeInequality[li[[#]]]&,Length[li]]
	}//StringJoin
	];

qedges:= Module[{li,St,qedge},
	qedge[i_]:= Array[qedge[i,#]&,Length[GBLregions[[i]]]  ];
	qedge[f_,r_]:= StringReplace[
 "- dih + 0.3257 y1 - 0.398 y2 -0.398 y3 -0.398 y5 -0.398 y6 < -4.14938\n" ,
	  {"dih"->StringVar["dih",f,r],
	 "y1"->StringVar["yG",f,r],
	 "y2"->StringVar["yG",f,r+1],
	 "y3"->StringVar["yG",f,r-1],
	 "y5"->StringVar["yG",f,r,r-1],
	 "y6"->StringVar["yG",f,r,r+1]
		}];
	li = Select[Range[1,Length[GBLregions]],Length[GBLregions[[#]]]>3&];
	{
	"\\ III.group 4.#6\n",
	Array[qedge[li[[#]] ]&,Length[li] ]
	}//StringJoin
	];


(******************** EXCEPTIONAL STUFF  *********************)

(* This is really rough stuff for the exceptions, used to
	obtain the SHORT/shortlist.m of 180 configurations. 
	Better inequalities are developed in MathToCplexExcept.m *)

(* ExceptionalStuff *)


SquanderFace= {0,0,0,0.1317,0.27113,0.41056,0.54999,0.6045};
ScoreFace  = {0,0,0,0,-0.05704,-0.11408,-0.1677,-0.1677};

fSquander[flist_List]:= 
	(Plus @@ Map[SquanderFace[[Length[GBLregions[[#]] ] ]]&,flist]) 

ExceptFaceText:= Module[{exc,g},
		exc = Select[Range[GBLregions//Lg],Length[GBLregions[[#]]]>4&];
		{"\n\n\\ Exceptional Face Text: \n\n",
		Map[
			 (g = GBLregions[[ # ]]//Length;
			 {StringVar["sig",#]," < ",S[ScoreFace[[g]] ],"\n",
			 StringVar["tau",#]," > ",S[SquanderFace[[g]] ],"\n"})&
		  ,exc]
		}//StringJoin
		];











(* list all subsets of {1,...,n} *)
subsets[0]:= {};
subsets[n_]:= subsets[n] =
	Union[subsets[n-1],Map[Join[{n},#]&,subsets[n-1]]]; /; n>1;
subsets[1]:= {{},{1}};

Edges[p_]:= Module[{EdgeOfCycle,t},
    EdgeOfCycle[t_]:= Array[Sort[{Period[t,#],Period[t,#+1]}]&,Length[t]];
    Flatten[Map[EdgeOfCycle,p],1]//Sort
    ];


	(* contributions of 1.4 or 1.5 at each vertex of types 
		(4,0,1),(3,1,1),(3,0,2)*)
VertexExcessText:= Module[{subn,edges,i,vertices,faces},
	(*build nonadjacVertices*) 
	subn=Module[{verticesOnF ,exVertices}, 
		verticesOnF = Select[GBLregions,Length[#]>4&]//Flatten//Union;
		exVertices=Select[verticesOnF,
			Inside[VertexType[#],{{4,0,1},{3,1,1},{3,0,2}}]&];
		Map[exVertices[[#]]&,
			Complement[subsets[Length[exVertices]],{{}}]]
		];
	(*build string*) 
	{
	 "\n\n\\ Vertex Excess Text\n\n", 
		edges = Edges[GBLregions]//Union;
		Table[
		vertices = subn[[i]];
		If[HasAdjacentPair[vertices,edges],"",
			faces = PSelect[GBLregions,Nonempty[#,vertices]&];
			{
			"vet",Map[("."<>S[#])&,vertices],": ",
			Array[{"+",StringVar["tau",faces[[#]]],CR[#,10]}&,Length[faces]],
			" > ",S[vSquander[vertices]+fSquander[faces]],"\n"
			}
		  ]
		,{i,1,Length[subn]}
			]
	}//StringJoin
	];

HasAdjacentPair[v_,edges_]:= Module[{e2,i,j},
	If[Length[v]<2,Return[False]];
	e2= 
	Flatten[Table[Sort[{v[[i]],v[[j]]}],
			{i,1,Length[v]},{j,i+1,Length[v]}],1];
	Nonempty[edges,e2]
	];

vSquander1[vNumber_Integer]:= Module[{pt},
	pt = 0.05537364566846414;
	Switch[VertexType[vNumber],
	{4,0,1},1.5 pt,
	{3,1,1},1.4 pt,
	{3,0,2},1.4 pt,
	_,0]];

vSquander[vlist_List]:= Plus @@ Map[vSquander1,vlist];


(*********************** QUAD CLUSTERS MATERIAL ********************)


(* division of quads into cases :
	  two flat quarters, split between vertices 1 & 3
	  two flat quarters, split between vertices 2 & 4,
	  octahedra, scored gamma/octavor on each upright quarter
	  no quarters, scored vorVc 
*)

quadInstall[basefile_]:= Module[{flist},
    flist=Select[Range[Length[GBLregions]],Length[GBLregions[[#]]]==4&];
	{"\n\n\\ quadInstall:\n",
    Array[quadInstallOne[basefile,flist[[#]]]&,Length[flist]]
		} //StringJoin
    ];

(* StringVarQuad *)
StringVar["yoct",f_]:= "y(F"<>S[f]<>")";
StringVar["yoct",f_,i_]:= "y(F"<>S[f]<>","<>S[GBLregions[[f,i]]]<>")";
StringVar["solF",f_,i_,j_]:= Module[{i1,j1},
    {i1,j1}=Sort[{GBLregions[[f,i]],GBLregions[[f,j]]}];
    "sol(F"<>S[f]<>","<>S[i1]<>","<>S[j1]<>")"];
StringVar["sigF",f_,i_,j_]:= Module[{i1,j1},
    {i1,j1}=Sort[{GBLregions[[f,i]],GBLregions[[f,j]]}];
    "sig(F"<>S[f]<>","<>S[i1]<>","<>S[j1]<>")"];
StringVar["dihF",f_,i_,j_]:= Module[{i1,j1},
    {i1,j1}=Sort[{GBLregions[[f,i]],GBLregions[[f,j]]}];
    "dih(F"<>S[f]<>","<>S[i1]<>","<>S[j1]<>")"];
StringVar["dihFx",i_,j_,f_]:="dih("<>S[GBLregions[[f,i]]]<>","<>
        S[GBLregions[[f,j]]]<>",F"<>S[f]<>")";
StringVar["slack",s_,f_]:= "slack"<>s<>S[f];
(* The following also appear in MathToCplexExcept.m *)
StringVar["pen",f_]:= "pen"<>S[f];
StringVar["cquo",f_,p1_,p2_]:= Module[{i,j},
    i=GBLregions[[f,p1]];
    j=GBLregions[[f,p2]];
    StringVar["cquoI",f,i,j]
    ];
StringVar["cquoI",f_,i_,j_]:= "cquo{"<>S[f]<>"}("<>S[Min[i,j]]<>","<>S[Max[i,j]]<>")";




quadInstallOne[basefile_,f_]:= 
    "\n\n\\ fInstall 412\n"<>
        fInstall[f,4,1,2]<>"\n\n\\ 234\n"<>fInstall[f,2,3,4]<>
    "\n\n\\ fInstall 123\n"<>
        fInstall[f,1,2,3]<>"\n\n\\ 341\n"<>fInstall[f,3,4,1]<>
    "\n\n\\ oInstall\n"<>
        oInstall[f,1,2]<>oInstall[f,2,3]<>oInstall[f,3,4]<>oInstall[f,4,1]<>
    "\n\n\\ vInstall\n"<> vInstall[basefile,f];

flatEquations =
    {
	"\\ flat Equations, flat quarters in the Q-system.\n",
	"\\ mu-scoring only.  Not to be used with erased flat quarters.\n",
    "- dih2 + 0.35 y2 - 0.15 y1 - 0.15 y3 + 0.7022 y5 - 0.17 y4 + slack > -0.0123",
    "- dih3 + 0.35 y3 - 0.15 y1 - 0.15 y2 + 0.7022 y6 - 0.17 y4 + slack > -0.0123",
    "  dih2 - 0.13 y2 + 0.631 y1 + "<>
               "0.31 y3 - 0.58 y5 + 0.413 y4 + 0.025 y6 + slack > 2.63363 ",
    "  dih3 - 0.13 y3 + 0.631 y1 + "<>
               "0.31 y2 - 0.58 y6 + 0.413 y4 + 0.025 y5 + slack > 2.63363 ",
    " -dih1 + 0.714 y1 - 0.221 y2 - 0.221 y3 + "<>
               "0.92 y4 - 0.221 y5 - 0.221 y6 + slack > 0.3482",
    "  dih1 - 0.315 y1 + 0.3972 y2 + 0.3972 y3 - "<>
               "0.715 y4 +  0.3972 y5 + 0.3972 y6 + slack > 2.37095",
    "- solid - 0.187 y1 - 0.187 y2 - "<>
               "0.187 y3 + 0.1185 y4 + 0.479 y5 + 0.479 y6 + slack > 0.437235 ",
    "+ solid + 0.488 y1 + 0.488 y2 + "<>
               "0.488 y3 - 0.334 y5 - 0.334 y6 + slack > 2.244 ",
    "- sigma - 0.159 y1 - 0.081 y2 - 0.081 y3 - "<>
               "0.133 y5 - 0.133 y6 + slack > -1.17401",
    "- sigma - 0.419351 solid + 0.0436 y5 + 0.0436 y6 + 0.079431 dih1 "<>
                " + slack > 0.0296 ", 
    " sigma + 0.197 y4 + 0.197 y5 + 0.197 y6 - slack < 1.34521 ", 
	" slack - Cqr1 - Cqr2 = 0",
    " y4 + slack > 2.51",
    " y4 - slack < 2.8284271247462"
    };
fEquations= Apply[StringJoin,Map[(#<>"\n")&,flatEquations]];

fInstall[f_,i_,j_,k_]:= (* quad f, vertices i,j,k in[1,2,3,4], *)
      StringReplace[fEquations,
        {"y2"->StringVar["yG",f,i],
         "y1"->StringVar["yG",f,j],
         "y3"->StringVar["yG",f,k],
         "y5"->StringVar["yG",f,j,k],
         "y4"->StringVar["yG",f,i,k],
         "y6"->StringVar["yG",f,i,j],
         "dih2"->StringVar["dih",f,i,j],
         "dih1"->StringVar["dih",f,j],
         "dih3"->StringVar["dih",f,k,j],
         "solid"->StringVar["sol",f,j],
		 "slack"->StringVar["slack",{"A","B"}[[1+Mod[i+j+k,2]]],f],
		 "Cqr1"->"Cqrl"<>S[f 100],
		 "Cqr2"->"Cqr"<>{"s","l"}[[1+Mod[i+j+k,2]]]<>S[f 100+1],
         "sigma"->StringVar["sig",f,j]}
    ];

OctahedralEquations:= {
	"\n\n\\ OctahedralEquations",
	"\n\\ Gamma/octavor scoring.  ",
	"\\ For use only on upright quarters in quad clusters",
    " y1 + slack > 2.51 ",
    " y1 - slack < 2.8284271247462 ",
    " y2 + slack > 2 ",
    " y3 + slack > 2 ",
    " y4 + slack > 2 ",
    " y5 + slack > 2 ",
    " y6 + slack > 2 ",
    " y2 - slack < 2.51 ",
    " y3 - slack < 2.51 ",
    " y4 - slack < 2.51 ",
    " y5 - slack < 2.51 ",
    " y6 - slack < 2.51 ",
    " dih1 - 0.636 y1 + 0.462 y2 + 0.462 y3 - 0.82 y4 + 0.462 y5 + "<>
        " 0.462 y6 + slack > 1.82419 ", (* case 20 *)
    " - dih1 + 0.55 y1 - 0.214 y2 - 0.214 y3 + 1.24 y4 - 0.214 y5 "<>
        " - 0.214 y6 + slack > 0.75281 ",  (* case21 *)
    " dih2 + 0.4 y1 - 0.15 y2 + 0.09 y3 + 0.631 y4 - 0.57 y5 + 0.23 y6 " <>
    "  + slack > 2.5481", (* case 22 *)
    " - dih2 - 0.454 y1 + 0.34 y2 + 0.154 y3 - 0.346 y4 + " <>
        "0.805 y5 + slack > -0.3429", (* case 23 *)
    " dih3 + 0.4 y1 - 0.15 y3 + 0.09 y2 + 0.631 y4 - 0.57 y6 + 0.23 y5 " <>
    "  + slack > 2.5481", (* case 22 *)
    " - dih3 - 0.454 y1 + 0.34 y3 + 0.154 y2 - 0.346 y4 + " <>
        "0.805 y6 + slack > -0.3429", (* case 23 *)
    " sol + 0.065 y2 + 0.065 y3 + 0.061 y4 - 0.115 y5 - "<>
        "0.115 y6 + slack > 0.2618", (* case 24 *)
    " - sol - 0.293 y1 - 0.03 y2 - 0.03 y3 + 0.12 y4 + " <>
        "0.325 y5 + 0.325 y6 + slack > 0.2514", (* case 25 *)
    " - sig - 0.054 y2 - 0.054 y3 - 0.083 y4 - 0.054 y5 - "<>
        "0.054 y6 + slack > -0.59834", (* case 26,27,28 *)
    " - sig - 0.419351 sol + 0.079431 dih2 - 0.0846 y1 + slack > -0.30592 " ,
                (* case 69,70,71, added 11/24/97 *)
    " - sig - 0.419351 sol + 0.079431 dih3 - 0.0846 y1 + slack > -0.30592 ",
                (* case 69,70,71, added 11/24/97 *)
	" slack - Cqrs - Cqrl =0",
	(* This equation holds if y2,y3 < 2.13 *)
	" - sig + 0.07 y1 - 0.133 y2 - 0.133 y3 - 0.135 y4 - "<>
	   "0.133 y5 - 0.133 y6 + ht213 + slack > -1.1583 "
    };


oEquations= Apply[StringJoin,Map[(#<>"\n")&,OctahedralEquations]];

oInstall[f_,i_,j_]:= (* quad f, vertices i,j,k in[1,2,3,4], *)
      StringReplace[oEquations,
        {"y1"->StringVar["yoct",f],
         "y2"->StringVar["yG",f,i],
         "y3"->StringVar["yG",f,j],
         "y4"->StringVar["yG",f,i,j],
         "y5"->StringVar["yoct",f,j],
         "y6"->StringVar["yoct",f,i],
         "dih1"->StringVar["dihF",f,i,j],
         "dih2"->StringVar["dihFx",i,j,f],
         "dih3"->StringVar["dihFx",j,i,f],
         "sol"->StringVar["solF",f,i,j],
		 "slack"->StringVar["slack","Oct",f],
		 "ht213"->StringVar["ht",f],
		 "Cqrs"->"Cqrs"<>S[f 100],
		 "Cqrl"->"Cqrl"<>S[f 100+1],
         "sig"->StringVar["sigF",f,i,j]}
    ];


(* vorVc stuff *)

QuoinText[f_]:= Module[{face,i,j},
    face = GBLregions[[f]];
	{"\n\n\\ Quoin Text (VI.4.10, and III.Appendix)\n",
    Table[
        j = IMod[i+1,Length[face]];
        {StringVar["cquo",f,i,j],
         SR[0.00758+0.0115]," ",StringVar["y",face[[i]] ],
         SR[0.00758+0.0115]," ",StringVar["y",face[[j]] ],
         SR[0.0115 +0.0115]," ",StringVar["y",face[[i]],face[[j]] ],
         " > ", S[2 0.06333],"\n"},
        {i,1,Length[face]}]}//StringJoin
    ];

QuadVcDefText[f_]:= StringReplace[
	StringJoin[{
	"\n\n\\ QuadVcDef The definition of truncated Voronoi\n",
	"vorVc -slack - Adih1 - Adih2 - Adih3 - Adih4 + negphi0 sol \n "<>
		" + 4doct cquo1 + 4doct cquo2 + 4doct cquo3 + 4doct cquo4 < 0\n",
	"Adih1 + slack>0\n",
	"Adih2 + slack>0\n",
	"Adih3 + slack>0\n",
	"Adih4 + slack>0\n"
	}],
		{"vorVc"->StringVar["sig",f],
		 "slack"->StringVar["slack","Vc",f],
		 "Adih1"->StringVar["Adih",f,1],
		 "Adih2"->StringVar["Adih",f,2],
		 "Adih3"->StringVar["Adih",f,3],
		 "Adih4"->StringVar["Adih",f,4],
		 "cquo1"->StringVar["cquo",f,1,2],
		 "cquo2"->StringVar["cquo",f,2,3],
		 "cquo3"->StringVar["cquo",f,3,4],
		 "cquo4"->StringVar["cquo",f,4,1],
		 "negphi0"-> "0.5666365478933329",
		 "4doct"-> "2.88361179806985"
		}
		];

vInstall[basefile_,f_]:= 
		{Array[VcSupplement[f,#]&,4],
		 QuadVcDefText[f], 
		 QuoinText[f],
		 If[basefile=="","",UpdateFace[basefile,f]]
		}//StringJoin


(* vorVc stuff *)
UpdateFace[basefile_,f_]:= Module[{face},
		face=GBLregions[[f]];
		 Array[LPmUpdateAdih[basefile,StringVar["y",face[[#]]],
			StringVar["dih",f,#],StringVar["Adih",f,#]]&,Length[face]]
		//StringJoin
        ];

(* first three arguments are text, last three the numerical bounds *)
AdihText[yvar_,dihvar_,Adihvar_,hmax_,dihmin_,dihmax_]:=
	Module[{A1,slope},
	A1 = (* A[1]= *) 0.109691511444153; 
	slope = 0.4806019663449766*(-1.493759012518814 + hmax)*
				(2.493759012518814 + hmax); (* (A[hmax]-A[1])/(hmax-1)*)
	{Adihvar,SR[-A1]," ",dihvar,SR[-slope dihmin/2]," ",yvar," < ",
		S[-slope dihmin],"\n",
	 Adihvar,SR[-slope dihmax/2]," ",yvar," < ",
		S[(A1-slope) dihmax],"\n"}//StringJoin
	];


VcSupplement[f_,i_]:= Module[{r,eqnlist},
	(* f is a face Number, 1..GBLregions//Length *)
	(* i = 1,2,3,4 *)
	(* unorthodox labels, dih4,quo42,etc. 4 refers to corner opposite 1,
		this allows standard labels 1,2,3 on the front simplex...,
		 *)
 
	eqnlist = {
		"\n\n\\ Supplementary Vc-inequalities. ["<>S[f]<>"] (See SPIII.A4)\n",
		" dih1 - 0.372 y1 + 0.465 y2 + 0.465 y3 + 0.465 y5 + "<>
			" 0.465 y6 + slack > 4.885" (* case 19 *),
		(* case 21 , assumes dih < 2.12 and heights <=2.26 *) 
		" - vor - 0.06 y2 - 0.06 y3 - 0.185 y5 - 0.185 y6 + "<>
			" slackDih+ slack > -0.9978 ",
		"  vor +0.419351 sol -slackDih - slack < 0.3072 ",
		" pen - slack =0",
		" slack - CqrsA - CqrsB =0"
		};
	eqnlist = Map[(#<>"\n")&,eqnlist]//StringJoin;
	r = RotateLeft[GBLregions[[f]],i-1];
	  StringReplace[eqnlist,
		{"y1"->StringVar["y",r[[1]]], 
		 "y2"->StringVar["y",r[[2]]],
		 "y3"->StringVar["y",r[[4]]],
		 "y4"->StringVar["y",r[[2]],r[[4]]],
		 "y5"->StringVar["y",r[[1]],r[[4]]],
		 "y6"->StringVar["y",r[[1]],r[[2]]],
		 "dih1"->StringVar["dih",f,i],
		 "Dih"->"Dih"<>S[r[[1]]],
		 "sol"->StringVar["sol",f,i],
		 "vor"->StringVar["sig",f,i],
		 "pen"->StringVar["pen",f],
		 "slack"->StringVar["slack","Vc",f],
		 "CqrsA"->"Cqrs"<>S[f 100],
		 "CqrsB"->"Cqrs"<>S[f 100+1]
		}
	]
	];

(* This is not part of WRITEOUT.  Used on archive PM(4,127). 6/1/98. *)
DihHalfVc[f_,i_]:= Module[{eqn,r},
	eqn = {
	"dih2 +0.59 y1 + 0.1 y2 + 0.1 y3 + 0.55 y4 -0.6 y5 -0.12y6>2.6506\n",
	"dih2 +0.35 y1 -0.24 y2 + 0.05y3 + 0.35 y4 -0.72y5 -0.18y6<0.47\n",
	"dih3 +0.59 y1 + 0.1 y3 + 0.1 y2 + 0.55 y4 -0.6 y6 -0.12y5>2.6506\n",
	"dih3 +0.35 y1 -0.24 y3 + 0.05y2 + 0.35 y4 -0.72y6 -0.18y5<0.47\n"
	}//StringJoin;
	r = RotateLeft[GBLregions[[f]],i-1];
	StringReplace[eqn,{
		"y1"-> StringVar["y",r[[1]]],
		"y2"-> StringVar["y",r[[2]]],
		"y3"-> StringVar["y",r[[4]]],
		"y4"-> StringVar["y",r[[2]],r[[4]]],
		"y5"-> StringVar["y",r[[1]],r[[4]]],
		"y6"-> StringVar["y",r[[1]],r[[2]]],
		"dih2"-> StringVar["dih",f,IMod[i+1,4],i],
		"dih3"-> StringVar["dih",f,IMod[i-1,4],i]
			}]
	];

QuadVarRelations:= Module[{fc},
	fc = Select[Range[Length[GBLregions]],Length[GBLregions[[#]]]==4&];
	Map[QuadVarOneRelation,fc]//StringJoin
	];

QuadVarOneRelation[f_]:= Module[{s},
    freeVar=freeVar~Join~Array[StringVar["sig",f,#]&,4];
	freeVar=freeVar~Join~Array[StringVar["sigF",f,#,AugMod[#,4]]&,4];
    freeVar=freeVar~Join~Array[StringVar["sol",f,#]&,4];
	freeVar=freeVar~Join~Array[StringVar["solF",f,#,AugMod[#,4]]&,4];
	freeVar=freeVar~Join~Array[StringVar["dihF",f,#,AugMod[#,4]]&,4];
	freeVar=freeVar~Join~Array[StringVar["dihFx",#,AugMod[#,4],f]&,4];
	freeVar=freeVar~Join~Array[StringVar["dihFx",AugMod[#,4],#,f]&,4];
	freeVar=freeVar~Join~Array[StringVar["dih",f,#,AugMod[#,4]]&,4];
	freeVar=freeVar~Join~Array[StringVar["dih",f,AugMod[#,4],#]&,4];
	freeVar=freeVar~Join~Array[StringVar["Adih",f,#]&,4];
	s = {
    "\n\n\\ QuadVarRelations"<>S[f]<>"\n",
	(* solid on flats relates to solid on quad *)
	StringVar["sol",f]<>"-"<>StringVar["sol",f,1]<>"-"<>
		StringVar["sol",f,3]<>" =0\n",
	StringVar["sol",f]<>"-"<>StringVar["sol",f,2]<>"-"<>
		StringVar["sol",f,4]<>" =0\n",
	StringVar["sig",f]<>"-"<>StringVar["sig",f,1]<>"-"<>
		StringVar["sig",f,3]<>" =0\n",
	StringVar["sig",f]<>"-"<>StringVar["sig",f,2]<>"-"<>
		StringVar["sig",f,4]<>" =0\n",
	(* solid on a flat relates to dihedrals *)
	StringVar["sol",f,1]<>"-",StringVar["dih",f,1],"-",
		StringVar["dih",f,2,1],"-",StringVar["dih",f,4,1],"+pi=0\n",
	StringVar["sol",f,2]<>"-",StringVar["dih",f,2],"-",
		StringVar["dih",f,1,2],"-",StringVar["dih",f,3,2],"+pi=0\n",
	StringVar["sol",f,3]<>"-",StringVar["dih",f,3],"-",
		StringVar["dih",f,2,3],"-",StringVar["dih",f,4,3],"+pi=0\n",
	StringVar["sol",f,4]<>"-",StringVar["dih",f,4],"-",
		StringVar["dih",f,1,4],"-",StringVar["dih",f,3,4],"+pi=0\n",	
	(* two dihedrals make the full angle *)
	StringVar["dih",f,1],"-",StringVar["dih",f,1,2],"-",
		StringVar["dih",f,1,4]," =0\n",
	StringVar["dih",f,2],"-",StringVar["dih",f,2,1],"-",
		StringVar["dih",f,2,3]," =0\n",
	StringVar["dih",f,3],"-",StringVar["dih",f,3,2],"-",
		StringVar["dih",f,3,4]," =0\n",
	StringVar["dih",f,4],"-",StringVar["dih",f,4,3],"-",
		StringVar["dih",f,4,1]," =0\n"}//StringJoin;
 
    (* octahedral constraints *)
    {s,
	"\n\\ octahedral relations:\n",
		StringVar["sigF",f,1,2]<>" +",
           StringVar["sigF",f,2,3]<>" +",
           StringVar["sigF",f,3,4]<>" +",
           StringVar["sigF",f,4,1]<>" -",
            StringVar["sig",f] <> " =0\n",
    	StringVar["dihF",f,1,2]<>" +"<>
           StringVar["dihF",f,2,3]<>" +"<>
           StringVar["dihF",f,3,4]<>" +"<>
           StringVar["dihF",f,4,1]<>" - 2 pi = 0\n",
    Array[StringVar["dihFx",#,AugMod[#,4],f]<>" +"<>
            StringVar["dihFx",#,IMod[#+3,4],f]<>" -"<>
            StringVar["dih",f,#]<>" =0\n"&,4],
    Array[StringVar["dihF",f,#,AugMod[#,4]]<>" +"<>
            StringVar["dihFx",#,AugMod[#,4],f]<>" +"<>
            StringVar["dihFx",AugMod[#,4],#,f]<>" -"<>
            StringVar["solF",f,#,AugMod[#,4]]<>" - pi = 0\n"&,4]
    }//StringJoin
];


(********************* TRIANGLE BRANCH AND BOUND ********************)


(*  Now add stuff for QRtets. 
	Cqrs, Cqrl are slack variables, making these
		equations useless unless the slacks are set to 0. 
	This is the QRtet branch & bound stuff.
	Constants are changed from MathToCplexQuadCase.m, because we do not assume
		y1+y2+y3<2.13. 
	Instead there is a slack variable ht213. Setting ht213=0 activates inequalities
		that hold if all heights are at most 2.13. 
*)
(* Now break qrtets into 2 cases, y4+y5+y6 vs. 2.25 *)
SmallQRtet= { (* partX.cc:cases 201--207 *)
	"\n\n\\ Inequalities for QR tets for which y4+y5+y6>6.25.",
	"\\ Set Cqrl slack variable=0 to activate. (References III.Appendix, VI (Kepler)",
	"\\ The last four assume that the heights are at most 2.13. Set ht213=0 to activate",
 
	" qrs.0 : y4 +y5 +y6 - Cqrs < 6.25 ",
    " qrs.1 : sol + 0.377076 y1 + 0.377076 y2 + 0.377076 y3 - 0.221 y4 - "<>
        " 0.221 y5 - 0.221 y6 + Cqrs > 1.487741 ",
    " qrs.2 : 0.221 y4 + 0.221 y5 + 0.221 y6 - sol + Cqrs > 0.76822 ",
	" qrs.3 : dih1 + 0.34 y2 + 0.34 y3 - 0.689 y4 + 0.27 y5 + 0.27 y6 + Cqrs > 2.29295 ",
	" qrs.4 : dih2 + 0.34 y1 + 0.34 y3 - 0.689 y5 + 0.27 y4 + 0.27 y6 + Cqrs > 2.29295 ",
	" qrs.5 : dih3 + 0.34 y1 + 0.34 y2 - 0.689 y6 + 0.27 y4 + 0.27 y5 + Cqrs > 2.29295 ",
    " qrs.6 : - dih1 + 0.498 y1 + 0.731 y4 - 0.212 y5 - 0.212 y6 + Cqrs > 0.37884 ",
    " qrs.7 : - dih2 + 0.498 y2 + 0.731 y5 - 0.212 y4 - 0.212 y6 + Cqrs > 0.37884 ",
    " qrs.8 : - dih3 + 0.498 y3 + 0.731 y6 - 0.212 y4 - 0.212 y5 + Cqrs > 0.37884 ",
	" qrs.9 : - sig - 0.109 y1 - 0.109 y2 - 0.109 y3 - 0.14135 y4 - "<>
		" 0.14135 y5 - 0.14135 y6 + Cqrs > -1.5574737 ",
    " qrs.10: - sig - 0.419351 sol - 0.2 y1 - 0.2 y2 - 0.2 y3 - 0.048 y4 - "<>
        " 0.048 y5 - 0.048 y6 + Cqrs > -1.77465 ",
    " qrs.11: tau - 0.0845696 y1 - 0.0845696 y2 - 0.0845696 y3 - 0.163 y4 - "<>
        " 0.163 y5 - 0.163 y6 + Cqrs > -1.48542 ",
	(* the next ones assume that the vertices have height at most 2.13 *)
	(* set ht213=0 to activate *)
	" qrs.12 : dih1 + 0.27 y2 + 0.27 y3 - 0.689 y4 + 0.27 y5 + 0.27 y6 +Cqrs+ht213> 2.01295 ",
	" qrs.13 : dih2 + 0.27 y1 + 0.27 y3 - 0.689 y5 + 0.27 y4 + 0.27 y6 +Cqrs+ht213> 2.01295 ",
	" qrs.14 : dih3 + 0.27 y1 + 0.27 y2 - 0.689 y6 + 0.27 y4 + 0.27 y5 +Cqrs+ht213> 2.01295 ",
	" qrs.15 : - sig - 0.14135 y1 - 0.14135 y2 - 0.14135 y3 - 0.14135 y4 - "<>
		" 0.14135 y5 - 0.14135 y6 +Cqrs+ht213 > -1.7515737 "
	};

LargeQRtet= { (* partX.cc:cases 208--212 *)
	"\n\n\\ Inequalities for QR tets for which y4+y5+y6>6.25.",
	"\\ Set Cqrl slack variable=0 to activate. (References III.Appendix, VI (Kepler)",
	"\\ The last three assume that the heights are at most 2.13. Set ht213=0 to activate",
 
    " qrl.0 : y4 +y5 +y6 + Cqrl > 6.25 ",
	" qrl.1 : sol + 0.378 y1 + 0.378 y2 + 0.378 y3 - 0.1781 y4 - "<>
   		" 0.1781 y5 - 0.1781 y6 +Cqrl > 1.761445 ",
	" qrl.2 : - sol - 0.171 y1 - 0.171 y2 - 0.171 y3 + 0.3405 y4 + " <>
   		" 0.3405 y5 + 0.3405 y6 +Cqrl > 0.489145 ",
	" qrl.3 : - sig - 0.1208 y1 - 0.1208 y2 - " <> (* fixed 5/30/98*)
		" 0.1208 y3 - 0.0781 y4 - 0.0781 y5 - 0.0781 y6 + Cqrl > -1.2436 ",
    " qrl.4 : - sig - 0.419351 sol - 0.2 y1 - 0.2 y2 - 0.2 y3 + 0.0106 y4 + "<>
        " 0.0106 y5 + 0.0106 y6 + Cqrl > - 1.40816 ",
	(* The next inequalities assume that the vertices have height at most 2.13.
		Set ht213=0  to activate. *)
    " qrl.5 : sol + 0.356 y1 + 0.356 y2 + 0.356 y3 - 0.1781 y4 - 0.1781 y5 - "<>
        " 0.1781 y6 +Cqrl+ht213 > 1.629445 ",
    " qrl.6 : - sol - 0.254 y1 - 0.254 y2 - 0.254 y3 + 0.3405 y4 + 0.3405 y5 + "<>
        " 0.3405 y6 +Cqrl+ht213 > - 0.008855 ",
    " qrl.7 : - sig - 0.167 y1 - 0.167 y2 - 0.167 y3 - 0.0781 y4 - 0.0781 y5 - "<>
        " 0.0781 y6 + Cqrl+ht213 > -1.51017 "
	};

(* branch & bound qrtet inequalities, activate with slacks *)
QRText:= Module[{g,x,Combine,QRTextOne,QRtetsub,f},
	Combine[x_]:= Apply[StringJoin,Map[(#<>"\n")&,x]];
	QRTextOne[f_]:=
		 StringReplace[
			Join[SmallQRtet,LargeQRtet],QRtetsub[f]]//Combine;
 
	QRtetsub[f_]:=
		{"qrl" -> "qrl"<>S[f],
		 "qrs" -> "qrs"<>S[f],
		 "sol" -> StringVar["sol",f],
		 "dih1"-> StringVar["dih",f,1],
		 "dih2"-> StringVar["dih",f,2],
		 "dih3"-> StringVar["dih",f,3],
		 "y1"->   StringVar["yG",f,1],
		 "y2"->   StringVar["yG",f,2],
		 "y3"->   StringVar["yG",f,3],
		 "y4"->   StringVar["yG",f,2,3],
		 "y5"->   StringVar["yG",f,1,3],
		 "y6"->   StringVar["yG",f,1,2],
		 "sig" -> StringVar["sig",f],
		 "ht213" -> StringVar["ht",f],
		 "tau" -> StringVar["tau",f]};
		
	{"\n\n\n\\ Branch and Bound Inequalities for QRTets,\n\n",
	 g= Select[Range[Length[GBLregions]],Length[GBLregions[[#]]]==3&]; 
	 Map[QRTextOne,g],
	 Array[StringVar["ht",#]<>" -ht213 < 0\n"&,Length[GBLregions]]
	}
	//StringJoin
	];

(************ GENERATE BRANCH AND BOUND CASES ******************)

optText[file_,added_]:=
	{
	"read ",file," lp\n",
	"add \n",added,
	"end\n",
	"opt\n"
	}//StringJoin;


Branch[list_,i_]:= Module[{digits,j},
	digits = Table[1+Mod[Floor[(i-1)/2^j],2],{j,0,Length[list]-1}];
	Array[{"Cqr",Part[{"s","l"},digits[[#]]],S[list[[#]]],"=0\n"}&,
		Length[list]]//StringJoin
	];

optBranch[file_,list_]:= 
	Array[optText[file,Branch[list,#]]&,2^Length[list]]//StringJoin;


