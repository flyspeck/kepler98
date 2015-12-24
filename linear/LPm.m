(* suppose we have an LP file "basefile" that starts out
	MAXIMIZE X ST ....
   This procedure adds the eqn 
		yvar - X = 0
	and maximizes with cplex.
*)
LPmFile=1;

(*
LPmax[basefile_String,yvar_String]:
LPmUpdateAdih[basefile_String,yvar_String,dihvar_String,Adihvar_String]:
LPmUpdateAdih[basefile_String,yvar_String,dihvar_String,Adihvar_String,
	add_String]:
*)

MakeRandom := (SeedRandom; MakeRandom := Random[Integer, {0, 10^9}]; 
    MakeRandom);

LPmax[basefile_String,yvar_String]:= LPmax[basefile,yvar,""];
LPmax[basefile_String,yvar_String,add_String]:= 
  Module[{stream,x,clopen,out,sedfile,ext},
	ext = S[MakeRandom];
	stream = OpenWrite["X/LPmaxExec"<>ext];
	clopen[x_]:= Close[OpenWrite[x]];
	out[x_]:= WriteString[stream,x<>"\n"];
	out["set logfile X/LPmax.log"<>ext];
	out["set preprocessing presolve no"];
	out["set output results n X/y.log"];
	out["set output dialog n X/x.log"];
	out["read "<>basefile<>" lp"];
	out["change delete constraints sigsumX"];
	out["add "];
	If[add!="",out[add]];
	out["sigsum > 0.4429"];
	out[yvar<>" - X =0"];
	out["end\nopt\nquit"];
	Close[stream];
	clopen["X/LPmax.log"<>ext]; 
	clopen["X/x.log"]; clopen["X/y.log"];
	Run["cplex < X/LPmaxExec"<>ext];
	clopen["X/LPmaxValues"<>ext];
	sedfile = "/afs/math.lsa.umich.edu/group/fac/hales/pub/html/countdown/"<>
        "linear/SHORT/getAddValues";
	Run["sed -f "<>sedfile<> " X/LPmax.log"<>ext<>" > X/LPmaxValues"<>ext];
	stream = OpenAppend["X/LPmaxValues"<>ext];
    WriteString[stream,"0}~Drop~-1;"];
    Close[stream];
	values={{"UNKNOWN"}};
	ToExpression["<< X/LPmaxValues"<>ext]; (* set values *)
	Run["rm X/LPmax.log"<>ext];
	Run["rm X/LPmaxValues"<>ext];
	Run["rm X/LPmaxExec"<>ext];
	Print[values//First//First];
	values//First//First
	];

LPdisplay[basefile_String]:= LPdisplay[basefile,""];
LPdisplay[basefile_String,add_String]:= 
  Module[{stream,x,clopen,out,sedfile,ext},
	ext = S[MakeRandom];
	stream = OpenWrite["X/LPmaxExec"<>ext];
	clopen[x_]:= Close[OpenWrite[x]];
	out[x_]:= WriteString[stream,x<>"\n"];
	out["set logfile X/LPmax.log"<>ext];
	out["set preprocessing presolve no"];
	out["set output results n X/y.log"];
	out["set output dialog n X/x.log"];
	out["read "<>basefile<>" lp"];
	out["add "];
	If[add!="",out[add]];
	out["end\nopt\nquit"];
	Close[stream];
	clopen["X/LPmax.log"<>ext]; 
	clopen["X/x.log"]; clopen["X/y.log"];
	Run["cplex < X/LPmaxExec"<>ext];
	clopen["X/LPmaxValues"<>ext];
	sedfile = "/afs/math.lsa.umich.edu/group/fac/hales/pub/html/countdown/"<>
        "linear/SHORT/getAddValues";
	Run["sed -f "<>sedfile<> " X/LPmax.log"<>ext<>" > X/LPmaxValues"<>ext];
	stream = OpenAppend["X/LPmaxValues"<>ext];
    WriteString[stream,"0}~Drop~-1;"];
    Close[stream];
	values={{"UNKNOWN"}};
	ToExpression["<< X/LPmaxValues"<>ext];
	Run["rm X/LPmax.log"<>ext];
	Run["rm X/LPmaxValues"<>ext];
	Run["rm X/LPmaxExec"<>ext];
	Print[values//First//First];
	values//First//First
	];


(* sample call. LPmUpdateAdih["SHORT/cplexE.lp174.F16.C2",
		StringVar["y",8],
		StringVar["dihE",16,8,4,12],
		StringVar["AdihE",16,8,4,12]] 
	This uses GBLregions[[16]]=={2, 7, 11, 12, 8, 4, 3},
		notices that 8 occurs between 4 and 12. *)


LPmUpdateAdih[basefile_String,yvar_String,dihvar_String,Adihvar_String]:= 
	Module[{hmax,dihmax,dihmin,},
	If[LPdisplay[basefile]<0.4429,
			Return["\n\\ UpdateAdih Infeasible: sigsum<0.4429\n"]];
	hmax = LPmax[basefile,yvar]/2. + 10.^-6; (* bug fixed 4/24/98 *)
	dihmin= -LPmax[basefile,"-"<>dihvar]- 10.^-6;
	dihmax= LPmax[basefile,dihvar] +10.^-6;
	"\\ hmax,dihmin,dihmax "<>S[hmax]<>" "<>S[dihmin]<>" "<>S[dihmax]<>"\n"<>
	AdihText[yvar,dihvar,Adihvar,hmax,dihmin,dihmax]
	];

LPmUpdateAdih[basefile_String,yvar_String,dihvar_String,Adihvar_String,
	add_String]:= 
	Module[{hmax,dihmax,dihmin,},
	If[LPdisplay[basefile,add]<0.4429,
			Return["\n\\ UpdateAdih Infeasible: sigsum<0.4429\n"]];
	hmax = LPmax[basefile,yvar,add]/2. + 10.^-6; (* bug fixed 4/24/98 *)
	dihmin= -LPmax[basefile,"-"<>dihvar,add]- 10.^-6;
	dihmax= LPmax[basefile,dihvar,add] +10.^-6;
	"\\ hmax,dihmin,dihmax "<>S[hmax]<>" "<>S[dihmin]<>" "<>S[dihmax]<>"\n"<>
	AdihText[yvar,dihvar,Adihvar,hmax,dihmin,dihmax]
	];

