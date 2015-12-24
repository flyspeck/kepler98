(* suppose we have an LP file "basefile" that starts out
	MAXIMIZE X ST ....
   This procedure adds the eqn 
		yvar - X = 0
	and maximizes with cplex.
*)
LPm`File=1;

(*
LPm`max[basefile_String,yvar_String]:
LPm`UpdateAdih[basefile_String,yvar_String,dihvar_String,Adihvar_String]:
*)

LPmax[basefile_,yvar_]:= Module[{stream,x,clopen,out,sedfile},
	stream = OpenWrite["/tmp/seanLPmaxExec"];
	clopen[x_]:= Close[OpenWrite[x]];
	out[x_]:= WriteString[stream,x<>"\n"];
	out["set logfile /tmp/seanLPmax.log"];
	out["set preprocessing presolve no"];
	out["set output results n /tmp/seanY.log"];
	out["set output dialog n /tmp/seanX.log"];
	out["read "<>basefile<>" lp"];
	out["add "];
	out[yvar<>" - X =0"];
	out["end\nopt\nquit"];
	Close[stream];
	clopen["/tmp/seanLPmax.log"]; clopen["/tmp/seanX.log"]; clopen["/tmp/seanY.log"];
	Run["cplex < /tmp/seanLPmaxExec"];
	clopen["/tmp/seanLPmaxValues"];
	sedfile = "/afs/math.lsa.umich.edu/group/fac/hales/pub/sphere/"<>
        "cplex/SHORT/getAddValues";
	Run["sed -f "<>sedfile<> " /tmp/seanLPmax.log > /tmp/seanLPmaxValues"];
	stream = OpenAppend["/tmp/seanLPmaxValues"];
    WriteString[stream,"0}~Drop~-1;"];
    Close[stream];
	<< /tmp/seanLPmaxValues;
	values
	];
LPm`max[basefile_String,yvar_String]:= LPmax[basefile,yvar]//First//First;


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

(* sample call. ReWorkAdih["SHORT/cplexE.lp174.F16.C2",
		StringVar["y",8],
		StringVar["dihE",16,8,4,12],
		StringVar["AdihE",16,8,4,12]] 
	This uses GBLregions[[16]]=={2, 7, 11, 12, 8, 4, 3},
		notices that 8 occurs between 4 and 12. *)


LPm`UpdateAdih[basefile_String,yvar_String,dihvar_String,Adihvar_String]:= 
	Module[{hmax,dihmax,dihmin,},
	Infeasible=0;
	If[LPm`max[basefile,"sigsum"]<0.4429,
			Return["\n\\ UpdateAdih Infeasible: sigsum<0.4429\n"]];
	hmax = LPm`max[basefile,yvar]/2. + 10.^-6; (* bug fixed 4/24/98 *)
	dihmin= -LPm`max[basefile,"-"<>dihvar]- 10.^-6;
	dihmax= LPm`max[basefile,dihvar] +10.^-6;
	"\\ hmax,dihmin,dihmax "<>S[hmax]<>" "<>S[dihmin]<>" "<>S[dihmax]<>"\n"<>
	AdihText[yvar,dihvar,Adihvar,hmax,dihmin,dihmax]
	];
