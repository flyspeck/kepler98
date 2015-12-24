

(*
<< QUAD18/lquad.m;
ConfigurationList = lpent;
*)

<< MathToCplex.m;


exInit:=(
Print[" -- Exceptional Cases -- "];
<< BASIC/lpent.m;
<< BASIC/lhex.m;
<< BASIC/lhept.m;
<< BASIC/loct.m;
ConfigurationList={};
<< LPm.m;
<< BASIC/l.sum; (* contains archive function *)
<< MathToCplexExcept.m;
);
(* exInit;  *)

pentInit:=(<< MathToCplexPent.m);
(* pentInit;  *)

quadInit:=(
<< QUAD18/lquad.m;
<< LPm.m;
ConfigurationList = lquad;
(* Quadgen:= Table[WRITEOUT[i]; Print[i],{i,1,Length[ConfigurationList]}]; *)
);

Print[" -- CPLEX Initialized -- "];
!seticon MToCplex
!settitle MathToCplex
Fix[x_] := Edit[Definition[x]]
SetAttributes[Fix,HoldFirst];
(Unprotect[In,Out]; Format[In]=MathToCplexIn;
    Format[Out]=MathToCplexOut; Protect[In,Out];)

