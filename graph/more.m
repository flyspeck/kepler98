(* java format conversion *)
set[i_] := jformat[i-1] = (Initialize[i]; arrx);
arrx := Flatten[Join[{0, Length[GBLregions]}, arr]];
arr := ({Length[#1], -1+ #1} & ) /@ GBLregions;


(* invariant testing 
relabel vertices for a trial
perm:={1->3,...}... 
f3 = finalX[[2,3]];
rewrite:= Table[{i/.perm,f3[[i,1]],f3[[i,2]]/.perm,f3[[i,3]]/.perm},{i,1,Length[f3]}]//Sort;

Map[Apply[ptcluster,#[[{2,3,4}]] ]&,%]; *)
  (* reinject into final... *)
