
Make[t_]:= (
	{v0,v1,v2}= SimplexCoordinates[t,2,2,t,2,2];
	v3 = FarFrom[v1,Vertex[{v0,2},{v2,2},{NullVector,t}]];
	v4 = FarFrom[v2,Vertex[{v1,t},{v3,2},{NullVector,2}]];
    v5 = FarFrom[v1,Vertex[{v2,t},{v4,t},{NullVector,2}]];
	dist[v0,v4]
	);

MakeQ[t_]:= (
	{v0,v1,v2} = SimplexCoordinates[2Sqrt[2],2,2,t,2,2];
	v3 = FarFrom[v1,Vertex[{v0,2},{v2,t},{NullVector,2}]];
	v4 = FarFrom[v2,Vertex[{v1,t},{v3,t},{NullVector,2}]];
	dist[v4,v0]
	);

Make[y1_,y2_,y3_,y4_,y5_]:= (
	{v0,v1,v2}= SimplexCoordinates[edgeMax,2,2,edgeMax,y2,y1];
	v3 = FarFrom[v1,Vertex[{v0,y3},{v2,edgeMax},{NullVector,2}]];
	v4 = FarFrom[v2,Vertex[{v1,edgeMax},{v3,edgeMax},{NullVector,2}]];
	v5 = CloseTo[v0,Vertex[{v1,y4},{v4,y5},{NullVector,edgeMax}]];
	{dist[v5,v0],dist[v5,v3]}
	);

MakeX[u_]:= (
	{v0,v1,v2}=SimplexCoordinates[2,2,2,u,edgeMax,edgeMax];
	v4 = FarFrom[v0,Vertex[{v2,edgeMax},{v1,edgeMax},{NullVector,2}]];
	v5 = CloseTo[v2,Vertex[{v4,2},{v1,2},{NullVector,edgeMax}]];
	v6 = CloseTo[v5,Vertex[{v0,edgeMax},{v2,2},{NullVector,edgeMax}]];
	{dist[v5,v6],dist[v5,v0],dist[v5,v2],dist[v6,v1],dist[v6,v4]}
	);

CloseTo[v_,pair_]:= Complement[pair,{FarFrom[v,pair]}][[1]];

tDef[x_]:= -2Pi + Dihedral[x,2,2,x,2,x]+Dihedral[x,2,2,x,2,2]+
	Dihedral[x,2,x,2,2,2]+Dihedral[x,2,x,x,2,x];

edgeMaxFn[x_]:= -2Pi + NewDihedral[x,x,2,2,x,2]+NewDihedral[x,2,2,x,x,2]+
	NewDihedral[x,2,2,x,2,2]+NewDihedral[x,x,2,x,2,2];

(* approx value of edgeMax =  2.5751835464 *)
edgeMax = 2.558913201377738 (* edgeMaxFn *)

(* This should be an excellent upper bound on Gamma for qrtet (with edgeMax) *)


PseudoGamma[y1_,y2_,y3_,y4_,y5_,y6_]:=
	pt - 0.17 (y1+y2+y3+y4+y5+y6-12) + 0.122 (
	(y1-2)^2 + (y2-2)^2 + (y3-2)^2 + (y4-2)^2 + (y5-2)^2 + (y6-2)^2)

NewVorVc[x__]:= NewTrunc[VorVc,x,edgeMax/2];
NewTauVc[x__]:= NewTrunc[tauVc,x,edgeMax/2];


MakeY[u1_,u2_]:= (
	u = Interval[{u1,u2}];
	{v0,v1,v2}= SimplexCoordinates[2,2,2,u,edgeMax,edgeMax];
	v3 = FarFrom[v0,Vertex[{v1,edgeMax},{v2,edgeMax},{NullVector,2}]];
	v4 = CloseTo[v0,Vertex[{v2,2},{v3,2},{NullVector,edgeMax}]];
	v5 = CloseTo[v4,Vertex[{v0,edgeMax},{v1,2},{NullVector,edgeMax}]];
	{dist[v4,v5],dist[v4,v1]}
	);
