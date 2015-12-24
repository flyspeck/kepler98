
(* cad example 1 *)
(* show there is no point in interior of tetrahedron with sides le sqrt8 *)

v1 = {x1,y1,z1};
v2 = {x2,y2,z2};
v3 = {x3,y3,z3};
v4 = {x4,y4,z4};

p1 = v1.v1 -4  (* ge 0 *)
p2 = v2.v2 - 4
p3 = v3.v3 - 4
p4 = v4.v4 -4
p5 = 8 - (v1-v2).(v1-v2);
p6 = 8 - (v1-v3).(v1-v3);
p7 = 8 - (v1-v4).(v1-v4);
p8 = 8 - (v2-v3).(v2-v3);
p9 = 8 - (v2-v4).(v2-v4);
p10= 8 - (v3-v4).(v3-v4);
z4 = 0; y4=0;

poly = {p1,p2,p3,p4,p5,p6,p7,p8,p9,p10};

Deg[A_,x_]:= Map[degree[#,x]&,A]//Max;
DegList[A_,y_]:= Map[{#,Deg[A,#]}&,y];
vars = {x1,y1,z1,x2,y2,z2,x3,y3,z3,x4};

calcz:=(
poly1 = PROJ[poly,z3];
poly2 = PROJ[poly1,z2];
poly3 = PROJ[poly2,z1];
Length[poly3]);

calcy:= (
poly4 = PROJ[poly3,y3];
poly5 = PROJ[poly4,y2];
poly6 = PROJ[poly6,y1];
Length[poly6]
);

