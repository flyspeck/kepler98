Ac := {0, 2/Sqrt[3]}
Bc := {1, -Sqrt[3]^(-1)}
cos[i_, j_] := Cos[2*Pi*{x, y} . (i*Ac + j*Bc)]
size[i_, j_] := 3*norm[i*Ac + j*Bc]
norm[x_] := x . x
ff[x1_, y1_] := sum /. {x -> x1, y -> y1}
derivx := D[sum, x]
derivy := D[sum, y]
ac := {1, 0}
bc := {1/2, Sqrt[3]/2}
fg[i_, j_] := ff @@ ((ac*i)/4 + (bc*j)/4)
fg[i_, j_] := ff @@ ((ac*i)/4 + (bc*j)/4)
 
fg[i_, j_, f_] := f @@ ((ac*i)/4 + (bc*j)/4)
eval[x1_, y1_, f_] := f /. {x -> x1, y -> y1}
cos[i_, j_] := Cos[2*Pi*{x, y} . (i*Ac + j*Bc)]
 
cos[i1_, j1_, i_, j_] := Cos[2*Pi*(i1*ac + j1*bc) . (i*Ac + j*Bc)]
sum[i1_, j1_, max_] := 
  Sum[v = isjs[k][[i]]; x[size @@ v]*cos[i1, j1, v[[1]], v[[2]]], 
    {k, 1, max}, {i, 1, Length[isjs[k]]}] + 1
cos[i_, j_] := Cos[2*Pi*{x, y} . (i*Ac + j*Bc)]
 
cos[i1_, j1_, i_, j_] := Cos[2*Pi*(i1*ac + j1*bc) . (i*Ac + j*Bc)]
isjs[0] = {{0, 0}}
 
isjs[1] = {{0, 1}, {1, 0}, {1, 1}}
 
isjs[2] = {}
 
isjs[3] = {{1, -1}, {1, 2}, {2, 1}}
 
isjs[4] = {{0, 2}, {2, 0}, {2, 2}}
 
isjs[5] = {}
 
isjs[6] = {}
 
isjs[7] = {{1, -2}, {1, 3}, {2, -1}, {2, 3}, {3, 1}, {3, 2}}
 
isjs[8] = {}
 
isjs[9] = {{0, 3}, {3, 0}, {3, 3}}
 
isjs[10] = {}
 
isjs[11] = {}
 
isjs[12] = {{2, -2}, {2, 4}, {4, 2}}
 
isjs[13] = {{1, -3}, {1, 4}, {3, -1}, {3, 4}, {4, 1}, {4, 3}}
 
isjs[14] = {}
 
isjs[15] = {}
 
isjs[16] = {{0, 4}, {4, 0}, {4, 4}}
 
isjs[17] = {}
 
isjs[18] = {}
 
isjs[19] = {{2, -3}, {2, 5}, {3, -2}, {3, 5}, {5, 2}, {5, 3}}
 
isjs[20] = {}
 
isjs[21] = {{1, -4}, {1, 5}, {4, -1}, {4, 5}, {5, 1}, {5, 4}}
 
isjs[22] = {}
 
isjs[23] = {}
 
isjs[24] = {}
 
isjs[25] = {{0, 5}, {5, 0}, {5, 5}}
 
isjs[26] = {}
 
isjs[27] = {{3, -3}, {3, 6}, {6, 3}}
 
isjs[28] = {{2, -4}, {2, 6}, {4, -2}, {4, 6}, {6, 2}, {6, 4}}
 
isjs[29] = {}
 
isjs[30] = {}
 
isjs[31] = {{1, -5}, {1, 6}, {5, -1}, {5, 6}, {6, 1}, {6, 5}}
 
isjs[32] = {}
 
isjs[33] = {}
 
isjs[34] = {}
 
isjs[35] = {}
 
isjs[36] = {{0, 6}, {6, 0}, {6, 6}}
 
isjs[37] = {{3, -4}, {3, 7}, {4, -3}, {4, 7}, {7, 3}, {7, 4}}
 

xx5sub = {x[2] -> 0, x[5] -> 0, x[6] -> 0, x[8] -> 0, x[10] -> 0, x[11] -> 0, 
   x[14] -> 0, x[15] -> 0, x[17] -> 0, x[18] -> 0, x[20] -> 0, x[21] -> 0, 
   x[22] -> 0, x[23] -> 0, x[24] -> 0, x[25] -> 0, x[26] -> 0, x[27] -> 0, 
   x[28] -> 0, x[29] -> 0, x[30] -> 0, x[31] -> 0, x[32] -> 0, x[33] -> 0, 
   x[34] -> 0, x[35] -> 0, x[36] -> 0}
x5sub = {x[1] -> 1.945731999999999, x[2] -> 0, x[3] -> 1.575004, 
   x[4] -> 1.394052999999999, x[5] -> 0, x[6] -> 0, x[7] -> 0.844747, 
   x[8] -> 0, x[9] -> 0.5655969999999999, x[10] -> 0, x[11] -> 0, 
   x[12] -> 0.30085, x[13] -> 0.2124979999999999, x[14] -> 0, x[15] -> 0, 
   x[16] -> 0.05426799999999999, x[17] -> 0, x[18] -> 0, x[19] -> 0.020175, 
   x[20] -> 0, x[21] -> 0, x[22] -> 0, x[23] -> 0, x[24] -> 0, x[25] -> 0, 
   x[26] -> 0, x[27] -> 0, x[28] -> 0, x[29] -> 0, x[30] -> 0, x[31] -> 0, 
   x[32] -> 0, x[33] -> 0, x[34] -> 0, x[35] -> 0, x[36] -> 0, 
   x[37] -> 0.004826999999999999}
sols = {{x[28] -> 0, x[31] -> 0, x[1] -> 2 - x[16], x[3] -> 2 - 2*x[13], 
    x[4] -> 2 - x[9] - 2*x[19], x[7] -> (2 - x[12] - 2*x[37])/2}}
sin[i1_, j1_, i_, j_] := Sin[2*Pi*(i1*ac + j1*bc) . (i*Ac + j*Bc)]
distR2[x_, y_, x1_, y1_] := 
  Module[{u, v}, u = x - x1 + (y - y1)/2; v = y - y1; u*u + (v*v*3)/4]
nearLattice[x_, y_, radius_] := 
  Module[{fr2 = 4*radius^2}, distR2[x, y, 0, 0] <= fr2 || 
    distR2[x, y, 1, 0] <= fr2 || distR2[x, y, 0, 1] <= fr2 || 
    distR2[x, y, 1, 1] <= fr2]
sumx[i1_, j1_, max_] := 
  Sum[v = isjs[k][[i]]; x[size @@ v]*sin[i1, j1, v[[1]], v[[2]]]*
     (v[[1]]*Ac + v[[2]]*Bc)[[1]], {k, 1, max}, {i, 1, Length[isjs[k]]}]
 
sumx[i1_, j1_, rad_, max_] := 
  If[nearLattice[i1, j1, rad], 0, sumx[i1, j1, max]]
sumy[i1_, j1_, max_] := 
  Sum[v = isjs[k][[i]]; x[size @@ v]*sin[i1, j1, v[[1]], v[[2]]]*
     (v[[1]]*Ac + v[[2]]*Bc)[[2]], {k, 1, max}, {i, 1, Length[isjs[k]]}]
 
sumy[i1_, j1_, rad_, max_] := 
  If[nearLattice[i1, j1, rad], 0, sumy[i1, j1, max]]
size[i_, j_] := (3*norm[i*Ac + j*Bc])/4
q[i_, j_] := If[i == 0 && j <= 0, 0, r[i, j]]
r[i_, j_] := If[i^2 + j^2 - i*j > 49, 0, size[i, j]]
isjs[0] = {}
 
isjs[1] = {{0, 1}, {1, 0}, {1, 1}}
 
isjs[2] = {}
 
isjs[3] = {{1, -1}, {1, 2}, {2, 1}}
 
isjs[4] = {{0, 2}, {2, 0}, {2, 2}}
 
isjs[5] = {}
 
isjs[6] = {}
 
isjs[7] = {{1, -2}, {1, 3}, {2, -1}, {2, 3}, {3, 1}, {3, 2}}
 
isjs[8] = {}
 
isjs[9] = {{0, 3}, {3, 0}, {3, 3}}
 
isjs[10] = {}
 
isjs[11] = {}
 
isjs[12] = {{2, -2}, {2, 4}, {4, 2}}
 
isjs[13] = {{1, -3}, {1, 4}, {3, -1}, {3, 4}, {4, 1}, {4, 3}}
 
isjs[14] = {}
 
isjs[15] = {}
 
isjs[16] = {{0, 4}, {4, 0}, {4, 4}}
 
isjs[17] = {}
 
isjs[18] = {}
 
isjs[19] = {{2, -3}, {2, 5}, {3, -2}, {3, 5}, {5, 2}, {5, 3}}
 
isjs[20] = {}
 
isjs[21] = {{1, -4}, {1, 5}, {4, -1}, {4, 5}, {5, 1}, {5, 4}}
 
isjs[22] = {}
 
isjs[23] = {}
 
isjs[24] = {}
 
isjs[25] = {{0, 5}, {5, 0}, {5, 5}}
 
isjs[26] = {}
 
isjs[27] = {{3, -3}, {3, 6}, {6, 3}}
 
isjs[28] = {{2, -4}, {2, 6}, {4, -2}, {4, 6}, {6, 2}, {6, 4}}
 
isjs[29] = {}
 
isjs[30] = {}
 
isjs[31] = {{1, -5}, {1, 6}, {5, -1}, {5, 6}, {6, 1}, {6, 5}}
 
isjs[32] = {}
 
isjs[33] = {}
 
isjs[34] = {}
 
isjs[35] = {}
 
isjs[36] = {{0, 6}, {6, 0}, {6, 6}}
 
isjs[37] = {{3, -4}, {3, 7}, {4, -3}, {4, 7}, {7, 3}, {7, 4}}
 
isjs[38] = {}
 
isjs[39] = {{2, -5}, {2, 7}, {5, -2}, {5, 7}, {7, 2}, {7, 5}}
 
isjs[40] = {}
 
isjs[41] = {}
 
isjs[42] = {}
 
isjs[43] = {{1, -6}, {1, 7}, {6, -1}, {6, 7}, {7, 1}, {7, 6}}
 
isjs[44] = {}
 
isjs[45] = {}
 
isjs[46] = {}
 
isjs[47] = {}
 
isjs[48] = {{4, -4}, {4, 8}, {8, 4}}
 
isjs[49] = {{0, 7}, {3, -5}, {3, 8}, {5, -3}, {5, 8}, {7, 0}, {7, 7}, {8, 3}, 
   {8, 5}}
size[i_, j_] := Simplify[(3*norm[i*Ac + j*Bc])/4]
AA[y_] := A[y/2]*2.*Pi
A3[y_, ymax_] := ((AA[ymax] - AA[2])*(y - 2))/(ymax - 2) + AA[2]
aDihedral[y1_, y2_, y3_, y4_, y5_, y6_] := 
  A[y1/2]*Dihedral[y1, y2, y3, y4, y5, y6]
seanT := 1.25841
seanM = 0.42755
