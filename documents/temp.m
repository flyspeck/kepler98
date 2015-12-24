f[y1_, y2_, y3_, y4_, y5_, y6_] :=
  Module[{cf = 0.08, slope = 0.4},
   (VolAn[y1, y2, y3, y4, y5, y6] - VolAn[2, 2, 2, 2, 2, 2] +
    (0.04*(6 - y1 - y2 - y3)))-slope*
     (Solid[y1, y2, y3, y4, y5, y6] - Solid[2, 2, 2, 2, 2, 2])]


fq[y1_, y2_, y3_, y4_, y5_, y6_] :=
  Module[{cf = 0.08, slope = 0.4, num, den},
   (VolAn[y1, y2, y3, y4, y5, y6] - VolAn[2, 2, 2, Sqrt[8], 2, 2] +
     cf*(y5 + y6 - 4) - 0.04*(y1 + y2 + y3 - 6))  -slope* 
    (Solid[y1, y2, y3, y4, y5, y6] - Solid[2, 2, 2, Sqrt[8], 2, 2])] 
