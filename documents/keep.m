Bisect[{v1_, v2_}, eps_] :=
  Module[{y1, y2, ym, x1, x2, xm, reps, f},
   {x1, x2} = Sort[N[{v1, v2}]]; f = Bisectf; xm := (x1 + x2)/2; y1 = f[x1];
    y2 = f[x2]; ym := f[xm]; If[y1*y2 > 0, Return[]];
    reps = Ceiling[Log[(x2 - x1)/eps]/Log[2]];
    Do[ym = f[xm]; If[ym*y1 < 0, x2 = xm, x1 = xm]; Null, {i, 1, reps}];
    {x1, x2}];
 

