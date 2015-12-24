s8 = Sqrt[8];
s3 = Sqrt[3];

w[0]= {0,0,0};

w[1]= {2,0,0};

w[2]= {1,s3,0};

w[3]= {-1,s3,0};

w[4]= {-2,0,0};

w[5]= {-1,-s3,0};

w[6]= {1,-s3,0};

w[7]= {1,s3/3, s8/s3};

w[8]= {-1,s3/3, s8/s3};

w[9]= {0,-2/s3, s8/s3};

w[10]= {0,2/s3, -s8/s3};

w[11]= {-1,-s3/3,-s8/s3};

w[12]={1,-s3/3,-s8/s3};


v[i_,theta_]:= Cos[theta] w[i] + Sin[theta] w[bijection[[i]]];

bijection = {10, 8, 11, 9, 12, 7, 6, 2, 4, 1, 3, 5};

$SpherePrecision=Infinity;

edge = dist[v[1,theta],v[2,theta]]//Simplify;

diag = dist[v[1,theta],v[10,theta]]//Simplify;

dodecsub = {theta -> 1/2 ArcSin[1/Sqrt[5]]};  (* diag=edge here *)

volume = VolAn[2,2,2,edge,edge,edge] 8 + 12 VolAn[2,2,2,diag,edge,edge];

cvolume = volume + 24 L[edge] 2;

ccvolume = cvolume + 0.006 6 (s8-diag);
