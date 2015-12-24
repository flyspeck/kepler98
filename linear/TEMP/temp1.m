(* The vorVc bound *)
AddMoreExcept:= Module[{multiplier,
	r,i,j,exc,phi0,pt,doct4,drastic,str,ln,Y},
	Y[f_,p1_]:= StringVar["y",GBLregions[[f,p1]]];
	Y[f_,p1_,p2_]:= StringVar["y",GBLregions[[f,p1]],GBLregions[[f,p2]]];
	exc = Select[Range[1,GBLregions//Lg],Length[GBLregions[[#]]]>4&];
	phi0 = -0.5666365478933329;
	pt = 0.05537364566846414;
	doct4 = 2.883611798069859;
	(* drastic = penalties for switching to VorVc *)
	drastic = {0,0,0,0,5*0.01561*pt,6*0.01561*pt,7*0.01561*pt,8*0.01561*pt};
	str = "\n\n";
 
	(* vorVc upper bound *)
	Do[str = str<>"\n";
	   r = exc[[i]];
		ln = GBLregions[[r]]//Length;
		str = str <> "-" <> StringVar["sig",r];
		str = str <> SR[phi0] <> " " <> StringVar["sol",r] <> " ";
		Do[str = str <> "+" <> StringVar["Adih",r,j],
			{j,1,ln}];
		Do[str = str <> "\n"<> S[-doct4] <> " "<> StringVar["quo",r,j,IMod[j+1,ln]];
			,{j,1,ln}];
		Do[str = str <> "\n"<> S[-doct4] <> " "<> StringVar["quo",r,IMod[j+1,ln],j];
			,{j,1,ln}];
		str = str<>"\n";
		str = str <> " > " <> S[-drastic[[ln]]] <> "\n";
		, {i,1,Length[exc]}
		];
 
	(* Quoin *)
	Do[r = exc[[i]];
		ln = GBLregions[[r]]//Length;
 
		Do[
		str = str <> StringVar["quo",exc[[i]],j,IMod[j+1,ln]] <>
		" +0.00758 "<> Y[r,j] <> 
		" +0.0115 " <> Y[r,IMod[j-1,ln]] <>
		" +0.0115 " <> Y[r,j,IMod[j-1,ln]] <>
		" > " <> S[0.00217 + 2 0.00758 + 4 0.0115] <> "\n";
		, {j,1,ln}];
 
		Do[
		str = str <> StringVar["quo",exc[[i]],j,IMod[j+1,ln]] <>
		" +0.00758 "<> Y[r,j] <> 
		" +0.0115 " <> Y[r,IMod[j+1,ln]] <>
		" +0.0115 " <> Y[r,j,IMod[j+1,ln]] <>
		" > " <> S[0.00217 + 2 0.00758 + 4 0.0115] <> "\n";
		, {j,1,ln}];
		, {i,1,Length[exc]}];
 
	(* Adih inequalities *)
	multiplier = 0.109691511444153; (* f[1];  
		(alpha f[h] is term in vorVc. ),
		 f[h_]:= (1 - h/t0)*(phi[h, t0] - phi[t0, t0]) *)
	Do[r = exc[[i]]; 
		ln = GBLregions[[r]]//Length;
		Do[
		str = str <> StringVar["Adih",r,j] <> S[-multiplier] <> " " <>
			StringVar["dih",r,j] <> " " <>
			SR[multiplier/(1.255-1.0) alphaMin[r,j]/2] <> " " <>
				StringVar["y",GBLregions[[r,j]]] <> " < " <>
			S[multiplier/(1.255-1.0) alphaMin[r,j]] <> "\n";
		str = str <> StringVar["Adih",r,j] <> " " <>
			SR[multiplier /(1.255-1)/2 alphaMax[r,j]] <> " " <>
			StringVar["y",GBLregions[[r,j]]] <> " < " <>
			S[multiplier/(1.255-1) 1.255 alphaMax[r,j]] <> "\n\n";
		(* finish *)
		,{j,1,ln}];
		, {i,1,Length[exc]}];
	str
	];
