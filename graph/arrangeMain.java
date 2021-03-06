import java.util.*;

public class arrangeMain extends arrangePlus {

	static { System.out.print(" -- arrangeMain $Revision: 1.3 $ -- \n"); }

	/* Push onto the stack all possible ngons constructed in face F
		 using the edge terminating at vertex A.  */
	 
static void generatePolygon(int ngon,arrangePlus a,vertex A,
	face F,Stack tentative,Stack terminal,parameterUse p)
        {
													debug.d1(ngon,a,A,F);
	int offset = F.indexOf(A);				       debug.d2(offset);
	int outer = F.size();							debug.d3(outer);
	boolean nix[][] = new boolean[outer][outer];
	init_nix(nix,offset,a,F);

	/* genpoly loop */
	iterator it = new iterator(ngon,outer);
	genpolyloop: while (it.next())
		{
		int poly[];
		int newCount[];

		/* GeneratePoly */ {
		int i;
		// look at whether edge (i,i+1) is nixed.
		for (i=0;i+1<it.list.length;i++)
		if ((it.list[i+1]>it.list[i])&&((i==0)||(it.list[i]>it.list[i-1])))
			if (nix[it.list[i]][it.list[i+1]]) { continue genpolyloop; }
		// convert to edge format.
		int count = 1;
		for (i=1;i<it.list.length;i++) if (it.list[i]!=it.list[i-1])
			count++;
		poly = new int[count];
		newCount = new int [count];
		poly[0]= offset;
		for (i=0;i<count;i++) newCount[i]=0;
		int index=0;
		for (i=1;i<it.list.length;i++)
			if (it.list[i]!=it.list[i-1])
			poly[++index] = misc.mod(it.list[i]+offset,outer);
			else newCount[index]++;
		}  /* GeneratePoly */

		if (!looksOK(a,F,poly,newCount,ngon,outer,p)) continue genpolyloop; 
		pushClone(a,F,poly,newCount,tentative,terminal,p);
		}
	}

private static boolean looksOK(arrangePlus a,face F,
		int poly[],int newCount[],int ngon,int outer,parameterUse p )  { 
		for (int i=0;i<poly.length;i++)
			{
			int t,q,e,temp;
			temp= -1;
			if (ngon>4) e = ngon; else e =0;
			if (ngon==4) q = 1; else q=0;
			if (ngon==3) t = 1; else t=0;
			/* forward */ if (i+1<poly.length) {
			int gon= 1+misc.mod(poly[i+1]-poly[i],outer)+newCount[i];
			if (gon==3) t += 1;
			else if (gon>3) temp += 1;
			} /* forward */
			/* back */ if (i>0) { 
			int gon = 1+misc.mod(poly[i]-poly[i-1],outer)+newCount[i-1];
			if (gon==3) t += 1;
			else if (gon>3) temp += 1;
			} /*back*/
			if (!a.LooksOKToAdd(t,q,e,temp,F.elementAt(poly[i]),p)) 
				 return false; 
			}
		return true;
		}

private static void pushClone(arrangePlus a,face F,int poly[],int newCount[],
		Stack tentative,Stack terminal, parameterUse p)
		{ 
		arrangePlus b = (arrangePlus)a.clone();
		int faceIndex = a.indexOf(F);
		face bF = b.faceAt(faceIndex);
		if (bF.size()!= F.size())
			{ misc.warn("clone error"); }
		if (bF.size()<4) { misc.warn("bad face"); }
		b.addPolygon(poly,newCount,bF);
		tentativePush(b,tentative,terminal,p);
		} 

private static void init_nix(boolean nix[][],int offset,arrangePlus a,face F) 
	{
	int outer = nix.length;
	int i,j;
	for (i=0;i<outer;i++) for (j=i+1;j<outer;j++)
		{
		int posi = misc.mod(i+offset,outer);
		int posj = misc.mod(j+offset,outer);
		vertex vi = F.elementAt(posi);
		vertex vj = F.elementAt(posj);
		if (a.divideFaceLooksOK(posi,posj,F)) nix[i][j]=false;
			else nix[i][j]=true;
		}
	for (i=0;i<outer;i++) for (j=0;j<i;j++)
			nix[i][j]=nix[j][i];
	for (i=0;i<outer;i++) nix[i][i]=false;
	}

static int highGon(int i,arrange a)
	{
	if (i<0) return 0;
	vertex v = a.vertexAt(i);
	int high = 0;
	for (int k=0;k<v.size();k++)
		{
		int temp = v.elementAt(k).size();
		if (temp>high) high = temp;
		}
	return high;
	}

	// we set a base point (a.getBasePosition()).  
	// The purpose of this function is to cut down on the number of
	// possibilities considered by factoring out the group action.
	// We make the base point the vertex at exceptional face with
	// the most sides.  Then we pick the vertex on such a face with
	// the most exceptionals, and if there is a tie, with the most quads,
	// and if there is a tie, that with the most triangles.
static boolean baseVertexIsMostCrowded(arrange a)
	{
	if (a.getBasePosition()<0) 
			{
			misc.throwIndex(a.getBasePosition());
			return true;
			}
	vertex bV = a.vertexAt(a.getBasePosition());
	if (bV.tempCount()>0) return true;
	int exC = bV.exCount();
	int quC = bV.quadCount();
	int trC = bV.triCount();
	int hiG = highGon(a.getBasePosition(),a);
	for (int j=0;j<a.vertexSize();j++)
		{
		vertex v = a.vertexAt(j);
		if (j==a.getBasePosition()) continue;
		if (v.tempCount()>0) continue;
		if (highGon(j,a)<hiG) continue;
		if (v.exCount() > exC) return false;
		if (v.exCount() < exC) continue;
		if (v.quadCount() > quC) return false;
		if (v.quadCount() < quC) continue;
		if (v.triCount() > trC) return false;
		}
	return true;
	}

static boolean passesScrutiny(arrangePlus a,parameterUse p)
	{
	if (a.vertexSize()> p.vertexCountMax) return false;
	if (a.vertexSize()< p.vertexCountMin) return false;
	for (int i=0;i<a.vertexSize();i++)
		{
		vertex v = a.vertexAt(i);
		if ((v.exCount()>0)&&
		    (v.size()>p.faceCountMaxAtExceptionalVertex))
			return false;
		if ((v.exCount()==0)&&
			(v.size()>p.faceCountMax)) return false;
		}
	if (a.ExcessNotAt(null,p)+a.faceSquanderLowerBound(p)> p.squanderTarget) 
		return false;
	if (a.ScoreUpperBound(p) < p.scoreTarget) return false;
	return true;
	}

private static boolean has11Type(int i,arrangePlus a)
	{
	// true if vertex i has 1 pentagon, 1 triangle and nothing else.
	vertex v = a.vertexAt(i);
	if (v.size()!=2) return false;
	if (highGon(i,a)!=5) return false;
	if (v.tempCount()>0) return false;
	if (v.quadCount()>0) return false;
	if (v.exCount()!=1) return false;
	if (v.triCount()!=1) return false;
	return true;
	}

static boolean passesTentativeScrutiny(arrangePlus a,parameterUse p)
	{
	if (!baseVertexIsMostCrowded(a)) return false;
	if (p.excludePentQRTet)
		for (int i=0;i<a.vertexSize();i++) if (has11Type(i,a)) 
			{ // System.out.print("rejected"+i+"\n"); 
				return false; }
	return true;
	}

static void terminalPush(arrangePlus a,Stack terminal,parameterUse p)
	{
	if (!passesScrutiny(a,p)) return;
	arrangeInvariant b = new arrangeInvariant(a);
	if (b.isContainedIn(terminal)) return;
	terminal.push(b);
	if (terminal.size()>=500) 
			{
			System.out.print("// file too big -- dumping...");
			arrange.dumpStack(terminal);
			terminal.removeAllElements();
			}
	}

static void tentativePush(arrangePlus a,Stack tentative,Stack terminal,
	parameterUse p)
	{
	face f;
	boolean hasTemp=false;
	for (int i=0;i<a.faceSize(); i++)
		{
		f = a.faceAt(i);
		if ((a.vertexSize()>3)&&(f.size() == 3)) 
			f.setTempState(false);
		if (f.getTempState()) { hasTemp = true; }
		}
	if (!passesTentativeScrutiny(a,p)) return;
	if (!hasTemp)  { terminalPush(a,terminal,p); return ; }
	tentative.push(a);
	}

static int PolyLimit(arrangePlus a,parameterUse p)
	{
	int polylimit = p.maxGon();
	int lb = a.faceSquanderLowerBound(p) + a.ExcessNotAt(null,p);
	while ((lb+p.squanderFace(polylimit)>=p.squanderTarget)&&
			(polylimit>0))
		polylimit --;
	return polylimit;
	}

static boolean niceQuad(arrangePlus a,parameterUse p)
	{
	int s = a.faceSquanderLowerBound(p)+ a.ExcessNotAt(null,p);
	boolean b = ((s +p.squanderFaceStartingAt(5)>=p.squanderTarget)&&
		(s+p.squanderVertex[2][1]>=p.squanderTarget)&&
		(s+p.squanderVertex[0][2]>=p.squanderTarget));
	if (b) discard.addLine("nQ passes"); else discard.addLine("nQ fails");
	return b;
	}

static void QuadCluster(arrangePlus a,face F,parameterUse p,Stack 
	tentative,Stack terminal,vertex v[])
	{
	boolean doquad = true;
	for (int i=0;i<4;i++)
		if (!a.LooksOKToAdd(0,1,0,-1,v[i],p)) 
			{ doquad=false; break; }
	debug.d5(a,F,doquad);
	if (doquad)
		{
		arrangePlus b = (arrangePlus)a.clone();
		int index = a.indexOf(F);
		b.faceAt(index).setTempState(false);
		tentativePush(b,tentative,terminal,p);
		}
	}

static void Do40(arrangePlus a,face F,parameterUse p,Stack 
	tentative,Stack terminal,vertex v[]) {
	boolean do40 = true;
	if (a.faceSquanderLowerBound(p) +a.ExcessNotAt(null,p)+
		p.forecast(4,0,0)>=
		p.squanderTarget) do40=false;
	else for (int i=0;i<4;i++)
		if (!a.LooksOKToAdd(2,0,0,-1,v[i],p)) { do40=false; break; }
	if (do40)
		{
		arrangePlus b = (arrangePlus)a.clone();
		int index = a.indexOf(F);
		face F40 = b.faceAt(index);
		vertex A1 = F40.elementAt(0);
		vertex A2 = F40.elementAt(1);
		vertex A3 = F40.elementAt(2);
		vertex A4 = F40.elementAt(3);
		b.divideFace(0,2,F40,1); 
		index  = b.vertexSize();
		if (b.vertexSize()!= a.vertexSize()+1)
			misc.warn("vertex was not added as expected");
		index = b.faceSize();
		face F40p = b.faceAt(index-1);
		if (F40p.size()!=4) misc.warn("quadrilateral expected");
	 	index = F40p.indexOf(A2);
		if (index<0) misc.warn("A2 not found on F40p");
		b.divideFace(index,misc.mod(index+2,4),F40p,0); 
		index = F40.indexOf(A4);
		if (F40.size()!=4) misc.warn("F40 quad expected");
		if (index<0) misc.warn("A4 not found on F40");
		b.divideFace(index,misc.mod(index+2,4),F40,0); 
		tentativePush(b,tentative,terminal,p);
		}
	}

static void Quad(arrangePlus a,face F,parameterUse p,Stack tentative,
		Stack terminal)
	{
	vertex v[] = new vertex[4];
	for (int i=0;i<4;i++) v[i]=F.elementAt(i);
	QuadCluster(a,F,p,tentative,terminal,v);
	Do40(a,F,p,tentative,terminal,v);

	/* tri2: */ {
	boolean dotri2=true;
	if (!a.divideFaceLooksOK(1,3,F)) dotri2=false;
	else for (int i=0;i<4;i++)
		if (!a.LooksOKToAdd(1+misc.mod(i,2),0,0,-1,v[i],p))
			{ dotri2=false; break; }
	if (dotri2)
		{
		arrangePlus b = (arrangePlus)a.clone();
		int index = a.indexOf(F);
		face F0 =  b.faceAt(index);
		b.divideFace(1,3,F0,0); 
		tentativePush(b,tentative,terminal,p);
		}
	}

	/* skewtri2: */ {
	boolean doskewtri2=true;
	if (!a.divideFaceLooksOK(0,2,F)) doskewtri2=false;
	for (int i=0;i<4;i++)
		if (!a.LooksOKToAdd(1+misc.mod(i+1,2),0,0,-1,v[i],p))
			{ doskewtri2=false; break; }
	if (doskewtri2)
		{
		arrangePlus b = (arrangePlus)a.clone();
		int index = a.indexOf(F);
		face F0 =  b.faceAt(index);
		b.divideFace(0,2,F0,0); 
		tentativePush(b,tentative,terminal,p);
		}
	}
	} // end quad.
	
private static int globalcount = 0;
private static void resetCount() { globalcount = 0; }

static boolean loop(Stack tentative,Stack terminal,parameterUse p)
	{
	arrangePlus a;
	if (tentative.size()>0) 
		{ a = (arrangePlus)tentative.pop();
		  globalcount++;
		  discard = new arrangeAnnotate(a);
		  if (0==misc.mod(globalcount,1000)) // adapt as needed 
			{
			System.out.println("// stack sizes = "
				+tentative.size()+" "+terminal.size());
			System.out.println("// cases considered= "+globalcount);
			}
		}
	else {  System.out.println("// Total pops: "+globalcount); return false; }

	/* ForcedTriangle: */ {
	forcedTriangle: for (int i=0;i<a.vertexSize();i++)
		{
		vertex v = a.vertexAt(i);
		if (v.tempCount()==0) continue forcedTriangle;
		if (!a.ForcedTriangleAt(v,p)) continue forcedTriangle;
		// if we get to here there is a forced triangle:
		face F = v.elementAt(0);
		int j=0;
		while (!F.getTempState()) F = v.elementAt(++j);
		int vpos = F.indexOf(v);
		int sz = F.size();
		vertex A = F.elementAt(misc.mod(vpos+1,sz));
		vertex B = F.elementAt(misc.mod(vpos-1,sz));
		discard.addLine("Forced triangle at vertex "+(i+1)+" face "+
				(1+a.indexOf(F))+" between vertices "+
				(1+a.indexOf(A))+" "+(1+a.indexOf(B)));
		if (a.divideFaceLooksOK(vpos+1,vpos-1,F)&&
			(a.LooksOKToAdd(1,0,0,-1,A,p))&&
			(a.LooksOKToAdd(1,0,0,-1,B,p)))
			{
			a.divideFace(vpos+1,vpos-1,F,0);
			tentativePush(a,tentative,terminal,p);
			return true;
			}
		discard.addLine("But forced triangle not added"+
				" because it didn't look ok");
		return true; /*true because it has been dealt with*/
		}
	}/*ForcedTriangle*/

	/* pick a face to use */
	face Fsmall = a.faceAt(0);
	int outer = p.vertexCountMax+1;
	for (int i=0;i<a.faceSize();i++)
		{
		face F = a.faceAt(i);
		if ((F.getTempState())&&(F.size()<outer))
			{ Fsmall = F; outer = F.size(); }
		}
	if (badData(outer,p,Fsmall)) 
		{ System.out.println("BAD DATA!"); return false; }

	/*polygons*/ {
	int startingIndex = a.selectVertex(Fsmall);
	if ((outer==4)&&(niceQuad(a,p)))  
		{ Quad(a,Fsmall,p,tentative,terminal); return true; }
	int polylimit = PolyLimit(a,p);
	discard.addLine("polylimit = "+polylimit);
	if ((outer==4)&&(a.vertexSize()>5)) 
		polylimit = Math.min(polylimit,5);
	vertex A = Fsmall.elementAt(startingIndex);
	for (int i=3;i<=polylimit;i++)
		{
		generatePolygon(i,a,A,Fsmall,tentative,terminal,p);
		}
	return true;
	}

	}

public static Stack generateAll(int NGON)
	{
	Stack tentative = new Stack();
	Stack terminal = new Stack();
	parameterUse p = new parameterUse(NGON);
	resetCount();
	tentative.push(new arrangePlus(arrange.makeNgon(NGON)));
	while ( loop(tentative,terminal,p) ) {  } 
	System.out.println("stack size = " + terminal.size());	
	System.out.println("count = " + globalcount);
	return terminal;
	}

public static Stack generateQuad(int casenum)
	{
	Stack tentative = new Stack();
	Stack terminal = new Stack();
	parameterUse p = new parameterUse(casenum,true);
	resetCount();
	tentative.push(new arrangePlus(arrange.seed(p.quadCases[casenum])));
	while ( loop(tentative,terminal,p) ) {  } 
	System.out.println("// stack size = " + terminal.size());	
	System.out.println("// count = " + globalcount);
	return terminal;
	}

private static boolean badData(int outer,parameterUse p,face Fsmall)
	{
	if ((outer> p.vertexCountMax) ||(!Fsmall.getTempState()))
			{ 
			String S = "Unable to locate face; outer = "+outer;
			discard.addLine(S);
			(new graphException(S)).report(); 
			return true; 
			} 
	return false;
	}

}  /* end of arrangeMain */


class iterator {
	static int inc(int i) { return i+1; }
	int maxelement;
	int list[];
	void print() { for (int i=0;i< list.length;i++) System.out.print(list[i]);
				System.out.println(" "); }
	iterator(int inner,int outer)
			{
			list = new int[inner];
			maxelement = outer-1;
			for (int i=0;i<list.length;i++) list[i]=0;
			list [list.length-1] = maxelement;
			list [list.length-2] = -1;
			}
	boolean next()
			{
			// always: list[list.length-1]==maxelement
			// always: list[0]=0;
			int r = list.length-2;
			int newvalue;
			while ((r>0)&&(list[r]==maxelement-1))
					{  r--; }
			if (r<=0) return false;
			newvalue = list[r]+1;
			for (int i=r;i<list.length-1;i++) { list[i]=newvalue; }
			return true;
			}
	}


class debug {

static void d1(int ngon,arrange a,arrange.vertex A,arrange.face F)
	{
	arrange.discard.addLine("generating "+ngon+
		" on face "+iterator.inc(a.indexOf(F)));
	arrange.discard.addLine("starting from vertex "
		+iterator.inc(a.indexOf(A)));
	}

static void d2(int offset) {
	if (offset<0) 
		{ (new graphException("generate Polygon A in F")).report();}}

static void d3(int outer) { if (outer<4) 
		misc.warn("illicit face"); }

static void d5(arrange a,arrange.face F,boolean doquad)
	{ int debug=a.indexOf(F); 
		if (doquad) arrange.discard.addLine("dQ on face "+(debug+1)); 
		else arrange.discard.addLine("dQ fails on face "+(debug+1)); }
	}
 
