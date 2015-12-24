import java.util.*;
import java.lang.*;

/* draw a clockwise circle around each vertex.  List the faces
	cyclicly going around in vertexCycle;
   draw a counterclockwise loop around the boundary of each face.
	List the vertices on the face going around.   */


class arrange {
	static { System.out.print(" -- arrange $Revision: 1.2 $ -- \n"); }
	static class vertex {
        int triCount()
                {
                int count = 0;
				face f;
                for (int i=0;i<size();i++)
				{
				f = elementAt(i);
				if ((!f.getTempState())&&
									(f.size() == 3)) count ++;
				}
                return count;
                }
        int quadCount()
                {
                int count = 0;
				face f;
                for (int i=0;i<size();i++)
				{
				f = elementAt(i);
				if ((!f.getTempState())&&
									(f.size() == 4)) count ++;
				}
                return count;
                }
        int exCount()
			{
			int count = 0;
			face f;
			for (int i=0;i<size();i++)
			{
			f = elementAt(i);
			if ((!f.getTempState())&&(f.size() > 4)) count ++;
			}
			return count;
			}
        int tempCount()
			{
			int count = 0;
			for (int i=0;i<size();i++)
				if (elementAt(i).getTempState()) count ++;
			return count;
			}
        private int height;
		int getHeight() { return height; }
		void setHeight(int h) { height = h; }
        Vector incidentFace;
		face elementAt(int i) { return (face)incidentFace.elementAt(i); }
		int size() { return incidentFace.size(); }
		int indexOf(face f) { return incidentFace.indexOf(f); }
        }

	static class face {
		Vector incidentVertex;
		int size() { return incidentVertex.size(); }
		int indexOf(vertex v) { return incidentVertex.indexOf(v); }
		vertex elementAt(int i) { return (vertex)incidentVertex.elementAt(i); }
		boolean getTempState() { return isTemp; }
		void setTempState(boolean state) { isTemp = state; }
		private boolean isTemp;
        }

	static arrange makeNgon(int ngon)
		{
		arrange a = new arrange();
		a.setBasePosition(0);
		a.vertexList = new Vector(ngon);
		a.faceList = new Vector(2);
		face F[] = {new face(),new face()};
		a.faceList.addElement(F[0]);
		a.faceList.addElement(F[1]);
		F[0].setTempState(false);
		F[1].setTempState(true);
		F[0].incidentVertex = new Vector();
		F[1].incidentVertex = new Vector();
		for (int i=0;i<ngon;i++) 
			{
			vertex v = new vertex();
			a.vertexList.addElement(v);
			v.incidentFace = new Vector();
			v.incidentFace.addElement(F[0]);
			v.incidentFace.addElement(F[1]);
			v.setHeight(0);
			F[0].incidentVertex.addElement(v);
			F[1].incidentVertex.insertElementAt(v,0);
			}
		return a;
		}

	vertex vertexAt(int i) { return (vertex)vertexList.elementAt(i); }
	int vertexSize() { return vertexList.size(); }
	int indexOf(vertex v) { return vertexList.indexOf(v); }
	face faceAt(int i) { return (face)faceList.elementAt(i); }
	int faceSize() { return faceList.size(); }
	int indexOf(face f) { return faceList.indexOf(f); }

	protected Vector vertexList;
	protected Vector faceList;

	private int basePosition =-1;
	void setBasePosition(int pos) { basePosition=pos; }
	int getBasePosition() { return basePosition; }

// break a face in two by drawing an edge between A and B;
// create newCount new vertices between A and B.
void divideFace(int posA,int posB,face F,int newCount) 
	{
	int ngon = F.size();
	if (misc.mod(posA,ngon)==misc.mod(posB,ngon)) 
		{ misc.warn("Exception coming: "+posA+" "+posB+" "+ngon); }
	posA = misc.mod(posA,ngon);
	posB = misc.mod(posB,ngon);
	vertex A = F.elementAt(posA);
	vertex B = F.elementAt(posB);
	if (!F.getTempState()) 
		{ (new graphException("F is not a temp")).report(); }
	if (posA==posB) 
		{ (new graphException("edge attempted between"+
				" vertex and itself "+posA+" "+posB+" "+ngon)).report(); }
	int i;
	face F2 = new face();
	F2.setTempState(true);
	Vector bdry1=  cyclicRange(posA,posB,F.incidentVertex);
	Vector bdry2= cyclicRange(posB,posA,F.incidentVertex);
	if (bdry1.size()+newCount<3) return;
	if (bdry2.size()+newCount<3) return;

	F.incidentVertex = bdry2;
	F2.incidentVertex = bdry1;

	// insert newCount more vertices between A and B:
	// update vertexList and faceList. and vertices vi.
	faceList.addElement(F2);
	vertex v[] = new vertex[newCount];
	for (i=0;i<newCount;i++) 
		{
		v[i]= new vertex();
		vertexList.addElement(v[i]);
		v[i].setHeight(
		  Math.min(A.getHeight()+1+i,B.getHeight()+newCount-i));
		v[i].incidentFace = new Vector(2);
		v[i].incidentFace.addElement(F);
		v[i].incidentFace.addElement(F2);
		}
	// update faces at A;
	int index = A.indexOf( F);
	if (index<0) 
		{ (new graphException("face not on vertex A")).report(); }
	A.incidentFace.insertElementAt(F2,index+1);
	// update faces at B;
	index = B.indexOf( F);
	if (index<0) 
		{ (new graphException("face not on vertex B")).report(); }
	B.incidentFace.insertElementAt(F2,index);
	// done with vertices

	// now update each face:
	index = F.indexOf( A);
	if (index<0) 
		{ (new graphException("vertex A not on face")).report(); }
	for (i=0;i<newCount;i++)
		F.incidentVertex.insertElementAt(v[i],index+i+1);
	// tell vertices around F2 that they are no longer on F.
	for (i=1;i<bdry1.size()-1;i++)
		{
		vertex w = (vertex) bdry1.elementAt(i);
		index = w.indexOf(F);
		if (index<0) 
		 	{ (new graphException("vertex w not on face")).report(); }
		w.incidentFace.setElementAt(F2,index);
		}
	index = F2.indexOf( B);
	if (index<0) 
		{ (new graphException("face(2) not on vertex B")).report(); }
	for (i=0;i<newCount;i++)
		F2.incidentVertex.insertElementAt(v[i],index+1);

	/*debug*/ if (false) {
	int j;
	System.out.print("F = "+indexOf(F)); 
	for (i=0;i<F.size();i++)
		{
		vertex iv = F.elementAt(i);
		System.out.print(" "+indexOf(iv));
		for (j=0;j<iv.size();j++)
			{
			face iff = iv.elementAt(j);
			System.out.print("*"+indexOf(iff));
			}
		}
	System.out.println("\nF2 = "+indexOf(F2));
	for (i=0;i<F.size();i++)
		{
		vertex iv =  F2.elementAt(i);
		System.out.print(" "+indexOf(iv));
		for (j=0;j<iv.size();j++)
			{
			face iff =  iv.elementAt(j);
			System.out.print("*"+indexOf(iff));
			}
		}
	System.out.println();
	} /*debug*/
	
	}

protected static vertex next(vertex A,int Findex) 
	{
	face F = A.elementAt(Findex);
	int index = F.indexOf(A);
	if (index<0)  { (new graphException
			("The vertex does not lie on the face")).report(); }
	return F.elementAt( misc.mod(index+1,F.size()));
	}

private boolean areAdjacent(vertex A,vertex B)
        {
        for (int i=0;i<A.size();i++)
		if (B==next(A,i)) return true; 
		return false;
        }

// newCount =0 in this, just straight from A to B.
boolean divideFaceLooksOK(int posA,int posB,face F) 
	{
	int nGon = F.size();
	if (nGon <5) return true;
	posA = misc.mod(posA,nGon);
	posB = misc.mod(posB,nGon);
	vertex A = F.elementAt(posA);
	vertex B = F.elementAt(posB);
	int i,j;
	int AtoB = misc.mod(posB-posA,nGon);
	int BtoA = misc.mod(posA-posB,nGon);
	int length = Math.min(AtoB,BtoA);

	if (length<2) return true;
	if (areAdjacent(A,B)) return false; // definitely a double join.

	if (length<3) return true;
	vertex vA[] = new vertex [A.size()];
	vertex vB[] = new vertex [B.size()];
	for (i=0;i<vA.length;i++)
		vA[i]= next(A,i);
	for (i=0;i<vB.length;i++)
		vB[i]= next(B,i);
	for (i=0;i<vA.length;i++) for (j=0;j<vB.length;j++)
	  if (vA[i]==vB[j]) return false; // definite bad triangle.
	
	// we'll check for bad quadrilaterals at the end.
	return true;
	}

void addPolygon(int poly[],int newCounts[],face F)
	{
	if (poly.length!=newCounts.length) 
		misc.warn("unequal lengths in addPolygon");
	vertex vpoly[] = new vertex[poly.length];
	for (int i=0;i<poly.length;i++) vpoly[i]=F.elementAt(poly[i]);
 
	/*debug*/{
	String S="adding on face "+(1+indexOf(F))+" at vertices";
	for (int i=0;i<poly.length;i++)
		{
		S = S + " "+ (indexOf(vpoly[i])+1);
		if (newCounts[i]>0) S = S+":"+newCounts[i];
		}
	discard.addLine(S);
	}/*debug*/
	
	
	for (int i=0;i<poly.length;i++)
		{
		int j = misc.mod(i+1,poly.length);
		int Indexi = F.indexOf(vpoly[i]);
		int Indexj = F.indexOf(vpoly[j]);
		if (Indexi<0) { misc.throwIndex(Indexi); }
		if (Indexj<0) { misc.throwIndex(Indexj); }
		divideFace(Indexi,Indexj,F,newCounts[i]); 
		}
	F.setTempState(false);
	}

// for initialization.
private void faceInit(int extra)
	{
	vertex v = vertexAt(0);
	face f=null;
	for (int i=0;i<v.size();i++)
		if (v.elementAt(i).getTempState()) f = v.elementAt(i);
	int index = f.indexOf(v);
	int poly[] = {misc.mod(index+1,f.size()),index};
	int count[] = {1+extra,0};
	addPolygon(poly,count,f);
	}

private void faceInitEnd(int extra)
	{
	vertex v = vertexAt(0);
	face f = null;
	for (int i=0;i<v.size();i++)
		if (v.elementAt(i).getTempState()) f = v.elementAt(i);
	int index = f.indexOf(v);
	int poly[] = {misc.mod(index+1,f.size()),misc.mod(index-1,f.size()),index};
	int count[] = {extra,0,0};
	addPolygon(poly,count,f);
	}

public static arrange seed(int array[])
	{
	arrange a = makeNgon(array[0]+3);
	for (int i=1;i<array.length-1;i++) a.faceInit(array[i]);
	a.faceInitEnd(array[array.length-1]);
	return a;
	}
	

private Vector cyclicRange(int indexA,int indexB,Vector c) 
	{
	indexA = misc.mod(indexA,c.size());
	indexB = misc.mod(indexB,c.size());
	int top = 1+misc.mod(indexB-indexA,c.size());
	Vector d = new Vector(top);
	for (int i=0;i<top;i++)
		d.addElement(c.elementAt(misc.mod(i+indexA,c.size())));
	return d;
	}

 public arrange(archive.GraphArrays gr)
	{
	int i,j;
	faceList = new Vector(gr.faceList.length);
	vertexList = new Vector(gr.adjacent.length);
	vertex v[] = new vertex[gr.vertexList.length];
	for (i=0;i<gr.vertexList.length;i++)
		{
		v[i] = new vertex();
		v[i].setHeight(0);
		v[i].incidentFace = new Vector(gr.vertexList[i].length);
		vertexList.addElement(v[i]);
		}
	face f[] = new face[gr.faceList.length];
	for (i=0;i<gr.faceList.length;i++)
		{
		f[i] = new face();
		f[i].setTempState(false);
		f[i].incidentVertex = new Vector(gr.faceList[i].length
);
		for (j=0;j<gr.faceList[i].length;j++)
			f[i].incidentVertex.addElement
			  (v[gr.faceList[i][j]]);
		faceList.addElement(f[i]);
		}
 
	for (i=0;i<gr.vertexList.length;i++)
	for (j=0;j<gr.vertexList[i].length;j++)
		{
		v[i].incidentFace.addElement
			(f[gr.vertexList[i][j]]);
		}
	for (i=0;i<gr.tempList.length;i++)
		{
		f[gr.tempList[i]].setTempState(true);
		}
	}

public archive.GraphArrays makeGr() // convert back
	{
	archive.GraphArrays gr = new archive.GraphArrays();
	int i,j,r;
	/* faceList: */ {
	gr.faceList= new int[faceSize()][];
	for (i=0;i<faceSize();i++)
		{
		face f = faceAt(i);
		gr.faceList[i]= new int[f.size()];
		for (j=0;j<f.size();j++)
			{
			gr.faceList[i][j]= indexOf(f.elementAt(j));
			}
		}
	}

	/* tempList: */ {
	r = 0;
	for (i=0;i<faceSize();i++) 
	   if (faceAt(i).getTempState()) { r++ ; }
	gr.tempList = new int[r];
	r = 0;
	for (i=0;i<faceSize();i++) 
	   if (faceAt(i).getTempState())
		{  gr.tempList[r++]= i; }
	}

	/* vertexList: */ {
	gr.vertexList = new int[vertexList.size()][];
	for (i=0;i<vertexSize();i++)
		{
		vertex v = vertexAt(i);
		//old: Vector f = vertexAt(i).incidentFace;
		gr.vertexList[i] = new int[v.size()];
		for (j=0;j<v.size();j++)
			{
			gr.vertexList[i][j]=indexOf(v.elementAt(j));
			}
		}
	}

	/* adjacent: */ {
	int k;
	gr.adjacent= new int[vertexSize()][];
	for (i=0;i<vertexSize();i++)
		{
		gr.adjacent[i]= new int[gr.vertexList[i].length];
		for (j=0;j<gr.adjacent[i].length;j++)
			{
			r = gr.vertexList[i][j];
			k=0;
			while (i!=gr.faceList[r][k]) { k++; } // hope to find
			k = misc.mod(k+1,gr.faceList[r].length);
			gr.adjacent[i][j]= gr.faceList[r][k];
			}
		}
	}
	return gr;

	}
	
arrange() { }

// to think that other languages do this with the equals sign!
protected synchronized Object clone() 
	{
	arrange a = new arrange();
	a.vertexList = new Vector(vertexSize());
	a.faceList = new Vector(faceSize());
	vertex[] avL = new vertex[vertexSize()];
	face[] afL = new face[faceSize()];
	int i,j,r;
	vertex u, v;
	face f, g;
	Object b,c;
	for (i=0;i<vertexSize();i++)
		{
		avL[i]= new vertex();
		a.vertexList.addElement(avL[i]);
		}
	for (i=0;i<faceSize();i++)
		{
		afL[i]= new face();
		a.faceList.addElement(afL[i]);
		}
	for (i=0;i<faceSize();i++)
		{
		f = faceAt(i);
		g = afL[i];
		g.setTempState(f.getTempState());
		g.incidentVertex = new Vector(f.size());
		for (j=0;j<f.size();j++)
			{ b = f.elementAt(j);
			  r = vertexList.indexOf(b);
			  if (r<0) { misc.warn("uncloned vertex"); }
			  c = avL[r];
			  g.incidentVertex.addElement(c);
		        }
			  
		}
	for (i=0;i<vertexSize();i++)
		{
		u = vertexAt(i);
		v = avL[i];
		v.incidentFace = new Vector(u.size());
		v.setHeight(u.getHeight());
		for (j=0;j<u.size();j++)
			{ b = u.elementAt(j);
			  r = faceList.indexOf(b);
			  if (r<0) { misc.warn("uncloned face"); }
			  c = afL[r];
			  v.incidentFace.addElement(c);
			}
		}
	a.basePosition = basePosition;
	return a;
	}

private static String faceToString(arrange a,face F)
    {
    String S = " ";
    S = S+F.size();
    for (int i=0;i<F.size();i++)
        S = S+" "+ a.indexOf(F.elementAt(i));
    return S;
    }

static String ToString(arrange a)
    {
    int numTemp = 0;
    String pre=" "+a.faceSize();
    String S = " ";
    for (int i=0;i<a.faceSize();i++) 
        {
        face F = a.faceAt(i);
        if (F.getTempState()) { numTemp++; pre = " "+i+pre; }
        S = S + faceToString(a,F);
        }
    return " "+numTemp+pre+S;
    }

static void dumpStack(Vector K)
	{
	if (K.size()==0) { System.out.println("-- no data--"); return; }
	System.out.println("public class graphUNKNOWN {");
	System.out.println("final static String data[] = {\n");
	for (int i=0;i<K.size();i++)
		{
		System.out.print("\""+ToString((arrange)K.elementAt(i)));
		if (i+1==K.size())
		System.out.println("\"};\n\n};");
		else System.out.println("\",\n" ) ;
		}
	}

private static String ptString(arrange a,vertex v)
	{
	String head = "{{"+v.triCount()+","+v.quadCount()+","+v.exCount()+"},{";
	for (int i=0;i<v.size();i++)
		{
		head = head + (1+a.indexOf(next(v,i)));
		if (i+1<v.size()) head = head+",";
		}
	head = head+"}}";
	return head;
	}
	

static String MathematicaString(arrange a)
	{
	// format {0,0,{p1,...,pr}}; pi = {{p,q,r},{adjacent}};
	String head = "{0,0,{";
	for (int i=0;i<a.vertexSize();i++)
			{
			head = head + ptString(a,a.vertexAt(i));
			if (i+1<a.vertexSize()) head = head+",\n";
			}
	head = head +"}}";
	return head;
	}


// for debugging:
public static arrangeAnnotate discard;



}

