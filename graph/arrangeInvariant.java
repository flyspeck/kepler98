import java.util.*;
import java.lang.*;

class arrangeInvariant extends arrange {
	static { System.out.print(" -- arrangeInvariant $Revision: 1.2 $ -- \n"); }
    private Hashtable vertexInvariant=null; 
    private Hashtable faceInvariant=null; 
    private long invariant=0;
    private static final long prime =15485863; //=Prime[10^6] //was=22801763489;

long getInvariant()
    {
	if (invariant==0) { misc.warn("nonzero invariant expected"); 
			setInvariant(); }
    return invariant;	
    }

arrangeInvariant(archive.GraphArrays gr) { super(gr); }

arrangeInvariant(arrange a)
	{
	arrange b = (arrange)a.clone();
	vertexList = b.vertexList;
	faceList = b.faceList;
	setBasePosition(-1);
	setInvariant();
	}

arrangeInvariant() {}

protected synchronized Object clone() 
	{
	arrangeInvariant b = new arrangeInvariant();
	arrange c = (arrange)super.clone();
	b.vertexList = c.vertexList;
	b.faceList = c.faceList;
	return b;
	}

private static Long symdih(long x[])
	{
	long total =0;
	for (int i=0;i<x.length;i++)
		total = (total + x[i]*x[(i+1) % x.length]) % prime;
	return new Long(total);
	}

private static long symcyc(long x[])
	{
	long total = 0;
	for (int i=0;i<x.length;i++)
		{
		long t = x[i]*x[i] % prime;
		t = t*x[(i+1) % x.length] % prime;
		total = (total + t) % prime;
		}
	return total;
	}

/* // Delete this.
private static Long symnone(long x[])
	{
	long total = 0;
	for (int i=0;i<x.length;i++)
		{ // 7919 a prime.
		total = (7919*total + x[i]) % prime;
		}
	return new Long (total);
	}
*/

private static Long vertexInvariant0(int tCount,int qCount,int exCount)
	{
	long r = tCount + 8*qCount + 64*exCount;
	return new Long(r*r);
	}

void setInvariant()
	{
	vertexInvariant = new Hashtable(vertexSize());
	faceInvariant = new Hashtable(faceSize());
	int i,j;
	invariant =0;
	for (i=0;i<vertexSize();i++)
		{
		vertex v = vertexAt(i);
		vertexInvariant.put(v,vertexInvariant0(
		  v.triCount(),v.quadCount(),v.exCount()));
		}
	for (i=0;i<faceSize();i++)
		{
		face f = faceAt(i);
		faceInvariant.put(f,new Long((long) 0));
		}
	for (int level=0;level<1;level++)
	{
	for (i=0;i<faceSize();i++)
		{
		face f = faceAt(i); 
		Vector vec = f.incidentVertex;
		long x[] = new long[f.size()];
		for (j=0;j<f.size();j++)
			x[j]=((Long) vertexInvariant.get
			   (vec.elementAt(j))).longValue();
		faceInvariant.remove(f);
		Long sd = symdih(x);
		faceInvariant.put(f,sd);
		invariant = (invariant + 
			sd.longValue()*sd.longValue()) % prime;
		}
	for (i=0;i<vertexSize();i++)
		{
		vertex v =  vertexAt(i);
		Vector vec = v.incidentFace;
		long x[] = new long[vec.size()];
		for (j=0;j<vec.size();j++)
			x[j]=((Long) faceInvariant.get
			   (vec.elementAt(j))).longValue();
		faceInvariant.remove(v);
		Long sd = symdih(x);
		vertexInvariant.put(v,sd);
		invariant = (invariant + 
			sd.longValue()*sd.longValue()) % prime;
		}
	}
	}

static void removeDuplication(Vector V /* of arrangeInvariants */)
	{
	int i,j;
	arrangeInvariant h[] = new arrangeInvariant[V.size()];
	for (i=0;i<h.length;i++) 
		{
		h[i]= (arrangeInvariant) V.elementAt(i);
		h[i].setInvariant();
		}
	arrangeInvariant temp;
	for (i=0;i<h.length;i++)
	for (j=0;j<h.length-1;j++)
		{
		if (h[j].invariant<h[j+1].invariant)
			{
			temp=h[j]; h[j]=h[j+1]; h[j+1]=temp;
			}
		}
	int newtop=0;
	for (i=0;i<h.length-1;i++)
		if (!Isomorphic(h[i],h[i+1])) h[newtop++]= h[i]; 
	V.removeAllElements();
	for (i=0;i<newtop;i++) V.addElement(h[i]);
	}

boolean isContainedIn(Vector V)
	{
	for (int i=0;i<V.size();i++)
			{
			arrangeInvariant a =(arrangeInvariant)V.elementAt(i);
			if ((a.getInvariant()==getInvariant())&&(Isomorphic(this,a)))
					return true;
			}
	return false;
	}

void addPolygon(int a[],int b[],face F)
	{
	super.addPolygon(a,b,F);
	invariant=0;
	}

void divideFace(int a,int b, face F,int c)
	{
	super.divideFace(a,b,F,c);
	invariant=0;
	}

private static Vector reverse(Vector V)
	{
	Vector W = new Vector (V.size());
	for (int i=V.size();i>0;i--) W.addElement(V.elementAt(i-1));
	return W;
	}

void flip()
	{
	int i;
	vertex v;
	face f;
	for (i=0;i<vertexSize();i++)
	  {
	  v = vertexAt(i);
	  v.incidentFace = reverse((Vector) v.incidentFace);
	  }
	for (i=0;i<faceSize();i++)
	  {
	  f = faceAt(i);
	  f.incidentVertex = reverse((Vector) f.incidentVertex);
	  }
	}

static boolean Isomorphic(arrangeInvariant a,arrangeInvariant b)
	{
	if (ProperlyIsomorphic(a,b)) return true;
	a.flip();
	if (ProperlyIsomorphic(a,b)) return true;
	return false;
	}

private long cyclicInvariant(vertex v)
	{
	long list[] = new long[v.size()];
	for (int i=0;i<list.length;i++)
		list[i]= ((Long)vertexInvariant.get(next(v,i))).longValue();
	return symcyc(list);
	}

private void SelectVertex0()
	{
	int index=-1;
	long minhash = prime;
	long symhash = prime;
	for (int i=0;i<vertexSize();i++)
		{
		vertex v = vertexAt(i);
		long r = ((Long) vertexInvariant.get(v))
				.longValue();
		if ((r<minhash)||
		    ((r==minhash)&&(cyclicInvariant(v)<symhash)))
			 { minhash = r; symhash = cyclicInvariant(v); index = i; }
		}
	if (index==0) return;
	
	vertex temp = vertexAt(index);
	vertexList.setElementAt(vertexAt(0),index);
	vertexList.setElementAt(temp,0);
	}

private void SelectVertex1()
	{
	vertex v0 = vertexAt(0);
	int index = -1;
	long minhash = prime;
	long symhash = prime;
	for (int i=0;i<v0.size();i++)
		{
		vertex v = next(v0,i);
		long r = ((Long) vertexInvariant.get(v))
				.longValue();
		if ((r<minhash)||
			((r==minhash)&&(cyclicInvariant(v)<symhash)))
			 {  minhash = r; symhash = cyclicInvariant(v); index = i; }
		}
	vertex temp = next(v0,index);
	index = indexOf(temp);
	vertexList.setElementAt(vertexAt(1),index);
	vertexList.setElementAt(temp,1);
	}

void canonify()
	{
	int lastset=1;
	int filllevel = 0;
	int i,j=0,k,index;
	vertex v,w,neighbors[];
	setInvariant();
	SelectVertex0(); SelectVertex1();
	while(lastset+1<vertexSize())
	{
	v = vertexAt(filllevel);
	neighbors = new vertex [v.size()];
	for (i=0;i<v.size();i++)
		neighbors[i]= next(v,i);
	int smallindex = vertexSize();
	for (i=0;i<v.size();i++)
		{
		index = indexOf(neighbors[i]);
		if (index<smallindex) {smallindex=index; j=i; }
		}
	for (i=j;i<j+v.size();i++)
		{
		index = indexOf(neighbors[misc.mod(i,neighbors.length)]);
		if (index>lastset) {
		lastset++;
		w = vertexAt(lastset);
		vertexList.setElementAt(w,index);
		vertexList.setElementAt(
			neighbors[misc.mod(i,neighbors.length)],lastset);
		}
		}
	filllevel++;
	} 
}


private static boolean ProperlyIsomorphic(arrangeInvariant a,arrangeInvariant b)
	{
	int i,j;
	if (a.vertexSize()!= b.vertexSize()) { return false; }
	if (a.faceSize()!= b.faceSize()) { return false; }
	a.canonify();
	b.canonify();
	for (i=0;i<a.vertexSize();i++)
		{
		vertex va = a.vertexAt(i);
		vertex vb = b.vertexAt(i);
		int ta[] = new int [va.size()];
		int tb[] = new int [vb.size()];
		if (ta.length!=tb.length) return false;
		for (j=0;j<va.size();j++)
			{
			ta[j]= a.indexOf(next(va,j));
			tb[j]= b.indexOf(next(vb,j));
			}
		int offset = -1;
		for (j=0;j<tb.length;j++) if (tb[j]==ta[0]) offset=j;
		if (offset<0) return false;
		for (j=0;j<ta.length;j++) 
			if (tb[misc.mod(j+offset,tb.length)]!=ta[j])
				return false;
		}
	// bijection of graphs establised!
	return true;
	}


}
