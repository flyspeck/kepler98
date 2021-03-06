import java.util.*;

class arrangePlus extends arrange {

	static { System.out.print(" -- arrangePlus $Revision: 1.2 $ -- \n"); }

arrangePlus(arrange a)
    {
    arrange b = (arrange)a.clone();
    vertexList = b.vertexList;
    faceList = b.faceList;
	setBasePosition(a.getBasePosition());
    }

arrangePlus() {}

protected synchronized Object clone()
    {
    arrangePlus b = new arrangePlus();
    arrange c = (arrange)super.clone();
    b.vertexList = c.vertexList;
    b.faceList = c.faceList;
	b.setBasePosition(c.getBasePosition());
    return b;
    }

int faceSquanderLowerBound(parameterUse p)
	{
	int i ,total =0;
	for (i=0;i<faceSize();i++)
		{
		face f = faceAt(i);
		if (!f.getTempState())
		total += p.squanderFace(f.size());
		}
	return total;
	}

// this is for finished vertices:
int ExcessAt(vertex v,parameterUse p)
	{
	if (v.tempCount()>0) return 0;
	int t = v.triCount();
	int q = v.quadCount();
	if (v.exCount()>0) return p.pqrExcess(t,q,v.exCount()); 
	if (q+t>p.faceCountMax) return p.squanderTarget; 
	int u = p.squanderVertex[t][q]- q*p.squanderFace(4)-t*p.squanderFace(3);
	if (u<0)  /* generate error: */ 
		{ 
		System.out.print("negative : " +t+" "+q+" "+p.squanderVertex[t][q]);
		misc.throwIndex(u); }
	return u;
	}

private void DeleteAround(vertex w,Vector List,Vector Excesses)
	{
	int index;
	index = List.indexOf(w);
	if (index>=0) { List.removeElementAt(index); 
			Excesses.removeElementAt(index); }
	for (int i=0;i<w.size();i++)
		{
		vertex v = next(w,i);
		index = List.indexOf(v);
		if (index>=0) { List.removeElementAt(index); 
			Excesses.removeElementAt(index); }
		}
	for (int i=0;i<w.size();i++)
		if (w.elementAt(i).size()==4)
		{
		face F = w.elementAt(i);
		for (int j=0;j<4;j++) 
			{
			vertex v = F.elementAt(j);
			index = List.indexOf(v);
			if (index>=0) { List.removeElementAt(index); 
				Excesses.removeElementAt(index); }
			}
		}
	}

int ExcessNotAt(vertex v,parameterUse p)
	{
	Vector UseList = new Vector(vertexSize());
	Vector SqList = new Vector(vertexSize());
	vertex c;
 
	/* BuildUseList */ {
	for (int i=0;i<vertexSize();i++)
		{
		vertex w = vertexAt(i);
		int ex = ExcessAt(w,p);
		if (ex<0)  /* generate error: */ { misc.throwIndex(ex); }
		if (ex>0)
			{ UseList.addElement(w); SqList.addElement(new Integer(ex)); }
		}
	if (v!=null) DeleteAround(v,UseList,SqList);
	}
	/* mainLoop: */ {	
	int total=0;
	while (UseList.size()>0)
		{
		int maxpos=0, maxsq = 0;
		if (UseList.size()!=SqList.size()) {
			(new graphException("unequal lists")).report(); }
		for (int i=0;i<SqList.size();i++) 
			{ int r = ((Integer)SqList.elementAt(i)).intValue(); 
			  if (r>maxsq) { maxsq = r; maxpos = i; }
			}
		if (maxpos>UseList.size()) {
			(new graphException("unequal listws "+maxpos+ " "
					+UseList.size()+" "+SqList.size())).report(); }
		vertex w = (vertex) UseList.elementAt(maxpos);
		total += maxsq;
		DeleteAround(w,UseList,SqList);
		}
	return total;
	}
	}


int ScoreUpperBound(parameterUse p) 
	{
	int temp = 0;
	for (int i=0;i<faceSize();i++) 
		{
		face F = faceAt(i);
		if (F.getTempState()) 
				misc.warn("\nsub called with temps");
		temp += p.fixedScoreFace[F.size()];
		}
	// Note:
	// you can get a better version by adding in the ScorePenalties.
	return temp;
	}



boolean LooksOKToAdd(int tri,int quad,int excep,	
	int temp,vertex v,parameterUse p)
	// Unusual parameter passing.
	// tri = number of additional triangles.
	// quad = number of additional quads,
	// excep=0 means no exceptionals added.
	// excep=n>0 means one exceptional, an n-gon added at v.
	// temp = number of additional temps.
	{
	int t = tri+v.triCount();
	int q = quad+v.quadCount();
	int tempX = temp + v.tempCount();
	int e = v.exCount();
	if (excep>0) { e++; }
 
	if ((e>0) && (t+q+tempX+e > p.faceCountMaxAtExceptionalVertex))
		 return false;
	if ((e==0) &&(t+q+tempX+e > p.faceCountMax)) return false; 

	int sq = faceSquanderLowerBound(p) + (excep>0 ? p.squanderFace(excep) : 0)+
			tri*p.squanderFace(3)+quad*p.squanderFace(4);
	int excess= p.forecast(t,q,tempX)-t*p.squanderFace(3)- q*p.squanderFace(4);

	if (sq >= p.squanderTarget) {  return false; }
	if (sq + ExcessNotAt(null,p) >= p.squanderTarget) { return false; }
	if ((e==0)&&( (sq+ExcessNotAt(v,p)+excess >=p.squanderTarget))) 
		{  return false; }
	return true;
	}

boolean ForcedTriangleAt(vertex v,parameterUse p)
	{
	int t = v.triCount();
	int q = v.quadCount();
	int tempX = v.tempCount();
	int e = v.exCount();
	int fsq = faceSquanderLowerBound(p);
	int fsqred = fsq - q*p.squanderFace(4) - t*p.squanderFace(3);
	int target = p.squanderTarget;
	int excessNot = ExcessNotAt(v,p);
	if (e==0)
		{
		if ((fsq+p.squanderFaceStartingAt(5)>target)&&
		 (p.forecast(t,q+1,tempX-1)+fsqred+excessNot>target)&&
		 (p.forecast(t+tempX+1,q,0)+fsqred+excessNot>target))
		{
		/*debug*/{ String S=" vars:t q tempX e fsq targ excessNot "+
				" svA svB = \n"+t+" "+
				q+" "+tempX+" "+e+" "+fsq+" "+target+" "+excessNot+" "+
				p.forecast(t,q+1,tempX-1)+" "+
				p.forecast(t+tempX+1,q,0);
				discard.addLine(S); }
		return true;
		}
		else return false;
		}
	// now v.exCount() >0
	int nextface = p.squanderFaceStartingAt(4);
	if ((fsq+excessNot+nextface >= target)&&
		(t+tempX+q+e+1 > p.faceCountMaxAtExceptionalVertex)) 
		// could add a pqrExcess condition here to improve things //
		return true;
	return false;
	}

int selectVertex(face f)  /**** returns local index ****/
        {
	int minheight = 100000;
	int index = -1;
	for (int i=0;i<f.size();i++)
		{
		int h = (f.elementAt(i)).getHeight();
		if (h<minheight) { minheight = h; index = i; }
		}
	if (index <0) /* generate error: */ { misc.throwIndex(index); }
	return index;
	}
} // end class arrangePlus



