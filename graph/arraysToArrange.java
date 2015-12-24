/**
  *  Conversion from GraphArrays format for a graph to
  *  an arrange.
  */

import java.util.*;

class arraysToArrange {
	static { System.out.print(" -- arraysToArrange $Revision: 1.2 $ -- \n"); }
	final static public arrange make(archive.GraphArrays gr)
		{
		int i,j;
		arrange a = new arrange();
		a.faceList = new Vector(gr.faceCycle.length);
		arrange.vertex v[] = new arrange.vertex[gr.faceCycle.length];
		for (i=0;i<gr.faceCycle.length;i++)
			{
			v[i] = new arrange.vertex();
			v[i].setHeight(0);
			v[i].incidentFace = new Vector(gr.faceCycle[i].length);
			}
		a.vertexList = new Vector(gr.GBLregions.length);
		arrange.face f[] = new arrange.face[gr.GBLregions.length];
		for (i=0;i<gr.GBLregions.length;i++)
			{
			f[i] = new arrange.face();
			f[i].isBoundary = false;
			f[i].incidentVertex = new Vector(gr.GBLregions[i].length);
			for (j=0;j<gr.GBLregions[i].length;j++)
				f[i].incidentVertex.addElement
				  (v[gr.GBLregions[i][j]]);
			}

		for (i=0;i<gr.faceCycle.length;i++)
		for (j=0;j<gr.faceCycle[i].length;j++)
			{
			v[i].incidentFace.addElement
				(f[gr.faceCycle[i][j]]);
			}
		for (i=0;i<gr.boundary.length;i++)
			{
			f[gr.boundary[i]].isBoundary=true;
			}
		return a;
		}
	}	
