package archive;
import java.util.*;

/**
  * These are the arrays needed to store graph incidence data for planar graphs.
  * faceList is a list of faces faceList[0]={v1,v2,v3,..,ngon} firstface
  *  faceList[1]={v'1,...} second face.
  *
  * The tempList is a list of face numbers of the faces that are temporary.
  *
  * adjacent is a field in the Mathematica.96 format for COMBIN.
  *   it lists the adjacent vertices around vertex i in cyclic order
  *
  * vertexList lists the faces around vertex i in cyclic order.
  *
**/

/** format of graph*.java:
	first group is (1+n) elements, where n is the first element.
	n = number of temporary faces, followed by the facenumber of each temp.
	the main group follows next m+f1+f2+...fm
	m = number of faces total, followed by data fi for each face.
	fi = (1+li) elements, li is the constant in the first element,
		it is the number of vertices on the polygon,
		followed by the vertex numbers of each element. **/


public class GraphArrays {
	public int faceList[][];
	public int adjacent[][];
	public int vertexList[][];
	public int tempList[];
	public GraphArrays() {};
	public GraphArrays(String S)

		{

		// Call String Tokenizer;
		StringTokenizer Token 
		  = new StringTokenizer(S);
		int L = Token.countTokens();
		int tokarray[] = new int[L];
		int i,j;
		for (i=0;i<L;i++)
			tokarray[i]=Integer.valueOf(Token.nextToken()).
					intValue();
		int offset=0;
		tempList = new int[tokarray[0]];
		offset=1;
		for (i=0;i<tempList.length;i++)
			tempList[i]= tokarray[offset+i];
		offset += tokarray[0];
		faceList = new int[tokarray[offset]][];
		
		offset++;
		for (i=0;i<faceList.length;i++)
			{
			faceList[i] = new int[tokarray[offset]];
			offset++;
			for (j=0;j<faceList[i].length;j++)
			  // BUGBUG -1 was temporary till strings
			  faceList[i][j]=tokarray[offset+j];
			offset += faceList[i].length;
			}

		int adjacentCount=-1;
		for (i=0;i<faceList.length;i++)
		for (j=0;j<faceList[i].length;j++)	
			if (faceList[i][j]>adjacentCount)
				adjacentCount = faceList[i][j];
		adjacentCount++;

		adjacent = new int[adjacentCount][];
		vertexList = new int[adjacentCount][];
		int adjacentSizes[] = new int[adjacentCount];
		for (i=0;i<adjacentCount;i++) adjacentSizes[i]=0;
		for (i=0;i<faceList.length;i++)
		for (j=0;j<faceList[i].length;j++)
			adjacentSizes[faceList[i][j]]++;
		for (i=0;i<adjacentCount;i++) 
			{
			adjacent[i]=new int[adjacentSizes[i]];
			vertexList[i]=new int[adjacentSizes[i]];
			}
		int d[][][] = new int[adjacentCount][][];
		int len[]= new int[adjacentCount];
		for (i=0;i<adjacentCount;i++) len[i]=0;
		int i1,i2,i3,r;
		for (i=0;i<adjacentCount;i++) 
			d[i]=new int[adjacentSizes[i]][3];
		for (i=0;i<faceList.length;i++)
		for (j=0;j<faceList[i].length;j++)
			{
			r = faceList[i].length;
			i1 = faceList[i][j%r];
			i2 = faceList[i][(j+1)%r];
			i3 = faceList[i][(j+2)%r];
			d[i2][len[i2]][0]=i1;
			d[i2][len[i2]][1]=i3;
			d[i2][len[i2]][2]=i; // face ID.
			len[i2]++;
			}
		// now sort;
		for (i=0;i<adjacentCount;i++)
		for (r=0;r+1<adjacentSizes[i];r++)
		for (j=0;j<len[i];j++)
			{
			if (d[i][r][1]== d[i][j][0])
				{
				i1= d[i][j][0]; 
				d[i][j][0]=d[i][r+1][0];
				d[i][r+1][0]=i1;
				i1=d[i][j][1];
				d[i][j][1]=d[i][r+1][1];
				d[i][r+1][1]=i1;
				i1=d[i][j][2];
				d[i][j][2]=d[i][r+1][2];
				d[i][r+1][2]=i1;
				}
			}
		// now install adjacent
		for (i=0;i<adjacentCount;i++)
		for (j=0;j<adjacentSizes[i];j++)
			{
			adjacent[i][j] = d[i][j][0];
			vertexList[i][j]=d[i][j][2];
			}
		}

	static public void printgraph (String S)
		{
		GraphArrays gr = new GraphArrays(S);
		int i,j;
		System.out.print("Boundary faces: ");
		for (i=0;i<gr.tempList.length;i++) 
			System.out.print(" "+gr.tempList[i]);
		if (gr.tempList.length==0)
			System.out.print("(none)");
		for (i=0;i<gr.faceList.length;i++) {
			System.out.print("\nFace "+i+": ");
			for (j=0;j<gr.faceList[i].length;j++) 
				System.out.print(" "+gr.faceList[i][j]);
			}
		System.out.println();
		// now print adjacents
		for (i=0;i<gr.adjacent.length;i++)
			{
			System.out.print("\nVertex "+i+": ");
			for (j=0;j<gr.adjacent[i].length;j++)
				System.out.print(" "+gr.adjacent[i][j]);
			}
		System.out.println();
		}
	


	
			
	
	}

	
