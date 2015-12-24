
class parameterUse extends parameters {

	static { System.out.print(" -- parameterUse $Revision: 1.2 $ -- \n"); }

private static int squanderFaceD[];
static int maxGon() { return squanderFaceD.length -1; }

int squanderVertex[][]; // update fixed squander vertex.

static int squanderFace(int ngon) { 
	if (ngon>=squanderFaceD.length) return squanderTarget;
	else return squanderFaceD[ngon]; }

static int squanderFaceStartingAt(int ngon)
	{
	int current = squanderTarget;
	for (int i=ngon;i<squanderFaceD.length;i++) if (current<squanderFaceD[i])
			{ current = squanderTarget; }
	return current;
	}

/* forecastD is accurate only at vertices that will eventually
	be completed without any exceptionals.  It gives the minimum
	squander given that there are already at least i,j,k 
	triangles, quads, and temps there */
int forecast(int tri,int quad,int temp) 
	{ 
	if ((tri>=forecastD.length)||
		(quad>=forecastD.length)||
		(temp>=forecastD.length))
	return squanderTarget;
	else return forecastD[tri][quad][temp]; }

private int forecastD[][][]; 

parameterUse(int ngon)
	{
	if (ngon<5) 
		{
		System.out.println("Error: the minimum number of sides is 5"+
			"\nfor the ngon construction.  Using ngon=5....");
		ngon=5;
		}
	if (ngon>fixedSquanderFace.length-1)
		{
		System.out.println("Error: the maximum number of sides is "+
			(fixedSquanderFace.length-1)+
			"\nUsing this for ngon...");
		ngon= fixedSquanderFace.length-1;
		}
	int Q = fixedSquanderVertex.length;
	int i,j,k,m;
	squanderVertex =new int[Q][Q];
	for (i=0;i<Q;i++) for (j=0;j<Q;j++)
		squanderVertex[i][j]=fixedSquanderVertex[i][j];
	squanderFaceD = new int[fixedSquanderFace.length];
	forecastD = new int[Q][Q][Q];
	for (i=0;i<squanderFaceD.length;i++) squanderFaceD[i] = 
		(i>ngon ? squanderTarget : fixedSquanderFace[i]);
	for (i=0;i<squanderFaceD.length;i++) if ((i!=ngon)&&(
		squanderFaceD[i]+squanderFaceD[ngon]>=squanderTarget))
		squanderFaceD[i]=squanderTarget;
	for (i=0;i<Q;i++)
	for (j=0;j<Q;j++)
		{
		forecastD[i][j][0]=fixedSquanderVertex[i][j];
		}
	for (i=Q-1;i>=0;i--)
	for (j=Q-1;j>=0;j--)
		{
		// make list monotonic.
		if ((i>0)&&(forecastD[i-1][j][0]> forecastD[i][j][0]))  
			forecastD[i-1][j][0]= forecastD[i][j][0];
		if ((j>0)&&(forecastD[i][j-1][0]>forecastD[i][j][0]))
			forecastD[i][j-1][0]=forecastD[i][j][0];
		}
	for (i=0;i<Q;i++)
	for (j=0;j<Q;j++)
	for (k=1;k<Q;k++)
		{
		forecastD[i][j][k]=squanderTarget;
		if (i+j+k<Q) for (m=0;m<=k;m++)
		if ((i+m<Q)&&(j+k-m<Q)&&(forecastD[i+m][j+k-m][0]<forecastD[i][j][k]))
			forecastD[i][j][k]=forecastD[i+m][j+k-m][0];
		}
	}

private int tCount(int gon)
	{
	int temp=0;
	for (int i=0;i<quadCases[gon].length;i++) 
	if (quadCases[gon][i]==0) temp++;
	return temp;
	}

private int qCount(int gon)
	{
	int temp=0;
	for (int i=0;i<quadCases[gon].length;i++) 
	if (quadCases[gon][i]==1) temp++;
	return temp;
	}

parameterUse(int ngon,boolean ignore)
	{
	squanderFaceD = new int[5]; // go only up to quads.
	for (int i=0;i<squanderFaceD.length;i++) 
		squanderFaceD[i]=fixedSquanderFace[i];
	int Q = fixedSquanderVertex.length;
	forecastD = new int[Q][Q][Q];
	int i,j,k,m;

	squanderVertex =new int[Q][Q];
	for (i=0;i<Q;i++) for (j=0;j<Q;j++)
		squanderVertex[i][j]=fixedSquanderVertex[i][j];
	int nTri = tCount(ngon),nQuad= qCount(ngon);
	for (i=0;i<ngon;i++) 
		{
		// we assume that cases (t,q) earlier on quadCase have been done.
		// so we exclude them from this pass.
		if ((nTri!= tCount(i))||(nQuad!=qCount(i)))
			squanderVertex[tCount(i)][qCount(i)]= squanderTarget;
		}

	for (i=0;i<Q;i++) for (j=0;j<Q;j++)
		forecastD[i][j][0]=squanderVertex[i][j];

	for (i=Q-1;i>=0;i--)
	for (j=Q-1;j>=0;j--)
		{
		// make list monotonic.
		if ((i>0)&&(forecastD[i-1][j][0]> forecastD[i][j][0]))  
			forecastD[i-1][j][0]= forecastD[i][j][0];
		if ((j>0)&&(forecastD[i][j-1][0]>forecastD[i][j][0]))
			forecastD[i][j-1][0]=forecastD[i][j][0];
		}
	for (i=0;i<Q;i++)
	for (j=0;j<Q;j++)
	for (k=1;k<Q;k++)
		{
		forecastD[i][j][k]=squanderTarget;
		if (i+j+k<Q) for (m=0;m<=k;m++)
		if ((i+m<Q)&&(j+k-m<Q)&&(forecastD[i+m][j+k-m][0]<forecastD[i][j][k]))
			forecastD[i][j][k]=forecastD[i+m][j+k-m][0];
		}
	}

int pqrExcess(int tCount,int qCount, int exCount)
	{
	if (tCount+qCount+exCount!= faceCountMaxAtExceptionalVertex)
			return 0;
	return excessTCount[tCount];
	}
}


