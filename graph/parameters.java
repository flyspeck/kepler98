// Put pass dependent constraints here.
// Kepler conjecture parameters

class parameters {
	static { System.out.print(" -- parameters $Revision: 1.5 $ -- \n"); }
	// true means that we are excluding a vertex with only two
	// standard clusters, a qrtet and a pentagonal.
	final static boolean excludePentQRTet = true;
	final static int vertexCountMin = 1; // There must be at least this many
	final static int vertexCountMax = 100; // And at most this many vertices.
	final static int faceCountMaxAtExceptionalVertex= 5; 
			// This is the maximum number
			// of faces that can be at a face with something
			// with an exceptional face.
	final static int faceCountMax = 6;

	// values multiplied by 1000 and rounded down to nearest integer;
	protected final static int 
		fixedSquanderFace[]={0,0,0,0,2378,4896,7414,9932,10916};
	protected final static int 
		fixedScoreFace[] ={1000,1000,1000,1000,0,-1030,-2060,-3030,-3030};
	final static int squanderTarget = 14800;
	private final static int x = 14800; // any value over squander Target.
	final static int scoreTarget = 8000;

	// TCount[count] is the excess around a vertex at an exceptional cluster
	// having count triangles, and faceCountMaxAtExceptionalVertex-count
	// nontriangles. The length of the array must be fCMAEVertex.
	protected final static int
		excessTCount[] = {x,x,x,1400,1500};
	protected final static int fixedSquanderVertex[][] = // type (i,j).
		// PartIII. table 5.1:
		// This must be a square array .
		{{x,x,x,7135,10649,x,x},
		 {x,x,6950,7135,x,x,x},
		 {x,8500,4756,12981,x,x,x},
		 {x,3642,8334,x,x,x,x},
		 {4139,3781,x,x,x,x,x},
		 {550,11220,x,x,x,x,x},
		 {6339,x,x,x,x,x,x}
		};

	final static int quadCases[][] = {
			{0,0,1,1,1}, //0
			{0,1,0,1,1}, //1
			{0,0,0,0,0,1}, //2
			{1,1,1,1}, //3
			{0,0,1}, //4
			{0,0,0,1,1},
			{0,0,1,0,1},
			{0,1,1},
			{0,1,1,1}, //8
			{1,1,1},
			{0,0,0,0,0,0},
			{0,0,1,1}, //11
			{1,0,1,0},
			{0,0,0,1},
			{0,0,0,0},
			{0,0,0,0,1},
			{0,0,0,0,0}
			};

	static void printParameters()
		{
		int i,j;
		System.out.println("***\n*  Parameters for graph generator");
		System.out.println("* vertexCountMin = "+vertexCountMin);
		System.out.println("* vertexCountMax = "+vertexCountMax);
		System.out.println("* faceCountMaxAtExceptionalVertex = "+
			faceCountMaxAtExceptionalVertex);
		System.out.print("* fixedSquanderFace = {");
		for (i=0;i<fixedSquanderFace.length;i++)
			System.out.print("("+i+","+fixedSquanderFace[i]+"),");
			System.out.println("}");
		System.out.println("* max # of edges on a face = "+
			fixedSquanderFace.length);

		System.out.print("* fixedScoreFace = {");
		for (i=0;i<fixedScoreFace.length;i++)
			System.out.print("("+i+","+fixedScoreFace[i]+"),");
			System.out.println("}");
		System.out.println("* scoreTarget = "+scoreTarget);
		System.out.print("* squanderTarget = "+squanderTarget);

		System.out.print("\n* fixedSquanderVertex = ");
			for (i=0;i<fixedSquanderVertex.length;i++)
			{
			for (j=0;j<fixedSquanderVertex[i].length;j++)
			if (fixedSquanderVertex[i][j]!=x)
			System.out.print("("+i+","+j+","
				+fixedSquanderVertex[i][j]+")");
			System.out.print("\n*  ");
			}
		/*
		System.out.print("\n* fixedScoreVertex = ");
			for (i=0;i<fixedScoreVertex.length;i++)
			{
			for (j=0;j<fixedScoreVertex[i].length;j++)
			if (fixedScoreVertex[i][j]!=0)
			System.out.print("("+i+","+j+","
				+fixedScoreVertex[i][j]+")");
			System.out.print("\n*  ");
			}
		*/

		System.out.println("\n***");
		}


	static {
                int i,j;
                int r = fixedSquanderVertex.length;
                //"There are at most 6 faces around each Vertex"
                if (r!=faceCountMax+1) System.out.print("faceCountMax required");
 
                for (i=0;i<r;i++)
                        if (r!= fixedSquanderVertex[i].length) 
			System.out.println("Error:"+i+" matrix is not square"+
					r+" "+fixedSquanderVertex[i].length);
			/*
                if (r!=fixedScoreVertex.length) 
			System.out.println("Error: matrices of unequal size");
                for (i=0;i<r;i++)
                        if (r!= fixedScoreVertex[i].length) 
			System.out.println("Error2: matrix is not square");
			*/
                if (vertexCountMin>vertexCountMax) 
			System.out.println("Error: min>max");
                if (vertexCountMin<0) 
			System.out.println("Error: min<0");
                if (fixedSquanderFace.length!=fixedScoreFace.length) 
			System.out.println("Error: unequal faces ");
                if (fixedSquanderFace.length>9) 
			System.out.println("Error: nonagons not allowed");
                for (i=0;i<fixedSquanderFace.length-1;i++)
                        if (fixedSquanderFace[i]>fixedSquanderFace[i+1]) 
			System.out.println("Error: need monotonicity");
                for (i=0;i<fixedScoreFace.length-1;i++)
                        if (fixedScoreFace[i]<fixedScoreFace[i+1]) 
			System.out.println("Error: need score monotonicity");
			if (fixedSquanderFace.length<=5)
				System.out.println("Error: need at least pentagons");
                }

	}
