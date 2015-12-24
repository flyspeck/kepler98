package archive;

public class graphDispatch {

		private static int s1 = 0,
						   s2 = 2000,
						   s3 = 6000,
						   s4 = 8000,
						   s5 =10000,
						   s6 =12000,
						   s7 =14000,
						   s8 =16000; 
		private static String descriptor  = 
			"QUAD: (1-"+graphQ.data.length+") "+
		    "\nPENT: ("+(s2+1)+"-"+(s2+graphPENT.data.length)+") "+
		    "\nHEX: ("+(s3+1)+"-"+(s3+graphHEX.data.length)+") "+
		    "\nHEPT: ("+(s4+1)+"-"+(s4+graphHEPT.data.length)+") "+
		    "\nOCT: ("+(s5+1)+"-"+(s5+graphOCT.data.length)+") " +
			"\nDODEC: ("+(s6+1)+"-"+(s6+seanA.data.length)+") ";

	static public String DataBaseDescription()
		{
		return descriptor;
		}

		

	// We advertize the lists as starting with one
	// but grab starts at zero.
	static public GraphArrays grab(int i)
		{
		if (i<0) return null;
		else if (i<s2)
			{
			i -= 0;
			if (i<graphQ.data.length)
				return new GraphArrays(graphQ.data[i]);
			else return null;
			}

		else if (i>=s7)  // Sean Experimental
			{
			i -= s7;
			if (i<seanQ.data.length)
				return new GraphArrays(seanQ.data[i]);
			else return null;
			}

		else if (i>=s6)  // Sean Experimental
			{
			i -= s6;
			if (i<seanA.data.length)
				return new GraphArrays(seanA.data[i]);
			else return null;
			}

/*
		else if (i>=s7)  // PARTIV experimental cases.
			{
			i -= s7;
			if (i<graphSOMEPENT2.data.length)
				return new GraphArrays(graphSOMEPENT2.data[i]);
			else return null;
			}

		else if (i>=s6)  // PARTIV stuff (based on loopM2.m of 1/96)
			{
			i -= s6;
			if (i<temp.data.length)
				return new GraphArrays(temp.data[i]);
			else return null;
			}
*/

		else if (i>=s5)
			{
			i -= s5;
			if (i<graphOCT.data.length)
				return new GraphArrays(graphOCT.data[i]);
			else return null;
			}

		else if (i>=s4)  // PARTIV stuff (8/3/97 generation)
			{
			i -= s4;
			if (i<graphHEPT.data.length)
				return new GraphArrays(graphHEPT.data[i]);
			else return null;
			}

		else if (i>=s3)  // PARTIV stuff (8/3/97 generation)
			{
			i -= s3;
			if (i<graphHEX.data.length)
				return new GraphArrays(graphHEX.data[i]);
			else return null;
			}

		else if (i>=s2)  // PARTIV stuff (8/3/97 generation)
			{
			i -= s2;
			if (i<graphPENT.data.length)
				return new GraphArrays(graphPENT.data[i]);
			else return null;
			}

		else return null;
		}

		}


