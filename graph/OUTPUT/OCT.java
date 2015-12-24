import java.util.*;

/* graphOCT.java in this directory is output from OCT case of
	renderGraph.  The output has been modified to combine all the
	cases into a single class.  We need to check that these all
	appear in the archive.  That is the purpose of this program. */

/* Program that identifies each octagon from graphOCT.java
with a configuration in the archive.  pos0 gives the array
positions in the archive.  pos0 was calculated with Mathematica
based on the list of invariants. */

/* The output shows that everything was matched with something
	in the archive.

maya:hales% javac OCT.java
maya:hales% java OCT
 -- sort $Revision: 1.3 $ -- 
 -- arrange $Revision: 1.2 $ -- 
 -- arrangeInvariant $Revision: 1.2 $ -- 
done!
maya:hales% date
Wed Feb 11 20:19:08 EST 1998


*/


class OCT {
static { System.out.print(" -- sort $Revision: 1.3 $ -- \n"); }
 
static public void main ( String ac[] )
 
    {
    final String dataArchive[] = archive.graphOCT.data;
    final String dataRecent[]  = graphOCT.data;
	Stack newstack = new Stack();

	    for (int i=0;i<dataRecent.length;i++)
        {
        archive.GraphArrays u = new archive.GraphArrays(dataRecent[i]);
        arrangeInvariant a = new arrangeInvariant(u);
        a.setInvariant();
		int j = pos0[i];
		archive.GraphArrays v = new archive.GraphArrays(dataArchive[j]);
		arrangeInvariant b = new arrangeInvariant(v);
		b.setInvariant();
		if (a.getInvariant() != b.getInvariant())
			System.out.println("invariant mismatch : "+ i);
		if (!arrangeInvariant.Isomorphic(b,a))
			{
			System.out.println(","+i+" compared to "+j); 
			System.out.println(" "+arrange.ToString((arrange)a));
			}
        }
	System.out.println("done!");
	}

 
 

static int pos0[] = 

 {32, 18, 0, 1, 2, 33, 21, 3, 25, 26, 4, 28, 5, 6, 7, 12, 29, 14, 30, 34, 
   11, 13, 35, 36, 15, 16, 37, 38, 10, 39, 40, 41, 42, 43};

	};
