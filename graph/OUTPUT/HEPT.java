import java.util.*;

/* graphHEPT.java in this directory is output from HEPT case of
	renderGraph.  The output has been modified to combine all the
	cases into a single class.  We need to check that these all
	appear in the archive.  That is the purpose of this program. */

/* Program that identifies each pentagon from graphHEPT.java
with a configuration in the archive.  pos0 gives the array
positions in the archive.  pos0 was calculated with Mathematica
based on the list of invariants. */

/* The output shows that everything was matched with something
	in the archive.

maya:hales% java HEPT
 -- sort $Revision: 1.3 $ -- 
 -- arrange $Revision: 1.2 $ -- 
 -- arrangeInvariant $Revision: 1.2 $ -- 
done!
maya:hales% date
Wed Feb 11 20:10:22 EST 1998


*/


class HEPT {
static { System.out.print(" -- sort $Revision: 1.3 $ -- \n"); }
 
static public void main ( String ac[] )
 
    {
    final String dataArchive[] = archive.graphHEPT.data;
    final String dataRecent[]  = graphHEPT.data;
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
  {127, 131, 125, 98, 157, 177, 221, 128, 211, 212, 14, 26, 35, 42, 45, 53, 
   306, 307, 69, 74, 79, 324, 309, 320, 326, 112};


	};
