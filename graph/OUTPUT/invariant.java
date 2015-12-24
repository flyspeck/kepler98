import java.util.*;

// print out the invariants in a list. //
class invariant {
static { System.out.print(" -- invariant $Revision$ -- \n"); }

static public void main ( String ac[] )

	{
	// adapt as needed...
	final String dataRecent[]  = graphOCT.data;



	Stack stack = new Stack();

	System.out.println("newlist size = "+dataRecent.length);
	for (int i=0;i<dataRecent.length;i++)
		{
		archive.GraphArrays u = new archive.GraphArrays(dataRecent[i]);
		arrangeInvariant a = new arrangeInvariant(u);
		a.setInvariant();
		System.out.print(", "+ a.getInvariant());
		if (0 == i % 15 ) System.out.println();
		}

	}
};
