import java.util.*;

// typical use:
/*
	1.  pick an archive file graphA.java for the first section and a new list
			of cases for the graphR.java in the second sectoin
	2.  run the program saving output to graphNEW.java
	3.  cat graphNEW.java >> graphA.java, remove the seam to get bigger
			archive file.
	4.  delete the files graphX.java graphNEW.java
*/

class sort {
static { System.out.print(" -- sort $Revision: 1.3 $ -- \n"); }

static public void main ( String ac[] )

	{
	final String dataArchive[] = archive.graphPENT.data;
	final String dataRecent[]  = graphPENT.data;


	// final String dataArchive2[] = archive.graphHEX.data;

	Stack stack = new Stack();

	System.out.println("archive size = "+dataArchive.length);
	// System.out.println("size2 = "+dataArchive2.length);
	System.out.println("newlist size = "+dataRecent.length);
	for (int i=0;i<dataArchive.length;i++)
		{
		archive.GraphArrays u = new archive.GraphArrays(dataArchive[i]);
		arrangeInvariant a = new arrangeInvariant(u);
		a.setInvariant();
		if (!a.isContainedIn(stack)) 
			stack.push(a);
		}

	/*
	for (int i=0;i<dataArchive2.length;i++)
		{
		archive.GraphArrays u = new archive.GraphArrays(dataArchive2[i]);
		arrangeInvariant a = new arrangeInvariant(u);
		a.setInvariant();
		if (!a.isContainedIn(stack)) 
			stack.push(a);
		}
	*/
	Stack newstack = new Stack();
	for (int i=0;i<dataRecent.length;i++)
		{
		arrangeInvariant a = new arrangeInvariant(new archive.GraphArrays(
				dataRecent[i]));
		a.setInvariant();
		if (!a.isContainedIn(stack)) 
			{
			newstack.push(a);
			stack.push(a);
			}
		}
	System.out.println("//// updated archive size = "+stack.size());
	System.out.println("//// new discoveries: = "+newstack.size());
	arrange.dumpStack(newstack);

	}
};
