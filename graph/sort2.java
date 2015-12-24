import java.util.*;

class sort2 {
static { System.out.print(" -- sort $Revision: 1.2 $ -- \n"); }

static public void main ( String ac[] )

	{
	final String dataArchive[] = graphPENT.data;
	final String dataRecent[]  = archive.graphPENT.data;

	Stack stack = new Stack();

	System.out.println("size = "+dataArchive.length);
	System.out.println("recent size = "+dataRecent.length);
	for (int i=0;i<dataArchive.length;i++)
		{
		archive.GraphArrays u = new archive.GraphArrays(dataArchive[i]);
		arrangeInvariant a = new arrangeInvariant(u);
		a.setInvariant();
		if (!a.isContainedIn(stack)) 
			stack.push(a);
		}

	Stack newstack = new Stack();
	for (int i=0;i<dataRecent.length;i++)
		{
		arrangeInvariant a = new arrangeInvariant(new archive.GraphArrays(
				dataRecent[i]));
		a.setInvariant();
		if (a.isContainedIn(stack)) 
			{
			System.out.println("\""+arrange.ToString((arrange)a)+"\",");
			}
		}
	
	System.out.println("done! ");
	}
};
