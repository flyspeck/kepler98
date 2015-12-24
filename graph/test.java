import java.util.*;

class test {


		static public void main(String ac[]) 
		{

	// misc.throwIndex(0);

	/* get list without repeats */  if (false) {
	Stack terminal = new Stack();
	for (int i=12000;i<13398;i++)
	{
	archive.GraphArrays gr = archive.graphDispatch.grab(i);
	arrange a = new arrange(gr);
	arrangeInvariant b = new arrangeInvariant(a);
	if (b.isContainedIn(terminal)) continue;
	terminal.push(b);
	}
	arrange.dumpStack(terminal);
	terminal.removeAllElements();
	}


	/* print invariants */ if (true) {
	//for (int i=12000;i<13066;i++)
	for (int i=12000;i<12001;i++)
	// for (int i=6000;i<6005;i++)
	{
	archive.GraphArrays gr = archive.graphDispatch.grab(i);
	arrange a = new arrange(gr);
	arrangeInvariant b = new arrangeInvariant(a);
	String S = arrange.MathematicaString(a);
	System.out.print(" {"+(i+1)+","+b.getInvariant()+"},\n");
	};
	}

	if (false)
	System.out.println(archive.seanQ.data.length);

	if (false)
	for (int i=0;i<parameters.quadCases[11].length;i++)
			System.out.println(parameters.quadCases[11][i]);

	/* parametertest*/ if (false) {
	for (int i=5;i<10;i++)
	{
	parameterUse p = new parameterUse(i);
	for (int j=0;j<10;j++)
		 System.out.println(" "+i+" "+j+" "+p.squanderFace(j));
	}}



	}

	  }
