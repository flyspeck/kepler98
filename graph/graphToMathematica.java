import java.util.*;

class graphToMathematica {

		static public void main(String ac[]) 
		{

	//for (int i=10000;i<10044;i++)
	for (int i=0;i<1749;i++)
	{
	archive.GraphArrays gr = archive.graphDispatch.grab(i);
	arrange a = new arrange(gr);
	String S = arrange.MathematicaString(a);
	System.out.print(S+",\n");
	};

	}
	  }
