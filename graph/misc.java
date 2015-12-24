
public class misc {
	public static int mod(int i,int j) 
		{ int r = i % j; return (r<0 ? r+j : r); }

	public static void warn(String s)
		{ System.err.println("Warning: "+s); }

	public static int throwIndex(int index)
		{  int a[] = new int[0]; 
			return a[index]; }
	}
