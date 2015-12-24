import java.util.*;

/* graphHEX.java in this directory is output from HEX case of
	renderGraph.  The output has been modified to combine all the
	cases into a single class.  We need to check that these all
	appear in the archive.  That is the purpose of this program. */

/* Program that identifies each pentagon from graphHEX.java
with a configuration in the archive.  pos0 gives the array
positions in the archive.  pos0 was calculated with Mathematica
based on the list of invariants. */

/* The output shows that everything was matched with something
	in the archive.

maya:hales% java HEX
 -- sort $Revision: 1.3 $ -- 
 -- arrange $Revision: 1.2 $ -- 
 -- arrangeInvariant $Revision: 1.2 $ -- 
done!
maya:hales% date
Wed Feb 11 20:01:28 EST 1998

*/


class HEX {
static { System.out.print(" -- sort $Revision: 1.3 $ -- \n"); }
 
static public void main ( String ac[] )
 
    {
    final String dataArchive[] = archive.graphHEX.data;
    final String dataRecent[]  = graphHEX.data;
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
  {178, 184, 185, 193, 196, 203, 214, 217, 225, 229, 230, 232, 233, 234, 236, 
   240, 5, 243, 245, 247, 248, 249, 250, 251, 252, 255, 256, 257, 259, 262, 
   263, 264, 265, 266, 267, 269, 271, 277, 280, 282, 284, 286, 289, 290, 291, 
   292, 293, 294, 295, 296, 297, 298, 300, 301, 302, 303, 304, 305, 306, 8, 
   308, 309, 310, 311, 312, 313, 11, 323, 13, 324, 325, 326, 327, 328, 329, 
   330, 331, 19, 20, 21, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 
   36, 39, 40, 41, 42, 43, 45, 46, 49, 50, 51, 52, 54, 55, 56, 57, 58, 59, 
   60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 333, 336, 
   338, 339, 80, 82, 84, 85, 348, 349, 86, 87, 88, 89, 352, 93, 94, 95, 96, 
   97, 98, 357, 101, 103, 104, 105, 106, 361, 107, 108, 110, 362, 363, 112, 
   113, 114, 115, 116, 117, 118, 119, 120, 121, 364, 122, 123, 124, 125, 126, 
   365, 367, 368, 127, 128, 129, 130, 369, 131, 132, 133, 370, 134, 135, 371, 
   372, 136, 137, 138, 139, 140, 141, 143, 374, 375, 376, 144, 145, 377, 146, 
   379, 380, 381, 383, 147, 148, 149, 384, 385, 386, 150, 151, 152, 153, 154, 
   387, 388, 389, 390, 155, 156, 391, 393, 45, 397, 398, 399, 400, 401, 402, 
   159, 160, 407, 409, 410, 411, 412, 161, 413, 414, 416, 166, 417, 167, 168, 
   169, 170, 418, 419, 171, 172, 173, 174, 422, 423, 424, 175, 427}
;

	};
