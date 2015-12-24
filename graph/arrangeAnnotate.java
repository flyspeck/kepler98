// this is for debugging:
class arrangeAnnotate extends arrangePlus {
	static { System.out.print(" -- arrangeAnnotate $Revision: 1.2 $ -- \n"); }
    String S;
    arrangeAnnotate(arrange a)
        {
        super(a);
        }
    void addLine(String S) { this.S = this.S +"\n"+ S; }
    }

