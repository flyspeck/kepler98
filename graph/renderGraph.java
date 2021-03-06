import java.util.*;
import java.applet.Applet;
import java.awt.Graphics;
import java.awt.Color;
import java.awt.*;

public class renderGraph extends Applet implements Runnable {
/***** stuff follows ******/

	static { System.out.print(" -- renderGraph $Revision: 1.3 $ -- \n"); }

	static boolean debugMode = true;

graph planar = new graph();

public String getAppletInfo()
	{
	return 
	"Combinatorics for the Kepler Conjecture\n"+
	"by Thomas C. Hales\n"+
	"July 1997";
	}

Thread thread;
public void start() { thread = new Thread(this);
	thread.start();
	}

public void stop() { thread.stop(); }

public void destroy() {}
public void run() {
	while (true) 
	   {
	  try {Thread.currentThread().sleep(100); }
	  catch (InterruptedException ie) {}
	  if (!freeze.getState()) 
		{ for (int i=0;i<5;i++) 
			if (planar.Move()) repaint();
			else freeze.setState(true);  
		}
	  }
	}

private TextField text1 = new TextField(9);
private TextField intext = new TextField(7);
private TextField text2 = new TextField(13);
private TextField inface = new TextField(7);
private TextArea programInfo = 
	new TextArea(archive.graphDispatch.DataBaseDescription(),5,20);
private Button loop = new Button("loop");
private Button pop = new Button("peek");
private Button test = new Button("test");
private Checkbox check_vertex = new Checkbox("vertices");
private Checkbox check_face = new Checkbox("faces");
private Checkbox check_flip = new Checkbox("flip");
Checkbox freeze = new Checkbox("freeze");
private Choice familyChoice = new Choice();
private Choice incrementChoice = new Choice();

public void init()
	{
	String incChoices[] = {"100","25","5","1","-1","-5","-25","-100"};
	for (int i=0;i<incChoices.length;i++) 
		incrementChoice.addItem(incChoices[i]);

	String famChoices[]={"III","quad","pent","hex","hept","oct"};
	for (int i=0;i<famChoices.length;i++)
		familyChoice.addItem(famChoices[i]);

	add(programInfo);
	add(incrementChoice); 	add(text1); 		add(intext);
	// add(text2); 			add(inface);		
	add(freeze);
	if (debugMode) { add(loop); add(pop);}	
	add(check_vertex); 	
	add(check_face); 		add(check_flip);	
	if (debugMode) { add(familyChoice);  }
	if (debugMode) { add(test); }

	text1.setEditable(false);
	text2.setEditable(false);
	programInfo.setEditable(false);

	planar.install(2);
	text1.setText("graph "+(planar.getSelectedIndex()+1));
	}

/* was in setcoords 
	text2.setText("missing face");
*/

public synchronized boolean action(Event e,Object o)
	{
	if (e.target==loop)
		{
		System.out.println("\nStarting new loop");
		System.out.println("Stacks = "+tentative.size()+" "+preterminal.size());
		if (tentative.size()==0) return true;
		arrangeMain.discard = new arrangeAnnotate((arrange)tentative.peek());
		planar.setCoords(arrangeMain.discard.makeGr(),true);
		freeze.setState(false);
		pop.setLabel("peek");
		repaint();
		arrangeMain.loop(tentative,preterminal,p);
		System.out.println(arrangeMain.discard.S);
		}
	if (e.target==test)
		{
		//parameterUse p = new parameterUse(0,true);
		parameterUse p = new parameterUse(5);
		arrange a = arrange.makeNgon(5);
		arrangeMain.discard = new arrangeAnnotate(a);
		Stack S = new Stack();
		/* // for quads
		for (int localvar = 0;localvar<p.quadCases.length;localvar++)
			{
			a = arrange.seed(p.quadCases[localvar]);
			S.push(a);
			}
		*/
		S.push(new arrangePlus(arrange.makeNgon(5)));
		planar.use(S);
		planar.install(0);
		freeze.setState(false);
		repaint();
		}
	if (e.target==pop)
		{
		if (tentative.size()<2) return true;
		System.out.println("\nPopping graph");
		if (pop.getLabel().equals("pop")) tentative.pop();
		pop.setLabel("pop");
		System.out.println("Stacks = "+tentative.size()+" "+preterminal.size());
		arrangeMain.discard = new arrangeAnnotate((arrange)tentative.peek());
		planar.setCoords(arrangeMain.discard.makeGr(),true);
		freeze.setState(false);
		repaint();
		}
	if (e.target==intext)
		{
		int i;
		try {
		i = Integer.valueOf(intext.getText()).intValue();
			}
		catch (NumberFormatException u) { return true; }
		planar.install(i-1);
		text1.setText("graph "+(planar.getSelectedIndex()+1));
		freeze.setState(false);
		repaint();
		}
	if (e.target==inface)
		{
		int i;
		try { i = Integer.valueOf(inface.getText()).intValue(); }
		catch (NumberFormatException u) { return true; }
		i--; // usr and java numbers differ;
		planar.setHiddenRegion(i);
		freeze.setState(false);
		repaint();
		}
	if (e.target==check_vertex)
		{
		repaint();
		}
	if (e.target==check_face)
		{
		repaint();
		}
	if (e.target==check_flip)
		{
		// flip
		arrangeInvariant q = new arrangeInvariant(planar.gr);
		q.flip();
		planar.gr = q.makeGr();
		planar.setCoords(planar.gr,true);
		freeze.setState(false);
		repaint();
		}
	if (e.target==incrementChoice)
		{
		String S = incrementChoice.getSelectedItem();
		int j=0;
		try { j = Integer.parseInt(S); }
		catch (NumberFormatException ee) { System.out.println(ee.toString()); }
		planar.install(planar.getSelectedIndex()+j); 
		text1.setText("graph "+(planar.getSelectedIndex()+1));
		freeze.setState(false);
		repaint();
		}
	if (e.target==familyChoice)
		{
		int i = familyChoice.getSelectedIndex();
		switch (i) {
			case 0 : planar.use(null); 
					 planar.install(2);
					 freeze.setState(false);
					 text1.setText("graph "+(planar.getSelectedIndex()+1));
					 break; 
			case 2 : case 3 : case 4 :
			case 5 : System.out.println("// generating..."
						+familyChoice.getSelectedItem());
					{
					generator gen = new generator(3+i);
					Thread t = new Thread(gen);
					t.setPriority(Thread.MIN_PRIORITY);
					t.start();
					}
					break;
			case 1 : System.out.println("// generating ..."+
					familyChoice.getSelectedItem());
					{
					arrange a = arrange.makeNgon(3);
					arrangeMain.discard = new arrangeAnnotate(a);
					quadGenerator gen = new quadGenerator();
					Thread t = new Thread(gen);
					t.setPriority(Thread.MIN_PRIORITY);
					t.start();
					}
					break;

			}
		repaint();
		}
	return true;
	}


class generator implements Runnable {
	int ngon;
	generator(int i) { ngon = i; }
	public void run() { 
		Stack terminal = new Stack();
		terminal = arrangeMain.generateAll(ngon);
		arrange.dumpStack(terminal);
		System.out.print("Generator complete. All "+ngon+"-gons determined\n");
		}
	}

class quadGenerator implements Runnable {
	public void run() {
	arrange a;
	for (int j =0;j<parameterUse.quadCases.length;j++)
		{
		Stack terminal = new Stack();
		System.out.println("// Running case "+j);
		Stack S = new Stack();
		a = arrange.seed(parameterUse.quadCases[j]);
		S.push(a);
		planar.use(S);
		planar.install(0);
		freeze.setState(false);
		repaint();
		terminal = arrangeMain.generateQuad(j);
		arrange.dumpStack(terminal);
		}
	System.out.print("Generator complete. All quad stuff determined\n");
	}
	}

// added 3/29/98 to make face numbering compatible with that in Mathematica.
public long magnitude(int region[])
	{
	int L = region.length;
	int min,minpos;
	min=10000; minpos=0;
	int i;
	for (i=0;i<L;i++) if (region[i]<min) { min=region[i]; minpos=i; }
	int r2[] = new int[L];
	for (i=0;i<L;i++) r2[i]= region[(L-i+minpos) % L];
	// I don't think Mathematica does reversals.
	//if (r2[1]>r2[L-1]) for (i=0;i<L;i++)
		//r2[i]= region[(minpos+L-i) % L];
	long t=0;
	for (i=0;i<L;i++) { t = 20*t + (1+r2[i]); }
	return t;
	}

// added 3/29/98...
public int rankIt(int i,long mags[])
	{
	long m = mags[i];
	int count =0;
	for (int j=0;j<mags.length;j++)
		{
		if (m>mags[j]) count++;
		}
	return count;
	}


public void paint( Graphics g)
	{
	int i,j;
	archive.GraphArrays gr = planar.gr;
	int b = planar.getHiddenRegion();
	for (i=0;i<gr.tempList.length;i++)
		{
		j = gr.tempList[i];
		if (j!=b) { render.drawFace(g,j,gr.faceList,planar.coords); }
		else render.markFace(g);
		}
	if (check_face.getState()) 
		{
		long mags[] = new long[gr.faceList.length];
		for (i=0;i<gr.faceList.length;i++)
			mags[i]=magnitude(gr.faceList[i]);
		for (i=0;i<gr.faceList.length;i++)
		render.numberFace(g,i,gr.faceList,planar.coords,b,rankIt(i,mags));
		}
	for (i=0;i<gr.adjacent.length;i++)
	for (j=0;j<gr.adjacent[i].length;j++)
		{
		if (i<gr.adjacent[i][j]) render.drawEdge(g,i,gr.adjacent[i][j],planar.coords);
		}
	for (i=0;i<gr.adjacent.length;i++)
		render.drawDot(g,i,planar.coords,check_vertex.getState());
	}



/* for viewing stacks stuff */
Stack tentative;
Stack preterminal;
parameterUse p = new parameterUse(11,true);
	{
	arrange a = arrange.makeNgon(3);
	arrangeMain.discard = new arrangeAnnotate(a);
	tentative = new Stack();
	preterminal = new Stack();
	tentative.push(new arrangePlus(arrange.seed(parameters.quadCases[11])));
	}

} /* end of renderGraph */






/**
  	* class render  
	*
**/


class render {

final static int centerx=170,centery=350,width=140;

static int xcord(double t) { return (int) (centerx + width*t); }
static int ycord(double t) { return (int) (centery + width*t); }

static void drawEdge(Graphics g,int i,int j,double xy[][])
	{
	g.setColor(Color.black);
	g.drawLine(xcord(xy[i][0]),ycord(xy[i][1]),
		xcord(xy[j][0]),ycord(xy[j][1]));
	}

static void drawDot(Graphics g,int i,double xy[][],boolean drawNumber)
	{
	g.setColor(Color.blue);
	int x = xcord(xy[i][0]), y = ycord(xy[i][1]);
	g.fillOval(x-5,y-5,10,10);
	if (drawNumber) g.drawString(" "+(i+1),x+3,y+3);
	}

static void drawFace(Graphics g,int i,int region[][],double xy[][])
	{
	g.setColor(Color.yellow);
	Polygon p = new Polygon();
	int r[] = region[i];
	for (int j=0;j<r.length;j++)
	  p.addPoint(xcord(xy[r[j]][0]),ycord(xy[r[j]][1]));
	g.fillPolygon(p);
	}

static void markFace(Graphics g)
	{
	g.setColor(Color.yellow);
	g.fillRect(xcord(1.0)-15,ycord(1.0)-15,15,15);
	}


static void numberFace(Graphics g, int i,int region[][],double xy[][],int bd,
	int rk)
	{
	g.setColor(Color.black);
	int r[]= region[i];
	double x=0,y=0;
	for (int j=0;j<r.length;j++)
		{ x += xy[r[j]][0]; y+= xy[r[j]][1]; } 
	x /= r.length; y /= r.length;
	if (i==bd) { x = 1.0; y= 1.0; }
	g.drawString(" "+(rk+1),xcord(x)-10,ycord(y));
	}
}

/**
	*
	* class graph
	*
**/

class graph {
	private Stack InUse = null;
	private int index = 0;
	archive.GraphArrays gr;
	double coords[][];
	int getSelectedIndex() { return index; }
	void use(Stack A) { InUse = A; }

	synchronized void install(int i)
		{
		if (InUse!=null) 
			{
			if ((i>=0)&&(i<InUse.size()))
				{
				arrange a = (arrange)InUse.elementAt(i);
				index =i;
				gr = a.makeGr();
				setCoords(gr,false);
				}
			else return;
			}
		else {
			archive.GraphArrays tempg = archive.graphDispatch.grab(i);
			if (tempg==null) return;
			gr = tempg;
			index =i;
			setCoords(gr,false);
			}
		}

	int getHiddenRegion() { return hiddenRegion; }

	void setHiddenRegion(int i)
		{
		if ((i>=0)&&(i<gr.faceList.length)) hiddenRegion=i;
		}
	private int hiddenRegion;

	private void setHiddenRegionWithTemp()
		{
		if (gr.tempList.length==0) { setHiddenRegion(); return; }
		int i,j,pos=0,t,temp=-1;
		for (j=0;j<gr.tempList.length;j++)
			{
			i = gr.tempList[j];
						t = gr.faceList[i].length;
						if (t >temp)
								{
								pos = i;
								temp = t;
								}
						}
				hiddenRegion = pos;
				}

	void setHiddenRegion()
		{
		int i,pos=0,t,temp=-1;
		for (i=0;i<gr.faceList.length;i++)
			{
			t = gr.faceList[i].length;
			if (t >temp)
				{
				pos = i;
				temp = t;
				}
			}
		hiddenRegion = pos;
		}

boolean IsHidden(int vertex)
	{
	for (int i=0;i<gr.faceList[hiddenRegion].length;i++)
		if (vertex==gr.faceList[hiddenRegion][i])
		return true;
	return false;
	}

static java.util.Random t = new java.util.Random();
synchronized void setCoords(archive.GraphArrays gr,boolean overlay)
	{
	if (gr==null) return;
	int i,j;
	this.gr = gr;
	if (!overlay)
	   {
		setHiddenRegion();
	   coords = new double[gr.adjacent.length][2];
	   for (i=0;i<gr.adjacent.length;i++)
		{ coords[i][0]=t.nextDouble() -0.5;
			coords[i][1]=t.nextDouble() -0.5; }
	   }
	else if (coords.length != gr.adjacent.length)
	   {
	   double [][] bak = coords;
	   coords = new double[gr.adjacent.length][2];
	   int b = Math.min(bak.length,coords.length);
	   for (i=0;i<b;i++) { 
		  coords[i][0]=bak[i][0]; coords[i][1]=bak[i][1]; }
	   for (i=b;i<coords.length;i++)
		{ coords[i][0]=t.nextDouble() -0.5;
			coords[i][1]=t.nextDouble() -0.5; }
	   }
	if (hiddenRegion>=gr.faceList.length) hiddenRegion=0;
	int r[] = gr.faceList[hiddenRegion];
	double N = r.length;
	for (i=0;i<r.length;i++)
		{
		coords[r[i]][0]=
			Math.cos(6.28319*i/N);
		coords[r[i]][1]=
			Math.sin(6.28319*i/N);
		}
	}

synchronized boolean Move()
	{
	double stepsize = 0.02;
	int i;
	double backup[][] = new double [coords.length][2];
	for (i=0;i<coords.length;i++) 
		{ backup[i][0]=coords[i][0];
		  backup[i][1]=coords[i][1]; 
		}
	double a[][],b[][],c[][];
	a = new double [gr.faceList.length][];
	b = new double [gr.faceList.length][];
	c = new double [gr.faceList.length][];
	for (i=0;i<gr.faceList.length;i++)
		{
		a[i] = new double[gr.faceList[i].length];
		b[i] = new double[gr.faceList[i].length];
		c[i] = new double[gr.faceList[i].length];
		}
	int k,j0,j1,j2;
	int e[];
	int count =0;
	double total=0.0;
	for (k=0;k<gr.faceList.length;k++)
	for (i=0;i<gr.faceList[k].length;i++)
		{
		e = gr.faceList[k];
		j2 = e[misc.mod(i+1 , e.length)];
		j1 = e[i];
		j0 = e[misc.mod(i+e.length-1, e.length)];
					// Use determinant to equalize areas:
					c[k][i] = (coords[j0][0]-coords[j1][0])*
						(coords[j2][1]-coords[j1][1]) -
						(coords[j2][0]-coords[j1][0])*
						(coords[j0][1]-coords[j1][1]);
		b[k][i] = coords[j0][1]-coords[j2][1];
		a[k][i] = coords[j2][0]-coords[j0][0];
		count ++;
		if (IsHidden(j1)) 
			{
			count --;
			c[k][i]=b[k][i]=a[k][i]=0.0;
			}
		total += c[k][i];
		}
	double average = total/Math.max(1,count);
	double dx[]=new double[gr.adjacent.length],
		dy[]= new double[gr.adjacent.length];
	for (i=0;i<gr.adjacent.length;i++)
		{ dx[i]=dy[i]=0.0; }
	for (k=0;k<gr.faceList.length;k++)
	for (i=0;i<gr.faceList[k].length;i++)
		{
		j1 = gr.faceList[k][i];
		dx[j1] += (c[k][i]-average)*b[k][i];
		dy[j1] += (c[k][i]-average)*a[k][i];
		// July 97 correction:
		if (gr.vertexList[j1].length==2) {
			for (int h=0;h<gr.vertexList[j1].length;h++)
			{ int hface = gr.vertexList[j1][h];
			for (int h2=0;h2<gr.faceList[hface].length;h2++)
			  {
			  int h3 = gr.faceList[hface][h2];
			  dx[j1] += (coords[j1][0]-coords[h3][0]);
			  dy[j1] += (coords[j1][1]-coords[h3][1]);
			  }
			}}
		}
	for (i=0;i<gr.adjacent.length;i++)
		{
		backup[i][0] = coords[i][0] - stepsize*dx[i];
		backup[i][1] = coords[i][1] - stepsize*dy[i];
		}
	// check vertices on the hidden face .
	int r[] = gr.faceList[hiddenRegion];
	double t1,t2;
	double N = Math.max(1,r.length);
	for (i=0;i<r.length;i++)
		{
		t1 = Math.cos(6.28319*i/N);
		t2 = Math.sin(6.28319*i/N);
		backup[r[i]][0]=t1;
		backup[r[i]][1]=t2;
		}

	boolean hasMoved=false;
	for (i=0;i<backup.length;i++) 
			if ((Math.abs(backup[i][0]-coords[i][0])> 0.001)||
			    ((Math.abs(backup[i][0]-coords[i][0])> 0.001))) // was0.001
				{ hasMoved=true; break; }
	// if (!hasMoved) return false;
	for (i=0;i<backup.length;i++)
			{ coords[i][0]=backup[i][0]; coords[i][1]=backup[i][1]; }
	return true;
	}
}
