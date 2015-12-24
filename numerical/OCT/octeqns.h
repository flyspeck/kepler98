


// break an octahedron into four upright quarters.

/* uprights are a,b,c,d,

	a = (0,1,2,3,4,5),
	b = (0,2,7,6,8,4),
	c = (0,7,11,9,10,8),
	d = (0,11,1,12,5,10),

	type[0] = face(0,1,5), between a and d
	type[1] = face(0,2,4), between a and b
	type[2] = face(0,7,8), between b and c
	type[3] = face(0,11,10), between c and d.

	type = 0 means face < sqrt(2),
	type = 1 means face > sqrt(2).

	ymin[13] labels 0..12 as above,
	ymax[13] as above.

	c[29]= outdata[29]... returned by separate, an array of doubles
		giving the coefficients of the inqualities.

	c[0] + c[1] y[0] + c[2] y[1] + c[3] y[2] + c[4] y[4] + c[5] y[5] +
				c[6] (dih(a)-pi/2) - c[28] < - f(a).
	c[7] + c[8] y[0] + c[9] y[2] + c[10] y[7] + c[11] y[8] + c[12] y[4] +
				c[13] (dih(b)-pi/2) - c[28] < - f(b).
	c[14] + c[15] y[0] + c[16] y[7] + c[17] y[11] + c[18] y[10] + c[19] y[8] +
				c[20] (dih(c)-pi/2) - c[28] < - f(c).
	c[21] + c[22] y[0] + c[23] y[11] + c[24] y[1] + c[25] y[5] + c[26] y[10] +
				c[27] (dih(d)-pi/2) - c[28] < - f(d).

	The outdata should always satisify.
	c[0]+c[7]+c[14]+c[14]>0,
	c[1]+c[8]+c[15]+c[22]>0,
	c[6]+c[13]+c[20]+c[27]=0,
	c[2]+c[24]>0, etc. (take coefficients of the same y[i], there are 8 such
						equations)

	fA[2] gives two doubles,
	For example if sigma < -9.494 + 3.0508 dih,
	fA[2] = {-3.0508,9.494};

	fa+fb+fc+fd = sigma + fA[0] d + fA[1],
	so if, for example, fa<0, fb<0, fc<0, fd<0, then the inequality holds.


*/

int separate(double ymin[13],double ymax[13],double outdata[29],
    double fA[2],int type[4]);

/* this version calls cfsqp first, it is much more reliable */
void separateChecked(double ymin[13],double ymax[13],double outdata[29],
    double fA[2],int type[4]);

// score of upright simplex a.
double fa(double y1,double y2,double y3,double y4,double y5,double y6,
    double fA[2],int type[4]);
double fb(double y1,double y2,double y3,double y4,double y5,double y6,
    double fA[2],int type[4]);
double fc(double y1,double y2,double y3,double y4,double y5,double y6,
    double fA[2],int type[4]);
double fd(double y1,double y2,double y3,double y4,double y5,double y6,
    double fA[2],int type[4]);
// dih-pi/2 : 
double dpi(double y1,double y2,double y3,double y4,double y5,double y6);


double edge(double y0,double z0); // eta(return,y0,z0)==sqrt(2).
