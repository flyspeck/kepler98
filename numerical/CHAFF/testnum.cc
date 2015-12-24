#include <iomanip.h>

class trace {
 
public:
    virtual double getY(int index)=0;         // 0..5
    virtual double getDih(int index)=0;       // 0..2
    virtual double getS()=0;
    virtual void setY(int index,double y)=0;
    virtual void setDih(int index,double d)=0;
    virtual void setS(double s)=0;
};

void numericallyAdjust(trace& t);

class traceD : public trace {
private:
    double y[6];
    double d[3];
    double s;
public:
    double getY(int i) { if ((i<0)||(i>5)) return 0; return y[i]; }
    double getDih(int i) { if ((i<0)||(i>2)) return 0; return d[i]; }
    double getS() { return s; }
    void setY(int i,double y0) { if ((i<0)||(i>5)) return; y[i]=y0; }
    void setDih(int i,double d0) { if ((i<0)||(i>2)) return; d[i]=d0; }
    void setS(double s0) { s = s0; }
    };

int main()
	{
	traceD t;
	int i;
	for (i=0;i<6;i++) t.setY(i,0);
	for (i=0;i<3;i++) t.setDih(i,0);
	t.setS(0);
	numericallyAdjust(t);
	for (i=0;i<6;i++) cout<< t.getY(i) << endl;
	for (i=0;i<3;i++) cout << t.getDih(i) << endl;
	cout << t.getS() << endl;
	}
	

