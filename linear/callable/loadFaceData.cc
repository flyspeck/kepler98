#include <string.h>
#include <fstream.h>
#include <iostream.h>
#include <iomanip.h>
#include <ctype.h>
#include <stdlib.h>
#include "cplexInit.h"

static void readThorpe(char* filename,char* buf,int sz)
	{
	ifstream reading(filename);
	while (!reading.eof())
	{
	reading.getline(buf,sz,'\n');
	if (buf[0]=='#') cout << buf << endl;
	}
	reading.close();
	}

static int getNext(char* buf,int& start)
	{
	while (isspace(buf[start])&&(buf[start]!=0)) start++;
	if (buf[start]==0) { cout << "no integer found " << endl; return 0; }
	char intbuf[10];
	int i=0;
	while ((!isspace(buf[start+i]))&&(buf[start+i]!=0)&&(i<10)) 
		{
		intbuf[i] = buf[i+start]; i++;
		}
	if (i>=10) return 0;
	intbuf[i]=0;
	start = start+i;
	return atoi(intbuf);
	}

static void setRow(char* buf,int& start,int outbuf[],int& r)
	{
	r = getNext(buf,start);
	for (int i=0;i<r;i++) outbuf[i]= getNext(buf,start);
	}

static void setMatrix(char* buf,int& start,int outbuf[50][10],int len[50],
	int& numFace)
	{
	numFace = getNext(buf,start);
	for (int i=0;i<numFace;i++) setRow(buf,start,outbuf[i],len[i]);
	}

class faceDataInstance : public faceData {
public:
	int M[50][10];
	int len[50];
	int numFace;
	int getVertexCount() const {
		int count = 0;
		for (int i=0;i<numFace;i++)
		for (int j=0;j<len[i];j++)
			if (count<M[i][j]) count = M[i][j];
		return count;
		}
	int getFaceCount() const { return numFace; }
	int getVertexCount(int face) const { return len[face]; }
	int getVertexAt(int face,int pos) const {
		if ((face<0)||(face>numFace)||(pos<0)||(pos>=len[face])) return 0;
		return M[face][pos]-1;
		}
	};

faceDataInstance fD;
	
	

faceData& loadFaceData(char* filename)
	{
	char buf[300];
	readThorpe(filename,buf,250);
	int start=1;
	setMatrix(buf,start,fD.M,fD.len,fD.numFace);
	return fD;
	}
