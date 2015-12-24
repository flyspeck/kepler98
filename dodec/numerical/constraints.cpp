#include<iostream.h>
#include<fstream.h>

int main()
{
  int i;
  ofstream outFile;
  outFile.open("case.dat");

  for(i=53;i<108;i++)
    {
      outFile<<"case "<<i+1<<" : ";
      
    }

  return 0;
}
	
