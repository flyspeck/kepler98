cplexInit.o : cplexInit.cc
	CC -I/usr/local/include -dalign -c cplexInit.cc

cplexVar : test.cc cplexInit.o loadFaceData.o
	CC -dalign test.cc cplexInit.o -I/usr/local/include -L/usr/local/cplex \
	/usr/local/cplex/libcplex.a loadFaceData.o \
	gradient.o numericallyAdjust.o numerical.o -o cplexVar -lm

loadFaceData.o : cplexInit.o loadFaceData.cc
	CC -dalign -c loadFaceData.cc  -I/usr/local/include -L/usr/local/cplex 
