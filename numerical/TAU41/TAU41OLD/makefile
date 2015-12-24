# for Silicon graphics cc compiler
CC = CC
CCOPTS1 =  -g -dalign

## for g++ on SUNS
#CC = g++
#CCOPTS1 = -g -O2


numerical.o : numerical.h numerical.cc
	$(CC) $(CCOPTS1) -c  numerical.cc
gradient.o : numerical.h gradient.cc
	$(CC) $(CCOPTS1) -c  gradient.cc   
morefn.o : morefn.h numerical.h gradient.cc
	$(CC) $(CCOPTS1) -c  morefn.cc   
bv.o : morefn.h numerical.h gradient.cc
	$(CC) $(CCOPTS1) -c  bv.cc   
constants.o : constants.h constants.cc
	$(CC) $(CCOPTS1) -c constants.cc
bv : morefn.h numerical.h gradient.cc
	$(CC) $(CCOPTS1)   bv.cc   -lm

octineq : numerical.h gradient.h numerical.h constants.h octeqns.h octineq.cc 
	$(CC) $(CCOPTS1) octineq.cc gradient.o numerical.o constants.o morefn.o octeqns.o constrained_min.o cfsqp.o qld.o -I/usr/local/include -L/usr/local/cplex \
    /usr/local/cplex/libcplex.a -lm -o octineq
gamma : numerical.h gradient.h numerical.h constants.h octeqns.h gamma.cc octeqns.cc
	$(CC) $(CCOPTS1) gamma.cc gradient.o numerical.o constants.o morefn.o octeqns.o octineq.o constrained_min.o cfsqp.o qld.o -I/usr/local/include -L/usr/local/cplex \
    /usr/local/cplex/libcplex.a -lm -o gamma

octineq.o : numerical.h gradient.h numerical.h constants.h octeqns.h octineq.cc 
	$(CC) $(CCOPTS1) -c octineq.cc -I/usr/local/include 

eqns : numerical.h gradient.h numerical.h morefn.h eqns.cc 
	$(CC) $(CCOPTS1) eqns.cc gradient.o numerical.o morefn.o constants.o -lm -o eqns
octeqns : numerical.h gradient.h numerical.h morefn.h octeqns.cc constants.h
	$(CC) $(CCOPTS1) octeqns.cc gradient.o numerical.o morefn.o constants.o \
	-I/usr/local/include -L/usr/local/cplex \
    /usr/local/cplex/libcplex.a -o octeqns
octeqns.o : numerical.h gradient.h numerical.h morefn.h octeqns.cc constants.h
	$(CC) $(CCOPTS1) -c octeqns.cc -I/usr/local/include 


partX.o : partX.cc
	$(CC) $(CCOPTS1) -c partX.cc


part4bis.o : part4bis.cc
	$(CC) $(CCOPTS1) -c part4bis.cc

ineq : numerical.h numerical.o morefn.o constants.o morefn.h constants.h part4.cc\
	gradient.o gradient.cc gradient.h ineq.cc constants.cc constrained_min.o \
	part4.o
	$(CC) $(CCOPTS1) ineq.cc gradient.o numerical.o morefn.o part4.o\
	constrained_min.o cfsqp.o qld.o \
	constants.o -lm -o ineq

numericallyAdjust.o : numerical.h morefn.h constants.h \
	gradient.cc gradient.h numericallyAdjust.cc constants.cc
	$(CC) $(CCOPTS1) -c numericallyAdjust.cc 
testnum : numericallyAdjust.o testnum.cc
	$(CC) $(CCOPTS) testnum.cc gradient.o numerical.o  \
    numericallyAdjust.o -lm -o testnum


checkineq : numerical.h numerical.o morefn.o morefn.h \
	gradient.o gradient.cc gradient.h checkineq.cc
	$(CC) $(CCOPTS1) checkineq.cc gradient.o numerical.o morefn.o -lm -o checkineq

sample: numerical.h 
	$(CC) $(CCOPTS1) sample.c gradient.o numerical.o morefn.o -lm -o sample
front : numerical.h numerical.o gradient.o gradient.cc gradient.h front.cc
	$(CC) $(CCOPTS1) front.cc gradient.o numerical.o -lm -o front

lpex1: lpex1.c octeqns.o numerical.h constants.h morefn.h 
	$(CC) $(CCOPTS1) lpex1.c -dalign -I/usr/local/include -L/usr/local/cplex \
	/usr/local/cplex/libcplex.a octeqns.o numerical.o constants.o morefn.o -o lpex1 -lm
