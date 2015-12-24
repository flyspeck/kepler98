# for Silicon graphics cc compiler
CC = CC
CCOPTS1 = -O2 -g -dalign

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
quoinfn.o : quoinfn.cc quoinfn.h
	$(CC) $(CCOPTS1) -c  quoinfn.cc   -lm
testquoin : morefn.h quoinfn.h quoinfn.cc quoinfn.o testquoin.cc numerical.h
	$(CC) $(CCOPTS1)  numerical.o quoinfn.o testquoin.cc  -lm -o testquoin

parts : numerical.h gradient.h numerical.h morefn.h partX.cc partX.o quoinfn.h
	$(CC) $(CCOPTS1) ineq.cc quoinfn.o gradient.o numerical.o morefn.o partX.o constrained_min.o cfsqp.o qld.o constants.o -lm -o parts 

partU : numerical.h gradient.h numerical.h morefn.h partU.cc partU.o quoinfn.h
	$(CC) $(CCOPTS1) ineq.cc quoinfn.o gradient.o numerical.o morefn.o partU.o constrained_min.o cfsqp.o qld.o constants.o -lm -o partU

part4may98 : numerical.h gradient.h numerical.h morefn.h part4may98.cc part4may98.o quoinfn.h
	$(CC) $(CCOPTS1) ineq.cc quoinfn.o gradient.o numerical.o morefn.o part4may98.o constrained_min.o cfsqp.o qld.o constants.o -lm -o part4may98

partY : numerical.h gradient.h numerical.h morefn.h partY.cc partY.o quoinfn.h
	$(CC) $(CCOPTS1) ineq.cc quoinfn.o gradient.o numerical.o morefn.o partY.o constrained_min.o cfsqp.o qld.o constants.o -lm -o partY 

part4C : numerical.h gradient.h numerical.h morefn.h part4C.cc part4C.o
	$(CC) $(CCOPTS1) ineq.cc gradient.o numerical.o morefn.o part4C.o constrained_min.o cfsqp.o qld.o constants.o -lm -o part4C

partK : numerical.h gradient.h numerical.h morefn.h partK.cc partK.o
	$(CC) $(CCOPTS1) ineq.cc gradient.o numerical.o morefn.o partK.o constrained_min.o cfsqp.o qld.o constants.o -lm -o partK
partM : numerical.h gradient.h numerical.h morefn.h partM.cc partM.o
	$(CC) $(CCOPTS1) ineq.cc gradient.o numerical.o morefn.o partM.o constrained_min.o cfsqp.o qld.o constants.o -lm -o partM
partS : numerical.h gradient.h numerical.h morefn.h partS.cc partS.o
	$(CC) $(CCOPTS1) ineq.cc gradient.o numerical.o morefn.o partS.o constrained_min.o cfsqp.o qld.o constants.o -lm -o partS

eqns : numerical.h gradient.h numerical.h morefn.h eqns.cc 
	$(CC) $(CCOPTS1) eqns.cc gradient.o numerical.o morefn.o constants.o -lm -o eqns
octeqns.o : numerical.h gradient.h numerical.h morefn.h octeqns.cc 
	$(CC) $(CCOPTS1) -c octeqns.cc 


partX.o : partX.cc
	$(CC) $(CCOPTS1) -c partX.cc
partU.o : partU.cc
	$(CC) $(CCOPTS1) -c partU.cc
partY.o : partY.cc
	$(CC) $(CCOPTS1) -c partY.cc
part4C.o : part4C.cc
	$(CC) $(CCOPTS1) -c part4C.cc
partK.o : partK.cc
	$(CC) $(CCOPTS1) -c partK.cc
partM.o : partM.cc
	$(CC) $(CCOPTS1) -c partM.cc
partS.o : partS.cc
	$(CC) $(CCOPTS1) -c partS.cc
part4may98.o : part4may98.cc
	$(CC) $(CCOPTS1) -c part4may98.cc


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
