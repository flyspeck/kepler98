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


quad : numerical.h numerical/gradient.h numerical.h numerical/morefn.h quadX.cc quadX.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o quadX.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o quad 


quadX.o : quadX.cc
	$(CC) $(CCOPTS1) -c quadX.cc



pent : numerical.h numerical/gradient.h numerical.h numerical/morefn.h pentX.cc pentX.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o pentX.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o pent


pentX.o : pentX.cc
	$(CC) $(CCOPTS1) -c pentX.cc


bad2 : numerical.h numerical/gradient.h numerical.h numerical/morefn.h bad2X.cc bad2X.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o bad2X.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o bad2


bad2X.o : bad2X.cc
	$(CC) $(CCOPTS1) -c bad2X.cc


bad2a : numerical.h numerical/gradient.h numerical.h numerical/morefn.h bad2aX.cc bad2aX.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o bad2aX.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o bad2a


bad2aX.o : bad2aX.cc
	$(CC) $(CCOPTS1) -c bad2aX.cc

exc : numerical.h numerical/gradient.h numerical.h numerical/morefn.h excX.cc excX.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o excX.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o exc


excX.o : excX.cc
	$(CC) $(CCOPTS1) -c excX.cc







bad4 : numerical.h numerical/gradient.h numerical.h numerical/morefn.h bad4X.cc bad4X.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o bad4X.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o bad4


bad4X.o : bad4X.cc
	$(CC) $(CCOPTS1) -c bad4X.cc


bad4a : numerical.h numerical/gradient.h numerical.h numerical/morefn.h bad4aX.cc bad4aX.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o bad4aX.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o bad4a


bad4aX.o : bad4aX.cc
	$(CC) $(CCOPTS1) -c bad4aX.cc












bad5 : numerical.h numerical/gradient.h numerical.h numerical/morefn.h bad5X.cc bad5X.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o bad5X.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o bad5


bad5X.o : bad5X.cc
	$(CC) $(CCOPTS1) -c bad5X.cc



bad6 : numerical.h numerical/gradient.h numerical.h numerical/morefn.h bad6X.cc bad6X.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o bad6X.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o bad6


bad6X.o : bad6X.cc
	$(CC) $(CCOPTS1) -c bad6X.cc


ceq0 : numerical.h numerical/gradient.h numerical.h numerical/morefn.h ceq0X.cc ceq0X.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o ceq0X.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o ceq0


ceq0X.o : ceq0X.cc
	$(CC) $(CCOPTS1) -c ceq0X.cc

ceq1 : numerical.h numerical/gradient.h numerical.h numerical/morefn.h ceq1X.cc ceq1X.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o ceq1X.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o ceq1


ceq1X.o : ceq1X.cc
	$(CC) $(CCOPTS1) -c ceq1X.cc


demo : numerical.h numerical/gradient.h numerical.h numerical/morefn.h demoX.cc demoX.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o demoX.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o demo


demoX.o : demoX.cc
	$(CC) $(CCOPTS1) -c demoX.cc




ceq3 : numerical.h numerical/gradient.h numerical.h numerical/morefn.h ceq3X.cc ceq3X.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o ceq3X.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o ceq3


ceq3X.o : ceq3X.cc
	$(CC) $(CCOPTS1) -c ceq3X.cc




pt2 : numerical.h numerical/gradient.h numerical.h numerical/morefn.h pt2X.cc pt2X.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o pt2X.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o pt2


pt2X.o : pt2X.cc
	$(CC) $(CCOPTS1) -c pt2X.cc











qe : numerical.h numerical/gradient.h numerical.h numerical/morefn.h qeX.cc qeX.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o qeX.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o qe


qeX.o : qeX.cc
	$(CC) $(CCOPTS1) -c qeX.cc

cep0 : numerical.h numerical/gradient.h numerical.h numerical/morefn.h cep0X.cc cep0X.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o cep0X.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o cep0


cep0X.o : cep0X.cc
	$(CC) $(CCOPTS1) -c cep0X.cc



cep1 : numerical.h numerical/gradient.h numerical.h numerical/morefn.h cep1X.cc cep1X.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o cep1X.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o cep1


cep1X.o : cep1X.cc
	$(CC) $(CCOPTS1) -c cep1X.cc




cepA : numerical.h numerical/gradient.h numerical.h numerical/morefn.h cepAX.cc cepAX.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o cepAX.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o cepA


cepAX.o : cepAX.cc
	$(CC) $(CCOPTS1) -c cepAX.cc


cepB : numerical.h numerical/gradient.h numerical.h numerical/morefn.h cepBX.cc cepBX.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o cepBX.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o cepB


cepBX.o : cepBX.cc
	$(CC) $(CCOPTS1) -c cepBX.cc

cepC : numerical.h numerical/gradient.h numerical.h numerical/morefn.h cepCX.cc cepCX.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o cepCX.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o cepC


cepCX.o : cepCX.cc
	$(CC) $(CCOPTS1) -c cepCX.cc

cepD : numerical.h numerical/gradient.h numerical.h numerical/morefn.h cepDX.cc cepDX.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o cepDX.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o cepD


cepDX.o : cepDX.cc
	$(CC) $(CCOPTS1) -c cepDX.cc








tcc : numerical.h numerical/gradient.h numerical.h numerical/morefn.h tccX.cc tccX.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o tccX.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o tcc


tccX.o : tccX.cc
	$(CC) $(CCOPTS1) -c tccX.cc

cep14 : numerical.h numerical/gradient.h numerical.h numerical/morefn.h cep14X.cc cep14X.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o cep14X.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o cep14


cep14X.o : cep14X.cc
	$(CC) $(CCOPTS1) -c cep14X.cc

cep13 : numerical.h numerical/gradient.h numerical.h numerical/morefn.h cep13X.cc cep13X.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o cep13X.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o cep13


cep13X.o : cep13X.cc
	$(CC) $(CCOPTS1) -c cep13X.cc


cep12 : numerical.h numerical/gradient.h numerical.h numerical/morefn.h cep12X.cc cep12X.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o cep12X.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o cep12


cep12X.o : cep12X.cc
	$(CC) $(CCOPTS1) -c cep12X.cc



cep11 : numerical.h numerical/gradient.h numerical.h numerical/morefn.h cep11X.cc cep11X.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o cep11X.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o cep11


cep11X.o : cep11X.cc
	$(CC) $(CCOPTS1) -c cep11X.cc


cep9 : numerical.h numerical/gradient.h numerical.h numerical/morefn.h cep9X.cc cep9X.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o cep9X.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o cep9


cep9X.o : cep9X.cc
	$(CC) $(CCOPTS1) -c cep9X.cc


cep7 : numerical.h numerical/gradient.h numerical.h numerical/morefn.h cep7X.cc cep7X.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o cep7X.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o cep7


cep7X.o : cep7X.cc
	$(CC) $(CCOPTS1) -c cep7X.cc


cep6 : numerical.h numerical/gradient.h numerical.h numerical/morefn.h cep6X.cc cep6X.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o cep6X.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o cep6


cep6X.o : cep6X.cc
	$(CC) $(CCOPTS1) -c cep6X.cc



cep4 : numerical.h numerical/gradient.h numerical.h numerical/morefn.h cep4X.cc cep4X.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o cep4X.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o cep4


cep4X.o : cep4X.cc
	$(CC) $(CCOPTS1) -c cep4X.cc














partU : numerical.h numerical/gradient.h numerical.h numerical/morefn.h partU.cc partU.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o partU.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o partU 


partU.o : partU.cc
	$(CC) $(CCOPTS1) -c partU.cc


parts : numerical.h numerical/gradient.h numerical.h numerical/morefn.h partX.cc partX.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o partX.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o parts 



parts3 : numerical.h numerical/gradient.h numerical.h numerical/morefn.h partXa.cc partXa.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o partXa.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o parts3 


partXa.o : partXa.cc
	$(CC) $(CCOPTS1) -c partXa.cc


parts4 : numerical.h numerical/gradient.h numerical.h numerical/morefn.h partX4.cc partX4.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o partX4.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o parts4 


partX4.o : partX4.cc
	$(CC) $(CCOPTS1) -c partX4.cc



parts4a : numerical.h numerical/gradient.h numerical.h numerical/morefn.h partX4a.cc partX4a.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o partX4a.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o parts4a 


partX4a.o : partX4a.cc
	$(CC) $(CCOPTS1) -c partX4a.cc


parts4b : numerical.h numerical/gradient.h numerical.h numerical/morefn.h partX4b.cc partX4b.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o partX4b.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o parts4b


partX4b.o : partX4b.cc
	$(CC) $(CCOPTS1) -c partX4b.cc








parts5 : numerical.h numerical/gradient.h numerical.h numerical/morefn.h partX5.cc partX5.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o partX5.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o parts5 


partX5.o : partX5.cc
	$(CC) $(CCOPTS1) -c partX5.cc



parts6 : numerical.h numerical/gradient.h numerical.h numerical/morefn.h partX6.cc partX6.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o partX6.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o parts6


partX6.o : partX6.cc
	$(CC) $(CCOPTS1) -c partX6.cc




parts7 : numerical.h numerical/gradient.h numerical.h numerical/morefn.h partX7.cc partX7.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o partX7.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o parts7


partX7.o : partX7.cc
	$(CC) $(CCOPTS1) -c partX7.cc




parts8 : numerical.h numerical/gradient.h numerical.h numerical/morefn.h partX8.cc partX8.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o partX8.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o parts8


partX8.o : partX8.cc
	$(CC) $(CCOPTS1) -c partX8.cc





parts9 : numerical.h numerical/gradient.h numerical.h numerical/morefn.h partX9.cc partX9.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o partX9.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o parts9


partX9.o : partX9.cc
	$(CC) $(CCOPTS1) -c partX9.cc



parts10 : numerical.h numerical/gradient.h numerical.h numerical/morefn.h partX10.cc partX10.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o partX10.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o parts10


partX10.o : partX10.cc
	$(CC) $(CCOPTS1) -c partX10.cc



parts11 : numerical.h numerical/gradient.h numerical.h numerical/morefn.h partX11.cc partX11.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o partX11.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o parts11


partX11.o : partX11.cc
	$(CC) $(CCOPTS1) -c partX11.cc




parts12 : numerical.h numerical/gradient.h numerical.h numerical/morefn.h partX12.cc partX12.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o partX12.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o parts12


partX12.o : partX12.cc
	$(CC) $(CCOPTS1) -c partX12.cc

edge : numerical.h numerical/gradient.h numerical.h numerical/morefn.h edgeX.cc edgeX.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o edgeX.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o edge


edgeX.o : edgeX.cc
	$(CC) $(CCOPTS1) -c edgeX.cc



xst2 : numerical.h numerical/gradient.h numerical.h numerical/morefn.h xstX2.cc xstX2.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o xstX2.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o xst2


xstX2.o : xstX2.cc
	$(CC) $(CCOPTS1) -c xstX2.cc







parts2 : numerical.h numerical/gradient.h numerical.h numerical/morefn.h quadX.cc quadX.o
	$(CC) $(CCOPTS1) numerical/ineq.cc numerical/gradient.o numerical.o numerical/morefn.o quadX.o constrained_min.o cfsqp.o qld.o numerical/constants.o -lm -o parts2 

eqns : numerical.h gradient.h numerical.h morefn.h eqns.cc 
	$(CC) $(CCOPTS1) eqns.cc gradient.o numerical.o morefn.o constants.o -lm -o eqns
octeqns.o : numerical.h gradient.h numerical.h morefn.h octeqns.cc 
	$(CC) $(CCOPTS1) -c octeqns.cc 


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
