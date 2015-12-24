# Makefile for dimred_quad, boundary_quad, etc.
#
#
#	make           to produce executables for all programs
#
#	make clean     to remove object files and core
#
#       BINDIR -- directory where to copy executables
#       CC --     ANSI C compiler
#       CCOPTS1 - options used to compile .c files
#       CCOPTS2 -- options used to link .o files
#
#       CFILES -- .c files for printing
#	FILES -- all other files
#

#for Gnu's gcc compiler -O2 for optimization, -g for debugging, -Wall for check
#
CC     = gcc
CCOPTS1 = -O2 -ansi 
BINDIR  = ../sun

# for Sun's cc compiler, -fast or O2 for optimization, -g for debugging
#CC = cc
#CCOPTS1 = -Xc -v -fast

# for Silicon Graphics cc compiler, -O2 for optimization, -g for debugging
#CC = cc
#CCOPTS1 = -ansi -O2 -OPT:IEEE_comparisons:IEEE_arithmetic=1:roundoff=0
#BINDIR  = ../sgi
#BINDIR  = bin

# for loader, ld
CCOPTS2 = $(CCOPTS1)

# OBJS is in library order, CFILES in alphabetical order after qhull.c
OBJS = 	interval.o i_sphere.o i_bounds.o sphere.o second_partials.o \
				i_taylor.o i_appendix.o i_voronoi.o
CFILES= interval.c i_sphere.c i_bounds.c dimred_quad.c sphere.c \
        testmain.c \
        second_partials.c \
        i_taylor.c i_appendix.c i_voronoi.c
        
FILES=  Makefile 

all: 		TestInterval


test: 	TestInterval

sec9:				Sec9_1 Sec9_2 Sec9_3a Sec9_3b Sec9_4 Sec9_5

sec12:			Sec12_1 Sec12_2

sec14:			Sec14_1 Sec14_2 Sec14_3 Sec14_4 Sec14_5 \
						Sec14_6 Sec14_7 Sec14_8 Sec14_9

sec15:			Sec15

sec17:			Sec17

sec18:			Sec18

sec19:			Sec19_1 Sec19_2

sec20:			Sec20

sec21:			Sec21

sec22:			Sec22_1 Sec22_2

sec23:			Sec23_1 Sec23_2 Sec23_2a

sec271:			Sec271 Sec271a Sec272 Sec273 Sec274 Sec275 Sec276 Sec277 \
						Sec278 Sec279 Sec2710 Sec2711 Sec2712

sec281:			Sec281 Sec281a Sec281b

sec282:			Sec282 Sec282a

sec283:			Sec283

sec284:			Sec284 Sec284a Sec284b

sec285:			Sec285

sec286:			Sec286 Sec286a

sec38:			Sec38_1 Sec38_2 Sec38_2a Sec38_2b Sec38_3 Sec38_4

fake:				Fake1a Fake1b Fake2

octa:				Octa

sec410:			Sec410_2 Sec410_3

sec44:			Sec442_1 Sec442_1a Sec442_2 Sec442_4 Sec442_5 Sec441_1 Sec441_2

sec45:			Sec45_1 Sec45_1a Sec45_7 Sec451_1 Sec451_2 Sec454_1 Sec454_2

sec461:			Sec46_1a Sec46_1b Sec46_1c Sec46_1d Sec46_1e

sec462:			Sec46_2a Sec46_2b Sec46_2c Sec46_2d Sec46_2e Sec46_2f Sec46_2g

sec463:			Sec46_3 Sec46_3a

sec468:			Sec46_8 Sec46_8a Sec46_10

sec412:			Sec412_4 Sec412_4a Sec412_4b Sec412_5 Sec412_6 \
						Sec412_8 Sec412_8b Sec412_9 Sec412_9b

sec49:			Sec49_1 Sec49_2

# all

interval.o:			system_headers.h macros.h
i_sphere.o:			system_headers.h i_sphere.h interval.h i_bounds.h macros.h
i_bounds.o:			system_headers.h i_sphere.h interval.h i_bounds.h macros.h
second_partials.o:			system_headers.h i_sphere.h interval.h i_bounds.h \
		macros.h i_taylor.h second_partials.h
i_taylor.o:			system_headers.h i_sphere.h interval.h i_bounds.h macros.h \
		i_taylor.h
sphere.o:				system_headers.h sphere.h
i_bounds.o:			system_headers.h i_sphere.h interval.h i_bounds.h macros.h
i_appendix.o:		system_headers.h i_sphere.h interval.h i_bounds.h macros.h
i_voronoi.o:		system_headers.h i_sphere.h interval.h i_bounds.h macros.h \
		second_partials.h i_voronoi.h

# test

testmain.o:		system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_taylor.h second_partials.h

section_a9_1.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h

section_a9_2.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h

section_a9_3a.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h

section_a9_3b.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h

section_a9_4.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h

section_a9_5.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h

fake_anchor_1a.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h

fake_anchor_1b.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h

fake_anchor_2.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h

octahedra.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h i_taylor.h

section_a15_1.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h

section_a17_1.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h

section_a18_1.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h

section_a19_1.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h

section_a19_2.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h

section_a20_1.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h

section_a21_1.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h

section_a22_1.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h

section_a22_2.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h

section_a23_1.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h

section_a23_2.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h

section_a23_2a.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h

section_a271.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a271a.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a272.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a273.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a274.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a275.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a276.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a277.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a278.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a279.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a2710.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a2711.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a2712.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a281.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a281a.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a281b.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a12_1.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a12_2.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a282.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a282a.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a283.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a284.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a284a.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a284b.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a285.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a286.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a286a.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a38_1.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a38_2.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a38_2a.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a38_2b.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a38_3.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a38_4.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a410_2.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a410_3.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a442_1.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a442_1a.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a442_2.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a442_4.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a442_5.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a441_1.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a441_2.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a45_1.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a45_1a.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a45_7.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a451_1.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a451_2.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a454_1.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a454_2.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a46_1a.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a46_1b.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a46_1c.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a46_1d.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a46_1e.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a46_2a.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a46_2b.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a46_2c.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a46_2d.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a46_2e.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a46_2f.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a46_2g.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a46_3.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a46_3a.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a46_8.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a46_8a.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a46_10.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a412_4.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a412_4a.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a412_4b.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a412_5.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a412_6.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a412_8.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a412_8b.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a412_9.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a412_9b.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a49_1.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

section_a49_2.o:	system_headers.h sphere.h interval.h i_sphere.h \
		i_bounds.h macros.h i_appendix.h second_partials.h i_voronoi.h \
		i_taylor.h

.c.o:
	$(CC) -c $(CCOPTS1) $<

clean:
	rm -f *.o core

move:
	mv Test* Sec* Fake* $(BINDIR); \
	chmod 744 do*; mv do* $(BINDIR)

# test

TestInterval:  testmain.o $(OBJS)
	$(CC) -o TestInterval $(CCOPTS2) testmain.o $(OBJS) -lm 

# sec9

Sec9_1:  section_a9_1.o $(OBJS)
	$(CC) -o Sec9_1 $(CCOPTS2) section_a9_1.o $(OBJS) -lm 

Sec9_2:  section_a9_2.o $(OBJS)
	$(CC) -o Sec9_2 $(CCOPTS2) section_a9_2.o $(OBJS) -lm 

Sec9_3a:  section_a9_3a.o $(OBJS)
	$(CC) -o Sec9_3a $(CCOPTS2) section_a9_3a.o $(OBJS) -lm 

Sec9_3b:  section_a9_3b.o $(OBJS)
	$(CC) -o Sec9_3b $(CCOPTS2) section_a9_3b.o $(OBJS) -lm 

Sec9_4:  section_a9_4.o $(OBJS)
	$(CC) -o Sec9_4 $(CCOPTS2) section_a9_4.o $(OBJS) -lm 

Sec9_5:  section_a9_5.o $(OBJS)
	$(CC) -o Sec9_5 $(CCOPTS2) section_a9_5.o $(OBJS) -lm 

# sec12

Sec12_1:  section_a12_1.o $(OBJS)
	$(CC) -o Sec12_1 $(CCOPTS2) section_a12_1.o $(OBJS) -lm 

Sec12_2:  section_a12_2.o $(OBJS)
	$(CC) -o Sec12_2 $(CCOPTS2) section_a12_2.o $(OBJS) -lm 

# sec15

Sec15:  section_a15_1.o $(OBJS)
	$(CC) -o Sec15 $(CCOPTS2) section_a15_1.o $(OBJS) -lm 

# sec17

Sec17:  section_a17_1.o $(OBJS)
	$(CC) -o Sec17 $(CCOPTS2) section_a17_1.o $(OBJS) -lm 

# sec18

Sec18:  section_a18_1.o $(OBJS)
	$(CC) -o Sec18 $(CCOPTS2) section_a18_1.o $(OBJS) -lm 

# sec19

Sec19_1:  section_a19_1.o $(OBJS)
	$(CC) -o Sec19_1 $(CCOPTS2) section_a19_1.o $(OBJS) -lm 

Sec19_2:  section_a19_2.o $(OBJS)
	$(CC) -o Sec19_2 $(CCOPTS2) section_a19_2.o $(OBJS) -lm 

# sec20

Sec20:  section_a20_1.o $(OBJS)
	$(CC) -o Sec20 $(CCOPTS2) section_a20_1.o $(OBJS) -lm 

# sec21

Sec21:  section_a21_1.o $(OBJS)
	$(CC) -o Sec21 $(CCOPTS2) section_a21_1.o $(OBJS) -lm 

# sec22

Sec22_1:  section_a22_1.o $(OBJS)
	$(CC) -o Sec22_1 $(CCOPTS2) section_a22_1.o $(OBJS) -lm 

Sec22_2:  section_a22_2.o $(OBJS)
	$(CC) -o Sec22_2 $(CCOPTS2) section_a22_2.o $(OBJS) -lm 

# sec23

Sec23_1:  section_a23_1.o $(OBJS)
	$(CC) -o Sec23_1 $(CCOPTS2) section_a23_1.o $(OBJS) -lm 

Sec23_2:  section_a23_2.o $(OBJS)
	$(CC) -o Sec23_2 $(CCOPTS2) section_a23_2.o $(OBJS) -lm 

Sec23_2a:  section_a23_2a.o $(OBJS)
	$(CC) -o Sec23_2a $(CCOPTS2) section_a23_2a.o $(OBJS) -lm 

# sec271

Sec271:  section_a271.o $(OBJS)
	$(CC) -o Sec271 $(CCOPTS2) section_a271.o $(OBJS) -lm 

Sec271a:  section_a271a.o $(OBJS)
	$(CC) -o Sec271a $(CCOPTS2) section_a271a.o $(OBJS) -lm 

Sec272:  section_a272.o $(OBJS)
	$(CC) -o Sec272 $(CCOPTS2) section_a272.o $(OBJS) -lm 

Sec273:  section_a273.o $(OBJS)
	$(CC) -o Sec273 $(CCOPTS2) section_a273.o $(OBJS) -lm 

Sec274:  section_a274.o $(OBJS)
	$(CC) -o Sec274 $(CCOPTS2) section_a274.o $(OBJS) -lm 

Sec275:  section_a275.o $(OBJS)
	$(CC) -o Sec275 $(CCOPTS2) section_a275.o $(OBJS) -lm 

Sec276:  section_a276.o $(OBJS)
	$(CC) -o Sec276 $(CCOPTS2) section_a276.o $(OBJS) -lm 

Sec277:  section_a277.o $(OBJS)
	$(CC) -o Sec277 $(CCOPTS2) section_a277.o $(OBJS) -lm 

Sec278:  section_a278.o $(OBJS)
	$(CC) -o Sec278 $(CCOPTS2) section_a278.o $(OBJS) -lm 

Sec279:  section_a279.o $(OBJS)
	$(CC) -o Sec279 $(CCOPTS2) section_a279.o $(OBJS) -lm 

Sec2710:  section_a2710.o $(OBJS)
	$(CC) -o Sec2710 $(CCOPTS2) section_a2710.o $(OBJS) -lm 

Sec2711:  section_a2711.o $(OBJS)
	$(CC) -o Sec2711 $(CCOPTS2) section_a2711.o $(OBJS) -lm 

Sec2712:  section_a2712.o $(OBJS)
	$(CC) -o Sec2712 $(CCOPTS2) section_a2712.o $(OBJS) -lm 

# sec281

Sec281:  section_a281.o $(OBJS)
	$(CC) -o Sec281 $(CCOPTS2) section_a281.o $(OBJS) -lm 

Sec281a:  section_a281a.o $(OBJS)
	$(CC) -o Sec281a $(CCOPTS2) section_a281a.o $(OBJS) -lm 

Sec281b:  section_a281b.o $(OBJS)
	$(CC) -o Sec281b $(CCOPTS2) section_a281b.o $(OBJS) -lm 

# sec282

Sec282:  section_a282.o $(OBJS)
	$(CC) -o Sec282 $(CCOPTS2) section_a282.o $(OBJS) -lm 

Sec282a:  section_a282a.o $(OBJS)
	$(CC) -o Sec282a $(CCOPTS2) section_a282a.o $(OBJS) -lm 

# sec283

Sec283:  section_a283.o $(OBJS)
	$(CC) -o Sec283 $(CCOPTS2) section_a283.o $(OBJS) -lm 

# sec284

Sec284:  section_a284.o $(OBJS)
	$(CC) -o Sec284 $(CCOPTS2) section_a284.o $(OBJS) -lm 

Sec284a:  section_a284a.o $(OBJS)
	$(CC) -o Sec284a $(CCOPTS2) section_a284a.o $(OBJS) -lm 

Sec284b:  section_a284b.o $(OBJS)
	$(CC) -o Sec284b $(CCOPTS2) section_a284b.o $(OBJS) -lm 

# sec285

Sec285:  section_a285.o $(OBJS)
	$(CC) -o Sec285 $(CCOPTS2) section_a285.o $(OBJS) -lm 

# sec286

Sec286:  section_a286.o $(OBJS)
	$(CC) -o Sec286 $(CCOPTS2) section_a286.o $(OBJS) -lm 

Sec286a:  section_a286a.o $(OBJS)
	$(CC) -o Sec286a $(CCOPTS2) section_a286a.o $(OBJS) -lm 

# sec38

Sec38_1:  section_a38_1.o $(OBJS)
	$(CC) -o Sec38_1 $(CCOPTS2) section_a38_1.o $(OBJS) -lm 

Sec38_2:  section_a38_2.o $(OBJS)
	$(CC) -o Sec38_2 $(CCOPTS2) section_a38_2.o $(OBJS) -lm 

Sec38_2a:  section_a38_2a.o $(OBJS)
	$(CC) -o Sec38_2a $(CCOPTS2) section_a38_2a.o $(OBJS) -lm 

Sec38_2b:  section_a38_2b.o $(OBJS)
	$(CC) -o Sec38_2b $(CCOPTS2) section_a38_2b.o $(OBJS) -lm 

Sec38_3:  section_a38_3.o $(OBJS)
	$(CC) -o Sec38_3 $(CCOPTS2) section_a38_3.o $(OBJS) -lm 

Sec38_4:  section_a38_4.o $(OBJS)
	$(CC) -o Sec38_4 $(CCOPTS2) section_a38_4.o $(OBJS) -lm 

# fake

Fake1a:  fake_anchor_1a.o $(OBJS)
	$(CC) -o Fake1a $(CCOPTS2) fake_anchor_1a.o $(OBJS) -lm 

Fake1b:  fake_anchor_1b.o $(OBJS)
	$(CC) -o Fake1b $(CCOPTS2) fake_anchor_1b.o $(OBJS) -lm 

Fake2:  fake_anchor_2.o $(OBJS)
	$(CC) -o Fake2 $(CCOPTS2) fake_anchor_2.o $(OBJS) -lm 

# octa

Octa:  octahedra.o $(OBJS)
	$(CC) -o Octa $(CCOPTS2) octahedra.o $(OBJS) -lm 

# sec410

Sec410_2:  section_a410_2.o $(OBJS)
	$(CC) -o Sec410_2 $(CCOPTS2) section_a410_2.o $(OBJS) -lm 

Sec410_3:  section_a410_3.o $(OBJS)
	$(CC) -o Sec410_3 $(CCOPTS2) section_a410_3.o $(OBJS) -lm 

# sec44

Sec442_1:  section_a442_1.o $(OBJS)
	$(CC) -o Sec442_1 $(CCOPTS2) section_a442_1.o $(OBJS) -lm 

Sec442_1a:  section_a442_1a.o $(OBJS)
	$(CC) -o Sec442_1a $(CCOPTS2) section_a442_1a.o $(OBJS) -lm 

Sec442_2:  section_a442_2.o $(OBJS)
	$(CC) -o Sec442_2 $(CCOPTS2) section_a442_2.o $(OBJS) -lm 

Sec442_4:  section_a442_4.o $(OBJS)
	$(CC) -o Sec442_4 $(CCOPTS2) section_a442_4.o $(OBJS) -lm 

Sec442_5:  section_a442_5.o $(OBJS)
	$(CC) -o Sec442_5 $(CCOPTS2) section_a442_5.o $(OBJS) -lm 

Sec441_1:  section_a441_1.o $(OBJS)
	$(CC) -o Sec441_1 $(CCOPTS2) section_a441_1.o $(OBJS) -lm 

Sec441_2:  section_a441_2.o $(OBJS)
	$(CC) -o Sec441_2 $(CCOPTS2) section_a441_2.o $(OBJS) -lm 

# sec45

Sec45_1:  section_a45_1.o $(OBJS)
	$(CC) -o Sec45_1 $(CCOPTS2) section_a45_1.o $(OBJS) -lm 

Sec45_1a:  section_a45_1a.o $(OBJS)
	$(CC) -o Sec45_1a $(CCOPTS2) section_a45_1a.o $(OBJS) -lm 

Sec45_7:  section_a45_7.o $(OBJS)
	$(CC) -o Sec45_7 $(CCOPTS2) section_a45_7.o $(OBJS) -lm 

Sec451_1:  section_a451_1.o $(OBJS)
	$(CC) -o Sec451_1 $(CCOPTS2) section_a451_1.o $(OBJS) -lm 

Sec451_2:  section_a451_2.o $(OBJS)
	$(CC) -o Sec451_2 $(CCOPTS2) section_a451_2.o $(OBJS) -lm 

Sec454_1:  section_a454_1.o $(OBJS)
	$(CC) -o Sec454_1 $(CCOPTS2) section_a454_1.o $(OBJS) -lm 

Sec454_2:  section_a454_2.o $(OBJS)
	$(CC) -o Sec454_2 $(CCOPTS2) section_a454_2.o $(OBJS) -lm 

# sec461

Sec46_1a:  section_a46_1a.o $(OBJS)
	$(CC) -o Sec46_1a $(CCOPTS2) section_a46_1a.o $(OBJS) -lm 

Sec46_1b:  section_a46_1b.o $(OBJS)
	$(CC) -o Sec46_1b $(CCOPTS2) section_a46_1b.o $(OBJS) -lm 

Sec46_1c:  section_a46_1c.o $(OBJS)
	$(CC) -o Sec46_1c $(CCOPTS2) section_a46_1c.o $(OBJS) -lm 

Sec46_1d:  section_a46_1d.o $(OBJS)
	$(CC) -o Sec46_1d $(CCOPTS2) section_a46_1d.o $(OBJS) -lm 

Sec46_1e:  section_a46_1e.o $(OBJS)
	$(CC) -o Sec46_1e $(CCOPTS2) section_a46_1e.o $(OBJS) -lm 

# sec462

Sec46_2a:  section_a46_2a.o $(OBJS)
	$(CC) -o Sec46_2a $(CCOPTS2) section_a46_2a.o $(OBJS) -lm 

Sec46_2b:  section_a46_2b.o $(OBJS)
	$(CC) -o Sec46_2b $(CCOPTS2) section_a46_2b.o $(OBJS) -lm 

Sec46_2c:  section_a46_2c.o $(OBJS)
	$(CC) -o Sec46_2c $(CCOPTS2) section_a46_2c.o $(OBJS) -lm 

Sec46_2d:  section_a46_2d.o $(OBJS)
	$(CC) -o Sec46_2d $(CCOPTS2) section_a46_2d.o $(OBJS) -lm 

Sec46_2e:  section_a46_2e.o $(OBJS)
	$(CC) -o Sec46_2e $(CCOPTS2) section_a46_2e.o $(OBJS) -lm 

Sec46_2f:  section_a46_2f.o $(OBJS)
	$(CC) -o Sec46_2f $(CCOPTS2) section_a46_2f.o $(OBJS) -lm 

Sec46_2g:  section_a46_2g.o $(OBJS)
	$(CC) -o Sec46_2g $(CCOPTS2) section_a46_2g.o $(OBJS) -lm 

# sec463

Sec46_3:  section_a46_3.o $(OBJS)
	$(CC) -o Sec46_3 $(CCOPTS2) section_a46_3.o $(OBJS) -lm 

Sec46_3a:  section_a46_3a.o $(OBJS)
	$(CC) -o Sec46_3a $(CCOPTS2) section_a46_3a.o $(OBJS) -lm 

# sec468

Sec46_8:  section_a46_8.o $(OBJS)
	$(CC) -o Sec46_8 $(CCOPTS2) section_a46_8.o $(OBJS) -lm 

Sec46_8a:  section_a46_8a.o $(OBJS)
	$(CC) -o Sec46_8a $(CCOPTS2) section_a46_8a.o $(OBJS) -lm 

Sec46_10:  section_a46_10.o $(OBJS)
	$(CC) -o Sec46_10 $(CCOPTS2) section_a46_10.o $(OBJS) -lm 

# sec412

Sec412_4:  section_a412_4.o $(OBJS)
	$(CC) -o Sec412_4 $(CCOPTS2) section_a412_4.o $(OBJS) -lm 

Sec412_4a:  section_a412_4a.o $(OBJS)
	$(CC) -o Sec412_4a $(CCOPTS2) section_a412_4a.o $(OBJS) -lm 

Sec412_4b:  section_a412_4b.o $(OBJS)
	$(CC) -o Sec412_4b $(CCOPTS2) section_a412_4b.o $(OBJS) -lm 

Sec412_5:  section_a412_5.o $(OBJS)
	$(CC) -o Sec412_5 $(CCOPTS2) section_a412_5.o $(OBJS) -lm 

Sec412_6:  section_a412_6.o $(OBJS)
	$(CC) -o Sec412_6 $(CCOPTS2) section_a412_6.o $(OBJS) -lm 

Sec412_8:  section_a412_8.o $(OBJS)
	$(CC) -o Sec412_8 $(CCOPTS2) section_a412_8.o $(OBJS) -lm 

Sec412_8b:  section_a412_8b.o $(OBJS)
	$(CC) -o Sec412_8b $(CCOPTS2) section_a412_8b.o $(OBJS) -lm 

Sec412_9:  section_a412_9.o $(OBJS)
	$(CC) -o Sec412_9 $(CCOPTS2) section_a412_9.o $(OBJS) -lm 

Sec412_9b:  section_a412_9b.o $(OBJS)
	$(CC) -o Sec412_9b $(CCOPTS2) section_a412_9b.o $(OBJS) -lm 

# sec49

Sec49_1:  section_a49_1.o $(OBJS)
	$(CC) -o Sec49_1 $(CCOPTS2) section_a49_1.o $(OBJS) -lm 

Sec49_2:  section_a49_2.o $(OBJS)
	$(CC) -o Sec49_2 $(CCOPTS2) section_a49_2.o $(OBJS) -lm 


# end of Makefile
