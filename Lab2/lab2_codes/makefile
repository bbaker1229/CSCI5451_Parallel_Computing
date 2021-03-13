.SUFFIXES: .c .o .f .F

##CC			=  gcc
CC			=  mpicc 
CFLAGS			= -g -Wall
##LIBS			=  -L/project/scicom/scicom00/SOFT/lib/linux/ -lblas
LIBS			=  -lblas 

FILES =  main.o auxil1.o
FILES2 = main_Stats.o auxil1.o

all: main.ex main_Stats.ex

main_Stats.ex: $(FILES2)
	${CC} ${CFLAGS} -o main_Stats.ex $(FILES2) -lm $(LIBS)

main.ex: $(FILES) 
	${CC} ${CFLAGS} -o main.ex -lm $(FILES) $(LIBS)

.c.o:
	${CC} ${CFLAGS} $< -c -o $@ $(LIBS) 

clean:
	rm *.o *.ex
