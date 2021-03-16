.SUFFIXES: .c .o .f .F

##CC			=  gcc
CC			=  mpicc 
CFLAGS			= -g -Wall
##LIBS			=  -L/project/scicom/scicom00/SOFT/lib/linux/ -lblas
LIBS			=  -lblas 

FILES =  main.o auxil1.o
FILES2 = main_Stats.o auxil1.o
FILES3 = main_Stats.o test_auxil1.o
FILES4 = main.o test_auxil1.o

all: main.ex main_Stats.ex fix_main_Stats.ex fix_main.ex

fix_main.ex: $(FILES4)
	${CC} ${CFLAGS} -o fix_main.ex $(FILES4) -lm $(LIBS)

fix_main_Stats.ex: $(FILES3)
	${CC} ${CFLAGS} -o fix_main_Stats.ex $(FILES3) -lm $(LIBS)

main_Stats.ex: $(FILES2)
	${CC} ${CFLAGS} -o main_Stats.ex $(FILES2) -lm $(LIBS)

main.ex: $(FILES) 
	${CC} ${CFLAGS} -o main.ex -lm $(FILES) $(LIBS)

.c.o:
	${CC} ${CFLAGS} $< -c -o $@ $(LIBS) 

clean:
	rm *.o *.ex
