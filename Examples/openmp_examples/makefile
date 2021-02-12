.SUFFIXES: .c .o .f .f

CC                      = gcc
CFLAGS                  = -fopenmp

FILES = testAxpy.o timer.o

dotp.ex: dotp.o
		${CC} ${CFLAGS} -o dotp.ex -lm dotp.o

dotp1.ex: dotp1.o
		${CC} ${CFLAGS} -o dotp1.ex -lm dotp1.o

dotp2.ex: dotp2.o
		${CC} ${CFLAGS} -o dotp2.ex -lm dotp2.o

Hello0.ex: Hello0.o
		${CC} ${CFLAGS} -o Hello0.ex -lm Hello0.o

Hello1.ex: Hello1.o
		${CC} ${CFLAGS} -o Hello1.ex -lm Hello1.o

racing.ex: racing.o
		${CC} ${CFLAGS} -o racing.ex -lm racing.o

racing1.ex: racing1.o
		${CC} ${CFLAGS} -o racing1.ex -lm racing1.o

racing2.ex: racing2.o
		${CC} ${CFLAGS} -o racing2.ex -lm racing2.o

reduc.ex: reduc.o
		${CC} ${CFLAGS} -o reduc.ex -lm reduc.o

axpy.ex: $(FILES)
		${CC} ${CFLAGS} -o axpy.ex -lm $(FILES)

.c.o:
		${CC} ${CFLAGS} $< -c -o $@

clean:
		rm *.o *.ex