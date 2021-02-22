CC	= gcc
FLAGS	= -fopenacc -lcuda -fcf-protection=none 

test.ex: VecAdd.o
	$(CC) $(FLAGS) VecAdd.o -o test.ex -lm 

clean:
	rm -f *.o *~ *.ex
