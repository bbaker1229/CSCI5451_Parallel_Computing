CC	= gcc
FLAGS	= -fopenacc -lcuda -fcf-protection=none 

test.ex: Mat_add.o
	$(CC) $(FLAGS) Mat_add.o -o test.ex -lm 

clean:
	rm -f *.o *~ *.ex
