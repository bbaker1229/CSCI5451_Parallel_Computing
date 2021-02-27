CC	= gcc
FLAGS	= -fopenacc -lcuda -fcf-protection=none 

test.ex: Mat_add.o timer.o
	$(CC) $(FLAGS) Mat_add.o timer.o -o test.ex -lm 

.c.o:
	$(CC) $(FLAGS) $< -c -o $@

clean:
	rm -f *.o *~ *.ex
