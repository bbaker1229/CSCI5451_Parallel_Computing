CC	= mpicc 
CFLAGS  = -g -Wall 

main.ex: main.o aux_1.o  
	$(CC) -o main.ex $(CFLAGS) -lm main.o aux_1.o

clean:
	'rm' *.o *.ex #*~ #*

.c.o:	
	$(CC) -c $(CFLAGS)  $< 
