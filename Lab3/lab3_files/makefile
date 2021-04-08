CC	= mpicc 
CFLAGS  = -g -Wall 

randw.ex: main.o aux.o  
	$(CC) -o randw.ex $(CFLAGS) -lm main.o aux.o

clean:
	'rm' *.o *.ex *~ #*

.c.o:	
	$(CC) -c $(CFLAGS)  $< 
