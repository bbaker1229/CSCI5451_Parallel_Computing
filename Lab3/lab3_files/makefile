CC	= mpicc 
CFLAGS  = -g -Wall 

randw.ex: main.o aux_1.o  
	$(CC) -o randw.ex $(CFLAGS) -lm main.o aux_1.o

clean:
	'rm' *.o *.ex #*~ #*

.c.o:	
	$(CC) -c $(CFLAGS)  $< 
