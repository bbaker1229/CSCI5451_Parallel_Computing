CC	= mpicc 
CFLAGS  = -g -Wall 

main.ex: main.o aux_1.o  
	$(CC) -o main.ex $(CFLAGS) -lm main.o aux_1.o

conv.ex: conv.o aux_1.o
	$(CC) -o conv.ex $(CFLAGS) -lm conv.o aux_1.o

fixed_grid.ex: fixed_grid.o aux_1.o
	$(CC) -o fixed_grid.ex $(CFLAGS) -lm fixed_grid.o aux_1.o

clean:
	'rm' *.o *.ex #*~ #*

.c.o:	
	$(CC) -c $(CFLAGS)  $< 
