CC = mpicc

  #  -g    adds debugging information to the executable file
  #  -Wall turns on most, but not all, compiler warnings 

CFLAGS  = -fopenmp

OBJS = matFact.o input.o factorization.o
EXEC = matFact

$(EXEC): $(OBJS)
	$(CC) -fopenmp -o $(EXEC) $(OBJS) -lm

matFact.o: matFact.c input.h
	$(CC) $(CFLAGS) -c matFact.c

input.o: input.c input.h
	$(CC) $(CFLAGS) -c input.c

factorization.o: factorization.c factorization.h
	$(CC) $(CFLAGS) -c factorization.c

clean:
		echo "Delete .o files and executable."
		rm -f *.o $(EXEC)
