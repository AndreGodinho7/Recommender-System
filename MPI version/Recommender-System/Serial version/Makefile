CC = gcc

  #  -g    adds debugging information to the executable file
  #  -Wall turns on most, but not all, compiler warnings 

CFLAGS  = -g -Wall

OBJS = matFact.o input.o factorization.o
EXEC = matFact

$(EXEC): $(OBJS)
	$(CC) $(CFFLAGS) -o $(EXEC) $(OBJS)

matFact.o: matFact.c input.h
	$(CC) $(CFLAGS) -c matFact.c

input.o: input.c input.h
	$(CC) $(CFLAGS) -c input.c

factorization.o: factorization.c factorization.h
	$(CC) $(CFLAGS) -c factorization.c

clean:
		echo "Delete .o files and executable."
		rm -f *.o $(EXEC)