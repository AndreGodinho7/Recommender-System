CC = gcc

  #  -g    adds debugging information to the executable file
  #  -Wall turns on most, but not all, compiler warnings 

CFLAGS  = -g -Wall

OBJS = main.o input.o factorization.o
EXEC = recsystem

$(EXEC): $(OBJS)
	$(CC) $(CFFLAGS) -o $(EXEC) $(OBJS)

main.o: main.c input.h
	$(CC) $(CFLAGS) -c main.c

input.o: input.c input.h
	$(CC) $(CFLAGS) -c input.c

factorization.o: factorization.c factorization.h
	$(CC) $(CFLAGS) -c factorization.c

clean:
		echo "Apagar todos os ficheiros objeto e executavel"
		rm -f *.o $(EXEC)