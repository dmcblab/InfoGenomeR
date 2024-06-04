CC = gcc
CFLAGS = -g -O2 -lm -Wall


PrepPois: DtaPrep.o ../lib/read.o ../lib/statlib.o
	$(CC) $(CFLAGS) DtaPrep.o ../lib/read.o ../lib/statlib.o -o PrepPois

DtaPrep.o: ../lib/read.o ../lib/statlib.o
	 $(CC) $(CFLAGS) -c DtaPrep.c ../lib/read.c ../lib/statlib.c

read.o: ../lib/read.h
	 $(CC) $(CFLAGS) -c ../lib/read.c

statlib.o: ../lib/statlib.c
	$(CC) $(CFLAGS) -c ../lib/statlib.c

clean:
	rm -rf *.o ../lib/*.o
