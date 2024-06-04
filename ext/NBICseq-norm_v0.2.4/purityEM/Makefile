CC = gcc
CFLAGS = -g -O2 -lm -Wall

purityEM: purityEM.o EM.o read.o gamma.o
	$(CC) $(CFLAGS) purityEM.o EM.o gamma.o read.o -o purityEM

purityEM.o: EM.h read.h gamma.h
	$(CC) $(CFLAGS) -c purityEM.c EM.c read.c gamma.c

EM.o: EM.h
	$(CC) $(CFLAGS) -c EM.c
gamma.o: gamma.h
	$(CC) $(CFLAGS) -c gamma.c

read.o: read.h
	$(CC) $(CFLAGS) -c read.c

clean:
	rm -rf *.o

