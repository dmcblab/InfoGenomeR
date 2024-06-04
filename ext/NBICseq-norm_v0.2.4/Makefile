CC = gcc
CFLAGS = -g -O2 -lm -Wall

all: combineFile PrepPois PrepPoisGAM BRS filterCNV purityEM random_split

combineFile: combineFile/combine.o lib/read.o
	$(CC) $(CFLAGS) combineFile/combine.o lib/read.o -o combineFile/combineFile

PrepPois: DataPrepare/DtaPrep.o lib/read.o lib/statlib.o
	$(CC) $(CFLAGS)  DataPrepare/DtaPrep.o lib/read.o lib/statlib.o -o DataPrepare/PrepPois

PrepPoisGAM: DataPrepareGAM/DtaPrepGAM.o lib/read.o lib/statlib.o
	$(CC) $(CFLAGS) DataPrepareGAM/DtaPrepGAM.o lib/read.o lib/statlib.o -o DataPrepareGAM/PrepPoisGAM

BRS: BRS_1bp/BRS.o BRS_1bp/pos_cnt_lst.o lib/read.o
	$(CC) $(CFLAGS) BRS_1bp/BRS.o BRS_1bp/pos_cnt_lst.o lib/read.o -o BRS_1bp/BRS

filterCNV: Filter/filter.o lib/read.o
	$(CC) $(CFLAGS) Filter/filter.o lib/read.o -o Filter/filterCNV

purityEM: purityEM/purityEM.o purityEM/EM.o lib/read.o purityEM/gamma.o
	 $(CC) $(CFLAGS) purityEM/purityEM.o purityEM/EM.o lib/read.o purityEM/gamma.o -o  purityEM/purityEM

random_split: randomSample/split.o
	$(CC) $(CFLAGS) randomSample/split.o -o randomSample/random_split


combine.o: combineFile/read.h
	 $(CC) $(CFLAGS) -c combineFile/combine.c lib/read.c

DtaPrep.o: lib/read.h lib/statlib.h
	 $(CC) $(CFLAGS) -c DataPrepare/DtaPrep.c lib/read.c lib/statlib.c

PrepPoisGAM.o: lib/read.h lib/statlib.h
	$(CC) $(CFLAGS) -c DataPrepareGAM/DtaPrepGAM.c lib/read.c lib/statlib.c

BRS.o: BRS_1bp/pos_cnt_lst.h lib/read.h
	$(CC) $(CFLAGS) -c BRS_1bp/pos_cnt_lst.c lib/read.c

filter.o: lib/read.h
	$(CC) $(CFLAGS) -c Filter/filter.c lib/read.c


purityEM.o: purityEM/EM.h lib/read.h purityEM/gamma.h
	$(CC) $(CFLAGS) -c purityEM/purityEM.c purityEM/EM.c lib/read.c purityEM/gamma.c

split.o:
	 $(CC) $(CFLAGS) -c randomSample/split.c



EM.o: purityEM/EM.h
	 $(CC) $(CFLAGS) -c purityEM/EM.c

gamma.o: purityEM/gamma.h
	$(CC) $(CFLAGS) -c purityEM/gamma.c


pos_cnt_lst.o: BRS_1bp/pos_cnt_lst.h
	$(CC) $(CFLAGS) -c BRS_1bp/pos_cnt_lst.c


statlib.o: lib/statlib.c
	$(CC) $(CFLAGS) -c lib/statlib.c

read.o: lib/read.h
	$(CC) $(CFLAGS) -c lib/read.c

clean:
	rm -rf combineFile/*.o DataPrepare/*.o lib/*.o DataPrepareGAM/*.o BRS_1bp/*.o Filter/*.o purityEM/*.o randomSample/*.o

