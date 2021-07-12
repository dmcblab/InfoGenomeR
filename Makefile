PROG1 = haplotype_graph/phasing/phasing
PROG2 = Eulerian/DAG_entropy/main

CC = g++
CPPFLAGS = --std=c++11
OBJS1 = haplotype_graph/phasing/dag_N.o
all: $(PROG1) $(PROG2)
$(PROG1) : $(OBJS1)
	$(CC) -o $(PROG1) $(OBJS1)
dag_N.o : 
	$(CC) $(CPPFLAGS) -c haplotype_graph/phasing/dag_N.cpp
$(PROG2) : 
	$(CC) $(CPPFLAGS) -pthread Eulerian/DAG_entropy/main.cpp Eulerian/DAG_entropy/breakpoint_graph_euler_fontier_merge_multiple_SV_dependency_precut_parallel_nodes_full.cpp Eulerian/DAG_entropy/set_cover_limited.cpp Eulerian/DAG_entropy/combination_precut_test.cpp -o Eulerian/DAG_entropy/main

clean:
	rm -f core $(PROG1) $(OBJS1) $(PROG2)

