#include "breakpoint_graph_euler_fontier_merge_multiple_SV_dependency_precut_parallel_nodes_full.h"
#include <iostream>
#include <string.h>
#include <algorithm>
#include <list>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <time.h>
#include <map>
#include <unordered_map>
#include <limits>
#include <thread>
#include <mutex>
#include <atomic>

//list<int> DAG_expand(int start_nodes_n);
int main(int argc, char* argv[]){
	string is_prerun = string(argv[1]);
	int target_node=atoi(argv[2]);
	int final_node=atoi(argv[3]);

	int total_nodes=atoi(argv[4]);
        int total_edges=atoi(argv[5]);
        int hidden_node=atoi(argv[6]);
        int hidden_edge=atoi(argv[7]);
 

	cout<<"stat"<<total_nodes<<"\t"<<total_edges<<"\t"<<hidden_node<<"\t"<<hidden_edge<<endl;
	if(string(argv[1])=="true"){
		vector<int> level_index_list = DAG_expand(true,0, 0 , total_nodes, total_edges, hidden_node, hidden_edge);
		cout<<"level_index"<<"\t"<<level_index_list[0]<<"\t"<<level_index_list[1]<<endl;
	}
	else{
                vector<int> level_index_list = DAG_expand(false,atoi(argv[2]),atoi(argv[3]), total_nodes, total_edges, hidden_node, hidden_edge);
	}
	return 1;
}
