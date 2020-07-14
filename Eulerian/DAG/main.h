// A C++ program print Eulerian Trail in a given Eulerian or Semi-Eulerian Graph
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
#include "combination_precut.h"
#include "set_cover_branch_limit.h"
using namespace std;
int sign = 1;
int branch_limit=99999;
int thread_n=1;
unordered_map<int,int> node_cantor_to_edge_index;

mutex mtx;

struct full_edges{
	int src;
	int dest;
	int degree1;
	int degree2;
	list<int> edges_index;
	int edge_by_genome_index;
	int multiplicity;
	int ***wtis;
	int type;
};

class Graph
{
public:
	int **edge_table;
	int **edge_by_g_adj_matrix;
	int **edge_adj_matrix;
        int hidden_node;
        int V; // No. of vertices
        int E;
	map<int,int> genome_index_to_type;
	map<int,vector<int>> genome_index_to_nodes;
        map<vector<int>,int> nodes_to_genome_index;
	multimap<vector<int>,int> nodes_to_edge_index;

	bool genome_edge_is_adjacent(int node, int genome_edge);
        list< vector<int> > *adj; // Pointer to an array containing adjacency lists
        Graph(int V, int E, int hn, int** edge_table); // Constructor
        ~Graph() { delete [] adj; for(int i=0; i<E; i++){delete [] edge_by_g_adj_matrix[i]; delete [] edge_adj_matrix[i];}  delete [] edge_by_g_adj_matrix; delete [] edge_adj_matrix;}
        void addEdge(int v, int w, int type, int index); // function to add an edge to graph
	void rmvEdge(int u, int v, int type);
        bool crossable(int type1, int type2);
        bool isConnected(int l_edge_table_n);
        void DFSUtil(int old_v, int v, bool visited[], bool edge_visited[]);
	void edge_adj_matrix_constr();
	void edge_adj_matrix_print();
	bool genome_edges_are_adjacent(int i, int j);
	list<int> BFS_except_hidden_by_edges(int s);
	void gtest();
	void DFSUtil_fromSV(int old_type, int v, bool visited[],bool edge_visited[], list<vector<int>> SV_adj[]);
	bool isConnected_fromSV(int l_edge_table_n, list<list<int>> pnodes, list<list<int>> ppaths,  list<vector<int>> SV_adj[]);
	bool isConnected_fromSV_src_to_dest(int l_edge_table_n, list<list<int>> pnodes, list<list<int>> ppaths, list<vector<int>> SV_adj[], int src, int dest);
};


void Graph::gtest(){
	cout<<"gtest"<<endl;
	for(int i=0; i<29; i++){
		for(int j=0; j<29; j++){
			cout<<edge_by_g_adj_matrix[i][j]<<"\t";
		}
		cout<<endl;
	}
	cout<<endl;

}

bool Graph::genome_edge_is_adjacent(int node, int genome_edge){
	for(list<vector<int>>::iterator it=adj[node].begin();it!=adj[node].end();it++){
		if(edge_table[(*it).at(2)][6] == genome_edge){
			return true;
		}
	}
	return false;
}

void Graph::edge_adj_matrix_print(){
	for(int i=0; i<5;i++){
		for(int j=0; j<5; j++){
			cout<<edge_adj_matrix[i][j]<<"\t";
		}
		cout<<endl;
	}
}
Graph::Graph(int V, int E,int hn, int** edge_table)
{
        hidden_node = hn;
        this->V = V;
        this->E = E;
        adj = new list< vector<int> >[V];

	edge_adj_matrix = new int*[E];
	for(int i=0; i<E; i++){
		edge_adj_matrix[i]=new int[E];
	}

 	this->edge_table = edge_table;
	edge_by_g_adj_matrix = new int*[E];
        for(int i=0; i<E; i++){
                edge_by_g_adj_matrix[i]=new int [E];
        }
}

void Graph::edge_adj_matrix_constr(){
	for(int i=0; i<edge_table[E-1][6];i++){
		for(int j=0; j<edge_table[E-1][6];j++){
			edge_by_g_adj_matrix[i][j]=0;
		}
	}
        for(int i=0; i<V; i++){
                for(list<vector<int>>::iterator it=adj[i].begin();it!=adj[i].end();it++){
                        for(list<vector<int>>::iterator jt=adj[i].begin(); jt!=adj[i].end();jt++){
                                edge_adj_matrix[(*it).at(2)][(*jt).at(2)]=1;
                                edge_by_g_adj_matrix[edge_table[(*it).at(2)][6]][edge_table[(*jt).at(2)][6]]=1;
	//			cout<<edge_table[(*it).at(2)][6]<<"\t"<<edge_table[(*jt).at(2)][6]<<endl;
			//	cout<<edge_by_g_adj_matrix[edge_table[(*it).at(2)][6]][edge_table[(*jt).at(2)][6]]<<endl;
                        }
                }
        }
}
bool Graph::genome_edges_are_adjacent(int i, int j){
        return edge_by_g_adj_matrix[i][j];
}


void Graph::addEdge(int v, int w, int type, int index)
{
        adj[w].push_back({v,type,index});
        adj[v].push_back({w,type,index}); // Add w to v’s list.
}

void Graph::rmvEdge(int u, int v, int type)
{
  for(list< vector<int> >::iterator iv=adj[u].begin(); iv!=adj[u].end();iv++){
	if((*iv).at(0)==v && (*iv).at(1)==type){
		(*iv).at(0)=-1;
		break;
	}
  }
  for(list< vector<int> >::iterator iu=adj[v].begin(); iu!=adj[v].end();iu++){
 	if((*iu).at(0)==u && (*iu).at(1)==type){
		(*iu).at(0)=-1;
		break;
	}
  }

}
list<int> Graph::BFS_except_hidden_by_edges(int s){
	
	list<int> BFS;
	list<int> hidden_edges;
	bool *visited = new bool[E];
	for(int i = 0; i < E; i++)
       		visited[i] = false;
	list<int> queue;
	//visited[s] = true;
	queue.push_back(s);
	list<vector<int>>::iterator i;
	while(!queue.empty()){
		s = queue.front();
		queue.pop_front();
	       for (i = adj[s].begin(); i != adj[s].end(); ++i){
		    if ( !visited[(*i).at(2)]){
			if((*i).at(0)==hidden_node){
                                visited[(*i).at(2)] = true;
				hidden_edges.push_back((*i).at(2));
			}else{
				visited[(*i).at(2)] = true;
				queue.push_back((*i).at(0));
	                        BFS.push_back((*i).at(2));
			}
		    }
		}
	}
	for(list<int>::iterator it=hidden_edges.begin();it!=hidden_edges.end();it++){
		BFS.push_back(*it);
	}
	return BFS;
}


void Graph::DFSUtil(int old_type, int v, bool visited[], bool edge_visited[])
{
    visited[v] = true;
    cout<<v<<"\t";
    list< vector<int> >::iterator i;
    for (i = adj[v].begin(); i != adj[v].end(); ++i){
 //       if ((*i).at(0) != -1 && !visited[(*i).at(0)] && (old_type==-1 || crossable(old_type,(*i).at(1) ))){
        if ((*i).at(0) != -1 && !edge_visited[(*i).at(2)] && (old_type==-1 || crossable(old_type,(*i).at(1) ))){
	    edge_visited[(*i).at(2)]=true;
            DFSUtil((*i).at(1),(*i).at(0), visited, edge_visited);
	}
    }
}
void Graph::DFSUtil_fromSV(int old_type, int v, bool visited[], bool edge_visited[] , list<vector<int>> SV_adj[])
{
   // cout<<v<<"\t";
    visited[v] = true;

    for (list<vector<int>>::iterator i = adj[v].begin(); i != adj[v].end(); ++i){
     //   cout<<"DFSUtil"<<v<<"\t"<<(*i).at(0)<<"\t"<<(*i).at(1)<<"\t"<<(*i).at(2)<<endl;
        if ((*i).at(0) != -1 && !edge_visited[(*i).at(2)] && (old_type==-1 || crossable(old_type,(*i).at(1) ))){
            edge_visited[(*i).at(2)]=true;
            DFSUtil_fromSV((*i).at(1),(*i).at(0), visited,edge_visited, SV_adj);
        }
    }

    for(list<vector<int>>::iterator i = SV_adj[v].begin(); i!=SV_adj[v].end(); i++){
//	cout<<"DFSUtilSV"<<v<<"\t"<<(*i).at(0)<<"\t"<<(*i).at(1)<<"\t"<<(*i).at(2)<<"\t"<<old_type<<endl;
        if ((*i).at(2) && (*i).at(0) != -1  && (old_type==-1 || crossable(old_type,(*i).at(1)))){
//           DFSUtil_fromSV((*it).at(1),(*i).at(0), visited,edge_visited, SV_adj);
	   (*i).at(2)=0;
           DFSUtil_fromSV((*i).at(1),(*i).at(0), visited,edge_visited, SV_adj);

	}
    }	
}
bool Graph::isConnected_fromSV_src_to_dest(int l_edge_table_n, list<list<int>> pnodes, list<list<int>> ppaths, list<vector<int>> SV_adj[], int src, int dest){
    bool visited[V];
    bool edge_visited[l_edge_table_n];
    int i;
    for (i = 0; i < V; i++){
        visited[i] = false;
    }
    for( i =0; i<l_edge_table_n;i++){
        edge_visited[i] = false;
   }


    for(list<list<int>>::iterator it = ppaths.begin(); it!=ppaths.end(); it++){
                for(list<int>::iterator kt= (*it).begin(); kt!=(*it).end(); kt++){
                        edge_visited[*kt]=true;
                }
    }

   for(list<list<int>>::iterator it = pnodes.begin(); it!=pnodes.end(); it++){
                for(list<int>::iterator kt= (*it).begin(); kt!=(*it).end(); kt++){
		if(kt!=(*it).begin() && kt!=(*it).end()){
		       visited[*kt]=true;
		}
	}
   }

cout<<"DFS:";
DFSUtil_fromSV(-1,src, visited,edge_visited, SV_adj);

cout<<endl;

if(visited[dest] == false){
	cout<<"novisited"<<endl;
	return false;
}

return true;
}


bool Graph::isConnected_fromSV(int l_edge_table_n, list<list<int>> pnodes, list<list<int>> ppaths, list<vector<int>> SV_adj[]){
    bool visited[V];
    bool edge_visited[l_edge_table_n];
    int i;
    for (i = 0; i < V; i++){
        visited[i] = false;
    }
    for( i =0; i<l_edge_table_n;i++){
        edge_visited[i] = false;
   }


    for(list<list<int>>::iterator it = ppaths.begin(); it!=ppaths.end(); it++){
		for(list<int>::iterator kt= (*it).begin(); kt!=(*it).end(); kt++){
			edge_visited[*kt]=true;
		}
    }

   for(list<list<int>>::iterator it = pnodes.begin(); it!=pnodes.end(); it++){
                for(list<int>::iterator kt= (*it).begin(); kt!=(*it).end(); kt++){
                        if(kt!=(*it).begin() && kt!=(*it).end()){
 	                       visited[*kt]=true;
			}
		}
    }

    for (i = 0; i < V; i++)
        if (adj[i].size() != 0 && visited[i]==false)
            break;

    if (i == V)
        return true;

   cout<<"DFS:";
    DFSUtil_fromSV(-1,i, visited,edge_visited, SV_adj);

  cout<<endl;

    for (i = 0; i < V; i++)
       if (visited[i] == false && adj[i].size() > 0){
           cout<<"novisited:"<<"\t"<<i<<endl;
            return false;
        }

    return true;
}



bool Graph::isConnected(int l_edge_table_n)
{ 
    bool visited[V];
    bool edge_visited[l_edge_table_n];
//    bool edge_visited[E];
    int i;
    for (i = 0; i < V; i++){
        visited[i] = false;
    }
    for( i =0; i<l_edge_table_n;i++){
        edge_visited[i] = false;
   }
    for( i=0; i<V;i++){
            if(adj[i].size()!=0){
		    int all_negative=1;
		    for(list< vector<int> >::iterator it=adj[i].begin();it!=adj[i].end();it++){
			if((*it).at(0)>=0){
				all_negative=0;
			}
		    }
		    if(all_negative==1){
			visited[i]=true;
		    }
	    }
   }


    for (i = 0; i < V; i++)
        if (adj[i].size() != 0 && visited[i]==false)
            break;


    if (i == V)
        return true;

   cout<<"DFS:";
    DFSUtil(-1,i, visited, edge_visited);
  cout<<endl;

    for (i = 0; i < V; i++)
       if (visited[i] == false && adj[i].size() > 0){
	   cout<<"novisited:"<<"\t"<<i<<endl;
            return false;
	}

    return true;
}

bool Graph::crossable(int type1, int type2){
	if(type1 == 3 || type2 ==3){
		return true;
	}else if((type1 == 0 && type2 > 0) || (type1 > 0 && type2 == 0)){
		return true;
	}
	return false;

}

struct vertex_info{
	int vertex_i;
	int** paired;
	int vertex_degrees;
};

struct AdjListNode
{
    int dest;
    int vertex_n;
   // struct vertex_info* edge_pair;
    struct AdjListNode* next;

    int path_n;
   // list<int>* path_frags;
//    list< list<int> > SV_paths;
    list< list<int> > gnodes;
    list< list<int> > frontier_types;
    list<list<int>> epaths;
    list<int> SVs_Ns;
    list<int> SVs_CNs;
    
  //  list< list<list<int>>* > equivalent_SV_paths;
    
  //  list<list<list<int>>> SV_paths_list;
    int SV_paths_magic_number;

};
/*
 *  * Adjacency List
 *   */  
struct AdjList
{
    struct AdjListNode *head;
};

class DAG
{
    private:
        int V;
        struct AdjList* array;
	struct AdjListNode** nodes;
	int level;
	int thres_level;
	int *level_index;
	int hidden_node;
	int hidden_edge;
	int *edge_by_genome;

	int** degree_table;
	int degree_table_n;
	int** edge_table;
	int edge_table_n;
	
    public:
	Graph* gp;
        DAG(int V, int hn, int he, int genome_edges, int** degree_table, int** edge_table, int degree_table_n, int edge_table_n, Graph* gp)
        {
            this->V = V;
            array = new AdjList [V];
            level_index = new int [V];
            for (int i = 0; i < V; ++i){
                array[i].head = NULL;
	    }
	   nodes= new AdjListNode*[V];
	   level=0;
	   thres_level=10;
	   level_index[0]=0;
	   hidden_node=hn;
	   hidden_edge=he;
	   edge_by_genome = new int [genome_edges];
	   this->degree_table = new int* [degree_table_n];
	   for(int i=0; i< degree_table_n ; i++){
		this->degree_table[i] = new int [3];
		for(int j=0; j<3; j++){
			this->degree_table[i][j] = degree_table[i][j];
		}
	   }
	   this->edge_table = edge_table;
	   this->degree_table_n = degree_table_n;
	   this->edge_table_n = edge_table_n;
	   this->gp = gp;
        }
        /*
 *          * Creating New Adjacency List Node
 *               */ 
	int call_level_index(int level){
		return level_index[level];
	}
	int call_level(){
		return level;
	}
	void set_target(int current_level, int target_index, int final_node){
		nodes[final_node-1] = nodes[target_index];
		level_index[current_level-1]=final_node-1;
	}

	void wtis_edit(int*** wtis, int wtis_size, int original_index, int original_orient, int new_index, int new_orient){
		for(int i=0; i<2; i++){
			for(int j=0; j<wtis_size; j++){
				if(wtis[i][j][0] == original_index && wtis[i][j][1] == original_orient){
					wtis[i][j][0] = new_index;
					wtis[i][j][1] = new_orient;
				}
			}
		}
	}


        AdjListNode* newAdjListNode(int parent, int offspring, list<full_edges> full_edges_list, int*** wtis)
        {
		AdjListNode* newNode = new AdjListNode;
		newNode->dest = offspring;
		newNode->next = NULL;
	//	list<list<int>> gnodes;
	//	list<list<int>> frontier_types;
		
		int wtis_size = 0;
                for(list<full_edges>::iterator flt = full_edges_list.begin(); flt!=full_edges_list.end(); flt++){
			wtis_size += (*flt).multiplicity;
		}
		cout<<"wtis_size"<<wtis_size<<endl;
		if(offspring!=0){
			newNode->gnodes = nodes[parent]->gnodes;
                        newNode->frontier_types=nodes[parent]->frontier_types;
			newNode->epaths = nodes[parent]->epaths;
			newNode->SVs_Ns = nodes[parent]-> SVs_Ns;
			newNode->SVs_CNs = nodes[parent]->SVs_CNs;
		}	
               if(newNode->gnodes.size()==0){
			for(list<full_edges>::iterator flt = full_edges_list.begin(); flt!=full_edges_list.end(); flt++){
				for(int i=0; i<(*flt).multiplicity; i++){	
					list<int>::iterator it = (*flt).edges_index.begin();
					list<int> gnode;
					list<int> frontier_type;
					list<int> epath;
					gnode.push_back((*flt).src);
					gnode.push_back((*flt).dest);
					newNode->gnodes.push_back(gnode);
					frontier_type.push_back(gp->genome_index_to_type.find((*flt).edge_by_genome_index)->second);
					newNode->frontier_types.push_back(frontier_type);
					epath.push_back(*it);
					it++;
					newNode->epaths.push_back(epath);
					if((*flt).type == 2){
						newNode->SVs_Ns.push_back(1);
						newNode->SVs_CNs.push_back((*flt).multiplicity);
					}else{
                                                newNode->SVs_Ns.push_back(0);
                                                newNode->SVs_CNs.push_back(0);
					}
				}
			}
		}else{
			      list<list<list<int>>::iterator> erases ;
			      list<list<list<int>>::iterator> erases_epaths;
			      list<list<list<int>>::iterator> erases_frontier_types;
                              list<list<int>::iterator> erases_SVs_Ns;
                              list<list<int>::iterator> erases_SVs_CNs;

				int i=0;
				for(list<full_edges>::iterator flt = full_edges_list.begin(); flt!=full_edges_list.end(); flt++){
				        int src = (*flt).src;
					int dest = (*flt).dest;
					int degree1 = (*flt).degree1;
					int degree2 = (*flt).degree2;
					list<int> edges_index = (*flt).edges_index;
					int edge_by_genome_index = (*flt).edge_by_genome_index;
					int multiplicity = (*flt).multiplicity;
 		                        int type = gp->genome_index_to_type.find((*flt).edge_by_genome_index)->second;

	                                list<int>::iterator eit = edges_index.begin();
	//				cout<<"dest"<<dest<<endl;
					cout<<"multiplicity"<<multiplicity<<endl;
					for(int j=0; j<multiplicity; j++){
						if((wtis[0][i][0] > 0 && wtis[0][i][0] >= newNode->gnodes.size()) || (wtis[0][i][0] <0 && wtis[1][i][0] <0)){
							list<int> gnode;
							list<int> frontier_type;
							list<int> epath;
							gnode.push_back(src);
							gnode.push_back(dest);
							newNode->gnodes.push_back(gnode);
							frontier_type.push_back(type);
							newNode->frontier_types.push_back(frontier_type);
							epath.push_back(*eit);
							eit++;
							newNode->epaths.push_back(epath);
							if((*flt).type == 2){
								newNode->SVs_Ns.push_back(1);
								newNode->SVs_CNs.push_back((*flt).multiplicity);
							}else{
							        newNode->SVs_Ns.push_back(0);
								newNode->SVs_CNs.push_back(0);
							}

						}else{
							int sod_start;
							int sod_end;
							if((wtis[0][i][0] >=0 &&  wtis[1][i][0]<0) || (wtis[0][i][0] < 0 && wtis[1][i][0]>=0)){
								if(wtis[0][i][0] >=0 &&  wtis[1][i][0]<0){
									sod_start=0;
									sod_end=0;
								}else if((wtis[0][i][0] < 0 && wtis[1][i][0]>=0)){

									sod_start=1;
									sod_end=1;
								}
								for(int sod = sod_start; sod<=sod_end; sod++){
									list<list<int>>::iterator jt = newNode->gnodes.begin();
									list<list<int>>::iterator kt = newNode->frontier_types.begin();
									list<list<int>>::iterator lt = newNode->epaths.begin();
									list<int>::iterator zt = newNode->SVs_Ns.begin();
									list<int>::iterator xt = newNode->SVs_CNs.begin();
									int wt=wtis[sod][i][0];
									while(wt>0){
										jt++;
										kt++;
										lt++;
										zt++;
										xt++;
										wt--;
									}
									if(wtis[sod][i][1]){
										(*jt).push_back(wtis[sod][i][2]);
										if((*kt).size()>1)
											(*kt).pop_back();
										(*kt).push_back(type);
										(*lt).push_back(*eit);
										eit++;
									}else{
										(*jt).push_front(wtis[sod][i][2]);
										if((*kt).size()>1)
											(*kt).pop_front();
										(*kt).push_front(type);
										(*lt).push_front(*eit);
										eit++;
									}
									if(type==2){
										(*zt)++;
										(*xt)+=multiplicity;
									}
								}
							}else if(wtis[0][i][0] >=0  && wtis[1][i][0]>=0){

								list<list<int>>::iterator path1_it = newNode->gnodes.begin();
								list<list<int>>::iterator path1_kt = newNode->frontier_types.begin();
								list<list<int>>::iterator path1_lt = newNode->epaths.begin();
								list<int>::iterator path1_zt = newNode->SVs_Ns.begin();
                                                                list<int>::iterator path1_xt = newNode->SVs_CNs.begin();
								int wt1=wtis[0][i][0];
								while(wt1>0){
									path1_it++;
									path1_kt++;
									path1_lt++;
									path1_zt++;
									path1_xt++;
									wt1--;
								}
							       list<list<int>>::iterator path2_it = newNode->gnodes.begin();
							       list<list<int>>::iterator path2_kt = newNode->frontier_types.begin();
							       list<list<int>>::iterator path2_lt = newNode->epaths.begin();
                                                               list<int>::iterator path2_zt = newNode->SVs_Ns.begin();
                                                               list<int>::iterator path2_xt = newNode->SVs_CNs.begin();

								int wt2=wtis[1][i][0];
								while(wt2>0){
									path2_it++;
									path2_kt++;
									path2_lt++;
									path2_zt++;
									path2_xt++;
									wt2--;
								}
		
								if(path1_it != path2_it){

									(*path1_zt) += (*path2_zt);
									(*path1_xt) += (*path2_xt);					

									if(wtis[0][i][1]){
										if(wtis[1][i][1]){
											for(list<int>::reverse_iterator pit=(*path2_it).rbegin();pit!=(*path2_it).rend(); pit++){
												(*path1_it).push_back(*pit);
											}
											(*path1_lt).push_back(*eit);
											eit++;
										       for(list<int>::reverse_iterator plt=(*path2_lt).rbegin();plt!=(*path2_lt).rend(); plt++){
												(*path1_lt).push_back(*plt);
											}
											if((*path1_kt).size()>1)
												(*path1_kt).pop_back();
											(*path1_kt).push_back((*path2_kt).front());

											wtis_edit ( wtis, wtis_size, wtis[1][i][0], 0, wtis[0][i][0], 1);
				
										}else{
											for(list<int>::iterator pit=(*path2_it).begin();pit!=(*path2_it).end(); pit++){
												(*path1_it).push_back(*pit);
											}
											(*path1_lt).push_back(*eit);
											eit++;
											for(list<int>::iterator plt=(*path2_lt).begin();plt!=(*path2_lt).end(); plt++){
												(*path1_lt).push_back(*plt);
											}
											if((*path1_kt).size()>1)
												(*path1_kt).pop_back();
											(*path1_kt).push_back((*path2_kt).back());

                                                                                        wtis_edit ( wtis, wtis_size, wtis[1][i][0], 1, wtis[0][i][0], 1);

										}

									}else{
									       if(wtis[1][i][1]){
											for(list<int>::reverse_iterator pit=(*path2_it).rbegin();pit!=(*path2_it).rend(); pit++){
												(*path1_it).push_front(*pit);
											}
											(*path1_lt).push_front(*eit);
											eit++;
											for(list<int>::reverse_iterator plt=(*path2_lt).rbegin();plt!=(*path2_lt).rend(); plt++){
												(*path1_lt).push_front(*plt);
											}
											if((*path1_kt).size()>1)
												(*path1_kt).pop_front();
											(*path1_kt).push_front((*path2_kt).front());

                                                                                        wtis_edit ( wtis, wtis_size, wtis[1][i][0], 0, wtis[0][i][0], 0);


										}else{
											for(list<int>::iterator pit=(*path2_it).begin();pit!=(*path2_it).end(); pit++){
												(*path1_it).push_front(*pit);
											}
											(*path1_lt).push_front(*eit);
											eit++;
											for(list<int>::iterator plt=(*path2_lt).begin();plt!=(*path2_lt).end(); plt++){
												(*path1_lt).push_front(*plt);
											}
											if((*path1_kt).size()>1)
												(*path1_kt).pop_front();
											(*path1_kt).push_front((*path2_kt).back());

                                                                                        wtis_edit ( wtis, wtis_size, wtis[1][i][0], 1, wtis[0][i][0], 0);

										}
									}
							//		cout<<"here"<<endl;
                                                                        if(erases.size() == 0 || find(erases.begin(), erases.end(),path2_it) == erases.end()){
                                                                                erases.push_back(path2_it);
                                                                        }
                                                                        if(erases_epaths.size() == 0 || find(erases_epaths.begin(), erases_epaths.end(),path2_lt) == erases_epaths.end()){
                                                                                erases_epaths.push_back(path2_lt);
                                                                        }
                                                                        if(erases_frontier_types.size() == 0 || find(erases_frontier_types.begin(), erases_frontier_types.end(),path2_kt) == erases_frontier_types.end()){
                                                                                erases_frontier_types.push_back(path2_kt);
                                                                        }
                                                                        if(erases_SVs_Ns.size() == 0 || find(erases_SVs_Ns.begin(), erases_SVs_Ns.end(),path2_zt) == erases_SVs_Ns.end()){
                                                                                erases_SVs_Ns.push_back(path2_zt);
                                                                        }
                                                                        if(erases_SVs_CNs.size() == 0 || find(erases_SVs_CNs.begin(), erases_SVs_CNs.end(),path2_xt) == erases_SVs_CNs.end()){
                                                                                erases_SVs_CNs.push_back(path2_xt);
                                                                        }

								//	(*it).erase(path2_it);
								}else{
									if(wtis[0][i][1]){
										(*path1_it).push_back(wtis[0][i][2]); ///////////////////////// cautious!! /////////////////////////////////
										 (*path1_lt).push_back(*eit);
										  eit++;
										if((*path1_kt).size()>1)
											(*path1_kt).pop_back();
									       (*path1_kt).push_back(type);
									}
									else{
										(*path1_it).push_front(wtis[0][i][2]);
										 (*path1_lt).push_front(*eit);
										  eit++;
										if((*path1_kt).size()>1)
											(*path1_kt).pop_front();
									       (*path1_kt).push_front(type);
									}
									(*path1_it).push_front(-1);
                                                                        (*path1_it).push_back(-1);
								}
								if(type==2){
									(*path1_zt)++;
									(*path1_xt)+=multiplicity;
								}
							}
							i++;
						}
					}
				}	
				//cout<<"test1"<<endl;

				list<list<list<list<int>>::iterator>::iterator> used;
				list<list<list<int>::iterator>::iterator> used2;
				for(list<list<list<int>>::iterator>::iterator erase_it=erases.begin(); erase_it!=erases.end(); erase_it++){
						if(used.size()==0 || find(used.begin(), used.end(), erase_it)==used.end()){
							cout<<"erased_check"<<endl;
							newNode->gnodes.erase(*erase_it);
							used.push_back(erase_it);
						}
				}
				used.clear();
                                for(list<list<list<int>>::iterator>::iterator erase_it=erases_epaths.begin(); erase_it!=erases_epaths.end(); erase_it++){
                                                if(used.size() ==0 || find(used.begin(), used.end(), erase_it)==used.end()){
                                                        newNode->epaths.erase(*erase_it);
                                                        used.push_back(erase_it);
                                                }
                                }
				used.clear();
                                for(list<list<list<int>>::iterator>::iterator erase_it=erases_frontier_types.begin(); erase_it!=erases_frontier_types.end(); erase_it++){
                                                if(used.size() == 0 || find(used.begin(), used.end(), erase_it)==used.end()){
                                                      newNode->frontier_types.erase(*erase_it);
                                                        used.push_back(erase_it);
                                                }
                                }
                                used.clear();
                               for(list<list<int>::iterator>::iterator erase_it=erases_SVs_Ns.begin(); erase_it!=erases_SVs_Ns.end(); erase_it++){
                                                if(used2.size() == 0 || find(used2.begin(), used2.end(), erase_it)==used2.end()){
                                                        newNode->SVs_Ns.erase(*erase_it);
                                                        used2.push_back(erase_it);
                                                }
                                }
                                used2.clear();
                                for(list<list<int>::iterator>::iterator erase_it=erases_SVs_CNs.begin(); erase_it!=erases_SVs_CNs.end(); erase_it++){
                                                if(used2.size() == 0 ||find(used2.begin(), used2.end(), erase_it)==used2.end()){
                                                      newNode->SVs_CNs.erase(*erase_it);
                                                        used2.push_back(erase_it);
                                                }
                                }
                                used2.clear();


		//		cout<<"test2"<<endl;
				

			}
//		newNode->frontier_types=frontier_types;

            int i=0; 
		cout<<"gnodes_list"<<endl;
		for(list<list<int>>::iterator jt=newNode->gnodes.begin(); jt!=newNode->gnodes.end();jt++){
			cout<<"gnodes:"<<"\t";
			for(list<int>::iterator kt=(*jt).begin(); kt!=(*jt).end(); kt++){
				cout<<*kt<<"\t";
			}
			cout<<endl;
		//	SV_path_to_inter_nodes(*jt);
		}
                cout<<"epaths"<<endl;
                for(list<list<int>>::iterator jt=newNode->epaths.begin(); jt!=newNode->epaths.end();jt++){
                        cout<<"epaths:"<<"\t";
                        for(list<int>::iterator kt=(*jt).begin(); kt!=(*jt).end(); kt++){
                                cout<<*kt<<"\t";
                        }
                        cout<<endl;
		}
		cout<<"frontier_types"<<endl;
                for(list<list<int>>::iterator jt=newNode->frontier_types.begin(); jt!=newNode->frontier_types.end();jt++){
                        cout<<"frontier_types:"<<"\t";
                        for(list<int>::iterator kt=(*jt).begin(); kt!=(*jt).end(); kt++){
				cout<<(*kt)<<"\t";
			}
			cout<<endl;
                }
                cout<<"SVs_Ns"<<endl;
                for(list<int>::iterator jt=newNode->SVs_Ns.begin(); jt!=newNode->SVs_Ns.end();jt++){
                        cout<<"SV_N:"<<(*jt)<<endl;
                }
                cout<<"SVs_CNs"<<endl;
                for(list<int>::iterator jt=newNode->SVs_CNs.begin(); jt!=newNode->SVs_CNs.end();jt++){
                        cout<<"SV_CN:"<<(*jt)<<endl;
                }


                //      SV_path_to_inter_nodes(*jt);
                //                      }
                //

	    cout<<endl;
            return newNode;
        }

	int cantor_pairing_function(int i, int j){
		return (i+j)*(i+j+1)/2+i;
	}
       void insert_SV_path_to_SV_paths(list<list<int>> *SV_paths, list<int> SV_path){
		list<list<int>>::iterator it=SV_paths->begin();
		it=SV_paths_util(SV_paths, SV_path, it, 1);
		SV_paths->insert(it,SV_path);
       }	
       list<list<int>>::iterator SV_paths_util(list<list<int>> *SV_paths, list<int> SV_path, list<list<int>>::iterator it_start, int which_element){
		
		if(it_start==SV_paths->end()){
			return it_start;
		}
		list<list<int>>::iterator it;
		if(SV_path.size()<which_element){
			return it_start;
		}

		it=it_start;

                
if((*it).size()<which_element){
			it=SV_paths_util(SV_paths, SV_path, ++it, which_element);
		}else{
                        list<int>::iterator compare_SV1=SV_path.begin();
                        list<int>::iterator compare_SV2=(*it).begin();

			for(int i=1; i<which_element;i++){
				if(*compare_SV1>=*compare_SV2){
					compare_SV1++;
					compare_SV2++;
				}else{
					return it;
				}
			}
	                if(*compare_SV1==*compare_SV2){
	 	               it=SV_paths_util(SV_paths, SV_path, it, ++which_element);
			}else if(*compare_SV1>*compare_SV2){
				it=SV_paths_util(SV_paths, SV_path, ++it, which_element);
			}
		}	
		return it;
        }

	bool SV_connectable(list<int> SV_path, bool is_front , int path_size, int edge_by_genome_index){
//		if(path_size == 1){
//			return true;
//		}
		list<int> nodes = SV_path_to_nodes(SV_path, true);
		if(is_front){
			if( nodes.front() == gp->genome_index_to_nodes.find(edge_by_genome_index)->second.at(0) || nodes.front() == gp->genome_index_to_nodes.find(edge_by_genome_index)->second.at(1)){
				return true;
			}
			return false;
		}else{
                        if( nodes.back() == gp->genome_index_to_nodes.find(edge_by_genome_index)->second.at(0) || nodes.back() == gp->genome_index_to_nodes.find(edge_by_genome_index)->second.at(1)){
                                return true;
                        }
                        return false;
		}

	}

	void remove_included_degree(int edge_by_genome_index, int src, int dest, int *degree1, int *degree2, list<list<int>>::iterator begin, list<list<int>>::iterator end){
		for(list<list<int>>::iterator it = begin; it !=end; it++){
			if((*it).size()>1){
				list<int> nodes=SV_path_to_nodes((*it), true);
				list<int> internodes= SV_path_to_inter_nodes((*it), true);
                                cout<<"remove_list:"<<src<<"\t"<<dest<<endl;
				for(list<int>::iterator jt=internodes.begin(); jt !=internodes.end(); jt++){
					cout<<*jt<<"\t";
					if((*jt)==src){
						(*degree1)--;
					}
					if((*jt)==dest){
						(*degree2)--;
					}

				}
				cout<<endl;
				cout<<"path_in_remove_included_degree"<<endl;
                                for(list<int>::iterator jt=nodes.begin(); jt !=nodes.end(); jt++){
					cout<<*jt<<"\t";
				}
				cout<<endl;
				if(nodes.front()==nodes.back() && ((gp->genome_index_to_type.find((*it).front())->second == 0 && gp->genome_index_to_type.find((*it).back())->second == 2) || (gp->genome_index_to_type.find((*it).front())->second == 2 && gp->genome_index_to_type.find((*it).back())->second == 0) ) ){
					if(nodes.front()==src){
                                                cout<<"found error"<<endl;

						(*degree1)--;
					}	
					if(nodes.front()==dest){
                                               (*degree2)--;
					}
				}
			}
		}
	}
       list<int> SV_path_to_nodes(list<int> SV_path, int return_nodes){
		list<int> nodes;
		list<int> types;
		if(return_nodes)
			nodes = SV_path_to_inter_nodes(SV_path, return_nodes);
		else{
			types = SV_path_to_inter_nodes(SV_path, return_nodes);
			return types;
		}

		int n1=gp->genome_index_to_nodes.find(SV_path.front())->second.at(0);
                int n2=gp->genome_index_to_nodes.find(SV_path.front())->second.at(1);
                int n3=gp->genome_index_to_nodes.find(SV_path.back())->second.at(0);
                int n4=gp->genome_index_to_nodes.find(SV_path.back())->second.at(1);
		if(n1==nodes.front()){
			nodes.push_front(n2);
		}else{
			nodes.push_front(n1);
		}
                if(n3==nodes.back()){
                        nodes.push_back(n4);
                }else{
                        nodes.push_back(n3);
                }
		return nodes;

	}


	list<int> SV_path_to_inter_nodes(list<int> SV_path, int return_nodes){
		list<int> nodes;
		list<int> types;
		if(SV_path.size()==1){
			types.push_back(gp->genome_index_to_type.find(SV_path.front())->second);
//			nodes.push_back(gp->genome_index_to_nodes.find(SV_path.front())->second.at(0));
  //                      nodes.push_back(gp->genome_index_to_nodes.find(SV_path.front())->second.at(1));
		}else{
			list<int>::iterator it = SV_path.begin();	
			int previous_edge = (*it);
                        int previous_edge_type=gp->genome_index_to_type.find(previous_edge)->second;

			for(it++; it!=SV_path.end(); it++){
					if(nodes.size()==0){
						types.push_back(gp->genome_index_to_type.find(previous_edge)->second);
					}

					int n1= gp->genome_index_to_nodes.find(previous_edge)->second.at(0);
                                        int n2= gp->genome_index_to_nodes.find(previous_edge)->second.at(1);
                                        int n3= gp->genome_index_to_nodes.find(*it)->second.at(0);
                                        int n4= gp->genome_index_to_nodes.find(*it)->second.at(1);
                                if(gp->genome_index_to_type.find(*it)->second <= 1 && gp->genome_index_to_type.find(previous_edge)->second <= 1){
					int min1=n1<n2? n1:n2;
					int max1=n3>n4 ? n3 : n4;
					if(min1 < max1){
						for(int i=min1+1; i<=max1-1; i++){
							nodes.push_back(i);
							int type= previous_edge_type > 0 ? 0 : 1;
							types.push_back(type);
							previous_edge_type=type;
						}		
					}
					 max1= n1>n2? n1:n2;
					 min1= n3<n4 ? n3:n4;
					if(max1 > min1){
						for( int i=max1-1; i>=min1+1; i--){
							nodes.push_back(i);
                                                        int type= previous_edge_type > 0 ? 0 : 1;
                                                        types.push_back(type);
                                                        previous_edge_type=type;
						}
					}
				}else{
					if((n1== n3 && n2 == n4) || ( n1==n4 && n2 == n3)){
						if(nodes.size()==0){
							nodes.push_back(n1);
						}else{
							if(n1==nodes.back()){
								nodes.push_back(n2);
							}else{
								nodes.push_back(n1);
							}
						}
					}else if(n1== n3 || n1 ==n4){
						nodes.push_back(n1);
					}else{
						nodes.push_back(n2);
					}
	                                types.push_back(gp->genome_index_to_type.find(*it)->second);
                                        previous_edge_type=gp->genome_index_to_type.find(*it)->second;

				}
				previous_edge= *it;
			}
		}
//		for(list<int>::iterator it=nodes.begin(); it !=nodes.end(); it++){
//			cout<<"nodes:";
//			cout<<*it<<"\t";
//		}
//		cout<<endl;
		if(return_nodes)
			return nodes;
		else
			return types;
	}

	bool is_cycle(list<int> SV_path){
		list<int> nodes= SV_path_to_nodes(SV_path,true);
                        cout<<"path_in_is_cycle:";
                        for(list<int>::iterator it =nodes.begin(); it!=nodes.end(); it++){
                                cout<<*it<<"\t";
                        }
		if((gp->genome_index_to_type.find(SV_path.front())->second == 0 && gp->genome_index_to_type.find(SV_path.back())->second == 2) || (gp->genome_index_to_type.find(SV_path.back())->second == 0 && gp->genome_index_to_type.find(SV_path.front())->second == 2)){
			if(nodes.front()==nodes.back()){
				cout<<"is_cycle"<<endl;
				return true;
			}
		}

		return false;
	}

	bool this_node_is_path_end(int src_or_dest, list<int> path, bool is_front_end){
	//	if(path.size()<=1){
	//		return true;
//		}else{
			list<int> nodes = SV_path_to_nodes(path, true);
			cout<<"path:";
			for(list<int>::iterator it =nodes.begin(); it!=nodes.end(); it++){
				cout<<*it<<"\t";
			}		
			cout<<endl;
			if(is_front_end){
				if(nodes.front()==src_or_dest){
                                       cout<<"this_node_is_path_end:"<<src_or_dest<<"\t"<<"true"<<endl;
					return true;
				}
				else{
                                       cout<<"this_node_is_path_end:"<<src_or_dest<<"\t"<<"false"<<endl;
					return false;
				}
			}else{
                               if(nodes.back()==src_or_dest){
                                       cout<<"this_node_is_path_end:"<<src_or_dest<<"\t"<<"true"<<endl;
                                        return true;
				}

                                else{
                                       cout<<"this_node_is_path_end:"<<src_or_dest<<"\t"<<"false"<<endl;
                                        return false;
				}
			}
			
	//	}
	}

        void addEdges(list<full_edges> full_edges_list)
        {

		cout<<"addEdges:"<<full_edges_list.front().dest<<endl;
	    int newNode_i=level_index[level];


	    if(newNode_i==0){
		    list<full_edges>::iterator it=  full_edges_list.begin();

		    int ***wtis;
		    wtis = new int**[2];
		    wtis[0] = new int*[(*it).degree1/2];
		    wtis[1] = new int*[(*it).degree2/2];
		    for(int i=0; i<(*it).degree1/2; i++){
			wtis[0][i] = new int[5];
			for(int j=0; j<5; j++)
				wtis[0][i][j]=-1;
		    }
		   for(int i=0; i<(*it).degree2/2; i++){
			wtis[1][i] = new int[5];
			for(int j=0; j<5; j++)
				wtis[1][i][j]=-1;
		    }

		     for(int i=0; i<(*it).multiplicity;i++){
			wtis[0][i][0]=i;
			wtis[0][i][1]=0;
			wtis[0][i][2]=0;
                        wtis[0][i][3]=0;
                        wtis[0][i][4]=0;


		     }
                     for(int i=0; i<(*it).multiplicity;i++){
                        wtis[1][i][0]=0;
                        wtis[1][i][1]=0;
                        wtis[1][i][2]=0;
                        wtis[1][i][3]=0;
                        wtis[0][i][4]=0;
		     }
		     (*it).wtis = wtis;
                     AdjListNode* newNode = newAdjListNode(0,0, full_edges_list, wtis);
                     nodes[0]=newNode;
                     newNode->next = NULL;
                     array[0].head = newNode;
                     newNode_i++;
	    }

	    int upto=level_index[level];

	     shared_ptr<int> p = make_shared<int>(newNode_i);
	if(thread_n ==1){
           cout<<"check"<< (level-1 <0 ? 0 : level_index[level-1]) <<"\t"<<upto<<endl;

		thread th(&DAG::addEdges_parallel,this, full_edges_list,  level-1 <0 ? 0 : level_index[level-1],  upto,  p);
		th.join();
	//	addEdges_parallel( src,  dest,  degree1,  degree2,  type,  edge_index,  edge_by_genome_index,  multiplicity,  level-1 <0 ? 0 : level_index[level-1],  upto,  &newNode_i);
	}else{
/*
	   int parent_start= level-1 <0 ? 0 : level_index[level-1];
	   int total_n= (upto - parent_start);
	   int n_for_thread= total_n /(thread_n) ;
	   int remain_n = total_n % thread_n;

	   if(n_for_thread == 0){
                thread th(&DAG::addEdges_parallel,this,  src,  dest,  degree1,  degree2,  type,  edge_index,  edge_by_genome_index,  multiplicity,  level-1 <0 ? 0 : level_index[level-1],  upto, p);
                th.join();
	   }else{
		   cout<<"here"<<n_for_thread<<"\t"<<remain_n<<endl;
                   int previous_start=parent_start;
                   int previous_end;

                   for(int thread_i  = 0 ; thread_i < thread_n ; thread_i++){
                        if(remain_n > 0){
                                previous_end=previous_start + n_for_thread + 1;
				remain_n--;
                        }else{
                                previous_end=previous_start + n_for_thread ;
                        }
			cout<<"here"<<parent_start<<"\t"<<upto<<"\t"<<previous_start<<"\t"<<previous_end<<endl;
                        previous_start = previous_end;

                   }
		   remain_n= total_n % thread_n;
		   thread threads[thread_n];
		   previous_start=parent_start;
		   for(int thread_i  = 0 ; thread_i < thread_n ; thread_i++){
			if(remain_n > 0){
				previous_end=previous_start + n_for_thread + 1;
				remain_n--;
			}else{
                                previous_end=previous_start + n_for_thread ;
			}
			threads[thread_i]=thread(&DAG::addEdges_parallel, this, src, dest, degree1,  degree2,  type,  edge_index,  edge_by_genome_index,  multiplicity, previous_start, previous_end, p);
			previous_start = previous_end; 
		   }
	
		   for(int thread_i = 0;  thread_i <thread_n ; thread_i++){
			threads[thread_i].join();
	          }
           }
*/
	}
	   cout<<"remove"<<(level-1 <0 ? 0 : level_index[level-1])<<"\t"<<upto<<endl;
	    for(int i=level-1 <0 ? 0 : level_index[level-1]; i<upto; i++){
		if(nodes[i]!=NULL){
			for(list<list<int>>::iterator jt=nodes[i]->gnodes.begin(); jt!=nodes[i]->gnodes.end(); jt++){
				(*jt).clear();
			}
			nodes[i]->gnodes.clear();
                        for(list<list<int>>::iterator jt=nodes[i]->frontier_types.begin(); jt!=nodes[i]->frontier_types.end(); jt++){
                                (*jt).clear();
                        }
                        nodes[i]->frontier_types.clear();
                        for(list<list<int>>::iterator jt=nodes[i]->epaths.begin(); jt!=nodes[i]->epaths.end(); jt++){
                                (*jt).clear();
                        }
                        nodes[i]->epaths.clear();
			nodes[i]->SVs_Ns.clear();
			nodes[i]->SVs_CNs.clear();

			nodes[i]=NULL;
		}
	    }

            newNode_i = *p;
            if(level > thres_level){
			thres_level+=4;
	    }
            level++;

            level_index[level]=newNode_i;

            cout<<"newNode_i:"<<newNode_i<<endl;
//       if(newNode_i==528){
  //             level_index[level-1]=527;
   //       }


	   for(list<full_edges>::iterator it = full_edges_list.begin(); it!=full_edges_list.end(); it++){
		    if(gp->genome_index_to_type.find((*it).edge_by_genome_index)->second==0){
			    degree_table[(*it).src][1]-=(*it).multiplicity;
			    degree_table[(*it).dest][1]-=(*it).multiplicity;
		    }else{
			    degree_table[(*it).src][2]-=(*it).multiplicity;
			    degree_table[(*it).dest][2]-=(*it).multiplicity;
		   }
	   }
	}

	void addEdges_parallel(list<full_edges> full_edges_list, int parent_start, int parent_end, shared_ptr<int> newNode_i){

            for(int i = parent_start ; i<parent_end; i++){
		if(nodes[i]==NULL)
			continue;

		list<full_edges_add> full_edges_add_list;
		int** wtis_src;
		int wtis_loop [2][128][4];
		int wtis_loop_i= 0;

		int total_multi=0;
		int total_degree1=0;
		int total_degree2=0;
		for(list<full_edges>::iterator flt = full_edges_list.begin(); flt!=full_edges_list.end(); flt++){
			int type = gp->genome_index_to_type.find((*flt).edge_by_genome_index)->second;		

			int ***wtis;
			wtis = new int**[2];
			wtis[0] = new int*[(*flt).degree1];
			wtis[1] = new int*[(*flt).degree2];
			for(int i=0; i<(*flt).degree1; i++){
			wtis[0][i] = new int[5];
			}
			for(int i=0; i<(*flt).degree2; i++){
			wtis[1][i] = new int[5];
			}

			int src = (*flt).src;
			int dest = (*flt).dest;
			int degree1 = (*flt).degree1;
			int degree2 = (*flt).degree2;
			list<int> edges_index = (*flt).edges_index;
			int edge_by_genome_index = (*flt).edge_by_genome_index;
			int multiplicity = (*flt).multiplicity;
			total_multi=total_multi+multiplicity;
			total_degree1=total_degree1 + degree1;
			total_degree2=total_degree2 + degree2;
			int wtis_i=0;
			int wtis_i_dest_adj=0;
			int wtis_i_src_adj=0;
			int it_i=0;
			int src_i=0;
			int dest_i=0;
			list<list<int>> l_SV_paths_src;
			list<list<int>> l_SV_paths_dest;
			int degree1_com;
			int degree2_com;

			if(type==0){
				degree1_com=degree_table[(*flt).src][1];
				degree2_com=degree_table[(*flt).dest][1];
			}else{
				degree1_com=degree_table[(*flt).src][2];
				degree2_com=degree_table[(*flt).dest][2];
			}
			cout<<"dest"<<degree1_com<<"\t"<<degree2_com<<endl;
			for(int i=0; i<degree1_com; i++){
				for(int j=0; j<4; j++)
					wtis[0][i][j]=-1;
				wtis[0][i][4]=0;
			}
			for(int i=0; i<degree2_com; i++){
				for(int j=0; j<4; j++)
					wtis[1][i][j]=-1;
				wtis[1][i][4]=0;
			}
			list<list<int>>::iterator ft=nodes[i]->frontier_types.begin();
			list<int>::iterator nt=nodes[i]->SVs_Ns.begin();
			list<int>::iterator cnt=nodes[i]->SVs_CNs.begin();
			for(list<list<int>>::iterator it=nodes[i]->gnodes.begin(); it !=nodes[i]->gnodes.end();it++){
				
				if(
					((src == (*it).front() && dest == (*it).back()) || (src == (*it).back() && dest == (*it).front())) 
					&& (((*ft).front() == 0 && type !=0)|| ((*ft).front() != 0 && type ==0)) 
					&& (((*ft).back() == 0 && type !=0)|| ((*ft).back() != 0 && type ==0)))	{
                                                wtis_loop[0][wtis_loop_i][0]=it_i;
                                                wtis_loop[0][wtis_loop_i][1]=0;
                                                wtis_loop[0][wtis_loop_i][2]= src == (*it).front() ? dest : src;
                                                wtis_loop[0][wtis_loop_i][3]=-1; /// loop
                                                wtis_loop[1][wtis_loop_i][0]=it_i;
                                                wtis_loop[1][wtis_loop_i][1]=1;
                                                wtis_loop[1][wtis_loop_i][2]= dest == (*it).back() ? src : dest;
                                                wtis_loop[1][wtis_loop_i][3]=-1; //// loop
                                              wtis_loop_i++;
                                                degree1_com--;
                                                degree2_com--;
                                                multiplicity--;
                                               it_i++;
                                                ft++;
                                                nt++;
                                                cnt++;
                                                continue;
				}



				if(( src == (*it).front() || dest == (*it).front())  && (((*ft).front() == 0 && type !=0)|| ((*ft).front() != 0 && type ==0) )  ){
	//				cout<<"destche"<<endl;
					int dest_adj=1;
					if(src==(*it).front()){
						dest_adj=0;
					}
					wtis_i = dest_adj ?  wtis_i_dest_adj : wtis_i_src_adj;

					wtis[dest_adj][wtis_i][0]=it_i;
					wtis[dest_adj][wtis_i][1]=0;
					wtis[dest_adj][wtis_i][2]=-1;
					wtis[dest_adj][wtis_i][4]=(*nt) > 0 ? round((double)(*cnt)/(*nt)) : 0;
	//				wtis[dest_adj][wtis_i][2]= dest_adj ? src : dest;

					int i_exist = dest_adj? SV_path_in_SV_paths(l_SV_paths_dest, *it) : SV_path_in_SV_paths(l_SV_paths_src, *it);
					if(i_exist>=0){
						wtis[dest_adj][wtis_i][3]=i_exist;

					}else{
						wtis[dest_adj][wtis_i][3]= dest_adj? dest_i : src_i;
						dest_adj? dest_i ++ : src_i++;
						dest_adj? l_SV_paths_dest.push_back(*it) : l_SV_paths_src.push_back(*it);
					}
					dest_adj? wtis_i_dest_adj++ : wtis_i_src_adj++;
		//			cout<<"destche"<<endl;
				}
				if((src == (*it).back() || dest == (*it).back())  && (((*ft).back() == 0 && type !=0)|| ((*ft).back() != 0 && type ==0) ) ){
 //                                       cout<<"destche"<<endl;

					int dest_adj=1;
					if(src==(*it).back()){
						dest_adj=0;
					}
					wtis_i = dest_adj ?  wtis_i_dest_adj : wtis_i_src_adj;

					wtis[dest_adj][wtis_i][0]=it_i;
					wtis[dest_adj][wtis_i][1]=1;
                                        wtis[dest_adj][wtis_i][2]=-1;
					wtis[dest_adj][wtis_i][4]=(*nt) > 0? round((double)(*cnt)/(*nt)) : 0;
					//wtis[dest_adj][wtis_i][2]= dest_adj ? src : dest;

					int i_exist = dest_adj? SV_path_in_SV_paths(l_SV_paths_dest, *it) : SV_path_in_SV_paths(l_SV_paths_src, *it);
					if(i_exist>=0){
						wtis[dest_adj][wtis_i][3]=i_exist;
					}else{
						wtis[dest_adj][wtis_i][3]= dest_adj? dest_i : src_i;
						dest_adj? dest_i ++ : src_i++;
						dest_adj? l_SV_paths_dest.push_back(*it) : l_SV_paths_src.push_back(*it);
					}
					dest_adj? wtis_i_dest_adj++ : wtis_i_src_adj++;
 //                                       cout<<"destche"<<endl;

				}
				it_i++;
				ft++;
				nt++;
				cnt++;
                	}
			for(int i=0; i<degree1_com; i++){
				if(wtis[0][i][3]<0){
					wtis[0][i][3]=src_i;
				}
			}
			for(int i=0; i<degree2_com; i++){
				if(wtis[1][i][3]<0){
					wtis[1][i][3]=dest_i;
				}
			}

			for(int i=0; i<degree1_com; i++){
				cout<<wtis[0][i][0]<<"\t"<<wtis[0][i][1]<<"\t"<<wtis[0][i][2]<<"\t"<<wtis[0][i][3]<<"\t"<<wtis[0][i][4]<<endl;
			}
			cout<<endl;
			for(int i=0; i<degree2_com; i++){
				cout<<wtis[1][i][0]<<"\t"<<wtis[1][i][1]<<"\t"<<wtis[1][i][2]<<"\t"<<wtis[1][i][3]<<"\t"<<wtis[1][i][4]<<endl;
			}
			cout<<endl;
			if(flt == full_edges_list.begin())
				wtis_src = wtis[0];
			else{
				for(int i=0; i<(*flt).degree1; i++){
					delete [] wtis[0][i];
				}
				delete [] wtis[0];
			}

			full_edges_add edges_add;
			edges_add.src = src;
			edges_add.dest = dest;
			edges_add.degree1_com=degree1_com;
                        edges_add.degree2_com=degree2_com;
			edges_add.multiplicity=multiplicity;
			edges_add.wtis_dest=wtis[1];
			full_edges_add_list.push_back(edges_add);
		}
	//	cout<<"whu_before"<<endl;
		list<int***> tcs = total_combinations_full (full_edges_add_list , wtis_src, full_edges_add_list.front().degree1_com, total_multi);
              //  list<int***> tcs = set_cover (full_edges_add_list , wtis_src, full_edges_add_list.front().degree1_com, branch_limit*10);

		cout<<"whu"<<tcs.size()<<endl;
		int tcs_branch=0;
//		int total_multi=0;
//	        for(list<full_edges_add>::iterator it = full_edges_add_list.begin(); it!=full_edges_add_list.end(); it++){
  //     	         	total_multi+=(*it).multiplicity;
    //  	 	 }

                for(list<int***>::iterator tcs_i= tcs.begin(); tcs_i!=tcs.end(); tcs_i++){
			mtx.lock();
			int l_newNode_i=*newNode_i;
			(*newNode_i)++;
			mtx.unlock();
                      clock_t tStart = clock();
			for( int i = 1; i <= wtis_loop_i; i++){
				for( int ii =0; ii <4; ii++){
					(*tcs_i)[0][total_multi  - i ][ii] = wtis_loop[0][i-1][ii];
					(*tcs_i)[1][total_multi - i ][ii] = wtis_loop[1][i-1][ii];
				}
			}
			
                        cout<<"tcs:";
                        for(int i=0; i<total_multi; i++){
                                cout<<(*tcs_i)[0][i][0]<<"\t"<<(*tcs_i)[0][i][1]<<"\t"<<(*tcs_i)[0][i][2]<<"\t"<<(*tcs_i)[0][i][3]<<endl;
                        }
                        cout<<endl;
                        for(int i=0; i<total_multi; i++){
                                cout<<(*tcs_i)[1][i][0]<<"\t"<<(*tcs_i)[1][i][1]<<"\t"<<(*tcs_i)[1][i][2]<<"\t"<<(*tcs_i)[1][i][3]<<endl;
                        }
                        cout<<endl;

               //         if( !tcs_validity(wtis, *tcs_i,  multiplicity,degree1_com, degree2_com)){
                        //if( !tcs_validity_error(wtis, *tcs_i,  multiplicity,degree1_com, degree2_com) || !tcs_validity(wtis, *tcs_i,  multiplicity,degree1_com, degree2_com)){
                 //               continue;
                   //     }
                        tcs_branch++;
                        if(tcs_branch > branch_limit){
                                break;
                        }
                       AdjListNode* newNode = newAdjListNode(i, l_newNode_i, full_edges_list, *tcs_i);
 //                       static int thres_level = 10;
                        bool connect_test_from;
                        cout<<"level"<<level<<"\t"<<"threslevel"<<thres_level<<"hidden_edge"<<hidden_edge<<endl;
                        if(level > thres_level || level == hidden_edge-1){
			//	connect_test_from= true;
				
                                connect_test_from = connect_test_from_SV_paths_speed(newNode);
                                //thres_level +=10;
                        }else{
				cout<<"always true?"<<endl;
                                connect_test_from = true;
                        }
                        if(connect_test_from){
                              int same_exist_index=-1;
                              int same_or_equivalent=0;
//                                cout<<"sum:"<<newNode->SV_paths_magic_number<<endl;
                              //  iter_pair = mm.equal_range(newNode->SV_paths_magic_number);
                                int node_i;

                                if( same_exist_index<0){
                                //        mm.insert(pair<int,int>(newNode->SV_paths_magic_number,newNode_i));
                                        nodes[l_newNode_i]=newNode;
                                     //    (*newNode_i)++;
                                }else{
					for(list<list<int>>::iterator jt=newNode->gnodes.begin(); jt!=newNode->gnodes.end(); jt++){
						(*jt).clear();
					}
					newNode->gnodes.clear();
					for(list<list<int>>::iterator jt=newNode->frontier_types.begin(); jt!=newNode->frontier_types.end(); jt++){
						(*jt).clear();
					}
					newNode->frontier_types.clear();
					for(list<list<int>>::iterator jt=newNode->epaths.begin(); jt!=newNode->epaths.end(); jt++){
						(*jt).clear();
					}
					newNode->epaths.clear();
					newNode->SVs_Ns.clear();
                                        newNode->SVs_CNs.clear();
					nodes[l_newNode_i]=NULL;
					delete newNode;
                                }

                        }else{
                                       for(list<list<int>>::iterator jt=newNode->gnodes.begin(); jt!=newNode->gnodes.end(); jt++){
                                                (*jt).clear();
                                        }
                                        newNode->gnodes.clear();
                                        for(list<list<int>>::iterator jt=newNode->frontier_types.begin(); jt!=newNode->frontier_types.end(); jt++){
                                                (*jt).clear();
                                        }
                                        newNode->frontier_types.clear();
                                        for(list<list<int>>::iterator jt=newNode->epaths.begin(); jt!=newNode->epaths.end(); jt++){
                                                (*jt).clear();
                                        }
                                        newNode->epaths.clear();
                                        newNode->SVs_Ns.clear();
                                        newNode->SVs_CNs.clear();
                                        nodes[l_newNode_i]=NULL;
                                        delete newNode;

                        }
                       printf("Time taken: %.8fs\t%d\n", (double)(clock() - tStart)/(CLOCKS_PER_SEC), *newNode_i);
/*			for(int i=0; i<2; i++){
				for(int j=0; j<total_multi; j++){
					delete[] (*tcs_i)[i][j];
				}
				delete[] (*tcs_i)[i];
			}
			delete[] (*tcs_i);
*/

              }
                for(list<int***>::iterator tcs_i= tcs.begin(); tcs_i!=tcs.end(); tcs_i++){
                       for(int i=0; i<2; i++){
                                for(int j=0; j<total_multi; j++){
                                        delete[] (*tcs_i)[i][j];
                                }
                                delete[] (*tcs_i)[i];
                        }
                        delete[] (*tcs_i);
		}
	      tcs.clear();
            
		for(list<full_edges_add>::iterator it = full_edges_add_list.begin(); it!=full_edges_add_list.end(); it++){
			for(int j=0; j<(*it).degree2_com; j++){
				delete[] (*it).wtis_dest[j];
			}
			delete[] (*it).wtis_dest;
		 }
	      }

	}

        bool tcs_validity_error(int ***wtis, int*** tcs_l, int multiplicity, int degree1_com, int degree2_com){
		
		vector<int> v1;
		vector<int> v2;
		int v1_count=0;
		int v2_count=0;
		for(int i=0; i<multiplicity; i++){
			if(tcs_l[0][i][0] >=0){
				v1.push_back(tcs_l[0][i][0]);
				v1_count++;
			}
		}

                for(int i=0; i<multiplicity; i++){
                        if(tcs_l[1][i][0] >=0){
                        	v2.push_back(tcs_l[1][i][0]);
				v2_count++;
			}
                }
		vector<int>::iterator new_end;
		sort(v1.begin(),v1.end());
		new_end=unique(v1.begin(),v1.end());
		v1.erase(new_end, v1.end());

		sort(v2.begin(),v2.end());
		new_end=unique(v2.begin(),v2.end());
                v2.erase(new_end, v2.end());

		for(vector<int>::iterator it=v1.begin(); it!=v1.end(); it++){
			cout<<*it<<"\t";
		}
		cout<<endl;
               for(vector<int>::iterator it=v2.begin(); it!=v2.end(); it++){
                        cout<<*it<<"\t";
                }
                cout<<endl;


		cout<<v1.size()<<"\t"<<v1_count<<"\t"<<v2.size()<<"\t"<<v2_count<<endl;
		if(v1.size() != v1_count || v2.size() != v2_count){
			cout<<"return false"<<endl;
			return false;
		}
		return true;
	}

	bool tcs_validity(int ***wtis, int*** tcs_l, int multiplicity, int degree1_com, int degree2_com){
		return true;

		//return true;
		int max_degree = degree1_com > degree2_com ? degree1_com : degree2_com;
                int how_many_index_in_tcs_l[max_degree];
                vector<int> tvs[max_degree];
                for(int i=0; i<max_degree; i++){
  	              tvs[i] = vector<int>(0);
                }
                for(int i=0; i<max_degree; i++){
                        how_many_index_in_tcs_l[i]=0;
                }

                for(int i=0; i<multiplicity; i++){
                        how_many_index_in_tcs_l[tcs_l[1][i][3]]++;
			tvs[tcs_l[1][i][3]].push_back(tcs_l[0][i][3]);
                }


		vector<int> arr;
		vector<int> max;

		for(int i=0; i<multiplicity; i++){
			vector<int>::iterator it = max.begin();
			vector<int>::iterator jt;
			for(jt = arr.begin(); jt!=arr.end(); jt++){
				if(*jt==tcs_l[0][i][3]){
					(*it)++;
					break;
				}
				it++;
			}
			if(jt == arr.end()){
				arr.push_back(tcs_l[0][i][3]);
				max.push_back(1);
			}
		}
		
		for(int i=0; i<max_degree; i++){
			if(how_many_index_in_tcs_l[i]!=0){
				for(int j=0; j<tvs[i].size(); j++){
					int m_i=0;
					for(int k=0; k<arr.size(); k++){
						if(arr.at(k)==tvs[i].at(j)){
							m_i=k;
							break;
						}
					}
				//	max.at(tvs[i].at(j))--;
					max.at(m_i)--;
				}	
				cout<<"tcs"<<endl;
				if(tv_validity(tvs[i], how_many_index_in_tcs_l[i], arr,max)){

				}else{
					return false;
				}
                              for(int j=0; j<tvs[i].size(); j++){
                                        int m_i=0;
                                        for(int k=0; k<arr.size(); k++){
                                                if(arr.at(k)==tvs[i].at(j)){
                                                        m_i=k;
                                                        break;
                                                }
                                        }
                                //      max.at(tvs[i].at(j))--;
                                        max.at(m_i)--;
				}
			}
		}

		return true;

/*
                for(int i=0; i<multiplicity; i++){
			how_many_index_in_tcs_l[i]=0;
		}
                for(int i=0; i<multiplicity; i++){
			how_many_index_in_tcs_l[tcs_l[0][i][3]]++;
		}

                int how_many_index_in_wtis[degree1_com];
               for(int i=0; i<degree1_com; i++){
			how_many_index_in_wtis[i]=0;
		}
               for(int i=0; i<degree1_com; i++){
                        how_many_index_in_wtis[wtis[0][i][3]]++;
                }

	        int how_many_non_zero=0;
		int how_many_max=0;
		for(int i=0; i<multiplicity; i++){
                        if(how_many_index_in_tcs_l[i]!=0){
				how_many_non_zero++;
				if(how_many_index_in_tcs_l[i] == how_many_index_in_wtis[i]){
					how_many_max++;
				}
			}
		}
		cout<<"tcs_validity1"<<how_many_non_zero<<"\t"<<how_many_max<<endl;
		if(how_many_non_zero == 1){
			
		}else if(how_many_non_zero > 1 && how_many_max >= how_many_non_zero -1 ){

		}else{
			return false;
		}


                for(int i=0; i<multiplicity; i++){
                        how_many_index_in_tcs_l[i]=0;
                }
                for(int i=0; i<multiplicity; i++){
                        how_many_index_in_tcs_l[tcs_l[1][i][3]]++;
                }
                int how_many_index_in_wtis2[degree2_com];

               for(int i=0; i<degree2_com; i++){
                        how_many_index_in_wtis2[i]=0;
                }
               for(int i=0; i<degree2_com; i++){
                        how_many_index_in_wtis2[wtis[1][i][3]]++;
                }

               how_many_non_zero=0;
                how_many_max=0;
                for(int i=0; i<multiplicity; i++){
                        if(how_many_index_in_tcs_l[i]!=0){
                                how_many_non_zero++;
                                if(how_many_index_in_tcs_l[i] == how_many_index_in_wtis2[i]){
                                        how_many_max++;
                                }
                        }
                }
                cout<<"tcs_validity2"<<how_many_non_zero<<"\t"<<how_many_max<<endl;
                if(how_many_non_zero == 1){
                }
                else if(how_many_non_zero > 1 && how_many_max >= how_many_non_zero -1 ){
                }else{
                        return false;
                }
		return true;
*/
	}

	void out_leafs(int target_node){
		ofstream outfile;
		outfile.open("euler_paths."+to_string(target_node));
		for(int i = level_index[level-1]; i < level_index[level]; i++){
			if(nodes[i]==NULL)
				continue;
			outfile<<"gnodes"<<i-level_index[level-1]+1<<endl;
			for(list<list<int>>::iterator it = nodes[i]->gnodes.begin(); it != nodes[i]->gnodes.end();it++){
				for(list<int>::iterator jt = (*it).begin(); jt!=(*it).end(); jt++)
					outfile<<*jt<<"\t";
				outfile<<endl;
			}
			outfile<<endl;
		}


	}

	int SV_path_in_SV_paths(list<list<int>> l_SV_paths, list<int> l_SV_path){
	//			cout<<"SV_PATH_IN_SV_PATHS"<<endl;
				int i=0;
                                for(list<list<int>>::iterator jt = l_SV_paths.begin(); jt!=l_SV_paths.end();jt++){
					cout<<endl;
					list<int>::iterator lt = l_SV_path.begin();
					if(l_SV_path.size() == (*jt).size()){
	                                        for(list<int>::iterator kt=(*jt).begin(); kt!=(*jt).end();kt++){
	                                                if((*kt)!=*lt)
								break;
							lt++;
	                                       	}
					}
					if(lt == l_SV_path.end()){
						//cout<<"right"<<i<<endl;
						return i;
					}
					i++;
                                }
				return -1;
	}

	bool connect_test_from_SV_paths_speed(AdjListNode* newNode){
		
		return true;//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////edit/////////////

		cout<<"connect_test_from_SV_paths_speed"<<endl;
		int check=0;
		int src[2];
		int dest[2];
		src[0]=0;
		dest[0]=19;
		src[1]=20;
		dest[1]=73;
		src[2]=74;
		dest[2]=101;
		for(int n=0; n<3; n++){
	                list<vector<int>> SV_adj[gp->V];
			list<list<int>>::iterator jt = newNode->frontier_types.begin();
			for(list<list<int>>::iterator it = newNode->gnodes.begin(); it!=newNode->gnodes.end(); it++){
				if((*it).size()>=2){
					cout<<(*it).front()<<"\t"<<(*it).back()<<"\t"<<(*jt).front()<<"\t"<<(*jt).back()<<endl;	
					SV_adj[(*it).front()].push_back({(*it).back(), (*jt).back(),1});
		                        SV_adj[(*it).back()].push_back({(*it).front(), (*jt).front(),1});
				}
				jt++;
			}
			if(gp->isConnected_fromSV_src_to_dest(gp->E, newNode->gnodes, newNode->epaths, SV_adj, src[n], dest[n])){
				check++;
			}
		}
		if(check ==3)
			return true;
		return false;
		
	}


	
	bool SV_paths_are_equivalent(list<list<int>>* paths1, int paths1_n, list<list<int>>* paths2, int paths2_n){
                if(paths1_n == paths2_n){
			list<list<int>>::reverse_iterator jt=paths2->rbegin();
			for(list<list<int>>::reverse_iterator it=paths1->rbegin(); it!=paths1->rend(); it++){
				if((*it).front() == (*jt).front() && (*it).back() == (*jt).back()){

				}else{
					return false;
				}
				jt++;
                      }

                }else{
                        return false;
                }
		return true;
	}

	bool single_fragment_same(list<int> fragment1, list<int> fragment2){
		list<int>::iterator it2=fragment2.begin();
	       for(list<int>::iterator it1=fragment1.begin(); it1!=fragment1.end();it1++){
			if(*it1 != *it2){
				break;
			}else{
				it2++;
			}
		}
		if(it2 == fragment2.end()){
			return true;
		}
		it2=fragment2.begin();
		for(list<int>::reverse_iterator it1=fragment1.rbegin(); it1!=fragment1.rend();it1++){
			if(*it1 != *it2){
				break;
			}else{
				it2++;
			}
		}
		if(it2 == fragment2.end()){
			return true;
		}
		return false;
	}

	bool paths_are_same(list<int>* paths1, int paths1_n, list<int>* paths2, int paths2_n){
		list<int> searched;
		if(paths1_n == paths2_n){
			for(int i=0; i<paths1_n;i++){
                                int where=-1;
                                for(int j=0; j<paths2_n;j++){
                                        if(!exist_in_list(j,searched)){
						if(paths1[i].size() == paths2[j].size() &&  fragment_is_same(paths1[i], paths2[j])){
							where=j;
							break;
						}
                                        }
                                }
                                if(where>=0){
                                        searched.push_front(where);
                                }else{
                                        return false;
                                }
                        }
                        return true;

                }else{
                        return false;
                }
        }

	bool fragment_is_same(list<int> fragment1, list<int> fragment2){
                list<int>* l_fragment1;
                list<int>* l_fragment2;
		l_fragment1= new list<int> [fragment1.size()];
                l_fragment2= new list<int> [fragment2.size()];
		int l_fragment1_size=1;
		int l_fragment2_size=1;
		
	        int previous_is_hidden=0;
		int its_initial=1;
		for(list<int>::iterator it=fragment1.begin(); it!=fragment1.end();it++){
			if(edge_by_genome[*it]!=hidden_edge){
				l_fragment1[l_fragment1_size-1].push_back(edge_by_genome[*it]);
				previous_is_hidden=0;
			}else{
				if(l_fragment1_size==1 && its_initial){

				}else if(previous_is_hidden == 1){
	
				}else{
					l_fragment1_size++;
				}
                                previous_is_hidden=1;
			}
			its_initial=0;
		}
		if(previous_is_hidden==0){
			l_fragment1_size++;
		}
		previous_is_hidden=0;
		its_initial=1;
                for(list<int>::iterator it=fragment2.begin(); it!=fragment2.end();it++){
                        if(edge_by_genome[*it]!=hidden_edge){
                                l_fragment2[l_fragment2_size-1].push_back(edge_by_genome[*it]);
                                previous_is_hidden=0;
                        }else{
                                if(l_fragment2_size==1 && its_initial){

                                }else if(previous_is_hidden == 1){

                                }else{
                                        l_fragment2_size++;
                                }
                                previous_is_hidden=1;
                        }
			its_initial=0;
                }
		if(previous_is_hidden==0){
			l_fragment2_size++;
		}
	/*	
                for(int i=0; i<l_fragment1_size-1; i++){
			cout<<"fragment:";
			for(list<int>::iterator it=l_fragment1[i].begin(); it!=l_fragment1[i].end();it++){
				cout<<*it<<"\t";
			}
			cout<<endl;
		}
		cout<<"here"<<endl;
                for(int i=0; i<l_fragment2_size-1; i++){
                        cout<<"fragment:";
                        for(list<int>::iterator it=l_fragment2[i].begin(); it!=l_fragment2[i].end();it++){
                                cout<<*it<<"\t";
                        }
                        cout<<endl;
                }
	*/



                list<int> searched;
		if(l_fragment1_size == l_fragment2_size){
			for(int i=0; i<l_fragment1_size-1; i++){
				int where=-1;
				for(int j=0;j<l_fragment2_size-1;j++){
					if(!exist_in_list(j,searched)){
						if(l_fragment1[i].size() == l_fragment2[j].size()){
							list<int>::iterator it2=l_fragment2[j].begin();
							for(list<int>::iterator it1=l_fragment1[i].begin(); it1!=l_fragment1[i].end();it1++){
								if(*it1 != *it2){
									break;
								}else{
									it2++;
								}
							}
							if(it2 == l_fragment2[j].end()){
								where=j;
								break;
							}
							it2=l_fragment2[j].begin();
                                                        for(list<int>::reverse_iterator it1=l_fragment1[i].rbegin(); it1!=l_fragment1[i].rend();it1++){
                                                                if(*it1 != *it2){
                                                                        break;
                                                                }else{
                                                                        it2++;
                                                                }
                                                        }
                                                        if(it2 == l_fragment2[j].end()){
                                                                where=j;
                                                                break;
                                                        }

						}
					}
				}
				if(where>=0){
					searched.push_front(where);
				}else{
					return false;
				}
			}
			return true;
		}else{
			return false;
		}

	}
	bool iterator_exist_in_list(list<list<int>>::iterator it, list<list<list<int>>::iterator> searched){
		for(list<list<list<int>>::iterator>::iterator search=searched.begin(); search!=searched.end(); search++){
			if(*search==it){
				return true;
			}
		}
		return false;
	}

	bool exist_in_list(int j, list<int> searched){
		for(list<int>::iterator it=searched.begin(); it!=searched.end(); it++){
			if(j== *it){
				return true;
			}
		}
		return false;
	}




};
