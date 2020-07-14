#include "main.h"
//list<int> DAG_expand(int start_nodes_n);

// Driver program to test above function
vector<int> DAG_expand(bool is_pre_expand, int target_node, int final_node, int total_nodes, int total_edges, int hidden_node, int hidden_edge){
//vector<int> DAG_expand(bool is_pre_expand, int target_node, int final_node)
  // Let us first create and test graphs shown in above figure
// printf("please check unknown node intiger value\n");a

       //  total_nodes=139;
         //total_edges=914;
       //  hidden_node=138;
       //  hidden_edge=138;

	int** degree_table=new int* [total_nodes];
	for(int i=0; i<total_nodes; i++){
		degree_table[i]=new int [3];
	}
	int** edge_table=new int* [total_edges];
        for(int i=0; i<total_edges; i++){
                edge_table[i]=new int [7];
        }

        ifstream degrees("degrees");
	string line;
	int i=0;
        while(getline(degrees,line,'\n')){
                stringstream f_iss(line);
                string val;
		int j=0;
		while(getline(f_iss,val,'\t')){
                        stringstream convertor(val);
                        int t;
                        convertor>>t;
			degree_table[i][j]=t;
			j++;
		}
		i++;
	}
	
        Graph g(total_nodes,total_edges,hidden_node, edge_table);
        Graph* gp;
        gp = &g;

	ifstream edges("edge_information.txt");
       i=0;
        while(getline(edges,line,'\n')){
                stringstream f_iss(line);
                string val;
                int j=0;
		int* temp;
                temp = new int [7];
                while(getline(f_iss,val,'\t')){
                        stringstream convertor(val);
                        int t;
                        convertor>>t;
			temp[j]=t;
                        edge_table[i][j]=t;
                        j++;
                }
		gp->addEdge(temp[0],temp[1],temp[4], temp[5]);
		cout<<"add"<<temp[0]<<"\t"<<temp[1]<<"\t"<<temp[4]<<"\t"<<temp[5]<<endl;
		map<int,int>::iterator mit=gp->genome_index_to_type.begin();
		map<int,vector<int>>::iterator mvit = gp->genome_index_to_nodes.begin();
                map<vector<int>,int>::iterator mvit2= gp->nodes_to_genome_index.begin();
                multimap<vector<int>,int>::iterator mvit3= gp->nodes_to_edge_index.begin();
                unordered_map<int,int>::iterator mvit4= node_cantor_to_edge_index.begin();


		vector<int> tv;
		tv.push_back(temp[0]);
		tv.push_back(temp[1]);
		gp->genome_index_to_type.insert(mit, pair<int,int>(temp[6],temp[4]));
		gp->genome_index_to_nodes.insert(mvit, pair<int,vector<int>>(temp[6], tv));
                gp->nodes_to_genome_index.insert(mvit2, pair<vector<int>,int>(tv, temp[6]));
                gp->nodes_to_edge_index.insert(mvit3, pair<vector<int>,int>(tv, temp[5]));
		node_cantor_to_edge_index.insert(mvit4, pair<int,int>((temp[0]+1+temp[1]+1)*(temp[0]+temp[1]+1+1+1)/2+temp[1], temp[5]));
                i++;
        }

/*
	list<int> BFS_by_edges=g.BFS_except_hidden_by_edges(0);
	cout<<"BFS:";
	for(list<int>::iterator it=BFS_by_edges.begin();it !=BFS_by_edges.end();it++){
		cout<<*it<<"\n";
	}
	cout<<endl;
*/
        g.edge_adj_matrix_constr();
	//g.gtest();

        DAG g1(10000000,hidden_node,hidden_edge,10000000, degree_table, edge_table, total_nodes, total_edges, gp);
        for(list<vector<int>>::iterator it = gp->adj[4].begin(); it!=gp->adj[4].end(); it++){
		cout<<(*it).at(0)<<endl;
	}

	bool visited[total_edges];

	for(i=0; i<total_edges; i++){
		visited[i]=false;
		
	}

        vector<int> level_index_list;

	for(i=0; i<gp->V; i++){
		list<full_edges> full_edges_list_seg;	
		list<full_edges> full_edges_list_rsv;
		for(list<vector<int>>::iterator it = gp->adj[i].begin(); it!=gp->adj[i].end(); it++){
			if(!visited[(*it).at(2)]){
				if((edge_table[(*it).at(2)][4] > 2)){
					break;
				}
	//			cout<<"here"<<(*it).at(0)<<endl;
				list<full_edges> &ltmp = edge_table[(*it).at(2)][4] == 0 ? full_edges_list_seg : full_edges_list_rsv;

				list<full_edges>::iterator jt;
				for(jt= ltmp.begin(); jt!=ltmp.end(); jt++){
					if( i  == (*jt).src  && (*it).at(0) ==   (*jt).dest){
						break;
					}
				}
				if(jt==ltmp.end()){

					full_edges full_edges_e;	
					full_edges_e.src=i;
					full_edges_e.dest=(*it).at(0);
					full_edges_e.degree1=degree_table[i][1]*2;
					full_edges_e.degree2=degree_table[(*it).at(0)][1]*2;
					full_edges_e.edges_index.push_back((*it).at(2));
					full_edges_e.edge_by_genome_index=(edge_table[(*it).at(2)][6]);
					full_edges_e.type = edge_table[(*it).at(2)][4];

					full_edges_e.multiplicity=1;

					ltmp.push_back(full_edges_e);
				}else{
					(*jt).edges_index.push_back((*it).at(2));
					(*jt).multiplicity++;
				}
				visited[(*it).at(2)]=true;
			}
		}	

		cout<<"seg"<<endl;
		for(list<full_edges>::iterator zt = full_edges_list_seg.begin(); zt!=full_edges_list_seg.end(); zt++){
			cout<<(*zt).src<<"\t"<<(*zt).dest<<"\t"<<(*zt).degree1<<"\t"<<(*zt).degree2<<endl;
			for(list<int>::iterator xt=(*zt).edges_index.begin(); xt!=(*zt).edges_index.end(); xt++){
				cout<<(*xt)<<"\t";
			}
			cout<<endl;
			cout<<(*zt).edge_by_genome_index<<"\t"<<(*zt).multiplicity<<endl;
		}
		cout<<"segend"<<endl;
                cout<<"rsv"<<endl;
                for(list<full_edges>::iterator zt = full_edges_list_rsv.begin(); zt!=full_edges_list_rsv.end(); zt++){
                        cout<<(*zt).src<<"\t"<<(*zt).dest<<"\t"<<(*zt).degree1<<"\t"<<(*zt).degree2<<endl;
                        for(list<int>::iterator xt=(*zt).edges_index.begin(); xt!=(*zt).edges_index.end(); xt++){
                                cout<<(*xt)<<"\t";
                        }
                        cout<<endl;
                        cout<<(*zt).edge_by_genome_index<<"\t"<<(*zt).multiplicity<<endl;
                }
                cout<<"rsvend"<<endl;

		if(i%2){
			if(full_edges_list_seg.size()!=0){
				g1.addEdges(full_edges_list_seg);
			}
			if(full_edges_list_rsv.size()!=0){
				g1.addEdges(full_edges_list_rsv);	
			}
		}else{
                       if(full_edges_list_rsv.size()!=0)
                                g1.addEdges(full_edges_list_rsv);
                        if(full_edges_list_seg.size()!=0)
                                g1.addEdges(full_edges_list_seg);
		}
		if(is_pre_expand && g1.call_level() > 1){
			cout<<"g1_call"<<g1.call_level_index(g1.call_level())<<"\t"<< g1.call_level_index(g1.call_level()-1)<<endl;
			if(g1.call_level_index(g1.call_level()) - g1.call_level_index(g1.call_level()-1) > 1000){
				level_index_list.push_back(g1.call_level_index(g1.call_level()-1));
                                level_index_list.push_back(g1.call_level_index(g1.call_level()));

				return level_index_list;
			}
		}else{
			if(g1.call_level_index(g1.call_level()) == final_node){
				cout<<"finalcheck"<<g1.call_level()<<"\t"<<final_node<<"\t"<<target_node<<endl;
 			        g1.set_target(g1.call_level(), target_node, final_node);///////////////////// edit////////////
			}
		}
	//	cout<<"call_level"<<g1.call_level_index()<<endl;
	}
                g1.out_leafs(target_node);
		return level_index_list;
}
