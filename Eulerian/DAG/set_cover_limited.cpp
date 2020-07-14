#include <iostream>
#include <string.h>
#include <algorithm>
#include <list>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <iterator>
#include "combination_precut.h"
#include "set_cover_branch_limit.h"
#include <time.h>
//int branch_limit=10;
using namespace std;
double time_sum=0;
list<int***> set_cover(list<full_edges_add> full_edges_add_list, int** wtis_src, int wtis_src_n, int branch_limit){
	list<int***> scs;

        vector<int> arr1;
        vector<int> max1;
	vector<int> SVs_CNs1;

        int total_multi=0;
        for(list<full_edges_add>::iterator it = full_edges_add_list.begin(); it!=full_edges_add_list.end(); it++){
                total_multi+=(*it).multiplicity;
        }
        int *** scs_e;
        scs_e = new int**[2];
        scs_e[0] = new int*[total_multi];
        scs_e[1] = new int*[total_multi];
        for(int i=0; i<total_multi; i++){
                scs_e[0][i] = new int[5];
                scs_e[1][i] = new int[5];
                for(int j=0; j<5; j++){
                        scs_e[0][i][j]=-1;
                        scs_e[1][i][j]=-1;
                }
        }
	
        for(int i=0; i<wtis_src_n; i++){
                vector<int>::iterator it = max1.begin();
                vector<int>::iterator jt;
                for(jt = arr1.begin(); jt!=arr1.end(); jt++){
                        if(*jt==wtis_src[i][3]){
                                (*it)++;
                                break;
                        }
                        it++;
                }
                if(jt == arr1.end()){
                        arr1.push_back(wtis_src[i][3]);
                        max1.push_back(1);
			SVs_CNs1.push_back(wtis_src[i][4]);
                }
        }

	vector<vector<int>> arr2_list;
	vector<vector<int>> max2_list;
	vector<vector<int>> SVs_CNs2_list;

	vector<int> multi;
	vector<bool> wtis_src_used;
	for(int i=0; i<wtis_src_n; i++){
		wtis_src_used.push_back(false);
	}
	vector<vector<bool>> wtis_dest_used_vector;

	for(list<full_edges_add>::iterator it = full_edges_add_list.begin(); it != full_edges_add_list.end(); it++){
		vector<int> arr2;
		vector<int> max2;
		vector<int> SVs_CNs2;
		vector<bool> wtis_dest_used;
		for(int i=0; i<(*it).degree2_com; i++){
			vector<int>::iterator jt = max2.begin();
			vector<int>::iterator kt;
			for(kt = arr2.begin(); kt!=arr2.end(); kt++){
				if(*kt==(*it).wtis_dest[i][3]){
					(*jt)++;
					break;
				}
				jt++;
			}
			if(kt == arr2.end()){
				arr2.push_back((*it).wtis_dest[i][3]);
				max2.push_back(1);
				SVs_CNs2.push_back((*it).wtis_dest[i][4]);
			}
			wtis_dest_used.push_back(false);
		}
		arr2_list.push_back(arr2);
		max2_list.push_back(max2);
		SVs_CNs2_list.push_back(SVs_CNs2);
		multi.push_back((*it).multiplicity);
		wtis_dest_used_vector.push_back(wtis_dest_used);
	}
	
	int to_hidden_edges = (wtis_src_n - total_multi > 0? wtis_src_n - total_multi : 0);	

	int selected_per_round[total_multi+to_hidden_edges];
	for(int i=0; i<total_multi+to_hidden_edges; i++){
		selected_per_round[i]=0;
	}
	int round =0;
	vector<int> selected_per_round_local;
      // clock_t tStart = clock();
	
        greedy_selection_recursive(&scs, scs_e ,0, full_edges_add_list.begin(), full_edges_add_list.end(), arr1, SVs_CNs1, max1, arr2_list, SVs_CNs2_list, max2_list, wtis_src, wtis_src_n,multi, total_multi, total_multi, wtis_src_used, wtis_dest_used_vector, round, selected_per_round, selected_per_round_local, 0, 0, 0, branch_limit, to_hidden_edges, to_hidden_edges);
 //       greedy_selection_recursive(&scs, scs_e ,0, full_edges_add_list.begin(), full_edges_add_list.end(), arr1, max1, arr2_list, max2_list, wtis_src, wtis_src_n,multi, total_multi, total_multi, wtis_src_used, wtis_dest_used_vector, round, selected_per_round, selected_per_round_local, 0, 0, 0, branch_limit, 0);

       // printf("greedy Time taken: %.8fs\t%d\n", (double)(clock() - tStart)/(CLOCKS_PER_SEC));
	for(int i=0; i<total_multi; i++){
	        delete[] scs_e[0][i];
	        delete[] scs_e[1][i];
	}
	delete[] scs_e[0];
	delete[] scs_e[1];
	delete[] scs_e;
       // cout<<"setCoverTime:"<<time_sum<<endl;
	return scs;
	//cout<<"setCoverTime:"<<time_sum<<endl;
}

void greedy_selection_recursive( list<int***>* scs, int*** scs_e, int scs_e_start, list<full_edges_add>::iterator start, list<full_edges_add>::iterator end, vector<int> arr1, vector<int> SVs_CNs1, vector<int> max1, vector<vector<int>> arr2_list, vector<vector<int>> SVs_CNs2_list, vector<vector<int>> max2_list, int** wtis_src, int wtis_src_n, vector<int> multi, int total_multi_initial, int total_multi_after, vector<bool> wtis_src_used, vector<vector<bool>> wtis_dest_used_vector, int round, int selected_per_round[], vector<int> selected_per_round_local, int max1_start, int max2_list_start, int max2_list_max2_start, int branch_limit, int to_hidden_edges_initial, int to_hidden_edges){
	if(scs->size() > branch_limit){
		return;
	}
	//cout<<"recursive call"<<total_multi_after<<"\t"<<to_hidden_edges<<endl;
	if(total_multi_after <= 0 && to_hidden_edges == 0 ){
			int i=0;
			bool more_optimal=false;
			for(i=0; i<selected_per_round_local.size(); i++){
				if(selected_per_round[i] > selected_per_round_local[i])
					break;
                                if(selected_per_round[i] < selected_per_round_local[i]){
                                         more_optimal=true;
				}
			}
                        if(more_optimal){
                                for(i=0; i<selected_per_round_local.size(); i++){
                                       selected_per_round[i]=selected_per_round_local[i];
                                }
                        }


			if(i== selected_per_round_local.size()){
				scs_e_sort(scs_e, total_multi_initial);
                                if(more_optimal){
                                        scs->clear();
                                        scs_push_back(scs, scs_e, total_multi_initial);
                                }else{
                                        scs_push_back(scs, scs_e, total_multi_initial);
                               }

			}else{
				for(int i=0; i<total_multi_initial; i++){
					delete[] scs_e[0][i];
					delete[] scs_e[1][i];
				}
				delete[] scs_e[0];
				delete[] scs_e[1];
				delete[] scs_e;
			}

		return;
	}
       clock_t tStart = clock();

	list<pair<vector<int>, list<full_edges_add>::iterator>> greedy_selected;
	list<int> greedy_selected_SVs_CN_diff;
	//cout<<"max1"<<max1.size()<<endl;
	for(int i=0; i<max1.size(); i++){
		list<full_edges_add>::iterator fit = start;
		int dest_i;
		for(dest_i=0; dest_i < max2_list.size(); dest_i++){
			for(int j=0; j<max2_list[dest_i].size(); j++){
                                vector<int> t = {min(multi.at(dest_i), max1.at(i), max2_list[dest_i].at(j)), i, dest_i,  j};
				if(greedy_selected.size() == 0 || greedy_selected.back().first.at(0) == min(multi.at(dest_i), max1.at(i), max2_list[dest_i].at(j))){
					greedy_selected_push_back(&greedy_selected, &greedy_selected_SVs_CN_diff, pair<vector<int>, list<full_edges_add>::iterator>(t, fit), abs(SVs_CNs1[i]-SVs_CNs2_list[dest_i][j]));
//					greedy_selected.push_back(pair<vector<int>, list<full_edges_add>::iterator>(t, fit));
				}else if(greedy_selected.back().first.at(0) < min(multi.at(dest_i), max1.at(i), max2_list[dest_i].at(j))){
					greedy_selected.clear();
					greedy_selected_SVs_CN_diff.clear();
                                        greedy_selected_push_back(&greedy_selected, &greedy_selected_SVs_CN_diff, pair<vector<int>, list<full_edges_add>::iterator>(t, fit), abs(SVs_CNs1[i]-SVs_CNs2_list[dest_i][j]));
				}
			}
			fit++;
		}
                if(to_hidden_edges > 0){
                        if(greedy_selected.size() == 0 || greedy_selected.back().first.at(0) == min(to_hidden_edges,to_hidden_edges, max1.at(i))){
                                vector<int> t = {min(to_hidden_edges,to_hidden_edges, max1.at(i)), i, -1, -1};
				greedy_selected_push_back(&greedy_selected, &greedy_selected_SVs_CN_diff, pair<vector<int>, list<full_edges_add>::iterator>(t, fit), to_hidden_edges_initial);
                                //greedy_selected.push_back(pair<vector<int>, list<full_edges_add>::iterator>(t, fit));
                        }else if(greedy_selected.back().first.at(0) < min(to_hidden_edges,to_hidden_edges, max1.at(i))){
                                vector<int> t = {min(to_hidden_edges,to_hidden_edges, max1.at(i)), i, -1, -1};
                                greedy_selected.clear();
				greedy_selected_SVs_CN_diff.clear();
				greedy_selected_push_back(&greedy_selected, &greedy_selected_SVs_CN_diff, pair<vector<int>, list<full_edges_add>::iterator>(t, fit), to_hidden_edges_initial);
 //                               greedy_selected.push_back(pair<vector<int>, list<full_edges_add>::iterator>(t, fit));
                        }
                }
	}

	time_sum += (double)(clock()-tStart)/CLOCKS_PER_SEC;
 //      printf("greedy Time taken: %.8fs\t%d\n", (double)(clock() - tStart)/(CLOCKS_PER_SEC));

	//cout<<"greedy"<<greedy_selected.size()<<endl;
	for(list<pair<vector<int>, list<full_edges_add>::iterator>>::iterator it = greedy_selected.begin(); it!=greedy_selected.end(); it++){
                if((*it).first.at(2) == -1){
 //                       cout<<"here"<<endl;
                        max1.at((*it).first.at(1)) -= (*it).first.at(0);
                        selected_per_round_local.push_back((*it).first.at(0));
			int *** scs_e_fill;
			scs_e_fill = new int **[2];
			scs_e_fill[0] = new int*[total_multi_initial];
			scs_e_fill[1] = new int*[total_multi_initial];
			for(int i=0; i<total_multi_initial; i++){
			       scs_e_fill[0][i] = new int[5];
				for(int j=0; j<5; j++){
					scs_e_fill[0][i][j] = scs_e[0][i][j];
				}
			}
			for(int i=0; i<total_multi_initial; i++){
				scs_e_fill[1][i] = new int[5];
				for(int j=0; j<5; j++){
					scs_e_fill[1][i][j] = scs_e[1][i][j];
				}
			}

			greedy_selection_recursive( scs, scs_e_fill, scs_e_start, start, end,arr1,SVs_CNs1,  max1, arr2_list,SVs_CNs2_list,  max2_list,  wtis_src, wtis_src_n,  multi, total_multi_initial, total_multi_after,wtis_src_used, wtis_dest_used_vector, round+1, selected_per_round, selected_per_round_local, (*it).first.at(1), (*it).first.at(2), (*it).first.at(3),branch_limit, to_hidden_edges_initial, to_hidden_edges-(*it).first.at(0));
			selected_per_round_local.pop_back();

			if(total_multi_after>0){
				for(int i=0; i<total_multi_initial; i++){
					delete[] scs_e_fill[0][i];
					delete[] scs_e_fill[1][i];
				}
				delete[] scs_e_fill[0];
				delete[] scs_e_fill[1];
				delete[] scs_e_fill;
			}

			max1.at((*it).first.at(1)) += (*it).first.at(0);
			continue;
		}

		max1.at((*it).first.at(1)) -= (*it).first.at(0);
		max2_list.at((*it).first.at(2)).at((*it).first.at(3)) -= (*it).first.at(0);
		multi.at((*it).first.at(2)) -= (*it).first.at(0);
		total_multi_after -= (*it).first.at(0);
	
	        int *** scs_e_fill;
                scs_e_fill = new int **[2];
                scs_e_fill[0] = new int*[total_multi_initial];
                scs_e_fill[1] = new int*[total_multi_initial];
                for(int i=0; i<total_multi_initial; i++){
                       scs_e_fill[0][i] = new int[5];
			for(int j=0; j<5; j++){
				scs_e_fill[0][i][j] = scs_e[0][i][j];
			}
                }
                for(int i=0; i<total_multi_initial; i++){
                        scs_e_fill[1][i] = new int[5];
                        for(int j=0; j<5; j++){
                                scs_e_fill[1][i][j] = scs_e[1][i][j];
			}
                }
		

		int fill_start = scs_e_start;
		vector<bool> wtis_src_used_l = wtis_src_used;
		vector<vector<bool>> wtis_dest_used_vector_l = wtis_dest_used_vector;
		for(int i=0; i<wtis_src_n && fill_start < scs_e_start + (*it).first.at(0) ; i++){
	//		cout<<"scheck"<<wtis_src[i][0]<<"\t"<< (scs_e_start + (*it).first.at(0))<<endl;
			if(wtis_src[i][3] == arr1.at((*it).first.at(1)) && !wtis_src_used_l.at(i)){
				for(int j=0; j<5; j++){
					scs_e_fill[0][fill_start][j]=wtis_src[i][j];
				}
				scs_e_fill[0][fill_start][2]=(*(*it).second).dest;
				wtis_src_used_l.at(i)=true;
				fill_start++;
			}
		}
		fill_start = scs_e_start;
                for(int i=0; i<(*(*it).second).degree2_com && fill_start < scs_e_start + (*it).first.at(0) ; i++){
                        if((*(*it).second).wtis_dest[i][3] == arr2_list.at((*it).first.at(2)).at((*it).first.at(3)) && !wtis_dest_used_vector_l.at((*it).first.at(2)).at(i)){
                                for(int j=0; j<5; j++){
                                        scs_e_fill[1][fill_start][j]=(*(*it).second).wtis_dest[i][j];
                                }
                                scs_e_fill[1][fill_start][2]=(*(*it).second).src;
				wtis_dest_used_vector_l.at((*it).first.at(2)).at(i)=true;
				fill_start++;
                        }
                }
		if(selected_per_round[round] <= (*it).first.at(0)){
	//		selected_per_round[round]=(*it).first.at(0);
			if(round > 0 && selected_per_round[round-1] == (*it).first.at(0)){
				if(max1_start < (*it).first.at(1) || (max1_start == (*it).first.at(1) &&  max2_list_start < (*it).first.at(2))  || (max1_start == (*it).first.at(1) &&  max2_list_start == (*it).first.at(2) && max2_list_max2_start <= (*it).first.at(3))){
					selected_per_round_local.push_back((*it).first.at(0));
					greedy_selection_recursive( scs, scs_e_fill, scs_e_start+ (*it).first.at(0), start, end,arr1,SVs_CNs1,  max1, arr2_list, SVs_CNs2_list,  max2_list,  wtis_src, wtis_src_n,  multi, total_multi_initial, total_multi_after,wtis_src_used_l, wtis_dest_used_vector_l, round+1, selected_per_round, selected_per_round_local, (*it).first.at(1), (*it).first.at(2), (*it).first.at(3),branch_limit, to_hidden_edges_initial, to_hidden_edges);
					selected_per_round_local.pop_back();
				}	
			}else{
				selected_per_round_local.push_back((*it).first.at(0));
				greedy_selection_recursive( scs, scs_e_fill, scs_e_start+ (*it).first.at(0), start, end,arr1,SVs_CNs1,  max1, arr2_list, SVs_CNs2_list, max2_list,  wtis_src, wtis_src_n,  multi, total_multi_initial, total_multi_after,wtis_src_used_l, wtis_dest_used_vector_l, round+1, selected_per_round, selected_per_round_local, (*it).first.at(1), (*it).first.at(2), (*it).first.at(3),branch_limit,to_hidden_edges_initial, to_hidden_edges);
				selected_per_round_local.pop_back();
			}
		}
		if(total_multi_after>0){
			for(int i=0; i<total_multi_initial; i++){
				delete[] scs_e_fill[0][i];
				delete[] scs_e_fill[1][i];
			}
			delete[] scs_e_fill[0];
			delete[] scs_e_fill[1];
			delete[] scs_e_fill;
		}
	
                max1.at((*it).first.at(1)) += (*it).first.at(0);
                max2_list.at((*it).first.at(2)).at((*it).first.at(3)) += (*it).first.at(0);
                multi.at((*it).first.at(2)) += (*it).first.at(0);
                total_multi_after += (*it).first.at(0);


	}
}



int min(int n1, int n2, int n3){
	int which_min=0;
	if(n1<n2)
		which_min=n1;
	else
		which_min=n2;
	if(which_min > n3)
		which_min=n3;
	return which_min;
}


void scs_e_sort(int*** scs_e, int total_multi_initial){
	vector<vector<int>> vec_vec;
	for(int i=0; i<total_multi_initial; i++){
		vector<int> vec;
		for(int j=0; j<4; j++){
			vec.push_back(scs_e[0][i][j]);
		}	
                for(int j=0; j<4; j++){
                        vec.push_back(scs_e[1][i][j]);
                }
		vec_vec.push_back(vec);
	}
	sort(vec_vec.begin(), vec_vec.end(), sortcol);
	for(int i=0; i<total_multi_initial; i++){
		for(int j=0; j<4; j++){
			scs_e[0][i][j]=vec_vec.at(i).at(j);
		}
               for(int j=0; j<4; j++){
                        scs_e[1][i][j]=vec_vec.at(i).at(j+4);
                }
	}
}

bool sortcol(const vector<int>& v1, const vector<int>& v2){
	if( v1[2]<v2[2] || (v1[2]==v2[2] && v1[3] <v2[3])  || (v1[2]==v2[2] && v1[3] ==v2[3] && v1[6] < v2[6]) || ((v1[2]==v2[2] && v1[3] ==v2[3] && v1[6] == v2[6] && v1[7] <v2[7]))){
		return true;
	}
//	if( v1[2]<v2[2] || (v1[2]==v2[2] && v1[3] <v2[3]) || (v1[2]==v2[2] && v1[3] ==v2[3] && v1[6] < v2[6]) || (v1[2]==v2[2] && v1[3] ==v2[3] && v1[6] == v2[6] && v1[7]<=v2[7])){
//		return true;
//	}
	return false;
}

void scs_push_back(list<int***>* scs, int*** scs_e, int total_multi_initial){

	list<int***>::iterator it ;
	for(it = scs->begin(); it!=scs->end(); it++){
		int i;
		for(i=0; i<total_multi_initial; i++){
			if(!((*it)[0][i][2] == scs_e[0][i][2] && (*it)[0][i][3] == scs_e[0][i][3] && (*it)[1][i][2] == scs_e[1][i][2] && (*it)[1][i][3] == scs_e[1][i][3]))
				break;
		}
		if(i== total_multi_initial)
			break;
	}
           cout<<"check"<<endl;

	if(scs->size()==0 || it==scs->end()){
              scs->push_back(scs_e);
	}else{
		for(int i=0; i<total_multi_initial; i++){
			delete[] scs_e[0][i];
			delete[] scs_e[1][i];
		}
		delete[] scs_e[0];
		delete[] scs_e[1];
		delete[] scs_e;

	}
}
void greedy_selected_push_back( list<pair<vector<int>, list<full_edges_add>::iterator>>* greedy_selected, list<int>* greedy_selected_SVs_CN_diff, pair<vector<int>, list<full_edges_add>::iterator> p, int  SVs_CN_diff){
	list<pair<vector<int>, list<full_edges_add>::iterator>>::iterator it  = greedy_selected->begin();
	list<int>::iterator jt;
	for(jt = greedy_selected_SVs_CN_diff->begin(); jt != greedy_selected_SVs_CN_diff->end(); jt++){
		if((*jt)>SVs_CN_diff)
			break;
		it++;
	}
	greedy_selected->insert(it, p);
	greedy_selected_SVs_CN_diff->insert(jt, SVs_CN_diff);
	return;
}
