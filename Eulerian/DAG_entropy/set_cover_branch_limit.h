#include <iostream>
#include <string.h>
#include <algorithm>
#include <list>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <iterator>

using namespace std;
list<int***> set_cover(list<full_edges_add> full_edges_add_list, int** wtis_src, int wtis_src_n, int branch_limit);
//void greedy_selection_recursive( list<int***>* scs, int*** scs_e, int scs_e_start, list<full_edges_add>::iterator start, list<full_edges_add>::iterator end, vector<int> arr1, vector<int> max1, vector<vector<int>> arr2_list, vector<vector<int>> max2_list, int** wtis_src, int wtis_src_n, vector<int> multi, int total_multi_initial, int total_multi_after, vector<bool> wtis_src_used, vector<vector<bool>> wtis_dest_used_vector);
int min(int n1, int n2, int n3);
void scs_e_sort(int*** scs_e, int total_multi_initial);
bool sortcol(const vector<int>& v1, const vector<int>& v2);
void scs_push_back(list<int***>* scs, int*** scs_e, int total_multi_initial);

//void greedy_selection_recursive( list<int***>* scs, int*** scs_e, int scs_e_start, list<full_edges_add>::iterator start, list<full_edges_add>::iterator end, vector<int> arr1, vector<int> max1, vector<vector<int>> arr2_list, vector<vector<int>> max2_list, int** wtis_src, int wtis_src_n, vector<int> multi, int total_multi_initial, int total_multi_after, vector<bool> wtis_src_used, vector<vector<bool>> wtis_dest_used_vector, int round, int selected_per_round[]);
void greedy_selection_recursive( list<int***>* scs, int*** scs_e, int scs_e_start, list<full_edges_add>::iterator start, list<full_edges_add>::iterator end, vector<int> arr1, vector<int> SVs_CNs1, vector<int> max1, vector<vector<int>> arr2_list, vector<vector<int>> SVs_CNs2_list, vector<vector<int>> max2_list, int** wtis_src, int wtis_src_n, vector<int> multi, int total_multi_initial, int total_multi_after, vector<bool> wtis_src_used, vector<vector<bool>> wtis_dest_used_vector, int round, int selected_per_round[], vector<int> selected_per_round_local, int max1_start, int max2_list_start, int max2_list_max2_start, int branch_limit, int to_hidden_edges_initial, int to_hidden_edges);
//void greedy_selection_recursive( list<int***>* scs, int*** scs_e, int scs_e_start, list<full_edges_add>::iterator start, list<full_edges_add>::iterator end, vector<int> arr1, vector<int> max1, vector<vector<int>> arr2_list, vector<vector<int>> max2_list, int** wtis_src, int wtis_src_n, vector<int> multi, int total_multi_initial, int total_multi_after, vector<bool> wtis_src_used, vector<vector<bool>> wtis_dest_used_vector, int round, int selected_per_round[], vector<int> selected_per_round_local, int max1_start, int max2_list_start, int max2_list_max2_start, int branch_limit, int to_hidden_edges);
void greedy_selected_push_back( list<pair<vector<int>, list<full_edges_add>::iterator>>* greedy_selected, list<int>* greedy_selected_SVs_CN_diff, pair<vector<int>, list<full_edges_add>::iterator> p, int  SVs_CN_diff);

