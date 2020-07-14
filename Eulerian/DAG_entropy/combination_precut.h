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
list<int***> total_combinations(int degree1, int degree2, int multiplicity, int ***wtis);
list<vector<int>> test(int m1[],int  m2[], int degree2, int multiplicity);
void combinationutil(int chosen[], vector<int> arr, vector<int> max, int index, int r, int start, int end, list<vector<int>> *vl, int max_sum);
list<vector<int>> combination(vector<int> arr, vector<int> max, int n, int r);
void recursive_combination(int m1_count[], list<vector<int>> *com, list<int> *com_entropy, int index, int m1_max, vector<int> arr, vector<int> max, vector<int> tv, int multiplicity);
bool tv_validity(vector<int> tv, int r, vector<int> arr, vector<int> max );
void tcs_insert ( list<int***>* tcs, list<int>* tcs_entropy, int*** wtis_com ,int entropy);
bool func(pair<int,int> i, pair<int,int> j);
bool imcomplete_occupied_exist(int chosen[], vector<int> max, int index);
double tv_entropy(vector<int> tv, int m1_count[], int multiplicity, vector<int> max);
struct full_edges_add{
	int src;
	int dest;
        int degree1_com;
        int degree2_com;
        int multiplicity;
        int **wtis_dest;
};
list<int***>  total_combinations_full (list<full_edges_add> full_edges_add_list , int** wtis_src, int wtis_src_n);
void total_combinations_full_recursive( list<int***> *tcs_p, int*** tcs_e, int tcs_e_start, list<full_edges_add>::iterator start, list<full_edges_add>::iterator end, int **wtis_src, int wtis_src_n, int total_multi);

