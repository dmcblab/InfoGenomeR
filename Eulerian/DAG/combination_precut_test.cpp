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
using namespace std;
int branch=99999;

bool func(pair<int,int> i, pair<int,int> j){
        return i.second > j.second;
}

bool func2(pair<int***,double> i, pair<int***,double> j){
        return i.second > j.second;
}

list<int***>  total_combinations_full (list<full_edges_add> full_edges_add_list , int** wtis_src, int wtis_src_n, int total_multi){
	list<int***> tcs;
//	int total_multi=0;
//	for(list<full_edges_add>::iterator it = full_edges_add_list.begin(); it!=full_edges_add_list.end(); it++){
//		total_multi+=(*it).multiplicity;
//	}
	int *** tcs_e;
	tcs_e = new int**[2];
	tcs_e[0] = new int*[total_multi];
        tcs_e[1] = new int*[total_multi];
	for(int i=0; i<total_multi; i++){
		tcs_e[0][i] = new int[4];
		tcs_e[1][i] = new int[4];
		for(int j=0; j<4; j++){
			tcs_e[0][i][j]=-1;
			tcs_e[1][i][j]=-1;
		}
	}

	int t_i=0;

	if(wtis_src_n == 0){
		tcs.push_back(tcs_e);
		return tcs;
	}

//	int wtis_src_used[wtis_src_n];
	total_combinations_full_recursive (&tcs, tcs_e,0, full_edges_add_list.begin(), full_edges_add_list.end(), wtis_src, wtis_src_n, total_multi);
	
	return tcs;
}

void total_combinations_full_recursive( list<int***> *tcs_p, int*** tcs_e, int tcs_e_start, list<full_edges_add>::iterator start, list<full_edges_add>::iterator end, int **wtis_src, int wtis_src_n, int total_multi){
        
	cout<<"heoo"<<endl;
	if(start == end){
		tcs_p->push_back(tcs_e);
		cout<<tcs_e[1][0][0]<<"\t"<<tcs_e[1][0][1]<<"\t"<<tcs_e[1][0][2]<<"\t"<<tcs_e[1][0][3]<<endl;
		return;
	}

	int count=0;
	for(int i=0; i<(*start).degree1_com; i++){
		if(wtis_src[i][2]<0)
			count++;
	}

	int *** temp;
	temp = new int **[2];
	temp[0] = new int* [count];
	for(int i=0; i<count; i++){
		temp[0][i]= new int[4];
	}
	temp[1] = new int* [(*start).degree2_com];
	for(int i=0; i<(*start).degree2_com; i++){
		temp[1][i] = new int [4];
	}


	int temp_i=0;
	for(int i=0; i<(*start).degree1_com; i++){
		if(wtis_src[i][2]<0){
			temp[0][temp_i][0]=wtis_src[i][0];
			temp[0][temp_i][1]=wtis_src[i][1];
			temp[0][temp_i][2]=(*start).dest;
			temp[0][temp_i][3]=wtis_src[i][3];
			temp_i++;
		}
	}
	for(int i=0; i<(*start).degree2_com; i++){
		temp[1][i][0]=(*start).wtis_dest[i][0];
		temp[1][i][1]=(*start).wtis_dest[i][1];
		temp[1][i][2]=(*start).src;
		temp[1][i][3]=(*start).wtis_dest[i][3];
	}

	list<int***> partial_combination = total_combinations(count, (*start).degree2_com, (*start).multiplicity, temp);
        cout<<"hi"<<count<<"\t"<<(*start).degree2_com<<"\t"<<(*start).multiplicity<<"\t"<<partial_combination.size()<<endl;

	for(list<int***>::iterator it = partial_combination.begin(); it!=partial_combination.end(); it++){
		int *** tcs_e_fill;
		tcs_e_fill = new int **[2];
	        tcs_e_fill[0] = new int*[total_multi];
		tcs_e_fill[1] = new int*[total_multi];
	        for(int i=0; i<total_multi; i++){
	 	       tcs_e_fill[0][i] = new int[4];
	        }
	        for(int i=0; i<total_multi; i++){
	        	tcs_e_fill[1][i] = new int[4];
	        }

		for(int i=0; i<total_multi; i++){
			for(int j=0; j<4; j++){
				tcs_e_fill[0][i][j]=tcs_e[0][i][j];
				tcs_e_fill[1][i][j]=tcs_e[1][i][j];
			}
		}	
		int k=0;
		for(int i=tcs_e_start; i<tcs_e_start+(*start).multiplicity; i++){
			for(int j=0; j<4; j++){
				tcs_e_fill[0][i][j]=(*it)[0][k][j];
				tcs_e_fill[1][i][j]=(*it)[1][k][j];
			}
			k++;
		}
	
		list<int> used;	
		for(int i=0; i<(*start).multiplicity; i++){
			for(int j=0; j<wtis_src_n; j++){
				if((*it)[0][i][3] == wtis_src[j][3] && wtis_src[j][2] < 0){
					wtis_src[j][2] = (*start).dest;
					used.push_back(j);
					break;
				}
			}
		}
		list<full_edges_add>::iterator start_new = start;
		start_new ++;
		total_combinations_full_recursive(tcs_p, tcs_e_fill, tcs_e_start+(*start).multiplicity, start_new ,end, wtis_src,  wtis_src_n, total_multi);


		for(list<int>::iterator jt=used.begin(); jt !=used.end(); jt++){
			wtis_src[(*jt)][2] = -1;
		}	
		used.clear();
	}
}


list<int***> total_combinations(int paths1_n, int paths2_n, int multiplicity, int ***wtis){
	list<int***> tcs;
	list<pair<int***, double>> tcs_with_entropy;
	list<int> tcs_entropy;	


	list<int> used;
	for(int i=0; i<paths1_n; i++){
		if(find(used.begin(), used.end(),wtis[0][i][0]) == used.end()){
			used.push_back(wtis[0][i][0]);
		}
	}
	int n1=used.size();
	used.clear();
        for(int i=0; i<paths2_n; i++){
                if(find(used.begin(), used.end(),wtis[1][i][0]) == used.end()){
                        used.push_back(wtis[1][i][0]);
                }
        }
	int n2=used.size();



	vector<int> arr;
        vector<int> max;
	
	
        for(int i=0; i<paths1_n; i++){
		vector<int>::iterator it = max.begin();
		vector<int>::iterator jt;
		for(jt = arr.begin(); jt!=arr.end(); jt++){
			if(*jt==wtis[0][i][3]){
				(*it)++;
				break;
			}
			it++;
		}
		if(jt == arr.end()){
			arr.push_back(wtis[0][i][3]);
			max.push_back(1);
		}
        }
//	cout<<"precombination"<<max.at(0)<<endl;
//        clock_t addStart = clock();


/*
	cout<<"here"<<endl;
	for(int i=0; i<arr.size(); i++){
	        cout<<arr.at(i)<<"\t";
	}
	cout<<endl;
        for(int i=0; i<max.size(); i++){
                cout<<max.at(i)<<"\t";
        }
        cout<<endl;
*/
/*
        vector<pair<int,int>> sort_l;

        vector<int>::iterator jt=max.begin();
        for(vector<int>::iterator it = arr.begin(); it!=arr.end(); it++){
                sort_l.push_back(pair<int,int>(*it, *jt));
                jt++;
        }
        sort(sort_l.begin(), sort_l.end(), func);
         jt=max.begin();
        vector<int>::iterator it = arr.begin();
        for(vector<pair<int,int>>::iterator vt=sort_l.begin(); vt!=sort_l.end(); vt++){
                *it=vt->first;
                *jt=vt->second;
                jt++;
                it++;
        }
        cout<<"here"<<endl;
        for(int i=0; i<arr.size(); i++){
                cout<<arr.at(i)<<"\t";
        }
        cout<<endl;
        for(int i=0; i<max.size(); i++){
                cout<<max.at(i)<<"\t";
        }
        cout<<endl;
*/


	list<vector<int>> com1= combination (arr, max, arr.size(), multiplicity);
	cout<<"comfinished"<<endl;
	cout<<"precombination"<<com1.size()<<endl;
       // printf("addTime taken: %.8fs\t%d\n", (double)(clock() - addStart)/(CLOCKS_PER_SEC), newNode_i);

/*
	for(list<vector<int>>::iterator it=com1.begin(); it!=com1.end(); it++){
		for(vector<int>::iterator jt=(*it).begin(); jt!=(*it).end(); jt++){
			cout<<*jt<<"\t";
		}
		cout<<endl;
	}	
	cout<<endl;
*/
	int m2[paths2_n];
	for(int i=0; i<paths2_n; i++){
		m2[i]=wtis[1][i][3];	
	}


	for(list<vector<int>>::iterator it=com1.begin(); it!=com1.end(); it++){
		int m1[multiplicity];
		int m1_com[multiplicity];

		list<int> used;
//		cout<<"here"<<endl;
  //              for(int i=0; i<(*it).size();i++){
//			cout<<(*it).at(i)<<"\t";
//		}
		for(int i=0; i<(*it).size();i++){
			for(int j=0; j<paths1_n; j++){
				if(wtis[0][j][3]== (*it).at(i) && find(used.begin(), used.end(), j)==used.end()){
					m1[i]=wtis[0][j][3];
					m1_com[i]=j;
					used.push_back(j);
					break;
				}
			}
		}
		used.clear();

		list<vector<int>> com2 = test(m1,m2,paths2_n, multiplicity);

		cout<<"aftercombination"<<com2.size()<<endl;

//		double min_entropy=9999;
		for(list<vector<int>>::iterator it = com2.begin(); it!=com2.end(); it++){
		   int ***wtis_com;
		    wtis_com = new int**[2];
		    wtis_com[0] = new int*[multiplicity];
		    wtis_com[1] = new int*[multiplicity];
		    for(int i=0; i<multiplicity; i++){
			wtis_com[0][i] = new int[4];
		    }
		   for(int i=0; i<multiplicity; i++){
			wtis_com[1][i] = new int[4];
		    }

			for(int i=0; i<multiplicity; i++){
				for(int j=0; j<4; j++){
					wtis_com[0][i][j]=wtis[0][m1_com[i]][j];
				}
			}

			int m2_com[multiplicity];
			used.clear();
			for(int i=0; i<(*it).size();i++){
				for(int j=0; j<paths2_n; j++){
					if(wtis[1][j][3]== (*it).at(i) && find(used.begin(), used.end(), j)==used.end()){
						m2_com[i]=j;
						used.push_back(j);
						break;
					}
				}
			}
			for(int i=0; i<multiplicity; i++){
				for(int j=0; j<4; j++){
					wtis_com[1][i][j]=wtis[1][m2_com[i]][j];
				}
    		        }
		
			list<vector<int>> unique_ele;
			list<int> number_ele;
			for(int i=0; i<multiplicity; i++){
				list<int>::iterator neit=number_ele.begin();
				list<vector<int>>::iterator eit;
				for(eit=unique_ele.begin(); eit!=unique_ele.end();eit++){
					if((*eit).at(0) == wtis_com[0][i][3] && (*eit).at(1) == wtis_com[1][i][3]){
						(*neit)++;
						break;
					}
					neit++;
				}
				if(eit == unique_ele.end()){
					unique_ele.push_back(vector<int>{wtis_com[0][i][3], wtis_com[1][i][3]});
					number_ele.push_back(1);
				}
			}
			double entropy=0;
			for(list<int>::iterator neit=number_ele.begin(); neit!=number_ele.end(); neit++){
				entropy-=((*neit)/(double)multiplicity)*log((*neit)/(double)multiplicity);
			}
/*
		cout<<"entriopy"<<entropy<<endl;
               for(int i=0; i<multiplicity; i++){
                        cout<<wtis_com[0][i][0]<<"\t"<<wtis_com[0][i][1]<<"\t"<<wtis_com[0][i][2]<<"\t"<<wtis_com[0][i][3]<<endl;
                }
                cout<<endl;
                for(int i=0; i<multiplicity; i++){
                        cout<<wtis_com[1][i][0]<<"\t"<<wtis_com[1][i][1]<<"\t"<<wtis_com[1][i][2]<<"\t"<<wtis_com[1][i][3]<<endl;
                }
*/	
//		if(min_entropy> entropy){
			tcs_with_entropy.push_back(pair<int***, double>(wtis_com, entropy));
//			min_entropy= entropy;
//		}

//			tcs_insert(&tcs, &tcs_entropy, wtis_com, entropy);
//			tcs.push_back(wtis_com);
		/*	for(int i=0; i<multiplicity; i++){
				for(int j=0; j<4; j++)
					cout<<wtis_com[0][i][j]<<"\t";
                                for(int j=0; j<4; j++)
                                        cout<<wtis_com[1][i][j]<<"\t";
				cout<<endl;
			}

			cout<<endl;
		*/
		}

		}
		
  	       tcs_with_entropy.sort(func2);
		int max_tcs_i=2;
		for(list<pair<int***,double>>::iterator it=tcs_with_entropy.begin(); it!=tcs_with_entropy.end(); it++){
			if(max_tcs_i > 0){
				tcs.push_back(it->first);
			}else{
  	                      for(int i=0; i<2; i++){
                               		for(int j=0; j<multiplicity; j++){
                                        delete[] (it->first)[i][j];
                                	}
                                delete[] (it->first)[i];
                        	}
                        	delete[] (it->first);

			}
			max_tcs_i--;
		}
		tcs_with_entropy.clear();	
		return tcs;

	}

list<vector<int>> test(int m1[],int  m2[], int paths2_n, int multiplicity){
	int index=0;
	list<vector<int>> com;
	list<int> com_entropy;
        int m1_max=0;
        for(int i=0; i<multiplicity; i++){
                if(m1[i]>m1_max){
                        m1_max=m1[i];
                }
        }

	int m1_count[multiplicity];
	for(int i=0; i<multiplicity; i++){
		m1_count[i]=0;
	}

	vector<int> arr_m1;
	vector<int> max_m1;
        for(int i=0; i<multiplicity; i++){
                vector<int>::iterator jt=max_m1.begin();
                vector<int>::iterator it;
                for(it=arr_m1.begin(); it!=arr_m1.end(); it++){
                        if(*it==m1[i]){
                                break;
                        }
                        jt++;
                }
                if(it==arr_m1.end()){
                        arr_m1.push_back(m1[i]);
                        max_m1.push_back(1);
                }else{
                        *jt+=1;
                }
        }


        vector<int> arr;
        vector<int> max;
        for(int i=0; i<paths2_n; i++){
                vector<int>::iterator jt=max.begin();
                vector<int>::iterator it;
                for(it=arr.begin(); it!=arr.end(); it++){
                        if(*it==m2[i]){
                                break;
                        }
                        jt++;
                }
                if(it==arr.end()){
                        arr.push_back(m2[i]);
                        max.push_back(1);
                }else{
                        *jt+=1;
                }
	}

        for(int i=0; i<max_m1.size(); i++){
		m1_count[i]=max_m1.at(i);
        }
	m1_max=max_m1.size()-1;

/*
        vector<pair<int,int>> sort_l;

        vector<int>::iterator jt=max.begin();
        for(vector<int>::iterator it = arr.begin(); it!=arr.end(); it++){
                sort_l.push_back(pair<int,int>(*it, *jt));
                jt++;
        }
//         sort(sort_l.begin(), sort_l.end(), func);
         jt=max.begin();
        vector<int>::iterator it = arr.begin();
        for(vector<pair<int,int>>::iterator vt=sort_l.begin(); vt!=sort_l.end(); vt++){
                *it=vt->first;
                *jt=vt->second;
                jt++;
                it++;
        }
*/

	vector<int> tv;
	cout<<"arr_size"<<arr.size()<<endl;
	recursive_combination(m1_count, &com, &com_entropy, index, m1_max,arr, max,tv, multiplicity);
	return com;

}

void tcs_insert ( list<int***>* tcs, list<int>* tcs_entropy, int*** wtis_com ,int entropy){
		if(tcs->size()==0){
			tcs->push_back(wtis_com);
			tcs_entropy->push_back(entropy);
		}else{
			list<int***>::iterator it=tcs->begin();
			list<int>::iterator jt;
			for(jt = tcs_entropy->begin(); jt!=tcs_entropy->end(); jt++){
				cout<<"tcs_insert"<<endl;
				if(entropy<=*jt){
					break;
				}
				it++;
			}
 	               tcs->insert(it, wtis_com);
        	        tcs_entropy->insert(jt, entropy);
		}
}

void com_insert( list<vector<int>> *com,  list<int> *com_entropy ,  vector<int> tv, int entropy){
		list<vector<int>>::iterator it = com->begin();
		list<int>::iterator jt = com_entropy->begin();
		for(it=com->begin(); it!=com->end(); it++){
			if(entropy<*jt){
				break;
			}
			jt++;
		}
		com->insert(it, tv);
		com_entropy->insert(jt, entropy);
		if(com->size() > branch){
			com->pop_back();
			com_entropy->pop_back();
		}
}

double tv_entropy(vector<int> tv, int m1_count[], int multiplicity, vector<int> max){
	cout<<"tv_entropy"<<endl;
	for(vector<int>::iterator it= tv.begin(); it!=tv.end(); it++){
		cout<<*it<<"\t";
	}
	cout<<endl;
	cout<<"m1_count"<<multiplicity<<endl;
	for(int i=0; i<multiplicity; i++){
		cout<<m1_count[i]<<"\t";
	}
	cout<<endl;
	cout<<"max"<<endl;
        for(vector<int>::iterator it= max.begin(); it!=max.end(); it++){
                cout<<*it<<"\t";
        }
        cout<<endl;


	vector<int> count;
	int tv_i=0;
	for(int i=0; i<multiplicity; i++){
		if(m1_count[i]!=0){
			int upto= tv_i + m1_count[i];
                        vector<int> used;
			
			int count_i=count.size();

			for(;tv_i<upto;tv_i++){
				int count_j=count_i;
				vector<int>::iterator it;
				for(it=used.begin(); it!=used.end(); it++){
					if(*it == tv.at(tv_i)){
						break;
					}
					count_j++;
				}
				if(it!=used.end()){
					count.at(count_j)++;
				}else{
					used.push_back(tv.at(tv_i));
					count.push_back(1);
				}
			}
		}
	}
	double entropy=0;
	int sum=0;
	for(vector<int>::iterator it=count.begin(); it!=count.end();it++){
		sum+=*it;
	}
	for(vector<int>::iterator it=max.begin(); it!=max.end();it++){
                sum+=*it;
	}
        for(vector<int>::iterator it=count.begin(); it!=count.end();it++){
		if(*it!=0)
	        	entropy-=((*it)/(double)sum)*log((*it)/(double)sum);
			cout<<((*it)/(double)sum)<<endl;
			cout<<((*it)/(double)sum)*log((*it)/(double)sum)<<endl;
	}
        for(vector<int>::iterator it=max.begin(); it!=max.end();it++){
                if(*it!=0)
                        entropy-=((*it)/(double)sum)*log((*it)/(double)sum);
        }
	cout<<"entropy"<<entropy<<endl;
	return entropy;

	
}

void recursive_combination(int m1_count[], list<vector<int>> *com, list<int> *com_entropy, int index, int m1_max, vector<int> arr, vector<int> max, vector<int> tv, int multiplicity){
	if( index >m1_max){
		com_insert(com, com_entropy, tv, tv_entropy(tv, m1_count, multiplicity, max));
		return;
	}
	cout<<"combinationenter"<<endl;	
	cout<<arr.size()<<"\t"<<index<<"\t"<<max.size()<<"\t"<<m1_count[index]<<"\t"<<m1_max<<endl;
        list<vector<int>> vl = combination( arr, max, arr.size(), m1_count[index]);
	cout<<"combinationend"<<endl;

	int z=99999;
	for(list<vector<int>>::iterator it=vl.begin(); it!=vl.end() && z>0;it++){
		vector<int> l_max=max;
		vector<int> l_tv=tv;
		for(vector<int>::iterator jt=(*it).begin(); jt!=(*it).end();jt++){
			l_tv.push_back(*jt);
			int m_i=0;
			for(int i=0; i<arr.size(); i++){
				if(arr.at(i)==*jt)
					m_i=i;
			}
                        l_max.at(m_i)--;
                      //  l_max.at((*jt))--;


		}
		cout<<"recur enter"<<endl;
		recursive_combination(m1_count, com, com_entropy, index+1, m1_max, arr,l_max, l_tv, multiplicity);
		cout<<"recur end"<<endl;
		z--;
	} 

}

list<vector<int>> combination(vector<int> arr, vector<int> max, int n, int r){

	list<vector<int>> vl;
	list<vector<int>>* vlp = &vl;
	int chosen[r+1];
	int max_sum=0;
	for(int i=0; i<max.size(); i++){
		max_sum+=max.at(i);
	}
	combinationutil(chosen, arr,max, 0, r, 0, n-1, vlp,max_sum);
	return vl;
}
void combinationutil(int chosen[], vector<int> arr, vector<int> max, int index, int r, int start, int end, list<vector<int>> *vlp, int  max_sum){
        max_sum=0;
        for(int i=start; i<max.size(); i++){
                max_sum+=max.at(i);
        }

        if(max_sum < r-index-1){
                return;
        }
        if(index==r){
                vector<int> tv;
                for( int i=0; i<r; i++){
                        tv.push_back(chosen[i]);
                }
                if(tv_validity(tv, r, arr,max)){
	             vlp->push_back(tv);
/*			cout<<"combinationtest:";
			for(int z=0; z<tv.size(); z++){
				cout<<tv.at(z)<<"\t";
			}
			cout<<endl;
*/
                }else{
                }
	//    cout<<"r"<<r<<endl;
 /*            cout<<"combinationutil:";
                for( int i=0; i<r; i++){
                        cout<<chosen[i]<<"\t";
                }
                cout<<endl;*/

                return;
        }
        for(int i=start; i<=end; i++){
                if(max.at(i)>0){
                       chosen[index]=arr[i];
                        max.at(i)--;
                        max_sum--;
                        combinationutil(chosen, arr,max, index+1, r, i ,end, vlp, max_sum);
                        max.at(i)++;
                }
        }
        return;
}



bool imcomplete_occupied_exist(int chosen[], vector<int> max, int index){
	for(int i=0; i<index; i++){
		if(max.at(chosen[i])!=0){
			return true;
		}
	}
	return false;
}


bool tv_validity(vector<int> tv, int r, vector<int> arr , vector<int> max){
/*
        for(int i=0; i<arr.size(); i++){
                cout<<arr.at(i)<<"\t";
        }
        cout<<endl;
        for(int i=0; i<max.size(); i++){
                cout<<max.at(i)<<"\t";
        }
        cout<<endl;
*/
        int how_many_index_in_tv[arr.size()];
        for(int i=0; i<arr.size(); i++){
		how_many_index_in_tv[i]=0;
	}
	for(int i=0; i<tv.size(); i++){
		how_many_index_in_tv[tv.at(i)]++;
	}
	
	int how_many_non_zero=0;
	int how_many_max_zero=0;
	for(int i=0; i<arr.size(); i++){
		if(how_many_index_in_tv[i]!=0){
			how_many_non_zero++;
			if(max.at(i)==0){
				how_many_max_zero++;
			}
		}
	}

	if(how_many_non_zero <1){
		return false;
	}
	return true;


	if(how_many_non_zero == 1){
		return true;
	}
	if(how_many_non_zero > 1 && how_many_max_zero >= how_many_non_zero -1 ){
		return true;
	}else{
		return false;
	}
/*
	if(max.size()>1){
	        cout<<"tv_validity"<<how_many_index_in_tv[0]<<"\t"<<max.at(0)<<endl;
	        cout<<"tv_validity"<<how_many_index_in_tv[1]<<"\t"<<max.at(1)<<endl;
	}


        for(int i=0; i<r; i++){
		if(how_many_index_in_tv[i]!=0){
                        if((r >= how_many_index_in_tv[i] && max.at(i)==0) || (r < max.at(i)+how_many_index_in_tv[i] && how_many_index_in_tv[i] == r)){
                        }else{
                                return false;
                        }
		}
	}
	return true;
*/
}
