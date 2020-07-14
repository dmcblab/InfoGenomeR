#include <iostream>
#include <string.h>
#include <algorithm>
#include <list>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <iterator>
#include <unordered_map>
#include <thread>
#include <mutex>

#include "combination_precut.h"
#include "set_cover_branch_limit.h"

using namespace std;
void cc(int*** t);

int main(){
	unsigned int a = 10000;
	double t = a;
	cout<<t<<endl;
	int degree1=6;
	int degree2=1;
	cout<<abs(((double)degree1/degree2)-100) << endl;

                    int ***wtis;
                    wtis = new int**[2];
                    wtis[0] = new int*[degree1];
                    wtis[1] = new int*[degree2];
                    for(int i=0; i<degree1; i++){
                        wtis[0][i] = new int[5];
                    }
                   for(int i=0; i<degree2; i++){
                        wtis[1][i] = new int[5];
                   }
                wtis[0][0][3]=0;
                wtis[0][1][3]=0;
                wtis[0][2][3]=0;
                wtis[0][3][3]=1;
                wtis[0][4][3]=1;
                wtis[0][5][3]=1;

                wtis[0][0][4]=3;
                wtis[0][1][4]=3;
                wtis[0][2][4]=3;
                wtis[0][3][4]=3;
                wtis[0][4][4]=3;
                wtis[0][5][4]=3;



                wtis[1][0][3]=0;
                wtis[1][0][4]=0;






                      clock_t tStart = clock();

                        full_edges_add edges_add;
                        edges_add.src = 2;
                        edges_add.dest = 3;
                        edges_add.degree1_com=degree1;
                        edges_add.degree2_com=degree2;
                        edges_add.multiplicity=1;
                        edges_add.wtis_dest=wtis[1];

		list<full_edges_add> edges_add_list;
		edges_add_list.push_back(edges_add);

                list<int***> tcs = set_cover(edges_add_list, wtis[0], degree1, 10);
 //               cout<<tcs.size()<<endl;
 		cout<<"tcssize"<<tcs.size()<<endl;
                       printf("Time taken: %.8fs\t%d\n", (double)(clock() - tStart)/(CLOCKS_PER_SEC));

		for(list<int***>::iterator it = tcs.begin(); it!=tcs.end(); it++){
			int*** n = (*it);
                       for(int i=0; i<1; i++){
                                cout<<n[0][i][0]<<"\t"<<n[0][i][1]<<"\t"<<n[0][i][2]<<"\t"<<n[0][i][3]<<"\t"<<n[0][i][4]<<endl;
                        }
                        cout<<endl;
                        for(int i=0; i<1; i++){
                                cout<<n[1][i][0]<<"\t"<<n[1][i][1]<<"\t"<<n[1][i][2]<<"\t"<<n[1][i][3]<<"\t"<<n[1][i][4]<<endl;
                        }
                        cout<<endl;

		}
}	
