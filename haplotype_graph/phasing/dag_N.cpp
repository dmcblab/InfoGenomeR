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
using namespace std;
int g_to_int (char g);
char int_to_g (int g);
vector<char> local_back ( int back_track, list<vector<int>**>*  viterbi_back_indices,  list<vector<char>**>*  viterbi_states, vector<int> which_indices);
int main(int argc, char* argv[]){

        ifstream dagf(argv[1]);
	//ifstream dagf("../chr11.chr11.bgl.dag.align");
	string line;
	getline(dagf,line,'\n');
	int level=-1;

	vector<int> v;
	vector<int> prev_v;
	v.push_back(5008);
	
	list< list<vector<double>> [4]>  DAG;
	list<vector<int>> DAG_node_n;

	list<vector<double>> DAG_e [4];
		
	int max_start=0;
	int max_end=0;
      while(getline(dagf,line,'\n')){
                stringstream f_iss(line);

		if(line != "					"){
			string val;
			int j=0;
		 	double  temp [6];
			char g;
			while(getline(f_iss,val,'\t')){
				stringstream convertor(val);
				double t;
				if(j != 4){
					convertor>>t;
				}else{
					convertor>>g;
				}
				temp[j]=t;
				j++;
			}

			int dis = temp[3] - v.size() + 1;
			while( dis > 0){
				v.push_back(0);
				dis--;
			}
			v[temp[3]]+=temp[5];
			DAG_e[g_to_int(g)].push_back({temp[2], temp[3], temp[5]/prev_v[temp[2]]});
			if(max_start < temp[2]){
				max_start=temp[2];
			}
			if(max_end < temp[3]){
				max_end = temp[3];
			}
	//		cout<<line<<"\t"<<g_to_int(g)<<"\t"<<temp[5]/prev_v[temp[2]]<<endl;
		}else{	
			if(level!=-1){
				DAG.push_back(DAG_e);
				DAG_node_n.push_back({max_start, max_end});
			}

			for(int i=0; i<4 ;i++){
				DAG_e[i].clear();	
			}
			prev_v=v;
			v.clear();
			level++;
			max_start=0;
			max_end=0;
		//	cout<<level<<endl;
		}
	}
	DAG.push_back(DAG_e);
	DAG_node_n.push_back({max_start, max_end});



//	cout<<"DAGsize"<<DAG.size()<<endl;
//	cout<<"DAG_nsize"<<DAG_node_n.size()<<endl;

	 list< list<vector<double>> [4]>::iterator it = DAG.begin();
	 list<vector<int>>::iterator jt = DAG_node_n.begin();
//	ifstream gf("genotype.format.11");
        ifstream gf(argv[2]);


	double** prev_viterbi_m;
	prev_viterbi_m = new double* [(*jt)[0] +1 ];
        for( int i=0; i<=(*jt)[0]; i++){
		prev_viterbi_m[i] = new double [(*jt)[0] + 1];
                for( int j=0; j<=(*jt)[0]; j++){
                        prev_viterbi_m[i][j]=0; 
                }
        }
	list<vector<int>**>  viterbi_back_indices;
	list<vector<char>**> viterbi_states;



        while(getline(gf,line,'\n')){
		double** viterbi_m;
		vector<int>** viterbi_index;
		vector<char>** viterbi_state;
		viterbi_m = new double* [(*jt)[1] + 1 ];
		viterbi_index = new vector<int>* [(*jt)[1] + 1];
		viterbi_state = new vector<char>* [(*jt)[1] + 1];
		for(int i=0; i<=(*jt)[1]; i++){
			viterbi_m[i] = new double [(*jt)[1] + 1 ];
			viterbi_index[i] = new vector<int> [(*jt)[1] + 1];
                        viterbi_state[i] = new vector<char> [(*jt)[1] + 1];
			for(int j=0; j<=(*jt)[1]; j++){
				viterbi_m[i][j]= - numeric_limits<double>::infinity();
				viterbi_index[i][j] = {-1,-1,-1,-1};
				viterbi_state[i][j] = {'N','N'};
			}
		}


                stringstream f_iss(line);
		string val;
		char g1;
		char g2;
		int j=0;
		int back_track=0;
		char p_g1, p_g2;
		int coordi;

		while(getline(f_iss, val, '\t')){
			stringstream convertor(val);
			if(j==5){
				convertor >> g1;
			}else if(j==6){
				convertor >> g2;
			}else if(j==8){
				convertor >> p_g1;
			}else if(j==9){
				convertor >> p_g2;
			}else if(j==10){
				convertor >> back_track;
			}else if(j==1){
				convertor >> coordi;
			}

			j++;
		}

		
		int genotype_missing  = 2;
		if( g1 != 'N' && g2 != 'N'){
			for( list<vector<double>> ::iterator lit = (*it)[g_to_int(g1)].begin();	lit != (*it)[g_to_int(g1)].end(); lit++){
				for( list <vector<double>>::iterator ljt = (*it)[g_to_int(g2)].begin(); ljt != (*it)[g_to_int(g2)].end(); ljt++){
					bool b = ( prev_viterbi_m [(int)(*lit)[0]][(int)(*ljt)[0]]  != - numeric_limits<double>::infinity());
					if (back_track !=0 && b){
						vector<char> genotypes = local_back (back_track,  &viterbi_back_indices, & viterbi_states, {(int)(*lit)[0],(int)(*ljt)[0]});
						b = b & (p_g1 == 'N' | genotypes[0] == p_g1) & (p_g2  == 'N' | genotypes[1] == p_g2);
						if( (int)(*lit)[0] == (int)(*ljt)[0] )
							b=true;
					}

					if(b){
						genotype_missing = 0;
						break;
					}
				}
			}
		}
		if(genotype_missing !=0){
			if( g1 !='N'){
				for( list<vector<double>> ::iterator lit = (*it)[g_to_int(g1)].begin(); lit != (*it)[g_to_int(g1)].end(); lit++){
					for( int k =0; k<= (*jt)[0] ; k++){
						bool b = ( prev_viterbi_m [(int)(*lit)[0]][k]  != - numeric_limits<double>::infinity());
						if(back_track !=0 && b){
						       vector<char> genotypes = local_back (back_track,  &viterbi_back_indices, & viterbi_states, {(int)(*lit)[0],k});
							if(g1 !='N')
								b = b & (p_g1 == 'N' | genotypes[0] == p_g1) & (p_g2  == 'N'| genotypes[1] == p_g2);
							else
								b= b & (p_g1 == 'N' | genotypes[0] == p_g1);
							if( (int)(*lit)[0] == k )
								b=true;
						}
						if(b){
							genotype_missing=1;
							break;
						}
					}
				}
			}
			if(genotype_missing != 1){
				if(g2 != 'N'){
					for( list<vector<double>> ::iterator lit = (*it)[g_to_int(g2)].begin(); lit != (*it)[g_to_int(g2)].end(); lit++){
						for( int k =0; k<= (*jt)[0] ; k++){
							bool b = ( prev_viterbi_m[k][(int)(*lit)[0]]  != - numeric_limits<double>::infinity());
							if(back_track !=0 && b){
							       vector<char> genotypes = local_back (back_track,  &viterbi_back_indices, & viterbi_states, {k,(int)(*lit)[0]});
								if(g2 != 'N')
									b = b & (p_g1 == 'N' | genotypes[0] == p_g1) & (p_g2  == 'N'| genotypes[1] == p_g2);
								else
 	                                                               b= b & (p_g2 == 'N' | genotypes[1] == p_g2);

								if( k == (int)(*lit)[0] )
									b=true;
							}
							if(b){
								genotype_missing=1;
								break;
							}
						}
					}
				}
			}
		}
		
//		cout<<"ge"<<genotype_missing<<endl;

		if(genotype_missing == 0){
                if( g1 != 'N' && g2 != 'N'){
			for( list<vector<double>> ::iterator lit = (*it)[g_to_int(g1)].begin(); lit != (*it)[g_to_int(g1)].end(); lit++){
				for( list <vector<double>>::iterator ljt = (*it)[g_to_int(g2)].begin(); ljt != (*it)[g_to_int(g2)].end(); ljt++){
					bool b = prev_viterbi_m [(int)(*lit)[0]][(int)(*ljt)[0]]  != - numeric_limits<double>::infinity() ;
					vector<char> genotypes;
					if(b){
                                                genotypes = local_back (back_track,  &viterbi_back_indices, & viterbi_states, {(int)(*lit)[0],(int)(*ljt)[0]});
						if( back_track !=0){
							b = b & (p_g1 == 'N' | genotypes[0] == p_g1) & (p_g2 == 'N' | genotypes[1] == p_g2);
							if( (int)(*lit)[0] == (int)(*ljt)[0] ){
								b=true;
							}
						}
					}
					if(b ){
						double val = prev_viterbi_m [(int)(*lit)[0]][(int)(*ljt)[0]] + log((*lit)[2])+log((*ljt)[2]);

						if( viterbi_m [(int)(*lit)[1]][(int)(*ljt)[1]] < val) {

							viterbi_m [(int)(*lit)[1]][(int)(*ljt)[1]] = val;
                                                        viterbi_m [(int)(*ljt)[1]][(int)(*lit)[1]] = val;

			
							int i1,i2;
							char s1,s2;
							i1=(int)(*lit)[0];
							i2=(int)(*ljt)[0];
						//	cout<<"gcheck"<<g1<<"\t"<<g2<<endl;
							s1=g1;
							s2=g2;
						
							int b1=0;
							int b2=1;
							//cout<<"here"<<genotypes[0]<<"\"t"<<p_g1<<endl;
					//		if(back_track != 0 && (p_g1 != 'N' & genotypes[0] != p_g1)){
                                                        if(back_track != 0 && ((p_g1 != 'N' & genotypes[0] != p_g1) | (p_g2!= 'N' & genotypes[1] !=p_g2))){
								b1=1;
								b2=0;
							}


							viterbi_index [(int)(*lit)[1]][(int)(*ljt)[1]] = {i1,i2,b1,b2};
                                                        viterbi_state [(int)(*lit)[1]][(int)(*ljt)[1]] = {s1,s2};
							if((int)(*lit)[1] != (int)(*ljt)[1]){
							       viterbi_index [(int)(*ljt)[1]][(int)(*lit)[1]] = {i1,i2,b2,b1};
								viterbi_state [(int)(*ljt)[1]][(int)(*lit)[1]] = {s2,s1};
							}
						}
					}
				}
			}
		}
		}else if (genotype_missing ==1){
			if( g1 != 'N'){
                        for( list<vector<double>> ::iterator lit = (*it)[g_to_int(g1)].begin(); lit != (*it)[g_to_int(g1)].end(); lit++){
				for( int k=0; k<4; k++){
				       if( g2 == 'N' || k != g_to_int(g2)){
                                //       if( 1){
	 	                               for( list <vector<double>>::iterator ljt = (*it)[k].begin(); ljt != (*it)[k].end(); ljt++){
							bool b = prev_viterbi_m [(int)(*lit)[0]][(int)(*ljt)[0]] != - numeric_limits<double>::infinity() ;
                                                        vector<char> genotypes;
						       if(b){
								genotypes = local_back (back_track,  &viterbi_back_indices, & viterbi_states, {(int)(*lit)[0],(int)(*ljt)[0]});
								if( back_track !=0){
									b = b & (p_g1 == 'N' | genotypes[0] == p_g1) & (p_g2 == 'N' | genotypes[1] == p_g2);
									if( (int)(*lit)[0] == (int)(*ljt)[0] ){
										b=true;
									}
								}
							}
							if(b){
							       double val = prev_viterbi_m [(int)(*lit)[0]][(int)(*ljt)[0]] + log((*lit)[2])+log((*ljt)[2]);
								if( viterbi_m [(int)(*lit)[1]][(int)(*ljt)[1]] < val) {
									viterbi_m [(int)(*lit)[1]][(int)(*ljt)[1]] = val;
									viterbi_m [(int)(*ljt)[1]][(int)(*lit)[1]] = val;

									int i1,i2;
									char s1,s2;
									i1=(int)(*lit)[0];
									i2=(int)(*ljt)[0];
									s1= g1;
									s2= g2;
//									s2= g2 != 'N' ? g2 : int_to_g(k);

									int b1=0;
									int b2=1;
								//	if(back_track != 0 && p_g1 != 'N' && genotypes[0] != p_g1){
									if(back_track != 0 && ((p_g1 != 'N' & genotypes[0] != p_g1) | (p_g2!= 'N' & genotypes[1] !=p_g2))){
										b1=1;
										b2=0;
									}

									viterbi_index [(int)(*lit)[1]][(int)(*ljt)[1]] = {i1,i2,b1,b2};
									viterbi_state [(int)(*lit)[1]][(int)(*ljt)[1]] = {s1,s2};
									if((int)(*lit)[1] != (int)(*ljt)[1]){
									       viterbi_index [(int)(*ljt)[1]][(int)(*lit)[1]] = {i1,i2,b2,b1};
										viterbi_state [(int)(*ljt)[1]][(int)(*lit)[1]] = {s2,s1};
									}
								}

							}

						}
					}
				} 
			}
			}
			if(g2 != 'N'){
                        for( list<vector<double>> ::iterator lit = (*it)[g_to_int(g2)].begin(); lit != (*it)[g_to_int(g2)].end(); lit++){
                                for( int k=0; k<4; k++){
                                      if( g1 == 'N' || k != g_to_int(g1) ){
	//				if(1){
                                               for( list <vector<double>>::iterator ljt = (*it)[k].begin(); ljt != (*it)[k].end(); ljt++){
                                                       bool b = prev_viterbi_m [(int)(*ljt)[0]][(int)(*lit)[0]] != - numeric_limits<double>::infinity() ;
                                                        vector<char> genotypes;
                                                       if(b){
                                                                genotypes = local_back (back_track,  &viterbi_back_indices, & viterbi_states, {(int)(*ljt)[0],(int)(*lit)[0]});
                                                                if( back_track !=0){
                                                                        b = b & (p_g1 == 'N' | genotypes[0] == p_g1) & (p_g2 == 'N' | genotypes[1] == p_g2);
                                                                        if( (int)(*lit)[0] == (int)(*ljt)[0] ){
                                                                                b=true;
                                                                        }
                                                                }
                                                        }
							if(b){
							       double val = prev_viterbi_m [(int)(*ljt)[0]][(int)(*lit)[0]] + log((*lit)[2])+log((*ljt)[2]);
								if( viterbi_m [(int)(*ljt)[1]][(int)(*lit)[1]] < val) {
									viterbi_m [(int)(*lit)[1]][(int)(*ljt)[1]] = val;
									viterbi_m [(int)(*ljt)[1]][(int)(*lit)[1]] = val;

									int i1,i2;
									char s1,s2;
									i1=(int)(*ljt)[0];
									i2=(int)(*lit)[0];
									//s1='N';
//									s1= g1 != 'N' ? g1 : int_to_g(k);
									s1=g1;
									s2=g2;

                                                                        int b1=0;
                                                                        int b2=1;
                                                                       // if(back_track != 0 && p_g1 != 'N' && genotypes[0] != p_g1){
									if(back_track != 0 && ((p_g1 != 'N' & genotypes[0] != p_g1) | (p_g2!= 'N' & genotypes[1] !=p_g2))){

                                                                                b1=1;
                                                                                b2=0;
                                                                        }

                                                                        viterbi_index [(int)(*ljt)[1]][(int)(*lit)[1]] = {i1,i2,b1,b2};
                                                                        viterbi_state [(int)(*ljt)[1]][(int)(*lit)[1]] = {s1,s2};
                                                                        if((int)(*lit)[1] != (int)(*ljt)[1]){
									       viterbi_index [(int)(*lit)[1]][(int)(*ljt)[1]] = {i1,i2,b2,b1};
										viterbi_state [(int)(*lit)[1]][(int)(*ljt)[1]] = {s2,s1};
                                                                        }

								}
							}
                                                }
                                        }
                                }
                        }
			}
		}else{
			for(int k = 0; k< 4; k++){
				for(int l = 0; l<4; l++){
					for( list<vector<double>> ::iterator lit = (*it)[k].begin(); lit != (*it)[k].end(); lit++){
						for( list <vector<double>>::iterator ljt = (*it)[l].begin(); ljt != (*it)[l].end(); ljt++){
                                                       bool b = prev_viterbi_m [(int)(*lit)[0]][(int)(*ljt)[0]] != - numeric_limits<double>::infinity() ;
                                                        vector<char> genotypes;
                                                       if(b){
                                                                genotypes = local_back (back_track,  &viterbi_back_indices, & viterbi_states, {(int)(*lit)[0],(int)(*ljt)[0]});
					//			cout<<"here0"<<genotypes[0]<<"\t"<<genotypes[1]<<endl;
							//	cout<<"cee"<<genotypes[0] <<"\t"<< p_g1<<endl;
                                                                if( back_track !=0){
                                                                        b = b & (p_g1 == 'N' | genotypes[0] == p_g1) & (p_g2 == 'N' | genotypes[1] == p_g2);
                                                                        if( (int)(*lit)[0] == (int)(*ljt)[0] ){
                                                                                b=true;
                                                                        }
                                                                  }

                                                        }
                                                        if(b){
                                                             double val = prev_viterbi_m [(int)(*lit)[0]][(int)(*ljt)[0]] + log((*lit)[2])+log((*ljt)[2]);
                                                                if( viterbi_m [(int)(*lit)[1]][(int)(*ljt)[1]] < val) {
                                                                        viterbi_m [(int)(*lit)[1]][(int)(*ljt)[1]] = val;
                                                                        viterbi_m [(int)(*ljt)[1]][(int)(*lit)[1]] = val;

                                                                        int i1,i2;
                                                                        char s1,s2;
									i1=(int)(*lit)[0];
									i2=(int)(*ljt)[0];
								//	cout<<"here"<<g1<<"\t"<<g2<<"\t"<<genotypes[0]<<"\t"<<genotypes[1]<<endl;
									s1=g1;
									s2=g2;
//									s1= g1 != 'N' ? g1 : int_to_g(k);
//									s2= g2 != 'N' ? g2 : int_to_g(l);
                                                                        int b1=0;
                                                                        int b2=1;
                                                                       // if(back_track != 0 && p_g1 != 'N' && genotypes[0] != p_g1){
									if(back_track != 0 && ((p_g1 != 'N' & genotypes[0] != p_g1) | (p_g2!= 'N' & genotypes[1] !=p_g2))){

                                                                                b1=1;
                                                                                b2=0;
                                                                        }

                                                                        viterbi_index [(int)(*lit)[1]][(int)(*ljt)[1]] = {i1,i2,b1,b2};
                                                                        viterbi_state [(int)(*lit)[1]][(int)(*ljt)[1]] = {s1,s2};
                                                                        if((int)(*lit)[1] != (int)(*ljt)[1]){
									       viterbi_index [(int)(*ljt)[1]][(int)(*lit)[1]] = {i1,i2,b2,b1};
										viterbi_state [(int)(*ljt)[1]][(int)(*lit)[1]] = {s2,s1};
                                                                        }
                                                                }
							}
						}
					}

				}
			}
		}


/*
		cout<<"viterbi_m"<<"\t"<<coordi<<"\t"<<"gm:"<<genotype_missing<<endl;
		for(int i = 0; i<=(*jt)[1] ; i++){
			for(int j =0 ; j<=(*jt)[1] ; j++){
				cout<<viterbi_m[i][j]<<"\t";
			}
			cout<<endl;
		}
                cout<<"viterbi_index"<<endl;
                for(int i = 0; i<=(*jt)[1] ; i++){
                        for(int j =0 ; j<=(*jt)[1] ; j++){
                                cout<<viterbi_index[i][j][0]<<","<<viterbi_index[i][j][1]<<","<<viterbi_index[i][j][2]<<","<<viterbi_index[i][j][3]<<"\t";
                        }
                        cout<<endl;
                }
                cout<<"viterbi_c"<<endl;
                for(int i = 0; i<=(*jt)[1] ; i++){
                        for(int j =0 ; j<=(*jt)[1] ; j++){
                                cout<<viterbi_state[i][j][0]<<","<<viterbi_state[i][j][1]<<"\t";
                        }
                        cout<<endl;
                }
*/


		it++;
		jt++;
		prev_viterbi_m = viterbi_m;
		viterbi_back_indices.push_back(viterbi_index);
		viterbi_states.push_back(viterbi_state);

	}
	
	double which_max= -numeric_limits<double>::infinity();
	vector<int> which_indices;
	for(int i=0; i<=DAG_node_n.back()[1]; i++){
		for( int j=0; j<= DAG_node_n.back()[1]; j++){
//			cout<<"viter"<<prev_viterbi_m[i][j]<<endl;
			if(which_max < prev_viterbi_m[i][j]){
				which_max=prev_viterbi_m[i][j];
				which_indices = {i,j};
			}
		}
	}
	
//	cout<<"start"<<which_indices[0]<<"\t"<<which_indices[1]<<endl;

	list<vector<int>**>::reverse_iterator rit=  viterbi_back_indices.rbegin(); 
        list<vector<char>**>::reverse_iterator rjt= viterbi_states.rbegin();

	list<vector<char>> reverse_states;

	int b1=0;
	int b2=1;
	for( rit = viterbi_back_indices.rbegin(); rit != viterbi_back_indices.rend(); rit++){
		int c_b1=(*rit)[which_indices[0]][which_indices[1]][ b1 == 1 ? 3: 2];
		int c_b2=(*rit)[which_indices[0]][which_indices[1]][b2 == 1? 3: 2];
		
		reverse_states.push_front({(*rjt)[which_indices[0]][which_indices[1]][b1], (*rjt)[which_indices[0]][which_indices[1]][b2]});
                //cout<<"debug"<<(*rjt)[which_indices[0]][which_indices[1]][b1]<<"\t"<<(*rjt)[which_indices[0]][which_indices[1]][b2]<<"\t"<<(*rit)[which_indices[0]][which_indices[1]][0]<<"\t"<<(*rit)[which_indices[0]][which_indices[1]][1]<<"\t"<<c_b1<<"\t"<<c_b2<<endl;
		which_indices = {(*rit)[which_indices[0]][which_indices[1]][0] ,(*rit)[which_indices[0]][which_indices[1]][1]};

		
		b1= c_b1;
		b2=c_b2;

		rjt++;
	}
	cout<<"hap1\thap2"<<endl;
	for ( list<vector<char>>::iterator it = reverse_states.begin(); it != reverse_states.end(); it++){
		cout<< (*it)[0]<<"\t" <<(*it)[1]<<endl;
	}

	return 0;
}

int g_to_int (char g){
	switch (g){
		case 'A':return 0;
		case 'C':return 1;
		case 'G':return 2;
		case 'T':return 3;
	}
}

char int_to_g (int g){
	switch (g){
               case 0:return 'A';
                case 1:return 'C';
                case 2:return 'G';
                case 3:return 'T';
	}
}


vector<char> local_back ( int back_track, list<vector<int>**>*  viterbi_back_indices,  list<vector<char>**>*  viterbi_states, vector<int> which_indices){
	
       list<vector<int>**>::reverse_iterator rit=  viterbi_back_indices->rbegin();
       list<vector<char>**>::reverse_iterator rjt= viterbi_states->rbegin();
	char a1;
	char a2;
	int b1=0;
	int b2=1;
	while(back_track > 0){
	//	reverse = (*rjt)[which_indices[0]][which_indices[1]][2];
		int c_b1 = (*rit)[which_indices[0]][which_indices[1]][b1 == 0 ? 2 : 3];
		int c_b2= (*rit)[which_indices[0]][which_indices[1]][b2 == 0 ? 2 : 3];

		a1=(*rjt)[which_indices[0]][which_indices[1]][b1];
		a2=(*rjt)[which_indices[0]][which_indices[1]][b2];

//                cout<<"local"<<a1<<"\t"<<which_indices[0]<<"\t"<< which_indices[1]<<"\t"<<b1<<endl;

		which_indices = {(*rit)[which_indices[0]][which_indices[1]][0] ,(*rit)[which_indices[0]][which_indices[1]][1]};
		b1= c_b1;
		b2= c_b2;
//		cout<<c_b1<<endl;

		
		rit++;
		rjt++;
		back_track--;
	}
	return {a1,a2};
}

