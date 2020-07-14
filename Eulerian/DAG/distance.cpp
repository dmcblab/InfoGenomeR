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
#include <regex>
int index(int gv_i , int gv_j, int gv_size);
int orientation ( int ori1, int ori2, bool swapped);
int ori_to_index (int ori);
using namespace std;
int main(){

       map<double,double> node_to_coordinate;
	map<double,double>::iterator mit= node_to_coordinate.begin();

        ifstream file;
	file.open("node_keys");
	//file.open("node_keys");
	
	string line;
	getline(file,line,'\n');
        while(getline(file,line,'\n')){
                stringstream f_iss(line);
                string val;
		double temp[4];
		int temp_j=0;
                while(getline(f_iss, val, '\t')){
                        stringstream convertor(val);
                        double t;
                        convertor>>t;
                        temp[temp_j]=t;
                        temp_j++;
                }
 	        node_to_coordinate.insert(mit, pair<double, double>(temp[3],temp[2]));

	}
	
	file.close();
	file.clear();


	file.open("node_keys_genes.recurrent.unique");
	vector<vector<int>> gv;
        while(getline(file,line,'\n')){
                stringstream f_iss(line);
                string val;
		int temp[3];
		int temp_j=0;
                while(getline(f_iss, val, '\t')){
                        stringstream convertor(val);
                        int t;
                        convertor>>t;
			temp[temp_j]=t;
			temp_j++;
		}
		gv.push_back({temp[1],temp[2]});
	}
	file.close();
	file.clear();


	file.open("euler_paths");
//	file.open("v");
	regex path("(gnodes)(.*)");
	bool distance_begin=false;

	vector<vector<double>> total_sum;
        vector<vector<int>> total_CN;
	for(int v=0; v<(gv.size()*(gv.size()-1)/2); v++){
		total_sum.push_back({0,0,0,0});
		total_CN.push_back({0,0,0,0});
	}
        while(getline(file,line,'\n')){
		double sum=0;
		distance_begin=false;

		if(regex_match(line,path)){
			cout<<"gnodes"<<endl;
		//	cout<<"onepath"<<endl;
		//	cout<<total_CN<<endl;
	//		if(total_CN[index(12,57,gv.size())][3]>0)
	//			cout<<((double)total_sum[index(12,57,gv.size())][3]/(double)total_CN[index(12,57,gv.size())][3])<<endl;
			for(int v=0; v<total_sum.size(); v++){
			//	total_sum[v]={0,0,0,0};
			//	total_CN[v]={0,0,0,0};
			}
		}else{
			stringstream f_iss(line);
			string val;
	//		unsigned int temp[2];
	//		unsigned int previous=0;
			vector<int> path;
			while(getline(f_iss, val, '\t')){
				stringstream convertor(val);
				int t;
				convertor>>t;
				path.push_back(t);
/*
				if(distance_begin && previous == 35 && t == 34){
					cout<<sum<<endl;
					total_sum+=sum;
					total_CN++;
					break;
				}

				if(distance_begin){
					if((previous%2)==0 && previous+1== t)G
						sum+=(node_to_coordinate.find(t)->second - node_to_coordinate.find(previous)->second);
					if((t%2)==0 && t+1 == previous)
                                                sum+=(node_to_coordinate.find(previous)->second - node_to_coordinate.find(t)->second);
				}
				if(previous== 26 && t==27){
	//				cout<<"who"<<endl;
					distance_begin=true;
					sum=0;
				}
				previous=t;
*/
			}
			vector<vector<int>> genes[gv.size()];
			int j=1;
			for(int i=0; i<path.size(); i+=2){
				for(int v=0; v<gv.size(); v++){
					if(path[i] == gv[v][0] && path[j] == gv[v][1]){
						genes[v].push_back({i,j,0});
					}else if(path[i] == gv[v][1] && path[j] == gv[v][0]){
						genes[v].push_back({i,j,1});
					}
				}
				j+=2;
			}

			double min_dis;
			double c_sum;
			int min_ori;
			//cout<<"size"<<genes[0].size()<<"\t"<<genes[2].size()<<endl;
			for(int gv_i=0; gv_i < gv.size(); gv_i++){
				for(int gv_j=gv_i+1; gv_j < gv.size(); gv_j++){
				c_sum=0;
				min_dis=numeric_limits<double>::infinity();

					for(int i =0; i<genes[gv_i].size(); i++){
						for(j=0; j<genes[gv_j].size(); j++){
							c_sum=0;
							int coor1= genes[gv_i][i][0];
							int coor2= genes[gv_j][j][0];
							int ori1=genes[gv_i][i][2];
							int ori2=genes[gv_j][j][2];
							int swapped= -1;
							if(coor1 > coor2){
								int swapped = coor1;
								coor1=coor2;
								coor2=swapped;
							}
							int l=coor1+3;
						//	cout<<"hs"<<c_sum<<endl;
							for(int k=coor1+2; k<coor2; k+=2){
								c_sum+=abs(node_to_coordinate.find(path[k])->second - node_to_coordinate.find(path[l])->second);
								l+=2;
							}

							if(min_dis > c_sum){
								min_dis=c_sum;
								min_ori = orientation(ori1, ori2, swapped> 0?  true : false);
							}
						}
					}

	//				if(genes[gv_i].size()>0 && genes[gv_j].size() >0){
                                        if(genes[gv_i].size()>0 && genes[gv_j].size() >0 && min_dis < 1e9){

						total_sum[index(gv_i,gv_j,gv.size())][ori_to_index(min_ori)]+=min_dis;
						total_CN[index(gv_i,gv_j,gv.size())][ori_to_index(min_ori)]++;
					}
				}
			}

		}


	}
 
        ofstream outfile;
        outfile.open("distance.output");

	for(int i=0; i<total_CN.size(); i++){
		for(int j=0; j<4;j++){
			if(total_sum[i][j]>0)
				total_sum[i][j]=total_sum[i][j]/total_CN[i][j];
			outfile<<total_sum[i][j]<<"\t";
		}
		outfile<<endl;
	}
	ofstream outfile2;
	outfile2.open("distance.output.CN");
        for(int i=0; i<total_CN.size(); i++){
                for(int j=0; j<4;j++){
                        outfile2<<total_CN[i][j]<<"\t";
                }
                outfile2<<endl;
        }


}
int ori_to_index (int ori){
	if(ori == 55)
		return 0;
	else if(ori == 53)
		return 1;
	else if(ori==35)
		return 2;
	else if(ori=33)
		return 3;

}
int orientation ( int ori1, int ori2, bool swapped){
	int min_ori;
	if(swapped){
                if(ori1 == 0 && ori2 == 0)
                        min_ori=53;
                else if(ori1 == 0 && ori2 ==1)
                        min_ori=55;
                else if(ori1 == 1 && ori2 == 0)
                        min_ori=33;
                else
                        min_ori=35;

	}else{
                if(ori1 == 0 && ori2 == 0)
			min_ori=35;
                else if(ori1 == 0 && ori2 ==1)
			min_ori=33;
		else if(ori1 == 1 && ori2 == 0)
			min_ori=55;
		else
			min_ori=53;
	}
	return min_ori;
}

int index(int gv_i , int gv_j, int gv_size){
	int sum=0;
	for(int i=0; i<gv_i; i++){
		sum+=gv_size-i-1;
	}
	sum+=gv_j-gv_i-1;
	return sum;
}
