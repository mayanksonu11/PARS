#include <bits/stdc++.h>
#include <stdlib.h>
#include <stdio.h>
#include <chrono>
#include <thread>
#include <fstream>
#include <random>
#include "header14.hpp"

using namespace std;
typedef long ll;
// #define double double
// #define MAX 100000
// #define R 100   // Connecitvity range definition
// int cost1,dist1,n;    // Global varaibles utilised in the code
// int flag_global=0,calls=0,max_length=0;
// double wt,throughput;
// double link_bw = 10000; // in MHz
// // double cpu_freq =8000; //in MHz 
// double capacity=8000;
// int k=8,maxVNFs=8; //k is types of vnf
// int numOfUsers;



//367910
// 1. Tangible Number of SFC request
// 2. VNF interference will be the part of core network and final throughput will be minimum of core network and RAN 

/*
	arg[1] is number of base stations
	arg[2] is seed 
	arg[3] is fraction to be reserved for Latency constrained users 
	arg[4] is number of users
*/

int main(int argc, char** argv){
	//=======Development of Core Archirtecture=============
	wt= 0; //to ensure execution of minimum utilised node in the algorithm
	numOfUsers = stoi(argv[4]); 
	n=stoi(argv[1]); //n represents number of nodes in the network 
    int seed=stoi(argv[2]); // Seed for the rand()
    srand(seed);
	double part = stod(argv[3]);  

	vector<pair<int,int>> coord(n);
	for(int i=0;i<n;i++){
		coord[i].first =rand()%250;
		coord[i].second =rand()%250;
	}

    vector<vector<int>> link(n,vector<int>(n,MAX));
	vector<vector<int>>dist(n,vector<int>(n,0));
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			dist[i][j] = ceil(sqrt(pow(coord[i].first-coord[j].first,2)+pow(coord[i].second-coord[j].second,2))); // Calculation of Distance of each node from every node
			// if(dist[i][j]<R){
			// 	link[i][j]=1;
			// 	link[j][i]=1;
			// }
		}		
	}
    int edge_count=0;
    // link[0][1]=1;link[0][2]=1;link[0][5]=1;
    // link[1][0]=1;link[1][2]=1;link[1][3]=1;
    // link[2][0]=1;link[2][1]=1;link[2][8]=1;
    // link[3][1]=1;link[3][6]=1;link[3][4]=1;link[3][13]=1;
    // link[4][3]=1;link[4][9]=1;
    // link[5][0]=1;link[5][10]=1;link[5][6]=1;
    // link[6][3]=1;link[6][5]=1;link[6][7]=1;
    // link[7][6]=1;link[7][8]=1;                                              //NSFNET Topology
    // link[8][7]=1;link[8][2]=1;link[8][9]=1;
    // link[9][8]=1;link[9][11]=1;link[9][12]=1;link[9][4]=1;
    // link[10][5]=1;link[10][11]=1;link[10][12]=1;
    // link[11][10]=1;link[11][9]=1;link[11][13]=1;
    // link[12][10]=1;link[12][13]=1;link[12][9]=1;
    // link[13][3]=1;link[13][11]=1;link[13][12]=1;

    link[0][1]=1;
	link[0][2]=1;
	link[1][0]=1;
	link[1][3]=1;
	link[1][6]=1;
	link[1][9]=1;
	link[2][4]=1;
	link[2][3]=1;
	link[2][0]=1;
	link[3][5]=1;
	link[3][6]=1;
	link[3][1]=1;
	link[3][2]=1;
	link[4][7]=1;
	link[4][2]=1;
	link[5][8]=1;
	link[5][3]=1;
	link[6][3]=1;
	link[6][8]=1;
	link[6][9]=1;
	link[6][1]=1;
	link[7][11]=1;
	link[7][8]=1;
	link[7][4]=1;
	link[8][20]=1;
	link[8][17]=1;
	link[8][18]=1;
	link[8][12]=1;
	link[8][6]=1;
	link[8][5]=1;
	link[8][7]=1;
	link[8][11]=1;
	link[9][6]=1;
	link[9][12]=1;
	link[9][13]=1;
	link[9][10]=1;
	link[9][1]=1;
	link[10][9]=1;
	link[10][13]=1;
	link[11][14]=1;
	link[11][20]=1;
	link[11][8]=1;
	link[11][7]=1;
	link[12][8]=1;
	link[12][21]=1;
	link[12][19]=1;
	link[12][13]=1;
	link[12][9]=1;
	link[13][12]=1;
	link[13][10]=1;
	link[13][9]=1;
	link[14][15]=1;
	link[14][11]=1;
	link[15][16]=1;
	link[15][14]=1;
	link[16][17]=1;
	link[16][15]=1;
	link[17][16]=1;
	link[17][18]=1;
	link[17][8]=1;
	link[18][17]=1;
	link[18][21]=1;
	link[18][8]=1;
	link[19][23]=1;
	link[19][12]=1;
	link[20][11]=1;
	link[20][8]=1;
	link[21][18]=1;
	link[21][22]=1;
	link[21][12]=1;
	link[22][23]=1;
	link[22][21]=1;
	link[23][19]=1;
	link[23][22]=1;

    for(int i=0;i<n;i++){
		for(int j=i+1;j<n;j++){
			if(link[i][j]==1){
				edge_count++;
			}
		}
	}
	cout << n << " "<< edge_count << endl;
	cout << (edge_count+4*n+edge_count*log(n)) << endl;
	// cin >> flag_global;

	vector<int> f[k];  // Contains information of nodes where this vnf could be found;
	vector<int> node[n];
	// cout << f[0].size() << " " << f[1].size() << " " << f[2].size() << " " << f[3].size() << " ";
	vector<vector<double>> per_node_vnf_resource (n,vector<double>(k,0));
	vector<vector<double>> per_node_vnf_resource_slice0 (n,vector<double>(k,0));
	vector<vector<double>> per_node_vnf_resource_slice1 (n,vector<double>(k,0));
	vector<vector<double>> link_ut(n,vector<double>(n,MAX));
	  
	default_random_engine generator;
  	poisson_distribution<int> distribution(0.5);


  	for(int i=0;i<n;i++){
  		int number = distribution(generator);
  		if(number>0){
  			node[i].push_back(0);
  			f[0].push_back(i);
  		}
  	}
	
	for(int i=0;i<n;i++){
		// for(int j=1;j<maxVNFs;j++){
		while(node[i].size()  < maxVNFs){
			int temp_vnf = rand()%k;
			if(find(node[i].begin(), node[i].end(),temp_vnf) == node[i].end()){
				node[i].push_back(temp_vnf);
				f[temp_vnf].push_back(i);
			}
		}
	}
	

	// Concept of resource sharing to be used with PARS

	for(int i=0;i<n;i++){
		int div = node[i].size();
		for(auto v:node[i]){
			if(v==0){
				per_node_vnf_resource_slice0[i][v]= part*1*capacity/div;
				per_node_vnf_resource_slice1[i][v]= (1-part)*1*capacity/div;
			}
			else{
				per_node_vnf_resource_slice0[i][v]= part*capacity/div; 
				per_node_vnf_resource_slice1[i][v]= (1-part)*capacity/div; 
			}
			per_node_vnf_resource[i][v] = per_node_vnf_resource_slice0[i][v]+per_node_vnf_resource_slice1[i][v];
		}
	}
	
	/*
		Uptil now 
			- We have placed points randomly in 250x250 area, these points act as nodes in the network
			- Based on the distance between them we have considered them to be connected or not. So if the distance between any two nodes is less than R (100) it is considered to be connected
			- We have placed VNF0 (VNF connnected to base station) randomly on to the nodes, it is performed using Posson distribution with lambda 0.5 
			- All the VNFs except VNF0 are placed on all the nodes
			- In the last step we equally divide the node resources among all the VNFs
		This completes the creation of an Edge Network Scenario, where there are nodes, which are having some resources, and these resources are divided equally among all the VNFs hosted on it. Links decides which nodes are directly connected and each link has some capacity.
	*/

		
	vector<vector<double>> node_vnf_resources= per_node_vnf_resource;
	vector<vector<double>> node_vnf_resources_slice0= per_node_vnf_resource_slice0;
	vector<vector<double>> node_vnf_resources_slice1= per_node_vnf_resource_slice1;
	vector<int> adj[n];
    for(int i=0,k;i<n;i++){
        for(int j=i;j<n;j++){
            if(link[i][j]==1) {
            	adj[i].push_back(j);
            	adj[j].push_back(i);
                link_ut[i][j] = 0;
                link_ut[j][i] = 0;
            }
        }
    }
    saveLinkInfo(link_ut,"link_res");
    displayNodeUt(per_node_vnf_resource);
    displayNodeUt(per_node_vnf_resource_slice0);
    displayNodeUt(per_node_vnf_resource_slice1);
	// cin >> flag_global;
    printGraph(adj);
    /*
		The above section of code contains creation of adjacency matrix for the graph. Along with that it also sets current utilisation value of all the links at zero. Moreover, the variable node_vnf_resources represents the data of amount of resources allocated to the VNFs in each nodes. Further it saves the current link utilisation in a file link_res.txt and also displays the current remaining resources at the node. Both in generalised sense and slice specific.
    */

    vector<vector<int>> reachable(n,vector<int>(n,0));
    vector<vector<vector<int>>> calculatedRoute(n,vector<vector<int>>(n)); //Contains all the routes from source to destination node
    for(int i=0;i<n;i++){  // Loop to generate route from one node to another using Dijkstra Algorithm based on physical location and connectivity
    	for(int j=i+1;j<n;j++){
            vector<int> visited(n,0);           
            reachable[i][j]= isReachable(i,j,adj);
            reachable[j][i] = reachable[i][j];
            if(reachable[i][j]){
            	calculatedRoute[i][j] = findPathDijkstra(j,i,link,visited);
            	calculatedRoute[j][i] = reverseArray(calculatedRoute[i][j]);
            	max_length = max(max_length,(int)calculatedRoute[i][j].size());
            }
    	}
    }   
    displayCost(reachable); 
	cout << "Max Path length:" << max_length << endl;
	// cin >> flag_global;
	/*
		Above section evaluates path from each node to every other node in the graph. 
	*/
    int numOfSlices=2;
	vector<pair<int,int>> service_type(numOfSlices);
	service_type[0] = {latency_slice0,rate_slice0};	
	service_type[1] = {latency_slice1,rate_slice1}; 	
	vector<double> service_charge(numOfSlices);
	service_charge[0] = 10000;
	service_charge[1] = 20000;
	vector<double> per_user_ut_node(k,0);
	for(int i=0;i<k;i++){
		per_user_ut_node[i]= 1*(rand()%4+1);
	}
	per_user_ut_node[0]= 0.1*(rand()%4+1);

	/*
		In the above section,
			- We have initialised the slice specific requirement parameters in terms of latency and data rate
			- We have also initialised resource requirement of VNFs but its not specific to the slice
	*/


	//=======Following section contains code for SFC generation for argv[4] number of users=========

	int  count1=0,count2=0;
	vector<pair<double,double>> per_user_req;
	vector<vector<int>> per_user_sfc;
	vector<int> per_user_slice;
	for(int key=0;key<stoi(argv[4]);key++){
		// Users' slice parameter
		int slice = rand()%30;
		if(slice<10) slice =1;
		else slice =0;
		per_user_slice.push_back(slice);
		
		// User specific data-rate and latency requirement
		double rate, latency;
		if(slice==0){ 
			count1++;
			latency=service_type[slice].first; 
			latency -= (rand()%50)/100.0;
			rate = service_type[slice].second;
		}
		else if(slice==1) {
			count2++;
			latency=service_type[slice].first;
			rate = service_type[slice].second;
			rate += rand()%(int)(rate/2);	
		}
		per_user_req.push_back({latency,rate});

		// SFC creation for individual users
		if(chain_length==0) chain_length++;
		vector<int> sfc_chain;
		// cout << "Iteration Count= " << key << endl;
		unordered_set<int> temp_chain;
		while(temp_chain.size() != chain_length){
			int t_ = rand()%k;
			t_ = max(t_,1);
			temp_chain.insert(t_);
		}
		sfc_chain.push_back(0);
		for(auto v:temp_chain){
		    sfc_chain.push_back(v);
		}
		per_user_sfc.push_back(sfc_chain);
	}
	// printf("Total SFCs:%ld\n",per_user_sfc.size());
	func1(per_user_req,per_user_sfc,per_user_slice,service_type,service_charge,per_node_vnf_resource,per_node_vnf_resource_slice0,per_node_vnf_resource_slice1,per_user_ut_node,reachable,calculatedRoute,link_ut,f,1,model);
}

	/*
	Whether user gets connected or not
	Amount of link utilisation for each user
	Exact of processing resource allocated 
	Reward associated with each user

	1. varied rate and latency users
	2. classification of users base on rate and latency --> we take it from user request itself
	3. reward functions for users of different slices
	4. Enforce maximum percentage of processor utilised by URLLC users for each VNFs in each Nodes
	5. Seperate VNFs per slice / without VNF sharing 
	6. Implement standard topology with respective
	7. Remove current VNF placement setup, VNF setup and resource allocation based on arrival of traffic
	8. Comparison with some baseline algorithm
	9. Algorithm pseudocode
	10. Plots:
		 =============== Analysis wrt baseline algorithm =======================
		* Execution time, -- line graph wrt varying number of sfc requests 
		* reward varying number of sfc requests
		* cost varying with number of sfc requests
		* acceptance ratio vs number of requests
		==================== VNF sharing and non-sharing difference ==============
		* fixing number of sfc requests plot the variation of reward with number of nodes for both the algorithms -- bar graph
		* Compare the above results with VNF sharing and without VNF sharing *** 
		* VNF utilisation improvement with addition of VNF sharing  vs  varied ratio of eMBB and URLLC traffic (embb=1 and vary urllc) --2 line graphs
	10. 
	*/