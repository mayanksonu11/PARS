#include <bits/stdc++.h>
#include <stdlib.h>
#include <stdio.h>
#include <chrono>
#include <thread>
#include <fstream>
#include <random>

using namespace std;
// #define double double
#define MAX 100000
#define R 70   // Connecitvity range definition
typedef long ll;
int cost1,dist1,n;    // Global varaibles utilised in the code
int flag_global=0,calls=0,max_length=0;
double wt;
double link_bw = 10000; // in MHz
// double cpu_freq =8000; //in MHz 
double capacity=236;


void findRoute(int s, int f, vector<vector<int>> dist,vector<vector<int>> routeCost,vector<vector<int>> adj[],vector<vector<int>>cost,vector<int> &route,vector<vector<vector<int>>> calculatedRoute);

template <typename t>
void display(vector<t> &a){
    for(int i=0;i<a.size();i++){
        cout << a[i] << " ";
    }
    cout << "\n";
}

void saveLinkInfo(vector<vector<double>> link_ut, string str){
	ofstream fout;
	fout.open(str+".txt");
	fout << 50 << " ";
	for(int i=0;i<link_ut.size();i++){
		fout << i+1 << " ";
	}
	fout << endl;
	for(int i=0;i<link_ut.size();i++){
		fout << i+1 << " ";
		for(int j=0;j<link_ut[0].size();j++){
			if(link_ut[i][j]<=link_bw)
				fout << link_ut[i][j]/link_bw << " ";
			else
				fout << "-1" << " ";	
		} 
		fout << endl;
	}
	fout.close();
}


void saveNodeInfo(vector<vector<double>> node_ut,vector<vector<double>> node_res, string str){
	ofstream fout;
	fout.open(str+".txt");
	fout << node_ut[0].size() << " ";
	for(int i=0;i<node_ut[0].size();i++){
		fout << i+1 << " ";
	}
	fout << endl;
	for(int i=0;i<node_ut.size();i++){
		fout << i+1 << " ";
		for(int j=0;j<node_ut[0].size();j++){
			if(node_res[i][j]>0.5)fout << (node_res[i][j]-node_ut[i][j])/node_res[i][j] << " ";
			else fout << 0 << " ";	
		} 
		fout << endl;
	}
	fout.close();
}

void displayLinkUt(vector<vector<double>> link_ut){ // To display distance matrix
	for(int i=0;i<link_ut.size();i++){
		for(int j=0;j<link_ut[0].size();j++){
			if(link_ut[i][j]>= link_bw){
                cout << " *  " ;
            } 
            else if(link_ut[i][j]==0)
            	cout << " = ";
			else cout<< link_ut[i][j] << " ";	
		} 
		cout << endl;
	}
	cout << endl;
}

void displayNodeUt(vector<vector<double>> node_ut){
	for(int i=0;i<node_ut.size();i++){
		for(int j=0;j<node_ut[0].size();j++){
			if(node_ut[i][j]==0)
            	cout << " = ";
			else cout<< node_ut[i][j] << " ";	
		} 
		cout << endl;
	}
	cout << endl;
}
void printGraph(vector<int> adj[]) //To display adjacency list 
{
    for (int v = 0; v < n; ++v)
    {
        cout << "\nAdjacency list of vertex "
             << v << "\nhead ";
        // for (auto x : adj[v])
        for(int i=0;i<adj[v].size();i++)		
           cout << "-> " << adj[v][i];
        printf("\n");
    }
}
void displayCost(vector<vector<int>> cost){ // To display cost matrix
	cout <<"Node   ";
	for(int i=0;i<cost.size();i++)
        cout << i << " ";
    cout << endl;
	for(int i=0;i<cost.size();i++){
        cout <<"Node "<< i << " ";
        for(int j=0;j<cost[0].size();j++){
            if(cost[i][j]>= MAX){
                cout << "* " ;
                continue; 
            } 
            cout << cost[i][j] << " ";
        }cout << endl;
    }
    cout << endl;
}
void displayRoute(vector<int> route){ // To display the route
	int i;
	for(i=0;i<route.size()-1;i++){
		cout << route[i] << " -> ";
	}
	cout << route[i] << endl;
}
double findCost(vector<int>route,vector<vector<double>> cost){ // Given a route in the form of vector, used to find out the total cost 
	double sum=0;
	for(int i=0;i<route.size()-1;i++){
		sum+= link_bw -cost[route[i]][route[i+1]];
		// cout  << cost[route[i]][route[i+1]]<< "--";
	}
	return sum;
}
double std_dev(vector<int>route, vector<vector<double>>&cost){
	double avg=0;
	int k =route.size()-1;
	for(int i=0;i<k;i++){
		avg+= cost[route[i]][route[i+1]];
	}
	avg/=(route.size());
	// cout<< avg << endl;
	double var=0;
	for(int i=0;i<k;i++){
	    var+= (cost[route[i]][route[i+1]]-avg)*(cost[route[i]][route[i+1]]-avg);
	}
	var/=k;
	double stdDev = sqrt(var);
	// cout << "StdDev "<< stdDev << endl;
	return stdDev;
}

void findBestRoute(vector<int>& route,vector<vector<double>> cost, vector<int>adj[]){
	//Used to find out the best route from given raw path based on distance
	//By changing and locking the nodes so as to find the optimal path
	calls++;
	vector<int> lock(route.size(),0);
	lock[0]=1; *(lock.rbegin())=1;
	// display(lock); 
	double cost_min = std_dev(route,cost);
	double routeCost= findCost(route,cost);
	// display(route);
	// cout << "cost_min=" << cost_min << " routeCost="<< findCost(route,cost) <<  endl;
	for(int i=0;i<lock.size();i++){
		int node_min=-1;
		if(lock[i]==0){
			int def= route[i];
			for(int j=0;j<adj[route[i-1]].size();j++){
				route[i]= adj[route[i-1]][j];
				double updatedCost = std_dev(route,cost);
				if(findCost(route,cost)<2*routeCost and updatedCost<cost_min){
					// cout<<"I have executed with " << updatedCost << endl; 
					flag_global++;
					node_min = adj[route[i-1]][j];
					cost_min = updatedCost;
				}
			}
			if(node_min>=0) {
				route[i]=node_min;
			}else{
				route[i]=def;
			}
			lock[i]=1;
		}
	}
	// display(route);
	// cout << "cost=" << findCost(route,cost) << endl;
	// cout << "StdDev "<< std_dev(route,cost) << endl;
	// cout << "end of adjustment" << endl;
}
int findPathDFS(int level, int path_len, int s, int f, vector<int> visited,vector<int>&parents ,vector<int> adj[]){
	int flag=0;
	if(level > path_len){
		return 0;
	}
	if(s==f){
		parents[f] = s;
		return 1;
	}

	for(int i=0;i<adj[s].size();i++){
		if(visited[adj[s][i]]==0){
			visited[adj[s][i]]=1;
			if(findPathDFS(level+1,path_len,adj[s][i],f,visited,parents,adj)){
				parents[adj[s][i]]=s;
				return 1;
			}
		}

	}

	return 0;

}
template <typename t>
vector<int> findPathDijkstra(int s,int f, vector<vector<t>> cost , vector <int> &visited){ // Implementation of Dijkstra Algorithm
    visited[s]=1;
    vector<t> dist = cost[s];
    vector<int> parents(dist.size(), s);
    int flag=0;
    vector<int> route;
    // cout << "source is:"<< s << endl;
    while (visited[f]!=1)
    {
        int k, minIndex=-1;
        t mini = MAX;
        for (int i = 0; i < dist.size(); i++)
        {
            if (visited[i] == 1)
                continue;
            if (dist[i] <= mini){
                mini = dist[i];
                minIndex = i;
            }
        }
        if(visited[f]==1){
        	break;
        }
        // cout << "I am here " << minIndex <<" "<< dist[minIndex] << " " << f<< endl;
        if(minIndex==-1){
        	// return route; 
        }
        
        for (int j = 0; j < dist.size(); j++)
        {
            if (visited[j] == 1)
                continue;
            if (dist[j] > dist[minIndex] + cost[minIndex][j])
                parents[j] = minIndex;
            flag = 1;
            dist[j] = min(dist[j], dist[minIndex] + cost[minIndex][j]);
            // cout << j << " " << dist[j] << " " << visited[j] <<  endl;
        }
        visited[minIndex] = 1;
    }
    if(flag==0){
    	return route;
    }
    for(int j=f;j!=s;){
        // cout << j << " <- ";
        route.push_back(j);
        j = parents[j];
    }
    route.push_back(s);
    // cout << s << " and cost= " << dist[f] << endl;
    return route;
}

vector<int> reverseArray(vector<int> arr){ // To reverse a given array
	int len = arr.size();
	vector<int> temp;
	for(int i=len-1;i>=0;i--){
		temp.push_back(arr[i]);
	}
	return temp;
}

// bool check(vector<vector<double>> &per_node_vnf_resource,vector<vector<double>> &link_ut, int k, int n){ //k is type of vnfs
// 	double node_ut=0;
// 	int count=0;
// 	for(int i=0;i<per_node_vnf_resource.size();i++){
// 		for(int j=0;j<k;j++){
// 			if(per_node_vnf_resource[i][j]>0.9){
// 				count++;
// 			}
// 		}
// 	}
// 	node_ut += count/((double)k*n);
// 	double node_ut = 
// 	int 
// }

double max_utilisation_node(vector<vector<double>> &mat,vector<vector<double>> &resources ){
	double val=-1;
	for(int i=0;i<mat.size();i++){
		for(int j=0;j<mat[0].size();j++){
			if(mat[i][j]!=0){
				val=max(val,(resources[i][j]-mat[i][j])/resources[i][j]*100);
			}
		}
	}
	cout << "Maximum utilised node status: " <<val << endl; 
	return val;
}

double max_utilisation(vector<vector<double>> &mat){
	double val=-1;
	for(int i=0;i<mat.size();i++){
		for(int j=0;j<mat[0].size();j++){
			if(mat[i][j]<MAX){
				val=max(val,mat[i][j]);
			}
		}
	}
	cout << "Max utilised link status: " <<val << endl; 
	return val;
}


bool notInRoute(vector<int> route, int k){
    for(int i=0;i<route.size();i++) 
    	if(route[i]==k) return false;
    return true;
}

	// node_map[i] = find_min_usage(sfc_chain[i],f[sfc_chain[i]],node_map,per_node_vnf_resource);


int find_min_usage(double req_res,int vnf, vector<int> f, vector<int> node_map ,vector<vector<double>> node_res,vector<vector<double>> res){
	int min_ut_node=-1;
	double max_remaining = -1e9;
	for(int i=0;i<f.size();i++){
		if(node_res[f[i]][vnf]*100/res[f[i]][vnf]>=max_remaining && node_res[f[i]][vnf]>req_res && notInRoute(node_map,f[i]) )
		{
			// cout << i << " ";
			max_remaining = node_res[f[i]][vnf];
			min_ut_node = f[i];
		}
	}
	// cout << endl;
	// display(node_map);
	// cout << "min_ut_node: " << min_ut_node << " "; 
	return min_ut_node;
}

// node_map[i] = find_good_node(req_res,sfc_chain[i],f[sfc_chain[i]],node_map[i+1],node_map,resource,per_node_vnf_resource,link_ut);

int find_good_node(double req_res, int vnf, vector<int> f,int next_node, vector<int> node_map ,vector<vector<double>> &avail_res,vector<vector<double>> &node_res, vector<vector<double>> link_ut){
	int best_node = -1;
	double evaluation = -1;
	for(int i=0;i<f.size();i++){
		if(avail_res[f[i]][vnf]<=req_res){
			continue;
		}
		if(!notInRoute(node_map,f[i])){
			continue;
		}
		// double temp_evaluation = (wt)*calculatedRoute[prev_node][f[i]].size()/(double)max_length+(1-wt)*avail_res[f[i]][vnf]/100.0;
		vector<int> visited(n,0);
		double t1 = findCost(findPathDijkstra(f[i],next_node,link_ut,visited),link_ut);
		// cout << "Cool"; cin >> flag_global;
		double t2 = 0;
		for(int j=0;j<avail_res[f[i]].size();j++){
			t2 += node_res[f[i]][vnf];
		}
		double t3 = avail_res[f[i]][vnf];
		// cout << t1 << " " << t2 << " " << t3 << endl;
		// cout << "Cool"; cin >> flag_global;
		double temp_evaluation = t1+t2+t3;
		// cout << temp_evaluation << "\n";
		if(temp_evaluation>evaluation){
			evaluation = temp_evaluation;
			best_node = f[i];
		}
	}

	return best_node;
}

int find_nearby_node(double req_res,int vnf, vector<int> f, int prev_node, vector<int> node_map,vector<vector<double>> &node_res, vector<vector<vector<int>>> &calculatedRoute){
	int min_dist_node=MAX;
	// display(f);
	// cout << vnf << " " << prev_node << endl; 
	for(int i=0;i<f.size();i++){
		if(calculatedRoute[prev_node][f[i]].size()<min_dist_node && notInRoute(node_map,f[i]) && node_res[f[i]][vnf]>=req_res){
			min_dist_node = f[i];
		}
	}
	return min_dist_node;
}

double average_bw(vector<int> route,vector<vector<double>> &link_ut){
	if(route.size()==0){
		return 0;
	}
	double avg_bw=0;
	for(int i=0;i<route.size()-1;i++){
		avg_bw+= link_ut[route[i]][route[i+1]];
	}
	return avg_bw/(route.size()-1);
}

int find_best_match(double req_res,int vnf, vector<int> f, int prev_node, vector<int> node_map,vector<vector<double>> &node_ut,vector<vector<double>> &link_ut ,vector<vector<vector<int>>> &calculatedRoute){
	vector<double> w(f.size(),0);
	vector<double> avg_bw(f.size(),0);
	double max_node_ut=-1;
	for(int i=0;i<f.size();i++){
		max_node_ut = max(max_node_ut,node_ut[f[i]][vnf]);
	}
	// cout << max_node_ut << " check-1" << endl;
	double max_link_ut=-1;
	for(int i=0;i<f.size();i++){
		avg_bw[i]=average_bw(calculatedRoute[prev_node][f[i]],link_ut);
		max_link_ut= max(max_link_ut,avg_bw[i]);
	}
	// display(avg_bw);
	// cout << max_link_ut << " check-2" << endl;
	for(int i=0;i<f.size();i++){
		if(calculatedRoute[prev_node][f[i]].size()==0)
			w[i]=-1;
		else if(max_link_ut==0)
			w[i]=0;
		else if(max_node_ut<0.05)
			return -1;
		else
			w[i]= (1/(double)calculatedRoute[prev_node][f[i]].size())*(avg_bw[i]/max_link_ut)*(node_ut[f[i]][vnf]/(double)max_node_ut);
		// cout << w[i] << " ";
	}
	// display(w);
	// cout << " check-3" << endl;
	double max_eval=-1;
	int node=-1;
	for(int i=0;i<f.size();i++){
		if(w[i]>max_eval){
			max_eval=w[i];
			node = f[i]; 
		}
	}
	// cout << "check-4" << endl;
	return node;
}

void display_utilisation(int vnf, vector<vector<double>> node_ut){
	for(int i=0;i<node_ut.size();i++){
		cout <<"Node " << i << " " << node_ut[i][vnf] << endl;
	}
	cout << endl;
}

double percentUtilNode(vector<vector<double>>&resources,vector<vector<double>> &node_ut){
	double value=0;
	int count=0;
	for(int i=0;i<node_ut.size();i++){
		for(int j=0;j<node_ut[0].size();j++){
			if((resources[i][j]-node_ut[i][j])>0)cout << "{" << (resources[i][j]-node_ut[i][j]) <<" "<< resources[i][j] << "} ";
			if(node_ut[i][j]!=0 and resources[i][j]!=0){
				value += (resources[i][j]-node_ut[i][j])/resources[i][j];
				count++;
			}
		}
		cout << endl;

	}
	cout << value << " ";
	cout << endl;
	return value/count*100;
}

double percentUtilLink(vector<vector<double>> link_ut){
	double val=0; int count=0;
	for(int i=0;i<link_ut.size();i++){
		for(int j=i+1;j<link_ut[0].size();j++){
			if(link_ut[i][j]<=link_bw){
				val+= link_ut[i][j];
				count++;
			}
		}
	}
	return val/(count*link_bw)*100;
}
int NoOfZeros(vector<int> recent_observation){
	int count_zero=0;
	for(int i=0;i<recent_observation.size();i++){
		if(recent_observation[i]==0)
			count_zero++;
	}
	return count_zero;
}

void addANode(vector<int>&route, vector<int> adj[]){
	int s = route[0]; int f =  route[1];
	vector<int> alt_route;
	alt_route.push_back(s);
	int flag=0;
	for(int i=0;i<adj[s].size();i++){
		for(int j=0;j<adj[s].size();j++){
			if(adj[s][i]==adj[s][j]){
				alt_route.push_back(adj[s][j]);
				flag=1;
				break;
			}
		}
		if(flag)
			break;
	}
	alt_route.push_back(f);
	route = alt_route;
}

double findLinkUt(int rate){
	double bw = rate/log2(1e6);
	// cout << bw << endl;
	return bw; // exact model of data rate to bandwidth translation has to be done for optical cable
	//bw is in MHz if rate is in Mbps
}


double findReqNodeResource(int latency,int rate){
	/* given a rate we need to come up with some heuristic which should tell how much of processor would be required to process that chunk of data every second
			Challenge would be find out how does the other process influence the upcoming one
	*/
	// double total_data = latency/3*rate*1e3; //assuming one-third of latency will be consumed at UE transmission

	return (rate*100)/(capacity*64); //information of VNFs in terms of number of machine cycles required for processing of 1MB of data
	/*instructions executed by VNF*/
}

double delayInNode(int node, double ut){
	// cout << ut << endl;
	// cin >> flag_global;
	double e=1;
	if(ut <=15) return 0.05*e;
	else if(ut >15 and ut <=25) return 1.0*e;
	else if(ut >25 and ut <=35) return 2.0*e;
	else if(ut >35 and ut <=45) return 2.25*e;
	else if(ut >45 and ut <=55) return 2.5*e;
	else if(ut >55 and ut <=65) return 3.0*e;
	else if(ut >65 and ut <=75) return 4.0*e;
	else if(ut >75 and ut <=85) return 5.5*e;
	else if(ut >85 and ut <=95) return 9.5*e;
	else if(ut >95 and ut <=100) return 15*e;
}

double findTotalLatency(vector<int> route, vector<int> sfc, vector<int> node_map,vector<vector<double>> &link_ut,vector<vector<double>> &node_ut,vector<vector<double>> &slice_res ,vector<vector<double>>& node_res){
	/*
		Given the route find out how much is the propagation delay and how much is the processing delay
		it is sum of the propagation delay in 
	*/
	// double latency = /*Transmission time in air and propagation time+ processing time at access points+/+/*SFC Processing time*/
	double tot_lat=0;
	for(int i=0;i<node_map.size();i++){
		if(node_res[node_map[i]][sfc[i]]==0) 
			return MAX;
		double ut = (slice_res[node_map[i]][sfc[i]]-node_ut[node_map[i]][sfc[i]])*100/node_res[node_map[i]][sfc[i]];
		cout << "Utilisation is :" << ut << endl;
		tot_lat+= delayInNode(node_map[i],ut);
		cout << tot_lat << " ";
	}
	cout <<"Latency is:" << tot_lat << endl;
	// latency in the core network links
	//delay in RAN
	//delay in UE transmission

	return tot_lat;
	
}

double findThreshold(vector<double> t){
	double val=0;
	sort(t.begin(),t.end());
	return t[t.size()/2]; // using median as the threshold
}


int findHighNode(vector<int> high_usage, set<int> vnf){
	for(auto v:vnf){
		if(find(high_usage.begin(), high_usage.end(),v)!= high_usage.end()){
			return v;
		}
	}
	return -1;
}

bool isReachable(int s, int d, vector<int> adj[]){
	if(s==d)
		return true;

	vector<bool> visited(n,0);
	queue<int> q;
	visited[s]= true;
	q.push(s);
	while(!q.empty()){
		s = q.front();
		q.pop();
		for(int i = 0;i<adj[s].size();i++){
			if(adj[s][i]==d)
				return true;

			if(!visited[adj[s][i]]){
				visited[adj[s][i]] = true;
				q.push(adj[s][i]);
			}
		}
	}
	return false;
}
//367910
// 1. Tangible Number of SFC request
// 2. VNF interference will be the part of core network and final throughput will be minimum of core network and RAN 

int main(int argc, char** argv){
	//=======Development of Core Archirtecture=============
	// 
	double av_reward=0,av_util=0,av_time=0,av_accRatio=0,av_cost=0,av_throughput=0;
	for(int iter=0;iter<stoi(argv[2]);iter++){
	n=stoi(argv[1]);
	int k=4,maxVNFs=4; //k is types of vnf
    int seed=iter;
	// cout << "Please enter the seed:";
	// cin >> seed;
    srand(seed);
	// Types of VNFs in the network
	// wt = stof(argv[3]);
	wt = 0;
	float part = stod(argv[3]);  
	// cout << wt << endl;
	ofstream fout;
	fout.open("sample.txt");
	// cout << "Please enter the number of nodes:";
	// cin >> n;
		fout << "# x y R" << endl; //Defines the coordinate of each node
	vector<pair<int,int>> coord(n);
	for(int i=0;i<n;i++){
		coord[i].first =rand()%100;
		coord[i].second =rand()%100;
		// cout << i << " " << coord[i].first << "," << coord[i].second << endl; //Defines the coordinate of each node
		fout  << coord[i].first << " " << coord[i].second << " " << i << endl; //Defines the coordinate of each node
	}
	fout << "e" << endl;

	fout.close();

	vector<vector<int>>dist(n,vector<int>(n,0));
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			dist[i][j] = ceil(sqrt(pow(coord[i].first-coord[j].first,2)+pow(coord[i].second-coord[j].second,2))); // Calculation of Distance of each node from every node
		}		
	}
    vector<vector<int>> link(n,vector<int>(n,MAX));
    int edge_count=0;
    for(int i=0;i<n;i++){
		for(int j=i+1;j<n;j++){
			if(dist[i][j]<=R){
				link[i][j]=1; // Defining available links between nodes in VM
				link[j][i]=1;
				edge_count++;
			}
		}
	}
    // link[0][1]=1;link[0][2]=1;link[0][5]=1;
    // link[1][0]=1;link[1][2]=1;link[1][3]=1;
    // link[2][0]=1;link[2][1]=1;link[2][8]=1;
    // link[3][1]=1;link[3][6]=1;link[3][4]=1;
    // link[4][3]=1;link[4][9]=1;link[3][13]=1;
    // link[5][0]=1;link[5][10]=1;link[5][6]=1;
    // link[6][3]=1;link[6][5]=1;link[6][7]=1;
    // link[7][6]=1;link[7][8]=1;link[9][4]=1;
    // link[8][7]=1;link[8][2]=1;link[8][9]=1;
    // link[9][8]=1;link[9][11]=1;link[9][12]=1;
    // link[10][5]=1;link[10][11]=1;link[10][12]=1;
    // link[11][10]=1;link[11][9]=1;link[11][13]=1;
    // link[12][10]=1;link[12][13]=1;link[12][9]=1;
    // link[13][3]=1;link[13][11]=1;link[13][12]=1;


	cout << n << " "<< edge_count << endl;
	cout << (edge_count+4*n+edge_count*log(n)) << endl;

	//

	vector<double> resource_demand(k); // this vector contains the resource demand by all the VNFs
	for(int i=0;i<k;i++){
		resource_demand[i]= rand()%(int)(capacity/k);
	}
	resource_demand[0]=0.25*(capacity/k);
	// we assume to be using assured network bandwidth to all the VMs resulting in negligible queuing delay from switch to server
	display(resource_demand);
	/*============Dividing the existing VNFs in 2 categories based on resource requirements=================*/
	double thres=findThreshold(resource_demand);
	vector<int> low_usage,high_usage;
	for(int i=0;i<k;i++){
		if(resource_demand[i]<thres){
			low_usage.push_back(i);
		}else{
			high_usage.push_back(i);
		}
	}
	//********************************************************************************************************

	// cin >> flag_global;
	vector<int> f[k];  // Contains information of nodes where this vnf could be found;
	vector<int> node[n];
	vector<double> available_resource(n,capacity);
	// cout << f[0].size() << " " << f[1].size() << " " << f[2].size() << " " << f[3].size() << " ";
	vector<vector<double>> per_node_vnf_resource (n,vector<double>(k,0));
	vector<vector<double>> per_node_vnf_resource_slice0 (n,vector<double>(k,0));
	vector<vector<double>> per_node_vnf_resource_slice1 (n,vector<double>(k,0));
	vector<vector<double>> link_ut(n,vector<double>(n,MAX));
	  
	default_random_engine generator;
  	poisson_distribution<int> distribution(10);


  	for(int i=0;i<n;i++){
  		int number = distribution(generator);
  		if(number>0){
  			node[i].push_back(0);
  			f[0].push_back(i);
  		}
  	}
	// int numOfBaseStation = n/10;
	// for(int i=0,c=0;i<n;i++){ //10-15 (use poisson for distributing the BS connection)
	// 	if(rand()%2 && c<numOfBaseStation){
	// 		node[i].push_back(0);
	// 		per_node_vnf_resource[i][0] += resource_demand[0];
	// 		available_resource[i]-= resource_demand[0];
	// 		f[0].push_back(i);
	// 		c++;
	// 	}
	// }
	/*
		As of now we are doing instatntiation of VNF randomly
		We can do the wise placement of the VNFs and finally include the VNF interference analysis

		low - high
		blank - high - low
		
		we should also ensure that all the VNFs are almost equally instantiated into the network
	*/


	for(int i=0;i<n;i++){
		for(int j=1;j<maxVNFs;j++){
			node[i].push_back(j);
			f[j].push_back(i);
		}
	}
	// greedily place all the low usage nodes 
	/*
82.042 98.5268 83.5151
0.65615
Maximum utilised node status: 99.2
Max utilised link status: 128.289
191 0 815 0
389 426
504 502
Total reward achieved is: 33.2136
Simulation executed in 0 seconds
	*/
	// for(int i=0;i<n;i++){
	// 	for(int j=0;j<low_usage.size();j++){
	// 		if(resource_demand[low_usage[j]]<=available_resource[i]){
	// 			int t3 = high_usage[j];
	// 			node[i].insert(t3);
	// 			per_node_vnf_resource[i][t3]+= resource_demand[t3];
	// 			available_resource[i]-=resource_demand[t3]; //t4 is the high usage VNF
	// 			f[t3].push_back(i);
	// 		}
	// 	}
	// 	int high_vnf = findHighNode(high_usage,node[i]);
	// 	if(high_vnf==-1) continue;
	// 	while(available_resource[i]>=resource_demand[high_vnf]){
	// 		per_node_vnf_resource[i][high_vnf]+= resource_demand[high_vnf];
	// 		available_resource[i]-=resource_demand[high_vnf];
	// 	}
	// }

	// for(int i=0;i<n;i++){
	// 	// display(node[i]);
	// 	cout <<node[i].size() << endl;
	// }
	// for(int i=0;i<n;i++){
	// 	cout << i<< " " <<available_resource[i] << endl; 
	// }

	for(int i=0;i<n;i++){
		int div = node[i].size();
		for(auto v:node[i]){
		
			per_node_vnf_resource_slice0[i][v]= part*capacity/div; 
			per_node_vnf_resource_slice1[i][v]= (1-part)*capacity/div; 
			
			per_node_vnf_resource[i][v] = per_node_vnf_resource_slice0[i][v]+per_node_vnf_resource_slice1[i][v];
		}
	}
		
		/*
			1. We here itself allocate the remaining resources to the nodes considering uniform frequency of occurrence of VNFs
			2. Keeping some resource reserved for later
				a. Prioritise latency constraint users 
				b. Go in the sequence
				
		*/
	
	// cin >> flag_global;
	vector<vector<double>> node_vnf_resources= per_node_vnf_resource;
	vector<vector<double>> node_vnf_resources_slice0= per_node_vnf_resource_slice0;
	vector<vector<double>> node_vnf_resources_slice1= per_node_vnf_resource_slice1;
	vector<int> adj[n];
	vector<vector<int>> cost(n,vector<int>(n,MAX));
    for(int i=0,k;i<n;i++){
        for(int j=i;j<n;j++){
            if(link[i][j]==1) {
            	adj[i].push_back(j);
            	adj[j].push_back(i);
                cost[i][j] = 1;
                cost[j][i] = 1;
                link_ut[i][j] = 0;
                link_ut[j][i] = 0;
            }
        }
    }
    saveLinkInfo(link_ut,"link_res");
    // printGraph(adj);
    // displayLinkUt(link_ut);
    // cout << "low: ";
    // for(int i=0;i<low_usage.size();i++){
    // 	cout << low_usage[i] << " ";
    // }
    // cout << endl;
    // cout << "high: ";
    // for(int i=0;i<high_usage.size();i++){
    // 	cout << high_usage[i] << " ";
    // }
    // cout << endl;
    displayNodeUt(per_node_vnf_resource);
    displayNodeUt(per_node_vnf_resource_slice0);
    displayNodeUt(per_node_vnf_resource_slice1);
    printGraph(adj);
    cin >> flag_global;
    vector<vector<int>> reachable(n,vector<int>(n,0));
    vector<vector<vector<int>>> calculatedRoute(n,vector<vector<int>>(n)); //Contains all the routes from source to destination node
    for(int i=0;i<n;i++){  // Loop to generate route from one node to another using Dijkstra Algorithm based on physical location and connectivity
    	for(int j=i+1;j<n;j++){
            vector<int> visited(n,0);           
            reachable[i][j]= isReachable(i,j,adj);
            reachable[j][i] = reachable[i][j];
            if(reachable[i][j]){
            	cout << " I am reachable from " << i << " to " << j << endl;
            	calculatedRoute[i][j] = findPathDijkstra(j,i,cost,visited);
            	displayRoute(calculatedRoute[i][j]);
            	calculatedRoute[j][i] = reverseArray(calculatedRoute[i][j]);
            	max_length = max(max_length,(int)calculatedRoute[i][j].size());
            }
    	}
    }   
    displayCost(reachable); 
	cout << "Max Path length:" << max_length << endl;
	// cin >> flag_global;
	
    //Service Type definition in terms of latency and data rate requirement
    int numOfSlices=2;
	vector<pair<int,int>> service_type(numOfSlices);
	service_type[0] = {1,10};	// latency constraint {3,70}
	service_type[1] = {100,100}; 	

	vector<double> per_user_ut_node(k,0);
	for(int i=0;i<k;i++){
		// per_user_ut_node[i]= 0.01; // find the user specific VNF util
		per_user_ut_node[i]= 1*(rand()%4+1);
	}
	// display(per_user_ut_node);
	// cin >> flag_global;

	int fail_node_map=0,success=0,fail_link_map=0,flag_completed=1,failed_latency=0,failed_s0=0,failed_s1=0;;
	int  count1=0,count2=0;
	double max_link_ut=0.99*link_bw;
	long double reward=0;
	int index=-1;
	auto start = std::chrono::high_resolution_clock::now();
	vector<int> recent_observation(100,-1);

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
		* fixing number of nodes plot the variation of reward with number of sfc request the algorithms -- bar graph
		* Compare the above results with VNF sharing and without VNF sharing *** 
		* VNF utilisation improvement with addition of VNF sharing  vs  varied ratio of eMBB and URLLC traffic (embb=1 and vary urllc) --2 line graphs
	10. Contribution of our work in terms of novelty
	11. Abstract 
	12. Title of the paper - On-demand,  resource aware, SFC mapping, Network Slice, VNF sharing 
	13. Name of the algorithm
	14. System model and problem formulation
	15. Table algorithm and problem formulation in over leaf
	16. Make a flow chart
	17. Total throughput with and without VNF sharing
	*/

	/*
	1. Variable capacity for 
	*/
	double tot_cost=0,throughput=0;
	// while(NoOfZeros(recent_observation)<=99){
	for(int key=0;key<stoi(argv[4]);key++){
		int slice = rand()%30;
		if(slice<10) slice =1;
		else slice =0;
		// int slice =1;
		vector<vector<double>> resource;
		double rate, latency;
		if(slice==0){ 
			count1++;
			latency=service_type[slice].first /*- (rand()%50)/100.0*/;
			rate = service_type[slice].second;
			resource = per_node_vnf_resource_slice0;
		}
		else if(slice==1) {
			count2++;
			latency=service_type[slice].first;
			rate = service_type[slice].second;
			// rate += rand()%(int)(rate/2);	
			resource = per_node_vnf_resource_slice1;
		}
		// double latency=service_type[slice].first - (rand()%50)/10.0;
		// double rate = service_type[slice].second + (rand()%10);
		double op_cost = 0;
		double per_user_ut_link=findLinkUt(rate);
		cout << per_user_ut_node[0] << " " <<per_user_ut_link << " key=" << key <<endl;
		index++;
		index = index%100;
		//recieve sfc request from the user
		// std::this_thread::sleep_for(std::chrono::milliseconds(500));
		flag_completed = 0;
		int chain_length = 3;
		if(chain_length==0) chain_length++;
		// SFC chain request creation: 
		vector<int> sfc_chain;
		// Have only 3 SFCs in the system
		cout << "Checkpoint-1 " << chain_length<< endl;
		unordered_set<int> temp_chain;
		while(temp_chain.size() != chain_length){
			int t_ = rand()%k;
			t_ = max(t_,1);
			temp_chain.insert(t_);
		}
		sfc_chain.push_back(0);
		// sfc_chain.push_back(1);
		// sfc_chain.push_back(3);
		// sfc_chain.push_back(5);
		// sfc_chain.push_back(7);
		for(auto v:temp_chain){
		    sfc_chain.push_back(v);
		}

		for(auto v:sfc_chain){
		    std::cout << v << " ";
		}
		cout << endl;
		//RES - res
		// cout << "well and good till here"; cin >> flag_global;
		//------------------SFC generation ends here------------------
		/*
		SFC to Node mapping starts here, options are:
			a. Choosing nodes with minimum utilised VNFs
			b. Choosing nodes closer to the 0th node having 
		*/
		vector<int> node_map(sfc_chain.size(),-1);
		// VNF 0  assignments begins here	
		//vnf 0  has to be assigned randomly as the BS connects to a particular AP in core network only
		// we can assume VNF0 doesn't gets filled up as rapidly as others with more and more people getting added to it
		// cout << "I am here" << endl;
		// cin >> flag_global;
		
		int temp = f[0].size();
		node_map[0] = f[0][rand()%temp];
		// node_map[0] = f[0][0];
		// cout << per_node_vnf_resource[node_map[0]][0] << endl;
		if(slice==0){ 
			if(resource[node_map[0]][0]<1){
				fail_node_map++;
				failed_s0++;
    			recent_observation[index]=0;
    			cout << "failed here" << endl;
				continue;
			}
		}
		else{ 
			if(resource[node_map[0]][0]<1){
				fail_node_map++;
				failed_s1++;
	    		recent_observation[index]=0;
	    		cout << "failed here" << endl;
				continue;
			}			
		}
		cout << "Checkpoint-3b" << endl;
		//fixing egress node into the system
		// int last_node = chain_length; 
		// int t4 = f[sfc_chain[last_node]].size();
		// node_map[last_node] = f[sfc_chain[last_node]][rand()%t4];
		// while(node_map[0]==node_map[last_node])
		// 	node_map[last_node] = f[sfc_chain[last_node]][rand()%t4];
		// cout << node_map[last_node] << " " << sfc_chain[last_node] << endl;
		// if(slice==0){ 
		// 	if(resource[node_map[last_node]][sfc_chain[last_node]]<per_user_ut_node[sfc_chain[last_node]]*rate){
		// 		fail_node_map++;
		// 		failed_s0++;
  //   			recent_observation[index]=0;
  //   			cout << "failed here" << endl;
		// 		continue;
		// 	}
		// }
		// else{ 
		// 	if(resource[node_map[last_node]][sfc_chain[last_node]]<per_user_ut_node[sfc_chain[last_node]]*rate){
		// 		fail_node_map++;
		// 		failed_s1++;
	 //    		recent_observation[index]=0;
	 //    		cout << "failed here" << endl;
		// 		continue;
		// 	}			
		// }
		
		
		// display(node_map);
		// cout << "cool";cin >> flag_global;

		//after fixing f0 we could have either path to f1 precalculated or we could calculate it in real time here. 
		// VNF1 onwards mapping using OPTION (A)
		int flag_loaded=0,flag_sos=0;
		double req_res=-1;
		cout << "Checkpoint-3c" << endl;
		for(int i=1;i<sfc_chain.size();i++){
		// for(int i=sfc_chain.size()-2;i>0;i--){
			// req_res = per_user_ut_node[sfc_chain[i]]*rate; 
			req_res = 1;
			node_map[i] = find_min_usage(req_res,sfc_chain[i],f[sfc_chain[i]],node_map,resource,node_vnf_resources);
			if(node_map[i]==-1){
				// cout << "Node mapping has failed for " << sfc_chain[i] <<  endl;
				// display_utilisation(sfc_chain[i],resource);
				// displayNodeUt(resource);
				flag_loaded=1;
				break;
			}

			// cout << per_node_vnf_resource[node_map[i]][sfc_chain[i]] << " ";
			
			// node_map[i] = find_nearby_node(req_res, sfc_chain[i],f[sfc_chain[i]],node_map[i-1],node_map,resource,calculatedRoute);
			// if(node_map[i]==MAX){
			// 	cout << "Node mapping has failed for " << sfc_chain[i] <<  endl;
			// 	display(f[sfc_chain[i]]);
			// 	display_utilisation(sfc_chain[i],resource);
			// 	display(node_map);
			// 	flag_loaded=1;
			// 	break;
			// }		
			// node_map[i] = find_good_node(req_res,sfc_chain[i],f[sfc_chain[i]],node_map[i+1],node_map,resource,per_node_vnf_resource,link_ut);
			// if(node_map[i]==-1 || resource[node_map[i]][sfc_chain[i]]<req_res){
			// 	cout << "Node mapping has failed for " << sfc_chain[i] <<  endl;
			// 	// display(f[sfc_chain[i]]);
			// 	// display_utilisation(sfc_chain[i],per_node_vnf_resource);
			// 	// display(node_map);
			// 	flag_loaded=1;
			// 	break;
			// }
			// node_map[i]= find_best_match(req_res,sfc_chain[i],f[sfc_chain[i]],node_map[i-1],node_map,resource,link_ut,calculatedRoute);
			// if(node_map[i]==-1|| resource[node_map[i]][sfc_chain[i]]<=req_res){
			// 	// cout << "Node mapping has failed for " << sfc_chain[i] <<  endl;
			// 	// display_utilisation(sfc_chain[i],per_node_vnf_resource);
			// 	// display_utilisation(per_node_vnf_resource);
			// 	flag_loaded=1;
			// 	break;
			// }

		}
		// cout << endl;
		cout << "Checkpoint-3d" << endl;
		if(flag_loaded){
			// fail_count++;
			if(slice==0) failed_s0++;
			if(slice==1) failed_s1++;
			// break;
			fail_node_map++;
    		recent_observation[index]=0;
			continue;
		}
		cout << "Checkpoint-4a" << endl;
		//======================End of Node Mapping============================


		//**********************Start of Node route========================
    	vector<int> route;
    	route.push_back(node_map[0]);
    	for(int i=0;i<node_map.size()-1;i++){
    		vector<int> visited(n,0);
    		vector<int> temp;
    		if(reachable[node_map[i]][node_map[i+1]]){
    			cout << " It is reachable" << endl;
    		}
    		else{
    			flag_loaded=1;
    			cout << "finding path failed " << node_map[i] << " " << node_map[i+1] << endl;
    		}
    		if(flag_loaded) break;
    		temp = findPathDijkstra(node_map[i],node_map[i+1],link_ut,visited);
    		temp = reverseArray(temp);
    		
    		// vector<int> temp = calculatedRoute[node_map[i]][node_map[i+1]];
    		// vector<int> temp = findPathDFS(node_map[i+1],node_map[i],cost,visited,adj);
    		// findBestRoute(temp,link_ut,adj);
    		// cout << i << " done" << endl;
    		for(auto v=temp.begin()+1;v!=temp.end();v++){
    			route.push_back(*v);
    		}
    		temp.clear();
    	}
    	if(flag_loaded){
    		if(slice==0) failed_s0++;
			if(slice==1) failed_s1++;
			// break;
			fail_link_map++;
    		recent_observation[index]=0;
			continue;
    	}
		cout << "Checkpoint-4b" << endl;

    	if(slice==0){
	    	if(findTotalLatency(route,sfc_chain,node_map,link_ut,resource,node_vnf_resources_slice0,node_vnf_resources)>1){
	    		cout << "failed in latency" << endl;
	    		failed_latency++;
	    		recent_observation[index]=0;
				if(slice==0) failed_s0++;
				if(slice==1) failed_s1++;
	    		continue;
	    	}
    	}else if(slice==1){
			for(int i=0;i<node_map.size();i++){
				if(per_node_vnf_resource_slice1[node_map[i]][sfc_chain[i]]<1){
					flag_loaded=1;
					fail_node_map++;
					break;
				}
			}
			if(flag_loaded){
    			recent_observation[index]=0;
				if(slice==0) failed_s0++;
				if(slice==1) failed_s1++;
				continue;
			}
    	}

		cout << "Checkpoint-5" << endl;
    	for(int i=1;i<route.size();i++){
    		if(link_ut[route[i-1]][route[i]]>=link_bw){
    			flag_loaded=1;
    			cout << "link overflow between "<< route[i-1] << " - "<< route[i]  << endl;
    			break;
    		}
    	}

    	//===================End of Node routing=============================

		cout << "Checkpoint-6" << endl;
    	//Check if anything failed in node placement as well as routing
    	if(flag_loaded){
    		// break;
			if(slice==0) failed_s0++;
			if(slice==1) failed_s1++;
    		fail_link_map++;
    		recent_observation[index]=0;
    		continue;
    	}
    	else{
    	//if everytihing went well we are updating the values of node and link as the user gets connected
			for(int j=0;j<node_map.size();j++){
				cout << node_map[j] << " " <<  sfc_chain[j]<< endl;
				// displayNodeUt(per_node_vnf_resource);
				// per_node_vnf_resource[node_map[j]][sfc_chain[j]]-= per_user_ut_node[sfc_chain[j]]*rate;
				per_node_vnf_resource[node_map[j]][sfc_chain[j]]-= 1;
				op_cost += per_user_ut_node[sfc_chain[j]]*rate;
				if(slice==0){
					// per_node_vnf_resource_slice0[node_map[j]][sfc_chain[j]]-= per_user_ut_node[sfc_chain[j]]*rate;
					per_node_vnf_resource_slice0[node_map[j]][sfc_chain[j]]-= 1;
					// cout<< slice << " "  << per_node_vnf_resource[node_map[j]][sfc_chain[j]] << " " << per_user_ut_node[sfc_chain[j]]/service_type[slice].first<< endl;
					// int t00; cin >> t00;
					// if(per_node_vnf_resource_slice0[node_map[j]][sfc_chain[j]]<0){
					// 	cout << slice << "Req Resource:" << per_user_ut_node[sfc_chain[j]] << " " << j << " " << rate << endl;
					// 	cout << node_map[j] << " " <<sfc_chain[j] << endl;
					// 	// displayNodeUt(per_node_vnf_resource);
					// 	displayNodeUt(per_node_vnf_resource_slice0);
					// 	// displayNodeUt(per_node_vnf_resource_slice1);
					// 	cin >> flag_global;
					// }
				}
				else if(slice==1){
					// per_node_vnf_resource_slice1[node_map[j]][sfc_chain[j]]-= per_user_ut_node[sfc_chain[j]]*rate;
					per_node_vnf_resource_slice1[node_map[j]][sfc_chain[j]]-= 1;
					// cout<< slice <<" " <<per_node_vnf_resource[node_map[j]][sfc_chain[j]] << " " << per_user_ut_node[sfc_chain[j]]*service_type[slice].second <<endl;
					// int t00; cin >> t00;
				}
			}
    		for(int i=1;i<route.size();i++){
    			link_ut[route[i-1]][route[i]]+=per_user_ut_link;
    			link_ut[route[i]][route[i-1]]+=per_user_ut_link;
    			op_cost += per_user_ut_link;
    			if(link_ut[route[i]][route[i-1]]>link_bw){
		    		flag_loaded=1;   					
    				for(int j=i;j>=1;j--){
		    			link_ut[route[i-1]][route[i]]-=per_user_ut_link;
		    			link_ut[route[i]][route[i-1]]-=per_user_ut_link; 
    				}
    				fail_link_map++;
    				break;
    			}
    		}
	    	if(flag_loaded){
	    		// break;
				if(slice==0) failed_s0++;
				if(slice==1) failed_s1++;
	    		fail_link_map++;
	    		recent_observation[index]=0;
	    		continue;
	    	}
    	}
    	//We increment the number of successful user
		cout << "Checkpoint-7" << endl;
    	success++;
    	recent_observation[index]=1;
    	tot_cost+=op_cost;
    	throughput += rate;
    	if(slice==0) reward+= service_type[slice].first/(latency*op_cost);
    	else if(slice==1) reward += rate/(op_cost*service_type[slice].second);
    	// reward = reward/;
		cout << success << " " << fail_link_map<<" "<<fail_node_map << " " << reward <<endl;
		flag_completed=1;
	}
	auto stop = std::chrono::high_resolution_clock::now();
    displayLinkUt(link_ut);
    saveLinkInfo(link_ut,"link"); saveNodeInfo(per_node_vnf_resource,node_vnf_resources,"node");
    saveNodeInfo(per_node_vnf_resource_slice0,node_vnf_resources_slice0,"slice0");
    saveNodeInfo(per_node_vnf_resource_slice1,node_vnf_resources_slice1,"slice1");
    cout << "Node resource" << endl;
    displayNodeUt(per_node_vnf_resource);
    // displayNodeUt(per_node_vnf_resource_slice0);
    cout << "Slice 1 resource" << endl;
    displayNodeUt(per_node_vnf_resource_slice1);
    double util =percentUtilNode(node_vnf_resources,per_node_vnf_resource); 
    cout  << util << " " << percentUtilNode(node_vnf_resources,per_node_vnf_resource_slice0) << " "<< percentUtilNode(node_vnf_resources,per_node_vnf_resource_slice1)<<endl;
    cout  <<percentUtilLink(link_ut) << endl;
    max_utilisation_node(per_node_vnf_resource,node_vnf_resources);
    max_utilisation(link_ut);
	cout << "success" << " " << "fail_link_map"<<" "<<"fail_node_map" <<" " << "failed_latency" <<endl;
	cout << success << " " << fail_link_map<<" "<<fail_node_map <<" " << failed_latency <<endl;
	cout <<failed_s0 << " " <<failed_s1 << endl;
	cout << count1 << " " << count2 << endl;
	// cout << success << endl;
	cout << "Total reward achieved is: " << reward << endl;
	double accRatio = success/(double)(count1+count2);
	cout << "Acceptance ratio: "  << accRatio  << endl;
	cout << "Total cost:" << tot_cost << endl;
	cout << "throughput:" << throughput << endl;
 	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);  
	cout << "Simulation executed in "<<duration.count() << " milliseconds" << endl;
	av_reward+= reward/stod(argv[2]);
	av_time += duration.count()/stod(argv[2]);
	av_util += util/stod(argv[2]);
	av_accRatio += accRatio/stod(argv[2]);
	av_cost += tot_cost/stod(argv[2]);
	av_throughput+= throughput/stod(argv[2]);
	cout << av_reward << endl;
	// cin >> flag_global;
	}

	cout << "\nPoint:" << stoi(argv[4]) << endl;
	cout << "av_time/1000.0 << " " << av_reward << " " << av_util << " " << av_accRatio << " " << av_cost << " " << av_throughput" << endl;
	cout << av_time/1000.0 << "\n" << av_reward << "\n" << av_util << "\n" << av_accRatio << "\n" << av_cost << "\n" << av_throughput << endl;

	ofstream fout;
	fout.open("data.txt");
	fout << "\n Point:" << stoi(argv[4]) << endl;
	fout << "Average Time: " << av_time << endl;
	fout << "Average Reward: " << av_reward << endl;
	fout << "Average Utilisation: " << av_util << endl;
	fout << "Acceptance ratio: " << av_accRatio << endl;
	fout << "Average cost: " << av_cost << endl;
	fout.close();
	// cout << '\a';
}