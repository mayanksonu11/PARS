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
#define R 100   // Connecitvity range definition
typedef long ll;
int cost1,dist1,n;    // Global varaibles utilised in the code
int flag_global=0,calls=0,max_length=0;
double wt,throughput,reward,total_latency;
double link_bw = 10000; // in MHz
// double cpu_freq =8000; //in MHz 
double capacity=40000;
int k=8,maxVNFs=8; //k is types of vnf
int numOfUsers,chain_length=k/2,iter_val;
int latency_slice0=1,latency_slice1=100;
int rate_slice0=10,rate_slice1=100;


void findRoute(int s, int f, vector<vector<int>> dist,vector<vector<int>> routeCost,vector<vector<int>> adj[],vector<vector<int>>cost,vector<int> &route,vector<vector<vector<int>>> calculatedRoute);

template <typename t>
void display(vector<t> &a){
    for(int i=0;i<a.size();i++){
        cout << a[i] << " ";
    }
    cout << "\n";
}

bool comp(const pair<double,vector<pair<int,vector<int>>>> &a, pair<double,vector<pair<int,vector<int>>>> &b)
{
	return a.first>b.first;
}


bool sortbysize(const vector<int> &a,
              const vector<int> &b)
{
    return (a.size() > b.size());
}
void saveLinkInfo(vector<vector<double>> link_ut, string str);
void saveNodeInfo(vector<vector<double>> node_ut,vector<vector<double>> avail_res, string str);
void displayLinkUt(vector<vector<double>> link_ut); // To display distance matrix
void displayNodeUt(vector<vector<double>> node_ut);
void printGraph(vector<int> adj[]); //To display adjacency list 
void displayCost(vector<vector<int>> cost); // To display cost matrix
void displayRoute(vector<int> route); // To display the route
double findCost(vector<int>route,vector<vector<double>> cost); // Given a route in the form of vector, used to find out the total cost 
double std_dev(vector<int>route, vector<vector<double>>&cost);
void findBestRoute(vector<int>& route,vector<vector<double>> cost, vector<int>adj[]);
int findPathDFS(int level, int path_len, int s, int f, vector<int> visited,vector<int>&parents ,vector<int> adj[]);

template <typename t>
vector<int> findPathDijkstra(int s,int f, vector<vector<t>> cost , vector <int> &visited){ // Verified-1
    visited[s]=1;
    vector<t> dist = cost[s];
    vector<int> parents(dist.size(), s);
    int flag=0;
    vector<int> route;
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
        
        for (int j = 0; j < dist.size(); j++)
        {
            if (visited[j] == 1)
                continue;
            if (dist[j] > dist[minIndex] + cost[minIndex][j])
                parents[j] = minIndex;
            flag = 1;
            dist[j] = min(dist[j], dist[minIndex] + cost[minIndex][j]);
        }
        visited[minIndex] = 1;
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

vector<int> reverseArray(vector<int> arr); // To reverse a given array


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

double max_utilisation_node(vector<vector<double>> &mat,vector<vector<double>> &resources );

double max_utilisation(vector<vector<double>> &mat);


bool notInRoute(vector<int> route, int k);

	// node_map[i] = find_min_usage(sfc_chain[i],f[sfc_chain[i]],node_map,per_node_vnf_resource);


int find_min_usage(double req_res,int vnf, vector<int> f, vector<int> node_map ,vector<vector<double>> avail_res,vector<vector<double>> res);
// node_map[i] = find_good_node(req_res,sfc_chain[i],f[sfc_chain[i]],node_map[i+1],node_map,resource,per_node_vnf_resource,link_ut);

int find_good_node(double req_res, int vnf, vector<int> f,int next_node, vector<int> node_map ,vector<vector<double>> &avail_res, vector<vector<double>> link_ut);

int find_nearby_node(double req_res,int vnf, vector<int> f, int prev_node, vector<int> node_map,vector<vector<double>> &avail_res, vector<vector<vector<int>>> &calculatedRoute);

int find_best_match(double req_res,int vnf, vector<int> f, int prev_node, vector<int> node_map,vector<vector<double>> &node_ut,vector<vector<double>> &link_ut ,vector<vector<vector<int>>> &calculatedRoute);

void display_utilisation(int vnf, vector<vector<double>> node_ut);

double percentUtilNode(vector<vector<double>>&resources,vector<vector<double>> &node_ut);

double percentUtilLink(vector<vector<double>> link_ut);

int NoOfZeros(vector<int> recent_observation);

void addANode(vector<int>&route, vector<int> adj[]);

double findLinkUt(int rate);

double findReqNodeResource(int latency,int rate);

double delayInNode(int node, double ut);

double findTotalLatency(vector<int> route, vector<int> sfc, vector<int> node_map,vector<vector<double>> &link_ut,vector<vector<double>> &avail_res,vector<vector<double>> &net_res,int slice);

double findThreshold(vector<double> t);


int findHighNode(vector<int> high_usage, set<int> vnf);

bool isReachable(int s, int d, vector<int> adj[]);

void display2d(vector<vector<int>> arr);
void displayPair(vector<pair<double,double>> per_user_req);

void modify(vector<vector<int>> & per_user_node_map){
	
}



vector<pair<int,vector<int>>> func1(vector<pair<double,double>> per_user_req,vector<vector<int>> per_user_sfc,vector<int> per_user_slice,vector<pair<int,int>> service_type,vector<double> service_charge,vector<vector<double>>per_node_vnf_resource,vector<vector<double>>per_node_vnf_resource_slice0,vector<vector<double>>per_node_vnf_resource_slice1,
	vector<double> per_user_ut_node,vector<vector<int>> reachable,vector<vector<vector<int>>>calculatedRoute,vector<vector<double>>link_ut,vector<int>f[],int type, string model);


// Verified - 1 	
int find_mixed_node(double req_res,int vnf, vector<int> f,int prev_node, vector<int> node_map ,vector<vector<double>> node_ut,vector<vector<vector<int>>> &calculatedRoute){
	int best_node = -1;
	double evaluation = -1;
	for(int i=0;i<f.size();i++){
		if(node_ut[f[i]][vnf]<req_res){
			continue;
		}if(calculatedRoute[prev_node][f[i]].size()==0){
			continue;
		}
		double temp_evaluation = (wt)*calculatedRoute[prev_node][f[i]].size()/(double)max_length+(1-wt)*node_ut[f[i]][vnf]/100.0;
		// cout << temp_evaluation << "\n";
		if(temp_evaluation>evaluation){
			evaluation = temp_evaluation;
			best_node = f[i];
		}
	}
	return best_node;
}

// Verified - 1
int findRandomNode(double req_res,int vnf, vector<int> f,int prev_node, vector<int> node_map ,vector<vector<double>> resource){
	int best_node = -1;
	int count=0;
	while(best_node==-1 and count<n){
		int temp =  f[rand()%f.size()];
		if(resource[temp][vnf]>req_res)
			best_node = temp;
		count++;
	}
	return best_node;
}


// Verified - 1
void displayPair(vector<pair<double,double>> per_user_req)
{
	for(int i=0;i<per_user_req.size();i++)
	{
		cout << per_user_req[i].first << " " << per_user_req[i].second << endl;
	}
	cout << per_user_req.size();
}

// Verified - 0
void saveLinkInfo(vector<vector<double>> link_ut, string str){
	ofstream fout;
	fout.open(str+".txt");
	fout << link_ut.size() << " ";
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

// Verified - 0
void saveNodeInfo(vector<vector<double>> node_ut,vector<vector<double>> avail_res, string str){
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
			if(avail_res[i][j]>0.5)fout << (avail_res[i][j]-node_ut[i][j])/avail_res[i][j] << " ";
			else fout << 0 << " ";	
		} 
		fout << endl;
	}
	fout.close();
}

//Verified - 0
void displayLinkUt(vector<vector<double>> link_ut){ //Verified-1
	for(int i=0;i<link_ut.size();i++){
		for(int j=0;j<link_ut[i].size();j++){
			if(link_ut[i][j]> link_bw){
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

// Verified - 0
void displayNodeUt(vector<vector<double>> node_ut){ //Verified-1
	for(int i=0;i<node_ut.size();i++){
		for(int j=0;j<node_ut[i].size();j++){
			if(node_ut[i][j]==0)
            	cout << " = ";
			else cout<< node_ut[i][j] << " ";	
		} 
		cout << endl;
	}
	cout << endl;
}

// Verified - 0
void printGraph(vector<int> adj[]) // Verified-1
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

// Verified - 0
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

// Verified - 0
void displayRoute(vector<int> route){ // To display the route
	int i;
	for(i=0;i<route.size()-1;i++){
		cout << route[i] << " -> ";
	}
	cout << route[i] << endl;
}

//Verified -1
double findCost(vector<int>route,vector<vector<double>> cost){ // Given a route in the form of vector, used to find out the total cost 
	double sum=0;
	for(int i=0;i<route.size()-1;i++){
		sum+= cost[route[i]][route[i+1]];
		// cout  << cost[route[i]][route[i+1]]<< "--";
	}
	return sum;
}

// Verified - 0 
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

// Crucial (NA)
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

// Crucial (NA)
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

// Verified - 1
vector<int> reverseArray(vector<int> arr){ // To reverse a given array
	int len = arr.size();
	vector<int> temp;
	for(int i=len-1;i>=0;i--){
		temp.push_back(arr[i]);
	}
	return temp;
}

// Verified - 1
double max_utilisation_node(vector<vector<double>> &avail_res,vector<vector<double>> &net_res ){
	double val=-1;
	for(int i=0;i<avail_res.size();i++){
		for(int j=0;j<avail_res[i].size();j++){
			if(avail_res[i][j]!=0){
				val=max(val,(net_res[i][j]-avail_res[i][j])/net_res[i][j]*100);
			}
		}
	}
	cout << "Maximum utilised node status: " <<val << endl; 
	return val;
}

// Verified - 1
double max_utilisation(vector<vector<double>> &avail_res){
	double val=-1;
	for(int i=0;i<avail_res.size();i++){
		for(int j=0;j<avail_res[i].size();j++){
			if(avail_res[i][j]<MAX){
				val=max(val,avail_res[i][j]);
			}
		}
	}
	cout << "Max utilised link status: " <<val << endl; 
	return val;
}

// Verified - 1 
bool notInRoute(vector<int> route, int k){
    for(int i=0;i<route.size();i++) 
    	if(route[i]==k) return false;
    return true;
}

// find_min_usage(req_res,sfc_chain[i],f[sfc_chain[i]],node_map,per_node_vnf_resource);

// Verified - 1
int find_min_usage(double req_res,int vnf, vector<int> f, vector<int> node_map ,vector<vector<double>> avail_res,vector<vector<double>> net_res){
	int min_ut_node=-1;
	double max_remaining = -1e9;
	for(int i=0;i<f.size();i++){
		if(avail_res[f[i]][vnf]*100/net_res[f[i]][vnf] >= max_remaining && avail_res[f[i]][vnf]>req_res && notInRoute(node_map,f[i]) )
		{
			// cout << i << " ";
			max_remaining = avail_res[f[i]][vnf];
			min_ut_node = f[i];
		}
	}
	return min_ut_node;
}

// Curcial 
int find_good_node(double req_res, int vnf, vector<int> f,int next_node, vector<int> node_map ,vector<vector<double>> &avail_res, vector<vector<double>> link_ut){ // 5GE
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
		vector<int> temp_route = findPathDijkstra(f[i],next_node,link_ut,visited);
		double t1 = link_bw*(temp_route.size()-1) - findCost(temp_route,link_ut);
		// cout << "Cool"; cin >> flag_global;
		double t2 = 0;
		for(int j=0;j<avail_res[f[i]].size();j++){
			t2 += avail_res[f[i]][j];
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

// Verified - 1 
int find_nearby_node(double req_res,int vnf, vector<int> f, int prev_node, vector<int> node_map,vector<vector<double>> &avail_res, vector<vector<vector<int>>> &calculatedRoute){
	int min_dist_node=MAX;
	// display(f);
	// cout << vnf << " " << prev_node << endl; 
	for(int i=0;i<f.size();i++){
		if(calculatedRoute[prev_node][f[i]].size()<min_dist_node && notInRoute(node_map,f[i]) && avail_res[f[i]][vnf]>=req_res){
			min_dist_node = f[i];
		}
	}
	return min_dist_node;
}

// Verified - 1
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

// Crucial 
int find_best_match(int slice,double req_res,int vnf, vector<int> f, int prev_node, vector<int> node_map,vector<vector<double>> &avail_res,vector<vector<double>> &link_ut ,vector<vector<vector<int>>> &calculatedRoute){
	vector<double> w(f.size(),0);
	vector<double> avg_bw(f.size(),0);
	double max_node_ut=-1;
	for(int i=0;i<f.size();i++){
		max_node_ut = max(max_node_ut,avail_res[f[i]][vnf]);
	}
	double max_node_delay = -1;
	for(int i=0;i<f.size();i++){
		max_node_delay = max(max_node_delay,delayInNode(f[i],avail_res[f[i]][vnf]));
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
		else if(slice==1)
			w[i]= (1/(double)calculatedRoute[prev_node][f[i]].size())*(avg_bw[i]/max_link_ut)*(avail_res[f[i]][vnf]/(double)max_node_ut);
		else if(slice==0)
			w[i]= (1- delayInNode(f[i],avail_res[f[i]][vnf])/max_node_delay)*(avail_res[f[i]][vnf]/(double)max_node_ut);
		// cout << w[i] << " ";
	}
	// display(w);
	// cout << " check-3" << endl;
	double max_eval=-1;
	int node=-1;
	for(int i=0;i<w.size();i++){
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
			if(node_ut[i][j]!=0){
				value += (resources[i][j]-node_ut[i][j])/resources[i][j];
				count++;
			}
		}
	}
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
	return 20*e;
}

	    // if(findTotalLatency(route,sfc_chain,node_map,link_ut,resource,node_res_slice0,node_res,slice)>latency){

double findTotalLatency(vector<int> route, vector<int> sfc, vector<int> node_map,vector<vector<double>> &link_ut,vector<vector<double>> &avail_res,vector<vector<double>> &net_res,int slice){
	/*
		Given the route find out how much is the propagation delay and how much is the processing delay
		it is sum of the propagation delay in 
	*/
	// double latency = /*Transmission time in air and propagation time+ processing time at access points+/+/*SFC Processing time*/
	double tot_lat=0;
	for(int i=0;i<node_map.size();i++){
		if(avail_res[node_map[i]][sfc[i]]==0) 
			return MAX;
		double ut = (net_res[node_map[i]][sfc[i]]-avail_res[node_map[i]][sfc[i]])*100/net_res[node_map[i]][sfc[i]];
		// cout << "Utilisation is :" << ut << endl;
		tot_lat+= delayInNode(node_map[i],ut);
		// cout << "tot_lat" <<tot_lat <<" ut=" << ut << endl ;
	}
	// cout <<"Latency is:" << tot_lat << endl;
	// latency in the core network links
	//delay in RAN
	//delay in UE transmission
	// 20*tot_lat; //Effective delay for embb Users

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


bool isReachable(int s, int d, vector<int> adj[]){ //Verified-1
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

void copyToClipboard(string b){
	string a = "echo -n ";
	// string b = "'Hey Sonu'";
	string c = " | xsel -b";
	string d = a+b+c;
	// cout << d << endl;
	char *e = &(d[0]);
	// printf("%s\n", e);
	// int exitStatus = 0;
	system(e);
}

void display2d(vector<vector<int>> arr){
	for(int i=0;i<arr.size();i++){
		for(int j=0;j<arr[i].size();j++)
			cout <<arr[i][j] << " ";
		cout << endl;
	}
}

double utility(double f1,double f2 )
{
	return (1/(1+exp(-f1/f2)));
}

vector<pair<int,vector<int>>> func1(vector<pair<double,double>> per_user_req,vector<vector<int>> per_user_sfc,vector<int> per_user_slice,vector<pair<int,int>> service_type,vector<double> service_charge,vector<vector<double>>per_node_vnf_resource,vector<vector<double>>per_node_vnf_resource_slice0,vector<vector<double>>per_node_vnf_resource_slice1,vector<double> per_user_ut_node,vector<vector<int>> reachable,vector<vector<vector<int>>>calculatedRoute,vector<vector<double>>link_ut,vector<int>f[],int type, string model){

	ofstream fout;
	fout.open("data.txt");
	throughput=0,reward=0;
	double tot_cost=0;
	int fail_node_map=0,success=0,fail_link_map=0,flag_completed=1,failed_latency=0,failed_s0=0,failed_s1=0,failed_latency_s0=0,failed_node_map_s0=0,failed_node_map_s1=0,failed_link_map_s0=0,failed_link_map_s1=0,failed_latency_s1=0;
	int index=0;
	int count1=0,count2=0;
	vector<int> recent_observation(100,-1);
	vector<vector<double>> node_res_slice0= per_node_vnf_resource_slice0;
	vector<vector<double>> node_res_slice1= per_node_vnf_resource_slice1;
	vector<vector<double>> node_res = per_node_vnf_resource;
	double av_reward=0,av_time=0,av_util=0,av_accRatio=0,av_cost=0,av_throughput=0;
	auto start = std::chrono::high_resolution_clock::now();
    cout  << percentUtilNode(node_res,per_node_vnf_resource) << " " << percentUtilNode(node_res_slice0,per_node_vnf_resource_slice0) << " "<< percentUtilNode(node_res_slice1,per_node_vnf_resource_slice1)<<endl;
	// cin >> flag_global;
	if(type==2){
		sort(per_user_sfc.begin(),per_user_sfc.end(),sortbysize);
	}
	vector<pair<int,vector<int>>> per_user_node_map;
	// while(NoOfZeros(recent_observation)<=99){
	for(int key=0;key<numOfUsers;key++){
		int slice = per_user_slice[key];
		// int slice =1;
		vector<vector<double>> avail_res,net_res;
		if(slice==0){ 
			count1++;
			avail_res = per_node_vnf_resource_slice0;
			net_res = node_res_slice0;
		}
		else if(slice==1) {
			count2++;	
			avail_res = per_node_vnf_resource_slice1;
			net_res = node_res_slice1;
		}
		double rate, latency;
		latency = per_user_req[key].first;
		rate = per_user_req[key].second;
		// double latency=service_type[slice].first - (rand()%50)/10.0;
		// double rate = service_type[slice].second + (rand()%10);
		double op_cost = 0;
		// double per_user_ut_link=findLinkUt(rate);
		double per_user_ut_link=rate;
		// cout << per_user_ut_node[0] << " " <<per_user_ut_link << " key=" << key <<endl;
		index++;
		index = index%100;
		//recieve sfc request from the user
		// std::this_thread::sleep_for(std::chrono::milliseconds(500));
		flag_completed = 0;
		// SFC chain request creation: 
		vector<int> sfc_chain = per_user_sfc[key];
		// Have only 3 SFCs in the system
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
		double req_res=-1;
		req_res = per_user_ut_node[sfc_chain[0]]*rate; 
				
		int temp = f[0].size();
		node_map[0] = f[0][rand()%temp];
		// node_map[0] = f[0][0];
		// cout << per_node_vnf_resource[node_map[0]][0] << endl;
		if(slice==0){ 
			if(avail_res[node_map[0]][0]<req_res){
				fail_node_map++;
				failed_node_map_s0++;
    			recent_observation[index]=0;
    			cout << "failed here1" << endl;
				continue;
			}
		}
		else{ 
			if(avail_res[node_map[0]][0]<req_res){
				fail_node_map++;
				failed_node_map_s1++;
	    		recent_observation[index]=0;
	    		cout << "failed here2" << endl;
				continue;
			}			
		}
		// cout << "Checkpoint-3b" << endl;
		int chain_length = sfc_chain.size()-1;
		//fixing egress node into the system
		if(type==2){
			int last_node = chain_length; 
			int t4 = f[sfc_chain[last_node]].size();
			// display(sfc_chain);
			// cout << sfc_chain[last_node] << "All good " << t4 << endl;
			// cin >> flag_global;
			node_map[last_node] = f[sfc_chain[last_node]][rand()%t4];
			while(node_map[0]==node_map[last_node])
				node_map[last_node] = f[sfc_chain[last_node]][rand()%t4];
			// cout << node_map[last_node] << " " << sfc_chain[last_node] << endl;
			if(slice==0){ 
				if(avail_res[node_map[last_node]][sfc_chain[last_node]]<per_user_ut_node[sfc_chain[last_node]]*rate){
					fail_node_map++;
					failed_node_map_s0++;
    				recent_observation[index]=0;
    				cout << "failed here3" << endl;
					continue;
				}
			}
			else{ 
				if(avail_res[node_map[last_node]][sfc_chain[last_node]]<per_user_ut_node[sfc_chain[last_node]]*rate){
					fail_node_map++;
					failed_node_map_s1++;
	    			recent_observation[index]=0;
	    			cout << "failed here4" << endl;
					continue;
				}			
			}
		}
		
		// display(node_map);
		// cout << "cool";cin >> flag_global;

		//after fixing f0 we could have either path to f1 precalculated or we could calculate it in real time here. 
		// VNF1 onwards mapping using OPTION (A)
		int flag_loaded=0,flag_sos=0;
		// cout << "Checkpoint-3c" << endl;
		for(int i=0;i<sfc_chain.size();i++){
		// for(int i=sfc_chain.size()-2;i>0;i--){
			if(node_map[i]!=-1)
				continue;
			req_res = per_user_ut_node[sfc_chain[i]]*rate; 
			// req_res = 1;
			// node_map[i] = findRandomNode(req_res,sfc_chain[i],f[sfc_chain[i]],node_map[i-1],node_map,resource);
			// node_map[i] = find_mixed_node(req_res,sfc_chain[i], f[sfc_chain[i]], node_map[i-1], node_map ,avail_res,calculatedRoute);
			// if(node_map[i]==-1){
			// 	// cout << "Node mapping has failed for " << sfc_chain[i] <<  endl;
			// 	// display_utilisation(sfc_chain[i],resource);
			// 	// displayNodeUt(resource);
			// 	flag_loaded=1;
			// 	break;
			// }
			// if(node_map[i]==-1){
			// 	// cout << "Node mapping has failed for " << sfc_chain[i] <<  endl;
			// 	// display_utilisation(sfc_chain[i],resource);
			// 	// displayNodeUt(resource);
			// 	flag_loaded=1;
			// 	break;
			// }
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
			if(type==2){		
				node_map[i] = find_good_node(req_res,sfc_chain[i],f[sfc_chain[i]],node_map[i+1],node_map,avail_res,link_ut);
				if(node_map[i]==-1 || avail_res[node_map[i]][sfc_chain[i]]<req_res){
					// cout << "Node mapping has failed for " << sfc_chain[i] <<  endl;
					// display(f[sfc_chain[i]]);
					// display_utilisation(sfc_chain[i],per_node_vnf_resource);
					// display(node_map);
					flag_loaded=1;
					break;
				}
			}
			/*PARS*/ 
			if(model=="PARS"){
				node_map[i] = find_min_usage(req_res,sfc_chain[i], f[sfc_chain[i]], node_map ,avail_res,net_res);
				if(node_map[i]==-1|| avail_res[node_map[i]][sfc_chain[i]]<=req_res){
					// cout << "Node mapping has failed for " << sfc_chain[i] <<  endl;
					// display_utilisation(sfc_chain[i],per_node_vnf_resource);
					// display_utilisation(per_node_vnf_resource);
					flag_loaded=1;
					break;
				}
			}

			/*DM*/ 
			if(model=="DM"){
				node_map[i]= find_best_match(slice,req_res,sfc_chain[i],f[sfc_chain[i]],node_map[i-1],node_map,avail_res,link_ut,calculatedRoute);
				if(node_map[i]==-1|| avail_res[node_map[i]][sfc_chain[i]]<=req_res){
					// cout << "Node mapping has failed for " << sfc_chain[i] <<  endl;
					// display_utilisation(sfc_chain[i],per_node_vnf_resource);
					// display_utilisation(per_node_vnf_resource);
					flag_loaded=1;
					break;
				}
			}

		}
		// cout << endl;
		// cout << "Checkpoint-3d" << endl;
		if(flag_loaded){
			// fail_count++;
			if(slice==0) failed_node_map_s0++;
			if(slice==1) failed_node_map_s1++;
			// break;
			// cout << "Node Mapping Failed!" << endl;
			fail_node_map++;
    		recent_observation[index]=0;
			continue;
		}
		// cout << "Checkpoint-4a" << endl;
		//======================End of Node Mapping============================

		// display(node_map);
		//**********************Start of Node route========================
    	vector<int> route;
    	route.push_back(node_map[0]);
    	for(int i=0;i<node_map.size()-1;i++){
    		vector<int> visited(n,0);
    		vector<int> temp;
    		if(reachable[node_map[i]][node_map[i+1]]){
    			// cout << " It is reachable" << endl;
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
    		if(slice==0) failed_link_map_s0++;
			if(slice==1) failed_link_map_s1++;
			// break;
			fail_link_map++;
    		recent_observation[index]=0;
			continue;
    	}
		// cout << "Checkpoint-4b" << endl;

    	// if(slice==0){
	    if(findTotalLatency(route,sfc_chain,node_map,link_ut,avail_res,net_res,slice)>latency){
	    	// cout << "failed in latency" << endl;
	    	failed_latency++;
	    	recent_observation[index]=0;
			if(slice==0) failed_latency_s0++;
			if(slice==1) failed_latency_s1++;
	    	continue;
	    }
    	// }
    	// else if(slice==1){
		for(int i=0;i<node_map.size();i++){
			if(avail_res[node_map[i]][sfc_chain[i]]<req_res){
				flag_loaded=1;
				fail_node_map++;
				break;
			}
		}
		if(flag_loaded){
    		recent_observation[index]=0;
			if(slice==0) failed_node_map_s0++;
			if(slice==1) failed_node_map_s1++;
			continue;
		}

		// cout << "Checkpoint-5" << endl;
    	for(int i=1;i<route.size();i++){
    		if(link_ut[route[i-1]][route[i]]>=link_bw){
    			flag_loaded=1;
    			cout << "link overflow between "<< route[i-1] << " - "<< route[i]  << endl;
    			break;
    		}
    	}

    	//===================End of Node routing=============================

		// cout << "Checkpoint-6" << endl;
    	//Check if anything failed in node placement as well as routing
    	if(flag_loaded){
    		// break;
			if(slice==0) failed_link_map_s0++;
			if(slice==1) failed_link_map_s1++;
    		fail_link_map++;
    		recent_observation[index]=0;
    		continue;
    	}
    	else{
    	//if everytihing went well we are updating the values of node and link as the user gets connected
			for(int j=0;j<node_map.size();j++){
				// cout << node_map[j] << " " <<  sfc_chain[j]<< endl;
				// displayNodeUt(per_node_vnf_resource);
				per_node_vnf_resource[node_map[j]][sfc_chain[j]]-= per_user_ut_node[sfc_chain[j]]*rate;
				// per_node_vnf_resource[node_map[j]][sfc_chain[j]]-= 1;
				op_cost += per_user_ut_node[sfc_chain[j]]*rate;
				if(slice==0){
					per_node_vnf_resource_slice0[node_map[j]][sfc_chain[j]]-= per_user_ut_node[sfc_chain[j]]*rate;
					// per_node_vnf_resource_slice0[node_map[j]][sfc_chain[j]]-= 1;
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
					per_node_vnf_resource_slice1[node_map[j]][sfc_chain[j]]-= per_user_ut_node[sfc_chain[j]]*rate;
					// per_node_vnf_resource_slice1[node_map[j]][sfc_chain[j]]-= 1;
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
    				// fail_link_map++;
    				break;
    			}
    		}
	    	if(flag_loaded){
	    		// break;
				if(slice==0) failed_link_map_s0++;
				if(slice==1) failed_link_map_s1++;
	    		fail_link_map++;
	    		recent_observation[index]=0;
	    		continue;
	    	}
    	}

    	//We increment the number of successful user
		// cout << "Checkpoint-7" << endl;
    	success++;
    	recent_observation[index]=1;
    	tot_cost+=op_cost;
    	throughput += rate;
    	// if(slice==0) reward+= utility(service_type[slice].first,latency)/(op_cost); //Change with current representation
    	// else if(slice==1) reward += utility(rate,service_type[slice].second)/(op_cost); 
    	reward += service_charge[slice] - op_cost;
    	fout << key << " " << percentUtilLink(link_ut) << " " << percentUtilNode(node_res_slice0,per_node_vnf_resource_slice0) << " " << percentUtilNode(node_res_slice1,per_node_vnf_resource_slice1) << " "<< percentUtilNode(node_res,per_node_vnf_resource) <<endl;

    	// reward += utility()/op_cost;
    	total_latency += findTotalLatency(route,sfc_chain,node_map,link_ut,avail_res,net_res,slice);
    	cout << total_latency << "\n"; /*cin >> flag_global;*/
		cout << success << " " << fail_link_map<<" "<<fail_node_map << " " << reward <<endl;
		flag_completed=1;
		per_user_node_map.push_back({key,node_map});
		// if(success==iter_val){
		// 	break;
		// }
	}
	// fitnessScores.push_back(throughput);
	failed_s0 = failed_node_map_s0+failed_link_map_s0+failed_latency_s0;
	failed_s1 = failed_node_map_s1+failed_link_map_s1+failed_latency_s1;
	auto stop = std::chrono::high_resolution_clock::now();
   	displayLinkUt(link_ut);
    saveLinkInfo(link_ut,"link"); saveNodeInfo(per_node_vnf_resource,node_res,"node");
    saveNodeInfo(per_node_vnf_resource_slice0,node_res_slice0,"slice0");
    saveNodeInfo(per_node_vnf_resource_slice1,node_res_slice1,"slice1");
    cout << "Node resource" << endl;
    displayNodeUt(per_node_vnf_resource);
    // displayNodeUt(per_node_vnf_resource_slice0);
    cout << "Slice 0 resource" << endl;
    displayNodeUt(per_node_vnf_resource_slice0);    
    cout << "Slice 1 resource" << endl;
    displayNodeUt(per_node_vnf_resource_slice1);
    double util =percentUtilNode(node_res,per_node_vnf_resource); 
    // cout  << util << " " << percentUtilNode(node_res,per_node_vnf_resource_slice0) << " "<< percentUtilNode(node_res,per_node_vnf_resource_slice1)<<endl;
    cout  << percentUtilNode(node_res,per_node_vnf_resource) << " " << percentUtilNode(node_res_slice0,per_node_vnf_resource_slice0) << " "<< percentUtilNode(node_res_slice1,per_node_vnf_resource_slice1)<<endl;
    cout  <<percentUtilLink(link_ut) << endl;
    // copyToClipboard(to_string(util)+","+to_string(percentUtilLink(link_ut)));
    max_utilisation_node(per_node_vnf_resource,node_res);
    max_utilisation(link_ut);
	cout << "success" << " " << "fail_link_map"<<" "<<"fail_node_map" <<" " << "failed_latency" <<endl;
	cout << success << " " << fail_link_map <<" "<<fail_node_map <<" " << failed_latency <<endl;	
	cout << "fail_link_map_s0"<<" "<<"fail_node_map_s0" <<" " << "failed_latency_s0" <<endl;
	cout << failed_link_map_s0<<" "<<failed_node_map_s0 <<" " << failed_latency_s0 <<endl;
	cout << "fail_link_map_s1"<<" "<<"fail_node_map_s1" <<" " << "failed_latency_s1" <<endl;
	cout << failed_link_map_s1<<" "<<failed_node_map_s1 <<" " << failed_latency_s1 <<endl;
	cout <<"Number of users failed: Slice0= " <<failed_s0 << " and Slice1= " <<failed_s1 << endl;
	cout << "Total number of users in slice0= "<< count1 << " and slice1= " << count2 << endl;
	// cout << success << endl;
	cout << "Total reward achieved is: " << reward << endl;
	double accRatio = success/(double)(count1+count2);
	cout << "Acceptance ratio (overall): "  << accRatio << " slice0:"<< 1-failed_s0/(double)count1 << " slice1:" << 1-failed_s1/(double)count2  << endl;
	cout << "Total cost:" << tot_cost << endl;
	cout << "throughput:" << throughput << endl;
	cout << "Average Latency:" << total_latency/(double)(count1+count2) << endl;
 	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);  
	cout << "Simulation executed in "<<duration.count() << " milliseconds" << endl;
	cout << accRatio << " " << 1-failed_s0/(double)count1 << " " << 1-failed_s1/(double)count2 << " " << reward << " " << throughput << endl;
	// copyToClipboard(to_string(numOfUsers)+","+to_string(success)+","+to_string(1-failed_s0/(double)count1) /*+" "+to_string(1-failed_s1/(double)count2)*/ + ","+to_string(throughput));
	av_reward+= reward/numOfUsers;
	av_time += duration.count()/numOfUsers;
	av_util += util/numOfUsers;
	av_accRatio += accRatio/numOfUsers;
	av_cost += tot_cost/numOfUsers;
	av_throughput+= throughput/numOfUsers;
	cout << av_reward << endl;
	// cin >> flag_global;
	cout << "\nNum of Users:" <<numOfUsers << endl;
	fout.close();
	// cout << " << " " << av_reward << " " << av_util << " " << av_accRatio << " " << av_cost << " " << av_throughput" << endl;
	// cout <<"av_time/1000.0: "<<av_time/1000.0 << "\n"<< "av_reward: " <<  av_reward << "\n"<< "av_util: " << av_util << "\n"<< "av_accRatio: "<< av_accRatio << "\n" << "av_cost: " << av_cost << "\n" << "av_throughput: " <<av_throughput << endl;
	cout << "Throughput (Mbps):" << throughput << "\n";
	cout << "Total Acceptance:" << success << "\n";
	cout << "Acceptance(T0):" << count1 - failed_s0 << "\n";
	cout << "Acceptance(T1):" << count2 - failed_s1 << "\n";
	cout << "Acceptance Ratio:" << accRatio << "\n";
	cout << "Acceptance Rate(T0):" << 1-failed_s0/(double)count1 << "\n";
	cout << "Acceptance Rate(T1):" << 1-failed_s1/(double)count2 << "\n";
	cout << "Execution time:" << duration.count() << "\n";

	return per_user_node_map;
}

// void generateNeigbours(vector<pair<double,vector<pair<int,vector<int>>>>> &neighbors,vector<pair<int,vector<int>>> currSolution,vector<pair<double,double>> per_user_req,vector<vector<int>> per_user_sfc,vector<int> per_user_slice,vector<pair<int,int>> service_type,vector<vector<double>>per_node_vnf_resource,vector<vector<double>>per_node_vnf_resource_slice0,vector<vector<double>>per_node_vnf_resource_slice1,vector<double> per_user_ut_node,vector<vector<int>> reachable,vector<vector<vector<int>>>calculatedRoute,vector<vector<double>>link_ut,vector<int>f[])
// {
// 	// vector<pair<double,vector<pair<int,vector<int>>>>> neigbours;
// 	for(int i=0;i<currSolution.size();i++){
// 		int user_id = currSolution[i].first;
// 		vector<int> node_map = currSolution[i].second;
// 		vector<int> chain = per_user_sfc[user_id];
// 		int index = rand()%(node_map.size());
// 		vector<int> adjacencies = f[chain[index]]; //
// 		int total_adjacencies= adjacencies.size();
// 		node_map[index] = adjacencies[rand()%total_adjacencies];
// 		currSolution[i].second = node_map;
// 		double fitness = fitnessEvaluation(currSolution,per_user_req,per_user_sfc,per_user_slice,service_type,per_node_vnf_resource,per_node_vnf_resource_slice0,per_node_vnf_resource_slice1,per_user_ut_node,reachable,calculatedRoute,link_ut,f);
// 		neighbors.push_back({fitness,currSolution});
// 	}
// }
























// double fitnessEvaluation(vector<pair<int,vector<int>>> per_user_node_map,vector<pair<double,double>> per_user_req,vector<vector<int>> per_user_sfc,vector<int> per_user_slice,vector<pair<int,int>> service_type,vector<vector<double>>per_node_vnf_resource,vector<vector<double>>per_node_vnf_resource_slice0,vector<vector<double>>per_node_vnf_resource_slice1,vector<double> per_user_ut_node,vector<vector<int>> reachable,vector<vector<vector<int>>>calculatedRoute,vector<vector<double>>link_ut,vector<int>f[])
// {
// 	throughput = 0;
// 	double tot_cost=0,reward=0;
// 	int fail_node_map=0,success=0,fail_link_map=0,flag_completed=1,failed_latency=0,failed_s0=0,failed_s1=0;;
// 	int index=0;
// 	int count1=0,count2=0;
// 	vector<int> recent_observation(100,-1);
// 	vector<vector<double>> node_res_slice0= per_node_vnf_resource_slice0;
// 	vector<vector<double>> node_res_slice1= per_node_vnf_resource_slice1;
// 	vector<vector<double>> node_res= per_node_vnf_resource;
// 	double av_reward=0,av_time=0,av_util=0,av_accRatio=0,av_cost=0,av_throughput=0;

// 	int users = per_user_node_map.size();
// 	for(int i=0;i<users;i++){
// 		int key = per_user_node_map[i].first;
// 		int slice = per_user_slice[key];
// 		// int slice =1;
// 		vector<vector<double>> avail_res,net_res;
// 		if(slice==0){ 
// 			count1++;
// 			avail_res = per_node_vnf_resource_slice0;
// 			net_res = node_res_slice0;
// 		}
// 		else if(slice==1) {
// 			count2++;	
// 			avail_res = per_node_vnf_resource_slice1;
// 			net_res = node_res_slice1;
// 		}
// 		double rate, latency;
// 		latency = per_user_req[key].first;
// 		rate = per_user_req[key].second;
// 		double op_cost = 0;
// 		// double per_user_ut_link=findLinkUt(rate);
// 		double per_user_ut_link=rate;
// 		index++;
// 		index = index%100;
// 		//recieve sfc request from the user
// 		// std::this_thread::sleep_for(std::chrono::milliseconds(500));
// 		flag_completed = 0;
// 		// SFC chain request creation: 
// 		vector<int> sfc_chain = per_user_sfc[key];
// 		// Have only 3 SFCs in the system
// 		//RES - res
// 		// cout << "well and good till here"; cin >> flag_global;
// 		//------------------SFC generation ends here------------------
// 		/*
// 		SFC to Node mapping starts here, options are:
// 			a. Choosing nodes with minimum utilised VNFs
// 			b. Choosing nodes closer to the 0th node having 
// 		*/
// 		vector<int> node_map;
// 		// VNF 0  assignments begins here	
// 		//vnf 0  has to be assigned randomly as the BS connects to a particular AP in core network only
// 		// we can assume VNF0 doesn't gets filled up as rapidly as others with more and more people getting added to it
// 		// cout << "I am here" << endl;
// 		// cin >> flag_global;
// 		node_map = per_user_node_map[i].second;
// 		req_res = per_user_ut_node[sfc_chain[0]]*rate; 
// 		// node_map[0] = f[0][0];
// 		// cout << per_node_vnf_resource[node_map[0]][0] << endl;
// 		if(slice==0){ 
// 			if(avail_res[node_map[0]][0]<req_res){
// 				fail_node_map++;
// 				failed_s0++;
//     			recent_observation[index]=0;
//     			// cout << "failed here0 " << key << endl;
// 				continue;
// 			}
// 		}
// 		else{ 
// 			if(avail_res[node_map[0]][0]<req_res){
// 				fail_node_map++;
// 				failed_s1++;
// 	    		recent_observation[index]=0;
// 	    		// cout << "failed here" << endl;
// 				continue;
// 			}			
// 		}
// 		// cout << "Checkpoint-3b" << endl;
// 		int chain_length = sfc_chain.size()-1;
// 		//fixing egress node into the system
		
// 		// display(node_map);
// 		// // cout << "cool";cin >> flag_global;

// 		//after fixing f0 we could have either path to f1 precalculated or we could calculate it in real time here. 
// 		// VNF1 onwards mapping using OPTION (A)
// 		int flag_loaded=0,flag_sos=0;
// 		// cout << "Checkpoint-3c" << endl;
// 		// for(int i=0;i<sfc_chain.size();i++){
// 		// // for(int i=sfc_chain.size()-2;i>0;i--){
// 		// 	if(node_map[i]!=-1)
// 		// 		continue;
// 		// 	req_res = per_user_ut_node[sfc_chain[i]]*rate; 
// 		// 	// req_res = 1;
// 		// 	if(type==1){
// 		// 		node_map[i] = find_mixed_node(sfc_chain[i], f[sfc_chain[i]], node_map[i-1], node_map ,resource,calculatedRoute);
// 		// 		if(node_map[i]==-1){
// 		// 			// // cout << "Node mapping has failed for " << sfc_chain[i] <<  endl;
// 		// 			// display_utilisation(sfc_chain[i],resource);
// 		// 			// displayNodeUt(resource);
// 		// 			flag_loaded=1;
// 		// 			break;
// 		// 		}
// 		// 	}
// 		// 	// // cout << per_node_vnf_resource[node_map[i]][sfc_chain[i]] << " ";
			
// 		// 	// node_map[i] = find_nearby_node(req_res, sfc_chain[i],f[sfc_chain[i]],node_map[i-1],node_map,resource,calculatedRoute);
// 		// 	// if(node_map[i]==MAX){
// 		// 	// 	// cout << "Node mapping has failed for " << sfc_chain[i] <<  endl;
// 		// 	// 	display(f[sfc_chain[i]]);
// 		// 	// 	display_utilisation(sfc_chain[i],resource);
// 		// 	// 	display(node_map);
// 		// 	// 	flag_loaded=1;
// 		// 	// 	break;
// 		// 	// }		
// 		// 	// node_map[i] = find_good_node(req_res,sfc_chain[i],f[sfc_chain[i]],node_map[i+1],node_map,resource,per_node_vnf_resource,link_ut);
// 		// 	// if(node_map[i]==-1 || resource[node_map[i]][sfc_chain[i]]<req_res){
// 		// 	// 	// cout << "Node mapping has failed for " << sfc_chain[i] <<  endl;
// 		// 	// 	// display(f[sfc_chain[i]]);
// 		// 	// 	// display_utilisation(sfc_chain[i],per_node_vnf_resource);
// 		// 	// 	// display(node_map);
// 		// 	// 	flag_loaded=1;
// 		// 	// 	break;
// 		// 	// }
// 		// 	if(type==2){
// 		// 		node_map[i]= find_best_match(req_res,sfc_chain[i],f[sfc_chain[i]],node_map[i-1],node_map,resource,link_ut,calculatedRoute);
// 		// 		if(node_map[i]==-1|| resource[node_map[i]][sfc_chain[i]]<=req_res){
// 		// 			// // cout << "Node mapping has failed for " << sfc_chain[i] <<  endl;
// 		// 			// display_utilisation(sfc_chain[i],per_node_vnf_resource);
// 		// 			// display_utilisation(per_node_vnf_resource);
// 		// 			flag_loaded=1;
// 		// 			break;
// 		// 		}
// 		// 	}

// 		// }
// 		// // cout << endl;
// 		// // cout << "Checkpoint-3d" << endl;
// 		// if(flag_loaded){
// 		// 	// fail_count++;
// 		// 	if(slice==0) failed_s0++;
// 		// 	if(slice==1) failed_s1++;
// 		// 	// break;
// 		// 	// cout << "Something wrong!" << endl;
// 		// 	fail_node_map++;
//   //   		recent_observation[index]=0;
// 		// 	continue;
// 		// }
// 		// // cout << "Checkpoint-4a" << endl;
// 		//======================End of Node Mapping============================


// 		//**********************Start of Node route========================
//     	vector<int> route;
//     	route.push_back(node_map[0]);
//     	for(int i=0;i<node_map.size()-1;i++){
//     		vector<int> visited(n,0);
//     		vector<int> temp;
//     		if(reachable[node_map[i]][node_map[i+1]]){
//     			// cout << " It is reachable" << endl;
//     		}
//     		else{
//     			flag_loaded=1;
//     			// cout << "finding path failed " << node_map[i] << " " << node_map[i+1] << endl;
//     		}
//     		if(flag_loaded) break;
//     		temp = findPathDijkstra(node_map[i],node_map[i+1],link_ut,visited);
//     		temp = reverseArray(temp);
    		
//     		// vector<int> temp = calculatedRoute[node_map[i]][node_map[i+1]];
//     		// vector<int> temp = findPathDFS(node_map[i+1],node_map[i],cost,visited,adj);
//     		// findBestRoute(temp,link_ut,adj);
//     		// // cout << i << " done" << endl;
//     		for(auto v=temp.begin()+1;v!=temp.end();v++){
//     			route.push_back(*v);
//     		}
//     		temp.clear();
//     	}
//     	if(flag_loaded){
//     		if(slice==0) failed_s0++;
// 			if(slice==1) failed_s1++;
// 			// break;
// 			fail_link_map++;
//     		recent_observation[index]=0;
// 			continue;
//     	}
// 		// cout << "Checkpoint-4b" << endl;

//     	if(slice==0){
// 	    	if(findTotalLatency(route,sfc_chain,node_map,link_ut,resource,node_res_slice0,slice)>latency){
// 	    		// cout << "failed in latency" << endl;
// 	    		failed_latency++;
// 	    		recent_observation[index]=0;
// 				if(slice==0) failed_s0++;
// 				if(slice==1) failed_s1++;
// 	    		continue;
// 	    	}
//     	}else if(slice==1){
// 			for(int i=0;i<node_map.size();i++){
// 				if(per_node_vnf_resource_slice1[node_map[i]][sfc_chain[i]]<req_res){
// 					flag_loaded=1;
// 					fail_node_map++;
// 					break;
// 				}
// 			}
// 			if(flag_loaded){
//     			recent_observation[index]=0;
// 				if(slice==0) failed_s0++;
// 				if(slice==1) failed_s1++;
// 				continue;
// 			}
//     	}

// 		// cout << "Checkpoint-5" << endl;
//     	for(int i=1;i<route.size();i++){
//     		if(link_ut[route[i-1]][route[i]]>=link_bw){
//     			flag_loaded=1;
//     			// cout << "link overflow between "<< route[i-1] << " - "<< route[i]  << endl;
//     			break;
//     		}
//     	}

//     	//===================End of Node routing=============================

// 		// cout << "Checkpoint-6" << endl;
//     	//Check if anything failed in node placement as well as routing
//     	if(flag_loaded){
//     		// break;
// 			if(slice==0) failed_s0++;
// 			if(slice==1) failed_s1++;
//     		fail_link_map++;
//     		recent_observation[index]=0;
//     		continue;
//     	}
//     	else{
//     	//if everytihing went well we are updating the values of node and link as the user gets connected
// 			for(int j=0;j<node_map.size();j++){
// 				// cout  << node_map[j] << " " <<  sfc_chain[j]<< endl;
// 				// displayNodeUt(per_node_vnf_resource);
// 				per_node_vnf_resource[node_map[j]][sfc_chain[j]]-= per_user_ut_node[sfc_chain[j]]*rate;
// 				// per_node_vnf_resource[node_map[j]][sfc_chain[j]]-= 1;
// 				op_cost += per_user_ut_node[sfc_chain[j]]*rate;
// 				if(slice==0){
// 					per_node_vnf_resource_slice0[node_map[j]][sfc_chain[j]]-= per_user_ut_node[sfc_chain[j]]*rate;
// 					// per_node_vnf_resource_slice0[node_map[j]][sfc_chain[j]]-= 1;
// 					// // cout << slice << " "  << per_node_vnf_resource[node_map[j]][sfc_chain[j]] << " " << per_user_ut_node[sfc_chain[j]]/service_type[slice].first<< endl;
// 					// int t00; cin >> t00;
// 					// if(per_node_vnf_resource_slice0[node_map[j]][sfc_chain[j]]<0){
// 					// 	// cout  << slice << "Req Resource:" << per_user_ut_node[sfc_chain[j]] << " " << j << " " << rate << endl;
// 					// 	// cout  << node_map[j] << " " <<sfc_chain[j] << endl;
// 					// 	// displayNodeUt(per_node_vnf_resource);
// 					// 	displayNodeUt(per_node_vnf_resource_slice0);
// 					// 	// displayNodeUt(per_node_vnf_resource_slice1);
// 					// 	cin >> flag_global;
// 					// }
// 				}
// 				else if(slice==1){
// 					per_node_vnf_resource_slice1[node_map[j]][sfc_chain[j]]-= per_user_ut_node[sfc_chain[j]]*rate;
// 					// per_node_vnf_resource_slice1[node_map[j]][sfc_chain[j]]-= 1;
// 					// // cout << slice <<" " <<per_node_vnf_resource[node_map[j]][sfc_chain[j]] << " " << per_user_ut_node[sfc_chain[j]]*service_type[slice].second <<endl;
// 					// int t00; cin >> t00;
// 				}
// 			}
//     		for(int i=1;i<route.size();i++){
//     			link_ut[route[i-1]][route[i]]+=per_user_ut_link;
//     			link_ut[route[i]][route[i-1]]+=per_user_ut_link;
//     			op_cost += per_user_ut_link;
//     			if(link_ut[route[i]][route[i-1]]>link_bw){
// 		    		flag_loaded=1;   					
//     				for(int j=i;j>=1;j--){
// 		    			link_ut[route[i-1]][route[i]]-=per_user_ut_link;
// 		    			link_ut[route[i]][route[i-1]]-=per_user_ut_link; 
//     				}
//     				fail_link_map++;
//     				break;
//     			}
//     		}
// 	    	if(flag_loaded){
// 	    		// break;
// 				if(slice==0) failed_s0++;
// 				if(slice==1) failed_s1++;
// 	    		fail_link_map++;
// 	    		recent_observation[index]=0;
// 	    		continue;
// 	    	}
//     	}

//     	//We increment the number of successful user
// 		// cout  << "Checkpoint-7" << endl;
//     	success++;
//     	recent_observation[index]=1;
//     	tot_cost+=op_cost;
//     	throughput += rate;
//     	// if(slice==0) reward+= service_type[slice].first/(latency*op_cost); //Change with current representation
//     	// else if(slice==1) reward += rate/(op_cost*service_type[slice].second); 
//     	if(slice==0) reward+= utility(service_type[slice].first,latency)/(op_cost); //Change with current representation
//     	else if(slice==1) reward += utility(rate,service_type[slice].second)/(op_cost); 
//     	// reward = reward/;
//     	// reward += rate/op_cost;
// 		// cout  << success << " " << fail_link_map<<" "<<fail_node_map << " " << reward <<endl;
// 		flag_completed=1;
// 	}
// 	return reward; //earlier I was sending throughput
// }


// void vertical_crossover(vector<pair<double,vector<pair<int,vector<int>>>>> &complete_node_map,vector<pair<int,vector<int>>>p1,vector<pair<int,vector<int>>>p2,vector<pair<double,double>> per_user_req,vector<vector<int>> per_user_sfc,vector<int> per_user_slice,vector<pair<int,int>> service_type,vector<vector<double>>per_node_vnf_resource,vector<vector<double>>per_node_vnf_resource_slice0,vector<vector<double>>per_node_vnf_resource_slice1,vector<double> per_user_ut_node,vector<vector<int>> reachable,vector<vector<vector<int>>>calculatedRoute,vector<vector<double>>link_ut,vector<int>f[])
// {
// 	int flag_crossover=0;
// 	int crossoverPoint = rand()%chain_length;
// 	// cout << "Inside crossover" << endl;
// 	// cout << p1.size() << " " << p2.size()  <<endl;
// 	// cin >> flag_global;
// 	for(int i=0;i<p1.size();i++){
// 		for(int j=0;j<p2.size();j++){
// 			// cout << i << " " <<  j << " a " << p1[i].first << " -> " << p2[j].first; cin >> flag_global;
// 			if(p1[i].first==p2[j].first){
// 				// cout << i << " " <<  j << " c "; cin >> flag_global;
// 				for(int t1=crossoverPoint-1;t1<chain_length;t1++){
// 					flag_crossover=1;
// 					swap((p1[i].second)[t1],(p2[j].second)[t1]);
// 				}
// 			}
// 			else if(p1[i].first<p2[j].first){
// 				break;
// 			}
// 			// cout << i << " " <<  j << " b "; cin >> flag_global;
// 		}
// 	}
// 	// cout << complete_node_map.size() << " ";
// 	if(flag_crossover){
// 		// cout << "I am here1"; cin >> flag_global;
// 		double fitness1 = fitnessEvaluation(p1,per_user_req,per_user_sfc,per_user_slice,service_type,per_node_vnf_resource,per_node_vnf_resource_slice0,per_node_vnf_resource_slice1,per_user_ut_node,reachable,calculatedRoute,link_ut,f);
// 		// cout << "I am here2"; cin >> flag_global;
// 		double fitness2 = fitnessEvaluation(p2,per_user_req,per_user_sfc,per_user_slice,service_type,per_node_vnf_resource,per_node_vnf_resource_slice0,per_node_vnf_resource_slice1,per_user_ut_node,reachable,calculatedRoute,link_ut,f);
// 		// cout << "I am here3"; cin >> flag_global;
// 		complete_node_map.push_back({fitness1,p1});
// 		// cout << "I am here4"; cin >> flag_global;
// 		complete_node_map.push_back({fitness2,p2}); 
// 	}
// 	// cout << complete_node_map.size() << endl;
// }

// void horizontal_crossover(vector<pair<double,vector<pair<int,vector<int>>>>> &complete_node_map,vector<pair<int,vector<int>>>p1,vector<pair<int,vector<int>>>p2,vector<pair<double,double>> per_user_req,vector<vector<int>> per_user_sfc,vector<int> per_user_slice,vector<pair<int,int>> service_type,vector<vector<double>>per_node_vnf_resource,vector<vector<double>>per_node_vnf_resource_slice0,vector<vector<double>>per_node_vnf_resource_slice1,vector<double> per_user_ut_node,vector<vector<int>> reachable,vector<vector<vector<int>>>calculatedRoute,vector<vector<double>>link_ut,vector<int>f[])
// {
// 	int flag_crossover=0;
// 	int size = min(p1.size(),p2.size());
// 	int crossoverPoint = rand()%size;
// 	// cout << "Inside crossover" << endl;
// 	// cout << p1.size() << " " << p2.size()  <<endl;
// 	// cin >> flag_global;
// 	for(int i=crossoverPoint;i<size;i++){
// 		swap(p1[i],p2[i]);
// 	}

// 	if(p1.size()>size){
// 		for(int i=0;i<p1.size();i++){
// 			p2.push_back(p1[i]);
// 		}
// 		while(p1.size()!=size){
// 			p1.pop_back();
// 		}
// 	}else{
// 		for(int i=0;i<p2.size();i++){
// 			p1.push_back(p2[i]);
// 		}
// 		while(p2.size()!=size){
// 			p2.pop_back();
// 		}
// 	}
// 	// cout << complete_node_map.size() << " ";
// 	if(flag_crossover){
// 		// cout << "I am here1"; cin >> flag_global;
// 		double fitness1 = fitnessEvaluation(p1,per_user_req,per_user_sfc,per_user_slice,service_type,per_node_vnf_resource,per_node_vnf_resource_slice0,per_node_vnf_resource_slice1,per_user_ut_node,reachable,calculatedRoute,link_ut,f);
// 		// cout << "I am here2"; cin >> flag_global;
// 		double fitness2 = fitnessEvaluation(p2,per_user_req,per_user_sfc,per_user_slice,service_type,per_node_vnf_resource,per_node_vnf_resource_slice0,per_node_vnf_resource_slice1,per_user_ut_node,reachable,calculatedRoute,link_ut,f);
// 		// cout << "I am here3"; cin >> flag_global;
// 		complete_node_map.push_back({fitness1,p1});
// 		// cout << "I am here4"; cin >> flag_global;
// 		complete_node_map.push_back({fitness2,p2}); 
// 	}
// 	// cout << complete_node_map.size() << endl;
// }

// void mutation(vector<pair<double,vector<pair<int,vector<int>>>>> &complete_node_map,vector<pair<int,vector<int>>>p1,vector<pair<double,double>> per_user_req,vector<vector<int>> per_user_sfc,vector<int> per_user_slice,vector<pair<int,int>> service_type,vector<vector<double>>per_node_vnf_resource,vector<vector<double>>per_node_vnf_resource_slice0,vector<vector<double>>per_node_vnf_resource_slice1,vector<double> per_user_ut_node,vector<vector<int>> reachable,vector<vector<vector<int>>>calculatedRoute,vector<vector<double>>link_ut,vector<int>f[])
// {
// 	int range1 = rand()%chain_length;
// 	int range2 = rand()%chain_length;
// 	range1 = min(range1,range2);
// 	range2 = max(range1,range2);
// 	if(range2==0){
// 		range1++;
// 		range2++;
// 	}
// 	// cout << "I am here in mutation with range1="<< range1<< " and range2="<< range2 << endl;
// 	// cout << "I am here"; cin >> flag_global;
// 	vector<pair<int,vector<int>>> node_map = p1;
// 	for(int j=0;j<node_map.size();j++){ // number of users mapped successfully using ith strategy
// 		vector<int> chain = per_user_sfc[node_map[j].first]; //returns user-id
// 		// for(int f1=0;f1<chain.size();f1++){
// 		for(int rng=range1;rng<=range2;rng++){
// 			vector<int> adjacencies = f[chain[rng]]; //
// 			int total_adjacencies= adjacencies.size();
// 			node_map[j].second[rng] = adjacencies[rand()%total_adjacencies];
// 		}
// 	}
// 	// cout << "I am done with mutation for " << i << endl;
// 	double fitness = fitnessEvaluation(node_map,per_user_req,per_user_sfc,per_user_slice,service_type,per_node_vnf_resource,per_node_vnf_resource_slice0,per_node_vnf_resource_slice1,per_user_ut_node,reachable,calculatedRoute,link_ut,f);
// 	complete_node_map.push_back({fitness,node_map});
// }





