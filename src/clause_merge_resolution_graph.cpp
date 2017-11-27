//============================================================================
// Name        : proof_graph_analyzer.cpp
// Author      : 
// Version     :
//============================================================================


/*
 * Creates a graph where each clause is a node,
 * and an edge exists iff 2 clauses can be merged + resolved.
 *
 *
 */

#include <algorithm>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <fstream>


using namespace std;





void print_vector(vector<int>& v){
	for(auto i: v){
		cout<<i<<" ";
	}
	cout<<endl;
}

void print_vector(vector<double>& v){
	for(auto i: v){
		cout<<i<<" ";
	}
	cout<<endl;
}

void growTo(vector<int>& v, int s){
	if(v.size() < s){
		for(int i = v.size(); i < s; i++){
			v.push_back(0);
		}
	}
}

void growTo(vector<double>& v, int s){
	if(v.size() < s){
		for(int i = v.size(); i < s; i++){
			v.push_back(0);
		}
	}
}


void read_dimacs(char* dimacs_file, vector< vector<int> >& clauses){
	vector<int>* curr_clause;

	std::fstream stream(dimacs_file, std::ios_base::in);

	int num;
	string s;

	while (stream >> s)
	{
		if(s == "cnf")
			break;

	}
	// read the rest of the header
	stream >> num;
	stream >> num;
	curr_clause = new vector<int>;
	while(stream >> num){
		if(num == 0){
			std::sort(curr_clause->begin(), curr_clause->end());
			clauses.push_back(*curr_clause);
			curr_clause = new vector<int>;
		}
		else
			curr_clause->push_back(num);
	}

}



void compute_merge_resolution_edges(vector< vector<int> >& clauses, vector< vector< pair<int, int> > >& edges){
	for(int i = 0; i < clauses.size() - 1; i++){
		edges.push_back(*(new vector< pair<int, int> >));
		vector<int> ci = clauses[i];
		for(int j = i + 1; j < clauses.size(); j++){
			vector<int> cj = clauses[j];
			int num_merges = 0;
			// assumes clauses are sorted by read_dimacs()
			// check first if there is a resolution between the clauses
			unsigned ni = 0;
			unsigned nj = 0;
			bool canResolve = false;
			vector<int> cj2(cj.size());
			for(int k = 0; k < cj2.size(); k++)
				cj2[k] = -cj[cj.size() - 1 - k];
			while(ni < ci.size() && nj < cj2.size()){
				if(ci[ni] == cj2[nj]){
					canResolve = true;
					break;
				}
				else if(ni == ci.size() || nj == cj2.size())
					break;
				else if(ci[ni] < cj2[nj])
					ni++;
				else
					nj++;
			}
			if(!canResolve)
				continue;
			ni = 0;
			nj = 0;
			while(ni < ci.size() && nj < cj.size()){
				if(ci[ni] == cj[nj]){
					num_merges++;
					ni++;
					nj++;
				}
				else if(ni == ci.size() || nj == cj.size())
					break;
				else if(ci[ni] < cj[nj])
					ni++;
				else
					nj++;
			}
			if(num_merges > 0){
				edges[i].push_back(std::pair<int,int>(j, num_merges));
			}
		}

	}
}

int num_isolated_clauses(vector< vector< pair<int, int> > >& edges){
	int n = 0;
	for(auto v : edges){
		if(v.empty())
			n++;
	}
	return n;
}


int main(int argc, char * argv[]) {
	vector< vector<int> > clauses;

	if(argc < 2){
		cout << "USAGE: ./proof_graph_analyzer cnf_file" << endl;
		return 1;
	}
	char* dimacs_file = argv[1];

	read_dimacs(dimacs_file, clauses);

	vector< vector<pair<int, int> > > edges; // only store the edge from the lower index to higher

	compute_merge_resolution_edges(clauses, edges);
	int numIsolated = num_isolated_clauses(edges);


	/*
	for(auto c : clauses){
		for(int i : c)
			cout<<i<<" ";
		cout<<endl;
	}
	*/
	/*
	int nVars = 0;
	for(auto c : clauses)
		for(auto l : c)
			if(abs(l) > nVars)
				nVars = abs(l);

	*/
	int max = 0;
	int ci = 0;
	int cj = 0;
	for(int i = 0; i < edges.size(); i++){
		for(auto p : edges[i]){

			cout<<i<<" "<<p.first<<" " << p.second<<" 0"<<endl;
			if(p.second > max)
			{
				max = p.second;
				ci = i;
				cj = p.first;

			}
			//print_vector(clauses[i]);
			//print_vector(clauses[p.first]);
		}
	}
	cout<<"NumIsolated,"<<numIsolated<<endl;
	print_vector(clauses[ci]);
	print_vector(clauses[cj]);

	return 0;
}
