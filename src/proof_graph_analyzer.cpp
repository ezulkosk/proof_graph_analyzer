//============================================================================
// Name        : proof_graph_analyzer.cpp
// Author      : 
// Version     :
//============================================================================

#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <assert.h>
#include <utility>

using namespace std;



void print_vector(vector<int>& v){
	for(auto i: v){
		cout<<i<<" ";
	}
	cout<<endl;
}

/*
 * Stores the deps of the final conflict in graph[0]
 */
void read_graph(char* file,
		vector< vector<int> >& graph,
		vector< vector<int> >& clauses,
		vector< pair<int, int> >& units){
	int num;
	vector<int> dummy;
	vector<int>* final_deps = new vector<int>;
	vector<int>* curr_deps;
	vector<int>* curr_clause;
	graph.push_back(*final_deps);
	clauses.push_back(dummy);

	cout << file << endl;

	std::fstream graph_file(file, std::ios_base::in);

	int MODE_ID = 0;
	int MODE_LITS = 1;
	int MODE_DEPS = 2;
	int MODE_FINAL_CONFLICT_UNITS = 3; // get the unit clauses that were useful for the final conflict
	int MODE_FINAL_CONFLICT_DEPS = 4; // get the clauses that propagated in the final conflict
	int mode = MODE_ID;

	while (graph_file >> num)
	{
		if(mode == MODE_ID){
			if(num == 0){
				// final conflicting state
				mode = MODE_FINAL_CONFLICT_UNITS;
			}
			else{
				assert(num == graph.size());
				curr_deps = new vector<int>;
				curr_clause = new vector<int>;
				mode++;
			}
		}
		else if(!num){
			if(mode == MODE_LITS){
				clauses.push_back(*curr_clause);
				if(curr_clause->size() == 1){
					// need size - 1 since we have dummys in index 0
					pair<int, int>* p = new pair<int, int>(curr_clause->at(0), clauses.size() - 1);
					units.push_back(*p);
				}
			}
			else if(mode == MODE_DEPS){
				// for each dep id, get the clause
				//   for each lit in the clause, if it's a unit, add the unit_id to the deps
				for(int i = 0; i < curr_deps->size(); i++){
					int cid = curr_deps->at(i);
					vector<int> clause = clauses[cid];
					for(auto p : units){
						for(int lit : clause){
							if(p.first == -lit){
								curr_deps->push_back(p.second);
								break;
							}
						}
					}
				}
				graph.push_back(*curr_deps);
			}
			if(mode == MODE_FINAL_CONFLICT_UNITS)
				mode++;
			else
				mode = (mode + 1) % 3;

		}
		else if(mode == MODE_LITS)
			curr_clause->push_back(num);
		else if(mode == MODE_DEPS)
			curr_deps->push_back(num);
		else if(mode == MODE_FINAL_CONFLICT_UNITS){
			for(auto p : units){
				if(num == p.first){
					graph[0].push_back(p.second);
					break;
				}
			}
		}
		else if(mode == MODE_FINAL_CONFLICT_DEPS){
			graph[0].push_back(num);
		}
	}

}



// given the full graph of clause dependencies, extract the proof containing only useful clauses
void graph_to_proof(vector< vector<int> >& graph){
	vector<int>* current_deps = new vector<int>;
	vector<int> workpool;
	set<int> seen;
	cout<<"Running graph to proof"<<endl;
	for(int i : graph[0]){
		if(seen.find(i) == seen.end()){
			workpool.push_back(i);
			seen.insert(i);
		}
	}


	// while the workpool is not empty, get the ids deps and add the edges to the proof
	// copy deps from the original graph rather than duplicating
	int count = 0;
	workpool.push_back(0); // push the root conflict node
	while(workpool.size() > 0){
		int cid = workpool.back();
		seen.insert(cid);
		count++;
		workpool.pop_back();
		for(int i : graph[cid]){
			if(seen.find(i) == seen.end()){
				workpool.push_back(i);
				seen.insert(i);
			}
		}
	}

	// delete entries from the graph that were never seen
	int delete_count = 0;
	for(int i = 0; i < graph.size(); i++){
		if(seen.find(i) == seen.end()){
			delete_count++;
			graph[i].clear();
		}
	}
	cout<<"Count: "<<count << "/" <<graph.size()<<endl;
	cout<<"Deleted: "<<delete_count << "/" <<graph.size()<<endl;
}

void convert_graph_to_dimacs_format(vector<vector<int> > graph, char* dimacs_out_file, int& nVars, int& nClauses){
	// for each edge in the graph, convert to a dimacs clause to run louvain
	ofstream file;
	file.open (dimacs_out_file);

	nVars = graph.size();
	for(int i = 1; i < graph.size(); i++){
		nClauses += graph[i].size();
	}

	file<<"p cnf " << nVars << " " << nClauses <<endl;
	for(int i = 1; i < graph.size(); i++){
		vector<int> adjacent_nodes = graph[i];
		for(auto j : adjacent_nodes){
			file<<i << " " << j<< " 0\n";
		}
	}
	file.close();
}


int main(int argc, char * argv[]) {
	vector< pair<int, int> > units; // (unit, unit_id) pairs
	vector< vector<int> > graph; // each index corresponds to a node, each element of graph[i] corresponds to deps[i]
	vector< vector<int> > proof;
	vector< vector<int> > clauses;


	if(argc < 3){
		cout << "USAGE: ./proof_graph_analyzer graph_file dimacs_out_file" << endl;
		return 1;
	}
	char* file = argv[1];
	char* dimacs_out_file = argv[2];

	read_graph(file, graph, clauses, units);
	graph_to_proof(graph);
	int nVars = 0;
	int nClauses = 0;
	convert_graph_to_dimacs_format(graph, dimacs_out_file, nVars, nClauses);
	/*
	for(auto c : clauses){
		for(int i : c)
			cout<<i<<" ";
		cout<<endl;
	}

	for(auto c : graph){
		for(int i : c)
			cout<<i<<" ";
		cout<<endl;
	}

	for(auto i : units)
		cout<<i.first<<" ";
	cout<<endl;
	*/

	return 0;
}
