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

// TODO UNITS
void read_graph(char* file,
		vector< vector<int> >& graph,
		vector< vector<int> >& clauses,
		vector< pair<int, int> >& units){
	int num;
	vector<int> dummy;
	vector<int>* curr_deps;
	vector<int>* curr_clause;
	graph.push_back(dummy);
	clauses.push_back(dummy);

	cout << file << endl;

	std::fstream graph_file(file, std::ios_base::in);

	int MODE_ID = 0;
	int MODE_LITS = 1;
	int MODE_DEPS = 2;
	int mode = MODE_ID;

	while (graph_file >> num)
	{
		if(mode == MODE_ID){
			assert(num == graph.size());
			curr_deps = new vector<int>;
			curr_clause = new vector<int>;
			mode++;
		}
		else if(!num){
			if(mode == MODE_LITS){
				clauses.push_back(*curr_clause);
				if(curr_clause->size() == 1){
					pair<int, int>* p = new pair<int, int>(curr_clause->at(0), clauses.size());
					units.push_back(*p);
				}
			}
			else if(mode == MODE_DEPS){
				// for each dep id, get the clause
				//   for each lit in the clause, if it's a unit, add the unit_id to the deps
				for(int i : *curr_deps){
					vector<int> clause = clauses[i];
					for(auto p : units){
						for(int lit : clause){
							if(p.first == -lit){
								//cout<<"Found unit: " << lit << endl;
								//for(int l : clause)
								//	cout<<l<<" ";
								//cout << endl;
								//cout<<"Lits in current clause: ";
								//for(int l : *curr_clause)
								//	cout<<l<<" ";
								//cout << endl;
								curr_deps->push_back(p.second);
								break;
							}
						}
					}
				}
				graph.push_back(*curr_deps);
			}
			mode = (mode + 1) % 3;
		}
		else if(mode == MODE_LITS)
			curr_clause->push_back(num);
		else if(mode == MODE_DEPS)
			curr_deps->push_back(num);
	}

}


int main(int argc, char * argv[]) {
	vector< pair<int, int> > units; // unit, unit_id pairs
	vector< vector<int> > graph; // each index corresponds to a node, each element of graph[i] corresponds to deps[i]
	vector< vector<int> > clauses;


	if(argc < 2){
		cout << "Must supply graph file" << endl;
		return 1;
	}
	char* file = argv[1];

	read_graph(file, graph, clauses, units);
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
	*/
	for(auto i : units)
		cout<<i.first<<" ";
	cout<<endl;

	return 0;
}
