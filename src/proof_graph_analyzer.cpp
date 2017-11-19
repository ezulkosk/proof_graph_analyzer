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
#include <map>
#include <algorithm>


using namespace std;





void print_vector(vector<int>& v){
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


static double gini(vector<double>& vals){
	// compute gini coefficient of normalized picks
	for(int i = 0; i < vals.size() - 1; i++){
		for(int j = i + 1; j < vals.size(); j++){
			if(vals[i] > vals[j]){
				double temp = vals[i];
				vals[i] = vals[j];
				vals[j] = temp;
			}
		}
	}
	double height = 0;
	double area = 0;
	double fair_area = 0;
	for(int i = 0; i < vals.size(); i++){
		height += vals[i];
		area += height - (vals[i] / 2);
	}
	fair_area = height * (float(vals.size()) / 2);
	if(fair_area != 0)
		return (fair_area - area) / fair_area;
	else
		return -1;
}


void read_cmtys(char* cmty_file,
		vector<int>& cmty,
		vector<int>& cmty_picks,
		vector<int>& cmty_size,
		vector<double>& cmty_clauses){
	FILE* cfile = fopen(cmty_file, "r");
	if (cfile == NULL)
		fprintf(stderr, "could not open file %s\n", (const char*) cfile), exit(1);
	int v;
	int c;
	// File is zero-based, vector should store starting at 1
	cmty.push_back(-1);
	int largest_cmty_index = -1;
	while (fscanf(cfile, "%d %d\n", &v, &c) == 2){
		if(cmty.size() <= v+1)
			growTo(cmty, v+2);
		cmty[v+1] = c;
		if(c > largest_cmty_index){
			largest_cmty_index = c;
			growTo(cmty_picks, largest_cmty_index + 1);
			growTo(cmty_size, largest_cmty_index + 1);
			growTo(cmty_clauses, largest_cmty_index + 1);
		}
		cmty_size[c] = cmty_size[c] + 1;
	}
}


/*
 * Reads in the proof-graph output of maplesat structure_logging.
 *
 * @param file: -proof-graph output of maplesat structure_logging
 * @param graph: the DAG of the proof clause IDs, where ID = 0 is the root conflict
 * @param clauses: stores the literals in each clause starting with ID = 1
 * @param units: maps each unit clause to its ID
 * @param nOrigVars: the number of variables in the original CNF formula
 * @param nOrigClauses: the number of clauses in the original formula
 *
 */
void read_graph(char* file,
		vector< vector<int> >& graph,
		vector< vector<int> >& clauses,
		vector< pair<int, int> >& units,
		int& nOrigVars,
		int& nOrigClauses){
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
				if(curr_deps->size() == 0)
					nOrigClauses++;
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
		else if(mode == MODE_LITS){
			curr_clause->push_back(num);
			if(abs(num) > nOrigVars)
				nOrigVars = abs(num);
		}
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



/*
 * Given the full graph of clause dependencies, extract the proof containing only useful clauses.
 * Clears any entries in graph and clauses that are not in the proof.
 *
 * @param graph: the DAG of the proof clause IDs, where ID = 0 is the root conflict
 * @param clauses: stores the literals in each clause starting with ID = 1
 */
void graph_to_proof(vector< vector<int> >& graph, vector< vector<int> >& clauses){
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
			clauses[i].clear();
		}
	}
	cout<<"Count: "<<count << "/" <<graph.size()<<endl;
	cout<<"Deleted: "<<delete_count << "/" <<graph.size()<<endl;
}

void convert_graph_to_dimacs_format(vector<vector<int> >& graph, vector<vector<int> >& clauses, char* dimacs_out_file,
		int& nProofNodes, int& nProofEdges){
	// for each edge in the graph, convert to a dimacs clause to run features_s
	ofstream file;
	file.open (dimacs_out_file);

	nProofNodes = 0;
	nProofEdges = 0;
	map<int, int> m; // maps cids from the original graph to new cids, ignoring unused clauses
	for(int i = 0; i < graph.size(); i++){
		// the i == 0 case needs an id but does not have a clause, since it's the root conflict
		if(clauses[i].size() > 0 || i == 0){
			nProofNodes++;
			m[i] = nProofNodes;
			nProofEdges += graph[i].size();
		}
	}

	file<<"p cnf " << nProofNodes << " " << nProofEdges <<endl;
	for(int i = 0; i < graph.size(); i++){
		vector<int> adjacent_nodes = graph[i];
		for(auto j : adjacent_nodes){
		  file<<m[i] << " " << m[j]<< " 0"<<endl;
		}
	}
	file.close();
}

void convert_proof_to_dimacs_format(vector<vector<int> >& graph, vector<vector<int> >& clauses,
		char* dimacs_out_file, int& nOrigVars, int& nProofNodes){
	// for each clause in the proof, output to a dimacs clause to run features_s
	ofstream file;
	file.open (dimacs_out_file);

	file<<"p cnf " << nOrigVars << " " << nProofNodes <<endl;
	for(int i = 0; i < graph.size(); i++){
		if(clauses[i].size() > 0 || i == 0){
			for(auto j: clauses[i])
				file<<j << " ";
			file<<"0"<<endl;
		}
	}
	file.close();
}


/*
 * For each clause, output how many times it was used to derive another clause.
 * Outputs a list of (cid, usages) pairs.
 *
 * @param graph: the DAG of the proof clause IDs, where ID = 0 is the root conflict
 * @param clauses: stores the literals in each clause starting with ID = 1
 * @param proofClauseUses: stores the number of times each clause_i was used to derive another clause
 * @param proof_clause_uses_file: the file to output the usage pairs
 * @param nProofNodes: the number of clauses used in the proof
 *
 */
void getProofClauseUses(vector<vector<int> >& graph, vector<vector<int> >& clauses, vector<int>& proofClauseUses,
		     char* proof_clause_uses_file, int& nProofNodes){
	ofstream file;
	for(int i = 0; i < graph.size(); i++)
		proofClauseUses.push_back(0);
	file.open (proof_clause_uses_file);
	for(unsigned i = 0; i < proofClauseUses.size(); i++){
		if(clauses[i].size() > 0 || i == 0){
			for(int j : graph[i]){
				proofClauseUses[j]++;
			}

		}
	}
	for(unsigned i = 0; i < proofClauseUses.size(); i++)
		file<<i<<" "<<proofClauseUses[i]<<endl;

	file.close();
}


// output the spatial locality of the proof clauses, as defined in PoCR
void proofSpatialLocality(vector< vector<int> >& clauses,
		vector<int>& cmty, vector<int>& cmty_picks, vector<int>& cmty_size, vector<double>& cmty_clauses){
	for(auto c: clauses){
		int s = c.size();
		for(int i = 0; i < s; i++){
			int ci = cmty[abs(c[i])];
			cmty_clauses[ci] = cmty_clauses[ci] + (double(1) / s);
			cmty_picks[ci] += 1;
		}
	}
}


/*
 * Performs multiple analyses of the properties of proof clauses
 *
 * @param proofClauseUses: stores the number of times each clause_i was used to derive another clause
 * @param clauses: stores the literals in each clause starting with ID = 1
 *
 */
void analyzeProofClauseUses(vector<int>& proofClauseUses, vector< vector<int> >& clauses,
		vector<int>& cmty, vector<int>& cmty_picks, vector<int>& cmty_size, vector<double>& cmty_clauses,
		char* proof_analyses_file){
	assert(proofClauseUses[0] == 0);
	// analyze the number of uses of each clause, bucketized by size
	vector<int> size_bucket_uses;
	vector<int> size_bucket_occs; // how often clauses of a given size occurred
	ofstream file;

	for(int i = 1; i < proofClauseUses.size(); i++){
		vector<int> c = clauses[i];
		int cs = c.size();
		if(c.size() == 0)
			continue;
		growTo(size_bucket_uses, cs + 1);
		growTo(size_bucket_occs, cs + 1);
		size_bucket_uses[cs] += proofClauseUses[i];
		size_bucket_occs[cs] += 1;
	}
	cout<<"Clauses uses bucketed by size\n";
	cout<<"i uses occs uses/occs\n";
	for(int i = 0; i < size_bucket_uses.size(); i++){
		printf("%3d %9d %9d %9.4f\n", i, size_bucket_uses[i], size_bucket_occs[i], size_bucket_occs[i] == 0 ? 0 : (float)size_bucket_uses[i]/size_bucket_occs[i]);
	}

	file.open(proof_analyses_file);
	file<<"SizeBucketUses";
	for(int i = 0; i < size_bucket_uses.size(); i++)
		file<<","<<size_bucket_uses[i];
	file<<endl;
	file<<"SizeBucketOccs";
	for(int i = 0; i < size_bucket_occs.size(); i++)
		file<<","<<size_bucket_occs[i];
	file<<endl;
	file.close();
}



int main(int argc, char * argv[]) {
	vector< pair<int, int> > units; // (unit, unit_id) pairs
	vector< vector<int> > graph; // each index corresponds to a node, each element of graph[i] corresponds to deps[i]
	vector< vector<int> > proof;
	vector<int> proofClauseUses; // for each index i, maps how many times clause_i was used to derive another clause
	vector< vector<int> > clauses;
	vector<int> clauseUses; // TODO use
	vector<int> cmty, cmty_picks, cmty_size; // cmty of each variable (one-based)
	vector<double> cmty_clauses;


	if(argc < 6){
		cout << "USAGE: ./proof_graph_analyzer graph_file cnf_cmty_file graph_cnf_out_file proof_out_file proof_clause_uses_file proof_analyses_file" << endl;
		return 1;
	}
	char* file = argv[1];
	char* cmty_file = argv[2]; // cmtys of the original formula (zero-based)
	char* graph_cnf_file = argv[3];
	char* proof_out_file = argv[4];
	char* proof_clause_uses_file = argv[5];
	char* proof_analyses_file = argv[6];


	int nOrigVars = 0;
	int nOrigClauses = 0;
	int nProofNodes = 0; // corresponds to the number of clauses in the proof
	int nProofEdges = 0;


	read_graph(file, graph, clauses, units, nOrigVars, nOrigClauses);
	read_cmtys(cmty_file, cmty, cmty_picks, cmty_size, cmty_clauses);
	graph_to_proof(graph, clauses);
	convert_graph_to_dimacs_format(graph, clauses, graph_cnf_file, nProofNodes, nProofEdges);
	convert_proof_to_dimacs_format(graph, clauses, proof_out_file, nOrigVars, nProofNodes);
	getProofClauseUses(graph, clauses, proofClauseUses, proof_clause_uses_file, nProofNodes);
	analyzeProofClauseUses(proofClauseUses, clauses, cmty, cmty_picks, cmty_size, cmty_clauses, proof_analyses_file);

	proofSpatialLocality(clauses, cmty, cmty_picks, cmty_size, cmty_clauses);

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

	for(int i = 0; i < proofClauseUses.size(); i++){
		cout<<clauses[i].size()<<" "<<proofClauseUses[i];
		cout<<endl;
	}
	*/


	//for(auto i : units)
	//	cout<<i.first<<" "<<i.second<<endl;


	return 0;
}
