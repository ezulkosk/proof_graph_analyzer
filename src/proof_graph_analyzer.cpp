//============================================================================
// Name        : proof_graph_analyzer.cpp
// Author      : 
// Version     :
//============================================================================

#include <assert.h>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <utility>
#include <vector>
#include <math.h>
#include <queue>
#include "tools.h"

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

// check if 2 clauses have a conflicting literal
bool canBeResolved(vector<int>& l, vector<int>& r){
	unsigned ni = 0;
	unsigned nj = 0;
	vector<int> r2(r.size());
	for(int i = 0; i < r2.size(); i++)
		r2[i] = -r[r2.size() - 1 - i];
	while(ni < l.size() && nj < r2.size()){
		if(l[ni] == r2[nj]){
			return true;
		}
		else if(ni == l.size() || nj == r2.size())
			break;
		else if(l[ni] < r2[nj])
			ni++;
		else
			nj++;
	}
	return false;
}

// counts the number of merges that occur between two clauses
// assumes clauses are sorted by read_graph()
int numMerges(vector<int>& l, vector<int>& r){
	unsigned ni = 0;
	unsigned nj = 0;
	int merges = 0;
	while(ni < l.size() && nj < r.size()){
		if(l[ni] == r[nj]){
			merges++;
			ni++;
			nj++;
		}
		else if(ni == l.size() || nj == r.size())
			break;
		else if(l[ni] < r[nj])
			ni++;
		else
			nj++;
	}
	return merges;
}



static double gini(vector<double> vals){
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
		vector<int>& lbd,
		vector<int>& conflict_merges,
		int& nOrigVars,
		int& nOrigClauses){
	int num;
	vector<int> dummy;
	vector<int>* final_deps = new vector<int>;
	vector<int>* curr_deps;
	vector<int>* curr_clause;
	graph.push_back(*final_deps);
	clauses.push_back(dummy);
	lbd.push_back(0);
	conflict_merges.push_back(0);

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
			graph_file >> num;
			lbd.push_back(num);
			graph_file >> num;
			conflict_merges.push_back(num);
			graph_file >> num;
			assert(num == 0);
		}
		else if(!num){
			if(mode == MODE_LITS){
				std::sort((*curr_clause).begin(), (*curr_clause).end());
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


/**
 * Records many properties of the clause at the given index
 */
void logClauseProperties(vector< vector<int> >& graph, vector< vector<int> >& clauses, vector<int>& lbd, vector< int >& cmty,
		set<int>& seen, int proofClauseUsesIndex, vector<double>& var_pop, vector<double>& lit_pop, int index, ofstream& log){
	int useful = 1;
	if(seen.find(index) == seen.end())
		useful = 0;

	// record the size of the clause
	vector<int> clause = clauses[index];
	int size = clause.size();

	// cmty span
	set<int> cmtySpanSet;
	for(auto l : clause){
		cmtySpanSet.insert(cmty[abs(l)]);
	}
	int cmtySpan = cmtySpanSet.size();

	// weighted cmty span
	vector<int> cmtySpanVec;
	for(auto l : clause){
		cmtySpanVec.push_back(cmty[abs(l)]);
	}
	std::sort(cmtySpanVec.begin(), cmtySpanVec.end());
	int prev = cmtySpanVec[0];
	vector<double> cmtyCounts;
	int count = 1;
	cmtyCounts.push_back(count);
	for(int i = 1; i < cmtySpanVec.size(); i++){
		if(cmtySpanVec[i] == prev)
			cmtyCounts[cmtyCounts.size()-1] += 1;
		else{
			cmtyCounts.push_back(1);
			prev = cmtySpanVec[i];
		}
	}

	double weightedCmtySpan = gini(cmtyCounts);

	// time span of oldest dep
	int minDep = graph.size();
	for(int i : graph[index])
		if(i < minDep)
			minDep = i;
	int timeSpanOfOldestDep;
	if(graph[index].size() > 0)
		timeSpanOfOldestDep = index - minDep;
	else
		timeSpanOfOldestDep = 0;

	// avg time span of deps
	int totalTimeSpanOfDeps = 0;
	int totalDeps = 0;
	for(int i : graph[index]){
		totalTimeSpanOfDeps += index - i;
		totalDeps++;
	}

	double avgTimeSpanOfDeps;
	if(totalDeps > 0)
		avgTimeSpanOfDeps = ((double)totalTimeSpanOfDeps) / totalDeps;
	else
		avgTimeSpanOfDeps = 0;

	// set of vars in deps
	// TODO compare to total vars in deps, and also do for lits
	set<int> varsInDeps;
	for(int i : graph[index]){
		vector<int> cl = clauses[i];
		for(int l : cl){
			varsInDeps.insert(abs(l));
		}
	}
	int totalVarsInDeps = varsInDeps.size();


	// merges of deps
	int merge_count = 0;
	int sizeOfAllDeps = 0;

	for(int i = 0; i < ((int)graph[index].size()) - 1; i++){

		vector<int> ci = clauses[i];
		for(unsigned j = i + 1; j < graph[index].size(); j++){
			vector<int> cj = clauses[j];

			sizeOfAllDeps += cj.size();
			if(canBeResolved(ci, cj)){
				// assumes clauses are sorted by read_graph()
				unsigned ni = 0;
				unsigned nj = 0;
				while(ni < ci.size() && nj < cj.size()){
					if(ci[ni] == cj[nj]){
						merge_count++;
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
			}
		}
	}
	int mergeRatio;
	if(sizeOfAllDeps > 0)
		mergeRatio = merge_count / sizeOfAllDeps;
	else
		mergeRatio = 0;

	// average popularity of clause vars and literals
	double clauseVarPop = 0;
	double clauseLitPop = 0;
	double worstLitPop = 1;
	double secondWorstLitPop = 1;
	double worstVarPop = 1;
	double secondWorstVarPop = 1;

	double v_temp = 0;
	double l_temp = 0;
	for(int l : clause){
		if(l < 0){
			l = abs(l);
			v_temp = var_pop[l];
			l_temp = lit_pop[(2*l) - 1];
		}
		else{
			v_temp = var_pop[l];
			l_temp = lit_pop[2*l];
		}
		if(v_temp < worstVarPop){
			secondWorstVarPop = worstVarPop;
			worstVarPop = v_temp;
		}
		else if(v_temp < secondWorstVarPop)
			secondWorstVarPop = v_temp;
		if(l_temp < worstLitPop){
			secondWorstLitPop = worstLitPop;
			worstLitPop = l_temp;
		}
		else if(l_temp < secondWorstLitPop)
			secondWorstLitPop = l_temp;
		clauseVarPop += v_temp;
		clauseLitPop += l_temp;
	}
	clauseLitPop /= clause.size();
	clauseVarPop /= clause.size();
	double worst2VarPops = worstVarPop + secondWorstVarPop;
	double worst2LitPops = worstLitPop + secondWorstLitPop;


	log<<index<<","
			<<useful<<","
			<<proofClauseUsesIndex<<","
			<<lbd[index]<<","
			<<size<<","
			<<cmtySpan<<","
			<<weightedCmtySpan<<","
			<<timeSpanOfOldestDep<<","
			<<avgTimeSpanOfDeps<<","
			<<totalVarsInDeps<<","
			<<totalDeps<<","
			<<clauseVarPop<<","
			<<clauseLitPop<<","
			<<worst2VarPops<<","
			<<worst2LitPops<<","
			<<mergeRatio<<endl;
}

/*
 * Normalize a vector of doubles such that they sum to 1.
 */
void normalize_vector(vector<double>& v){
	double weight = 0;
	for(double d : v)
		weight += d;
	for(int i = 0; i < v.size(); i++){
		v[i] = v[i] / weight;
	}
}


void get_var_and_lit_popularity(vector< vector<int> >& graph, vector< vector<int> >& clauses,
		vector<double>& var_pop, vector<double>& lit_pop)
{
	for(int i = 0; i < graph.size(); i++){
		if(graph[i].size() == 0){
			// no deps, must be an input clause
			vector<int> clause = clauses[i];
			for(auto l : clause){
				if(l < 0){
					l = abs(l);
					lit_pop[(2*l)-1]++;
					var_pop[l]++;
				}
				else{
					lit_pop[2*l]++;
					var_pop[l]++;
				}
			}
		}
	}
	normalize_vector(var_pop);
	normalize_vector(lit_pop);
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
void getProofClauseUses(vector<vector<int> >& graph, vector<vector<int> >& clauses, set<int>& seen,
		vector<int>& proofClauseUses){
	for(int i = 0; i < graph.size(); i++)
		proofClauseUses.push_back(0);
	for(unsigned i = 0; i < proofClauseUses.size(); i++){
		if(seen.find(i) != seen.end() || i == 0){
			for(int j : graph[i]){
				proofClauseUses[j]++;
			}

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
void graph_to_proof(vector< vector<int> >& graph, vector< vector<int> >& clauses, vector<int>& lbd, vector<int>& cmty,
		vector<int>& proofClauseUses, vector<double>& var_pop, vector<double>& lit_pop, ofstream& log){
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

	//getProofClauseUses(graph, clauses, seen, proofClauseUses);

	/*
	clausePropertiesLog.open(clausePropertiesFile);
	clausePropertiesLog<<"id,useful,uses,lbd,size,cmty_span,weighted_cmty_span,time_span_of_oldest_dep,"<<
			"avg_time_span_of_deps,set_of_vars_in_deps_size,total_deps,popularity_of_vars,"<<
			"popularity_of_literals,var_pop_worst2,lit_pop_worst2,merge_ratio"<<endl;

	for(int i = 1; i < graph.size(); i++){
		logClauseProperties(graph, clauses, lbd, cmty, seen, proofClauseUses[i], var_pop, lit_pop, i, clausePropertiesLog);
	}
	clausePropertiesLog.close();
	*/

	// compute 1) the number of useless input clauses;
	// 2) the number of useless clauses with no useless deps but have useful deps (but not input clauses)
	int num_useless_input = 0;
	int num_useless_fringe = 0; //if we cut these, the other useless clauses are disconnected from the proofs

	// delete entries from the graph that were never seen
	int delete_count = 0;
	for(int i = 0; i < graph.size(); i++){
		if(seen.find(i) == seen.end()){
			delete_count++;
			if(graph[i].empty()){
				num_useless_input++;
			}
			else{
				for(int j = 0; j < graph[i].size(); j++){
					if(seen.find(graph[i][j]) != seen.end()){
						num_useless_fringe++;
						break;
					}
				}
			}


			graph[i].clear();
			clauses[i].clear();
		}
	}

	log<<"Useful: "<<count << " / " <<graph.size()<<endl;
	log<<"Deleted: "<<delete_count << " / " <<graph.size()<<endl;
	log<<"UselessInput: "<<num_useless_input << " / " <<delete_count<<endl;
	log<<"FringeUseless: "<<num_useless_fringe << " / " <<delete_count<<endl;
	log<<"CuttableUseless: "<<(delete_count - num_useless_fringe - num_useless_input) << " / " <<delete_count<<endl;


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

	nProofNodes = 0;
	for(int i = 1; i < graph.size(); i++){
			if(clauses[i].size() > 0)
				nProofNodes++;
	}

	file<<"p cnf " << nOrigVars << " " << nProofNodes <<endl;
	for(int i = 1; i < graph.size(); i++){
		if(clauses[i].size() > 0 || i == 0){
			for(auto j: clauses[i])
				file<<j << " ";
			file<<"0"<<endl;
		}
	}


	file.close();
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
		char* proof_analyses_file, ofstream& proofAnalysesFile){
	assert(proofClauseUses[0] == 0);
	// analyze the number of uses of each clause, bucketized by size
	vector<int> size_bucket_uses;
	vector<int> size_bucket_occs; // how often clauses of a given size occurred

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

	proofAnalysesFile<<"SizeBucketUses";
	for(int i = 0; i < size_bucket_uses.size(); i++)
		proofAnalysesFile<<","<<size_bucket_uses[i];
	proofAnalysesFile<<endl;
	proofAnalysesFile<<"SizeBucketOccs";
	for(int i = 0; i < size_bucket_occs.size(); i++)
		proofAnalysesFile<<","<<size_bucket_occs[i];
	proofAnalysesFile<<endl;

	// analyze the number of uses of each clause, bucketized by cmty span (not normalized by size)
	vector<int> cmtySpanBucketUses;
	vector<int> cmtySpanBucketOccs; // how often clauses of a given size occurred

	for(int i = 1; i < proofClauseUses.size(); i++){
		vector<int> c = clauses[i];
		// get the number of cmtys that the clause spans
		if(c.size() == 0)
			continue;
		int cs = 0;
		set<int> cmtysInClause;
		cmtysInClause.clear();

		for(auto l : c){
			int curr_cmty = cmty[abs(l)];
			cmtysInClause.insert(curr_cmty);
		}
		cs = cmtysInClause.size();

		growTo(cmtySpanBucketUses, cs + 1);
		growTo(cmtySpanBucketOccs, cs + 1);
		cmtySpanBucketUses[cs] += proofClauseUses[i];
		cmtySpanBucketOccs[cs] += 1;
	}
	cout<<"Clauses uses bucketed by cmty span\n";
	cout<<"i uses occs uses/occs\n";
	for(int i = 0; i < cmtySpanBucketUses.size(); i++){
		printf("%3d %9d %9d %9.4f\n", i, cmtySpanBucketUses[i], cmtySpanBucketOccs[i], cmtySpanBucketOccs[i] == 0 ? 0 : (float)cmtySpanBucketUses[i]/cmtySpanBucketOccs[i]);
	}

	proofAnalysesFile<<"CmtySpanBucketUses";
	for(int i = 0; i < cmtySpanBucketUses.size(); i++)
		proofAnalysesFile<<","<<cmtySpanBucketUses[i];
	proofAnalysesFile<<endl;
	proofAnalysesFile<<"CmtySpanBucketOccs";
	for(int i = 0; i < cmtySpanBucketOccs.size(); i++)
		proofAnalysesFile<<","<<cmtySpanBucketOccs[i];
	proofAnalysesFile<<endl;
}


// output the spatial locality of the proof clauses, as defined in PoCR
void proofSpatialLocality(vector< vector<int> >& clauses,
		vector<int>& cmty, vector<int>& cmty_picks, vector<int>& cmty_size, vector<double>& cmty_clauses,
		ofstream& proofAnalysesFile){
	for(auto c: clauses){
		int s = c.size();
		for(int i = 0; i < s; i++){
			int ci = cmty[abs(c[i])];
			cmty_clauses[ci] = cmty_clauses[ci] + (double(1) / s);
			cmty_picks[ci] += 1;
		}
	}


	vector<double> cmtyOccs;
	for(auto e : cmty_picks)
		cmtyOccs.push_back(e);
	// non-normalized by size
	proofAnalysesFile<<"GiniCmtyOccs,"<<gini(cmtyOccs)<<endl;
	proofAnalysesFile<<"GiniCmtyClauses,"<<gini(cmty_clauses)<<endl;

	// normalized by size
	vector<double> cmtyOccsNormalized;
	for(int i = 0; i < cmtyOccs.size(); i++){
		if(cmty_size[i] == 0)
			cmtyOccsNormalized.push_back(0);
		else
			cmtyOccsNormalized.push_back(cmtyOccs[i] / cmty_size[i]);
	}
	proofAnalysesFile<<"GiniCmtyOccsNormalizedBySize,"<<gini(cmtyOccsNormalized)<<endl;

	vector<double> cmtyClausesNormalized;
	for(int i = 0; i < cmty_clauses.size(); i++){
		if(cmty_size[i] == 0)
			cmtyClausesNormalized.push_back(0);
		else
			cmtyClausesNormalized.push_back((double)cmty_clauses[i] / cmty_size[i]);
	}
	proofAnalysesFile<<"GiniCmtyClausesNormalizedBySize,"<<gini(cmtyClausesNormalized)<<endl;
}


// output the merge locality of the proof clauses
void proofMergeLocality(vector< vector<int> >& graph, vector< vector<int> >& clauses, vector<int>& conflict_merges, ofstream& proofAnalysesFile){
	// for each clause in the proof, output how many merges occur between resolvable clauses of its deps
	long merges = 0;
	double mergesNormalizedByNumDeps = 0;
	double mergesNormalizedByNumDepsSquared = 0;
	long pf_size = 0;
	long numDeps = 0;
	for(int i = 0; i < clauses.size(); i++){
		vector<int> learnt = clauses[i];
		vector<int> deps = graph[i];
		if(learnt.size() == 0 || i == 0) // OLD: && i != 0)
			continue;
		pf_size++;
		/*
		for(int j = 0; j < ((int)deps.size()) - 1; j++){
			vector<int> dj = clauses[deps[j]];
			for(int k = j + 1; k < ((int) deps.size()); k++){
				vector<int> dk = clauses[deps[k]];
				if(canBeResolved(dj, dk)){
					int m = numMerges(dj, dk);
					merges += m;
					curr_merges += m;
				}
			}
		}
		*/
		merges += conflict_merges[i];


		if(deps.size() > 0){
			mergesNormalizedByNumDeps += (double) conflict_merges[i] / deps.size();
			if(deps.size() > 1)
				mergesNormalizedByNumDepsSquared += (double) conflict_merges[i] / ((deps.size() * (deps.size()-1))/2);

			numDeps += deps.size();
		}
	}

	proofAnalysesFile<<"Merges," << merges<<endl;
	proofAnalysesFile<<"MergeLocalityAverage,"<<((double) merges) / pf_size<<endl;
	proofAnalysesFile<<"MergeLocalityNormalizedByNumDeps,"<<mergesNormalizedByNumDeps / pf_size<<endl;
	proofAnalysesFile<<"MergeLocalityNormalizedByNumDepsSquared,"<<mergesNormalizedByNumDepsSquared / pf_size<<endl;
	proofAnalysesFile<<"AverageDeps,"<<(double)numDeps / pf_size<<endl;


}


void analyzeProofLSRPercentiles(vector< vector<int> >& graph,
		vector< vector<int> >& clauses,
		ofstream& proofAnalysesFile){
	set<int> lsr_vars;
	vector<int> lsr_over_time;
	for(int i = 1; i < graph.size(); i++){
		if(clauses[i].size() != 0){
			if(graph[i].size() == 0)
				continue;
			else{
				for(int j : clauses[i])
					lsr_vars.insert(abs(j));
				lsr_over_time.push_back(lsr_vars.size());
			}
		}
	}
	double percentile_inc = 0.25 * lsr_vars.size();
	double percentile = percentile_inc;
	cout<<"LSR25Checkpoints";
	for(int i = 0; i < lsr_over_time.size(); i++){
		if(percentile > lsr_vars.size()+1){
			cout<<endl;
			break;
		}
		if(lsr_over_time[i] >= percentile){
			printf(",%0.7f", (double)i / lsr_over_time.size());
			percentile += percentile_inc;

		}
	}
	//cout<<"Percentile:"<<percentile;
	cout<<"LSR,"<<lsr_vars.size()<<endl;
}

void analyzeProofLSRAvg(vector< vector<int> >& graph,
		vector< vector<int> >& clauses,
		ofstream& proofAnalysesFile){

	/*
	 * Optimization: for each clause, find it's last direct usage, record in vector
	 *
	 * When i exceeds the last use of a clause, delete it's associated lsr_set.
	 */

	int total_learnts = 0;
	long total_lsr = 0;
	long total_clause_proof_size = 0;

	set<int> lsr_vars;
	long curr_clause_proof_size = 0;

	// for each clause, find it's last usage
	vector<int> last_usage_vec;
	map<int, vector<int> > dead_clauses;
	for(int i = 0; i < clauses.size(); i++)
		last_usage_vec.push_back(0);
	for(int i = 1; i < graph.size(); i++)
		for(int j: graph[i])
			last_usage_vec[j] = i;
	// store the last usages in a map
	for(int i = 1; i < graph.size(); i++){
		int usage = last_usage_vec[i];
		if(dead_clauses.find(usage + 1) == dead_clauses.end()){
			dead_clauses[usage + 1] = *(new vector<int>);
		}
		dead_clauses[usage + 1].push_back(i);
	}
	last_usage_vec.clear();

	map<int, set<int> > lsr_sets;




	for(int i = 1; i < graph.size(); i++){
		if(dead_clauses.find(i) != dead_clauses.end()){
			for(int j: dead_clauses[i]){
				lsr_sets[j].clear();
			}
		}
		if(clauses[i].size() != 0){
			if(graph[i].size() == 0)
				continue;
			else{
				//cout<<"PfLSR: "<<i<<endl;
				total_learnts++;
				//curr_clause_proof_size = 0;
				lsr_vars.clear();
				for(int j: clauses[i])
					lsr_vars.insert(abs(j));
				for(int cid: graph[i]){
					// INV: lsr_vars are always positive, don't need to do abs() on deps.
					for(int j : lsr_sets[cid])
						lsr_vars.insert(j);
				}
				total_lsr += lsr_vars.size();
				set<int>* curr_set = new set<int>(lsr_vars);
				lsr_sets[i] = *curr_set;
				//total_clause_proof_size += curr_clause_proof_size;

			}
		}
	}
	cout<<"TotalLSR,"<<total_lsr<<endl;
	cout<<"LSRAvg,"<<((double) total_lsr)/total_learnts<<endl;
	//cout<<"AverageProofSizeInClauses,"<<((double)total_clause_proof_size)/total_learnts<<endl;
}




void compareInputPopWithProofPop(vector< vector<int> >& graph,
		vector< vector<int> >& clauses,
		vector<double>& var_pop,
		vector<double>& lit_pop,
		ofstream& proofAnalysesFile){
	vector<double> proof_pop;
	vector<double> proof_lit_pop;
	growTo(proof_pop, var_pop.size());
	growTo(proof_lit_pop, lit_pop.size());
	cout<<"IN comp\n";
	// i = 0 is the final conflict
	for(int i = 1; i < graph.size(); i++){
		if(clauses[i].size() != 0){
			// the clause is involved in the proof
			vector<int> clause = clauses[i];
			for(auto l : clause){
				if(l < 0){
					l = abs(l);
					proof_pop[l]++;
					proof_lit_pop[(2*l)-1]++;
				}
				else{
					proof_pop[l]++;
					proof_lit_pop[2*l]++;
				}
			}
		}
	}
	normalize_vector(proof_pop);
	normalize_vector(proof_lit_pop);
	double var_amt = 0;
	for(int i = 0; i < var_pop.size(); i++){
		var_amt += fabs((var_pop[i] - proof_pop[i]) * (var_pop[i] - proof_pop[i]));
	}
	var_amt /= var_pop.size();
	var_amt = sqrt(var_amt);

	double lit_amt = 0;
	for(int i = 0; i < lit_pop.size(); i++){
		lit_amt += fabs((lit_pop[i] - proof_lit_pop[i]) * (lit_pop[i] - proof_lit_pop[i]));
	}
	lit_amt /= lit_pop.size();
	lit_amt = sqrt(lit_amt);

	//proofAnalysesFile<<"NormalizedInputAndProofPopDiff,"<<total_diff/var_pop.size()<<endl;
	proofAnalysesFile<<"VarPopProofRMSE,"<<var_amt<<endl;
	proofAnalysesFile<<"LitPopProofRMSE,"<<lit_amt<<endl;


}

double fixed_cmty_modularity(Graph* g, vector <int>& n2c) {
//-------------------------------------------------------------------------------------------
// Given the graph "g" and the partition "n2c" computes the modularity
//-------------------------------------------------------------------------------------------
	double w = 0;
	double arity = 0;
	vector <double> aritym(g->size(), 0);

	//std::cout<<"arity: " << g->size()<<std::endl;
	for (Graph::EdgeIter it=g->begin(); it != g->end(); it++) {
		//assert(it->orig >= 0 && it->orig < n2c.size());
		//assert(it->dest >= 0 && it->dest < n2c.size());
		//std::cout<<it->orig<<std::endl;
		//std::cout<< it->orig<< " " << it->dest <<std::endl;
		if (n2c[it->orig] == n2c[it->dest])
			w += it->weight;
	}

	for (int i=0; i<g->size(); i++)
		aritym[n2c[i]] += g->arity(i);

	for (int i=0; i<g->size(); i++)
		arity += aritym[i] * aritym[i] / g->arity() / g->arity();


//  cerr <<" Modularity = "<<w / g->arity()<<" - "<<arity<<" = "<<w / g->arity() - arity<<endl;
	return 2*w / g->arity() - arity;
}


int main(int argc, char * argv[]) {
	vector< pair<int, int> > units; // (unit, unit_id) pairs
	vector< vector<int> > graph; // each index corresponds to a node, each element of graph[i] corresponds to deps[i]
	vector< vector<int> > proof;
	vector<int> lbd;
	vector<int> conflict_merges;
	vector<int> proofClauseUses; // for each index i, maps how many times clause_i was used to derive another clause
	vector< vector<int> > clauses;
	vector<int> clauseUses; // TODO use
	vector<int> cmty, cmty_picks, cmty_size; // cmty of each variable (one-based)
	vector<double> cmty_clauses;
	vector<double> var_pop;
	vector<double> lit_pop; // +i at index 2i, -i at index 2i - 1
	set<int> seen;

	if(argc < 3){
		cout << "USAGE: ./proof_graph_analyzer graph_file cnf_cmty_file proof_analyses_file clause_properties file" << endl;
		return 1;
	}
	char* graph_file = argv[1];
	char* cmty_file = argv[2]; // cmtys of the original formula (zero-based)
	char* proof_analyses_file = argv[3];
	//char* clause_properties_file = argv[4];
	//char* proof_out_file = argv[5];

	int nOrigVars = 0;
	int nOrigClauses = 0;
	int nProofNodes = 0; // corresponds to the number of clauses in the proof
	int nProofEdges = 0;


	read_graph(graph_file, graph, clauses, units, lbd, conflict_merges, nOrigVars, nOrigClauses);
	growTo(var_pop, nOrigVars+1);
	growTo(lit_pop, (nOrigVars+1) * 2);
	get_var_and_lit_popularity(graph, clauses, var_pop, lit_pop);

	read_cmtys(cmty_file, cmty, cmty_picks, cmty_size, cmty_clauses);
	//<Graph*,Graph*> p = readFormula("/home/ezulkosk/cp2017_benchmarks/agile/simp/cnf/bench_8784.smt2.cnf", 400);
	//Graph* vig = p.first;

	vector<int> n2c;
	for(int i = 1; i < cmty.size(); i++){
		n2c.push_back(cmty[i]);
	}

	//double q = fixed_cmty_modularity(vig, n2c);
	//cout<<q<<endl;

	//cout<<"after\n";

	ofstream proofAnalysesFile;
	proofAnalysesFile.open(proof_analyses_file);

	graph_to_proof(graph, clauses, lbd, cmty, proofClauseUses, var_pop, lit_pop, proofAnalysesFile);

	//convert_graph_to_dimacs_format(graph, clauses, graph_cnf_file, nProofNodes, nProofEdges);

	// convert_proof_to_dimacs_format(graph, clauses, proof_out_file, nOrigVars, nProofNodes);
	getProofClauseUses(graph, clauses, seen, proofClauseUses);


	cout<<"after\n";




	// 	analyzeProofLSRPercentiles(graph, clauses, proofAnalysesFile);
	//analyzeProofLSRAvg(graph, clauses, proofAnalysesFile);

	analyzeProofClauseUses(proofClauseUses, clauses, cmty, cmty_picks, cmty_size, cmty_clauses, proof_analyses_file, proofAnalysesFile);
	proofSpatialLocality(clauses, cmty, cmty_picks, cmty_size, cmty_clauses, proofAnalysesFile);

	proofMergeLocality(graph, clauses, conflict_merges, proofAnalysesFile);

	compareInputPopWithProofPop(graph, clauses, var_pop, lit_pop, proofAnalysesFile);

	/*
	pair<Graph*,Graph*> p = readFormula(proof_out_file, 400);
	Graph* vig = p.first;
	double q = fixed_cmty_modularity(vig, n2c);
	cout<<q<<endl;
	proofAnalysesFile<<"ProofModularityOriginalPartition,"<<q<<endl;
	proofAnalysesFile.close();
	*/

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
