//============================================================================
// Name        : proof_graph_analyzer.cpp
// Author      : 
// Version     :
//============================================================================


/*
 * Checks the merge resolution statistics for a dimacs file.
 *
 *
 */

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


// checks if resolvable clauses are merges in a cmty-based fashion
// if 1/n of a clause is in a cmty, count that fraction toward the community, both for clause count and merges

long compute_num_merges(vector< vector<int> >& clauses,
		vector<int>& cmty,
		vector<double>& cmty_merges,
		vector<double>& cmty_resolutions){
	long num_merges = 0;
	for(int i = 0; i < clauses.size() - 1; i++){
		vector<int> ci = clauses[i];

		for(int j = i + 1; j < clauses.size(); j++){
			vector<int> cj = clauses[j];

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
			// add to cmty_resolutions
			for(auto l: ci){
				cmty_resolutions[cmty[abs(l)]] += 1 / ((double) ci.size());
			}
			for(auto l: cj){
				cmty_resolutions[cmty[abs(l)]] += 1 / ((double) cj.size());
			}


			ni = 0;
			nj = 0;
			while(ni < ci.size() && nj < cj.size()){
				if(ci[ni] == cj[nj]){
					num_merges++;
					ni++;
					nj++;
					for(auto l: ci){
						cmty_merges[cmty[abs(l)]] += 1 / ((double) ci.size());
					}
					for(auto l: cj){
						cmty_merges[cmty[abs(l)]] += 1 / ((double) cj.size());
					}
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


	return num_merges;
}


int main(int argc, char * argv[]) {
	vector< vector<int> > clauses;
	vector<int> cmty;
	vector<int> cmty_picks;
	vector<int> cmty_size;
	vector<double> cmty_clauses;
	vector<double> cmty_merges;
	vector<double> cmty_resolutions;

	if(argc < 3){
		cout << "USAGE: ./proof_graph_analyzer cnf_file cmty_file" << endl;
		return 1;
	}
	char* dimacs_file = argv[1];
	char* cmty_file = argv[2];

	read_dimacs(dimacs_file, clauses);
	read_cmtys(cmty_file, cmty, cmty_picks, cmty_size, cmty_clauses);

	// add to cmty_clauses
	for(auto c: clauses){
		for(auto l: c){
			cmty_clauses[cmty[abs(l)]] += 1 / ((double) c.size());
		}
	}


	growTo(cmty_merges, cmty_clauses.size());
	growTo(cmty_resolutions, cmty_clauses.size());

	long num_pairwise_merges = compute_num_merges(clauses, cmty, cmty_merges, cmty_resolutions);
	printf("cmty clauses resolutions merges\n");
	for(int i = 0; i < cmty_clauses.size(); i++){
		printf("%d %8.4f %8.4f %8.4f\n", i, cmty_clauses[i], cmty_resolutions[i], cmty_merges[i]);
	}

	/*
	for(auto c : clauses){
		for(int i : c)
			cout<<i<<" ";
		cout<<endl;
	}
	*/
	int nVars = 0;
	for(auto c : clauses)
		for(auto l : c)
			if(abs(l) > nVars)
				nVars = abs(l);

	long clause_pairs = (long)clauses.size() * (clauses.size() - 1) / 2;
	cout<<dimacs_file<<","<<nVars<<","<<clauses.size()<<","<<clause_pairs<<","<<num_pairwise_merges<<","<<(double)num_pairwise_merges/clause_pairs<<endl;

	return 0;
}
