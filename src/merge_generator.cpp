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

void growTo(vector<vector <vector<int>*> >& v, int s){
	if(v.size() < s){
		for(int i = v.size(); i < s; i++){
			vector< vector<int>* >* t = new vector<vector<int>* >;
			v.push_back(*t);
		}
	}
}


int toLit(int i){
	if(i > 0)
		return 2*i;
	else
		return 2*-i+1;
}


void is_res_merge(vector<int>& c1, vector<int>& c2, bool& res, int& num_merge){
	num_merge = 0;
	res = false;
	for(int i = 0; i < c1.size(); i++){
		int l1 = c1[i];
		for(int j = 0; j < c2.size(); j++){
			int l2 = c2[j];
			if(l1 == l2){
				num_merge++;
			}
			else if(l1 == -l2){
				res = true;
			}
		}
	}
}





void read_cmtys(char* cmty_file,
		vector<int>& cmty){
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
			for(int i = 0; i < curr_clause->size() - 1; i++){
				if(curr_clause->at(i) == curr_clause->at(i+1))
					cout<<"HIT!\n";
			}
			int low = 0;
			int high = curr_clause->size()-1;
			while(low < high){
				if(curr_clause->at(low) == curr_clause->at(high)){
					cout<<"HIT2!\n";
				}
				else if(curr_clause->at(low) > 0 || curr_clause->at(high) < 0)
					break;
				else if(abs(curr_clause->at(low)) > curr_clause->at(high))
					low++;
				else
					high--;
			}
			for(int i = 0; i < curr_clause->size() - 1; i++){
				if(curr_clause->at(i) == curr_clause->at(i+1))
					cout<<"HIT!\n";
			}

			clauses.push_back(*curr_clause);
			curr_clause = new vector<int>;
		}
		else
			curr_clause->push_back(num);
	}

}


/*
 * increases enabled and merge counters for appropriate literals
 */
bool res_merge_counter(vector<int>& c1, vector<int>& c2,
		vector<int>* e1, vector<int>* e2,
		vector<int>* m1, vector<int>* m2){

	// can assume that a literal will never be both a merge and a res for the same clause pair
	vector<int> res_indices1;
	vector<int> res_indices2;
	vector<int> merge_indices1;
	vector<int> merge_indices2;

	int num_merges = 0;

	for(int i = 0; i < c1.size(); i++){
		int l1 = c1[i];
		for(int j = 0; j < c2.size(); j++){
			int l2 = c2[j];
			if(l1 == l2){
				num_merges++;
				merge_indices1.push_back(i);
				merge_indices2.push_back(j);
			}
			else if(l1 == -l2){
				res_indices1.push_back(i);
				res_indices2.push_back(j);
			}
		}
	}
	if(res_indices1.size() > 0 && num_merges > 0){
		// print_vector(*m1);
		for(int i = 0; i < res_indices1.size(); i++){
			e1->at(res_indices1[i]) += num_merges;
		}

		for(int i = 0; i < res_indices2.size(); i++){
			e2->at(res_indices2[i]) += num_merges;
		}

		for(int i = 0; i < merge_indices1.size(); i++){
			m1->at(merge_indices1[i])++;
		}

		for(int i = 0; i < merge_indices2.size(); i++){
			m2->at(merge_indices2[i])++;
		}
		// print_vector(*m1);
		return true;
	}
	return false;
}


void compute_enabled_and_merges_on_lits(vector< vector<int> >& clauses,
		vector< vector<int> >& enabled,
		vector< vector<int> >& merges){

	for(int i = 0; i < (int) clauses.size() - 1; i++){
		vector<int> c1 = clauses[i];
		vector<int> e1 = enabled[i];
		vector<int> m1 = merges[i];
		for(int j = i + 1; j < (int) clauses.size(); j++){
			vector<int> c2 = clauses[j];
			vector<int> e2 = enabled[j];
			vector<int> m2 = merges[j];
			bool updated = res_merge_counter(c1, c2, &e1, &e2, &m1, &m2);
			if(updated){
				// print_vector(m1);
				merges[i] = m1;
				merges[j] = m2;
				enabled[i] = e1;
				enabled[j] = e2;
				// cout<<"===="<<endl;
			}
		}
	}
	cout<<endl;
	for(int i = 0; i < clauses.size(); i++){
			vector<int> e = enabled[i];
			vector<int> m = merges[i];
			for(int j = 0; j < e.size(); j++){
				if(e[j] + m[j] == 0)
					printf("%10d %10d %10d %10d\n", i, j, e[j], m[j]);
			}
		}
}


int main(int argc, char * argv[]) {
	vector< vector<int> > clauses;
	vector<int> cmty;

	if(argc < 4){
		cout << "USAGE: ./merge_generator cnf_file cmty_file out_file" << endl;
		return 1;
	}

	char* dimacs_file = argv[1];
	char* cmty_file = argv[2];
	char* out_file = argv[3];

	read_dimacs(dimacs_file, clauses);
	read_cmtys(cmty_file, cmty);

	vector< vector<int> > enabled; // vector to store the number of enabled merges of every individual literal
	vector< vector<int> > merges; // vector to store the number of actual merges on every individual literal
	for(auto c : clauses){
		vector<int> vec = *(new vector<int>);
		vector<int> vec2 = *(new vector<int>);
		growTo(vec, c.size());
		growTo(vec2, c.size());
		enabled.push_back(vec);
		merges.push_back(vec);
	}

	compute_enabled_and_merges_on_lits(clauses, enabled, merges);



	ofstream outFile;
	outFile.open(out_file);

	long num_merges = 0;
	long num_resolutions = 0;

	outFile.close();
	return 0;
}
