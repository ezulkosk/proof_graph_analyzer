//============================================================================
// Name        : proof_graph_analyzer.cpp
// Author      : 
// Version     :
//============================================================================


/*
 * Takes a cnf and greedily increases merges in the formula.
 *
 * Defintions:
 *   enabler literal - for a clause pair, resolution over the literal allows a merge to happen over some other literal
 *   merge literal - there exists a clause pair that has a merge over the literal
 *   double resolution - a clause pair resolves on two literals, and merges on a third: (a v b v c) ^ (!a v !b v c)
 *
 * The program takes an input CNF. The clause "skeleton" is "locked", i.e., the set of clauses remain the same,
 * but the polarities of literals can be changed. The number of each literal is also locked, so in order to
 * swap v to !v, !v must be swapped to v elsewhere. This ensures that the community structure and literal popularity
 * distribution of the original instance is maintained. It also ensures that if 2 clauses are already merge+res, they
 * will still be.
 *
 * If a literal is an enabler or a merge, it is also locked. If there is a double resolution, one may flip.
 *
 */

#include <ctime>
#include <assert.h>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>      // std::stringstream

#include <set>
#include <utility>
#include <vector>


using namespace std;



void dump_clauses(vector< vector<int> >& clauses, int& nVars, const char* out_file){
	ofstream out;
	out.open(out_file);

	out<<"p cnf "<<nVars<<" "<<clauses.size()<<endl;
	for(auto c : clauses){
		for(int l: c){
			out<<l<<" ";
		}
		out<<"0"<<endl;
	}
	out.close();
}

bool my_comparator (int i, int j) {
	return abs(i) < abs(j) || (abs(i) == abs(j) && i < j);
}


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

void growTo(vector<bool>& v, int s){
	if(v.size() < s){
		for(int i = v.size(); i < s; i++){
			v.push_back(false);
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


void read_dimacs(char* dimacs_file, vector< vector<int> >& clauses, int& nVars){
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
	stream >> nVars;
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
 * Increases enabled and merge counters for appropriate literals
 * Condition for positive change: the clause has no merges but has a double-res with another clause.
 * 		   If we flip one literal, it will have a res+merge with that clause.
 */
bool res_merge_counter(vector<int>& c1, vector<int>& c2,
		vector<int>* e1, vector<int>* e2,
		vector<int>* m1, vector<int>* m2,
		vector<bool>* lock1, vector<bool>* lock2,
		vector< pair<int,int> >& curr_flip_pairs){

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
	bool updated = false;
	if(res_indices1.size() > 0 && num_merges > 0){
		// print_vector(*m1);
		for(int i = 0; i < res_indices1.size(); i++){
			e1->at(res_indices1[i]) += num_merges;
		}

		for(int i = 0; i < res_indices2.size(); i++){
			e2->at(res_indices2[i]) += num_merges;
		}

		if(res_indices1.size() == 1){
			lock1->at(res_indices1[0]) = true;
			lock2->at(res_indices2[0]) = true;
		}

		for(int i = 0; i < merge_indices1.size(); i++){
			m1->at(merge_indices1[i])++;
			lock1->at(merge_indices1[i]) = true;
		}

		for(int i = 0; i < merge_indices2.size(); i++){
			m2->at(merge_indices2[i])++;
			lock2->at(merge_indices2[i]) = true;

		}
		// print_vector(*m1);
		updated = true;
	}
	if(res_indices1.size() > 1){
		// create flip pairs
		for(int i = 0; i < res_indices1.size(); i++){
			pair<int, int> p = std::make_pair(res_indices1[i], res_indices2[i]);
			curr_flip_pairs.push_back(p);
		}
		updated = true;
	}
	return updated;
}




void compute_enabled_and_merges_on_lits(vector< vector<int> >& clauses,
		vector< vector<int> >& enabled,
		vector< vector<int> >& merges,
		vector< vector<bool> >& locked,
		vector< pair< pair<int,int>, pair<int,int> > >& flip_pairs,
		int& num_flips,
		int num_flips_per_dump,
		int max_flips,
		int& nVars,
		char* out_file){

	for(int i = 0; i < (int) clauses.size() - 1; i++){
		vector<int> c1 = clauses[i];
		vector<int> e1 = enabled[i];
		vector<int> m1 = merges[i];
		vector<bool> lock1 = locked[i];
		for(int j = i + 1; j < (int) clauses.size(); j++){
			vector<int> c2 = clauses[j];
			vector<int> e2 = enabled[j];
			vector<int> m2 = merges[j];
			vector<bool> lock2 = locked[j];
			vector< pair<int, int> > curr_flip_pairs; // does not contain clause_ids. Incorporate those after.
			bool updated = res_merge_counter(c1, c2, &e1, &e2, &m1, &m2, &lock1, &lock2, curr_flip_pairs);


			if(updated){
				// print_vector(m1);
				merges[i] = m1;
				merges[j] = m2;
				enabled[i] = e1;
				enabled[j] = e2;
				locked[i] = lock1;
				locked[j] = lock2;
				// cout<<"===="<<endl;
				for(auto p: curr_flip_pairs){
					pair<int, int> p1 = std::make_pair(i, p.first);
					pair<int, int> p2 = std::make_pair(j, p.second);
					pair< pair<int,int>, pair<int,int> > p3 = std::make_pair(p1, p2);
					flip_pairs.push_back(p3);
				}
			}
		}
	}
	cout<<endl;

	map<int, vector< pair<int, int> > > unlocked_lit_locs;
	for(int i = 0; i < clauses.size(); i++){
		vector<int> c = clauses[i];
		vector<bool> lock = locked[i];

		for(int j = 0; j < c.size(); j++){
			if(!lock[j]){
				// get the current set associated with c[j]
				vector< pair<int, int> >* v;
				if(unlocked_lit_locs.find(c[j]) != unlocked_lit_locs.end()){
					v = &unlocked_lit_locs[c[j]];
				}
				else{
					v = new vector< pair<int, int> >;
				}
				v->push_back(std::make_pair(i, j));
				unlocked_lit_locs[c[j]] = *v;
			}
		}
	}

	set<int> locked_clauses;

	// randomize
	std::random_shuffle (flip_pairs.begin(), flip_pairs.end());

	for(auto p: flip_pairs){
		int c1 = p.first.first;
		int l1 = p.first.second;
		int lit1 = clauses[c1][l1];
		int c2 = p.second.first;
		int l2 = p.second.second;
		int lit2 = clauses[c2][l2];

		vector<int> clause1 = clauses[c1];
		vector<int> clause2 = clauses[c2];

		// either of the clauses in the pair have been changed
		if(locked_clauses.find(c1) != locked_clauses.end() || locked_clauses.find(c2) != locked_clauses.end())
			continue;

		bool actual_flip = false;
		if(!locked[c1][l1] && unlocked_lit_locs.find(-clauses[c1][l1]) != unlocked_lit_locs.end()){
			// swap for c1[l1]
			vector< pair<int,int> >* vec = &(unlocked_lit_locs[-clauses[c1][l1]]);
			for(auto pp : (*vec)){
				if(locked_clauses.find(pp.first) == locked_clauses.end() && pp.first != c1 && pp.first != c2){
					//cout<<"======"<<endl;
					print_vector(clauses[c1]);
					print_vector(clauses[c2]);
					print_vector(clauses[pp.first]);
					clauses[c1][l1] *= -1;
					clauses[pp.first][pp.second] *= -1;

					//print_vector(clauses[c1]);
					//print_vector(clauses[c2]);
					//cout<<"======"<<endl;

					// lock the two involved clauses, and the swapper clause
					locked_clauses.insert(c1);
					locked_clauses.insert(c2);
					locked_clauses.insert(pp.first);
					num_flips++;
					if(num_flips % num_flips_per_dump == 0){
						stringstream ss;
						ss << out_file << num_flips;
						cout<<"Dump: "<<ss.str()<<endl;
						dump_clauses(clauses, nVars, ss.str().c_str());
					}
					actual_flip = true;
					break;
				}
			}
		}
		else if(!actual_flip && !locked[c2][l2] && unlocked_lit_locs.find(-clauses[c2][l2]) != unlocked_lit_locs.end()){

			// swap for c2[l2]
			vector< pair<int,int> >* vec = &(unlocked_lit_locs[-clauses[c2][l2]]);
			for(auto pp : (*vec)){
				if(locked_clauses.find(pp.first) == locked_clauses.end() && pp.first != c1 && pp.first != c2){

					//cout<<"======"<<endl;
					//print_vector(clauses[c1]);
					//print_vector(clauses[c2]);

					clauses[c2][l2] *= -1;
					clauses[pp.first][pp.second] *= -1;

					//print_vector(clauses[c1]);
					//print_vector(clauses[c2]);
					//cout<<"======"<<endl;

					// lock the two involved clauses, and the swapper clause
					locked_clauses.insert(c1);
					locked_clauses.insert(c2);
					locked_clauses.insert(pp.first);
					num_flips++;
					if(num_flips % num_flips_per_dump == 0){
						stringstream ss;
						ss << out_file << num_flips;
						cout<<"Dump: "<<ss.str()<<endl;
						dump_clauses(clauses, nVars, ss.str().c_str());
					}
					break;
				}
			}
		}
		if(num_flips >= max_flips)
			break;

	}
}


int main(int argc, char * argv[]) {

	std::srand ( unsigned ( std::time(0) ) );

	vector< vector<int> > clauses;
	vector<int> cmty;
	int nVars = 0;

	if(argc < 4){
		cout << "USAGE: ./merge_generator cnf_file cmty_file out_cnf_file" << endl;
		return 1;
	}

	char* dimacs_file = argv[1];
	//char* cmty_file = argv[2];
	char* out_file = argv[2];
	int num_flips_per_dump = atoi(argv[3]);
	int max_flips = atoi(argv[4]);


	read_dimacs(dimacs_file, clauses, nVars);
	//read_cmtys(cmty_file, cmty);

	int old_num_flips = -1;
	int num_flips = 0;
	while(old_num_flips != num_flips && num_flips < max_flips){

		vector< vector<int> > enabled; // vector to store the number of enabled merges of every individual literal
		vector< vector<int> > merges; // vector to store the number of actual merges on every individual literal
		vector< vector<bool> > locked; // if a literal is locked (single enabler/merge), it cannot be flipped
		// stores ((clause_id, lit_id), (clause_id, lit_id)) pairs which can be flipped,
		// but must ensure that the 2 involved lits are not locked.
		vector< pair< pair<int, int>, pair<int, int> > > flip_pairs;

		for(auto c : clauses){
			vector<int> vec = *(new vector<int>);
			vector<int> vec2 = *(new vector<int>);
			vector<bool> vec3 = *(new vector<bool>);
			growTo(vec, c.size());
			growTo(vec2, c.size());
			growTo(vec3, c.size());
			enabled.push_back(vec);
			merges.push_back(vec2);
			locked.push_back(vec3);
		}
		old_num_flips = num_flips;
		compute_enabled_and_merges_on_lits(clauses, enabled, merges, locked, flip_pairs, num_flips, num_flips_per_dump, max_flips, nVars, out_file);
		cout<<"Num Flips: "<<num_flips<<endl;
		//dump_clauses(clauses, nVars, out_file);
	}
	stringstream ss;
	ss << out_file << num_flips;
	cout<<"Dump: "<<ss.str();
	dump_clauses(clauses, nVars, ss.str().c_str());

	return 0;
}
