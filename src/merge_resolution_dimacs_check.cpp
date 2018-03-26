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





void print_vector(vector<int>& v, bool newline=true){
	for(auto i: v){
		cout<<i<<" ";
	}
	if(newline)
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


// checks if resolvable clauses are merges in a cmty-based fashion
// if 1/n of a clause is in a cmty, count that fraction toward the community, both for clause count and merges
void compute_num_merges(vector< vector<int> >& clauses,
		vector<int>& cmty,
		vector<double>& cmty_merges,
		vector<double>& cmty_resolutions,
		long& num_merges,
		long& num_resolutions){
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
				if(ci[ni] == cj2[nj] && !canResolve){
					canResolve = true;
					ni++;
					nj++;
				}
				else if(ci[ni] == cj2[nj] && !canResolve){
					canResolve = false;
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
			num_resolutions++;
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

}


void intra_community_merge_res(vector< vector<int> >& clauses,
		vector<int>& cmty,
		int num_cmtys,
		ofstream& out){

	vector< vector<int> > clauses_in_cmty; // stores the indices of clauses in each cmty
	vector<int> intra_cmty_res;
	vector<int> intra_cmty_merge;
	vector<int> intra_cmty_clauses;
	vector<int> clauses_outside_cmtys;
	long inter_cmty_res = 0;
	long inter_cmty_merge = 0;

	growTo(intra_cmty_res, num_cmtys);
	growTo(intra_cmty_merge, num_cmtys);
	growTo(intra_cmty_clauses, num_cmtys);

	for(int i = 0; i <= num_cmtys; i++){
		vector<int>* vec = new vector<int>;
		clauses_in_cmty.push_back(*(vec));
	}

	for(int i = 0; i < clauses.size(); i++){
		vector<int> c = clauses[i];
		int curr_cmty = -1;
		bool is_intra = true;
		for(auto l: c){
			int v = abs(l);
			int v_cmty = cmty[v];
			if(curr_cmty == -1)
				curr_cmty = v_cmty;
			else if(curr_cmty != v_cmty){
				is_intra = false;
				break;
			}
		}
		if(is_intra && curr_cmty >= 0){
			clauses_in_cmty[curr_cmty].push_back(i);
			intra_cmty_clauses[curr_cmty] += 1;
		}
		else{
			clauses_outside_cmtys.push_back(i);
		}
	}


	for(int i = 0; i < (int)clauses_in_cmty.size(); i++){
		for(int j = 0; j < (int)clauses_in_cmty[i].size()-1; j++){
			vector<int> c1 = clauses[clauses_in_cmty[i][j]];
			for(int k = j + 1; k < clauses_in_cmty[i].size(); k++){
				vector<int> c2 = clauses[clauses_in_cmty[i][k]];
				int num_merge;
				bool res;
				is_res_merge(c1, c2, res, num_merge);

				if(res){
					intra_cmty_res[i] += 1;
					intra_cmty_merge[i] += num_merge;
				}
			}
		}
	}

	for(int i = 0; i < (int)clauses_outside_cmtys.size()-1; i++){
		vector<int> c1 = clauses[clauses_outside_cmtys[i]];
		for(int j = i + 1; j < clauses_outside_cmtys.size(); j++){
			vector<int> c2 = clauses[clauses_outside_cmtys[j]];
			int num_merge;
			bool res;
			is_res_merge(c1, c2, res, num_merge);

			if(res){
				inter_cmty_res += 1;
				inter_cmty_merge += num_merge;
			}
		}
	}


	out<<"i clauses res merges\n";
	for(int i = 0; i < num_cmtys; i++){
		out<<"CM "<<i<<" "<<intra_cmty_clauses[i]<<" "<<intra_cmty_res[i]<<" "<<intra_cmty_merge[i]<<endl;
	}
	out<<"CMINTER "<<-1<<" "<<clauses_outside_cmtys.size()<<" "<<inter_cmty_res<<" "<<inter_cmty_merge<<endl;

}


void compute_num_merges2(vector< vector<int> >& clauses,
		vector<int>& cmty,
		vector< vector< vector<int>* > >& var_inclusion,
		vector<double>& cmty_merges,
		vector<double>& cmty_resolutions,
		vector<int>& var_merges,
		vector<int>& var_resolutions,
		long& num_merges,
		long& num_resolutions){

	for(int i = 2; i < var_inclusion.size(); i+=2){
		num_resolutions += var_inclusion[i].size() * var_inclusion[i+1].size();
		var_resolutions[i/2] += var_inclusion[i].size() * var_inclusion[i+1].size();
		for(int j = 0; j < var_inclusion[i].size(); j++){
			for(auto l: *(var_inclusion[i][j])){
				cmty_resolutions[cmty[abs(l)]] += var_inclusion[i+1].size() / ((double) (*var_inclusion[i][j]).size());
			}
		}
		for(int j = 0; j < var_inclusion[i+1].size(); j++){
			for(auto l: *(var_inclusion[i+1][j])){
				cmty_resolutions[cmty[abs(l)]] += var_inclusion[i].size() / ((double) (*var_inclusion[i+1][j]).size());
			}
		}

		int ni = 0;
		int nj = 0;
		for(int j = 0; j < var_inclusion[i].size(); j++){
			for(int k = 0; k < var_inclusion[i+1].size(); k++){
				vector<int> ci = *var_inclusion[i][j];
				vector<int> cj = *var_inclusion[i+1][k];
				ni = 0;
				nj = 0;
				while(ni < ci.size() && nj < cj.size()){
					if(ci[ni] == cj[nj]){
						var_merges[abs(ci[ni])]++;
						num_merges++;
						//cout<<"CCC ";
						//print_vector(ci, false);
						//cout<<" CCC ";
						//print_vector(cj);
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


	}

	cout<<num_resolutions<<endl;
	return;
}


double remap_highest(vector< vector<int> >& clauses, int highest){
	ofstream out;
	out.open("/home/ezulkosk/git/proof_graph_analyzer/remap.cnf");
	for(auto c: clauses){
		for(auto l : c){
			if(abs(l) == highest){
				l = abs(l);
			}
			out<<l<<" ";
		}
		out<<"0\n";
	}
	out.close();
}



int main(int argc, char * argv[]) {
	vector< vector<int> > clauses;
	vector<int> cmty;
	vector<int> cmty_picks;
	vector<int> cmty_size;
	vector<int> var_merges;
	vector<int> var_resolutions;
	vector<double> cmty_clauses;
	vector<double> cmty_merges;
	vector<double> cmty_resolutions;
	vector< vector< vector<int>* > > var_inclusion; // for each literal, include a pointer to all clauses containing it


	if(argc < 4){
		cout << "USAGE: ./proof_graph_analyzer cnf_file cmty_file out_file" << endl;
		return 1;
	}
	char* dimacs_file = argv[1];
	char* cmty_file = argv[2];
	char* out_file = argv[3];


	read_dimacs(dimacs_file, clauses);
	read_cmtys(cmty_file, cmty, cmty_picks, cmty_size, cmty_clauses);

	// add to cmty_clauses
	for(auto c: clauses){
		for(auto l: c)
		{
			cmty_clauses[cmty[abs(l)]] += 1 / ((double) c.size());
		}
	}

	// create maps from each literal to the clauses that contain that literal
	growTo(var_inclusion, (cmty.size()+2)*2);
	growTo(var_resolutions, cmty.size());
	growTo(var_merges, cmty.size());
	for(int i = 0; i < clauses.size(); i++){
		vector<int>* c = &(clauses[i]);
		for(auto l: *c){
			//cout<<toLit(l)<<endl;
			assert(toLit(l) < var_inclusion.size());
			//cout<<var_inclusion[toLit(l)].size()<< " "<<var_inclusion[toLit(-l)].size()<<endl;

			var_inclusion[toLit(l)].push_back(c);
			//print_vector(*var_inclusion[toLit(l)][0]);
			//print_vector(*var_inclusion[toLit(l)][0]);
			//print_vector(*var_inclusion[toLit(l)][0]);


		}
	}

	//assert(toLit(431) == 862);
	//print_vector(*var_inclusion[toLit(431)][0]);
	//exit(1);

	//exit(1);
	/*
	for(int i = 2; i < var_inclusion.size(); i+=2){
		cout<<"A "<<var_inclusion[i].size()<< " "<<var_inclusion[i+1].size()<<endl;
		if(var_inclusion[i].size() > 0)
			print_vector(*var_inclusion[i][0]);
		if(var_inclusion[i+1].size() > 0)
			print_vector(*var_inclusion[i+1][0]);

		for(int j = 0; j < var_inclusion[i].size(); j++){
			for(int k = 0; k < var_inclusion[i+1].size(); k++){
				//cout<<"-------\n";
				//cout<<var_inclusion[i].size()<< " "<<var_inclusion[i+1].size()<<endl;
				//print_vector(*var_inclusion[i][j]);
				//print_vector(*var_inclusion[i+1][k]);
			}
		}
	}

	exit(1);
	*/




	growTo(cmty_merges, cmty_clauses.size());
	growTo(cmty_resolutions, cmty_clauses.size());

	ofstream outFile;
	outFile.open(out_file);

	long num_merges = 0;
	long num_resolutions = 0;

	compute_num_merges(clauses, cmty, cmty_merges, cmty_resolutions, num_merges, num_resolutions);
	//compute_num_merges2(clauses, cmty, var_inclusion, cmty_merges, cmty_resolutions, var_merges, var_resolutions, num_merges, num_resolutions);

	intra_community_merge_res(clauses, cmty, cmty_clauses.size(), outFile);

	int highest = 0;
	int high_var = 1;
	for(int i = 1; i < var_resolutions.size(); i++){
		if(abs(int(var_inclusion[2*i].size()) - int(var_inclusion[2*i+1].size())) > highest)
		{
			highest = abs(int(var_inclusion[2*i].size()) - int(var_inclusion[2*i+1].size()));
			high_var = i;
		}
		// cout<<"V "<<i<<" "<<abs(int(var_inclusion[2*i].size()) - int(var_inclusion[2*i+1].size()))<<" "<<var_resolutions[i]<<" "<<var_merges[i]<<endl;
	}
	cout<<"Highest "<<high_var<<" "<<highest<<endl;

	remap_highest(clauses, 16); // high_var);


	// discuss instance of all positive literals ("counterexample of mergability")
	// sc2 initiative

	printf("cmty clauses resolutions merges\n");
	for(int i = 0; i < cmty_clauses.size(); i++){
		printf("%d %8.4f %8.4f %8.4f\n", i, cmty_clauses[i], cmty_resolutions[i], cmty_merges[i]);
	}

	double avg_cmty_merge_resolutions = 0;
	for(int i = 0; i < cmty_clauses.size(); i++){
		avg_cmty_merge_resolutions += cmty_merges[i] / cmty_resolutions[i];
	}
	avg_cmty_merge_resolutions /= cmty_clauses.size();


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
	cout<<dimacs_file<<","<<nVars<<","<<clauses.size()<<","<<clause_pairs<<","<<num_merges<<","<<(double)num_merges/clause_pairs<<endl;
	outFile<<"ClausePairs,"<<clause_pairs<<endl;
	outFile<<"NumResolutions,"<<num_resolutions<<endl;
	outFile<<"NumMerges,"<<num_merges<<endl;
	outFile<<"AvgCmtyMergeOverResolutions,"<<avg_cmty_merge_resolutions<<endl;
	for(int i = 0; i < cmty_clauses.size(); i++){
		outFile<<"c,"<<i<<","<<cmty_clauses[i]<<","<<cmty_resolutions[i]<<","<<cmty_merges[i]<<endl;
	}
	return 0;
}
