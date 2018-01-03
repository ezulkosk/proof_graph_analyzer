/*
Graph Features Computation for SAT instances.
Version 2.2
Authors:
  - Carlos Ansótegui (DIEI - UdL)
  - María Luisa Bonet (LSI - UPC)
  - Jesús Giráldez-Cru (IIIA-CSIC)
  - Jordi Levy (IIIA-CSIC)

Contact: jgiraldez@iiia.csic.es

    Copyright (C) 2014  C. Ansótegui, M.L. Bonet, J. Giráldez-Cru, J. Levy

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <vector>
#include <assert.h>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <set>

#ifndef Graph_H
#define Graph_H

using namespace std;

struct classcomp {
  bool operator() (const pair<int,double>& a, const pair<int,double>& b) const
  {return a.first<b.first;}
};

class Graph {
private:
	int nnodes;                        // Number of nodes
	int typeA;
	double tarity;                      // Sum of the arities
	vector <double> narity;             // Arity of each node
	vector<set<pair<int,double>,classcomp > > neigh; // For each node, a list of its neighbors

public:

	Graph() : tarity(0), narity(0), nnodes(0), neigh(0) {}
	Graph(int n, int m) {
		tarity = 0;
		typeA = n;
		narity.resize(n+m, 0); 
		neigh.resize(n+m); 
		nnodes = n+m;
	}

	~Graph(){
		narity.clear();
		neigh.clear();	
	}

	double arity(int x) {     
		assert(x>=0 && x<= nnodes-1); 
		return narity[x]; 
	}
	
	int nNeighs(int x){
		return neigh[x].size();
	}
	
	int getTypeA(){
		return typeA;
	}

	int size() { return nnodes; }

	double arity() { return tarity; }

	void add_edge(int x, int y) { add_edge(x,y,1); };

	void add_edge(int x, int y, double w) {
		assert(x>=0 && x<= nnodes-1); 
		assert(y>=0 && y<= nnodes-1);
		//std::cout<<x<<" " << y<<std::endl;
		narity[x] += w;
		narity[y] += w;
		tarity += 2 * w;

		set<pair<int,double> >::iterator it=neigh[x].find(pair<int,double>(y,0));
		if(y>=typeA || it==neigh[x].end()){
			neigh[x].insert(pair<int,double>(y,w));
			if(x!=y) neigh[y].insert(pair<int,double>(x,w));
		}else{
			pair<int,double> p;
			p=make_pair(y,(it->second)+w);
			neigh[x].erase(it);
			neigh[x].insert(p);
			if(x!=y){
				it=neigh[y].find(pair<int,double>(x,0));
				p=make_pair(x,(it->second)+w);
				neigh[y].erase(it);
				neigh[y].insert(p);
			}
		}
	}

	double connected(int x, int y) {
		assert(x>=0 && x<= nnodes-1); 
		assert(y>=0 && y<= nnodes-1);

		set<pair<int,double> >::iterator it=neigh[x].find(pair<int,double>(y,0));
		if(it == neigh[x].end()) return (0);
		else return it->second;
	}

	void print() {
		//cerr <<"------ GRAPH ------\n";
		//cerr <<"NNodes = "<<nnodes<<endl;
		//cerr <<"TArity = "<<tarity<<endl;
		for(int i=0; i<nnodes; i++){
			cerr << i << " --> ";
			for(set<pair<int,double> >::iterator it=neigh[i].begin(); it!=neigh[i].end(); ++it)
				//if (i <= it->first) 
					cerr << it->first << " (" << it->second << ") ";
			cerr << endl;
		}
	}

//--------------- ITERATOR ON EDGES --------------------------------------------------

	typedef struct {int orig, dest; double weight;} edge;

	class EdgeIter;                          // Iterator to traverse all EDGES of the Graph
	friend class EdgeIter;

	class EdgeIter : public iterator<input_iterator_tag, edge> {
		Graph &g;
		set <pair <int,double> >::iterator it;
		int node;
		Graph::edge e;

	public:

		EdgeIter(Graph &x) : g(x){}
		EdgeIter(Graph &x, set<pair<int,double> >::iterator y, int n): g(x), it(y), node(n) { };

		EdgeIter (const EdgeIter &x) : g(x.g), it(x.it), node(x.node) {}

		EdgeIter& operator++() {
			assert(node >= 0 && node < g.nnodes);
			assert(node < g.nnodes-1 || it != g.neigh[node].end()); 
			do {
				if (++it == g.neigh[node].end() && node < g.nnodes-1) {
					do {
						node++; 
						it = g.neigh[node].begin(); 
					}
					while (it == g.neigh[node].end() && node < g.nnodes-1); //skip nodes without neighs
				}
			}
			while (!(it == g.neigh[node].end() && node == g.nnodes-1)   //skip when orig > dest
				&& node > it->first);
			return *this;
		}
		EdgeIter& operator++(int) {
			EdgeIter clone(*this);
			assert(node >= 0 && node < g.nnodes);
			assert(node < g.nnodes-1 || it != g.neigh[node].end()); 
			do {
				if (++it == g.neigh[node].end() && node < g.nnodes-1) {
					do {
						node++; 
						it = g.neigh[node].begin(); 
					}
					while (it == g.neigh[node].end() && node < g.nnodes-1); //skip nodes without neighs
				}
			}
			while (!(it == g.neigh[node].end() && node == g.nnodes-1)   //skip when orig > dest
				&& node > it->first);
			return *this;
		}

		bool operator==(const EdgeIter &rhs) {return it==rhs.it;}
		bool operator!=(const EdgeIter &rhs) {return it!=rhs.it;}

		Graph::edge& operator*() {
			cerr<<"ENTER *\n";
			e.orig = node; 
			e.dest = it->first; 
			e.weight = it->second; 
			return e;
		}
		Graph::edge *operator->() {
			e.orig = node; 
			e.dest = it->first; 
			e.weight = it->second; 
			return &e;
		}
	};

	EdgeIter begin() {
		int node = 0;
		while (neigh[node].size()==0 && node < nnodes) node++;
		return EdgeIter(*this, neigh[node].begin(), node);  
	// First neight of first node
	}

	EdgeIter end() {
		return EdgeIter(*this, neigh[nnodes-1].end(), nnodes);  
	// Last neight of last node
	}

//--------------- ITERATOR ON NEIGHBORS --------------------------------------------------

	class NeighIter : public set<pair<int,double> >::iterator {
		set<pair <int,double> >::iterator it;
		Graph::edge e;
	public:

		NeighIter(set<pair<int,double> >::iterator x) : it(x){}
		Graph::edge *operator->() {
			e.dest = it->first; 
			e.weight = it->second; 
			return &e;
		}
		NeighIter& operator++() {
			++it;
			return *this;
		}
		NeighIter& operator++(int) {
			it++;
			return *this;
		}
		bool operator==(const NeighIter& x) {return it==x.it;}
		bool operator!=(const NeighIter& x) {return it!=x.it;}
	};

	NeighIter begin(int x) { 
		assert(x>=0 && x<= nnodes-1); 
		return (NeighIter)neigh[x].begin(); 
	}

	NeighIter end(int x) { 
		assert(x>=0 && x<= nnodes-1); 
		return (NeighIter)neigh[x].end(); 
	}


};

#endif
