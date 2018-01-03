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
#include <stack>
#ifndef VECTOR
#include "graph_set.h"
#else
#include "graph_vector.h"
#endif
#include <algorithm>

#ifndef DIMENSION_H
#define DIMENSION_H

using namespace std;

typedef int diameter;

extern bool verbose;
extern int maxx;

//-------------------------------------------------------------------------------------------
int components(Graph *g) {
//-------------------------------------------------------------------------------------------
// Computes the number of disconected components of a graph
//-------------------------------------------------------------------------------------------
	stack <int> s;
	vector <bool> covered(g->size(), false);
	int comp = 0;

	for (int c=0; c<g->size(); c++) {
		if (!covered[c]) {
			covered[c] = true;
			s.push(c);
			while(!s.empty()) {
				int v = s.top();
				s.pop();
				for (Graph::NeighIter it=g->begin(v); it != g->end(v); it++)
				if (!covered[it->dest]) {
					covered[it->dest] = true;
					s.push(it->dest);
				}
			}
			comp++;
		}
	}
	return comp;
}

//-----------------------------------------------------------------------------
int tile(int c, diameter d, vector <diameter> &cover, Graph *g)
//-----------------------------------------------------------------------------
//Given a node c and a distance d, returns the number of nodes c2 reachable from c
//at distance d and not marked as cover[c2]==-1
//-----------------------------------------------------------------------------
{
	int ncover = 0;

	stack <pair <int, diameter> > s;

	s.push(make_pair(c,d));
// cerr <<"center "<<c<<" "<<d<<endl;
	while(!s.empty()) {
		int c = s.top().first;
		diameter d = s.top().second;
		s.pop();
	// cerr <<c<<" "<<d<<endl;

		if (d > cover[c]) {   //we've found a best path till node c
		if(cover[c] == -1) { //node c is not covered
			ncover++;
		}

		cover[c] = d;

	//push neighs of c to the stack with d-1
		for (Graph::NeighIter it=g->begin(c); it != g->end(c); it++)
			if (d > 1) 
			s.push(make_pair(it->dest, d-1));
	}
}
return ncover;
}

//-----------------------------------------------------------------------------
int needed(Graph *g, diameter d, vector <int> &centers) {  
//-----------------------------------------------------------------------------
// Given a graph g, and a diameter d, computes how many tiles of diameter d are
// needed for covering the graph, trying tile centers as centers
//-----------------------------------------------------------------------------


	if(d==1)
		return g->size();

	int ncover = 0;
	vector <diameter> cover(g->size(), -1);
	int needed = 0;       // needed[l] = 0;

	int i=0;
	while (ncover < g->size()) {
		int c = centers[i++];
		while (cover[c] != -1) 
			c = centers[i++];
		int tc = tile(c, d, cover, g);
	//cerr<<(double)i*100/g->size()<<" "<<(double)ncover*100/g->size()<<endl;
		if (tc > 0) {
			needed++;
			ncover += tc;
		}
	}
	return needed;
}


//-----------------------------------------------------------------------------
bool comparesecond(pair <int,int> a, pair <int,int> b) {
	return a.second > b.second;
}

//-----------------------------------------------------------------------------
vector <int> computeNeeded(Graph *g) {  
//-----------------------------------------------------------------------------
// Given a (weighted) graph g, computes needed[i] as the number of tiles of diameter i
// needed for covering the graph
//-----------------------------------------------------------------------------

	vector <int> v(1, g->size());     // v[0] = g->size();

	int comp = components(g);
	if(verbose)
		cerr << "\tComponents: " << comp << endl;

	vector <pair <int,int> > aux(g->size());
	for (int i=0; i<g->size(); i++) {
		aux[i].first = i;
		//aux[i].second = (int)g->arity(i);
		aux[i].second = (int)g->nNeighs(i);
	}
	sort(aux.begin(), aux.end(), comparesecond);

	vector <int> centers(g->size());
	for (int i=0; i<g->size(); i++) 
		centers[i] = aux[i].first;
//shuffle(centers);

	for (int d=1; d<=maxx && v[d-1]>comp; d++) {
		v.push_back(needed(g,(diameter)d,centers)); // v[d] = needed(g,d,centers);
		if(verbose)
			cerr << "\t" << d << " => " << v[d] <<endl;
	}
	return v;
}

#endif
