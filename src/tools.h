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

#include <math.h>
#include <vector>
#include <map>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#ifndef VECTOR
#include "graph_set.h"
#else
#include "graph_vector.h"
#endif

#ifndef TOOLS_H
#define TOOLS_H

using namespace std;

extern bool verbose = true; 


//-----------------------------------------------------------------------------
/*double abs(double x) {
//-----------------------------------------------------------------------------
  if (x < 0) return -x;
  else return x;
}
*/


//-----------------------------------------------------------------------------
double powlawc(int x, int xmin, double alpha) {
//-----------------------------------------------------------------------------
// Computes sum_{i=x}^{\infty} x^{alpha} / sum_{i=xmin}^{\infty} x^{alpha} 
// or approximates it as (x/xmin)^(alpha+1)
//-----------------------------------------------------------------------------
  assert(alpha < -1);
  assert(xmin <= x);
#define MAXITER 10000
  double num = 0, den = 0;
  int i;
  
  if (xmin < 25) {
    
    for (i=xmin; i<x; i++)
      den += pow(i, alpha);
    
    double pold = -2, p = -1;
    int n = 0;
    
    while (abs(p - pold) > 0.00000001 && n < MAXITER) {
      den += pow((double)i, alpha);
      num += pow((double)i, alpha);
      i++;
      n++;
      pold = p;
      p = num/den;
    }
    if (n < MAXITER)
      return p;
  }
  return pow((double)x/xmin, alpha + 1);
}


//-----------------------------------------------------------------------------
double exponc(int x, int xmin, double beta) {
//-----------------------------------------------------------------------------
  return exp(beta*(xmin - x)) ;
}


//-----------------------------------------------------------------------------
double mostlikely(vector <pair <int,int> > v, int maxxmin, char* fileout, char *nint, char *nplot, bool var) {
//-----------------------------------------------------------------------------

  //---- Compute vectors x, y, sxy, sylogx ------------------

	int n=v.size();
	vector <double> x(n), y(n+1), syx(n+1), sylogx(n+1);

	double Sy = 0;
	for (int i=0; i<n; i++) Sy += v[i].second;

	syx[n] = 0;
	sylogx[n] = 0;
	y[n] = 0;
	for(int i=n-1; i>=0; i--) {
		x[i]      = v[i].first;
		y[i]      = y[i+1] + v[i].second / Sy;
		sylogx[i] = sylogx[i+1] + v[i].second / Sy * log(x[i]);
		syx[i]    = syx[i+1] + v[i].second / Sy *     x[i];
    }

  //------ Compute, for powerlaw (a) and exponential (b), 
  //       the best alpha, xmin, dif and where is located--------------------

	double bestalpha, bestbeta;
	int bestxmina=0, bestxminb=0;
	double bestdifa = 1, bestdifb = 1;
	int bestinda, bestindb;
	int wherea, whereb;
  
	int ind = 0;
	int xmin;

	for (int ind=1; ind<=maxxmin && ind<n-3; ind++) {
		xmin = (int)x[ind];
    
		double alpha = -1 - 1 / (sylogx[ind] / y[ind] - log((double)xmin - 0.5));
		double beta = log(1 / (syx[ind] / y[ind] - xmin) + 1);

    //------------- Model powerlaw -----------------------------------------

		double worstdif = -1;
		int worstx = -1;
    
		for (int j=ind+1; j<n; j++) {
			double aux;
			aux = abs((double)y[j]/y[ind] - powlawc((int)x[j],xmin,alpha));
			if (aux >= bestdifa) {
				worstdif = aux;
				worstx = (int)x[j]; 
				j = n;  //Finish search of worst diff
			} else if (aux >= worstdif) {
				worstdif = aux;
				worstx = (int)x[j];
			}
		}
		for (int j=ind; j<n; j++) {
			if (x[j] + 1 < x[j+1]) {
				double aux;
				aux = abs((double)y[j+1]/y[ind] - powlawc((int)x[j]+1,xmin,alpha));
				if (aux >= bestdifa) {
					worstdif = aux;
					worstx = (int)x[j]+1; 
					j = n;  //Finish search of worst diff
				} else if (aux >= worstdif) {
					worstdif = aux;
					worstx = (int)x[j]+1;
				}
			}
		}
		if(worstdif < bestdifa) { 
			bestalpha = alpha; 
			bestxmina = xmin;
			bestdifa = worstdif;
			bestinda = ind;
			wherea = worstx;
		}

		if(fileout != NULL){	
		
		//------------- Model exponential -----------------------------------------
			worstdif = -1;
			worstx = -1;

			for (int j=ind+1; j<n; j++) {
				double aux;
				aux = abs((double)y[j]/y[ind] - exponc(x[j],xmin,beta));
				if (aux >= bestdifb) {
					worstdif = aux;
					worstx = x[j]; 
					j = n;  //Finish search of worst diff
				} else if (aux >= worstdif) {
					worstdif = aux;
					worstx = x[j];
				}
			}
			for (int j=ind; j<n; j++) {
				if (x[j] + 1 < x[j+1]) {
					double aux;
					aux = abs((double)y[j+1]/y[ind] - exponc(x[j]+1,xmin,beta));
					if (aux >= bestdifb) {
						worstdif = aux;
						worstx = x[j]+1; 
						j = n;  //Finish search of worst diff
					} else if (aux >= worstdif) {
						worstdif = aux;
						worstx = x[j]+1;
					}
				}
			}
			if(worstdif < bestdifb) { 
				bestbeta = beta;
				bestxminb = xmin;
				bestdifb = worstdif;
				bestindb = ind;
				whereb = worstx;
			}
		}
	}
	
	if(verbose){
	  cerr << "POWERLAW:\n";
	  cerr << "  alpha = "<<-bestalpha<<endl;
	  cerr << "  min = "<<bestxmina<<endl;
	  cerr << "  error = "<<bestdifa<<" in "<<wherea<<endl;
	  cerr << "EXPONENTIAL:\n";
	  cerr << "  beta = "<<bestbeta<<endl;
	  cerr << "  min = "<<bestxminb<<endl;
	  cerr << "  error = "<<bestdifb<<" in "<<whereb<<endl;
	}
	
	if(fileout != NULL){
		//------------ Print results ----------------------------------------------
		FILE* file = fopen(fileout, "w");
		fprintf(file, "POWERLAW:\n");
		fprintf(file, "alpha = %f\n",-bestalpha);
		fprintf(file, "min = %i\n",bestxmina);
		fprintf(file, "error = %f in %i\n",bestdifa,wherea);
		fprintf(file, "EXPONENTIAL:\n");
		fprintf(file, "beta = %f\n",bestbeta);
		fprintf(file, "min = %i\n",bestxminb);
		fprintf(file, "error = %f in %i\n",bestdifb,whereb);
		fclose(file);
	}
	
	if( (nint!=NULL && nplot==NULL) || (nint==NULL && nplot!=NULL) ){
		cerr << "[WARNING]: Both scale-free files for normalized distribution and plot should be defined." << endl;
	}else if(!(nint==NULL && nplot==NULL)){
		FILE *fint = fopen(nint, "w");
		if(fint == NULL)
			cerr << "Unable to open file: " << nint << endl;
		FILE *fplot = fopen(nplot, "w");
		if(fplot == NULL)
			cerr << "Unable to open file: " << nplot << endl;
		if(fint != NULL && fplot != NULL){
		    for (int i=0; i<n; i++)
		      fprintf(fint,"%f %f\n", x[i], y[i]);
		    fclose(fint);
			
			//char *eps = new char[strlen(nplot)];
			//strcpy(eps, nplot, sizeof(nplot)-4);
			//strcat(eps, ".eps");
			
			string eps(nplot);
			eps = eps.substr(0,strlen(nplot)-4);
			eps.append(".eps");
		
			fprintf(fplot,"set logscale xy\nset term postscript eps enhanced color\nset size 0.7,0.7\nset yrange [0.00001:1.1]\nset output \"%s\"\n", eps.c_str());
			if(var){
				fprintf(fplot, "plot \"%s\" ti \"variables\" lt 1 pt 7",nint, y[0]);
			    fprintf(fplot, ",%lf * x**%lf lt 1 ti \"{/Symbol a}=%0.2f\"", (double)y[bestinda]/ pow(bestxmina,bestalpha+1), bestalpha+1, -bestalpha);
				fprintf(fplot, ",%lf*exp(-%lf*(x - %d)) lt 2 ti \"{/Symbol b}=%0.3f\"",  (double)y[bestindb], bestbeta, bestxminb, bestbeta);
			}else{
				fprintf(fplot, "plot \"%s\" ti \"clauses\" lt 3 pt 7", nint, y[0]);
				fprintf(fplot, ",%lf * x**%lf lt 3 ti \"{/Symbol a}=%0.2f\"", (double)y[bestinda]/ pow(bestxmina,bestalpha+1), bestalpha+1, -bestalpha); 
				fprintf(fplot, ",%lf*exp(-%lf*(x - %d)) lt 4 ti \"{/Symbol b}=%0.3f\"",  (double)y[bestindb], bestbeta, bestxminb, bestbeta);
			}
			fprintf(fplot, ",\"< echo '%f %f'\" ti \"\"\n", (double)x[bestinda], y[bestinda]);
			fprintf(fplot, "quit\n");
			fclose(fplot);	
		}
	}

return -bestalpha;
}

//-----------------------------------------------------------------------------
pair <double,double> regresion(vector <pair <double,double> > &v) {
//-----------------------------------------------------------------------------
//Given a vector of points, computes the alpha and beta of a regression
//-----------------------------------------------------------------------------
  double Sx = 0, Sy = 0, Sxx = 0, Syy = 0, Sxy = 0;
  for (vector<pair <double,double> >::iterator it=v.begin(); it != v.end(); it++) {
    double x = it->first;
    double y = it->second;
    Sx += x;
    Sy += y;
    Sxx += x * x;
    Syy += y * y;
    Sxy += x * y;
  }
  
  double alpha = (Sx * Sy - v.size() * Sxy)/( Sx * Sx - v.size() * Sxx);
  double beta = Sy / v.size() - alpha * Sx / v.size();
  return pair <double,double>(alpha,beta);
}

/*
//-------------------------------------------------------------------------------------------
void shuffle(vector <int> &x) {
//-------------------------------------------------------------------------------------------
// Randomly re-order elements of a vector
//-------------------------------------------------------------------------------------------
	for (int i=0 ; i<x.size()-1 ; i++) {
		int j = rand()%(x.size()-i)+i;
		int aux = x[i];
		x[i] = x[j];
		x[j] = aux;
	}
}
*/


//--------------------------------------------------------------------------------
pair<Graph*,Graph*> readFormula(char* filename, int MAXCLAUSE){
//--------------------------------------------------------------------------------
// Given a CNF file (filename) in DIMACS format, creates it correspondent formula
//  disregarding clauses of size greater than MAXCLAUSE
//--------------------------------------------------------------------------------	
	FILE *source;
	source = fopen(filename, "r");
	if(!source){
		cerr << "Unable to read CNF file " << filename << endl;
		exit(-1);
	}

	int totVars=0, totClauses=0;
	int var=0;
	int aux=-1;

	// Skip comments
	while((aux=getc(source))=='c'){
		while (getc(source)!='\n')
			;
	}
	ungetc(aux,source);

	// File Head
	if( !fscanf(source, "p cnf %i %i", &totVars, &totClauses)) {
		cerr << "Invalid CNF file\n";
		exit(-1);
	}

	Graph* vig = new Graph(totVars, 0);
	Graph* cvig = new Graph(totVars,totClauses);
	
	// Read the clauses
	vector<int> clause;
	int nclauses=0;

	while(fscanf(source, "%i", &var)==1) {
		if (var==0) {
			if (clause.size() <= MAXCLAUSE && clause.size()>0) {	
				double weight_vig = 0;
				double weight_cvig = 1.0/clause.size();
				if(clause.size()>1){
					weight_vig = 2.0 / (clause.size() * (clause.size()-1) );
					for (int i=0; i<clause.size()-1; i++){
						cvig->add_edge(clause[i], totVars+nclauses, weight_cvig);
						for (int j=i+1; j<clause.size(); j++){
							vig->add_edge(clause[i], clause[j], weight_vig);
						}
					}
					cvig->add_edge(clause[clause.size()-1], totVars+nclauses, weight_cvig);
				}else if(clause.size()==1){
					cvig->add_edge(clause[0], totVars+nclauses, weight_cvig);	
				}
			} else {
				if(verbose)
					cerr << "\tDisregarded clause of size " << clause.size() << endl;
			}
			nclauses++;
			clause.clear();
		} else {
			if (abs(var) > totVars) {
				cerr << "Unvalid variable number " << abs(var) << endl;
				exit(-1);
			}
			clause.push_back(abs(var)-1);
		}
	}
	
	fclose(source);
	
	return make_pair(vig,cvig);
}

//--------------------------------------------------------------------------------
Graph* readVIG(char* filename, int MAXCLAUSE){
//--------------------------------------------------------------------------------
// Given a CNF file (filename) in DIMACS format, creates it correspondent formula
//  disregarding clauses of size greater than MAXCLAUSE
//--------------------------------------------------------------------------------	
	FILE *source;
	source = fopen(filename, "r");
	if(!source){
		cerr << "Unable to read CNF file " << filename << endl;
		exit(-1);
	}

	int totVars=0, totClauses=0;
	int var=0;
	int aux=-1;

	// Skip comments
	while((aux=getc(source))=='c'){
		while (getc(source)!='\n')
			;
	}
	ungetc(aux,source);

	// File Head
	if( !fscanf(source, "p cnf %i %i", &totVars, &totClauses)) {
		cerr << "Invalid CNF file\n";
		exit(-1);
	}

	Graph* vig = new Graph(totVars, 0);
	
	// Read the clauses
	vector<int> clause;

	while(fscanf(source, "%i", &var)==1) {
		if (var==0) {
			if (clause.size() <= MAXCLAUSE && clause.size()>1) {	
				double weight_vig = 2.0 / (clause.size() * (clause.size()-1) );
				for (int i=0; i<clause.size()-1; i++){
					for (int j=i+1; j<clause.size(); j++){
					  //std::cout<<"ASDF "<<clause[i]<<std::endl;
						vig->add_edge(clause[i], clause[j], weight_vig);
					}
				}
			} else {
				if(verbose && clause.size()>1)
					cerr << "\tDisregarded clause of size " << clause.size() << endl;
			}
			clause.clear();
		} else {
			if (abs(var) > totVars) {
				cerr << "Unvalid variable number " << abs(var) << endl;
				exit(-1);
			}
			clause.push_back(abs(var)-1);
		}
	}
	
	fclose(source);
	
	return vig;
}

//--------------------------------------------------------------------------------
Graph* readCVIG(char* filename, int MAXCLAUSE){
//--------------------------------------------------------------------------------
// Given a CNF file (filename) in DIMACS format, creates it correspondent formula
//  disregarding clauses of size greater than MAXCLAUSE
//--------------------------------------------------------------------------------	
	FILE *source;
	source = fopen(filename, "r");
	if(!source){
		cerr << "Unable to read CNF file " << filename << endl;
		exit(-1);
	}

	int totVars=0, totClauses=0;
	int var=0;
	int aux=-1;

	// Skip comments
	while((aux=getc(source))=='c'){
		while (getc(source)!='\n')
			;
	}
	ungetc(aux,source);

	// File Head
	if( !fscanf(source, "p cnf %i %i", &totVars, &totClauses)) {
		cerr << "Invalid CNF file\n";
		exit(-1);
	}

	Graph* cvig = new Graph(totVars,totClauses);
	
	// Read the clauses
	vector<int> clause;
	int nclauses=0;

	while(fscanf(source, "%i", &var)==1) {
		if (var==0) {
			if (clause.size() <= MAXCLAUSE && clause.size()>0) {	
				double weight_cvig = 1.0/clause.size();
				for (int i=0; i<clause.size(); i++){
					cvig->add_edge(clause[i], totVars+nclauses, weight_cvig);
				}
			} else {
				if(verbose)
					cerr << "\tDisregarded clause of size " << clause.size() << endl;
			}
			nclauses++;
			clause.clear();
		} else {
			if (abs(var) > totVars) {
				cerr << "Unvalid variable number " << abs(var) << endl;
				exit(-1);
			}
			clause.push_back(abs(var)-1);
		}
	}
	
	fclose(source);
	
	return cvig;
}



bool mycompare(pair<double,double> x, pair<double,double> y) {
  return x.first <= y.first;
}

//-----------------------------------------------------------------------------
vector<pair<int,int> > readPoints(char *filein) {
//-----------------------------------------------------------------------------
// Read a list of points. When two points share the same first component
// as (x, y1) (x, y2), interpret it as (x, y1+y2)
//-----------------------------------------------------------------------------
	FILE *fin = fopen(filein, "r");
	if(fin == NULL){
		cerr << "Unable to open file: " << filein << endl;
	}
	  map<int,int> m;
	  vector<pair<int,int> > v;
	  int x,y;

	  int read = fscanf(fin, "%d %d\n", &x, &y);
	  while (read == 2) {
	    if (m.find(x) != m.end()) m[x] += y;
	    else m[x] = y;
	    read = fscanf(fin, "%d %d\n", &x, &y);
	  }

	  for (map<int,int>::iterator it=m.begin(); it!=m.end(); it++) 
	    v.push_back(pair<int,int>(it->first, it->second));

	  sort(v.begin(), v.end(), mycompare);
	  
	fclose(fin);

	  return v;
}



#endif
