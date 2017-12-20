all:
	g++ -o proof_graph_analyzer -std=c++11 src/proof_graph_analyzer.cpp
	g++ -g -o clause_merge_analyzer -std=c++11 src/merge_resolution_dimacs_check.cpp
