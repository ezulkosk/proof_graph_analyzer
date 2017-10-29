#!/bin/bash

#maplesat="/home/ezulkosk/git/maplesat/maplesat/simp/maplesat_static"
#graph_features="/home/ezulkosk/git/proof_graph_analyzer/features_s"

base_dir=$1
maplesat=$2
proof_graph=$3
graph_features=$4
cnf_file=$5

for i in graph cmty graph_cmty graph_cnf graph_scale_free
do
	mkdir -p ${base_dir}/${i}
done	

name=`basename $cnf_file .cnf`
echo $name

# get the cmty structure data for the cnf
#${graph_features} -5 -q ${base_dir}/cmty/${name}.cmty -y ${base_dir}/cmty/${name}.q ${cnf_file}

# get the proof graph certificate from maplesat
#${maplesat} -proof-graph=${base_dir}/graph/${name}.graph ${cnf_file}

# dump the graph_cnf
${proof_graph} ${base_dir}/graph/${name}.graph ${base_dir}/graph_cnf/${name}.graph_cnf

# get the cmty structure and power-law data for the graph_cnf
${graph_features} -1 -5 -k ${base_dir}/graph_scale_free/${name}.normalized_var_dist -g ${base_dir}/graph_scale_free/${name}.var_dist_plot -t ${base_dir}/graph_scale_free/${name}.var_dist -l ${base_dir}/graph_scale_free/${name}.graph_scale_free -q ${base_dir}/graph_cmty/${name}.graph_cmty -y ${base_dir}/graph_cmty/${name}.graph_q ${base_dir}/graph_cnf/${name}.graph_cnf


# graph clause-size|lbd vs utility 

# consider time
