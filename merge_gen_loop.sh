#!/bin/bash

cnf=$1
cmty=$2

cp $cnf temp0.cnf

count=0
old_num_merges="None"
echo $old_num_merges


./clause_merge_analyzer temp0.cnf $cmty merge.r0 > analyzer.r0
curr_num_merges=`cat merge.r0 | grep NumMerges`
echo "Iteration 0"
echo $curr_num_merges

while [ "$old_num_merges" != "$curr_num_merges" ]
do
    #echo "IN"
    new_count=$(( count + 1 ))
    echo "Iteration ${new_count}"
    ./merge_generator temp${count}.cnf $cmty temp${new_count}.cnf > gen_out.r${new_count}
    ./clause_merge_analyzer temp${new_count}.cnf ${cmty} merge.r${new_count} > analyzer.r${new_count}
    count=$new_count
    old_num_merges=$curr_num_merges
    curr_num_merges=`cat merge.r${new_count} | grep NumMerges`
    echo $curr_num_merges

done


