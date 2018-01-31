#!/bin/bash

cnf=$1
cmty=$2

name=`basename $cnf .cnf`

prefix=$name

cp $cnf ${prefix}_temp0.cnf

count=0
old_num_merges="None"
echo $old_num_merges


./clause_merge_analyzer ${prefix}_temp0.cnf $cmty ${prefix}_merge.r0 > ${prefix}_analyzer.r0
curr_num_merges=`cat ${prefix}_merge.r0 | grep NumMerges`

echo "Iteration 0"
echo $curr_num_merges
echo $curr_num_merges > ${prefix}_merge_main.out

while [ "$old_num_merges" != "$curr_num_merges" ]
do
    #echo "IN"
    new_count=$(( count + 1 ))
    echo "Iteration ${new_count}"
    ./merge_generator ${prefix}_temp${count}.cnf $cmty ${prefix}_temp${new_count}.cnf > ${prefix}_gen_out.r${new_count}
    ./clause_merge_analyzer ${prefix}_temp${new_count}.cnf ${cmty} ${prefix}_merge.r${new_count} > ${prefix}_analyzer.r${new_count}
    count=$new_count
    old_num_merges=$curr_num_merges
    curr_num_merges=`cat ${prefix}_merge.r${new_count} | grep NumMerges`
    echo $curr_num_merges
    echo $curr_num_merges >> ${prefix}_merge_main.out

done


