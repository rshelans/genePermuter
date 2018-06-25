#!/bin/bash
perm=$1
file=$2

>&2 echo `date`
>&2 echo "run_experiment.sh" ${perm} ${file}
 #export PATH=$PATH:~/drive/manny/bin/
for ((i=1; i<=${perm};i++));
do
	cat $2 |\
	permute_genome -r -b ${i}|\
	find_motif -f GTAYGT TACTAAC -s "(60-400)" 
	#cat tmp |find_motif -f GTNNGN NRYTRAY -s "(40-600)"
done
