#!/bin/bash
perm=$1
file=$2
tmpfile=`strings /dev/urandom | grep -o '[[:alnum:]]' | head -n 30 | tr -d '\n'; echo`
>&2 echo `date`
>&2 echo "run_experiment.sh" ${perm} ${file}
 #export PATH=$PATH:~/drive/manny/bin/
echo -e "sense\tantisence"
for ((i=1; i<=${perm};i++));
do
	cat $2 |\
	permute_genome -b ${i} >${tmpfile}
	sence=`cat ${tmpfile}|find_motif -f GTAYGT TACTAAC -s "(60-400)" `
	antisence=`cat ${tmpfile} |find_motif -rf  GTAYGT TACTAAC -s "(40-600)"`
	echo -e ${sence} "\t" ${antisence}
	rm ${tmpfile}
done
