#!/bin/bash
perm=$1
ORF=/home/rshelans/drive/manny/data/ALL/ORF
>&2 echo `date`
>&2 echo "kmer_experiment.sh" ${perm}

echo "HELLO" | kmer -d > perm_cer_kmer.tab
echo "HELLO" | kmer -d > perm_mik_kmer.tab
echo "HELLO" | kmer -d > perm_bay_kmer.tab

for ((i=1; i<=${perm};i++));
do
	cat ${ORF}/Scer.fasta|\
	permute_genome -b ${i}|\
	kmer >>perm_cer_kmer.tab
	cat ${ORF}/Sbay.fasta|\
	permute_genome -b ${i}|\
	kmer >>perm_bay_kmer.tab
	cat ${ORF}/Smik.fasta|\
	permute_genome -b ${i}|\
	kmer >>perm_mik_kmer.tab
done
