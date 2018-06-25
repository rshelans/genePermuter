#!/bin/bash
perm=$1

inter=/home/rshelans/drive/manny/data/ALL/INTERGENIC
orf=/home/rshelans/drive/manny/data/ALL/ORF

>&2 echo `date`
>&2 echo "kmer_experiment.sh" ${perm}

echo ""  > perm_cer_orf_weak.tab
echo ""  > perm_mik_orf_weak.tab
echo ""  > perm_bay_orf_weak.tab

echo -e "SENSE\tANTISENSE"  > perm_cer_inter_weak.tab
echo -e "SENSE\tANTISENSE"  > perm_mik_inter_weak.tab
echo -e "SENSE\tANTISENSE"  > perm_bay_inter_weak.tab

for ((i=1; i<=${perm};i++));
do
	cat ${orf}/Scer.fasta|permute_genome -b ${i}|find_motif -f GYWHGN NRYTRAY -s "(60-400)" >> perm_cer_orf_weak.tab
	cat ${orf}/Sbay.fasta|permute_genome -b ${i}|find_motif -f GYWHGN NRYTRAY -s "(60-400)" >> perm_bay_orf_weak.tab
	cat ${orf}/Smik.fasta|permute_genome -b ${i}|find_motif -f GYWHGN NRYTRAY -s "(60-400)" >> perm_mik_orf_weak.tab
	
	cat ${inter}/scer.fasta|permute_genome -b ${i} > tmp
	sence=`cat tmp|find_motif -f GYWHGN NRYTRAY -s "(60-400)"`
	antisence=`cat tmp |find_motif -rf GYWHGN NRYTRAY -s "(60-400)"`
	echo -e ${sence} "\t" ${antisence} >> perm_cer_inter_weak.tab

	rm tmp

	cat ${inter}/smik.fasta|permute_genome -b ${i} > tmp
	sence=`cat tmp|find_motif -f GYWHGN NRYTRAY -s "(60-400)"`
	antisence=`cat tmp |find_motif -rf GYWHGN NRYTRAY -s "(60-400)"`
	echo -e ${sence} "\t" ${antisence} >> perm_mik_inter_weak.tab

	rm tmp

	cat ${inter}/sbay.fasta|permute_genome -b ${i} > tmp
	sence=`cat tmp|find_motif -f GYWHGN NRYTRAY -s "(60-400)"`
	antisence=`cat tmp |find_motif -rf GYWHGN NRYTRAY -s "(60-400)"`
	echo -e ${sence} "\t" ${antisence} >> perm_bay_inter_weak.tab

	rm tmp
done
