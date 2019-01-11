annotated_introns.bed:
	04/25/15
	Set of Annotated Introns from the Ares Lab and previous Literature.
	orf_genomic.fasta.gz:
	Pulled from SGD 04/25/15
	ORF sequences from the initial ATG to the stop codon, including
	intron sequences and any bases not translated due to translational
	frameshifting, but not including 5'-UTR or 3'-UTR sequences, for all 
	"Verified" and "Uncharacterized" ORFs, and transposable element genes. 
	Does NOT include sequences for ORFs classified as "Dubious" or "pseudogene".



FILTERS:::

introns.tab:
	From annotated_introns.bed: However only one occurance of a gene is allowed.
	cut -f5 annotated_introns.bed|grep -v "chrM"|uniq|awk -F"chr" '{print $1}'|uniq|sed 's/\(.*\)_/\1/'|tail -n +2 >introns.tab 
	
	+YNL194C, YNL194C
transposons.tab:
	List of transposons from orf_genomic
	grep "transposon" orf_genomic.fasta | cut -d' ' -f2 >transposons.tab

mito.tab
	List of Mitochondrial genes in orf_genomic
	grep "Chr M" orf_genomic.fasta|cut -d' ' -f2 >mito.tab

improper.tab
	List of Improper Orfs from orf_genomic

remove.tab
	List to remove from all curated sources



curated_orf.fasta
	##Removed ORFS from orf set
	cat orf_genomic.fasta |select.py -r remove.tab >curated_orf.fasta 
	remove.tab
	Removed 413 from 5917 expected 428 Missed 15.
	DB639721
	YDR535C
	DB639731
	DB641776
	SNR17A
	YOR318C
	DB639229
	YLR202C
	SNR17B
	Q0070
	Q0075
	Q0105
	Q0255
	Q0060
	Q0045










