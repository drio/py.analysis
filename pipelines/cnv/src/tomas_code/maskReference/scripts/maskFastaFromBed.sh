#!/bin/bash


# Script that hard-masks FASTA files for the regions given in the BED file


ASSEMBLY="rheMac2"


# Define paths
CURRENT_DIR="/scratch/primate/Assemblies/rhesusMac2/CNV/"
FASTA_FILES=/scratch/primate/Assemblies/rhesusMac2/CNV/hardMask
BED_FILES=/scratch/primate/Assemblies/rhesusMac2/CNV/intervalsToBeMasked
OUTDIR=/scratch/primate/Assemblies/rhesusMac2/CNV/chromFaMasked_kmerMasked
MASKFASTAFROMBED=/aplic/bedtools/maskFastaFromBed
mkdir -p $OUTDIR/qu_cmd/out_mn


for file in `find $BED_FILES/*.bed -type f`; do
	
	chr=`basename $file`;
	chr=${chr/.intervalsToBeMasked.bed/};	
	
	#continue
	BED=$BED_FILES/$chr.intervalsToBeMasked.bed
	#LINES_IN_BED=`cat $BED | wc -l`
	INFASTA=$FASTA_FILES/$chr.fa.masked
	OUTFASTA=$OUTDIR/$chr.fa.kmer.masked
	
	#-s If file has content (is not zero size)
	if [ -s $BED ]; then
		# Case A: regions to be masked
		# Job building	
		jobName=$OUTDIR/qu_cmd/$chr.maskFastaFromBed.sh
		echo "#!/bin/bash" > $jobName
		jobCmd="$MASKFASTAFROMBED -fi $INFASTA -bed $BED -fo $OUTFASTA"
		echo $jobCmd >> $jobName
		# Queue submission
		cd $OUTDIR/qu_cmd
		chmod a+x $jobName
		qsub -o $OUTDIR/qu_cmd/out_mn -e $OUTDIR/qu_cmd/out_mn -l h_vmem=1G -N $chr.maskFastaFromBed -cwd $jobName
		cd $CURRENT_DIR
	else
		# Case B: no regions to be masked; simply copy the fasta file	
		echo "$chr being copied to ../chromFaMasked_kmerMasked"
		cp $INFASTA $OUTFASTA
	fi	

done
