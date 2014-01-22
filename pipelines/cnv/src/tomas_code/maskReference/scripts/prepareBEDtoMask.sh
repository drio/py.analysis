#!/bin/bash

# Script that for a single input BED file with regions to be masked makes, for each chromosome, BED files with merged regions
# The majority of the regions in the input file correspond to chrUn so when looping over chromosomes, bedtools will fail because of file size
# In the case of chrUn, comment accordingly at chromosomes definition and use bedops 


# Define paths
CWD="/scratch/primate/Assemblies/rhesusMac2/CNV/"
CHROMOSOME_FILE="/scratch/primate/Assemblies/rhesusMac2/CNV/chromInfo.txt"

INPUT=$CWD/intervalsToBeMasked.bed
INPUT_TMP="\${TMP}/intervalsToBeMasked.bed"
OUTDIR=$CWD/intervalsToBeMasked
OUTDIR_TMP="\${TMP}"
mkdir -p $OUTDIR/qu_cmd/out_mn


# Define chromosomes
while read line;
do
	chr=`echo $line | awk '{print $1}'`;
	# Job building
	jobName=$OUTDIR/qu_cmd/$chr.prepareBEDtoBeMasked.sh
	echo "#!/bin/bash" > $jobName
        echo "cp $INPUT $INPUT_TMP" >> $jobName
        echo "mkdir -p $OUTDIR_TMP" >> $jobName
	#jobCmd="cat $INPUT | grep $chr | /aplic/bedops/sort-bed - | /aplic/bedops/bedops --merge - > $OUTDIR/$chr.intervalsToBeMasked.bed"
	jobCmd="cat $INPUT_TMP | grep -P '$chr\t' | /aplic/bedops/sort-bed - | /aplic/bedops/bedops --merge - > $OUTDIR_TMP/$chr.intervalsToBeMasked.bed"
        echo $jobCmd >> $jobName
	echo "cp $OUTDIR_TMP/$chr.intervalsToBeMasked.bed $OUTDIR" >> $jobName
	
	# Queue submission
	cd $OUTDIR/qu_cmd
	chmod a+x $jobName		
	# All chromosomes but chrUn
	qsub -o $OUTDIR/qu_cmd/out_mn -e $OUTDIR/qu_cmd/out_mn -cwd -l h_vmem=10G,local_disk=800M -N vmem=10G.$chr.prepare.bed.to.mask $jobName
	# chrUn
	#qsub -o $OUTDIR/qu_cmd/out_mn -e $OUTDIR/qu_cmd/out_mn -cwd -l h_vmem=20G -N vmem=20G.$chr.prepare.bed.to.mask $jobName

	cd $CWD
done < $CHROMOSOME_FILE











