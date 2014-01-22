#!/bin/bash


# Script that:
# (1) Makes BED files with chromosome intervals
# (2) Gets FASTA sequences for the intervals


#print "***** Check assembly paths! *****\n";
#exit;

CURRENT_DIR=/scratch/primate/Assemblies/rhesusMac2/CNV/
mkdir -p $CURRENT_DIR/intervals/qu_cmd/out_mn
mkdir -p $CURRENT_DIR/intervalsFasta/qu_cmd/out_mn
ASSEMBLY="rheMac2"
ASSEMBLY_lc="rheMac2"

while read line
do

	# Get chromosome name and length
	cName=$(echo -e $line | cut -f1 -d' ')
	cLen=$(echo -e $line | cut -f2 -d' ')
	# Skip header	
	if [ "$cName" == "#chrom" ]; then
		continue
	fi
	#echo $cName $cLen
	#continue

	# make intervals: Job building
	jobName=$CURRENT_DIR/intervals/qu_cmd/$cName"_makeIntervals.cmd"
	echo "perl /scratch/primate/Assemblies/rhesusMac2/CNV/scripts/makeIntervalsBED_v2.pl $cName $cLen $ASSEMBLY" > $jobName

	# make intervals: Queue submission
	cd $CURRENT_DIR/intervals/qu_cmd/
	chmod a+x $jobName
	jobId=`qsub -o $CURRENT_DIR/intervals/qu_cmd/out_mn -e $CURRENT_DIR/intervals/qu_cmd/out_mn -l h_vmem=500M -N make.intervals.$cName -cwd -b y $jobName | awk '{print $3}'`
	cd $CURRENT_DIR

	# get fasta from bed: job building
	jobName=$CURRENT_DIR/intervals/qu_cmd/$cName"_fastaFromBed.sh"
	echo "#!/bin/bash" > $jobName
	inputFasta="/scratch/primate/Assemblies/rhesusMac2/CNV/hardMask/"$cName.fa.masked
	inputBed="/scratch/primate/Assemblies/rhesusMac2/CNV/intervals/"$cName"_k36_step5_intervals.bed"
	outputFasta="/scratch/primate/Assemblies/rhesusMac2/CNV/intervalsFasta/"$cName"_k36_step5_intervals.fa.masked"
	jobCmd="/aplic/bedtools/fastaFromBed -fi $inputFasta -bed $inputBed -fo $outputFasta"
	echo $jobCmd > $jobName

	# get fasta from bed: queue submission
	cd $CURRENT_DIR/intervalsFasta/qu_cmd
	chmod a+x $jobName
	qsub -hold_jid $jobId -o $CURRENT_DIR/intervalsFasta/qu_cmd/out_mn -e $CURRENT_DIR/intervalsFasta/qu_cmd/out_mn -l h_vmem=4G -N fasta.from.bed.$cName -cwd $jobName
	cd $CURRENT_DIR

done <"/scratch/primate/Assemblies/rhesusMac2/CNV/chromInfo.txt"
