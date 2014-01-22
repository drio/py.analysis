#!/bin/bash


# Script that:
# (1) Make index files for each chromosome
# (2) Map kmers to chromosomes


# Define paths
CURRENT_DIR=/scratch/primate/Assemblies/rhesusMac2/CNV/
ASSEMBLY="rheMac2"
ASSEMBLY_lc="rheMac2"
MRSFAST="/aplic/mrsfast/mrsfast"


mkdir -p $CURRENT_DIR/mrsFast/qu_cmd/out_mn


# Index
# Job building
jobName=$CURRENT_DIR/mrsFast/qu_cmd/index.sh
echo "#!/bin/bash" > $jobName
# window size calculated as in http://mrfast.sourceforge.net/manual.html
# length=36
# e=2
jobCmd="$MRSFAST --index /scratch/primate/Assemblies/rhesusMac2/CNV/rheMac2.masked.fa --ws 12"
echo $jobCmd >> $jobName
# Queue submission
cd $CURRENT_DIR/mrsFast/qu_cmd
chmod a+x $jobName
jobId=`qsub -o $CURRENT_DIR/mrsFast/qu_cmd/out_mn -e $CURRENT_DIR/mrsFast/qu_cmd/out_mn -l h_vmem=4G -N vmem=4G.index.mr.fast -cwd $jobName | awk '{print $3}'`
cd $CURRENT_DIR


# Mapping
for i in `ls /scratch/primate/Assemblies/rhesusMac2/CNV/intervalsFasta/*k36_step5_intervals.fa.masked`; do
	chrName0=`echo $i | cut -f8 -d'/'`
	chrName=${chrName0//"_k36_step5_intervals.fa.masked"/""}
	#echo $chrName	
	#continue
	# Job building
	jobName=$CURRENT_DIR/mrsFast/qu_cmd/$chrName"_mapping.sh"
	echo "#!/bin/bash" > $jobName
	jobCmd="$MRSFAST --search /scratch/primate/Assemblies/rhesusMac2/CNV/rheMac2.masked.fa --seq $i -o $CURRENT_DIR/mrsFast/$chrName"_k36_step5_intervals.map" -e 2"
	echo $jobCmd >> $jobName
	# Queue submission
	cd $CURRENT_DIR/mrsFast/qu_cmd
	chmod a+x $jobName
	qsub -hold_jid $jobId -o $CURRENT_DIR/mrsFast/qu_cmd/out_mn -e $CURRENT_DIR/mrsFast/qu_cmd/out_mn -l h_vmem=8G -N vmem=8G.mapping.$chrName.mrs.fast -cwd $jobName
	cd $CURRENT_DIR
done





