#!/bin/bash


# Specify parameters
assembly="rheMac2"
#version="chromFa"
version="hardMask"
#version="chromFaMasked_kmerMasked"


# Define paths
ASSEMBLY_PATH=/scratch/primate/Assemblies/rhesusMac2/CNV/$version
OUTDIR=$ASSEMBLY_PATH/AssemblyStats
mkdir -p $OUTDIR/out_mn
mkdir -p $OUTDIR/tmp

# Output file
outfile=$OUTDIR/assemblyStats.txt
logfile=$OUTDIR/assemblyStats.log
echo -e "chromosome\tTotal\tA\tC\tG\tT\tN\tlower_case\tother" > $outfile


jobsIds=""


# Calculate statistics for each chromosome
for fasta in `ls $ASSEMBLY_PATH/chr* | grep -v fai `; do
	chr=`basename $fasta`;
	chr=${chr/.fa.masked/};

	echo $chr
	#continue
	jobId=`qsub -e $OUTDIR/out_mn -o $OUTDIR/out_mn -l h_vmem=1G -cwd -N $chr.calculate.assembly.stats -b y perl calculateAssemblyStats_v2.pl $fasta $assembly $version $chr $assembly | awk '{print $3}'`
	jobsIds=$jobsIds","$jobId
done


# Combine results when jobs finished
echo '#!/bin/bash/' > $OUTDIR/combine.files.sh
echo "cat $OUTDIR/tmp/tmp.chr* >> $OUTDIR/assemblyStats.txt" >> $OUTDIR/combine.files.sh
echo "rm -f $OUTDIR/tmp/tmp*" >> $OUTDIR/combine.files.sh
qsub -hold_jid $jobsIds -e $OUTDIR/out_mn -o $OUTDIR/out_mn -l h_vmem=500M -cwd -N combine.assembly.stats $OUTDIR/combine.files.sh


