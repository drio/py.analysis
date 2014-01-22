INPUT: fasta genome (chrm)

getKmers.sh*
        makeIntervalsBED_v2.pl
        /aplic/bedtools/fastaFromBed -fi $inputFasta -bed $inputBed -fo $outputFasta"

mrsFast_reference_is_multiFASTA.sh

countMappings_allKmers_skip_scaffolds_wo_mrsFast_output.pl

calculateAssemblyStats.sh*
makeHistogramCounts.R*

maskFastaFromBed.sh*


----
#!/bin/bash
scripts=src/tomas_code/maskReference/scripts
fasta=data/chr17.fa
size=64391591
chrm=chr17

$ perl src/step1/makeIntervalsBED_v2.pl $chrm $size > intervals.bed
$ fastaFromBed -fi $inputFasta -bed intervals.bed -fo kmers.fa
