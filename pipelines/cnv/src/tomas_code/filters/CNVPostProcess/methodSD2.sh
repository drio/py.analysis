#!/bin/bash


# Samples
samples="sample1 sample2 sample3 sample4 sample5 sample6"


# Define paths
PROJECT=/CNV/
DUPS_CALLS=$PROJECT/Final_calls_M2
BED_RMASK=rMask.bed
BED_TRF=trf.sorted.bed
BED_GAPS=gaps.bed

sdT=4

for sample in $samples; do

        CN_CALLS=$PROJECT/mrCNVR/$sample/calls

        avg=`cat $CN_CALLS/$sample.calls.log | grep 'Depth' |cut -f1 -d "," |cut -f2 -d ":"`
        sd=`cat $CN_CALLS/$sample.calls.log | grep 'Depth' |cut -f2 -d "," |cut -f2 -d ":"`

	avgLWnoX=`echo $avg |cut -f1 -d " "`
	avgSWnoX=`echo $avg |cut -f2 -d " "`
	sdLWnoX=`echo $sd |cut -f1 -d " "`
	sdSWnoX=`echo $sd |cut -f2 -d " "`

        WSSD=$DUPS_CALLS/$sample
        mkdir -p $WSSD

        grep -v chrX $CN_CALLS/$sample.calls.sw_norm.bed | grep -v chrM |grep -v chrY | grep -e "^$" -v > $WSSD/$sample.calls.sw_norm_wochrXMY.bed
        grep -v chrX $CN_CALLS/$sample.calls.lw_norm.bed | grep -v chrM |grep -v chrY | grep -e "^$" -v > $WSSD/$sample.calls.lw_norm_wochrXMY.bed
        grep -v chrX $CN_CALLS/$sample.calls.cw_norm.bed | grep -v chrM |grep -v chrY | grep -e "^$" -v > $WSSD/$sample.calls.cw_norm_wochrXMY.bed

        mkdir -p $WSSD/wssd

        thr5=$(echo "scale=9; $avgLWnoX + $sdT * $sdLWnoX" |bc)
        thr1=$(echo "scale=9; $avgSWnoX + $sdT * $sdSWnoX" |bc)

        perl ./wssd_picker.pl -f $WSSD/$sample.calls.lw_norm_wochrXMY.bed -w 7 -s 6 -b 4 -k $WSSD/$sample.calls.sw_norm_wochrXMY.bed -n 5 -i 1 -c $thr5 -t $thr1 -o $WSSD/wssd/$sample.wssd_"$sdT"sd_woChrXMY.tab

        perl ./mergeCoord.pl -f1 $WSSD/wssd/$sample.wssd_"$sdT"sd_woChrXMY.tab -h1 -outF $WSSD/wssd/$sample.wssd_"$sdT"sd_woChrXMY.merged

        cat $WSSD/wssd/$sample.wssd_"$sdT"sd_woChrXMY.merged| awk '{ if($3-$2>=10000) print $0; }' > $WSSD/wssd/$sample.wssd_"$sdT"sd_woChrXMY_10K.merged

        perl ./twoPartOvp_mgsrt.pl -i $WSSD/wssd/$sample.wssd_"$sdT"sd_woChrXMY_10K.merged -f -j $BED_GAPS -t  -L -o $WSSD/wssd/$sample.wssd_"$sdT"sd_woChrXMY_10K_woGaps.merged


done




