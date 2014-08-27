#!/bin/bash
#
input_truth="truth/originals/Homo_sapiens-French_HGDP00521_1000bp_simple_per_1kb.HM.bedgraph.bed"
input_pcalls="pipeline.calls.bed"
threshold_min_truth=-1
threshold_max_truth=3
threshold_min_pcalls=-1
threshold_max_pcalls=3
sort="sort -k1,1 -k2,2n"
merge="bedtools merge -i stdin"

compute_overlaps() {
  local a=$1
  local b=$2

  n_events_a=$(cat $a | wc -l)
  bp_events_a=$(cat $a | awk 'BEGIN{a=0}; {a+=$3-$2}; END{print a}')
  n_events_b=$(cat $b | wc -l)
  bp_events_b=$(cat $b | awk 'BEGIN{a=0}; {a+=$3-$2}; END{print a}')

  n_events_in_a_overlapping_with_b=$(bedtools intersect -a $a -b $b -wa -u | wc -l)
  bp_events_in_a_overlapping_with_b=$(bedtools intersect -a $a -b $b -wa -u | awk 'BEGIN{a=0}; {a+=$3-$2}; END{print a}')

  pro=$(echo "scale=2;$n_events_in_a_overlapping_with_b/$n_events_a"| bc)
  bp_pro=$(echo "scale=2;$bp_events_in_a_overlapping_with_b/$bp_events_a"| bc)
  echo $pro $bp_pro $n_events_a $bp_events_a $n_events_b $bp_events_b $n_events_in_a_overlapping_with_b $bp_events_in_a_overlapping_with_b
}

echo "$input_truth ($threshold_min_truth, $threshold_max_truth)"
echo "$input_pcalls ($threshold_min_pcalls, $threshold_max_pcalls)"
echo ""

cat $input_truth | sed 's/10+/10/' | sed 's/chr//'  | \
  awk -v min=$threshold_min_truth -v max=$threshold_max_truth '{if ($4<= min || $4>= max) print $1"\t"$2"\t"$3}' | $sort | $merge > ./A.bed &
sed '1,2d' $input_pcalls | \
  awk -v min=$threshold_min_pcalls -v max=$threshold_max_pcalls '{if ($5<= min || $5>= max) print $1"\t"$2"\t"$3}' | $sort | $merge > ./B.bed &

wait

echo "overlap bp_overlap num_events_a bp_events_a num_events_b bp_events_b num_overlaps_a_b bp_overalps_a_b"
echo "--"
echo -ne "proportion of events in Truth that overlap with Our calls: "
compute_overlaps A.bed B.bed
echo -ne "proportion of events in Our calls that overlap with A Truth: "
compute_overlaps B.bed A.bed
