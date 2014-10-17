#!/bin/bash
#
# vim:set ts=2 ft=sh:
#
# This is an example on how to run the pipeline
#
src="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

usage() {
  local msg=$1
  echo $1
  echo "Usage: `basename $0` <sample_id> <path_to_vcf>"
  exit 1
}

#id=34601
#vcf=/stornext/snfs5/rogers/Fawcett/CRV/phaseI/exomes/vcfs/34601.vcf.anno.gz

[ $# -ne 2 ] && usage "Wrong number of parameters."
id=$1
vcf=$2

[ $id == "." ] && usage "Need sample id."
[ $vcf == "." ] && usage "Need vcf file."
[ ! -f $vcf ] && usage "I can't find vcf file."

[ ! -f inputs/pp.txt ] && usage "I can't find inputs/pp.txt"
[ ! -f inputs/sift.txt ] && usage "I can't find inputs/sift.txt"
[ ! -f inputs/lo.txt.gz ] && usage "I can't find inputs/lo.txt.gz"

output_dir=output
mkdir -p $output_dir

skipping() {
  echo "Skipping, already there."
}

#
# First, we generate a combined version of the predictions
# Notice that we have to do some fixing because some of the input
# records have inconsistent chrmosomal names.
#
echo "Combining predictions ..."
if [ ! -f "./$output_dir/pre.bed" ];then
  $src/damages-predictions.py inputs/pp.txt inputs/sift.txt |\
  sed '1d' |\
  ruby -ane 'if ($_[0] == "c") then puts $_ else puts "chr" + $_ end' > $output_dir/pre.bed
else
  skipping
fi

#
# Now we have to generate the link between the human coordinates and
# hsap.
# But before doing so, we want to generate a more standard version of
# the lift overs, namely a bed format. In that new file we have the
# human coordinates first and then the liftover coordinates for rhesus.
# In 6 columns
#
echo "Fixing lift overs ..."
if [ ! -f "$output_dir/lo.fixed.bed.gz" ];then
  gzip -cd inputs/lo.txt.gz | \
  ruby -ane '
    chrm, start = $F[3].split(":")
    puts $F[0..2].join("\t") + "\t" + chrm + "\t" + start + "\t" + start
  ' | gzip -c - > $output_dir/lo.fixed.bed.gz
else
  skipping
fi

echo "Computing lift overs ..."
cat $output_dir/pre.bed | $src/damages-lift.py $output_dir/lo.fixed.bed.gz > $output_dir/pre-lift.bed

#
# Finally we generate the results by only displaying the SNV that
# have interesting predictions
#
echo "Generating results ..."
extension="${vcf##*.}"
[ "$extension" == "gz" ] && _cat="gzip -cd " || _cat="cat"
cat $output_dir/pre-lift.bed | $src/damages-vcf.py $id <($_cat $vcf | sed 's/Chr/chr/') > $output_dir/$id.bed

echo "Done!"
