#!/bin/bash
#
src_dir=$(dirname $(readlink -f $0))
th_over=$1

[ ".$th_over" == "." ] && echo "Need threshold." && exit 1

cat counts* | xvfb-run R CMD BATCH $src_dir/makeHistogramCounts.R
mv cumu.dist.png cumu.before.png
mv histograms.png hist.before.png

cat counts* | awk -v th=$1 '{ if ($5>th) print;}' > all.$th_over.counts
cat all.$th_over.counts | xvfb-run R CMD BATCH $src_dir/makeHistogramCounts.R
mv cumu.dist.png cumu.after.png
mv histograms.png hist.after.png
