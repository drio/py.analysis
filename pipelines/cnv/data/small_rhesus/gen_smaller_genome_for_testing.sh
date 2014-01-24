#!/bin/bash

genome=rhemac2
nbp=$(echo "100000" | bc)
n_lines=$(echo "($nbp/50)+1" | bc)

out_seed=small.$genome.$nbp
rm -f $out_seed*

for f in softMask/*.fa
do
  chrm=$(echo $f | sed 's/softMask\///g')
  chrm=$(echo $chrm | sed 's/.fa//g')
  head -$n_lines $f >> $out_seed.fa
  echo -e "$chrm\t$nbp\t$nbp" >> $out_seed.sizes.bed
done


