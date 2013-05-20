#!/bin/bash
#
# snpn qual > 200
# both samples have to have more than 10 reads with MAPQ>=1 at locus
#
gzip -cd merged.vcf.anno.rdp.gz | \
  awk '{if (substr($0,0,0) || $6>200) print;}' |\
  ruby -ane 'puts $_; next if $_[0]=="#"; o,t = $F[6].match(/RDP=([\d,]+)/)[1].split(","); print $_ if o.to_i > 10 and t.to_i > 10' | \
  gzip -c - > papio.vcf.anno.rdp.q200.rdp10.mapq10.gz
