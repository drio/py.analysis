#!/bin/bash

(
cat <<EOF
Chr1	0	10000	2
Chr1	10001	20000	4
Chr2	0	10000	2
Chr2	10001	20000	2
Chr2	20001	30000	2
Chr3	0	10000	0
Chr3	10001	20000	3
EOF
) > truth.bed
cat test.vcf | ./snp_density.py 10000  > foo.bed
diff foo.bed truth.bed > /dev/null
if [ $? -eq 0 ];then
  echo "PASS"
  rm -f foo.bed truth.bed
else
  echo "ERROR"
fi
