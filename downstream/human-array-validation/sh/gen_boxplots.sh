#!/bin/bash

rm -rf boxplots
rm -f box*.png mapq*.png

for s in Omni affy6 exon
do
  # Generate the boxplot graphs for the eg cov values.
  echo ">> eg $s"
  ./src/boxplots.py "*$s*.txt" 2> /dev/null
  mv boxplot.png boxplot.$s.png

  # Generate the boxplot for the mapq values
  #echo ">> mapq $s"
  #grep -vP "^@" *$s*.sam  | awk '{print $5}' > _mapq.$s.txt
  #./src/boxplots.py "_mapq.*$s*.txt"
  #mv boxplot.png mapq.$s.png
  #rm -f _mapq*.txt
done

mkdir boxplots
mv boxplot*.png boxplots/
#mv mapq*.png boxplots/
echo "
<img src='boxplot.Omni.png'><br>
<img src='boxplot.affy6.png'><br>
<img src='boxplot.exon.png'><br>
<-- <img src='mapq.Omni.png'><br> -->
<-- <img src='mapq.affy6.png'><br> -->
<-- <img src='mapq.exon.png'><br> -->
" > boxplots/index.html
chmod 755 boxplots/*

scp -r boxplots drio@is04607.com:public_html



