#!/bin/bash
set -e

usage()
{
  local msg=$1
  local tool_name=$( basename ${BASH_SOURCE[0]} )
  [ ".$msg" != "." ] && echo "Ups!: $msg"
  cat <<EOF
Usage: $tool_name -e <vcf.wes.gz> -g <vcf.wgs.gz> [-f func_cons] [-v]

OPTIONS:
  -h  help
  -g  [STRING] wgs vcf file (gz)
  -e  [STRING] wes vcf file (gz)
  -f  [STRING] only report snps with this functional consequence
  -v  dump vcf entries instead gene report
  -t  enable test mode
EOF
  exit 1
}

vcf_wes=""
vcf_wgs=""
only_this_fc="cat -"
vcf_output=0
enable_test="cat -"
while getopts “e:g:f:tvh” OPTION
do
 case $OPTION in
   h)
     usage
     exit 1
     ;;
   e)
     vcf_wes=$OPTARG
     ;;
   g)
     vcf_wgs=$OPTARG
     ;;
   f)
     only_this_fc="grep EFF=$OPTARG"
     ;;
   v)
     vcf_output=1
     ;;
   t)
    enable_test="head -10000"
     ;;
   ?)
     usage
     exit
     ;;
   esac
done

[ ! -f "$vcf_wes" ] && usage "Cannot find wes file."
[ ! -f "$vcf_wgs" ] && usage "Cannot find wgs file."

DIR=$( dirname "${BASH_SOURCE[0]}" )
SRC_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/.."
source $SRC_DIR/common.sh

##INFO=<ID=EFF,Number=.,Type=String,Description="Predicted effects for this variant.Format:
# 'Effect ( Effect_Impact 0 | Functional_Class 1 | Codon_Change 2 | Amino_Acid_change 3 | Amino_A
# cid_length | Gene_Name 5 | Gene_BioType | Coding | Transcript | Exon [ | ERRORS | WARNINGS ] )'">

(
gzip -cd $vcf_wes | $enable_test | $SRC_DIR/population.filter.sh "wes"
gzip -cd $vcf_wgs | $enable_test | $SRC_DIR/population.filter.sh "wgs"
) | sort -T/tmp -S8G -k1,1 -k2,2n |\
  ruby -ane '
    # Collapse snp present in both wes and wgs
    BEGIN {@prev=[]}
    if @prev.size == 0
      @prev = $F
    else
      ch, coor = $F[0], $F[1];
      if @prev[0] == ch && @prev[1] == coor
        puts $F.join("\t")
        @prev=[]
      else
        puts @prev.join("\t")
        @prev = $F
      end
    end
  ' | $only_this_fc |\
if [ $vcf_output == 0 ];then
  awk -F\| '{if ($6!="") print}' |\
  awk -F\| '{print $6}' | sort | uniq -c | sort -k1,1nr
else
  cat -
fi
