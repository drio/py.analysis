#!/bin/bash
set -e

usage()
{
  local msg=$1
  local tool_name=$( basename ${BASH_SOURCE[0]} )
  [ ".$msg" != "." ] && echo "Ups!: $msg"
  cat <<EOF
Usage: $tool_name -e <vcf.wes.gz> -g <vcf.wgs.gz>

OPTIONS:
  -h  help
  -g  wgs vcf file (gz)
  -e  wes vcf file (gz)
EOF
  exit 1
}

vcf_wes=""
vcf_wgs=""
while getopts “e:g:h” OPTION
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

(
gzip -cd $vcf_wes | $SRC_DIR/population.filter.sh "wes" | awk '{print $1" "$2" e"}'
gzip -cd $vcf_wgs | $SRC_DIR/population.filter.sh "wgs" | awk '{print $1" "$2" g"}'
) | sort -T/tmp -S8G -k1,1 -k2,2n |\
  ruby -ane '
    BEGIN {@prev=[]}
    if @prev.size == 0
      @prev = $F
    else
      ch, coor, what = $F;
      if @prev[0] == ch && @prev[1] == coor
        puts "eg"
        @prev=[]
      else
        puts @prev[2]
        @prev = $F
      end
    end
  ' | sort -k1,1 | uniq -c
