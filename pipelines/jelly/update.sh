#!/bin/bash

error() {
  local msg=$1
  echo "$msg"
  exit 1
}

[ ! -f ./stats/jelly.stats.txt ] && error "jelly stats not found."
[ ! -f ./stats/input.stats.txt ] && error "input stats not found."

echo -ne "Genome:"
read genome
echo -ne "Subdir name:"
read subdir

main="/home/drio/public_html/primates/jelly"
dir=$main/$genome/$subdir
ssh drio@is04607.com mkdir -p $dir
scp stats/*stats.txt drio@is04607.com:$dir/

echo "http://is04607.com/primates/jelly/$genome/$subdir"

