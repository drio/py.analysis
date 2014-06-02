#!/bin/bash
#
# vim: ts=2 et:
#
norm="/stornext/snfs6/rogers/drio_scratch/dev/py.analysis/pipelines/cnv/src/validation/normalize_windows.py"

usage() {
  local msg=$1
  echo "ERROR: $msg" >&2
  echo "Usage: $0 <bed_file1 (truth)> <bed_file2 (raw_canavar_out)>"
  exit 1
}

truth=$1
cvar=$2
[ $# -ne 2 ] && usage "Wrong num of arguments ($# instead of 2)"
for f in $cvar $truth
do
  [ ! -f "$f" ] && usage "Input files do not exists"
done

mkdir -p all
for chrm in `seq 1 22`
do
  #cmd="cut -f1,2,3,5 $cvar | sed '1,2d' | $norm $truth - chr$chrm > all/$chrm.bed"
  cmd="$norm $truth $cvar chr$chrm > all/$chrm.bed"
  echo $cmd | submit -s norm.$chrm -m8G "$cmd"
done

