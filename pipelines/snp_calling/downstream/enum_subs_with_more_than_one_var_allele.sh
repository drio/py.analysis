#!/bin/bash
#
cat - | \
  grep -v INDEL | \
  awk '{if ($6>50) print;}' | \
  ruby -ane 'puts $_ if $F[4].chomp.gsub(/,/,"").size > 1 || $_[0] == "#"'
