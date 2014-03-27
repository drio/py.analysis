#!/bin/bash
#
# vim: set ts=2 expandtab:
#

cat - | \
  cut -f 8  | \
  tr ";" "\n"  | \
  ruby -ane '
    BEGIN { @h = {}; @id = 1; @g = {}; @prev_hash = false}

    if $_[0] == "#"
      @prev_hash = true
      next
    else
      if @prev_hash
        @prev_hash = false
        if @g.length > 0
          print @id, " ", @g.length, "\n"
          @id += 1
        end
      end
    end

    if ($_[0..2] == "EFF" and $_.split("|")[5] != "")
        gene = $_.split("|")[5]
        @g[gene] = 1
    end

    END { print @id, " ", @g.length, "\n"}
  '
