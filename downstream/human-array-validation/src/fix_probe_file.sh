#!/bin/bash
#
cat - | ruby -ane 'puts $_ if $F[3].match(/^[ACGT]+$/) and $F[4].match(/^[ACGT]+$/)'
