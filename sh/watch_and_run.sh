#!/bin/bash
#
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $script_dir/common.sh

WATCH_THIS="tests/*.py \
drio.py/*.py \
downstream/diff_allele_freq/*.py "

check_requirements
export PYTHONPATH=$PYTHONPATH:$script_dir/../drio.py
cd $script_dir/..
filewatcher "$WATCH_THIS" 'nosetests tests/*.py;echo "\n\n\n"'
cd ..
