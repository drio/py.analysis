#!/bin/bash
#
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $script_dir/common.sh

WATCH_THIS=`find tests -type f; find downstream -type f; find drio.py -type f`

check_requirements
export PYTHONPATH=$PYTHONPATH:$script_dir/../drio.py
cd $script_dir/..
filewatcher "$WATCH_THIS" './sh/sync_ardmore.sh'
cd ..
