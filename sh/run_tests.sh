#!/bin/bash

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $script_dir/common.sh

check_requirements
export PYTHONPATH=$PYTHONPATH:$script_dir/../drio.py
cd $script_dir/..
nosetests tests/
cd ..
