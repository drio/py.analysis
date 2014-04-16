#!/bin/bash
# vim: set ts=2 et:
#
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

help() {
  local msg=$1
  echo "UPS!: $msg"
  exit 1
}

which sapi.py &>/dev/null; [ $? -ne 0 ] && help "sapi.py not found."
which submit &>/dev/null; [ $? -ne 0 ] && help "submit not found."

cat $DIR/bootstrap/config.sh > config.sh
cat $DIR/bootstrap/second_part.sh > second_part.sh
cat $DIR/bootstrap/run.sh > run.sh
chmod 755 second_part.sh run.sh

echo "All good! Make sure you edit ./config.sh"
