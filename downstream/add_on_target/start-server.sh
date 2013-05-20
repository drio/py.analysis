#!/bin/bash

#time gzip -cd rhemac.lifted_over.human.vchrome.targets.bed.gz | sed -e 's/:/\t/' -e 's/-/\t/' -e 's/Chr//g' | go run bed_target_server.go
echo "Example of curl query for testing:"
echo 'curl -X POST -H 'Content-Type: application/json' -d {"Sites":[{"Chrm":"1","Start":100}]} http://localhost:8080'
cat small.bed | \
sed -e 's/:/\t/' -e 's/-/\t/' |\
/usr/local/go/bin/go run bed_target_server.go
