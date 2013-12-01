#!/bin/bash

mkdir -p fastqc
cd fastqc
fastqc -o . $1
touch done.txt
