#!/bin/bash

mkdir fastqc
cd fastqc
fastqc -o . $1
touch done.txt
