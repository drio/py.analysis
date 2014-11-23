#!/usr/bin/env python
#
# vim: set ts=4 sw=4 ft=python:
import sys
import os
from glob import glob
from itertools import chain

def trimed_lines_no_comments(f):
    return chain(l.rstrip() for l in open(f) if l[0] != "#")

def get_fasta(common="../common.sh"):
    for l in trimed_lines_no_comments(common):
        var, val = l.split("=")
        if var == "fasta":
            return val[1:-1]
    return None

def submit(cmd, name, cores=1, mem="8G"):
    print "echo '%s' | submit -s %s -m %s -c %s" % (cmd, name, mem, cores)

def bwa(config, fasta):
    # RG fa fa fq1 fa fq2 fq1 fq2 out_name.bam
    tmpl = 'bwa sampe -r "`cat %s`" %s <(bwa aln -t4 -1 %s %s) <(bwa aln -t4 -2 %s %s) %s %s |  samtools view -Shb /dev/stdin > %s.bam'
    for idx, l in enumerate(open(config)):
        fq1, fq2, lane, lib, sample = l.rstrip().split(' ')
        _id = "%s_%s_%s" % (sample, lane, lib)
        rg_f_name = "rg_%s" % _id
        with open(rg_f_name, 'w') as rg_f:
            rg = "@RG\tID:%s\tPL:ILLUMINA\tLB:%s\tSM:%s" % (_id, lib, sample)
            rg_f.write(rg)
            submit(tmpl % (rg_f_name, fasta, fasta, fq1, fasta, fq2, fq1, fq2, _id), "bwa" + _id, 8, "16G")

def merge():
    bams = " ".join(["INPUT=%s" % b for b in glob("*.bam")])
    cmd = 'java -Xmx12G -jar $PICARD/MergeSamFiles.jar \
        TMP_DIR=/space1/tmp \
        SORT_ORDER=coordinate \
        USE_THREADING=true \
        VALIDATION_STRINGENCY=LENIENT \
        %s \
        OUTPUT=merged.sorted.bam' % (bams)
    submit(cmd, "merge", 3, "16G")

def dups(_id):
    tmpl = 'java -Xmx12G -jar $PICARD/MarkDuplicates.jar \
        TMP_DIR=/space1/tmp \
        METRICS_FILE=/dev/null \
        VALIDATION_STRINGENCY=LENIENT \
        INPUT=merged.sorted.bam \
        REMOVE_DUPLICATES=True \
        OUTPUT=%s.merged.sorted.dups.bam' % (_id)
    submit(tmpl, "dups", 3, "16G")

if __name__ == "__main__":
    config = "config.txt"
    if os.path.isfile(config):
        #bwa("config.txt", get_fasta())
        #merge()
        dups("french")
    else:
        sys.stderr.write("config.txt not found, generating ...\n")
        with open(config, 'w') as f:
            f.write('fq_read1 fq_read2 lane_number library_id sample_id' + "\n")
            f.write('reads/raw/French_1_1.fastq.gz reads/raw/French_1_2.fastq.gz 1 libXXXX French' + "\n")
            f.write('reads/raw/French_2_1.fastq.gz reads/raw/French_2_2.fastq.gz 2 libBBBB French' + "\n")
            f.write('...' + "\n")
