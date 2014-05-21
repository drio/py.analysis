#
SHELL=/bin/bash
THIS_MAKE:=$(lastword $(MAKEFILE_LIST))
MAKE_DIR:=$(dir $(THIS_MAKE))
TMP:=$(shell [ -d /space1/tmp ] && echo "/space1/tmp" || echo "/tmp")
HTS:=/stornext/snfs6/rogers/drio_scratch/local/bin/htscmd
JAVA_OPTS=-Xmx4G
PICARD_OPTS=MAX_RECORDS_IN_RAM=5000000 TMP_DIR=$(TMP) VALIDATION_STRINGENCY=LENIENT
ID=$(shell ls *.bam |sed -s 's/.bam//')

define HELP
First thing first:
	_ ln -s $(THIS_MAKE)

# For testing:
  _ samtools view -h /stornext/snfs5/rogers/drio_scratch/bams/Cyno_OConnor/stornext/snfs5/next-gen/solid/p-merge-work/merges/CIP_CYNO_WG_ILM/11200-03/11200-03.dups.bam | \
	head -100000 | samtools view -bS - > test.bam
  _ make test.fq.tar.bz2

# For real data:
  _ mkdir workdir; cd workdir
  _ ln -s /stornext/snfs5/rogers/drio_scratch/bams/Cyno_OConnor/stornext/snfs5/next-gen/solid/p-merge-work/merges/CIP_CYNO_WG_ILM/11200-03/11200-03.dups.bam ./11200-03.bam
  _ make 11200-03.fq.tar.bz2
endef

export HELP
help:
	@echo "$$HELP"

$(ID).fq.tar.bz2: $(ID).1.fq $(ID).2.fq
	tar -jcvf $@ $<
	rm -f *_tmp *.fq

%.1.fq:
	$(HTS) bamshuf -O -u $*.bam _tmp | \
	java $(JAVA_OPTS) -jar $$PICARD/SamToFastq.jar I=/dev/stdin F=$*.1.fq F2=$*.2.fq $(PICARD_OPTS)

%.2.fq:
	@echo "noop" > /dev/null


.PHONY: help


