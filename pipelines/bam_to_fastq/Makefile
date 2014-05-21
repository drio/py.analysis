#
SHELL=/bin/bash
THIS_MAKE:=$(lastword $(MAKEFILE_LIST))
MAKE_DIR:=$(dir $(THIS_MAKE))
TMP:=$(shell [ -d /space1/tmp ] && echo "/space1/tmp" || echo "/tmp")
HTS:=/stornext/snfs6/rogers/drio_scratch/local/bin/htscmd
JAVA_OPTS=-Xmx4G
PICARD_OPTS=MAX_RECORDS_IN_RAM=5000000 TMP_DIR=$(TMP) VALIDATION_STRINGENCY=LENIENT
ID:=$(shell ls *.bam 2>/dev/null |sed -s 's/.bam//')

define HELP

Your ID is:=$(ID)

First thing first:
	_ ln -s $(THIS_MAKE)

# For testing:
  _ samtools view -h /stornext/snfs5/rogers/drio_scratch/bams/Cyno_OConnor/stornext/snfs5/next-gen/solid/p-merge-work/merges/CIP_CYNO_WG_ILM/11200-03/11200-03.dups.bam | \
	head -100000 | samtools view -h -bS - > test.bam
  _ make test.fq.tar.bz2

# For real data:
  _ mkdir workdir; cd workdir
  _ ln -s /stornext/snfs5/rogers/drio_scratch/bams/Cyno_OConnor/stornext/snfs5/next-gen/solid/p-merge-work/merges/CIP_CYNO_WG_ILM/11200-03/11200-03.dups.bam ./11200-03.bam
  _ make 11200-03.fq.tar.bz2

# The pipeline generates a file called check.txt. Use the metrics inside
  to determine the correctness of the results.
endef

export HELP
help:
	@echo "$$HELP"


$(ID).fq.tar.bz2: check.txt
	tar -jcvf $@ $<
	rm -f *_tmp *.fq

check.txt: $(ID).1.fq $(ID).2.fq unpaired.fq reads_in_bam.txt reads.fq.1 reads.fq.2
	@(_1=`cat reads.fq.1`;\
	_2=`cat reads.fq.2`;\
	_B=`cat reads_in_bam.txt`;\
	_DF=`echo "$$_1-$$_2" | bc`;\
	_DBF=`echo "$$_B-($$_1+$$_2)" | bc`;\
	echo "id n_reads_bam n_reads_fq_1 n_reads_fq2 diff_fq diff_bam_fq" > $@;\
	echo "$(ID) $$_B $$_1 $$_2 $$_DF $$_DBF" >> $@)

%.1.fq:
	$(HTS) bamshuf -O -u $*.bam _tmp | \
	java $(JAVA_OPTS) -jar $$PICARD/SamToFastq.jar I=/dev/stdin F=$*.1.fq F2=$*.2.fq $(PICARD_OPTS)

%.2.fq:
	@echo "noop" > /dev/null

unpaired.fq:
	@echo "noop" > /dev/null

reads_in_bam.txt: $(ID).bam
	samtools view $< | wc -l > $@

reads.fq.%:
	_wc=`wc -l $(ID).$*.fq | awk '{print $$1}'`; echo "$$_wc/4" | bc > $@

test.bam:
	samtools view -h /stornext/snfs5/rogers/drio_scratch/bams/Cyno_OConnor/stornext/snfs5/next-gen/solid/p-merge-work/merges/CIP_CYNO_WG_ILM/11200-03/11200-03.dups.bam | head -100000 | samtools view -hSb - > test.bam

clean:
	rm -f reads* check.txt *.fq _tmp*

.PHONY: help clean

