#
SHELL=/bin/bash
THIS_MAKE := $(lastword $(MAKEFILE_LIST))
MAKE_DIR := $(dir $(THIS_MAKE))
SRC:=$(MAKE_DIR)/src
# Change this as necessary
GTF="ftp://ftp.ensembl.org/pub/release-73/gtf/macaca_mulatta/Macaca_mulatta.MMUL_1.73.gtf.gz"
SOT=/stornext/snfs6/rogers/drio_scratch/dev/py.analysis/downstream/target_coverage/stats_on_target.py
SHELL_NUM_SAMPLES=$(shell echo $$NUM_SAMPLES)

# PARTI: individual filtering
$(ID).out.gz $(ID).pass.gz: $(ID).stats.gz
	(gzip -cd $(ID).stats.gz | $(SRC)/filter.py | gzip -c > $(ID).pass.gz) 2>&1 | gzip -c > $(ID).out.gz

$(ID).stats.gz: $(ID).depth.gz
	$(SOT) exons.bed $< | gzip -c > $(ID).stats.gz

$(ID).depth.gz: exons.bed
	samtools depth -Q 20 -b $< $(BAM) | gzip -c > $(ID).depth.gz

exons.bed: gtf.gz
	gzip -cd $< | $(SRC)/gtf2bed.py | sed 's/chr/Chr/g' > $@

gtf.gz:
	curl $(GTF) > $@

# PARTII: aggregated filtering
all.bed:
	$(SRC)/matrix.py > all.bed

all.pass.%.bed: all.bed
	awk '{if ($$7 >= $*) print;}' $< > $@

all.out.%.bed: all.bed
	awk '{if ($$7 < $*) print;}' $< > $@
