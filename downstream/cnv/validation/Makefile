#
# Change for your purposes
THIS_MAKE := $(lastword $(MAKEFILE_LIST))
MAKE_DIR := $(dir $(THIS_MAKE))
SHELL=/bin/bash

PICARD:=$(shell echo $$PICARD)
CNV_PIPE=$(MAKE_DIR)/..
SRC=$(MAKE_DIR)/src
BED_EVENTS=data/chr21.bed
REF=data/chr21.small.fa.gz
X_COV=30
READ_LENGTH=100
CHRM=chr21
WIN_SIZE=1000
THRESHOLD=50
MIN_EVENT_SIZE=100
PADDING=1000
SEED=$(CHRM)_$(WIN_SIZE)_$(THRESHOLD)_$(MIN_EVENT_SIZE)
#
CMD_DEC_REF=gzip -cd $(REF)
G_SIZE:=$(shell $(CMD_DEC_REF) | grep -v "^>" | wc -c)
N_READ_PAIRS:=$(shell echo "($(X_COV)*$(G_SIZE))/$(READ_LENGTH)" | bc)
READS_WIN:=$(shell echo "$(X_COV) * ($(WIN_SIZE)/$(READ_LENGTH))" | bc)

NORMAL=tput sgr0
GREEN=tput setaf 2
YELLOW=tput setaf 3
PURPLE=tput setaf 5
RED=tput setaf 1
OK=$(GREEN); echo ">> [ok]"; $(NORMAL)
FAILURE=($(RED); echo ">> [FAILURE]"; $(NORMAL); exit 1)

define LOG
@$(PURPLE); echo ">>$1"; $(NORMAL)
endef

define HELP_MSG
HELP:
	$$ make clean
	$$ make -f PATH_TO_MAKEFILE BED_EVENTS=data/chr21.bed REF=data/chr21.small.fa.gz X_COV=30 CHRM=chrm21 WIN_SIZE=1000 THRESHOLD=50 MIN_EVENT_SIZE=100 PADDING=1000
endef

export HELP_MSG

all: normal.reads.fq.bam  mutated.reads.fq.bam calls.txt roc.data.$(SEED).txt

# ROC
# min_event_size threshold sensitivity fdr
roc.data.%.txt: calls.txt
	@(echo -ne "$(MIN_EVENT_SIZE) $(THRESHOLD) " ;\
		cat $< | grep -v coor | sort -k2,2n | awk -v c=$(CHRM) '{print c" "$$2" "$$2+$$3" "$$1" "$$4}' | \
		$(SRC)/roc.py -c - -t $(BED_EVENTS) -m $(MIN_EVENT_SIZE) -w $(PADDING) | tail -1 ) > $@

# Call CNV events with our pipeline (3)
# chrm start end size type
calls.txt: mutated.ratios.txt
	$(call LOG, "Calling events")
	$(CNV_PIPE)/join_ratios_file.py $< "([a-z]+).ratios.txt" . | sort -k2,2n | ruby -ne 'puts $$_.gsub(/\s$$/, "")' |\
		$(CNV_PIPE)/downstream.py -i - -n $(READS_WIN) -t $(THRESHOLD) > $@
	@$(OK)

mutated.ratios.txt: bins.txt
	$(call LOG, "Computing Ratios")
	samtools index mutated.reads.fq.bam
	cat $< | $(CNV_PIPE)/compute_log_ratio.py mutated.reads.fq.bam > $@
	@$(OK)

bins.txt: normal.reads.fq.bam
	$(call LOG, "Binning")
	samtools index $<
	$(CNV_PIPE)/binner.py $(READS_WIN) $< $(CHRM) > $@
	@$(OK)

# Generate mutated ref, generate reads, align reads
%.reads.fq.bam: %.reads.fq.sam
	$(call LOG, "Sorting $*")
	java -Xmx512m -jar $(PICARD)/SortSam.jar SORT_ORDER=coordinate INPUT=$< OUTPUT=$@
	@$(OK)

%.reads.fq.sam: %.fa.sa %.reads.fq
	$(call LOG, "Mapping $*")
	bwa mem -t4 normal.fa $*.reads.fq > $@ || $(FAILURE)
	@$(OK)

%.fa.sa: %.fa
	$(call LOG, "Indexing $*")
	bwa index $<

%.reads.fq: %.fa
	$(call LOG, "Generating reads for $*")
	@$(YELLOW)
	@echo "G_SIZE       $(G_SIZE)"
	@echo "N_READ_PAIRS $(N_READ_PAIRS)"
	@echo "READS_WIN    $(READS_WIN)"
	@$(NORMAL)
	wgsim -N $(N_READ_PAIRS) -1$(READ_LENGTH) -d0 -S11 -e0.050 -r0.0020 $< $@ /dev/null > /dev/null
	@$(OK)

mutated.fa: $(BED_EVENTS) normal.fa
	$(call LOG, "Creating new mutated reference")
	@$(SRC)/introduce_events.py -e $(BED_EVENTS) -r normal.fa | $(SRC)/normalize_ref.py > $@
	@$(OK)

normal.fa:
	$(call LOG, "Preparing normal reference.")
	@([[ $(REF) == *.gz ]] && gzip -cd $(REF) | $(SRC)/normalize_ref.py > $@ || cat $(REF) | $(SRC)/normalize_ref.py > $@)
	@$(OK)

wgsim:
	$(call LOG, "Cloning dwsim ...")
	git clone https://github.com/lh3/wgsim.git
	cd wgsim; gcc -g -O2 -Wall -o wgsim wgsim.c -lz -lm
	mv wgsim/wgsim ./w
	rm -rf wgsim
	mv w wgsim
	@$(OK)

test:
	$(SRC)/introduce_events.py -e data/test.bed -r data/test.fa | $(SRC)/normalize_ref.py > _tmp
	@diff _tmp data/test.ok.fa && echo "test ok" || echo "test FAILED!"
	rm _tmp

rebuild_test:
	$(SRC)/introduce_events.py -e data/test.bed -r data/test.fa | $(SRC)/normalize_ref.py > data/test.ok.fa
	make test

clean:
	$(call LOG, "Cleaning up.")
	@rm -f wgsim mutated* reads* normal* *.txt
	@$(OK)

help:
	@echo "make dir: $(dir $(THIS_MAKE))"
	@echo "$$HELP_MSG"

.PHONY: clean test rebuild_test help foo
.PRECIOUS: %.sa ratios.txt
