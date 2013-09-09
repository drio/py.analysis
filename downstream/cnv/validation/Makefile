#
# (*) Change for your purposes
SHELL=/bin/bash
SRC=src
PICARD:=$(shell echo $$PICARD)
CNV_PIPE=..
BED_EVENTS=data/test.bed
REF=data/chr21.small.fa.gz
X_COV=30
READ_LENGTH=100
CHRM=chr21
WIN_SIZE=1000
THRESHOLD=20
#
CMD_DEC_REF=gzip -cd $(REF)
G_SIZE:=$(shell $(CMD_DEC_REF) | grep -v "^>" | wc -c)
N_READ_PAIRS:=$(shell echo "($(X_COV)*$(G_SIZE))/$(READ_LENGTH)" | bc)
READS_WIN:=$(shell echo "$(X_COV) * ($(WIN_SIZE)/$(READ_LENGTH))" | bc)

all: mutated.reads.fq.bam normal.reads.fq.bam mutated.calls.txt

# Call CNV events with our pipeline (3)
%.calls.txt: %.ratios.txt
	$(CNV_PIPE)/join_ratios_file.py $< "([a-z]+).ratios.txt" . | sort -k2,2n | ruby -ne 'puts $$_.gsub(/\s$$/, "")' |\
		$(CNV_PIPE)/downstream.py -i - -n $(READS_WIN) -t $(THRESHOLD) > $@

%.ratios.txt: %.bins.txt
	cat $< | $(CNV_PIPE)/compute_log_ratio.py $*.reads.fq.bam > $@

%.bins.txt: %.reads.fq.bam
	samtools index $<
	$(CNV_PIPE)/binner.py $(READS_WIN) $< $(CHRM) > $@

# Generate mutated ref, generate reads, align reads
%.reads.fq.bam: %.reads.fq.sam
	java -Xmx512m -jar $(PICARD)/SortSam.jar SORT_ORDER=coordinate INPUT=$< OUTPUT=$@

%.reads.fq.sam: %.fa.sa %.reads.fq
	bwa mem -t4 $(basename $<) $*.reads.fq > $@

%.fa.sa: %.fa
	bwa index $<

%.reads.fq: %.fa
	@echo "G_SIZE       $(G_SIZE)"
	@echo "N_READ_PAIRS $(N_READ_PAIRS)"
	@echo "READS_WIN    $(READS_WIN)"
	wgsim -N $(N_READ_PAIRS) -1$(READ_LENGTH) -d0 -S11 -e0.050 -r0.0020 $< $@ /dev/null > /dev/null

mutated.fa: $(BED_EVENTS) normal.fa
	@echo ">> Creating new mutated reference"
	$(SRC)/introduce_events.py -e $(BED_EVENTS) -r normal.fa | ./$(SRC)/normalize_ref.py > $@

normal.fa:
	@echo ">> Preparing normal reference."
	([[ $(REF) == *.gz ]] && gzip -cd $(REF) > $@ || ln -s $(REF) $@)

wgsim:
	git clone https://github.com/lh3/wgsim.git
	cd wgsim; gcc -g -O2 -Wall -o wgsim wgsim.c -lz -lm
	mv wgsim/wgsim ./w
	rm -rf wgsim
	mv w wgsim

test:
	$(SRC)/introduce_events.py -e data/test.bed -r data/test.fa | ./$(SRC)/normalize_ref.py > _tmp
	@diff _tmp data/test.ok.fa && echo "test ok" || echo "test FAILED!"
	rm _tmp

rebuild_test:
	$(SRC)/introduce_events.py -e data/test.bed -r data/test.fa | ./$(SRC)/normalize_ref.py > data/test.ok.fa
	make test

clean:
	rm -f wgsim mutated* reads* normal*

.PHONY: clean test rebuild_test
.PRECIOUS: %.sa
