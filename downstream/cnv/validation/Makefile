#
SRC=src
BED_EVENTS=data/test.bed
REF=data/phix.fa
N_READ_PAIRS=1000
READ_LENGTH=100

all: mutated.reads.fq.bam normal.reads.fq.bam

%.reads.fq.bam: %.reads.fq.sam
	samtools view -bS $< > $@

%.reads.fq.sam: %.fa.sa %.reads.fq
	bwa mem $(basename $<) $*.reads.fq > $@

%.fa.sa: %.fa
	bwa index $<

%.reads.fq: %.fa
	wgsim -N $(N_READ_PAIRS) -1$(READ_LENGTH) -d0 -S11 -e0.050 -r0.0020 $< $@ /dev/null

mutated.fa: $(BED_EVENTS) normal.fa
	$(SRC)/introduce_events.py -e $(BED_EVENTS) -r normal.fa | ./$(SRC)/normalize_ref.py > $@

normal.fa:
	ln -s $(REF) $@

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


