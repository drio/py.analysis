#
SRC=src
BED_EVENTS=data/small_events.bed
REF=data/phix.fa
N_READ_PAIRS=1000
READ_LENGTH=100

reads.fq.bam: reads.fq.sam
	samtools view -bS $< > $@

reads.fq.sam: reads.fq mutated_ref.fa.sa
	bwa mem mutated_ref.fa $< > $@

mutated_ref.fa.sa: mutated_ref.fa
	bwa index $<

reads.fq: wgsim mutated_ref.fa
	wgsim -N $(N_READ_PAIRS) -1$(READ_LENGTH) -d0 -S11 -e0.050 -r0.0020 mutated_ref.fa $@ /dev/null

mutated_ref.fa: $(BED_EVENTS) $(REF)
	$(SRC)/introduce_events.py -e $(BED_EVENTS) -r $(REF) | ./$(SRC)/normalize_ref.py > $@

wgsim:
	git clone https://github.com/lh3/wgsim.git
	cd wgsim; gcc -g -O2 -Wall -o wgsim wgsim.c -lz -lm
	mv wgsim/wgsim ./w
	rm -rf wgsim
	mv w wgsim

clean:
	rm -f wgsim mutated* reads*

.PHONY: clean $(BIN_DEPS)
