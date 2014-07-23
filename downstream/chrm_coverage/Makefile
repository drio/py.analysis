#
SRC := $(dir $(shell readlink -f Makefile))
BAM := input.bam
BT  := /hgsc_software/BEDTools/latest/bin/bedtools

plot.%.png: coverage.bed.gz
	gzip -cd $< | grep -P "^$*\t" | $(SRC)/plot.py > $@

coverage.bed.gz: genome.bed
	$(BT) genomecov -g $< -ibam $(BAM) -d | gzip -c > $@

genome.bed: $(BAM)
	samtools view -H $< $(BAM) | grep SQ | awk -F: '{print $$2"1\t"$$3}' | sed 's/LN//g' > $@


.PHONY: coverage
