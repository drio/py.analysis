#
SHELL=/bin/bash
THIS_MAKE := $(lastword $(MAKEFILE_LIST))
MAKE_DIR := $(dir $(THIS_MAKE))
SRC:=$(MAKE_DIR)/src
MF:=$(shell pwd)/bin/mrsfast

MRSF_URL="http://downloads.sourceforge.net/project/mrsfast/2.6.0.4/mrsfast-2.6.0.4.zip?r=http%3A%2F%2Fmrsfast.sourceforge.net%2FHome&ts=1388504359&use_mirror=hivelocity"

all: bin/mrsfast

bin/mrsfast: mrsf.tar.gz
	mkdir bin
	(cd bin;\
	mv ../*.tar.gz . ;\
	tar zxvf mrsf.tar.gz;\
	cd mrsfast-2.6.0.4;\
	cat Makefile | grep -v static > m;\
	mv m Makefile;\
	make;\
	mv mrsfast ../;\
	cd ..;  rm -rf mrsfast-* *.gz;)

mrsf.tar.gz:
	curl -L $(MRSF_URL) > $@

sample_genome:
	curl -L "http://downloads.sourceforge.net/project/mrfast/mrfast/sample_genome.tar.gz?r=http%3A%2F%2Fmrfast.sourceforge.net%2F&ts=1388506275&use_mirror=hivelocity" | tar -zxvf -

sample_genome_two: sample_genome
	(mkdir $@;\
		cd $@;\
		cp -r ../sample_genome/*.fasta .;\
		cp -r ../sample_genome/reads*  .;\
		$(MF) --index genome.fasta;\
		samtools faidx genome.fasta;\
	  $(MF) --search genome.fasta --seq reads1.sim -e 2 -o test_out.sam -u test_unmapped.fastq)
	  #$(MF) --search genome.fasta --seq1 reads1.sim --seq2 reads2.sim -e 5 --min 100 --max 400 --discordant-vh --pe -o test_out.sam -u test_unmapped.fastq)
		# -e: num of mismatches
		# -min/max: distances between ends
		# -discordant-vh: report discordant pairs
	  #$(MF) --search genome.fasta --seq1 reads1.sim --seq2 reads2.sim -e 3 --min 300 --max 400 --discordant-vh --pe -o test_out.sam -u test_unmapped.fastq)

paper:
	curl -L http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3115707/pdf/nihms297017.pdf > paper.nihms297017.pdf

clean:
	rm -rf bin *.pdf
