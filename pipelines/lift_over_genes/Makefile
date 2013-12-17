# Given a list of gene names (GDOC), find its locations using a gtf file from ensembl
# (GTF) and liftover them to another genome (CHAIN) using a minimun
# edit distance of (MIN_MATCH)
#
SHELL=/bin/bash
THIS_MAKE := $(lastword $(MAKEFILE_LIST))
MAKE_DIR := $(dir $(THIS_MAKE))
SRC:=$(MAKE_DIR)/src
URL_LO="http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver"
LO_BIN:=$(shell echo $(URL_LO) | awk -F/ '{print $$7}')

GDOC="https://docs.google.com/spreadsheet/pub?key=0AlwwLefqWuS8dGp0REp5eHNlLU1mWUUzbmVfRUJuc0E&single=true&gid=0&output=txt"
GTF="ftp://ftp.ensembl.org/pub/release-73/gtf/macaca_mulatta/Macaca_mulatta.MMUL_1.73.gtf.gz"
CHAIN:="http://hgdownload.cse.ucsc.edu/goldenPath/rheMac2/liftOver/rheMac2ToHg19.over.chain.gz"
CHAIN_BASE:=$(shell echo $(CHAIN) | awk -F/ '{print $$7}')
MIN_MATCH=0.9

all: stats.txt

# Report stats
stats.txt: lifted.bed
	@( \
		echo "n_genes n_exons_with_coor n_genes_lifted min_match"; \
		echo -n "`wc -l gene_list.txt | awk '{print $$1}'` "; \
  	echo -n "`wc -l exon.coor.bed | awk '{print $$1}'` "; \
		echo -n "`wc -l $< | awk '{print $$1}'` " ; \
		echo -n "$(MIN_MATCH)" ) > $@

# Liftover coordinates
lifted.bed: exon.coor.bed $(CHAIN_BASE) $(LO_BIN)
	./$(LO_BIN) -tab -bedPlus=3 minMatch=$(MIN_MATCH) exon.coor.bed $(CHAIN_BASE) $@ unMapped

# Get genes coordinates from gene names
exon.coor.bed: gtf.gz gene_list.txt
	gzip -cd gtf.gz | $(SRC)/get_coordinates.py gene_list.txt | sort -k4,4 -k6,6 -k5,5n > $@

$(CHAIN_BASE):
	curl $(CHAIN) > $@

gtf.gz:
	curl $(GTF) > $@

gene_list.txt:
	curl $(GDOC) > $@

$(LO_BIN):
	curl $(URL_LO) > $@
	chmod 755 $@
