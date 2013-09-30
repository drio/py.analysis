# Given a list of gene names (GDOC), find its locations using a refgene file from
# uscs (URL_REFGENE) and liftover them to another genome (URL_LO) using a minimun
# edit distance of (MIN_MATCH)
#
THIS_MAKE := $(lastword $(MAKEFILE_LIST))
MAKE_DIR := $(dir $(THIS_MAKE))
SRC:=$(MAKE_DIR)/src
LO_BIN:="$(MAKE_DIR)/bin/liftOver"

GDOC="https://docs.google.com/spreadsheet/pub?key=0AlwwLefqWuS8dGp0REp5eHNlLU1mWUUzbmVfRUJuc0E&single=true&gid=0&output=txt"
URL_REFGENE="http://hgdownload.cse.ucsc.edu/goldenPath/rheMac3/database/refGene.txt.gz"
URL_LO="http://hgdownload.cse.ucsc.edu/goldenPath/rheMac3/liftOver/rheMac3ToHg19.over.chain.gz"
MIN_MATCH=0.9

all: stats.txt

# Report stats
stats.txt: lifted.bed
	@( \
		echo "n_genes n_genes_with_coor n_genes_lifted min_match"; \
		echo -n "`wc -l gene_list.txt | awk '{print $$1}'` "; \
  	echo -n "`wc -l genes.coor.txt | awk '{print $$1}'` "; \
		echo -n "`wc -l $< | awk '{print $$1}'` " ; \
		echo -n "$(MIN_MATCH)" ) > $@

# Liftover coordinates
lifted.bed: genes.coor.txt rheMac3ToHg19.over.chain.gz
	$(LO_BIN) -minMatch=$(MIN_MATCH) genes.coor.txt rheMac3ToHg19.over.chain.gz $@ unMapped

# Get genes coordinates from gene names
genes.coor.txt: refGene.txt.gz gene_list.txt
	gzip -cd refGene.txt.gz | $(SRC)/get_coordinates.py gene_list.txt > $@

rheMac3ToHg19.over.chain.gz:
	curl $(URL_LO) > $@

refGene.txt.gz:
	curl $(URL_REFGENE) > $@

gene_list.txt:
	curl $(GDOC) > $@
