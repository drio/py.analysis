# This is to generate a small genome (for testing) from rheMac2
#
url:="http://hgdownload.cse.ucsc.edu/goldenPath/rheMac2/bigZips/chromFa.tar.gz"
genome:=rhemac2
nbp:=$(shell echo "100000" | bc)

all: small.rhemac2.$(nbp).fa small.rhemac2.$(nbp).sizes.bed

chromFa.tar.gz:
	curl http://hgdownload.cse.ucsc.edu/goldenPath/rheMac2/bigZips/chromFa.tar.gz > $@

softMask: chromFa.tar.gz
	tar zxvf $<

%.fa %.bed:
	./gen_smaller_genome_for_testing.sh $(genome) $(nbp)
