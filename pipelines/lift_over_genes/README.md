## Lifting over genomic regions between genomes

### Description

This pipeline will help you asses how a set of genomic regions are conserved
between genomes.

### Usage

```sh
$ mkdir foo; cd foo
$ makefile -f PATH_TO_MAKEFILE
```

A more realistic example will be:

```sh
$ makefile -f MAKE GDOC="https://docs.google.com/spreadsheet/pub?key=0AlwwL" \
GTF="ftp://ftp.ensembl.org/pub/release-73/gtf/macaca_mulatta/Macaca_mulatta.MMUL_1.73.gtf.gz" \
URL_LO="http://hgdownload.cse.ucsc.edu/goldenPath/rheMac3/liftOver/rheMac3ToHg19.over.chain.gz" \
```

```GDOC``` points to the list of genes.
```URL_REFGENE``` is the refGene file (ucsc) to use to extract the gene coordinates.
```URL_LO```  points to the lift over database.


