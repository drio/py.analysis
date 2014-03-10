bams.bed:
	curl -s "https://docs.google.com/spreadsheet/pub?key=0Ao0CXEbTWQu2dHJTcUxIOEFtNy1qWGlXQ0Z4RGlOLVE&single=true&gid=0&output=txt" | awk -F"\t" '{print $11"\t"$2}' > bams.bed

coverage:
	for b in *.bam;do submit -s cov.`basename $$b` ./extract_cov.sh $$b ; done

links: bams.bed
	cat bams.bed | awk -F"\t" '{print "ln -s "$$1" "$$2".bam"}' | grep -v DNA | bash

kill:
	qstat -a | grep "cov." | awk '{print $$1}' |xargs -I{} bkill {}

.PHONY: links coverage kill
