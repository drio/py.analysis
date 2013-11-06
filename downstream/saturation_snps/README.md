## Computing the saturation of snps as a function of new samples added.

### What

This code should be used when we want to find how we saturate the number of
new snps found when we add more samples.

### Usage examples

1. saturation for substitutions:

```sh
$ gzip -cd *.gz | saturation.py | awk '{print $1" "$2}' | plot_sat.py "title: subs" "x: samples " "y: # of events" > subs.plot.py
```

2. gene saturation

```sh
$ gzip -cd *.gz | saturation.py | awk '{print $1" "$4}' | plot_sat.py "title: subs" "x: samples " "y: # of genes" > subs.plot.py
```

### Advance execution

The input stream (stdin) contains data from vcf files. The tool detects changes of samples using
the header of the vcf. If you understand this, you can refine the input stream to feed to the saturation
script anything you want. For example, here I am feeding only NON_SYNONYMOUS snps:

```sh
( export _FC="NON_SYNONYMOUS_CODING"; \
gzip -cd *.gz | \
  ruby -ne 'if $_[0] == "#" then puts $_ else m = $_.match(/;EFF=(\w+)/); puts $_ if m[1] == ENV["_FC"] end' \
) | saturation.py | awk '{print $1" "$4}' | grep -v num | plot_sat.py > output.png
```

As usual, we don't want to run this in the cluster so we can create a script and send it to the cluster:

```sh
$ cat run.sh
#!/bin/bash

gzip -cd /stornext/snfs6/rogers/drio_scratch/projects/rhesus/crv.v2/wgs/snps/*.vcf.anno.gz | /stornext/snfs6/rogers/drio_scratch/dev/py.analysis/downstream/saturation_snps/saturation.py | awk '{print $1" "$2}' | /stornext/snfs6/rogers/drio_scratch/dev/py.analysis/downstream/saturation_snps/plot_sat.py "title: subs" "x: samples " "y: # of events" > Ravi_subs.plot.png

$ submit -s saturation -m12G -c2 "./run.sh"
```

Notice that I am requesting quite a lot of memory. Use 16G to play save (for 55 wgs 30x samples).
The walltime time was around 10 hours.
