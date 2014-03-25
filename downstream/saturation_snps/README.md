## Computing the saturation of snps as a function of new samples added.

### What

This code should be used when we want to find how we saturate the number of
new snps found when we add more samples.

### Usage examples

1. saturation for substitutions:

```sh
$ gzip -cd *.gz | saturation.py | tee data.txt | grep -v samples | awk '{print $1" "$2}' | plot_sat.py "title: subs" "x: samples " "y: # of events" > subs.plot.py
```

2. gene saturation

```sh
$ gzip -cd *.gz | saturation.py | tee data.txt | grep -v samples | awk '{print $1" "$4}' | plot_sat.py "title: subs" "x: samples " "y: # of genes" > subs.plot.py
```

### Advance execution

The input stream (stdin) contains data from vcf files. The tool detects changes of samples using
the header of the vcf. If you understand this, you can refine the input stream to feed to the saturation
script anything you want. For example, here I am feeding only NON_SYNONYMOUS snps:

```sh
( export _FC="NON_SYNONYMOUS_CODING"; \
gzip -cd *.gz | \
  ruby -ne 'if $_[0] == "#" then puts $_ else m = $_.match(/;EFF=(\w+)/); puts $_ if m[1] == ENV["_FC"] end' \
) | saturation.py | tee data.txt |grep -v samples | awk '{print $1" "$4}' | grep -v num | plot_sat.py > output.png
```

As usual, we don't want to run this in the cluster so we can create a script and send it to the cluster:

```sh
$ cat run.sh
#!/bin/bash

gzip -cd /stornext/snfs6/rogers/drio_scratch/projects/rhesus/crv.v2/wgs/snps/*.vcf.anno.gz | \
/stornext/snfs6/rogers/drio_scratch/dev/py.analysis/downstream/saturation_snps/saturation.py | \
tee data.txt | \
grep -v samples | \
awk '{print $1" "$2}' | /stornext/snfs6/rogers/drio_scratch/dev/py.analysis/downstream/saturation_snps/plot_sat.py "title: subs" "x: samples " "y: # of events" > Ravi_subs.plot.png

$ submit -s saturation -m12G -c2 "./run.sh"
```

Notice that I am requesting quite a lot of memory. Use 16G to play save (for 55 wgs 30x samples).
The walltime time was around 10 hours.

#### Using the `-a <num>` parameter to specify when to count a snp

You can use the `-a <num>` parameter to tell the saturation tool in how many samples
we have to see the snp in order to start counting it. By default this value is 2. Use `-a 1`
if you want a snp to be counted as soon as is seen in one animal.

#### How can I control the order in which the samples are computed?

In all the previous sections we have been using `gzip -cd *.gz` to generate the stream that we feed
into the saturation tool. Let's imagine you have two subgroups of samples and you want to feed first
one of the groups to the tool. Here is one possible approach to acomplish this:

```sh
$ ls group1
file1.vcf.gz file2.vcf.gz
$ ls group2
file1.vcf.gz file2.vcf.gz
$ (gzip -cd group1/*.gz; gzip -cd group2/*.gz) | saturation.py | tee data.txt | grep -v samples | awk '{print $1" "$2}' | plot_sat.py "title: subs" "x: samples " "y: # of events" > subs.plot.py
```

### Other questions:


#### Just to understand the curve well (There are some sudden bumps in the curve and want to know which specific samples through this) is it possible to include the sample names either in the sat.e or data.txt?

We don't have the sample names available when the script is running, just the
vcf data. So generated the sample id cannot be done at that level. That being
said, the order in which the vcfs are fed into the script is determined by the
shell. In your script your are using:

`cat /stornext/snfs5/rogers/Ravi/saturation/freez1/Indian/30XWGS/*.recode.vcf`

To find the actual files (and therefore the samples ids), you can just run a
previous ls command with the same argument and store it to a file:

`ls /stornext/snfs5/rogers/Ravi/saturation/freez1/Indian/30XWGS/*.recode.vcf > order.txt`

So in order.txt you have now the actual files and the line number is the order
in which they are run.


#### Can I check the saturation for genes instead of snps?

Yes, if your vcf files are annotated with snp effect. Example:

```sh
$ src="/stornext/snfs6/rogers/drio_scratch/dev/py.analysis/downstream/saturation_snps"
$ gzip -cd *.anno*.gz  | $src/gene_saturation_snpeff.sh | $src/plot_sat.py "title: subs" "x: samples " "y: # of events" > plot.png
```



