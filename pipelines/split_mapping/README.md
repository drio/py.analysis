## Split mapping pipeline

### What is this?

The amount of reads generated in a single illumina lane is increasing constantly. Mapping all those reads is quite straight forward thanks to aligners like BWA, doing so in an efficient way (both time and space) is more complicated.

### Mapping steps

1. fastqc
2. sai generation (finding alignment locations for the reads)
3. to sam (Computing local alignment if necessary and generate SAM records)
4. sam to bam (+ coordinate sorting)
5. mark duplicates
6. generate stats
7. clean up

### Approaches: A high level view

#### Approach 1: The naive, simple way.

One possible first approach is to compute steps 1 and 2 for reads 1 and reads 2. Then continue sequentially with the other steps. To accomplish this we need to leverage the job dependency facilities that your cluster scheduler provides.

Depending the amounts of data you have and how fast you want the alignments this approach may be fine. Maybe you don't mind waiting some extra days to get your results. That makes sense if you run these pipeline from time to time and you don't have tied deadlines to return the results

But remember, we are talking about 1 Billion reads for a 30x whole genome experiment. Processing 500M reads in a decent 8 core machine can take between 5 and 7 days.

#### Approach 2: Split reads

We can do better by exploiting the embarrassing parallel nature of the problem. We could split the reads in multiple chunks and compute step 2 in parallel. 

The obvious benefit is the reduction of running time. Also, we could potentially avoid recomputations if a split has already been processed. That would add up complexity to the pipeline but it maybe worth it for your particular situation.

#### What do I need?

This is just a brainstorming on my end to think about what problem we have in my group and what's the best way to solve it.

Approach 1 is no viable for various reasons: First, we have tons of data coming. Second, our cluster is taken down for maintenance once a month (don't ask), that means we could loose a lot of computations if step 2 does not complete prior to the downtime.

#### The new pipeline steps: 

1. fastqc (1)
2. sai generation (N)
3. to sam (N)
4. clean up sais
5. sam to bam (1)
6. clean up sams
7. mark duplicates (1)
8. clean up sorted bam
9. generate stats (1)

### Implementation

This is the fun part.

We have to make an strategic decision here: do we control dependencies and fire up jobs automatically as previous steps are done? You are probably thinking this should be mandatory. This is where the cluster reliability comes into play. I cannot assume my cluster's software will properly signal job completions. Yes, it shouldn't be like that but it is: If the cluster gets flooded with jobs it will stop controlling job decencies so your pipeline could get stuck at any point. I have heart some schedulers handle this much better (LSF and SGE) compared to others (PBS/moab). 

So based on these, I am not going to bother controlling dependencies. I will let the user take care of that if they want to. Good luck with that.

Now that we have made some decisions, let's add some more assumptions: We will have a bam as input. 

So, how would the execution of the pipeline look like?

```
$ sapi 
Available actions:
    help
    version
    fastqc
    sai
    sam
    rm_sais
    sam2bam
    ....

$ sapi fastqc
java -jar path.jar -xmx4G input.bam

$ sapi split -s 10 input.bam
split -s 10 input.bam

$ sapi sai
bwa sai -1 input.1.bam -2 input.1.bam > sai.1.1.sai
bwa sai -1 input.2.bam -2 input.2.bam > sai.1.2.sai
....

$ sapi sam
bwa sampe -1 sai.1.1.sai -2 sai.1.2.sai > output.1.sam
....

$ sapi rmsai
rm -f *.sai

$ sapi merge
java -jar xxxx -xmx12G .... INPUT OUTPUT:merged.sorted.bam
```

I hope you get the point.

How do we send this to the cluster?

```
$ sapi sam | to_cluster 
sapi merge | \
qsub 
-N one 
-q analysis 
-d `pwd` 
-o moab_logs/one.o 
-e moab_logs/one.e 
-l nodes=1:ppn=1,mem=4Gb -V
```

This is very flexible but we have to tell ```to_cluster``` a bunch of options that will change per each job:

1. job name
2. queue
3. output log
4. stderr log
5. resources
6. others?

In order to avoid coupling between sapi and to_cluster, we could pass those options in the stderr? but what if we have an actual problem ? Ok, it would be better to have text files with the options required per each step. And we could then remove the to_cluste all together:

```
$ sapi merge
java -jar XXXX .... 
$ sapi merge -c pbs
echo 'java -jar XXXXX' | qsub -c XXXX ..... 
```

By default, we generate just the shell command, if you use the -c parameter (and provide the scheduler type as value) you 
will get the shell command ready to be submitted to your cluster.

Ok.. hacking time... let's build a prototype.



    
