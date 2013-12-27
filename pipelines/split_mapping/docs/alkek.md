### Alleviating pain by reusing old hardware

I would like to share some ideas about efficient usage of hardware.

I have written a generic mapping pipeline from scratch which I use to map WGS and WES samples using our ardmore hardware. In doing so, once again I have experienced the pains associated with our HPC setup: the poor performance and dubious interface of the PBS system and the bottlenecks caused by the network storage and the analysis heterogeneity. Moreover, the PBS software does not seem to react very well when under heavy load (high number of jobs). This gets accentuated when many dependencies exists among the jobs.

The PBS interface is inferior to other solutions, e.g. SGE or LSF. It makes it more difficult for the developer. One particular painpoint is that we have to express dependencies between jobs using jobs IDs instead of the job names. This is not a major problem as this is something we can overcome with patience and experience. However, the behavior of the system under heavy load is a far more important and concerning.

We all enjoy the benefits of network storage such as huge amounts of storage available anywhere within the cluster. But this comes at a price of I/O bottlenecks. This is particularly evident in our environment where we have a diverse set of users with different levels of programming and systems skills running analysis with different I/O profiles.

I believe building a "one-size-fits-all" solution --something we have been trying to do for a while now, is not possible. Setting constraints and mixing environments is the only way to go. However, that is not the objective of this discussion.

In the rest of this document I am going to talk about an implementation which should help mitigate some of the issues outlined above. 

Some premises first:

1) I suspect a good part of the cpu cycles in the cluster are burned in mapping activities. 

2) I think we have a pretty decent pull of hardware in Alkek not being used right now. By hardware I mean commodity rackable machines, not storage. Probably those machines are not very powerful but we should have plenty of them. They should have local drives of decent size (250Gb or more). We should also have more powerful machines with more local disk (1T or more), but we just need a few of those (at least one).

My proposal:

The idea is very simple: The sysadmins will load these Alkek machines with a vanilla Linux installation (ideally this process should be automatic) and give root access to the user (or users of a group) so they can use it any way they want.

One concrete usage for this environment would be setting up a parallel mapping pipeline. The head node pulls bams from ardmore, splits them and performs the mapping on the rest of the nodes. Once completed the results are moved back to the head node where the merge and other pipeline steps are completed. Finally, the resulting bam can be sent back to ardmore. If properly designed, this setup can be easily reused by different groups --including production-- as they have picks of load.

With this setup the sysadmin involvement is minimal so they can focus on more important tasks. I/O performance will be better and, most importantly, we will not encounter unexpected drops of performance since I/O is done against local disk.

This approach has a great advantage. As I said earlier, much of the CPU cycles are spent in mapping the reads, moving mapping to Alkek substantially frees up the nodes on Ardmore. Moreover, if mapping pipeline can be made parallel (e.g. split the reads and use multiple machines per sample for mapping), the mapping jobs will get completed much faster. Thus, we might be able to reduce the overall load on the Ardmore infrastructure.

It is obvious that this approach requires writing additional code to implement the functionality I discussed. However, once written, it can be used by multiple groups. One more aspect to consider is machine failures. The code must be designed to be resilient against such machine failures.








