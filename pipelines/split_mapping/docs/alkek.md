### Alleviating pain by reusing old hardware

Here are some thoughts about an idea that has been going over my mind lately.

I have written a generic mapping pipeline from scratch and using it to map wgs
and wes samples in our ardmore hardware. In doing so I have experienced, once again, the pains associated with our HPC setup, namely: the poor performance and dubious interface of the PBS system and the bottlenecks caused by the network storage and the heterogeneity of the analysis.

The PBS software does not seem to react very well when the number of jobs in
the cluster increases. This gets accentuated when many dependencies exists among the jobs.

The interface against the PBS system is inferior to other solutions --particularly if we compare it against SGE or LSF-- and makes it a little bit more difficult for the developer. One particular annoyance is the way we have to express dependencies, using the job id instead of the job name. But this is something we can overcome with patience and experience. The behavior of the system under heavy load is a far more important and concerning.

We all enjoy the benefits of network storage: Huge amounts of storage available everywhere within the cluster. But this comes at a price: I/O bottlenecks. This is particularly evident in our environment where we have a decent number of users with different levels of programming and systems skills running analysis with different I/O profiles.

I believe building a "one-size-fits-all" solution --something we have been
trying to do for a while now-- is not possible. Setting constrains and mixing
environments is the only way to go. But that's another discussion.

In the rest of this document I am going to talk about an implementation which should help mitigate some of the issues outlined above. Some premises first:

I suspect a good part of the cpu cycles in the cluster are burned in mapping
activities. 

I think we have a pretty decent pull of hardware in Alkek not being
used right now.  When I say hardware I am talking about commodity rackable
machines, not storage. Probably those machines are not very powerful but we
should have plenty of them. They should have local drives of decent size (250Gb or more). We should also have more powerful machines with more local
disk (1T or more), but we just need a few of those (1 at least).

The idea is very simple: The sysadmins will load these  machines with a vanilla
Linux installation (ideally this process should be automatic) and give root access to the user (or users of a group) so they can use it any way they want.

One particular usage for this environment would be setting up a parallel mapping pipeline. The head node pulls bams from ardmore, splits them and performs the mapping on the rest of the nodes. Once completed the results are moved back to the head node where the merge and other pipeline steps are completed. Finally, the resulting bam can be sent back to ardmore. If properly designed, this setup can be easily reused by different groups --including production-- as they have picks of load.

With this setup the sysadmin involvement is minimal so they can focus in more
important tasks. I/O performance will be better and, most importantly, we will not encounter unexpected drops of performance since I/O is done against local disk.

Of course, code has to be written to implement all these, but once done it will be reused by multiple groups. Also, machines may fail, code will have to be written to be resilient against that.







