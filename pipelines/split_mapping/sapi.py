#!/usr/bin/env python
"""
    sapi: split mapping pipeline
    TODO:
        1. add pbs
        2. add steps in -h
        3. remove unnecessary files

"""

import sys
import os
import argparse
from argparse import RawTextHelpFormatter
import re
import inspect


def error(msg):
    sys.stderr.write(msg + "\n")
    exit(1)


def match_this(reg_exp, string):
    match = re.search(reg_exp, string)
    if match:
        return match.group(1)
    else:
        error("Failing matching re (" + reg_exp + ") against" + string)


def check_bam(bam):
    if not bam:
        error("Need input bam.")


def check_num_reads(nreads):
    if not nreads:
        error("Need the number of reads you want to put per split.")


def check_done(prev_step, _file):
    if not os.path.isfile(prev_step + "/" + _file):
        error("Previous step (%s) hasn't completed yet." % prev_step)


def check_file(_file, msg):
    if not _file:
        error(msg)
    if not os.path.lexists(_file):
        error(msg + " (" + _file + ")")


def gen_job_name(sample_id, step, index):
    return "%s_%s_%s" % (sample_id, step, index)


def extra_mem(o_mem):
    if o_mem[-1] == 'g' or o_mem[-1] == 'G':
        num = int(match_this("^(\d+)", o_mem))
        num += 2
        return "%sg" % num
    else:
        return o_mem


def cmd_to_pbs(cmd, sample_id, queue, step, index, mem, cores, tmp):
    t = "mkdir -p logs; echo '_CMD_' | qsub -N _NAME_ -q _QUEUE_ -d `pwd` "
    t += "-o logs/_NAME_.o -e logs/_NAME_.e "
    t += "-l nodes=1:ppn=_CORES_,mem=_MEM_ -V"
    t = t.replace('_QUEUE_', queue)
    t = t.replace('_NAME_', gen_job_name(sample_id, step, index))
    t = t.replace('_CORES_', cores)
    # The JVM uses more memory than specified by the -X parameters (check
    # (http://blogs.vmware.com/apps/2011/06/taking-a-closer-look-at-sizing-the-java-process.html) for
    # more details. If the cluster software has hard memory limits jobs may
    # get killed. To avoid this, we request some extra memory when sending
    # the job to the cluster.
    t = t.replace('_MEM_', extra_mem(mem))
    t = t.replace('_CMD_', cmd)
    return t


def cmd_to_csv(cmd, sample_id, queue, step, index, mem, cores, tmp):
    t = "_NAME_,_CMD_,_CORES_,_MEM_"
    t = t.replace('_CMD_', cmd)
    t = t.replace('_NAME_', gen_job_name(sample_id, step, index))
    t = t.replace('_CORES_', cores)
    t = t.replace('_MEM_', mem)
    return t


def schedulify(cmds, scheduler, sample_id, queue, step, mem, cores, tmp):
    """ If user provides a scheduler, we have to wrap the cmds with
        the appropiate scheduler command
    """
    out = ""
    for idx, c in enumerate(cmds):
        if scheduler == 'pbs':
            out += cmd_to_pbs(c, sample_id, queue, step, idx, mem, cores, tmp)
        elif scheduler == 'csv':
            out += cmd_to_csv(c, sample_id, queue, step, idx, mem, cores, tmp)
        else:
            out += c
        out += "\n"
    return out


class Action(object):
    def __init__(self, args):
        self.args = args
        self.scripts_dir = os.path.dirname(os.path.realpath(__file__)) + "/sh"

    def cmds(self):
        # take care of the action
        try:
            getattr(self, self.args.step)
        except:
            error("[%s] is not a valid step." % self.args.step)
        else:
            cmd = getattr(self, self.args.step)(
                self.scripts_dir, self.args.step,
                bam=self.args.bam, nreads=self.args.num_reads,
                fasta=self.args.fasta, n_threads=self.args.n_threads,
                curr_dir=os.getcwd(), tmp=self.args.tmp,
                mem=self.args.mem, sample_id=self.args.sample_id,
                url=self.args.url)
            return schedulify(cmd, self.args.scheduler,
                              self.args.sample_id, self.args.queue,
                              self.args.step,
                              self.args.mem, self.args.n_threads,
                              self.args.tmp)

    def init(self, sdir, act, **kwargs):
        bam, nreads = kwargs["bam"], kwargs["nreads"]
        url = kwargs["url"]
        _id = kwargs["sample_id"]
        if not nreads:
            error("Need number of reads per split.")
        check_bam(bam)
        return ["%s/%s.sh %s %s %s %s" % (sdir, act, bam, nreads, url, _id)]

    def validate(self, sdir, act, **kwargs):
        tmp_dir = kwargs["tmp"]
        mem = kwargs["mem"]
        bam = kwargs["bam"]
        url = kwargs["url"]
        _id = kwargs["sample_id"]
        check_bam(bam)
        return ["%s/%s.sh %s %s %s %s %s" % (sdir, act, tmp_dir, mem, bam, url, _id)]

    def fastqc(self, sdir, act, **kwargs):
        bam = kwargs["bam"]
        url = kwargs["url"]
        _id = kwargs["sample_id"]
        check_bam(bam)
        return ["%s/%s.sh %s %s %s" % (sdir, act, bam, url, _id)]

    def splits(self, sdir, act, **kwargs):
        tmp_dir = kwargs["tmp"]
        mem = kwargs["mem"]
        url = kwargs["url"]
        bam, nreads = kwargs["bam"], kwargs["nreads"]
        _id = kwargs["sample_id"]
        check_bam(bam)
        check_num_reads(nreads)
        return ["%s/%s.sh %s %s %s %s %s %s" % (sdir, act, bam,
                                          nreads, mem, tmp_dir, url, _id)]

    def sais(self, sdir, act, **kwargs):
        fasta, n_threads = kwargs["fasta"], kwargs["n_threads"]
        check_file(fasta, "Need fasta file.")
        _id = kwargs["sample_id"]
        cmd = []
        url = kwargs["url"]
        n_of_splits = int(open('init/nsplits.txt').read().rstrip())
        for sp_num in range(0, n_of_splits):
            for ot in [1, 2]:
                f = "../splits/split.%2.2d.bam" % sp_num
                c = (("%s/%s.sh " + "%s " * 7) %
                    (sdir, act, fasta, ot, n_threads, f, "%2.2d" % sp_num, url, _id))
                cmd.append(c)
        return cmd

    def sampe(self, sdir, act, **kwargs):
        fasta, bam = kwargs["fasta"], kwargs["bam"]
        curr_dir = kwargs["curr_dir"]
        tmp_dir = kwargs["tmp"]
        _id = kwargs["sample_id"]
        url = kwargs["url"]
        mem = kwargs["mem"]
        check_bam(bam)
        check_file(fasta, "Need fasta file.")
        cmd = []
        n_of_splits = int(open('init/nsplits.txt').read().rstrip())
        for sp_num in range(0, n_of_splits):
            bam = curr_dir + "/splits/split.%2.2d.bam" % sp_num
            one = "../sais/1.%2.2d.sai" % sp_num
            two = "../sais/2.%2.2d.sai" % sp_num
            c = (("%s/%s.sh " + "%s " * 9) %
                (sdir, act, fasta, one, two, sp_num, bam, mem, tmp_dir, url, _id))
            cmd.append(c)
        return cmd

    def merge(self, sdir, act, **kwargs):
        url = kwargs["url"]
        fasta, bam = kwargs["fasta"], kwargs["bam"]
        tmp_dir = kwargs["tmp"]
        mem = kwargs["mem"]
        _id = kwargs["sample_id"]
        check_bam(bam)
        check_file(fasta, "Need fasta file.")
        return ["%s/%s.sh %s %s %s %s %s %s" % (sdir, act, tmp_dir, mem, bam, fasta, url, _id)]

    def dups(self, sdir, act, **kwargs):
        tmp_dir = kwargs["tmp"]
        mem = kwargs["mem"]
        url = kwargs["url"]
        _id = kwargs["sample_id"]
        return ["%s/%s.sh %s %s %s %s" % (sdir, act, tmp_dir, mem, _id, url)]

    def stats(self, sdir, act, **kwargs):
        url = kwargs["url"]
        sample_id = kwargs["sample_id"]
        return ["%s/%s.sh %s %s" % (sdir, act, url, sample_id)]


def list_steps():
    order = ['fastqc', 'validate', 'splits', 'sais',
             'sampe', 'merge', 'dups', 'stats']
    desc = {
        'fastqc':   "Run fastQC on bam",
        'validate': "Validate a bam.",
        'init':     "Perform various computations necessary in other steps",
        'splits':   "Split a bam so it has <n_reads> per split",
        'sais':     "Find candidate aligment locations for reads in splitted input",
        'sampe':    "Generate sam records from sai files",
        'merge':    "Merge splitted bams",
        'dups':     "Mark duplicates",
        'stats':    "Compute some metrics for given bam",
    }
    methods = set({})
    for name, _ in inspect.getmembers(Action, predicate=inspect.ismethod):
        methods.add(name)

    s = "+ Available steps:\n"
    for name in order:
        s = s + "* " + name + ": " + desc[name] + "\n"

    s += "\n"
    s += "Use check_env to verify your environment"

    return s


def check_env():
    def which(program):
        def is_exe(fpath):
            return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

        fpath, fname = os.path.split(program)
        if fpath:
            if is_exe(program):
                return program
        else:
            for path in os.environ["PATH"].split(os.pathsep):
                path = path.strip('"')
                exe_file = os.path.join(path, program)
                if is_exe(exe_file):
                    return exe_file

        return None

    for p in ["bwa", "samtools", "bash", "fastqc", "sapi.py", "signal.py"]:
        where = which(p)
        msg = "found" if where else "!! NOT found"
        print "[%s] %s (%s)" % (msg, p, where)


def process_args():
    if len(sys.argv) == 2 and sys.argv[1] == "check_env":
        check_env()
        exit(0)

    parser = argparse.ArgumentParser(description="Split mapping pipeline",
                                     formatter_class=RawTextHelpFormatter)

    parser.add_argument(dest='step', action='store', default='single',
                        help="Step you want to compute\n" + list_steps())

    parser.add_argument('-n', dest='num_reads', action='store',
                        help='number of reads per split')
    parser.add_argument('-b', dest='bam', action='store',
                        help='input bam file.')
    parser.add_argument('-f', dest='fasta', action='store',
                        help='input fasta file')

    parser.add_argument('-t', dest='n_threads', action='store', default='1',
                        help='Number of threads')
    parser.add_argument('-m', dest='mem', action='store', default='8g',
                        help='Amount of mem to use (Use g for gibabytes)')

    parser.add_argument('-q', dest='queue', action='store',
                        help='Cluster queue')
    parser.add_argument('-i', dest='sample_id', required=True,
                        action='store', help='sample id')

    parser.add_argument('-x', dest='execute', action='store_true', default=False,
                        required=False,
                        help='Do not show the cmd, execute it')

    parser.add_argument('-u', dest='url', action='store', required=False,
                        help='URL for the rest service.')

    parser.add_argument('-c', dest='check_env', action='store_true',
                        help='Check binary dependencies')

    tmp_default = "/tmp"
    if os.path.isdir("/space1/tmp"):
        tmp_default = "/space1/tmp"
    parser.add_argument('-p', dest='tmp', action='store', default=tmp_default,
                        help='temporary directory to use')

    parser.add_argument('-s', dest='scheduler', action='store',
                        choices={'pbs', 'single', 'csv'}, default='single',
                        help='Specify what scheduler to use')

    p = parser.parse_args()

    if p.bam and not os.path.isabs(p.bam):
        p.bam = os.path.realpath(p.bam)
    if p.fasta:
        p.fasta = os.path.realpath(p.fasta)
    if p.scheduler == 'pbs' and not(p.queue):
        error("What queue do you want me to use?")

    return p


def do_work():
    args = process_args()
    cmd = Action(args).cmds()
    if args.execute:
        os.system(cmd)
    else:
        print cmd

if __name__ == "__main__":
    do_work()
