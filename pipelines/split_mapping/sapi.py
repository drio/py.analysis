#!/usr/bin/env python
"""
    sapi: split mapping pipeline
    TODO:
        1. add pbs
        2. add steps in -h

"""

import sys
import os
import argparse
from argparse import RawTextHelpFormatter
import glob
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


def schedulify(cmds, scheduler):
    out = ""
    for c in cmds:
        if scheduler == 'pbs':
            out += ("TODO")
        out += c + "\n"
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
                curr_dir=os.getcwd())
            return schedulify(cmd, self.args.scheduler)

    def fastqc(self, sdir, act, **kwargs):
        bam = kwargs["bam"]
        check_bam(bam)
        return ["%s/%s.sh %s" % (sdir, act, bam)]

    def splits(self, sdir, act, **kwargs):
        bam, nreads = kwargs["bam"], kwargs["nreads"]
        check_bam(bam)
        check_num_reads(nreads)
        return ["%s/%s.sh %s %s" % (sdir, act, bam, nreads)]

    def sais(self, sdir, act, **kwargs):
        fasta, n_threads = kwargs["fasta"], kwargs["n_threads"]
        curr_dir = kwargs["curr_dir"]
        check_done("splits", "done.txt")
        check_file(fasta, "Need fasta file.")
        cmd = []
        for f in glob.glob(curr_dir + "/splits/*.bam"):
            sp_num = match_this("\.(\d+)\.", f)
            for ot in [1, 2]:
                c = (("%s/%s.sh " + "%s " * 5) %
                    (sdir, act, fasta, ot, n_threads, f, sp_num))
                cmd.append(c)
        return cmd

    def sampe(self, sdir, act, **kwargs):
        fasta, bam = kwargs["fasta"], kwargs["bam"]
        curr_dir = kwargs["curr_dir"]
        check_bam(bam)
        check_file(curr_dir + "/splits/done.txt",
                   "Splits not completed. bailing out.")
        check_file(fasta, "Need fasta file.")
        actual_n_splits = int(open(curr_dir +
                              "/splits/done.txt").read().strip())
        if len(glob.glob("./sais/*.sai"))/2 != actual_n_splits:
            error("Number of splits does not match number of sai files")
        ones = glob.glob(curr_dir + "/sais/1.*.sai")
        twos = glob.glob(curr_dir + "/sais/2.*.sai")
        cmd = []
        for one, two in zip(ones, twos):
            sp_num = match_this("1.(\d+)\.", one)
            bam = glob.glob(curr_dir + "/splits/split.%s.bam" % sp_num)[0]
            c = (("%s/%s.sh " + "%s " * 5) %
                (sdir, act, fasta, one, two, sp_num, bam))
            cmd.append(c)
        return cmd

    def merge(self, sdir, act, **kwargs):
        curr_dir = kwargs["curr_dir"]
        actual_n_splits = int(open(curr_dir +
                              "/splits/done.txt").read().strip())
        if len(glob.glob("./sampe/*.bam")) != actual_n_splits:
            error("Number of sampe bams does not match number of splits")
        return ["%s/%s.sh" % (sdir, act)]

    def dups(self, sdir, act, **kwargs):
        if len(glob.glob("./merge/*.bam")) != 1:
            error("No merge bam found. Run the merge step first")
        return ["%s/%s.sh" % (sdir, act)]

    def stats(self, sdir, act, **kwargs):
        if len(glob.glob("./dups/*.bam")) != 1:
            error("No dups bam found. Run the dups step first")
        return ["%s/%s.sh" % (sdir, act)]


def list_steps():
    methods = set({})
    for name, _ in inspect.getmembers(Action, predicate=inspect.ismethod):
        methods.add(name)

    s = "Available steps:\n"
    for name in methods.difference({"__init__", "cmds"}):
        s = s + name + "\n"
    return s


def process_args():
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
    parser.add_argument('--scheduler', dest='scheduler', action='store',
                        choices={'pbs', 'single'}, default='single',
                        help='Specify what scheduler to use')
    return parser.parse_args()


def do_work():
    args = process_args()
    print Action(args).cmds()

if __name__ == "__main__":
    do_work()
