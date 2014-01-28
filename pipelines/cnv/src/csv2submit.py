#!/usr/bin/env python
#
# Given a csv/tab file with the list of jobs and its dependencies,
# generate target code to submit those to the cluster. I am
# assuming a tool exists (submit) for sending jobs to the cluster.
# Take a look to my implementation of that: https://github.com/drio/drd.bio.toolbox/blob/master/lua/submit.lua
# The reason for using this approach is to delegate the details of
# job submission for a particular cluster to the submit tool.
#
# Notice I am simplifying things here since I don't implement multiple
# dependencies. This is good for now as enough Parallelization levels
# are achieved with this simple approach.
#
import sys
import random
import drdcommon as common # https://github.com/drio/py.analysis/tree/master/drio.py

TOOL_NAME = "cmds2submit.py"
MAIN_HELP = """Usage: %s file.tsv\n""" % (TOOL_NAME)


def cmd2submit(lines, prev, sep="\t"):
    output_dep = "./deps.%s.txt" % int(random.random()*10000)
    i = 1
    for idx, l in enumerate(lines):
        _cmd, _name, _dep = l.split(sep)
        _name += "." + str(i)
        i += 1

        dep_input, dep_flag = '', ''
        if prev:
            dep_input = 'cat %s | ' % prev
            dep_flag = '-f -'

        redi = None
        if idx == 0:
            redi = '>'
        else:
            redi = '>>'

        print "%s submit %s -s %s -m %s -c %s '%s' | bash %s %s" % \
              (dep_input, dep_flag, _name, "6G", "1", _cmd, redi, output_dep)

    return output_dep


def main_help(msg=None, main_help=False):
    if msg:
        sys.stderr.write('ERROR: ' + msg + "\n")
    if main_help:
        sys.stderr.write(MAIN_HELP)


def main():
    if len(sys.argv) == 2:
        lines = []
        prev_dep_file = None
        for idx, line in enumerate(common.xopen(sys.argv[1])):
            if line in ['\n', '\r\n']:
                print ""
                prev_dep_file = cmd2submit(lines, prev_dep_file)
                lines = []
            else:
                lines.append(line.rstrip())
        print ""
        cmd2submit(lines, prev_dep_file)
    else:
        main_help('Need input file (use - for stdin).', main_help=True)


if __name__ == "__main__":
    main()
