#!/usr/bin/env python

import sys
import re
import random
import drdcommon as common

TOOL_NAME = "cmds2submit.py"
MAIN_HELP = """Usage: %s file.tsv\n""" % (TOOL_NAME)


def cmd2submit(lines, prev):
    output_dep = "./deps.%s.txt" % int(random.random()*10000)
    for idx, l in enumerate(lines):
        name, _cmd, cores, mem = l.split(",")

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
              (dep_input, dep_flag, name, mem, cores, _cmd, redi, output_dep)

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
                print "#--------"
                prev_dep_file = cmd2submit(lines, prev_dep_file)
                lines = []
            else:
                lines.append(line.rstrip())
    else:
        main_help('Need input file (use - for stdin).', main_help=True)

if __name__ == "__main__":
    main()
