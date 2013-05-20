#!/usr/bin/env python

import sys, os, tempfile
import pandas as pd
from random import randint
__lib_dir = os.path.dirname( os.path.realpath(__file__) ) + "/../drio.py"
sys.path.append(__lib_dir)
import drdcommon as common
import ngsproject as np

TOOL_NAME = "cmds2submit.py"
MAIN_HELP = """Usage: %s file.tsv\n""" % (TOOL_NAME)

class Job(object):
  def __init__(self, series, dep_fn, first=False):
    #self.cmd="submit -s fmi.{} -m 8G -c 8 'sga index --no-reverse -d 5000000 -t 8 {}'
    self.dep_fn = dep_fn
    self.series = series
    self.first  = first

  def __str__(self):
    s = self.series
    cmd = "submit -s %s -m %s -c %s " % (s["id"], s["ram"], s["threads"])
    redirect = ">" if self.first else ">>"
    if s["deps"] != "-":
      cmd += "-f %s " % self.dep_fn
    cmd += "'%s' " % s["cmd"]
    cmd += " | bash "
    cmd += "%s %s " % (redirect, self.dep_fn)
    return cmd

def main_help(msg=None, main_help=False):
  if msg: sys.stderr.write('ERROR: ' + msg + "\n")
  if main_help: sys.stderr.write(MAIN_HELP)

def main():
  dep_fn = "deps." + str(randint(1,1000000))
  if len(sys.argv) == 2:
    # Dirty hack since I don't know how to make pandas.read_table work
    # off of stdin.
    data = common.xopen(sys.argv[1]).read()
    f = tempfile.NamedTemporaryFile(delete=False)
    f.write(data)
    f.close()
    for i,s in pd.read_table(f.name).iterrows(): # index, pandas series (line)
      print Job(s, dep_fn, i == 0)
  else:
    main_help('Need input file (use - for stdin).', main_help=True)

if __name__ == "__main__":
  main()
