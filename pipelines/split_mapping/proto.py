#!/usr/bin/env python

import sys, os, tempfile
from random import randint
__lib_dir = os.path.dirname( os.path.realpath(__file__) ) + "/../drio.py"
sys.path.append(__lib_dir)
import drdcommon as common
import ngsproject as np

TOOL_NAME = "split_mapping_pipe.py"
MAIN_HELP = """Usage: %s <ref.fa> <input.bam> <sample_name> <#_of_reads_per_split>\n""" % (TOOL_NAME)

class Pipeline(object):
  def __init__(self):
    pass

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

class Cmd(object):
  def __init__(self, name):
    self.name = name
    self.cd = "."

class Picard(Cmd):
  def __init__(self):
    self.common   = "cd %s; java -jar $PICARD/%s TMP_DIR=/tmp VALIDATION_STRINGENCY=LENIENT "

class SamToFq(Picard):
  def __init__(self, input_fn, oseed):
    Picard.__init__(self)
    Cmd.__init__(self, "Picard-Sam2Fq")
    self.jar       = "SamToFastq.jar"
    self.input_fn  = input_fn
    self.output_e1 = "%s.e1.fq" % oseed
    self.output_e2 = "%s.e2.fq" % oseed

  def __str__(self):
    return (self.common + "I=%s F=%s F2=%s") % \
      (self.cd, self.jar, self.input_fn, self.output_e1, self.output_e2)

class Unix(Cmd):
  def __init__(self, cmd=""):
    self.cmd = cmd

  def __str__(self):
    return self.cmd

class Args(object):
  def __init__(self):
    self.process()

  def process(self):
    if len(sys.argv) == 5:
      # TODO: common, check file exists....
      self.ref_fn, self.inputb_fn, self.sample, self.rp_split = sys.argv[1:]
      self.rp_split = int(self.rp_split)
    else:
      self.main_help('Incorrect num of params.', main_help=True)
      sys.exit(1)

  def main_help(self, msg=None, main_help=False):
    if msg: sys.stderr.write('ERROR: ' + msg + "\n")
    if main_help: sys.stderr.write(MAIN_HELP)

def do_work():
  lc   = []
  a    = Args()
  lc.append(Unix("mkdir " + a.sample))
  lc.append(SamToFq(a.inputb_fn, a.sample))

  for c in lc:
    print c

if __name__ == "__main__":
  do_work()
