#!/usr/bin/env python
#
# (rm -rf plots ; /stornext/snfs6/rogers/drio_scratch/py.analysis/downstream/cnv/parallel_log_plots.py | xargs -t -I {} -n 1 -P 32 sh -c "{}" &> /dev/null ) &
#
from collections import defaultdict
import sys
import os
import fnmatch
import re
from drdcommon import mkdir_p

RE_PATH = re.compile(r"/(\d+)_vs_(\d+)_Chr(\w+).txt.gz")
def extract_info(filename):
  # ./log2ratios/5000/30158/35087_vs_30158_Chr15.txt.gz
  match = RE_PATH.search(filename)
  if match:
    id_control, sample_id, chrm = match.group(1), match.group(2), match.group(3)
    return id_control, sample_id, chrm
  else:
    raise Exception('Unexpected fname provided: %s' % filename)

def main():
  current_dir = os.path.dirname(os.path.abspath(__file__))
  single_plot_py = current_dir + "/plot_log2ratios.py"
  pngs_by_chrm = defaultdict(list)
  for root, dirnames, filenames in os.walk('.'):
    for filename in fnmatch.filter(filenames, '*.txt.gz'):
      full_path      = os.path.join(root, filename)
      idc, ids, chrm = extract_info(full_path)
      odir           = "plots/%s" % chrm
      png            = "%s_%s_%s.png" % (idc, ids, chrm)
      pngs_by_chrm[chrm].append(png)
      mkdir_p(odir)
      cmd  =  "gzip -cd %s | %s '%s_vs_%s' %s" %\
        (full_path, single_plot_py, idc, ids, odir + "/" + png)
      #print "echo '%s' | submit -s%s -c1 -m1" % (cmd, png)
      print cmd

  # The index files
  for chrm, list_files in pngs_by_chrm.items():
    fd = open("plots/%s/index.html" % chrm, 'w')
    for png in list_files:
      fd.write("<img src='%s'><br>\n" % png)
    fd.close

if __name__ == "__main__":
  main()
