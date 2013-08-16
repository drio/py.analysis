import sys, gzip, bz2
import os, errno
import glob
from time import gmtime, strftime
import re

_human_genome_sizes = {
  '1': 249250621, '2': 243199373, '3': 198022430, '4': 191154276, '5': 180915260,
  '6': 171115067, '7': 159138663, '8': 146364022, '9': 141213431, '10': 135534747,
  '11': 135006516, '12': 133851895, '13': 115169878, '14': 107349540, '15': 102531392,
  '16': 90354753, '17': 81195210, '18': 78077248, '19': 59128983, '20': 63025520,
  '21': 48129895, '22': 51304566, 'X': 155270560, 'Y': 59373566, 'M': 16571,
}

def ds_chrm_coor(cell_val=0, log=True):
  '''
    d[chrm] : list(each cell is a coordinate)
    Cells values will be initialized to cell_val
  '''
  d_to_list = {}
  for chrm, size in _human_genome_sizes.iteritems():
    if log:
      print >> sys.stderr, "ds_chrm_coor():" , chrm , size , "mem:" , memory_usage()
    d_to_list[chrm] = [cell_val] * size

  return d_to_list

def setup_logging(logging):
  logging.basicConfig(format='%(asctime)s >> %(message)s ', level=logging.DEBUG)

def log(msg, output_to=sys.stderr):
  now = strftime("%Y-%m-%d %H:%M:%S", gmtime())
  output_to.write(now + ">> " + msg + "\n")

def files_in_dir(directory, pattern='*'):
  """ return a list of all files in directory that match pattern """
  #files = [directory + "/" + f for f in os.listdir(directory) if os.path.isfile(directory + "/" + f)]
  return glob.glob(directory + "/" + pattern)
  #return fnmatch.filter(files, pattern)

def rec_find_files(directory, pattern):
  """ Find recursively all files that match pattern in directory """
  for root, dirs, files in os.walk(directory):
    for basename in files:
      if fnmatch.fnmatch(basename, pattern):
        filename = os.path.join(root, basename)
        yield filename

def error(msg, usage=None):
  sys.stderr.write("Ups!: " + msg + "\n")
  if usage:
    sys.stderr.write(usage)
  sys.exit(1)

# Are we running this under ipython?
def in_ipython():
  try:
    __IPYTHON__
  except NameError:
    return False
  else:
    return True

SUBMIT = 'submit'
def cmd(root_dir, o_dir, seed, b):
  d   = root_dir + "/" + o_dir
  cmd = seed % (d, b.path, d + "/" + b.id)
  print "%s -s %s '%s'" % (SUBMIT, o_dir + "_" + b.id, cmd)

# Do we have data in stdin?
def data_in_stdin():
  return not sys.stdin.isatty()

def xopen(fn):
  if fn == '-':
    return sys.stdin
  else:
    ext = fn[-3:]
    try:
      if ext == '.gz':
        return gzip.open(fn, 'rb')
      elif ext == '.bz2':
        return bz2.open(fn, 'rb')
      else:
        return open(fn, 'r')
    except:
      error("Problem opening input file: %s\n" % (fn,))

def memory_usage():
  """Memory usage of the current process in Megabytes."""
  status = None
  result = {'peak': 0, 'rss': 0}
  try:
    # This will only work on systems with a /proc file system
    # (like Linux).
    status = open('/proc/self/status')
    for line in status:
      parts = line.split()
      key = parts[0][2:-1].lower()
      if key in result:
        result[key] = int(parts[1])
  except:
    pass
  finally:
    if status is not None:
      status.close()
  return float(result['peak']) / float(1000)

def mkdir_p(path):
  try:
    os.makedirs(path)
  except OSError as exc: # Python >2.5
    if exc.errno == errno.EEXIST and os.path.isdir(path):
      pass
    else:
      raise

def canonic_chrm(s):
  return re.sub(r'(^[cC]hrm?)', '', s)
