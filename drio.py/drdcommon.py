import sys, gzip, bz2
import os, errno
import fnmatch
from time import gmtime, strftime

def setup_logging(logging):
  logging.basicConfig(format='%(asctime)s >> %(message)s ', level=logging.DEBUG)

def log(msg, output_to=sys.stderr):
  now = strftime("%Y-%m-%d %H:%M:%S", gmtime())
  output_to.write(now + ">> " + msg + "\n")

def files_in_dir(directory, pattern='*'):
  """ return a list of all files in directory that match pattern """
  files = [f for f in os.listdir(directory) if os.path.isfile(f)]
  return fnmatch.filter(files, pattern)

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
