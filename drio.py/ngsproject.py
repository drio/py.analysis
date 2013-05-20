import pandas as pd
import drdcommon as common

class SeqEvent():
  def __init__(self, h):
    self.set_attr(h)
    if self.id == '':
      common.error('I cannot create a bam without id ')
    if not self.valid_path():
      common.error('Invalid path [%s] for bam with id = [%s]' % (self.path, self.id))

  def set_attr(self, h):
    self.library_name = h['library_name']
    self.id           = str(h['DNAid']).rstrip()
    self.merge        = h['merge']
    self.path         = h['bam_path'].rstrip()
    self.project      = h['project']

  def valid_path(self):
    try:
      with open(self.path) as f:
        return True
    except IOError as e:
      return False

class NgsPrj:
  MANDATORY_COLS = set(['library_name', 'DNAid', 'merge', 'bam_path', 'project'])

  def __init__(self, fn):
    self.df = pd.read_table(fn)
    for c in self.MANDATORY_COLS:
      if c not in self.df.columns:
        common.error("I couldn't find column [%s] in project tsv" % c)

  def seq_events(self):
    events = []
    for e_line in self.df.values: # iterate over the lines
      h = {}
      for k, v in zip(self.df.columns, e_line): h[k] = v
      events.append(SeqEvent(h))
    return events
