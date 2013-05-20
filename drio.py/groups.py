import pandas as pd
from collections import defaultdict

class Group:
  """Captures all the logic related to grouping the different genotypes of a
  snp per sample.
  """
  def __init__(self, fd):
    self.df = pd.read_table(fd)
    self.load_ids_to_groups()

  def load_ids_to_groups(self):
    """Set the number of groups detected in a vcf file. And also a dict that
    tells you to what group a sample belongs to.
    """
    self.groups = set([])
    self.h_group_ids = defaultdict(lambda: set([]))
    self.h_id_to_group = defaultdict(lambda: set([]))
    for i,g in zip(self.df.sample_id, self.df.group_name):
      self.h_group_ids[g].add(i)
      self.h_id_to_group[i] = g
      self.groups.add(g)

  def what_is(self, _id):
    """What's the group for this sample id?
    """
    for g in self.groups:
      if _id in self.h_group_ids[g]:
        return g
    return None

  def id_for_group(self, g):
    return self.h_id_to_group[g]

  def indices_for_grp(self, grp_name):
    return self.h_group_ids[grp_name]

  def num(self):
    return len(self.groups)

  def fresh_hash(self):
    """ return a hash where the keys are the groups and the values set to 0
    """
    _h = defaultdict(lambda: 0)
    very_small = 0.000000000001
    for g in self.groups: _h[g] = { "total": very_small, "var_all": 0 }
    return _h

