import sys, unittest
import binner
import nose
from nose.tools import assert_raises

class TestDrdVcf:
  def setUp(self):
    pass

  def test_basic(self):
    n_max_reads = 2
    chrms = ['chr1', 'chr2']
    lengths = {'chr1':100, 'chr2':200}
    b = binner.Binner(n_max_reads, chrms, lengths)

    assert not b.add('chr1', 10)
    win = b.add('chr1', 11)
    assert win == 'chr1 1 11'

    assert not b.add('chr1', 50)
    win = b.add('chr1', 60)
    assert win == 'chr1 11 60'

    assert not b.add('chr1', 90)
    win = b.add('chr2', 10)
    print win
    assert win == 'chr1 60 100'

    assert not b.add('chr2', 70)
    win = b.add('chr2', 80)
    assert win == 'chr2 1 80'

    assert not b.add('chr2', 120)
    win = b.add('chr2', 130)
    assert win == 'chr2 80 130'

if __name__ == '__main__':
  nose.run(defaultTest=__name__)
