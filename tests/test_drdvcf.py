import sys, unittest
from drdvcf import VcfSnp

class TestDrdVcf:
  def setUp(self):
    self.fixtures = {
      "non_syn"       : VcfSnp("Chr1    627540  .       A       G       11.30   .       AC1=1;AC=1;AF1=0.5;AN=2;DP4=3,2,2,0;DP=14;EFF=NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|Tcc/Ccc|S8P|668|RGS12|protein_coding|CODING|ENSMMUT00000009994|exon_1_625681_627664);FQ=14.2;MQ=53;PV4=1,0.31,0.26,1;SF=5;VDB=0.0279 "),
      "syn_coding"    : VcfSnp("Chr1    24428   .       T       G       7.59    .       AC1=2;AC=2;AF1=1;AN=2;DP4=0,0,2,0;DP=2;EFF=SYNONYMOUS_CODING(LOW|SILENT|ggA/ggC|G145|368|HMX1|protein_coding|CODING|ENSMMUT00000019076|exon_1_24407_24477);FQ=-33;MQ=22;SF=4;VDB=0.0133 GT:GQ:SP:PL     .       ."),
      "intron"        : VcfSnp("Chr1    26208   .       N       T       11.10   .       AC1=2;AC=2;AF1=1;AN=2;DP4=0,0,0,2;DP=2;EFF=INTRON(MODIFIER||||368|HMX1|protein_coding|CODING|ENSMMUT00000019076|);FQ=-33;MQ=25;SF=14;VDB=0.0099 GT:GQ:SP:PL     .       .       .       .       .       .       .       .       .       .       .       .       .    ."),
      "intergenic"    : VcfSnp("Chr1    2986    .       G       T       9.52    .       AC1=1;AC=1;AF1=0.5;AN=2;DP4=6,0,0,1;DP=10;EFF=INTERGENIC(MODIFIER|||||||||);FQ=12.3;MQ=60;PV4=0.14,1,1,1;SF=5   GT:GQ:SP:PL     .       .       .       .       .       0/1:41:0:39,0,96        .       .       .       .       .       .       .       .       .    ."),
      "not_annotated" : VcfSnp("20     1110696 rs6040355 A      G,T     67   PASS   NS=2;DP=10;AF=0.333,0.667;AA=T;DB GT:GQ:DP:HQ 1|2:21:6:23,27 2|1:2:0:18,2   2/2:35:4"),
      "full_snp"      : VcfSnp("Chr1    7053    .       A       G       203.90  .       AC1=1;AC=12;AF1=0.5;AN=20;DP4=78,96,122,120;DP=442;EFF=INTERGENIC(MODIFIER|||||||||);FQ=225;MQ=45;PV4=0.75,1,1.6e-09,1;SF=0,1,2,7,9,11,14,15,16,19;VDB=0.0399   GT:GQ:SP:PL     0/1:99:1:255,0,255      1/1:99:0:255,111,0      0/1:99:1:255,0,255      .       .       .      .       0/1:99:4:218,0,255      .       1/1:99:0:255,132,0      .       0/1:99:5:255,0,255      .       .       0/1:99:0:234,0,255      0/1:99:0:108,0,126      0/1:99:1:255,0,254      .       .       0/1:99:2:255,0,255"),
      "all_same_gtype": VcfSnp("Chr1    2222    .       A       G       203.90  .       AC1=1;AC=12;AF1=0.5;AN=20;DP4=78,96,122,120;DP=442;EFF=INTERGENIC(MODIFIER|||||||||);FQ=225;MQ=45;PV4=0.75,1,1.6e-09,1;SF=0,1,2,7,9,11,14,15,16,19;VDB=0.0399   GT:GQ:SP:PL     0/1:99:1:255,0,255      0/1:99:0:255,111,0      0/1:99:1:255,0,255      .       .       .      .       0/1:99:4:218,0,255"),
      "second_vars"   : VcfSnp("Chr1    8239092 .       G       T,C     222.00  AC1=2;AC=12,6;AF1=1;AN=18;DP4=0,0,186,177;DP=385;EFF=INTERGENIC(MODIFIER|||||||||);FQ=-126;MQ=40;SF=3,4,7,8,9,10,12,14,15;VDB=0.0384;RDP=37,55,44,34,43,39,48,49,45,56,40,61,43,45,43,28,51,52,45,46   GT:GQ:SP:PL     .     .     .       2/2:99:0:255,.,.,99,.,0    1/1:99:0:255,123,0,.,.,.        .       .       2/2:99:0:255,.,.,138,.,0        1/1:99:0:255,132,0,.,., .       1/1:99:0:255,166,0,.,.,.        1/1:99:0:255,111,0,.,.,.        .       2/2:99:0:255,.,.,126,.,0        .       1/1:99:0:255,123,0,.,.,.       1/1:99:0:255,66,0,.,.,. .       .       .       .")
    }

  def test_all_samples_same_genotype(self):
    assert self.fixtures['all_same_gtype'].species_snp()
    assert not self.fixtures['full_snp'].species_snp()
    assert not self.fixtures['second_vars'].species_snp()
    pass

  def test_if_we_correctly_detect_annotated_snps(self):
    for ft in ['non_syn', 'syn_coding', 'intron']:
      assert self.fixtures[ft].annotated
    assert self.fixtures['not_annotated'].annotated == False

  def test_we_properly_extract_the_annotated_gene(self):
    assert self.fixtures['non_syn'].gene    == 'RGS12'
    assert self.fixtures['syn_coding'].gene == 'HMX1'
    assert self.fixtures['intergenic'].gene == ''

  def test_checking_the_functional_consequence_of_the_snp(self):
    assert self.fixtures['non_syn'].func_cons == "NON_SYNONYMOUS_CODING"
    assert self.fixtures['not_annotated'].func_cons == ""

  def test_checking_the_impact_of_the_snp(self):
    assert self.fixtures['non_syn'].impact == "MODERATE"
    assert self.fixtures['not_annotated'].impact == ""

  def test_genotypes(self):
    gtypes = self.fixtures['full_snp'].gtypes()
    assert len(gtypes) == 10
    # Each genotype
    assert gtypes[0][0] == 0; assert gtypes[0][1] == 1
    assert gtypes[1][0] == 1; assert gtypes[1][1] == 1
    assert gtypes[2][0] == 0; assert gtypes[2][1] == 1
    assert gtypes[7][0] == 0; assert gtypes[7][1] == 1
    assert gtypes[9][0] == 1; assert gtypes[9][1] == 1
    assert gtypes[11][0] == 0; assert gtypes[11][1] == 1
    assert gtypes[14][0] == 0; assert gtypes[14][1] == 1
    assert gtypes[15][0] == 0; assert gtypes[15][1] == 1
    assert gtypes[16][0] == 0; assert gtypes[16][1] == 1
    assert gtypes[19][0] == 0; assert gtypes[19][1] == 1
    for i in [3, 4, 5, 6, 8, 10, 12, 13, 17, 18]:
      assert not gtypes.has_key(i)

if __name__ == '__main__':
  import nose
  nose.run(defaultTest=__name__)
