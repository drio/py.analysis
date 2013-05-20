
from drdvcf import *

v = VcfSnp("Chr1    232980  .       T       A,G     62.75   .       AC1=2;AC=1,2;AF1=1;AN=4;DP4=1,1,6,2;DP=10;EFF=INTERGENIC(MODIFIER|||||||||);FQ=-39;MQ=60;PV4=1,1,1,0.42;SF=19,36;VDB=0.0216       GT:GQ:PL        .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       1/2:21:133,30,18,110,0,107        .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       .       0/2:50:73,.,.,0,.,47      .       .       .       .       .       .       .       .       .       .       .       .       .")
v.add_info('OT')
print str(v)
