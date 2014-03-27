▾ CNVPostProcess/
      copyNumberDistribution.R
      mergeCoord.pl
      methodSD1.sh
      methodSD2.sh
      twoPartOvp_mgsrt.pl
      wssd_picker.pl


Attached in the present mail you are going to find a set of scripts we are
using once mrCaNaVar step is complete. Main scripts are
copyNumberDistribution.R, methodSD1.sh and methodSD2.

copyNumberDistribution.R Is an R script which calculates copy number
distribution genome-wide, control regions, non control regions and non-control
region in autosomes. Pay attention, you must edit this script to adapt to your
project environment. As you will see this script is though to be run in 6
samples, so you must redefine R vectors.

methodSD1.sh This script performs Method 1. As in the first case, this script
must be adapted to your project environment. Please do not forget to change
binary paths for bedtools according to your computer cluster.

methodSD2.sh And the final script for Method 2. This script calls
mergeCoord.pl, toPartOvp_mgsrt.pl and wssd_picker.pl perl scripts. Please
review all of them in order to adapt it to your environment.
