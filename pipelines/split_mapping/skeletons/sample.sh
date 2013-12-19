#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo -ne "Path to your bam (make sure we have rw/permissions unless index exists): "
read bam
echo -ne "Path to fasta file: "
read fasta
echo -ne "Sample id: "
read id
echo -ne "# of max cores: "
read cores
echo -ne "# of reads per split: "
read nrps
echo -ne "path to local tmp scratch: "
read tmp
echo -ne "what's the rest server? (hit enter to skip this)"
read host
if [ ".$rest_host" != "." ];then
    echo -ne "what's the rest port? "
    read port
    echo -ne "what's the http auth user? "
    read user
    echo -ne "what's the http auth password? "
    read pass
fi


rm -f ./`basename $bam`
ln -s $bam

cat $DIR/config.sh | \
    sed -e "s!_BAM_!$bam!; s!_FA_!$fasta!; \
            s!_ID_!$id!; s!_CR_!$cores!; \
            s!_NR_!$nrps!; s!_TMP_!$tmp!; \
            s!_USER_!$user!; s!_PWD_!$pass!; \
            s!_URL_!$URL!" > config.sh

cp $DIR/run.sh .

cat <<EOF

Example of executing the pipeline in HPC a cluster. submit knows
the details of how to send the jobs to the cluster:

export SAPI_USER=$user
export SAPI_PWD=$pass
$ (submit -s _one 'source config.sh; sapi.py -u \$url -i \$id -b \$bam -n \$nrps init -x' | bash) | submit -s _two -f - "./run.sh | bash" | bash

EOF
