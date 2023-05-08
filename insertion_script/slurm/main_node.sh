#! /bin/bash

DATADIR='/gpfs/commons/groups/gursoy_lab/ubaymuradov/main_node'
TESTDATA='/gpfs/commons/groups/gursoy_lab/ubaymuradov/PrecisionChain/data'
CODE='/gpfs/commons/groups/gursoy_lab/ubaymuradov/PrecisionChain/code'

export NTASKS=4

# clean up $DATADIR
rm -rf $DATADIR
mkdir $DATADIR
# create permissions file
touch $DATADIR/permissions.txt
chmod 666 $DATADIR/permissions.txt

echo $NTASKS >$DATADIR/ntasks.txt

# provide file for nohup output
nohup bash /gpfs/commons/groups/gursoy_lab/ubaymuradov/PrecisionChain/insertion_script/slurm/permissions_watch.sh >permissions_watch.out 2>&1 &
# keep track of status
echo "started" >$DATADIR/status.txt
python3 $CODE/buildChain.py -cn=combchain -dr=$DATADIR
# write stream names to a file
python3 $CODE/createStream-OMOP-Domain.py -cn=combchain --datadir=$DATADIR -hp=$TESTDATA/mapping_vocab/CONCEPT_ANCESTOR.csv -dp=$TESTDATA/clinical/
python3 $CODE/createStream-OMOP-Person.py -cn=combchain --datadir=$DATADIR
python3 $CODE/createStream-variants.py -cn=combchain --datadir=$DATADIR
python3 $CODE/createStream-gtf.py -cn=combchain -dr=$DATADIR

sbatch $CODE/../insertion_script/slurm/nodes.slurm

python3 $CODE/insertData-OMOP-Person.py -cn=combchain --datadir=$DATADIR -dp=$TESTDATA/clinical/ --personPath=$TESTDATA/clinical/person.csv
python3 $CODE/insertData-OMOP-Domain.py -cn=combchain --datadir=$DATADIR -hp=$TESTDATA/mapping_vocab/CONCEPT_ANCESTOR.csv -dp=$TESTDATA/clinical/
python3 $CODE/insertData-variantPerson.py -cn=combchain --datadir=$DATADIR -mf=$TESTDATA/mapping_vocab/person_sample_mapping.txt -dp=$TESTDATA/vcf
python3 $CODE/insertData-gtf.py -cn=combchain -dr=$DATADIR -gp=$TESTDATA/gtf/ -vp=$TESTDATA/vcf/

while true; do
    if [ $(cat $DATADIR/ntasks.txt) -eq 1 ]; then
        echo "about to run last task"
        nohup bash /gpfs/commons/groups/gursoy_lab/ubaymuradov/PrecisionChain/insertion_script/slurm/last_task.sh >last_task.out 2>&1 &
        break
    else
        sleep 10
    fi
done
