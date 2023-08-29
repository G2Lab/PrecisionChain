#!/bin/bash
#SBATCH --job-name=multichain-nodes
#SBATCH --partition=pe2
#SBATCH --mail-user=ubaymuradov@nygenome.org
#SBATCH --ntasks-per-node=1
#SBATCH --array=1-3
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=stdout_%j_%a.log
set -e

source /gpfs/commons/home/ubaymuradov/.bashrc

CN=query_4k
DR='/gpfs/commons/groups/gursoy_lab/ubaymuradov/main_node'
PT=/gpfs/commons/groups/gursoy_lab/aelhussein/blockchain/public
AT=/gpfs/commons/groups/gursoy_lab/aelhussein/blockchain/data/preprocessing
CODE='/gpfs/commons/groups/gursoy_lab/ubaymuradov/PrecisionChain/code'

rm -rf $DR/$CN

#SAMPLE NUMBERS
PPL=100
PPL_GS=100
PPL_ES=0
PPL_AF=0
PPL_BC=0

export NTASKS=4

# create permissions file
touch $DR/permissions.txt
chmod 666 $DR/permissions.txt

echo $NTASKS >$DR/ntasks.txt

# provide file for nohup output
nohup bash /gpfs/commons/groups/gursoy_lab/ubaymuradov/PrecisionChain/insertion_script/slurm/permissions_watch.sh >permissions_watch.out 2>&1 &
# keep track of status
echo "started" >$DR/status.txt
python3 $CODE/buildChain.py -cn=$CN -dr=$DR
# write stream names to a file
python3 $CODE/createStream-OMOP-Domain.py -cn=$CN --datadir=$DR -hp=$PT/data/mapping_vocab/CONCEPT_ANCESTOR.csv -dp=$PT/data/clinical/
python3 $CODE/createStream-OMOP-Person.py -cn=$CN --datadir=$DR
python3 $CODE/createStream-variants.py -cn=$CN --datadir=$DR
python3 $CODE/createStream-gtf.py -cn=$CN -dr=$DR

sbatch $CODE/../insertion_script/slurm/nodes_analysis.slurm

python3 $CODE/insertData-OMOP-Person.py -cn=$CN --datadir=$DR -dp=$PT/data/clinical/ --personPath=$PT/data/clinical/person.csv -np=$PPL -mf=$PT/data/samples/metadata.csv
python3 $CODE/insertData-OMOP-Domain.py -cn=$CN --datadir=$DR -hp=$PT/data/mapping_vocab/CONCEPT_ANCESTOR.csv -dp=$PT/data/clinical/ -np=$PPL -mf=$PT/data/samples/metadata.csv
# WGS
python3 $CODE/insertData-variantPerson.py -cn=$CN --datadir=$DR -dp=$PT/data/vcf/wgs -mf=$PT/data/samples/metadata.csv -np=$PPL_GS -sq=WGS
# WES
python3 $CODE/insertData-variantPerson.py -cn=$CN --datadir=$DR -dp=$PT/data/vcf/wes -mf=$PT/data/samples/metadata.csv -np=$PPL_ES -sq=WES
# AFFYMETRIX
python3 $CODE/insertData-variantPerson.py -cn=$CN --datadir=$DR -dp=$PT/data/vcf/affymetrix -mf=$PT/data/samples/metadata.csv -np=$PPL_AF -sq=Affymetrix
# BEADCHIP
python3 $CODE/insertData-variantPerson.py -cn=$CN --datadir=$DR -dp=$PT/data/vcf/beadchip -mf=$PT/data/samples/metadata.csv -np=$PPL_BC -sq=BeadChip
# Analysis
python3 $CODE/insertData-analysis.py -cn=$CN --datadir=$DR -dp=$PT/data/vcf/wgs -mf=$PT/data/samples/metadata.csv -pc=$AT/ -rf=$AT/relatedness -np=$PPL_GS

#python3 $CODE/insertData-gtf.py -cn=$CN -dr=$DR -gp=$PT/gtf/ -vp=$PT/vcf/

while true; do
    if [ $(cat $DR/ntasks.txt) -eq 1 ]; then
        echo "about to run last task"
        nohup bash /gpfs/commons/groups/gursoy_lab/ubaymuradov/PrecisionChain/insertion_script/slurm/last_task.sh >last_task.out 2>&1 &
        break
    else
        sleep 10
    fi
done
