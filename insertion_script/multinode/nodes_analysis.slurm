#!/bin/bash
#SBATCH --job-name=multichain-nodes
#SBATCH --partition=pe2
#SBATCH --mail-user=aelhussein@nygenome.org
#SBATCH --ntasks-per-node=1
#SBATCH --array=1-3
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=4:00:00
#SBATCH --output=stdout_%j_%a.log

module load multichain
module load bcftools
source ~/anaconda3/bin/activate ~/anaconda3/envs/blockchain_env/
conda activate blockchain_env

export NTASKS=4

CN=4_node
#SAMPLE NUMBERS
PPL=100
PPL_GS=100
PPL_ES=0
PPL_AF=0
PPL_BC=0

MULTINODE='/gpfs/commons/groups/gursoy_lab/aelhussein/blockchain/multichain/multinode'
NODEADDRESS=$(cat $MULTINODE/main_node/nodeaddress.txt)
# Seed node port
PORT=$(echo $NODEADDRESS | cut -d':' -f2)
# Node port
NODEPORT=$(($PORT + $SLURM_ARRAY_TASK_ID * 2))
RPCPORT=$(($PORT + $SLURM_ARRAY_TASK_ID * 2 - 1))

CODE='/gpfs/commons/groups/gursoy_lab/aelhussein/blockchain/public/code'
DR=$MULTINODE/node$SLURM_ARRAY_TASK_ID
MAINDATADIR=$MULTINODE/main_node
PT=/gpfs/commons/groups/gursoy_lab/aelhussein/blockchain/public
AT=/gpfs/commons/groups/gursoy_lab/aelhussein/blockchain/data/preprocessing
BENCHMARKS='/gpfs/commons/groups/gursoy_lab/aelhussein/blockchain/multichain/multinode/benchmarks'

# Re-create the data directories
rm -rf $MULTINODE/node$SLURM_ARRAY_TASK_ID
mkdir $MULTINODE/node$SLURM_ARRAY_TASK_ID

# Add node ports to multichain.conf
echo "port=$NODEPORT" >>$MULTINODE/node$SLURM_ARRAY_TASK_ID/multichain.conf
echo "rpcport=$RPCPORT" >>$MULTINODE/node$SLURM_ARRAY_TASK_ID/multichain.conf
echo "about to start"
# capture the output of multichaind command to a file
multichaind $NODEADDRESS autosubscribe=streams -port=$NODEPORT -rpcport=$RPCPORT -datadir=$DR -daemon &>start$SLURM_ARRAY_TASK_ID.out
cat start$SLURM_ARRAY_TASK_ID.out | grep -m 1 receive | sed "s/$CN/$CN -datadir=${MAINDATADIR//\//\\/}/" >>$MULTINODE/main_node/permissions.txt
# stop multichaind
sleep 20
echo "finished sleeping, about to start with permissions"
multichaind $NODEADDRESS -autosubscribe=streams -port=$NODEPORT -rpcport=$RPCPORT -datadir=$DR -daemon
# Subscribe to all streams
echo "python scripts starting"
python3 $CODE/insertData-OMOP-Person.py -cn=$CN --datadir=$DR -dp=$PT/data/clinical/ --personPath=$PT/data/clinical/person.csv -np=$PPL -mf=$PT/data/samples/metadata.csv
python3 $CODE/insertData-OMOP-Domain.py -cn=$CN --datadir=$DR -hp=$PT/data/mapping_vocab/CONCEPT_ANCESTOR.csv -dp=$PT/data/clinical/ -np=$PPL -mf=$PT/data/samples/metadata.csv
# WGS
python3 $CODE/insertData-variantPerson.py -cn=$CN --datadir=$DR -dp=$PT/data/vcf/wgs -mf=$PT/data/samples/metadata.csv -np=$PPL_GS -sq=WGS
# Analysis
python3 $CODE/insertData-analysis.py -cn=$CN --datadir=$DR -dp=$PT/data/vcf/wgs -mf=$PT/data/samples/metadata.csv -pc=$AT/ -rf=$AT/relatedness -np=$PPL_GS
echo "python scripts ran"

awk '{print $1 - 1}' $MAINDATADIR/ntasks.txt >$MAINDATADIR/ntasks_copy.txt
mv $MAINDATADIR/ntasks_copy.txt $MAINDATADIR/ntasks.txt

# read file, if it contains "finished", then run the command else sleep
while true; do
    if grep -q "finished" $MAINDATADIR/status.txt; then
        du -h $DR/$CN >$BENCHMARKS/node$SLURM_ARRAY_TASK_ID-space-usage-stats-$NTASKS-tasks.txt
        multichain-cli $CN -datadir=$DR liststreams | jq -c '.[] | {name: .name, items: .items, keys: .keys, publishers: .publishers}' >$BENCHMARKS/node$SLURM_ARRAY_TASK_ID-stats-$NTASKS-tasks.txt
        break
    else
        sleep 10
    fi
done

sleep 3600
