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

source /gpfs/commons/home/ubaymuradov/.bashrc

export NTASKS=4

USERHOME=/gpfs/commons/groups/gursoy_lab/ubaymuradov
NODEADDRESS=$(cat $USERHOME/main_node/nodeaddress.txt)
# Seed node port
PORT=$(echo $NODEADDRESS | cut -d':' -f2)
# Node port
NODEPORT=$(($PORT + $SLURM_ARRAY_TASK_ID * 2))
RPCPORT=$(($PORT + $SLURM_ARRAY_TASK_ID * 2 - 1))

CODE='/gpfs/commons/groups/gursoy_lab/ubaymuradov/PrecisionChain/code'
DATADIR=$USERHOME/node$SLURM_ARRAY_TASK_ID
MAINDATADIR=$USERHOME/main_node
PT=/gpfs/commons/groups/gursoy_lab/aelhussein/blockchain/public
AT=/gpfs/commons/groups/gursoy_lab/aelhussein/blockchain/data/preprocessing
BENCHMARKS='/gpfs/commons/groups/gursoy_lab/ubaymuradov/benchmarks'

# Re-create the data directories
rm -rf $USERHOME/node$SLURM_ARRAY_TASK_ID
mkdir $USERHOME/node$SLURM_ARRAY_TASK_ID

#SAMPLE NUMBERS
PPL=4000
PPL_GS=2000
PPL_ES=1200
PPL_AF=400
PPL_BC=400

# Add node ports to multichain.conf
echo "port=$NODEPORT" >>$USERHOME/node$SLURM_ARRAY_TASK_ID/multichain.conf
echo "rpcport=$RPCPORT" >>$USERHOME/node$SLURM_ARRAY_TASK_ID/multichain.conf
echo "about to start"
# capture the output of multichaind command to a file
multichaind $NODEADDRESS autosubscribe=streams -port=$NODEPORT -rpcport=$RPCPORT -datadir=$DATADIR -daemon &>start$SLURM_ARRAY_TASK_ID.out
cat start$SLURM_ARRAY_TASK_ID.out | grep -m 1 receive | sed "s/combchain/combchain -datadir=${MAINDATADIR//\//\\/}/" >>$USERHOME/main_node/permissions.txt
# stop multichaind
sleep 20
echo "finished sleeping, about to start with permissions"
multichaind $NODEADDRESS -autosubscribe=streams -port=$NODEPORT -rpcport=$RPCPORT -datadir=$DATADIR -daemon
# Subscribe to all streams

echo "python scripts starting"
python3 $CODE/insertData-OMOP-Person.py -cn=combchain --datadir=$DATADIR -dp=$PT/clinical/ --personPath=$PT/clinical/person.csv
python3 $CODE/insertData-OMOP-Domain.py -cn=combchain --datadir=$DATADIR -hp=$PT/mapping_vocab/CONCEPT_ANCESTOR.csv -dp=$PT/clinical/
python3 $CODE/insertData-variantPerson.py -cn=combchain --datadir=$DATADIR -mf=$PT/mapping_vocab/person_sample_mapping.txt -dp=$PT/vcf
echo "python scripts ran"

awk '{print $1 - 1}' $MAINDATADIR/ntasks.txt >$MAINDATADIR/ntasks_copy.txt
mv $MAINDATADIR/ntasks_copy.txt $MAINDATADIR/ntasks.txt

# read file, if it contains "finished", then run the command else sleep
while true; do
    if grep -q "finished" $MAINDATADIR/status.txt; then
        du -h $DATADIR/combchain >$BENCHMARKS/node$SLURM_ARRAY_TASK_ID-space-usage-stats-$NTASKS-tasks.txt
        multichain-cli combchain -datadir=$DATADIR liststreams | jq -c '.[] | {name: .name, items: .items, keys: .keys, publishers: .publishers}' >$BENCHMARKS/node$SLURM_ARRAY_TASK_ID-stats-$NTASKS-tasks.txt
        break
    else
        sleep 10
    fi
done

sleep 3600
