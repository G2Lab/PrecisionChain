#! /bin/bash

DATADIR='/gpfs/commons/groups/gursoy_lab/ubaymuradov/main_node'
TESTDATA='/gpfs/commons/groups/gursoy_lab/ubaymuradov/PrecisionChain/data'
BENCHMARKS='/gpfs/commons/groups/gursoy_lab/ubaymuradov/benchmarks'
CODE='/gpfs/commons/groups/gursoy_lab/ubaymuradov/PrecisionChain/code'
export NTASKS=4

python3 $CODE/insertData-variant.py -cn=combchain --datadir=$DATADIR -mf=$TESTDATA/mapping_vocab/person_sample_mapping.txt -dp=$TESTDATA/vcf
sleep 300 # give nodes time to transfer data
du -h $DATADIR/combchain >$BENCHMARKS/main-node-space-usage-stats-$NTASKS-tasks.txt
multichain-cli combchain -datadir=$DATADIR liststreams | jq -c '.[] | {name: .name, items: .items, keys: .keys, publishers: .publishers}' >$BENCHMARKS/main-node-stats-$NTASKS-tasks.txt
echo "finished" >/gpfs/commons/groups/gursoy_lab/ubaymuradov/main_node/status.txt
