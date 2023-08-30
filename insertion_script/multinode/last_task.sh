#!/bin/bash

CN=4_node
DR='/gpfs/commons/groups/gursoy_lab/aelhussein/blockchain/multichain/multinode/main_node'
PT=/gpfs/commons/groups/gursoy_lab/aelhussein/blockchain/public
AT=/gpfs/commons/groups/gursoy_lab/aelhussein/blockchain/data/preprocessing
CODE='/gpfs/commons/groups/gursoy_lab/aelhussein/blockchain/public/code'
BENCHMARKS='/gpfs/commons/groups/gursoy_lab/aelhussein/blockchain/multichain/multinode/benchmarks'
export NTASKS=4

python3 $CODE/insertData-variant.py -cn=$CN --datadir=$DR -mf=--datadir=$DR -dp=$PT/data/vcf/wgs -mf=$PT/data/samples/metadata.csv -np=$PPL_GS -sq=WGS
sleep 300 # give nodes time to transfer data
du -h $DR/main_node/$CN >$BENCHMARKS/main-node-space-usage-stats-$NTASKS-tasks.txt
multichain-cli $CN -datadir=$DR liststreams | jq -c '.[] | {name: .name, items: .items, keys: .keys, publishers: .publishers}' >$BENCHMARKS/main-node-stats-$NTASKS-tasks.txt
echo "finished" >$DR/main_node/status.txt
