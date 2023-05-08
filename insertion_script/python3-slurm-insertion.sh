#! /bin/bash

DATADIR='/gpfs/commons/groups/gursoy_lab/ubaymuradov/main_node'
TESTDATA='/gpfs/commons/groups/gursoy_lab/ubaymuradov/PrecisionChain/data'
BENCHMARKS='/gpfs/commons/groups/gursoy_lab/ubaymuradov/benchmarks'
python3 buildChain.py -cn=combchain -dr=$DATADIR

python3 createStream-OMOP-Domain.py -cn=combchain --datadir=$DATADIR -hp=$TESTDATA/mapping_vocab/CONCEPT_ANCESTOR.csv -dp=$TESTDATA/clinical/

python3 insertData-OMOP-Domain.py -cn=combchain --datadir=$DATADIR -hp=$TESTDATA/mapping_vocab/CONCEPT_ANCESTOR.csv -dp=$TESTDATA/clinical/

python3 createStream-OMOP-Person.py -cn=combchain --datadir=$DATADIR

python3 insertData-OMOP-Person.py -cn=combchain --datadir=$DATADIR -dp=$TESTDATA/clinical/ --personPath=$TESTDATA/clinical/person.csv

python3 createStream-variants.py -cn=combchain --datadir=$DATADIR

# run this last
python3 insertData-variant.py -cn=combchain --datadir=$DATADIR -mf=$TESTDATA/mapping_vocab/person_sample_mapping.txt -dp=$TESTDATA/vcf

python3 insertData-variantPerson.py -cn=combchain --datadir=$DATADIR -mf=$TESTDATA/mapping_vocab/person_sample_mapping.txt -dp=$TESTDATA/vcf

python3 createStream-gtf.py -cn=combchain -dr=$DATADIR

python3 insertData-gtf.py -cn=combchain -dr=$DATADIR -gp=$TESTDATA/gtf/ -vp=$TESTDATA/vcf/

du -h $DATADIR/combchain >$BENCHMARKS/single-node-space-usage-stats.txt
multichain-cli combchain -datadir=$DATADIR liststreams | jq -c '.[] | {name: .name, items: .items, keys: .keys, publishers: .publishers}' >$BENCHMARKS/single-node-stats.txt
