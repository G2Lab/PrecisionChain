#! /bin/bash
set -e
CN=query_4k
DR='/gpfs/commons/groups/gursoy_lab/ubaymuradov/main_node'
PT=/gpfs/commons/groups/gursoy_lab/aelhussein/blockchain/public
AT=/gpfs/commons/groups/gursoy_lab/aelhussein/blockchain/data/preprocessing
BENCHMARKS='/gpfs/commons/groups/gursoy_lab/ubaymuradov/benchmarks'

#SAMPLE NUMBERS
PPL=4000
PPL_GS=2000
PPL_ES=1200
PPL_AF=400
PPL_BC=400

python buildChain.py -cn=$CN --datadir=$DR

#CLINICAL

python3 createStream-OMOP-Domain.py -cn=$CN --datadir=$DR -hp=$PT/data/mapping_vocab/CONCEPT_ANCESTOR.csv -dp=$PT/data/clinical/

python3 insertData-OMOP-Domain.py -cn=$CN --datadir=$DR -hp=$PT/data/mapping_vocab/CONCEPT_ANCESTOR.csv -dp=$PT/data/clinical/ -np=$PPL -mf=$PT/data/samples/metadata.csv

python3 createStream-OMOP-Person.py -cn=$CN --datadir=$DR

python3 insertData-OMOP-Person.py -cn=$CN --datadir=$DR -dp=$PT/data/clinical/ --personPath=$PT/data/samples/metadata.csv -np=$PPL -mf=$PT/data/samples/metadata.csv

#GENETIC
python3 createStream-variants.py -cn=$CN --datadir=$DR

##WGS
python3 insertData-variant.py -cn=$CN --datadir=$DR -dp=$PT/data/vcf/wgs -mf=$PT/data/samples/metadata.csv -np=$PPL_GS -sq=WGS

python3 insertData-variantPerson.py -cn=$CN --datadir=$DR -dp=$PT/data/vcf/wgs -mf=$PT/data/samples/metadata.csv -np=$PPL_GS -sq=WGS

##WES
python3 insertData-variant.py -cn=$CN --datadir=$DR -dp=$PT/data/vcf/wes -mf=$PT/data/samples/metadata.csv -np=$PPL_ES -sq=WES

python3 insertData-variantPerson.py -cn=$CN --datadir=$DR -dp=$PT/data/vcf/wes -mf=$PT/data/samples/metadata.csv -np=$PPL_ES -sq=WES

##AFFYMETRIX
python3 insertData-variant.py -cn=$CN --datadir=$DR -dp=$PT/data/vcf/affymetrix -mf=$PT/data/samples/metadata.csv -np=$PPL_AF -sq=Affymetrix

python3 insertData-variantPerson.py -cn=$CN --datadir=$DR -dp=$PT/data/vcf/affymetrix -mf=$PT/data/samples/metadata.csv -np=$PPL_AF -sq=Affymetrix

##BEADCHIP
python3 insertData-variant.py -cn=$CN --datadir=$DR -dp=$PT/data/vcf/beadchip -mf=$PT/data/samples/metadata.csv -np=$PPL_BC -sq=BeadChip

python3 insertData-variantPerson.py -cn=$CN --datadir=$DR -dp=$PT/data/vcf/beadchip -mf=$PT/data/samples/metadata.csv -np=$PPL_BC -sq=BeadChip

##Analysis
python3 insertData-analysis.py -cn=$CN --datadir=$DR -dp=$PT/data/vcf/wgs -mf=$PT/data/samples/metadata.csv -pc=$AT/ -rf=$AT/relatedness -np=$PPL_GS

##GTF
python3 createStream-gtf.py -cn=$CN -dr=$DR

python3 insertData-gtf.py -cn=$CN -dr=$DR -gp=$PT/data/gtf -vp=$PT/data/vcf/wgs -ap=$PT/data/annotations

du -h $DR/$CN >$BENCHMARKS/single-node-space-usage-stats.txt
multichain-cli $CN -datadir=$DR liststreams | jq -c '.[] | {name: .name, items: .items, keys: .keys, publishers: .publishers}' >$BENCHMARKS/single-node-stats.txt
