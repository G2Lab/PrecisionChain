#! /bin/bash
set -e
CN=query_4k
DR=/mnt/disks/data/multichain_data
PT=/mnt/disks/data/test-data/data
AT=/mnt/disks/data/test-data/preprocessing

#SAMPLE NUMBERS
PPL=100
PPL_GS=50
PPL_ES=30
PPL_AF=10
PPL_BC=10

python3 buildChain.py -cn=$CN --datadir=$DR

#CLINICAL

python3 createStream-OMOP-Domain.py -cn=$CN --datadir=$DR -hp=$PT/mapping_vocab/CONCEPT_ANCESTOR.csv -dp=$PT/clinical/

python3 insertData-OMOP-Domain.py -cn=$CN --datadir=$DR -hp=$PT/mapping_vocab/CONCEPT_ANCESTOR.csv -dp=$PT/clinical/ -np=$PPL -mf=$PT/samples/metadata.csv

python3 createStream-OMOP-Person.py -cn=$CN --datadir=$DR

python3 insertData-OMOP-Person.py -cn=$CN --datadir=$DR -dp=$PT/clinical/ --personPath=$PT/clinical/person.csv -np=$PPL -mf=$PT/samples/metadata.csv

#GENETIC
python3 createStream-variants.py -cn=$CN --datadir=$DR

##WGS
python3 insertData-variant.py -cn=$CN --datadir=$DR -dp=$PT/vcf/wgs -mf=$PT/samples/metadata.csv -np=$PPL_GS -sq=WGS

python3 insertData-variantPerson.py -cn=$CN --datadir=$DR -dp=$PT/vcf/wgs -mf=$PT/samples/metadata.csv -np=$PPL_GS -sq=WGS

##WES
python3 insertData-variant.py -cn=$CN --datadir=$DR -dp=$PT/vcf/wes -mf=$PT/samples/metadata.csv -np=$PPL_ES -sq=WES

python3 insertData-variantPerson.py -cn=$CN --datadir=$DR -dp=$PT/vcf/wes -mf=$PT/samples/metadata.csv -np=$PPL_ES -sq=WES

##AFFYMETRIX
python3 insertData-variant.py -cn=$CN --datadir=$DR -dp=$PT/vcf/affymetrix -mf=$PT/samples/metadata.csv -np=$PPL_AF -sq=Affymetrix

python3 insertData-variantPerson.py -cn=$CN --datadir=$DR -dp=$PT/vcf/affymetrix -mf=$PT/samples/metadata.csv -np=$PPL_AF -sq=Affymetrix

##BEADCHIP
python3 insertData-variant.py -cn=$CN --datadir=$DR -dp=$PT/vcf/beadchip -mf=$PT/samples/metadata.csv -np=$PPL_BC -sq=BeadChip

python3 insertData-variantPerson.py -cn=$CN --datadir=$DR -dp=$PT/vcf/beadchip -mf=$PT/samples/metadata.csv -np=$PPL_BC -sq=BeadChip

##Analysis
python3 insertData-analysis.py -cn=$CN --datadir=$DR -dp=$PT/vcf/wgs -mf=$PT/samples/metadata.csv -pc=$AT/ -rf=$AT/relatedness -np=$PPL_GS

##GTF
python3 createStream-gtf.py -cn=$CN -dr=$DR

python3 insertData-gtf.py -cn=$CN -dr=$DR -gp=$PT/gtf -vp=$PT/vcf/wgs -ap=$PT/annotations
