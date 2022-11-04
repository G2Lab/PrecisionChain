#! /bin/bash

python3 buildChain.py -cn=combchain -dr=/mnt/disks/data/multichain_data/

python3 createStream-OMOP-Domain.py -cn=combchain --datadir=/mnt/disks/data/multichain_data/ -hp=/home/bek/PrecisionChain/data/mapping_vocab/CONCEPT_ANCESTOR.csv -dp=/home/bek/PrecisionChain/data/clinical/

python3 insertData-OMOP-Domain.py -cn=combchain --datadir=/mnt/disks/data/multichain_data/ -hp=/home/bek/PrecisionChain/data/mapping_vocab/CONCEPT_ANCESTOR.csv -dp=/home/bek/PrecisionChain/data/clinical/

python3 createStream-OMOP-Person.py -cn=combchain --datadir=/mnt/disks/data/multichain_data/

python3 insertData-OMOP-Person.py -cn=combchain --datadir=/mnt/disks/data/multichain_data/ -dp=/home/bek/PrecisionChain/data/clinical/ --personPath=/home/bek/PrecisionChain/data/clinical/person.csv

python3 createStream-variants.py -cn=combchain --datadir=/mnt/disks/data/multichain_data/

python3 insertData-variant.py -cn=combchain --datadir=/mnt/disks/data/multichain_data/ -mf=/home/bek/PrecisionChain/data/mapping_vocab/person_sample_mapping.txt -dp=/home/bek/PrecisionChain/data/vcf

python3 insertData-variantPerson.py -cn=combchain --datadir=/mnt/disks/data/multichain_data/ -mf=/home/bek/PrecisionChain/data/mapping_vocab/person_sample_mapping.txt -dp=/home/bek/PrecisionChain/data/vcf

python3 createStream-gtf.py -cn=combchain -dr=/mnt/disks/data/multichain_data/

python3 insertData-gtf.py -cn=combchain -dr=/mnt/disks/data/multichain_data/ -gp=/home/bek/PrecisionChain/data/gtf/ -vp=/home/bek/PrecisionChain/data/vcf/
