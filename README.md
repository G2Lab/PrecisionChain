# PrecisionChain
This ReadMe accompanies the paper Realizing the potential of secure and decentralized harmonization of clinical and genetic data for precision medicine. The contribution of the paper is to create a decentralized blockchain based data-sharing framework that unifies storage and access clinical and genetic data for precision medicine purposes. The ReadMe is a step by step guide to replicate the publicly available network available at LINK. Note this code supports the set-up, insertion and querying functionality but not the front-end site.

All data is available in the data folder. If changing the structure of the directories, please be sure to update the directory paths in the scripts. For clinical data we use MIMIC-IV in OMOP CDM format (https://physionet.org/content/mimic-iv-demo-omop/0.9/). For genetic data we use a sample of 100 patients from the 1000 GenomesProject in VCF format (https://www.internationalgenome.org/category/vcf/). 



**Requirements:**
We use the Multichain blockchain API for the platform and bcftools to process VCF files.
- Multichain Blockchain API (https://www.multichain.com/getting-started/)
To set up Multichain, please download and install MultiChain Community. Once installed, keep a note of the directory in which MultiChain will store blockchain data. Note users do not need to create a chain, chain creation is automated by the scripts. However, before continuing we reccommend creating a dummy chain to familiarize yourself with setup.

- bcftools (https://samtools.github.io/bcftools/bcftools.html)
This is a tool developed by Samtools to support the processing of VCF files. To install please follow the instructions provided at: https://samtools.github.io/bcftools/howtos/install.html. Note, scripts will automatically call bcftools during genetic data insertion. Users are not expected to use the tool themselves.

**Scripts: Chain creation and insertion**
This is a list of scripts to call for chain creation and data insertion. Please follow chain creation and data insertion scripts in the specified order.

## Arguments used in scripts
[CHAIN NAME] = Name of the chain
[MULTICHAIN DIR] = Directory where multichain data is being stored in
[CONCEPT HIERARCHY DIR] = Directory where OMOP clinical concepts hierarch is stored (provided in same directory as clinical data)
[CLINICAL DATA DIR] = Directory where OMOP formatted clinical data is stored
[PATIENT TABLE DIR] = File where patient demographic data is stored (provided in same directory as clinical data)
[PATIENT:GENETIC SAMPLE MAPPING FILE DIR] =  Dictionary with a mapping between Patient IDs and VCF sample IDs (this is necessary if your lab samples have different ids to your hospital records)
[VCF FILE DIR] = Directory where VCF files are stored
[GTF FILE DIR] = Directory where GTF files are stored


#### Builds the empty chain
python buildChain.py -cn=[CHAIN NAME] -dr=[MULTICHAIN DIR]

#### Creates concept streams for clinical concepts in OMOP
python createStream-OMOP-Domain.py -cn=[CHAIN NAME] --datadir=[MULTICHAIN DIR] -hp=[CONCEPT HIERARCHY DIR] -dp=[CLINICAL DATA DIR]

#### Inserts OMOP concepts in domain-view
python insertData-OMOP-Domain.py -cn=[CHAIN NAME] --datadir=[MULTICHAIN DIR] -hp=[CONCEPT HIERARCHY DIR] -dp=[CLINICAL DATA DIR]

#### Creates person streams for clinical concepts 

python createStream-OMOP-Person.py -cn=[CHAIN NAME] --datadir=[MULTICHAIN DIR]

#### Inserts OMOP concepts in person-view

python insertData-OMOP-Person.py -cn=[CHAIN NAME] --datadir=[MULTICHAIN DIR] -dp=[CLINICAL DATA DIR] --personPath=[PATIENT TABLE DIR]

#### Creates chromosome streams for VCF files

python createStream-variants.py -cn=[CHAIN NAME] --datadir=[MULTICHAIN DIR]

#### Inserts variant data from VCF files, inserted per position

python insertData-variant.py -cn=[CHAIN NAME] --datadir=[MULTICHAIN DIR] -mf=[PATIENT:GENETIC SAMPLE MAPPING FILE DIR] -dp=[VCF FILE DIR]

#### Inserts variant data from VCF files, inserted per sample

python insertData-variantPerson.py -cn=[CHAIN NAME] --datadir=[MULTICHAIN DIR] -mf=[PATIENT:GENETIC SAMPLE MAPPING FILE DIR] -dp=[VCF FILE DIR]

#### Creates chromosome streams for GTF files

python createStream-gtf.py -cn=[CHAIN NAME] -dr=[MULTICHAIN DIR]

#### Inserts genetic data from GTF files

python insertData-gtf.py -cn=[CHAIN NAME] -dr=[MULTICHAIN DIR] -gp=[GTF FILE DIR] -vp=[VCF FILE DIR]

**Scripts: Querying**






