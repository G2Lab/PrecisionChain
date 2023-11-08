# PrecisionChain
This README accompanies the paper '_Realizing the potential of secure and decentralized harmonization of clinical and genetic data for precision medicine'_. The contribution of the paper is to create a decentralized blockchain based data-sharing framework that unifies storage and access of clinical and genetic data for precision medicine purposes. The README is a step by step guide to replicate the publicly available network available at  https://precisionchain.g2lab.org/. Note this code supports the set-up, insertion and querying functionality but not the front-end site.

Previews of the data are available in the data folder. For the complete data please following te links. If changing the structure of the directories, please be sure to update the directory paths in the scripts. For clinical data we use Synthea OMOP CDM format accessed via: (https://registry.opendata.aws/synthea-omop/) (Walonoski et al., JAMIA, 2018). For genetic data we simulate data using 1000 GenomesProject in VCF format (https://www.internationalgenome.org/category/vcf/). A link to the code used to generate the simulated data can be found at https://github.com/aelhussein/blockchain. 


## Requirements:
We use the Multichain blockchain API for the platform and bcftools to process VCF files. We use v2.3.2.
- Multichain Blockchain API (https://www.multichain.com/getting-started/)
To set up Multichain, please download and install MultiChain Community. Once installed, keep a note of the directory in which MultiChain will store blockchain data. Note users do not need to create a chain, chain creation is automated by the scripts. However, before continuing we reccommend creating a dummy chain to familiarize yourself with setup.

- bcftools (https://samtools.github.io/bcftools/bcftools.html)
This is a tool developed by Samtools to support the processing of VCF files. To install please follow the instructions provided at: https://samtools.github.io/bcftools/howtos/install.html. Note, scripts will automatically call bcftools during genetic data insertion. Users are not expected to use the tool themselves.
- Python version 3.6.13. Please see requirements.txt file for the packages used.

## Scripts: Chain creation and insertion
This is a list of scripts to call for chain creation and data insertion. Please follow chain creation and data insertion scripts in the specified order.

#### Arguments used in scripts

- [CHAIN NAME] = Name of the chain <br/>
- [MULTICHAIN DIR] = Directory where multichain data is being stored in <br/>
- [CONCEPT HIERARCHY DIR] = Directory where OMOP clinical concepts hierarch is stored (provided in same directory as clinical data) <br/>
- [CLINICAL DATA DIR] = Directory where OMOP formatted clinical data is stored <br/>
- [PATIENT TABLE DIR] = File where patient demographic data is stored (provided in same directory as clinical data) <br/>
- [PATIENT:GENETIC SAMPLE MAPPING FILE DIR] =  Dictionary with a mapping between Patient IDs and VCF sample IDs (this is necessary if your lab samples have different ids to your hospital records) <br/>
- [VCF FILE DIR] = Directory where VCF files are stored <br/>
- [GTF FILE DIR] = Directory where GTF files are stored <br/>
- [METADATA FILE DIR] = Directory where VCF files are stored <br/>
- [PRINCIPLE COMPONENTS FILE DIR] = Directory where GTF files are stored <br/>
- [RELATEDNESS SNPS FILE DIR] = Directory where GTF files are stored <br/>
- [PPL] = Number of people to add <br/>
- [SEQUENCING] = Genetic sample sequencing type <br/>

#### Builds the empty chain
```
python buildChain.py -cn=[CHAIN NAME] -dr=[MULTICHAIN DIR]
```
#### Creates concept streams for clinical concepts in OMOP
```
python createStream-OMOP-Domain.py -cn=[CHAIN NAME] --datadir=[MULTICHAIN DIR] -hp=[CONCEPT HIERARCHY DIR] -dp=[CLINICAL DATA DIR]
```
#### Inserts OMOP concepts in domain-view
```
python insertData-OMOP-Domain.py -cn=[CHAIN NAME] --datadir=[MULTICHAIN DIR] -hp=[CONCEPT HIERARCHY DIR] -dp=[CLINICAL DATA DIR] -np=[PPL]
```
#### Creates person streams for clinical concepts 
```
python createStream-OMOP-Person.py -cn=[CHAIN NAME] --datadir=[MULTICHAIN DIR]
```
#### Inserts OMOP concepts in person-view
```
python insertData-OMOP-Person.py -cn=[CHAIN NAME] --datadir=[MULTICHAIN DIR] -dp=[CLINICAL DATA DIR] --personPath=[PATIENT TABLE DIR] -np=[PPL]
```
#### Creates chromosome streams for VCF files
```
python createStream-variants.py -cn=[CHAIN NAME] --datadir=[MULTICHAIN DIR]
```
#### Inserts variant data from VCF files, inserted per position
```
python insertData-variant.py -cn=[CHAIN NAME] --datadir=[MULTICHAIN DIR] -mf=[PATIENT:GENETIC SAMPLE MAPPING FILE DIR] -dp=[VCF FILE DIR] -np=[PPL] -sq=[SEQUENCING]
```
#### Inserts variant data from VCF files, inserted per sample
```
python insertData-variantPerson.py -cn=[CHAIN NAME] --datadir=[MULTICHAIN DIR] -mf=[PATIENT:GENETIC SAMPLE MAPPING FILE DIR] -dp=[VCF FILE DIR] -np=[PPL] -sq=[SEQUENCING]
```
#### Inserts analysis data
```
python insertData-analysis.pycn=[CHAIN NAME] --datadir=[MULTICHAIN DIR] -dp=[VCF FILE DIR] -mf=[METADATA FILE DIR] -pc=[PRINCIPLE COMPONENTS FILE DIR] -rf=[RELATEDNESS SNPS FILE DIR] -np=[PPL]
```
#### Creates chromosome streams for GTF files
```
python createStream-gtf.py -cn=[CHAIN NAME] -dr=[MULTICHAIN DIR]
```
#### Inserts genetic data from GTF files
```
python insertData-gtf.py -cn=[CHAIN NAME] -dr=[MULTICHAIN DIR] -gp=[GTF FILE DIR] -vp=[VCF FILE DIR]
```
## Scripts: Querying

#### Arguments used in scripts
- [CHAIN NAME] = Name of the chain <br/>
- [MULTICHAIN DIR] = Directory where multichain data is being stored in <br/>
- [COHORT KEYS] = OMOP concept codes used to define the cohort e.g. OMOP code 201826 to select patients with Type II diabetes diagnosis. Multiple keys can be provided. Please separate keys with ','. <br/>
- [SEARCH KEYS] = OMOP concept codes to extract within the cohort e.g. OMOP code 44790340 extract all drugs taken by cohort. Multiple keys can be provided. Please separate keys with ','. <br/>
- [VIEW_CLINICAL] QueryClinical = Option of domain or person. Domain view takes a cohortKey and returns data for that cohort. Person view extracts all clinical data for a given set of patients. Both results can be filtered by SearchKey which filters the clinical information returned. <br/>
- [VIEW_GENETIC]QueryVariant = Option of variant, person, MAF or annotation. Variant view extracts all data for a given set of positions. Person view extracts all genotypes for a given set of patients. MAF view extracts all variants within a certain MAF range. Annotation view extracts relevant annotations for a list of variants. <br/>
- [VIEW_COMBINATION] QueryCombination = Option of variant or clinical. Variant view extracts clinical information for patients with the searched variant (position and genotype).Clinical view extracts variant information in a given gene for patients in a specified clinical cohort <br/>
- [VIEW_ANALYSIS] QueryAnalysis = Option of pca, kin or meta. pca view extracts extracts the principle components for the samples, kin view extracts information on whether samples are related and meta view extracts technical sequencing metadata <br/>
- [CHROMOSOMES] = Chromsomes to search. Multiple chromosomes can be provided. Please separate with ','. Only necessary if [VIEW] = Variant OR Person <br/>
- [GENOTYPES] = Genotypes to extract from each variant i.e. '0/0', '1/0', '1/1'. Only necessary if [VIEW] = Variant <br/>
- [PERSON_IDS] = Person IDs to search. Multiple IDs can be provided. Please separate with ','. Only necessary if [VIEW] = Person <br/>
- [INPUT RANGE] = MAF range to search. Input values between 0-1 in 'X-Y' format. Only necessary if [VIEW] = MAF <br/>
- [SAMPLE SEARCH] = Person IDs to extract data for, default is all samples. Only necessary if [VIEW] = pca or kin<br/>
- [K SEARCH] = Number of PCs to extract per user min is 1 and max is 20. Only necessary if [VIEW] = pca <br/>
- [METADATA] = Sequencing metadata to search for. Only necessary if [VIEW] = meta <br/>
- [ANNOTATIONS] = Annotation type to extract can be any or vepp, clinvar or cadd. Only necessary if [VIEW] = annot <br/>


#### Query clinical data
```
python QueryClinical.py -cn=[CHAIN NAME] -dr=[MULTICHAIN DIR] --view=[VIEW_CLINICAL] -ck=[COHORT KEYS] -sk=[SEARCH KEY] -pi=[PERSON_IDS]
```
#### Query genetic data
```
python QueryVariant.py -cn=[CHAIN NAME] -dr=[MULTICHAIN DIR] --view=[VIEW_GENETIC] -ch=[CHROMOSOMES] -ps=[POSITIONS]  -gt=[GENOTYPES] -pi=[PERSON_IDS] -ir=[INPUT RANGE] -md=[METADATA] -at=[ANNOTATIONS]
```
#### Query combination of clinical and genetic data
```
python QueryCombination.py -cn=[CHAIN NAME] -dr=[MULTICHAIN DIR] --view=[VIEW_COMBINATION] -ch=[CHROMOSOMES] -gn=[GENE]  -ir=[INPUT RANGE] -ck=[COHORT KEYS] -gt=[GENOTYPES]
```
#### Query analysis data
```
python QueryAnalysis.py -cn=[CHAIN NAME] -dr=[MULTICHAIN DIR] --view=[VIEW_ANALYSIS] -ss=[SAMPLE SEARCH] -ks=[K SEARCH]  -md=[METADATA]
```



