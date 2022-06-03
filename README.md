# multidomain_chain
software to automate the set-up, insertion and querying of clinical and genetic data. For clinical data we use OMOP CDM. We provide 2 views domain view and person view. In domain view the data is structured as in an SQL table with each entry representing a row. In the Person view, data pertaining to a person is aggregated and inserted together, much like a noSQL view. This makes it easier to get complete 360 views of a patient.

**Requirements:**
- multichain blockchain API (
https://www.multichain.com/getting-started/)
- bcftools (https://samtools.github.io/bcftools/bcftools.html)

**Script:**
##Builds the empty chain
python buildChain.py -cn=combchain -dr=[MULTICHAIN DIR]

##Creates concept streams for clinical concepts in OMOP
python createStream-OMOP-Domain.py -cn=combchain --datadir=[MULTICHAIN DIR] -hp=[CONCEPT HIERARCHY DIR] -dp=[CLINICAL DATA DIR]

##Inserts OMOP concepts in domain-view
python insertData-OMOP-Domain.py -cn=combchain --datadir=[MULTICHAIN DIR] -hp=[CONCEPT HIERARCHY DIR] -dp=[CLINICAL DATA DIR]

##Creates person streams for clinical concepts 
python createStream-OMOP-Person.py -cn=combchain --datadir=[MULTICHAIN DIR]

##Inserts OMOP concepts in person-view
python insertData-OMOP-Person.py -cn=combchain --datadir=[MULTICHAIN DIR] -dp=[CLINICAL DATA DIR] --personPath=[PATIENT TABLE DIR]

##Creates chromosome streams for VCF files
python createStream-variants.py -cn=combchain --datadir=[MULTICHAIN DIR]

##Inserts variant data from VCF files, inserted per position
python insertData-variant.py -cn=combchain --datadir=[MULTICHAIN DIR] -mf=[PATIENT:GENETIC SAMPLE MAPPING FILE DIR] -dp=[VCF FILE DIR]

##Inserts variant data from VCF files, inserted per sample
python insertData-variantPerson.py -cn=combchain --datadir=[MULTICHAIN DIR] -mf=[PATIENT:GENETIC SAMPLE MAPPING FILE DIR] -dp=[VCF FILE DIR]

##Creates chromosome streams for GTF files
python createStream-gtf.py -cn=combchain -dr=[MULTICHAIN DIR]

##Insersts genetic data from GTF files
python insertData-gtf.py -cn=combchain -dr=[MULTICHAIN DIR] -gp=[GTF FILE DIR] -vp=[VCF FILE DIR]

**Args**
[MULTICHAIN DIR] = Directory where multichain data is being stored in
[CONCEPT HIERARCHY DIR] = Directory where OMOP clinical concepts hierarch is stored (provided here too)
[CLINICAL DATA DIR] = Directory where OMOP formatted clinical data is stored
[PATIENT TABLE DIR] = File where patient demographic data is stored (this may be a table with your other clinical data)
[PATIENT:GENETIC SAMPLE MAPPING FILE DIR] =  Dictionary with a mapping between Patient IDs and VCF sample IDs (this is necessary if your lab samples have different ids to your hospital records)
[VCF FILE DIR] = Directory where VCF files are stored
[GTF FILE DIR] = Directory where GTF files are stored



