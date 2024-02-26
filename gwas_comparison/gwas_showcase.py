#!/usr/bin/env python
# coding: utf-8

""" 
This notebook compares the cohort definition and GWAS analysis 
using our blockchain implmenetation and a standard software

Key differences include:
- Simpler process to create cohort
- Integrated process to harmonize genetic and clinical data
- No requirment to use multiple different packages to handle different data types
"""

# # STANDARD IMPLEMENTATION
import warnings
warnings.simplefilter(action='ignore')
import pandas as pd
import numpy as np
import  copy

# ## BUILD COHORT - Clinical

# ### Load Clinical data
#lOAD THE OMOP TABLES
ROOT_DIR = 'ROOT_DIR'
hierarchy = pd.read_csv(f'{ROOT_DIR}/omop_concept_ancestor.csv', sep = '\t')
condition = pd.read_csv(f'{ROOT_DIR}/omop_condition_occurrence.txt', sep = '\t')
drug = pd.read_csv(f'{ROOT_DIR}/omop_drug_exposure.txt', sep = '\t')
procedure = pd.read_csv(f'{ROOT_DIR}/omop_procedure_occurrence.txt', sep = '\t')

#lOAD THE SELF_REPORT
#MUST BE DONE THIS WAY DUE TO ISSUES WITH TOKENIZATION
file_path = f'{ROOT_DIR}/anewbury/ADO/self_report_data/ukb675190.txt'
with open(file_path, 'r') as file:
    lines = []
    for line in file:
        line = line.strip().split('\t')
        if len(line) < 137:
            pad_size = 137-len(line)
            line = line + [np.nan for i in range(pad_size)]
        lines.append(line)
self_report = pd.DataFrame(lines[1:], columns = lines[0])
self_report.set_index('eid', inplace = True)

# All DM condition codes that are not gestional 
dm_cond_codes = hierarchy['descendant_concept_id'][hierarchy['ancestor_concept_id'] == 201820].unique() # We have to traverse the hierarchy. This is automatically done by the platform
gestational_dm_codes = hierarchy['descendant_concept_id'][hierarchy['ancestor_concept_id'] == 4024659].unique() 
dm_cond_codes = list(set(dm_cond_codes) - set(gestational_dm_codes)) + [201820]
dm_concept_count = condition[condition['condition_concept_id'].isin(dm_cond_codes)]['eid'].value_counts() 
dm_concept = list(dm_concept_count.index)

#All oral DM drugs that are not metformin
dm_drug_codes = hierarchy['descendant_concept_id'][hierarchy['ancestor_concept_id'] == 21600744].unique() # We have to traverse the hierarchy. This is automatically done by the platform
metformin_codes = hierarchy['descendant_concept_id'][hierarchy['ancestor_concept_id'] == 21600745].unique()
dm_drug_codes = list(set(dm_drug_codes) - set(metformin_codes)) + [21600744]
dm_drug_count = drug[drug['drug_concept_id'].isin(dm_drug_codes)]['eid'].value_counts() 
dm_drugs = list(dm_drug_count.index)

# Self report
values_of_interest = ['1220', '1223', '1222']
self_report_counts = self_report.isin(values_of_interest).sum(axis = 1)
filtered_ids = self_report_counts[self_report_counts >= 1].index.tolist()
dm_self_report = [int(i) for i in filtered_ids]

# Get DM ids
dm_ids = list(set(dm_concept + dm_drugs+ dm_self_report))

# All CADD codes that is not congenital or radiation
cadd_cond_codes = hierarchy['descendant_concept_id'][hierarchy['ancestor_concept_id'] == 317576].unique()
exc_cadd_codes = hierarchy['descendant_concept_id'][hierarchy['ancestor_concept_id'].isin([4178321, 4175846, 4119951])].unique()
cadd_cond_codes = list(set(cadd_cond_codes) - set(exc_cadd_codes)) + [317576]
cadd_concept_count = condition[condition['condition_concept_id'].isin(cadd_cond_codes)]['eid'].value_counts() 
cadd_concept = list(cadd_concept_count.index)

# All CADD procedure
cadd_proc_codes = hierarchy['descendant_concept_id'][hierarchy['ancestor_concept_id'] == 4336464].unique()
cadd_proc_codes = list(set(cadd_proc_codes)) + [4336464]
cadd_proc_count = procedure[procedure['procedure_concept_id'].isin(cadd_proc_codes)]['eid'].value_counts() 
cadd_proc = list(cadd_proc_count.index)

# Self report
values_of_interest = ['1074','1075']
self_report_counts = self_report.isin(values_of_interest).sum(axis = 1)
filtered_ids = self_report_counts[self_report_counts >= 1].index.tolist()
cadd_self_report = [int(i) for i in filtered_ids]

#CADD ids
cadd_ids = set(cadd_concept + cadd_proc + cadd_self_report)
cadd_ids = list(cadd_ids.intersection(set(dm_ids)))

#Create phenotype table
dm_df = pd.DataFrame(dm_ids, columns=['eid']) 
dm_df['phenotype'] = np.where(dm_df['eid'].isin(cadd_ids), 2, 1)
phenotype = dm_df

#Extract demographic information
person = pd.read_csv(f'{ROOT_DIR}/omop/omop_person.txt', sep = '\t')
person_cohort = person[person['eid'].isin(phenotype['eid'])]
person_cohort = person_cohort[['eid', 'gender_concept_id', 'year_of_birth', 'month_of_birth', 'race_concept_id']]

gender_mapping = {8507:1, 8532:2} #0=male, 1=female
race_mapping = {8527:'White',38003574:'Asian_indian',38003600:'African', 38003575:'Bangladeshi', 8515:'Asian',38003589:'Pakistani', 38003579:'Chinese', 38003598:'Black'}
person_cohort['sex'] = person_cohort['gender_concept_id'].replace(gender_mapping)
person_cohort['race'] = person_cohort['race_concept_id'].replace(race_mapping)

from datetime import datetime
person_cohort['birthdate'] = pd.to_datetime(person_cohort['year_of_birth'].astype(str) + '-' + person_cohort['month_of_birth'].astype(str))
reference_date = datetime(2011, 1, 1)
person_cohort['age'] = ((reference_date - person_cohort['birthdate']).dt.days / 365.25).astype(int)

#filter for white patients only
person_cohort = person_cohort[person_cohort['race'] == 'White']

# ## BUILD COHORT - GENETIC
# ### Harmonize genetic data
# We first must identify all patients with sequencing data and the machine used
samples = pd.read_csv(f'{ROOT_DIR}/Chr1/chr1.fam', sep= ' ', names = ['eid', 'fid', '0', '1', '2', 'batch'], usecols=['eid', 'batch'])

def map_values(value):
    if value.startswith('UKBiLEVEAX'):
        return 0
    elif value.startswith('Batch_b'):
        batch_number = int(value.split('_')[1][1:]) 
        if 1 <= batch_number <= 22:
            return 1
        elif 23 <= batch_number:
            return 2
    elif value == 'redacted3':
        return 
    return -1 
samples = samples[samples['batch'].isin([0,2])]

#Intersect those that meet all criteria
samples_used = list(set(samples['eid']).intersection(set(person_cohort['eid'])))
person_cohort = person_cohort[person_cohort['eid'].isin(samples_used)]

#Create metadata table
metadata = phenotype.merge(person_cohort[['eid', 'age', 'sex']], on = 'eid')

# Wrangle metadata table for PLINK
metadata.insert(0, 'family_id', 0)
metadata['family_id'] = metadata['eid']

# ### Extract PC's
# We use preloaded ones here
pca = pd.read_csv(f'{ROOT_DIR}/principal_components.csv')
pca = pca.apply(lambda x:x.fillna(x.mean()), axis = 0)
pca.columns = ['eid', 'PC1', 'PC2', 'PC3', 'PC4']

# ### Create covariates and phenotypes tables for PLINK
#covariates
covars = metadata[['family_id','eid', 'age', 'sex']].merge(pca, on = 'eid')
covars.to_csv(f'{ROOT_DIR}/covars.txt', index = False, sep = '\t')

#phenotypes
pheno = metadata[['family_id', 'eid','phenotype']]
pheno.to_csv(f'{ROOT_DIR}/phenotypes.txt', index = False, sep = '\t', header = False)

# ### PLINK SCRIPT
# we show the script required to run PLINK which is done on an external tool
#  Note we have removed QC steps as this was done prior analysis steps
""" 
#First convert VCF files to bed file.
process() {
    file="chr_$1".vcf.gz
    bed_file=$(echo "$file" | sed 's/.*\(chr_[0-9]*\).*/\1/')
    plink --vcf "$ROOT_DIR/$file" --make-bed --out "$ROOT_DIR/$bed_file" --threads 2
}
export -f process

parallel -j $SLURM_CPUS_PER_TASK process ::: $(seq 1 22)

# Merge the datasets
plink --merge-list $ROOT_DIR/merge.txt --out $ROOT_DIR/merged

# Run GWAS
plink --bfile $ROOT_DIR/merged --covar $ROOT_DIR/covars.txt --glm  --out $ROOT_DIR/results 

#Filter results for the variant data only
files="${ROOT_DIR}/results.PHENO1.glm.logistic"
awk_cmd="awk -F'\t' '\$7 == \"ADD\"' ${files} > ${files}_filtered"
sed_cmd="sed -i '1s/^/#CHROM\\tPOS\\tID\\tREF\\tALT\\tA1\\tTEST\\tOBS_CT\\tBETA\\tSE\\tT_STAT\\tP\\n/' ${files}_filtered"
eval $awk_cmd
eval $sed_cmd

"""

# # BLOCKCHAIN IMPLEMENTATION
ROOT_DIR = 'ROOT_DIR'
multichainLoc = ''
chainName = 'CHAIN'
datadir = f'{ROOT_DIR}/multichain'
querydir = f'{ROOT_DIR}/code/chain_code'

# Standard libaries
import pandas as pd
import json
import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf
import sys
sys.path.append(f'{querydir}')
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
#Network functions
from QueryParse import (
                        harmonizeMetadata,
                        getPCs,
                        getAgeGenderRace,
                        getPhenotype,
                        getVariantDF2,
                        runGwas2
                          )

# ## BUILD COHORT
# ### Harmonize genetic data
#Search for patients with Affymetrix
metadata = 'Affymetrix'
meta_ids = harmonizeMetadata(metadata)

# ### EXTRACT PC'S
kSearch = 20
pc_df = getPCs(meta_ids, kSearch)

# ### Get phenotypes
# ####  Age and gender
# Get the age and gender and filter for patients who are race = white
demos = getAgeGenderRace(chainName, multichainLoc, datadir)
demos = demos[demos['race_concept_id'] == '8527']

# #### Phenotype of interest
## INCLUSION ##
#DM diagnosis
dm_cond = '201820'
dm_id_cond =  getPhenotype(dm_cond)

#DM medication
dm_drug = '21600744'
dm_id_drug=  getPhenotype(dm_drug)

#DM self reprot
dm_self = '123456789'
dm_id_self=  getPhenotype(dm_self)

## EXCLUSION ##
# Gestational DM
gest_dm_cond = '4024659'
gest_dm_id_cond =  getPhenotype(gest_dm_cond)

# Metformin
dm_drug_met = '21600745'
dm_id_met=  getPhenotype(dm_drug_met)

## CADD PHENOTYPE ##
# CADD diagnosis
cadd_cond = '317576'
cadd_id_cond =  getPhenotype(cadd_cond)

# CADD Procedure
cadd_proc = '4336464'
cadd_id_proc =  getPhenotype(cadd_proc)

#CADD self report
cadd_self = '23456789'
cadd_id_self=  getPhenotype(cadd_self)

# Phenotyping logic
all_dm_ids = list(set(list(dm_id_drug['person_id']) + list(dm_id_cond['person_id']) + list(dm_id_self['person_id'])) - set(list(gest_dm_id_cond['person_id']) + list(dm_id_met['person_id'])))
all_cadd_ids = list(set(list(cadd_id_cond['person_id']) + list(cadd_id_proc['person_id']) + list(cadd_id_self['person_id'])))
only_dm_ids = list(set(all_dm_ids) - set(all_cadd_ids))
phenotype = pd.DataFrame(all_dm_ids, columns=['eid']) 
phenotype['phenotype'] = np.where(phenotype['eid'].isin(all_cadd_ids), 2, 1)
phenos = phenotype.set_index('eid')

# ### EXTRACT GENOTYPE INFORMATION
chroms = [i for i in range(1,23)]
variants_df = pd.DataFrame()
for chrom in chroms:
    genotype = 'all'
    metadata = None
    variants_df_chrom = getVariantDF2(chrom, variants = 'all', genotype = 'all', metadata = None)
    variants_df  = pd.concat([variants_df, variants_df_chrom],  axis = 1)

# ### RUN GWAS
results = runGwas2(pc_df, phenos, demos, variants_df)


