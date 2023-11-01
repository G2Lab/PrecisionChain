#!/usr/bin/env python
# coding: utf-8

# In[2]:


# %load /gpfs/commons/groups/gursoy_lab/aelhussein/Code/SAMChain/SAMchain/buildChain.py
"""
Helper functions to abstract parsing of returned data
"""

ROOT_DIR = '/gpfs/commons/groups/gursoy_lab/aelhussein/blockchain'
multichainLoc = ''
chainName = 'public_access_2'
datadir = f'{ROOT_DIR}/multichain'
querydir = f'{ROOT_DIR}/public/code'
metafile = f'{ROOT_DIR}/public/data/samples/metadata.csv'
annotation_path = f'{ROOT_DIR}/public/data/annotations'
personPath = f'{ROOT_DIR}/public/data/clinical/person.csv'
dataPath = f'{ROOT_DIR}/public/data/clinical/'

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
from QueryClinical import ( 
                            queryGroupDemographics,
                            domainQuery,
                            personQuery,
                            parseKeys
                          )

from QueryVariant import ( 
                            queryVariants,
                            queryPersonsChroms,
                            getPatientVariantAnnotation
                         )

from QueryCombination import (
                                extractGeneVariants,
                                queryClinicalGeneVariantRange,
                                queryVariantClinical,
                             )

from QueryAnalysis import ( 
                            queryMetadata, 
                            queryKinship,
                            querySamplePCA
                          )


def harmonizeMetadata(metadata):
    """ Returns list of harmonized IDs """
    response = queryMetadata(chainName, datadir, metadata)
    metadata_list = metadata.split(',')
    ids = [value for meta in metadata_list for value in list(response[meta].values())[0]]
    meta_ids = list(set(ids))
    print(f'{len(meta_ids)} patients meet sequencing metadata criteria')
    return meta_ids

def removeRelated(ids):
    """ Checks and removes related patients """
    if isinstance(ids, list):
        ids = ','.join(ids)
    response = queryKinship(chainName, datadir, ids)
    response_json = json.loads(response)
    kin_df = pd.DataFrame(response_json)
    unrelated = kin_df.apply(lambda col: (col == 'UR').sum() == kin_df.shape[0] - 1)
    unrelated_ids = unrelated.index.tolist()
    print(f'{len(unrelated_ids)} remaining after removing related samples')
    return unrelated_ids

def getPCs(ids, kSearch):
    """ Get PCs for list of samples """
    kSearch = 20
    if isinstance(ids, list):
        ids = ','.join(ids)
    response = querySamplePCA(chainName, datadir, ids, kSearch)
    pc_df = pd.DataFrame(json.loads(response))
    return pc_df

def getAgeGender(chainName, multichainLoc, datadir, ids):
    """ extract the age and gender of all patients """
    searchKeys = 'birth_datetime,gender_concept_id'
    demo_data = queryGroupDemographics(chainName, multichainLoc, datadir, searchKeys)
    demographics = pd.DataFrame(demo_data).T
    keys = searchKeys.split(',')
    demographics.columns = keys

    demographics['birth_datetime'] = pd.to_datetime(demographics['birth_datetime'])
    current_date = pd.to_datetime('2023-07-06')
    demographics['age'] = (current_date - demographics['birth_datetime']).dt.days / 365.25
    demo_processed = demographics[['gender_concept_id', 'age']]
    demo_processed['gender'] = demo_processed['gender_concept_id'].replace({8507:0, 8532:1})
    demo_processed.drop(columns = 'gender_concept_id', inplace = True)

    if isinstance(ids,str):
        ids = ids.split(',')

    return demo_processed.loc[ids]

def getPhenotype(pheno_id, demos):
    """ Get phenotype of interest """
    searchKeys = 'demographics' # returns basic information for patients with disease. can be changed if more complex info needed
    response = domainQuery(chainName, multichainLoc, datadir, pheno_id, searchKeys)
    data = [r['data']['json'] for r in response]
    df = pd.DataFrame(data)
    disease_ids = list(df['person_id'].unique())
    phenos = demos.copy()
    phenos['phenotype'] = 0
    phenos.loc[phenos.index.isin(disease_ids), 'phenotype']  = 1
    return phenos

def getVariantDF(chrom, variants, genotype = 'all', metadata = None):
    """ Get variant in DF format """
    response = queryVariants(chainName, multichainLoc, datadir, chrom, variants, genotype, metadata)
    variants_dict = json.loads(response)
    data_for_df = []
    for variant, genotypes in variants_dict.items():
        for genotype, ids in genotypes.items():
            for id_ in ids:
                data_for_df.append({'variant': variant, 'genotype': genotype, 'id': id_})
    df = pd.DataFrame(data_for_df)
    variants_df = df.pivot(index='id', columns='variant', values='genotype').reset_index()
    variants_df.set_index('id', inplace=True)
    variants_df.replace({"0|0":0, "1|0":1, "1|1":1}, inplace = True)
    variants_df.columns = ['variant']
    variants_df.index = variants_df.index.astype(str)
    return variants_df

def runGwas(pc_df, phenos, variants_df):
    #Linear mixed model with age, gender and phenotype
    covariates = pc_df.merge(phenos, left_index =True, right_index=True)
    data = covariates.merge(variants_df, left_index =True, right_index=True)
    formula = f"variant ~ phenotype + age + gender + " + ' + '.join([f'V{i}' for i in range(1,kSearch+1)])
    md = smf.ols(formula, data)
    mdf = md.fit()
    return mdf