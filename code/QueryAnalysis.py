#!/usr/bin/env python
# coding: utf-8

# In[3]:


'''
QueryAnalysis.py
Queries analysis data - mostly pca and kinship
modified by AE 07/2023
'''
import sys
import time
import math
import binascii
import argparse
import subprocess
from subprocess import Popen, PIPE
import os
import psutil
import time
import multiprocessing
import glob
import pandas as pd
import json
import warnings
import multiprocessing
import numpy as np
import dask.array as da
import gc
from datetime import datetime
from scipy.stats import multivariate_normal

warnings.simplefilter(action='ignore')


# In[ ]:


#Given a chain name subscribe to audit log stream to ensure query is recorded
def subscribeToStream(chainName, multichainLoc, datadir):
    subscribeStreamCommand=multichainLoc+'multichain-cli {} -datadir={} subscribe audit_log'.format(chainName,datadir)
    #subscribe to the stream
    procSubscribe = subprocess.Popen(subscribeStreamCommand.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    procSubscribe.wait()
    return

def querySamplePCA(chainName, datadir, sampleSearch, kSearch):
    '''
    Query PC's of samples on chain
    Input:
        sampleSearch - sampled ids to extract
        kSearch - number of PCs to get
    '''
    #Extract data
    queryCommand = 'multichain-cli {} -datadir={} liststreamkeyitems analysis {} false 9999999'.format(chainName, datadir, "PCA")
    items = subprocess.check_output(queryCommand.split())
    matches = json.loads(items, parse_int= int)
    publishToAuditstream(chainName, datadir, queryCommand)
    #Create DF
    sample_pcs = [(match_['keys'][1], match_['data']['json']) for match_ in matches]
    pc_df = pd.DataFrame.from_dict(dict(sample_pcs)).T
    #Filter based on search
    kSearch = int(kSearch)
    if sampleSearch:
        sampleSearch = sampleSearch.split(',')
        valid_samples = pc_df.index.intersection(sampleSearch)
        pc_df = pc_df.loc[valid_samples]
        pc_df = pc_df.drop_duplicates()
    if kSearch != 20:
        pc_df = pc_df.iloc[:,:kSearch]
    print(pc_df.to_json())
    return pc_df.to_json()
# In[ ]:

def querySampleRelatedness(chainName, datadir, sampleSearch):
    '''
    Query relatedness snps from chain
    Input:
        sampleSearch - sampled ids to extract
    '''
    #Extract data
    queryCommand = 'multichain-cli {} -datadir={} liststreamkeyitems analysis {} false 99999999'.format(chainName, datadir, 'Relatedness')
    items = subprocess.check_output(queryCommand.split())
    matches = json.loads(items, parse_int= int)
    publishToAuditstream(chainName, datadir, queryCommand)
    #Create DF
    sample_snps = [(match_['keys'][1], match_['data']['json']) for match_ in matches if match_['keys'][1].upper() != 'AF' ]
    rl_df = pd.DataFrame.from_dict(dict(sample_snps)).T
    #Filter based on search
    if sampleSearch:
        sampleSearch = sampleSearch.split(',')
        valid_samples = rl_df.index.intersection(sampleSearch)
        rl_df = rl_df.loc[valid_samples]
    #AF
    queryCommand = 'multichain-cli {} -datadir={} liststreamkeyitems analysis {} false 99999999'.format(chainName, datadir, 'AF')
    items = subprocess.check_output(queryCommand.split())
    matches = json.loads(items, parse_int= int)
    txid = matches[0]['data']['txid']
    queryCommand = 'multichain-cli {} -datadir={} gettxoutdata {} 0'.format(chainName, datadir, txid)
    publishToAuditstream(chainName, datadir, queryCommand)
    items = subprocess.check_output(queryCommand.split())
    matches = json.loads(items, parse_int= int)['json']
    af = [float(match_) for match_ in matches]
    return rl_df, af

def calculate_mismatch(df):
    """ Cacluate hgmr and agmr matrices """
    #Convert to int to save space
    df_array = df.to_numpy(dtype=np.int8) 
    # Convert to Dask Array
    df_dask_array = da.from_array(df.to_numpy(dtype=np.int8), chunks=(100, 1000))  # Chunk sizes are examples

    # AGMR Calculation
    agmr_dask_matrix = da.not_equal(df_dask_array[:, None, :], df_dask_array).mean(axis=2)
    agmr_matrix = agmr_dask_matrix.compute()  # This triggers the actual computation
    agmr_df = pd.DataFrame(agmr_matrix, index=df.index, columns=df.index)
    del agmr_matrix, agmr_dask_matrix
    gc.collect()

    # HGMR Calculation
    homozygous = df_dask_array % 2 == 0
    homozygous_pairs = homozygous[:, None, :] & homozygous
    num_homozygous_pairs = homozygous_pairs.sum(axis=2)
    no_homozygous_pairs = num_homozygous_pairs == 0
    num_homozygous_pairs = da.where(no_homozygous_pairs, 1, num_homozygous_pairs)
    mismatches = da.not_equal(df_dask_array[:, None, :], df_dask_array)
    homozygous_mismatches = mismatches * homozygous_pairs
    sum_homozygous_mismatches = homozygous_mismatches.sum(axis=2)
    hgmr_matrix = da.divide(sum_homozygous_mismatches, num_homozygous_pairs)
    hgmr_matrix = da.where(num_homozygous_pairs != 0, hgmr_matrix, 0)
    hgmr_matrix_computed = hgmr_matrix.compute()
    hgmr_df = pd.DataFrame(hgmr_matrix_computed, index=df.index, columns=df.index)
    del hgmr_matrix, sum_homozygous_mismatches, num_homozygous_pairs, mismatches, homozygous_pairs, homozygous, df_array, df_dask_array
    gc.collect()
    return agmr_df, hgmr_df

def determine_relationships(agmr_df, hgmr_df, relationship_parameters):
    ''' Use MLE to estimate relationships. See original paper: 
    https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0179106'''
    agmr_array = agmr_df.values
    hgmr_array = hgmr_df.values
    relationships = np.empty(agmr_array.shape, dtype=object)

    # Calculate likelihoods for all relationships
    likelihoods = np.empty((*agmr_array.shape, len(relationship_parameters)))
    for i, (rel, params) in enumerate(relationship_parameters.items()):
        likelihoods[:, :, i] = multivariate_normal(mean=params['mean'], cov=params['cov']).pdf(np.stack([agmr_array, hgmr_array], axis=-1))
    max_indices = np.argmax(likelihoods, axis=-1)
    mapping = {i: rel for i, rel in enumerate(relationship_parameters.keys())}
    max_relationships = np.vectorize(mapping.get)(max_indices)

    relationships = pd.DataFrame(max_relationships, index=agmr_df.index, columns=agmr_df.columns)

    return relationships

def assess_kinship(rl_df):
    relationship_parameters = {
        'ID': {'mean': [0.0, 0.0], 'cov': [[0.01**2, 0.01*0.38*0.317], [0.01*0.38*0.317, 0.38**2]]},
        'PO': {'mean': [0.04, 0.02], 'cov': [[0.07**2, 0.07*1.47*0.104], [ 0.07*1.47*0.104, 1.47**2]]},
        'FS': {'mean': [0.36, 0.1], 'cov': [[1.02**2,  1.02*2.32*0.784], [1.02*2.32*0.784, 2.32**2]]},
        'D2': {'mean': [0.40, 0.15], 'cov': [[1.35**2, 1.35*2.32*0.784], [0, 1.10**2]]},
        'D3': {'mean': [0.44, 0.19], 'cov': [[1.35**2, 1.35*0.96*0.830], [1.35*0.96*0.830, 0.96**2]]},
        'UR': {'mean': [0.54, 0.22], 'cov': [[1.35**2, 1.35*0.96*0.830], [1.35*0.96*0.830, 0.96**2]]}
    }
    agmr_df, hgmr_df  = calculate_mismatch(rl_df)
    relationships = determine_relationships(agmr_df, hgmr_df, relationship_parameters)
    return relationships

def queryKinship(chainName, datadir, sampleSearch):
    rl_df, af = querySampleRelatedness(chainName, datadir, sampleSearch)
    relationships = assess_kinship(rl_df)
    #print(relationships.to_json())
    return relationships.to_json()

# Metadata queries
def queryMetadata(chainName, datadir, search_values):
    '''
    Extract metadata for patients
    Inputs:
        search values - list of metadata types to query. use of this query alone is intended with general metadata searches e.g. variant calling or seq_machine
    Output:
        dictionary with {search_key:{sub_key: list}}, sub_key is the specific value associate with each general metadata search 
    '''
    #Metadata query
    queryCommand = 'multichain-cli {} -datadir={} liststreamitems mappingData_metadata false 9999999'.format(chainName, datadir)
    items = subprocess.check_output(queryCommand.split())
    matches = json.loads(items, parse_int= int)
    #Parse the search values specified
    publishToAuditstream(chainName, datadir, queryCommand)
    all_patient_ids = {}
    if search_values:
        search_values = search_values.split(',')
        for search_value in search_values:
            filtered_dicts = [d for d in matches if all(key in d['keys'] for key in [search_value])]
            return_value = [[x for x in d["keys"] if x!=search_value][0] for d in filtered_dicts]
            patient_ids = {key:d['data']['json'] for d, key in zip(filtered_dicts, return_value)}
            all_patient_ids[search_value] = patient_ids
    #print(all_patient_ids)
    return all_patient_ids


def publishToAuditstream(chainName, datadir, queryCommand):
    #load wallet
    walletCommand='multichain-cli {} -datadir={} listaddresses'.format(chainName, datadir)
    items = subprocess.check_output(walletCommand.split())
    matches = json.loads(items, parse_int= int)
    wallets = pd.DataFrame(matches)
    wallet = wallets['address'][wallets['ismine'] == True][0]
    
    #load time and parse query conducted
    time = str(datetime.utcnow())
    query = ' '.join(queryCommand.split(' ')[3:])
    
    ##publish to auditstream
    streamKeys = (wallet, time, query)
    publishCommand = ['multichain-cli', 
            str('{}'.format(chainName)), 
            str('-datadir={}'.format(datadir)),
            'publish',
            str('audit_log'), 
            str('["{}", "{}", "{}"]'.format(streamKeys[0], streamKeys[1], streamKeys[2])),
            '{'+'"json":{}'+'}']
    procPublish = subprocess.Popen(publishCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    procPublish.wait()
    return


# ## Run queries


def main():
    parser = argparse.ArgumentParser()
    action_choices = ['pca', 'kin', 'meta']
    parser.add_argument('--view', choices=action_choices)
    parser.add_argument("-cn", "--chainName", help = "the name of the chain to store data", default = "chain1")
    parser.add_argument("-ml", "--multichainLoc", help = "path to multichain commands", default = "")
    parser.add_argument("-dr", "--datadir", help = "path to store the chain")
    parser.add_argument("-ss", "--sampleSearch", required=(action_choices[0:2]in sys.argv), help = "samples to search", default = None)
    parser.add_argument("-ks", "--kSearch", required=(action_choices[0]in sys.argv), help = "k loadings to search", default = 20)
    parser.add_argument("-md", "--metadata", required=(action_choices[0] in sys.argv), help = "metadata to search or filter on", default = None) 

    args = parser.parse_args()
    start = time.time()
    try:
        subscribeToStream(args.chainName, args.multichainLoc, args.datadir)
        if args.view == action_choices[0]:
            #PCA
            querySamplePCA(args.chainName,args.datadir, args.sampleSearch, args.kSearch)
        elif args.view == action_choices[1]:
            #Kinship
            queryKinship(args.chainName, args.datadir, args.sampleSearch)
        elif args.view == action_choices[2]:
            #Metadata
            queryMetadata(args.chainName, args.datadir, args.metadata)
        
        end = time.time()
        e = int(end - start)
        print('\n\n Time elapsed:\n\n')
        print( '{:02d}:{:02d}:{:02d}'.format(e // 3600, (e % 3600 // 60), e % 60))
    
    except Exception as e:
        print(e)
        sys.stderr.write("\nERROR: Failed query. Please try again.\n")
        quit()
        

if __name__ == "__main__":
    main()

