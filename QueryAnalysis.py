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
from itertools import compress
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
    sampleSearch = sampleSearch.split(',')
    kSearch = int(kSearch)
    if sampleSearch != ['none']:
        pc_df = pc_df.loc[sampleSearch]
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
    sampleSearch = sampleSearch.split(',')
    if sampleSearch != ['none']:
        rl_df = rl_df.loc[sampleSearch]
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
    n = df.shape[0] 
    agmr_matrix = np.zeros((n, n)) 
    hgmr_matrix = np.zeros((n, n)) 

    for i in range(n):
        for j in range(i+1, n): 
            x = df.iloc[i]  
            y = df.iloc[j]  
            
            # AGMR: percentage of SNPs on which the two genotypes are not identical
            agmr = np.mean(x != y)
            agmr_matrix[i, j] = agmr
            agmr_matrix[j, i] = agmr  
            
            # HGMR: genotype mismatch rate when only the SNPs with homozygous calls for both samples are considered
            homozygous_indices = (x % 2 == 0) & (y % 2 == 0)  
            hgmr = np.mean(x[homozygous_indices] != y[homozygous_indices]) if np.sum(homozygous_indices) > 0 else 0
            hgmr_matrix[i, j] = hgmr
            hgmr_matrix[j, i] = hgmr 

    
    agmr_df = pd.DataFrame(agmr_matrix, index=df.index, columns=df.index)
    hgmr_df = pd.DataFrame(hgmr_matrix, index=df.index, columns=df.index)

    return agmr_df, hgmr_df

def determine_relationships(agmr_df, hgmr_df, relationship_parameters):
    relationships = pd.DataFrame(index=agmr_df.index, columns=agmr_df.columns)
    
    for i in agmr_df.index:
        for j in agmr_df.columns:
            agmr = agmr_df.loc[i, j]
            hgmr = hgmr_df.loc[i, j]

            max_likelihood = -np.inf
            max_rel = None
            for rel, params in relationship_parameters.items():
                likelihood = multivariate_normal(mean=params['mean'], cov=params['cov']).pdf([agmr, hgmr])
                if likelihood > max_likelihood:
                    max_likelihood = likelihood
                    max_rel = rel
            relationships.loc[i, j] = max_rel
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
    print(relationships.to_json())
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
    queryCommand = 'multichain-cli {} -datadir={} liststreamitems mappingData_metadata'.format(chainName, datadir)
    items = subprocess.check_output(queryCommand.split())
    matches = json.loads(items, parse_int= int)
    #Parse the search values specified
    publishToAuditstream(chainName, datadir, queryCommand)
    all_patient_ids = {}
    search_values = search_values.split(',')
    for search_value in search_values:
        filtered_dicts = [d for d in matches if all(key in d['keys'] for key in [search_value])]
        return_value = [[x for x in d["keys"] if x!=search_value][0] for d in filtered_dicts]
        patient_ids = {key:d['data']['json'] for d, key in zip(filtered_dicts, return_value)}
        all_patient_ids[search_value] = patient_ids
    print(all_patient_ids)
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
    parser.add_argument("-ss", "--sampleSearch", required=(action_choices[0:2]in sys.argv), help = "samples to search")
    parser.add_argument("-ks", "--kSearch", required=(action_choices[0]in sys.argv), help = "k loadings to search", default = 20)
    parser.add_argument("-md", "--metadata", required=(action_choices[0] in sys.argv), help = "metadata to search or filter on", default="none") 

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

