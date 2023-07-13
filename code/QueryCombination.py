#!/usr/bin/env python
# coding: utf-8
 
# In[538]:


'''
QueryCombination.py
Queries data across the 3 data inputs - clinical, gtf and vcf
Usage: $ python QueryCombination.py -cn=<chain name> -dr=<Chain path> -ch=<Chromosome> -ps=<Variant position> -gn=<Gene> -ir=<MAF range> -ck<OMOP keys to define Cohort>
modified by AE 05/2022
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
import ast
import warnings
import multiprocessing
import numpy as np
from itertools import compress
from datetime import datetime
from pprint import pprint as pp
import json
import pdb
import itertools
from json.decoder import JSONDecodeError
warnings.simplefilter(action='ignore')


# In[540]:


#Given a chain name subscribe to audit log stream to ensure query is recorded
def subscribeToStream(chainName, multichainLoc, datadir):
    subscribeStreamCommand=multichainLoc+'multichain-cli {} -datadir={} subscribe audit_log'.format(chainName,datadir)
    #subscribe to the stream
    procSubscribe = subprocess.Popen(subscribeStreamCommand.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    procSubscribe.wait()
    return


# ## variant -> gene queries

# In[541]:


def extractVariantGenes(chainName, multichainLoc, datadir, variant, chrom):
    '''
    For a given variant position, extract the gene associated with it
    Inputs:
        variant - position of the variant of interest
        chrom - which chromosome the variant is on
    '''
    ##Multichain query
    queryCommand = 'multichain-cli {} -datadir={} liststreamkeyitems gene_variant_chrom_{} {}'.format(chainName, datadir,
                                                                                                     chrom, variant)
    items = subprocess.check_output(queryCommand.split())
    matches = json.loads(items, parse_int= int)
    publishToAuditstream(chainName, multichainLoc, datadir, queryCommand)
    ##Extract gene id from the query
    if matches != []:
        gene = matches[0]['keys'][1]
        return gene
    return None
    
    
def exractGeneData(chainName, multichainLoc, datadir, gene, variant, chrom):
    '''
    For the extracted gene, get its information e.g. function, type etc.
    As GTF has multiple genes under the same ID (depending on the source), filter for Genes that actually contain variant
    Inputs:
        gene - gene of interest (extracted using extractVariantGenes)
        variant - position of the variant of interest
        chrom - which chromosome the variant is on
    '''
    ##Multichain query
    queryCommand = 'multichain-cli {} -datadir={} liststreamkeyitems gene_chrom_{} {} false 9999999999'.format(chainName, datadir,
                                                                                                     chrom, gene)
    items = subprocess.check_output(queryCommand.split())
    matches = json.loads(items, parse_int= int)
    
    ##wrangle the data json object and the keys object to a DF
    df = pd.DataFrame(matches)[['keys', 'data']]
    keys = pd.DataFrame(df["keys"].to_list(), columns=['gene_id', 'name', 'feature'])

    def createrow(row):
        return pd.Series(row['data']['json'])
    
    gene_info = df.apply(createrow, axis =1)
    gene_df = keys.merge(gene_info, left_index = True, right_index=True)
    
    ##filter the genes that actually contain the variant (needed as multiple genes have same ID depending on sequencing source)
    gene_df_filtered = gene_df[(gene_df['start'] <= int(variant)) & (gene_df['end'] >= int(variant))]
    gene_df_filtered['variant'] = variant
    
    return gene_df_filtered


# In[543]:


def queryVariantGene(chainName, multichainLoc, datadir, variants, chrom):
    '''
    Full query that takes in variants of interest and extracts gene information
    Inputs:
        variant - position of the variant of interest
        chrom - which chromosome the variant is on
    '''
    if type(variants) is str:
        variants = variants.split(',')
    for variant in variants:
        gene = extractVariantGenes(chainName, multichainLoc, datadir, variant, chrom)
        if gene:
            gene_data = exractGeneData(chainName, multichainLoc, datadir, gene, variant, chrom)
            return gene_data
        else:
            print('No gene associated with variant position {}'.format(variant))


# ## gene -> variant query

# In[544]:


def extractGeneVariants(chainName, multichainLoc, datadir, gene, chrom):
    '''
    Given a gene of interest, extract all variant positions associated with it
    Inputs:
        gene - gene_id of gene of interest
        chrom - which chromosome the gene is on
    '''
    ##multichain command
    queryCommand = 'multichain-cli {} -datadir={} liststreamkeyitems gene_variant_chrom_{} {} false 99'.format(chainName, datadir,
                                                                                                     chrom, gene)
    items = subprocess.check_output(queryCommand.split())
    matches = json.loads(items, parse_int= int)
    publishToAuditstream(chainName, multichainLoc, datadir, queryCommand)
    ##parse returned json object to get the matches (will be multiple)
    if matches != []:
        variants = [match['keys'][0] for match in matches]
        return variants
    return None


# ## MAF -> variant -> gene queries

# In[545]:


def queryRangeParser(inputRange):
    '''
    Parse the MAF range that the user inputs to format necessary to extract data from MAF stream
    This requires finding the MAF buckets (0-0.1 etc) that the input range covers (variant MAFs are stored as data entries in these buckets)
    Inputs:
        inputRange - the MAF range user has inputted
    '''
    ##parse string input to numeric (*100 as want integers not floats)
    inputRange_split = list(map(float, inputRange.split('-')))
    inputRange_int = [split*100 for split in inputRange_split]
    ##identify which buckets the input range covers
    streamRanges_int = [(0,5), (5,10), (10, 15), (15,20), (20,30), (30,40), (40,50), (50,100) ]
    desiredRanges = [1 if ((inputRange_int[0] in range(r[0], r[1])) or (inputRange_int[1] in range(r[0],r[1]))) else 0 for r in streamRanges_int]
    #if range goes across a whole stream entirely
    desiredRanges2 = [1 if ((desiredRanges[i] ==1) or (sum(desiredRanges[i+1:]) >= 1 & sum(desiredRanges[:i]) >= 1)) else 0 for i in range(len(desiredRanges))]
    ##convert back to string ranges used to store data on multichain
    rangesStreams_str = ['0-0.05','0.05-0.1', '0.1-0.15', '0.15-0.2','0.2-0.3', '0.3-0.4','0.4-0.5', '0.5-1']
    desiredRanges_streams = list(compress(rangesStreams_str, desiredRanges2))
    return desiredRanges_streams, np.array(inputRange_int) / 100


# In[546]:


def MAFquery(chainName, multichainLoc, datadir, chrom, streamRange, numericRange):
    '''
    Given a MAF range user has inputted, extract the variants from the MAF stream that correspond to it
    Inputs:
        streamRange - the MAF stream entries that are included in the user range
        numericRange - the numeric values of the user inputted range
    '''
    ##multichain command to extract positions from MAF stream using streamRange
    queryCommand = 'multichain-cli {} -datadir={} liststreamkeyitems MAF_chrom_{} {}'.format(chainName, datadir, chrom, streamRange)
    items = subprocess.check_output(queryCommand.split())
    matches = json.loads(items, parse_int= int)
    MAF_variants = matches[0]['data']['json']
    try:
        MAF_variants = json.loads(MAF_variants.replace("(",'"(').replace(")",')"'))
    except:
        MAF_variants = MAF_variants.replace("[", '(').replace("]",")").replace('A', '"A"').replace('C', '"C"').replace('G', '"G"').replace('T', '"T"').replace('0|0', '"0|0"').replace('1|0', '"1|0"'
                                    ).replace('1|1', '"1|1"').replace('2|0', '"2|0"').replace('2|2', '"2|2"').replace('2|1', '"2|1"').replace('3|0', '"3|0"').replace('3|1', '"3|1"').replace(
                                        '3|2', '"3|2"').replace('3|3', '"3|3"')
        MAF_variants = eval(MAF_variants)
    ##create DF from returned values and filter only for those within actual numeric range
    MAF_df = pd.DataFrame.from_dict(MAF_variants, orient = 'index', columns = ['MAF'])
    MAF_df = MAF_df[(MAF_df['MAF'] >= numericRange[0]) & (MAF_df['MAF'] <= numericRange[1])]
    MAF_df.reset_index(inplace = True)
    ##function to parse the returned JSON object with variant information
    def parseVariantInfo(row):
        try:
            info = row[0].replace(' ','').replace('(','').replace(')','').split(',')
            pos, ref, alt, gt = info[0], info[1], info[2:-1], info[-1]
        except:
            info = row[0]
            pos, ref, alt, gt = info[0], info[1], info[2:-1], info[-1]
        return [pos, ref, alt, gt, row['MAF']]
    ##loop through the df and parse through every entry to create cleaned df
    if not MAF_df.empty:
        MAF_parsed = MAF_df.apply(parseVariantInfo, axis = 1)
        MAF_parsed_df = pd.DataFrame(np.vstack(MAF_parsed), columns = ['pos','ref','alt','gt','maf'])
        publishToAuditstream(chainName, multichainLoc, datadir, queryCommand)
        return MAF_parsed_df
    
    publishToAuditstream(chainName, multichainLoc, datadir, queryCommand)
    return


# In[547]:


def MAFqueries(chainName, multichainLoc, datadir, chrom, inputRange):
    '''
    Full MAF query that uses queryRangeParser and MAFquery to extract relevant positions
    Inputs:
        chrom - chromosome of interest
        inputRange - the MAF range user has inputted
    '''
    streamRanges, numericRanges = queryRangeParser(inputRange)
    MAFquery_df = pd.DataFrame(columns = ['pos','ref','alt','gt','maf'])
    for streamRange in streamRanges:
        MAFquery_df = MAFquery_df.append(MAFquery(chainName, multichainLoc, datadir, chrom, streamRange, numericRanges))
    return MAFquery_df


# In[548]:


def extractMAFPersonIDs(chainName, multichainLoc, datadir, chrom, inputRange):
    '''
    Builds on MAF queries to extract the PersonIDs of samples that contain a variant-position with a MAF within the desired range
    Inputs:
        chrom - chromosome of interest
        inputRange - the MAF range user has inputted
    '''
    ##extract the variants within MAF range
    MAF_df = MAFqueries(chainName, multichainLoc, datadir, chrom, inputRange)
    variant_dict = {}
    ##loop through each variant and search the variant stream for the person IDs that have that variant
    for _ , row in MAF_df.iterrows():
        pos = row['pos']
        gt = row['gt']
        ##parses gt if heterozygous (unlikely to be needed now)
        gt = str(gt).replace('[','(').replace(']',')') if len(gt) > 1  else gt[0]
        queryCommand = 'multichain-cli {} -datadir={} liststreamqueryitems chrom_{} {{"keys":["{}","{}"]}}'.format(
            chainName, datadir, chrom, pos, gt)
        items = subprocess.check_output(queryCommand.split())
        matches = json.loads(items, parse_int= int)
        variant_dict[(pos, gt)] = matches[0]['data']['json']

    return variant_dict


# In[549]:


def queryMAFVariantGene(chainName, multichainLoc, datadir, chrom, inputRange):
    '''
    Takes the variants with a specific MAF and and identifies genes that are associated. Also returns the personIDs with said variant
    Inputs:
        chrom - chromosome of interest
        inputRange - the MAF range user has inputted
    '''
    ##dictionary of variant, MAF and associated personIDs
    variant_dict = extractMAFPersonIDs(chainName, multichainLoc, datadir, chrom, inputRange)
    ##loop through each variant and query the gene stream for associated genes
    results = []
    for variant, person_ids in variant_dict.items():
        gene_df = queryVariantGene(chainName, multichainLoc, datadir, [variant[0]], chrom)
        try:
            gene_df['gt'] = variant[1]
            gene_df['person_ids'] = [person_ids for _ in range(len(gene_df))]
            results.append(gene_df)
            # Avoid printing for now
            # print(gene_df)
        except:
            pass
    return results



# ## disease -> gene -> variant

# In[550]:


def queryMappingStream(chainName, multichainLoc, datadir, keys):
    '''
    Find the streams where the OMOP cohort keys are inserted from the mapping stream
    Input:
        keys - OMOP keys to build cohort
    '''
    concepts = []
    for key in keys: 
        queryCommand=multichainLoc+'multichain-cli {} -datadir={} liststreamkeyitems mappingData_clinical {}'.format(chainName, datadir, key)
        items = subprocess.check_output(queryCommand.split())
        info = tuple(json.loads(items, parse_int= int)[0]['keys'])
        concepts.append(info)
    return list(set(concepts))


# In[551]:


def extractDataStream(chainName, multichainLoc, datadir, cohortKeys):
    '''
    Using the cohort keys, extract the stream and bucket that the data is stored under
    This is mostly parsing through the queryMappingStream function
    Inputs:
        cohortKeys - OMOP keys to build cohort
    '''
    matches = queryMappingStream(chainName, multichainLoc, datadir, cohortKeys)
    concept_domain, concept_stream, concept_bucket = matches[0][1], matches[0][2], ast.literal_eval(matches[0][3])
    return concept_domain, concept_stream, concept_bucket


# In[552]:


def extractPersonIDs(chainName, multichainLoc, datadir, cohortKeys):
    '''
    Given cohortKeys, get the personIDs who have that feature (i.e. diagnosis, medication)
    Inputs:
        cohortKeys - OMOP keys to build cohort
    '''
    ##parse inputted cohort keys
    cohortKeys = str.split(cohortKeys.replace(" ",""), ',')
    ##extract which stream and bucket they are stored under
    concept_domain, concept_stream, concept_bucket = extractDataStream(chainName, multichainLoc, datadir, cohortKeys)
    matches = []
    ##for every stream bucket, query for the concept
    for bucket in concept_bucket:
        ##potentially change to liststreamkeyitems
        queryCommand=multichainLoc+'multichain-cli {} -datadir={} liststreamkeyitems {}_id_{}_bucket_{} {}'.format(chainName, datadir,
                                                                                                               concept_domain, concept_stream,
                                                                                                               bucket+1, cohortKeys[0])
        items = subprocess.check_output(queryCommand.split())
        matches.extend(json.loads(items, parse_int= int))
    
    ids = []
    ##parse through the match to get the personIDs
    for match in matches:
        ids.append(int(match['keys'][1]))
    return(list(set(ids)))


# In[553]:


def calculateMAF(chainName, multichainLoc, datadir, chrom, variant, gt):
    '''
    Given a variant and genotype, calculate the MAF (on the fly)
    Takes the number of samples stored and the number of samples with that variant to calculate the MAD
    Inputs:
        chrom - chromosome of variant
        variant - position of the variant
        gt - geneotype e.g. 1|0, 0|0
    '''

    ##Search mapping stream for all the samples added to the chain
    queryCommand=multichainLoc+'multichain-cli {} -datadir={} liststreamkeyitems mappingData_variants samples'.format(chainName, datadir)

    items = subprocess.check_output(queryCommand.split())
    matches = json.loads(items, parse_int= int)
    ##count the number of samples (not there will be multiple matches from the query as each time a batch of samples is added a new entry is created)
    samples = []
    #BEGIN_NEW#
    for match in matches:
            ##track mapping files for all unique samples added
                try:
                    samples.extend(match['data']['json'])
                except:
                    items = get_json_payload_from_txid(match['data'].get('txid'), chainName, datadir)
                    match = json.loads(items, parse_int= int)
                    samples.extend(match['json'])
    #END_NEW#
    samples = len(set(samples))

    ##if not a homozygous gt then carry out search and count number of samples that match
    if gt != '0|0':
        ##Search mapping stream for all the samples added to the chain
        queryCommand='multichain-cli {} -datadir={} liststreamqueryitems chrom_{} {{"keys":["{}","{}"]}}'.format(chainName, datadir,
                                                                                                    chrom, variant, gt)

        items = subprocess.check_output(queryCommand.split())
        matches = json.loads(items, parse_int= int)
        ##count the number of samples (not there will be multiple matches from the query as each time a batch of samples is added a new entry is created)
        alleleMatch = []
        #BEGIN_NEW#
        for match in matches:
            ##track mapping files for all unique samples added
                try:
                    alleleMatch.extend(match['data']['json'])
                except:
                    items = get_json_payload_from_txid(match['data'].get('txid'), chainName, datadir)
                    match = json.loads(items, parse_int= int)
                    alleleMatch.extend(match['json'])
        #END_NEW#
        alleleMatch = len(set( alleleMatch ))
    ##if a homozygous gt then add up all the matches and take #full samples - result (this is because 0|0 is not stored on chain)
    else:
        queryCommand='multichain-cli {} -datadir={} liststreamkeyitems chrom_{} {}'.format(chainName, datadir,
                                                                                                    chrom, variant)
        items = subprocess.check_output(queryCommand.split())
        matches = json.loads(items, parse_int= int)

        alleleMatch = []
        #BEGIN_NEW#
        for match in matches:
            ##track mapping files for all unique samples added
                try:
                    alleleMatch.extend(match['data']['json'])
                except:
                    items = get_json_payload_from_txid(match['data'].get('txid'), chainName, datadir)
                    match = json.loads(items, parse_int= int)
                    alleleMatch.extend(match['json'])
        #END_NEW#
        alleleMatch = len(set( alleleMatch ))
        alleleMatch = samples - alleleMatch
    return round(alleleMatch / samples,3)


# In[554]:

def queryClinicalGeneVariant(chainName, multichainLoc, datadir, cohortKeys, gene, chrom):
    '''
    Given the gene of interest and clinical cohort definition, get the variants of those patients in that gene
    Inputs:
        cohortKeys - OMOP keys that define the cohort
        gene - gene_id of gene of interest
        chrom - chromosome of variant
    '''

    ##extract personIDs for cohort and variants associated with the gene of interest
    person_ids = extractPersonIDs(chainName, multichainLoc, datadir, cohortKeys)
    variants = extractGeneVariants(chainName, multichainLoc, datadir, gene, chrom)

    ##build dictionary that loops through the personIDs and variants and extracts the patients alleles i.e. what genotype do they hold in that position
    variants_dict = {}
    for person_id in person_ids:
        queryCommand=multichainLoc+'multichain-cli {} -datadir={} liststreamkeyitems person_chrom_{} {} false 999999'.format(chainName, datadir,
                                                                                                        chrom, person_id)

        items = subprocess.check_output(queryCommand.split())
        matches = json.loads(items, parse_int= int)
        #BEGIN_NEW#
        variants_person = []
        for match_ in matches:
            try:
                variants_person.append(match_['data']['json'])
            except:
                items = get_json_payload_from_txid(match_['data'].get('txid'), chainName, datadir)
                match_ = json.loads(items, parse_int= int)
                variants_person.append(match_['json'])
        variants_person = dict(itertools.chain.from_iterable(d.items() for d in variants_person))
        variants_person_filtered = {k:v for k,v in variants_person.items() if k in variants}
        for variant in variants_person_filtered:
            try:
                gt = variants_person_filtered[variant][-1]
            except:
                ##if patient has no variant stored here it must be 0|0 which is not stored on chain
                gt = '0|0'
            key = (variant, gt)
            if key in variants_dict:
                variants_dict[key]['person_id'].append(person_id)
            else:
                variants_dict[key] = {'person_id':[person_id]}
                ##for each variant, calculate the MAF and also add this to the dictionary
                MAF = calculateMAF(chainName, multichainLoc, datadir, chrom, variant, gt)
                variants_dict[key]['MAF'] = MAF

    
    publishToAuditstream(chainName, multichainLoc, datadir, queryCommand)
    
    return variants_dict

# In[556]:


def queryClinicalGeneVariantRange(chainName, multichainLoc, datadir, cohortKeys, gene, chrom, MAFRange):
    '''
    Full query that gives the alleles of variants associated with a gene of interest in a cohort
    Can filter the results for only variants within MAF range i.e. rare variants
    Inputs:
        cohortKeys - OMOP keys that define the cohort
        gene - gene_id of gene of interest
        chrom - chromosome of variant
        MAFRange - inputted MAF range to filter variants
    '''
    if MAFRange != "":
        _, numericRanges = queryRangeParser(MAFRange)
    variants_dict = queryClinicalGeneVariant(chainName, multichainLoc, datadir, cohortKeys, gene, chrom)
    variants_df = pd.DataFrame.from_dict(variants_dict, orient = 'index')
    if MAFRange != "":
        variants_df_filtered = variants_df[(variants_df['MAF'] >= numericRanges[0]) & (variants_df['MAF'] <= numericRanges[1])]
    else:
        variants_df_filtered = variants_df
    # Avoid printing for now
    #print(variants_df_filtered)
    return variants_df_filtered

# ## Variant to clinical queries
def extractVariantPersonIDs(chainName, datadir, chrom, variant):
    '''
    Get the IDs of patients with specific variant
    '''
    variant_dict = {}
    queryCommand = 'multichain-cli {} -datadir={} liststreamkeyitems chrom_{} {} false 9999999999999999'.format(chainName, datadir, chrom, variant)
    items = subprocess.check_output(queryCommand.split())
    matches = json.loads(items, parse_int= int)
    persons_gt = {}
    for match in matches:
        pos, gt = match['keys'][0], match['keys'][3]
        persons_gt[(pos, gt)] =  match['data']['json']
    return persons_gt

def queryDemographics(chainName, multichainLoc, datadir, person_ids):
    '''
    Get demographic data of returned patients
    '''
    ##check if search from running query or person specific query
    if isinstance(person_ids,str):
        personids = [int(id) for id in cohortKeys.split(',')]
    else:
        personids = person_ids
    persons_df = pd.DataFrame()
    
    for personid in personids:
        queryCommand=multichainLoc+'multichain-cli {} -datadir={} liststreamkeyitems person_demographics {}'.format(chainName, datadir, personid)
        items = subprocess.check_output(queryCommand.split())
        json_item = json.loads(items, parse_int= int)[0]['data']['json']
        person_df = pd.DataFrame.from_dict(json_item, orient = 'index').T
        persons_df = pd.concat([persons_df,person_df])
    persons_df.set_index('person_id', inplace = True)
    publishToAuditstream(chainName, multichainLoc, datadir, queryCommand)
    persons_json = persons_df.to_json(orient = 'index')
    print(persons_json)
    return persons_json


def extractPersonStreams(chainName, multichainLoc, datadir, person_ids):
    '''
    Get the perso stream that patients are stored in
    '''
    person_streams = {}
    for person_id in person_ids:
        queryCommand=multichainLoc+'multichain-cli {} -datadir={} liststreamkeyitems mappingData_person {}'.format(chainName, datadir, person_id)
        items = subprocess.check_output(queryCommand.split())
        matches = json.loads(items, parse_int= int) 
        person_streams[person_id] = matches[0]['data']['json']
    return person_streams

def queryPersonStreams(chainName, multichainLoc, datadir, person_ids, searchKeys):
    '''
    Get the clinical data for patients with variant
    '''
    person_streams = extractPersonStreams(chainName, multichainLoc, datadir, person_ids)
    data = {}
    for person_id in person_streams.keys():
        queryCommand = multichainLoc+'multichain-cli {} -datadir={}  liststreamkeyitems person_stream_{} {} false 150'.format(chainName, datadir,
                                                                                person_streams[person_id], person_id)
        items = subprocess.check_output(queryCommand.split())
        matches = json.loads(items, parse_int= int)
        for match in matches:
            key = match['keys'][0]
            if ('all' in searchKeys) | (key in searchKeys): 
                value = match['data']['json']
                d = pd.DataFrame.from_dict(value, orient = 'index').T
                if key in data:
                    data[key] = pd.concat([data[key],d])
                else:
                    data[key] = pd.DataFrame()
                    data[key] = pd.concat([data[key],d])
    
    for key in data:
        print(data[key])
        
    return data

def queryVariantClinical(chainName, multichainLoc, datadir, searchKeys, chrom, variant, gt):
    '''
    Overall query takes variant of interest (position and gt) and returns specific clinical data on patients with variant
    Note only 1 variant per search
    '''
    ##extract person_ids
    persons_gt = extractVariantPersonIDs(chainName, datadir, chrom, variant)
    searchKeys = searchKeys.split(',')
    ##filter for gt of interest
    for key in persons_gt:
        if key[1] == gt:
            person_ids = persons_gt[key]
    ##get the info of interest
    if searchKeys[0] == 'demographics':
        data = queryDemographics(chainName, multichainLoc, datadir, person_ids)
    else:
        data = queryPersonStreams(chainName, multichainLoc, datadir, person_ids, searchKeys)

    return data

# ## log queries

# In[ ]:
#BEGIN_NEW#
def get_json_payload_from_txid(txid, chainName, datadir):
    queryCommand = 'multichain-cli {} -datadir={} gettxoutdata {} 0'.format(chainName, datadir, txid)
    items = subprocess.check_output(queryCommand.split())
    return items
#END_NEW#

def publishToAuditstream(chainName, multichainLoc, datadir, queryCommand):
    #load wallet
    walletCommand=multichainLoc+'multichain-cli {} -datadir={} listaddresses'.format(chainName, datadir)
    items = subprocess.check_output(walletCommand.split())
    matches = json.loads(items, parse_int= int)
    wallets = pd.DataFrame(matches)
    wallet = wallets['address'][wallets['ismine'] == True][0]
    
    #load time and parse query conducted
    time = str(datetime.utcnow())
    query = ' '.join(queryCommand.split(' ')[3:])
    
    ##publish to auditstream
    streamKeys = (wallet, time, query)
    publishCommand = [multichainLoc+'multichain-cli', 
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

# In[ ]:


def main():
    parser = argparse.ArgumentParser()
    action_choices = ['variant', 'gene', 'maf', 'clinical']
    parser.add_argument('--query', choices=action_choices)
    parser.add_argument("-cn", "--chainName", help = "the name of the chain to store data", default = "chain1")
    parser.add_argument("-ml", "--multichainLoc", help = "path to multichain commands", default = "")
    parser.add_argument("-dr", "--datadir", help = "path to store the chain")
    parser.add_argument("-ch", "--chromosome", help = "chromosome to search")
    parser.add_argument("-ps", "--positions",required=(action_choices[0] in sys.argv), help = "variant positions to search")
    parser.add_argument("-gt", "--genotype",required=(action_choices[0] in sys.argv), help = "genotype to search")
    parser.add_argument("-sk", "--searchKeys", required=(action_choices[0] in sys.argv), default = 'all', help = "OMOP keys to extract clinical data")
    parser.add_argument("-gn", "--gene", required=(action_choices[1:4:2] in sys.argv), help = "genes to search")
    parser.add_argument("-ir", "--inputRange", required=(action_choices[2:4] in sys.argv), help = "MAF range to search", default = "")
    parser.add_argument("-ck", "--cohortKeys", required=(action_choices[3] in sys.argv), help = "OMOP keys to define cohort by")
    
    

    args = parser.parse_args()
    start = time.time()
    try:
        subscribeToStream(args.chainName, args.multichainLoc, args.datadir)
        if args.query == action_choices[0]:
            queryVariantClinical(args.chainName, args.multichainLoc, args.datadir, args.searchKeys, args.chromosome, args.positions, args.genotype)
        
        ##deprecated
        elif args.query == action_choices[1]:
            extractGeneVariants(args.chainName, args.multichainLoc, args.datadir, args.gene, args.chromosome)
        
        ##deprecated
        elif args.query== action_choices[2]:
            result = queryMAFVariantGene(args.chainName, args.multichainLoc, args.datadir, args.chromosome, args.inputRange)
        
        elif args.query == action_choices[3]:
            result, result_dict = queryClinicalGeneVariantRange(args.chainName, args.multichainLoc, args.datadir, args.cohortKeys, args.gene, args.chromosome, args.inputRange)
        
        filtered_results = []
        if args.query == action_choices[3] and result_dict:
            result = {str(k):v for k, v in result_dict.items()}
        elif type(result) is not list:
            result = json.loads(result.to_json(orient="records"))
        else:
            for item in result:
                if type(item) is pd.DataFrame and item.empty:
                    continue
                elif type(item) is pd.DataFrame:
                    filtered_results.append(json.loads(item.to_json(orient="records")))
                elif item:
                    filtered_results.append(item)
            result = filtered_results
        current_search = {"{}".format(" ".join(sys.argv[3:])): result}

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

