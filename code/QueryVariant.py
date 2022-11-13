#!/usr/bin/env python
# coding: utf-8

# In[3]:


'''
QueryVariant.py
Queries VCF data either by searching for specific variants, sepcific samples or MAF ranges
Usage: $ python QueryCombination.py -cn=<chain name> -dr=<Chain path> -ch=<Chromosome> -ps=<Variant position> -gn=<Gene> -pi=<person ids> -ir=<MAF range> 
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
import warnings
import multiprocessing
import numpy as np
from itertools import compress
from datetime import datetime
warnings.simplefilter(action='ignore')


# In[ ]:


#Given a chain name subscribe to audit log stream to ensure query is recorded
def subscribeToStream(chainName, multichainLoc, datadir):
    subscribeStreamCommand=multichainLoc+'multichain-cli {} -datadir={} subscribe audit_log'.format(chainName,datadir)
    #subscribe to the stream
    procSubscribe = subprocess.Popen(subscribeStreamCommand.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    procSubscribe.wait()
    return


# ## variant queries

# In[1]:


def extractVariantsGenotypes(position, genotype):
    '''
    Parse the user inputted positions and genotypes of interest
    Input:
        position - the user-inputted variant positions
        genotype - the user-inputted genotype positions
    '''
    position = position.replace(' ','').split(',')
    genotype = genotype.replace(' ','').split(',')
    ##if heterozygous gt queried, parse gt to right format stored in blockchain (i.e. 1|0 )
    for gt in genotype:
        if gt[0] < gt[2]:
            genotype.remove(gt)
            genotype.append('{}|{}'.format(gt[2], gt[0]))
    genotype = list(set(genotype))
    return position, genotype
    


# In[4]:


def homozgyousPersons(chainName, multichainLoc, datadir, chrom, variant):
    '''
    Extract the person ids for those who have homozygous ref (0|0) alleles.
    Necessary as 0|0 alleles are not stored on chain
    Input:
        chrom - chromosome the variant is in
        variant - dictionary with person_ids for samples with non-reference homozygous alles
    '''
    ##command to extract all the samples added
    queryCommand = 'multichain-cli {} -datadir={} liststreamitems mappingData_variants'.format(chainName, datadir)
    items = subprocess.check_output(queryCommand.split())
    matches = json.loads(items, parse_int= int)
    all_persons = matches[0]['data']['json']
    ##extract the non-reference personIDs
    non_ref = []
    for persons in variant.values():
        non_ref.extend(persons)

    homozygous = list(set(all_persons) - set(non_ref))
    return homozygous


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

def extractGeneData(chainName, multichainLoc, datadir, gene, variant, chrom):
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

def queryVariantGene(chainName, multichainLoc, datadir, variants, chrom):
    '''
    Full query that takes in variants of interest and extracts gene information
    Inputs:
        variant - position of the variant of interest
        chrom - which chromosome the variant is on
    '''
    gene_data = pd.DataFrame()
    for variant in variants:
        gene = extractVariantGenes(chainName, multichainLoc, datadir, variant, chrom)
        if gene:
            gene_info = extractGeneData(chainName, multichainLoc, datadir, gene, variant, chrom)
            gene_info = gene_info.iloc[0:1]
            
        else:
            gene_info = pd.DataFrame(columns = ['gene_id', 'name', 'feature', 'start', 'end', 'gene_type', 'strand','variant'], 
                                data = np.array([('None', 'None', 'NA', 'NA', 'NA', 'NA','NA', variant)]))
        gene_data = pd.concat([gene_data, gene_info])
    return gene_data


# In[5]:

def queryVariant(chainName, multichainLoc, datadir, chrom, variant, genotype):
    '''
    For given genotype and variant, extract the person ids who have that variant
    Input:
        chrom - chromosome the variant is in
        variant - position of the variant of interest
        genotype - allele of interest i.e. 0|0 1|1 
    '''
    variant_dict = {}
    queryCommand = 'multichain-cli {} -datadir={} liststreamkeyitems chrom_{} {} false 9999999999999999'.format(chainName, datadir, chrom, variant)
    items = subprocess.check_output(queryCommand.split())
    matches = json.loads(items, parse_int= int)
    ##parse the matches normally if not homo-ref, if homo-ref then use specific function
    for match in matches:
        gt = match['keys'][3]
        if gt in genotype:
            variant_dict[gt] = match['data']['json']
    if '0|0' in genotype:
        variant_dict['0|0'] = homozgyousPersons(chainName, multichainLoc, datadir, chrom, variant_dict)
    publishToAuditstream(chainName, multichainLoc, datadir, queryCommand)
    return variant_dict


# In[6]:


def queryVariants(chainName, multichainLoc, datadir, chrom, variants, genotype):
    '''
    Loop through all variants inputted by user and extract person_ids
    Input:
        chrom - chromosome the variant is in
        variants - positions of the variant of interest
        genotype - allele of interest i.e. 0|0 1|1 
    '''
    variants, genotype = extractVariantsGenotypes(variants, genotype)
    variants_dict = {}
    for variant in variants:
        variants_dict[variant] = queryVariant(chainName, multichainLoc, datadir, chrom, variant, genotype)
    #get gene info
    gene_df = queryVariantGene(chainName, multichainLoc, datadir, variants, chrom)
    ##merge gene info and person variant info
    variants_df = pd.DataFrame.from_dict(variants_dict, orient = 'index')
    variants_df = variants_df.merge(gene_df, left_index = True, right_on = 'variant')
    ##make into json (similar format)
    variants_df.set_index('variant', inplace = True)
    variants_json = variants_df.to_json(orient = 'index')

    print(variants_json)
    return variants_json
        


# ## person queries

# In[7]:


def extractPersonIDsChromsPos(person_ids, chroms, pos):
    '''
    Parse through user-inputted person ids, chromosomes and positions
    Input:
        person_ids - the ids of samples to search
        chroms - chromosomes the variant is in
        pos - positions to search
    '''
    chroms = chroms.split(',')
    person_ids = person_ids.split(',')
    if pos != 'all':
        pos = pos.split(',')
    return chroms, person_ids, pos


# In[18]:


def queryPersonChrom(chainName, multichainLoc, datadir, chrom, person_id):
    '''
    For a given personID, extract the variants the person has
    Input:
        person_id - the id of sample to search
        chrom - chromosome the variant is in
    '''
    queryCommand = 'multichain-cli {} -datadir={} liststreamkeyitems person_chrom_{} {} false 9999999999999999'.format(chainName, datadir, chrom, person_id)
    items = subprocess.check_output(queryCommand.split())
    matches = json.loads(items, parse_int= int)
    person_variants_df = pd.DataFrame()
    for i, match in enumerate(matches):
        try:
            person_variants = match['data']['json']
            person_variants_df_ = pd.DataFrame.from_dict(person_variants, orient = 'index', columns = [
            'ref_allele_{}'.format(person_id),
            'alt_allele_{}'.format(person_id),
            'gt_{}'.format(person_id)])
            person_variants_df = pd.concat([person_variants_df, person_variants_df_], axis = 0)

        except:
            pass
    person_variants_df.drop_duplicates(inplace = True)
    
    publishToAuditstream(chainName, multichainLoc, datadir, queryCommand)
    
    return person_variants_df


# In[9]:


def extractAllPosition(chainName, multichainLoc, datadir, chrom):
    '''
    If user wants all variant data, extract every position stored on chain
    Input:
        chrom - chromosome the variants are in
    '''
    queryCommand = 'multichain-cli {} -datadir={} liststreamkeyitems mappingData_variants chrom_{}'.format(chainName, datadir, chrom)
    items = subprocess.check_output(queryCommand.split())
    matches = json.loads(items, parse_int= int)
    return matches[0]['data']['json']


# In[61]:


def queryPersonsChrom(chainName, multichainLoc, datadir, chrom, person_ids, pos):
    '''
    create dataframe of all variants (of interest) for inputted person_ids
    Input:
        person_ids - the ids of samples to search
        chrom - chromosome the variant is in
        pos - the position of the variants
    '''
    ##create empty dataframe
    person_df = pd.DataFrame(columns = ['ref_allele','alt_allele'])
    ##for every sample, extract their relevant positions and merge to the dataframe
    for person_id in person_ids:
        query_df = queryPersonChrom(chainName, multichainLoc, datadir, chrom, person_id)
        query_df.index = pd.to_numeric(query_df.index)
        person_df = person_df.merge(query_df, right_index = True, left_index = True, how = 'outer')
        person_df['ref_allele'].update(person_df.pop('ref_allele_{}'.format(person_id)))
        person_df['alt_allele'].update(person_df.pop('alt_allele_{}'.format(person_id)))
        ##where variant data for person doesnt exist, must be 0|0 genotype - this only fills in the 0|0 if at least
        ##one sample has non-ref allele
        person_df.iloc[:,1:] = person_df.iloc[:,1:].fillna('0|0')
        ##this extract all positions and identifies which ones are missing from dataframe and so must be 0|0 for all patients
        allPositions = extractAllPosition(chainName, multichainLoc, datadir, chrom)
        homo_pos = list(set(allPositions) - set(person_df.index))
        homo = pd.DataFrame(index = homo_pos, columns = person_df.columns)
        homo.iloc[:,2:] = homo.iloc[:,2:].fillna('0|0')
        person_full_df  = pd.concat([person_df, homo])
    ##filter out positions not of interest
    if pos != 'all':
        person_full_df = person_df.loc[pos]
    return person_full_df


# In[63]:


def queryPersonsChroms(chainName, multichainLoc, datadir, chroms, person_ids, pos):
    '''
    Full person query
    Input:
        person_ids - the ids of samples to search
        chroms - chromosomes the variants are in
        pos - the position of the variants
    '''
    chroms, person_ids, pos = extractPersonIDsChromsPos(person_ids, chroms, pos)
    queryReturn = {}
    for chrom in chroms:
        queryReturn[int(chrom)] = queryPersonsChrom(chainName, multichainLoc, datadir, chrom, person_ids, pos)
    print(queryReturn)
    return queryReturn


# ## MAF queries

# In[67]:


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

# ## variant-gene queries

def extractGeneVariants(chainName, multichainLoc, datadir, gene, chrom):
    '''
    Given a gene of interest, extract all variant positions associated with it
    Inputs:
        gene - gene_id of gene of interest
        chrom - which chromosome the gene is on
    '''
    ##multichain command
    queryCommand = 'multichain-cli {} -datadir={} liststreamkeyitems gene_variant_chrom_{} {} false 999'.format(chainName, datadir,
                                                                                                     chrom, gene)
    items = subprocess.check_output(queryCommand.split())
    matches = json.loads(items, parse_int= int)
    publishToAuditstream(chainName, multichainLoc, datadir, queryCommand)
    ##parse returned json object to get the matches (will be multiple)
    if matches != []:
        variants = [match['keys'][0] for match in matches]
        return variants
    return None


# ## log queries

# In[ ]:


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
    action_choices = ['variant', 'person', 'gene', 'maf']
    parser.add_argument('--view', choices=action_choices)
    parser.add_argument("-cn", "--chainName", help = "the name of the chain to store data", default = "chain1")
    parser.add_argument("-ml", "--multichainLoc", help = "path to multichain commands", default = "")
    parser.add_argument("-dr", "--datadir", help = "path to store the chain")
    parser.add_argument("-ch", "--chromosomes", help = "chromosome to search")
    parser.add_argument("-ps", "--positions",required=(action_choices[0:2] in sys.argv), help = "positions to search", default = "all")
    parser.add_argument("-gt", "--genotypes", required=(action_choices[0] in sys.argv), help = "genotypes to search")
    parser.add_argument("-pi", "--person_ids", required=(action_choices[1] in sys.argv), help = "person_ids to search")
    parser.add_argument("-gn", "--gene", required=(action_choices[2] in sys.argv), help = "genes to search")
    parser.add_argument("-ir", "--inputRange", required=(action_choices[3] in sys.argv), help = "MAF range to search")

    args = parser.parse_args()
    start = time.time()
    try:
        subscribeToStream(args.chainName, args.multichainLoc, args.datadir)
        if args.view == action_choices[0]:
            queryVariants(args.chainName, args.multichainLoc, args.datadir, args.chromosomes, args.positions, args.genotypes)
        
        elif args.view == action_choices[1]:
            queryPersonsChroms(args.chainName, args.multichainLoc, args.datadir, args.chromosomes, args.person_ids, args.positions)
        
        elif args.view == action_choices[2]:
             extractGeneVariants(args.chainName, args.multichainLoc, args.datadir, args.gene, args.chromosomes)

        elif args.view == action_choices[3]:
            MAFqueries(args.chainName, args.datadir, args.chromosomes, args.inputRange)
        
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

