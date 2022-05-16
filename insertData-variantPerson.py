#!/usr/bin/env python
# coding: utf-8

# In[1]:


'''
insertData-variantPerson.py
Inserts data from VCF files (variant data)
Usage: $ python insertData.py -cn=<chain name> -dr=<Chain path> -vf=<VCF path> -mf=<Person_id:VCF mapping file path>
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
import numpy as np
import json
import warnings
import random
import itertools
from io import BytesIO
warnings.simplefilter("ignore")


# In[2]:


#Given a chain name and the name of the new stream, subscribe to that stream
def subscribeToStream(chainName, streamName, multichainLoc, datadir):
    subscribeStreamCommand=multichainLoc+'multichain-cli {} -datadir={} subscribe {}'.format(chainName,datadir,streamName)
    #subscribe to the stream
    procSubscribe = subprocess.Popen(subscribeStreamCommand.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    procSubscribe.wait()
    return


# In[3]:


#subscribe to person-variant streams
def subscribeToStreams(chainName, multichainLoc, datadir):
    for i in range(1, 23):
        subscribeToStream(chainName, "person_chrom_{}".format(i), multichainLoc, datadir)
        
    subscribeToStream(chainName, "mappingData_variants", multichainLoc, datadir)
    return


# In[4]:


def loadFilePaths(dataPath, variantFiles):
    '''
    Load the VCF file paths that will be inserted
    Input:
        dataPath - path where the VCF files are stored
        variantFiles - files to be added (one per chromosome)
    '''
    #parse the files that are part of user input
    if variantFiles == 'all':
        files = [x for x in range(1,23)]
    else:
        files = str.split(files.replace(" ",""), ',')
    paths = []    
    for file in files:
        filePath = 'ALL.chr{}.phase3_shapeit2_mvncall_integrated_v3plus_nounphased.rsID.genotypes.GRCh38_dbSNP_no_SVs.vcf.gz'.format(file)
        paths.append(dataPath+filePath)
    random.shuffle(paths)
    return paths


# In[5]:


def mappingSamplePerson(mappingFile, variantFile):
    '''
    Mapping of person_ids to sample_ids (necessary as using dummy data so the ids are not aligned)
    Input:
        mappingFile - path where the person:vcf ID mapping is held
        variantFiles - files to be added (one per chromosome)
    '''
    #mapping of person_ids to sample_ids, needed as using fake sample genes
    with open('{}.txt'.format(mappingFile)) as f:
        mapping = f.read()    
    sample_person = json.loads(mapping)
    #extract samples from vcf file
    request = 'bcftools query -l {}| head -n 100'.format(variantFile)
    output = subprocess.check_output(request, shell = True)
    samples = output.decode().split()
    #get dict in the same order as the vcf file
    sample_mapping = {x:sample_person[x] for x in samples}
    return sample_mapping


# In[25]:


def publishMappingPerson(chainName, multichainLoc, datadir, sample_mapping):
    '''
    Publish the samples that were added during this insertion
    Input:
        sample_mapping: dictionary converting person_ids:vcf sample_ids
    '''
    streamName = 'mappingData_variants'
    streamKeys = 'samples'
    streamValues = '{'+'"json":{}'.format(list(sample_mapping.values())) + '}'
    publishCommand = [multichainLoc+'multichain-cli', 
        str('{}'.format(chainName)), 
        str('-datadir={}'.format(datadir)),
        'publish',
        str('{}'.format(streamName)), 
        str('{}'.format(streamKeys)),
        str('{}'.format(streamValues))]
    procPublish = subprocess.Popen(publishCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    procPublish.wait()
    return
    


# In[7]:


def extractPersonVariants(file, sample_id, sample_col):
    '''
    Extract all variants for patient and wrangle to format for blockchain data entry
    Input:
        sample_id: id of the sample being added
        sample_col: the column in the VCF file that sample data is stored
    '''
    ##BCFtools request (CHANGE FROM JUST TAKING THE HEAD)
    request = 'bcftools view -s {}  -e \'GT[{}]="RR"\' -H  {} | head -n 100'.format(sample_id, sample_col, file)
    output = subprocess.check_output(request, shell = True)
    ##extract the data from the output
    df = pd.read_csv(BytesIO(output), sep='\t', usecols = [0,1,3,4,9], names= ['chrom', 'pos', 'ref', 'alt', 'gt'], index_col = 'pos')
    df.index = df.index.map(str)
    #streamname
    chrom = df.chrom.iloc[0]
    #dict format to add
    values = df[['ref','alt', 'gt']].T.to_dict(orient = 'list')
    return chrom, values


# In[8]:


def publishToDataStream(chainName, multichainLoc, datadir, streamName, streamKeys, streamValues):
    '''
    Publish person entry to field
    Input:
        streamName: chromosome being added
        streamKeys: person_id
        streamValues: all the variants:genotype for that person 
    '''
    publishCommand = [multichainLoc+'multichain-cli', 
        str('{}'.format(chainName)), 
        str('-datadir={}'.format(datadir)),
        'publish',
        str('person_chrom_{}'.format(streamName)), 
        str('{}'.format(streamKeys)),
        str('{}'.format(streamValues))]
    procPublish = subprocess.Popen(publishCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    procPublish.wait()
    return


# In[9]:


def publishToDataStreams(chainName, multichainLoc, datadir, mappingFile, paths):
    '''
    loop through all samples and add to multichain
    Input:
        mappingFile - path where the person:vcf ID mapping is held
        paths - path for the VCF files being added
    '''
    for path in paths:
        #load mapping dictionary and file paths
        sample_mapping = mappingSamplePerson(mappingFile, path)
        publishMappingPerson(chainName, multichainLoc, datadir, sample_mapping)
        for i, sample_id in enumerate(sample_mapping.keys()):
            streamName, streamValues = extractPersonVariants(path, sample_id, i)
            streamKeys = sample_mapping[sample_id]
            streamValues ='{'+'"json":{}'.format(str(streamValues)) +'}'#create JSON data object
            streamValues = streamValues.replace("'",'"')
            publishToDataStream(chainName, multichainLoc, datadir, streamName, streamKeys, streamValues)
    return


# In[ ]:


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-dp", "--dataPath", type = str, help = "path of data to add to chain")
    parser.add_argument("-cn", "--chainName", help = "the name of the chain to store data", default = "chain1")
    parser.add_argument("-ml", "--multichainLoc", help = "path to multichain commands", default = "")
    parser.add_argument("-dr", "--datadir", help = "path to store the chain")
    parser.add_argument("-mf", "--mappingfile", help = "path to sample mapping file")
    parser.add_argument("-vf", "--variantfile", help = "variant files to add", default = "all")
    args = parser.parse_args()

    start = time.time()
    
    cpu = multiprocessing.cpu_count()
    print('CPUs available: {}'.format(cpu))
    
    try:
        subscribeToStreams(args.chainName, args.multichainLoc, args.datadir)
        print('Subscribed to streams') 
        
        paths = loadFilePaths(args.dataPath, args.variantfile)
        paths_split = np.array_split(paths, cpu)
        processes = []
        for i in range(cpu):
            paths_split_ins = paths_split[i]
            p_ins = multiprocessing.Process(target=publishToDataStreams, args = (args.chainName, args.multichainLoc,
                                                                                     args.datadir, args.mappingfile, paths_split_ins))
            processes.append(p_ins)
            p_ins.start()
            
            
            for process in processes:
                process.join()
            
            print('Inserted {}'.format(paths_split_ins))
            end = time.time()
            e = int(end - start)
            print('\n\n Time elapsed:\n\n')
            print( '{:02d}:{:02d}:{:02d}'.format(e // 3600, (e % 3600 // 60), e % 60))

        
        end = time.time()
        e = int(end - start)
        print('\n\n Time elapsed:\n\n')
        print( '{:02d}:{:02d}:{:02d}'.format(e // 3600, (e % 3600 // 60), e % 60))

        process = psutil.Process(os.getpid())
        print('\n\n Total memory in bytes:\n\n')
        print(process.memory_info().rss)
    
    except Exception as e:
        print(e)
        sys.stderr.write("\nERROR: Failed stream publishing. Please try again.\n")
        quit()
        

if __name__ == "__main__":
    main()

