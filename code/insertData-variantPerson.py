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
from itertools import islice
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
        filePath= f'{dataPath}/chr_{file}.vcf.gz'
        paths.append(filePath)
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


def extractPersonVariants(file, sample_id):
    '''
    Extract all variants for patient and wrangle to format for blockchain data entry
    Input:
        sample_id: id of the sample being added
    '''
    ##BCFtools request (CHANGE FROM JUST TAKING THE HEAD)
    filter_stmt = '''awk '{if($10 != "0|0") { print } }' '''
    request = f'''bcftools view -s {sample_id}  -H  {file}| {filter_stmt} '''
    output = subprocess.check_output(request, shell = True)
    ##extract the data from the output
    df = pd.read_csv(BytesIO(output), sep='\t', usecols = [0,1,3,4,9], names= ['chrom', 'pos', 'ref', 'alt', 'gt'], index_col = 'pos')
    ##clean the ./. genotype
    df['gt'] = df['gt'].str[:3]
    df['gt'] =  df['gt'].str.replace('.','0')
    df.index = df.index.map(str)
    try:
        #streamname
        chrom = df['chrom'].iloc[0]
        #remove homoz genotypes
        df = df[df['gt'] != '0|0']
        #dict format to add
        values = df[['ref','alt', 'gt']].T.to_dict(orient = 'list')
        return chrom, values
    except:
        chrom = file.split('/')[-1].split('.')[0]
        chrom, False


# In[26]:


def chunkDictionary(values_dict, SIZE=2000):
    '''
    chunk the person variant dictionary as its too large for a single entry
    Input:
        values_dict: person variant dictionary from extractPersonVariants
    '''
    def chunks(values_dict, SIZE=2000):
        it = iter(values_dict)
        for i in range(0, len(values_dict), SIZE):
            yield {k:values_dict[k] for k in islice(it, SIZE)}
    
    split_variants = {}
    for i, chunk in enumerate(chunks(values_dict, SIZE = 2000)):
        split_variants[i] = chunk
    return split_variants

# In[27]:


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


# In[28]:


def publishMappingPerson(chainName, multichainLoc, datadir, samples):
    '''
    Publish the samples that were added during this insertion
    Input:
        samples: List of all sample ID in file
    '''
    streamName = 'mappingData_variants'
    streamKeys = 'samples'
    streamValues = '{'+'"json":{}'.format(json.dumps(samples)) + '}'
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


# In[29]:


def publishToDataStreams(fields):
    chainName, multichainLoc, datadir, mappingFile, paths = fields
    '''
    loop through all samples and add to multichain
    Input:
        mappingFile - path where the person:vcf ID mapping is held
        paths - path for the VCF files being added
    '''
    for variantFile in paths:
        #load mapping dictionary and file paths
        sample_person = mappingSamplePerson(mappingFile, variantFile)
        publishMappingPerson(chainName, multichainLoc, datadir, list(sample_person.values()))
        for sample_id in sample_person:
            streamName, streamValues = extractPersonVariants(variantFile, sample_id)
            ## only submit if have non ./. alleles
            if isinstance(streamValues, dict):
                ##chunk the dictionary as too large for one entry
                split_variants = chunkDictionary(streamValues, SIZE=200)
                for v in split_variants.values():
                    streamKeys = sample_person[sample_id]
                    streamValues ='{'+'"json":{}'.format(json.dumps(v)) +'}'#create JSON data object
                    publishToDataStream(chainName, multichainLoc, datadir, streamName, streamKeys, streamValues)
            else:
                pass
                    
        print('Inserted {}'.format(variantFile))

    return


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
        arguments = []
        for paths_split_ins in paths_split:
            arguments.append((args.chainName, args.multichainLoc, args.datadir, args.mappingfile, paths_split_ins))
        pool = multiprocessing.Pool(cpu)
        pool.map(publishToDataStreams, arguments)
        pool.close()
        pool.join()


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

