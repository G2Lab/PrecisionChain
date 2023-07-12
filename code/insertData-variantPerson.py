'''
insertData-variantPerson.py
Inserts data from VCF files (variant data)
Usage: $ python insertData.py -cn=<chain name> -dr=<Chain path> -vf=<VCF path>
modified by AE 03/2023
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
        files = str.split(variantFiles.replace(" ",""), ',')
    paths = []    
    for file in files:
        filePath= f'{dataPath}/chr_{file}.vcf.gz'
        paths.append(filePath)
    random.shuffle(paths)
    return paths


# In[5]:


#BEGIN_NEW#
def metadataPerson(metaFile, sequencing, people):
    '''
    load metadata associated with the samples in the variant file
    Input:
        metaFile - path where the metadata is held
        sequencing - the sequencing type
    '''

    #metadata
    meta = pd.read_csv(f'{metaFile}')
    meta_seq = meta[meta['sequence'] == sequencing].iloc[:people]
    meta_seq['id'] = meta_seq['id'].astype(str)
    samples = meta_seq['id'].values
    return meta_seq, samples
#END_NEW#


# In[25]:


def extractPersonVariants(file, sample_id):
    '''
    Extract all variants for patient and wrangle to format for blockchain data entry
    Input:
        sample_id: id of the sample being added
    '''
    ##BCFtools request
    filter_stmt = '''awk '{if($10 != "0|0" && $10 != "0/0") { print } }' '''
    request = f'''bcftools view -s {sample_id}  -H  {file} | head -10000 | {filter_stmt} '''
    output = subprocess.check_output(request, shell = True)
    ##extract the data from the output
    df = pd.read_csv(BytesIO(output), sep='\t', usecols = [0,1,3,4,9], skiprows = 1, names= ['chrom', 'pos', 'ref', 'alt', 'gt'], index_col = 'pos')
    ##clean the ./. genotype
    df['gt'] = df['gt'].str[:3]
    df['gt'] =  df['gt'].str.replace('.','0')
    df.index = df.index.map(str)
    try:
        #streamname
        chrom = df['chrom'].iloc[0]
        #remove homoz genotypes
        df = df[df['gt'] != '0|0']
        ##reorder the genotypes
        def gt_parser(row):
            if int(row['gt'][0]) < int(row['gt'][2]):
                row['gt'] = f"{row['gt'][2]}|{row['gt'][0]}"
            return row
        df = df.apply(gt_parser, axis = 1)
        #dict format to add
        values = df[['ref','alt', 'gt']].T.to_dict(orient = 'list')
        return chrom, values
    except:
        chrom = file.split('/')[-1].split('.')[0]
        return chrom , False


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


def publishMappingPerson(chainName, multichainLoc, datadir, meta):
    '''
    Publish the samples that were added during this insertion
    Input:
        samples: List of all sample ID in file
    '''
    samples = list(meta['id'].values)
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
    
    ## person mapping
    for i, row in meta.iterrows():
        streamName = 'mappingData_variants'
        streamKeys = f"{row['id']}, {row['company']}, {row['seq_machine']}, {row['seq_protocol']}, {str(row['coverage'])}, {row['alignment_protocol']}, {row['variant_calling']}"
        streamValues = '{'+'"json":{}'+ '}'
        publishCommand = [multichainLoc+'multichain-cli', 
            str('{}'.format(chainName)), 
            str('-datadir={}'.format(datadir)),
            'publish',
            str('{}'.format(streamName)), 
            str('{}'.format(streamKeys)),
            str('{}'.format(streamValues))]
        procPublish = subprocess.Popen(publishCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        procPublish.wait()

    ## by sequence
    for seq in meta['sequence'].unique():
        s = list(meta['id'][meta['sequence'] == seq])
        streamName = 'mappingData_variants'
        streamKeys = f'{seq}'
        streamValues = '{'+'"json":{}'.format(json.dumps(s)) + '}'
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
    chainName, multichainLoc, datadir, metaFile, paths, people, sequencing = fields #NEW_LINE#
    '''
    loop through all samples and add to multichain
    Input:
        mappingFile - path where the person:vcf ID mapping is held
        paths - path for the VCF files being added
    '''
    for variantFile in paths:
        #load mapping dictionary and file paths
        meta, samples = metadataPerson(metaFile, sequencing, people) #NEW_LINE#
        publishMappingPerson(chainName, multichainLoc, datadir, meta)
        for sample_id in samples:
            streamName, streamValues = extractPersonVariants(variantFile, sample_id)
            ## only submit if have non ./. alleles
            if isinstance(streamValues, dict):
                ##chunk the dictionary as too large for one entry
                split_variants = chunkDictionary(streamValues, SIZE=200)
                for v in split_variants.values():
                    streamKeys = sample_id
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
    parser.add_argument("-mf", "--metafile", help = "path to sample metadata file")
    parser.add_argument("-vf", "--variantfile", help = "variant files to add", default = "all")
    parser.add_argument("-np", "--numberPeople", help = "number of people to add", default = "100")
    parser.add_argument("-sq", "--sequencing", help = "sequencing type") #NEWLINE
    
    args = parser.parse_args()

    start = time.time()
    
    cpu = multiprocessing.cpu_count()
    people = int(args.numberPeople)
    cpu = min(cpu, people)
    print('CPUs available: {}'.format(cpu))
    
    try:
        subscribeToStreams(args.chainName, args.multichainLoc, args.datadir)
        print('Subscribed to streams') 
        
        paths = loadFilePaths(args.dataPath, args.variantfile)
        paths_split = np.array_split(paths, cpu)
        arguments = []
        for paths_split_ins in paths_split:
            arguments.append((args.chainName, args.multichainLoc, args.datadir, args.metafile, paths_split_ins, people, args.sequencing))
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
        sys.stderr.write("\nERROR: Failed stream publishing - variantPerson. Please try again.\n")
        quit()




if __name__ == "__main__":
    main()

