#!/usr/bin/env python
# coding: utf-8

# In[1]:


'''
insertData-anaylsis.py
Inserts analysis data including metdadata, pc loadings, and kinship
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
import time
import multiprocessing
import glob
import pandas as pd
import numpy as np
import json
import warnings
import random
import itertools
import psutil
from io import BytesIO, StringIO
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
    subscribeToStream(chainName, "mappingData_metadata", multichainLoc, datadir)
    subscribeToStream(chainName, "analysis", multichainLoc, datadir)
    return

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

def metadataPerson(metaFile, variantFile, numPeople):
    '''
    load metadata associated with the samples in the variant file
    Input:
        metaFile - path where the metadata is held
        variantFiles - files to be added (one per chromosome)
    '''
    
    #extract samples from vcf file
    request = 'bcftools query -l {} | head -n {}'.format(variantFile, numPeople)
    output = subprocess.check_output(request, shell = True)
    samples = output.decode().split()
    #metadata
    meta = pd.read_csv(f'{metaFile}')
    meta['id'] = meta['id'].astype(str)
    #get dict in the same order as the vcf file
    meta_seq = meta[meta['id'].isin(samples)]
    return meta_seq, samples

def publishMetadata(chainName, multichainLoc, datadir, meta):
    '''
    Publish the positions added into mapping stream
    Input:
        positions: list of positions added from vcf file
        chrom: chromosome positions come from
    '''
    ## metadata insert
    for _, row in meta.iterrows():
        row['id']
        streamName = 'mappingData_metadata'
        streamKeys = '{}'.format(row['id'])
        streamValues = '{'+'"json":{}'.format(json.dumps(list(row.values[1:]))) + '}'
        publishCommand = [multichainLoc+'multichain-cli', 
            str('{}'.format(chainName)), 
            str('-datadir={}'.format(datadir)),
            'publish',
            str('{}'.format(streamName)), 
            str('{}'.format(streamKeys)),
            str('{}'.format(streamValues))]
        procPublish = subprocess.Popen(publishCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        procPublish.wait()

    for col in meta.columns[1:]:
        for values in meta[col].unique():
            ids = list(meta.loc[:,'id'][meta[col] ==values].values)
            streamName = 'mappingData_metadata'
            streamKeys = '["{}", "{}"]'.format(col, values)
            streamValues = '{'+'"json":{}'.format(json.dumps(ids)) + '}'
            publishCommand = [multichainLoc+'multichain-cli', 
                str('{}'.format(chainName)), 
                str('-datadir={}'.format(datadir)),
                'publish',
                str('{}'.format(streamName)), 
                str('{}'.format(streamKeys)),
                str('{}'.format(streamValues))]
            procPublish = subprocess.Popen(publishCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            procPublish.wait()


### PCA
def loadData(pcaFile, samples):
    '''
    Load sample principal components
    Input:
        file paths for pca data
    '''
    #Sample PCs
    sample_pcs = pd.read_csv(f'{pcaFile}pca/samples_pcs.csv')
    sample_pcs = sample_pcs.iloc[:len(samples)]
    #Set sample IDs
    sample_pcs.index = samples
    return sample_pcs

def publishSamplePC(chainName, multichainLoc, datadir, pcaFile, samples):
    '''
    Insert samples principal components
    Input:
        file paths for pca data and sample ids
    '''
    sample_pcs = loadData(pcaFile, samples)
    for i, row in sample_pcs.iterrows():
        streamName = 'analysis'
        streamKeys = ['PCA', i]
        streamValues ='{"json":'+row.to_json()+'}' #create JSON and remove brackets

        publishCommand = [multichainLoc+'multichain-cli', 
            str('{}'.format(chainName)), 
            str('-datadir={}'.format(datadir)),
            'publish',
            str('{}'.format(streamName)), 
            str('["{}", "{}"]'.format(streamKeys[0], streamKeys[1])),
            str('{}'.format(streamValues))]
        procPublish = subprocess.Popen(publishCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        procPublish.wait()

### Kinship
def getRelatednessSNP(sampleids, relatedFile):
    '''
    Extract the related snps used by GRAF
    Input:
        sampleids - list of sampleids to insert
        relatedFile - location of the vcf file with snps
    '''
    #SNPS
    gt_final = {}
    for sampleid in sampleids:
        queryCommand = f"bcftools query -f '[%GT]\\n' -s {sampleid} {relatedFile}/merged.vcf.gz"
        #Load genotypes
        process = subprocess.Popen(queryCommand, shell=True, stdout=subprocess.PIPE)
        gt = []
        while True:
            line = process.stdout.readline().decode().strip()
            if not line:
                break
            gt.append(line)

        process.wait()

        gt_mapper = {'0|0':0, '0|1':1, '1|0':1, '1|1':2}
        #Update gt
        modified_gt = []
        for item in gt:
            x, y = item.split('|')
            x = min(int(x), 1)
            y = min(int(y), 1)
            modified_gt.append(f'{x}|{y}')
        #Create final list per sample
        gt_final[sampleid] = [gt_mapper[c] for c in modified_gt]

    #AF
    queryCommand = f"bcftools query -f '%AF\n'  {relatedFile}/merged.vcf.gz"
    #Load AF's
    process = subprocess.Popen(queryCommand, shell=True, stdout=subprocess.PIPE)
    af = []
    while True:
        line = process.stdout.readline().decode().strip()
        if not line:
            break
        af.append(line)
    process.wait()
    af = [a[:4] for a in af]
    return gt_final, af

def publishRelatednessSNP(chainName, datadir, sampleids, relatedFile):
    '''
    Publish the related snp gentoypes per sample onto the chain
    Input:
        sampleids - list of sampleids to insert
        relatedFile - location of the vcf file with snps
    '''
    gt_final, af = getRelatednessSNP(sampleids, relatedFile)
    #SNPs
    for sample in gt_final:
        streamName = 'analysis'
        streamKeys = ['Relatedness', sample]
        streamValues ='{"json":'+json.dumps(gt_final[sample])+'}' #create JSON and remove brackets
        publishCommand = ['multichain-cli', 
            str('{}'.format(chainName)), 
            str('-datadir={}'.format(datadir)),
            'publish',
            str('{}'.format(streamName)), 
            str('["{}","{}"]'.format(streamKeys[0], streamKeys[1])),
            str('{}'.format(streamValues))]
        procPublish = subprocess.Popen(publishCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        procPublish.wait()

    #AF
    streamName = 'analysis'
    streamKeys = ['Relatedness', 'AF']
    streamValues ='{"json":'+json.dumps(af)+'}' #create JSON and remove brackets
    publishCommand = ['multichain-cli', 
        str('{}'.format(chainName)), 
        str('-datadir={}'.format(datadir)),
        'publish',
        str('{}'.format(streamName)), 
        str('["{}","{}"]'.format(streamKeys[0], streamKeys[1])),
        str('{}'.format(streamValues))]
    procPublish = subprocess.Popen(publishCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-dp", "--dataPath", type = str, help = "path of data to add to chain")
    parser.add_argument("-cn", "--chainName", help = "the name of the chain to store data", default = "chain1")
    parser.add_argument("-ml", "--multichainLoc", help = "path to multichain commands", default = "")
    parser.add_argument("-dr", "--datadir", help = "path to store the chain")
    parser.add_argument("-mf", "--metafile", help = "path to sample metadata file")
    parser.add_argument("-vf", "--variantfile", help = "variant files to add", default = "all")
    parser.add_argument("-pc", "--pcafile", help = "pca loadings files for samples")
    parser.add_argument("-rf", "--relatedfile", help = "10k snps used to assess kinship")
    parser.add_argument("-np", "--numberPeople", help = "number of people to add", default = "100")
    args = parser.parse_args()

    start = time.time()
    
    cpu = multiprocessing.cpu_count()
    person = int(args.numberPeople)
    cpu = min(cpu, person)
    print('CPUs available: {}'.format(cpu))
    
    try:
        subscribeToStreams(args.chainName, args.multichainLoc, args.datadir)
        print('Subscribed to streams') 
        
        paths = loadFilePaths(args.dataPath, args.variantfile)

        #publish metadata
        meta, samples = metadataPerson(args.metafile, paths[0], int(args.numberPeople))
        publishMetadata(args.chainName, args.multichainLoc, args.datadir, meta)
        
        #publish pca
        publishSamplePC(args.chainName, args.multichainLoc, args.datadir, args.pcafile, samples)

        #publish kinship
        publishRelatednessSNP(args.chainName, args.datadir, samples, args.relatedfile)

                
        end = time.time()
        e = int(end - start)
        print('\n\n Time elapsed:\n\n')
        print( '{:02d}:{:02d}:{:02d}'.format(e // 3600, (e % 3600 // 60), e % 60))

        process = psutil.Process(os.getpid())
        print('\n\n Total memory in bytes:\n\n')
        print(process.memory_info().rss)
    
    except Exception as e:
        print(e)
        sys.stderr.write("\nERROR: Failed stream publishing analysis. Please try again.\n")
        quit()
        

if __name__ == "__main__":
    main()