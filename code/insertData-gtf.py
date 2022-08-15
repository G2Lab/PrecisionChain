#!/usr/bin/env python
# coding: utf-8

# In[1]:


'''
insertData-gtf.py
Inserts data from gtf files (genetic data)
Usage: $ python insertData.py -cn=<chain name> -dr=<Chain path> -gp=<GTF path> -ch=<Chromosomes>
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
from gtfparse import read_gtf
from io import BytesIO, StringIO
import numpy as np
import random
warnings.simplefilter(action='ignore', category=FutureWarning)


# In[37]:


#Given a chain name and the name of the new stream, subscribe to that stream
def subscribeToStream(chainName, streamName, multichainLoc, datadir):
    subscribeStreamCommand=multichainLoc+'multichain-cli {} -datadir={} subscribe {}'.format(chainName,datadir,streamName)
    #subscribe to the stream
    procSubscribe = subprocess.Popen(subscribeStreamCommand.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    procSubscribe.wait()
    return


# In[38]:


#subscribe to person-variant streams
def subscribeToStreams(chainName, multichainLoc, datadir):
    for i in range(1, 23):
        subscribeToStream(chainName, "gene_chrom_{}".format(i), multichainLoc, datadir)
        subscribeToStream(chainName, "gene_variant_chrom_{}".format(i), multichainLoc, datadir)
    return


# In[39]:


def loadFilePaths(genePath, geneFiles, variantPath):
    '''
    Load the GTF file paths that will be inserted
    Input:
        genePath - path where the files are stored
        geneFiles - files to be added (one per chromosome)
    '''
    #parse the files that are part of user input
    if geneFiles == 'all':
        files = [x for x in range(1,23)]
    else:
        files = str.split(geneFiles.replace(" ",""), ',')
        
    #variant paths includes the VCF file of the patients that are being added - this is fixed based on the VCF file insertion
    variantPaths = [] 
    genePaths = []
    for file in files:
        variantPath_ = f'{variantPath}/chr{file}.vcf.gz'
        variantPaths.append(variantPath_)
        gene = '{}/gtf_chr{}.txt'.format(genePath, file)
        genePaths.append(gene)
    paths = list(zip(genePaths, variantPaths, files))
    #shuffle insertion for multiprocessing as some files are much larger than others
    random.shuffle(paths)
    return paths


# In[40]:


def publishToVariantStream(chainName, multichainLoc, datadir, streamName, streamKeys):
    '''
    Publish the data to the variant stream (keys are the variant position and gene it is associated with).
    Note that the data entry is empty
    Input:
        streamName - chromosome that the variant is in
        streamKeys - keys of the entry (position, gene)
    '''
    publishCommand = [multichainLoc+'multichain-cli', 
            str('{}'.format(chainName)), 
            str('-datadir={}'.format(datadir)),
            'publish',
            str('gene_variant_chrom_{}'.format(streamName)), 
            str('["{}", "{}", "{}"]'.format(streamKeys[0], streamKeys[1], streamKeys[2])),
            '{'+'"json":{}'+'}']
    procPublish = subprocess.Popen(publishCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    procPublish.wait()
    return


# In[41]:


def publishToGeneStream(chainName, multichainLoc, datadir, streamName, streamKeys, streamValues):
    '''
    Publish the data to the gene stream
    Input:
        streamName - chromosome that the variant is in
        streamKeys - keys of the entry (gene id, gene feature)
        streamValues - gene information (start, end, type, strand)
    '''
    
    publishCommand = [multichainLoc+'multichain-cli', 
            str('{}'.format(chainName)), 
            str('-datadir={}'.format(datadir)),
            'publish',
            str('gene_chrom_{}'.format(streamName)), 
            str('["{}", "{}", "{}"]'.format(streamKeys[0], streamKeys[1], streamKeys[2])),
            str('{}'.format(streamValues))]
    procPublish = subprocess.Popen(publishCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    procPublish.wait()
    return


# In[42]:


###CONSIDER MAKING IT ALL POSITIONS IN 1 GENE AS ENTRY TO SAVE ON INSERTION SPACE
def publishToStreams(gene, chainName, multichainLoc, datadir, chrom, variantFile):
    '''
    Extract relevant data for the gene including position, feature, type etc
    Input:
        chrom - chromosome
        variantfile - the VCF file associated with the same chromosome as the GTF file
    '''
    ##extract relevant gene info
    gene_id = gene['gene_id']
    gene_name = gene['gene_name']
    gene_feature = gene['feature']
    position = gene['start'], gene['end']
    ##create data entry information and publish to gene stream
    streamName = chrom
    streamKeys = gene_id, gene_name, gene_feature
    values = gene[['start', 'end', 'gene_type','strand']].to_dict()
    values = str(values).replace("\'", '"')
    streamValues = '{'+'"json":{}'.format(values) +'}'
    publishToGeneStream(chainName, multichainLoc, datadir, streamName, streamKeys, streamValues)
    
    ##extract variant info from VCF file related to that gene (i.e. all positions within start and end of gene)
    request = 'bcftools query -r {}:{}-{} -f \'%POS \' {}'.format(chrom, position[0], position[1], variantFile)
    output = subprocess.check_output(request, shell = True)  
    variants = output.decode('utf-8').split(' ')[:-1]
    ##for every variant position create entry with gene 
    for variant in variants:
        streamName = chrom
        streamKeys = [variant, gene_id, gene_name]
        publishToVariantStream(chainName, multichainLoc, datadir, streamName, streamKeys)
    return


# In[43]:


def publishGTF(chainName, multichainLoc, datadir, paths):
    '''
    For all GTF files - insert data
    Input:
        paths - the path to GTF files
    '''
    for path in paths:
        geneFile, variantFile, chrom = path
        #read in gtf file using data
        df = read_gtf(geneFile, 
                      usecols=['seqname','gene_id','feature','start','end', 'gene_type', 'gene_name','strand'])
        df.apply(publishToStreams, axis =1, args= (chainName, multichainLoc, datadir, chrom, variantFile))
    return


# In[ ]:


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-dp", "--dataPath", type = str, help = "path of data to add to chain")
    parser.add_argument("-cn", "--chainName", help = "the name of the chain to store data", default = "chain1")
    parser.add_argument("-ml", "--multichainLoc", help = "path to multichain commands", default = "")
    parser.add_argument("-dr", "--datadir", help = "path to store the chain")
    parser.add_argument("-gp", "--genepaths", help = "path to GTF files")
    parser.add_argument("-vp", "--variantpaths", help = "path to VCF files")
    parser.add_argument("-ch", "--chromosomes", help = "chromosomes to add", default = "all")
    args = parser.parse_args()

    start = time.time()
    
    cpu = multiprocessing.cpu_count()
    print('CPUs available: {}'.format(cpu))
    
    try:
        subscribeToStreams(args.chainName, args.multichainLoc, args.datadir)
        print('Subscribed to streams') 
        
        paths = loadFilePaths(args.genepaths, args.chromosomes, args.variantpaths)
        cpu = min(cpu, len(paths))
        paths_split = np.array_split(paths, cpu)
        processes = []
        for i in range(cpu):
            paths_split_ins = paths_split[i]
            p_ins = multiprocessing.Process(target=publishGTF, args = (args.chainName, args.multichainLoc,
                                                                                     args.datadir, paths_split_ins))
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
        sys.stderr.write("\nERROR: Failed stream insertion. Please try again.\n")
        quit()
        

if __name__ == "__main__":
    main()

