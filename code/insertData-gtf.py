#!/usr/bin/env python
# coding: utf-8

# In[1]:


'''
insertData-gtf.py
Inserts data from gtf files (genetic data)
Usage: $ python insertData.py -cn=<chain name> -dr=<Chain path> -gp=<GTF path> -ch=<Chromosomes>
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
        variantPath_ = f'{variantPath}/chr_{file}.vcf.gz'
        variantPaths.append(variantPath_)
        gene = '{}/gtf_variant_chr{}.txt'.format(genePath, file)
        genePaths.append(gene)
    paths = list(zip(genePaths, variantPaths, files))
    #shuffle insertion for multiprocessing as some files are much larger than others
    random.shuffle(paths)
    return paths


# In[40]:

def publishToVariantStream(chainName, multichainLoc, datadir, streamName, streamKeys, streamValues):
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
            '["{}"]'.format('", "'.join(streamKeys)),
            streamValues
            ]
    procPublish = subprocess.Popen(publishCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    procPublish.wait()
    return




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
    #publishToGeneStream(chainName, multichainLoc, datadir, streamName, streamKeys, streamValues)
    
    ##extract variant info from VCF file related to that gene (i.e. all positions within start and end of gene)
    request = 'bcftools query -r {}:{}-{} -f \'%POS \' {}'.format(chrom, position[0], position[1], variantFile)
    output = subprocess.check_output(request, shell = True)  
    variants = output.decode('utf-8').split(' ')[:-1]
    variants = [int(v) for v in variants]
    #BEGIN_NEW#
    ##get annotations
    clinvar =pd.read_csv(f'{annotation_path}/clinvar_annot.txt')
    vep = pd.read_csv(f'{annotation_path}/vep_annot.txt')
    columns = ['#Chrom', 'Pos', 'AnnoType', 'Consequence', 'ConsScore', 'ConsDetail', 'RawScore', 'PHRED']
    cadd =pd.read_csv(f'{annotation_path}/cadd_{chrom}.csv', skiprows = 1, sep = '\t', usecols = columns)

    for variant in variants:
        streamKeys = [str(variant), gene_id, gene_name]
        streamValues = {}
        if str(variant) in vep['position'][vep['chromosome'] == chrom].values:
            streamKeys.append('vep')
            vep_annot = list(vep[['#Uploaded_variation', 'Consequence']][(vep['position'] == str(variant)) & (vep['chromosome'] == 1)].values[0])
            streamValues['vep'] = vep_annot
        #clinvar
        if variant in clinvar['PositionVCF'][clinvar['Chromosome'] == chrom].values:
            streamKeys.append('clinvar')
            clinvar_annot = list(clinvar[['ClinicalSignificance', 'ClinSigSimple', 'PhenotypeList']][(clinvar['PositionVCF'] == variant) & (clinvar['Chromosome'] == str(chrom))].values[0])
            streamValues['clinvar'] = clinvar_annot
        if int(variant) in cadd['Pos'].values:
            streamKeys.append('cadd')
            cadd_annot = list(cadd[cadd['Pos'] == variant].iloc[:,2:].values[0])
            streamValues['cadd'] = cadd_annot

        streamValues = json.dumps({'json': streamValues})

        publishToVariantStream(chainName, multichainLoc, datadir, streamName, streamKeys, streamValues)
    #END_NEW#

def publishGTF(arguments):
    '''
    For all GTF files - insert data
    Input:
        paths - the path to GTF files
    '''
    chainName, multichainLoc, datadir, paths = arguments
    for path in paths:
        geneFile, variantFile, chrom = path
        #read in gtf file
        df = pd.read_csv(geneFile, usecols=['seqname','gene_id','feature','start','end', 'gene_type', 'gene_name','strand'])
        df = df.head(5)
        df.apply(publishToStreams, axis =1, args= (chainName, multichainLoc, datadir, chrom, variantFile))
    return


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-dp", "--dataPath", type = str, help = "path of data to add to chain")
    parser.add_argument("-cn", "--chainName", help = "the name of the chain to store data", default = "chain1")
    parser.add_argument("-ml", "--multichainLoc", help = "path to multichain commands", default = "")
    parser.add_argument("-dr", "--datadir", help = "path to store the chain")
    parser.add_argument("-gp", "--genepaths", help = "path to GTF files")
    parser.add_argument("-vp", "--variantpaths", help = "path to VCF files")
    parser.add_argument("-ap", "--annotationpath", help = "path to annotations") #NEW_LINE
    parser.add_argument("-ch", "--chromosomes", help = "chromosomes to add", default = "all")
    args = parser.parse_args()

    start = time.time()
    
    #cpu = multiprocessing.cpu_count() * 2
    cpu = 4
    print('CPUs available: {}'.format(cpu))
    
    global annotation_path
    annotation_path = args.annotationpath
    
    try:
        subscribeToStreams(args.chainName, args.multichainLoc, args.datadir)
        print('Subscribed to streams') 
        
        paths = loadFilePaths(args.genepaths, args.chromosomes, args.variantpaths)
        cpu = min(cpu, len(paths))
        paths_split = np.array_split(paths, cpu)
        arguments = []
        for paths_split_ins in paths_split:
            arguments.append((args.chainName, args.multichainLoc, args.datadir, paths_split_ins))
        pool = multiprocessing.Pool(cpu)
        pool.map(publishGTF, arguments)
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
        sys.stderr.write("\nERROR: Failed stream insertion. Please try again.\n")
        quit()
        

if __name__ == "__main__":
    main()