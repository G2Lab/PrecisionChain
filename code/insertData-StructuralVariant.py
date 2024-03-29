#!/usr/bin/env python
# coding: utf-8

# In[1]:


'''
insertData-variant.py
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
from io import BytesIO, StringIO
warnings.simplefilter("ignore")


# # Variant view

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
        subscribeToStream(chainName, "structural_chrom_{}".format(i), multichainLoc, datadir)
        subscribeToStream(chainName, "MAF_structural_chrom_{}".format(i), multichainLoc, datadir)
    subscribeToStream(chainName, "mappingData_variants_structural", multichainLoc, datadir)
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
        filePath = 'ALL.chr{}.phase3_shapeit2_mvncall_integrated_v3plus_nounphased.rsID.genotypes.GRCh38_dbSNP_no_SVs.vcf.gz'.format(file)
        paths.append(dataPath+filePath)
    random.shuffle(paths)
    return paths


# In[5]:


def extractSamples(variantFile):
    '''
    Extract samples in the VCF file to be added
    Input:
        variantFiles - files to be added (one per chromosome)
    '''
    #extract samples from vcf file
    request = 'bcftools query -l {}'.format(variantFile)
    output = subprocess.check_output(request, shell = True)
    samples = output.decode().split()
    return samples


# In[6]:


def extractPositions(variantFile):
    '''
    Extract the positions in the VCF file -> this is needed to efficiently extract data using
    bcftools query while specifying the start:end positions, without this bcftools would load 
    all positions at once which would be too large to fit in memory
    Input:
        variantFiles - files to be added (one per chromosome)
    '''
    #get chrom
    request = 'bcftools query -f \'%CHROM\n\' {} | head -n 1'.format(variantFile)
    output = subprocess.check_output(request, shell = True) 
    chrom = output.decode('utf-8').split('\n')[0].split('chr')[1]
    
    #get first position
    request = 'bcftools query -f \'%POS\n\' {} | head -n 1'.format(variantFile)
    output = subprocess.check_output(request, shell = True)
    positions = [output.decode('utf-8').split('\n')[0]]
    
    #get every 3rd position (MUST CHANGE IN FULL IMPLEMENTATION, CURRENTLY USING HEAD TO SPEED UP TESTING)
    request = 'bcftools query -f \'%POS\n\' {} | head -n 1000'.format(variantFile)
    output = subprocess.check_output(request, shell = True)
    positions.extend(output.decode('utf-8').split('\n'))
    
    #chunk the data to be uploaded
    positions_split = np.array_split(positions, 10)
    pos_regions = {}
    
    #create file that links together the chunked regions
    last = positions[0]
    for i, split in enumerate(positions_split):
        pos_regions[i] = [last, split[-1]]
        last = split[-1]
    return (chrom, pos_regions)


# In[7]:


def getAllelesPosition(row):
    '''
    Extract the alleles (genotypes) for each position. This is based on the number of alt alleles for that position.
    Heterozygous alleles are all grouped under one bucket e.g. 1|0 includes 1|0 and 0|1
    Input:
        row - row from the VCF file which includes the position, ref and alt alleles and sample data
    '''
    ##get the number of alternative alleles in each position
    split = row['alt'].split(',')
    alts = [str(i) for i in range(len(split)+1)]
    
    ##enumerate the possible genotypes by combining the different alleles
    combinations = []
    for subset in itertools.product(alts, repeat = 2):
        combinations.append("{}/{}".format(subset[0], subset[1]))
    
    #remove homo ref allele
    combinations.remove('0/0')
    
    #group genotypes if heterozygous (always done to larger number first then smaller so 1|0 or 2|1)
    alleles = []
    for gt in combinations:
        if gt[0] != gt[2]:
            alleles.append(gt) if int(gt[0]) < int(gt[2]) else alleles.append('{}/{}'.format(gt[2], gt[0]))
        ##add homozygous
        else:
            alleles.append(gt)
    alleles = list(set(alleles))
    return alleles


# In[8]:


def extractRelevantGenotypes(row):
    '''
    For each variant position extract the genotypes and samples associated with that genotype. Note does not get data for 0|0 samples
    Store in a dictionary that maps genotype:sample_ids
    Input:
        row - row from the VCF file which includes the position, ref and alt alleles and sample data
    '''
    ##get the possible alleles
    alleles = getAllelesPosition(row)
    relevantGenotypes = {}
    ##for every allele filter for that genotype
    for allele in alleles:
        genotypes = list(row[(row == allele)].index)
        ##restructure heterozygous genotypes to match that in rest of chain i.e. a|b wtih a > b
        if int(allele[0]) < int(allele[2]):
            x = list(allele)
            x[0], x[2] = x[2], x[0]
            allele = ''.join(x)
        relevantGenotypes[(row['ref'], row['alt'], allele.replace('/', '|'))] = genotypes
    return relevantGenotypes


# In[9]:


def extractVariant(variantFile, positions, pos_region, colnames, chrom):
    '''
    Within a position range in the VCF file, extract data on the genotypes and samples
    Input:
        positions - stores all the positions included in the VCF to be added to mapping chain
        pos_region - range of positions for bcftools to query (done to avoid overloading ram)
        colnames - person_ids included (and ref, alt columns)
        
    '''
    #get variant data ADD IN POSITIONS
    request = 'bcftools query -i \'FILTER="PASS"\' -f \'%POS %REF %ALT [ %GT]\n\' {}'.format(variantFile)
    output = subprocess.check_output(request, shell = True) 
    ##CHANGE THE NUMBER OF COLS USED WHEN HAVE FULL FILES
    df = pd.read_csv(BytesIO(output), delim_whitespace=True, usecols=[x for x in range(4)], names = colnames, index_col = 'pos')
    #keep track of positions added (for mapping)
    positions.extend(df.index)
    #get dictionary of samples associated with each genotype
    alt_genotypes = pd.DataFrame(df.apply(extractRelevantGenotypes, axis = 1))
    #wrangle the dataset to be genotype:sample_ids
    def filterows(row):
        return {k:v for k,v in row[0].items() if v}
    alt_genotypes = alt_genotypes.apply(filterows, axis = 1)
    ##filter out any genotypes that have no samples
    return alt_genotypes[alt_genotypes != {}], positions


# In[24]:


def publishToDataStream(chainName, multichainLoc, datadir, streamName, streamKeys, streamValues, publishVariant):
    '''
    Request to add data to multichain
    Input:
        streamName: chromosome of variant data being added
        streamKeys: position, ref, alt, genotype, # of samples with this gt, # of samples, # MAF
        streamValues: the person_ids that have this genotype
        publishVariant: Boolean on whether publishing variant data or MAF data
        
    '''
    ##if publishing variant data
    if publishVariant == True:
        publishCommand = [multichainLoc+'multichain-cli', 
            str('{}'.format(chainName)), 
            str('-datadir={}'.format(datadir)),
            'publish',
            str('structural_chrom_{}'.format(streamName)), 
            str('["{}", "{}", "{}", "{}", "{}", "{}", "{}"]'.format(streamKeys[0], streamKeys[1], streamKeys[2],
                                                                    streamKeys[3], streamKeys[4], streamKeys[5],
                                                                   streamKeys[6])),
            str('{}'.format(streamValues))]
    
    ##if publishing MAF data
    else:
        publishCommand = [multichainLoc+'multichain-cli', 
            str('{}'.format(chainName)), 
            str('-datadir={}'.format(datadir)),
            'publish',
            str('MAF_structural_chrom_{}'.format(streamName)), 
            str('["{}"]'.format(streamKeys)),
            str('{}'.format(streamValues))]

    procPublish = subprocess.Popen(publishCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    procPublish.wait()
    return


# In[25]:


def publishToDataStreams(chainName, multichainLoc, datadir, alt_genotypes, chrom, MAF, sample_size,  prevMAF, publishVariant):
    '''
    For a given position, get the genotypes associated and samples for each genotype and process data for insertion
    Every genotype is a separate insertion
    Input:
        alt_genotypes: The genotypes and samples attached to each genotype to be added (excludes the 0|0 gt)
        chrom: The chromosome VCF file belongs to
        MAF: the current MAF for that position and genotype
        sample_size: the # of samples being added
        prevMAF: the MAF prior to the new samples being added
        publishVariant: Boolean on whether publishing variant data or MAF data  
    '''
    ##for every position loop through data insertion
    for i, row in enumerate(alt_genotypes):
       ##for every genotype for that position loop through to process data for insertion
        for j, gt in enumerate(row):
            pos = str(alt_genotypes.index[i])
            ref,alt = gt[0], gt[1]
            allele = str(gt[2])
            ##if there is no previous MAF (as first data inserstion then just use current sample information)
            if isinstance(prevMAF, type(False)):
                count, total = len(row[gt]), sample_size
            ##else add in current and previous sample information together
            else:
                if (pos,ref,alt,allele) in prevMAF.index:
                    prevgt_count, prev_sample_size = int(prevMAF.loc[(pos,ref,alt,allele)]['count']), int(prevMAF.loc[(pos,ref,alt,allele)]['total'])
                    count, total = len(row[gt])+ prevgt_count, sample_size + prev_sample_size
                else:
                    count, total = len(row[gt]), sample_size

            ##keep track of the MAF for adding to the MAF streams later
            MAF[(pos,ref,alt,allele)]= [count, total, round(count/total,2)]
            ##process data for insertion
            streamName = chrom
            streamKeys = [pos, ref, alt, allele, MAF[(pos,ref,alt, allele)][0], MAF[(pos,ref,alt, allele)][1], MAF[(pos,ref,alt, allele)][2]]
            streamKeys[1] = streamKeys[1].replace("'",'')
            streamValues ='{'+'"json":{}'.format(row[gt]) +'}'#create JSON data object
            streamValues = streamValues.replace("'",'"')
            publishToDataStream(chainName, multichainLoc, datadir, streamName, streamKeys, streamValues, publishVariant)
    return 


# In[26]:


def extractPreviousMAF(chainName, multichainLoc, datadir, chrom):
    '''
    Get the existing MAF for all positions being added
    Input:
        chrom: The chromosome VCF file belongs to
    '''
    ##get the data for each position stored in the stream keys
    request = 'multichain-cli {} -datadir={} liststreamitems structural_chrom_{} false 9999999999999999'.format(chainName, datadir, chrom)
    output = subprocess.check_output(request, shell = True)
    ##parse the output (FIND A BETTER WAY HERE)
    output = output.decode(
        'utf-8').replace('\n', '').replace(
        '[','').replace(']','').replace('},    {','},, {').replace(
        ' ','').replace('        ','').replace(
        '    ','').replace('     ','').replace(
        '"keys":','"keys":[').replace(
        ',"offchain"','],"offchain"').replace(
        '"json":','"json":[').replace('},"confirmations"',']},"confirmations"').split(
        ',,')
    
    ##if there is no output then means no samples added for this position
    if output != ['']:
        output = [json.loads(x) for x in output]

        #extract old counts and get the latest one using blocktime
        df = pd.DataFrame(output).loc[:,['keys', 'blocktime']]
        #create dataframe of positions and MAF information  
        df[['position', 'ref', 'alt', 'allele', 'count', 'total', 'freq']] = pd.DataFrame(df['keys'].tolist())
        df = df.iloc[df.groupby(['position', 'allele'])['blocktime'].idxmax(),2:].set_index(['position', 'ref','alt','allele'])
        return df
    else:
        return False


# In[27]:


def publishMAF(chainName, multichainLoc, datadir, MAF, chrom):
    '''
    Publish the MAF taking into account all previous samples and newly added samples
    Input:
        MAF: dictionary that contains the MAF for every position-genotype
    '''
    ##convert dictionary to dataframe and create multi-index using position-genotype
    MAF_df = pd.DataFrame.from_dict(MAF, orient='index', columns = ['count', 'total', 'freq'])
    MAF_df.index = pd.MultiIndex.from_tuples(MAF_df.index)
    ##calculate the new MAF
    MAF_df['freq'] = MAF_df['count'].div(MAF_df['total'])
    ##bucket the MAFs into ranges and insert into stream
    ranges = [(0,0.05), (0.05,0.1), (0.1, 0.15), (0.15,0.2), (0.2,0.3), (0.3,0.4), (0.4,0.5), (0.5,1) ]
    for range in ranges:
        MAF_group = MAF_df[(MAF_df['freq'] > range[0]) &  (MAF_df['freq'] <= range[1])]
        streamName = chrom
        streamKeys = "{}-{}".format(range[0], range[1])
        streamValues ='{'+'"json":"{}"'.format(MAF_group['freq'].to_json().replace('"','').replace("'","")) +'}'#create JSON data object
        streamValues = streamValues.replace("'",'"')
        publishToDataStream(chainName, multichainLoc, datadir, streamName, streamKeys, streamValues, publishVariant = False)
    
    return


# In[28]:


def publishPositions(chainName, multichainLoc, datadir, positions, chrom):
    '''
    Publish the positions added into mapping stream
    Input:
        positions: list of positions added from vcf file
        chrom: chromosome positions come from
    '''
    streamName = 'mappingData_variants_structural'
    streamKeys = 'chrom_{}'.format(chrom)
    streamValues = '{'+'"json":{}'.format(positions) + '}'
    publishCommand = [multichainLoc+'multichain-cli', 
        str('{}'.format(chainName)), 
        str('-datadir={}'.format(datadir)),
        'publish',
        str('{}'.format(streamName)), 
        str('{}'.format(streamKeys)),
        str('{}'.format(streamValues))]
    procPublish = subprocess.Popen(publishCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    procPublish.wait()


# In[29]:


def publishVariants(chainName, multichainLoc, datadir, variantFiles):
    '''
    Publish the variant data
    Input:
        variantFiles - files to be added (one per chromosome)
    '''
    for variantFile in variantFiles:
        #extract samples
        samples = extractSamples(variantFile)
        colnames= ['pos', 'ref', 'alt']
        colnames.extend(samples)
        sample_size = len(colnames) - 3
        #split into groups and extract relevant variants
        chrom, pos_regions = extractPositions(variantFile)
        #store the MAF for all variants and extract old MAFs
        MAF = {}
        prevMAF_df = extractPreviousMAF(chainName, multichainLoc, datadir, chrom)
        ##store all positions
        positions = []
        for pos_region in list(pos_regions.values())[:len(pos_regions.values())-1]:
            ##insert variant data
            alt_genotypes, positions = extractVariant(variantFile, positions, pos_region, colnames, chrom)
            publishToDataStreams(chainName, multichainLoc, datadir, alt_genotypes, chrom, MAF, sample_size, prevMAF_df, publishVariant = True )
        ##publish MAF data
        publishMAF(chainName, multichainLoc, datadir, MAF, chrom)
        ##publish mapping of positions added
        publishPositions(chainName, multichainLoc, datadir, positions, chrom)
    return


# In[ ]:


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-dp", "--dataPath", type = str, help = "path of data to add to chain")
    parser.add_argument("-cn", "--chainName", help = "the name of the chain to store data", default = "chain1")
    parser.add_argument("-ml", "--multichainLoc", help = "path to multichain commands", default = "")
    parser.add_argument("-dr", "--datadir", help = "path to store the chain")
    parser.add_argument("-vf", "--variantfile", help = "variant files to add", default = "all")
    args = parser.parse_args()

    start = time.time()
    
    cpu = multiprocessing.cpu_count() * 2
    # cpu = 2
    print('CPUs available: {}'.format(cpu))
    
    try:
        subscribeToStreams(args.chainName, args.multichainLoc, args.datadir)
        print('Subscribed to streams') 
        
        paths = loadFilePaths(args.dataPath, args.variantfile)
        cpu = min(cpu, len(paths))
        paths_split = np.array_split(paths, cpu)
        processes = []
        for i in range(cpu):
            paths_split_ins = paths_split[i]
            p_ins = multiprocessing.Process(target=publishVariants, args = (args.chainName, args.multichainLoc,
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
        sys.stderr.write("\nERROR: Failed stream publishing. Please try again.\n")
        quit()
        

if __name__ == "__main__":
    main()

