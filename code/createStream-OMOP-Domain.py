#!/usr/bin/env python
# coding: utf-8

# In[2]:


'''
createStream-mimic-Domain.py
Creates streams for the domain view of the clinical data in the blockchain
Usage: $ python createStream-mimic-Domain.py -cn=<chain name> -dr=<Chain path> -hp=<Vocabulary path> -dp=<Data tables path> 
modified by AE 02/2022
'''

import sys
import time
import math
import binascii
import argparse
import subprocess
from subprocess import Popen, PIPE
import json
import os
import psutil
import time
import pandas as pd
import warnings
import multiprocessing
import random
warnings.simplefilter("ignore")


# In[ ]:


def parseTables(tables):
    '''
    takes the input tables and parses them to a list for use in code
    '''

    #if user wants to input to all tables
    if tables == 'all':
        tables = ['condition_occurrence', 'drug_exposure',
              'measurement', 'observation', 'procedure_occurrence', 'visit_occurrence']
    
    #parse the tables that are part of user input
    else:
        tables = str.split(tables.replace(" ",""), ',')
    return tables


# In[ ]:


def loadConcepts(dataPath, table):
    '''
    Load the data and extract the unique OMOP codes used in the table
    This is necessary to create streams relevant to the dataset
    
    Inputs:
        dataPath - path for the data that is being added
        table - the specific table being added
    '''

    dataPath = '{}{}.csv'.format(dataPath, table) 
    df = pd.read_csv(dataPath)
    unique_codes = list(df.iloc[:,2].unique())
    concept_type = df.columns[2]
    concept_type = concept_type.split('_')[0]
    return unique_codes, concept_type, df


# In[8]:


def loadSuperConcepts(hierarchyPath, dataPath, table):
    '''
    Load the super concepts in the vocabulary hierarchy that will make up the streams
    Concepts are then assigned to the relevant super concept based on the hierarchy
    
    Inputs:
        hierarchyPath - path for the vocabulary hierarchy (used to assign concepts to super concepts)
        dataPath - path for the data that is being added
        table - the specific table being added
    '''
    
    #load the list of super concepts for each domain from the dictionary
    dataPath = dataPath + "table_super_concepts.txt"
    with open(dataPath) as f:
        load = f.read()
    super_concept_dict = json.loads(load)
    
    #parse each super concept for the specific table being uploaded
    super_concept = str.split(super_concept_dict[table].replace(" ",""), ',')
    super_concept = [int(x) for x in super_concept]
    
    #load the hierarchy tables and clean
    df = pd.read_csv(hierarchyPath, sep = '|')
    df.drop_duplicates(inplace = True)
    df_clean = df[df['min_levels_of_separation'] != 0]
    
    #filter for the concepts one level below the super concept
    df_hierarchy = df_clean[(df_clean['ancestor_concept_id'].isin(super_concept))]
    df_stream_concepts = df_hierarchy[(df_hierarchy['min_levels_of_separation']  == 1)]
   
    #use these concepts to group all other concepts
    df_stream_ancestor = df_clean[df_clean['ancestor_concept_id'].isin(df_stream_concepts['descendant_concept_id'])]
    
    return df_stream_ancestor.sort_values('ancestor_concept_id')


# In[6]:


def conceptStreamMapping(hierarchyPath, dataPath, table):
    '''
    Create a dictionary which indicates which stream each concept in the data could go
    Note there may be multiple streams that a concept is eligible to go to 
    Inputs:
        hierarchyPath - path for the vocabulary hierarchy (used to assign concepts to super concepts)
        dataPath - path for the data that is being added
        table - the specific table being added
    '''
    stream_mapping = {}
    unique_concepts, _, _ = loadConcepts(dataPath, table)
    df_stream_ancestor = loadSuperConcepts(hierarchyPath, dataPath, table)
    
    #for every concept in the dataset assign the potential super concepts that it could go to
    for concept in unique_concepts:
            super_concepts = list(df_stream_ancestor['ancestor_concept_id'][df_stream_ancestor['descendant_concept_id'] == concept].values)
            stream_mapping[concept] = super_concepts
    return stream_mapping    


# In[11]:


def assignConceptStream(hierarchyPath, dataPath, table):
    '''
    Assign a specific stream for each concept, taking a random stream to help with stream size balancing
    Inputs:
        hierarchyPath - path for the vocabulary hierarchy (used to assign concepts to super concepts)
        dataPath - path for the data that is being added
        table - the specific table being added
    '''
    _, _, data_df = loadConcepts(dataPath, table)
    stream_mapping = conceptStreamMapping(hierarchyPath, dataPath, table)
    
    ##assign concept stream randomly from those that concept could go to; if no super concepts attached to concept then place in stream '0'
    concept_stream = []
    random.seed(1)
    for concept in data_df.iloc[:,2]:
        try:
            concept_stream.append(random.choice(stream_mapping[concept]))
        except:
            concept_stream.append(0)
    return concept_stream 


# In[1]:


def numStreamBuckets(hierarchyPath, dataPath, table):
    '''
    Determine the number of buckets that a super concept needs, more than 1 stream will be needed for a concept if 1.2gb+ data is in a stream 
    Inputs:
        hierarchyPath - path for the vocabulary hierarchy (used to assign concepts to super concepts)
        dataPath - path for the data that is being added
        table - the specific table being added
    '''
    ##get assignments per stream
    stream_assignments = assignConceptStream(hierarchyPath, dataPath, table)
    ##get the size of each stream (based on size of a single entry being 900kb)
    stream_size = {i:stream_assignments.count(i) for i in list(set(stream_assignments))}
    stream_max_count = round(1.2e+9 / 900)
    num_stream_buckets = {i:math.ceil(stream_size[i]/stream_max_count) for i in stream_size}
    ##assign stream 0 one bucket (typcially just a few concepts that are rare and without a parent)
    num_stream_buckets[0] = 1
    return num_stream_buckets


# In[12]:


#Given a chain name and the name of the new stream, makes that stream
def makeStream(chainName, streamName, multichainLoc, datadir):
    createStreamCommand=multichainLoc+'multichain-cli {} -datadir={} create stream {} true'.format(chainName,datadir,streamName)
    subscribeStreamCommand=multichainLoc+'multichain-cli {} -datadir={} subscribe {}'.format(chainName,datadir,streamName)

    #make stream of name StreamName
    procCreate = subprocess.Popen(createStreamCommand.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    procCreate.wait()

    #subscribe to the stream
    procSubscribe = subprocess.Popen(subscribeStreamCommand.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    procCreate.wait()
    return


# In[22]:


#Making the streams that will be used to store the file using concept id
def createStreams(chainName, multichainLoc, datadir, hierarchyPath, dataPath, table):
    #extract concepts
    _, concept_type, _ = loadConcepts(dataPath, table)
    num_stream_buckets = numStreamBuckets(hierarchyPath, dataPath, table)
    #data stream by for mapping
    makeStream(chainName, "mappingData_clinical", multichainLoc, datadir) #header, other things for the whole file
    #data stream for concepts
    for concept in num_stream_buckets.keys():
        for num_buckets in range(num_stream_buckets[concept]):
            makeStream(chainName, "{}_id_{}_bucket_{}".format(concept_type, concept, num_buckets+1), multichainLoc, datadir)
    print("Created {} streams".format(table))
    return


# In[1]:


######################################################################################################
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-cn", "--chainName", help = "the name of the chain to store data", default = "chain1")
    parser.add_argument("-ml", "--multichainLoc", help = "path to multichain commands", default = "")
    parser.add_argument("-dr", "--datadir", help = "path to store the chain")
    parser.add_argument("-hp", "--hierarchyPath", help = "path to vocabulary file")
    parser.add_argument("-dp", "--dataPath", type = str, help = "path of data to add to chain")
    parser.add_argument("-tb", "--tables", help = "tables to add", default = "all")
    args = parser.parse_args()

    start = time.time()
    tables = parseTables(args.tables)
    processes = []
    try:
        print("--STREAM CONSTRUCTION--")
        for table in tables:
            p = multiprocessing.Process(target=createStreams, args = (args.chainName, args.multichainLoc,
                                                                      args.datadir, args.hierarchyPath, args.dataPath,  table))
            processes.append(p)
            p.start()
            
        for process in processes:
            process.join()
  
    except:
        sys.stderr.write("\nERROR: Failed stream creation. Please try again.\n")
        quit()
    
            
    end = time.time()
    e = int(end - start)
    print('\n\n Time elapsed:\n\n')
    print( '{:02d}:{:02d}:{:02d}'.format(e // 3600, (e % 3600 // 60), e % 60))

    process = psutil.Process(os.getpid())
    print('\n\n Total memory in bytes:\n\n')
    print(process.memory_info().rss)



        


# In[ ]:


if __name__ == "__main__":
    main()

