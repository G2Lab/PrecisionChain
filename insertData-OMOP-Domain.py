#!/usr/bin/env python
# coding: utf-8

# In[2]:


'''
insertData-mimic-Domain.py
Loads tab-separated text file data onto an existing data stream on an existing multichain
Usage: $ python insertData-mimic-Domain.py -cn=<Chain name> -dr=<Chain path> hp=<Vocabulary path> -dp=<Data tables path> 
modified by AE 02/2022
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
warnings.simplefilter("ignore")


# In[3]:


def parseTables(tables):
    '''
    takes the input tables and parses them to a list for use in code
    '''
    #parse the tables that are part of user input (ADD IN MEASUREMENT LATER REMOVED TO SPEED UP INSERTION)
    if tables == 'all':
        tables = ['condition_occurrence', 'drug_exposure', 'device_exposure', 'observation', 'procedure_occurrence', 'specimen', 'visit_occurrence']
    else:
        tables = str.split(tables.replace(" ",""), ',')
    return tables


# In[4]:


def loadConcepts(dataPath, table):
    '''
    Load the data and extract the unique OMOP codes used in the table
    This is necessary to create streams relevant to the dataset
    
    Inputs:
        dataPath - path for the data that is being added
        table - the specific table being added
    '''
    dataPath = '{}/{}.csv'.format(dataPath, table) 
    df = pd.read_csv(dataPath)
    unique_codes = list(df.iloc[:,2].unique())
    concept_type = df.columns[2]
    concept_type = concept_type.split('_')[0]
    return unique_codes, concept_type, df


# In[5]:


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


# In[7]:


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


# In[8]:


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


# In[9]:


def processKeys(dataPath, table):
    '''
    For a given table load the name of the keys that are used to insert data , this is fixed to ensure consistency across data insertions
    Keys will include the concept id and date
    Inputs:
        dataPath - path for the data that is being added
        table - the specific table being added
    '''
    dataPath = dataPath + 'table_keys.txt'
    with open(dataPath) as f:
        load = f.read()
    keys_dict = json.loads(load)
    keys_processed = str.split(keys_dict[table].replace(" ",""), ',')
    return keys_processed


# In[10]:


def processDateKey(keys_df):
    '''
    Process the date key to have full date, month-year and years (necessary for querying as cant do range queries)
    Inputs:
        keys_df: the dataframe of the extracted keys
    '''
    date_df = keys_df.filter(like='date')
    keys_df['year'] = pd.to_datetime(date_df.iloc[:,0]).dt.to_period('Y')
    keys_df['month_year'] = pd.to_datetime(date_df.iloc[:,0]).dt.to_period('M')
    return keys_df
    


# In[11]:


#load the data and organise the data into the stream_concept, keys and values
def loadData(dataPath, table, keys):
    '''
    For a given table, load the data and extract the concepts, keys, values
    Inputs:
        dataPath - path for the data that is being added
        table - the specific table being added
    '''
    #load data and extract relevant data
    dataPath = '{}/{}.csv'.format(dataPath, table) 
    df = pd.read_csv(dataPath)
    #get concept
    concept_type = df.columns[2]
    concept_type = concept_type.split('_')[0]
    #get keys
    keys_df = df.loc[:,keys]
    keys_df = processDateKey(keys_df)
    #get values
    values_df = df
    return concept_type, keys_df, values_df


# In[12]:


#Given a chain name and the name of the new stream, subscribe to that stream
def subscribeToStream(chainName, streamName, multichainLoc, datadir):
    subscribeStreamCommand=multichainLoc+'multichain-cli {} -datadir={} subscribe {}'.format(chainName,datadir,streamName)

    #subscribe to the stream
    procSubscribe = subprocess.Popen(subscribeStreamCommand.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    procSubscribe.wait()
    return


# In[13]:


#Making the streams that will be used to store the file
#making the stream bins
#hash for location
def subscribeToStreams(chainName, multichainLoc, datadir, hierarchyPath, dataPath, table):
    _, concept_type, _ = loadConcepts(dataPath, table)
    #data stream by for mapping
    subscribeToStream(chainName, "mappingData_clinical", multichainLoc, datadir) #header, other things for the whole file
    #subscribe to concept streams
    stream_buckets = numStreamBuckets(hierarchyPath, dataPath, table)
    for concept in stream_buckets.keys():
        for num_buckets in range(stream_buckets[concept]):
            subscribeToStream(chainName, "{}_id_{}_bucket_{}".format(concept_type, concept, num_buckets+1), multichainLoc, datadir)
    return


# In[14]:


def processTable(chainName, multichainLoc, datadir, hierarchyPath, dataPath, table):
    '''
    Process the data table to get ready for insertion - includes extracting the streams and stream buckets, keys, and values
    Inputs:
        hierarchyPath - path for the vocabulary hierarchy (used to assign concepts to super concepts)
        dataPath - path for the data that is being added
        table - the specific table being added
    '''
    #get # of buckets for each stream
    stream_buckets = numStreamBuckets(hierarchyPath, dataPath, table)
    #Get keys and values to publish
    keys = processKeys(dataPath, table)
    concept_type, keys_df, values_df = loadData(dataPath, table, keys)
    #Assign each concept to a super concept stream
    concept_stream = assignConceptStream(hierarchyPath, dataPath, table)
    values_df['concept_stream'] = concept_stream
    #Assign a bucket to each entry
    df = pd.DataFrame()
    for stream_concept in stream_buckets.keys():
        stream_df = values_df[values_df['concept_stream'] == stream_concept].sort_values(by = values_df.columns[2])
        stream_df['bucket'] = pd.qcut(stream_df.iloc[:,2], stream_buckets[stream_concept], labels = [x for x in range(stream_buckets[stream_concept])])
        df = df.append(stream_df)
    return concept_type, df, keys_df   


# In[15]:


def assignedStreamDictionary(df):
    '''
    create dictionary of the streams created and buckets used for each stream (super concept)
    create dictionary of the concepts:stream/bucket mapping that were actually used during data insertion
    Inputs:
        df = dataset being inserted
    '''
    #create stream mapper for what buckets each super concept has
    stream_dictionary ={}
    #create mapping of each concept to super-concept stream
    stream_concept_dictionary = {}
    #for each concept in the dataset find which super concept it was assigned to and which buckets were used, add these to a dictionary (to be added to the mapping stream)
    for concept in df.iloc[:,2].unique():
        concept = int(concept)
        stream_concept = int(df[df.iloc[:,2] == concept].concept_stream.iloc[0])
        stream_dictionary[str(stream_concept)] = list(df[df.iloc[:,2] == concept].bucket.unique())
        stream_concept_dictionary[concept] = {}
        stream_concept_dictionary[concept][stream_concept] = list(df[df.iloc[:,2] == concept].bucket.unique())
    return stream_dictionary, stream_concept_dictionary


# In[16]:


#Publish mapping data to the mapping stream
def publishToMappingStream(chainName, multichainLoc, datadir, concept, concept_type, stream):
    '''
    Publish data to the clinical mapping stream - indicates what stream and bucket each concept is stored under
    Inputs:
        concept - OMOP concept recorded in the data table
        concept_type - which data table the concept is from
        stream - what stream (super concept it is under)
    '''
    streamName = "mappingData_clinical"
    streamKeys = concept, concept_type, list(stream.keys())[0], list(stream.values())[0]
    streamValues ='{"json":'+json.dumps(stream)+'}' #create JSON 

    publishCommand = [multichainLoc+'multichain-cli', 
        str('{}'.format(chainName)), 
        str('-datadir={}'.format(datadir)),
        'publish',
        str('{}'.format(streamName)), 
        str('["{}","{}","{}","{}"]'.format(streamKeys[0],streamKeys[1],streamKeys[2],streamKeys[3])),
        str('{}'.format(streamValues))
            ]

    procPublish = subprocess.Popen(publishCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    procPublish.wait()
    return


# In[17]:


def publishToDataStream(row, chainName, multichainLoc, datadir, concept_type, keys_df):
    '''
    Publish data row to the appropriate stream - streamKeys include concept_id and dates
    Inputs:
        row - the recorded value for that data row
        concept_type - which data table the concept is from
        keys_df - the keys datframe for key insertion
    '''
    streamName = row['concept_stream']
    streamBucket = row['bucket']
    streamKeys = keys_df.loc[row.name]
    streamValues ='{"json":'+row.to_json().strip('[').strip(']')+'}' #create JSON and remove brackets
    
    publishCommand = [multichainLoc+'multichain-cli', 
        str('{}'.format(chainName)), 
        str('-datadir={}'.format(datadir)),
        'publish',
        str('{}_id_{}_bucket_{}'.format(concept_type,streamName, streamBucket+1)), 
        str('["{}", "{}", "{}", "{}", "{}"]'.format(streamKeys[0], streamKeys[1], streamKeys[2], streamKeys[3], streamKeys[4])),
        str('{}'.format(streamValues))]
    procPublish = subprocess.Popen(publishCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    procPublish.wait()
    return


# In[18]:


def publishToDataStreams(chainName, multichainLoc, datadir, data_ins):
    '''
    Publish data for the whole dataset
    Inputs:
        data_ins - the dataset to be inserted
        keys_df - the keys datframe for key insertion
    '''
    concept_type, df_split_ins, keys_df = data_ins 
    df_split_ins.apply(publishToDataStream, args = (chainName, multichainLoc, datadir,concept_type, keys_df), axis = 1)
    return


# In[ ]:


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
    cpu = multiprocessing.cpu_count()
    print('CPUs available: {}'.format(cpu))
    
    tables = parseTables(args.tables)
    
    processes_sub = []
    processes_map = []
    processes_ins = []
    try:
        for table in tables:
            p_sub = multiprocessing.Process(target=subscribeToStreams, args = (args.chainName, 
                                                                               args.multichainLoc, args.datadir, args.hierarchyPath, args.dataPath, table))
            processes_sub.append((p_sub, table))
            p_sub.start()
            
        for process, table in processes_sub:
            process.join()
            print('Subscribed to {} table'.format(table))
        
        for table in tables:
            #process the data to be inserted into the relevant pieces of information
            concept_type, df, keys_df = processTable(args.chainName, args.multichainLoc, args.datadir, args.hierarchyPath, args.dataPath, table)
            #record the streams used
            stream_dictionary, stream_concept_dictionary = assignedStreamDictionary(df)
            #publish the mapping stream with the streams and buckets created
            publishToMappingStream(args.chainName, args.multichainLoc, args.datadir, 'StreamsUsed', concept_type, stream_dictionary)
            #for each concept insert data to mapping stream
            for concept in stream_concept_dictionary.keys():
                stream = stream_concept_dictionary[concept]
                p_map = multiprocessing.Process(target=publishToMappingStream, args = (args.chainName, args.multichainLoc, args.datadir,
                                                                               concept, concept_type, stream))
                processes_map.append(p_map)
                p_map.start()
            
            for process in processes_map:
                process.join()
            
            print('Published mapping stream for {}'.format(table))
            
            ##split the dataset based on number of CPUs available and insert
            df_split = np.array_split(df, cpu)  
            for i in range(cpu):
                df_split_ins = df_split[i]
                data_ins = (concept_type, df_split_ins, keys_df)
                p_ins = multiprocessing.Process(target=publishToDataStreams, args = (args.chainName,
                                                                                    args.multichainLoc, args.datadir, data_ins))
                processes_ins.append(p_ins)
                p_ins.start()
            
            for process in processes_ins:
                p_ins.join()
            print('Inserted to {}'.format(table))
            end = time.time()
            e = int(end - start)
            print('\n\n Time elapsed:\n\n')
            print( '{:02d}:{:02d}:{:02d}'.format(e // 3600, (e % 3600 // 60), e % 60))
    
    except:
        sys.stderr.write("\nERROR: Failed stream publishing. Please try again.\n")
        quit()

    end = time.time()
    e = int(end - start)
    print('\n\n Time elapsed:\n\n')
    print( '{:02d}:{:02d}:{:02d}'.format(e // 3600, (e % 3600 // 60), e % 60))

    process = psutil.Process(os.getpid())
    print('\n\n Total memory in bytes:\n\n')
    print(process.memory_info().rss)

if __name__ == "__main__":
    main()


# tables = parseTables('all')
# chainName, multichainLoc, datadir = 'mimic2', '', '/gpfs/commons/groups/gursoy_lab/aelhussein/multichain' 
# hierarchyPath = '/gpfs/commons/groups/gursoy_lab/aelhussein/Data/Vocabulary/CONCEPT_ANCESTOR.csv'
# dataPath = '/gpfs/commons/groups/gursoy_lab/aelhussein/Data/mimic_iv_omop/Clinical_data/'
# concept_type, df, keys_df  = processTable(chainName, multichainLoc, datadir, hierarchyPath, dataPath, tables[0])
# stream_dictionary, stream_concept_dictionary = streamDictionary(df)
