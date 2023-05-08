#!/usr/bin/env python
# coding: utf-8

# In[1]:


'''
insertData-mimic-Person.py
Inserts data for the person view of the clinical data in the blockchain
Usage: $ python insertData-mimic-Person.py -cn=<chain name> -dr=<Chain path> -dp=<Data path> -pp=<Person data path>
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
warnings.simplefilter(action='ignore')
import persons


# In[2]:


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


# In[3]:


def loadPatients(path):
    '''
    load the patients being inserted and split into 20 groups based on person_id
    '''
    df = pd.read_csv(path)
    df['stream'] = pd.qcut(df['person_id'], q = 20, labels =np.arange(1,21,1))
    df['person_id'] = df['person_id'].astype(str)
    # Loads only the people/NTASKS number of people
    df = df[df.person_id.isin(persons.getChunk())]
    return df


# In[4]:

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
    


# In[5]:


def loadData(path, table):
    '''
    For a given table, load the data and extract the concepts, keys, values
    Inputs:
        path - path for the data that is being added
        table - the specific table being added
    '''
    df = pd.read_csv(path + table +'.csv')
    df['person_id'] = df['person_id'].astype(str)
    # Loads only the people/NTASKS number of people
    df = df[df.person_id.isin(persons.getChunk())]
    concept_type = df.columns[2]
    concept_type = concept_type.split('_')[0]
    
    #take the patient_id, concept_id, date as keys
    keys_df = df.iloc[:,1:3]
    keys_df['date'] = df.filter(like='date').iloc[:,0]
    keys_df = processDateKey(keys_df)
    return (concept_type, keys_df, df)


# In[6]:


def loadDataStream(stream, data_df, person_df):
    '''
    For a given stream, filter the patients in the stream in the dataset
    Inputs:
        stream - the stream (1-20) being inserted into
        data_df - the dataset to be inserted
        person_df - the table of person_ids and assigned streams
    '''
    ##get the concept, keys dataframe and dataset from the tuple of data_df
    concept_type, keys_df, df = data_df
    ##filter for the patients in the person stream
    person_df_stream = person_df[person_df['stream'] == stream]
    person_id = person_df_stream['person_id'].unique()
    ##filter for the dataset belonging to those patients
    df_stream = df[df['person_id'].isin(person_id)]
    keys_df_stream = keys_df[df['person_id'].isin(person_id)]
    rows = df_stream.shape[0]
    ##save the relevant data in a tuple
    data_df_stream = (rows, concept_type,keys_df_stream, df_stream)
    return person_df_stream, data_df_stream


# In[7]:


#Given a chain name and the name of the new stream, subscribe to that stream
def subscribeToStream(chainName, streamName, multichainLoc, datadir):
    subscribeStreamCommand=multichainLoc+'multichain-cli {} -datadir={} subscribe {}'.format(chainName,datadir,streamName)
    #subscribe to the stream
    procSubscribe = subprocess.Popen(subscribeStreamCommand.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    procSubscribe.wait()
    return


# In[8]:


#subscribe to person streams
def subscribeToStreams(chainName, multichainLoc, datadir):
    subscribeToStream(chainName, "mappingData_person", multichainLoc, datadir) 
    subscribeToStream(chainName, "person_demographics", multichainLoc, datadir)
    for i in range(1, 21):
        subscribeToStream(chainName, "person_stream_{}".format(i), multichainLoc, datadir)
    return


# In[9]:


def publishToMappingStreams(chainName, multichainLoc, datadir, person_df): 
    '''
    The mapping stream records which stream each person is added to
    The demographic stream contains demographic data on the patient (based on demographic OMOP table)
    Inputs:
        person_df - the table of person_ids and assigned streams
    '''
    #insert into the mapping and demographics stream
    for person_id in person_df['person_id']:
            #mappingstream
            streamName = "mappingData_person"
            streamKeys = person_id
            streamValues ='{"json":'+json.dumps(int(person_df['stream'][person_df['person_id'] == person_id].iloc[0]))+'}' #create JSON 
            publishCommand = [multichainLoc+'multichain-cli', 
                str('{}'.format(chainName)), 
                str('-datadir={}'.format(datadir)),
                'publish',
                str('{}'.format(streamName)), 
                str('["{}"]'.format(streamKeys)),
                str('{}'.format(streamValues))]
            procPublish = subprocess.Popen(publishCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            procPublish.wait()

            #demographics stream
            streamName = "person_demographics"
            row = person_df[person_df['person_id']==person_id]
            if row['race_concept_id'].iloc[0] == 0:
                race = row['ethnicity_source_value'].iloc[0]
            else:
                race = row['race_source_value'].iloc[0]
            streamKeys = row['person_id'].iloc[0], row['gender_source_value'].iloc[0], race

            streamValues ='{"json":'+row.to_json(orient = 'records').strip('[').strip(']')+'}'
            publishCommand = [multichainLoc+'multichain-cli', 
                str('{}'.format(chainName)), 
                str('-datadir={}'.format(datadir)),
                'publish',
                str('{}'.format(streamName)), 
                str('["{}","{}","{}"]'.format(streamKeys[0],streamKeys[1],streamKeys[2])),
                str('{}'.format(streamValues))
                    ]
            print(" ".join(publishCommand))
            procPublish = subprocess.run(publishCommand, capture_output=True, text=True)
            print(procPublish.stdout)
            print(procPublish.stderr)
    return


# In[15]:


def publishToDataStream(row, chainName, multichainLoc, datadir, concept_type, keys_df, person_df):
    '''
    Publish data row to the appropriate stream - streamKeys include concept_id and dates
    Inputs:
        row - the recorded value for that data row
        concept_type - which data table the concept is from
        keys_df - the keys datframe for key insertion
    '''
    streamName = person_df['stream'][person_df['person_id'] == row.person_id].iloc[0]
    streamKeys = keys_df.loc[row.name]
    streamValues ='{"json":'+row.to_json().strip('[').strip(']')+'}' #create JSON and remove brackets
    
    publishCommand = [multichainLoc+'multichain-cli', 
        str('{}'.format(chainName)), 
        str('-datadir={}'.format(datadir)),
        'publish',
        str('person_stream_{}'.format(streamName)), 
        str('["{}","{}", "{}", "{}", "{}", "{}"]'.format(concept_type, streamKeys[0], streamKeys[1], streamKeys[2], streamKeys[3], streamKeys[4])),
        str('{}'.format(streamValues))]
    procPublish = subprocess.Popen(publishCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    procPublish.wait()
    return


# In[11]:


def publishToDataStreams(chainName, multichainLoc, datadir, data_ins, person_df):
    '''
    Publish data for the whole dataset
    Inputs:
        data_ins - the dataset to be inserted
        person_df - dataset of person_ids and stream assignment
    '''
    concept_type, df_split_ins, keys_df = data_ins 
    df_split_ins.apply(publishToDataStream, args = (chainName, multichainLoc, datadir,concept_type, keys_df, person_df), axis = 1)
    return



# In[ ]:

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-dp", "--dataPath", type = str, help = "path of data to add to chain")
    parser.add_argument("-cn", "--chainName", help = "the name of the chain to store data", default = "chain1")
    parser.add_argument("-ml", "--multichainLoc", help = "path to multichain commands", default = "")
    parser.add_argument("-dr", "--datadir", help = "path to store the chain")
    parser.add_argument("-pp", "--personPath", help = "path to patient demographics file")
    parser.add_argument("-tb", "--tables", help = "tables to add", default = "all")
    args = parser.parse_args()

    start = time.time()
    
    tables = parseTables(args.tables)
    processes = []
    # cpu = multiprocessing.cpu_count()
    cpu = 16
    print('CPUs available: {}'.format(cpu))
    
    try:
        subscribeToStreams(args.chainName, args.multichainLoc, args.datadir)
        print('Subscribed to streams') 
        
        person_df = loadPatients(args.personPath)
        publishToMappingStreams(args.chainName, args.multichainLoc, args.datadir, person_df)
        print('Published to Mapping streams') 
        
        for table in tables:
            concept_type, keys_df, data_df = loadData(args.dataPath, table)
            ##split the dataset based on number of CPUs
            data_df_split = np.array_split(data_df, cpu)
            keys_df_split = np.array_split(keys_df,cpu)
            ##insert via multiprocessing
            for i in range(cpu):
                data_df_split_ins = data_df_split[i]
                keys_df_split_ins = keys_df_split[i]
                rows = data_df_split_ins.shape[0]
                data_ins = (concept_type, data_df_split_ins, keys_df_split_ins)
                p_ins = multiprocessing.Process(target=publishToDataStreams, args = (args.chainName, args.multichainLoc,
                                                                                         args.datadir, data_ins, person_df))
                processes.append(p_ins)
                p_ins.start()
            
            
            for process in processes:
                process.join()
            print('Inserted {}'.format(table))
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

