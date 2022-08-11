#!/usr/bin/env python
# coding: utf-8

# In[2]:


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
import ast
import multiprocessing
import glob
import pandas as pd
import json
import warnings
import multiprocessing
from datetime import datetime 
warnings.simplefilter(action='ignore', category=FutureWarning)


# In[ ]:


def parseKeys(cohortKeys, searchKeys):
    cohortKeys = str.split(cohortKeys.replace(" ",""), ',')
    searchKeys = str.split(searchKeys.replace(" ",""), ',')
    return cohortKeys, searchKeys


# In[ ]:


#Given a chain name subscribe to audit log stream to ensure query is recorded
def subscribeToStream(chainName, multichainLoc, datadir):
    subscribeStreamCommand=multichainLoc+'multichain-cli {} -datadir={} subscribe audit_log'.format(chainName,datadir)
    #subscribe to the stream
    procSubscribe = subprocess.Popen(subscribeStreamCommand.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    procSubscribe.wait()
    return


# ## Cohort building

# In[ ]:


def queryMappingStream(chainName, multichainLoc, datadir, keys):
    concepts = []
    for key in keys: 
        queryCommand=multichainLoc+'multichain-cli {} -datadir={} liststreamkeyitems mappingData_clinical {}'.format(chainName, datadir, key)
        items = subprocess.check_output(queryCommand.split())
        info = tuple(json.loads(items, parse_int= int)[0]['keys'])
        concepts.append(info)
    return list(set(concepts))


# In[ ]:


def extractDataStream(chainName, multichainLoc, datadir, cohortKeys):
    matches = queryMappingStream(chainName, multichainLoc, datadir, cohortKeys)
    concept_domain, concept_stream, concept_bucket = matches[0][1], matches[0][2], ast.literal_eval(matches[0][3])
    return concept_domain, concept_stream, concept_bucket


# In[ ]:


def extractPersonIDs(chainName, multichainLoc, datadir, cohortKeys):
    concept_domain, concept_stream, concept_bucket = extractDataStream(chainName, multichainLoc, datadir, cohortKeys)
    matches = []
    for bucket in concept_bucket:
        ##potentially change to liststreamkeyitems
        queryCommand=multichainLoc+'multichain-cli {} -datadir={} liststreamkeyitems {}_id_{}_bucket_{} {}'.format(chainName, datadir,
                                                                                                               concept_domain, concept_stream,
                                                                                                               bucket+1, cohortKeys[0])
        items = subprocess.check_output(queryCommand.split())
        matches.extend(json.loads(items, parse_int= int))
    ids = []
    for match in matches:
        ids.append(int(match['keys'][0]))
    return(list(set(ids)))


# In[ ]:


def extractPersonStreams(chainName, multichainLoc, datadir, cohortKeys):
    person_ids = extractPersonIDs(chainName, multichainLoc, datadir, cohortKeys)
    person_streams = {}
    for person_id in person_ids:
        queryCommand=multichainLoc+'multichain-cli {} -datadir={} liststreamkeyitems mappingData_person {}'.format(chainName, datadir, person_id)
        items = subprocess.check_output(queryCommand.split())
        matches = json.loads(items, parse_int= int) 
        person_streams[person_id] = matches[0]['data']['json']
    return person_streams


# ## queries

# In[ ]:

def queryDemographics(chainName, multichainLoc, datadir, cohortKeys):
    ##check if search from running query or person specific query
    if isinstance(cohortKeys,str):
        personids = [int(id) for id in cohortKeys.split(',')]
    else:
        personids = cohortKeys
    persons_df = pd.DataFrame()
    
    for personid in personids:
        queryCommand=multichainLoc+'multichain-cli {} -datadir={} liststreamkeyitems person_demographics {}'.format(chainName, datadir, personid)
        items = subprocess.check_output(queryCommand.split())
        json_item = json.loads(items, parse_int= int)[1]['data']['json']
        person_df = pd.DataFrame.from_dict(json_item, orient = 'index').T
        persons_df = pd.concat([persons_df,person_df])
    persons_df.set_index('person_id', inplace = True)
    publishToAuditstream(chainName, multichainLoc, datadir, queryCommand)
    persons_json = persons_df.to_json(orient = 'index')
    print(persons_json)
    return persons_json


def queryDomainStream(chainName, multichainLoc, datadir, cohortKeys, searchKeys):
    person_ids = extractPersonIDs(chainName, multichainLoc, datadir, cohortKeys)
    searchStreams = queryMappingStream(chainName, multichainLoc, datadir, searchKeys)

    for stream in searchStreams:
        buckets = ast.literal_eval(stream[3])
        for bucket in buckets:
            for person_id in person_ids:
                try:
                    queryCommand = multichainLoc+'multichain-cli {} -datadir={} liststreamqueryitems {}_id_{}_bucket_{} {{"keys":["{}","{}"]}}'.format(chainName, datadir,
                                                                                                stream[1], stream[2], bucket+1, person_id, stream[0])
                    items = subprocess.check_output(queryCommand.split())
                    matches = json.loads(items, parse_int= int)
                    print(matches)
                    publishToAuditstream(chainName, multichainLoc, datadir, queryCommand)
                except:
                    pass


# In[ ]:


def queryPersonStreams(chainName, multichainLoc, datadir, cohortKeys, searchKeys):
    person_streams = extractPersonStreams(chainName, multichainLoc, datadir, cohortKeys)
    for person_id in person_streams.keys():
        queryCommand = multichainLoc+'multichain-cli {} -datadir={} liststreamkeyitems person_stream_{} {}'.format(chainName, datadir,
                                                                                    person_streams[person_id], person_id)
        items = subprocess.check_output(queryCommand.split())
        matches = json.loads(items, parse_int= int)
        print(matches)
        publishToAuditstream(chainName, multichainLoc, datadir, queryCommand)
    return


# In[ ]:


def query(chainName, multichainLoc, datadir, cohortKeys, searchKeys):
    
    cohortKeys, searchKeys = parseKeys(cohortKeys, searchKeys)
    if searchKeys[0] == 'demographics':
        queryDemographics(chainName, multichainLoc, datadir, cohortKeys)
    elif searchKeys[0] == 'all' :
        queryPersonStreams(chainName, multichainLoc, datadir, cohortKeys, searchKeys)
    else:
        queryDomainStream(chainName, multichainLoc, datadir, cohortKeys, searchKeys)
       


# ## audit log

# In[35]:


def publishToAuditstream(chainName, multichainLoc, datadir, queryCommand):
    #load wallet
    walletCommand=multichainLoc+'multichain-cli {} -datadir={} listaddresses'.format(chainName, datadir)
    items = subprocess.check_output(walletCommand.split())
    matches = json.loads(items, parse_int= int)
    wallets = pd.DataFrame(matches)
    wallet = wallets['address'][wallets['ismine'] == True][0]
    
    #load time and parse query conducted
    time = str(datetime.utcnow())
    query = ' '.join(queryCommand.split(' ')[3:])
    
    ##publish to auditstream
    streamKeys = (wallet, time, query)
    publishCommand = [multichainLoc+'multichain-cli', 
            str('{}'.format(chainName)), 
            str('-datadir={}'.format(datadir)),
            'publish',
            str('audit_log'), 
            str('["{}", "{}", "{}"]'.format(streamKeys[0], streamKeys[1], streamKeys[2])),
            '{'+'"json":{}'+'}']
    procPublish = subprocess.Popen(publishCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    procPublish.wait()
    return


# In[ ]:


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-cn", "--chainName", help = "the name of the chain to store data", default = "chain1")
    parser.add_argument("-ml", "--multichainLoc", help = "path to multichain commands", default = "")
    parser.add_argument("-dr", "--datadir", help = "path to store the chain")
    parser.add_argument("-ck", "--cohortKeys", help = "concepts to build cohort from e.g. diagnosis code")
    parser.add_argument("-sk", "--searchKeys", type = str, help = "concepts/domain to search for cohort e.g. drugs taken")
    args = parser.parse_args()

    start = time.time()
    try:
        print("--QUERYING--")
        subscribeToStream(args.chainName, args.multichainLoc, args.datadir)
        query(args.chainName, args.multichainLoc, args.datadir, args.cohortKeys, args.searchKeys)
        
        end = time.time()
        
        e = int(end - start)
        print('\n\n Time elapsed:\n\n')
        print( '{:02d}:{:02d}:{:02d}'.format(e // 3600, (e % 3600 // 60), e % 60))

    
    except:
        sys.stderr.write("\nERROR: Failed querying. Please try again.\n")
        quit()
        

if __name__ == "__main__":
    main()

