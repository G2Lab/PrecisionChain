#!/usr/bin/env python
# coding: utf-8

# In[1]:


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
warnings.simplefilter(action='ignore', category=FutureWarning)


# In[ ]:


def queryMappingStream(chainName, multichainLoc, datadir, cohortDomain, key ):
    queryCommand=multichainLoc+'multichain-cli {} -datadir={} liststreamkeyitems mappingData_{} {}'.format(chainName, datadir, cohortDomain, key)
    items = subprocess.check_output(queryCommand.split())
    matches = json.loads(items, parse_int= int)
    return matches[0]


# In[ ]:


def extractDataStream(chainName, multichainLoc, datadir, cohortDomain, key ):
    matches = queryMappingStream(chainName, multichainLoc, datadir, cohortDomain, key)
    data = matches['data']['json']
    concept_stream = list(data.keys())
    concept_bucket = list(data.values())[0]
    return concept_stream, concept_bucket


# In[ ]:


def extractPersonIDs(chainName, multichainLoc, datadir, cohortDomain, key):
    concept_stream, concept_bucket = extractDataStream(chainName, multichainLoc, datadir, cohortDomain, key)
    matches = []
    for bucket in concept_bucket:
        queryCommand=multichainLoc+'multichain-cli {} -datadir={} liststreamqueryitems {}_id_{}_bucket_{} {{"keys":["{}"]}}'.format(chainName, datadir,
                                                                                                               cohortDomain, concept_stream[0],
                                                                                                               bucket+1, key)
        items = subprocess.check_output(queryCommand.split())
        matches.extend(json.loads(items, parse_int= int))
    
    ids = []
    for match in matches:
        ids.append(int(match['keys'][0]))
    return(list(set(ids)))


# In[ ]:


def queryDomainStream(chainName, multichainLoc, datadir, cohortDomain, key, searchDomain, streamsPath):
    person_ids = extractPersonIDs(chainName, multichainLoc, datadir, cohortDomain, key)

    with open(streamsPath + 'streams.txt') as f:
        streams = json.loads(f.read())
    
    for stream in streams[searchDomain]:
        for person_id in person_ids:
            try:
                queryCommand = multichainLoc+'multichain-cli {} -datadir={} liststreamkeyitems condition_id_{}_bucket_1 {}'.format(chainName, datadir,
                                                                                            stream, person_id)
                items = subprocess.check_output(queryCommand.split())
                matches = json.loads(items, parse_int= int)
            except:
                pass


# In[ ]:


def extractPersonStreams(chainName, multichainLoc, datadir, cohortDomain, key):
    person_ids = extractPersonIDs(chainName, multichainLoc, datadir, cohortDomain, key)
    person_streams = {}
    for person_id in person_ids:
        queryCommand=multichainLoc+'multichain-cli {} -datadir={} liststreamkeyitems mappingData_person {}'.format(chainName, datadir, person_id)
        items = subprocess.check_output(queryCommand.split())
        matches = json.loads(items, parse_int= int) 
        person_streams[person_id] = matches[0]['data']['json']
    return person_streams


# In[ ]:


def queryPersonStreams(chainName, multichainLoc, datadir, cohortDomain, key, searchDomain):
    person_streams = extractPersonStreams(chainName, multichainLoc, datadir, cohortDomain, key)
    for person_id in person_streams.keys():
        queryCommand = multichainLoc+'multichain-cli {} -datadir={} liststreamkeyitems patient_stream_{} {}'.format(chainName, datadir,
                                                                                    person_streams[person_id], person_id)
        items = subprocess.check_output(queryCommand.split())
        matches = json.loads(items, parse_int= int)
        print(matches)
    return matches


# In[ ]:


def query(chainName, multichainLoc, datadir, cohortDomain, key, searchDomain, streamsPath, view):
    if view == 'Domain':
        queryDomainStream(chainName, multichainLoc, datadir, cohortDomain, key, searchDomain, streamsPath)
    else:
        queryPersonStreams(chainName, multichainLoc, datadir, cohortDomain, key, searchDomain)


# In[ ]:


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-cn", "--chainName", help = "the name of the chain to store data", default = "chain1")
    parser.add_argument("-ml", "--multichainLoc", help = "path to multichain commands", default = "")
    parser.add_argument("-dr", "--datadir", help = "path to store the chain")
    parser.add_argument("-cd", "--cohortdomain", help = "cohort domain")
    parser.add_argument("-ks", "--keys", type = str, help = "keys for the stream")
    parser.add_argument("-sd", "--searchdomain", help = "search domain")
    parser.add_argument("-sp", "--streamspath", help = "streams path")
    parser.add_argument("-vw", "--view", help = "view")
    args = parser.parse_args()

    start = time.time()
    try:
        print("--QUERYING--")
        query(args.chainName, args.multichainLoc, args.datadir, args.cohortdomain, args.keys, args.searchdomain, args.streamspath, args.view)
        
        end = time.time()
        
        e = int(end - start)
        print('\n\n Time elapsed:\n\n')
        print( '{:02d}:{:02d}:{:02d}'.format(e // 3600, (e % 3600 // 60), e % 60))

    
    except:
        sys.stderr.write("\nERROR: Failed querying. Please try again.\n")
        quit()
        

if __name__ == "__main__":
    main()

