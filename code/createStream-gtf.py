#!/usr/bin/env python
# coding: utf-8

# In[1]:


'''
createStream-gtf.py
Creates streams for the gtf files
Usage: $ python createStream-gtf.py -cn=<chain name> -dr=<Chain path> 
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
import pandas as pd


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


#Name the streams that will be used to store the GTF data - there will be one stream per chromosome
def createStreams(chainName, multichainLoc, datadir):

    #data stream  split by chrom
    for i in range(1, 23):
        makeStream(chainName, "gene_chrom_{}".format(i), multichainLoc, datadir)  
        makeStream(chainName, "gene_variant_chrom_{}".format(i), multichainLoc, datadir) 
    return


# In[ ]:


######################################################################################################
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-cn", "--chainName", help = "the name of the chain to store data", default = "chain1")
    parser.add_argument("-ml", "--multichainLoc", help = "path to multichain commands", default = "")
    parser.add_argument("-dr", "--datadir", help = "path to store the chain")
    args = parser.parse_args()

    start = time.time()
    try:
        #make the streams
        print("--STREAM CREATION--")
        createStreams(args.chainName, args.multichainLoc, args.datadir)

        print("Stream construction complete!")
        end = time.time()
        
        e = int(end - start)
        print('\n\n Time elapsed:\n\n')
        print( '{:02d}:{:02d}:{:02d}'.format(e // 3600, (e % 3600 // 60), e % 60))

        process = psutil.Process(os.getpid())
        print('\n\n Total memory in bytes:\n\n')
        print(process.memory_info().rss)
    
    except:
        sys.stderr.write("\nERROR: Failed stream creation. Please try again.\n")
        quit()
        


# In[ ]:


if __name__ == "__main__":
    main()

