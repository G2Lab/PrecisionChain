#!/usr/bin/env python
# coding: utf-8

# In[10]:


import sys
import time
import argparse
import subprocess
from subprocess import Popen, PIPE
import os
import pandas as pd
import warnings
import multiprocessing
import random
import numpy as np
import itertools
from itertools import permutations
warnings.simplefilter(action='ignore', category=FutureWarning)


# In[30]:


##load what files still need to be done
path = '/gpfs/commons/groups/gursoy_lab/aelhussein/Data/VCF_ALS/'
completed = pd.read_csv(path + 'completed.txt', sep = ' ', names = ['files', 'chrom'])
files = pd.read_csv(path + 'files.txt', names = ['files'])
all_files = list(itertools.product(list(files.files.values), [i for i in range(1,23)]))
all_files= pd.DataFrame(all_files, columns = ['files', 'chrom'])
left = pd.concat([all_files,completed],axis=0)
left = df.drop_duplicates(keep=False)


# In[ ]:


def ViewVCFConvert(left):

    for i, row in left.iterrows():
        file = row['files']
        chrom = row['chrom']
        request = '''if [[ -f "chr{}.vcf.gz" ]];
                then
                    bcftools view -Oz {} --regions chr{} > add_{}.vcf.gz 
                    tabix -p vcf add_{}.vcf.gz
                    bcftools merge -Oz chr{}.vcf.gz add_{}.vcf.gz > merged_chr{}.vcf.gz
                    mv merged_chr{}.vcf.gz chr{}.vcf.gz
                    tabix -p vcf chr{}.vcf.gz;
                    echo "{} {}" >> completed.txt;
                else
                    bcftools view -Oz {} --regions chr{} > chr{}.vcf.gz
                    tabix -p vcf chr{}.vcf.gz
                    echo "{} {}" >> completed.txt;
                fi;'''.format(chrom,
                              file, chrom, chrom, chrom, chrom, chrom, chrom, chrom, chrom, chrom, file, chrom,
                              file, chrom, chrom, chrom, file, chrom)
        output = subprocess.run(request, shell = True)


def main():
    args = [i for i in range(1,23)]

    number_of_cores = 22

    with Pool(number_of_cores) as pool:

        results = pool.map(ViewVCFConvert, args)

if __name__ == "__main__":
    main()

