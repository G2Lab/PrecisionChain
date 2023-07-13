#!/usr/bin/env python
# coding: utf-8

# In[2]:

"""
preprocessVCF.py
Preprocesses files using plink. Filters for missingness, Sample call rate and HWE
Usage: $ python preprocessVCF.py <dataPath>
modified by AE 02/2023
"""


import argparse
import subprocess
from subprocess import Popen, PIPE
import os
import multiprocessing
import pandas as pd
import warnings
import multiprocessing
warnings.simplefilter(action='ignore', category=FutureWarning)


def preprocessVCF(input):
    vcf_file_path, file, chrom = input
    commands = [
                f'plink --vcf {vcf_file_path}/original/{file} --geno 0.05 --make-bed --out {vcf_file_path}/intermediate1_{chrom}',
                f'plink --bfile {vcf_file_path}/intermediate1_{chrom} --hwe 1e-6 --make-bed --out {vcf_file_path}/intermediate2_{chrom}',
                f'rm {vcf_file_path}/intermediate1_{chrom}.bed {vcf_file_path}/intermediate1_{chrom}.bim {vcf_file_path}/intermediate1_{chrom}.fam',
                f'plink --bfile {vcf_file_path}/intermediate2_{chrom} --mind 0.05 --recode vcf --out {vcf_file_path}/final_{file}',
                f'rm {vcf_file_path}/intermediate2_{chrom}.bed {vcf_file_path}/intermediate2_{chrom}.bim {vcf_file_path}/intermediate2_{chrom}.fam',
                f'bgzip {vcf_file_path}/final_{file}.vcf',
                f'mv {vcf_file_path}/final_{file}.vcf.gz {vcf_file_path}/{file}',
                f'tabix -p vcf {vcf_file_path}/{file}',
                f'rm {vcf_file_path}/*log'
            ]
    for cmd in commands:
        subprocess.run(cmd, shell=True, check=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-dp", "--dataPath", type = str, help = "path to vcf files for preprocessing")
    args = parser.parse_args()

    cpu = multiprocessing.cpu_count()
    vcf_file_path = args.dataPath
    files = [file for file in  os.listdir(f'{vcf_file_path}/original') if file.endswith('.vcf.gz')]
    chroms = [chrom.split('.')[0][3:] for chrom in files]
    inputs = [(vcf_file_path, file, chrom) for file, chrom in zip(files, chroms)]
    
    pool =  multiprocessing.Pool(cpu)
    pool.map(preprocessVCF, inputs)
    pool.close()
    pool.join()

if __name__ == '__main__':
     main()
