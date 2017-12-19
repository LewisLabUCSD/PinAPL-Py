#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 10:22:56 2016

@author: philipp
"""
# Build library index
# =======================================================================
# Imports
import pandas as pd
import yaml
import os
import time
import sys
from Bowtie2 import BuildIndex

def BuildBowtieIndex():

    # ------------------------------------------------
    # Print header
    # ------------------------------------------------
    print('++++++++++++++++++++++++++++++++++++++++++++++++')  
    start_total = time.time()  

    # ------------------------------------------------
    # Get parameters
    # ------------------------------------------------
    configFile = open('configuration.yaml','r')
    config = yaml.load(configFile)
    configFile.close()
    ScriptsDir = config['ScriptsDir']
    bw2Dir = config['bw2Dir']    
    LibDir = config['LibDir']
    IndexDir = config['IndexDir']
    LibFilename = config['LibFilename']
    LibFormat = LibFilename[-3:]
    if LibFormat == 'tsv':
        libsep = '\t'
    elif LibFormat == 'csv':
        libsep = ','

    # ----------------------------------
    # Convert library to fasta
    # ----------------------------------
    print('Converting library to fasta format ...') 
    os.chdir(LibDir)
    LibCols = ['gene','ID','seq']
    LibFile = pd.read_table(LibFilename, sep = libsep, skiprows = 1, names = LibCols)
    seq = list(LibFile['seq'])
    IDs = list(LibFile['ID'])    
    with open('library.fasta','w') as library_fasta:
        for k in range(len(IDs)):
            library_fasta.write('>'+str(IDs[k])+'\n')
            library_fasta.write(str(seq[k])+'\n')
    library_fasta.close()
    
    # ----------------------------------
    # Call bowtie2 index builder
    # ----------------------------------   
    print('Building library index ...')     
    if not os.path.exists(IndexDir):
        os.makedirs(IndexDir)
    os.rename('library.fasta',IndexDir+'library.fasta')
    BuildIndex('library.fasta',IndexDir,bw2Dir)

    # --------------------------------------
    # Final time stamp
    # --------------------------------------        
    os.chdir(ScriptsDir)    
    end_total = time.time()
    # Final time stamp
    print('------------------------------------------------')
    print('Script completed.')    
    sec_elapsed = end_total - start_total
    if sec_elapsed < 60:
        time_elapsed = sec_elapsed
        print('Time elapsed (Total) [secs]: ' + '%.3f' % time_elapsed +'\n')
    elif sec_elapsed < 3600:
        time_elapsed = sec_elapsed/60
        print('Time elapsed (Total) [mins]: ' + '%.3f' % time_elapsed +'\n')
    else:
        time_elapsed = sec_elapsed/3600
        print('Time elapsed (Total) [hours]: ' + '%.3f' % time_elapsed +'\n')            


   
if __name__ == "__main__":
    BuildBowtieIndex()
