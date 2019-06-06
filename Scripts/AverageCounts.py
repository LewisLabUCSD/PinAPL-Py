#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 17:49:59 2017

@author: philipp
"""

# Average counts across replicates
# =======================================================================
# Imports
from __future__ import division # floating point division by default
import pandas
import numpy
import os
import glob
import time
import yaml
import sys


def AverageReadCounts(treatment):
    # ------------------------------------------------
    # Get parameters
    # ------------------------------------------------
    configFile = open('configuration.yaml','r')
    config = yaml.load(configFile)
    configFile.close()
    ScriptsDir = config['ScriptsDir']
    sgRNAReadCountDir = config['sgRNAReadCountDir']
    GeneReadCountDir = config['GeneReadCountDir']    
    repl_avg = config['repl_avg']
    
    # ------------------------------------------------
    # Get counts from each replicate 
    # ------------------------------------------------ 
    start = time.time()
    print('Processing '+treatment+' ...')
    colnames_s = ['sgRNA','gene','counts']
    colnames_g = ['gene','counts']    
    
    # sgRNA counts
    os.chdir(sgRNAReadCountDir)
    ReplFiles = [d for d in os.listdir(sgRNAReadCountDir) if treatment in d and '_avg' not in d]
    R = len(ReplFiles)
    if R >= 2:
        print('Averaging read counts over '+str(R)+' replicates ...')
        AllGuideCounts = pandas.DataFrame()  
        for filename in ReplFiles:            
            CountsFile = pandas.read_table(filename, sep='\t',names=colnames_s)
            CountsFile = CountsFile.sort_values(['gene','sgRNA'])
            sgIDs = list(CountsFile['sgRNA'])        
            genes = list(CountsFile['gene'])                
            counts = list(CountsFile['counts'])    
            AllGuideCounts['sgRNA'] = sgIDs 
            AllGuideCounts['gene'] = genes         
            AllGuideCounts[filename] = counts
        # Compute averages  
        repl_counts = AllGuideCounts.iloc[:,2:]
        if repl_avg == 'median':
            avg_counts = repl_counts.median(axis=1)
        elif repl_avg == 'mean':
            avg_counts = repl_counts.mean(axis=1)
        AllGuideCounts[treatment+'_avg'] = avg_counts
        del_columns = range(2,2+R)
        AllGuideCounts.drop(AllGuideCounts.columns[del_columns],axis=1,inplace=True)
        # Write result dataframe
        AllGuideCounts.to_csv(treatment+'_avg_GuideCounts.txt', sep = '\t', index = False, header = False)
    else:
        print('(No filenames found)')               
       
    # gene counts
    os.chdir(GeneReadCountDir)
    ReplFiles = [d for d in os.listdir(GeneReadCountDir) if treatment in d and '_avg' not in d]
    R = len(ReplFiles)
    if R >= 2:
        AllGeneCounts = pandas.DataFrame()
        for filename in ReplFiles:
            CountsFile = pandas.read_table(filename, sep='\t',names=colnames_g)
            CountsFile = CountsFile.sort_values(['gene'])
            genes = list(CountsFile['gene'])                
            counts = list(CountsFile['counts'])    
            AllGeneCounts['gene'] = genes         
            AllGeneCounts[filename] = counts                
        # Compute averages 
        repl_counts = AllGeneCounts.iloc[:,1:]
        if repl_avg == 'median':
            avg_counts = repl_counts.median(axis=1)
        elif repl_avg == 'mean':
            avg_counts = repl_counts.mean(axis=1)
        AllGeneCounts[treatment+'_avg'] = avg_counts
        del_columns = range(1,1+R)
        AllGeneCounts.drop(AllGeneCounts.columns[del_columns],axis=1,inplace=True)
        # Write result dataframes
        os.chdir(GeneReadCountDir)
        AllGeneCounts.to_csv(treatment+'_avg_GeneCounts.txt', sep = '\t', index = False, header = False)


    # --------------------------------------
    # Time stamp
    # --------------------------------------        
    os.chdir(ScriptsDir)
    end = time.time()
    # Final time stamp
    print('------------------------------------------------')
    print('Script completed.')      
    sec_elapsed = end - start
    if sec_elapsed < 60:
        time_elapsed = sec_elapsed
        print('Time elapsed [secs]: ' + '%.3f' % time_elapsed +'\n')
    elif sec_elapsed < 3600:
        time_elapsed = sec_elapsed/60
        print('Time elapsed [mins]: ' + '%.3f' % time_elapsed +'\n')
    else:
        time_elapsed = sec_elapsed/3600
        print('Time elapsed [hours]: ' + '%.3f' % time_elapsed +'\n')    
    
if __name__ == "__main__":
    input1 = sys.argv[1]    
    AverageReadCounts(input1)     