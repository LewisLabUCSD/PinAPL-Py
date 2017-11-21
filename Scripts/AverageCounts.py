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
    AlnQCDir = config['AlnQCDir']
    repl_avg = config['repl_avg']
    AvgDir = treatment+'_avg'
    
    # ------------------------------------------------
    # Get counts from each replicate 
    # ------------------------------------------------ 
    start = time.time()
    print('Processing '+treatment+' ...')
    colnames_s = ['sgRNA','gene','counts']
    colnames_g = ['gene','counts']    
    os.chdir(AlnQCDir)
    ReplDirs = [d for d in os.listdir(AlnQCDir) if treatment in d and '_avg' not in d]
    R = len(ReplDirs)
    if R >= 2:
        print('Averaging read counts over '+str(R)+' replicates ...')
        AllGuideCounts = pandas.DataFrame()  
        AllGuideCounts0 = pandas.DataFrame()
        AllGeneCounts = pandas.DataFrame()
        AllGeneCounts0 = pandas.DataFrame()    
        for replicate in ReplDirs:
            os.chdir(replicate)
            # sgRNA counts
            filename = glob.glob('*GuideCounts.txt')[0]
            CountsFile = pandas.read_table(filename, sep='\t',names=colnames_s)
            CountsFile = CountsFile.sort_values(['gene','sgRNA'])
            sgIDs = list(CountsFile['sgRNA'])        
            genes = list(CountsFile['gene'])                
            counts = list(CountsFile['counts'])    
            AllGuideCounts['sgRNA'] = sgIDs 
            AllGuideCounts['gene'] = genes         
            AllGuideCounts[replicate] = counts
            # normalized sgRNA counts
            filename = glob.glob('*GuideCounts_0.txt')[0]
            CountsFile = pandas.read_table(filename, sep='\t',names=colnames_s)
            CountsFile = CountsFile.sort_values(['gene','sgRNA'])        
            sgIDs = list(CountsFile['sgRNA'])        
            genes = list(CountsFile['gene'])                
            counts = list(CountsFile['counts'])    
            AllGuideCounts0['sgRNA'] = sgIDs 
            AllGuideCounts0['gene'] = genes         
            AllGuideCounts0[replicate] = counts        
            # gene counts
            filename = glob.glob('*GeneCounts.txt')[0]
            CountsFile = pandas.read_table(filename, sep='\t',names=colnames_g)
            CountsFile = CountsFile.sort_values(['gene'])
            genes = list(CountsFile['gene'])                
            counts = list(CountsFile['counts'])    
            AllGeneCounts['gene'] = genes         
            AllGeneCounts[replicate] = counts        
            # normalized gene counts
            filename = glob.glob('*GeneCounts_0.txt')[0]
            CountsFile = pandas.read_table(filename, sep='\t',names=colnames_g)
            CountsFile = CountsFile.sort_values(['gene'])
            genes = list(CountsFile['gene'])                
            counts = list(CountsFile['counts'])    
            AllGeneCounts0['gene'] = genes     
            AllGeneCounts0[replicate] = counts
            os.chdir(AlnQCDir)            
        # ------------------------------------------------
        # Compute averages
        # ------------------------------------------------ 
        # sgRNA counts    
        repl_counts = AllGuideCounts.iloc[:,2:]
        if repl_avg == 'median':
            avg_counts = repl_counts.median(axis=1)
        elif repl_avg == 'mean':
            avg_counts = repl_counts.mean(axis=1)
        AllGuideCounts[treatment+'_avg'] = avg_counts
        del_columns = range(2,2+R)
        AllGuideCounts.drop(AllGuideCounts.columns[del_columns],axis=1,inplace=True) 
        # normalized sgRNA counts    
        repl_counts = AllGuideCounts0.iloc[:,2:]
        if repl_avg == 'median':
            avg_counts = repl_counts.median(axis=1)
        elif repl_avg == 'mean':
            avg_counts = repl_counts.mean(axis=1)
        AllGuideCounts0[treatment+'_avg'] = avg_counts
        del_columns = range(2,2+R)
        AllGuideCounts0.drop(AllGuideCounts0.columns[del_columns],axis=1,inplace=True)
        # gene counts    
        repl_counts = AllGeneCounts.iloc[:,1:]
        if repl_avg == 'median':
            avg_counts = repl_counts.median(axis=1)
        elif repl_avg == 'mean':
            avg_counts = repl_counts.mean(axis=1)
        AllGeneCounts[treatment+'_avg'] = avg_counts
        del_columns = range(1,1+R)
        AllGeneCounts.drop(AllGeneCounts.columns[del_columns],axis=1,inplace=True)
        # normalized gene counts    
        repl_counts = AllGeneCounts0.iloc[:,1:]
        if repl_avg == 'median':
            avg_counts = repl_counts.median(axis=1)
        elif repl_avg == 'mean':
            avg_counts = repl_counts.mean(axis=1)
        AllGeneCounts0[treatment+'_avg'] = avg_counts
        del_columns = range(1,1+R)
        AllGeneCounts0.drop(AllGeneCounts0.columns[del_columns],axis=1,inplace=True)    
        # ------------------------------------------------
        # Write result dataframes
        # ------------------------------------------------     
        os.chdir(AlnQCDir)
        if not os.path.exists(AvgDir):
            os.makedirs(AvgDir)
        os.chdir(AvgDir)
        AllGuideCounts.to_csv(treatment+'_avg_GuideCounts.txt', sep = '\t', index = False, header = False)
        AllGuideCounts0.to_csv(treatment+'_avg_GuideCounts_0.txt', sep = '\t', index = False, header = False)
        AllGeneCounts.to_csv(treatment+'_avg_GeneCounts.txt', sep = '\t', index = False, header = False)
        AllGeneCounts0.to_csv(treatment+'_avg_GeneCounts_0.txt', sep = '\t', index = False, header = False)
    else:
        print('(No replicates found)')


    # --------------------------------------
    # Time stamp
    # --------------------------------------        
    os.chdir(ScriptsDir)
    end = time.time()
    # Final time stamp
    sec_elapsed = end - start
    if sec_elapsed < 60:
        time_elapsed = sec_elapsed
        print('Time elapsed [secs]: ' + '%.3f' % time_elapsed)
    elif sec_elapsed < 3600:
        time_elapsed = sec_elapsed/60
        print('Time elapsed [mins]: ' + '%.3f' % time_elapsed)
    else:
        time_elapsed = sec_elapsed/3600
        print('Time elapsed [hours]: ' + '%.3f' % time_elapsed)    
    
if __name__ == "__main__":
    input1 = sys.argv[1]    
    AverageReadCounts(input1)     