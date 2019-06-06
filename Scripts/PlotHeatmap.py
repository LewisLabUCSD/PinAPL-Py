#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 24 17:43:58 2016

@author: philipp
"""
# Clustering of top N enriched sgRNAs
# ======================================================================================
import os
import glob
import pandas as pd
import yaml
import sys
import time
import subprocess
import re
import numpy

def TopN_Clustering():
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
    WorkingDir = config['WorkingDir']
    AnalysisDir = config['AnalysisDir']
    sgRNAReadCountDir = config['sgRNAReadCountDir']
    ClusterDir = config['HeatDir']
    ScriptsDir = config['ScriptsDir']
    ClusterBy = config['ClusterBy']
    N = config['TopN']
    delta = config['delta']
    width = config['width_p']
    height = config['height_p']
    fontsize = config['fontsize_p']
    marginsize = config['marginsize']
    ScreenType = config['ScreenType']
    svg = config['svg']
    
    # ------------------------------------------------
    # Dataframe containing all counts 
    # ------------------------------------------------
    print('Loading read counts ...')
    os.chdir(sgRNAReadCountDir)
    Filenames = [d for d in os.listdir(sgRNAReadCountDir) if 'normalized' in d \
        and '_avg' not in d]
    AllCounts = pd.DataFrame()    
    colnames = ['sgRNA','gene','counts']
    for filename in Filenames:
        sample = filename[0:-27]    
        CountsFile = pd.read_table(filename, sep='\t',names=colnames)
        sgIDs = list(CountsFile['sgRNA'])
        genes = list(CountsFile['gene'])                
        counts = list(CountsFile['counts'])    
        AllCounts['sgRNA'] = sgIDs 
        AllCounts['gene'] = genes         
        AllCounts[sample] = counts

    # ------------------------------------------------
    # Sort read counts
    # ------------------------------------------------    
    if ClusterBy == 'variance':
        # Compute variance
        print('Computing variance ...')
        L = len(sgIDs)
        AC = AllCounts.iloc[:,2:].as_matrix()
        Var = [numpy.var(AC[k,]) for k in range(L)]
        # Insert into dataframe
        AllCounts['Variance'] = Var
        # Sort by highest variance
        print('Extracting top '+str(N)+' most variable sgRNAs ...')
        Q = AllCounts.sort_values('Variance',ascending=False)    
        Q = Q[0:(N-1)]
        # Write dataframe
        print('Writing dataframe ...')
        OutputSheetname = 'Top'+str(N)+'_Var.txt'
        if not os.path.exists(ClusterDir):
           os.makedirs(ClusterDir)
        os.chdir(ClusterDir)
        Q.to_csv(OutputSheetname,sep='\t',index=False)     
        HeatmapFilename = 'Heatmap_Top'+str(N)+'_most_variable_sgRNAs'
    elif ClusterBy == 'counts':
        TopGuides = set()
        # Assembling top N guides  
        Samples = [filename[0:-27] for filename in Filenames]
        for sample in Samples:
            print('Reading top '+str(N)+' sgRNAs in sample '+sample+' ...')     
            SampleCounts = AllCounts[['sgRNA',sample]]
            if ScreenType == 'enrichment':
                SampleCounts = SampleCounts.sort_values(sample,ascending=0)
            elif ScreenType == 'depletion':
                SampleCounts = SampleCounts.sort_values(sample,ascending=1)
            SampleIDs = list(SampleCounts['sgRNA'])
            TopGuides_sample = SampleIDs[0:int(N-1)]
            TopGuides = set.union(TopGuides,TopGuides_sample)
        TopGuides = list(TopGuides)    
        T = len(TopGuides)      
        # Establish data frame
        print('Assembling data frame ...')
        TopIndex = [sgIDs.index(TopGuides[k]) for k in range(T)]
        TopGenes = [genes[TopIndex[k]] for k in range(T)]
        TopGuides_df = pd.DataFrame(data = {'sgRNA': [TopGuides[k] for k in range(T)],
                                            'gene': [TopGenes[k] for k in range(T)]},
                                    columns = ['sgRNA','gene'])  
        # Extract counts for top guides and write into data frame
        for sample in Samples:   
            SampleCounts = AllCounts[['sgRNA',sample]]            
            if ScreenType == 'enrichment':
                SampleCounts = SampleCounts.sort_values([sample],ascending=[0])
            elif ScreenType == 'depletion':
                SampleCounts = SampleCounts.sort_values(sample,ascending=[1])           
            SampleIDs = list(SampleCounts['sgRNA'])
            TopIndex_sample = [SampleIDs.index(TopGuides[k]) for k in range(T)]
            counts = list(SampleCounts[sample])
            TopCounts = [counts[TopIndex_sample[k]] for k in range(T)]
            TopGuides_df[sample] = TopCounts
        # Write dataframe  
        TopGuides_df[' '] = ' '
        OutputSheetname = 'Top'+str(N)+'_Counts.txt'
        if not os.path.exists(ClusterDir):
            os.makedirs(ClusterDir)  
        os.chdir(ClusterDir)
        TopGuides_df.to_csv(OutputSheetname,sep='\t',index=False)
        HeatmapFilename = 'Heatmap_combined_Top'+str(N)+'_highest_sgRNA_counts'
   

    # ------------------------------------------------
    # Calling R for clustering analysis
    # ------------------------------------------------
    print('Plotting heatmap ...')
    print('--------------------------------------------')    
    # Define command and arguments
    command = 'Rscript'
    path2script = ScriptsDir+'ClusterAnalysis_Heatmap.r'
    args = [ClusterDir,OutputSheetname,str(delta),HeatmapFilename,\
            str(width),str(height),str(fontsize),str(marginsize),str(svg)]  
    cmd = [command, path2script] + args    
    # Run R script
    subprocess.check_output(cmd)
    print('--------------------------------------------')    
    
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
    TopN_Clustering() 