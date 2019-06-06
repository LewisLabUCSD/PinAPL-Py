#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 09:23:51 2017

@author: philipp
"""
# Analyze count distribution
# =======================================================================
# Imports
from __future__ import division # floating point division by default
import sys
import yaml
import os
import glob
import pandas
import scipy.stats.mstats as sc
import numpy
import time


def Normalization():
    # ------------------------------------------------
    # Print header
    # ------------------------------------------------
    print('++++++++++++++++++++++++++++++++++++++++++++++++')
    start = time.time() 

    # ------------------------------------------------
    # Get parameters
    # ------------------------------------------------
    configFile = open('configuration.yaml','r')
    config = yaml.load(configFile)
    configFile.close()
    ScriptsDir = config['ScriptsDir']
    sgRNAReadCountDir = config['sgRNAReadCountDir']
    GeneReadCountDir = config['GeneReadCountDir']      
    delta = config['delta']
    norm = config['Normalization']
    RoundCount = config['RoundCount']
    NormSuffix = '_normalized.txt'
    N0 = 1000000
    eps = 0.001 
    
    # ------------------------------------------------
    # Get files
    # ------------------------------------------------  
    os.chdir(sgRNAReadCountDir)
    FileNames_u = glob.glob('*_GuideCounts.txt')
    colnames_u = ['sgRNA','gene','counts']
    os.chdir(GeneReadCountDir)
    FileNames_g = glob.glob('*_GeneCounts.txt')     
    colnames_g = ['gene','counts'] 

    # ------------------------------------------------
    # Normalization to counts per million
    # ------------------------------------------------               
    if norm == 'cpm':    
        print('Normalizing to counts per million reads ...')   
        # sgRNA counts   
        os.chdir(sgRNAReadCountDir) 
        for filename in FileNames_u:
            print('Processing file '+filename+' ...') 
            GuideCounts = pandas.read_table(filename,sep='\t',names=colnames_u)
            L = len(GuideCounts)
            sgIDs = list(GuideCounts['sgRNA'])        
            geneIDs = list(GuideCounts['gene'])                            
            ReadsPerGuide = list(GuideCounts['counts'])
            N = sum(ReadsPerGuide)
            if RoundCount:
                ReadsPerGuide_0 = [int(numpy.round(ReadsPerGuide[k]/N * N0)) for k in range(L)]
            else:
                ReadsPerGuide_0 = [ReadsPerGuide[k]/N * N0 for k in range(L)]
            GuideCounts0_Filename = filename[0:-4] + NormSuffix            
            GuideCounts0 = pandas.DataFrame()
            GuideCounts0['sgID'] = sgIDs
            GuideCounts0['geneID'] = geneIDs
            GuideCounts0['Norm. Read Counts'] = ReadsPerGuide_0
            GuideCounts0.to_csv(GuideCounts0_Filename, sep = '\t', index = False, header = False)
        # gene counts
        os.chdir(GeneReadCountDir) 
        for filename in FileNames_g:
            print('Processing file '+filename+' ...')
            GeneCounts = pandas.read_table(filename,sep='\t',names=colnames_g)
            G = len(GeneCounts)            
            geneIDs = list(GeneCounts['gene'])
            ReadsPerGene = list(GeneCounts['counts'])
            N = sum(ReadsPerGene)
            if RoundCount:
                ReadsPerGene_0 = [int(numpy.round(ReadsPerGene[j]/N * N0)) for j in range(G)]
            else:
                ReadsPerGene_0 = [ReadsPerGene[j]/N * N0 for j in range(G)]
            GeneCounts0_Filename = filename[0:-4] + NormSuffix
            GeneCounts0 = pandas.DataFrame()
            GeneCounts0['geneID'] = geneIDs
            GeneCounts0['Norm. Read Counts'] = ReadsPerGene_0
            GeneCounts0.to_csv(GeneCounts0_Filename, sep = '\t', index = False, header = False)
    # ------------------------------------------------------------
    # Normalization to mean total read count across replicates
    # ------------------------------------------------------------    
    elif norm == 'total':
        print('Normalizing to mean total read count ...')
        os.chdir(sgRNAReadCountDir)       
        TotalCounts = list()
        for filename in FileNames_u:
            SampleFile = pandas.read_table(filename, sep='\t',names=colnames_u)
            x = list(SampleFile['counts'])
            TotalCounts.append(numpy.sum(x))
        MeanCount = numpy.mean(TotalCounts)
        # sgRNA counts 
        os.chdir(sgRNAReadCountDir)
        for filename in FileNames_u:
            print('Processing file '+filename+' ...')             
            GuideCounts = pandas.read_table(filename,sep='\t',names=colnames_u)
            L = len(GuideCounts)
            sgIDs = list(GuideCounts['sgRNA'])        
            geneIDs = list(GuideCounts['gene'])                            
            ReadsPerGuide = list(GuideCounts['counts'])
            N = sum(ReadsPerGuide)
            if RoundCount:            
                ReadsPerGuide_0 = [int(numpy.round(ReadsPerGuide[k]/N * MeanCount)) for k in range(L)]             
            else:
                ReadsPerGuide_0 = [ReadsPerGuide[k]/N * MeanCount for k in range(L)]                       
            GuideCounts0_Filename = filename[0:-4] + NormSuffix
            GuideCounts0 = pandas.DataFrame()
            GuideCounts0['sgID'] = sgIDs
            GuideCounts0['geneID'] = geneIDs
            GuideCounts0['Norm. Read Counts'] = ReadsPerGuide_0
            GuideCounts0.to_csv(GuideCounts0_Filename, sep = '\t', index = False, header = False)
        # gene counts
        os.chdir(GeneReadCountDir)
        for filename in FileNames_g:
            print('Processing file '+filename+' ...')
            GeneCounts = pandas.read_table(filename,sep='\t',names=colnames_g)
            G = len(GeneCounts)
            geneIDs = list(GeneCounts['gene'])                            
            ReadsPerGene = list(GeneCounts['counts'])
            N = sum(ReadsPerGene)
            if RoundCount:
                ReadsPerGene_0 = [int(numpy.round(ReadsPerGene[j]/N * MeanCount)) for j in range(G)]
            else:
                ReadsPerGene_0 = [ReadsPerGene[j]/N * MeanCount for j in range(G)]
            GeneCounts0_Filename = filename[0:-4] + NormSuffix
            GeneCounts0 = pandas.DataFrame()
            GeneCounts0['geneID'] = geneIDs
            GeneCounts0['Norm. Read Counts'] = ReadsPerGene_0
            GeneCounts0.to_csv(GeneCounts0_Filename, sep = '\t', index = False, header = False)
    # ------------------------------------------------------------
    # Normalization by size-factor (Love et al., Genome Biol 2014)
    # ------------------------------------------------------------     
    elif norm == 'size':
        print('Normalizing by size-factors ...')       
        # Establish data frame
        os.chdir(sgRNAReadCountDir)
        filename = FileNames_u[0]
        SampleFile = pandas.read_table(filename, sep='\t',names=colnames_u)
        sgIDs = list(SampleFile['sgRNA'])
        geneIDs = list(SampleFile['gene'])
        L = len(sgIDs)        
        RawCounts = pandas.DataFrame(data = {'sgRNA': [sgIDs[k] for k in range(L)],
                                         'gene': [geneIDs[k] for k in range(L)]},
                                columns = ['sgRNA','gene'])
        SizeFactors = pandas.DataFrame(data = {'sgRNA': [sgIDs[k] for k in range(L)],
                                         'gene': [geneIDs[k] for k in range(L)]},
                                columns = ['sgRNA','gene'])
        # Compute geometric means for all sgRNAs
        print('Computing geometric means ...')
        for filename in FileNames_u:
            sample = filename[0:-16]
            SampleFile = pandas.read_table(filename, sep='\t',names=colnames_u)
            x = list(SampleFile['counts'])
            RawCounts[sample] = x
            SizeFactors[sample] = [x[k] if x[k]>0 else x[k]+eps for k in range(L)]
        geomean = [sc.gmean(list(SizeFactors.iloc[k,2:])) for k in range(L)]
        SizeFactors['Geom mean'] = geomean
        # Compute size-factors for each sgRNA and each sample   
        print('Computing sgRNA size-factors ...')
        for filename in FileNames_u:
            sample = filename[0:-16]
            x = SizeFactors[sample]
            g0 = SizeFactors['Geom mean']
            x0_k = [x[k]/g0[k] for k in range(L)]
            SizeFactors[sample+' sgRNA size-factors'] = [x0_k[k] for k in range(L)]
        # Compute size-factor for each sample
        print('Computing sample size-factors ...')
        for filename in FileNames_u:
            sample = filename[0:-16]
            SizeFactors[sample+' size-factor'] = numpy.median(SizeFactors[sample+' sgRNA size-factors'])
        # Write size-factor dataframe
        SizeFactors.to_csv('Size-factors.txt',sep='\t',index=False)
        # Write normalized counts dataframe
        print('Writing normalized read counts ...')
        # sgRNA counts         
        for filename in FileNames_u:
            sample = filename[0:-16]
            if RoundCount:
                ReadsPerGuide_0 = [int(numpy.round(RawCounts[sample][k]/SizeFactors[sample+' size-factor'][k])) \
                        for k in range(L)]
            else:
                ReadsPerGuide_0 = [RawCounts[sample][k]/SizeFactors[sample+' size-factor'][k] for k in range(L)]
            GuideCounts0_Filename = filename[0:-4] + NormSuffix
            GuideCounts0 = pandas.DataFrame()
            GuideCounts0['sgID'] = sgIDs
            GuideCounts0['geneID'] = geneIDs
            GuideCounts0['Norm. Read Counts'] = ReadsPerGuide_0
            GuideCounts0.to_csv(GuideCounts0_Filename, sep = '\t', index = False, header = False)
        # gene counts        
        os.chdir(GeneReadCountDir)  
        for filename in FileNames_g:
            sample = filename[0:-15]
            GeneCounts = pandas.read_table(filename,sep='\t',names=colnames_g)
            G = len(GeneCounts)
            geneIDs = list(GeneCounts['gene'])
            ReadsPerGene = list(GeneCounts['counts'])                    
            if RoundCount:
                ReadsPerGene_0 = [int(numpy.round(ReadsPerGene[j]/SizeFactors[sample+' size-factor'][j])) \
                    for j in range(G)]
            else:
                ReadsPerGene_0 = [ReadsPerGene[j]/SizeFactors[sample+' size-factor'][j] for j in range(G)] 
            GeneCounts0_Filename = filename[0:-4] + NormSuffix
            GeneCounts0 = pandas.DataFrame()
            GeneCounts0['geneID'] = geneIDs
            GeneCounts0['Norm. Read Counts'] = ReadsPerGene_0
            GeneCounts0.to_csv(GeneCounts0_Filename, sep = '\t', index = False, header = False)                            
    # ------------------------------------------------------------
    # Spelling error catch
    # ------------------------------------------------------------  
    else:
        print('### ERROR: Check spelling of Normalization parameter in configuration file! ###')
    
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
        print('Time elapsed (Total) [secs]: ' + '%.3f' % time_elapsed +'\n')
    elif sec_elapsed < 3600:
        time_elapsed = sec_elapsed/60
        print('Time elapsed (Total) [mins]: ' + '%.3f' % time_elapsed +'\n')
    else:
        time_elapsed = sec_elapsed/3600
        print('Time elapsed (Total) [hours]: ' + '%.3f' % time_elapsed +'\n')
            


if __name__ == "__main__":
    Normalization()
