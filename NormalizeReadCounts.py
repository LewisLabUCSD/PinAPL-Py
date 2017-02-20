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
    print('++++++++++++++++++++++++++++++++++++')
    print('PinAPL-Py: Read Count Normalization')
    print('++++++++++++++++++++++++++++++++++++')  
    start = time.time() 

    # ------------------------------------------------
    # Get parameters
    # ------------------------------------------------
    configFile = open('../configuration.yaml','r')
    config = yaml.load(configFile)
    configFile.close()
    ScriptsDir = config['ScriptsDir']
    AlnQCDir = config['AlnQCDir']     
    DepthDir = config['DepthDir']
    delta = config['delta_d']
    norm = config['Normalization']
    NormSuffix = '_0.tsv'
    N0 = 1000000
    
    # ------------------------------------------------
    # Normalization
    # ------------------------------------------------        
    os.chdir(AlnQCDir)
    SampleNames = [d for d in os.listdir(AlnQCDir) if os.path.isdir(d)]
    colnames_u = ['sgRNA','gene','counts']
    colnames_g = ['gene','counts']    
    if norm == 'cpm':
        print('Normalizing to counts per million reads ...')   
        for sample in SampleNames:
            print('Processing '+sample+' ...') 
            os.chdir(sample)
            # sgRNA counts            
            GuideCountsFilename = glob.glob('*GuideCounts.tsv')[0]
            GuideCounts = pandas.read_table(GuideCountsFilename,sep='\t',names=colnames_u)
            sgIDs = list(GuideCounts['sgRNA'].values)        
            geneIDs = list(GuideCounts['gene'].values)                            
            ReadsPerGuide = list(GuideCounts['counts'].values)
            N = sum(ReadsPerGuide)
            GuideCounts0_Filename = GuideCountsFilename[0:-4] + NormSuffix
            GuideCounts0 = open(GuideCounts0_Filename,'w')
            ReadsPerGuide_0 = list()
            for k in range(len(sgIDs)):      
                ReadsPerGuide_0 = int(numpy.ceil(ReadsPerGuide[k]/N * N0))
                GuideCounts0.write(str(sgIDs[k]) + '\t' + str(geneIDs[k]) + '\t' + \
                    str(ReadsPerGuide_0) + '\n')
            GuideCounts0.close()
            # gene counts
            GeneCountsFilename = glob.glob('*GeneCounts.tsv')[0]
            GeneCounts = pandas.read_table(GeneCountsFilename,sep='\t',names=colnames_g)            
            geneIDs = list(GeneCounts['gene'].values)                            
            ReadsPerGene = list(GeneCounts['counts'].values)    
            N = sum(ReadsPerGene)                    
            GeneCounts0_Filename = GeneCountsFilename[0:-4] + NormSuffix
            GeneCounts0 = open(GeneCounts0_Filename,'w')
            ReadsPerGene_0 = list()
            for j in range(len(geneIDs)):    
                ReadsPerGene_0 = int(numpy.ceil(ReadsPerGene[j]/N * N0))
                GeneCounts0.write(str(geneIDs[j]) + '\t' + str(ReadsPerGene_0) + '\n')
            GeneCounts0.close()            
            os.chdir(AlnQCDir)
    elif norm == 'size':
        print('Normalizing by size-factors ...')       
        # Establish data frame
        os.chdir(SampleNames[0])
        filename = glob.glob('*GuideCounts.tsv')[0]
        SampleFile = pandas.read_table(filename, sep='\t',names=colnames_u)
        sgIDs = list(SampleFile['sgRNA'].values)
        genes = list(SampleFile['gene'].values)
        L = len(sgIDs)
        Counts = pandas.DataFrame(data = {'sgRNA': [sgIDs[k] for k in range(L)],
                                         'gene': [genes[k] for k in range(L)]},
                                columns = ['sgRNA','gene'])
        # Compute geometric means for all sgRNAs
        print('Computing geometric means ...')
        os.chdir(AlnQCDir)
        for sample in SampleNames:
            os.chdir(sample)
            filename = glob.glob('*GuideCounts.tsv')[0]
            SampleFile = pandas.read_table(filename, sep='\t',names=colnames_u)
            x = list(SampleFile['counts'].values)
            Counts[sample] = [x[k]+delta for k in range(L)]
            os.chdir(AlnQCDir)
        geomean = [sc.gmean(list(Counts.iloc[k,2:])) for k in range(L)]
        Counts['Geom mean'] = geomean
        # Compute size-factors for each sgRNA and each sample   
        print('Computing sgRNA size-factors ...')
        for sample in SampleNames:
            x = Counts[sample]
            g0 = Counts['Geom mean']
            x0_k = [x[k]/g0[k] for k in range(L)]
            Counts[sample+' sgRNA size-factors'] = [x0_k[k] for k in range(L)]
        # Compute size-factor for each sample
        print('Computing sample size-factors ...')
        for sample in SampleNames:
            Counts[sample+' size-factor'] = numpy.median(Counts[sample+' sgRNA size-factors'])
        # Write size-factor dataframe
        os.chdir(DepthDir)
        Counts.to_csv('Size-factors.tsv',sep='\t',index=False)
        # Write normalized counts dataframe
        print('Writing normalized read counts ...')
        os.chdir(AlnQCDir)        
        Counts0 = pandas.DataFrame(data = {'sgRNA': [sgIDs[k] for k in range(L)],
                                         'gene': [genes[k] for k in range(L)]},
                                columns = ['sgRNA','gene'])
        G = len(list(set(genes)))                                
        GeneCounts0 = pandas.DataFrame(data = {'gene': [genes[j] for j in range(G)]},
                                columns = ['gene'])
        for sample in SampleNames:
            os.chdir(sample)
            Counts0[sample] = [int(numpy.ceil(Counts[sample][k]/Counts[sample+' size-factor'][k])) \
                for k in range(L)]            
            Counts0.to_csv(sample+'_GuideCounts'+NormSuffix,sep='\t',columns=['sgRNA','gene',sample],\
                header=False,index=False)            
            GeneCounts = pandas.read_table(glob.glob('*GeneCounts.tsv')[0],sep='\t',names=colnames_g) 
            ReadsPerGene = list(GeneCounts['counts'].values)
            GeneCounts0[sample] = [int(numpy.ceil(ReadsPerGene[j]/Counts[sample+' size-factor'][0])) \
                for j in range(G)]
            GeneCounts0.to_csv(sample+'_GeneCounts'+NormSuffix,sep='\t',columns=['gene',sample],\
                header=False,index=False)                        
            os.chdir(AlnQCDir)
    else:
        print('ERROR: Check spelling of Normalization parameter in configuration file!')
    
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