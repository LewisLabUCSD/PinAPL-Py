#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 15:54:46 2018

@author: philipp
"""

# Apply (optional) read count cutoff and obtain counts per gene
# =======================================================================
# Imports
import yaml
import os
import sys
import pandas
from joblib import Parallel, delayed
import multiprocessing

def CountReadsPerGene(g):       # faster than slicing
    gene = GeneList[g]        
    I = geneIDs.index(gene)
    j = I
    g_counts = 0
    terminate = False
    while geneIDs[j] == gene and terminate == False:
        g_counts = g_counts + sgRNA_counts[j]
        if j <= L-2:
            j+=1  
        else:
            terminate = True    
    return g_counts  


def ReadCountCutoff(sample):
    # ------------------------------------------------
    # Get parameters
    # ------------------------------------------------
    configFile = open('configuration.yaml','r')
    config = yaml.load(configFile)
    configFile.close()    
    sgRNAReadCountDir = config['sgRNAReadCountDir']
    GeneReadCountDir = config['GeneReadCountDir']
    LibDir = config['LibDir']    
    LibFilename = config['LibFilename']
    LibFormat = LibFilename[-3:]
    if LibFormat == 'tsv':
        libsep = '\t'
    elif LibFormat == 'csv':
        libsep = ','   
    N0 = 1000000
    minN = config['Cutoff']
    GuideCount_Suffix = '_GuideCounts.txt'
    GeneCount_Suffix = '_GeneCounts.txt'    

    # ------------------------------------------------
    # Read library
    # ------------------------------------------------  
    os.chdir(LibDir)
    LibCols = ['gene','ID','seq']
    LibFile = pandas.read_table(LibFilename, sep = libsep, skiprows = 1, names = LibCols)
    LibFile = LibFile.sort_values(['gene','ID'])    
    sgIDs = list(LibFile['ID'])
    global L; L = len(sgIDs)
    global geneIDs; geneIDs = list(LibFile['gene'])
    G = len(set(geneIDs)) 

        
    # Load sgRNA counts
    os.chdir(sgRNAReadCountDir)
    colnames = ['ID','gene','counts']
    FileName = sample+'_GuideCounts.txt'
    ReadCountTable = pandas.read_table(FileName, sep='\t', names=colnames)
    ReadsPerGuide = list(ReadCountTable['counts'])    
    NReads = sum(ReadsPerGuide)
    
    # Apply read count cut-off
    if minN > 0:    
        print('Applying sgRNA minimal count cutoff...')
    ReadSel = [ReadsPerGuide[k] >= NReads/N0*minN for k in range(L)]
    global sgRNA_counts
    sgRNA_counts = [ReadSel[k]*ReadsPerGuide[k] for k in range(L)]    
    
    # Read counts per gene in library       
    global GeneList; GeneList = list(set(geneIDs))
    num_cores = multiprocessing.cpu_count()
    gene_counts = Parallel(n_jobs=num_cores)(delayed(CountReadsPerGene)(g) for g in range(G))  
    
    # Write output  
    # sgRNAs
    GuideCountsFilename = sample + GuideCount_Suffix
    GuideCounts = pandas.DataFrame()
    GuideCounts['sgIDs'] = sgIDs
    GuideCounts['geneIDs'] = geneIDs
    GuideCounts['sgRNA_counts'] = sgRNA_counts
    GuideCounts.to_csv(GuideCountsFilename, sep = '\t', index = False, header = False)
    # genes
    if not os.path.exists(GeneReadCountDir):
        os.makedirs(GeneReadCountDir)
    os.chdir(GeneReadCountDir)
    GeneCountsFilename = sample + GeneCount_Suffix
    GeneCounts = pandas.DataFrame()
    GeneCounts['gene'] = GeneList
    GeneCounts['gene_counts'] = gene_counts
    GeneCounts = GeneCounts.sort_values('gene')
    GeneCounts.to_csv(GeneCountsFilename, sep = '\t', index = False, header = False)
    
    
if __name__ == "__main__":
    input1 = sys.argv[1]
    ReadCountCutoff(input1)