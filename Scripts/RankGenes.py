#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 14:46:37 2016

@author: philipp
"""

# Find canidate genes
# =======================================================================
# Imports
from __future__ import division # floating point division by default
import yaml
import os
import time
import pandas
import glob
import numpy
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
import multiprocessing
from collections import Counter
import sys
from pvalPlots import *
from RankGenes_SigmaFC import *
from RankGenes_AvgLogFC import *
from RankGenes_aRRA import *
from RankGenes_STARS import *


def CountGuidesPerGene(g):
    gene = geneList[g]       
    i0 = genes_bundled.index(gene)
    i = i0
    nG = 0
    terminate = False
    while genes_bundled[i] == gene and terminate == False:
        nG += 1
        if i <= L-2:
            i+=1
        else:
            terminate = True    
    return nG
     


def GeneRankingAnalysis(sample):
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
    sgRNARanksDir = config['sgRNARanksDir']
    EffDir = config['EffDir']
    GeneDir = config['GeneDir']
    alpha_g = config['alpha_g'] 
    padj = config['padj']
    screentype = config['ScreenType']    
    num_cores = multiprocessing.cpu_count()
    GeneMetric = config['GeneMetric']
    SheetFormat = config['HitListFormat']
    pvalDir = config['pvalDir_genes'] 
    res = config['dpi']
    svg = config['svg']
    r = config['NumGuidesPerGene']
 
   
    # ------------------------------------------------
    # Read sgRNA enrichment/depletion table
    # ------------------------------------------------
    os.chdir(sgRNARanksDir)
    print('Loading sgRNA '+screentype+' table ...')
    filename = glob.glob(sample+'_*sgRNAList.txt')[0]
    sgRNARanking = pandas.read_table(filename, sep='\t')
    if screentype == 'enrichment':
        sgRNARanking = sgRNARanking.sort_values(['significant','p-value','fold change','sgRNA'],ascending=[0,1,0,1])
    elif screentype == 'depletion':
        sgRNARanking = sgRNARanking.sort_values(['significant','p-value','fold change','sgRNA'],ascending=[0,1,1,1])        
    global L; L = len(sgRNARanking)
    genes = list(sgRNARanking['gene'])
    global geneList; geneList = list(set(genes))
    fc = list(sgRNARanking['fold change'])
    NB_pval = list(sgRNARanking['p-value'])
    G = len(geneList)
    NB_sig = list(sgRNARanking['significant'])


    # ------------------------------------------------
    # Find number of sgRNAs per gene (some genes have less than r sgRNAs)
    # ------------------------------------------------
    Aux_sgRNARanking = sgRNARanking.sort_values(['gene'])
    global genes_bundled; genes_bundled = list(Aux_sgRNARanking['gene'])
    nGuides = Parallel(n_jobs=num_cores)(delayed(CountGuidesPerGene)(g) for g in range(G))

 
    # ------------------------------------------
    # Find number of significant sgRNAs per gene
    # ------------------------------------------  
    # Find genes with at least 1 signif. sgRNA
    print('Looking for sgRNAs with significant fold change ...')
    Temp_DF = pandas.DataFrame(data = {'gene': genes,
                                       'NB_significant': NB_sig},
                                       columns = ['gene','NB_significant'])        
    Temp_DF = Temp_DF.sort_values(['NB_significant'],ascending=False)
    NB_sig = list(Temp_DF['NB_significant'])
    genes0 = list(Temp_DF['gene'])
    n0 = NB_sig.index(False)
    sigGenes = [genes0[k] for k in range(n0)]       # genes with at least 1 signif. sgRNA
    sigGenesList = list(set(sigGenes))
    # count significant sgRNAs for all genes with at least 1 sign. sgRNA        
    guidesPerGene = list()
    GeneCounts = Counter()
    for gene in sigGenes:
        GeneCounts[gene] += 1
    for gene in sigGenesList:
        guidesPerGene.append(GeneCounts[gene])
    # Save number of significant sgRNAs per gene
    sigGuides = list()
    for gene in geneList:
        if gene not in sigGenesList:
            sigGuides.append(0)
        else:
            sigGuides.append(guidesPerGene[sigGenesList.index(gene)])  
    # Plot histogram
    if not os.path.exists(EffDir):
        os.makedirs(EffDir)
    os.chdir(EffDir)
    fig, ax = plt.subplots(figsize=(3.5,2.9))
    if len(guidesPerGene) > 0:
        plt.hist(guidesPerGene, bins = range(1,r+2), align = 'left',color='#42f4a1',edgecolor='black')
    else:
        plt.figtext(0.5,0.5,'N/A')
    plt.title('sgRNA on-Target Efficacy', fontsize=12)
    plt.xlabel('# Significant sgRNAs', fontsize=11)
    plt.ylabel('Number of Genes', fontsize=11)
    plt.xticks(range(1,r+1))
    plt.tick_params(labelsize=11)
    plt.tight_layout()
    plt.savefig(sample+'_sgRNA_Efficacy.png',dpi=res)   
    if svg:
        plt.savefig(sample+'_sgRNA_Efficacy.svg')                   
    
    
    # ------------------------------------------
    # Rank genes according to specified method
    # ------------------------------------------     
    os.chdir(ScriptsDir)    
    if GeneMetric == 'SigmaFC':
        metric, metric_pval, metric_sig = compute_SigmaFC(sgRNARanking)
        SortFlag = False if screentype=='enrichment' else True   # metric based on fold-change
    elif GeneMetric == 'aRRA':
        metric, metric_pval, metric_sig = compute_aRRA(sgRNARanking)
        SortFlag = True    # metric based on p-val     
    elif GeneMetric == 'STARS':
        metric,metric_pval,metric_sig = compute_STARS(sgRNARanking)        
        SortFlag = False      # metric always from high to low  
    elif GeneMetric == 'AvgLogFC':
        metric, metric_pval, metric_sig = compute_AvgLogFC(sgRNARanking,nGuides)
        SortFlag = False if screentype=='enrichment' else True  # metric based on fold-change
    else:
        print('### ERROR: Cannot find gene ranking method! ###')
    # Correcting cdf artifact in case of no significant sgRNAs
    if False not in metric_sig:
        metric_sig = [False for g in range(G)]
        

    # -------------------------------------------------  
    # Plotting p-value distribution
    # -------------------------------------------------  
    if set(['N/A']) != set(NB_pval):    
        print('Plotting gene metric p-value distribution...')
        PlotTitle = 'Gene '+screentype.capitalize()+' ('+GeneMetric+')'
        pvalHist(metric_pval,pvalDir,sample,res,svg,'#c8d1ca',PlotTitle)           

           
    # -------------------------------------------------  
    # Output list
    # -------------------------------------------------  
    if not os.path.exists(GeneDir):
        os.makedirs(GeneDir)
    os.chdir(GeneDir)
    print('Writing results dataframe ...')
    Results_df = pandas.DataFrame(data = {'gene': [geneList[g] for g in range(G)],
                                    GeneMetric: [metric[g] for g in range(G)],
                                     'p_value': [metric_pval[g] for g in range(G)],
                                     'significant': [str(metric_sig[g]) for g in range(G)],
                                     '# sgRNAs': [nGuides[g] for g in range(G)],                
                                     '# signif. sgRNAs': [sigGuides[g] for g in range(G)]},
                            columns = ['gene',GeneMetric,'p_value',\
                            'significant','# sgRNAs','# signif. sgRNAs'])
    Results_df_0 = Results_df.sort_values(['significant',GeneMetric],ascending=[False,SortFlag])
    GeneListFilename = filename[0:-14]+'_'+GeneMetric+'_GeneList.txt'
    Results_df_0.to_csv(GeneListFilename, sep = '\t', index = False)      
    if SheetFormat == 'xlsx':
        print('Converting to xlsx ...')
        GeneListFilename = filename[0:-14]+'_'+GeneMetric+'_GeneList.xlsx'
        Results_df_0.to_excel(GeneListFilename)              


    # -------------------------------------------------  
    # Time Stamp
    # -------------------------------------------------      
    os.chdir(ScriptsDir)
    end_total = time.time()
    print('------------------------------------------------')
    print('Script completed.')
    sec_elapsed = end_total - start_total    
    if sec_elapsed < 60: 
        time_elapsed = sec_elapsed
        print('Time elapsed [secs]: ' + '%.3f' %time_elapsed +'\n')
    elif sec_elapsed < 3600:
        time_elapsed = sec_elapsed/60
        print('Time elapsed [mins]: ' + '%.3f' % time_elapsed +'\n')
    else:
        time_elapsed = sec_elapsed/3600
        print('Time elapsed [hours]: ' + '%.3f' % time_elapsed +'\n')      
    


if __name__ == "__main__":
    input1 = sys.argv[1]
    GeneRankingAnalysis(input1) 