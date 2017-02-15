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
import pandas as pd
import glob
import numpy
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
import multiprocessing
from scipy.stats import beta
from scipy.stats import norm
from statsmodels.distributions.empirical_distribution import ECDF
from statsmodels.sandbox.stats.multicomp import multipletests
from collections import Counter
from decimal import *
import sys


def computeES(g):
    GOI = geneList[g]
    GOI_ranks = list()
    noGOI_ranks = list()
    for k in range(L):
        gene = genes[k]
        if GOI == gene:
            GOI_ranks.append(ranks[k])
            noGOI_ranks.append(0)
        else:
            GOI_ranks.append(0)
            noGOI_ranks.append(1/(G-r))
    RankSum = sum(GOI_ranks)
    GOI_cumsum = numpy.cumsum(GOI_ranks).tolist()
    noGOI_cumsum = numpy.cumsum(noGOI_ranks).tolist()
    # compute GSEA score
    p_hit = list()
    p_miss = list()
    p = list()
    for i in range(L):
        p_hit.append(GOI_cumsum[i]/RankSum)
        p_miss.append(noGOI_cumsum[i])
        p.append(p_hit[i] - p_miss[i])
    ES = max(p)
    return ES

def computeESnull(I):
    I_ranks = list()
    noI_ranks = list()
    for k in range(L):
        if k in I:
            I_ranks.append(ranks[k])
            noI_ranks.append(0)
        else:
            I_ranks.append(0)
            noI_ranks.append(1/(G-r))
    RankSum = sum(I_ranks)
    I_cumsum = numpy.cumsum(I_ranks).tolist()
    noI_cumsum = numpy.cumsum(noI_ranks).tolist()
    # compute enrichment score
    p_hit = list()
    p_miss = list()
    p = list()
    for i in range(L):
        p_hit.append(I_cumsum[i]/RankSum)
        p_miss.append(noI_cumsum[i])
        p.append(p_hit[i] - p_miss[i])
    ES_null = max(p)
    return ES_null  

def compute_aRRA(g):
    GOI = geneList[g]    
    GOI_ranks_sig = list()    
    for i in range(L):
        if GOI == genes[i] and (NB_pval[i] < alpha):
            GOI_ranks_sig.append(ranks[i])
    j = len(GOI_ranks_sig)
    if j>0:
        U = [GOI_ranks_sig[i]/L for i in range(j)]
        rho_k = list()
        for k in range(j):
            kth_smallest = U[j-1-k]
            pval = 1 - beta.cdf(kth_smallest,k+1,j-k)
            rho_k.append(pval)
        rho = min(rho_k)
    else:
        rho = 1
    return rho

def compute_aRRA_null(I):
    I_ranks_sig = list()
    for i in range(L):
        if (i in I) and (NB_pval[i] < alpha):
            I_ranks_sig.append(ranks[i])
    j = len(I_ranks_sig)
    if j>0:
        U = [I_ranks_sig[i]/L for i in range(j)]
        rho_k = list()
        for k in range(j):
            kth_smallest = U[j-1-k]
            pval = 1 - beta.cdf(kth_smallest,k+1,j-k)
            rho_k.append(pval)
        rho_null = min(rho_k)
    else:
        rho_null = 1
    return rho_null    

def TimeStamp(sec_elapsed,ProcessName):
    if sec_elapsed < 60: 
        time_elapsed = sec_elapsed
        print('Time elapsed ('+ProcessName+') [secs]: ' + '%.3f' %time_elapsed)
    elif sec_elapsed < 3600:
        time_elapsed = sec_elapsed/60
        print('Time elapsed ('+ProcessName+') [mins]: ' + '%.3f' % time_elapsed)
    else:
        time_elapsed = sec_elapsed/3600
        print('Time elapsed ('+ProcessName+') [hours]: ' + '%.3f' % time_elapsed)             




def PrepareGeneList(sample):
    # ------------------------------------------------
    # Print header
    # ------------------------------------------------
    print('++++++++++++++++++++++++++++++++')
    print('PinAPL-Py: Gene Ranking Analysis')
    print('++++++++++++++++++++++++++++++++')
    start_total = time.time()    

    # ------------------------------------------------
    # Get parameters
    # ------------------------------------------------
    configFile = open('configuration.yaml','r')
    config = yaml.load(configFile)
    configFile.close()    
    AnalysisDir = config['AnalysisDir']
    ListDir = config['HitDir']
    EffDir = config['EffDir']
    GeneDir = config['GeneDir']
    global alpha
    alpha = config['alpha']            
    pcorr = config['pcorr']
    num_cores = multiprocessing.cpu_count()
    global r
    r = config['sgRNAsPerGene']
    GeneMetric = config['GeneMetric']
    SheetFormat = config['HitListFormat']
    
    # --------------------------        
    # Read data from sgRNA list
    # --------------------------   
    os.chdir(ListDir)
    print('Loading sgRNA counts ...')        
    if SheetFormat == 'tsv':
        filename = glob.glob(sample+'_*sgRNAList.tsv')[0]
        HitList = pd.read_table(filename, sep='\t')
    elif SheetFormat == 'xlsx':
        filename = glob.glob(sample+'_*sgRNAList.xlsx')[0]
        HitList = pd.read_excel(filename)    
    sgIDs = list(HitList['sgRNA'].values)
    global genes
    genes = list(HitList['gene'].values)
    counts = list(HitList['counts [norm.]'].values)
    global NB_pval
    NB_pval = list(HitList['NB_pval'].values)
    sig = list(HitList['significant'].values)
    global L        
    L = len(sgIDs)
    global geneList
    geneList = list(set(genes))
    global G
    G = len(geneList)
    
    # ------------------------------------------
    # Find number of significant sgRNAs per gene
    # ------------------------------------------
    if min(NB_pval) > -1:    
        print('Looking for sgRNAs with significant fold change ...')
        n0 = sig.index(False)
        sigGenes = [genes[k] for k in range(n0)]
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
        plt.figure()
        plt.hist(guidesPerGene, bins = range(1,r+2), align = 'left',color="g")
        plt.suptitle(sample+': on-Target Efficiency', fontsize=14, fontweight='bold')
        plt.xlabel('# sign. sgRNAs', fontsize=14)
        plt.ylabel('Number of Genes', fontsize=14)
        plt.savefig(sample+'_sgRNA_Efficiency.png')
    else: # no control replicates
        print('WARNING: No control replicates found! No significant sgRNAs counted.')
        sigGuides = ['N/A' for k in range(G)]
        # Plot histogram
        if not os.path.exists(EffDir):
            os.makedirs(EffDir)
        os.chdir(EffDir)
        plt.figure()
        plt.figtext(0.5,0.5,'N/A')
        plt.suptitle(sample+': on-Target Efficiency', fontsize=14, fontweight='bold')
        plt.xlabel('# sign. sgRNAs', fontsize=14)
        plt.ylabel('Number of Genes', fontsize=14)
        plt.savefig(sample+'_sgRNA_Efficiency.png')        
    
    # -------------------------------------------        
    # Compute sgRNA count ranks
    # -------------------------------------------
    countsDF = pd.DataFrame(data = {'counts': [counts[k] for k in range(L)]},
                                columns = ['counts'])
    countsDF_Ranked = countsDF.rank(ascending=True)
    global ranks
    ranks = countsDF_Ranked['counts'].tolist()    
    
    if GeneMetric == 'ES':
        # -------------------------------------------------        
        # compute ES 
        # -------------------------------------------------
        start = time.time()
        SortFlag = True        
        print('Computing gene enrichment scores ...')
        metric = Parallel(n_jobs=num_cores)(delayed(computeES)(g) for g in range(G))
        end = time.time()
        sec_elapsed = end - start
        TimeStamp(sec_elapsed,'ES Computation')
        # Permutation
        start = time.time()
        Np = config['Np_ES']
        print('Estimating ES null distribution ('+str(Np)+' Permutations) ...')
        I_perm = numpy.random.choice(L,size=(Np,r),replace=False)          
        metric_null = Parallel(n_jobs=num_cores)(delayed(computeESnull)(I) for I in I_perm)
        end = time.time()
        sec_elapsed = end-start
        TimeStamp(sec_elapsed,'ES Permutation')
        # p-value
        print('Computing ES p-values ...')        
        ecdf = ECDF(metric_null)
        pval_list= list()
        for g in range(G):
            GOI = geneList[g]
            pval = 1 - ecdf(metric[g])
            pval_list.append(pval)
        # Determine critical p value (FDR correction)
        pval_corr = multipletests(pval_list,alpha,pcorr)
        metric_sig = pval_corr[0]
    elif GeneMetric == 'aRRA':            
        # -------------------------------------------------        
        # compute aRRA 
        # -------------------------------------------------
        if min(NB_pval) > -1:
            start = time.time()
            SortFlag = True
            print('Computing a-RRA scores ...')
            metric = Parallel(n_jobs=num_cores)(delayed(compute_aRRA)(g) for g in range(G))        
            end = time.time()
            sec_elapsed = end - start
            TimeStamp(sec_elapsed,'a-RRA Computation')        
            # Permutation
            start = time.time()
            Np = config['Np_aRRA']
            print('Estimating a-RRA null distribution ('+str(Np)+' Permutations)...')
            I_perm = numpy.random.choice(L,size=(Np,r),replace=False)
            metric_null = Parallel(n_jobs=num_cores)(delayed(compute_aRRA_null)(I) for I in I_perm)
            end = time.time()
            sec_elapsed = end - start
            TimeStamp(sec_elapsed,'a-RRA Permutation')
            # p-value
            print('Computing a-RRA p-values ...')
            ecdf = ECDF(metric_null)
            pval_list = list()
            for g in range(G):    
                GOI = geneList[g]
                pval = ecdf(metric[g])
                pval_list.append(pval)
            # Determine critical p value (FDR correction)
            pval_corr = multipletests(pval_list,alpha,pcorr)
            metric_sig = pval_corr[0]
        else: # no control replicates
            print('ERROR: Cannot compute a-RRA scores without control replicates!')
            SortFlag = True
            Np = config['Np_aRRA']
            metric = [-1 for k in range(G)]
            pval_list = [-1 for k in range(G)]
            metric_sig = ['N/A' for k in range(G)]
    elif GeneMetric == 'STARS':
        # -------------------------------------------------        
        # compute STARS score 
        # -------------------------------------------------            
        start = time.time()
        STARSDir = config['STARSDir']
        Np = config['Np_STARS']
        thr = config['thr_STARS']
        screentype = config['ScreenType']
        os.chdir(ListDir)
        SortFlag = False
        # Running null distribution
        print('Estimating STARS null distribution ('+str(Np)+' Permutations)...')
        STARS_input = open('STARS_input.txt','w')
        STARS_input.write('sgID\tcounts\n')
        STARS_chip = open('STARS_chip.txt','w')
        STARS_chip.write('sgID\tgene\n')
        for k in range(L):
            STARS_input.write(sgIDs[k]+'\t'+str(counts[k])+'\n')
            STARS_chip.write(sgIDs[k]+'\t'+genes[k]+'\n')                
        STARS_input.close()    
        STARS_chip.close()
        os.system('mv STARS_input.txt STARS_chip.txt '+STARSDir)
        os.chdir(STARSDir)
        STARS_null_cmd = 'python stars_null_v1.2.py --input-file \
            STARS_input.txt --chip-file STARS_chip.txt --thr '+str(thr)+' --num-ite '+str(Np)
        os.system(STARS_null_cmd)
        # Computing STARS score
        print('Computing STARS scores ...')
        if screentype == 'enrichment':
            d = 'P'
        elif screentype == 'depletion':
            d = 'N'
        STARS_cmd = 'python stars_v1.2.py --input-file STARS_input.txt \
            --chip-file STARS_chip.txt --thr '+str(thr)+' --dir '+d+' \
            --null Null_STARSOutput8_'+str(thr)+'.txt'
        os.system(STARS_cmd)
        # Extracting statistics
        STARS_output = glob.glob('counts_STARSOutput*.txt')[0]
        STARS = pd.read_table(STARS_output, sep='\t')
        geneList_s = list(STARS['Gene Symbol'].values)
        G = len(geneList_s)
        s_index = [geneList.index(geneList_s[k]) for k in range(G)]
        geneList = geneList_s
        sigGuides_s = [sigGuides[s_index[k]] for k in range(G)]
        sigGuides = sigGuides_s
        metric = list(STARS['STARS Score'].values)
        pval_list = list(STARS['p-value'].values)
        pval_corr = multipletests(pval_list,alpha,pcorr)
        metric_sig = pval_corr[0]                     
        # Deleting STARS files            
        os.system('rm STARS_input.txt STARS_chip.txt')
        os.system('rm Null_STARSOutput8_'+str(thr)+'.txt')
        os.system('rm '+STARS_output)
            
            
    # -------------------------------------------------  
    # Output list
    # -------------------------------------------------  
    if not os.path.exists(GeneDir):
        os.makedirs(GeneDir)
    os.chdir(GeneDir)
    print('Writing results dataframe ...')
    Results_df = pd.DataFrame(data = {'gene': [geneList[g] for g in range(G)],
                                    GeneMetric: [metric[g] for g in range(G)],
                                     GeneMetric+' p_value': ['%.2E' % Decimal(pval_list[g]) for g in range(G)],
                                     'significant': [metric_sig[g] for g in range(G)],                                                    
                                     'signif. sgRNAs': [sigGuides[g] for g in range(G)]},
                            columns = ['gene',GeneMetric,GeneMetric+' p_value','significant',\
                            'signif. sgRNAs'])
    Results_df_0 = Results_df.sort_values(['significant',GeneMetric],ascending=[False,SortFlag])
    if SheetFormat == 'tsv':
        GeneListFilename = filename[0:-14]+'_'+GeneMetric+'_'+'P'+str(Np)+'_GeneList.tsv'
        Results_df_0.to_csv(GeneListFilename, sep = '\t', index = False)      
    elif SheetFormat == 'xlsx':
        GeneListFilename = filename[0:-14]+'_'+GeneMetric+'_'+'P'+str(Np)+'_GeneList.xlsx'
        Results_df_0.to_excel(GeneListFilename)              
    
    end_total = time.time()
    print('------------------------------------------------')
    print('Script completed.')
    sec_elapsed = end_total - start_total    
    TimeStamp(sec_elapsed,'Total')
    print('\n')

if __name__ == "__main__":
    input1 = sys.argv[1]
    PrepareGeneList(input1) 
