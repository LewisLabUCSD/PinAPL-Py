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
from scipy.stats import beta
from scipy.stats import norm
from statsmodels.distributions.empirical_distribution import ECDF
from statsmodels.sandbox.stats.multicomp import multipletests
from collections import Counter
from decimal import *
import sys
from pvalPlots import *


def compute_aRRAx(g):
    GOI = geneList[g]    
    GOI_ranks_sig = list() 
    I = genes_x.index(GOI)
    i = I
    terminate = False
    while genes_x[i] == GOI and terminate == False:
        if NB_pval_x[i] < P_0:
            GOI_ranks_sig.append(ranks_x[i])
        if i <= L-2:
            i+=1
        else:
            terminate = True
    GOI_ranks_sig.sort()
    J = len(GOI_ranks_sig)
    if J>0:
        U = [GOI_ranks_sig[i]/L for i in range(J)]
        rho_k = list()
        for k in range(J):
            kth_smallest = U[k]
            pval = beta.cdf(kth_smallest,k+1,J-k)
            rho_k.append(pval)
        rho = min(rho_k)
    else:
        rho = 1
    return rho

def compute_aRRA_nullx(I):
    I_ranks_sig = list()
    for i in I:
        if NB_pval[i] < P_0:
            I_ranks_sig.append(ranks[i])
    I_ranks_sig.sort()
    J = len(I_ranks_sig)
    if J>0:
        U = [I_ranks_sig[i]/L for i in range(J)]
        rho_k = list()
        for k in range(J):
            kth_smallest = U[k]
            pval = beta.cdf(kth_smallest,k+1,J-k)
            rho_k.append(pval)
        rho_null = min(rho_k)
    else:
        rho_null = 1
    return rho_null 

    
def AverageLogFC(g):
    gene = geneList[g]        
    geneIndex = [i for i,x in enumerate(genes) if x==pl-pygene]
    nG = len(geneIndex)
    logfcs = [numpy.log2(fc[i]) for i in geneIndex]
    avglfc = numpy.mean(logfcs)
    return nG,avglfc

def AverageLogFCx(g):
    gene = geneList[g]       
    lfc_list = list()
    I = genes_X.index(gene)
    i = I
    terminate = False
    while genes_X[i] == gene and terminate == False:
        lfc_list.append(lfc_X[i])
        if i <= L-2:
            i+=1
        else:
            terminate = True    
    nG = len(lfc_list)
    avglfc = numpy.mean(lfc_list)
    return nG,avglfc


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


def GeneRankingAnalysis(sample):
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
    ScriptsDir = config['ScriptsDir']
    AnalysisDir = config['AnalysisDir']
    ListDir = config['HitDir']
    EffDir = config['EffDir']
    GeneDir = config['GeneDir']
    alpha = config['alpha']            
    padj = config['padj']
    screentype = config['ScreenType']    
    num_cores = multiprocessing.cpu_count()
    GeneMetric = config['GeneMetric']
    Np = config['Np']
    SheetFormat = config['HitListFormat']
    pvalDir = config['pvalDir']
    res = config['dpi']
    svg = config['svg']
    r = config['NumGuidesPerGene']
    global P_0; P_0 = config['P_0']
    
    # ------------------------------------------------
    # Read sgRNA enrichment/depletion table
    # ------------------------------------------------
    os.chdir(ListDir)
    print('Loading sgRNA '+screentype+' table ...')    
    filename = glob.glob(sample+'_*sgRNAList.tsv')[0]
    HitList = pandas.read_table(filename, sep='\t')
    if screentype == 'enrichment':
        HitList = HitList.sort_values(['significant','p-value','fold change','sgRNA'],ascending=[0,1,0,1])
    elif screentype == 'depletion':
        HitList = HitList.sort_values(['significant','p-value','fold change','sgRNA'],ascending=[0,1,1,1])        
    sgIDs = list(HitList['sgRNA'])
    global genes; genes = list(HitList['gene'])
    global geneList; geneList = list(set(genes))
    fc = list(HitList['fold change'])
    global NB_pval; NB_pval = list(HitList['p-value'])
    global L; L = len(sgIDs)
    global G; G = len(geneList)
    NB_sig = list(HitList['significant'].values)
    
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
    sigGenes = [genes0[k] for k in range(n0)]
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
    plt.figure(figsize=(4,3.5))
    if len(guidesPerGene) > 0:
        plt.hist(guidesPerGene, bins = range(1,r+2), align = 'left',color='#42f4a1',edgecolor='black')
    else:
        plt.figtext(0.5,0.5,'N/A')
    plt.title('on-Target Efficacy', fontsize=14)
    plt.xlabel('# Significant sgRNAs', fontsize=12)
    plt.ylabel('Number of Genes', fontsize=12)
    plt.xticks(range(1,r+1))
    plt.tick_params(labelsize=12)
    plt.tight_layout()
    plt.savefig(sample+'_sgRNA_Efficacy.png',dpi=res)   
    if svg:
        plt.savefig(sample+'_sgRNA_Efficacy.svg')                   
    
    # -------------------------------------------        
    # Rank data sheet after fold change
    # -------------------------------------------
    fc_DF = pandas.DataFrame(data = {'fc': list(HitList['fold change'])},
                                 columns = ['fc'])
    if screentype == 'enrichment':
        fc_DF = fc_DF.rank(ascending = False)
    elif screentype == 'depletion':
        fc_DF = fc_DF.rank(ascending = True)
    global ranks; ranks = list(fc_DF['fc'])    
        
    # -------------------------------------------        
    # Carry out gene ranking analysis
    # -------------------------------------------    
    if GeneMetric == 'aRRA':            
        # -------------------------------------------------        
        # compute aRRA 
        # -------------------------------------------------
        if min(NB_pval) < 1:
            start = time.time()
            SortFlag = True
            print('Computing aRRA scores ...')
            aRRA_DF = pandas.DataFrame(data = {'gene': [genes[i] for i in range(L)],
                                            'ranks': [ranks[i] for i in range(L)],
                                            'NB_pval': [NB_pval[i] for i in range(L)]},
                                      columns = ['gene','ranks','NB_pval'])
            aRRA_DF = aRRA_DF.sort_values(['gene'])
            global genes_x; genes_x = list(aRRA_DF['gene'])
            global ranks_x; ranks_x = list(aRRA_DF['ranks'])
            global NB_pval_x; NB_pval_x = list(aRRA_DF['NB_pval'])
            metric = Parallel(n_jobs=num_cores)(delayed(compute_aRRAx)(g) for g in range(G))        
            end = time.time()
            sec_elapsed = end - start
            TimeStamp(sec_elapsed,'aRRA Computation')        
            # Permutation
            start = time.time()
            print('Estimating aRRA null distribution ('+str(Np)+' Permutations)...')
            I_perm = numpy.random.choice(L,size=(Np,r),replace=False)
            metric_null = Parallel(n_jobs=num_cores)(delayed(compute_aRRA_nullx)(I) for I in I_perm)
            end = time.time()
            sec_elapsed = end - start
            TimeStamp(sec_elapsed,'aRRA Permutation')
            # p-value
            print('Computing aRRA p-values ...')
            ecdf = ECDF(metric_null)
            metric_pval = list()
            for g in range(G):    
                GOI = geneList[g]
                pval = ecdf(metric[g])
                metric_pval.append(pval)
            # Determine critical p value (FDR correction)
            multTest = multipletests(metric_pval,alpha,padj)
            metric_sig = multTest[0]
            metric_pval0 = multTest[1]
        else: # no control replicates
            print('ERROR: Cannot compute aRRA scores without significant sgRNAs!')
            SortFlag = True
            metric = [-1 for k in range(G)]
            metric_pval = [-1 for k in range(G)]
            metric_pval0 = [-1 for k in range(G)]            
            metric_sig = ['N/A' for k in range(G)]
    elif GeneMetric == 'STARS':
        # -------------------------------------------------        
        # compute STARS score 
        # -------------------------------------------------            
        start = time.time()
        STARSDir = config['STARSDir']
        thr = config['thr_STARS']
        os.chdir(ListDir)
        SortFlag = False
        # Running null distribution
        print('Estimating STARS null distribution ('+str(Np)+' Permutations)...')
        STARS_input = open('STARS_input.txt','w')
        STARS_input.write('sgID\tcounts\n')
        STARS_chip = open('STARS_chip.txt','w')
        STARS_chip.write('sgID\tgene\n')
        for k in range(L):
            STARS_input.write(sgIDs[k]+'\t'+str(fc[k])+'\n')
            STARS_chip.write(sgIDs[k]+'\t'+genes[k]+'\n')                
        STARS_input.close()    
        STARS_chip.close()
        os.system('mv STARS_input.txt STARS_chip.txt '+STARSDir)
        os.chdir(STARSDir)
        STARS_null_cmd = 'python -u stars_null_v1.2.py --input-file \
            STARS_input.txt --chip-file STARS_chip.txt --thr '+str(thr)+' --num-ite '+str(Np)
        os.system(STARS_null_cmd)
        # Computing STARS score
        print('Computing STARS scores ...')
        if screentype == 'enrichment':
            d = 'P'
        elif screentype == 'depletion':
            d = 'N'
        STARS_cmd = 'python -u stars_v1.2.py --input-file STARS_input.txt \
            --chip-file STARS_chip.txt --thr '+str(thr)+' --dir '+d+' \
            --null Null_STARSOutput8_'+str(thr)+'.txt'
        os.system(STARS_cmd)
        # Extracting statistics
        STARS_output = glob.glob('counts_STARSOutput*.txt')[0]
        STARS = pandas.read_table(STARS_output, sep='\t')
        geneList_s = list(STARS['Gene Symbol'].values)
        G = len(geneList_s)
        s_index = [geneList.index(geneList_s[k]) for k in range(G)]
        geneList = geneList_s
        sigGuides_s = [sigGuides[s_index[k]] for k in range(G)]
        sigGuides = sigGuides_s
        metric = list(STARS['STARS Score'].values)
        metric_pval = list(STARS['p-value'].values)
        multTest = multipletests(metric_pval,alpha,padj)
        metric_sig = multTest[0]  
        metric_pval0 = multTest[1]                   
        # Deleting STARS files            
        os.system('rm STARS_input.txt STARS_chip.txt')
        os.system('rm Null_STARSOutput8_'+str(thr)+'.txt')
        os.system('rm '+STARS_output)


    # -------------------------------------------------  
    # p-value plots
    # -------------------------------------------------  
    if min(NB_pval) < 1:    
        print('Plotting p-values ...')
        pvalHist_metric(metric_pval,metric_pval0,GeneMetric,pvalDir,sample,res,svg)           

    # -------------------------------------------------  
    # Compute average fold change across sgRNAs
    # -------------------------------------------------  
    print('Computing average fold change across sgRNAs ...')
    avglfc_DF = pandas.DataFrame(data = {'gene': [genes[i] for i in range(L)],
                                    'lfc': [numpy.log2(fc[i]) for i in range(L)]},
                              columns = ['gene','lfc'])
    avglfc_DF = avglfc_DF.sort_values(['gene'])
    global genes_X; genes_X = list(avglfc_DF['gene'])
    global lfc_X; lfc_X = list(avglfc_DF['lfc'])
    parjob = Parallel(n_jobs=num_cores)(delayed(AverageLogFCx)(g) for g in range(G))      
    nGuides = [parjob[g][0] for g in range(G)]
    AvgLogFCs = [parjob[g][1] for g in range(G)]
           
    # -------------------------------------------------  
    # Output list
    # -------------------------------------------------  
    if not os.path.exists(GeneDir):
        os.makedirs(GeneDir)
    os.chdir(GeneDir)
    print('Writing results dataframe ...')
    Results_df = pandas.DataFrame(data = {'gene': [geneList[g] for g in range(G)],
                                    GeneMetric: [metric[g] for g in range(G)],
                                     GeneMetric+' p_value': [metric_pval[g] for g in range(G)],
                                    GeneMetric+' FDR': [metric_pval0[g] for g in range(G)],
                                     'significant': [str(metric_sig[g]) for g in range(G)],
                                     '# sgRNAs': [nGuides[g] for g in range(G)],                
                                     '# signif. sgRNAs': [sigGuides[g] for g in range(G)],
                                    'avg. logFC': [AvgLogFCs[g] for g in range(G)]},
                            columns = ['gene',GeneMetric,GeneMetric+' p_value',GeneMetric+' FDR',\
                            'significant','# sgRNAs','# signif. sgRNAs','avg. logFC'])
    Results_df_0 = Results_df.sort_values(['significant',GeneMetric],ascending=[False,SortFlag])
    GeneListFilename = filename[0:-14]+'_'+GeneMetric+'_'+'P'+str(Np)+'_GeneList.tsv'
    Results_df_0.to_csv(GeneListFilename, sep = '\t', index = False)      
    if SheetFormat == 'xlsx':
        print('Converting to xlsx ...')
        GeneListFilename = filename[0:-14]+'_'+GeneMetric+'_'+'P'+str(Np)+'_GeneList.xlsx'
        Results_df_0.to_excel(GeneListFilename)              
    
    end_total = time.time()
    print('------------------------------------------------')
    print('Script completed.')
    sec_elapsed = end_total - start_total    
    TimeStamp(sec_elapsed,'Total')
    print('\n')
    os.chdir(ScriptsDir)

if __name__ == "__main__":
    input1 = sys.argv[1]
    GeneRankingAnalysis(input1) 
