#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 17:15:19 2019

@author: philipp
"""

# =======================================================================
# Sort genes by average Log FC 
# =======================================================================
import numpy
import yaml
from joblib import Parallel, delayed
import multiprocessing
from statsmodels.distributions.empirical_distribution import ECDF
import pandas


def AverageLogFC(g):
    gene = geneList[g]   
    lfc_list = list()
    i0 = GeneBundle.index(gene)
    i = i0
    terminate = False
    while GeneBundle[i] == gene and terminate == False:
        lfc_list.append(lfc[i])
        if i <= L-2:
            i+=1
        else:
            terminate = True
    if repl_avg == 'mean':
        avglfc = numpy.mean(lfc_list)
    elif repl_avg == 'median':
        avglfc = numpy.median(lfc_list)
    return avglfc


def AvgLogFC_null(I):
    logFC_I = [lfc[i] for i in I]
    if repl_avg == 'mean':
        avglfc_I = numpy.mean(logFC_I)
    elif repl_avg == 'median':
        avglfc_I = numpy.median(logFC_I)     
    return avglfc_I


def compute_AvgLogFC(HitList,nGuides):    
    configFile = open('configuration.yaml','r')
    config = yaml.load(configFile)
    configFile.close()  
    r = config['NumGuidesPerGene']
    screentype = config['ScreenType']
    num_cores = multiprocessing.cpu_count()
    Np = config['Np']
    alpha_g = config['alpha_g']
    global repl_avg; repl_avg = config['repl_avg']
    # -------------------------------------------------  
    # Compute average log fold change across sgRNAs
    # ------------------------------------------------- 
    print('Computing average fold change across sgRNAs ...')
    genes = list(HitList['gene'])    
    fc = list(HitList['fold change'])
    global L; L = len(genes)
    global geneList; geneList = list(set(genes))
    G = len(geneList)
    Aux_DF = pandas.DataFrame(data = {'gene': [genes[i] for i in range(L)],
                                    'lfc': [numpy.log2(fc[i]) for i in range(L)]},
                              columns = ['gene','lfc'])
    Aux_DF = Aux_DF.sort_values(['gene'])
    global GeneBundle; GeneBundle = list(Aux_DF['gene'])      # genes bundled 
    global lfc; lfc = list(Aux_DF['lfc'])
    AvgLogFC = Parallel(n_jobs=num_cores)(delayed(AverageLogFC)(g) for g in range(G))
    # ------------------------------------------------- 
    # Compute permutations
    # -------------------------------------------------     
    I_perm = numpy.random.choice(L,size=(Np,r),replace=True)
    metric_null = Parallel(n_jobs=num_cores)(delayed(AvgLogFC_null)(I) for I in I_perm)
    ecdf = ECDF(metric_null)
    metric_pval = list()
    for g in range(G):
        if nGuides[g] == 1:
            pval = 'N/A'        # exclude p-value for miRNAs etc 
            metric_pval.append(pval)            
        elif screentype == 'enrichment':
            pval = 1 - ecdf(AvgLogFC[g])
            metric_pval.append(pval)   
        elif screentype == 'depletion':
            pval = ecdf(AvgLogFC[g])
            metric_pval.append(pval)            
        else:
            print('### ERROR: Check spelling of ScreenType in configuration file! ###')
    metric_sig = [True if metric_pval[g] < alpha_g else False for g in range(G)]
    return AvgLogFC, metric_pval, metric_sig