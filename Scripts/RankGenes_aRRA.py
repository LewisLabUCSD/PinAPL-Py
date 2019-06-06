#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 17:15:19 2019

@author: philipp
"""

# =======================================================================
# Sort genes by adjusted robust rank aggregation (Li et al., Genome Biology 2014)
# =======================================================================
from joblib import Parallel, delayed
import multiprocessing
from scipy.stats import beta
from statsmodels.distributions.empirical_distribution import ECDF
import yaml
import pandas
import numpy


def aRRA(g):
    gene = geneList[g]    
    gene_ranks_sig = list() 
    I = genes_x.index(gene)
    i = I
    terminate = False
    while genes_x[i] == gene and terminate == False:
        if NB_pval_x[i] < P0_aRRA:
            gene_ranks_sig.append(ranks_x[i])
        if i <= L-2:
            i+=1
        else:
            terminate = True
    gene_ranks_sig.sort()
    J = len(gene_ranks_sig)
    if J>0:
        U = [gene_ranks_sig[i]/L for i in range(J)]
        rho_k = list()
        for k in range(J):
            kth_smallest = U[k]
            pval = beta.cdf(kth_smallest,k+1,J-k)
            rho_k.append(pval)
        rho = min(rho_k)
    else:           # no sgRNAs are significant
        rho = 1
    return rho


def aRRA_null(I):
    I_ranks_sig = list()
    for i in I:
        if NB_pval_x[i] < P0_aRRA:
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


def compute_aRRA(HitList):    
    configFile = open('configuration.yaml','r')
    config = yaml.load(configFile)
    configFile.close()          
    Np = config['Np']
    r = config['NumGuidesPerGene']
    alpha_g = config['alpha_g']    
    screentype = config['ScreenType']  
    num_cores = multiprocessing.cpu_count()
    fc = list(HitList['fold change'])
    NB_pval = list(HitList['p-value']) 
    sig = list(HitList['significant'])
    genes = list(HitList['gene'])       
    global L; L = len(genes)
    global geneList; geneList = list(set(genes))
    G = len(geneList)
    # -------------------------------------------        
    # Compute fold change ranks
    # -------------------------------------------
    fc_DF = pandas.DataFrame(data = {'fc': fc},
                                 columns = ['fc'])
    if screentype == 'enrichment':
        fc_DF = fc_DF.rank(ascending = False)
    elif screentype == 'depletion':
        fc_DF = fc_DF.rank(ascending = True)
    global ranks; ranks = list(fc_DF['fc'])    
    # -------------------------------------------        
    # Read out maximal sgRNA p-value that's still significant
    # -------------------------------------------    
    I0 = sig.index(False)
    global P0_aRRA; P0_aRRA = NB_pval[I0-1]
    # -------------------------------------------        
    # Compute aRRA metric
    # -------------------------------------------    
    if set(['N/A']) != set(NB_pval):
        print('Computing aRRA scores ...')
        aRRA_DF = pandas.DataFrame(data = {'gene': [genes[i] for i in range(L)],
                                        'ranks': [ranks[i] for i in range(L)],
                                        'NB_pval': [NB_pval[i] for i in range(L)]},
                                  columns = ['gene','ranks','NB_pval'])
        aRRA_DF = aRRA_DF.sort_values(['gene'])
        global genes_x; genes_x = list(aRRA_DF['gene'])
        global ranks_x; ranks_x = list(aRRA_DF['ranks'])
        global NB_pval_x; NB_pval_x = list(aRRA_DF['NB_pval'])
        metric = Parallel(n_jobs=num_cores)(delayed(aRRA)(g) for g in range(G))        
        # Permutation
        print('Estimating aRRA null distribution ('+str(Np)+' permutations)...')
        I_perm = numpy.random.choice(L,size=(Np,r),replace=True)
        metric_null = Parallel(n_jobs=num_cores)(delayed(aRRA_null)(I) for I in I_perm)
        # p-value
        print('Estimating aRRA p-values ...')
        ecdf = ECDF(metric_null)
        metric_pval = list()
        for g in range(G):    
            pval = ecdf(metric[g])
            metric_pval.append(pval)
        sig_metric = [True if metric_pval[g] < alpha_g else False for g in range(G)]
    else: # no control replicates
        print('### ERROR: Cannot compute aRRA scores without significant sgRNAs! ###')
        metric = [-1 for k in range(G)]
        metric_pval = [-1 for k in range(G)]
        metric_sig = ['N/A' for k in range(G)]
    # Return output
    return metric, metric_pval, sig_metric