#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 17:15:19 2019

@author: philipp
"""

# =======================================================================
# Sort genes by accumulated fold-change
# =======================================================================
import yaml
import pandas
import numpy as np
import os
import multiprocessing
from joblib import Parallel, delayed
from statsmodels.distributions.empirical_distribution import ECDF



def SigmaFC_Permutation(P):  
    P = list(P)
    P_data = temp_df.iloc[P]
    sg_table = P_data[P_data['sgRNA_signif']==True]
    if len(sg_table)>0:     
        k = len(sg_table)         # number of signif. sgRNAs
        Nx = list(sg_table['counts']); 
        N0 = list(sg_table['control mean']);               
        pi = 0        
        for j in range(k):
            # compute fold-change weight   (discontinued)
#            y = np.log10((Nx[j]+delta)/(N0[j]+delta)) # log fold-change
#            if y < FCmin:
#                u = np.exp(Lambda*(y-FCmin))
#            else:
            u = 1
            pi_j = u * np.log10((Nx[j]+delta)/(N0[j]+delta))
            pi = pi + pi_j    
        #Pi = w[k]*pi       # only w = id supported
        Pi = k*pi   # Avoid "out of range" error in case of an unexpected number of sgRNAs per gene
    else:
        Pi = 0
    return Pi


def compute_SigmaFC(sgRNAList):
    # ------------------------------------------------
    # Get parameters
    # ------------------------------------------------    
    configFile = open('configuration.yaml','r')
    config = yaml.load(configFile)
    configFile.close()
    ScriptsDir = config['ScriptsDir']
    GeneDir = config['GeneDir']
    global delta; delta = config['delta']
    r = config['NumGuidesPerGene']
    Np = config['Np']
    padj = config['padj']
    alpha_g = config['alpha_g']
    ScreenType = config['ScreenType']
    num_cores = multiprocessing.cpu_count()    
#    global Lambda; Lambda = -np.log(0.1)   # deprecation rate
#    global FCmin;           # deprecation cutoff (log10 fold change). DISCONTINUED
#    FCmin = config['FCmin_SigmaFC']
    global w;                # reward function for numbers of signif. sgRNAs
    w = config['w_SigmaFC']


    # ------------------------------------------------
    # Reading counts from sample and control
    # ------------------------------------------------    
    print('Reading sgRNA counts ...')       
    L = len(sgRNAList)
    genes = list(sgRNAList['gene']) 
    sig = list(sgRNAList['significant'])        
    counts = list(sgRNAList['counts'])
    control = list(sgRNAList['control mean'])        
    gene_list = list(set(genes))
    G = len(gene_list)
       
    # ------------------------------------------------
    # Make temporary dataframe (necessary since the original order must not change!)
    # ------------------------------------------------   
    global temp_df
    temp_df = pandas.DataFrame(data = {'gene': genes,
                                       'sgRNA_signif': sig,
                                       'counts': counts,
                                       'control mean': control},
                                       columns = ['gene','sgRNA_signif',\
                                       'counts','control mean'])                                               
    temp_df = temp_df.sort_values(['sgRNA_signif'],ascending=False)
    sig = list(temp_df['sgRNA_signif'])
    genes = list(temp_df['gene'])       # genes sorted (genes with no signif. sgRNAs last)    
    T0 = sig.index(False)       # index of first non-significant sgRNA     
    Genes = [genes[i] for i in range(T0)]      # genes with at least 1 signif. sgRNA
    GeneList = list(set(Genes))    
    

    # ------------------------------------------------    
    # Compute SigmaFC Score (non-parallel)
    # ------------------------------------------------          
    print('Computing accumulated fold-change scores ...')    
    metric = list()
    for g in range(G):
        gene = gene_list[g]
        if gene in GeneList: 
            I = [i for i in range(T0) if genes[i] == gene] 
            sg_table = temp_df.iloc[I]
            k = len(sg_table)         # number of signif. sgRNAs      
            Nx = list(sg_table['counts']) 
            N0 = list(sg_table['control mean'])
            pi = 0
            for j in range(k):                
                # compute fold-change weight   (discontinued)             
#                y = np.log10((Nx[j]+delta)/(N0[j]+delta)) # log fold-change
#                if y < FCmin:
#                    u = np.exp(Lambda*(y-FCmin))
#                else:
                u = 1
                pi_j = u * np.log10((Nx[j]+delta)/(N0[j]+delta))
                pi = pi + pi_j
        else:
            pi = 0; k = 0
        #metric.append(w[k]*pi)
        metric.append(k*pi)     # Avoid "out of range" error 

    # ------------------------------------------------    
    # Estimate SigmaFC Score distribution (Permutation)
    # ------------------------------------------------   
    print('Estimating null distribution ('+str(Np)+' permutations)...')
    P_set = np.random.choice(L,size=(Np,r),replace=True)
    SigmaFC_null = Parallel(n_jobs=num_cores)(delayed(SigmaFC_Permutation)(P) for P in P_set)
    ecdf = ECDF(SigmaFC_null,side='left')
    if ScreenType == 'enrichment':
        metric_pval = [1 - ecdf(metric[g]) for g in range(G)]    
    elif ScreenType == 'depletion':
        metric_pval = [ecdf(metric[g]) for g in range(G)]            
    else:
        print('### ERROR: Check spelling of ScreenType in configuration file! ###')
    metric_sig = [True if metric_pval[g] < alpha_g else False for g in range(G)]
    
    # ------------------------------------------------    
    # Return output
    # ------------------------------------------------   
    return metric, metric_pval, metric_sig