#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 17:15:19 2019

@author: philipp
"""

# =======================================================================
# Sort genes by STARS (Doench et al., Nat Biotech 2016)
# =======================================================================
from joblib import Parallel, delayed
import multiprocessing
from scipy.stats import beta
from statsmodels.distributions.empirical_distribution import ECDF
import yaml
import pandas
import numpy
import os
import glob


def compute_STARS(HitList): 
    configFile = open('configuration.yaml','r')
    config = yaml.load(configFile)
    configFile.close()      
    STARSDir = config['STARSDir']
    thr = config['thr_STARS']
    screentype = config['ScreenType'] 
    alpha_g = config['alpha_g']
    Np = 10
    L = len(HitList)
    sgIDs = list(HitList['sgRNA'])
    genes = list(HitList['gene'])
    geneList = list(set(genes))
    G = len(geneList)    
    fc = list(HitList['fold change'])
    # Estimating null distribution
    print('Estimating STARS null distribution ('+str(Np)+' permutations)...')
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
    geneList_s = list(STARS['Gene Symbol'])
    metric_s = list(STARS['STARS Score'])
    metric_pval_s = list(STARS['p-value'])               
    # Add N/A for genes not reported by STARS    
    metric = list()
    metric_pval = list()
    for g in range(G):
        gene = geneList[g]
        if gene in geneList_s:
            gene_index = geneList_s.index(gene)
            metric.append(metric_s[gene_index])
            metric_pval.append(metric_pval_s[gene_index])
        else:
            metric.append(-1)
            metric_pval.append('N/A')
    sig = [True if metric_pval[g] < alpha_g else False for g in range(G)]
    # Deleting STARS files            
    os.system('rm STARS_input.txt STARS_chip.txt')
    os.system('rm Null_STARSOutput8_'+str(thr)+'.txt')
    os.system('rm '+STARS_output)
    # Return output
    return metric, metric_pval, sig
    return 