#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 19:15:32 2017

@author: philipp
"""

# Combine gene ranks across replicates
# =======================================================================
# Imports
import sys
import yaml
import os
import time
import pandas
import numpy
import scipy
from scipy import stats

def GeneRankCombination(treatment):
    print('++++++++++++++++++++++++++++++++++++++++++++++++')   
    
    # Get parameters    
    start = time.time()
    configFile = open('configuration.yaml','r')
    config = yaml.load(configFile)
    configFile.close()    
    alpha = config['alpha_g']
    padj = config['padj']
    metric = config['GeneMetric']
    Np = config['Np']
    GeneDir = config['GeneDir']    
    eps = 1e-16    
   
    # Read replicate files and compute Fisher's statistic
    os.chdir(GeneDir)
    treatment_files = [f for f in os.listdir(GeneDir) if treatment in f\
        and metric in f and 'combined' not in f]
    if len(treatment_files) > 1:
        treatment_files.sort()
        K = len(treatment_files)
        ResultTable = pandas.DataFrame()
        X1 = pandas.read_table(treatment_files[0], sep='\t')
        G = len(X1)        
        # Read replicates
        chi = list(numpy.zeros(G))    
        k = 0    
        for treatment_file in treatment_files:
            k+=1   
            print('Reading '+treatment+' replicate '+str(k)+' ...')            
            X = pandas.read_table(treatment_file, sep='\t')
            X0 = X.sort_values('gene',ascending=1)    
            genes = list(X0['gene'])
            ResultTable['gene'] = genes
            pval = list(X0['p_value'])
            ResultTable['p-value Repl. '+str(k)] = pval
            ln_pval = [numpy.log(pval[i]+eps) for i in range(G)]
            chi = numpy.add(chi,ln_pval)         
        
        # Combine p-values [REF 1]
        print('Computing Fisher statistic ...')
        chi = [-2*chi[i] for i in range(G)]
        ResultTable['Fisher Statistic'] = chi
        PVal = [1 - scipy.stats.chi2.cdf(chi[i],2*K) for i in range(G)]
        ResultTable['p-value combined'] = PVal
        significant = [PVal[i] < alpha for i in range(G)]
        ResultTable['significant'] = significant
        ResultTable = ResultTable.sort_values(['significant','p-value combined'],ascending=[0,1])
        print('Writing results dataframe ...')
        ResultFilename = treatment+'_combined_'+str(alpha)+'_'+str(padj)+'_'+str(metric)\
            +'_P'+str(Np)+'_GeneList.txt'
        ResultTable.to_csv(ResultFilename, sep = '\t', index = False)  
   
   # Time stamp
    end = time.time()
    sec_elapsed = end - start
    if sec_elapsed < 60:
        time_elapsed = sec_elapsed
        print('Time elapsed [secs]: ' + '%.3f' % time_elapsed)
    elif sec_elapsed < 3600:
        time_elapsed = sec_elapsed/60
        print('Time elapsed [mins]: ' + '%.3f' % time_elapsed)
    else:
        time_elapsed = sec_elapsed/3600
        print('Time elapsed [hours]: ' + '%.3f' % time_elapsed)    


        
if __name__ == "__main__":
    input1 = sys.argv[1]
    GeneRankCombination(input1)
    
# REFERENCES:
# [REF 1]: 'Fisher's Method' (wikipedia)    
    