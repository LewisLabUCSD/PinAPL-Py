#!/usr/bin/python
#-*- coding: utf-8 -*-
"""
Created on Sat Apr 30 14:22:06 2016

@author: philipp
"""
# Find candidate sgRNAs
# =======================================================================
# Imports
from __future__ import division # floating point division by default
import pandas as pd
import numpy as np
import scipy
from scipy import stats
from decimal import *
import os
import glob
import sys
from statsmodels.sandbox.stats.multicomp import multipletests
import time
import yaml
import warnings
import sys



def ListCandidateGuides(sample):
    # ------------------------------------------------
    # Print header
    # ------------------------------------------------
    print('********************************************************')
    print('PinAPL-Py: List Candidate sgRNAs')
    print('********************************************************')
    start_total = time.time()

    # ------------------------------------------------
    # Get parameters
    # ------------------------------------------------
    os.chdir('/workingdir/')    
    configFile = open('configuration.yaml','r')
    config = yaml.load(configFile)
    configFile.close()
    WorkingDir = config['WorkingDir']
    AnalysisDir = config['AnalysisDir']
    QCDir = config['QCDir']
    CtrlDir = AnalysisDir + 'Control/'
    ListDir = AnalysisDir + 'Candidate_Lists/'
    CtrlPrefix = config['CtrlPrefix']
    CtrlCounts_Filename = 'Control_GuideCounts_0.tsv'
    ScreenType = config['ScreenType']
    alpha = config['alpha']    
    pcorr = config['pcorr']
    VarEst = config['VarEst']  

    
    # --------------------------------
    # Read control and sample data
    # --------------------------------   
    print('Preparing candidate sgRNA list for '+sample+' ...')     
    os.chdir(CtrlDir)
    Ctrl_File = pd.read_table(CtrlCounts_Filename, sep='\t')
    sgIDs = list(Ctrl_File['sgID'].values)
    genes = list(Ctrl_File['gene'].values)
    mu = list(Ctrl_File['Mean'].values)
    if VarEst == 'model':
        sigma2 = list(Ctrl_File['Model Variance'].values)
    elif VarEst == 'sample':
        sigma2 = list(Ctrl_File['Sample Variance'].values)   
    L = len(sgIDs)
    os.chdir(QCDir+sample)
    colnames = ['sgID','gene','counts']
    filename = glob.glob('*GuideCounts_0.tsv')[0]
    SampleFile = pd.read_table(filename, sep='\t',names=colnames)
    x = list(SampleFile['counts'].values)
     
    # -----------------------------------------------
    # Compute p-values and fold change
    # -----------------------------------------------
    # Compute fold change (compared to control)
    print('Computing fold changes ...')
    fc = list()
    for k in range(L):
        if mu[k] > 0:
            fc.append(x[k]/mu[k])
        else:
            fc.append(0)               
    # Compute negative binomial p-values
    print('Computing p-values ...')
    NBpval = list()
    if ScreenType == 'enrichment':
        for i in range(L):
            if mu[i]==0 and x[i]==0:
                NBpval.append(1)
            elif mu[i]==0 and x[i]>0:
                NBpval.append(1)
            elif mu[i]>0 and x[i]<=mu[i]:
                NBpval.append(1)
            elif x[i]>mu[i]:
                n = (mu[i]**2/sigma2[i])/(1-mu[i]/sigma2[i])
                p = mu[i]/sigma2[i]
                NBpval.append(1 - scipy.stats.nbinom.cdf(x[i],n,p))
            else:
                NBpval.append(1)
    elif ScreenType == 'depletion':
        for i in range(L):
            if mu[i]==0 and x[i]==0:
                NBpval.append(0)
            elif mu[i]==0 and x[i]>0:
                NBpval.append(1)
            elif mu[i]>0 and x[i]>=mu[i]:
                NBpval.append(1)
            elif x[i]<mu[i]:
                n = (mu[i]**2/sigma2[i])/(1-mu[i]/sigma2[i])
                p = mu[i]/sigma2[i]
                NBpval.append(scipy.stats.nbinom.cdf(x[i],n,p))
            else:
                NBpval.append(1)                 
    # Determine critical threshold
    NBpval_corr = multipletests(NBpval,alpha,pcorr)
    significant = NBpval_corr[0]
             
    # -----------------------------------------------
    # Save sgRNA dataframe
    # -----------------------------------------------
    # Write dataframe                     
    print('Writing results dataframe ...')
    Results_df = pd.DataFrame(data = {'sgRNA': [sgIDs[k] for k in range(L)],
                                     'gene': [genes[k] for k in range(L)],
                                     'counts': [x[k] for k in range(L)],
                                     'control mean': [np.rint(mu[k]) for k in range(L)],
                                     'control stdev': [np.rint(np.sqrt(sigma2[k])) for k in range(L)],
                                     'fold change': [round(fc[k]*100)/100 for k in range(L)],   
                                     'NB_pval': ['%.2E' % Decimal(NBpval[k]) for k in range(L)],
                                     'significant': [significant[k] for k in range(L)]},
                            columns = ['sgRNA','gene','counts','control mean','control stdev',\
                            'fold change','NB_pval','significant'])
    if ScreenType == 'enrichment':
        Results_df_0 = Results_df.sort_values(['significant','counts'],ascending=[0,0])                
    else:
        Results_df_0 = Results_df.sort_values(['significant','counts'],ascending=[0,1])            
    if not os.path.exists(ListDir):
        os.makedirs(ListDir)
    os.chdir(ListDir) 
    ListFilename = sample+'_'+str(alpha)+'_'+pcorr+'_sgRNAList.tsv'
    Results_df_0.to_csv(ListFilename, sep = '\t', index = False)
    os.chdir(QCDir)

    # --------------------------------------
    # Final time stamp
    # --------------------------------------        
    end_total = time.time()
    # Final time stamp
    print('------------------------------------------------')
    print('Script completed.')    
    sec_elapsed = end_total - start_total
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
    input1 = sys.argv[1]    
    ListCandidateGuides(input1) 