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
import pandas
import numpy
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
from pvalPlots import *


def PrepareHitList(sample):
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
    WorkingDir = config['WorkingDir']
    AnalysisDir = config['AnalysisDir']
    InputDir = config['AlnQCDir']+sample
    CtrlDir = AnalysisDir + 'Control/'
    ListDir = config['HitDir']
    CtrlCounts_Filename = 'Control_GuideCounts_0.tsv'
    ScreenType = config['ScreenType']
    alpha = config['alpha']    
    padj = config['padj']
    SheetFormat = config['HitListFormat']
    delta = config['delta']
    pvalDir = config['pvalDir']
    res = config['dpi']
    svg = config['svg']
    
    # --------------------------------
    # Read control and sample data
    # --------------------------------   
    print('Loading read counts ...')     
    os.chdir(CtrlDir)
    Ctrl_File = pandas.read_table(CtrlCounts_Filename, sep='\t')
    sgIDs = list(Ctrl_File['sgID'])
    genes = list(Ctrl_File['gene'])
    mu = list(Ctrl_File['Mean'])
    L = len(sgIDs)
    sigma2 = list(Ctrl_File['Model Variance'])
    n = list(Ctrl_File['n'])
    p = list(Ctrl_File['p'])    
    os.chdir(InputDir)
    colnames = ['sgID','gene','counts']
    filename = glob.glob('*GuideCounts_0.tsv')[0]
    SampleFile = pandas.read_table(filename, sep='\t',names=colnames)
    x = list(SampleFile['counts'])
     
    # -----------------------------------------------
    # Compute fold change and p-values
    # -----------------------------------------------    
    if max(sigma2) == 0:        # check for control replicates
    # -----------------------------------------------------------
        print('WARNING: No control replicates! No p-values computed...')
        print('Computing fold changes ...')
        fc = list()
        for k in range(L):
            if x[k]==0 or mu[k]==0:
                fc.append((x[k]+delta)/(mu[k]+delta))
            else:
                fc.append(x[k]/mu[k])        
        NBpval = [1 for k in range(L)]
        NBpval_0 = [1 for k in range(L)]
        significant = [False for k in range(L)]     
    # -----------------------------------------------------------
    elif ScreenType == 'enrichment':       # enrichment screen
    # -----------------------------------------------------------
        fc = list(); NBpval = list(); NBpval2 = list()
        print('Computing fold-changes and p-values...')
        for k in range(L):
            # fold-change            
            if x[k]==0 or mu[k]==0:
                fc.append((x[k]+delta)/(mu[k]+delta))
            else:
                fc.append(x[k]/mu[k])
            # one-sided p-value
            if mu[k]==0 and x[k]==0:
                  NBpval.append(1)
            elif x[k]<=mu[k]:
                  NBpval.append(1)
            else: 
                  NBpval.append(1 - scipy.stats.nbinom.cdf(x[k],n[k],p[k]))    
    # -----------------------------------------------------------                 
    elif ScreenType == 'depletion':       # depletion screen        
    # -----------------------------------------------------------
        fc = list(); NBpval = list(); NBpval2 = list()
        print('Computing fold-changes and p-values...')
        for k in range(L):
            # fold-change
            if x[k]==0 or mu[k]==0:
               fc.append((x[k]+delta)/(mu[k]+delta))
            else:
                fc.append(x[k]/mu[k])
            # one-sided p-value
            if mu[k]==0 and x[k]==0:
                NBpval.append(1)
            elif x[k]>=mu[k]:
                NBpval.append(1)
            else:
                NBpval.append(scipy.stats.nbinom.cdf(x[k],n[k],p[k]))                    
    # -----------------------------------------------------------                  
    else:                           # error in scree type
    # -----------------------------------------------------------   
        print('ERROR: Check spelling of ScreenType in configuration file!')

    # -----------------------------------------------
    # p-value Correction and Plots
    # -----------------------------------------------  
    if max(sigma2) > 0:
        # p-value correction for multiple tests
        print('p-value correction ...')
        multTest = multipletests(NBpval,alpha,padj)
        significant = multTest[0]
        NBpval_0 = multTest[1]
        # Plots
        print('Plotting p-values ...')
        pvalHist(NBpval,NBpval_0,pvalDir,sample,res,svg)
        VolcanoPlot(fc,NBpval,significant,pvalDir,ScreenType,sample,res,svg,alpha)
        QQPlot(NBpval,significant,pvalDir,sample,res,svg,alpha)
        zScorePlot(fc,significant,pvalDir,ScreenType,sample,res,svg,alpha)

               
    # -----------------------------------------------
    # Save sgRNA dataframe
    # -----------------------------------------------
    # Write dataframe                     
    if not os.path.exists(ListDir):
        os.makedirs(ListDir)
    os.chdir(ListDir)             
    print('Writing results dataframe ...')
    Results_df = pandas.DataFrame(data = {'sgRNA': [sgIDs[k] for k in range(L)],
                                     'gene': [genes[k] for k in range(L)],
                                     'counts': [x[k] for k in range(L)],
                                     'control mean': [mu[k] for k in range(L)],
                                     'control stdev': [numpy.sqrt(sigma2[k]) for k in range(L)],
                                     'fold change': [fc[k] for k in range(L)],   
                                     'p-value': [NBpval[k] for k in range(L)],
                                     'FDR': [NBpval_0[k] for k in range(L)],                                                 
                                     'significant': [str(significant[k]) for k in range(L)]},
                            columns = ['sgRNA','gene','counts','control mean',\
                            'control stdev','fold change','p-value','FDR','significant'])
    if ScreenType == 'enrichment': 
        Results_df_0 = Results_df.sort_values(['significant','p-value','fold change','sgRNA'],ascending=[0,1,0,1])
    elif ScreenType == 'depletion': 
        Results_df_0 = Results_df.sort_values(['significant','p-value','fold change','sgRNA'],ascending=[0,1,1,1])
    ListFilename = sample+'_'+str(alpha)+'_'+padj+'_sgRNAList.tsv'
    Results_df_0.to_csv(ListFilename, sep = '\t', index = False)
    if SheetFormat == 'xlsx':
        print('Converting to xlsx ...')
        ListFilename = sample+'_'+str(alpha)+'_'+padj+'_sgRNAList.xlsx'
        Results_df_0.to_excel(ListFilename)

    # --------------------------------------
    # Final time stamp
    # --------------------------------------        
    os.chdir(ScriptsDir)
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
    PrepareHitList(input1) 
