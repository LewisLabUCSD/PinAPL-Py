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


def PrepareGuideRanking(sample):
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
    sgRNAReadCountDir = config['sgRNAReadCountDir']
    CtrlDir = config['ControlDir']
    ListDir = config['sgRNARanksDir']
    CtrlCounts_Filename = 'Control_GuideCounts_normalized.txt'
    ScreenType = config['ScreenType']
    alpha = config['alpha_s']
    padj = config['padj']
    SheetFormat = config['HitListFormat']
    delta = config['delta']
    pvalDir = config['pvalDir_sgRNA']
    zScoreDir = config['zScoreDir_sgRNA']
    res = config['dpi']
    svg = config['svg']
    
    # --------------------------------
    # Read control and sample data
    # --------------------------------   
    print('Loading sgRNA read counts ...')     
    os.chdir(CtrlDir)
    Ctrl_File = pandas.read_table(CtrlCounts_Filename, sep='\t')
    Model = Ctrl_File['Model'][0]
    sgIDs = list(Ctrl_File['sgID'])
    genes = list(Ctrl_File['gene'])
    mu = list(Ctrl_File['Mean'])
    L = len(sgIDs)
    SampleVar = list(Ctrl_File['Sample Variance'])
    sigma2 = list(Ctrl_File['Model Variance'])
    n = list(Ctrl_File['n'])
    p = list(Ctrl_File['p'])    
    os.chdir(sgRNAReadCountDir)
    colnames = ['sgID','gene','counts']
    filename = sample+'_GuideCounts_normalized.txt'
    SampleFile = pandas.read_table(filename, sep='\t',names=colnames)
    x = list(SampleFile['counts'])
     
    # -----------------------------------------------
    # Compute fold change 
    # -----------------------------------------------    
    print('Computing fold-changes ...')
    fc = list()
    for k in range(L):
        if x[k]==0 or mu[k]==0:
            fc.append((x[k]+delta)/(mu[k]+delta))
        else:
            fc.append(x[k]/mu[k])     
     
    # -----------------------------------------------
    # Compute p-values 
    # -----------------------------------------------              
    if Model == 'none':        
    # -----------------------------------------------------------
        print('WARNING: Zero variance or no control replicates! Cannot compute p-values ...')   
        pval = ['N/A' for k in range(L)]
        pval0 = ['N/A' for k in range(L)]
        significant = ['N/A' for k in range(L)]     
    # -----------------------------------------------------------
    elif ScreenType == 'enrichment':       # enrichment screen
    # -----------------------------------------------------------
        pval = list(); 
        print('Computing p-values ('+Model+' model) ...')
        # one-sided p-value
        if Model == 'Neg. Binomial':
            for k in range(L):  
                pval.append(1 - scipy.stats.nbinom.cdf(x[k],n[k],p[k]))
        elif Model == 'Poisson':
            for k in range(L):
                pval.append(1 - scipy.stats.poisson.cdf(x[k],sigma2[k]))
    # -----------------------------------------------------------                 
    elif ScreenType == 'depletion':       # depletion screen        
    # -----------------------------------------------------------
        pval = list();
        print('Computing p-values ('+Model+' model) ...')
        # one-sided p-value
        if Model == 'Neg. Binomial':
            for k in range(L):         
                pval.append(scipy.stats.nbinom.cdf(x[k],n[k],p[k]))
        elif Model == 'Poisson':
                pval.append(scipy.stats.poisson.cdf(x[k],sigma2[k]))        
    # -----------------------------------------------------------                  
    else:                           # error in scree type
    # -----------------------------------------------------------   
        print('### ERROR: Check spelling of ScreenType in configuration file! ###')

    # -----------------------------------------------
    # p-value Correction and Plots
    # -----------------------------------------------  
    if max(SampleVar) > 0:
        # p-value correction for multiple tests
        print('Multiple tests correction ...')
        pval_0 = list()
        for k in range(L):
            condition1 = x[k]!=0 or mu[k]!=0        # prevent artificial peak at p=0.5 because of fc=1
            condition2 = ~numpy.isnan(pval[k])      # remove nan (shouldn't have any though...)
            if condition1 and condition2:            
                pval_0.append(pval[k])
        multTest = multipletests(pval_0,alpha,padj)
        significant_0 = multTest[0]    
        significant = list()        # map back to full list
        j=-1
        for k in range(L):
            condition1 = x[k]!=0 or mu[k]!=0 
            condition2 = ~numpy.isnan(pval[k])
            if condition1 and condition2:
                j += 1
                significant.append(significant_0[j])
            else:
                significant.append(False)
        # Plots
        print('Plotting p-value histogram ...')
        PlotTitle = 'sgRNA '+ScreenType.capitalize()
        pvalHist(pval_0,pvalDir,sample,res,svg,'#0be52c',PlotTitle)

               
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
                                     'p-value': [pval[k] for k in range(L)],                                                
                                     'significant': [str(significant[k]) for k in range(L)]},
                            columns = ['sgRNA','gene','counts','control mean',\
                            'control stdev','fold change','p-value','significant'])
    if ScreenType == 'enrichment': 
        Results_df_0 = Results_df.sort_values(['significant','p-value','fold change','sgRNA'],ascending=[0,1,0,1])
    elif ScreenType == 'depletion': 
        Results_df_0 = Results_df.sort_values(['significant','p-value','fold change','sgRNA'],ascending=[0,1,1,1])
    ListFilename = sample+'_'+str(alpha)+'_'+padj+'_sgRNAList.txt'
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
        print('Time elapsed [secs]: ' + '%.3f' % time_elapsed +'\n')
    elif sec_elapsed < 3600:
        time_elapsed = sec_elapsed/60
        print('Time elapsed [mins]: ' + '%.3f' % time_elapsed +'\n')
    else:
        time_elapsed = sec_elapsed/3600
        print('Time elapsed [hours]: ' + '%.3f' % time_elapsed +'\n')

    
if __name__ == "__main__":
    input1 = sys.argv[1]    
    PrepareGuideRanking(input1) 