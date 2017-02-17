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
    print('++++++++++++++++++++++++++++++++++++++++++++++')
    print('PinAPL-Py: sgRNA Enrichment/Depletion Analysis')
    print('++++++++++++++++++++++++++++++++++++++++++++++')
    start_total = time.time()

    # ------------------------------------------------
    # Get parameters
    # ------------------------------------------------
    configFile = open('configuration.yaml','r')
    config = yaml.load(configFile)
    configFile.close()
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
    delta = config['delta_2']
    pvalDir = config['pvalDir']
    res = config['dpi']
    svg = config['svg']
    
    # --------------------------------
    # Read control and sample data
    # --------------------------------   
    print('Loading read counts ...')     
    os.chdir(CtrlDir)
    Ctrl_File = pandas.read_table(CtrlCounts_Filename, sep='\t')
    sgIDs = list(Ctrl_File['sgID'].values)
    genes = list(Ctrl_File['gene'].values)
    mu = list(Ctrl_File['Mean'].values)
    L = len(sgIDs)
    sigma2 = list(Ctrl_File['Model Variance'].values)
    if numpy.isnan(max(sigma2)): # no control replicates
        sigma2 = [0 for k in range(L)]
    os.chdir(InputDir)
    colnames = ['sgID','gene','counts']
    filename = glob.glob('*GuideCounts_0.tsv')[0]
    SampleFile = pandas.read_table(filename, sep='\t',names=colnames)
    x = list(SampleFile['counts'].values)
     
    # -----------------------------------------------
    # Compute p-values and fold change
    # -----------------------------------------------
    # Compute fold change (compared to control)
    print('Computing fold changes ...')
    fc = list()
    for k in range(L):
        if x[k]==0 or mu[k]==0:
            fc.append((x[k]+delta)/(mu[k]+delta))
        else:
            fc.append(x[k]/mu[k])
    # Compute negative binomial p-values    
    if max(sigma2) > 0: 
        print('Computing p-values ...')
        # Neg. Binom. Parameters  n: number of failures, p: probability of failure
        n = list(); p = list()
        for i in range(L):
            if mu[i]==0 or sigma2[i]==0:
                n.append(((mu[i]+delta)**2/(sigma2[i]+2*delta))/(1-(mu[i]+delta)/(sigma2[i]+2*delta)))
                p.append((mu[i]+delta)/(sigma2[i]+2*delta))
            else:
                n.append((mu[i]**2/sigma2[i])/(1-mu[i]/sigma2[i]))
                p.append(mu[i]/sigma2[i])
        NBpval = list(); 
        if ScreenType == 'enrichment':
            for i in range(L):
                if mu[i]==0 and x[i]==0:
                    NBpval.append(1)
                elif x[i]<=mu[i]:
                    NBpval.append(1)
                else: 
                    NBpval.append(1 - scipy.stats.nbinom.cdf(x[i],n[i],p[i]))
        elif ScreenType == 'depletion':
            for i in range(L):
                if mu[i]==0 and x[i]==0:
                    NBpval.append(1)
                elif x[i]>=mu[i]:
                    NBpval.append(1)                                                           
                else:
                    NBpval.append(scipy.stats.nbinom.cdf(x[i],n[i],p[i]))
        else:
            print('ERROR: Check spelling of ScreenType in configuration file!')
        # Compute two-sided pvalues (for volcano plot only!)
        NBpval2 = list()
        for i in range(L):
            if x[i]<=mu[i]:
                NBpval2.append(scipy.stats.nbinom.cdf(x[i],n[i],p[i]))
            else:
                NBpval2.append(1 - scipy.stats.nbinom.cdf(x[i],n[i],p[i]))                        
        # p-value correction for multiple tests
        print('p-value correction ...')
        multTest = multipletests(NBpval,alpha,padj)
        significant = multTest[0]
        NBpval_0 = multTest[1]
        # Plots
        print('Plotting p-values ...')
        pvalHist(NBpval,NBpval_0,pvalDir,sample,res,svg)
        VolcanoPlot(fc,NBpval2,significant,pvalDir,ScreenType,sample,res,svg)
        QQPlot(NBpval,significant,pvalDir,sample,res,svg)
        zScorePlot(fc,significant,pvalDir,ScreenType,sample,res,svg)
    else:         # no control replicates
        print('WARNING: No control replicates! No p-values computed...')
        NBpval = [-1 for k in range(L)]
        NBpval_0 = [-1 for k in range(L)]
        significant = ['N/A' for k in range(L)]
            
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
                                     'counts [norm.]': [x[k] for k in range(L)],
                                     'control mean [norm.]': [numpy.rint(mu[k]) for k in range(L)],
                                     'control stdev [norm.]': [numpy.rint(numpy.sqrt(sigma2[k])) for k in range(L)],
                                     'fold change': [fc[k] for k in range(L)],   
                                     'p-value': ['%.2E' % Decimal(NBpval[k]) for k in range(L)],
                                     'adj. p-value': ['%.2E' % Decimal(NBpval_0[k]) for k in range(L)],                                                 
                                     'significant': [str(significant[k]) for k in range(L)]},
                            columns = ['sgRNA','gene','counts [norm.]','control mean [norm.]',\
                            'control stdev [norm.]','fold change','p-value','adj. p-value','significant'])
    if ScreenType == 'enrichment':
        Results_df_0 = Results_df.sort_values(['significant','fold change'],ascending=[0,0])                
    elif ScreenType == 'depletion':
        Results_df_0 = Results_df.sort_values(['significant','fold change'],ascending=[0,1])            
    if SheetFormat == 'tsv':
        ListFilename = sample+'_'+str(alpha)+'_'+padj+'_sgRNAList.tsv'
        Results_df_0.to_csv(ListFilename, sep = '\t', index = False)
    elif SheetFormat == 'xlsx':
        print('Converting to xlsx ...')
        ListFilename = sample+'_'+str(alpha)+'_'+padj+'_sgRNAList.xlsx'
        Results_df_0.to_excel(ListFilename)

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
    PrepareHitList(input1) 
