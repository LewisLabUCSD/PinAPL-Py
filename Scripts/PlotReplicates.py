#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 09:41:39 2016

@author: philipp
"""
# Scatterplot with two replicates
# =======================================================================
# Imports 
from __future__ import division # floating point division by default
import os
import pandas as pd
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy
import glob
import yaml
import sys
import time
import scipy.stats

def Repl_Scatterplot(Repl1,Repl2):
    # ------------------------------------------------
    # Print header
    # ------------------------------------------------
    print('++++++++++++++++++++++++++++++++')
    print('PinAPL-Py: Replicate Scatterplot')
    print('++++++++++++++++++++++++++++++++')  
    print('Sample 1: '+Repl1+' | Sample 2: '+Repl2)
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
    AlnQCDir = config['AlnQCDir']
    PlotDir = config['CorrelDir']
    alpha = config['alpha']
    delta = config['delta']
    NonTPrefix = config['NonTargetPrefix']
    res = config['dpi']
    svg = config['svg']    
    dotsize = config['dotsize']
    logbase = config['logbase']    
        
    # ------------------------------------------------
    # Reading counts from sample and control
    # ------------------------------------------------
    print('Reading counts ...')    
    colnames = ['sgRNA','gene','counts']
    os.chdir(AlnQCDir+Repl1)
    filename1 = glob.glob('*_GuideCounts_0.tsv')[0]
    ListFile1 = pd.read_table(filename1, sep='\t',low_memory=False,names=colnames)    
    ListFile1 = ListFile1.sort_values('sgRNA')    
    os.chdir(AlnQCDir+Repl2)    
    filename2 = glob.glob('*_GuideCounts_0.tsv')[0]
    ListFile2 = pd.read_table(filename2, sep='\t',low_memory=False,names=colnames)
    ListFile2 = ListFile2.sort_values('sgRNA')    
    sgIDs = list(ListFile1['sgRNA'].values)
    genes = list(ListFile1['gene'].values)  
    L = len(sgIDs)
    repl1_counts = list(ListFile1['counts'].values)
    repl2_counts = list(ListFile2['counts'].values)
    # Log transformation
    print('Log'+str(logbase)+' transformation ...')    
    if logbase == 2:
        repl1_log = [numpy.log2(repl1_counts[k]+delta) for k in range(L)]
        repl2_log = [numpy.log2(repl2_counts[k]+delta) for k in range(L)]
    elif logbase == 10:
        repl1_log = [numpy.log10(repl1_counts[k]+delta) for k in range(L)]
        repl2_log = [numpy.log10(repl2_counts[k]+delta) for k in range(L)]        
    
    # ------------------------------------------------
    # Creating gene subsets
    # ------------------------------------------------    
    K_nonT = [k for k in range(L) if NonTPrefix in genes[k]]
    K_rest = list(set(range(L)) - set(K_nonT)) 
    repl1_nonT = [repl1_log[k] for k in K_nonT]
    repl2_nonT = [repl2_log[k] for k in K_nonT]
    repl1_rest = [repl1_log[k] for k in K_rest]
    repl2_rest = [repl2_log[k] for k in K_rest]

    # ----------------------------------------    
    # Compute correlation
    # ----------------------------------------
    print('Computing correlation ...')
    CorrCoeffP = numpy.corrcoef(repl1_log,repl2_log)[0,1]
    CorrCoeffS = scipy.stats.spearmanr(repl1_log,repl2_log)[0]           

    # ------------------------------------------------
    # Plotting
    # ------------------------------------------------       
    print('Generating scatterplot ...')
    if not os.path.exists(PlotDir):
        os.makedirs(PlotDir)      
    os.chdir(PlotDir)   
    plt.figure(figsize=(6,5))
    plt.scatter(repl1_rest,repl2_rest,s=dotsize,facecolor='black',lw=0,alpha=0.35)
    if len(K_nonT)>0:
        plt.scatter(repl1_nonT,repl2_nonT,s=dotsize,facecolor='orange',lw=0,alpha=0.75,\
            label='Non Targeting')
    axes = plt.gca()
    x0 = axes.get_xlim()  
    plt.plot((x0[0],x0[1]), (x0[0],x0[1]), ls="--", color=(51/255,153/255,1))
    plt.title('Correlation '+Repl1+' '+Repl2, fontsize=14)
    plt.xlabel(Repl1+' log'+str(logbase)+' counts [norm.]', fontsize=12)    
    plt.ylabel(Repl2+' log'+str(logbase)+' counts [norm.]', fontsize=12)
    plt.legend(loc='upper left', prop={'size':10})
    plt.text(.6,.2,'Corr (Pearson) = '+str(round(CorrCoeffP*1000)/1000),transform=axes.transAxes,\
        fontsize=10) 
    plt.text(.6,.15,'Corr (Spearman) = '+str(round(CorrCoeffS*1000)/1000),transform=axes.transAxes,\
        fontsize=10)    
    plt.tight_layout()  
    plt.savefig(Repl1+' '+Repl2+' correlation.png', dpi=res)  
    if svg:
        plt.savefig(Repl1+' '+Repl2+' correlation.svg')        
    plt.close()

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
    input2 = sys.argv[2]
    Repl_Scatterplot(input1,input2)    
