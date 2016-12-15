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
    print('***********************************************************')
    print('PinAPL-Py: Replicate Scatterplot')
    print('***********************************************************')  
    start_total = time.time()  
    
    # ------------------------------------------------
    # Get parameters
    # ------------------------------------------------
    os.chdir('/workingdir/')    
    configFile = open('configuration.yaml','r')
    config = yaml.load(configFile)
    configFile.close()
    WorkingDir = config['WorkingDir'] 
    AnalysisDir = WorkingDir + 'Analysis/'
    ListDir = AnalysisDir + 'Candidate_Lists/'
    PlotDir = AnalysisDir + 'ScatterPlots/'
    alpha = config['alpha']
    pcorr = config['pcorr']
    delta = config['delta_s']
    NonTPrefix = config['NonTargetPrefix']
    res = config['dpi']
    dotsize = config['dotsize']
   
    # ------------------------------------------------
    # Reading counts from sample and control
    # ------------------------------------------------
    print('Reading counts ...')    
    os.chdir(ListDir)
    filename1 = glob.glob(Repl1+'_'+str(alpha)+'_'+str(pcorr)+'_sgRNAList.tsv')[0]
    filename2 = glob.glob(Repl2+'_'+str(alpha)+'_'+str(pcorr)+'_sgRNAList.tsv')[0]    
    ListFile1 = pd.read_table(filename1, sep='\t',low_memory=False)
    ListFile2 = pd.read_table(filename2, sep='\t',low_memory=False)
    ListFile1 = ListFile1.sort_values('sgRNA')
    ListFile2 = ListFile2.sort_values('sgRNA')
    sgIDs = list(ListFile1['sgRNA'].values)
    genes = list(ListFile1['gene'].values)  
    L = len(sgIDs)
    repl1_counts = list(ListFile1['counts'].values)
    repl2_counts = list(ListFile2['counts'].values)
    # Log transformation
    repl1_log = [numpy.log(repl1_counts[k]+delta) for k in range(len(repl1_counts))]
    repl2_log = [numpy.log(repl2_counts[k]+delta) for k in range(len(repl2_counts))]
    
    # ------------------------------------------------
    # Creating gene subsets
    # ------------------------------------------------    
    print('Creating gene subsets ...') 
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
    Corr_p = scipy.stats.spearmanr(repl1_log,repl2_log)[1]          

    # ------------------------------------------------
    # Plotting
    # ------------------------------------------------       
    print('Generating scatterplot ...')
    if not os.path.exists(PlotDir):
        os.makedirs(PlotDir)      
    os.chdir(PlotDir)   
    plt.figure()
    plt.scatter(repl1_rest,repl2_rest,s=dotsize,facecolor='black',lw=0)
    if len(K_nonT)>0:
        plt.scatter(repl1_nonT,repl2_nonT,s=dotsize,facecolor=(255/255,0,255/255),lw=0,label='Non Targeting')
    axes = plt.gca()
    x0 = axes.get_xlim()  
    plt.plot((0,x0[1]-1), (0,x0[1]-1), ls="--", color=(51/255,153/255,1))
    plt.suptitle('Correlation '+Repl1+' '+Repl2, fontsize=12, fontweight='bold')
    plt.xlabel('log counts '+Repl1, fontsize=14)    
    plt.ylabel('log counts '+Repl2, fontsize=14)
    plt.legend(loc='upper left', prop={'size':10})
    plt.text(.6,.2,'Corr (Pearson) = '+str(round(CorrCoeffP*1000)/1000),transform=axes.transAxes) 
    plt.text(.6,.15,'Corr (Spearman) = '+str(round(CorrCoeffS*1000)/1000),transform=axes.transAxes) 
    plt.text(.6,.1,'p-val (Spearman) = '+str(round(Corr_p*1000)/1000),transform=axes.transAxes)      
    plt.savefig(Repl1+' '+Repl2+' correlation.png', dpi=res)        
    plt.close()

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
    input2 = sys.argv[2]
    Repl_Scatterplot(input1,input2)    