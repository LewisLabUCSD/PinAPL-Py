#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 09:41:39 2016

@author: philipp
"""
# Scatterplot highlighting a particular gene of interest
# =======================================================================
# Imports 
from __future__ import division # floating point division by default
import os
import pandas as pd
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import glob
import yaml
import sys
import time

def GOI_Scatterplot(sample,GOI='None'):
    # ------------------------------------------------
    # Print header
    # ------------------------------------------------
    print('+++++++++++++++++++++++++++++++++++++++++++')
    print('PinAPL-Py: Read Counts Scatterplot')
    print('+++++++++++++++++++++++++++++++++++++++++++')  
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
    ListDir = config['HitDir']
    PlotDir = config['ScatterDir']
    annotate = config['scatter_annotate']
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
    os.chdir(ListDir)
    filename = glob.glob(sample+'_*sgRNAList.tsv')[0]
    HitList = pd.read_table(filename, sep='\t')       
    sgIDs = list(HitList['sgRNA'].values)
    genes = list(HitList['gene'].values)  
    L = len(sgIDs)
    sample_counts = list(HitList['counts [norm.]'].values)
    control_counts = list(HitList['control mean [norm.]'].values)
    sig = list(HitList['significant'].values)  
    # Log transformation
    print('Log'+str(logbase)+' transformation ...')
    if logbase == 2:
        sample_log = [np.log2(sample_counts[k]+delta) for k in range(len(sample_counts))]
        control_log = [np.log2(control_counts[k]+delta) for k in range(len(control_counts))]
    elif logbase == 10:
        sample_log = [np.log10(sample_counts[k]+delta) for k in range(len(sample_counts))]
        control_log = [np.log10(control_counts[k]+delta) for k in range(len(control_counts))]
    
    # ------------------------------------------------
    # Creating gene subsets
    # ------------------------------------------------    
    K_nonT = [k for k in range(L) if NonTPrefix in genes[k]]
    K_sig = [k for k in range(L) if sig[k]==True]
    K_goi = [k for k in range(L) if genes[k] == GOI]   
    K_rest = list(set(range(L)) - set.union(set(K_nonT),set(K_sig),set(K_goi)))   
    sample_nonT = [sample_log[k] for k in K_nonT]
    control_nonT = [control_log[k] for k in K_nonT]
    sample_sig = [sample_log[k] for k in K_sig]
    control_sig = [control_log[k] for k in K_sig]
    sample_goi = [sample_log[k] for k in K_goi]
    control_goi = [control_log[k] for k in K_goi]
    sample_rest = [sample_log[k] for k in K_rest]
    control_rest = [control_log[k] for k in K_rest]
    goi_sgIDs = [sgIDs[k] for k in K_goi]

    # ------------------------------------------------
    # Plotting
    # ------------------------------------------------       
    print('Generating scatterplot ...')
    if not os.path.exists(PlotDir):
        os.makedirs(PlotDir)      
    os.chdir(PlotDir)   
    plt.figure(figsize=(6,5))
    plt.scatter(control_rest,sample_rest,s=dotsize,facecolor='black',lw=0,alpha=0.35)
    plt.scatter(control_sig,sample_sig,s=dotsize,facecolor='green',lw=0,alpha=0.35,label='Significant')
    if GOI != 'None':
        plt.scatter(control_goi,sample_goi,s=1.5*dotsize,facecolor='red',lw=0,alpha=0.35,label=GOI)
    if len(K_nonT)>0:
        plt.scatter(control_nonT,sample_nonT,s=dotsize,facecolor='orange',lw=0,alpha=0.75,\
            label='Non Targeting')
    axes = plt.gca()
    x0 = axes.get_xlim()  
    plt.plot((x0[0],x0[1]), (x0[0],x0[1]), ls="--", color=(51/255,153/255,1))
    plt.title(sample+' log'+str(logbase)+' counts [norm.]', fontsize=14)
    plt.xlabel('Control (avg.)', fontsize=12)    
    plt.ylabel(sample, fontsize=12)
    plt.legend(loc='upper left', prop={'size':10})
    if annotate:
        for label, x, y in zip(goi_sgIDs,control_goi,sample_goi):
            plt.annotate(label,xy=(x,y),color='red',fontsize=8)  
    plt.tight_layout()
    if GOI != 'None':
        plt.savefig(sample+' '+GOI+' counts.png', dpi=res)
        if svg:
            plt.savefig(sample+' '+GOI+' counts.svg')
    else:
        plt.savefig(sample+' '+' counts.png', dpi=res)
        if svg:
            plt.savefig(sample+' '+' counts.svg')
    plt.close()

    # ------------------------------------------------
    # Printing 
    # ------------------------------------------------
    if GOI != 'None':
        print('-----------------------------------------------')     
        print('sgID\t\tCounts\tControl\tSignificant')    
        print('-----------------------------------------------')       
        if not K_goi:
            print('ERROR: Gene name not found!')
        else:            
            for k in K_goi:        
                println = str(sgIDs[k])+'\t'+str(int(sample_counts[k]))+'\t'+ \
                    str(int(control_counts[k]))+'\t'+str(sig[k])
                print(println)

    # --------------------------------------
    # Time stamp
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
    if len(sys.argv) == 3:
        input1 = sys.argv[1]
        input2 = sys.argv[2]
        GOI_Scatterplot(input1,input2)    
    else:
        input1 = sys.argv[1]
        GOI_Scatterplot(input1)            
