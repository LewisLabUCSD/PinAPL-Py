#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu May 30 17:26:55 2019

@author: philipp
"""
# Volcano plot of fold change
# =======================================================================
# Imports 
from __future__ import division # floating point division by default
import sys
import time
import os
import pandas
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy
import yaml
import glob
from matplotlib.ticker import FuncFormatter

def VolcanoPlot(sample,GOI='none',Annot='none',NonT='none'):
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
    sgRNARanksDir = config['sgRNARanksDir']
    ScriptsDir = config['ScriptsDir']
    outputDir = config['VolcanoDir_sgRNA']
    ScreenType = config['ScreenType']
    res = config['dpi']
    svg = config['svg']
    NonTPrefix = config['NonTargetPrefix']
    PrintHighlights = config['PrintHighlights']
    # Show non-targeting controls?
    if NonT == 'none':
        ShowNonTargets = config['ShowNonTargets']
    elif NonT == 'False':
        ShowNonTargets = False
    elif NonT == 'True':
        ShowNonTargets = True   
    # Annotate sgRNAs ?
    if Annot == 'none':
        annotate = config['scatter_annotate']    
    elif Annot == 'False':
        annotate = False
    elif Annot == 'True':
        annotate = True        

    # ------------------------------------------------
    # Reading fold-change data
    # ------------------------------------------------
    print('Reading sgRNA read counts ...')    
    os.chdir(sgRNARanksDir)
    filename = glob.glob(sample+'_*sgRNAList.txt')[0]
    sgRNARanking = pandas.read_table(filename, sep='\t')
    L = len(sgRNARanking)    
    if ScreenType == 'enrichment':
        sgRNARanking = sgRNARanking.sort_values('fold change',ascending=True)
    elif ScreenType == 'depletion':
        sgRNARanking = sgRNARanking.sort_values('fold change',ascending=False)            
    fc = list(sgRNARanking['fold change'])
    sig = list(sgRNARanking['significant'])    
    genes = list(sgRNARanking['gene'])
    sgIDs = list(sgRNARanking['sgRNA'])
    pval = list(sgRNARanking['p-value'])
    
    # ------------------------------------------------
    # Log transformation
    # ------------------------------------------------  
    print('Log transformation ...')    
    eps = 1e-16
    logfc = [numpy.log2(fc[k]) for k in range(L)]
    neglogp = [-numpy.log10(pval[k]+eps) for k in range(L)]

    # ------------------------------------------------
    # Creating subsets
    # ------------------------------------------------ 
    K_nonT = [k for k in range(L) if NonTPrefix in genes[k]]
    K_sig = [k for k in range(L) if sig[k]==True]
    K_goi = [k for k in range(L) if genes[k] == GOI]
    K_rest = list(set(range(L)) - set.union(set(K_nonT),set(K_sig),set(K_goi)))
    FC_nonT = [logfc[k] for k in K_nonT]
    p_nonT = [neglogp[k] for k in K_nonT]
    FC_sig = [logfc[k] for k in K_sig]
    p_sig = [neglogp[k] for k in K_sig]
    FC_goi = [logfc[k] for k in K_goi]
    p_goi = [neglogp[k] for k in K_goi]
    FC_rest = [logfc[k] for k in K_rest]
    p_rest = [neglogp[k] for k in K_rest]  
    sgIDs_goi = [sgIDs[k] for k in K_goi]

    # ------------------------------------------------
    # Plot
    # ------------------------------------------------    
    print('Generating volcano plot ...')    
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    os.chdir(outputDir)  
    fig, ax = plt.subplots(figsize=(3.5,2.9))
    if len(FC_sig)>100:
        tpcy = 0.35
    else:
        tpcy = 1
    plt.scatter(FC_rest,p_rest,s=6,facecolor='grey',lw=0,alpha=0.25,rasterized=True)
    plt.scatter(FC_sig,p_sig,s=6,facecolor='green',lw=0,alpha=tpcy,label='significant',rasterized=True) 
    if len(K_nonT)>0 and ShowNonTargets:
        plt.scatter(FC_nonT,p_nonT,s=6,lw=0,color='orange',alpha=0.15,label='non-targeting',rasterized=True)
    if GOI != 'none':
        plt.scatter(FC_goi,p_goi,s=6,lw=0,color='red',label=GOI,rasterized=True)   
    if annotate:
        for label, x, y in zip(sgIDs_goi,FC_goi,p_goi):
            plt.annotate(label,xy=(x,y),color='red',fontsize=4,fontweight='bold')          
    plt.xlabel('sgRNA log2 Fold-Change', fontsize=11)
    plt.ylabel('-log10 p-value', fontsize=11)
    plt.tick_params(labelsize=11)
    plt.title('sgRNA '+ScreenType.capitalize(), fontsize=12) 
    if ScreenType == 'enrichment':
        leg = plt.legend(loc='upper left', prop={'size':6})
        for lh in leg.legendHandles: lh.set_alpha(1)
    elif ScreenType == 'depletion':
        leg = plt.legend(loc='upper right', prop={'size':6})        
        for lh in leg.legendHandles: lh.set_alpha(1)
    plt.tight_layout()
    # Define file name
    figurename = sample+'_'+'sgRNA_volcano.png'
    if GOI != 'none':    
        figurename = figurename[:-4]+'_'+GOI+'.png'
    if annotate:
        figurename = figurename[:-4]+'_IDs.png'        
    if ShowNonTargets:
        figurename = figurename[:-4]+'_nonT.png'        
    # save figure    
    if GOI != 'none':       
        if not os.path.exists(outputDir+'/'+sample+'_Highlighted_Genes'):
        	os.makedirs(outputDir+'/'+sample+'_Highlighted_Genes')         
        os.chdir(outputDir+'/'+sample+'_Highlighted_Genes')  
    plt.savefig(figurename, dpi=res)  
    if svg:
        plt.savefig(figurename[:-4]+'.svg') 


    # ------------------------------------------------
    # Printing 
    # ------------------------------------------------
    if GOI != 'none' and PrintHighlights:
        print('-----------------------------------------------------------')     
        print('sgID\t\tFold-Change\tp-value\t\tSignificant')    
        print('-----------------------------------------------------------')
        if not K_goi:
            print('### ERROR: Gene name not found! ###')
        else:            
            for k in K_goi:        
                println = str(sgIDs[k])+'\t'+str(fc[k])+'\t'+ \
                    str(pval[k])+'\t'+str(sig[k])
                print(println)


    # Final time stamp
    os.chdir(ScriptsDir)
    end_total = time.time()
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
    if len(sys.argv) == 2:
        input1 = sys.argv[1]
        VolcanoPlot(input1)            
    elif len(sys.argv) == 3:
        input1 = sys.argv[1]
        input2 = sys.argv[2]
        VolcanoPlot(input1,input2)    
    elif len(sys.argv) == 4:
        input1 = sys.argv[1]
        input2 = sys.argv[2]
        input3 = sys.argv[3]    
        VolcanoPlot(input1,input2,input3) 
    elif len(sys.argv) == 5:
        input1 = sys.argv[1]
        input2 = sys.argv[2]
        input3 = sys.argv[3]
        input4 = sys.argv[4]
        VolcanoPlot(input1,input2,input3,input4)        