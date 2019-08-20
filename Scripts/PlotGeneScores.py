#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 16:49:51 2019

@author: philipp
"""


# Scatterplot of gene scores (optional: highlighting a particular set of genes)
# =======================================================================
# Imports 
import sys
import yaml
import pandas
import numpy
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import time
import os
import glob
from matplotlib.ticker import FuncFormatter

def kilos(x, pos):
    return '%1.0fk' % (x*1e-3)

def GeneScoreScatterplot(sample,GOI='none'):
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
    GeneMetric = config['GeneMetric']
    RankingsDir = config['GeneDir']
    PlotDir = config['GenePlotDir']    
    ScriptsDir = config['ScriptsDir']
    res = config['dpi']
    svg = config['svg']
    NonTPrefix = config['NonTargetPrefix']
    
    # ------------------------------------------------
    # Read gene ranking input
    # ------------------------------------------------    
    os.chdir(RankingsDir)
    print('Reading gene ranking data....')
    filename = glob.glob(sample+'*_'+GeneMetric+'_GeneList.txt')[0]
    X = pandas.read_table(filename, sep='\t')
    X0 = X.sort_values(['gene'])        # sort gene table alphabetically
    genes = list(X0['gene'])    
    G = len(genes)
    metric = list(X0[GeneMetric])
    sig = list(X0['significant'])
    if GeneMetric == 'AvgLogFC':
        V = metric
        yAxisLabel = 'Avg. sgRNA LogFC'
    elif GeneMetric == 'STARS':
        V = metric
        yAxisLabel = 'STARS Score'
    elif GeneMetric == 'aRRA':
        V = [-numpy.log10(metric[g]) for g in range(G)]
        yAxisLabel = '-log(aRRA) Score'
    elif GeneMetric == 'SigmaFC':
        V = metric
        yAxisLabel = 'SigmaFC Score'
        

    # ------------------------------------------------      
    # Define sets
    # ------------------------------------------------  
    GeneIndex = range(G)
    SigIndex = [g for g in range(G) if sig[g]==True]
    NonTIndex = [g for g in range(G) if NonTPrefix in genes[g]]
    GOIIndex = [g for g in range(G) if genes[g]==GOI] 
    RestIndex = list(set(GeneIndex)-set(SigIndex)-set(NonTIndex)-set(GOIIndex))
    V_sig = [V[i] for i in SigIndex]
    V_nonT = [V[i] for i in NonTIndex]
    V_goi = [V[i] for i in GOIIndex]
    V_rest = [V[i] for i in RestIndex]    

    # ------------------------------------------------      
    # Plot
    # ------------------------------------------------  
    print('Plotting '+str(GeneMetric)+' scores ...')
    fig, ax = plt.subplots(figsize=(3.5,2.9))
    plt.scatter(RestIndex,V_rest,s=15,color='#99a399',lw=0,alpha=0.2)
    plt.scatter(SigIndex,V_sig,s=15,color='#54e84c',lw=0,alpha=0.2,label='significant',\
        rasterized=True)
    plt.scatter(NonTIndex,V_nonT,s=15,color='#ffd575',lw=0,alpha=0.2,label='non-targeting',\
        rasterized=True)
    if GOI != 'none':
        GOILabel = GOI
        plt.scatter(GOIIndex,V_goi,s=15,color='red',lw=0,alpha=1,label=GOILabel,rasterized=True)    
    formatter = FuncFormatter(kilos)
    ax.xaxis.set_major_formatter(formatter)  
    plt.xlabel('Gene Index', fontsize=11)
    plt.ylabel(yAxisLabel, fontsize=11)
    plt.title('Gene Ranking', fontsize=12)
    if len(V_sig)>0 and len(V_rest)>0:
        ymax = 1.75*max([max(V_sig),max(V_rest)])
        plt.ylim(top=ymax)
    elif len(V_rest)>0:
        ymax = 1.75*max(V_rest)
        plt.ylim(top=ymax)
    leg = plt.legend(loc='upper center', prop={'size':5})
    for lh in leg.legendHandles: lh.set_alpha(1)    
    plt.tight_layout()
    # Save figure
    if not os.path.exists(PlotDir):
    	os.makedirs(PlotDir)
    os.chdir(PlotDir)
    figurename = sample+'_GeneScores.png'
    if GOI != 'none':    
        figurename = figurename[:-4]+'_'+GOI+'.png'    
        if not os.path.exists(PlotDir+'/'+sample+'_Highlighted_Genes'):
        	os.makedirs(PlotDir+'/'+sample+'_Highlighted_Genes')         
        os.chdir(PlotDir+'/'+sample+'_Highlighted_Genes')  
    plt.savefig(figurename, dpi=res)  
    if svg:
        plt.savefig(figurename[:-4]+'.svg')     
    
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
        GeneScoreScatterplot(input1)            
    elif len(sys.argv) == 3:
        input1 = sys.argv[1]
        input2 = sys.argv[2]
        GeneScoreScatterplot(input1,input2)     