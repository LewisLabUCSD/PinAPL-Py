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
from matplotlib.ticker import FormatStrFormatter


def GOI_Scatterplot(sample,GOI='none',Annot='none',NonT='none',Transp='none'):
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
    sgRNARanksDir = config['sgRNARanksDir']
    PlotDir = config['ScatterDir']
    alpha = config['alpha_s']
    delta = config['delta']
    NonTPrefix = config['NonTargetPrefix']
    res = config['dpi']
    svg = config['svg']
    dotsize = config['dotsize']
    logbase = config['logbase']
    ScreenType = config['ScreenType']
    PrintHighlights = config['PrintHighlights']
    if Annot == 'none':
        annotate = config['scatter_annotate']    
    elif Annot == 'False':
        annotate = False
    elif Annot == 'True':
        annotate = True
    if NonT == 'none':
        ShowNonTargets = config['ShowNonTargets']
    elif NonT == 'False':
        ShowNonTargets = False
    elif NonT == 'True':
        ShowNonTargets = True
    if Transp == 'none':        
        TransparencyLevel = config['TransparencyLevel']
    else:
        TransparencyLevel = float(Transp)          
    
    # ------------------------------------------------
    # Reading counts from sample and control
    # ------------------------------------------------
    print('Reading sgRNA read counts ...')    
    os.chdir(sgRNARanksDir)
    filename = glob.glob(sample+'_*sgRNAList.txt')[0]
    sgRNARanking = pd.read_table(filename, sep='\t')       
    sgIDs = list(sgRNARanking['sgRNA'].values)
    genes = list(sgRNARanking['gene'].values)  
    L = len(sgIDs)
    sample_counts = list(sgRNARanking['counts'].values)
    control_counts = list(sgRNARanking['control mean'].values)
    sig = list(sgRNARanking['significant'].values)  
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
    # Plots
    # ------------------------------------------------   
    if len(sample_sig) > 100:
        tpcy = TransparencyLevel
    else:
        tpcy = 1
    print('Generating scatterplot ...')
    if not os.path.exists(PlotDir):
        os.makedirs(PlotDir)      
    os.chdir(PlotDir)   
    fig,ax = plt.subplots(figsize=(3.5,3.5))
    plt.scatter(control_rest,sample_rest,s=dotsize,facecolor='black',lw=0,alpha=TransparencyLevel,\
        rasterized=True)
    plt.scatter(control_sig,sample_sig,s=dotsize,facecolor='green',lw=0,alpha=tpcy,\
        label='significant',rasterized=True)
    if len(K_nonT)>0 and ShowNonTargets:
        plt.scatter(control_nonT,sample_nonT,s=dotsize,facecolor='orange',lw=0,alpha=0.15,\
            label='non-targeting', rasterized=True)
    if GOI != 'none':
        plt.scatter(control_goi,sample_goi,s=dotsize,facecolor='red',lw=0,alpha=1.00,label=GOI,\
            rasterized=True)
    if len(K_sig)>0:
        xmax = 1.05*max([max(control_rest),max(control_sig)])
        ymax = 1.25*max([max(sample_rest),max(sample_sig)])  
        xmin = -0.1*max([max(control_rest),max(control_sig)])
        ymin = -0.1*max([max(sample_rest),max(sample_sig)])
    else:
        xmax = 1.05*max(control_rest)
        ymax = 1.25*max(sample_rest)
        xmin = -0.1*max(control_rest)
        ymin = -0.1*max(sample_rest)
    plt.xlim([xmin,xmax]); plt.ylim([ymin,ymax])
    plt.tick_params(labelsize=11)
    axes = plt.gca()
    x0 = axes.get_xlim()  
    plt.plot((x0[0],x0[1]), (x0[0],x0[1]), ls="--", color=(51/255,153/255,1), alpha=0.75)    
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))    
    plt.xlabel('log'+str(logbase)+' counts (Control)', fontsize=11)    
    plt.ylabel('log'+str(logbase)+' counts ('+sample+')', fontsize=11)
    plt.title('sgRNA '+ScreenType.capitalize(), fontsize=12)
    leg = plt.legend(loc='upper left', prop={'size':6})
    for lh in leg.legendHandles: lh.set_alpha(1)
    if annotate:
        for label, x, y in zip(goi_sgIDs,control_goi,sample_goi):
            plt.annotate(label,xy=(x,y),color='red',fontsize=4,fontweight='bold')  
    plt.tight_layout()
    # Define file name
    figurename = 'counts_'+sample+'_scatterplot.png'
    if GOI != 'none':    
        figurename = figurename[:-4]+'_'+GOI+'.png'
    if annotate:
        figurename = figurename[:-4]+'_IDs.png'        
    if ShowNonTargets:
        figurename = figurename[:-4]+'_nonT.png'        
    # Save figure
    if GOI != 'none':        
        if not os.path.exists(PlotDir+'/'+sample+'_Highlighted_Genes/'):
        	os.makedirs(PlotDir+'/'+sample+'_Highlighted_Genes/')         
        os.chdir(PlotDir+'/'+sample+'_Highlighted_Genes/')        
        plt.savefig(figurename, dpi=res)
        if svg:
            plt.savefig(figurename[:-4]+'.svg')
        os.chdir(PlotDir)        
    else:
        plt.savefig(figurename, dpi=res)
        if svg:
            plt.savefig(figurename[:-4]+'.svg')
    plt.close()

    # ------------------------------------------------
    # Printing 
    # ------------------------------------------------
    if GOI != 'none' and PrintHighlights:
        print('-----------------------------------------------')     
        print('sgID\t\tCounts\tControl\tSignificant')    
        print('-----------------------------------------------')       
        if not K_goi:
            print('### ERROR: Gene name not found! ###')
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
    if len(sys.argv) == 2:
        input1 = sys.argv[1]
        GOI_Scatterplot(input1)            
    elif len(sys.argv) == 3:
        input1 = sys.argv[1]
        input2 = sys.argv[2]
        GOI_Scatterplot(input1,input2)    
    elif len(sys.argv) == 4:
        input1 = sys.argv[1]
        input2 = sys.argv[2]
        input3 = sys.argv[3]
        GOI_Scatterplot(input1,input2,input3)            
    elif len(sys.argv) == 5:
        input1 = sys.argv[1]
        input2 = sys.argv[2]
        input3 = sys.argv[3]
        input4 = sys.argv[4]
        GOI_Scatterplot(input1,input2,input3,input4)
    elif len(sys.argv) == 6:
        input1 = sys.argv[1]
        input2 = sys.argv[2]
        input3 = sys.argv[3]
        input4 = sys.argv[4]
        input5 = sys.argv[5]
        GOI_Scatterplot(input1,input2,input3,input4,input5)