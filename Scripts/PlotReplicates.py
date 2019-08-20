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

def Repl_Scatterplot(Repl1,Repl2,GOI='none',Annot='none',NonT='none',Transp='none'):
    # ------------------------------------------------
    # Print header
    # ------------------------------------------------
    print('++++++++++++++++++++++++++++++++++++++++++++++++')
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
    sgRNAReadCountDir = config['sgRNAReadCountDir']
    PlotDir = config['CorrelDir']
    HiLiteDir2 = config['HiLiteDir2']
    delta = config['delta']
    NonTPrefix = config['NonTargetPrefix']
    res = config['dpi']
    svg = config['svg']
    dotsize = config['dotsize']
    logbase = config['logbase'] 
    ShowNonTargets = config['ShowNonTargets']
    TransparencyLevel = config['TransparencyLevel']
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
    print('Reading sgRNA counts ...')    
    os.chdir(sgRNAReadCountDir)
    colnames = ['sgRNA','gene','counts']    
    filename1 = Repl1+'_GuideCounts_normalized.txt'
    ListFile1 = pd.read_table(filename1, sep='\t',low_memory=False,names=colnames)    
    ListFile1 = ListFile1.sort_values('sgRNA')    
    filename2 = Repl2+'_GuideCounts_normalized.txt'
    ListFile2 = pd.read_table(filename2, sep='\t',low_memory=False,names=colnames)
    ListFile2 = ListFile2.sort_values('sgRNA')    
    sgIDs = list(ListFile1['sgRNA'].values)
    genes = list(ListFile1['gene'].values)  
    L = len(sgIDs)
    repl1_counts = list(ListFile1['counts'])
    repl2_counts = list(ListFile2['counts'])
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
    K_goi = [k for k in range(L) if genes[k] == GOI] 
    K_nonT = [k for k in range(L) if NonTPrefix in genes[k]]
    K_rest = list(set(range(L)) - set.union(set(K_goi),set(K_nonT))) 
    repl1_goi = [repl1_log[k] for k in K_goi]
    repl2_goi = [repl2_log[k] for k in K_goi]
    repl1_nonT = [repl1_log[k] for k in K_nonT]
    repl2_nonT = [repl2_log[k] for k in K_nonT]
    repl1_rest = [repl1_log[k] for k in K_rest]
    repl2_rest = [repl2_log[k] for k in K_rest]
    goi_sgIDs = [sgIDs[k] for k in K_goi]

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
    goi_highlight = False; nonT_highlight = False
    os.chdir(PlotDir)   
    fig,ax = plt.subplots(figsize=(3.5,3.5))
    plt.scatter(repl1_rest,repl2_rest,s=dotsize,facecolor='black',lw=0,alpha=TransparencyLevel,\
        rasterized=True)
    if len(K_nonT)>0 and ShowNonTargets:
        plt.scatter(repl1_nonT,repl2_nonT,s=dotsize,facecolor='orange',lw=0,alpha=0.15,\
            label='non-targeting',rasterized=True)
        nonT_highlight = True
    if GOI != 'none':
        plt.scatter(repl1_goi,repl2_goi,s=2*dotsize,facecolor='red',lw=0,alpha=1.00,label=GOI,\
            rasterized=True)
        goi_highlight = True
    if goi_highlight or nonT_highlight:
        leg = plt.legend(loc='upper left', prop={'size':8})
        for lh in leg.legendHandles: lh.set_alpha(1)
    axes = plt.gca()
    x0 = axes.get_xlim()  
    plt.plot((x0[0],x0[1]), (x0[0],x0[1]), ls="--", color=(51/255,153/255,1), alpha=0.75)
    plt.title('Correlation '+Repl1+' '+Repl2, fontsize=11)
    plt.xlabel(Repl1+' log'+str(logbase)+' counts', fontsize=11)    
    plt.ylabel(Repl2+' log'+str(logbase)+' counts', fontsize=11)    
    plt.text(.45,.15,'Corr (Pearson) = '+str(round(CorrCoeffP*1000)/1000),transform=axes.transAxes,\
        fontsize=7, color='blue') 
    plt.text(.45,.1,'Corr (Spearman) = '+str(round(CorrCoeffS*1000)/1000),transform=axes.transAxes,\
        fontsize=7, color='blue') 
    if len(repl1_nonT) != 0:
        xmin = -0.1*(max([max(repl1_rest),max(repl1_nonT)]))
        xmax = 1.05*(max([max(repl1_rest),max(repl1_nonT)]))
    else:        
        xmin = -0.1*max(repl1_rest)
        xmax = 1.05*max(repl1_rest)
    if len(repl2_nonT) != 0:
        ymin = -0.1*(max([max(repl2_rest),max(repl2_nonT)]))
        ymax = 1.05*(max([max(repl2_rest),max(repl2_nonT)]))    
    else:
        ymin = -0.1*max(repl2_rest)
        ymax = 1.05*max(repl2_rest)    
    plt.xlim([xmin,xmax]); plt.ylim([ymin,ymax])
    plt.tick_params(labelsize=11)
    if annotate:
        for label, x, y in zip(goi_sgIDs,repl1_goi,repl2_goi):
            plt.annotate(label,xy=(x,y),color='red',fontsize=5,fontweight='bold')     
    plt.tight_layout()  
    # Define file name
    figurename = 'counts_'+Repl1+'_'+Repl2+'.png'    
    if GOI != 'none':    
        figurename = figurename[:-4]+'_'+GOI+'.png'
    if annotate:
        figurename = figurename[:-4]+'_IDs.png'        
    if ShowNonTargets:
        figurename = figurename[:-4]+'_nonT.png'  
    # Save figure        
    if GOI != 'none':
        if not os.path.exists(HiLiteDir2):
        	os.makedirs(HiLiteDir2)         
        os.chdir(HiLiteDir2)
        plt.savefig(figurename, dpi=res)
        if svg:
            plt.savefig(figurename[:-4]+'.svg')
        os.chdir(PlotDir) 
    else:
        plt.savefig(figurename, dpi=res)
        if svg:
            plt.savefig(figurename[:-4]+'.svg')
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
    if len(sys.argv) == 3:
        input1 = sys.argv[1]
        input2 = sys.argv[2]
        Repl_Scatterplot(input1,input2)    
    elif len(sys.argv) == 4:
        input1 = sys.argv[1]
        input2 = sys.argv[2]
        input3 = sys.argv[3]
        Repl_Scatterplot(input1,input2,input3)
    elif len(sys.argv) == 5:
        input1 = sys.argv[1]
        input2 = sys.argv[2]
        input3 = sys.argv[3]
        input4 = sys.argv[4]
        Repl_Scatterplot(input1,input2,input3,input4)
    elif len(sys.argv) == 6:
        input1 = sys.argv[1]
        input2 = sys.argv[2]
        input3 = sys.argv[3]
        input4 = sys.argv[4]
        input5 = sys.argv[5]        
        Repl_Scatterplot(input1,input2,input3,input4,input5)
    elif len(sys.argv) == 7:
        input1 = sys.argv[1]
        input2 = sys.argv[2]
        input3 = sys.argv[3]
        input4 = sys.argv[4]
        input5 = sys.argv[5]        
        input6 = sys.argv[6] 
        Repl_Scatterplot(input1,input2,input3,input4,input5,input6)