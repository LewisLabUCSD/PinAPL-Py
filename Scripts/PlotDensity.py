#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 16:29:08 2019

@author: philipp
"""

# Density plot of sgRNA counts
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
import seaborn


def DensityPlot(sample):
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
    DensityDir = config['DensityDir']
    ScreenType = config['ScreenType']        
    delta = config['delta']
    res = config['dpi']
    svg = config['svg']
    logbase = config['logbase']
    
    # ------------------------------------------------
    # Reading counts from sample and control
    # ------------------------------------------------
    print('Reading sgRNA read counts ...')    
    os.chdir(sgRNARanksDir)
    filename = glob.glob(sample+'_*sgRNAList.txt')[0]
    HitList = pd.read_table(filename, sep='\t')       
    sgIDs = list(HitList['sgRNA'].values)
    genes = list(HitList['gene'].values)  
    L = len(sgIDs)
    sample_counts = list(HitList['counts'].values)
    control_counts = list(HitList['control mean'].values)
    sig = list(HitList['significant'].values)  
    # Log transformation
    print('Log'+str(logbase)+' transformation ...')
    if logbase == 2:
        sample_log = [np.log2(sample_counts[k]+delta) for k in range(len(sample_counts))]
        control_log = [np.log2(control_counts[k]+delta) for k in range(len(control_counts))]
    elif logbase == 10:
        sample_log = [np.log10(sample_counts[k]+delta) for k in range(len(sample_counts))]
        control_log = [np.log10(control_counts[k]+delta) for k in range(len(control_counts))]    
    
    # --------------
    # DENSITY PLOT  
    # --------------
    print('Generating density plot ...')
    if not os.path.exists(DensityDir):
        os.makedirs(DensityDir)    
    os.chdir(DensityDir)   
    fig,ax = plt.subplots(figsize=(3.5,3.7)) 
    seaborn.set_style("white")
    seaborn.kdeplot(control_log,sample_log,cmap='Reds',shade=True,bw=.15,shade_lowest=False)
    xmax = 1.05*max(control_log)
    ymax = 1.25*max(sample_log)
    xmin = -0.1*max(control_log)
    ymin = -0.1*max(sample_log)    
    plt.xlim([xmin,xmax]); plt.ylim([ymin,ymax])
    plt.tick_params(reset = True,labelsize=11)
    axes = plt.gca()
    x0 = axes.get_xlim()  
    plt.plot((x0[0],x0[1]), (x0[0],x0[1]), ls="--", color='#808080', alpha=0.75)    
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))    
    plt.xlabel('log'+str(logbase)+' counts (Control)', fontsize=11)    
    plt.ylabel('log'+str(logbase)+' counts ('+sample+')', fontsize=11)    
    plt.title('sgRNA Density', fontsize=12)
    plt.tight_layout()
    figurename = 'counts_'+sample+'_density.png'
    plt.savefig(figurename, dpi=res)
    if svg:
        plt.savefig(figurename[:-4]+'.svg')
    plt.close()    

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
    input1 = sys.argv[1]
    DensityPlot(input1) 