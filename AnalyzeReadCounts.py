#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 13:13:27 2016

@author: philipp
"""

# Analyze count distribution
# =======================================================================
# Imports
from __future__ import division # floating point division by default
import glob
import pandas as pd
from Lorenz import gini
import yaml
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy
import os
import sys
import time
import warnings
warnings.simplefilter(action = "ignore", category = FutureWarning)


def AnalyzeCounts(sample):
    # ------------------------------------------------
    # Print header
    # ------------------------------------------------
    print('++++++++++++++++++++++++++++++++++++++')
    print('PinAPL-Py: Exploratory Counts Analysis')
    print('++++++++++++++++++++++++++++++++++++++')    
    start_total = time.time()   

    # ------------------------------------------------
    # Get parameters
    # ------------------------------------------------
    configFile = open('../configuration.yaml','r')
    config = yaml.load(configFile)
    configFile.close()
    ScriptsDir = config['ScriptsDir']    
    DataDir = config['DataDir']
    AnalysisDir = config['AnalysisDir']
    ScreenType = config['ScreenType']
    max_q = config['max_q']
    InputDir = config['AlnQCDir']+sample
    OutputDir = config['CountQCDir']+sample
    res = config['dpi']    
    svg = config['svg']
    logfilename = sample+'_ReadStatistics.txt'
   
    # --------------------------------------
    # Load counts
    # --------------------------------------
    os.chdir(InputDir)
    colnames = ['ID','gene','counts']
    GuideFileName = glob.glob('*_GuideCounts_0.tsv')[0]
    GuideFile = pd.read_table(GuideFileName, sep='\t', names=colnames)
    ReadsPerGuide = list(GuideFile['counts'].values)
    L = len(ReadsPerGuide)   
    colnames = ['gene','counts']    
    GeneFileName = glob.glob('*_GeneCounts_0.tsv')[0]
    GeneFile = pd.read_table(GeneFileName, sep='\t', names=colnames)
    ReadsPerGene = list(GeneFile['counts'].values)
    sgID = list(GuideFile['ID'].values)    
    gene = list(GuideFile['gene'].values)    
    if not os.path.exists(OutputDir):
        os.makedirs(OutputDir)
    os.chdir(OutputDir)

    # --------------------------------------
    # Lorenz Curve
    # --------------------------------------
    print('Computing Gini coefficients ... ')    
    GiniIndex_u,xu,yu = gini(ReadsPerGuide)
    GiniIndex_g,xg,yg = gini(ReadsPerGene)
    GiniIndex_u = round(GiniIndex_u*1000)/1000
    GiniIndex_g = round(GiniIndex_g*1000)/1000
    print('Gini Index (sgRNAs): ' + str(round(GiniIndex_u*1000)/1000))
    print('Gini Index (genes): ' + str(round(GiniIndex_g*1000)/1000))
    # Plot Lorenz curves
    print('Plotting Lorenz curves ...')
    plt.figure(figsize=(3.5,7))
    gs = gridspec.GridSpec(2, 1)
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])   
    ax0.plot(xu,yu, linewidth=2)   
    ax0.plot(xu,xu, 'r--')
    ax0.set_ylim([0,1])    
    ax0.set_xlabel('Cumulative Fraction of sgRNAs', fontsize=10)    
    ax0.set_ylabel('Cumulative Fraction of Reads', fontsize=10)
    ax0.set_title('Read Disparity (sgRNAs)', fontsize=14)
    ax0.text(.1,.8,'Gini coefficient: '+str((round(GiniIndex_u*1000)/1000)),fontsize=9)     
    ax1.plot(xg,yg, linewidth=2)   
    ax1.plot(xg,xg, 'r--')
    ax1.set_ylim([0,1])
    ax1.set_xlabel('Cumulative Fraction of Genes', fontsize=10)    
    ax1.set_ylabel('Cumulative Fraction of Reads', fontsize=10)
    ax1.set_title('Read Disparity (Genes)', fontsize=14)     
    ax1.text(.1,.8,'Gini coefficient: '+str((round(GiniIndex_g*1000)/1000)),fontsize=9)    
    plt.tight_layout()
    plt.savefig(sample+'_LorenzCurves.png',dpi=res)
    if svg:
        plt.savefig(sample+'_LorenzCurves.svg')        
   
    # --------------------------------------
    # Boxplots & Histograms
    # --------------------------------------
    plt.figure(figsize=(7,5))
    gs = gridspec.GridSpec(2, 2, width_ratios=[1, 2])
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])
    ax2 = plt.subplot(gs[2])
    ax3 = plt.subplot(gs[3])    
    print('Generating boxplots...')
    # Reads per guide: Boxplot
    bp = ax0.boxplot(ReadsPerGuide, showfliers = False) # No outliers
    plt.setp(bp['boxes'], color='black')
    ax0.set_title('Reads per sgRNA', fontsize=12)
    ax0.set_xticks([''])
    ax0.set_ylabel('counts [norm.]', fontsize=10)
    # Reads per gene: Boxplot
    bp = ax2.boxplot(ReadsPerGene, showfliers = False) # No outliers
    plt.setp(bp['boxes'], color='black')
    ax2.set_title('Reads per Gene', fontsize=12)
    ax2.set_xticks([''])
    ax2.set_ylabel('counts [norm.]', fontsize=10)
    print('Generating histograms...')
    # Reads per guide: Histogram
    Counts_noFliers = list()
    max_count = int(numpy.percentile(ReadsPerGuide,max_q))
    for count in ReadsPerGuide:
        if count <= max_count:
            Counts_noFliers.append(count)
    ax1.hist(Counts_noFliers, bins = range(max_count+2), align = 'left')
    ax1.set_title(sample+' Read Distribution', fontsize=12)
    ax1.set_xlabel('counts [norm.] per sgRNA', fontsize=10)
    ax1.set_ylabel('Number of sgRNAs', fontsize=10)
    # Reads per gene: Histogram
    Counts_noFliers = list()
    max_count = int(numpy.percentile(ReadsPerGene,max_q))
    for count in ReadsPerGene:
        if count <= max_count:
            Counts_noFliers.append(count)
    ax3.hist(Counts_noFliers, bins = range(max_count+2), align = 'left')
    ax3.set_title(sample+' Read Distribution', fontsize=12)
    ax3.set_xlabel('counts [norm.] per Gene', fontsize=10)
    ax3.set_ylabel('Number of Genes', fontsize=10)
    plt.tight_layout()
    plt.savefig(sample+'_ReadCount_Distribution.png',dpi=res)
    if svg:
        plt.savefig(sample+'_ReadCount_Distribution.svg')        
    
    # --------------------------------------------
    # Counts distribution
    # --------------------------------------------
    print('Writing distribution file ... ')
    N_Guides = len(ReadsPerGuide)
    N_Genes = len(ReadsPerGene)
    guide_m = int(numpy.median(ReadsPerGuide))
    guide_sd = int(numpy.std(ReadsPerGuide))
    guide_q25 = int(numpy.percentile(ReadsPerGuide,25))
    guide_q75 = int(numpy.percentile(ReadsPerGuide,75))
    guide_min = int(min(ReadsPerGuide))    
    guide_max = int(max(ReadsPerGuide))
    guide_pres = [1 for n in ReadsPerGuide if n > 0]
    guide_pres = len(guide_pres)
    guide_pres100 = round((guide_pres/N_Guides)*1000)/10
    gene_m = int(numpy.median(ReadsPerGene))
    gene_sd = int(numpy.std(ReadsPerGene))
    gene_q25 = int(numpy.percentile(ReadsPerGene,25))
    gene_q75 = int(numpy.percentile(ReadsPerGene,75))
    gene_min = int(min(ReadsPerGene))    
    gene_max = int(max(ReadsPerGene))
    gene_pres = [1 for n in ReadsPerGene if n > 0]
    gene_pres = len(gene_pres)
    gene_pres100 = round((gene_pres/N_Genes)*1000)/10
    # Write log file
    LogFile = open(logfilename,'w')     
    LogFile.write(sample+' Read Counts Distribution:\n')
    LogFile.write('***********************************\n')
    LogFile.write('\n')
    LogFile.write('Read Counts per sgRNA [norm.]\n')    
    LogFile.write('------------------------------------\n')
    LogFile.write('Median:\t\t\t'+str(guide_m)+'\n')
    LogFile.write('Standard Deviation:\t'+str(guide_sd)+'\n')    
    LogFile.write('25% Quantile:\t\t'+str(guide_q25)+'\n')
    LogFile.write('75% Quantile:\t\t'+str(guide_q75)+'\n')
    LogFile.write('Minimum:\t\t'+str(guide_min)+'\n')
    LogFile.write('Maximum:\t\t'+str(guide_max)+'\n')
    LogFile.write('sgRNA Representation:\t'+str(guide_pres)+' ('+str(guide_pres100)+'%)\n')
    LogFile.write('Gini coefficient:\t'+str(GiniIndex_u)+'\n')
    LogFile.write('\n')
    LogFile.write('Read Counts per Gene [norm.]\n')    
    LogFile.write('------------------------------------\n')
    LogFile.write('Median:\t\t\t'+str(gene_m)+'\n')
    LogFile.write('Standard Deviation:\t'+str(gene_sd)+'\n')    
    LogFile.write('25% Quantile:\t\t'+str(gene_q25)+'\n')
    LogFile.write('75% Quantile:\t\t'+str(gene_q75)+'\n')
    LogFile.write('Minimum:\t\t'+str(gene_min)+'\n')
    LogFile.write('Maximum:\t\t'+str(gene_max)+'\n')
    LogFile.write('Gene Representation:\t'+str(gene_pres)+' ('+str(gene_pres100)+'%)\n')
    LogFile.write('Gini coefficient:\t'+str(GiniIndex_g)+'\n')    
    LogFile.close()

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
    AnalyzeCounts(input1)
