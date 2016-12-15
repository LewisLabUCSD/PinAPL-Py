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
    print('****************************************************************')
    print('PinAPL-Py: Exploratory Counts Analysis')
    print('****************************************************************')    
    start_total = time.time()   

    # ------------------------------------------------
    # Get parameters
    # ------------------------------------------------
    os.chdir('/workingdir/')
    configFile = open('configuration.yaml','r')
    config = yaml.load(configFile)
    configFile.close()
    DataDir = config['DataDir']
    AnalysisDir = config['AnalysisDir']
    ScreenType = config['ScreenType']
    max_q = config['max_q']
    QCDir = config['QCDir']
    res = config['dpi']
    LogDir = QCDir+sample+'/' 
    logfilename = 'ReadCounts_Statistics.txt'
   
    # --------------------------------------
    # Load counts
    # --------------------------------------
    os.chdir(LogDir)
    colnames = ['ID','gene','counts']
    GuideFileName = glob.glob(sample+'*'+'_GuideCounts_0.tsv')[0]
    GuideFile = pd.read_table(GuideFileName, sep='\t', names=colnames)
    ReadsPerGuide = list(GuideFile['counts'].values)
    L = len(ReadsPerGuide)   
    colnames = ['gene','counts']    
    GeneFileName = glob.glob(sample+'*'+'_GeneCounts_0.tsv')[0]
    GeneFile = pd.read_table(GeneFileName, sep='\t', names=colnames)
    ReadsPerGene = list(GeneFile['counts'].values)
    sgID = list(GuideFile['ID'].values)    
    gene = list(GuideFile['gene'].values)    
 

    # --------------------------------------
    # Rank counts
    # --------------------------------------
    print('Ranking sgRNA counts ...')    
    Results_df = pd.DataFrame(data = {'sgRNA': [sgID[k] for k in range(L)],
                                     'gene': [gene[k] for k in range(L)],
                                     'counts': [ReadsPerGuide[k] for k in range(L)]},
                                    columns = ['sgRNA','gene','counts'])
    if ScreenType == 'enrichment':
        Results_df_0 = Results_df.sort_values(['counts'],ascending=[0])
    else:
        Results_df_0 = Results_df.sort_values(['counts'],ascending=[1])
    Results_Filename = sample+'_GuideCounts_0_sorted.tsv'
    Results_df_0.to_csv(Results_Filename,sep='\t',index=False)
    
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
    plt.figure(figsize=(5,10))
    gs = gridspec.GridSpec(2, 1)
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])   
    ax0.plot(xu,yu, linewidth=2)   
    ax0.plot(xu,xu, 'r--')
    ax0.set_ylim([0,1])    
    ax0.set_xlabel('Cumulative Fraction of sgRNAs', fontsize=12)    
    ax0.set_ylabel('Cumulative Fraction of Reads', fontsize=12)
    ax0.set_title('Read Disparity (sgRNAs)', fontsize=12, fontweight='bold')
    ax0.text(.1,.8,'Gini coefficient: '+str((round(GiniIndex_u*1000)/1000)),fontsize=14)     
    ax1.plot(xg,yg, linewidth=2)   
    ax1.plot(xg,xg, 'r--')
    ax1.set_ylim([0,1])
    ax1.set_xlabel('Cumulative Fraction of Genes', fontsize=12)    
    ax1.set_ylabel('Cumulative Fraction of Reads', fontsize=12)
    ax1.set_title('Read Disparity (Genes)', fontsize=12, fontweight='bold')     
    ax1.text(.1,.8,'Gini coefficient: '+str((round(GiniIndex_g*1000)/1000)),fontsize=14)    
    plt.savefig('ReadCounts_LorenzCurves.png',dpi=res)
   
    # --------------------------------------
    # Boxplots & Histograms
    # --------------------------------------
    plt.figure()
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
    ax0.set_ylabel('Counts', fontsize=12)
    # Reads per gene: Boxplot
    bp = ax2.boxplot(ReadsPerGene, showfliers = False) # No outliers
    plt.setp(bp['boxes'], color='black')
    ax2.set_title('Reads per Gene', fontsize=12)
    ax2.set_xticks([''])
    ax2.set_ylabel('Counts', fontsize=12)
    print('Generating histograms...')
    # Reads per guide: Histogram
    Counts_noFliers = list()
    max_count = int(numpy.percentile(ReadsPerGuide,max_q))
    for count in ReadsPerGuide:
        if count <= max_count:
            Counts_noFliers.append(count)
    ax1.hist(Counts_noFliers, bins = range(max_count+2), align = 'left')
    ax1.set_title('Read Distribution', fontsize=12)
    ax1.set_xlabel('Reads per sgRNA', fontsize=12)
    ax1.set_ylabel('Number of sgRNAs', fontsize=12)
    # Reads per gene: Histogram
    Counts_noFliers = list()
    max_count = int(numpy.percentile(ReadsPerGene,max_q))
    for count in ReadsPerGene:
        if count <= max_count:
            Counts_noFliers.append(count)
    ax3.hist(Counts_noFliers, bins = range(max_count+2), align = 'left')
    ax3.set_title('Read Distribution', fontsize=12)
    ax3.set_xlabel('Reads per Gene', fontsize=12)
    ax3.set_ylabel('Number of Genes', fontsize=12)
    plt.tight_layout()
    plt.savefig('ReadCounts_Distribution.png',dpi=res)
    
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
    guide_miss = [1 for n in ReadsPerGuide if n == 0]
    guide_miss = len(guide_miss)
    guide_miss100 = round((guide_miss/N_Guides)*1000)/10
    gene_m = int(numpy.median(ReadsPerGene))
    gene_sd = int(numpy.std(ReadsPerGene))
    gene_q25 = int(numpy.percentile(ReadsPerGene,25))
    gene_q75 = int(numpy.percentile(ReadsPerGene,75))
    gene_min = int(min(ReadsPerGene))    
    gene_max = int(max(ReadsPerGene))
    gene_miss = [1 for n in ReadsPerGene if n == 0]
    gene_miss = len(gene_miss)
    gene_miss100 = round((gene_miss/N_Genes)*1000)/10
    # Write log file
    LogFile = open(logfilename,'w')     
    LogFile.write('Read Counts Distribution:\n')
    LogFile.write('***********************************\n')
    LogFile.write('\n')
    LogFile.write('Read Counts per sgRNA\n')    
    LogFile.write('------------------------------------\n')
    LogFile.write('Median:\t\t\t'+str(guide_m)+'\n')
    LogFile.write('Standard Deviation:\t'+str(guide_sd)+'\n')    
    LogFile.write('25% Quantile:\t\t'+str(guide_q25)+'\n')
    LogFile.write('75% Quantile:\t\t'+str(guide_q75)+'\n')
    LogFile.write('Minimum:\t\t'+str(guide_min)+'\n')
    LogFile.write('Maximum:\t\t'+str(guide_max)+'\n')
    LogFile.write('Missing sgRNAs:\t\t'+str(guide_miss)+' ('+str(guide_miss100)+'%)\n')
    LogFile.write('Gini coefficient:\t'+str(GiniIndex_u)+'\n')
    LogFile.write('\n')
    LogFile.write('Read Counts per Gene\n')    
    LogFile.write('------------------------------------\n')
    LogFile.write('Median:\t\t\t'+str(gene_m)+'\n')
    LogFile.write('Standard Deviation:\t'+str(gene_sd)+'\n')    
    LogFile.write('25% Quantile:\t\t'+str(gene_q25)+'\n')
    LogFile.write('75% Quantile:\t\t'+str(gene_q75)+'\n')
    LogFile.write('Minimum:\t\t'+str(gene_min)+'\n')
    LogFile.write('Maximum:\t\t'+str(gene_max)+'\n')
    LogFile.write('Missing genes:\t\t'+str(gene_miss)+' ('+str(gene_miss100)+'%)\n')
    LogFile.write('Gini coefficient:\t'+str(GiniIndex_g)+'\n')    
    LogFile.close()

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
    AnalyzeCounts(input1)
