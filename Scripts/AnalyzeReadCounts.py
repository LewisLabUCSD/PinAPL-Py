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
from matplotlib.ticker import FuncFormatter

def kilos(x, pos):
    return '%1.0fk' % (x*1e-3)

def kilos1(x, pos):
    return '%1.1fk' % (x*1e-3)

def AnalyzeCounts(sample):
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
    DataDir = config['DataDir']
    AnalysisDir = config['AnalysisDir']
    sgRNAReadCountDir = config['sgRNAReadCountDir']
    GeneReadCountDir = config['GeneReadCountDir']
    OutputDir = config['CountQCDir']+sample
    res = config['dpi']    
    svg = config['svg']
    logfilename = sample+'_ReadCount_Statistics.txt'
   
    # --------------------------------------
    # Load counts
    # --------------------------------------
    os.chdir(sgRNAReadCountDir)
    colnames = ['ID','gene','counts']
    GuideFileName = sample+'_GuideCounts_normalized.txt'
    GuideFile = pd.read_table(GuideFileName, sep='\t', names=colnames)
    ReadsPerGuide = list(GuideFile['counts'])
    L = len(ReadsPerGuide)   
    os.chdir(GeneReadCountDir)
    colnames = ['gene','counts']    
    GeneFileName = sample+'_GeneCounts_normalized.txt'
    GeneFile = pd.read_table(GeneFileName, sep='\t', names=colnames)
    ReadsPerGene = list(GeneFile['counts'])
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
    plt.figure(figsize=(3,6))
    gs = gridspec.GridSpec(2, 1)
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])   
    ax0.plot(xu,yu, linewidth=2, color='green')   
    ax0.plot(xu,xu, '--', color='#dbdcdd')
    ax0.set_ylim([0,1])    
    ax0.set_xlabel('Cumulative Fraction of sgRNAs', fontsize=10)
    ax0.set_ylabel('Cumulative Fraction of Reads', fontsize=10)
    ax0.tick_params(labelsize=11)
    ax0.set_title('Read Count Inequality (sgRNAs)', fontsize=10)
    ax0.text(.05,.8,'Gini coefficient: '+str((round(GiniIndex_u*1000)/1000)),fontsize=9)
    ax1.plot(xg,yg, linewidth=2, color='blue')   
    ax1.plot(xg,xg, '--', color='#dbdcdd')
    ax1.set_ylim([0,1])
    ax1.set_xlabel('Cumulative Fraction of Genes', fontsize=10)    
    ax1.set_ylabel('Cumulative Fraction of Reads', fontsize=10)
    ax1.tick_params(labelsize=11)
    ax1.set_title('Read Count Inequality (Genes)', fontsize=10)
    ax1.text(.05,.8,'Gini coefficient: '+str((round(GiniIndex_g*1000)/1000)),fontsize=9)    
    plt.tight_layout()
    plt.savefig(sample+'_LorenzCurves.png',dpi=res)
    if svg:
        plt.savefig(sample+'_LorenzCurves.svg')        
   
    # --------------------------------------
    # Boxplots & Histograms
    # --------------------------------------
    fig = plt.figure(figsize=(6,5))   
    gs = gridspec.GridSpec(2, 2, width_ratios=[1, 2])
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])
    ax2 = plt.subplot(gs[2])
    ax3 = plt.subplot(gs[3])    
    print('Generating boxplots...')
    # Reads per guide: Boxplot
    bp = ax0.boxplot(ReadsPerGuide, showfliers = False, patch_artist=True) # No outliers
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['medians'], color='red') 
    plt.setp(bp['whiskers'], color='black')
    for patch in bp['boxes']:
        patch.set(facecolor='#92fcae') 
    ax0.set_xticks([''])
    ax0.set_ylabel('Counts per sgRNA', fontsize=11)
    ax0.tick_params(labelsize=11)
    # Reads per gene: Boxplot
    bp = ax2.boxplot(ReadsPerGene, showfliers = False, patch_artist=True) # No outliers
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['medians'], color='red') 
    plt.setp(bp['whiskers'], color='black')
    for patch in bp['boxes']:
        patch.set(facecolor='#9de4f9') 
    ax2.set_xticks([''])
    ax2.set_ylabel('Counts per Gene', fontsize=11)
    ax2.tick_params(labelsize=11)
    print('Generating histograms...')
    # Reads per guide: Histogram
    ax1.set_title('EMPTY', color='white', fontsize=14)
    fig.text(.3,.95,'Read Distribution (sgRNAs)', fontsize=12)
    Counts_noFliers = list()
    max_count = max(40,int(numpy.percentile(ReadsPerGuide,99)))
    max_count = min(max_count,60)    
    for count in ReadsPerGuide:
        if count <= max_count:
            Counts_noFliers.append(count)
    ax1.hist(Counts_noFliers, color='green', bins=range(max_count+2), align = 'left')   
    ax1.set_xlabel('Counts per sgRNA', fontsize=11)
    ax1.set_ylabel('Number of sgRNAs', fontsize=11)    
    ax1.tick_params(labelsize=11)
    formatter = FuncFormatter(kilos)
    ax1.yaxis.set_major_formatter(formatter)
    ax1.set_xlim([-10,max_count])
    # Reads per gene: Histogram
    ax3.set_title('EMPTY', color='white', fontsize=14)
    fig.text(.3,.47,'Read Distribution (Genes)', fontsize=12)
    Counts_noFliers = list()
    max_count = max(40,int(numpy.percentile(ReadsPerGene,99)))
    max_count = min(max_count,120)
    for count in ReadsPerGene:
        if count <= max_count:
            Counts_noFliers.append(count)
    ax3.hist(Counts_noFliers, color='blue', bins=range(max_count+2), align = 'left')
    ax3.set_xlabel('Counts per Gene', fontsize=11)
    ax3.set_ylabel('Number of Genes', fontsize=11)
    ax3.tick_params(labelsize=11)
    formatter = FuncFormatter(kilos1)
    ax3.yaxis.set_major_formatter(formatter)        
    ax3.set_xlim([-5,max_count])
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
    LogFile.write('Read Counts per sgRNA (normalized)\n')    
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
    LogFile.write('Read Counts per Gene (normalized)\n')    
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
        print('Time elapsed [secs]: ' + '%.3f' % time_elapsed +'\n')
    elif sec_elapsed < 3600:
        time_elapsed = sec_elapsed/60
        print('Time elapsed [mins]: ' + '%.3f' % time_elapsed +'\n')
    else:
        time_elapsed = sec_elapsed/3600
        print('Time elapsed [hours]: ' + '%.3f' % time_elapsed +'\n')            


if __name__ == "__main__":
    input1 = sys.argv[1]
    AnalyzeCounts(input1)
