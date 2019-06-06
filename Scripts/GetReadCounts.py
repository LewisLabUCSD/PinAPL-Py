#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 09:57:05 2018

@author: philipp
"""

# Perform Read counting, classification and normalization
# =======================================================================
# Imports
from __future__ import division # floating point division by default
import sys
import yaml
import time
import os
from os import path
import pandas
from collections import Counter
from mpl_toolkits.mplot3d import axes3d
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import numpy
import pysam
from matplotlib.ticker import FuncFormatter



def millions(x, pos):
    return '%1.1fM' % (x*1e-6)

def millions2(x, pos):
    return '%1.2fM' % (x*1e-6)    


def CountReads(sample):  
    # ------------------------------------------------
    # Print header
    # ------------------------------------------------
    print('++++++++++++++++++++++++++++++++++++++++++++++++')  
    start = time.time()

    # ------------------------------------------------
    # Get parameters
    # ------------------------------------------------
    configFile = open('configuration.yaml','r')
    config = yaml.load(configFile)
    configFile.close()
    WorkingDir = config['WorkingDir']    
    ScriptsDir = config['ScriptsDir']
    LibDir = config['LibDir']    
    LibFilename = config['LibFilename']
    LibFormat = LibFilename[-3:]
    if LibFormat == 'tsv':
        libsep = '\t'
    elif LibFormat == 'csv':
        libsep = ','    
    OutputDir = config['AlnQCDir']+sample   
    AlnStemDir = config['AlignDir']
    AlnDir = AlnStemDir+sample+'/'    
    minN = config['Cutoff']
    Theta = config['Theta']
    AS_min = config['AS_min']    
    N0 = 1000000
    res = config['dpi']
    svg = config['svg']
    AlnOutput = config['AlnOutput']
    keepCutReads = config['keepCutReads']
    GuideCount_Suffix = '_GuideCounts.txt'
    logfilename = sample+'_AlignmentResults.txt'
    L_bw = config['L_bw']
    N_bw = config['N_bw']
    i_bw = config['i_bw']    
    ReadCountDir = config['sgRNAReadCountDir']

    # ------------------------------------------------
    # Read library
    # ------------------------------------------------  
    os.chdir(LibDir)
    LibCols = ['gene','ID','seq']
    LibFile = pandas.read_table(LibFilename, sep = libsep, skiprows = 1, names = LibCols)
    LibFile = LibFile.sort_values(['gene','ID'])    
    sgIDs = list(LibFile['ID'])
    global L; L = len(sgIDs)
    global geneIDs; geneIDs = list(LibFile['gene'])
    G = len(set(geneIDs))

    # ------------------------------------------------
    # Get sample read file
    # ------------------------------------------------      
    os.chdir(WorkingDir)
    DataSheet = pandas.read_excel('DataSheet.xlsx')
    FileNames = list(DataSheet['FILENAME'].values)
    n = len(FileNames)
    Samples = list(DataSheet['SAMPLE NAME'].values)
    ReadsFilename = [FileNames[j] for j in range(n) if Samples[j] == sample][0] 
    ReadsFilename0 = 'Trim_'+ReadsFilename          
    
    # ------------------------------------------
    # Extract and analyze alignments
    # ------------------------------------------     
    print('Loading alignment ...')      
    # CLASSIFY ALIGNMENTS 
    os.chdir(AlnDir)
    sam_file = ReadsFilename0 + '_bw2output.sam' 
    sam_file_present = path.exists(sam_file)
    bam_file = ReadsFilename0 + '_bw2output.bam'
    bam_file_present = path.exists(bam_file)    
    if sam_file_present:
        sam_file = sam_file
    elif bam_file_present: 
        os.system('samtools view -h '+bam_file+' > '+sam_file)
    else:
        print('### ERROR: No alignment file present ###')
    bw2sam = pysam.AlignmentFile(sam_file,'rb')
    print('Applying matching threshold ...')
    print('Applying ambiguity threshold ...')       
    NFail = 0; NUnique = 0; NTol = 0; NAmb = 0
    mapQ = list()
    primScore = list()
    secScore = list()
    AlnStatus = list()
    sgRNA_Hitlist = list()
    for read in bw2sam.fetch():
        mapQ.append(read.mapping_quality)
        # read with primary and seconday alignment
        if read.has_tag('AS'):
            AS = read.get_tag('AS')            
            if AS >= AS_min:                
                if read.has_tag('XS'):
                    XS = read.get_tag('XS')       
                    primScore.append(AS)
                    secScore.append(XS)            
                    if XS <= AS - Theta:
                        NTol += 1
                        AlnStatus.append('Tolerate')
                        sgRNA_Hitlist.append(read.reference_name)
                    else:
                        NAmb += 1  
                        AlnStatus.append('Ambiguous')
                # read with only primary alignment
                else:
                    primScore.append(AS)
                    secScore.append(0)
                    NUnique += 1
                    AlnStatus.append('Unique')
                    sgRNA_Hitlist.append(read.reference_name)                
            else:
                # read with insufficient primary score 
                primScore.append(AS)                
                if read.has_tag('XS'):
                    XS = read.get_tag('XS')
                    secScore.append(XS)
                else:
                    secScore.append(0)          
                NFail += 1
                AlnStatus.append('Fail')
        # read with failed alignment
        else:
            primScore.append(0)
            secScore.append(0)
            NFail += 1
            AlnStatus.append('Fail')
    bw2sam.close();          
    NReads = NTol + NAmb + NUnique + NFail
    FracUnique = round(NUnique/NReads*1000)/10    
    FracTol = round(NTol/NReads*1000)/10
    FracAmb = round(NAmb/NReads*1000)/10
    FracFail = round(NFail/NReads*1000)/10         
    print('*** Successfully mapped reads: '\
        +str(NUnique+NTol)+' ('+str(FracUnique+FracTol)+'%) ***')
        
    # ------------------------------------------
    # Mapping statistics
    # ------------------------------------------ 
    print('Writing alignment logfile ...')  
    if not os.path.exists(OutputDir):
        os.makedirs(OutputDir)        
    os.chdir(OutputDir)    
    LogFile = open(logfilename,'w')         
    LogFile.write(sample+' Alignment Results\n')
    LogFile.write('**************************************\n')
    LogFile.write('\n')
    LogFile.write('Total Number of Reads: \t'+str(NReads)+'\n')    
    LogFile.write('\n')    
    LogFile.write('Number of Reads with unique Alignments: \t'+str(NUnique)+' ('+str(FracUnique)+'%)\n')
    LogFile.write('Number of Reads above Ambiguity Tolerance: \t'+str(NTol)+' ('+str(FracTol)+'%)\n')
    LogFile.write('---------------------------------------------------------------\n')
    LogFile.write('Total Number of Reads matched: \t\t\t'+str(NUnique+NTol)+' ('+str(FracUnique+FracTol)+'%)\n')
    LogFile.write('\n\n')    
    LogFile.write('Number of Reads below Ambiguity Tolerance: \t'+str(NAmb)+' ('+str(FracAmb)+'%)\n')
    LogFile.write('Number of Reads with failed Alignment: \t\t'+str(NFail)+' ('+str(FracFail)+'%)\n')
    LogFile.write('---------------------------------------------------------------\n')
    LogFile.write('Total Number of Reads discarded: \t\t'+str(NFail+NAmb)+' ('+str(FracFail+FracAmb)+'%)\n')    
    LogFile.write('\n\n\n')    
    LogFile.write('Parameter settings\n')
    LogFile.write('L_bw (Seed Length): '+str(L_bw)+'\n')    
    LogFile.write('N_bw (Seed Mismatch): '+str(N_bw)+'\n')    
    LogFile.write('i_bw (Interval Function): '+str(i_bw)+'\n')        
    LogFile.write('Theta (Ambiguity Threshold): '+str(Theta)+'\n')    
    LogFile.write('AS_min (Matching Threshold): '+str(AS_min)+'\n')        
    LogFile.close()              

    # ------------------------------------------
    # MAPPING QUALITY HISTOGRAM
    # ------------------------------------------
    os.chdir(OutputDir)
    print('Plotting mapping quality ...')
    maxQuality = max(mapQ)
    fig, ax = plt.subplots(figsize=(3.5,2.9))
    plt.hist(mapQ, bins = range(maxQuality+1), width = 1, align = 'left')
    plt.xlim([-1,maxQuality+1])
    plt.title('Read Alignment Quality', fontsize=12)
    plt.xlabel('Bowtie2 Mapping Quality', fontsize=10)
    plt.ylabel('Number of Reads', fontsize=10)
    plt.tick_params(labelsize=10)
    formatter = FuncFormatter(millions2)
    ax.yaxis.set_major_formatter(formatter) 
    plt.tight_layout()
    plt.savefig(sample+'_MappingQuality.png',dpi=res)  
    if svg:
        plt.savefig(sample+'_MappingQuality.svg')  
        
    # ------------------------------------------        
    # ALIGNMENT SCORE BARPLOT
    # ------------------------------------------        
    print('Plotting alignment scores ...')
    primScoreKeep = [primScore[k] for k in range(NReads) if AlnStatus[k] in ['Unique','Tolerate']]
    secScoreKeep = [secScore[k] for k in range(NReads) if AlnStatus[k] in ['Unique','Tolerate']]    
    primScoreToss = [primScore[k] for k in range(NReads) if AlnStatus[k] in ['Fail','Ambiguous']]
    secScoreToss = [secScore[k] for k in range(NReads) if AlnStatus[k] in ['Fail','Ambiguous']]
    # Plot bars (matched reads)
    xedges = range(42)
    yedges = range(42)
    H, xedges, yedges = numpy.histogram2d(secScoreKeep, primScoreKeep, bins=(xedges, yedges))
    Y, X = numpy.nonzero(H)
    Z = H[Y,X]
    Z_off = numpy.zeros(len(Y))
    dX = numpy.ones(len(X))
    dY = numpy.ones(len(Y))
    # Plot bars (discarded reads)
    h, xedges, yedges = numpy.histogram2d(secScoreToss, primScoreToss, bins=(xedges, yedges))
    y, x = numpy.nonzero(h)
    z = h[y,x]
    z_off = numpy.zeros(len(y))
    dx = numpy.ones(len(x))
    dy = numpy.ones(len(y))
    # Show plot    
    fig = plt.figure(figsize=(5,4.2))
    ax = fig.gca(projection='3d')
    ax.bar3d(X,Y,Z_off,dX,dY,Z, color=(97/255, 252/255, 80/255))  # bar3d has a bug with hex code
    green_proxy = plt.Rectangle((0, 0), 1, 1, fc=(97/255, 252/255, 80/255))
    ax.bar3d(x,y,z_off,dx,dy,z,color=(255/255, 51/255, 51/255))
    red_proxy = plt.Rectangle((0, 0), 1, 1, fc=(255/255, 51/255, 51/255))
    ax.set_title('Alignment Scores',fontsize=14)
    ax.legend([green_proxy,red_proxy],['Reads Accepted','Reads Discarded'],loc='upper left',prop={'size':7})
    ax.set_xticklabels([0,10,20,30,40],fontsize=8)
    ax.set_yticklabels([0,10,20,30,40],fontsize=8)  
    ax.set_xlabel('Prim. Alignment Score', fontsize=8)
    ax.set_ylabel('Sec. Alignment Score', fontsize=8)    
    formatter = FuncFormatter(millions)
    ax.zaxis.set_major_formatter(formatter)    
    ax.tick_params(axis="z", labelsize=8)
    ax.set_zlabel('Number of Reads',fontsize=8)
    plt.tight_layout()    
    plt.savefig(sample+'_AlignmentScores.png',dpi=res)   
    if svg:
        plt.savefig(sample+'_AlignmentScores.svg',dpi=res)
 
    # --------------------------------------
    # Read count acquisition
    # -------------------------------------- 
    # Read counts per sgRNA
    ReadCounts = Counter()
    for sgRNA in sgRNA_Hitlist:
        ReadCounts[sgRNA] += 1
    global ReadsPerGuide
    ReadsPerGuide = list()
    for sgRNA in sgIDs:
        ReadsPerGuide.append(ReadCounts[sgRNA])      
    # Write read counts 
    if not os.path.exists(ReadCountDir):
        os.makedirs(ReadCountDir)
    os.chdir(ReadCountDir)         
    GuideCountsFilename = sample + GuideCount_Suffix
    GuideCounts = open(GuideCountsFilename,'w')
    for k in range(L):
        GuideCounts.write(str(sgIDs[k]) + '\t'+ str(geneIDs[k]) + '\t' + str(ReadsPerGuide[k]) + '\n')
    GuideCounts.close()     
    # No-mapping warning
    if sum(ReadsPerGuide) == 0:
        print('### ERROR: Zero read counts! Check library and alignment ###')    
    end = time.time()
    # Time stamp    
    print('------------------------------------------------')
    print('Script completed.')      
    sec_elapsed = end - start
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
    CountReads(input1)