#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon May 16 19:11:28 2016

@author: philipp
"""

# Perform alignment, count & normalize reads
# =======================================================================
# Imports
from __future__ import division # floating point division by default
import pandas
from Bowtie2 import RunBowtie2
import yaml
import os
import glob
import time
import sys
from collections import Counter
from joblib import Parallel, delayed
import multiprocessing
from mpl_toolkits.mplot3d import axes3d
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import numpy
import pysam
from matplotlib.ticker import FuncFormatter

def kilos(x, pos):
    return '%1.0fk' % (x*1e-3)

def CountReadsPerGene(g):
    global GeneList
    global geneIDs
    global ReadsPerGuide
    gene = GeneList[g]        
    geneIndex = [i for i,x in enumerate(geneIDs) if x==gene]
    sgCounts = [ReadsPerGuide[i] for i in geneIndex]
    g_counts = sum(sgCounts)
    return g_counts  


def MapAndCount(sample):    
    # ------------------------------------------------
    # Print header
    # ------------------------------------------------
    print('++++++++++++++++++++++++++++++++++++')
    print('PinAPL-Py: Alignment & Read Counting')
    print('++++++++++++++++++++++++++++++++++++')  
    start_total = time.time()   

    # ------------------------------------------------
    # Get parameters
    # ------------------------------------------------
    configFile = open('configuration.yaml','r')
    config = yaml.load(configFile)
    configFile.close()
    ScriptsDir = config['ScriptsDir']
    WorkingDir = config['WorkingDir']
    DataDir = config['DataDir']
    AnalysisDir = config['AnalysisDir']
    CutAdaptDir = config['CutAdaptDir']
    bw2Dir = config['bw2Dir']
    IndexDir = config['IndexDir']
    LibDir = config['LibDir']
    AlnStemDir = config['AlignDir']
    AlnDir = AlnStemDir+sample+'/'
    OutputDir = config['AlnQCDir']+sample    
    seq_5_end = config['seq_5_end']
    seq_3_end = config['seq_3_end']
    CutErrorTol = config['CutErrorTol']
    R_min = config['R_min']
    minN = config['Cutoff']
    LibFilename = config['LibFilename']
    LibFormat = LibFilename[-3:]
    if LibFormat == 'tsv':
        libsep = '\t'
    elif LibFormat == 'csv':
        libsep = ','
    Theta = config['Theta']
    L_bw = config['L_bw']
    N_bw = config['N_bw']
    i_bw = config['i_bw']
    N0 = 1000000
    res = config['dpi']
    svg = config['svg']
    AlnOutput = config['AlnOutput']
    keepCutReads = config['keepCutReads']
    AlnFileSuffix = '_bw2Aln.tsv'
    GuideCount_Suffix = '_GuideCounts.tsv'
    GeneCount_Suffix = '_GeneCounts.tsv'
    cutadaptLog = sample+'_cutadapt_log.txt'
    logfilename = sample+'_AlignmentResults.txt'
    
    
    # ------------------------------------------------
    # Get sample read file
    # ------------------------------------------------  
    start = time.time()   
    os.chdir(WorkingDir)
    DataSheet = pandas.read_excel('DataSheet.xlsx')
    FileNames = list(DataSheet['FILENAME'].values)
    n = len(FileNames)
    Samples = list(DataSheet['SAMPLE NAME'].values)
    ReadsFilename = [FileNames[j] for j in range(n) if Samples[j] == sample][0] 

               
    # ------------------------------------------------
    # Apply cutadapt to trim off the constant region
    # ------------------------------------------------           
    print('Read clipping in progress ...')        
    os.chdir(DataDir) 
    Reads_filename = ReadsFilename[0:-9]
    ReadsCutFilename = Reads_filename + '_cut.fastq.gz'
    ReadsCut_filename = ReadsCutFilename[0:-9]
    CutAdaptCmdLine = CutAdaptDir+'cutadapt -a '+seq_5_end+'...'+seq_3_end  \
                        +' '+ReadsFilename+' -o '+ReadsCutFilename \
                        +' -e '+str(CutErrorTol) \
                        +' --minimum-length '+str(R_min)+' > '+cutadaptLog
    os.system(CutAdaptCmdLine)
    if not os.path.exists(OutputDir):
        os.makedirs(OutputDir)
    mv_cmdline = 'mv '+cutadaptLog+' '+OutputDir
    os.system(mv_cmdline)
    print('Read clipping completed.')
    end = time.time()    
    # Time stamp
    sec_elapsed = end-start
    if sec_elapsed < 60:
        time_elapsed = sec_elapsed
        print('Time elapsed (Read clipping) [secs]: ' + '%.3f' % time_elapsed +'\n')
    elif sec_elapsed < 3600:
        time_elapsed = sec_elapsed/60
        print('Time elapsed (Read clipping) [mins]: ' + '%.3f' % time_elapsed +'\n')
    else:
        time_elapsed = sec_elapsed/3600
        print('Time elapsed (Read clipping) [hours]: ' + '%.3f' % time_elapsed +'\n')    


    # ----------------------------------------------
    # Run alignment
    # ----------------------------------------------                  
    start = time.time()  
    print('Alignment in progress ...')        
    RunBowtie2(ReadsCut_filename,DataDir,AlnDir,bw2Dir,IndexDir,L_bw,N_bw,i_bw)
    print('Alignment completed.')
    end = time.time()
    # Time stamp
    aln_time = end-start
    if aln_time < 60: 
        time_elapsed = aln_time
        print('Time elapsed (Alignment) [secs]: ' + '%.3f' % time_elapsed +'\n')
    elif aln_time < 3600:
        time_elapsed = aln_time/60
        print('Time elapsed (Alignment) [mins]: ' + '%.3f' % time_elapsed +'\n')
    else:
        time_elapsed = aln_time/3600
        print('Time elapsed (Alignment) [hours]: ' + '%.3f' % time_elapsed +'\n')


    # ------------------------------------------
    # Extract and analyze alignments
    # ------------------------------------------ 
    start = time.time()
    print('Analyzing alignment ...') 
    # CLASSIFY ALIGNMENTS 
    os.chdir(AlnDir)
    bw2outputFilename = ReadsCut_filename + '_bw2output.sam'
    bw2sam = pysam.AlignmentFile(bw2outputFilename,'rb')
    fail = pysam.Samfile('failed_alignments.bam','wb',template=bw2sam)
    unique = pysam.Samfile('unique_alignments.bam','wb',template=bw2sam)
    keep = pysam.Samfile('tolerated_alignments.bam','wb',template=bw2sam)
    toss = pysam.Samfile('ambiguous_alignments.bam','wb',template=bw2sam)
    NFail = 0; NUnique = 0; NTol = 0; NAmb = 0
    mapQ = list()
    primScore = list()
    secScore = list()
    AlnStatus = list()
    sgRNA_Hitlist = list()
    for read in bw2sam.fetch():        
        mapQ.append(read.mapping_quality)
        # read with primary and seconday alignment
        if read.has_tag('AS') and read.has_tag('XS'):
            AS = read.get_tag('AS')
            XS = read.get_tag('XS')       
            primScore.append(AS)
            secScore.append(XS)            
            if XS <= AS - Theta:
                keep.write(read); NTol += 1
                AlnStatus.append('Tolerate')
                sgRNA_Hitlist.append(read.reference_name)
            else:
                toss.write(read); NAmb += 1  
                AlnStatus.append('Ambiguous')
        # read with only primary alignment
        elif read.has_tag('AS'):
            AS = read.get_tag('AS')
            primScore.append(AS)
            secScore.append(0)
            unique.write(read); NUnique += 1
            AlnStatus.append('Unique')
            sgRNA_Hitlist.append(read.reference_name)
        # read with failed alignment
        else:
            primScore.append(0)
            secScore.append(0)
            fail.write(read); NFail += 1
            AlnStatus.append('Fail')
    bw2sam.close(); fail.close(); unique.close(); keep.close(); toss.close()   
    os.system('samtools view failed_alignments.bam -o failed_alignments.txt')
    os.system('samtools view unique_alignments.bam -o unique_alignments.txt')
    os.system('samtools view tolerated_alignments.bam -o tolerated_alignments.txt')  
    os.system('samtools view ambiguous_alignments.bam -o ambiguous_alignments.txt')  
    os.system('rm failed_alignments.bam')
    os.system('rm unique_alignments.bam')
    os.system('rm tolerated_alignments.bam')
    os.system('rm ambiguous_alignments.bam')              
    # PRINT ANALYSIS LOGFILE
    print('Writing alignment logfile ...')
    NReads = NTol + NAmb + NUnique + NFail
    FracUnique = round(NUnique/NReads*1000)/10    
    FracTol = round(NTol/NReads*1000)/10
    FracAmb = round(NAmb/NReads*1000)/10
    FracFail = round(NFail/NReads*1000)/10    
    if aln_time < 60:
        time_elapsed = aln_time
        time_unit = ' [secs]'
    elif aln_time < 3600:
        time_elapsed = aln_time/60
        time_unit = ' [mins]'
    else:
        aln_time = aln_time/3600
        time_unit = ' [hours]'   
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
    LogFile.write('---- Alignment completed in %.2f' % time_elapsed + time_unit+' ----')
    LogFile.write('\n\n\n')        
    LogFile.write('Parameter settings\n')
    LogFile.write('Seed Length: '+str(L_bw)+'\n')    
    LogFile.write('Seed Mismatch: '+str(N_bw)+'\n')    
    LogFile.write('Interval Function: '+str(i_bw)+'\n')        
    LogFile.write('Ambiguity Cutoff: '+str(Theta)+'\n')    
    LogFile.close()              
    # DRAW PLOTS
    # MAPPING QUALITY HISTOGRAM
    os.chdir(OutputDir)
    print('Plotting mapping quality ...')
    maxQuality = max(mapQ)
    plt.figure(figsize=(5,4))
    plt.hist(mapQ, bins = range(maxQuality+1), align = 'left')
    plt.title(sample+' Mapping Quality', fontsize=14)
    plt.xlabel('Mapping Quality', fontsize=12)
    plt.ylabel('Number of Reads', fontsize=12)
    plt.tight_layout()
    plt.savefig(sample+'_MappingQuality.png',dpi=res)  
    if svg:
        plt.savefig(sample+'_MappingQuality.svg')  
    # ALIGNMENT SCORE BARPLOT
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
    fig = plt.figure(figsize=(6,5))
    ax = fig.gca(projection='3d')
    ax.bar3d(X,Y,Z_off,dX,dY,Z,color='#00FF00')
    green_proxy = plt.Rectangle((0, 0), 1, 1, fc='#00FF00')
    ax.bar3d(x,y,z_off,dx,dy,z,color='#FF3333')
    red_proxy = plt.Rectangle((0, 0), 1, 1, fc='#FF3333')
    ax.legend([green_proxy,red_proxy],['Matched','Discarded'],loc='upper left',prop={'size':10})
    plt.suptitle(sample+' Alignment Scores', fontsize=14)
    ax.set_xlabel('Prim. Alignment Score', fontsize=12)
    ax.set_ylabel('Sec. Alignment Score', fontsize=12)
    formatter = FuncFormatter(kilos)
    ax.zaxis.set_major_formatter(formatter)    
    ax.set_zlabel('Number of Reads',fontsize=12)
    plt.tight_layout()
    plt.savefig(sample+'_AlignmentScores.png',dpi=res)   
    if svg:
        plt.savefig(sample+'_AlignmentScores.svg',dpi=res)
    print('Alignment analysis completed.')  
    end = time.time()
    # Time stamp
    sec_elapsed = end-start
    if sec_elapsed < 60: 
        time_elapsed = sec_elapsed
        print('Time elapsed (Alignment analysis) [secs]: ' + '%.3f' % time_elapsed +'\n')
    elif sec_elapsed < 3600:
        time_elapsed = sec_elapsed/60
        print('Time elapsed (Alignment analysis) [mins]: ' + '%.3f' % time_elapsed +'\n')
    else:
        time_elapsed = sec_elapsed/3600
        print('Time elapsed (Alignment analysis) [hours]: ' + '%.3f' % time_elapsed +'\n')
        
    
    # --------------------------------------
    # Get read counts
    # --------------------------------------
    start = time.time()    
    print('Read counting in progress ...')
    # Read library
    os.chdir(LibDir)
    LibCols = ['gene','ID','seq']
    LibFile = pandas.read_table(LibFilename, sep = libsep, skiprows = 1, names = LibCols)
    sgIDs = list(LibFile['ID'].values)
    L = len(sgIDs)
    global geneIDs
    geneIDs = list(LibFile['gene'].values)
    G = len(set(geneIDs))
    # Read counts per sgRNA
    print('Counting reads per sgRNA ...')
    ReadCounts = Counter()
    for sgRNA in sgRNA_Hitlist:
        ReadCounts[sgRNA] += 1
    global ReadsPerGuide
    ReadsPerGuide = list()
    for sgRNA in sgIDs:
        ReadsPerGuide.append(ReadCounts[sgRNA])      
    # Apply read count cut-off
    ReadSel = [ReadsPerGuide[k] >= NReads/N0*minN for k in range(L)]
    ReadsPerGuide = [ReadSel[k]*ReadsPerGuide[k] for k in range(L)]
    # Write read counts   
    os.chdir(OutputDir)         
    GuideCountsFilename = sample + GuideCount_Suffix
    GuideCounts = open(GuideCountsFilename,'w')
    for k in range(L):
        GuideCounts.write(str(sgIDs[k]) + '\t'+ str(geneIDs[k]) + '\t' + str(ReadsPerGuide[k]) + '\n')
    GuideCounts.close()
    # No-mapping warning
    if sum(ReadsPerGuide) == 0:
        print('!! ERROR: Zero total read counts! Check library file and index !!')
    # Read counts per gene in library       
    print('Counting reads per gene ...')   
    global GeneList
    GeneList = list(set(geneIDs))   
    G = len(GeneList)  
    num_cores = multiprocessing.cpu_count()
    ReadsPerGene = Parallel(n_jobs=num_cores)(delayed(CountReadsPerGene)(g) for g in range(G))  
    GeneCountsFilename = sample + GeneCount_Suffix
    GeneCounts = open(GeneCountsFilename,'w')
    for g in range(G):
        GeneCounts.write(str(GeneList[g]) + '\t' + str(ReadsPerGene[g]) + '\n')
    GeneCounts.close()        
    end = time.time()
    print('Read counting completed.')
    # Time stamp
    sec_elapsed = end-start
    if sec_elapsed < 60:
        time_elapsed = sec_elapsed
        print('Time elapsed (Read Counting) [secs]: ' + '%.3f' % time_elapsed +'\n')
    elif sec_elapsed < 3600:
        time_elapsed = sec_elapsed/60
        print('Time elapsed (Read Counting) [mins]: ' + '%.3f' % time_elapsed + '\n')
    else:
        time_elapsed = sec_elapsed/3600
        print('Time elapsed (Read Counting) [hours]: ' + '%.3f' % time_elapsed + '\n')


    # --------------------------------------
    # Cleaning up...
    # --------------------------------------
    start = time.time()  
    # clipped reads file    
    if not keepCutReads:
        os.chdir(DataDir)
        os.system('rm '+ReadsCutFilename)    
    # alignment output        
    if AlnOutput == 'Compress':
        print('Compressing raw alignment output...')
        os.chdir(AlnDir)
        # converting SAM to BAM
        SAM_output = glob.glob('*bw2output.sam')[0]
        BAM_output = SAM_output[:-3] + 'bam'
        os.system('samtools view -buSH '+SAM_output+' > '+BAM_output)
        os.system('rm '+SAM_output) 
        # zipping up txt outputs
        os.system('tar -cvf - unique_alignments.txt | gzip -5 - > unique_alignments.tar.gz')
        os.system('tar -cvf - failed_alignments.txt | gzip -5 - > failed_alignments.tar.gz')
        os.system('tar -cvf - tolerated_alignments.txt | gzip -5 - > tolerated_alignments.tar.gz')
        os.system('tar -cvf - ambiguous_alignments.txt | gzip -5 - > ambiguous_alignments.tar.gz')
        os.system('rm *.txt')
    elif AlnOutput == 'Delete':
        print('Removing raw alignment output...')
        os.chdir(AlnDir)
        os.system('rm *')
    elif AlnOutput == 'Keep':
        print('Keeping raw alignment output ...')
    end = time.time()
    # Time stamp
    sec_elapsed = end-start
    if sec_elapsed < 60: 
        time_elapsed = sec_elapsed
        print('Time elapsed (Clean-up) [secs]: ' + '%.3f' % time_elapsed +'\n')
    elif sec_elapsed < 3600:
        time_elapsed = sec_elapsed/60
        print('Time elapsed (Clean-up) [mins]: ' + '%.3f' % time_elapsed +'\n')
    else:
        time_elapsed = sec_elapsed/3600
        print('Time elapsed (Clean-up) [hours]: ' + '%.3f' % time_elapsed +'\n')     

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
    MapAndCount(input1)
