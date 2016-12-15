#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 17:53:39 2016

@author: philipp
"""

# Bowtie2 alignment
# =======================================================================
# Imports
from __future__ import division # floating point division by default
import os
import pandas as pd
import time
import matplotlib
matplotlib.use('Agg') 
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import matplotlib.pyplot as plt
import subprocess


def BuildIndex(LibFastA,IndexDir,bw2Dir):
    os.chdir(IndexDir)
    bw2_cmdline = bw2Dir+'bowtie2-build -f library.fasta Library'
    os.system(bw2_cmdline)

def RunBowtie2(ReadsCut_filename,ReadsDir,AlnDir,LogDir,bw2Dir,IndexDir,Theta,L_bw,N_bw,i_bw,res):   
    # ------------------------------------------    
    # Set output parameters
    # ------------------------------------------    
    ReadsCutFilename = ReadsCut_filename + '.fastq.gz'    
    SAMoutput = ReadsCut_filename + '_bw2output.sam'
    SAM_AlnFilename = ReadsCut_filename + '_SAM_alignments.txt'
    AlnOutputFile = ReadsCut_filename + '_bw2Aln.tsv'
    logfilename = 'Bowtie2_Results.txt'    
       
    # ------------------------------------------
    # Run bowtie alignment
    # ------------------------------------------ 
    os.chdir(ReadsDir)
    cp_cmdline = 'cp ' + ReadsCutFilename + ' ' + IndexDir   
    os.system(cp_cmdline)
    os.chdir(IndexDir)
    bw2_cmdline1 = bw2Dir+'bowtie2 --local -L '+str(L_bw)+' -N '+str(N_bw)+' -i '+str(i_bw) 
    bw2_cmdline2 = ' -q -x Library -U ' + ReadsCutFilename + ' -S ' + SAMoutput
    bw2_cmdline = bw2_cmdline1 + bw2_cmdline2
    start = time.time()
    os.system(bw2_cmdline)
    end = time.time()
    sec_elapsed = end - start
    rm_cmdline = 'rm ' + ReadsCutFilename
    os.system(rm_cmdline)
    print('Bowtie2 Alignment completed.')
    
    # ------------------------------------------
    # Format output
    # ------------------------------------------ 
    if not os.path.exists(AlnDir):
        os.makedirs(AlnDir)
    os.chdir(IndexDir)
    mv_cmdline = 'mv ' + SAMoutput + ' ' + AlnDir
    os.system(mv_cmdline)
    # Read bowtie2 output
    os.chdir(AlnDir)
    sam_file = open(SAMoutput)
    all_text = sam_file.read()
    sam_file.close()
    alltext = all_text.splitlines()
    # Extract alignments
    AlnStartFound = False
    k = -1
    while AlnStartFound == False:
        k += 1
        currline = alltext[k]
        if currline.find('ID:bowtie2') != -1:
            AlnStartFound = True
    m1 = k+1                     # Alignment start in SAM file
    m2 = len(alltext)            # Alignment end in SAM file    
    # Write alignments-only file
    with open(SAM_AlnFilename,'w') as outfile:
        for k in range(m1,m2):
            outfile.write(alltext[k]+'\n')
    outfile.close()


    # ------------------------------------------
    # Analyze alignment results
    # ------------------------------------------   
    # Determine Number of Reads
    os.chdir(ReadsDir)
    NLines = int(subprocess.check_output('gunzip -c '+ReadsCutFilename+' | wc -l', shell=True))
    NReads = int(NLines/4)    
    # Extract alignment results
    os.chdir(AlnDir)
    colnames = ['ReadIndex','flag','libID','offset','MappingQuality','cigar','MateRef', \
                'MateOffset','MateFrag','ReadSeq','Quality','Opt01','Opt02','Opt03', \
                'Opt04','Opt05','Opt06','Opt07','Opt08','Opt09','Opt10','Opt11','Opt12']
    sam_results = pd.read_table(SAM_AlnFilename, sep = '\t', names = colnames)
    bw_libIDs = sam_results.libID.tolist()
    bw_Opt01 = sam_results.Opt01.tolist()
    bw_Opt02 = sam_results.Opt02.tolist()
    bw_quality = sam_results.MappingQuality.tolist()
    # Extract primary (AS) and secondary (XS) alignment scores    
    AlnPrimScore = list()
    AlnSecScore = list()
    AlnStatus = list()
    for k in range(NReads):
        if bw_libIDs[k] == '*':
            AlnPrimScore.append(0)
            AlnSecScore.append(0)
            AlnStatus.append('Fail')
        else:
            if 'XS' not in bw_Opt02[k]:
                AlnPrimScore.append(int(bw_Opt01[k][5:]))
                AlnSecScore.append(0)
                AlnStatus.append('Unique')
            else:
                AlnPrimScore.append(int(bw_Opt01[k][5:]))
                AlnSecScore.append(int(bw_Opt02[k][5:]))
                if AlnSecScore[k] <= AlnPrimScore[k] - Theta:
                    AlnStatus.append('Tolerate')
                else:
                    AlnStatus.append('Ambiguous')     
               
    # ------------------------------------------
    # Print output files
    # ------------------------------------------  
    NFail = int(AlnStatus.count('Fail'))
    NUnique = int(AlnStatus.count('Unique'))
    NTol = int(AlnStatus.count('Tolerate'))     
    NAmb = int(AlnStatus.count('Ambiguous'))
    FracUnique = round(NUnique/NReads*1000)/10    
    FracTol = round(NTol/NReads*1000)/10
    FracAmb = round(NAmb/NReads*1000)/10
    FracFail = round(NFail/NReads*1000)/10    
    if sec_elapsed < 60:
        time_elapsed = sec_elapsed
        time_unit = ' [secs]'
    elif sec_elapsed < 3600:
        time_elapsed = sec_elapsed/60
        time_unit = ' [mins]'
    else:
        time_elapsed = sec_elapsed/3600
        time_unit = ' [hours]'
    # Write logfile
    if not os.path.exists(LogDir):
        os.makedirs(LogDir)      
    os.chdir(LogDir)    
    LogFile = open(logfilename,'w')         
    LogFile.write('Bowtie2 Alignment Results\n')
    LogFile.write('**************************************\n')
    LogFile.write('Call: '+bw2_cmdline+'\n')
    LogFile.write('Ambiguity Cutoff: '+str(Theta)+'\n')
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
    LogFile.close()              
    # Write alignment results file
    os.chdir(AlnDir)
    print('Generating alignment file ...')
    AlnFile = open(AlnOutputFile,'w')  
    for k in range(NReads):
        if AlnStatus[k] not in ['Unique','Tolerate']:
            bw_libIDs[k] = '*'
        AlnFile.write(bw_libIDs[k] + '\t' + str(bw_quality[k]) + '\n')
    AlnFile.close()
        
    # ------------------------------------------
    # Plots
    # ------------------------------------------      
    # Plot mapping quality histogram
    os.chdir(LogDir)
    print('Generating mapping quality histogram...')
    maxQuality = max(bw_quality)
    plt.hist(bw_quality, bins = range(maxQuality+1), align = 'left')
    plt.suptitle('Bowtie2 Mapping Quality', fontsize=16, fontweight='bold')
    plt.xlabel('Mapping Quality', fontsize=14)
    plt.ylabel('Number of Reads', fontsize=14)
    plt.savefig('Bowtie2_MappingQuality_Hist.png',dpi=res)  
    # Alignment score barplot
    print('Generating alignment score barplot...')
    AlnPrimAcc = [AlnPrimScore[k] for k in range(NReads) if AlnStatus[k] in ['Unique','Tolerate']]
    AlnSecAcc = [AlnSecScore[k] for k in range(NReads) if AlnStatus[k] in ['Unique','Tolerate']]    
    AlnPrimDisc = [AlnPrimScore[k] for k in range(NReads) if AlnStatus[k] in ['Fail','Ambiguous']]
    AlnSecDisc = [AlnSecScore[k] for k in range(NReads) if AlnStatus[k] in ['Fail','Ambiguous']]
    # Plot bars (matched reads)
    xedges = range(42)
    yedges = range(42)
    H, xedges, yedges = np.histogram2d(AlnSecAcc, AlnPrimAcc, bins=(xedges, yedges))
    Y, X = np.nonzero(H)
    Z = H[Y,X]
    Z_off = np.zeros(len(Y))
    dX = np.ones(len(X))
    dY = np.ones(len(Y))
    # Plot bars (non-matched reads)
    h, xedges, yedges = np.histogram2d(AlnSecDisc, AlnPrimDisc, bins=(xedges, yedges))
    y, x = np.nonzero(h)
    z = h[y,x]
    z_off = np.zeros(len(y))
    dx = np.ones(len(x))
    dy = np.ones(len(y))
    # Show plot    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.bar3d(X,Y,Z_off,dX,dY,Z,color='#00FF00')
    green_proxy = plt.Rectangle((0, 0), 1, 1, fc='#00FF00')
    ax.bar3d(x,y,z_off,dx,dy,z,color='#FF3333')
    red_proxy = plt.Rectangle((0, 0), 1, 1, fc='#FF3333')
    ax.legend([green_proxy,red_proxy],['Matched','Discarded'],loc='upper left')
    plt.suptitle('Bowtie2 Alignment Scores', fontsize=16, fontweight='bold')
    ax.set_xlabel('Primary Alignment Score')
    ax.set_ylabel('Secondary Alignment Score')
    ax.set_zlabel('Number of Reads')
    plt.savefig('Bowtie2_AlignmentScores.png',dpi=res)
