#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 17:53:39 2016

@author: philipp
"""

# Bowtie2 alignment
# =======================================================================
# Imports
import os


def BuildIndex(LibFastA,IndexDir,bw2Dir):
    os.chdir(IndexDir)
    bw2_cmdline = bw2Dir+'bowtie2-build -f library.fasta Library'
    os.system(bw2_cmdline)

def RunBowtie2(ReadsCut_filename,DataDir,AlnDir,bw2Dir,IndexDir,L_bw,N_bw,i_bw):   
    # ------------------------------------------    
    # Set output parameters
    # ------------------------------------------    
    ReadsCutFilename = ReadsCut_filename + '.fastq.gz'    
    bw2output = ReadsCut_filename + '_bw2output.sam'
               
    # ------------------------------------------
    # Run bowtie alignment
    # ------------------------------------------ 
    os.chdir(DataDir) 
    cp_cmdline = 'cp ' + ReadsCutFilename + ' ' + IndexDir   
    os.system(cp_cmdline)
    os.chdir(IndexDir)
    bw2_cmdline1 = bw2Dir+'bowtie2 --local -L '+str(L_bw)+' -N '+str(N_bw)+' -i '+str(i_bw) 
    bw2_cmdline2 = ' -q -x Library -U ' + ReadsCutFilename + ' -S ' + bw2output
    bw2_cmdline = bw2_cmdline1 + bw2_cmdline2
    os.system(bw2_cmdline)
    rm_cmdline = 'rm ' + ReadsCutFilename
    os.system(rm_cmdline)
    if not os.path.exists(AlnDir):
        os.makedirs(AlnDir)
    os.chdir(IndexDir)
    mv_cmdline = 'mv ' + bw2output + ' ' + AlnDir
    os.system(mv_cmdline)