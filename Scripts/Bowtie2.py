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
import multiprocessing

def BuildIndex(LibFastA,IndexDir,bw2Dir):
    os.chdir(IndexDir)
    bw2_cmdline = bw2Dir+'bowtie2-build -f library.fasta Library'
    os.system(bw2_cmdline)

def RunBowtie2(ReadsFilename,DataDir,AlnDir,bw2Dir,IndexDir,L_bw,N_bw,i_bw):   
    bw2output = ReadsFilename + '_bw2output.sam'
    os.chdir(DataDir) 
    cp_cmdline = 'cp ' + ReadsFilename + ' ' + IndexDir   
    os.system(cp_cmdline)
    os.chdir(IndexDir)
    num_cores = multiprocessing.cpu_count()
    bw2_cmdline = bw2Dir+'bowtie2 --local -L '+str(L_bw)+' -N '+str(N_bw)+' -i '+str(i_bw) +\
        ' -q -x Library -U ' + ReadsFilename + ' -S ' + bw2output + ' -p '+str(num_cores)
    os.system(bw2_cmdline)
    rm_cmdline = 'rm ' + ReadsFilename
    os.system(rm_cmdline)
    if not os.path.exists(AlnDir):
        os.makedirs(AlnDir)
    os.chdir(IndexDir)
    mv_cmdline = 'mv ' + bw2output + ' ' + AlnDir
    os.system(mv_cmdline)