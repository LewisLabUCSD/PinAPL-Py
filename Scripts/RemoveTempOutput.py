#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 11:24:09 2018

@author: philipp
"""
# Remove temporary trimming and alignment output
# =======================================================================
# Imports
import sys
import yaml
import os
import glob
import time

def CleanUp():  
    # ------------------------------------------------
    # Get parameters
    # ------------------------------------------------
    configFile = open('configuration.yaml','r')
    config = yaml.load(configFile)
    configFile.close()
    AlnOutput = config['AlnOutput']
    keepCutReads = config['keepCutReads']    
    TempDataDir = config['TempDataDir']
    AlnStemDir = config['AlignDir']
    ScriptsDir = config['ScriptsDir']
    
    # --------------------------------------
    # Cleaning up...
    # --------------------------------------
    start = time.time()  
    # clipped reads file    
    if not keepCutReads:        
        os.chdir(TempDataDir)
        os.system('rm *.*')
    # alignment output        
    if AlnOutput == 'Compress':
        print('Compressing raw alignment output...')
        os.chdir(AlnStemDir)
        SubDirs = [d for d in os.listdir(AlnStemDir)]
        for Dir in SubDirs:
            os.chdir(Dir)
            # converting SAM to BAM
            SAM_output = glob.glob('*bw2output.sam')[0]
            BAM_output = SAM_output[:-3] + 'bam'
            os.system('samtools view -bS '+SAM_output+' > '+BAM_output)
            os.system('rm '+SAM_output) 
            os.chdir(AlnStemDir)
    elif AlnOutput == 'Delete':
        print('Removing raw alignment output...')
        os.chdir(AlnStemDir)
        SubDirs = [d for d in os.listdir(AlnStemDir)]
        for Dir in SubDirs:
            os.chdir(Dir)
            os.system('rm *')
            os.chdir(AlnStemDir)
    elif AlnOutput == 'Keep':
        print('Keeping raw alignment output ...')
    end = time.time()
    os.chdir(ScriptsDir)
    # Time stamp    
    sec_elapsed = end-start
    if sec_elapsed < 60: 
        time_elapsed = sec_elapsed
        print('Time elapsed [secs]: ' + '%.3f' % time_elapsed)
    elif sec_elapsed < 3600:
        time_elapsed = sec_elapsed/60
        print('Time elapsed [mins]: ' + '%.3f' % time_elapsed)
    else:
        time_elapsed = sec_elapsed/3600
        print('Time elapsed [hours]: ' + '%.3f' % time_elapsed)    
        

if __name__ == "__main__":
    CleanUp()