#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 18:34:42 2017

@author: philipp
"""

import sys
import os
import time
import yaml


def RunSeqQC():
    start_total = time.time()
    # Get parameters    
    configFile = open('configuration.yaml','r')
    config = yaml.load(configFile)
    configFile.close()
    WorkingDir = config['WorkingDir']
    DataDir = config['DataDir']
    SeqQCDir = config['SeqQCDir']
    ScriptsDir = config['ScriptsDir']

    # Run fastqc
    if not os.path.exists(SeqQCDir):
        os.makedirs(SeqQCDir)
    os.chdir(DataDir)
    FileNames = [d for d in os.listdir(DataDir)]
    for filename in FileNames:
        if filename[-8:] == 'fastq.gz':
            os.system('fastqc -o '+SeqQCDir+' --extract '+filename)
        elif filename[-5:] == 'fastq':
            os.system('fastqc -o '+SeqQCDir+' '+filename)
    os.chdir(SeqQCDir)    
    os.system('rm *.zip')
    os.chdir(ScriptsDir)
    end = time.time()    
    
    end_total = time.time()
    # Final time stamp
    print('------------------------------------------------')
    print('Script completed.')    
    sec_elapsed = end_total - start_total
    if sec_elapsed < 60:
        time_elapsed = sec_elapsed
        print('Time elapsed (Quality control) [secs]: ' + '%.3f' % time_elapsed +'\n')
    elif sec_elapsed < 3600:
        time_elapsed = sec_elapsed/60
        print('Time elapsed (Quality control) [mins]: ' + '%.3f' % time_elapsed +'\n')
    else:
        time_elapsed = sec_elapsed/3600
        print('Time elapsed (Quality control) [hours]: ' + '%.3f' % time_elapsed +'\n')


if __name__ == "__main__":
    RunSeqQC() 