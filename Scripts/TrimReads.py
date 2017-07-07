#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 14:43:38 2017

@author: philipp
"""

import sys
import os
import yaml
import time

def RunCutadapt():
    start_total = time.time()  
    # Get parameters    
    configFile = open('configuration.yaml','r')
    config = yaml.load(configFile)
    configFile.close()    
    DataDir = config['DataDir']
    TempDataDir = config['TempDataDir']
    CutAdaptDir = config['CutAdaptDir']
    TrimLogDir = config['TrimLogDir']
    seq_5_end = config['seq_5_end']
    CutErrorTol = config['CutErrorTol']
    R_min = config['R_min']
    
    # Trim 5' adapters
    print('Trimming 5\' adapters ...') 
    if not os.path.exists(TempDataDir):
        os.makedirs(TempDataDir)
    if not os.path.exists(TrimLogDir):
        os.makedirs(TrimLogDir)        
    os.chdir(DataDir)
    FileNames = [d for d in os.listdir(DataDir)]
    command5 = str()    
    for ReadsFilename in FileNames:
        ReadsFilename5 = 'Trim5_'+ReadsFilename
        LogFilename = 'cutadapt_'+ReadsFilename+'.txt'
        command5 = command5 + CutAdaptDir+'cutadapt -g '+seq_5_end\
                            +' '+ReadsFilename+' -o '+TempDataDir+ReadsFilename5\
                            +' -e '+str(CutErrorTol)+' -m '+str(R_min)\
                            +' > '+LogFilename+' & '
    command5 = command5[:-3]
    os.system(command5)
    os.system('mv *.txt '+TrimLogDir)
    
    # Trim 3' adapters
    print('Trimming 3\' adapters ...') 
    os.chdir(TempDataDir)
    command3 = str()
    for ReadsFilename in FileNames:
        ReadsFilename5 = 'Trim5_'+ReadsFilename
        ReadsFilename53 = 'Trim53_'+ReadsFilename
        DumpFilename = ReadsFilename+'_dump.log'
        command3 = command3 + CutAdaptDir+'cutadapt -l 20 '\
                            +ReadsFilename5+' -o '+ReadsFilename53+' >'+DumpFilename+' & '
    command3 = command3[:-3]
    os.system(command3)    
    os.system('rm *.log')
    end = time.time()
    
    # Time stamp
    end_total = time.time()
    # Final time stamp
    print('------------------------------------------------')
    print('Script completed.')    
    sec_elapsed = end_total - start_total
    if sec_elapsed < 60:
        time_elapsed = sec_elapsed
        print('Time elapsed (Read trimming) [secs]: ' + '%.3f' % time_elapsed +'\n')
    elif sec_elapsed < 3600:
        time_elapsed = sec_elapsed/60
        print('Time elapsed (Read trimming) [mins]: ' + '%.3f' % time_elapsed +'\n')
    else:
        time_elapsed = sec_elapsed/3600
        print('Time elapsed (Read trimming) [hours]: ' + '%.3f' % time_elapsed +'\n')





if __name__ == "__main__":
    RunCutadapt() 