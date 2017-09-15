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
import subprocess
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
    sgLength = config['sgRNALength']
    R_min = sgLength
    ScriptsDir = config['ScriptsDir']
    RunInBack = 'RunInBackground.sh'
    
    # Trim 5' adapters
    if not os.path.exists(TempDataDir):
        os.makedirs(TempDataDir)
    if not os.path.exists(TrimLogDir):
        os.makedirs(TrimLogDir)        
    os.chdir(DataDir)
    FileNames = [d for d in os.listdir(DataDir)]
    os.chdir(TempDataDir)
    filesNotDone = []
    # Load processes into the background
    for ReadsFilename in FileNames:
        ReadsFilename0 = 'Trim_'+ReadsFilename
        LogFilename = 'cutadapt_'+ReadsFilename+'.txt'
        command = CutAdaptDir+'cutadapt -g '+seq_5_end\
                                +' '+DataDir+ReadsFilename+' -o '+ReadsFilename0\
                                +' -e '+str(CutErrorTol)+' -m '+str(R_min)+' -l '+str(sgLength)\
                                +' 2>&1 > '+LogFilename
        subprocess.call(ScriptsDir+RunInBack+' "'+command+'" '+ReadsFilename+' cutadapt_status.log &',shell=True,\
            stdin=None, stdout=None, stderr=None, close_fds=True)
        filesNotDone.append(ReadsFilename+'_cutadapt_status.log')
        print('Loading '+ReadsFilename)
    print('Removing adapters ...')
    print('Extracting '+str(sgLength)+' bp sgRNA sequences ...')
    # Check for completion
    while len(filesNotDone) > 0:
        time.sleep(.1)
        statusfile = filesNotDone[0]
        status = open(statusfile,'r').readline().rstrip('\n')
        if status == 'done':
            filesNotDone.pop(0)
            os.system('rm '+statusfile)        
    print('Adapter removal completed.')
    print('Writing logfiles...')    
    os.system('mv *.txt '+TrimLogDir)    
       
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