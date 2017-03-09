#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 18:34:42 2017

@author: philipp
"""

import sys
import os
import yaml


def RunSeqQC():
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
        os.system('fastqc -o '+SeqQCDir+' '+filename)
    os.chdir(SeqQCDir)    
    os.system('rm *.zip')
    os.chdir(ScriptsDir)

if __name__ == "__main__":
    RunSeqQC() 