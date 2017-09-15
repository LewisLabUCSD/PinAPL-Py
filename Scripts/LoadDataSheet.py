#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 11 14:57:08 2017

@author: philipp
"""
# Load excel sheet with read filenames
# =======================================================================
import sys
import os
import yaml
import pandas

def LoadExcelDataSheet():
    # Get parameters
    configFile = open('configuration.yaml','r')
    config = yaml.load(configFile)
    configFile.close()
    WorkingDir = config['WorkingDir']    
    ScriptsDir = config['ScriptsDir']

    # Get read files
    os.chdir(WorkingDir)
    DataSheet = pandas.read_excel('DataSheet.xlsx')
    FileNames = list(DataSheet['FILENAME'].values)
    TreatmentList = list(DataSheet['TREATMENT'])
    Treatments = list(set(TreatmentList))
    if 'Control' in Treatments:
        N = len(FileNames)
        SampleNames = ['' for n in range(N)]
        for treatment in Treatments:    
            n = 0
            for j in range(N):    
               if TreatmentList[j] == treatment:
                   n += 1
                   SampleNames[j] = treatment+'_'+str(n)
        for j in range(N):
            print('Found sample '+SampleNames[j])
        DataSheet['SAMPLE NAME'] = SampleNames           
        DataSheet.to_excel('DataSheet.xlsx',columns=['FILENAME','TREATMENT','SAMPLE NAME'])
        os.chdir(ScriptsDir)
    else:
        print('### ERROR: No control treatment defined! ###')
    
    
if __name__ == "__main__":
    LoadExcelDataSheet()
