# -*- coding: utf-8 -*-
"""
Created on Sat Feb 11 14:57:08 2017

@author: philipp
"""
# Load sample names, treatments and replicates
# =======================================================================
import sys
import os
import yaml
import pandas

def GetSamples():
    configFile = open('configuration.yaml','r')
    config = yaml.load(configFile)
    configFile.close()
    WorkingDir = config['WorkingDir']     
    ScriptsDir = config['ScriptsDir']
    os.chdir(WorkingDir)
    DataSheet = pandas.read_excel('DataSheet.xlsx')
    FileNames = list(DataSheet['FILENAME'].values)
    TreatmentList = list(DataSheet['TREATMENT'].values)
    SampleNames = list(DataSheet['SAMPLE NAME'].values)
    TreatmentSamples = list()
    for sample in SampleNames:
        if 'Control' not in sample:
            TreatmentSamples.append(sample)
    Treatments = list(set(TreatmentList))
    Replicates = dict()
    n = len(SampleNames)
    for treatment in Treatments:
        replist = [SampleNames[j] for j in range(n) if TreatmentList[j] == treatment]
        Replicates[treatment] = replist
    os.chdir(ScriptsDir)
    return SampleNames, Treatments, TreatmentSamples, Replicates
