#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 11 12:23:50 2017

@author: philipp
"""

# Number of reads barplot
# =======================================================================
# Imports
import matplotlib.pyplot as plt
import os
import subprocess
import yaml
import sys
import pandas
from matplotlib.ticker import FuncFormatter

def millions(x, pos):
    return '%1.1fM' % (x*1e-6)

def PlotReadDepth():    
    # Get parameters    
    configFile = open('configuration.yaml','r')
    config = yaml.load(configFile)
    configFile.close()
    DataDir = config['DataDir']
    WorkingDir = config['WorkingDir']
    AnalysisDir = config['AnalysisDir']
    DepthDir = config['DepthDir']
    res = config['dpi']
    
    # Get data
    os.chdir(WorkingDir)
    DataSheet = pandas.read_excel('DataSheet.xlsx')
    FileNames = list(DataSheet['FILENAME'].values)
    SampleNames = list(DataSheet['SAMPLE NAME'].values)
    
    # Get reads    
    n = len(SampleNames)
    NumReads = list()
    os.chdir(DataDir)
    print('\nPlotting read depth ... ')
    for k in range(n):
        filename = FileNames[k]
        sample = SampleNames[k]
        NLines = int(subprocess.check_output('gunzip -c '+filename+' | wc -l', shell=True))
        NReads = int(NLines/4)
        print(sample+': '+str(NReads)+' Reads')
        NumReads.append(NReads)
    
    # Plot    
    fig, ax = plt.subplots()
    plt.bar(range(n), NumReads, align='center', color=(0.2588,0.4433,1.0))
    plt.xticks(range(n), SampleNames, rotation='vertical')
    formatter = FuncFormatter(millions)
    ax.yaxis.set_major_formatter(formatter)
    plt.title('Read Depth', fontsize=18, fontweight='bold')
    plt.xlabel('Sample', fontsize=14)  
    plt.ylabel('Number of Reads', fontsize=14)
    plt.tight_layout()
    if not os.path.exists(DepthDir):
        os.makedirs(DepthDir)
    os.chdir(DepthDir)
    plt.savefig('Read_Depth.png',dpi=res)


    
if __name__ == "__main__":
    PlotReadDepth()     