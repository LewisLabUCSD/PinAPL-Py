#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 11 12:23:50 2017

@author: philipp
"""

# Number of reads barplot
# =======================================================================
# Imports
from __future__ import division # floating point division by default
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import os
import yaml
import sys
import re
import pandas
import numpy
from matplotlib.ticker import FuncFormatter
from matplotlib.font_manager import FontProperties

def millions(x, pos):
    return '%1.1fM' % (x*1e-6)

def PlotReadDepth():    
    # Get parameters    
    configFile = open('configuration.yaml','r')
    config = yaml.load(configFile)
    configFile.close()
    AlnQCDir = config['AlnQCDir']
    WorkingDir = config['WorkingDir']
    DepthDir = config['DepthDir']
    ScriptsDir = config['ScriptsDir']
    res = config['dpi']
    svg = config['svg']
    
    # Get data
    os.chdir(WorkingDir)
    DataSheet = pandas.read_excel('DataSheet.xlsx')
    FileNames = list(DataSheet['FILENAME'].values)
    SampleNames = list(DataSheet['SAMPLE NAME'].values)
    
    # Get number of reads per category        
    os.chdir(AlnQCDir)    
    S = len(SampleNames)
    N1 = list(); N2 = list()
    n1 = list(); n2 = list()    
    for sample in SampleNames:
        os.chdir(sample)
        filename = sample+'_AlignmentResults.txt'
        AlnStats = open(filename,'r')
        Lines = AlnStats.read()
        p = re.search('(?<=Number of Reads with unique Alignments: \t)\d+',Lines)
        N1.append(int(p.group(0)))
        p = re.search('(?<=Number of Reads above Ambiguity Tolerance: \t)\d+',Lines)
        N2.append(int(p.group(0)))
        p = re.search('(?<=Number of Reads below Ambiguity Tolerance: \t)\d+',Lines)
        n1.append(int(p.group(0)))
        p = re.search('(?<=Number of Reads with failed Alignment: \t\t)\d+',Lines)
        n2.append(int(p.group(0)))
        os.chdir(AlnQCDir)
    N2N1 = numpy.add(N2,N1)
    n1N2N1 = numpy.add(N2N1,n1)
    
    # Plot    
    S = len(SampleNames)
    fig, ax = plt.subplots(figsize=(4,3.5))
    plt.bar(range(S), N1, align='center', color=(66/255,128/255,244/255),\
        label='Alignment Unique')      
    plt.bar(range(S), N2, bottom = N1, align='center', color=(66/255,244/255,217/255),\
        label='Alignment Tolerated')          
    plt.bar(range(S), n1, bottom = N2N1, align='center', color=(244/255,188/255,66/255),\
        label='Alignment Ambiguous')
    plt.bar(range(S), n2, bottom = n1N2N1, align='center', color=(244/255,66/255,66/255),\
        label='Alignment Failed')
    plt.xticks(range(S), SampleNames, rotation='vertical', fontsize=11)
    formatter = FuncFormatter(millions)
    ax.yaxis.set_major_formatter(formatter)
    plt.title('Read Depth', fontsize=12)
    plt.ylabel('Number of Reads', fontsize=11)
    plt.xlim([-1,S])    
    plt.tick_params(labelsize=11)
    fontP = FontProperties()
    fontP.set_size('x-small')
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], prop=fontP,loc='lower center',bbox_to_anchor=(1.3, 0.5))
    if not os.path.exists(DepthDir):
        os.makedirs(DepthDir)
    os.chdir(DepthDir)
    plt.savefig('Read_Depth.png',dpi=res,bbox_inches="tight")
    if svg:
        plt.savefig('Read_Depth.svg',bbox_inches="tight")
    os.chdir(ScriptsDir)


    
if __name__ == "__main__":
    PlotReadDepth()     