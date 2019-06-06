#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 15:10:27 2016

@author: philipp
"""
# p-value plots
# =======================================================================
# Imports
from __future__ import division # floating point division by default
import numpy
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import sys
import pandas
import os
from matplotlib.ticker import FuncFormatter

def kilos(x, pos):
    return '%1.0fk' % (x*1e-3)

def millions(x, pos):
    return '%1.1fM' % (x*1e-6)

def pvalHist(pval0,pvalDir,sample,res,svg,bar_color,PlotTitle):
    if not os.path.exists(pvalDir):
        os.makedirs(pvalDir)
    os.chdir(pvalDir)   
    if 'N/A' in pval0:       # eliminate N/A
        pval = [pval0[i] for i in range(len(pval0)) if pval0[i]!='N/A']
    else:        
        pval = pval0            
    bin_size = 0.1; min_edge = 0; max_edge = 1
    n = (max_edge-min_edge)/bin_size; Nplus1 = n + 1
    bin_list = numpy.linspace(min_edge, max_edge, Nplus1)    
    fig, ax = plt.subplots(figsize=(3.5,2.9))
    plt.hist(pval,bin_list,color=bar_color,edgecolor='#4e4f51')
    plt.xticks([0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1])   
    for label in ax.xaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    plt.xlabel('p-value', fontsize=11)    
    plt.ylabel('Frequency', fontsize=11) 
    if len(pval)>10000:
        formatter = FuncFormatter(kilos)
        ax.yaxis.set_major_formatter(formatter)
    plt.tick_params(labelsize=11)    
    plt.title(PlotTitle, fontsize=11) 
    plt.tight_layout()
    plt.savefig(sample+' '+PlotTitle+'.png', dpi=res)
    if svg:
        plt.savefig(sample+' '+PlotTitle+'.svg')            
