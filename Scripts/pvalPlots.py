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

def pvalHist(NBpval,NBpval_0,pvalDir,sample,res,svg):
    sampleDir = pvalDir+sample    
    if not os.path.exists(sampleDir):
        os.makedirs(sampleDir)
    os.chdir(sampleDir)      
    bin_size = 0.1; min_edge = 0; max_edge = 1
    n = (max_edge-min_edge)/bin_size; Nplus1 = n + 1
    bin_list = numpy.linspace(min_edge, max_edge, Nplus1)    
    plt.figure(figsize=(5,4))
    plt.hist(NBpval,bin_list,color=(106/255,240/255,247/255),label='Unadjusted')
    plt.hist(NBpval_0,bin_list,color=(232/255,27/255,16/255),rwidth=0.8,alpha=0.75,label='Adjusted')
    plt.xticks([0,.1,.2,.3,.4,.5,0.6,.7,.8,.9,1])
    plt.xlabel('p-value', fontsize=12)    
    plt.ylabel('Frequency', fontsize=12) 
    plt.title(sample+' sgRNA Significance', fontsize=14) 
    plt.legend(loc='upper center', prop={'size':10})    
    plt.tight_layout()
    plt.savefig(sample+'_sgRNA_Significance.png', dpi=res)
    if svg:
        plt.savefig(sample+'_sgRNA_Significance.svg')        


def pvalHist_metric(pval_list,pval_list0,GeneMetric,pvalDir,sample,res,svg):
    sampleDir = pvalDir+sample    
    if not os.path.exists(sampleDir):
        os.makedirs(sampleDir)
    os.chdir(sampleDir)      
    bin_size = 0.1; min_edge = 0; max_edge = 1
    n = (max_edge-min_edge)/bin_size; Nplus1 = n + 1
    bin_list = numpy.linspace(min_edge, max_edge, Nplus1)    
    plt.figure(figsize=(5,4))
    plt.hist(pval_list,bin_list,color=(150/255,255/255,94/255),label=GeneMetric+' Unadjusted')
    plt.hist(pval_list0,bin_list,color=(222/255,32/255,247/255),rwidth=0.8,alpha=0.75,\
        label=GeneMetric+' Adjusted')
    plt.xticks([0,.1,.2,.3,.4,.5,0.6,.7,.8,.9,1])
    plt.xlabel('p-value', fontsize=12)    
    plt.ylabel('Frequency', fontsize=12) 
    plt.title(sample+' Gene Significance', fontsize=14) 
    plt.legend(loc='upper center', prop={'size':8})    
    plt.tight_layout()
    plt.savefig(sample+'_Gene_Significance.png', dpi=res)
    if svg:    
        plt.savefig(sample+'_Gene_Significance.svg')


def HalfVolcanoPlot(metric,pval_list,metric_sig,GeneMetric,pvalDir,ScreenType,sample,res,svg): 
    sampleDir = pvalDir+sample    
    if not os.path.exists(sampleDir):
        os.makedirs(sampleDir)
    os.chdir(sampleDir) 
    eps = 1e-17
    L = len(metric)
    neglogm = [-numpy.log10(metric[k]+eps) for k in range(L) if metric_sig[k]==False]
    neglogp = [-numpy.log10(pval_list[k]+eps) for k in range(L) if metric_sig[k]==False]
    neglogm_sig = [-numpy.log10(metric[k]+eps) for k in range(L) if metric_sig[k]==True]
    neglogp_sig = [-numpy.log10(pval_list[k]+eps) for k in range(L) if metric_sig[k]==True]      
    plt.figure(figsize=(5,4))
    plt.scatter(neglogm,neglogp,s=3,facecolor='black',lw=0,alpha=0.35)  
    plt.scatter(neglogm_sig,neglogp_sig,s=3,facecolor='green',lw=0,alpha=0.35,label='Significant')    
    xmax = max(max(neglogm_sig),max(neglogm))
    ymax = max(max(neglogp_sig),max(neglogp))
    plt.xlim([0,1.05*xmax]); plt.ylim([0,1.05*ymax])
    plt.xlabel('-log10 '+GeneMetric, fontsize=12)    
    plt.ylabel('-log10 p-value', fontsize=12) 
    plt.title(sample+' '+GeneMetric+' Metric', fontsize=14) 
    plt.legend(loc='lower right', prop={'size':10})
    plt.tight_layout()
    plt.savefig(sample+'_'+GeneMetric+'_Metric.png', dpi=res)
    if svg:
        plt.savefig(sample+'_'+GeneMetric+'_Metric.svg')        
    
    
def VolcanoPlot(fc,NBpval2,significant,pvalDir,ScreenType,sample,res,svg): 
    sampleDir = pvalDir+sample    
    if not os.path.exists(sampleDir):
        os.makedirs(sampleDir)
    os.chdir(sampleDir)  
    L = len(fc)
    logfc = [numpy.log10(fc[k]) for k in range(L) if significant[k]==False]
    neglogp2 = [-numpy.log10(NBpval2[k]) for k in range(L) if significant[k]==False]
    logfc_sig = [numpy.log10(fc[k]) for k in range(L) if significant[k]==True]
    neglogp2_sig = [-numpy.log10(NBpval2[k]) for k in range(L) if significant[k]==True]    
    plt.figure(figsize=(5,4))
    plt.scatter(logfc,neglogp2,s=3,facecolor='grey',lw=0,alpha=0.35)
    plt.scatter(logfc_sig,neglogp2_sig,s=3,facecolor='green',lw=0,alpha=0.35,label='Significant')    
    xmin = min(min(logfc_sig),min(logfc))        
    xmax = max(max(logfc_sig),max(logfc))
    ymax = max(max(neglogp2_sig),max(neglogp2))
    plt.xlim([0.95*xmin,1.05*xmax]); plt.ylim([0,1.05*ymax])
    plt.xlabel('log10 fold change', fontsize=12)    
    plt.ylabel('-log10 p-value (two-sided)', fontsize=12) 
    plt.title(sample+' sgRNA '+ScreenType, fontsize=14) 
    plt.legend(loc='lower right', prop={'size':10})
    plt.tight_layout()
    plt.savefig(sample+'_'+'sgRNA_volcano.png', dpi=res)
    if svg:
        plt.savefig(sample+'_'+'sgRNA_volcano.svg')


def QQPlot(NBpval,significant,pvalDir,sample,res,svg): 
    sampleDir = pvalDir+sample    
    if not os.path.exists(sampleDir):
        os.makedirs(sampleDir)
    os.chdir(sampleDir)  
    L = len(NBpval)
    neglogp = [-numpy.log10(NBpval[k]) for k in range(L)]
    neglogp.sort()    
    pExp = numpy.linspace(1,1/L,L)
    neglogpExp = [-numpy.log10(pExp[k]) for k in range(L)]
    neglogpExp.sort() 
    sig = [significant[k] for k in range(L)]
    sig.sort()
    S = list(sig).index(True)
    plt.figure(figsize=(5,4))
    plt.scatter(neglogpExp[0:S],neglogp[0:S],s=8,facecolor='grey',lw=0,alpha=0.35)
    plt.scatter(neglogpExp[S:],neglogp[S:],s=8,facecolor='green',lw=0,alpha=0.35,label='Significant')
    xmax = 1.05*max(neglogpExp)
    ymax = 1.05*max(neglogp)
    plt.xlim([0,xmax]); plt.ylim([0,ymax])
    plt.plot((0,xmax), (0,xmax), ls="--", color=(51/255,153/255,1))
    plt.xlabel('-log10 expected p-value', fontsize=12)    
    plt.ylabel('-log10 observed p-value', fontsize=12) 
    plt.title(sample+' sgRNA QQ', fontsize=14) 
    plt.legend(loc='upper left', prop={'size':10})
    plt.tight_layout()
    plt.savefig(sample+'_'+'sgRNA_QQ.png', dpi=res)    
    if svg:
        plt.savefig(sample+'_'+'sgRNA_QQ.svg') 

def zScorePlot(fc,significant,pvalDir,ScreenType,sample,res,svg): 
    sampleDir = pvalDir+sample    
    if not os.path.exists(sampleDir):
        os.makedirs(sampleDir)
    os.chdir(sampleDir)  
    L = len(fc)
    logfc = [numpy.log10(fc[k]) for k in range(L)]
    m = numpy.mean(logfc)
    std = numpy.std(logfc)
    zScores = [(logfc[k]-m)/std for k in range(L)]
    sig = [significant[k] for k in range(L)]
    sig.sort()
    S = list(sig).index(True)
    if ScreenType == 'enrichment':    
        zScores.sort();
    elif ScreenType == 'depletion':
        zScores.sort(reverse=True);
    fig, ax = plt.subplots(figsize=(5,4))
    plt.scatter(range(1,S+1),zScores[0:S],s=4,lw=0,color='grey')
    plt.scatter(range(S+1,L+1),zScores[S:],s=4,lw=0,color='green',label='Significant')
    plt.plot((0,L), (0,0), ls="--", color=(51/255,153/255,1))
    ymax = 1.05*max(zScores); ymin = 1.05*min(zScores)    
    plt.ylim([ymin,ymax]) 
    formatter = FuncFormatter(kilos)
    ax.xaxis.set_major_formatter(formatter)        
    plt.xlabel('ranked sgRNAs', fontsize=12)
    plt.ylabel('z-Score', fontsize=12)
    plt.title(sample+' '+ScreenType+' '+' (fold change)',fontsize=14)
    plt.legend(loc='upper left', prop={'size':10})
    plt.tight_layout()
    plt.savefig(sample+'_'+'sgRNA_zScores.png', dpi=res)    
    if svg:
        plt.savefig(sample+'_'+'sgRNA_zScores.svg')     