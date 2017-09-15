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

def pvalHist(NBpval,NBpval_0,pvalDir,sample,res,svg):
    sampleDir = pvalDir+sample    
    if not os.path.exists(sampleDir):
        os.makedirs(sampleDir)
    os.chdir(sampleDir)      
    bin_size = 0.1; min_edge = 0; max_edge = 1
    n = (max_edge-min_edge)/bin_size; Nplus1 = n + 1
    bin_list = numpy.linspace(min_edge, max_edge, Nplus1)    
    fig, ax = plt.subplots(figsize=(4.2,3.5))
    plt.hist(NBpval,bin_list,color='#c8d1ca',label='Unadjusted')
    plt.hist(NBpval_0,bin_list,color='#0be52c',rwidth=0.8,alpha=0.75,label='Adjusted')
    plt.xticks([0,.2,.4,.6,.8,1])
    plt.xlabel('p-value', fontsize=12)    
    plt.ylabel('Frequency', fontsize=12) 
    formatter = FuncFormatter(kilos)
    ax.yaxis.set_major_formatter(formatter)
    plt.tick_params(labelsize=12)    
    plt.title('sgRNA Significance', fontsize=14) 
    plt.legend(loc='upper center', prop={'size':9})    
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
    fig, ax = plt.subplots(figsize=(4.2,3.5))
    plt.hist(pval_list,bin_list,color='#c8d1ca',label='Unadjusted')
    plt.hist(pval_list0,bin_list,color='#0be52c',rwidth=0.75,alpha=0.75,label='Adjusted')
    plt.xticks([0,.2,.4,.6,.8,1])
    plt.xlabel('p-value', fontsize=12)    
    plt.ylabel('Frequency', fontsize=12) 
    formatter = FuncFormatter(kilos)
    ax.yaxis.set_major_formatter(formatter)    
    plt.tick_params(labelsize=12)
    plt.title('Gene Significance', fontsize=14) 
    plt.legend(loc='upper center', prop={'size':9})    
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
    plt.figure(figsize=(4,3.5))
    if min(len(neglogm_sig),len(neglogp_sig))>100:
        tpcy = 0.35
    else:
        tpcy = 1    
    plt.scatter(neglogm,neglogp,s=8,facecolor='black',lw=0,alpha=0.25)  
    plt.scatter(neglogm_sig,neglogp_sig,s=8,facecolor='green',lw=0,alpha=tpcy,label='Significant')    
    xmax = max(max(neglogm_sig),max(neglogm))
    ymax = max(max(neglogp_sig),max(neglogp))
    plt.xlim([0,1.05*xmax]); plt.ylim([0,1.05*ymax])
    plt.xlabel('-log10 '+GeneMetric, fontsize=12)    
    plt.ylabel('-log10 p-value', fontsize=12) 
    plt.tick_params(labelsize=12)
    plt.title(GeneMetric+' Metric', fontsize=14) 
    leg = plt.legend(loc='lower right', prop={'size':9})
    for lh in leg.legendHandles: lh.set_alpha(1)
    plt.tight_layout()
    plt.savefig(sample+'_'+GeneMetric+'_Metric.png', dpi=res)   
    
    
def VolcanoPlot(fc,NBpval,significant,pvalDir,ScreenType,sample,res,svg,alpha): 
    sampleDir = pvalDir+sample    
    if not os.path.exists(sampleDir):
        os.makedirs(sampleDir)
    os.chdir(sampleDir)  
    L = len(fc)
    logfc = [numpy.log2(fc[k]) for k in range(L) if significant[k]==False]
    neglogp2 = [-numpy.log10(NBpval[k]) for k in range(L) if significant[k]==False]
    logfc_sig = [numpy.log2(fc[k]) for k in range(L) if significant[k]==True]
    neglogp2_sig = [-numpy.log10(NBpval[k]) for k in range(L) if significant[k]==True]    
    plt.figure(figsize=(4,3.5))
    if len(logfc_sig)>100:
        tpcy = 0.35
    else:
        tpcy = 1
    plt.scatter(logfc,neglogp2,s=8,facecolor='grey',lw=0,alpha=0.25)
    plt.scatter(logfc_sig,neglogp2_sig,s=8,facecolor='green',lw=0,alpha=tpcy,label='p < '+str(alpha))    
    if len(logfc_sig)>0 and len(neglogp2_sig)>0:
        xmin = min(min(logfc_sig),min(logfc))        
        xmax = max(max(logfc_sig),max(logfc))
        ymax = max(max(neglogp2_sig),max(neglogp2))
        plt.xlim([0.95*xmin,1.05*xmax]); plt.ylim([0,1.05*ymax])
    plt.xlabel('log2 Fold-Change', fontsize=12)    
    plt.ylabel('-log10 p-value', fontsize=12) 
    plt.tick_params(labelsize=12)
    plt.title('sgRNA Volcano Plot', fontsize=14) 
    if len(logfc_sig)>0:
        leg = plt.legend(loc='lower right', prop={'size':9})
        for lh in leg.legendHandles: lh.set_alpha(1)
    plt.tight_layout()
    plt.savefig(sample+'_'+'sgRNA_volcano.png', dpi=res)


def QQPlot(NBpval,significant,pvalDir,sample,res,svg,alpha): 
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
    if True in sig:
        S = list(sig).index(True)
    else:
        S = L
    if len(neglogp[S:])>100:
        tpcy = 0.35
    else:
        tpcy = 1    
    plt.figure(figsize=(4,3.5))
    plt.scatter(neglogpExp[0:S],neglogp[0:S],s=8,facecolor='grey',lw=0,alpha=0.25)
    plt.scatter(neglogpExp[S:],neglogp[S:],s=8,facecolor='green',lw=0,alpha=tpcy,label='p < '+str(alpha))
    xmax = 1.05*max(neglogpExp)
    ymax = 1.05*max(neglogp)
    plt.xlim([0,xmax]); plt.ylim([0,ymax])
    plt.plot((0,xmax), (0,xmax), ls="--", color=(51/255,153/255,1))
    plt.xlabel('-log10 Expected p-value', fontsize=12)    
    plt.ylabel('-log10 Observed p-value', fontsize=12) 
    plt.tick_params(labelsize=12)
    plt.title('sgRNA QQ Plot', fontsize=14) 
    if True in sig:    
        leg = plt.legend(loc='lower right', prop={'size':9})
        for lh in leg.legendHandles: lh.set_alpha(1)
    plt.tight_layout()
    plt.savefig(sample+'_'+'sgRNA_QQ.png', dpi=res)    

def zScorePlot(fc,significant,pvalDir,ScreenType,sample,res,svg,alpha): 
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
    if True in sig:
        S = list(sig).index(True)
    else:
        S = L    
    if ScreenType == 'enrichment':    
        zScores.sort();
    elif ScreenType == 'depletion':
        zScores.sort(reverse=True);
    fig, ax = plt.subplots(figsize=(5,3.5))
    plt.scatter(range(1,S+1),zScores[0:S],s=8,lw=0,color='#d3d3d3')
    plt.scatter(range(S+1,L+1),zScores[S:],s=8,lw=0,color='green',label='p < '+str(alpha))
    plt.plot((0,L), (0,0), ls="--", color=(51/255,153/255,1))
    ymax = 1.05*max(zScores); ymin = 1.05*min(zScores)    
    plt.ylim([ymin,ymax]) 
    formatter = FuncFormatter(kilos)
    ax.xaxis.set_major_formatter(formatter)        
    plt.xlabel('Ranked sgRNAs', fontsize=12)
    plt.ylabel('z-Score', fontsize=12)
    plt.tick_params(labelsize=12)
    plt.title('sgRNA '+ScreenType+' z-Scores',fontsize=14)
    if True in sig:
        if ScreenType == 'enrichment':
            leg = plt.legend(loc='upper left', prop={'size':9})
            for lh in leg.legendHandles: lh.set_alpha(1)
        elif ScreenType == 'depletion':
            leg = plt.legend(loc='upper right', prop={'size':9})        
            for lh in leg.legendHandles: lh.set_alpha(1)
    plt.tight_layout()
    plt.savefig(sample+'_'+'sgRNA_zScores.png', dpi=res)     