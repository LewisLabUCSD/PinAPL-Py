#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 15:10:27 2016

@author: philipp
"""
# Estimate mean counts and variance from control samples
# =======================================================================
# Imports
from __future__ import division # floating point division by default
import numpy
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import pandas as pd
import os
import re
import pandas
import sys
import time
import yaml
import scipy
from scipy import stats

def EstimateControlCounts(): 
    # ------------------------------------------------
    # Print header
    # ------------------------------------------------
    print('++++++++++++++++++++++++++++++++++++++++++++++++')
    start_total = time.time()    
       
    # ------------------------------------------------
    # Get parameters
    # ------------------------------------------------
    configFile = open('configuration.yaml','r')
    config = yaml.load(configFile)
    configFile.close()
    delta = config['delta']
    ScriptsDir = config['ScriptsDir']
    WorkingDir = config['WorkingDir']
    AlnQCDir = config['AlnQCDir']
    sgRNAReadCountDir = config['sgRNAReadCountDir']
    ControlDir = config['ControlDir']
    res = config['dpi']
    svg = config['svg']
    p_overdisp = config['p_overdisp']
    CtrlCounts_Filename = 'Control_GuideCounts_normalized.txt'
   
    # --------------------------------    
    # Generate table of control counts
    # --------------------------------    
    print('Reading sgRNA counts ...')    
    os.chdir(sgRNAReadCountDir)
    ControlSamples = [d for d in os.listdir(sgRNAReadCountDir) if 'Control' in d and \
        'normalized' in d and '_avg' not in d]
    CountFile = pd.read_table(ControlSamples[0], sep = '\t', names=['sgID','gene','counts'])
    sgIDs = list(CountFile['sgID'])
    genes = list(CountFile['gene'])
    L = len(CountFile)
    CtrlCounts_df = pd.DataFrame()
    CtrlCounts_df['sgID'] = sgIDs
    CtrlCounts_df['gene'] = genes
    if len(ControlSamples) == 0:
        print('### ERROR: No control samples found! ###')
    else:
        for filename in ControlSamples:
            controlsample = filename[0:-27]
            CountFile = pd.read_table(filename, sep='\t',names=['sgID','gene','counts'])
            counts = list(CountFile['counts'])
            CtrlCounts_df[controlsample] = counts
    
    # --------------------------------------------    
    # Compute sample means and sample variance
    # --------------------------------------------  
    print('Estimating mean and variance ...')    
    CtrlCounts_array = CtrlCounts_df.iloc[0:,2:]
    CtrlCounts_matrix = numpy.asmatrix(CtrlCounts_array)    
    Mean_matrix = numpy.mean(CtrlCounts_matrix,axis=1)
    Mean_array = numpy.array(Mean_matrix.T)[0]
    Mean = list(Mean_array)   
    Var_matrix = numpy.var(CtrlCounts_matrix,axis=1)
    Var_array = numpy.array(Var_matrix.T)[0]
    SampleVar = list(Var_array)
    
    # --------------------------------------------------------------    
    # Determine if the variance equals the mean (Poisson distribution)
    # --------------------------------------------------------------       
    if len(ControlSamples) == 1:
        Model = 'none'
        print('WARNING: No control replicates! Cannot choose statistical model.')
    else:
        I = [i for i in range(L) if Mean[i]>0]        
        Mean0 = [Mean[i] for i in I]
        Var0 = [SampleVar[i] for i in I]
        TestStat = scipy.stats.wilcoxon(Var0,Mean0)
        if TestStat[1] >= p_overdisp:
            Model = 'Poisson'
            print('Cannot reject equality of read count mean and variance (p='+str(p_overdisp)+'). Choosing Poisson model ...')
        else:
            # compute rank sums manually (** scipy does not allow one-sided Wilcoxon tests **)
            I = [i for i in range(len(Mean0)) if Var0[i]!=Mean0[i]]
            Mean00 = [Mean0[i] for i in I]
            Var00 = [Var0[i] for i in I]
            Delta = [numpy.abs(Var00[i]-Mean00[i]) for i in range(len(Mean00))]
            sig = [1 if Var00[i]>Mean00[i] else -1 for i in range(len(Mean00))]
            Ranks = scipy.stats.mstats.rankdata(Delta)
            Ranks_pos = [Ranks[i] for i in range(len(Mean00)) if sig[i]>0]
            Ranks_neg = [Ranks[i] for i in range(len(Mean00)) if sig[i]<0]
            W_pos = sum(Ranks_pos)
            W_neg = sum(Ranks_neg)
            if W_pos > W_neg:
                Model = 'Neg. Binomial'
                print('Overdispersion detected at p='+str(TestStat[1])+'. Choosing negative binomial model ...')
            else:
                Model = 'Neg. Binomial'             # for lack of better choice...
                print('WARNING: Low variance in control samples (underdispersion)!')            


    # -----------------------------------------------    
    # Model variance
    # -----------------------------------------------
    if Model == 'none':
        Var = [0 for k in range(len(SampleVar))]    
        n = 'N/A'
        p = 'N/A'        
    elif Model == 'Neg. Binomial':
        x = [numpy.log(Mean[k]) for k in range(L) if Mean[k]>0 and SampleVar[k]>Mean[k]]
        y = [numpy.log(SampleVar[k]-Mean[k]) for k in range(L) if Mean[k]>0 and SampleVar[k]>Mean[k]]    
        c = [y[k]-2*x[k] for k in range(len(x))]
        c_0 = numpy.mean(c)
        D = numpy.exp(c_0)
        Var = [Mean[k] + D*Mean[k]**2 for k in range(L)]  
        # -----------------------------------------------    
        # Compute parameters for neg. binom. distribution 
        # n: number of failures, p: probability of failure
        # -----------------------------------------------
        print('Computing distribution parameters ...')
        n = list(); p = list()
        for k in range(L):
            if Mean[k]==0 or Var[k]==0 :
                n.append(((Mean[k]+delta)**2/(Var[k]+2*delta))/(1-(Mean[k]+delta)/(Var[k]+2*delta)))
                p.append((Mean[k]+delta)/(Var[k]+2*delta))
            else:
                n.append((Mean[k]**2/Var[k])/(1-Mean[k]/Var[k]))
                p.append(Mean[k]/Var[k])
    elif Model == 'Poisson':
        Var = [Mean[k] if Mean[k]>0 else 1 for k in range(L)]
        n = 'N/A'
        p = 'N/A'
                
    # --------------------------------    
    # Write data frame
    # --------------------------------     
    if not os.path.exists(ControlDir):
        os.makedirs(ControlDir)         
    os.chdir(ControlDir)    
    print('Writing dataframe ...')
    CtrlCounts_df['Model'] = Model
    CtrlCounts_df['Mean'] = Mean
    CtrlCounts_df['Sample Variance'] = SampleVar
    CtrlCounts_df['Model Variance'] = Var
    CtrlCounts_df['n'] = n
    CtrlCounts_df['p'] = p
    CtrlCounts_df.to_csv(CtrlCounts_Filename,sep='\t')    

    # --------------------------------    
    # Plots
    # -------------------------------- 
    plt.figure(figsize=(6.5,2.9))
    # Mean/Variance plot
    print('Generating dispersion plot ...')
    plt.subplot(121)        
    if max(SampleVar) > 0:
        Mmax = numpy.percentile(Mean_array,99)
        x = [Mean[k] for k in range(L) if Mean[k] < Mmax]
        y = [SampleVar[k] for k in range(L) if Mean[k] < Mmax]
        plt.scatter(x,y,s=4,lw=0,alpha=0.25,rasterized=True)
        plt.plot(x,x,'--',color='orange',label='Mean = Variance')
        leg = plt.legend(loc='upper left', prop={'size':8})
        for lh in leg.legendHandles: lh.set_alpha(1)
    else: # no control replicates
        plt.figtext(0.25,0.5,'N/A')
    plt.xlabel('Mean', fontsize=11)    
    plt.ylabel('Variance', fontsize=11)    
    plt.tick_params(labelsize=11)
    plt.title('Read Count Overdispersion', fontsize=12)
    # Log Plot with Regression
    print('Generating log regression plot ...')
    plt.subplot(122)
    if Model == 'Neg. Binomial' and max(SampleVar)>0:
        logx = [numpy.log(Mean[k]) for k in range(L) if Mean[k]>0 and SampleVar[k]>Mean[k]]
        logy = [numpy.log(SampleVar[k]-Mean[k]) for k in range(L) if Mean[k]>0 and SampleVar[k]>Mean[k]]
        plt.scatter(logx,logy,s=4,lw=0,alpha=0.25,rasterized=True)  
        logy_0 = [2*logx[k] + c_0 for k in range(len(logx))]
        plt.plot(logx,logy_0,'r--')
        Disp = '%.2f' % D
        plt.figtext(0.62,0.8,'Var = Mean+'+Disp+'Mean2',color='red',fontsize=8)
    else: # no control replicates
        plt.figtext(0.75,0.5,'N/A')
    plt.xlabel('log (Mean)', fontsize=11)    
    plt.ylabel('log (Variance - Mean)', fontsize=11)     
    plt.title('Mean/Variance Model', fontsize=12)
    plt.tick_params(labelsize=11)
    plt.tight_layout()
    plt.savefig('Control_MeanVariance.png',dpi=res)
    if svg:
        plt.savefig('Control_MeanVariance.svg')      
    

    # --------------------------------------
    # Final time stamp
    # --------------------------------------        
    os.chdir(ScriptsDir)
    end_total = time.time()
    # Final time stamp
    print('------------------------------------------------')
    print('Script completed.')    
    sec_elapsed = end_total - start_total
    if sec_elapsed < 60:
        time_elapsed = sec_elapsed
        print('Time elapsed (Total) [secs]: ' + '%.3f' % time_elapsed +'\n')
    elif sec_elapsed < 3600:
        time_elapsed = sec_elapsed/60
        print('Time elapsed (Total) [mins]: ' + '%.3f' % time_elapsed +'\n')
    else:
        time_elapsed = sec_elapsed/3600
        print('Time elapsed (Total) [hours]: ' + '%.3f' % time_elapsed +'\n')            


if __name__ == "__main__":
    EstimateControlCounts()
