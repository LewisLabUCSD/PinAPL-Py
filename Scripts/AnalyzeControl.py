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
import glob
import sys
import time
import yaml

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
    ControlDir = config['ControlDir']
    res = config['dpi']
    thr_overdisp = config['thr_overdisp']
    CtrlCounts_Filename = 'Control_GuideCounts_0.tsv'
   
    # --------------------------------    
    # Generate table of control counts
    # --------------------------------    
    print('Reading control counts ...')    
    os.chdir(AlnQCDir)
    ControlSamples = [d for d in os.listdir(AlnQCDir) if 'Control' in d]
    os.chdir(ControlSamples[0])
    colnames = ['sgID','gene','counts']                      
    CountFile = pd.read_table(glob.glob('*GuideCounts_0.tsv')[0], sep='\t',names=colnames)
    sgIDs = list(CountFile['sgID'].values)
    genes = list(CountFile['gene'].values)
    L = len(sgIDs)
    CtrlCounts_df = pd.DataFrame(data = {'sgID': sgIDs,
                                    'gene': genes},
                            columns = ['sgID','gene'])        
    if len(ControlSamples) == 0:
        print('### ERROR: No control samples found! ###')
    else:
        os.chdir(AlnQCDir)
        for controlsample in ControlSamples:
            os.chdir(controlsample)
            filename = glob.glob('*GuideCounts_0.tsv')[0]                          
            CountFile = pd.read_table(filename, sep='\t',names=colnames)
            counts = list(CountFile['counts'].values)
            CtrlCounts_df[controlsample] = counts
            os.chdir(AlnQCDir)
    
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
    Svar0 = numpy.mean(SampleVar)
    if Svar0 == 0:
        Model = 'none'
        print('WARNING: Zero variance or no control replicates! Cannot choose statistical model.')
    else:
        L0_list = [1 if Mean[k]>0 else 0 for k in range(L)]
        overdisp_list = [1 if Mean[k]>0 and SampleVar[k]>Mean[k] else 0 for k in range(L)]
        overdisp = sum(overdisp_list)/sum(L0_list)
        print('Overdispersion fraction: '+str(overdisp))
        if overdisp >= thr_overdisp:
            Model = 'Neg. Binomial'
            print('Choosing negative binomial model ...')
        else:
            Model = 'Poisson'
            print('Choosing Poisson model ...')

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
        print('Computing parameters of negative binomial distribution ...')
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
    plt.figure(figsize=(8,3.33))
    # Mean/Variance plot
    print('Generating dispersion plot ...')
    plt.subplot(121)        
    if max(SampleVar) > 0:
        Mmax = numpy.percentile(Mean_array,99)
        x = [Mean[k] for k in range(L) if Mean[k] < Mmax]
        y = [SampleVar[k] for k in range(L) if Mean[k] < Mmax]
        plt.scatter(x,y,s=4,lw=0,alpha=0.25)
        plt.plot(x,x,'--',color='orange',label='Mean = Variance')
        leg = plt.legend(loc='upper left', prop={'size':8})
        for lh in leg.legendHandles: lh.set_alpha(1)
    else: # no control replicates
        plt.figtext(0.25,0.5,'N/A')
    plt.xlabel('Mean', fontsize=12)    
    plt.ylabel('Variance', fontsize=12)    
    plt.tick_params(labelsize=12)
    plt.title('Read Count Overdispersion', fontsize=14)
    # Log Plot with Regression
    print('Generating log regression plot ...')
    plt.subplot(122)
    if Model == 'Neg. Binomial' and max(SampleVar)>0:
        logx = [numpy.log(Mean[k]) for k in range(L) if Mean[k]>0 and SampleVar[k]>Mean[k]]
        logy = [numpy.log(SampleVar[k]-Mean[k]) for k in range(L) if Mean[k]>0 and SampleVar[k]>Mean[k]]
        plt.scatter(logx,logy,s=4,lw=0,alpha=0.25)  
        logy_0 = [2*logx[k] + c_0 for k in range(len(logx))]
        plt.plot(logx,logy_0,'r--')
        Disp = '%.2f' % D
        plt.figtext(0.62,0.8,'Var = Mean+'+Disp+'Mean2',color='red',fontsize=8)
    else: # no control replicates
        plt.figtext(0.75,0.5,'N/A')
    plt.xlabel('log (Mean)', fontsize=12)    
    plt.ylabel('log (Variance - Mean)', fontsize=12)     
    plt.title('Mean/Variance Model', fontsize=14)
    plt.tick_params(labelsize=12)
    plt.tight_layout()
    plt.savefig('Control_MeanVariance.png',dpi=res)
    

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
