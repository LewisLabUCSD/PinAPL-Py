#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 15:10:27 2016

@author: philipp
"""
# Estimate mean counts and variance from control samples
# =======================================================================
# Imports
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
    print('*************************************************************')
    print('PinAPL-Py: Control Count Estimation')
    print('*************************************************************')
    start_total = time.time()    
       
    # ------------------------------------------------
    # Get parameters
    # ------------------------------------------------
    os.chdir('/workingdir/')
    configFile = open('configuration.yaml','r')
    config = yaml.load(configFile)
    configFile.close()
    WorkingDir = config['WorkingDir']
    LibDir = config['LibDir']
    LibFilename = config['LibFilename']    
    CtrlPrefix = config['CtrlPrefix']
    AnalysisDir = config['AnalysisDir']
    QCDir = config['QCDir']
    ControlDir = AnalysisDir + 'Control/'
    CtrlCounts_Filename = 'Control_GuideCounts_0.tsv'

    # ------------------------------------------------
    # Get library entries
    # ------------------------------------------------
    os.chdir(LibDir)
    LibCols = ['gene','ID','seq']
    LibFile = pd.read_table(LibFilename, sep = '\t', skiprows = 1, names = LibCols)
    sgIDs = list(LibFile['ID'].values)
    genes = list(LibFile['gene'].values)
    L = len(sgIDs)
    CtrlCounts_df = pd.DataFrame(data = {'sgID': sgIDs,
                                    'gene': genes},
                            columns = ['sgID','gene'])
    
    # --------------------------------    
    # Generate table of control counts
    # --------------------------------    
    print('Reading control counts ...')    
    os.chdir(QCDir)
    ControlSamples = [d for d in os.listdir(QCDir) if CtrlPrefix in d]
    if len(ControlSamples) == 0:
        print('ERROR: No control sample directories found!')
    else:
        colnames = ['sgID','gene','counts']        
        for controlsample in ControlSamples:
            os.chdir(controlsample)
            filename = glob.glob('*GuideCounts_0.tsv')[0]                          
            CountFile = pd.read_table(filename, sep='\t',names=colnames)
            counts = list(CountFile['counts'].values)
            CtrlCounts_df[controlsample] = counts
            os.chdir(QCDir)
    
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
    Var = list(Var_array)    

    # -----------------------------------------------    
    # Estimate variance from negative binomial model
    # -------------------------------- ---------------
    x = [numpy.log(Mean[k]) for k in range(L) if Mean[k]>0 and Var[k]>Mean[k]]
    y = [numpy.log(Var[k]-Mean[k]) for k in range(L) if Mean[k]>0 and Var[k]>Mean[k]]    
    c = [y[k]-2*x[k] for k in range(len(x))]
    c_0 = numpy.mean(c)
    D = numpy.exp(c_0)
    Var_Model = [Mean[k] + D*Mean[k]**2 for k in range(L)]    
    
    # --------------------------------    
    # Write data frame
    # --------------------------------     
    if not os.path.exists(ControlDir):
        os.makedirs(ControlDir)         
    os.chdir(ControlDir)    
    print('Writing dataframe ...')
    CtrlCounts_df['Mean'] = Mean
    CtrlCounts_df['Sample Variance'] = Var
    CtrlCounts_df['Model Variance'] = Var_Model 
    CtrlCounts_df.to_csv(CtrlCounts_Filename,sep='\t')    

    # --------------------------------    
    # Plots
    # -------------------------------- 
    plt.figure(figsize=(12,5))
    # Mean/Variance plot
    print('Generating dispersion plot ...')
    Mmax = numpy.percentile(Mean_array,99)
    x = [Mean[k] for k in range(L) if Mean[k] < Mmax]
    y = [Var[k] for k in range(L) if Mean[k] < Mmax]
    plt.subplot(121)
    plt.scatter(x,y,s=2)
    plt.xlabel('Mean', fontsize=14)    
    plt.ylabel('Variance', fontsize=14)
    plt.plot(x,x,'g--')
    plt.title('Read Counts (Control Samples)', fontsize=16)
    # Log Plot with Regression
    print('Generating log regression plot ...')
    plt.subplot(122)
    logx = [numpy.log(Mean[k]) for k in range(L) if Mean[k]>0 and Var[k]>Mean[k]]
    logy = [numpy.log(Var[k]-Mean[k]) for k in range(L) if Mean[k]>0 and Var[k]>Mean[k]]
    plt.scatter(logx,logy,s=2)
    plt.xlabel('log (Mean)', fontsize=14)    
    plt.ylabel('log (Variance - Mean)', fontsize=14)    
    logy_0 = [2*logx[k] + c_0 for k in range(len(logx))]
    plt.plot(logx,logy_0,'r--')
    plt.title('log Regression', fontsize=16)
    plt.savefig('Control_MeanVariance.png')

    # --------------------------------------
    # Final time stamp
    # --------------------------------------        
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
