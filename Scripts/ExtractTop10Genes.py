#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 18:34:09 2019

@author: philipp
"""
import sys
import os
import yaml
import pandas
import glob

def printTop10Genes(sample):
    configFile = open('configuration.yaml','r')
    config = yaml.load(configFile)
    configFile.close()
    GeneDir = config['GeneDir']
    ScriptsDir = config['ScriptsDir']
    os.chdir(GeneDir)
    filename = glob.glob(sample+'*_GeneList.txt')[0]
    GeneRanking = pandas.read_table(filename, sep='\t')
    sig = list(GeneRanking['significant'])
    if True in sig:
        S = sig.index(False)
        S0 = min(S,10)        
    else:
        S0 = 0
    GeneRankingTop10 = GeneRanking[0:S0] 
    print('-------------------------------------------------------------------------------------')    
    print(sample + ': Gene Ranking Results')
    print('-------------------------------------------------------------------------------------')    
    if S0>0:
        print(GeneRankingTop10)
    else:
        print('### No significant genes found ###')
    print('\n')
    os.chdir(ScriptsDir)

def GetTop10Genes(sample):
    configFile = open('configuration.yaml','r')
    config = yaml.load(configFile)
    configFile.close()
    GeneDir = config['GeneDir']
    ScriptsDir = config['ScriptsDir']
    os.chdir(GeneDir)
    filename = glob.glob(sample+'*_GeneList.txt')[0]
    GeneRanking = pandas.read_table(filename, sep='\t')
    sig = list(GeneRanking['significant'])
    if True in sig:
        S = sig.index(False)
        S0 = min(S,10)        
    else:
        S0 = 0
    GeneRankingTop10 = GeneRanking[0:S0]     
    Top10Genes = list(GeneRankingTop10['gene'])
    os.chdir(ScriptsDir)
    return Top10Genes


if __name__ == "__main__":
    input1 = sys.argv[1]    
    printTop10Genes(input1)    