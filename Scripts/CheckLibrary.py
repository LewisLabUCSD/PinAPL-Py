#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 10:22:56 2016

@author: philipp
"""
# Library sanity check
# =======================================================================
# Imports
import yaml
import sys
import os
import pandas

def LibrarySanityCheck():
    # ------------------------------------------------
    # Get parameters
    # ------------------------------------------------
    configFile = open('configuration.yaml','r')
    config = yaml.load(configFile)
    configFile.close()
    LibDir = config['LibDir']
    LibFilename = config['LibFilename']
    LibFormat = LibFilename[-3:]
    if LibFormat == 'tsv':
        libsep = '\t'
    elif LibFormat == 'csv':
        libsep = ','

    # ----------------------------------
    # Read library
    # ----------------------------------
    os.chdir(LibDir)
    LibCols = ['gene','ID','seq']
    LibFile = pandas.read_table(LibFilename, sep = libsep, skiprows = 1, names = LibCols)
    GeneNames = list(LibFile['gene'])
    ID = list(LibFile['ID'])
    seq = list(LibFile['seq'])

    # ----------------------------------
    # Replace non-printable characters (...these cause problems in PlotCount.py)
    # ----------------------------------    
    GeneNames0 = []
    ID0 = []
    for gene in GeneNames:
        gene = gene.replace('|','_')
        gene = gene.replace('(','_')
        gene = gene.replace(')','_')
        gene = gene.replace(';','_')
        gene = gene.replace('"','')
        gene = gene.replace('/','_')
        gene = gene.replace('\\','_')
        gene = gene.replace(' ','_')        
        GeneNames0.append(gene)
    for sgRNA in ID:
        sgRNA = sgRNA.replace('|','_')
        sgRNA = sgRNA.replace('(','_')
        sgRNA = sgRNA.replace(')','_')
        sgRNA = sgRNA.replace(';','_')
        sgRNA = sgRNA.replace('"','')
        sgRNA = sgRNA.replace('/','_')
        sgRNA = sgRNA.replace('\\','_')
        sgRNA = sgRNA.replace(' ','_')        
        ID0.append(sgRNA)    
    if GeneNames != GeneNames0 or ID != ID0:
            LibFile0 = pandas.DataFrame(data = {'gene': [gene for gene in GeneNames0],
                                     'ID': [sgRNA for sgRNA in ID0],
                                     'seq': [s for s in seq]},
                            columns = ['gene','ID','seq'])
            LibFile0.to_csv(LibFilename, sep = libsep, index = False)
            print("WARNING: Found non-printable characters in library file. Replaced by '_' ")


   
if __name__ == "__main__":
    LibrarySanityCheck()
