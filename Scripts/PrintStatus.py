# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 13:50:52 2017

@author: philipp
"""

# Print status messages
import sys


def PrintStatus_Header():
    print('**************************************************')
    print('Launching PinAPL-Py ...')
    print('P. Spahn et al., UC San Diego (10/2016)')
    print('**************************************************')
    
def PrintStatus_SubHeader(msg):
    print('**************************************************')
    print(msg)
    print('**************************************************')

def PrintStatus_Done(msg):
    print(msg)
    print('\n\n')

def PrintStatus_ProcessSample(sample):
    print('\nProcessing sample '+sample+' ... ')

def PrintStatus_AllDone():
    print('**************************************************')
    print('PinAPL-Py completed.')  

def PrintStatus_TimeStamp(msg):
    print(msg)


   
if __name__ == "__main__":
    input1 = sys.argv[1]
    input2 = sys.argv[2]    
    if input1 == 'Header':
        PrintStatus_Header()
    elif input1 == 'SubHeader':
        PrintStatus_SubHeader(input2)
    elif input1 == 'Done':
        PrintStatus_Done(input2)
    elif input1 == 'ProcessSample':
        PrintStatus_ProcessSample(input2)
    elif input1 == 'AllDone':
        PrintStatus_AllDone()
    elif input1 == 'TimeStamp':
        PrintStatus_TimeStamp(input2)        