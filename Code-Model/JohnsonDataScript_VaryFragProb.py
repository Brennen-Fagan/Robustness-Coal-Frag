# -*- coding: utf-8 -*-
"""
Created on Thursday July 7 2018

@author: btf500
"""

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("start")
parser.add_argument("end")
args = parser.parse_args()

from counting import random_partition_large as StartPosition
from os import makedirs
from os import getcwd
from os.path import isdir
from datetime import date
from scipy.io import savemat
import numpy as np

Slices = 10

#Expecting j's of 0 -> 9, which means task IDs of 1 -> 10
for j in range(int(args.start), int(args.end)):
    from JohnsonOnePopModels import MixedKernelModel as Run
    Label = 'MixedKernel_MC'
    FKernel = 'Constant'
    FDistri = 'Total'

    if(j%10 == 0):
        Tag = 'Frag_025'
        pfrag = 0.025
    elif(j%10  == 1):
        Tag = 'Frag_050'
        pfrag = 0.05
    elif(j%10 == 2):
        Tag = 'Frag_100'
        pfrag = 0.10
    elif(j%10 == 3):
        Tag = 'Frag_125'
        pfrag = 0.125
    elif(j%10 == 4):
        Tag = 'Frag_150'
        pfrag = 0.150
    elif(j%10 == 5):
        Tag = 'Frag_175'
        pfrag = 0.175
    elif(j%10 == 6):
        Tag = 'Frag_200'
        pfrag = 0.200
    elif(j%10 == 7):
        Tag = 'Frag_225'
        pfrag = 0.225
    elif(j%10 == 8):
        Tag = 'Frag_250'
        pfrag = 0.250
    elif(j%10 == 9):
        Tag = 'Frag_275'
        pfrag = 0.275

    Pops = list()
    #if((j//4)%3 == 0):         
    #        Pops.append(20000)
    #        Pops.append(100)
    #        Pops.append(30000)
    #elif((j//4)%3 == 1):
    #        Pops.append(10000)
    #        Pops.append(200)
    #        Pops.append(300)
    #        Pops.append(70000)
    #elif((j//4)%3 == 2):
    #        Pops.append(7000)
    #        Pops.append(3000)
    #        Pops.append(1000)
    #        Pops.append(2000)
    Pops.append(100000)

    DirectoryName = getcwd() + '/Data' + date.today().strftime('%y-%m-%d') + '/' + Label + Tag + '/'
    if not isdir(DirectoryName):
        makedirs(DirectoryName)

    for i in Pops:
        Start = StartPosition(i)
        #Get Start significantly closer to the equilibrium "state" before recording.
        Start = Run(Population = Start,
                          EndTime = 1000000,
                          PCoal = 1 - pfrag, PFrag = -1,
                          History = 0,
                          FragmentKernel = FKernel,
                          FragmentDistribution = FDistri,
                          HistorySteps = 1)

        FileName = 'RandomPartition_' + Label + Tag + '_' +str(int(i)) + '_Run_' + str(int(j)) + '.mat'

        for s in range(Slices):          
            History = Run(Population = Start,
                              EndTime = (200000)//Slices,
                              PCoal = 1 - pfrag, PFrag = -1,
                              History = 1,
                              FragmentKernel = FKernel,
                              FragmentDistribution = FDistri,
                              HistorySteps = 10)

            #Following: https://stackoverflow.com/a/29042033
            if(s == 0):
                #Create .mat file
                savemat(DirectoryName + FileName, 
                        mdict = {Label + Tag + '_' + str(int(s)) : History}, 
                        do_compression = True)
            else:
                #append .mat file
                with open(DirectoryName + FileName, 'ab') as f:
                    savemat(f, 
                            mdict = {Label + Tag + '_' + str(int(s)) : np.array(History)}, 
                            do_compression = True)
            
            Start = History[-1]
                

