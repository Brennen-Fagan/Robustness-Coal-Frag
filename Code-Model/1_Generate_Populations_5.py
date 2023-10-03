# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 10:37:11 2019

@author: btf500
"""

#Python3 Code to generate populations. 
#Second Iteration (different populations).

#Terminal Parsing.
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("start")
parser.add_argument("end")
args = parser.parse_args()

#System imports
from counting import random_partition_large as StartPosition
from os import makedirs
from os import getcwd
from os.path import isdir
from datetime import date
from scipy.io import savemat
import numpy as np

#System Settings
from JohnsonOnePopModels import MixedKernelModel as Run
Slices = 10      #Number of sub-files.
EndTime = 200000 #Total run time.
Sampling = 10    #Once per _ steps.

FKernel = 'Constant'
Barrier = 0     #Fragmentation Barrier. Default = 0.

#Expecting j's of 0 -> 11, since 6 sigma, 2 frag probs, ID: (1->12)*4
population_sets = 4
system_sets = 12

for j in range(int(args.start), int(args.end)):

#    if(j%system_sets == 0):
#        Tag = 'CC_Total'
#        Label = 'MixedKernel_'
#        FDistri = 'Total'
    if(j%system_sets == 0):
        Tag = 'CC_PL195_v001'
        Label = 'MixedKernel_'
        alph = 1.95
        FDistri = 'PL'
        pfrag = 0.01
    elif(j%system_sets == 1):
        Tag = 'CC_PL175_v001'
        Label = 'MixedKernel_'
        alph = 1.75
        FDistri = 'PL'
        pfrag = 0.01
    elif(j%system_sets == 2):
        Tag = 'CC_PL150_v001'
        Label = 'MixedKernel_'
        alph = 1.50
        FDistri = 'PL'
        pfrag = 0.01
    elif(j%system_sets == 3):
        Tag = 'CC_PL125_v001'
        Label = 'MixedKernel_'
        alph = 1.25
        FDistri = 'PL'
        pfrag = 0.01
    elif(j%system_sets == 4):
        Tag = 'CC_PL105_v001'
        Label = 'MixedKernel_'
        alph = 1.05
        FDistri = 'PL'
        pfrag = 0.01
    elif(j%system_sets == 5):
        Tag = 'CC_Total_v001'
        Label = 'MixedKernel_'
        FDistri = 'Total'
        pfrag = 0.01
    elif(j%system_sets == 6):
        Tag = 'CC_PL195_v0001'
        Label = 'MixedKernel_'
        alph = 1.95
        FDistri = 'PL'
        pfrag = 0.001
    elif(j%system_sets == 7):
        Tag = 'CC_PL175_v0001'
        Label = 'MixedKernel_'
        alph = 1.75
        FDistri = 'PL'
        pfrag = 0.001
    elif(j%system_sets == 8):
        Tag = 'CC_PL150_v0001'
        Label = 'MixedKernel_'
        alph = 1.50
        FDistri = 'PL'
        pfrag = 0.001
    elif(j%system_sets == 9):
        Tag = 'CC_PL125_v0001'
        Label = 'MixedKernel_'
        alph = 1.25
        FDistri = 'PL'
        pfrag = 0.001
    elif(j%system_sets == 10):
        Tag = 'CC_PL105_v0001'
        Label = 'MixedKernel_'
        alph = 1.05
        FDistri = 'PL'
        pfrag = 0.001
    elif(j%system_sets == 11):
        Tag = 'CC_Total_v0001'
        Label = 'MixedKernel_'
        FDistri = 'Total'
        pfrag = 0.001
		
    Pops = list()
    if((j//system_sets)%population_sets == 0):         
        Pops.append(30000)
        Pops.append(100)
    elif((j//system_sets)%population_sets == 1):
        Pops.append(10000)
        Pops.append(300)
    elif((j//system_sets)%population_sets == 2):
        Pops.append(3000)
        Pops.append(1000)
    elif((j//system_sets)%population_sets == 3):
        Pops.append(100000)
        
    #Pops.append(100000)

    DirectoryName = getcwd() + '/Data' + date.today().strftime('%y-%m-%d') + '/' + Label + Tag + '/'
    if not isdir(DirectoryName):
        makedirs(DirectoryName)

    for i in Pops:
        Start = StartPosition(i)
        #Get Start significantly closer to the equilibrium "state" before recording.
        if 'alph' not in locals():
            Start = Run(Population = Start,
                              EndTime = 1000000,
                              PCoal = 1 - pfrag, PFrag = -1,
                              History = 0,
                              FragmentKernel = FKernel,
                              FragmentDistribution = FDistri,
                              FragmentBarrier = Barrier,
                              CoalescenceKernel = FKernel,
                              HistorySteps = 1)
        else:
            Start = Run(Population = Start,
                                  EndTime = 1000000,
                                  PCoal = 1 - pfrag, PFrag = -1,
                                  History = 0, alpha = alph,
                                  FragmentKernel = FKernel,
                                  FragmentDistribution = FDistri,
                                  FragmentBarrier = Barrier,
                                  CoalescenceKernel = FKernel,
                                  HistorySteps = 1)
        FileName = 'RandomPartition_' + Label + Tag + '_' +str(int(i)) + '_Run_' + str(int(j)) + '.mat'

        for s in range(Slices):          
            if 'alph' not in locals():
                History = Run(Population = Start,
                                  EndTime = (EndTime)//Slices,
                                  PCoal = 1 - pfrag, PFrag = -1,
                                  History = 1,
                                  FragmentKernel = FKernel,
                                  FragmentDistribution = FDistri,
                                  FragmentBarrier = Barrier,
                                  CoalescenceKernel = FKernel,
                                  HistorySteps = Sampling)
            else:
                History = Run(Population = Start,
                                      EndTime = (EndTime)//Slices,
                                      PCoal = 1 - pfrag, PFrag = -1,
                                      History = 1, alpha = alph,
                                      FragmentKernel = FKernel,
                                      FragmentDistribution = FDistri,
                                      FragmentBarrier = Barrier,
                                      CoalescenceKernel = FKernel,
                                      HistorySteps = Sampling)

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
