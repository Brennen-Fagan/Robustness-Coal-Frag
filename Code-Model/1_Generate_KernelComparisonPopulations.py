# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 15:21:39 2023

@author: btf500
"""

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("start")
parser.add_argument("end")
args = parser.parse_args()

from counting import random_partition_Probabilistic_Divide_Conquer_1 as StartPosition
from os import makedirs
from os import getcwd
from os.path import isdir
from datetime import date
from scipy.io import savemat
import numpy as np
from numpy.random import seed as np_seed
from random import seed as py_seed

#System Settings
from JohnsonOnePopModels import MixedKernelModel as Run
Slices = 10      #Number of sub-files.
EndTime = 200000 #Total run time.
BrnTime = 10**6  #Amount of time to burn-in for.
Sampling = 10    #Once per _ steps.
Slices = 10
Barrier = 0     #Fragmentation Barrier. Default = 0.

population_sets = 1
system_sets = 4

seed_val = int(date.today().toordinal())
seed_val = 10000000 * np.random.rand(population_sets*system_sets*4) # 4 reps.
    
#Expecting j's of 0 -> 3, which means task IDs of (1 -> 4)*4
for j in range(int(args.start), int(args.end)):
    
    Tag = 'PFragv03'
    pfrag = 0.3 # 0.3 * 1e4 / 0.7 should be no cycles region.
    
    if (j%system_sets == 0):
        Label = 'MixedKernel_CC'
        CKernel = 'Constant'
        FKernel = 'Constant'
        FDistri = 'Total'
    elif (j%system_sets == 1):
        Label = 'MixedKernel_MC'
        CKernel = 'Multiplicative'
        FKernel = 'Constant'
        FDistri = 'Total'
    elif (j%system_sets == 2):
        Label = 'MixedKernel_CM'
        CKernel = 'Constant'
        FKernel = 'Multiplicative'
        FDistri = 'Total'
    elif (j%system_sets == 3):
        Label = 'MixedKernel_MM'
        CKernel = 'Multiplicative'
        FKernel = 'Multiplicative'
        FDistri = 'Total'

    Pops = list()
    Pops.append(10000)

    DirectoryName = (getcwd() + '/Data' + date.today().strftime('%y-%m-%d') +
                        '/' + Label + Tag + '/')
    if not isdir(DirectoryName):
        makedirs(DirectoryName)
        
    print(int(seed_val[j]))
    py_seed(int(seed_val[j]))
    np_seed(int(seed_val[j]))

    for i in Pops:
        Start = StartPosition(i)
        #Get Start significantly closer to the equilibrium "state" before recording.
        if 'alph' not in locals():
            Start = Run(Population = Start,
                              EndTime = BrnTime,
                              PCoal = 1 - pfrag, PFrag = -1,
                              History = 0,
                              FragmentKernel = FKernel,
                              FragmentDistribution = FDistri,
                              FragmentBarrier = Barrier,
                              CoalescenceKernel = FKernel,
                              Flux = False,
                              HistorySteps = 1)
        else:
            Start = Run(Population = Start,
                                  EndTime = BrnTime,
                                  PCoal = 1 - pfrag, PFrag = -1,
                                  History = 0, alpha = alph,
                                  FragmentKernel = FKernel,
                                  FragmentDistribution = FDistri,
                                  FragmentBarrier = Barrier,
                                  CoalescenceKernel = FKernel,
                                  Flux = False,
                                  HistorySteps = 1)
        FileName = ('RandomPartition_' + Label + Tag + 
                    '_' +str(int(i)) + '_Run_' + str(int(j)) + '.mat')

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

            Fluxes = list(zip(*History[1]))
            #Following: https://stackoverflow.com/a/29042033
            if(s == 0):
                #Create .mat file
                savemat(DirectoryName + FileName, 
                        mdict = {Label + Tag + '_' + str(int(s)) : np.array(History[0]),
                                 Label + Tag + '_' + str(int(s)) + '_CFlux': np.array(Fluxes[0]),
                                 Label + Tag + '_' + str(int(s)) + '_FFlux': np.array(Fluxes[1])}, 
                        do_compression = True)
            else:
                #append .mat file
                with open(DirectoryName + FileName, 'ab') as f:
                    savemat(f, 
                            mdict = {Label + Tag + '_' + str(int(s)) : np.array(History[0]),
                                     Label + Tag + '_' + str(int(s)) + '_CFlux': np.array(Fluxes[0]),
                                     Label + Tag + '_' + str(int(s)) + '_FFlux': np.array(Fluxes[1])}, 
                            do_compression = True)
            
            Start = History[0][-1]
