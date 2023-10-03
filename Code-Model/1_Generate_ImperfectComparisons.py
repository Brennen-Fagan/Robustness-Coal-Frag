# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 12:27 2023

@author: btf500
"""

#Python3 Code to generate populations.

#%% Terminal Parsing.
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("start")
parser.add_argument("end")
args = parser.parse_args()
# class Dummy:
#     def __init__(self, start, end):
#         self.start = start
#         self.end = end
# args = Dummy(0, 1)

#%% System imports
from counting import random_partition_large as StartPosition
from os import makedirs
from os import getcwd
from os.path import isdir
from datetime import date
from scipy.io import savemat
import numpy as np

from random import seed as rseed
from numpy.random import seed as npseed
from functools import partial # Trying to make this easier to run with variety of functions.

#%% System Settings
# Importing different model variations.
# from JohnsonOnePopModels import MixedKernelModel as Run
Slices = 10      #Number of sub-files.
EndTime = 200000 #Total run time.
Sampling = 10    #Once per _ steps.

# New: need for 0.30 and 10000, already have for 0.01, 10000.
# See: YARCC_CoFr_Analysis_19-02-04_19-02-20.7z\YARCC_CoFr_Analysis_19-02-04_19-02-20\YARCC_CoFr_Analysis_19_02_20\Data19-02-20\
#      for 7 below at 0.01 at 10,000.
# See: YARCC-Scratch-Data.7z\YARCC-Scratch-Data\YARCC_Python_18_04_28\
#      for an assortment of the below at 0.01 with various populations.
pfrag = 0.30
Pops = [10000]

# 6 fragmentations + 1 coalescence, but 1 fragmentation has two parameters
# So 9 combinations in total ((2,4), (3,3), and (4, 2)).
# Then there is the combined coal-frag (imperfect all).
# Then we add on the 110, 120, 130, 140, 160, 170, 180, and 190 high pFrag
# cases. So 9 + 1 + 8 = 18
#  1. imperfectFrag BetaBinomial ((2,4), (3,3), and (4, 2))
#  2. imperfectFrag Partitioning
#  3. imperfectFrag Halving
#  4. imperfectFrag Halving+Rel.
#  5. imperfectFrag Decimating
#  6. imperfectFrag Antidecimating
#  7. imperfectFrag CRP (110, 120, 130, 140, 160, 170, 180, and 190)
#  8. imperfectCoal (Additive)
#  9. imperfectAll  (Additive)


system_sets = 18

#%% Main
for j in range(int(args.start), int(args.end)):
    #%% Initialise Randomness
    theSeed = np.random.randint(1e8)
    rseed(theSeed); npseed(theSeed)

    #%% Setup the run
    Label = 'SteadyState_'
    if(j%system_sets == 0):
        Tag = 'BBa2b4_v030'
        from JohnsonOnePopModels import OnePopImperfectFrag as Engine
        Run = partial(Engine, a = 2, b = 4, PCoal = 1 - pfrag, PFrag = pfrag)
    elif(j%system_sets == 1):
        Tag = 'BBa3b3_v030'
        from JohnsonOnePopModels import OnePopImperfectFrag as Engine
        Run = partial(Engine, a = 3, b = 3, PCoal = 1 - pfrag, PFrag = pfrag)
    elif(j%system_sets == 2):
        Tag = 'BBa4b2_v030'
        from JohnsonOnePopModels import OnePopImperfectFrag as Engine
        Run = partial(Engine, a = 4, b = 2, PCoal = 1 - pfrag, PFrag = pfrag)
    elif(j%system_sets == 3):
        Tag = 'Partitioning_v030'
        from JohnsonOnePopModels import OnePopImperfectFrag2 as Engine
        Run = partial(Engine, PCoal = 1 - pfrag, PFrag = pfrag)
    elif(j%system_sets == 4):
        Tag = 'Splitn1d2_v030' # Numerator and Denominator
        from JohnsonOnePopModels import OnePopImperfectFrag3 as Engine
        Run = partial(Engine, PCoal = 1 - pfrag, PFrag = pfrag)
    elif(j%system_sets == 5):
        Tag = 'Fragn1d2_v030' # Numerator and Denominator
        from JohnsonOnePopModels import OnePopImperfectFrag4 as Engine
        Run = partial(Engine, PCoal = 1 - pfrag, PFrag = pfrag)
    elif(j%system_sets == 6):
        Tag = 'Fragn1d10_v030' # Numerator and Denominator
        from JohnsonOnePopModels import OnePopImperfectFrag5 as Engine
        Run = partial(Engine, PCoal = 1 - pfrag, PFrag = pfrag)
    elif(j%system_sets == 7):
        Tag = 'Fragn9d10_v030' # Numerator and Denominator
        from JohnsonOnePopModels import OnePopImperfectFrag6 as Engine
        Run = partial(Engine, PCoal = 1 - pfrag, PFrag = pfrag)
    elif(j%system_sets == 8):
        Tag = 'PL110_v030' # Note output PL is 1 + alpha
        from JohnsonOnePopModels import OnePopImperfectFrag7 as Engine
        Run = partial(Engine, alpha = 0.10, PCoal = 1 - pfrag, PFrag = pfrag)
    elif(j%system_sets == 9):
        Tag = 'PL120_v030' # Note output PL is 1 + alpha
        from JohnsonOnePopModels import OnePopImperfectFrag7 as Engine
        Run = partial(Engine, alpha = 0.20, PCoal = 1 - pfrag, PFrag = pfrag)
    elif(j%system_sets == 10):
        Tag = 'PL130_v030' # Note output PL is 1 + alpha
        from JohnsonOnePopModels import OnePopImperfectFrag7 as Engine
        Run = partial(Engine, alpha = 0.30, PCoal = 1 - pfrag, PFrag = pfrag)
    elif(j%system_sets == 11):
        Tag = 'PL140_v030' # Note output PL is 1 + alpha
        from JohnsonOnePopModels import OnePopImperfectFrag7 as Engine
        Run = partial(Engine, alpha = 0.40, PCoal = 1 - pfrag, PFrag = pfrag)
    elif(j%system_sets == 12):
        Tag = 'PL160_v030' # Note output PL is 1 + alpha
        from JohnsonOnePopModels import OnePopImperfectFrag7 as Engine
        Run = partial(Engine, alpha = 0.60, PCoal = 1 - pfrag, PFrag = pfrag)
    elif(j%system_sets == 13):
        Tag = 'PL170_v030' # Note output PL is 1 + alpha
        from JohnsonOnePopModels import OnePopImperfectFrag7 as Engine
        Run = partial(Engine, alpha = 0.70, PCoal = 1 - pfrag, PFrag = pfrag)
    elif(j%system_sets == 14):
        Tag = 'PL180_v030' # Note output PL is 1 + alpha
        from JohnsonOnePopModels import OnePopImperfectFrag7 as Engine
        Run = partial(Engine, alpha = 0.80, PCoal = 1 - pfrag, PFrag = pfrag)
    elif(j%system_sets == 15):
        Tag = 'PL190_v030' # Note output PL is 1 + alpha
        from JohnsonOnePopModels import OnePopImperfectFrag7 as Engine
        Run = partial(Engine, alpha = 0.90, PCoal = 1 - pfrag, PFrag = pfrag)
    elif(j%system_sets == 16):
        Tag = 'ImpComb_v030' 
        from JohnsonOnePopModels import OnePopImperfectComb as Engine
        Run = partial(Engine, PCoal = 1 - pfrag, PFrag = pfrag)
    elif(j%system_sets == 17):
        Tag = 'ImpAll_v030' 
        from JohnsonOnePopModels import OnePopImperfect as Engine
        Run = partial(Engine, a = 1, b = 1, PCoal = 1 - pfrag, PFrag = pfrag)
        
    #%% Initialise Directory Structure
    DirectoryName = getcwd() + '/Data' + date.today().strftime('%y-%m-%d') + '/' + Label + Tag + '/'
    if not isdir(DirectoryName):
        makedirs(DirectoryName)

    #%% Proceed to Runs
    for i in Pops:
        #%% "Burn-in"
        Start = StartPosition(i)
        #Get Start significantly closer to the equilibrium "state" before recording.
        Start = Run(Population = Start,
                          EndTime = 1000000,
                          History = 0,
                          HistorySteps = 1)
        FileName = 'RandomPartition_' + Label + Tag + '_' +str(int(i)) + '_Run_' + str(int(j)) + '.mat'

        #%% Record seed and begin main run.
        savemat(DirectoryName + FileName, 
                mdict = {'RandomSeed' : theSeed}, 
                do_compression = True)
            
        for s in range(Slices):          
            History = Run(Population = Start,
                              EndTime = (EndTime)//Slices,
                              History = 1,
                              HistorySteps = Sampling)

            #Following: https://stackoverflow.com/a/29042033
            #append .mat file
            with open(DirectoryName + FileName, 'ab') as f:
                savemat(f, 
                        mdict = {Label + Tag + '_' + str(int(s)) : np.array(History)}, 
                        do_compression = True)
            
            Start = History[-1]
