# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 13:45:23 2016

@author: btf500
"""
from scipy.sparse import issparse as IsSparse
from numpy.random import uniform as Pick

def PickProportionalToSize(NumberOfGroupsOfSize, PopSize = -1):
    '''Picks an agent from NUMBEROFGROUPSOFSIZE
    and Returns that agent's group's size.
    '''
    if(IsSparse(NumberOfGroupsOfSize)):
        return SparsePickProportionalToSize(NumberOfGroupsOfSize, PopSize)
    
    lengthOfNumberOfGroupsOfSize = len(NumberOfGroupsOfSize)
    
    #In Case the population's size is not actually its length
    if(PopSize == -1):
        PopSize = lengthOfNumberOfGroupsOfSize
    
    AgentNumber = int(Pick(0,PopSize))
    
    AgentsSoFar = 0
    for Size in range(0, lengthOfNumberOfGroupsOfSize):
        AgentsSoFar += NumberOfGroupsOfSize[Size]*(Size+1)
        if(AgentNumber < AgentsSoFar):
            return Size + 1
            
def PickProportionalToSize_opt(NumberOfGroupsOfSize, PopSize = -1):
    '''Picks an agent from NUMBEROFGROUPSOFSIZE
    and Returns that agent's group's size.
    '''
    lengthOfNumberOfGroupsOfSize = len(NumberOfGroupsOfSize)
    
    #In Case the population's size is not actually its length
    if(PopSize == -1):
        PopSize = lengthOfNumberOfGroupsOfSize
    
    AgentNumber = int(Pick(0, PopSize))
    
    AgentsSoFar = 0
    for Size, NumberOfGroups in enumerate(NumberOfGroupsOfSize):
        AgentsSoFar += NumberOfGroups*(Size+1)
        if(AgentNumber < AgentsSoFar):
            return Size + 1
            
def PickUniformlyAtRandom(NumberOfGroupsOfSize):
    '''Picks a group from NUMBEROFGROUPSOFSIZE
    and returns that group's size.
    '''
    if(IsSparse(NumberOfGroupsOfSize)):
        return SparsePickUniformlyAtRandom(NumberOfGroupsOfSize)
    
    NumberOfGroups = sum(NumberOfGroupsOfSize)
    GroupNumber = int(Pick(0, NumberOfGroups))
    
    GroupsSoFar = 0
    for Size in range(0,len(NumberOfGroupsOfSize)):
        GroupsSoFar += NumberOfGroupsOfSize[Size]
        if(GroupNumber < GroupsSoFar):
            return Size + 1
            
def PickFromNGroupsPropToSize(ListOfPopulations, ListOfPopSizes):
    '''Picks an agent from LISTOFPOPULATIONS,
    which have matching LISTOFPOPSIZES, and
    Returns that agent's group's Population
    and Size as a list
    '''
    AgentNumber = int(Pick(0, sum(ListOfPopSizes)))
    
    AgentsSoFar = 0
    
    #We traverse by Population
    for PopulationNumber in range(0, len(ListOfPopulations)):
        for Size in range(0, len(ListOfPopulations[PopulationNumber])):
            AgentsSoFar += ListOfPopulations[PopulationNumber][Size]*(Size+1)
            if(AgentNumber < AgentsSoFar):
                return([PopulationNumber, Size + 1])
    
    return None
    
def PickFromNGroupsUnifAtRand(ListOfPopulations):
    '''Picks a group from LISTOFPOPULATIONS and
    Returns that group's Population and Size
    as a list.
    '''
    NumberOfGroups = sum([sum(Pop) for Pop in ListOfPopulations])
    GroupNumber = int(Pick(0, NumberOfGroups))
    
    GroupsSoFar = 0
    for PopulationNumber in range(0, len(ListOfPopulations)):
        for Size in range(0, len(ListOfPopulations[PopulationNumber])):
            GroupsSoFar += ListOfPopulations[PopulationNumber][Size]
            if(GroupNumber < GroupsSoFar):
                return ([PopulationNumber, Size + 1])

    return None
    
def PickNTimesProportionalToSize(NumberOfGroupsOfSize, 
                                 N, PopulationSize = -1,
                                 Unique = 1):
    '''Picks an agent from NUMBEROFGROUPSOFSIZE
    N times. If N < 1 (and N > 0) then N * 
    POPULATIONSIZE agents are picked.
    Groups picked will be unique if UNIQUE =/= 0
    and Returns those agent's group's sizes.
    WARNING: DOES NOT WORK SPARSELY.
    '''
    from counting import isum
    
    #Retrieve total Population Size.
    PopSize = PopulationSize
    if(PopulationSize == -1):
        PopSize = isum(NumberOfGroupsOfSize)
    
    #Create a Copy of the NumberOfGroupsOfSize
    CopyPop = list(NumberOfGroupsOfSize)
    
    #Create a list to return:
    if(N >= 1):
        RetVal = [-1] * N
        FullRange = range(N)
    else:
        RetVal = [-1] * int(N * PopSize)
        FullRange = range(int(N * PopSize))
        
    for i in FullRange:
        AgentNumber = int(Pick(0, PopSize))
    
        AgentsSoFar = 0
        for Size in range(0, len(CopyPop)):
            AgentsSoFar += CopyPop[Size] * (Size + 1)
            if(AgentNumber < AgentsSoFar):
                RetVal[i] = Size + 1
                break
            
        if(Unique):
            #Update Population so as not to choose same group
            CopyPop[RetVal[i] - 1] += -1
            PopSize += -RetVal[i]
    return RetVal
            
def PickNTimesUniformlyAtRandom(NumberOfGroupsOfSize, N):
    '''Picks an agent from NUMBEROFGROUPSOFSIZE
    N times. Groups picked will be unique
    and Returns those agent's group's sizes.
    WARNING: DOES NOT WORK SPARSELY.
    '''
    NumberOfGroups = sum(NumberOfGroupsOfSize)
    
    #Create a Copy of the NumberOfGroupsOfSize
    CopyPop = list(NumberOfGroupsOfSize)
    
    #Create a list to return:
    RetVal = [-1] * N

    for i in range(N):
        GroupNumber = int(Pick(0, NumberOfGroups))
    
        GroupsSoFar = 0
        for Size in range(0,len(CopyPop)):
            GroupsSoFar += CopyPop[Size]
            if(GroupNumber < GroupsSoFar):
                RetVal[i] = Size + 1
                break
        #Update Population so as not to choose same group
        CopyPop[RetVal[i]-1] += -1
        NumberOfGroups += -1
    return RetVal
                
                
    
############################################SPARSE:
def SparsePickProportionalToSize(NumberOfGroupsOfSize, PopSize = -1):
    '''Picks an agent from NUMBEROFGROUPSOFSIZE
    and Returns that agent's group's size.
    WARNING: Treats Input as a ROW
    '''
    #Retrieve elements list
    from scipy.sparse import find
    Elements = find(NumberOfGroupsOfSize)
    
    #Pick an Individual, uniformly at random:
    from numpy.random import uniform as PickIndividual
    from counting import isum
    
    #In Case the population's size is not actually its length
    if(PopSize == -1):
        PopSize = isum(NumberOfGroupsOfSize)
    
    AgentNumber = int(PickIndividual(0,PopSize))
    
    
    AgentsSoFar = 0
    for ElementNumber in range(len(Elements[0])):
        #Elements = [ListOfRows, ListOfColumns, ListOfValues]
        #           Should be 0,    Size - 1,   Number of Groups for a Size
        AgentsSoFar += (Elements[1][ElementNumber]+1) * Elements[2][ElementNumber]
        if(AgentNumber < AgentsSoFar):
            return Elements[1][ElementNumber]+1
            
def SparsePickUniformlyAtRandom(NumberOfGroupsOfSize):
    '''Picks a group from NUMBEROFGROUPSOFSIZE
    and returns that group's size.
    WARNING: Treats Input as a ROW
    '''
    #Retrieve elements list
    from scipy.sparse import find
    Elements = find(NumberOfGroupsOfSize)
        
    from numpy.random import uniform as PickGroup
    
    NumberOfGroups = sum(Elements[2])
    GroupNumber = int(PickGroup(0, NumberOfGroups))
    
    GroupsSoFar = 0
    for ElementNumber in range(len(Elements[0])):
        #Elements = [ListOfRows, ListOfColumns, ListOfValues]
        #           Should be 0,    Size - 1,   Number of Groups for a Size
        GroupsSoFar += Elements[2][ElementNumber]
        if(GroupNumber < GroupsSoFar):
            return Elements[1][ElementNumber]+1
    
#def SparsePickFromNGroupsPropToSize(ListOfPopulations, ListOfPopSizes):
    
#def SparsePickFromNGroupsUnifAtRand(ListOfPopulations):
    