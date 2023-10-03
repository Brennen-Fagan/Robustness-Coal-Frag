# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 10:25:20 2016

@author: btf500
"""

#Model According to Johnson et al. (2009, 2016):
#Common ecology quantifies human insurgency and
#New online ecology of adversarial aggregates: ISIS and beyond.

#(image of process may be found in New's Supplemental Materials, p. 15-18, 30)
#Process description from Common's Supplemental Materials:
#   ''On each timestep in the model, a[n agent] is selected... 
#   With probability [vfrag], the [agent's] group fragments into single agents.
#   Otherwise, with probability [vcoal], a second [agent from another group] is
#   selected and the two [agents'] groups coalesce...''

#Design note: 
#   We provide two versions of this model. The primary difference
#   is whether the coalescence step can end in failure due to lack of another
#   group. It appears from the code written by Johnson et al. in 2016 that 
#   they thought that this failure should not happen, but they neglect to 
#   mention this as a part of their model from what I can tell.

#Comparison with CoalFragDiscrete01 series:
#   The biggest difference here is whether we have one action per timestep,
#   or whether all groups can interact per timestep. CoalFragDiscrete01 takes
#   the latter approach. Here, we take the former. Note that the differential
#   equations appear to take the latter approach, which seems ignored by
#   Johnson et al. as a major difference between their discrete and continuous
#   models.

#Input assumptions from Johnson et al. (2016):
#   The initial conditions of this numerical simulation are such 
#   that all aggregates have size 1.
#   vfrag = .01, vcoal = .99, N = 10000, endtime from code = 10000000


#UNIVERSAL UTILITY FUNCTIONS AND DEFAULTS:
#from ReadDefaultOnePop import OnePopDefault as DefPop
#DefaultPopulation = DefPop()
from counting     import random_partition_Probabilistic_Divide_Conquer_1 as part
DefaultPopulation = part(10**4)

from numpy.random import uniform    as ActionValue
from numpy        import array # Used for adding numbers to lists.
from counting     import MinSum     as MinimumCount
from counting     import isum
from counting     import IndexAbove
from math         import ceil
from itertools    import zip_longest

def JohnsonOnePopMustSucceed(Population = DefaultPopulation, EndTime = 10000, 
                             PCoal = 0.99, PFrag = -1, verbose = 0,
                             History = 0, HistorySteps = 1):
    '''Take a POPULATION and, until the ENDTIME, Combine 
    (with probability PCOAL) and Fragment (with 
    probability PFRAG) it. Once the final timestep has
    occurred, we Return the final population. Note that 
    combining can only be chosen if there are two groups.
    Additionally, will Return the population at each timestep
    if HISTORY is not 0.
    '''
    
    #Some Python fixing:
    #If PFrag not assigned, assign it as the complement of PCoal
    if(verbose):
        print('Performing Fixing of Inputs to match expected style.')
    if(PFrag<0):
        PFrag = 1 - PCoal
    
    #Fix The received population so it is some form of list to work on.
    PopAsList = []
    if(isinstance(Population, int)):
        PopAsList = [Population]
    else:
        PopAsList = list(Population)
        
    #Transform the Population as a list into 
    #a list of the number of groups by size.
    #Use isum to make this transition easy.
     			    #Lazy Scaled Summation
    NumberOfGroupsOfSize = [0] * isum(PopAsList)
    for size in range(0, len(PopAsList)):
        NumberOfGroupsOfSize[size] = PopAsList[size]

    #Import Procedure Related Functions:
    from JohnsonAggregator import SingleAggregate as Combine
    from JohnsonDisaggregator import SingleDisaggregate as Fragment

    #import pdb; pdb.set_trace()
    
    #Procedure:Record History if applicable
    if(History):
        returnHistory = list()
    #On each timestep
    if(verbose):
        print('Beginning Simulation, Time = 0')
    NowTime = 0;
    if(History):
        HistorySinceLastSave = 0
        returnHistory.append(NumberOfGroupsOfSize)
    for NowTime in range(EndTime):
        AtLeast2Clusters = MinimumCount(NumberOfGroupsOfSize, 2)
    
        #Call appropriate action
        #Fragmentation if no other choice or Probability < PFrag
        #Coalesce otherwise
        if(not(AtLeast2Clusters) or ActionValue() < PFrag):
            if(verbose):
                print('Performing Fragment')
            NumberOfGroupsOfSize = Fragment(NumberOfGroupsOfSize)[0]
            
        else:
            if(verbose):
                print('Performing Combine')
            NumberOfGroupsOfSize = Combine(NumberOfGroupsOfSize)[0]
            
        #Record History if applicable
        if(History):
            if(HistorySinceLastSave + 1 == HistorySteps):
                returnHistory.append(NumberOfGroupsOfSize)
                HistorySinceLastSave = 0
            else:
                HistorySinceLastSave += 1
        if(verbose):
            print('Ending Timestep:', NowTime)

    if(History):
        return returnHistory
    return NumberOfGroupsOfSize
    
def PagettOnePop(Population = DefaultPopulation, EndTime = 10000, 
                            PCoal = 0.99, PFrag = -1, verbose = 0,
                            History = 0, HistorySteps = 1):
    '''Take a POPULATION and, until the ENDTIME, Combine 
    (with probability PCOAL) and Fragment (with 
    probability PFRAG) it. Once the final timestep has
    occurred, we Return the final population. Note that 
    combining can only be chosen if there are two groups.
    Additionally, will Return the population at each timestep
    if HISTORY is not 0. Picks groups without size bias.
    '''
    
    #Some Python fixing:
    #If PFrag not assigned, assign it as the complement of PCoal
    if(verbose):
        print('Performing Fixing of Inputs to match expected style.')
    if(PFrag<0):
        PFrag = 1 - PCoal
    
    #Fix The received population so it is some form of list to work on.
    PopAsList = []
    if(isinstance(Population, int)):
        PopAsList = [Population]
    else:
        PopAsList = list(Population)
        
    #Transform the Population as a list into 
    #a list of the number of groups by size.
    #Use isum to make this transition easy.
     			    #Lazy Scaled Summation
    NumberOfGroupsOfSize = [0] * isum(PopAsList)
    for size in range(0, len(PopAsList)):
        NumberOfGroupsOfSize[size] = PopAsList[size]

    #Import Procedure Related Functions:
    from JohnsonAggregator import PagettAggregate as Combine
    from JohnsonDisaggregator import PagettDisaggregate as Fragment

    #import pdb; pdb.set_trace()
    
    #Procedure:Record History if applicable
    if(History):
        returnHistory = list()
    #On each timestep
    if(verbose):
        print('Beginning Simulation, Time = 0')
    NowTime = 0;
    if(History):
        HistorySinceLastSave = 0
        returnHistory.append(NumberOfGroupsOfSize)
    for NowTime in range(EndTime):
        AtLeast2Clusters = MinimumCount(NumberOfGroupsOfSize, 2)
    
        #Call appropriate action
        #Fragmentation if no other choice or Probability < PFrag
        #Coalesce otherwise
        if(not(AtLeast2Clusters) or ActionValue() < PFrag):
            if(verbose):
                print('Performing Fragment')
            NumberOfGroupsOfSize = Fragment(NumberOfGroupsOfSize)[0]
            
        else:
            if(verbose):
                print('Performing Combine')
            NumberOfGroupsOfSize = Combine(NumberOfGroupsOfSize)[0]
            
        #Record History if applicable
        if(History):
            if(HistorySinceLastSave + 1 == HistorySteps):
                returnHistory.append(NumberOfGroupsOfSize)
                HistorySinceLastSave = 0
            else:
                HistorySinceLastSave += 1
        if(verbose):
            print('Ending Timestep:', NowTime)

    if(History):
        return returnHistory
    return NumberOfGroupsOfSize
    
def JohnsonOnePopMaySucceed(Population = DefaultPopulation, EndTime = 10000, 
                            PCoal = 0.99, PFrag = -1, verbose = 0,
                            History = 0, HistorySteps = 1):
    '''Take a POPULATION and, until the ENDTIME, Combine 
    (with probability PCOAL) and Fragment (with 
    probability PFRAG) it. Once the final timestep has
    occurred, we Return the final population. Note that 
    if there are not at least two groups, combining
    performs no action.
    Additionally, will Return the population at each timestep
    if HISTORY is not 0.
    '''
    
    #Some Python fixing:
    #If PFrag not assigned, assign it as the complement of PCoal
    if(verbose):
        print('Performing Fixing of Inputs to match expected style.')
    if(PFrag<0):
        PFrag = 1 - PCoal
    
    #Fix The received population so it is some form of list to work on.
    PopAsList = []
    if(isinstance(Population, int)):
        PopAsList = [Population]
    else:
        PopAsList = list(Population)
        
    #Transform the Population as a list into 
    #a list of the number of groups by size.
    #Use isum to make this transition easy.
     			    #Lazy Scaled Summation
    NumberOfGroupsOfSize = [0] * isum(PopAsList)
    for size in range(0, len(PopAsList)):
        NumberOfGroupsOfSize[size] = PopAsList[size]

    #Import Procedure Related Functions:
    from JohnsonAggregator import SingleAggregate as Combine
    from JohnsonDisaggregator import SingleDisaggregate as Fragment

    #Procedure: Record History if applicable
    if(History):
        returnHistory = list()
    #On each timestep
    if(verbose):
        print('Beginning Simulation, Time = 0')
    NowTime = 0;
    if(History):
        HistorySinceLastSave = 0
        returnHistory.append(NumberOfGroupsOfSize)
    for NowTime in range(EndTime):
        AtLeast2Clusters = MinimumCount(NumberOfGroupsOfSize, 2)
    
        #Call appropriate action
        #Fragmentation if Probability < PFrag
        #Coalesce otherwise
        if(ActionValue() < PFrag):
            if(verbose):
                print('Performing Fragment')
            NumberOfGroupsOfSize = Fragment(NumberOfGroupsOfSize)[0]
            
        else:
            if(AtLeast2Clusters):
                if(verbose):
                    print('Performing Combine')
                NumberOfGroupsOfSize = Combine(NumberOfGroupsOfSize)[0]
            else:
                if(verbose):
                    print('Combine failed, insufficient number of groups')
            
        #Record History if applicable
        if(History):
            if(HistorySinceLastSave + 1 == HistorySteps):
                returnHistory.append(NumberOfGroupsOfSize)
                HistorySinceLastSave = 0
            else:
                HistorySinceLastSave += 1
        if(verbose):
            print('Ending Timestep:', NowTime)

    if(History):
        return returnHistory        
    return NumberOfGroupsOfSize
    
def OnePopWithAttrition(Population = DefaultPopulation, EndTime = 10000, 
                             PCoal = 0.99, PFrag = -1, Attrition = 1,
                             Unique = 0, History = 0, HistorySteps = 1,
                             verbose = 0, debug = 0):
    '''Take a POPULATION and, until the ENDTIME, Combine 
    (with probability PCOAL) and Fragment (with probability
    PFRAG) it. Every timestep, up to ATTRITION unique
    groups lose an individual. Once the final timestep has
    occurred, we Return the final population. Note that 
    combining can only be chosen if there are two groups.
    Additionally, will Return the population at each timestep
    if HISTORY is not 0. 
    '''
    
    #Some Python fixing:
    #If PFrag not assigned, assign it as the complement of PCoal
    if(verbose):
        print('Performing Fixing of Inputs to match expected style.')
    if(PFrag<0):
        PFrag = 1 - PCoal
    
    #Fix The received population so it is some form of list to work on.
    PopAsList = []
    if(isinstance(Population, int)):
        PopAsList = [Population]
    else:
        PopAsList = list(Population)
        
    #Transform the Population as a list into 
    #a list of the number of groups by size.
    #Use isum to make this transition easy.
     			    #Lazy Scaled Summation
    PopulationSize = isum(PopAsList)
    NumberOfGroupsOfSize = [0] * PopulationSize
    for size in range(0, len(PopAsList)):
        NumberOfGroupsOfSize[size] = PopAsList[size]

    #Import Procedure Related Functions:
    from JohnsonAggregator import SingleAggregate as Combine
    from JohnsonDisaggregator import SingleDisaggregate as Fragment
    from JohnsonDisaggregator import IndividualsDisaggregate as Attrit
    
    #Procedure:Record History if applicable
    if(History):
        returnHistory = list()
    #On each timestep
    if(verbose):
        print('Beginning Simulation, Time = 0')
    NowTime = 0;
    if(History):
        HistorySinceLastSave = 0
        returnHistory.append(NumberOfGroupsOfSize)
        
    if(debug):
        import pdb; pdb.set_trace()

    for NowTime in range(EndTime):
        AtLeast2Clusters = MinimumCount(NumberOfGroupsOfSize, 2)
    
        #Call appropriate action
        #Fragmentation if no other choice or Probability < PFrag
        #Coalesce otherwise
        if(not(AtLeast2Clusters) or ActionValue() < PFrag):
            if(verbose):
                print('Performing Fragment')
            NumberOfGroupsOfSize = Fragment(NumberOfGroupsOfSize)[0]
            
        else:
            if(verbose):
                print('Performing Combine')
            NumberOfGroupsOfSize = Combine(NumberOfGroupsOfSize)[0]
            
        if(verbose):
            print('Performing Attrition')
        NumberOfGroupsOfSize = Attrit(NumberOfGroupsOfSize, Attrition, 
                                      PopulationSize, Unique)[0]
        
        #Record History if applicable
        if(History):
            if(HistorySinceLastSave + 1 == HistorySteps):
                returnHistory.append(NumberOfGroupsOfSize)
                HistorySinceLastSave = 0
            else:
                HistorySinceLastSave += 1
        if(verbose):
            print('Ending Timestep:', NowTime)

    if(History):
        return returnHistory
    return NumberOfGroupsOfSize   
    
def OnePopWithAccretion(Population = DefaultPopulation, EndTime = 10000, 
                             PCoal = 0.99, PFrag = -1, Accretion = 1,
                             Unique = 0, History = 0, HistorySteps = 1,
                             verbose = 0, debug = 0):
    '''Take a POPULATION and, until the ENDTIME, Combine 
    (with probability PCOAL) and Fragment (with probability
    PFRAG) it. Every timestep, up to ACCRETION unique
    groups gain an individual. Once the final timestep has
    occurred, we Return the final population. Note that 
    combining can only be chosen if there are two groups.
    Additionally, will Return the population at each timestep
    if HISTORY is not 0. 
    '''
    
    #Some Python fixing:
    #If PFrag not assigned, assign it as the complement of PCoal
    if(verbose):
        print('Performing Fixing of Inputs to match expected style.')
    if(PFrag<0):
        PFrag = 1 - PCoal
    
    #Fix The received population so it is some form of list to work on.
    PopAsList = []
    if(isinstance(Population, int)):
        PopAsList = [Population]
    else:
        PopAsList = list(Population)
        
    #Transform the Population as a list into 
    #a list of the number of groups by size.
    #Use isum to make this transition easy.
     			    #Lazy Scaled Summation
    PopulationSize = isum(PopAsList)
    NumberOfGroupsOfSize = [0] * PopulationSize
    for size in range(0, len(PopAsList)):
        NumberOfGroupsOfSize[size] = PopAsList[size]

    #Import Procedure Related Functions:
    from JohnsonAggregator import SingleAggregate as Combine
    from JohnsonDisaggregator import SingleDisaggregate as Fragment
    from JohnsonAggregator import IndividualsAggregate as Accrete
    
    #Procedure:Record History if applicable
    if(History):
        returnHistory = list()
    #On each timestep
    if(verbose):
        print('Beginning Simulation, Time = 0')
    NowTime = 0;
    if(History):
        HistorySinceLastSave = 0
        returnHistory.append(NumberOfGroupsOfSize)
        
    if(debug):
        import pdb; pdb.set_trace()

    for NowTime in range(EndTime):
        AtLeast2Clusters = MinimumCount(NumberOfGroupsOfSize, 2)
    
        #Call appropriate action
        #Fragmentation if no other choice or Probability < PFrag
        #Coalesce otherwise
        if(not(AtLeast2Clusters) or ActionValue() < PFrag):
            if(verbose):
                print('Performing Fragment')
            NumberOfGroupsOfSize = Fragment(NumberOfGroupsOfSize)[0]
            
        else:
            if(verbose):
                print('Performing Combine')
            NumberOfGroupsOfSize = Combine(NumberOfGroupsOfSize)[0]
            
        if(verbose):
            print('Performing Accretion')
        NumberOfGroupsOfSize = Accrete(NumberOfGroupsOfSize, Accretion, 
                                       PopulationSize, Unique)[0]
        
        #Record History if applicable
        if(History):
            if(HistorySinceLastSave + 1 == HistorySteps):
                returnHistory.append(NumberOfGroupsOfSize)
                HistorySinceLastSave = 0
            else:
                HistorySinceLastSave += 1
        if(verbose):
            print('Ending Timestep:', NowTime)

    if(History):
        return returnHistory
    return NumberOfGroupsOfSize
    
def OnePopWithAccrAttr(Population = DefaultPopulation, EndTime = 10000, 
                             PCoal = 0.99, PFrag = -1, Accretion = 1,
                             Attrition = 1,
                             Unique = 0, History = 0, HistorySteps = 1,
                             verbose = 0, debug = 0):
    '''Take a POPULATION and, until the ENDTIME, Combine 
    (with probability PCOAL) and Fragment (with probability
    PFRAG) it. Every timestep, up to ACCRETION unique
    groups gain an individual and up to ATTRITION unique
    groups lose an individual. Once the final timestep has
    occurred, we Return the final population. Note that 
    combining can only be chosen if there are two groups.
    Additionally, will Return the population at each timestep
    if HISTORY is not 0. 
    '''
    
    #Some Python fixing:
    #If PFrag not assigned, assign it as the complement of PCoal
    if(verbose):
        print('Performing Fixing of Inputs to match expected style.')
    if(PFrag<0):
        PFrag = 1 - PCoal
    
    #Fix The received population so it is some form of list to work on.
    PopAsList = []
    if(isinstance(Population, int)):
        PopAsList = [Population]
    else:
        PopAsList = list(Population)
        
    #Transform the Population as a list into 
    #a list of the number of groups by size.
    #Use isum to make this transition easy.
     			    #Lazy Scaled Summation
    PopulationSize = isum(PopAsList)
    NumberOfGroupsOfSize = [0] * PopulationSize
    for size in range(0, len(PopAsList)):
        NumberOfGroupsOfSize[size] = PopAsList[size]

    #Import Procedure Related Functions:
    from JohnsonAggregator import SingleAggregate as Combine
    from JohnsonAggregator import IndividualsAggregate as Accrete
    from JohnsonDisaggregator import SingleDisaggregate as Fragment
    from JohnsonDisaggregator import IndividualsDisaggregate as Attrit
    
    #Procedure:Record History if applicable
    if(History):
        returnHistory = list()
    #On each timestep
    if(verbose):
        print('Beginning Simulation, Time = 0')
    NowTime = 0;
    if(History):
        HistorySinceLastSave = 0
        returnHistory.append(NumberOfGroupsOfSize)
        
    if(debug):
        import pdb; pdb.set_trace()

        
    for NowTime in range(EndTime):
        AtLeast2Clusters = MinimumCount(NumberOfGroupsOfSize, 2)
    
        #Call appropriate action
        #Fragmentation if no other choice or Probability < PFrag
        #Coalesce otherwise
        if(not(AtLeast2Clusters) or ActionValue() < PFrag):
            if(verbose):
                print('Performing Fragment')
            NumberOfGroupsOfSize = Fragment(NumberOfGroupsOfSize)[0]
            
        else:
            if(verbose):
                print('Performing Combine')
            NumberOfGroupsOfSize = Combine(NumberOfGroupsOfSize)[0]
            
        if(verbose):
            print('Performing Accretion')
        NumberOfGroupsOfSize = Accrete(NumberOfGroupsOfSize, Accretion, 
                                       PopulationSize, Unique)[0]
        
        if(verbose):
            print('Performing Attrition')
        NumberOfGroupsOfSize = Attrit(NumberOfGroupsOfSize, Attrition, 
                                      PopulationSize, Unique)[0]
        
        #Record History if applicable
        if(History):
            if(HistorySinceLastSave + 1 == HistorySteps):
                returnHistory.append(NumberOfGroupsOfSize)
                HistorySinceLastSave = 0
            else:
                HistorySinceLastSave += 1
        if(verbose):
            print('Ending Timestep:', NowTime)

    if(History):
        return returnHistory
    return NumberOfGroupsOfSize
 
def OnePopWithAttrAccr(Population = DefaultPopulation, EndTime = 10000, 
                             PCoal = 0.99, PFrag = -1, Accretion = 1,
                             Attrition = 1,
                             Unique = 0, History = 0, HistorySteps = 1,
                             verbose = 0, debug = 0):
    '''Take a POPULATION and, until the ENDTIME, Combine 
    (with probability PCOAL) and Fragment (with probability
    PFRAG) it. Every timestep, up to ACCRETION unique
    groups gain an individual and up to ATTRITION unique
    groups lose an individual. Once the final timestep has
    occurred, we Return the final population. Note that 
    combining can only be chosen if there are two groups.
    Additionally, will Return the population at each timestep
    if HISTORY is not 0. 
    '''
    
    #Some Python fixing:
    #If PFrag not assigned, assign it as the complement of PCoal
    if(verbose):
        print('Performing Fixing of Inputs to match expected style.')
    if(PFrag<0):
        PFrag = 1 - PCoal
    
    #Fix The received population so it is some form of list to work on.
    PopAsList = []
    if(isinstance(Population, int)):
        PopAsList = [Population]
    else:
        PopAsList = list(Population)
        
    #Transform the Population as a list into 
    #a list of the number of groups by size.
    #Use isum to make this transition easy.
     			    #Lazy Scaled Summation
    PopulationSize = isum(PopAsList)
    NumberOfGroupsOfSize = [0] * PopulationSize
    for size in range(0, len(PopAsList)):
        NumberOfGroupsOfSize[size] = PopAsList[size]

    #Import Procedure Related Functions:
    from JohnsonAggregator import SingleAggregate as Combine
    from JohnsonAggregator import IndividualsAggregate as Accrete
    from JohnsonDisaggregator import SingleDisaggregate as Fragment
    from JohnsonDisaggregator import IndividualsDisaggregate as Attrit
    
    #Procedure:Record History if applicable
    if(History):
        returnHistory = list()
    #On each timestep
    if(verbose):
        print('Beginning Simulation, Time = 0')
    NowTime = 0;
    if(History):
        HistorySinceLastSave = 0
        returnHistory.append(NumberOfGroupsOfSize)
        
    if(debug):
        import pdb; pdb.set_trace()

        
    for NowTime in range(EndTime):
        AtLeast2Clusters = MinimumCount(NumberOfGroupsOfSize, 2)
    
        #Call appropriate action
        #Fragmentation if no other choice or Probability < PFrag
        #Coalesce otherwise
        if(not(AtLeast2Clusters) or ActionValue() < PFrag):
            if(verbose):
                print('Performing Fragment')
            NumberOfGroupsOfSize = Fragment(NumberOfGroupsOfSize)[0]
            
        else:
            if(verbose):
                print('Performing Combine')
            NumberOfGroupsOfSize = Combine(NumberOfGroupsOfSize)[0]
        
        if(verbose):
            print('Performing Attrition')
        NumberOfGroupsOfSize = Attrit(NumberOfGroupsOfSize, Attrition, 
                                      PopulationSize, Unique)[0]
            
        if(verbose):
            print('Performing Accretion')
        NumberOfGroupsOfSize = Accrete(NumberOfGroupsOfSize, Accretion, 
                                       PopulationSize, Unique)[0]
        
        #Record History if applicable
        if(History):
            if(HistorySinceLastSave + 1 == HistorySteps):
                returnHistory.append(NumberOfGroupsOfSize)
                HistorySinceLastSave = 0
            else:
                HistorySinceLastSave += 1
        if(verbose):
            print('Ending Timestep:', NowTime)

    if(History):
        return returnHistory
    return NumberOfGroupsOfSize    

def OnePopImperfect(Population = DefaultPopulation, EndTime = 10000, 
                             PCoal = 0.99, PFrag = -1, History = 0, 
                             a = 1, b = 1, HistorySteps = 1,
                             verbose = 0, debug = 0):
    '''Take a POPULATION and, until the ENDTIME, Combine 
    (with probability PCOAL) and Fragment (with probability
    PFRAG) it. Combination combines a group and part of
    another group, while fragmentation breaks groups into
    smaller groups. Once the final timestep has
    occurred, we Return the final population. Note that 
    combining can only be chosen if there are two groups.
    Additionally, will Return the population at each timestep
    if HISTORY is not 0. 
    '''
    
    #Some Python fixing:
    #If PFrag not assigned, assign it as the complement of PCoal
    if(verbose):
        print('Performing Fixing of Inputs to match expected style.')
    if(PFrag<0):
        PFrag = 1 - PCoal
    
    #Fix The received population so it is some form of list to work on.
    PopAsList = []
    if(isinstance(Population, int)):
        PopAsList = [Population]
    else:
        PopAsList = list(Population)
        
    #Transform the Population as a list into 
    #a list of the number of groups by size.
    #Use isum to make this transition easy.
     			    #Lazy Scaled Summation
    NumberOfGroupsOfSize = [0] * isum(PopAsList)
    for size in range(0, len(PopAsList)):
        NumberOfGroupsOfSize[size] = PopAsList[size]

    #Import Procedure Related Functions:
    from JohnsonAggregator import ImperfectAggregate as Combine
    from JohnsonDisaggregator import ImperfectFragmentation1 as Fragment
    
    #Procedure:Record History if applicable
    if(History):
        returnHistory = list()
    #On each timestep
    if(verbose):
        print('Beginning Simulation, Time = 0')
    NowTime = 0;
    if(History):
        HistorySinceLastSave = 0
        returnHistory.append(NumberOfGroupsOfSize)
        
    if(debug):
        import pdb; pdb.set_trace()

        
    for NowTime in range(EndTime):
        AtLeast2Clusters = MinimumCount(NumberOfGroupsOfSize, 2)
    
        #Call appropriate action
        #Fragmentation if no other choice or Probability < PFrag
        #Coalesce otherwise
        if(not(AtLeast2Clusters) or ActionValue() < PFrag):
            if(verbose):
                print('Performing Imperfect Fragment')
            NumberOfGroupsOfSize = Fragment(NumberOfGroupsOfSize, a, b)[0]
            
        else:
            if(verbose):
                print('Performing Imperfect Combine')
            NumberOfGroupsOfSize = Combine(NumberOfGroupsOfSize)[0]
            
        
        #Record History if applicable
        if(History):
            if(HistorySinceLastSave + 1 == HistorySteps):
                returnHistory.append(NumberOfGroupsOfSize)
                HistorySinceLastSave = 0
            else:
                HistorySinceLastSave += 1
        if(verbose):
            print('Ending Timestep:', NowTime)

    if(History):
        return returnHistory
    return NumberOfGroupsOfSize
    
def OnePopImperfect2(Population = DefaultPopulation, EndTime = 10000, 
                             PCoal = 0.99, PFrag = -1, History = 0,
                             HistorySteps = 1,
                             verbose = 0, debug = 0):
    '''Take a POPULATION and, until the ENDTIME, Combine 
    (with probability PCOAL) and Fragment (with probability
    PFRAG) it. Combination combines a group and part of
    another group, while fragmentation breaks groups into
    smaller groups. Once the final timestep has
    occurred, we Return the final population. Note that 
    combining can only be chosen if there are two groups.
    Additionally, will Return the population at each timestep
    if HISTORY is not 0. 
    '''
    
    #Some Python fixing:
    #If PFrag not assigned, assign it as the complement of PCoal
    if(verbose):
        print('Performing Fixing of Inputs to match expected style.')
    if(PFrag<0):
        PFrag = 1 - PCoal
    
    #Fix The received population so it is some form of list to work on.
    PopAsList = []
    if(isinstance(Population, int)):
        PopAsList = [Population]
    else:
        PopAsList = list(Population)
        
    #Transform the Population as a list into 
    #a list of the number of groups by size.
    #Use isum to make this transition easy.
     			    #Lazy Scaled Summation
    NumberOfGroupsOfSize = [0] * isum(PopAsList)
    for size in range(0, len(PopAsList)):
        NumberOfGroupsOfSize[size] = PopAsList[size]

    #Import Procedure Related Functions:
    from JohnsonAggregator import ImperfectAggregate as Combine
    from JohnsonDisaggregator import ImperfectFragmentation2 as Fragment
    
    #Procedure:Record History if applicable
    if(History):
        returnHistory = list()
    #On each timestep
    if(verbose):
        print('Beginning Simulation, Time = 0')
    NowTime = 0;
    if(History):
        HistorySinceLastSave = 0
        returnHistory.append(NumberOfGroupsOfSize)
        
    if(debug):
        import pdb; pdb.set_trace()

        
    for NowTime in range(EndTime):
        AtLeast2Clusters = MinimumCount(NumberOfGroupsOfSize, 2)
    
        #Call appropriate action
        #Fragmentation if no other choice or Probability < PFrag
        #Coalesce otherwise
        if(not(AtLeast2Clusters) or ActionValue() < PFrag):
            if(verbose):
                print('Performing Imperfect Fragment')
            NumberOfGroupsOfSize = Fragment(NumberOfGroupsOfSize)[0]
            
        else:
            if(verbose):
                print('Performing Imperfect Combine')
            NumberOfGroupsOfSize = Combine(NumberOfGroupsOfSize)[0]
            
        
        #Record History if applicable
        if(History):
            if(HistorySinceLastSave + 1 == HistorySteps):
                returnHistory.append(NumberOfGroupsOfSize)
                HistorySinceLastSave = 0
            else:
                HistorySinceLastSave += 1
        if(verbose):
            print('Ending Timestep:', NowTime)

    if(History):
        return returnHistory
    return NumberOfGroupsOfSize
    
def OnePopImperfect3(Population = DefaultPopulation, EndTime = 10000, 
                             PCoal = 0.99, PFrag = -1, History = 0,
                             HistorySteps = 1,
                             verbose = 0, debug = 0):
    '''Take a POPULATION and, until the ENDTIME, Combine 
    (with probability PCOAL) and Fragment (with probability
    PFRAG) it. Combination combines a group and part of
    another group, while fragmentation breaks groups into
    smaller groups. Once the final timestep has
    occurred, we Return the final population. Note that 
    combining can only be chosen if there are two groups.
    Additionally, will Return the population at each timestep
    if HISTORY is not 0. 
    '''
    
    #Some Python fixing:
    #If PFrag not assigned, assign it as the complement of PCoal
    if(verbose):
        print('Performing Fixing of Inputs to match expected style.')
    if(PFrag<0):
        PFrag = 1 - PCoal
    
    #Fix The received population so it is some form of list to work on.
    PopAsList = []
    if(isinstance(Population, int)):
        PopAsList = [Population]
    else:
        PopAsList = list(Population)
        
    #Transform the Population as a list into 
    #a list of the number of groups by size.
    #Use isum to make this transition easy.
     			    #Lazy Scaled Summation
    NumberOfGroupsOfSize = [0] * isum(PopAsList)
    for size in range(0, len(PopAsList)):
        NumberOfGroupsOfSize[size] = PopAsList[size]

    #Import Procedure Related Functions:
    from JohnsonAggregator import ImperfectAggregate as Combine
    from JohnsonDisaggregator import ImperfectFragmentation3 as Fragment
    
    #Procedure:Record History if applicable
    if(History):
        returnHistory = list()
    #On each timestep
    if(verbose):
        print('Beginning Simulation, Time = 0')
    NowTime = 0;
    if(History):
        HistorySinceLastSave = 0
        returnHistory.append(NumberOfGroupsOfSize)
        
    if(debug):
        import pdb; pdb.set_trace()

        
    for NowTime in range(EndTime):
        AtLeast2Clusters = MinimumCount(NumberOfGroupsOfSize, 2)
    
        #Call appropriate action
        #Fragmentation if no other choice or Probability < PFrag
        #Coalesce otherwise
        if(not(AtLeast2Clusters) or ActionValue() < PFrag):
            if(verbose):
                print('Performing Imperfect Fragment')
            NumberOfGroupsOfSize = Fragment(NumberOfGroupsOfSize)[0]
            
        else:
            if(verbose):
                print('Performing Imperfect Combine')
            NumberOfGroupsOfSize = Combine(NumberOfGroupsOfSize)[0]
            
        
        #Record History if applicable
        if(History):
            if(HistorySinceLastSave + 1 == HistorySteps):
                returnHistory.append(NumberOfGroupsOfSize)
                HistorySinceLastSave = 0
            else:
                HistorySinceLastSave += 1
        if(verbose):
            print('Ending Timestep:', NowTime)

    if(History):
        return returnHistory
    return NumberOfGroupsOfSize
    
def OnePopImperfect4(Population = DefaultPopulation, EndTime = 10000, 
                             PCoal = 0.99, PFrag = -1, History = 0,
                             HistorySteps = 1,
                             verbose = 0, debug = 0):
    '''Take a POPULATION and, until the ENDTIME, Combine 
    (with probability PCOAL) and Fragment (with probability
    PFRAG) it. Combination combines a group and part of
    another group, while fragmentation breaks groups into
    smaller groups. Once the final timestep has
    occurred, we Return the final population. Note that 
    combining can only be chosen if there are two groups.
    Additionally, will Return the population at each timestep
    if HISTORY is not 0. 
    '''
    
    #Some Python fixing:
    #If PFrag not assigned, assign it as the complement of PCoal
    if(verbose):
        print('Performing Fixing of Inputs to match expected style.')
    if(PFrag<0):
        PFrag = 1 - PCoal
    
    #Fix The received population so it is some form of list to work on.
    PopAsList = []
    if(isinstance(Population, int)):
        PopAsList = [Population]
    else:
        PopAsList = list(Population)
        
    #Transform the Population as a list into 
    #a list of the number of groups by size.
    #Use isum to make this transition easy.
     			    #Lazy Scaled Summation
    NumberOfGroupsOfSize = [0] * isum(PopAsList)
    for size in range(0, len(PopAsList)):
        NumberOfGroupsOfSize[size] = PopAsList[size]

    #Import Procedure Related Functions:
    from JohnsonAggregator import ImperfectAggregate as Combine
    from JohnsonDisaggregator import ImperfectFragmentation4 as Fragment
    
    #Procedure:Record History if applicable
    if(History):
        returnHistory = list()
    #On each timestep
    if(verbose):
        print('Beginning Simulation, Time = 0')
    NowTime = 0;
    if(History):
        HistorySinceLastSave = 0
        returnHistory.append(NumberOfGroupsOfSize)
        
    if(debug):
        import pdb; pdb.set_trace()

        
    for NowTime in range(EndTime):
        AtLeast2Clusters = MinimumCount(NumberOfGroupsOfSize, 2)
    
        #Call appropriate action
        #Fragmentation if no other choice or Probability < PFrag
        #Coalesce otherwise
        if(not(AtLeast2Clusters) or ActionValue() < PFrag):
            if(verbose):
                print('Performing Imperfect Fragment')
            NumberOfGroupsOfSize = Fragment(NumberOfGroupsOfSize)[0]
            
        else:
            if(verbose):
                print('Performing Imperfect Combine')
            NumberOfGroupsOfSize = Combine(NumberOfGroupsOfSize)[0]
            
        
        #Record History if applicable
        if(History):
            if(HistorySinceLastSave + 1 == HistorySteps):
                returnHistory.append(NumberOfGroupsOfSize)
                HistorySinceLastSave = 0
            else:
                HistorySinceLastSave += 1
        if(verbose):
            print('Ending Timestep:', NowTime)

    if(History):
        return returnHistory
    return NumberOfGroupsOfSize
    
def OnePopImperfect5(Population = DefaultPopulation, EndTime = 10000, 
                             PCoal = 0.99, PFrag = -1, History = 0,
                             HistorySteps = 1,
                             verbose = 0, debug = 0):
    '''Take a POPULATION and, until the ENDTIME, Combine 
    (with probability PCOAL) and Fragment (with probability
    PFRAG) it. Combination combines a group and part of
    another group, while fragmentation breaks groups into
    smaller groups. Once the final timestep has
    occurred, we Return the final population. Note that 
    combining can only be chosen if there are two groups.
    Additionally, will Return the population at each timestep
    if HISTORY is not 0. 
    '''
    
    #Some Python fixing:
    #If PFrag not assigned, assign it as the complement of PCoal
    if(verbose):
        print('Performing Fixing of Inputs to match expected style.')
    if(PFrag<0):
        PFrag = 1 - PCoal
    
    #Fix The received population so it is some form of list to work on.
    PopAsList = []
    if(isinstance(Population, int)):
        PopAsList = [Population]
    else:
        PopAsList = list(Population)
        
    #Transform the Population as a list into 
    #a list of the number of groups by size.
    #Use isum to make this transition easy.
     			    #Lazy Scaled Summation
    NumberOfGroupsOfSize = [0] * isum(PopAsList)
    for size in range(0, len(PopAsList)):
        NumberOfGroupsOfSize[size] = PopAsList[size]

    #Import Procedure Related Functions:
    from JohnsonAggregator import ImperfectAggregate as Combine
    from JohnsonDisaggregator import ImperfectFragmentation5 as Fragment
    
    #Procedure:Record History if applicable
    if(History):
        returnHistory = list()
    #On each timestep
    if(verbose):
        print('Beginning Simulation, Time = 0')
    NowTime = 0;
    if(History):
        HistorySinceLastSave = 0
        returnHistory.append(NumberOfGroupsOfSize)
        
    if(debug):
        import pdb; pdb.set_trace()

        
    for NowTime in range(EndTime):
        AtLeast2Clusters = MinimumCount(NumberOfGroupsOfSize, 2)
    
        #Call appropriate action
        #Fragmentation if no other choice or Probability < PFrag
        #Coalesce otherwise
        if(not(AtLeast2Clusters) or ActionValue() < PFrag):
            if(verbose):
                print('Performing Imperfect Fragment')
            NumberOfGroupsOfSize = Fragment(NumberOfGroupsOfSize)[0]
            
        else:
            if(verbose):
                print('Performing Imperfect Combine')
            NumberOfGroupsOfSize = Combine(NumberOfGroupsOfSize)[0]
            
        
        #Record History if applicable
        if(History):
            if(HistorySinceLastSave + 1 == HistorySteps):
                returnHistory.append(NumberOfGroupsOfSize)
                HistorySinceLastSave = 0
            else:
                HistorySinceLastSave += 1
        if(verbose):
            print('Ending Timestep:', NowTime)

    if(History):
        return returnHistory
    return NumberOfGroupsOfSize
    
def OnePopImperfect6(Population = DefaultPopulation, EndTime = 10000, 
                             PCoal = 0.99, PFrag = -1, History = 0,
                             HistorySteps = 1,
                             verbose = 0, debug = 0):
    '''Take a POPULATION and, until the ENDTIME, Combine 
    (with probability PCOAL) and Fragment (with probability
    PFRAG) it. Combination combines a group and part of
    another group, while fragmentation breaks groups into
    smaller groups. Once the final timestep has
    occurred, we Return the final population. Note that 
    combining can only be chosen if there are two groups.
    Additionally, will Return the population at each timestep
    if HISTORY is not 0. 
    '''
    
    #Some Python fixing:
    #If PFrag not assigned, assign it as the complement of PCoal
    if(verbose):
        print('Performing Fixing of Inputs to match expected style.')
    if(PFrag<0):
        PFrag = 1 - PCoal
    
    #Fix The received population so it is some form of list to work on.
    PopAsList = []
    if(isinstance(Population, int)):
        PopAsList = [Population]
    else:
        PopAsList = list(Population)
        
    #Transform the Population as a list into 
    #a list of the number of groups by size.
    #Use isum to make this transition easy.
     			    #Lazy Scaled Summation
    NumberOfGroupsOfSize = [0] * isum(PopAsList)
    for size in range(0, len(PopAsList)):
        NumberOfGroupsOfSize[size] = PopAsList[size]

    #Import Procedure Related Functions:
    from JohnsonAggregator import ImperfectAggregate as Combine
    from JohnsonDisaggregator import ImperfectFragmentation6 as Fragment
    
    #Procedure:Record History if applicable
    if(History):
        returnHistory = list()
    #On each timestep
    if(verbose):
        print('Beginning Simulation, Time = 0')
    NowTime = 0;
    if(History):
        HistorySinceLastSave = 0
        returnHistory.append(NumberOfGroupsOfSize)
        
    if(debug):
        import pdb; pdb.set_trace()

        
    for NowTime in range(EndTime):
        AtLeast2Clusters = MinimumCount(NumberOfGroupsOfSize, 2)
    
        #Call appropriate action
        #Fragmentation if no other choice or Probability < PFrag
        #Coalesce otherwise
        if(not(AtLeast2Clusters) or ActionValue() < PFrag):
            if(verbose):
                print('Performing Imperfect Fragment')
            NumberOfGroupsOfSize = Fragment(NumberOfGroupsOfSize)[0]
            
        else:
            if(verbose):
                print('Performing Imperfect Combine')
            NumberOfGroupsOfSize = Combine(NumberOfGroupsOfSize)[0]
            
        
        #Record History if applicable
        if(History):
            if(HistorySinceLastSave + 1 == HistorySteps):
                returnHistory.append(NumberOfGroupsOfSize)
                HistorySinceLastSave = 0
            else:
                HistorySinceLastSave += 1
        if(verbose):
            print('Ending Timestep:', NowTime)

    if(History):
        return returnHistory
    return NumberOfGroupsOfSize
    
def OnePopImperfectFrag(Population = DefaultPopulation, EndTime = 10000, 
                             PCoal = 0.99, PFrag = -1, History = 0, 
                             a = 1, b = 1, HistorySteps = 1,
                             verbose = 0, debug = 0):
    '''Take a POPULATION and, until the ENDTIME, Combine 
    (with probability PCOAL) and Fragment (with probability
    PFRAG) it. Fragmentation breaks groups into
    smaller groups. Once the final timestep has
    occurred, we Return the final population. Note that 
    combining can only be chosen if there are two groups.
    Additionally, will Return the population at each timestep
    if HISTORY is not 0. 
    '''
    
    #Some Python fixing:
    #If PFrag not assigned, assign it as the complement of PCoal
    if(verbose):
        print('Performing Fixing of Inputs to match expected style.')
    if(PFrag<0):
        PFrag = 1 - PCoal
    
    #Fix The received population so it is some form of list to work on.
    PopAsList = []
    if(isinstance(Population, int)):
        PopAsList = [Population]
    else:
        PopAsList = list(Population)
        
    #Transform the Population as a list into 
    #a list of the number of groups by size.
    #Use isum to make this transition easy.
     			    #Lazy Scaled Summation
    NumberOfGroupsOfSize = [0] * isum(PopAsList)
    for size in range(0, len(PopAsList)):
        NumberOfGroupsOfSize[size] = PopAsList[size]

    #Import Procedure Related Functions:
    from JohnsonAggregator import SingleAggregate as Combine
    from JohnsonDisaggregator import ImperfectFragmentation1 as Fragment
    
    #Procedure:Record History if applicable
    if(History):
        returnHistory = list()
    #On each timestep
    if(verbose):
        print('Beginning Simulation, Time = 0')
    NowTime = 0;
    if(History):
        HistorySinceLastSave = 0
        returnHistory.append(NumberOfGroupsOfSize)
        
    if(debug):
        import pdb; pdb.set_trace()

    for NowTime in range(EndTime):
        AtLeast2Clusters = MinimumCount(NumberOfGroupsOfSize, 2)
    
        #Call appropriate action
        #Fragmentation if no other choice or Probability < PFrag
        #Coalesce otherwise
        if(not(AtLeast2Clusters) or ActionValue() < PFrag):
            if(verbose):
                print('Performing Imperfect Fragment')
            NumberOfGroupsOfSize = Fragment(NumberOfGroupsOfSize, a, b)[0]
            
        else:
            if(verbose):
                print('Performing Combine')
            NumberOfGroupsOfSize = Combine(NumberOfGroupsOfSize)[0]
            
        #Record History if applicable
        if(History):
            if(HistorySinceLastSave + 1 == HistorySteps):
                returnHistory.append(NumberOfGroupsOfSize)
                HistorySinceLastSave = 0
            else:
                HistorySinceLastSave += 1
        if(verbose):
            print('Ending Timestep:', NowTime)

    if(History):
        return returnHistory
    return NumberOfGroupsOfSize
    
def OnePopImperfectFrag2(Population = DefaultPopulation, EndTime = 10000, 
                             PCoal = 0.99, PFrag = -1, History = 0,
                             HistorySteps = 1,
                             verbose = 0, debug = 0):
    '''Take a POPULATION and, until the ENDTIME, Combine 
    (with probability PCOAL) and Fragment (with probability
    PFRAG) it. Fragmentation breaks groups into
    smaller groups. Once the final timestep has
    occurred, we Return the final population. Note that 
    combining can only be chosen if there are two groups.
    Additionally, will Return the population at each timestep
    if HISTORY is not 0. 
    '''
    
    #Some Python fixing:
    #If PFrag not assigned, assign it as the complement of PCoal
    if(verbose):
        print('Performing Fixing of Inputs to match expected style.')
    if(PFrag<0):
        PFrag = 1 - PCoal
    
    #Fix The received population so it is some form of list to work on.
    PopAsList = []
    if(isinstance(Population, int)):
        PopAsList = [Population]
    else:
        PopAsList = list(Population)
        
    #Transform the Population as a list into 
    #a list of the number of groups by size.
    #Use isum to make this transition easy.
     			    #Lazy Scaled Summation
    NumberOfGroupsOfSize = [0] * isum(PopAsList)
    for size in range(0, len(PopAsList)):
        NumberOfGroupsOfSize[size] = PopAsList[size]

    #Import Procedure Related Functions:
    from JohnsonAggregator import SingleAggregate as Combine
    from JohnsonDisaggregator import ImperfectFragmentation2 as Fragment
    
    #Procedure:Record History if applicable
    if(History):
        returnHistory = list()
    #On each timestep
    if(verbose):
        print('Beginning Simulation, Time = 0')
    NowTime = 0;
    if(History):
        HistorySinceLastSave = 0
        returnHistory.append(NumberOfGroupsOfSize)
        
    if(debug):
        import pdb; pdb.set_trace()

    for NowTime in range(EndTime):
        AtLeast2Clusters = MinimumCount(NumberOfGroupsOfSize, 2)
    
        #Call appropriate action
        #Fragmentation if no other choice or Probability < PFrag
        #Coalesce otherwise
        if(not(AtLeast2Clusters) or ActionValue() < PFrag):
            if(verbose):
                print('Performing Imperfect Fragment')
            NumberOfGroupsOfSize = Fragment(NumberOfGroupsOfSize)[0]
            
        else:
            if(verbose):
                print('Performing Combine')
            NumberOfGroupsOfSize = Combine(NumberOfGroupsOfSize)[0]
            
        #Record History if applicable
        if(History):
            if(HistorySinceLastSave + 1 == HistorySteps):
                returnHistory.append(NumberOfGroupsOfSize)
                HistorySinceLastSave = 0
            else:
                HistorySinceLastSave += 1
        if(verbose):
            print('Ending Timestep:', NowTime)

    if(History):
        return returnHistory
    return NumberOfGroupsOfSize
    
def OnePopImperfectFrag3(Population = DefaultPopulation, EndTime = 10000, 
                             PCoal = 0.99, PFrag = -1, History = 0,
                             HistorySteps = 1,
                             verbose = 0, debug = 0):
    '''Take a POPULATION and, until the ENDTIME, Combine 
    (with probability PCOAL) and Fragment (with probability
    PFRAG) it. Fragmentation breaks groups into
    smaller groups. Once the final timestep has
    occurred, we Return the final population. Note that 
    combining can only be chosen if there are two groups.
    Additionally, will Return the population at each timestep
    if HISTORY is not 0. 
    '''
    
    #Some Python fixing:
    #If PFrag not assigned, assign it as the complement of PCoal
    if(verbose):
        print('Performing Fixing of Inputs to match expected style.')
    if(PFrag<0):
        PFrag = 1 - PCoal
    
    #Fix The received population so it is some form of list to work on.
    PopAsList = []
    if(isinstance(Population, int)):
        PopAsList = [Population]
    else:
        PopAsList = list(Population)
        
    #Transform the Population as a list into 
    #a list of the number of groups by size.
    #Use isum to make this transition easy.
     			    #Lazy Scaled Summation
    NumberOfGroupsOfSize = [0] * isum(PopAsList)
    for size in range(0, len(PopAsList)):
        NumberOfGroupsOfSize[size] = PopAsList[size]

    #Import Procedure Related Functions:
    from JohnsonAggregator import SingleAggregate as Combine
    from JohnsonDisaggregator import ImperfectFragmentation3 as Fragment
    
    #Procedure:Record History if applicable
    if(History):
        returnHistory = list()
    #On each timestep
    if(verbose):
        print('Beginning Simulation, Time = 0')
    NowTime = 0;
    if(History):
        HistorySinceLastSave = 0
        returnHistory.append(NumberOfGroupsOfSize)
        
    if(debug):
        import pdb; pdb.set_trace()

    for NowTime in range(EndTime):
        AtLeast2Clusters = MinimumCount(NumberOfGroupsOfSize, 2)
    
        #Call appropriate action
        #Fragmentation if no other choice or Probability < PFrag
        #Coalesce otherwise
        if(not(AtLeast2Clusters) or ActionValue() < PFrag):
            if(verbose):
                print('Performing Imperfect Fragment')
            NumberOfGroupsOfSize = Fragment(NumberOfGroupsOfSize)[0]
            
        else:
            if(verbose):
                print('Performing Combine')
            NumberOfGroupsOfSize = Combine(NumberOfGroupsOfSize)[0]
            
        #Record History if applicable
        if(History):
            if(HistorySinceLastSave + 1 == HistorySteps):
                returnHistory.append(NumberOfGroupsOfSize)
                HistorySinceLastSave = 0
            else:
                HistorySinceLastSave += 1
        if(verbose):
            print('Ending Timestep:', NowTime)

    if(History):
        return returnHistory
    return NumberOfGroupsOfSize
    
def OnePopImperfectFrag4(Population = DefaultPopulation, EndTime = 10000, 
                             PCoal = 0.99, PFrag = -1, History = 0,
                             HistorySteps = 1,
                             verbose = 0, debug = 0):
    '''Take a POPULATION and, until the ENDTIME, Combine 
    (with probability PCOAL) and Fragment (with probability
    PFRAG) it. Fragmentation breaks groups into
    smaller groups. Once the final timestep has
    occurred, we Return the final population. Note that 
    combining can only be chosen if there are two groups.
    Additionally, will Return the population at each timestep
    if HISTORY is not 0. 
    '''
    
    #Some Python fixing:
    #If PFrag not assigned, assign it as the complement of PCoal
    if(verbose):
        print('Performing Fixing of Inputs to match expected style.')
    if(PFrag<0):
        PFrag = 1 - PCoal
    
    #Fix The received population so it is some form of list to work on.
    PopAsList = []
    if(isinstance(Population, int)):
        PopAsList = [Population]
    else:
        PopAsList = list(Population)
        
    #Transform the Population as a list into 
    #a list of the number of groups by size.
    #Use isum to make this transition easy.
     			    #Lazy Scaled Summation
    NumberOfGroupsOfSize = [0] * isum(PopAsList)
    for size in range(0, len(PopAsList)):
        NumberOfGroupsOfSize[size] = PopAsList[size]

    #Import Procedure Related Functions:
    from JohnsonAggregator import SingleAggregate as Combine
    from JohnsonDisaggregator import ImperfectFragmentation4 as Fragment
    
    #Procedure:Record History if applicable
    if(History):
        returnHistory = list()
    #On each timestep
    if(verbose):
        print('Beginning Simulation, Time = 0')
    NowTime = 0;
    if(History):
        HistorySinceLastSave = 0
        returnHistory.append(NumberOfGroupsOfSize)
        
    if(debug):
        import pdb; pdb.set_trace()

    for NowTime in range(EndTime):
        AtLeast2Clusters = MinimumCount(NumberOfGroupsOfSize, 2)
    
        #Call appropriate action
        #Fragmentation if no other choice or Probability < PFrag
        #Coalesce otherwise
        if(not(AtLeast2Clusters) or ActionValue() < PFrag):
            if(verbose):
                print('Performing Imperfect Fragment')
            NumberOfGroupsOfSize = Fragment(NumberOfGroupsOfSize)[0]
            
        else:
            if(verbose):
                print('Performing Combine')
            NumberOfGroupsOfSize = Combine(NumberOfGroupsOfSize)[0]
            
        #Record History if applicable
        if(History):
            if(HistorySinceLastSave + 1 == HistorySteps):
                returnHistory.append(NumberOfGroupsOfSize)
                HistorySinceLastSave = 0
            else:
                HistorySinceLastSave += 1
        if(verbose):
            print('Ending Timestep:', NowTime)

    if(History):
        return returnHistory
    return NumberOfGroupsOfSize
    
def OnePopImperfectFrag5(Population = DefaultPopulation, EndTime = 10000, 
                             PCoal = 0.99, PFrag = -1, History = 0,
                             HistorySteps = 1,
                             verbose = 0, debug = 0):
    '''Take a POPULATION and, until the ENDTIME, Combine 
    (with probability PCOAL) and Fragment (with probability
    PFRAG) it. Fragmentation breaks groups into
    smaller groups. Once the final timestep has
    occurred, we Return the final population. Note that 
    combining can only be chosen if there are two groups.
    Additionally, will Return the population at each timestep
    if HISTORY is not 0. 
    '''
    
    #Some Python fixing:
    #If PFrag not assigned, assign it as the complement of PCoal
    if(verbose):
        print('Performing Fixing of Inputs to match expected style.')
    if(PFrag<0):
        PFrag = 1 - PCoal
    
    #Fix The received population so it is some form of list to work on.
    PopAsList = []
    if(isinstance(Population, int)):
        PopAsList = [Population]
    else:
        PopAsList = list(Population)
        
    #Transform the Population as a list into 
    #a list of the number of groups by size.
    #Use isum to make this transition easy.
     			    #Lazy Scaled Summation
    NumberOfGroupsOfSize = [0] * isum(PopAsList)
    for size in range(0, len(PopAsList)):
        NumberOfGroupsOfSize[size] = PopAsList[size]

    #Import Procedure Related Functions:
    from JohnsonAggregator import SingleAggregate as Combine
    from JohnsonDisaggregator import ImperfectFragmentation5 as Fragment
    
    #Procedure:Record History if applicable
    if(History):
        returnHistory = list()
    #On each timestep
    if(verbose):
        print('Beginning Simulation, Time = 0')
    NowTime = 0;
    if(History):
        HistorySinceLastSave = 0
        returnHistory.append(NumberOfGroupsOfSize)
        
    if(debug):
        import pdb; pdb.set_trace()

    for NowTime in range(EndTime):
        AtLeast2Clusters = MinimumCount(NumberOfGroupsOfSize, 2)
    
        #Call appropriate action
        #Fragmentation if no other choice or Probability < PFrag
        #Coalesce otherwise
        if(not(AtLeast2Clusters) or ActionValue() < PFrag):
            if(verbose):
                print('Performing Imperfect Fragment')
            NumberOfGroupsOfSize = Fragment(NumberOfGroupsOfSize)[0]
            
        else:
            if(verbose):
                print('Performing Combine')
            NumberOfGroupsOfSize = Combine(NumberOfGroupsOfSize)[0]
            
        #Record History if applicable
        if(History):
            if(HistorySinceLastSave + 1 == HistorySteps):
                returnHistory.append(NumberOfGroupsOfSize)
                HistorySinceLastSave = 0
            else:
                HistorySinceLastSave += 1
        if(verbose):
            print('Ending Timestep:', NowTime)

    if(History):
        return returnHistory
    return NumberOfGroupsOfSize
    
def OnePopImperfectFrag6(Population = DefaultPopulation, EndTime = 10000, 
                             PCoal = 0.99, PFrag = -1, History = 0,
                             HistorySteps = 1,
                             verbose = 0, debug = 0):
    '''Take a POPULATION and, until the ENDTIME, Combine 
    (with probability PCOAL) and Fragment (with probability
    PFRAG) it. Fragmentation breaks groups into
    smaller groups. Once the final timestep has
    occurred, we Return the final population. Note that 
    combining can only be chosen if there are two groups.
    Additionally, will Return the population at each timestep
    if HISTORY is not 0. 
    '''
    
    #Some Python fixing:
    #If PFrag not assigned, assign it as the complement of PCoal
    if(verbose):
        print('Performing Fixing of Inputs to match expected style.')
    if(PFrag<0):
        PFrag = 1 - PCoal
    
    #Fix The received population so it is some form of list to work on.
    PopAsList = []
    if(isinstance(Population, int)):
        PopAsList = [Population]
    else:
        PopAsList = list(Population)
        
    #Transform the Population as a list into 
    #a list of the number of groups by size.
    #Use isum to make this transition easy.
     			    #Lazy Scaled Summation
    NumberOfGroupsOfSize = [0] * isum(PopAsList)
    for size in range(0, len(PopAsList)):
        NumberOfGroupsOfSize[size] = PopAsList[size]

    #Import Procedure Related Functions:
    from JohnsonAggregator import SingleAggregate as Combine
    from JohnsonDisaggregator import ImperfectFragmentation6 as Fragment
    
    #Procedure:Record History if applicable
    if(History):
        returnHistory = list()
    #On each timestep
    if(verbose):
        print('Beginning Simulation, Time = 0')
    NowTime = 0;
    if(History):
        HistorySinceLastSave = 0
        returnHistory.append(NumberOfGroupsOfSize)
        
    if(debug):
        import pdb; pdb.set_trace()

    for NowTime in range(EndTime):
        AtLeast2Clusters = MinimumCount(NumberOfGroupsOfSize, 2)
    
        #Call appropriate action
        #Fragmentation if no other choice or Probability < PFrag
        #Coalesce otherwise
        if(not(AtLeast2Clusters) or ActionValue() < PFrag):
            if(verbose):
                print('Performing Imperfect Fragment')
            NumberOfGroupsOfSize = Fragment(NumberOfGroupsOfSize)[0]
            
        else:
            if(verbose):
                print('Performing Combine')
            NumberOfGroupsOfSize = Combine(NumberOfGroupsOfSize)[0]
            
        #Record History if applicable
        if(History):
            if(HistorySinceLastSave + 1 == HistorySteps):
                returnHistory.append(NumberOfGroupsOfSize)
                HistorySinceLastSave = 0
            else:
                HistorySinceLastSave += 1
        if(verbose):
            print('Ending Timestep:', NowTime)

    if(History):
        return returnHistory
    return NumberOfGroupsOfSize
    
def OnePopImperfectFrag7(Population = DefaultPopulation, EndTime = 10000, 
                             PCoal = 0.99, PFrag = -1, History = 0,
                             alpha = 0,
                             HistorySteps = 1,
                             verbose = 0, debug = 0):
    '''Take a POPULATION and, until the ENDTIME, Combine 
    (with probability PCOAL) and Fragment (with probability
    PFRAG) it. Fragmentation breaks groups into
    smaller groups. Once the final timestep has
    occurred, we Return the final population. Note that 
    combining can only be chosen if there are two groups.
    Additionally, will Return the population at each timestep
    if HISTORY is not 0. 
    '''
    
    #Some Python fixing:
    #If PFrag not assigned, assign it as the complement of PCoal
    if(verbose):
        print('Performing Fixing of Inputs to match expected style.')
    if(PFrag<0):
        PFrag = 1 - PCoal
    
    #Fix The received population so it is some form of list to work on.
    PopAsList = []
    if(isinstance(Population, int)):
        PopAsList = [Population]
    else:
        PopAsList = list(Population)
        
    #Transform the Population as a list into 
    #a list of the number of groups by size.
    #Use isum to make this transition easy.
     			    #Lazy Scaled Summation
    NumberOfGroupsOfSize = [0] * isum(PopAsList)
    for size in range(0, len(PopAsList)):
        NumberOfGroupsOfSize[size] = PopAsList[size]

    #Import Procedure Related Functions:
    from JohnsonAggregator import SingleAggregate as Combine
    from JohnsonDisaggregator import ImperfectFragmentation7 as Fragment
    
    #Procedure:Record History if applicable
    if(History):
        returnHistory = list()
    #On each timestep
    if(verbose):
        print('Beginning Simulation, Time = 0')
    NowTime = 0;
    if(History):
        HistorySinceLastSave = 0
        returnHistory.append(NumberOfGroupsOfSize)
        
    if(debug):
        import pdb; pdb.set_trace()

    for NowTime in range(EndTime):
        AtLeast2Clusters = MinimumCount(NumberOfGroupsOfSize, 2)
    
        #Call appropriate action
        #Fragmentation if no other choice or Probability < PFrag
        #Coalesce otherwise
        if(not(AtLeast2Clusters) or ActionValue() < PFrag):
            if(verbose):
                print('Performing Imperfect Fragment')
            NumberOfGroupsOfSize = Fragment(NumberOfGroupsOfSize, alpha)[0]
            
        else:
            if(verbose):
                print('Performing Combine')
            NumberOfGroupsOfSize = Combine(NumberOfGroupsOfSize)[0]
            
        #Record History if applicable
        if(History):
            if(HistorySinceLastSave + 1 == HistorySteps):
                returnHistory.append(NumberOfGroupsOfSize)
                HistorySinceLastSave = 0
            else:
                HistorySinceLastSave += 1
        if(verbose):
            print('Ending Timestep:', NowTime)

    if(History):
        return returnHistory
    return NumberOfGroupsOfSize
    
def OnePopImperfectFrag8(Population = DefaultPopulation, EndTime = 10000, 
                             PCoal = 0.99, PFrag = -1, History = 0,
                             alpha = 2,
                             HistorySteps = 1,
                             verbose = 0, debug = 0):
    '''Take a POPULATION and, until the ENDTIME, Combine 
    (with probability PCOAL) and Fragment (with probability
    PFRAG) it. Fragmentation breaks groups into
    smaller groups. Once the final timestep has
    occurred, we Return the final population. Note that 
    combining can only be chosen if there are two groups.
    Additionally, will Return the population at each timestep
    if HISTORY is not 0. 
    Alpha in (2, inf), as this version uses Price's Network
    simulations for determining fragmentation (the degree
    of his network simulations is powerlaw of exponent at
    least 2).
    '''
    
    #Some Python fixing:
    #If PFrag not assigned, assign it as the complement of PCoal
    if(verbose):
        print('Performing Fixing of Inputs to match expected style.')
    if(PFrag<0):
        PFrag = 1 - PCoal
    
    #Fix The received population so it is some form of list to work on.
    PopAsList = []
    if(isinstance(Population, int)):
        PopAsList = [Population]
    else:
        PopAsList = list(Population)
        
    #Transform the Population as a list into 
    #a list of the number of groups by size.
    #Use isum to make this transition easy.
     			    #Lazy Scaled Summation
    NumberOfGroupsOfSize = [0] * isum(PopAsList)
    for size in range(0, len(PopAsList)):
        NumberOfGroupsOfSize[size] = PopAsList[size]

    #Import Procedure Related Functions:
    from JohnsonAggregator import SingleAggregate as Combine
    from JohnsonDisaggregator import ImperfectFragmentation8 as Fragment
    
    #Procedure:Record History if applicable
    if(History):
        returnHistory = list()
    #On each timestep
    if(verbose):
        print('Beginning Simulation, Time = 0')
    NowTime = 0;
    if(History):
        HistorySinceLastSave = 0
        returnHistory.append(NumberOfGroupsOfSize)
        
    if(debug):
        import pdb; pdb.set_trace()

    for NowTime in range(EndTime):
        AtLeast2Clusters = MinimumCount(NumberOfGroupsOfSize, 2)
    
        #Call appropriate action
        #Fragmentation if no other choice or Probability < PFrag
        #Coalesce otherwise
        if(not(AtLeast2Clusters) or ActionValue() < PFrag):
            if(verbose):
                print('Performing Imperfect Fragment')
            NumberOfGroupsOfSize = Fragment(NumberOfGroupsOfSize, alpha)[0]
            
        else:
            if(verbose):
                print('Performing Combine')
            NumberOfGroupsOfSize = Combine(NumberOfGroupsOfSize)[0]
            
        #Record History if applicable
        if(History):
            if(HistorySinceLastSave + 1 == HistorySteps):
                returnHistory.append(NumberOfGroupsOfSize)
                HistorySinceLastSave = 0
            else:
                HistorySinceLastSave += 1
        if(verbose):
            print('Ending Timestep:', NowTime)

    if(History):
        return returnHistory
    return NumberOfGroupsOfSize
    
def OnePopImperfectComb(Population = DefaultPopulation, EndTime = 10000, 
                             PCoal = 0.99, PFrag = -1, History = 0,
                             HistorySteps = 1,
                             verbose = 0, debug = 0):
    '''Take a POPULATION and, until the ENDTIME, Combine 
    (with probability PCOAL) and Fragment (with probability
    PFRAG) it. Combination combines a group and part of
    another group. Once the final timestep has
    occurred, we Return the final population. Note that 
    combining can only be chosen if there are two groups.
    Additionally, will Return the population at each timestep
    if HISTORY is not 0. 
    '''
    
    #Some Python fixing:
    #If PFrag not assigned, assign it as the complement of PCoal
    if(verbose):
        print('Performing Fixing of Inputs to match expected style.')
    if(PFrag<0):
        PFrag = 1 - PCoal
    
    #Fix The received population so it is some form of list to work on.
    PopAsList = []
    if(isinstance(Population, int)):
        PopAsList = [Population]
    else:
        PopAsList = list(Population)
        
    #Transform the Population as a list into 
    #a list of the number of groups by size.
    #Use isum to make this transition easy.
     			    #Lazy Scaled Summation
    NumberOfGroupsOfSize = [0] * isum(PopAsList)
    for size in range(0, len(PopAsList)):
        NumberOfGroupsOfSize[size] = PopAsList[size]

    #Import Procedure Related Functions:
    from JohnsonAggregator import ImperfectAggregate as Combine
    from JohnsonDisaggregator import SingleDisaggregate as Fragment
    
    #Procedure:Record History if applicable
    if(History):
        returnHistory = list()
    #On each timestep
    if(verbose):
        print('Beginning Simulation, Time = 0')
    NowTime = 0;
    if(History):
        HistorySinceLastSave = 0
        returnHistory.append(NumberOfGroupsOfSize)
        
    if(debug):
        import pdb; pdb.set_trace()

        
    for NowTime in range(EndTime):
        AtLeast2Clusters = MinimumCount(NumberOfGroupsOfSize, 2)
    
        #Call appropriate action
        #Fragmentation if no other choice or Probability < PFrag
        #Coalesce otherwise
        if(not(AtLeast2Clusters) or ActionValue() < PFrag):
            if(verbose):
                print('Performing Fragment')
            NumberOfGroupsOfSize = Fragment(NumberOfGroupsOfSize)[0]
            
        else:
            if(verbose):
                print('Performing Imperfect Combine')
            NumberOfGroupsOfSize = Combine(NumberOfGroupsOfSize)[0]
            
        
        #Record History if applicable
        if(History):
            if(HistorySinceLastSave + 1 == HistorySteps):
                returnHistory.append(NumberOfGroupsOfSize)
                HistorySinceLastSave = 0
            else:
                HistorySinceLastSave += 1
        if(verbose):
            print('Ending Timestep:', NowTime)

    if(History):
        return returnHistory
    return NumberOfGroupsOfSize

def OnePopMultiStep(Population = DefaultPopulation, EndTime = 10000, 
                             PCoal = 0.99, PFrag = -1, Steps = 4,
                             History = 0, HistorySteps = 1,
                             verbose = 0, debug = 0):
    '''Take a POPULATION and, until the ENDTIME, Combine 
    (with probability PCOAL) and Fragment (with probability
    PFRAG) STEPS groups per timestep. Once the final timestep 
    has occurred, we Return the final population. Note that 
    combining can only be chosen if there are STEPS groups.
    Additionally, will Return the population at each timestep
    if HISTORY is not 0. 
    '''
    
    #Some Python fixing:
    #If PFrag not assigned, assign it as the complement of PCoal
    if(verbose):
        print('Performing Fixing of Inputs to match expected style.')
    if(PFrag<0):
        PFrag = 1 - PCoal
    
    #Fix The received population so it is some form of list to work on.
    PopAsList = []
    if(isinstance(Population, int)):
        PopAsList = [Population]
    else:
        PopAsList = list(Population)
        
    #Transform the Population as a list into 
    #a list of the number of groups by size.
    #Use isum to make this transition easy.
     			    #Lazy Scaled Summation
    PopulationSize = isum(PopAsList)
    NumberOfGroupsOfSize = [0] * PopulationSize
    for size in range(0, len(PopAsList)):
        NumberOfGroupsOfSize[size] = PopAsList[size]

    #Import Procedure Related Functions:
    from JohnsonAggregator import MultipleAggregate as Combine
    from JohnsonDisaggregator import MultipleDisaggregate as Fragment
    
    #Procedure:Record History if applicable
    if(History):
        returnHistory = list()
    #On each timestep
    if(verbose):
        print('Beginning Simulation, Time = 0')
    NowTime = 0;
    if(History):
        HistorySinceLastSave = 0
        returnHistory.append(NumberOfGroupsOfSize)
        
    if(debug):
        import pdb; pdb.set_trace()

        
    for NowTime in range(EndTime):
        AtLeastSTEPSClusters = MinimumCount(NumberOfGroupsOfSize, Steps)
    
        #Call appropriate action
        #Fragmentation if no other choice or Probability < PFrag
        #Coalesce otherwise
        if(not(AtLeastSTEPSClusters) or ActionValue() < PFrag):
            if(verbose):
                print('Performing Fragment ', Steps, ' times')
            NumberOfGroupsOfSize = Fragment(NumberOfGroupsOfSize, Steps)[0]
            
        else:
            if(verbose):
                print('Performing  Combine', Steps, ' times')
            NumberOfGroupsOfSize = Combine(NumberOfGroupsOfSize, Steps, PopulationSize)[0]
            
        
        #Record History if applicable
        if(History):
            if(HistorySinceLastSave + 1 == HistorySteps):
                returnHistory.append(NumberOfGroupsOfSize)
                HistorySinceLastSave = 0
            else:
                HistorySinceLastSave += 1
        if(verbose):
            print('Ending Timestep:', NowTime)

    if(History):
        return returnHistory
    return NumberOfGroupsOfSize
    

DefaultEndTime = 27*365
def ExamplePopFunc(Time, TotalPopulation, AveragePopulation = -1, YearPeriod = 365):
    '''An example population function which fluctuates
    a population periodically. Takes TIME as a non-constant
    input variable. It takes constants TOTALPOPULATION,
    AVERAGEPOPULATION, and YEARPERIOD to rescale a cosine
    function such that the resulting functions period is
    the YearPeriod, the average is the AveragePopulation,
    and the maximum is the TotalPopulation. It returns
    the current population as an integer.
    '''
    from numpy import cos, pi
    if(AveragePopulation == -1):
        AveragePopulation = 3/4 * TotalPopulation
    #Average Population at 0, Total Population at Pi/2
    #YearPeriod is scaling so that the period is
    #an appropriate value.
    Result = AveragePopulation + (TotalPopulation-AveragePopulation) * cos(Time * 2*pi/YearPeriod)
    return int(round(Result))
    
def DecreasingPopFunc(Time, TotalPopulation, Period = DefaultEndTime):
    '''Decreases the population until it reaches the
    value originally specified. Starts at 10 times that
    number, and decreases via the function
    f(t) = floor(9*TotalPopulation*exp(-t*ln(Period)/Period)+TotalPopulation).
    This function starts at 10 times and at approximately Period,
    it reaches 1 times.
    Moves quickly at first, then slows down.
    '''
    from numpy import exp
    from numpy import log as ln
    from math import floor
    
    return floor(9*TotalPopulation*exp(-1*Time*ln(Period)/Period)+TotalPopulation)
    
def IncreasingPopFunc(Time, TotalPopulation, Period = DefaultEndTime): 
    '''Increases the population until it reaches the value
    originally specified. Starts at 1 individual and
    increases via the function
    f(t) = floor(TotalPopulation*1.442*ln(t/(Period)+1)+1)
    Moves more quickly at first, but not quickly overall.
    '''
    from numpy import log as ln
    from math import floor
    
    return floor(TotalPopulation*1.442*ln(Time/Period+1)+1)
    
    
def UnstablePop(Population = DefaultPopulation, EndTime = DefaultEndTime, 
                             PCoal = 0.99, PFrag = -1, PopFunc = ExamplePopFunc,
                             History = 0, HistorySteps = 1,
                             verbose = 0, debug = 0):
    '''Take a POPULATION and, until the ENDTIME, Combine 
    (with probability PCOAL) and Fragment (with probability
    PFRAG) each timestep. Additionally varies the Current
    Population via the function POPFUNC, which takes the 
    current time and total population as inputs and 
    returns the current population as an integer. Once 
    the final timestep has occurred, we Return the final 
    population. Note that combining can only be chosen 
    if there are two groups.
    Additionally, will Return the population at each timestep
    if HISTORY is not 0. 
    '''
    
    #Some Python fixing:
    #If PFrag not assigned, assign it as the complement of PCoal
    if(verbose):
        print('Performing Fixing of Inputs to match expected style.')
    if(PFrag<0):
        PFrag = 1 - PCoal
    
    #Fix The received population so it is some form of list to work on.
    PopAsList = []
    if(isinstance(Population, int)):
        PopAsList = [Population]
    else:
        PopAsList = list(Population)
        
    #Transform the Population as a list into 
    #a list of the number of groups by size.
    #Use isum to make this transition easy.
     			    #Lazy Scaled Summation
    TotalPopulation = isum(PopAsList)
    NumberOfGroupsOfSize = [0] * TotalPopulation
    for size in range(0, len(PopAsList)):
        NumberOfGroupsOfSize[size] = PopAsList[size]

    #Import Procedure Related Functions:
    from JohnsonAggregator import SingleAggregate2 as Combine
    from JohnsonDisaggregator import SingleDisaggregate2 as Fragment
    from JohnsonDisaggregator import RemoveAnIndividualUAR as IndivLoss
    
    #Procedure:Record History if applicable
    if(History):
        returnHistory = list()
    #On each timestep
    if(verbose):
        print('Beginning Simulation, Time = 0')
    NowTime = 0;
    if(History):
        HistorySinceLastSave = 0
        returnHistory.append(NumberOfGroupsOfSize)
        
    if(debug):
        import pdb; pdb.set_trace()
        
    #Pre-allocate so as to not continuously have to reperform calculations.
    PreviousPopulationLevel = TotalPopulation
    CurrentPopulationLevels = [PopFunc(x, TotalPopulation) for x in range(EndTime)]
    for NowTime in range(EndTime):
        
        if(verbose):
            print('Adjusting Population Level: ', CurrentPopulationLevels[NowTime])
        
        if(PreviousPopulationLevel < CurrentPopulationLevels[NowTime]):
            if(verbose):
                print('Adding Individuals')
            NumberOfGroupsOfSize[0] += CurrentPopulationLevels[NowTime] - PreviousPopulationLevel
            PreviousPopulationLevel += CurrentPopulationLevels[NowTime] - PreviousPopulationLevel
        
        if(verbose):
            print('Removing Individuals from groups non-uniquely, proportional to size')
        while(PreviousPopulationLevel > CurrentPopulationLevels[NowTime]):
            #Move an agent from a group to the individuals
            NumberOfGroupsOfSize = IndivLoss(NumberOfGroupsOfSize, 1)
            #Remove an individual
            NumberOfGroupsOfSize[0] += -1
            #Note the current population level
            PreviousPopulationLevel += -1
        
        AtLeast2Clusters = MinimumCount(NumberOfGroupsOfSize, 2)
    
        #Call appropriate action
        #Fragmentation if no other choice or Probability < PFrag
        #Coalesce otherwise
        if(not(AtLeast2Clusters) or ActionValue() < PFrag):
            if(verbose):
                print('Performing Fragment')
            NumberOfGroupsOfSize = Fragment(NumberOfGroupsOfSize, PreviousPopulationLevel)[0]
            
        else:
            if(verbose):
                print('Performing Combine')
            NumberOfGroupsOfSize = Combine(NumberOfGroupsOfSize, PreviousPopulationLevel)[0]
            
        #Record History if applicable
        if(History):
            if(HistorySinceLastSave + 1 == HistorySteps):
                returnHistory.append(NumberOfGroupsOfSize)
                HistorySinceLastSave = 0
            else:
                HistorySinceLastSave += 1
            
        if(verbose):
            print('Ending Timestep:', NowTime)

    if(History):
        return returnHistory
    return NumberOfGroupsOfSize

def MixedKernelModel(
        Population = DefaultPopulation, EndTime = 10000, 
        PCoal = 0.99, PFrag = -1, 
        FragmentKernel = 'Constant',
        FragmentDistribution = 'Total',
        FragmentBarrier = 0,
        CoalescenceKernel = 'Multiplicative',
        verbose = 0, History = 0, HistorySteps = 1,
        Flux = True,
        debug = 0, **kwargs):
    '''
    Inputs:
        Population = DefaultPopulation, 
            A number designating the number of individuals, or a list where
            the ith entry corresponds to the number of groups of (i + 1)th 
            size. Note that the list need not be full, e.g. [0, 4] has a
            population of 8 organized in 4 pairs.
        EndTime = 10000, 
            The number of steps taken during the simulation.
        PCoal = 0.99, 
            The probability of coalescence/combination/coagulation/aggregation 
            on a given time step.
        PFrag = -1, 
            If -1, 1 - PCoal.
            Else, if PCoal + PFrag > 1, rescale until == 1.
            Else, 1 - PCoal - PFrag == P(Nothing) on a time step.
        FragmentKernel = 'Constant',
            Options: 
                'Constant' or 'Pagett' or '0' == Pick groups Unif.@Rand.
                'Multiplicative' or 'Mult' or '1' == Pick groups with
                probability proportional to their size.
        FragmentDistribution = 'Total',
            Options:
                'Total' or 'Shattering' or 'Perfect' or 'Atomize' or 'Atomise'
                Group to be fragmented converted to individuals.
                'CRP'
                'PN'
        FragmentBarrier = 0,
            Number within [0, 1]. This is a percentage of the population
            required for fragmentation to occur. If 0, then any size group can 
            fragment; if 1, only a group with the entire population in it can 
            fragment. Since we are interested in fragmentation as an exogeneous 
            process (e.g. a government rooting out terrorists) as opposed to an
            endogeneous process (e.g. infighting amongst groups of terrorists), 
            we apply the barrier before determining if fragmentation occurs, as
            opposed to selecting a group, and then determining if it can 
            fragment.
        CoalescenceKernel = 'Multiplicative',
            Options: 
                'Constant' or 'Pagett' or '0' == Pick groups Unif.@Rand.
                'Multiplicative' or 'Mult' or '1' == Pick groups with
                probability proportional to their size.
        verbose = 0, 
            If 1, comment on each step taken.
        History = 0, 
            If 1, output the configuration once every HistorySteps steps.
        HistorySteps = 1
            Controls the number of steps in between saved configurations, i.e.
            2 would mean saving every other step.
        **kwargs,
            Can contain Accretion, Attrition, Exponent, Exp, Alpha, or alpha,
            with the latter four considered synonyms. Make sure to choose
            an appropriate Distribution or Kernel above with an appropriate
            keyword. If multiple versions are provided, the last key is used.
    Outputs:
        if(History):
            Rows of the Number of Groups of Specific Sizes at a given time.
        if(!History):
            Rows of the Number of Groups of Specific Sizes at the end.
    '''
    
    #Some Python fixing:
    #If PFrag not assigned, assign it as the complement of PCoal
    if(verbose):
        print('Performing Fixing of Inputs to match expected style.')
    if(PFrag < 0):
        PFrag = 1 - PCoal
    if(PCoal + PFrag > 1):
        PTotal = PCoal + PFrag
        PCoal = PCoal / PTotal
        PFrag = PFrag / PTotal
    
    #Fix The received population so it is some form of list to work on.
    PopAsList = []
    if(isinstance(Population, int)):
        PopAsList = [Population]
    else:
        PopAsList = list(Population)
        
    #Transform the Population as a list into 
    #a list of the number of groups by size.
    #Use isum to make this transition easy.
                         #Lazy Scaled Summation
    NumberOfGroupsOfSize = [0] * isum(PopAsList)
    for size in range(0, len(PopAsList)):
        NumberOfGroupsOfSize[size] = PopAsList[size]

    #Handle Keyworded Arguments.
    alph = None
    accretio = None
    attritio = None
    for key in kwargs:
        if(key in set(['Exponent', 'Exp', 'Alpha', 
                       'alpha', 'exp', 'exponent'])):
            alph = kwargs[key]
        elif(key in set(['Accretion', 'accretion', 'accr', 'Accr'])):
            accretio = kwargs[key]
        elif(key in set(['Attrition', 'attrition', 'attr', 'Attr'])):
            attritio = kwargs[key]

    #Import Procedure Related Functions:
    if(CoalescenceKernel in set(['Constant', 'Pagett', '0', 0])):
        from JohnsonAggregator import PagettAggregate as Combine
    elif(CoalescenceKernel in set(['Multiplicative', 'Mult', '1', 1])):
        from JohnsonAggregator import SingleAggregate as Combine
    
    if(FragmentDistribution in 
       set(['Total', 'Shattering', 'Perfect', 'Atomize', 'Atomise'])
       ):
        
        if(FragmentKernel in set(['Constant', 'Pagett', '0', 0])):
            from JohnsonDisaggregator import PagettDisaggregate as Fragment
        elif(FragmentKernel in set(['Multiplicative', 'Mult', '1', 1])):
            from JohnsonDisaggregator import SingleDisaggregate as Fragment
            
    elif(FragmentDistribution in 
         set(['Chinese', 'ChineseRestaurant', 'Chinese_Restaurant', 'CRP'])
         or (FragmentDistribution in 
             set(['Power', 'PowerLaw', 'Power_Law', 'PL'])
             and alph != None and 1 <= alph and alph < 2)
         ):
        
        if(FragmentKernel in set(['Constant', 'Pagett', '0', 0])):
            from JohnsonDisaggregator import (
                    ImperfectFragmentation7_ConstantKernel as tempF
                    )
            from functools import partial
            Fragment = partial(tempF, alpha = alph - 1)
        elif(FragmentKernel in set(['Multiplicative', 'Mult', '1', 1])):
            from JohnsonDisaggregator import (
                    ImperfectFragmentation7 as tempF
                    )
            from functools import partial
            Fragment = partial(tempF, alpha = alph - 1)
            
    elif(FragmentDistribution in 
         set(['Price', 'PriceNetwork', 'Price_Network', 'PN'])
         or (FragmentDistribution in 
             set(['Power', 'PowerLaw', 'Power_Law', 'PL'])
             and alph != None and 2 < alph)
         ):
        
        if(FragmentKernel in set(['Constant', 'Pagett', '0', 0])):
            from JohnsonDisaggregator import (
                    ImperfectFragmentation8_ConstantKernel as tempF
                    )
            from functools import partial
            Fragment = partial(tempF, alpha = alph)
        elif(FragmentKernel in set(['Multiplicative', 'Mult', '1', 1])):
            from JohnsonDisaggregator import (
                    ImperfectFragmentation8 as tempF
                    )
            from functools import partial
            Fragment = partial(tempF, alpha = alph)
        
    #Determine the Fragmentation Threshold as a raw number.
    FragThreshold = ceil(len(NumberOfGroupsOfSize) * FragmentBarrier)
    
    if(debug):
        import pdb; pdb.set_trace()
    
    #Procedure:Record History if applicable
    if(History):
        returnHistory = list()
    if(Flux is not None and Flux is not False):
        CoalFlux = array([0] * (isum(PopAsList) - 1))
        FragFlux = array([0] * (isum(PopAsList) - 1))
    #On each timestep
    if(verbose):
        print('Beginning Simulation, Time = 0')
    NowTime = 0;
    if(History):
        HistorySinceLastSave = 0
        returnHistory.append(NumberOfGroupsOfSize)
        if(Flux is not None and Flux is not False):
            returnFlux = list()
            returnFlux.append((array(CoalFlux), array(FragFlux)))
    
    for NowTime in range(EndTime):
        #Check Coalescence and Fragmentaiton conditions.
        #Note in Fragmentation condition: We deduct one from the threshold
        #due to zero-indexing.
        AtLeast2Clusters = MinimumCount(NumberOfGroupsOfSize, 2)
        if(FragThreshold > 1):
            CanFragment = IndexAbove(FragThreshold - 1, NumberOfGroupsOfSize)
        else:
            CanFragment = True
    
        #Call appropriate action
        #Fragmentation nothing to coalesce or Probability < PFrag
        #Coalesce if there is nothing capable of fragmenting or 'otherwise'.
        #(Unless P(Nothing) is non-zero.)
        Act = ActionValue()
        if(not(AtLeast2Clusters) or (CanFragment and Act < PFrag)):
            if(verbose):
                print('Performing Fragment')
                
            if(FragThreshold > 1):
                NonCandidates = NumberOfGroupsOfSize[:(FragThreshold - 1)]
                #To work with existing infrastructure, we need to know how much
                #is in the applicable population and pass an appropriate list.
                NumCandidates = len(NumberOfGroupsOfSize) - isum(NonCandidates)
                FragCandidates = ([0] * (FragThreshold - 1) + 
                                  NumberOfGroupsOfSize[
                                          (FragThreshold - 1):NumCandidates])
                #Note, we can just trim off the end because it must be empty;
                #in order to be in the end, we would need mass from those in
                #the beginning of the list!
                FragResults = Fragment(FragCandidates)
                
                GroupsAfterFragment = (FragResults[0] + 
                                       [0] * (len(NumberOfGroupsOfSize) 
                                               - NumCandidates))
                NumberOfGroupsOfSize = [x + y for (x, y) in zip_longest(
                                        NonCandidates, GroupsAfterFragment,
                                        fillvalue = 0)]
            else:
                FragResults = Fragment(NumberOfGroupsOfSize)
                NumberOfGroupsOfSize = FragResults[0]
                
            if (Flux is not None and Flux is not False):
                Fragmented = sum(i*j for (i,j) in FragResults[1])
                if(verbose):
                    print(FragResults[1])
                if Fragmented != 1:
                    for (i, j) in FragResults[1]:
                        FragFlux[(i - 1) : (Fragmented - 1)] += i * j
            
        elif(not(CanFragment) or PFrag <= Act < PFrag + PCoal):
            if(verbose):
                print('Performing Combine')
            CoalResults = Combine(NumberOfGroupsOfSize)
            NumberOfGroupsOfSize = CoalResults[0]
            
            if (Flux is not None and Flux is not False):
                Combined = sum(i for i in CoalResults[1])
                if(verbose):
                    print(CoalResults[1])
                for i in CoalResults[1]:
                    CoalFlux[(i - 1) : (Combined - 1)] += i
            
        #Record History if applicable
        if(History):
            if(HistorySinceLastSave + 1 == HistorySteps):
                returnHistory.append(NumberOfGroupsOfSize)
                HistorySinceLastSave = 0
                if(Flux is not None and Flux is not False):
                    returnFlux.append((array(CoalFlux), array(FragFlux)))
                    if Flux != 'Sum':
                        FragFlux = array([0] * len(FragFlux))
                        CoalFlux = array([0] * len(CoalFlux))
            else:
                HistorySinceLastSave += 1
        if(verbose):
            print('Ending Timestep:', NowTime)

    if(History and not Flux):
        return returnHistory
    if(Flux and not History):
        return (CoalFlux, FragFlux)
    if(History and Flux):
        return (returnHistory, returnFlux)
    return NumberOfGroupsOfSize