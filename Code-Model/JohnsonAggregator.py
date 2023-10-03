# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 13:37:41 2016

@author: btf500
"""

def SingleAggregate(NumberOfGroupsOfSize):
    '''Combines two groups from
    NUMBEROFGROUPSOFSIZE. This process is 
    random with probability proportional 
    to the size of the groups.
    Assumes length of NUMBEROFGROUPSOFSIZE
    is the number of total agents in the
    population.
    '''
    
    from picker import PickProportionalToSize as PickCluster
    
    RetVal = list(NumberOfGroupsOfSize)
    
    #PickProportionalToSize Returns the Size of the Group. We want its index,
    #which is one less.
    ClusterPick1 = PickCluster(RetVal) - 1
    RetVal[ClusterPick1] += -1

    #Need to inform PickCluster of the Actual Population Size
    ClusterPick2 = PickCluster(RetVal, len(RetVal) - ClusterPick1 - 1) - 1
    RetVal[ClusterPick2] += -1

    RetVal[ClusterPick1 + ClusterPick2 + 1] += 1
    
    return (RetVal, (ClusterPick1 + 1, ClusterPick2 + 1))
    
def PagettAggregate(NumberOfGroupsOfSize):
    '''Combines two groups from
    NUMBEROFGROUPSOFSIZE. Groups are picked
    uniformly at random irrespective of size.
    Assumes length of NUMBEROFGROUPSOFSIZE
    is the number of total agents in the
    population.
    '''
    
    from picker import PickUniformlyAtRandom as PickCluster
    
    RetVal = list(NumberOfGroupsOfSize)
    
    #PickProportionalToSize Returns the Size of the Group. We want its index,
    #which is one less.
    ClusterPick1 = PickCluster(RetVal) - 1
    RetVal[ClusterPick1] += -1

    #Need to inform PickCluster of the Actual Population Size
    ClusterPick2 = PickCluster(RetVal) - 1
    RetVal[ClusterPick2] += -1

    RetVal[ClusterPick1 + ClusterPick2 + 1] += 1
    
    return (RetVal, (ClusterPick1 + 1, ClusterPick2 + 1))
    
def SparseSingleAggregate(NumberOfGroupsOfSize):
    '''Combines two groups from
    NUMBEROFGROUPSOFSIZE. This process is 
    random with probability proportional 
    to the size of the groups.
    Assumes length of NUMBEROFGROUPSOFSIZE
    is the number of total agents in the
    population.
    '''
    #import pdb; pdb.set_trace()
    from picker import PickProportionalToSize as PickCluster
    from scipy.sparse import lil_matrix
    
    RetVal = lil_matrix(NumberOfGroupsOfSize.copy())
    
    #PickProportionalToSize Returns the Size of the Group. We want its index,
    #which is one less.
    ClusterPick1 = PickCluster(RetVal) - 1
    RetVal[(0,ClusterPick1)] += -1

    #Automatically updated Population Value
    ClusterPick2 = PickCluster(RetVal) - 1
    RetVal[(0,ClusterPick2)] += -1

    RetVal[(0,ClusterPick1 + ClusterPick2 + 1)] += 1
    
    return (RetVal, (ClusterPick1 + 1, ClusterPick2 + 1))
    
def IndividualsAggregate(NumberOfGroupsOfSize, Accretion, 
                         PopulationSize, Unique = 0):
    '''From NUMBEROFGROUPSOFSIZE combine
    ACCRETION unique groups with individuals.
    This process is performed in the same manner
    as the SingleAggregate process; i.e.
    This process is random with probability 
    proportional to the size of the groups.
    Assumes length of NUMBEROFGROUPSOFSIZE
    is the number of total agents in the
    population.
    '''
    
    from picker import PickNTimesProportionalToSize as PickClusters
        
    RetVal = list(NumberOfGroupsOfSize)
    if(RetVal[0] == 0):
        return RetVal
    
    if(Unique):
        #PickProportionalToSize Returns the Size of the Group. We want its index,
        #which is one less.
        ClustersPicked = PickClusters(RetVal, Accretion, PopulationSize, Unique) 
        
        for Size in ClustersPicked:
            
            if(Size < 0):
                continue
            
            #If at least two individuals, we're safe
            if(RetVal[0]>2):
                RetVal[Size-1] += -1
                RetVal[0] += -1
                RetVal[Size] += 1
            #If one individual, we need to check to make sure we're not combining
            #the individual with his/herself.
            elif(RetVal[0]==1):
                #If group has size at least 2
                if(Size>1):
                    RetVal[Size-1] += -1
                    RetVal[0] += -1
                    RetVal[Size] += 1
                #else, trying to combine the singleton with his/herself. Do nothing
            else:
                break
            
    else:
        #Need to do the same idea for disaggregation, but in reverse.
        from numpy.random import choice as Pick
        
        PopSize = PopulationSize
        if(PopulationSize == -1):
            from counting import isum
            PopSize = isum(NumberOfGroupsOfSize)
        
        if(Accretion >= 1):
            NumberAccreted = int(Accretion)
        elif(Accretion < 1 and Accretion > 0):
            NumberAccreted = int(Accretion * PopSize)
        elif(Accretion <= 0):
            return RetVal
        
        NumberAccreted = min(NumberAccreted, RetVal[0])
        
        #Disallow self-aggregation: would not appear in chemical representation.
        #=> PopSize - 1 agents to accrete with.
        #Remove Active Agent
        RetVal[0] += -1
        #Calculate weights for choices from different groups.
        Weights = [(Size + 1) * num / (PopSize - 1) for (Size, num) in zip(range(PopSize), RetVal)]
        
                   
        for i in range(NumberAccreted):
            if(RetVal[0] == -1):
                #No Individuals, not even an active one.
                break
            
            GroupSizePicked = Pick(range(PopSize), p = Weights)

            #Coalesce, automatically switch active individual.
            RetVal[0] += -1
            RetVal[GroupSizePicked] += -1
            RetVal[GroupSizePicked + 1] += 1

            #Recalculate Weights.
            #(Note if Weights[0] < 0, we'll break out on next iteration)
            #(Note that we have to deal with round-off error about 0.)
            if(RetVal[0] == 0):
                Weights[0] = 0
            else:
                Weights[0] += -1 / (PopSize - 1)
            if(RetVal[GroupSizePicked] == 0):
                Weights[GroupSizePicked] = 0
            else:
                Weights[GroupSizePicked] += -(GroupSizePicked + 1) / (PopSize - 1)
            Weights[GroupSizePicked + 1] += (GroupSizePicked + 1 + 1) / (PopSize - 1)

        RetVal[0] += 1
    
    return (RetVal, )
    
def MultipleAggregate(NumberOfGroupsOfSize, NumberToAggregate, PopulationSize):
    '''From NUMBEROFGROUPSOFSIZE, combines
    NUMBERTOAGGREGATE groups into one group.
    This process is performed in the same manner
    as the SingleAggregate process; i.e.
    this process is random with probability
    proportional to the size of the group.
    Assumes length of NUMBEROFGROUPSOFSIZE
    is the number of total agents in the
    population.
    '''
    from picker import PickNTimesProportionalToSize as PickClusters
        
    RetVal = list(NumberOfGroupsOfSize)
    
    #PickProportionalToSize Returns the Size of the Group. 
    ClustersPicked = PickClusters(RetVal, NumberToAggregate, PopulationSize) 
    
    for Size in ClustersPicked:
        RetVal[Size-1] += -1

    RetVal[sum(ClustersPicked) -1] += 1
    
    return (RetVal, (ClustersPicked))
    
def ImperfectAggregate(NumberOfGroupsOfSize):
    '''Combines two groups from
    NUMBEROFGROUPSOFSIZE imperfectly. The
    larger group stays intact and gains at least
    one and the remainder are additively fragmented.
    This process is random with probability 
    proportional to the size of the groups.
    Assumes length of NUMBEROFGROUPSOFSIZE
    is the number of total agents in the
    population.
    '''
    
    from picker import PickProportionalToSize as PickCluster
    from numpy.random import uniform as Gain
    from counting import AdditiveFactors as BreakUp
    
    RetVal = list(NumberOfGroupsOfSize)
    
    #PickProportionalToSize Returns the Size of the Group. 
    ClusterPick1 = PickCluster(RetVal)
    RetVal[ClusterPick1-1] += -1

    #Need to inform PickCluster of the Actual Population Size
    ClusterPick2 = PickCluster(RetVal, len(RetVal) - ClusterPick1)
    RetVal[ClusterPick2-1] += -1

    #Determine Size
    LargePick = ClusterPick1
    SmallPick = ClusterPick2
    if(LargePick<ClusterPick2):
        LargePick = ClusterPick2
        SmallPick = ClusterPick1
    
    #Determine Large's Gain. Note that int will round down and
    #we always gain 1. So our range is between (1, SmallPick+1)
    #=> [1, SmallPick]
    LargeGain = int(Gain(1,SmallPick+1))
    
    #Assign LargeGroup and its gain
    RetVal[LargePick+LargeGain - 1] += 1

    #Assign the small groups remaining
    SmallPick += - LargeGain
    if(SmallPick):
        SmallGroups = BreakUp(SmallPick)
        for i in SmallGroups:
            RetVal[i-1] += 1
    
    return (RetVal, )
    
    
def SafeAggregate(NumberOfGroupsOfSize):
    '''Combines two groups from
    NUMBEROFGROUPSOFSIZE. This process is 
    random with probability proportional 
    to the size of the groups. Additionally,
    this variant is more computationally 
    expensive but more guarded than normal.
    Assumes length of NUMBEROFGROUPSOFSIZE
    is the number of total agents in the
    population.
    '''
    
    from picker import PickProportionalToSize as PickCluster
    from counting import isum
    
    RetVal = list(NumberOfGroupsOfSize)
    NumPop = isum(RetVal)
    
    #PickProportionalToSize Returns the Size of the Group.
    ClusterPick1 = PickCluster(RetVal, NumPop) 
    RetVal[ClusterPick1 - 1] += -1
    NumPop += -ClusterPick1

    #Need to inform PickCluster of the Actual Population Size
    ClusterPick2 = PickCluster(RetVal, NumPop)
    RetVal[ClusterPick2 - 1] += -1

    RetVal[ClusterPick1 + ClusterPick2 - 1] += 1
    
    return (RetVal, (ClusterPick1 + 1, ClusterPick2 + 1))
    
def SingleAggregate2(NumberOfGroupsOfSize, PopulationSize):
    '''Combines two groups from
    NUMBEROFGROUPSOFSIZE. Requires the
    POPULATIONSIZE for safety. This process is 
    random with probability proportional 
    to the size of the groups.
    Assumes length of NUMBEROFGROUPSOFSIZE
    is the number of total agents in the
    population.
    '''
    
    from picker import PickProportionalToSize as PickCluster
    
    RetVal = list(NumberOfGroupsOfSize)
    
    #PickProportionalToSize Returns the Size of the Group. 
    ClusterPick1 = PickCluster(RetVal, PopulationSize)
    RetVal[ClusterPick1 - 1] += -1

    #Need to inform PickCluster of the Actual Population Size
    ClusterPick2 = PickCluster(RetVal, PopulationSize - ClusterPick1)
    RetVal[ClusterPick2 - 1] += -1

    RetVal[ClusterPick1 + ClusterPick2 -1] += 1
    
    return (RetVal, (ClusterPick1 + 1, ClusterPick2 + 1))