# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 13:40:59 2016

@author: btf500
"""

def SingleDisaggregate(NumberOfGroupsOfSize):
    '''Fragments a single group from
    NUMBEROFGROUPSOFSIZE into individuals.
    This process is random with probability
    proportional to the size of the group.
    Assumes length of NUMBEROFGROUPSOFSIZE
    is the number of total agents in the
    population.
    '''
    
    from picker import PickProportionalToSize as PickCluster
    
    RetVal = list(NumberOfGroupsOfSize)
    
    #PickProportionalToSize Returns the Size of the Group. We want its index,
    #which is one less.
    ClusterPick = PickCluster(RetVal) - 1

    RetVal[ClusterPick] += -1
    RetVal[0] += ClusterPick + 1
    
    return (RetVal, ((1, ClusterPick + 1), )) 
    #       (Population, ((Size, Multiplicity), ...))
    
def PagettDisaggregate(NumberOfGroupsOfSize):
    '''Fragments a single group from
    NUMBEROFGROUPSOFSIZE into individuals.
    The Group is picked uniformly at random 
    irrespective of size.
    Assumes length of NUMBEROFGROUPSOFSIZE
    is the number of total agents in the
    population.
    '''
    
    from picker import PickUniformlyAtRandom as PickCluster
    
    RetVal = list(NumberOfGroupsOfSize)
    
    #PickProportionalToSize Returns the Size of the Group. We want its index,
    #which is one less.
    ClusterPick = PickCluster(RetVal) - 1

    RetVal[ClusterPick] += -1
    RetVal[0] += ClusterPick + 1
    
    return (RetVal, ((1, ClusterPick + 1), )) 
    #       (Population, ((Size, Multiplicity), ...))
    
def SparseSingleDisaggregate(NumberOfGroupsOfSize):
    '''Fragments a single group from
    NUMBEROFGROUPSOFSIZE into individuals.
    This process is random with probability
    proportional to the size of the group.
    Assumes length of NUMBEROFGROUPSOFSIZE
    is the number of total agents in the
    population.
    '''
    
    from picker import PickProportionalToSize as PickCluster
    from scipy.sparse import lil_matrix
    
    RetVal = lil_matrix(NumberOfGroupsOfSize.copy())
    
    #PickProportionalToSize Returns the Size of the Group. We want its index,
    #which is one less.
    ClusterPick = PickCluster(RetVal) - 1

    RetVal[(0,ClusterPick)] += -1
    RetVal[(0,0)] += ClusterPick + 1
    
    return (RetVal, ((1, ClusterPick + 1), )) 
    #       (Population, ((Size, Multiplicity), ...))

def IndividualsDisaggregate(NumberOfGroupsOfSize, Attrition, 
                            PopulationSize = -1, Unique = 1):
    '''From NUMBEROFGROUPSOFSIZE, removes
    an individual from ATTRITION unique groups.
    This process is performed in the same manner
    as the SingleDisaggregate process; i.e.
    this process is random with probability
    proportional to the size of the group.
    Assumes length of NUMBEROFGROUPSOFSIZE
    is the number of total agents in the
    population.
    '''
    RetVal = list(NumberOfGroupsOfSize)
    
    if(Unique):
        from picker import PickNTimesProportionalToSize as PickClusters
        ClustersPicked = PickClusters(RetVal, Attrition, PopulationSize, Unique) 
        
        #PickProportionalToSize Returns the Size of the Group. 
        for Size in ClustersPicked:
            if(Size>1):
                RetVal[Size-1] += -1
                RetVal[0] += 1
                RetVal[Size-1-1] += 1 
            else:
                break
    else:
        from numpy.random import randint as Pick
        from numpy import unique
        
        PopSize = PopulationSize
        if(PopulationSize == -1):
            from counting import isum
            PopSize = isum(NumberOfGroupsOfSize)
        
        if(Attrition >= 1):
            NumberAttrited = int(Attrition)
        elif(Attrition < 1 and Attrition > 0):
            NumberAttrited = int(Attrition * PopSize)
        elif(Attrition <= 0):
            return (RetVal,)
            
        AgentsPicked = Pick(PopSize, size = NumberAttrited)
        
        #By Sorting Agents, we can then deal with non-uniqueness:
        #   i.e. Sort the agents, then traverse.
        #       for each agent found, check how many agents are
        #       within a group size of that agent.
        #       Deduct that many agents from the group, then
        #       place the group appropriately. This won't
        #       change the number of agents counted so far.
        #       Then just need to be careful to iterate until
        #       we've "caught up".
        #       Duplicate agents are removed once, twice has no
        #       effect. Hence why we use unique instead of sort.
        
        AgentsPicked = unique(AgentsPicked)
        AgentsPickedNumber = 0
        
        AgentsSoFar = 0
        
        for Size in range(0, len(RetVal)):
            AgentsInGroupsOfSize = RetVal[Size] * (Size + 1)
            
            while(AgentsPicked[AgentsPickedNumber] < 
                  AgentsSoFar + AgentsInGroupsOfSize):
                #Identify Group + AgentsPicked in Group
                GroupFirstAgent = AgentsSoFar + (Size + 1) * int((AgentsPicked[AgentsPickedNumber] - AgentsSoFar) / (Size + 1)) 
                AgentsInGroup = 0
                while(GroupFirstAgent + (Size + 1) > AgentsPicked[AgentsPickedNumber + AgentsInGroup]):
                    AgentsInGroup += 1
                    if(AgentsPickedNumber + AgentsInGroup >= len(AgentsPicked)):
                        break
                
                RetVal[Size] += -1
                if(Size + 1 == AgentsInGroup):
                    RetVal[0] += AgentsInGroup
                else:
                    RetVal[Size - AgentsInGroup] += 1
                    RetVal[0] += AgentsInGroup
                
                AgentsPickedNumber += AgentsInGroup
                
                if(AgentsPickedNumber >= len(AgentsPicked)):
                    break
            
            if(AgentsPickedNumber >= len(AgentsPicked)):
                break
            
            AgentsSoFar += AgentsInGroupsOfSize

    return (RetVal, )
    
def ImperfectFragmentation1(NumberOfGroupsOfSize, a = 1, b = 1):
    '''Fragments a single group from
    NUMBEROFGROUPSOFSIZE into smaller groups,
    at least one of which is an individual.
    This process is random with probability
    proportional to the size of the group.
    Assumes length of NUMBEROFGROUPSOFSIZE
    is the number of total agents in the
    population.
    ALIAS: BranchingFragmentation
           SingleDisaggregateAdditively
    '''
    
    from picker import PickProportionalToSize as PickCluster
    if(a == 1 and b == 1):
        #Faster, but specialized
        from counting import AdditiveFactors as BreakUp
    else:
        #More generalized but approx. 3x slower
        from counting import AdditiveFactors_BB as BreakUp
    
    RetVal = list(NumberOfGroupsOfSize)
    
    #PickProportionalToSize Returns the Size of the Group. 
    ClusterPick = PickCluster(RetVal)
    #If we pick an individual, then fragmenting returns just an individual;
    if(ClusterPick == 1):
        return (RetVal,)
    RetVal[ClusterPick-1] += -1

    #At Least one of which is an individual
    ClusterPick += -1
    RetVal[0] += 1
    
    #The BreakUp Function returns the sizes of the new groups.
    if (a == 1 and b == 1):
        NewGroupSizes = BreakUp(ClusterPick)
    else:
        NewGroupSizes = BreakUp(ClusterPick, a, b)
        
    for i in NewGroupSizes:
        RetVal[i-1] += 1
    
    return (RetVal, ((i, 1) for i in NewGroupSizes))

SingleDisaggregateAdditively = ImperfectFragmentation1 
BranchingFragmentation = SingleDisaggregateAdditively   

def ImperfectFragmentation2(NumberOfGroupsOfSize):
    '''Fragments a single group from
    NUMBEROFGROUPSOFSIZE, into smaller groups,
    although it can fail.
    The group selected is random with probability
    proportional to the size of the group.
    The fragmentation picks a partition
    uniformly at random from the set of all 
    partitions of the size of the selected group.
    ALIAS: PartitioningFragmentation
    '''
    from picker import PickProportionalToSize as PickCluster
    from counting import random_partition_Probabilistic_Divide_Conquer_1 as BreakUp
    
    RetVal = list(NumberOfGroupsOfSize)
    
    #PickProportionalToSize Returns the Size of the Group. 
    ClusterPick = PickCluster(RetVal)
    #If we pick an individual, then fragmenting returns just an individual;
    if(ClusterPick == 1):
        return (RetVal,)
    RetVal[ClusterPick-1] += -1
        
    #The BreakUp Function returns the number of individuals, pairs, etc.
    NewGroupsOfSize = BreakUp(ClusterPick)

    LastHalfOfRetVal = RetVal[ClusterPick:]
    FirstHalfOfRetVal = [x+y for x, y in zip(NewGroupsOfSize, RetVal)]
    RetVal = FirstHalfOfRetVal + LastHalfOfRetVal
    
    return (RetVal,)
    
PartitioningFragmentation = ImperfectFragmentation2

def ImperfectFragmentation3(NumberOfGroupsOfSize):
    '''Fragments a single group from
    NUMBEROFGROUPSOFSIZE, into smaller groups,
    although it can fail.
    The group selected is random with probability
    proportional to the size of the group.
    The fragmentation simply splits the group
    into two groups. If odd, clearly one resultant
    group is larger by one than the other.
    ALIAS: HalvingFragmentation
    '''
    from picker import PickProportionalToSize as PickCluster
    from math import ceil
    
    RetVal = list(NumberOfGroupsOfSize)
    
    #PickProportionalToSize Returns the Size of the Group. 
    ClusterPick = PickCluster(RetVal)
    #If we pick an individual, then fragmenting returns just an individual;
    if(ClusterPick == 1):
        return (RetVal,)
    RetVal[ClusterPick-1] += -1

    Bigger=int(ceil(ClusterPick/2))
    Smaller = ClusterPick-Bigger
    
    RetVal[Bigger-1] += 1
    RetVal[Smaller-1] += 1
    return (RetVal,)
    
HalvingFragmentation = ImperfectFragmentation3

def ImperfectFragmentation4(NumberOfGroupsOfSize):
    '''Fragments a single group from
    NUMBEROFGROUPSOFSIZE, into smaller groups,
    although it can fail.
    The group selected is random with probability
    proportional to the size of the group.
    The fragmentation simply splits the group
    into two. The smaller split then is released
    as individuals.
    ALIAS: HalfAndReleaseFragmentation
    '''
    from picker import PickProportionalToSize as PickCluster
    from math import ceil
    
    RetVal = list(NumberOfGroupsOfSize)
    
    #PickProportionalToSize Returns the Size of the Group. 
    ClusterPick = PickCluster(RetVal)
    #If we pick an individual, then fragmenting returns just an individual;
    if(ClusterPick == 1):
        return (RetVal,)
    RetVal[ClusterPick-1] += -1

    Bigger=int(ceil(ClusterPick/2))
    Smaller = ClusterPick-Bigger
    
    RetVal[Bigger-1] += 1
    RetVal[0] += Smaller
    return (RetVal,)
    
HalfAndReleaseFragmentation = ImperfectFragmentation4

def ImperfectFragmentation5(NumberOfGroupsOfSize):
    '''Fragments a single group from
    NUMBEROFGROUPSOFSIZE, into smaller groups,
    although it can fail.
    The group selected is random with probability
    proportional to the size of the group.
    The fragmentation splits the group by removing
    the ceiling of one-tenth of its size and 
    returning those removed as individuals.
    ALIAS:  DecimatingFragmentation,
            Roman_Fragmentation
    '''
    from picker import PickProportionalToSize as PickCluster
    from math import ceil
    
    RetVal = list(NumberOfGroupsOfSize)
    
    #PickProportionalToSize Returns the Size of the Group. 
    ClusterPick = PickCluster(RetVal)
    #If we pick an individual, then fragmenting returns just an individual;
    if(ClusterPick == 1):
        return (RetVal,)
    RetVal[ClusterPick-1] += -1

    Removed=int(ceil(ClusterPick/10))
    Staying = ClusterPick-Removed
    
    if(Staying):
        RetVal[Staying-1] += 1
    RetVal[0] += Removed
    return (RetVal,)
    
DecimatingFragmentation = ImperfectFragmentation5
Roman_Fragmentation = ImperfectFragmentation5

def ImperfectFragmentation6(NumberOfGroupsOfSize):
    '''Fragments a single group from
    NUMBEROFGROUPSOFSIZE, into smaller groups,
    although it can fail.
    The group selected is random with probability
    proportional to the size of the group.
    The fragmentation splits the group by removing
    the ceiling of nine-tenths of its size and 
    returning those removed as individuals.
    ALIAS:  AntiDecimatingFragmentation,
            Anti_Roman_Fragmentation
    '''
    from picker import PickProportionalToSize as PickCluster
    from math import ceil
    
    RetVal = list(NumberOfGroupsOfSize)
    
    #PickProportionalToSize Returns the Size of the Group. 
    ClusterPick = PickCluster(RetVal)
    #If we pick an individual, then fragmenting returns just an individual;
    if(ClusterPick == 1):
        return (RetVal,)
    RetVal[ClusterPick-1] += -1

    Removed=int(ceil(ClusterPick/10*9))
    Staying = ClusterPick-Removed
    
    if(Staying):
        RetVal[Staying-1] += 1
    RetVal[0] += Removed
    return (RetVal,)
    
AntiDecimatingFragmentation = ImperfectFragmentation6
Anti_Roman_Fragmentation = ImperfectFragmentation6

def ImperfectFragmentation7(NumberOfGroupsOfSize, alpha = 0):
    '''Fragments a single group from
    NUMBEROFGROUPSOFSIZE into smaller groups,
    at least one of which is an individual.
    This process is random with probability
    proportional to the size of the group.
    Assumes length of NUMBEROFGROUPSOFSIZE
    is the number of total agents in the
    population.
    ALIAS: ChineseRestaurantFragmentation
    '''
    
    from picker import PickProportionalToSize as PickCluster
    from ChineseRestaurant import Chinese_Restaurant_opt as BreakUp
    
    RetVal = list(NumberOfGroupsOfSize)
    
    #PickProportionalToSize Returns the Size of the Group. 
    ClusterPick = PickCluster(RetVal)
    #If we pick an individual, then fragmenting returns just an individual;
    if(ClusterPick == 1):
        return (RetVal, ((1, 1),))
    RetVal[ClusterPick-1] += -1

    #At Least one of which is an individual
    ClusterPick += -1
    RetVal[0] += 1
    
    NewGroupSizes = BreakUp(ClusterPick, alpha, 1)

    for i in NewGroupSizes:
        RetVal[i-1] += 1
        
    # Don't forget the first guy.
    NewGroupSizes = NewGroupSizes + [1]
    
    return (RetVal, tuple((i, 1) for i in NewGroupSizes))

ChineseRestaurantFragmentation = ImperfectFragmentation7

def ImperfectFragmentation7_ConstantKernel(NumberOfGroupsOfSize, alpha = 0):
    '''Fragments a single group from
    NUMBEROFGROUPSOFSIZE into smaller groups,
    at least one of which is an individual.
    This process is done uniformly at random.
    Assumes length of NUMBEROFGROUPSOFSIZE
    is the number of total agents in the
    population.
    '''
    
    from picker import PickUniformlyAtRandom as PickCluster
    from ChineseRestaurant import Chinese_Restaurant_opt as BreakUp
    
    RetVal = list(NumberOfGroupsOfSize)
    
    #Returns the Size of the Group. 
    ClusterPick = PickCluster(RetVal)
    #If we pick an individual, then fragmenting returns just an individual;
    if(ClusterPick == 1):
        return (RetVal, ((1, 1),))
    RetVal[ClusterPick-1] += -1

    #At Least one of which is an individual
    ClusterPick += -1
    RetVal[0] += 1
    
    NewGroupSizes = BreakUp(ClusterPick, alpha, 1)

    for i in NewGroupSizes:
        RetVal[i-1] += 1
        
    # Don't forget the first guy.
    NewGroupSizes = NewGroupSizes + [1]
    
    return (RetVal, tuple((i, 1) for i in NewGroupSizes))

def ImperfectFragmentation8(NumberOfGroupsOfSize, alpha = 2.01):
    '''Fragments a single group from
    NUMBEROFGROUPSOFSIZE into smaller groups,
    at least one of which is an individual.
    This process is random with probability
    proportional to the size of the group.
    Assumes length of NUMBEROFGROUPSOFSIZE
    is the number of total agents in the
    population.
    ALIAS: PriceFragmentation
           PriceNetworkFragmentation
    '''
    
    from picker import PickProportionalToSize as PickCluster
    from ChineseRestaurant import Price_Network as BreakUp
    if(alpha<2):
        return -1
    
    RetVal = list(NumberOfGroupsOfSize)
    
    #PickProportionalToSize Returns the Size of the Group. 
    ClusterPick = PickCluster(RetVal)
    #If we pick an individual, then fragmenting returns just an individual;
    if(ClusterPick == 1):
        return (RetVal, ((1, 1),))
    if(ClusterPick == 2):
        RetVal[0] += 2
        RetVal[1] -= 1
        return (RetVal, ((1, 2),))
    RetVal[ClusterPick-1] += -1

    #At Least one of which is an individual
    ClusterPick += -1
    RetVal[0] += 1
    
    NewGroupSizes = BreakUp(ClusterPick, c = 1, a = alpha-2)

    for i in NewGroupSizes:
        if(i):
            RetVal[i-1] += 1
        
    # Don't forget the first guy.
    NewGroupSizes = NewGroupSizes + [1]
    
    return (RetVal, tuple((i, 1) for i in NewGroupSizes))

PriceFragmentation = ImperfectFragmentation8
PriceNetworkFragmentation = ImperfectFragmentation8

def ImperfectFragmentation8_ConstantKernel(NumberOfGroupsOfSize, alpha = 2.01):
    '''Fragments a single group from
    NUMBEROFGROUPSOFSIZE into smaller groups,
    at least one of which is an individual.
    This process is done uniformly at random.
    Assumes length of NUMBEROFGROUPSOFSIZE
    is the number of total agents in the
    population.
    '''
    
    from picker import PickUniformlyAtRandom as PickCluster
    from ChineseRestaurant import Price_Network as BreakUp
    if(alpha<2):
        return -1
    
    RetVal = list(NumberOfGroupsOfSize)
    
    #Returns the Size of the Group. 
    ClusterPick = PickCluster(RetVal)
    #If we pick an individual, then fragmenting returns just an individual;
    if(ClusterPick == 1):
        return (RetVal, ((1, 1),))
    if(ClusterPick == 2):
        RetVal[0] += 2
        RetVal[1] -= 1
        return (RetVal, ((1, 2),))
    RetVal[ClusterPick-1] += -1

    #At Least one of which is an individual
    ClusterPick += -1
    RetVal[0] += 1
    
    NewGroupSizes = BreakUp(ClusterPick, c = 1, a = alpha-2)

    for i in NewGroupSizes:
        if(i):
            RetVal[i-1] += 1
        
    # Don't forget the first guy.
    NewGroupSizes = NewGroupSizes + [1]
    
    return (RetVal, tuple((i, 1) for i in NewGroupSizes))
    
def MultipleDisaggregate(NumberOfGroupsOfSize, NumberToDisaggregate):
    '''From NUMBEROFGROUPSOFSIZE, fragments
    NUMBERTODISAGGREGATE groups into individuals.
    This process is performed in the same manner
    as the SingleDisaggregate process; i.e.
    this process is random with probability
    proportional to the size of the group.
    Assumes length of NUMBEROFGROUPSOFSIZE
    is the number of total agents in the
    population.
    '''
    from picker import PickNTimesProportionalToSize as PickClusters
        
    RetVal = list(NumberOfGroupsOfSize)
    
    #PickProportionalToSize Returns the Size of the Group. 
    ClustersPicked = PickClusters(RetVal, NumberToDisaggregate) 
    
    for Size in ClustersPicked:
        if(Size>1):
            RetVal[Size-1] += -1
            RetVal[0] += Size
    
    return (RetVal,)
    
def SafeDisaggregate(NumberOfGroupsOfSize):
    '''Fragments a single group from
    NUMBEROFGROUPSOFSIZE into individuals.
    This process is random with probability
    proportional to the size of the group.
    Additionally, this variant is more 
    computationally expensive and more guarded.
    Assumes length of NUMBEROFGROUPSOFSIZE
    is the number of total agents in the
    population.
    '''
    
    from picker import PickProportionalToSize as PickCluster
    from counting import isum
        
    
    RetVal = list(NumberOfGroupsOfSize)
    NumPop = isum(RetVal)
    
    #PickProportionalToSize Returns the Size of the Group.
    ClusterPick = PickCluster(RetVal, NumPop)

    RetVal[ClusterPick - 1] += -1
    RetVal[0] += ClusterPick
    
    return (RetVal, ((1, ClusterPick + 1), )) 
    
def SingleDisaggregate2(NumberOfGroupsOfSize, PopulationSize):
    '''Fragments a single group from
    NUMBEROFGROUPSOFSIZE into individuals.
    Additionally requires POPULATIONSIZE
    for safety.
    This process is random with probability
    proportional to the size of the group.
    Assumes length of NUMBEROFGROUPSOFSIZE
    is the number of total agents in the
    population.
    '''
    
    from picker import PickProportionalToSize as PickCluster
    
    RetVal = list(NumberOfGroupsOfSize)
    
    #PickProportionalToSize Returns the Size of the Group.
    ClusterPick = PickCluster(RetVal, PopulationSize)

    RetVal[ClusterPick - 1] += -1
    RetVal[0] += ClusterPick
    
    return (RetVal, ((1, ClusterPick + 1), )) 
    
def RemoveAnIndividualUAR(NumberOfGroupsOfSize, PopulationSize):
    '''Similar to IndividualsDisaggregate, but
    performs the singleton case. Takes
    NUMBEROFGROUPSOFSIZE and the total
    POPULATIONSIZE to determine which
    group to remove an individual from.
    The process is Uniformly At Random
    with respect to the agents in the population.
    '''
    
    from picker import PickProportionalToSize as PickCluster
    
    RetVal = list(NumberOfGroupsOfSize)
    
    #PickProportionalToSize Returns the Size of the Group.
    ClusterPick = PickCluster(RetVal, PopulationSize)

    if(ClusterPick > 1):
        RetVal[ClusterPick - 1] += -1
        RetVal[0] += 1
        RetVal[ClusterPick - 2] += 1
    
    return (RetVal,)