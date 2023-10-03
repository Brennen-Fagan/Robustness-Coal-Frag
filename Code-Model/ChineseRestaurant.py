# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 16:16:58 2017

Chinese Restaurant Dirichlet Process

@author: btf500
"""

def Dirichlet_Process(n,alpha,theta):
    '''Constructs an exchangeable 
    (Alpha, Theta) partition of n
    by the Chinese Restaurant method
    of Dong, Goldschmidt, and Martin
    in Coagulation-Fragmentation
    Duality, (2005).'''
    return Chinese_Restaurant(n,alpha,theta)

from numpy.random import uniform as RandomNumber
from numpy.random import seed as seed
    
def Chinese_Restaurant(Customers, Asocial, Table_making, debug = 0, Seed = 0):
    '''
    Returns power law with exponent
    1 + Asocial.
    '''
    if(debug):
        import pdb; pdb.set_trace()
    if(Seed):
        seed(Seed)
    if(Asocial >= 1):
        print('Alpha/Asocial >= 1')
        return False
    if(Asocial < 0):
        print('Alpha/Asocial < 0')
        return False
    if(Table_making < -Asocial):
        print('Theta/Table_making < - Alpha/Asocial')
        return False
    
    #Put the first person at the first table
    Tables = list()
    Tables.append(1)
    
    #Add the next people repetitively.
    ps = RandomNumber(0,1,Customers-1)
    for NewCustomer in range(Customers-1):
        p = ps[NewCustomer]
        
        NormalizingConstant = NewCustomer+1+Table_making
        NewTableValue = Table_making + len(Tables)*Asocial
        
        p = p - NewTableValue/NormalizingConstant
        
        if(p<0):
            Tables.append(1)
            
        else:
            i = -1
            while(p>0):
                i = i + 1
                p = p - (Tables[i]-Asocial)/NormalizingConstant
            Tables[i] = Tables[i] + 1
    if(Seed):
        seed()
    return Tables
    
def Chinese_Restaurant_opt(Customers, Asocial, Table_making, debug = 0, bypass = 1, Seed = 0):
    if(~bypass):
        if(debug):
            import pdb; pdb.set_trace()
        if(Asocial >= 1):
            print('Alpha/Asocial >= 1')
            return False
        if(Asocial < 0):
            print('Alpha/Asocial < 0')
            return False
        if(Table_making < -Asocial):
            print('Theta/Table_making < - Alpha/Asocial')
            return False
    if(Seed):
        seed(Seed)
    #Put the first person at the first table
    Tables = [1]
    
    #Add the next people repetitively.
    ps = RandomNumber(0,1,Customers-1)
    NormalizingConstant = Table_making
    NewTableValue = Table_making + Asocial
    for NewCustomer in range(Customers-1):
        NormalizingConstant += 1
        
        p = ps[NewCustomer]*NormalizingConstant
        
        p = p - NewTableValue
        
        if(p<0):
            Tables.append(1)
            NewTableValue += Asocial
            
        else:
            i = -1
            while(p>0):
                i = i + 1
                p = p - (Tables[i] - Asocial)
            Tables[i] = Tables[i] + 1
    if(Seed):
        seed()
    return Tables
    
#For some bizarre reason that I am unable to pin down, this version is broken
#The next version does work it appears.
def Stick_Breaking(length, alpha, theta, debug = 0):
    if(debug):
        import pdb; pdb.set_trace()
    #For 0<=alpha < 1 and theta > -alpha
    #let B1, B2,... be independent random variables such that
    #    Bn ~ Beta(1-alpha, theta + n alpha)
    #for all n >= 1.
    #Let X_1 = B1 and X_n = (1-B1)...(1-B(n-1))Bn for n>=2.
    
    #We instead want integer breaks. So
    
    from numpy.random import beta as RandomNumber
    from numpy import prod
    
    #Perform one cycle to initialize
    Return = [] 
    complements = []
    counter = 1
    Return.append(RandomNumber(1-alpha, theta + counter*alpha))
    complements.append(1-Return[-1])
    while(sum(Return)<(length-1)/length):
        counter = counter + 1
        Return.append(RandomNumber(1-alpha, theta + counter*alpha)*prod(complements))
        complements.append(1-Return[-1])
    Return = [int(x*length) for x in Return if int(x*length)]
    if(length-sum(Return)):
        Return.append(length-sum(Return))

    return Return
    
def Stick_Breaking_as_written(length, alpha, theta, debug = 0):
    if(debug):
        import pdb; pdb.set_trace()
    from numpy.random import beta as RandomNumber
    from numpy import prod
        
    draws = []
    complements = []
    for i in range(length):
        draws.append(RandomNumber(1-alpha, theta + (i+1)*alpha))
        complements.append(1-draws[-1])
        
    Return = [draws[0]]
    for i in range(1,len(draws)):
        Return.append(draws[i]*prod(complements[0:i]))
        
    Return = [int(x*length) for x in Return if int(x*length)]
    if(length-sum(Return)):
        Return.append(length-sum(Return))

    return Return
    
def Stick_Breaking_2(length, alpha, theta, debug = 0):
    if(debug):
        import pdb; pdb.set_trace()
    from numpy.random import beta as RandomNumber
    from numpy import log2
    from numpy import prod
        
    draws = []
    complements = []
    for i in range(int(round(log2(length)*2))):
        draws.append(RandomNumber(1-alpha, theta + (i+1)*alpha))
        complements.append(1-draws[-1])
        
    Return = [draws[0]]
    for i in range(1,len(draws)):
        Return.append(draws[i]*prod(complements[0:i]))
        
    Return = [int(x*length) for x in Return if int(x*length)]
    if(length-sum(Return)):
        Return.append(length-sum(Return))

    return Return
    
def Price_Network_opt_4(N, c, a):
    '''Starts with a pair of connected nodes.
    Returns a network with expected power law
    exponent of 2+a/c. (a > 0, c > 0, N Natural)
    '''
    from numpy.random import uniform as rng
    from numpy.random import randint as pick
    In_Degree = [0 for x in range(N)]
    Edges = [-1 for x in range(N*c)]
    
    In_Degree[0] = 1
    Edges[0] = 0
    In_Degree[1] = 1
    Edges[1] = 1
    EdgeNum = 2
    
    allrs = rng(size = (N-2, c)) < (c/(c+a))
    
    for i in range(2,N):
        rs = allrs[i-2]
        for r in rs:
            if(r):
                new_edge = Edges[pick(EdgeNum)]
            else:
                new_edge = pick(i)
                 
            In_Degree[new_edge] += 1
            Edges[EdgeNum] = new_edge
            EdgeNum += 1
            
    return In_Degree

def Price_Network(N, c, a):
    '''Starts with a node.
    Returns a network with expected power law
    exponent of 2+a/c. (a > 0, c > 0, N Natural)
    Implemented based upon Newman's Networks.
    '''
    from numpy.random import uniform as rng
    from numpy.random import randint as pick
    In_Degree = [0 for x in range(N)]
    Edges = [-1 for x in range(N*c)]
    
    In_Degree[0] = 1
    Edges[0] = 0
    EdgeNum = 1
    
    allrs = rng(size = (N-1, c)) < (c/(c+a))
    
    for i in range(1,N):
        rs = allrs[i-1]
        for r in rs:
            if(r):
                new_edge = Edges[pick(EdgeNum)]
            else:
                new_edge = pick(i)
                 
            In_Degree[new_edge] += 1
            Edges[EdgeNum] = new_edge
            EdgeNum += 1
            
    return In_Degree