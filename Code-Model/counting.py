# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 13:25:34 2016

@author: btf500
"""

#Return True if there are AtLeast or more non-zero elements in Vector.
def MinCount(Vector, AtLeast):
    '''Returns True if VECTOR has ATLEAST
    or more elements.
    '''
    #Guard against sparse matrices.
    from scipy.sparse import issparse
    if(issparse(Vector)):
        return SparseMinCount(Vector, AtLeast)
        
    if(AtLeast < 1):
        return True
    
    #Guard In case of Vector not a list
    if(isinstance(Vector, int)):
        if(AtLeast > 1):
            return False
        if(AtLeast == 1):
            return bool(Vector)
        
    
    #Trivial Cases
    if(AtLeast == 1):
        return any(Vector)
    if(AtLeast == len(Vector)):
        return all(Vector)
    if(AtLeast > len(Vector)):
        return False
        
    currentCount = 0
    for Item in Vector:
        if(Item):
            currentCount += 1
        if(currentCount >= AtLeast):
            return True
    return False
    
#Return True if the elements in Vector sum to AtLeast or more.
def MinSum(Vector, AtLeast):
    '''Returns True if VECTOR's elements
    sum to ATLEAST or more.
    '''
    #Guard against sparse matrices.
    from scipy.sparse import issparse
    if(issparse(Vector)):
        return SparseMinSum(Vector, AtLeast)
    
    #Guard In case of Vector not a list
    if(isinstance(Vector, int)):
        return Vector>=AtLeast
        
    currentCount = 0
    for Item in Vector:
        if(Item):
            currentCount += Item
        if(currentCount >= AtLeast):
            return True
    return False
    
def AdditiveFactors(Number):
    '''Randomly generate values in a list
    (with magnitude <= input)
    that will add up to the integer NUMBER
    and Return said list.
    '''
    from numpy.random import uniform as RandomNumber
    
    ActingNumber = int(Number)
    
    if(Number == 0):
        return [0]
    if(Number < 0):
        ActingNumber = -ActingNumber
    
    Return = []
    while(sum(Return) < ActingNumber):
        Return.append(int(RandomNumber(1, ActingNumber - sum(Return) + 1)))

    if(Number < 0):
        return [-x for x in Return]
    else:
        return Return
    
def AdditiveFactors_no_round(Number):
    '''Randomly generate values in a list
    (with magnitude <= input)
    that will add up to the integer NUMBER
    and Return said list.
    '''
    from numpy.random import randint as RandomNumber
    
    ActingNumber = int(Number)
    
    if(Number == 0):
        return [0]
    if(Number < 0):
        ActingNumber = -ActingNumber
    
    Return = []
    while(sum(Return) < ActingNumber):
        Return.append(RandomNumber(1, ActingNumber - sum(Return) + 1))

    if(Number < 0):
        return [-x for x in Return]
    else:
        return Return
    
def BetaBinomial(n, a = 1, b = 1, size = None):
    '''Draws p from a Beta(a, b) distribution
    in order to draw from a Binomial(n, p). '''
    #Speed boost by only drawing one random number
    from numpy.random import beta, binomial, randint
    if(a == 1 and b == 1):
        return randint(0, n + 1)
    else:
        return binomial(n, beta(a, b, size))
    
def AdditiveFactors_BB(Number, alpha = 1, beta = 1):
    '''Randomly generate values in a list
    (with magnitude <= input)
    that will add up to the integer NUMBER
    and Return said list.
    '''
    ActingNumber = int(Number)
    
    if(Number == 0):
        return [0]
    if(Number < 0):
        ActingNumber = -ActingNumber
    
    Return = []
    sumReturn = sum(Return)
    while(sumReturn < ActingNumber):
        Return.append(BetaBinomial(ActingNumber - sumReturn - 1, alpha, beta) + 1)
        sumReturn = sum(Return)

    if(Number < 0):
        return [-x for x in Return]
    else:
        return Return    

def isum(X, iScale = 1):
    '''Sums from i = 1 to (end of X) of X_i*i**iScale, iScale defaults to 1
    '''
    from scipy.sparse import issparse
    if(issparse(X)):
        return sparsesum(X, iScale)
        
    if(isinstance(X,int)):
        return X
	
    summand = 0
    i = 1
    for x_i in X:
        summand += x_i*i**iScale
        i += 1
    return summand

def IndexAbove(index, targetlist):
    '''
    Checks for a nonzero entry in targetlist in a position at least index.
    E.g. 3, [0,0,0,0,0,1,0,0,2] returns True.
    E.g. 3, [1,2,3,0,0,0,0,0,0] returns False. (Remember we zero-index.)
    '''
    if(index < 0):
        return(any([x > 0 for x in targetlist]))
    elif(index >= len(targetlist)):
        return(False) 
    else:
        return any([x > 0 for x in targetlist[index:]])
    

#FOLLOWING NOT TO BE CALLED DIRECTLY!
def SparseMinSum(Matrix, AtLeast):
    from scipy.sparse import find
    return sum(find(Matrix)[2]) >= AtLeast
        
    
def SparseMinCount(Matrix, AtLeast):
    return Matrix.size >= AtLeast
    
 
def ijsum(X, iScale = 1, jScale = 1):
    '''Sums from i,j = 1 to (end of X) of 
    X[i][j]* (i**iScale) * (j**jScale), i,Scale defaults to 1
    '''
    from scipy.sparse import issparse
    if(issparse(X)):
        return sparsesum(X, iScale, jScale)
        
    if(isinstance(X,int)):
        return X
	
    summand = 0
    for i in range(len(X)):
        for j in range(len(X[i])):
            summand += X[i][j] * (i**iScale) * (j**jScale)
    
    return summand
    
    
def sparsesum(X, iScale = 1, jScale = 1):
    '''Sums from i,j = 1 to (end of X) of 
    X[i][j]* (i**iScale) * (j**jScale), i,Scale defaults to 1
    '''
    from scipy.sparse import find
    XElements = find(X)
    
    summand = 0
    #For each elements
    for i in range(len(XElements[0])):
        summand += ((XElements[0][i]+1)**iScale) * ((XElements[1][i]+1)**jScale) * XElements[2][i]
    return summand

    
#http://stackoverflow.com/questions/2161406/how-do-i-generate-a-uniform-random-integer-partition
# answered Jan 29 '10 at 17:24 by Jason Orendorff
# site accessed at 13:40 on Jan 23 '17
# Issue with this solution is the maximum recursion depth.
import random

cache = {}

def count_partitions(n, limit):
    if n == 0:
        return 1
    if (n, limit) in cache:
        return cache[n, limit]
    x = cache[n, limit] = sum(count_partitions(n-k, k) for k in range(1, min(limit, n) + 1))
    return x

def random_partition(n):
    a = []
    limit = n
    total = count_partitions(n, limit)
    which = random.randrange(total)
    while n:
        for k in range(1, min(limit, n) + 1):
            count = count_partitions(n-k, k)
            if which < count:
                break
            which -= count
        a.append(k)
        limit = k
        n -= k
    return a
    
#From the same source:
# answered Nov 7 '13 at 6:46 by Stephen DeSalvo

from numpy import exp, pi, sqrt
from numpy.random import geometric
def random_partition_large(n):    
    '''Fristedt's algorithm. Found on
    https://stackoverflow.com/a/19829615
    '''
    
    if(n == 1):
        return [1]
    GeometricParam =  exp(-pi/sqrt(6*n))
    GeometricVars = [0]*n
    while(not(isum(GeometricVars)==n)):
        #Geometric support in numpy is N\{0}, and we want N+{0}
        #Python traverses from 0:(n-1), we need from 1:n.
        GeometricVars = [geometric(1-GeometricParam**(i+1))-1 for i in range(n)]
    return GeometricVars

def random_partition_Probabilistic_Divide_Conquer_1(n):
    '''Based on Algorithm 5 from 
    https://arxiv.org/pdf/1110.3856.pdf
    '''
    
    if(n == 1):
        return [1]
    from numpy.random import rand
    GeometricParam =  exp(-pi/sqrt(6*n))
    GeometricVars = [0]*n
    while(1):
        #Geometric support in numpy is N\{0}, and we want N+{0}
        #Python traverses from 0:(n-1), we need from 1:n.
        GeometricVars = [0] + [geometric(1-GeometricParam**(i+1))-1 
                                for i in range(n) if i]
        k = n - isum(GeometricVars)
        if(k < 0 or rand() >= exp(-k * pi / sqrt(6 * n))):
            continue
        else:
            GeometricVars[0] = k
            break
    return GeometricVars

#Claims that the following is inefficient and you only need
#to accept with probability P(X = k)/max_l P(X = l).
#def random_partition_Probabilistic_Divide_Conquer_Fix(n):
#    '''Based on Algorithm 5 from 
#    https://arxiv.org/pdf/1110.3856.pdf
#    Fixed for proper P value?
#    '''
#    
#    if(n == 1):
#        return [1]
#    from numpy.random import rand
#    GeometricParam =  exp(-pi/sqrt(6*n))
#    GeometricVars = [0]*n
#    while(1):
#        #Geometric support in numpy is N\{0}, and we want N+{0}
#        #Python traverses from 0:(n-1), we need from 1:n.
#        GeometricVars = [0] + [geometric(1-GeometricParam**(i+1))-1 
#                                for i in range(n) if i]
#        k = n - isum(GeometricVars)
#        if(k < 0 or rand() >= exp(-k * pi / sqrt(6 * n))*(1-exp(-pi/sqrt(6*n)))):
#            continue
#        else:
#            GeometricVars[0] = k
#            break
#    return GeometricVars

#def DeSalvoEqn15(n):
#    return exp((-pi/sqrt(6*n))/sqrt(n))
#
#def DeSalvoEqn11(x, k):
#    return (1 - x**i)*((x**i)**k)
#
#def DeSalvoEqn5(n, j):
#    
#    
#def random_partition_Probabilistic_Divide_Conquer_2(n):
#    '''Based on Algorithm 6 from 
#    https://arxiv.org/pdf/1110.3856.pdf
#    '''
#    if(n == 1):
#        return [1]
#    
#    from numpy.random import rand
#    GeometricParam =  exp(-pi/sqrt(6*n)) #Eqn 15: x(n) = ...
#    GeometricVars = [0]*n
#    
#    k = -1
#    #Formality to specify the leave conditions up here.
#    while(k < 0 or k % 2):
#        #Compute Odd values.
#        #Geometric support in numpy is N\{0}, and we want N+{0}
#        #Python traverses from 0:(n-1), we need from 1:n.
#        GeometricVars = [geometric(1-GeometricParam**(i+1))-1 
#                         if (i + 1) % 2 else 0 for i in range(n) ]
#        k = n - isum(GeometricVars)
#        #Don't bother rejection calculation if k not right.
#        #Rejection calculation is probability that
#        #the even groups times their group sizes sum to k
#        #divided by the maximum probability of even 
#        #group sizes summing to an (even) value.
#        if(k < 0 or k % 2 or rand() < DeSalvoEqn5(n, k/2) / ):
#            continue
#    
#    Evens = random_partition_Probabilistic_Divide_Conquer_2(k/2)
#    GeometricVars = [GeometricVars[i] if (i + 1) % 2
#                     else Evens[i] for i in range(n)]
#    return GeometricVars

def random_partition_large_sizes_indexed_by_group(n):
    from numpy import repeat
    GeometricParam = exp(-pi/sqrt(6*n))
    sumvalue = 0
    RetVal = list()
    while(not(sumvalue == n)):
        sumvalue = 0
        candidates = list()
        for i in range(n):
            temp = geometric(1-GeometricParam**(i+1))-1
            if(temp):
                candidates.append(((i + 1), temp))
                sumvalue += temp*(i+1)
                
    for j in range(len(candidates)):
        RetVal.extend(repeat(candidates[j][0], candidates[j][1]))
    return RetVal

def Mathematica_Partition(n):
    from subprocess import check_output as Command
    from platform import system as osname
    #Run the Mathematica code, then remove extraneous feedback
    #and then convert the csv's to a list of ints.
    if(osname() == 'Linux'):
        execommand = 'math'
    elif(osname() == 'Windows'):
        execommand = 'math'
    temp = Command(
            [execommand, '-noprompt', '-run', 
             'Needs["Combinatorica`"];' +
             'Print[Combinatorica`RandomPartition[' + 
             str(int(n)) + 
             ']]; Exit[]']
            ).decode().split('{')[1].split('}')[0]
    #Recommendation from https://stackoverflow.com/a/19334399
    RetVal = [int(s) for s in temp.split(',')]
    return RetVal
