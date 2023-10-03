# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 12:16:23 2016

@author: btf500
"""

def RunBirnstiel(ModelNumber, FragAlpha = 1, Barrier = 0, 
                 Kernel = 'Multiplicative', Population = 10000, Time = 10000, 
                 PCoal = 0.99, HistorySteps = 1, VersionNumber = -1):
    #Select Model:
    if(FragAlpha == 2):
        from Birnstiel import CoFr_with_Barrier as Model
    else:
        from Birnstiel import CRP_CoFr_with_Barrier as Model

    from os import makedirs
    from os import getcwd
    from os.path import isdir
    from datetime import date
    from scipy.io import savemat
    from counting import isum

    DirectoryName = getcwd() + '/Data' + date.today().strftime('%y-%m-%d') + '/Birnstiel/CRPVal' + str(int(FragAlpha*100))+'/'
    if not isdir(DirectoryName):
        makedirs(DirectoryName)
        
    FileName = 'RandomPartition_Birnstiel_' + str(ModelNumber) +  '_alpha' + str(int(FragAlpha*100)) + '_barrier' + str(int(Barrier*100))

    if isinstance(Population, list):
        #From: https://stackoverflow.com/questions/19502378/python-find-first-instance-of-non-zero-number-in-list
        if(VersionNumber == -1):
            LastZero = '_Pos' + str(next((i for i, x in enumerate(Population) if x), None)+1)
        else:
            LastZero = '_Ver' + str(VersionNumber)
        PopName = isum(Population)
    else:
        if(VersionNumber == -1):
            LastZero = '_Pos' + str(0+1)
        else:
            LastZero = '_Ver' + str(VersionNumber)
        PopName = Population    
    FileName = FileName + '_' + str(PopName)+'_'+str(int(100*(1-PCoal)))+'percent'+LastZero + '_' + Kernel
    if(FragAlpha == 2): 
        if(isinstance(Population, list) or Population):
            ResultsFromRun = Model(Population = Population, PCoal=PCoal, EndTime = Time, History = 1, HistorySteps = HistorySteps, FragBarrier = Barrier, Kernel = Kernel)
        else:
            ResultsFromRun = Model(EndTime = Time, PCoal=PCoal, History = 1, HistorySteps = HistorySteps, FragBarrier = Barrier, Kernel = Kernel)
    else:
        if(isinstance(Population, list) or Population):
            ResultsFromRun = Model(Population = Population, PCoal=PCoal, EndTime = Time, History = 1, HistorySteps = HistorySteps, FragExponent = FragAlpha, FragBarrier = Barrier, Kernel = Kernel)
        else:
            ResultsFromRun = Model(PCoal=PCoal, EndTime = Time, History = 1, HistorySteps = HistorySteps, FragExponent = FragAlpha, FragBarrier = Barrier, Kernel = Kernel)

    savemat(DirectoryName + FileName, mdict = {FileName: ResultsFromRun}, do_compression = True)
    print(FileName + ' is Done!')

    #Return the last result in the Run
    return ResultsFromRun[-1]

def RunPagett(ModelNumber, Population = 10000, Time = 10000, 
             PCoal = 0.99, HistorySteps = 1):
        #Select Model:
    from JohnsonOnePopModels import PagettOnePop as Model

    from os import makedirs
    from os import getcwd
    from os.path import isdir
    from datetime import date
    from scipy.io import savemat
    from counting import isum

    DirectoryName = getcwd() + '/Data' + date.today().strftime('%y-%m-%d') + '/Pagett'
    if not isdir(DirectoryName):
        makedirs(DirectoryName)
        
    FileName = 'RandomPartition_Pag_'+str(ModelNumber)

    if isinstance(Population, list):
        #From: https://stackoverflow.com/questions/19502378/python-find-first-instance-of-non-zero-number-in-list
        #if(VersionNumber == -1):
        #    LastZero = '_Pos' + str(next((i for i, x in enumerate(Population) if x), None)+1)
        #else:
        #    LastZero = '_Ver' + str(Version)
        PopName = isum(Population)
    else:
        #if(VersionNumber == -1):
        #    LastZero = '_Pos' + str(0+1)
        #else:
        #    LastZero = '_Ver' + str(VersionNumber)
        PopName = Population    
    FileName = FileName + '_' + str(PopName)+'_'+str(int(100*(1-PCoal)))+'percent'#+LastZero

    if(isinstance(Population,list) or Population):
        ResultsFromRun = Model(Population = Population, PCoal=PCoal, EndTime = Time, History = 1, HistorySteps = HistorySteps)
    else:
        ResultsFromRun = Model(EndTime = Time, PCoal=PCoal, History = 1, HistorySteps = HistorySteps)

    savemat(DirectoryName + FileName, mdict = {FileName: ResultsFromRun}, do_compression = True)
    print(FileName + ' is Done!')

    #Return the last result in the Run
    return ResultsFromRun[-1]


def RunJ1CRP(ModelNumber, CRPalpha = 0, Population = 10000, Time = 10000, 
             PCoal = 0.99, HistorySteps = 1, VersionNumber = -1):
        #Select Model:
    from JohnsonOnePopModels import OnePopImperfectFrag7 as Model

    from os import makedirs
    from os import getcwd
    from os.path import isdir
    from datetime import date
    from scipy.io import savemat
    from counting import isum

    DirectoryName = getcwd() + '/Data' + date.today().strftime('%y-%m-%d') + '/CRP/Trial' + str(ModelNumber)+'/'
    if not isdir(DirectoryName):
        makedirs(DirectoryName)
        
    FileName = 'RandomPartition_CRP_'+str(ModelNumber) + '_alpha' + str(int(CRPalpha*100))

    if isinstance(Population, list):
        #From: https://stackoverflow.com/questions/19502378/python-find-first-instance-of-non-zero-number-in-list
        if(VersionNumber == -1):
            LastZero = '_Pos' + str(next((i for i, x in enumerate(Population) if x), None)+1)
        else:
            LastZero = '_Ver' + str(VersionNumber)
        PopName = isum(Population)
    else:
        if(VersionNumber == -1):
            LastZero = '_Pos' + str(0+1)
        else:
            LastZero = '_Ver' + str(VersionNumber)
        PopName = Population    
    FileName = FileName + str(PopName)+'_'+str(int(100*(1-PCoal)))+'percent'+LastZero

    if(isinstance(Population,list) or Population):
        ResultsFromRun = Model(Population = Population, PCoal=PCoal, EndTime = Time, History = 1, HistorySteps = HistorySteps, alpha = CRPalpha)
    else:
        ResultsFromRun = Model(EndTime = Time, PCoal=PCoal, History = 1, HistorySteps = HistorySteps, alpha = CRPalpha)

    savemat(DirectoryName + FileName, mdict = {FileName: ResultsFromRun}, do_compression = True)
    print(FileName + ' is Done!')

    #Return the last result in the Run
    return ResultsFromRun[-1]

def RunJ1Price(ModelNumber, Pricealpha = 2, Population = 10000, Time = 10000, 
             PCoal = 0.99, HistorySteps = 1, VersionNumber = -1):
        #Select Model:
    from JohnsonOnePopModels import OnePopImperfectFrag8 as Model

    from os import makedirs
    from os import getcwd
    from os.path import isdir
    from datetime import date
    from scipy.io import savemat
    from counting import isum

    DirectoryName = getcwd() + '/Data' + date.today().strftime('%y-%m-%d') + '/Price/Trial' + str(ModelNumber)+'/'
    if not isdir(DirectoryName):
        makedirs(DirectoryName)
        
    FileName = 'RandomPartition_Price_'+str(ModelNumber) + '_alpha' + str(int(Pricealpha*100))

    if isinstance(Population, list):
        #From: https://stackoverflow.com/questions/19502378/python-find-first-instance-of-non-zero-number-in-list
        if(VersionNumber == -1):
            LastZero = '_Pos' + str(next((i for i, x in enumerate(Population) if x), None)+1)
        else:
            LastZero = '_Ver' + str(VersionNumber)
        PopName = isum(Population)
    else:
        if(VersionNumber == -1):
            LastZero = '_Pos' + str(0+1)
        else:
            LastZero = '_Ver' + str(VersionNumber)
        PopName = Population    
    FileName = FileName + str(PopName)+'_'+str(int(100*(1-PCoal)))+'percent'+LastZero

    if(isinstance(Population,list) or Population):
        ResultsFromRun = Model(Population = Population, PCoal=PCoal, EndTime = Time, History = 1, HistorySteps = HistorySteps, alpha = Pricealpha)
    else:
        ResultsFromRun = Model(EndTime = Time, PCoal=PCoal, History = 1, HistorySteps = HistorySteps, alpha = Pricealpha)

    savemat(DirectoryName + FileName, mdict = {FileName: ResultsFromRun}, do_compression = True)
    print(FileName + ' is Done!')

    #Return the last result in the Run
    return ResultsFromRun[-1]

def RunJ1withStandardModelSelection(ModelNumber, Population = 10000, Time = 10000, 
                                    PCoal = 0.99, HistorySteps = 1, VersionNumber = -1):
    '''Chooses a model from the standard list
    of models according to MODELNUMBER, and
    then runs a !single! run of that model with
    the POPULATION (=10000), TIME (=10000),
    PCOAL (=0.99), and HISTORYSTEPS (=1) provided.
    Then saves it with the provided VERSIONNUMBER;
    if this is not provided, the first nonzero
    entry in the population is used instead.
    Standard Models:
        0 = Standard Must
        1 = Accretion (Ac = 3)
        2 = Attrition (At = 3)
        3 = Ac. + At. (Both = 3)
        4 = Imperfect1:Both
        5 = Imperfect Frag. 1
        6 = Imperfect Coal.
    '''
    
    #Select Model:
    if(ModelNumber%7==0):
        from JohnsonOnePopModels import JohnsonOnePopMustSucceed as Model
        FileName = ('RandomPartition_Must_'+str(ModelNumber)+'_')
    elif(ModelNumber%7 == 1):
        from JohnsonOnePopModels import OnePopWithAccretion as Model
        FileName = ('RandomPartition_Accr3_'+str(ModelNumber)+'_')
    elif(ModelNumber%7 == 2):
        from JohnsonOnePopModels import OnePopWithAttrition as Model
        FileName = ('RandomPartition_Attr3_'+str(ModelNumber)+'_')
    elif(ModelNumber%7 == 3):
        from JohnsonOnePopModels import OnePopWithAccrAttr as Model
        FileName = ('RandomPartition_AcAt3_'+str(ModelNumber)+'_')
    elif(ModelNumber%7 == 4):
        from JohnsonOnePopModels import OnePopImperfect as Model
        FileName = ('RandomPartition_ImpAl_'+str(ModelNumber)+'_')
    elif(ModelNumber%7 == 5):
        from JohnsonOnePopModels import OnePopImperfectFrag as Model
        FileName = ('RandomPartition_ImpFr_'+str(ModelNumber)+'_')
    elif(ModelNumber%7 == 6):
        from JohnsonOnePopModels import OnePopImperfectComb as Model
        FileName = ('RandomPartition_ImpCo_'+str(ModelNumber)+'_')
    else:
        print('ModelNumber%7 did not return a result!')
        return 1

    from os import makedirs
    from os import getcwd
    from os.path import isdir
    from datetime import date
    from scipy.io import savemat
    from counting import isum

    DirectoryName = getcwd() + '/Data' + date.today().strftime('%y-%m-%d') + '/StandardSuite/Trial' + str(ModelNumber)+'/'
    if not isdir(DirectoryName):
        makedirs(DirectoryName)

    if isinstance(Population, list):
        #From: https://stackoverflow.com/questions/19502378/python-find-first-instance-of-non-zero-number-in-list
        if(VersionNumber == -1):
            LastZero = '_Pos' + str(next((i for i, x in enumerate(Population) if x), None)+1)
        else:
            LastZero = '_Ver' + str(VersionNumber)
        PopName = isum(Population)
    else:
        if(VersionNumber == -1):
            LastZero = '_Pos' + str(0+1)
        else:
            LastZero = '_Ver' + str(VersionNumber)
        PopName = Population    
    FileName = FileName + str(PopName)+'_'+str(int(100*(1-PCoal)))+'percent'+LastZero

    if(isinstance(Population,list) or Population):
        ResultsFromRun = Model(Population = Population, PCoal=PCoal, EndTime = Time, History = 1, HistorySteps = HistorySteps)
    else:
        ResultsFromRun = Model(EndTime = Time, PCoal=PCoal, History = 1, HistorySteps = HistorySteps)

    savemat(DirectoryName + FileName, mdict = {FileName: ResultsFromRun}, do_compression = True)
    print(FileName + ' is Done!')

    #Return the last result in the Run
    return ResultsFromRun[-1]


def RunTrialsOfJ1PopModels(TrialNumber, Pop = 0, Time = 10000):
    #from JohnsonOnePopModels import JohnsonOnePopMustSucceed as Model1
    #from JohnsonOnePopModels import JohnsonOnePopMaySucceed as Model2
    #from JohnsonOnePopModels import OnePopWithAttrition as Model3
    #from JohnsonOnePopModels import OnePopWithAccretion as Model4
    #from JohnsonOnePopModels import OnePopWithAccrAttr as Model5
    from JohnsonOnePopModels import OnePopImperfect as Model6
    from JohnsonOnePopModels import OnePopImperfectFrag as Model7
    #from JohnsonOnePopModels import OnePopImperfectComb as Model8
    #from JohnsonOnePopModels import OnePopMultiStep as Model9
    #from JohnsonOnePopModels import UnstablePop as Model10
    #Storage List
    List =    [[]]*1000
    
    from os import makedirs
    from os import getcwd
    from os.path import isdir
    from datetime import date
    from scipy.io import savemat
    from counting import isum
    #C:\User Files\Local Data
    #H:\Python\Johnson03
    #DirectoryName = 'C:\\User Files\\Local Data\\Data'+date.today().strftime('%y-%m-%d')+'\\Trial'+str(TrialNumber)+'\\'
    #DirectoryName = 'H:\\Python\\Johnson03\\Data'+date.today().strftime('%y-%m-%d')+'\\Trial'+str(TrialNumber)+'\\'
    #DirectoryName = '~\\w2k\\Python\\Johnson03\\Data'+date.today().strftime('%y-%m-%d')+'\\Trial'+str(TrialNumber)+'\\'
    DirectoryName = getcwd() + '/Data' + date.today().strftime('%y-%m-%d') + '/Trial' + str(TrialNumber)+'/'
    if not isdir(DirectoryName):
        makedirs(DirectoryName)

    if isinstance(Pop, list):
        PopName = isum(Pop)
    else:
        PopName = Pop
    
    #for i in range(1000):
    #    if(Pop):
    #        List[i] = Model1(Population = Pop, EndTime = Time)
    #    else:
    #        List[i] = Model1(EndTime = Time)
    #    
    #FileName = 'StandardMust' + str(TrialNumber)+'_'+str(Pop)
    #savemat(DirectoryName + FileName, mdict = {FileName: List}, do_compression = True)
    #print(FileName + ' is Done!')
    
#    for i in range(1000):
#        List[i] = Model2()
#        
#    FileName = 'StandardMay'+ str(TrialNumber)
#    savemat(DirectoryName + FileName, mdict = {FileName: List}, do_compression = True)
#    print(FileName + ' is Done!')
    
    #for i in range(1000):
    #    if(Pop):
    #        List[i] = Model3(Population = Pop, EndTime = Time)
    #    else:
    #        List[i] = Model3(EndTime = Time)
    #    
    #FileName = 'Attrition'+ str(TrialNumber)+'_'+str(Pop)
    #savemat(DirectoryName + FileName, mdict = {FileName: List}, do_compression = True)
    #print(FileName + ' is Done!')
    
    #for i in range(1000):
    #    if(Pop):
    #        List[i] = Model4(Population = Pop, EndTime = Time)
    #    else:
    #        List[i] = Model4(EndTime = Time)
    #    
    #FileName = 'Accretion'+ str(TrialNumber)+'_'+str(Pop)
    #savemat(DirectoryName + FileName, mdict = {FileName: List}, do_compression = True)
    #print(FileName + ' is Done!')
    
    #for i in range(1000):
    #    if(Pop):
    #        List[i] = Model5(Population = Pop, EndTime = Time)
    #    else:
    #        List[i] = Model5(EndTime = Time)
    #    
    #FileName = 'AccrAttr'+ str(TrialNumber)+'_'+str(Pop)
    #savemat(DirectoryName + FileName, mdict = {FileName: List}, do_compression = True)
    #print(FileName + ' is Done!')
    
    for i in range(1000):
        if(Pop):
            List[i] = Model6(Population = Pop, EndTime = Time)
        else:
            List[i] = Model6(EndTime = Time)
        
    FileName = 'ImperfectAll_'+ str(TrialNumber)+'_'+str(PopName)
    savemat(DirectoryName + FileName, mdict = {FileName: List}, do_compression = True)
    print(FileName + ' is Done!')
    
    for i in range(1000):
        if(Pop):
            List[i] = Model7(Population = Pop, EndTime = Time)
        else:
            List[i] = Model7(EndTime = Time)
        
    FileName = 'ImperfectFrag_'+ str(TrialNumber)+'_'+str(PopName)
    savemat(DirectoryName + FileName, mdict = {FileName: List}, do_compression = True)
    print(FileName + ' is Done!')
    
    #for i in range(1000):
    #    if(Pop):
    #        List[i] = Model8(Population = Pop, EndTime = Time)
    #    else:
    #        List[i] = Model8(EndTime = Time)
    #    
    #FileName = 'ImperfectComb'+ str(TrialNumber)+'_'+str(Pop)
    #savemat(DirectoryName + FileName, mdict = {FileName: List}, do_compression = True)
   # print(FileName + ' is Done!')
#    
#    for i in range(1000):
#        if(Pop):
#            List[i] = Model9(Population = Pop, EndTime = Time)
#        else:
#            List[i] = Model9(EndTime = Time)
#        
#    FileName = 'MultiSteps'+ str(TrialNumber)+'_'+str(Pop)
#    savemat(DirectoryName + FileName, mdict = {FileName: List}, do_compression = True)
#    print(FileName + ' is Done!')
#    
#    for i in range(1000):
#        if(Pop):
#            List[i] = Model10(Population = Pop, EndTime = Time)
#        else:
#            List[i] = Model10(EndTime = Time)
#        
#    FileName = 'UnstablePop'+ str(TrialNumber)+'_'+str(Pop)
#    savemat(DirectoryName + FileName, mdict = {FileName: List}, do_compression = True)
#    print(FileName + ' is Done!')
    
def RunTrialsOfImpFrag2Models(TrialNumber, Pop = 0, Time = 10000):
    from JohnsonOnePopModels import OnePopImperfect2 as Model7
    from JohnsonOnePopModels import OnePopImperfectFrag2 as Model6
    List =    [[]]*1000

    from counting import isum
    from os import makedirs
    from os import getcwd
    from os.path import isdir
    from datetime import date
    from scipy.io import savemat

    #C:\User Files\Local Data
    #H:\Python\Johnson03
    #DirectoryName = 'C:\\User Files\\Local Data\\Data'+date.today().strftime('$
    #DirectoryName = 'H:\\Python\\Johnson03\\Data'+date.today().strftime('%y-%m$
    #DirectoryName = '~\\w2k\\Python\\Johnson03\\Data'+date.today().strftime('%$
    DirectoryName = getcwd() + '/Data'+date.today().strftime('%y-%m-%d') + '/Trial'+str(TrialNumber)+'/'
    if not isdir(DirectoryName):
        makedirs(DirectoryName)

    if isinstance(Pop, list):
        PopName = isum(Pop)
    else:
        PopName = Pop

    for i in range(1000):
        if(Pop):
            List[i] = Model6(Population = Pop, EndTime = Time)
        else:
            List[i] = Model6(EndTime = Time)
    FileName = 'ImperfectFrag2_'+ str(TrialNumber)+'_'+str(PopName)
    savemat(DirectoryName + FileName, mdict = {FileName: List}, do_compression = True)
    print(FileName + ' is Done!')

    for i in range(1000):
        if(Pop):
            List[i] = Model7(Population = Pop, EndTime = Time)
        else:
            List[i] = Model7(EndTime = Time)
    FileName = 'ImperfectAll2_'+ str(TrialNumber)+'_'+str(PopName)
    savemat(DirectoryName + FileName, mdict = {FileName: List}, do_compression = True)
    print(FileName + ' is Done!')


def RunTrialsOfImpFrag3Models(TrialNumber, Pop = 0, Time = 10000):
    from JohnsonOnePopModels import OnePopImperfect3 as Model6
    from JohnsonOnePopModels import OnePopImperfectFrag3 as Model7
    List =    [[]]*1000
    
    from counting import isum
    from os import makedirs
    from os import getcwd
    from os.path import isdir
    from datetime import date
    from scipy.io import savemat
    #C:\User Files\Local Data
    #H:\Python\Johnson03
    #DirectoryName = 'C:\\User Files\\Local Data\\Data'+date.today().strftime('%y-%m-%d')+'\\Trial'+str(TrialNumber)+'\\'
    #DirectoryName = 'H:\\Python\\Johnson03\\Data'+date.today().strftime('%y-%m-%d')+'\\Trial'+str(TrialNumber)+'\\'
    #DirectoryName = '~\\w2k\\Python\\Johnson03\\Data'+date.today().strftime('%y-%m-%d')+'\\Trial'+str(TrialNumber)+'\\'
    DirectoryName = getcwd() + '/Data'+date.today().strftime('%y-%m-%d') + '/Trial'+str(TrialNumber)+'/'
    if not isdir(DirectoryName):
        makedirs(DirectoryName)
    
    if isinstance(Pop, list):
        PopName = isum(Pop)
    else:
        PopName = Pop

    for i in range(1000):
        if(Pop):
            List[i] = Model6(Population = Pop, EndTime = Time)
        else:
            List[i] = Model6(EndTime = Time)
        
    FileName = 'ImperfectAll3_'+ str(TrialNumber)+'_'+str(PopName)
    savemat(DirectoryName + FileName, mdict = {FileName: List}, do_compression = True)
    print(FileName + ' is Done!')
    
    for i in range(1000):
        if(Pop):
            List[i] = Model7(Population = Pop, EndTime = Time)
        else:
            List[i] = Model7(EndTime = Time)
        
    FileName = 'ImperfectFrag3_'+ str(TrialNumber)+'_'+str(PopName)
    savemat(DirectoryName + FileName, mdict = {FileName: List}, do_compression = True)
    print(FileName + ' is Done!')
    
def RunTrialsOfImpFrag4Models(TrialNumber, Pop = 0, Time = 10000):
    from JohnsonOnePopModels import OnePopImperfect4 as Model6
    from JohnsonOnePopModels import OnePopImperfectFrag4 as Model7
    List =    [[]]*1000
    
    from counting import isum
    from os import makedirs
    from os import getcwd
    from os.path import isdir
    from datetime import date
    from scipy.io import savemat
    #C:\User Files\Local Data
    #H:\Python\Johnson03
    #DirectoryName = 'C:\\User Files\\Local Data\\Data'+date.today().strftime('%y-%m-%d')+'\\Trial'+str(TrialNumber)+'\\'
    #DirectoryName = 'H:\\Python\\Johnson03\\Data'+date.today().strftime('%y-%m-%d')+'\\Trial'+str(TrialNumber)+'\\'
    #DirectoryName = '~\\w2k\\Python\\Johnson03\\Data'+date.today().strftime('%y-%m-%d')+'\\Trial'+str(TrialNumber)+'\\'
    DirectoryName = getcwd() + '/Data'+date.today().strftime('%y-%m-%d') + '/Trial'+str(TrialNumber)+'/'
    if not isdir(DirectoryName):
        makedirs(DirectoryName)
    if isinstance(Pop, list):
        PopName = isum(Pop)
    else:
        PopName = Pop
    for i in range(1000):
        if(Pop):
            List[i] = Model6(Population = Pop, EndTime = Time)
        else:
            List[i] = Model6(EndTime = Time)
        
    FileName = 'ImperfectAll4_'+ str(TrialNumber)+'_'+str(PopName)
    savemat(DirectoryName + FileName, mdict = {FileName: List}, do_compression = True)
    print(FileName + ' is Done!')
    
    for i in range(1000):
        if(Pop):
            List[i] = Model7(Population = Pop, EndTime = Time)
        else:
            List[i] = Model7(EndTime = Time)
        
    FileName = 'ImperfectFrag4_'+ str(TrialNumber)+'_'+str(PopName)
    savemat(DirectoryName + FileName, mdict = {FileName: List}, do_compression = True)
    print(FileName + ' is Done!')
    
def RunTrialsOfImpFrag5Models(TrialNumber, Pop = 0, Time = 10000):
    from JohnsonOnePopModels import OnePopImperfect5 as Model6
    from JohnsonOnePopModels import OnePopImperfectFrag5 as Model7
    List =    [[]]*1000
    
    from counting import isum
    from os import makedirs
    from os import getcwd
    from os.path import isdir
    from datetime import date
    from scipy.io import savemat
    #C:\User Files\Local Data
    #H:\Python\Johnson03
    #DirectoryName = 'C:\\User Files\\Local Data\\Data'+date.today().strftime('%y-%m-%d')+'\\Trial'+str(TrialNumber)+'\\'
    #DirectoryName = 'H:\\Python\\Johnson03\\Data'+date.today().strftime('%y-%m-%d')+'\\Trial'+str(TrialNumber)+'\\'
    #DirectoryName = '~\\w2k\\Python\\Johnson03\\Data'+date.today().strftime('%y-%m-%d')+'\\Trial'+str(TrialNumber)+'\\'
    DirectoryName = getcwd() + '/Data'+date.today().strftime('%y-%m-%d') + '/Trial'+str(TrialNumber)+'/'
    if not isdir(DirectoryName):
        makedirs(DirectoryName)
    if isinstance(Pop, list):
        PopName = isum(Pop)
    else:
        PopName = Pop
    
    for i in range(1000):
        if(Pop):
            List[i] = Model6(Population = Pop, EndTime = Time)
        else:
            List[i] = Model6(EndTime = Time)
        
    FileName = 'ImperfectAll5_'+ str(TrialNumber)+'_'+str(PopName)
    savemat(DirectoryName + FileName, mdict = {FileName: List}, do_compression = True)
    print(FileName + ' is Done!')
    
    for i in range(1000):
        if(Pop):
            List[i] = Model7(Population = Pop, EndTime = Time)
        else:
            List[i] = Model7(EndTime = Time)
        
    FileName = 'ImperfectFrag5_'+ str(TrialNumber)+'_'+str(PopName)
    savemat(DirectoryName + FileName, mdict = {FileName: List}, do_compression = True)
    print(FileName + ' is Done!')
    
def RunTrialsOfImpFrag6Models(TrialNumber, Pop = 0, Time = 10000):
    from JohnsonOnePopModels import OnePopImperfect6 as Model6
    from JohnsonOnePopModels import OnePopImperfectFrag6 as Model7
    List =    [[]]*1000
    
    from counting import isum
    from os import makedirs
    from os import getcwd
    from os.path import isdir
    from datetime import date
    from scipy.io import savemat
    #C:\User Files\Local Data
    #H:\Python\Johnson03
    #DirectoryName = 'C:\\User Files\\Local Data\\Data'+date.today().strftime('%y-%m-%d')+'\\Trial'+str(TrialNumber)+'\\'
    #DirectoryName = 'H:\\Python\\Johnson03\\Data'+date.today().strftime('%y-%m-%d')+'\\Trial'+str(TrialNumber)+'\\'
    #DirectoryName = '~\\w2k\\Python\\Johnson03\\Data'+date.today().strftime('%y-%m-%d')+'\\Trial'+str(TrialNumber)+'\\'
    DirectoryName = getcwd() + '/Data'+date.today().strftime('%y-%m-%d') + '/Trial'+str(TrialNumber)+'/'
    if not isdir(DirectoryName):
        makedirs(DirectoryName)
    if isinstance(Pop, list):
        PopName = isum(Pop)
    else:
        PopName = Pop
    
    for i in range(1000):
        if(Pop):
            List[i] = Model6(Population = Pop, EndTime = Time)
        else:
            List[i] = Model6(EndTime = Time)
        
    FileName = 'ImperfectAll6_'+ str(TrialNumber)+'_'+str(PopName)
    savemat(DirectoryName + FileName, mdict = {FileName: List}, do_compression = True)
    print(FileName + ' is Done!')
    
    for i in range(1000):
        if(Pop):
            List[i] = Model7(Population = Pop, EndTime = Time)
        else:
            List[i] = Model7(EndTime = Time)
        
    FileName = 'ImperfectFrag6_'+ str(TrialNumber)+'_'+str(PopName)
    savemat(DirectoryName + FileName, mdict = {FileName: List}, do_compression = True)
    print(FileName + ' is Done!')

def RunTrialsOfJ1PopImperfects(NumberOfRuns, DirectoryNumber):
    #from JohnsonOnePopModels import OnePopImperfect as Model6
    from JohnsonOnePopModels import OnePopImperfectFrag2 as Model7
    #from JohnsonOnePopModels import OnePopImperfectComb as Model8
    #from JohnsonOnePopModels import OnePopMultiStep as Model9
    #Storage List
    List =    [[]]*1000
    
    from os import makedirs
    from os.path import isdir
    from datetime import date
    from scipy.io import savemat
    
    #C:\User Files\Local Data
    #H:\Python\Johnson03
    #I:\research\combatmodelling\Johnson One Population Model\One Population Model Data
    #DirectoryName = 'C:\\User Files\\Local Data\\Data'+date.today().strftime('%y-%m-%d')+'\\Trial'+str(DirectoryNumber)+'\\'
    #DirectoryName = 'H:\\Python\\Johnson03\\Data'+date.today().strftime('%y-%m-%d')+'\\Trial'+str(DirectoryNumber)+'\\'
    DirectoryName = 'I:\\research\\combatmodelling\\Johnson One Population Model\\One Population Model Data\\Data'+date.today().strftime('%y-%m-%d')+'\\Trial'+str(DirectoryNumber)+'\\'
    if not isdir(DirectoryName):
        makedirs(DirectoryName)
    
#    for RunNumber in range(NumberOfRuns):
#        for i in range(1000):
#            List[i] = Model6()
#        
#        FileName = 'ImperfectAll'+ str(RunNumber)
#        savemat(DirectoryName + FileName, mdict = {FileName: List})
#        print(FileName + ' is Done!')
        
    for RunNumber in range(NumberOfRuns):
        for i in range(1000):
            List[i] = Model7()
        
        FileName = 'ImperfectFrag2'+ str(RunNumber)
        savemat(DirectoryName + FileName, mdict = {FileName: List}, do_compression = True)
        print(FileName + ' is Done!')
        
#    for RunNumber in range(NumberOfRuns):
#        for i in range(1000):
#            List[i] = Model8()
#        
#        FileName = 'ImperfectComb'+ str(RunNumber)
#        savemat(DirectoryName + FileName, mdict = {FileName: List})
#        print(FileName + ' is Done!')
#        
#    for RunNumber in range(NumberOfRuns):
#        for i in range(1000):
#            List[i] = Model9()
#        
#        FileName = 'MultiSteps'+ str(RunNumber)
#        savemat(DirectoryName + FileName, mdict = {FileName: List})
#        print(FileName + ' is Done!')
        
            
def RunTrialsOfJ1PopAcAts(NumberOfRuns, DirectoryNumber, AcAts):
#    from JohnsonOnePopModels import OnePopWithAttrition as Model3
    from JohnsonOnePopModels import OnePopWithAccretion as Model4
    from JohnsonOnePopModels import OnePopWithAccrAttr as Model5
    #Storage List
    List =    [[]]*1000
    
    from os import makedirs
    from os.path import isdir
    from datetime import date
    from scipy.io import savemat
    
    #C:\User Files\Local Data
    #H:\Python\Johnson03
    #I:\research\combatmodelling\Johnson One Population Model\One Population Model Data
    #DirectoryName = 'C:\\User Files\\Local Data\\Data'+date.today().strftime('%y-%m-%d')+'\\Trial'+str(DirectoryNumber)+'\\'
    #DirectoryName = 'H:\\Python\\Johnson03\\Data'+date.today().strftime('%y-%m-%d')+'\\Trial'+str(DirectoryNumber)+'\\'
    DirectoryName = 'I:\\research\\combatmodelling\\Johnson One Population Model\\One Population Model Data\\Data'+date.today().strftime('%y-%m-%d')+'\\Trial'+str(DirectoryNumber)+'\\'
    if not isdir(DirectoryName):
        makedirs(DirectoryName)
    
#    for RunNumber in range(NumberOfRuns):
#        for i in range(1000):
#            List[i] = Model3(Attrition = AcAts)
#        
#        FileName = 'Attrition'+ str(RunNumber)
#        savemat(DirectoryName + FileName, mdict = {FileName: List})
#        print(FileName + ' is Done!')
        
    for RunNumber in range(NumberOfRuns):
        for i in range(1000):
            List[i] = Model4(Accretion = AcAts)
        
        FileName = 'Accretion'+ str(RunNumber)
        savemat(DirectoryName + FileName, mdict = {FileName: List})
        print(FileName + ' is Done!')
        
    for RunNumber in range(NumberOfRuns):
        for i in range(1000):
            List[i] = Model5(Accretion = AcAts, Attrition = AcAts)
        
        FileName = 'AccrAttr'+ str(RunNumber)
        savemat(DirectoryName + FileName, mdict = {FileName: List})
        print(FileName + ' is Done!')
        
def RunTrialsOfJ1PopAcAt(NumberOfRuns, DirectoryNumber, AcAts):
#    from JohnsonOnePopModels import OnePopWithAttrition as Model3
#    from JohnsonOnePopModels import OnePopWithAccretion as Model4
    from JohnsonOnePopModels import OnePopWithAccrAttr as Model5
    #Storage List
    List =    [[]]*1000
    
    from os import makedirs
    from os.path import isdir
    from datetime import date
    from scipy.io import savemat
    
    #C:\User Files\Local Data
    #H:\Python\Johnson03
    #I:\research\combatmodelling\Johnson One Population Model\One Population Model Data
    #DirectoryName = 'C:\\User Files\\Local Data\\Data'+date.today().strftime('%y-%m-%d')+'\\Trial'+str(DirectoryNumber)+'\\'
    #DirectoryName = 'H:\\Python\\Johnson03\\Data'+date.today().strftime('%y-%m-%d')+'\\Trial'+str(DirectoryNumber)+'\\'
    DirectoryName = 'I:\\research\\combatmodelling\\Johnson One Population Model\\One Population Model Data\\Data'+date.today().strftime('%y-%m-%d')+'\\Trial'+str(DirectoryNumber)+'\\'
    if not isdir(DirectoryName):
        makedirs(DirectoryName)
    
#    for RunNumber in range(NumberOfRuns):
#        for i in range(1000):
#            List[i] = Model3(Attrition = AcAts)
#        
#        FileName = 'Attrition'+ str(RunNumber)
#        savemat(DirectoryName + FileName, mdict = {FileName: List})
#        print(FileName + ' is Done!')
        
#    for RunNumber in range(NumberOfRuns):
#        for i in range(1000):
#            List[i] = Model4(Accretion = AcAts)
#        
#        FileName = 'Accretion'+ str(RunNumber)
#        savemat(DirectoryName + FileName, mdict = {FileName: List})
#        print(FileName + ' is Done!')
        
    for RunNumber in range(NumberOfRuns):
        for i in range(1000):
            List[i] = Model5(Accretion = AcAts, Attrition = AcAts)
        
        FileName = 'AccrAttr'+ str(RunNumber)
        savemat(DirectoryName + FileName, mdict = {FileName: List})
        print(FileName + ' is Done!')

def RunTrialsOfJ1PopAc(NumberOfRuns, DirectoryNumber, AcAts):
#    from JohnsonOnePopModels import OnePopWithAttrition as Model3
    from JohnsonOnePopModels import OnePopWithAccretion as Model4
#    from JohnsonOnePopModels import OnePopWithAccrAttr as Model5
    #Storage List
    List =    [[]]*1000
    
    from os import makedirs
    from os.path import isdir
    from datetime import date
    from scipy.io import savemat
    
    #C:\User Files\Local Data
    #H:\Python\Johnson03
    #I:\research\combatmodelling\Johnson One Population Model\One Population Model Data
    #DirectoryName = 'C:\\User Files\\Local Data\\Data'+date.today().strftime('%y-%m-%d')+'\\Trial'+str(DirectoryNumber)+'\\'
    #DirectoryName = 'H:\\Python\\Johnson03\\Data'+date.today().strftime('%y-%m-%d')+'\\Trial'+str(DirectoryNumber)+'\\'
    DirectoryName = 'I:\\research\\combatmodelling\\Johnson One Population Model\\One Population Model Data\\Data'+date.today().strftime('%y-%m-%d')+'\\Trial'+str(DirectoryNumber)+'\\'
    if not isdir(DirectoryName):
        makedirs(DirectoryName)
    
#    for RunNumber in range(NumberOfRuns):
#        for i in range(1000):
#            List[i] = Model3(Attrition = AcAts)
#        
#        FileName = 'Attrition'+ str(RunNumber)
#        savemat(DirectoryName + FileName, mdict = {FileName: List})
#        print(FileName + ' is Done!')
        
    for RunNumber in range(NumberOfRuns):
        for i in range(1000):
            List[i] = Model4(Accretion = AcAts)
        
        FileName = 'Accretion'+ str(RunNumber)
        savemat(DirectoryName + FileName, mdict = {FileName: List})
        print(FileName + ' is Done!')
        
#    for RunNumber in range(NumberOfRuns):
#        for i in range(1000):
#            List[i] = Model5(Accretion = AcAts, Attrition = AcAts)
#        
#        FileName = 'AccrAttr'+ str(RunNumber)
#        savemat(DirectoryName + FileName, mdict = {FileName: List})
#        print(FileName + ' is Done!')

def RunTrialsOfJ1PopUnstable(NumberOfRuns, DirectoryNumber):
    from JohnsonOnePopModels import UnstablePop as Model10
    #Storage List
    List =    [[]]*1000
    
    from os import makedirs
    from os.path import isdir
    from datetime import date
    from scipy.io import savemat
    
    #C:\User Files\Local Data
    #H:\Python\Johnson03
    #I:\research\combatmodelling\Johnson One Population Model\One Population Model Data
    #DirectoryName = 'C:\\User Files\\Local Data\\Data'+date.today().strftime('%y-%m-%d')+'\\Trial'+str(DirectoryNumber)+'\\'
    #DirectoryName = 'H:\\Python\\Johnson03\\Data'+date.today().strftime('%y-%m-%d')+'\\Trial'+str(DirectoryNumber)+'\\'
    DirectoryName = 'I:\\research\\combatmodelling\\Johnson One Population Model\\One Population Model Data\\Data'+date.today().strftime('%y-%m-%d')+'\\Trial'+str(DirectoryNumber)+'\\'
    if not isdir(DirectoryName):
        makedirs(DirectoryName)
    
    for RunNumber in range(NumberOfRuns):
        for i in range(1000):
            List[i] = Model10()
        FileName = 'UnstablePop'+ str(RunNumber)
        savemat(DirectoryName + FileName, mdict = {FileName: List})
        print(FileName + ' is Done!')
        
def RunTrialsOfJ1PopMultiStep(NumberOfRuns, DirectoryNumber):
    from JohnsonOnePopModels import OnePopMultiStep as Model9
    #Storage List
    List =    [[]]*1000
    
    from os import makedirs
    from os.path import isdir
    from datetime import date
    from scipy.io import savemat
    
    #C:\User Files\Local Data
    #H:\Python\Johnson03
    #I:\research\combatmodelling\Johnson One Population Model\One Population Model Data
    #DirectoryName = 'C:\\User Files\\Local Data\\Data'+date.today().strftime('%y-%m-%d')+'\\Trial'+str(DirectoryNumber)+'\\'
    #DirectoryName = 'H:\\Python\\Johnson03\\Data'+date.today().strftime('%y-%m-%d')+'\\Trial'+str(DirectoryNumber)+'\\'
    DirectoryName = 'I:\\research\\combatmodelling\\Johnson One Population Model\\One Population Model Data\\Data'+date.today().strftime('%y-%m-%d')+'\\Trial'+str(DirectoryNumber)+'\\'
    if not isdir(DirectoryName):
        makedirs(DirectoryName)

    for RunNumber in range(NumberOfRuns):
        for i in range(1000):
            List[i] = Model9(Population = 10000, EndTime = 100000)
        
        FileName = 'MultiSteps'+ str(RunNumber)
        savemat(DirectoryName + FileName, mdict = {FileName: List})
        print(FileName + ' is Done!')
        
def RunTrialsOfJ1PopMust(NumberOfRuns, DirectoryNumber, Pop = 0, Time = 10000):
    from JohnsonOnePopModels import JohnsonOnePopMustSucceed as Model1
    #Storage List
    List =    [[]]*1000
    
    from os import makedirs
    from os.path import isdir
    from datetime import date
    from scipy.io import savemat
    
    #C:\User Files\Local Data
    #H:\Python\Johnson03
    #I:\research\combatmodelling\Johnson One Population Model\One Population Model Data
    DirectoryName = 'C:\\User Files\\Local Data\\Data'+date.today().strftime('%y-%m-%d')+'\\Trial'+str(DirectoryNumber)+'\\'
    #DirectoryName = 'H:\\Python\\Johnson03\\Data'+date.today().strftime('%y-%m-%d')+'\\Trial'+str(DirectoryNumber)+'\\'
    #DirectoryName = 'I:\\research\\combatmodelling\\Johnson One Population Model\\One Population Model Data\\Data'+date.today().strftime('%y-%m-%d')+'\\Trial'+str(DirectoryNumber)+'\\'
    if not isdir(DirectoryName):
        makedirs(DirectoryName)

    for RunNumber in range(NumberOfRuns):
        for i in range(1000):
            if(Pop):
                List[i] = Model1(Pop, Time)
            else:
                List[i] = Model1(Time)
                
        
        FileName = 'Must'+ str(RunNumber)
        savemat(DirectoryName + FileName, mdict = {FileName: List}, do_compression = True)
        print(FileName + ' is Done!')
        
def HistoryTrialsJ1Must(PopVector, TimeSteps,DirectoryNumber,  HistSteps = 1):
    from JohnsonOnePopModels import JohnsonOnePopMustSucceed as Model1
    from os import makedirs
    from os.path import isdir
    from datetime import date
    from scipy.io import savemat
    
    DirectoryName = 'I:\\research\\combatmodelling\\Johnson One Population Model\\One Population Model Data\\Data'+date.today().strftime('%y-%m-%d')+'\\Trial'+str(DirectoryNumber)+'\\'
    if not isdir(DirectoryName):
        makedirs(DirectoryName)
        
    ReportedHistory = Model1(Population = PopVector, EndTime = TimeSteps, History = 1, HistorySteps = HistSteps)
    
    FileName = 'MustHistory'+str(TimeSteps)
    savemat(DirectoryName + FileName, mdict = {FileName: ReportedHistory}, do_compression = True)
    print(FileName + ' is Done!')
