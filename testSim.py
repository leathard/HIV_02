import time,os
import initSim
from Sim import onesim,get_param_set,buildPopulation
from loadConfig import *
from Persons import *
from Scheduler import *

import unittest
from unittest import suite
import sys
import importlib
import multiprocessing as mp
import shutil
import glob

class TestLoaderWithKwargs(unittest.TestLoader):
    """A test loader which allows to parse keyword arguments to the
       test case class."""
    def loadTestsFromTestCase(self, testCaseClass, **kwargs):
        """Return a suite of all tests cases contained in 
           testCaseClass."""
        if issubclass(testCaseClass, suite.TestSuite):
            raise TypeError("Inherit unittest.TestCase for test class")
        testCaseNames = self.getTestCaseNames(testCaseClass)
        if not testCaseNames and hasattr(testCaseClass, 'runTest'):
            testCaseNames = ['runTest']

        # Modification here: parse keyword arguments to testCaseClass.
        test_cases = []
        for test_case_name in testCaseNames:
            test_cases.append(testCaseClass(test_case_name, **kwargs))
        loaded_suite = self.suiteClass(test_cases)

        return loaded_suite 

def insertDefaultParameters(params):
    model_params={}
    model_name = 'test_model'
    outfolder = 'temp01'
    model_params['outfile_pattern'] = 'eps??sim??.out'
    model_params['nSim'] = 10 
    model_params['nMF'] = 500,500
    model_params['nSimDays'] = 365 * 3
    model_params['nBurnDays'] = 80
    model_params['nOutputInterval'] = 20 
    model_params['outfolder'] = outfolder
    model_params['model_name'] = 'experiment'

    for k,v in model_params.items():
        params[k]=params.get(k,model_params[k])
    return params



class testSimParts(unittest.TestCase):
    def __init__(self, testName, params=None):
        super(testSimParts, self).__init__(testName)
        #inserting default dummy parameters
        self.params = insertDefaultParameters(params.copy())
        
    def setUp(self):
        print "building population for test"
        for params in get_param_set(self.params):
            break
        self.params=params
        self.males,self.females=buildPopulation(params)
        
    def tearDown(self):
        pass

    def test_buildPopulation(self):
        print "\nPopulation composed of: %d males, %d females"%(len(self.males),len(self.females))        
        self.assertTrue(abs(len(self.males)-len(self.females))<2)
        #computing expected ratios
        n_sexworkers=params['p_nsexworkers']*len(self.females)
        nM,nF=params['nMF']
        pop_size=nM+nF
        n_prep=params['p_PREP']*pop_size
        n_f_art=params['p_nM_ART']*len(self.males)
        n_m_art=params['p_nM_ART']*len(self.females)
        n_art=n_f_art+n_m_art
        #counting ratios achieved
        nsw=len(filter(lambda x:x.comsex,self.females))
        nprep=len(filter(lambda x:x.PrEP,self.females+self.males))
        nart=len(filter(lambda x:x.ART,self.females+self.males))
        
        print "\t sexworkers found: %d, expected: %f error:%f"%(nsw,n_sexworkers,float(nsw-n_sexworkers)/float(n_sexworkers))
        print "\t art found: %d, expected: %f error:%f"%(nart,n_art,float(nart-n_art)/float(n_art))
        print "\t prep found: %d, expected: %f error:%f"%(nprep,n_prep,float(nprep-n_prep)/float(n_prep))
        self.assertTrue(abs(n_sexworkers-nsw)/float(n_sexworkers)<0.05)
        self.assertTrue(abs(n_f_art+n_m_art-nart)/float(n_f_art+n_m_art)<0.05)
        self.assertTrue(abs(n_prep-nprep)/float(n_prep)<0.05)

    def test_formAndEvolvePartnerships(self):
        # use the previously build population to evolve it a little and check the results of the evolution
        schedule = Scheduler(params=self.params) 
        for m in self.males: schedule.register_person(m)
        for f in self.females: schedule.register_person(f)
        print "forming partnerships"
        for day in xrange(params['nBurnDays']):
            schedule.coresim(males=self.males,females=self.females,day=day,params=params)
        print "n_partnerships:%d"%schedule.n_partnerships
        self.assertTrue(schedule.n_partnerships>1)

        print "seeding infections"
        diseases = seed_infections(self.males, self.females, params['nBurnDays'], schedule=schedule, params=self.params)
        initial_f_infections=sum(f.is_infected for f in self.females)
        initial_m_infections=sum(m.is_infected for m in self.males)
        print "evolving partnerships"
        n_transmissions_max=0
        for day in xrange(params['nBurnDays'],2*params['nBurnDays']):
            schedule.coresim(males=self.males,females=self.females,day=day,params=self.params)            
            n_transmissions_max=max(schedule.count_transmissions(),n_transmissions_max)
        self.nTransmissions=n_transmissions_max            
        self.assertTrue(n_transmissions_max>1)
        
        final_f_infections=sum(f.is_infected for f in self.females)
        final_m_infections=sum(m.is_infected for m in self.males)
        
        self.assertTrue(sum(f.is_infected for f in self.females) > 0)
        self.assertTrue(sum(m.is_infected for m in self.males) > 0 )
        print "initial infections:%d,%d"%(initial_f_infections,initial_m_infections)
        print "final infections:%d,%d"%(final_f_infections,final_m_infections)
        # not strictly true, but highly unlikely
        # cannot guarantee final  will be higher than initial
        # nor the reverse
        self.assertTrue(initial_f_infections!=final_f_infections)
        self.assertTrue(initial_m_infections!=final_m_infections)
        # more criteria could be inserted here 

def load_file(fname):
    "loader for data files obtained during simulation"
    data=np.genfromtxt(fname,dtype=None,delimiter=',',names=True)
    return data


class testSim(unittest.TestCase):
    def __init__(self, testName, params=None):
        super(testSim, self).__init__(testName)
        self.testPassed=False
        try:
            shutil.rmtree("test")
        except OSError,e:
            #doesn't exist, nothing to clean up
            pass
        self.params = insertDefaultParameters(params.copy())
        #override the length of the simulation
        self.params['nSimDays']=self.params['nBurnDays'] 

        #make up alternative realities to compare the simulation results in the end
        neutralParams = self.params.copy()
        neutralParams.update({'p_nF_ART':0.0,'p_nM_ART':0.0,'p_PREP':0.0})
        neutralParams['outfolder']="test/neutral_"+neutralParams['outfolder']
        neutralParams['sScenarioName']=neutralParams['outfolder']
        
        infectParams = self.params.copy()
        infectParams.update({'p_nF_ART':0.0,'p_nM_ART':0.0,'p_PREP':0.0,'pSeedHIV':map(lambda x:2*x,neutralParams['pSeedHIV'])})
        infectParams['outfolder']="test/infected_"+infectParams['outfolder']
        infectParams['sScenarioName']=infectParams['outfolder']
        
        artParams = self.params.copy()
        artParams.update({'p_PREP':0.0})
        artParams['outfolder']="test/art_"+artParams['outfolder']
        artParams['sScenarioName']=artParams['outfolder']
        
        prepParams = self.params.copy()
        prepParams.update({'p_ART':0.0})        
        prepParams['outfolder']="test/prep_"+ prepParams['outfolder']
        prepParams['sScenarioName']=prepParams['outfolder']
                
        noComWorkersParams=self.params.copy()
        noComWorkersParams.update({'p_nclients': 0.0,'p_nsexworkers': 0.0})
        noComWorkersParams['outfolder']="test/nocsw_"+noComWorkersParams['outfolder']
        noComWorkersParams['sScenarioName']=noComWorkersParams['outfolder']
                
        self.paramDict={'neutral':neutralParams,'art':artParams,'prep':prepParams,
                            'noCSW':noComWorkersParams,'infected':infectParams}
        # upgrade parameter set with the simulation parameters
        for runType in self.paramDict.keys():
            for params in get_param_set(self.paramDict[runType]):
                self.paramDict[runType]=params.copy()
                break
        
    def runSimulations(self):
        if os.path.exists('test'): return
        print "runSimulations"
        for p in self.paramDict.values():
            initSim.prepare(p['outfolder'], logPath=p['outfolder']+'_out/temp.log')
        #run the alternative simulations
        n_processors = mp.cpu_count()        
        workerpool = mp.Pool(n_processors)
        workerpool.map(onesim, self.paramDict.values())
        #for p in self.paramDict.values():
        #    onesim(p)
        workerpool.close()
        workerpool.join()
        
    def setUp(self):
        self.runSimulations()
            
    def testCompareSimulationsResults(self):
        self.results={}
        for runType,params in self.paramDict.items():
            datalist=glob.glob(params['outfolder']+"/"+params['outfile_pattern'])
            # expecting only one
            self.results[runType]={'file':datalist[0]}
            data=load_file(self.results[runType]['file'])
            self.results[runType]['nInfect']=data['nMinfect']+data['nFinfect']
        # expect the last known data as obtained by evolving simulation to obey some conditions
        # no treatment is better than treatment, no sex workers is better than with sexworkers
        # higher initial infection ratio will propagate through the end
        # these conditions should hold at least on average
        print
        print "neutral \t art \t prep \t noCSW \t infected"
        print "%4.1f \t \t %4.1f \t %4.1f \t %4.1f \t %4.1f"%\
          (np.mean(self.results['neutral']['nInfect']),np.mean(self.results['art']['nInfect']),
                                                np.mean(self.results['prep']['nInfect']),np.mean(self.results['noCSW']['nInfect']),
                                                np.mean(self.results['infected']['nInfect']))
        self.assertTrue(np.mean(self.results['neutral']['nInfect'])>np.mean(self.results['art']['nInfect']))
        self.assertTrue(np.mean(self.results['neutral']['nInfect'])>np.mean(self.results['prep']['nInfect']))
        self.assertTrue(np.mean(self.results['neutral']['nInfect'])>np.mean(self.results['noCSW']['nInfect']))
        self.assertTrue(np.mean(self.results['neutral']['nInfect'])<np.mean(self.results['infected']['nInfect']))
        self.testPassed=True

    def testCFSW(self):
        self.results={}
        for runType,params in self.paramDict.items():
            datalist=glob.glob(params['outfolder']+"/"+params['outfile_pattern'])
            self.results[runType]={'file':datalist[0]}
            data=load_file(self.results[runType]['file'])
            self.results[runType]['nFinfect']=data['nFinfect']
            self.results[runType]['nMinfect']=data['nMinfect']
            self.results[runType]['nFSWinf']=data['nFSWinf']
            self.results[runType]['nCSWinf']=data['nCSWinf']
            
        nFSWinfectEstimate=self.results['neutral']['nFinfect']*self.paramDict['neutral']['p_nsexworkers']
        nCSWinfectEstimate=self.results['neutral']['nMinfect']*self.paramDict['neutral']['p_nclients']
        nFSWinfectSim=self.results['neutral']['nFSWinf']
        nCSWinfectSim=self.results['neutral']['nCSWinf']
        print
        print "FSW infected simulation:%r,\testimate:%r"%(nFSWinfectSim,nFSWinfectEstimate)
        print "CSW infected simulation:%r,\testimate:%r"%(nCSWinfectSim,nCSWinfectEstimate)
        
        self.assertTrue(np.mean(abs(nFSWinfectSim-nFSWinfectEstimate)/nFSWinfectSim)<0.3)        
        self.assertTrue(np.mean(abs(nCSWinfectSim-nCSWinfectEstimate)/nCSWinfectSim)<0.3)        
        self.testPassed=True

    def testMiners(self):
        self.testPassed = False # won't clean up logs if test failed
        self.results={}
        for runType,params in self.paramDict.items():
            datalist=glob.glob(params['outfolder']+"/"+params['outfile_pattern'])
            self.results[runType]={'file':datalist[0]}
            data=load_file(self.results[runType]['file'])
            self.results[runType]['nMinfect']=data['nMinfect']
            self.results[runType]['nMinerInf']=data['nMinerInf']
        nMinersInfEstimate=self.results['neutral']['nMinfect']*self.paramDict['neutral']['p_miners'] 
        nMinersInfSim=self.results['neutral']['nMinerInf']
        print
        print "Miners infected simulation:%r,\testimate:%r"%( nMinersInfSim,nMinersInfEstimate)    
        self.assertTrue(np.mean(abs(nMinersInfSim-nMinersInfEstimate)/nMinersInfSim)<0.4)
        self.testPassed=True
        
    def testARTPREP(self):
        self.results={}
        for runType,params in self.paramDict.items():
            datalist=glob.glob(params['outfolder']+"/"+params['outfile_pattern'])
            self.results[runType]={'file':datalist[0]}
            data=load_file(self.results[runType]['file'])
            self.results[runType]['nFinfect']=data['nFinfect']
            self.results[runType]['nMinfect']=data['nMinfect']
            self.results[runType]['nARTinf']=data['nARTinf']
            self.results[runType]['nPREPinf']=data['nPREPinf']

        nArtInfectEstimate = self.results['art']['nFinfect']*self.paramDict['art']['p_nF_ART']+\
          self.results['art']['nMinfect']*self.paramDict['art']['p_nM_ART']
        nArtInfectSim      = self.results['art']['nARTinf']
        
        nInfect=self.results['prep']['nFinfect']+self.results['prep']['nMinfect']
        nPREPInfectEstimate = nInfect*self.paramDict['prep']['p_PREP']
        nPREPInfectSim      = self.results['prep']['nPREPinf']
        print
        print "ART infected simulation:%r,\testimate:%r"%(nArtInfectSim,nArtInfectEstimate)
        print "PREP infected simulation:%r,\testimate:%r"%(nPREPInfectSim,nPREPInfectEstimate)
        
        self.assertTrue(np.mean(abs(nArtInfectSim-nArtInfectEstimate)/nArtInfectSim)<0.4) 
        self.assertTrue(np.mean(abs(nPREPInfectSim-nPREPInfectEstimate)/nPREPInfectSim)<0.3) 
        self.testPassed=True
        
    def testCFSWMiners_ART(self):
        self.results={}
        for runType,params in self.paramDict.items():
            datalist=glob.glob(params['outfolder']+"/"+params['outfile_pattern'])
            self.results[runType]={'file':datalist[0]}
            data=load_file(self.results[runType]['file'])
            self.results[runType]['nFinfect']=data['nFinfect']
            self.results[runType]['nMinfect']=data['nMinfect']
            self.results[runType]['nFSWARTinf']=data['nFSWARTinf']
            self.results[runType]['nCSWARTinf']=data['nCSWARTinf']            
            self.results[runType]['nMinARTinf']=data['nMinARTinf']
        
        nFSWARTInfectEstimate   = self.results['art']['nFinfect']*self.paramDict['art']['p_nsexworkers']*\
                               self.paramDict['art']['p_nF_ART']
        nFSWARTInfectSim        = self.results['art']['nFSWARTinf']                               
        nCSWARTInfectEstimate   = self.results['art']['nMinfect']*self.paramDict['art']['p_nclients']*\
                               self.paramDict['art']['p_nM_ART']
        nCSWARTInfectSim        = self.results['art']['nCSWARTinf']
        
        nMinerARTInfectEstimate = self.results['art']['nMinfect']*self.paramDict['art']['p_miners']*\
                               self.paramDict['art']['p_nM_ART']
        nMinerARTInfectSim      = self.results['art']['nMinARTinf']
        print
        print "FSW(ART) infected simulation:%r,\testimate:%r"%(nFSWARTInfectSim,nFSWARTInfectEstimate)
        print "CSW(ART) infected simulation:%r,\testimate:%r"%(nCSWARTInfectSim ,nCSWARTInfectEstimate)
        print "Miners(ART) infected simulation:%r,\testimate:%r"%(nMinerARTInfectSim,nMinerARTInfectEstimate)

        
        self.assertTrue(any(nFSWARTInfectSim==0) or \
                        np.mean(abs(nFSWARTInfectSim-nFSWARTInfectEstimate)/nFSWARTInfectSim)<0.35) 
        self.assertTrue(any(nCSWARTInfectSim==0) or \
                        np.mean(abs(nCSWARTInfectSim-nCSWARTInfectEstimate)/nCSWARTInfectSim)<0.35) 
        self.assertTrue(any(nMinerARTInfectSim==0) or\
                        np.mean(abs(nMinerARTInfectSim-nMinerARTInfectEstimate)/nMinerARTInfectSim)<0.3)
        self.testPassed=True

        
    def testCFSWMiners_PREP(self):

        self.results={}
        for runType,params in self.paramDict.items():
            datalist=glob.glob(params['outfolder']+"/"+params['outfile_pattern'])
            self.results[runType]={'file':datalist[0]}
            data=load_file(self.results[runType]['file'])
            self.results[runType]['nFinfect']=data['nFinfect']
            self.results[runType]['nMinfect']=data['nMinfect']
            self.results[runType]['nFSWPREPinf']=data['nFSWPREPinf']
            self.results[runType]['nCSWPREPinf']=data['nCSWPREPinf']            
            self.results[runType]['nMinPREPinf']=data['nMinPREPinf']
        
        nFSWPREPInfectEstimate = self.results['prep']['nFinfect']*self.paramDict['prep']['p_nsexworkers']*\
                               self.paramDict['prep']['p_PREP']
        nFSWPREPInfectSim      = self.results['prep']['nFSWPREPinf']
        
        nCSWPREPInfectEstimate=self.results['prep']['nMinfect']*self.paramDict['prep']['p_nclients']*\
                               self.paramDict['prep']['p_PREP']
        nCSWPREPInfectSim     = self.results['prep']['nCSWPREPinf']
        
        nMinerPREPInfectEstimate = self.results['prep']['nMinfect']*self.paramDict['prep']['p_miners']*\
                               self.paramDict['prep']['p_PREP']
        nMinerPREPInfectSim      = self.results['prep']['nMinPREPinf']
        print        
        print "FSW(PREP) infected simulation:%r,\testimate:%r"%(nFSWPREPInfectSim,nFSWPREPInfectEstimate)
        print "CSW(PREP) infected simulation:%r,\testimate:%r"%(nCSWPREPInfectSim,nCSWPREPInfectEstimate)
        print "Miners(PREP) infected simulation:%r,\testimate:%r"%( nMinerPREPInfectSim,nMinerPREPInfectEstimate)

        self.assertTrue(any(nFSWPREPInfectSim==0) or\
                         np.mean(abs(nFSWPREPInfectSim-nFSWPREPInfectEstimate)/nFSWPREPInfectSim)<0.35)
        self.assertTrue(np.mean(abs(nCSWPREPInfectSim-nCSWPREPInfectEstimate)/nCSWPREPInfectSim)<0.7) 
        self.assertTrue(any(nMinerPREPInfectSim==0) or\
                         np.mean(abs(nMinerPREPInfectSim-nMinerPREPInfectEstimate)/nMinerPREPInfectSim)<0.35)        
        self.testPassed=True
        
if __name__=='__main__':
    baseline = 'tests.ini'
    params=loadConfig(baseline)
    params['model_name']='experiment'
    updateComputableParts(params)

    # stand-alone run for interactive debug
    #test=testSimParts("test_buildPopulation",params=params)
    #test.setUp()
    #test.test_buildPopulation()
    #test.test_formAndEvolvePartnerships()
    # call your test
    #testS=testSim("testCompareSimulationsResults",params=params)
    #testS.setUp()
    #testS.testCompareSimulationsResults()
    #testS.testCFSW()
    #testS.testMiners()
    #testS.testARTPREP()
    #testS.testCFSWMiners_ART()
    #testS.testCFSWMiners_PREP()
    
    # test using framework

    loader = TestLoaderWithKwargs()
    suiteParts = loader.loadTestsFromTestCase(testSimParts, params=params)    
    unittest.TextTestRunner(verbosity=2).run(suiteParts)
    suiteSimulation = loader.loadTestsFromTestCase(testSim, params=params)    
    unittest.TextTestRunner(verbosity=2).run(suiteSimulation)
    # cleanup the output directory if all tests passed
    if all(map(lambda x:x.testPassed,suiteSimulation._tests)):
        shutil.rmtree("test")

