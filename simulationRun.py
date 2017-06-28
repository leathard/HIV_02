import time
import sys
import initSim
from Sim import onesim
from loadConfig import *
from runSim import runSimulationsMP,runSimulations



if __name__=='__main__':
    baseline = 'baseline.ini'
    experimentsFile='experiments.ini'
    experiment=None
    rParams=loadConfig(baseline)
    rParams['model_name']='baseline'
    #sys.argv.append('test')
    if len(sys.argv)>1:
        experiment=sys.argv[1]
        xParams=loadConfig(experimentsFile,defsection=experiment)
        rParams.update(xParams)
        rParams['model_name']=experiment
    updateComputableParts(rParams)
    outfolder=rParams['sScenarioName']
    initSim.prepare(outfolder, logPath='out/temp.log')
    
    T0 = time.time()
    runSimulationsMP(onesim, rParams)
    print "simulation took {minutes:.2f} minutes".format(minutes=(time.time() - T0) / 60.0)
    print


