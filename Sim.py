"""Run a single concurrency simulation. 

Provides function ``get_param_set`` needed to build population and run sim. 
Builds population in ``buildPopulation`` from given parameters.
Provides function ``onesim``, which can run a single concurrency simulation.

This code is an interation of Sim.py first created by Alan G. Isaac by adding
sub-populations: CSW, client, miner, ART, PrEP to buildPopulation and onesim. 

"""

from collections import defaultdict
import copy
import gc
import pprint
import shutil
import sys
import tempfile, time
from numpy.random import RandomState
from Infection import stagedHIVfactory
#import Persons
from Persons import *
from Scheduler import *


def get_param_set(params):
    """Get a specific set of parameters.

    Yield specific set of parameters as dictionary.
    There are n_sim (100) replicates for each epsilon.
    Each epsilon has its own random seed, sim_name, and outfilename.

    This generator function should yield dicts with **pickleable**
    elements only, so that it can be used with multiprocessing.

    Vars:
        eps:    global constant ehg_epsilon
        params:

    Returns:
        simparams:
    """

    from os import path
    import glob
    outfolder = params['sScenarioName']
    outfile_pattern = 'eps??sim??.out'
    search_pattern = path.normpath(path.join(outfolder, outfile_pattern))
    previously_done = glob.glob(search_pattern)
    for n in range(params['nSim']):
        for i, eps in enumerate(params['epsilon']):
            assert type(eps) is float
            sim_name = 'eps{0:02d}sim{1:02d}'.format(i, n)
            outfilename = path.normpath(path.join(outfolder, sim_name + '.out'))
            if path.normpath(outfilename) in previously_done:
                print 'Skip simulation ' + sim_name + ' (already completed).'
                continue
            else:
                """fresh copy of params used for each replicate.
                Get phi for this simulation.
                phi needs sim_epsilon=eps
                """
                seed = params['nRndSeed'], n, int(eps * 10 ** 5)
                simparams = params.copy()
                simparams.update(
                    prng=RandomState(seed=seed),
                    sim_epsilon=eps,
                    sim_name=sim_name,
                    outfilename=outfilename,
                )
                yield simparams

def buildPopulation(params):
    """Add docstring
    """
    nM, nF = params["nMF"]
    sim_days = params['nSimDays']
    burn_days = params['nBurnDays']
    out_interval = params['nOutputInterval']  # usu. 365 (i.e., once a year)
    p_nclients = params.get('p_nclients',0.0)
    p_nsexworkers = params.get('p_nsexworkers',0.0)
    p_nM_ART = params.get('p_nM_ART',0.0) 
    p_nF_ART = params.get('p_nF_ART',0.0)
    p_miners = params.get('p_miners',0.0)
    p_PREP = params.get('p_PREP',0.0)
    beta_PREP = params.get("beta_PREP",(0.,0.,0.,0.))
    params['Disease'] = stagedHIVfactory(durations=params['durationsHIV'],
                                         transM2F=params['beta_M2F'],
                                         transF2M=params['beta_F2M'],
                                         )

    params['sim_phi'] = lambda male, female: params['model_phi'](male, female, params['sim_epsilon'])
    params['counters'] = counters = defaultdict(int)
    params["DiseasePREP"] = stagedHIVfactory(durations=params['durationsHIV'],
                                              transM2F=beta_PREP,
                                              transF2M=beta_PREP,
                                              )
    params["DiseaseART"] = stagedHIVfactory(durations=params['durationsHIV'],
                                             transM2F=params['beta_M2F_ART'] if p_nM_ART>0 else (0.,0.,0.,0.),
                                             transF2M=params['beta_F2M_ART'] if p_nF_ART>0 else (0.,0.,0.,0.)
                                             )
    
    mPopComb = [('Person', {'sex': 'M', 'registry': None, 'params': params}, {'fraction': 1 - p_miners - p_nclients})]
    if p_nclients>0:
        mPopComb +=[('PersonCSWclient',{'sex':'M', 'registry': None, 'params':params},{'fraction':p_nclients})]
    if p_miners>0:
        mPopComb +=[('PersonMiner', {'sex': 'M', 'registry': None, 'params': params}, {'fraction': p_miners})]

    fPopComb = [('Person', {'sex': 'F', 'registry': None, 'params': params}, {'fraction': 1 - p_nsexworkers})]
    if p_nsexworkers>0:
        fPopComb +=[('PersonCSW', {'sex': 'F', 'registry': None, 'params': params}, {'fraction': p_nsexworkers})]

    males = typesCombinations(copy.deepcopy(mPopComb), nM)
    females = typesCombinations(copy.deepcopy(fPopComb), nF)
    return males,females

def onesim(params):
    """Setup and run one replicate.

    TODO: explain in detail why are the following:
        * what are the issues in multiprocessing:
          - that prevent pickling of model_params with factory or phi?
          - That prevent counters be global?
        * Why is params['nOutputInterval'] usually 365 (i.e., once a year) ?
        * Why we seed infections only once after burn_days have passed?
        * Why do we process/have only fatal disease (for now)?

    Vars:
        tempfh: temp file to hold the resulst of one simulation
        params:

    Returns:
        None.

    Raises:
        KeyboardInterrupt: Used in code development.
    """

    t0 = time.time()
    try:
        sim_name = '{0}: {1}'.format(params['model_name'], params['sScenarioName'])
        print '\n\nBegin ' + sim_name + '\n' + pprint.pformat(params)

        
        # create a temp file to hold the results of one simulation
        outfolder = params['sScenarioName']
        tempfh = tempfile.NamedTemporaryFile(mode='w', suffix='.out', dir=outfolder, delete=False)
        assert params.get('fout', None) is None
        params['fout'] = tempfh  # used by `record_output`        
        tempfname = tempfh.name  # we'll rename this if the simulation runs to completion
        # prepare outfile by writing header
        header = 'nMinfect,nFinfect,'
        header += 'MPrimaryTrans,FPrimaryTrans,MAsymptomaticTrans,FAsymptomaticTrans,MSymptomaticTrans,FSymptomaticTrans,'
        header += 'MPships,,,,,,,,,,,,FPships,,,,,,,,,,,,iMPships,,,,,,,,,,,,iFPships,,,,,,,,,,,,'
        header += 'MPrimary,FPrimary'
        # Add header components one by one
        listHeader  = [header]
        listHeader += ['nFSWinf','nCSWinf','nMinerInf','nARTinf','nPREPinf']
        stages = ['primary','asymptomatic','sympotmatic']
        persons = ['FSW','CSW', 'ART', 'Miner', 'PREP'] 
        for pers in persons:
            for stage in stages:
                listHeader += ['n%s_%s'%(pers,stage)]
                    
        listHeader += ['nFSWARTinf','nCSWARTinf','nMinARTinf','nFSWPREPinf','nCSWPREPinf','nMinPREPinf']
        listHeader += ['nFSWprim','nCSWprim','nMinprim','nARTprim','nPREPprim']
        # then join them together comma-delimited
        header=",".join(listHeader)        
        tempfh.write(header)

        nM,nF=params['nMF']
        sim_days = params['nSimDays']
        burn_days = params['nBurnDays']
        out_interval = params['nOutputInterval']  # usu. 365 (i.e., once a year)
        p_nclients = params.get('p_nclients',0)
        p_nsexworkers = params.get('p_nsexworkers',0)
        params['counters'] = counters = defaultdict(int)

        for _type in (nM, nF, sim_days, burn_days, out_interval):
            assert_type(_type, int)

        males,females=buildPopulation(params)
        schedule = Scheduler(params=params)
        for m in males: schedule.register_person(m)
        for f in females: schedule.register_person(f)

        prng = params['prng'] #JK: explain this assignment of prng here.

        nclients = int(p_nclients * len(males))
        nsexworkers = int(p_nsexworkers * len(females))
        clients = prng.choice(males, nclients)

        # begin simulation loop /* Do the simulations */
        for day in range(sim_days + burn_days):
            logging.info('\nBEGIN ITERATION for day {0}.'.format(day))
            logging.debug(schedule.show_one_day(day))

           # seed_infections_after_burn_days_passed(day, burn_days, deathday,
           #         schedule, params)

            # Seed infections after burn_days have passed (ONCE)
            if (day == burn_days):
                assert schedule.count_scheduled_deaths() == 0
                diseases = seed_infections(males, females, day, schedule=schedule, params=params)
                assert schedule.count_scheduled_deaths() == len(diseases)  # for now only have fatal disease
                assert all(deathday >= day for deathday in schedule.deaths)  # equal only possible on day==burn_days
            
            start_treatment(males, females, day, params) #administer treatment

            """Run the core of the simulation (runs even during burn days).
            """
            schedule.coresim(
                males=males,
                females=females,
                day=day,
                params=params
            )

            """Record the output once a "period" (i.e., every out_interval days)
                :note: this won't record after last year is run (it needs one more day to pass the test).
                We keep it this way just to match EHG.
            """
            if (day >= burn_days and (day - burn_days) % out_interval == 0):
                print '.',
                sys.stdout.flush()
                outIndex = (day - burn_days) / out_interval
                record_output(males, females, params,day,params['nOutGroups'])
                # reset counters
                counters.clear()

            gc.collect()

        """END of simulation. 
        Clean up: reset static elements.
        Prepare classes for reuse. 
        Clear partnerships and transmissions multimaps.
        """
        schedule.clear_partnerships()  
        schedule.deaths.clear()

        tempfh.close()
        dt = time.time() - t0
        outfilename = params['outfilename']
        shutil.move(tempfname, outfilename)
        msg = """
        {sim_name} completed successfully in {minutes:.2f} minutes.
        {sim_name} output written to {outfilename}.
        """.format(sim_name=sim_name, minutes=dt / 60., outfilename=outfilename)
        logging.info(msg)
    except KeyboardInterrupt:
        logging.exception("Interrupted")
        sys.exit(0)
        raise

