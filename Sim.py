"""
Provides function ``onesim``, which can run a single concurrency simulation.
Compare to EHG's sim.cpp.
"""

from Infection import stagedHIVfactory
from Scheduler import *
from Persons import *
#from constants import *
from numpy.random import RandomState
import pprint
import tempfile, time
from collections import defaultdict
import gc
import shutil


def get_param_set(params):
    """Yield dict, a replicate specific set of parameters.
    There are n_sim (e.g., 100) replicates for each epsilon,
    and each one has its own random seed, sim_name, and outfilename.
    This generator function should yield dicts with **pickleable** elements only,
    so that it can be used with multiprocessing.
    """
    from os import path
    import glob
    outfolder = params['sScenarioName']
    outfile_pattern = 'eps??sim??.out'
    search_pattern = path.normpath(path.join(outfolder, outfile_pattern))
    previously_done = glob.glob(search_pattern)
    for n in range(params['nSim']):
        for i, eps in enumerate(params['epsilon']):  # ehg_epsilons is a global constant
            assert type(eps) is float
            sim_name = 'eps{0:02d}sim{1:02d}'.format(i, n)
            outfilename = path.normpath(path.join(outfolder, sim_name + '.out'))
            if path.normpath(outfilename) in previously_done:
                print 'Skip simulation ' + sim_name + ' (already completed).'
                continue
            else:
                seed = params['nRndSeed'], n, int(eps * 10 ** 5)  # replicate(n)-specific seed
                simparams = params.copy()  # fresh copy for each replicate (IMPORTANT!)
                simparams.update(
                    prng=RandomState(seed=seed),
                    # get phi for this simulation (to be passed to `random_pairings` via `coresim`)
                    sim_epsilon=eps,  # `phi` needs this!
                    sim_name=sim_name,
                    outfilename=outfilename,
                )
                yield simparams

def buildPopulation(params):
    #after updating params set in get_params_set
    # time the simulation
    nM,nF = params["nMF"]
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
    # note: don't add factory to model_params as it won't pickle -> MP problems
    params['Disease'] = stagedHIVfactory(durations=params['durationsHIV'],
                                         transM2F=params['beta_M2F'],
                                         transF2M=params['beta_F2M'],
                                         )

    # create the phi for this model
    # note: don't add to model_params as it won't pickle -> MP problems

    params['sim_phi'] = lambda male, female: params['model_phi'](male, female, params['sim_epsilon'])

    # counters used to tally incidence transmission by stage of infection
    # note: don't make this global (MP!)
    params['counters'] = counters = defaultdict(int)
    prep_params = params.copy()
    prep_params["Disease"] = stagedHIVfactory(durations=params['durationsHIV'],
                                              transM2F=beta_PREP,
                                              transF2M=beta_PREP,
                                              )
    params_ART = params.copy()
    params_ART["Disease"] = stagedHIVfactory(durations=params['durationsHIV'],
                                             transM2F=params['beta_M2F_ART'] if p_nM_ART>0 else (0.,0.,0.,0.),
                                             transF2M=params['beta_F2M_ART'] if p_nF_ART>0 else (0.,0.,0.,0.)
                                             )
    
    mPopComb = [(Person, {'sex': 'M', 'params': params}, {'fraction': 1 - p_PREP - p_miners - p_nM_ART - p_nclients})]
    if p_nclients>0:
        mPopComb +=[(Person01,{'sex':'M','params':params},{'fraction':p_nclients})]
    if p_nM_ART>0:
        mPopComb +=[(Person02, {'sex': 'M', 'params': params_ART}, {'fraction': p_nM_ART})]
    if p_PREP>0:
        mPopComb +=[(Person04, {'sex': 'M', 'params': prep_params}, {'fraction': p_PREP})]
    if p_miners>0:
        mPopComb +=[(Person03, {'sex': 'M', 'params': params}, {'fraction': p_miners})]
    
    fPopComb = [(Person, {'sex': 'F', 'params': params}, {'fraction': 1 - p_nF_ART - p_nsexworkers - p_PREP})]
    if p_nF_ART>0:
        fPopComb +=[(Person02, {'sex': 'F', 'params': params_ART}, {'fraction': p_nF_ART})]
    if p_nsexworkers>0:
        fPopComb +=[(Person01, {'sex': 'F', 'params': params}, {'fraction': p_nsexworkers})]
    if p_PREP>0:
        fPopComb +=[(Person04, {'sex': 'F', 'params': prep_params}, {'fraction': p_PREP})]
    import copy
    males = typesCombinations(copy.deepcopy(mPopComb), nM)
    females = typesCombinations(copy.deepcopy(fPopComb), nF)
    return males,females

def onesim(params):
    """Return None; setup and run one replicate
    (e.g., all daily iterations for 1 out of 100 replicates, for 1 value of epsilon)
    """
    t0 = time.time()
    import sys
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
        # type check parameters

        assert type(nM) is int
        assert type(nF) is int
        assert type(sim_days) is int
        assert type(burn_days) is int
        assert type(out_interval) is int
        
        males,females=buildPopulation(params)                
        schedule = Scheduler(params=params)  # event scheduler (transmissions, deaths, dissolutions)
        for m in males: schedule.register_person(m)
        for f in females: schedule.register_person(f)

        prng = params['prng']

        # ai: look for data!
        # jk: vandepitte (2006); carael (2006) data        
        nclients = int(p_nclients * len(males))
        nsexworkers = int(p_nsexworkers * len(females))

        clients = prng.choice(males, nclients)

        # begin simulation loop /* Do the simulations */
        for day in range(sim_days + burn_days):
            logging.info('\nBEGIN ITERATION for day {0}.'.format(day))
            logging.debug(schedule.show_one_day(day))

            # Seed infections after burn_days have passed (ONCE)
            if (day == burn_days):
                assert schedule.count_scheduled_deaths() == 0
                diseases = seed_infections(males, females, day, schedule=schedule, params=params)
                assert schedule.count_scheduled_deaths() == len(diseases)  # for now only have fatal disease
                assert all(deathday >= day for deathday in schedule.deaths)  # equal only possible on day==burn_days

            # run the core of the simulation (runs even during burn days)
            # params holds counters and fout
            schedule.coresim(
                males=males,
                females=females,
                day=day,
                params=params
            )

            # Record the output once a "period" (i.e., every out_interval days)
            # :note: this won't record after last year is run (it needs one more day to pass the test).
            #        We keep it this way just to match EHG.
            if (day >= burn_days and (day - burn_days) % out_interval == 0):
                print '.',
                sys.stdout.flush()
                outIndex = (day - burn_days) / out_interval
                # ai: recording strategy
                #    params holds counters and fout
                record_output(males, females, params,day,params['nOutGroups'])
                # reset counters
                counters.clear()

            gc.collect()
        # END of simulation; just need to clean up: reset static elements
        # ai: mostly don't need to clean up in Python unless we intend to reuse objects
        # (and EHG don't even reuse the persons, so the rest of object creation is a trivial cost)
        # prepare classes for reuse
        schedule.clear_partnerships()  # clears the `partnerships` and `transmissions` multimaps
        schedule.deaths.clear()

        tempfh.close()
        dt = time.time() - t0
        # since the simulation ran to completion, we can use the output file
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
        # except TypeError,e:
        #    logging.exception("onesim")
        #    raise
