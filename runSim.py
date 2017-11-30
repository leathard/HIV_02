"""Runner routines for the simulation. 

This code provides multiprocessing and sequential run options for running
simulation defined in Sim.py
Sequential run can be useful for debugging problems.
"""
import multiprocessing as mp
from mpi4py import MPI
import logging
from Sim import get_param_set


def runSimulations(oneSim, params):
    """Run all pending simulation for this model.

    Input is the 1-simulation routine and  a dict of model-specific parameters.
    Use it for debugging, because it's a simpler, sequential run where errors can only occurs one at a time.

    Returns: 
        None
    """
    for simparams in get_param_set(params):
        oneSim(simparams)


def runSimulationsMP(oneSim, params, workerpool=None):
    """Runs sumulations utilizing multiprocessing. 
    
    This is the same as `runSimulations`.
    Expected input is 1-simulation routine and the dictionary of input parameters for simulation run.
    It uses multiprocessing to distribute work to more processors.

    Allows for console interruption on linux (Zorro) and darwin (Mac). Relevant
    for development. 

    Returns: 
        None
    """
    import sys, os
    if 'linux' in sys.platform or 'darwin' in sys.platform:
        logging.info("Installing interrupt handler")
        import signal
        def signal_handler(signum, frame):
            logging.warn('Killing worker %d' % os.getpid())
            sys.exit(0)

        # only for *nix-based
        signal.signal(signal.SIGINT, signal_handler)

    if workerpool is None:
        n_processors = mp.cpu_count()
        print '{0} cores available'.format(n_processors)
        workerpool = mp.Pool(n_processors)

    if 'linux' in sys.platform or 'darwin' in sys.platform:
        signal.signal(signal.SIGINT, signal.SIG_DFL)

    try:
        workerpool.map(oneSim, [p for p in get_param_set(params)])
    except KeyboardInterrupt:
        logging.exception("INTERRUPTED")
        pass
    workerpool.close()
    workerpool.join()

#Use default communicator. 
#COMM = MPI.COMM_WORLD 

#Collect things to be done in a list. 
#if COMM.rank == 0: 
#    jobs = list(range(100))
#    jobs = split(jobs, COMM.size)

#else: 
#    jobs = None

#Scatter jobs across cores. 
#jobs = COMM.scatter(jobs, root=0)

#results = []
#for job in jobs: 
    
