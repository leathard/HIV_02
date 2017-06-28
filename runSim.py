"""
Runner routines for the simulation.
Both multiprocessing run and sequential run are available.
Sequential run can be useful for debugging problems.
"""
import multiprocessing as mp
import logging
from Sim import get_param_set


def runSimulations(oneSim, params):
    """Return None. Run all pending simulation for this model.
    Input is the 1-simulation routine and  a dict of model-specific parameters.
    Use it for debugging, because it's a simpler, sequential run where errors can only occurs one at a time.
    """
    for simparams in get_param_set(params):
        oneSim(simparams)


def runSimulationsMP(oneSim, params, workerpool=None):
    """Doesn't return.
    Runs sumulations, same as runSimulation.
    Expected input is 1-simulation routine and the dictionary of input parameters for simulation run.
    It uses multiprocessing to distribute work to more processors.
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
