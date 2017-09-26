"""
Provides an intial set of params that can be altered to add specific data.
"""
from constants import *

"""
11 different values of epsilon (the concurrency parameter)
EHG's R version: epsilon = 1-1.02*(exp(3.9*(seq(0,1,.1)-1))-.02)  
:note: lowest index is biggest (approx 1) epsilon
"""
nOutGroups = 12  # based on EHG's Rsim.cpp; otherwise arbitrary
ehg_epsilons = tuple(1 - 1.02 * (math.exp(3.9 * (xi / 10. - 1)) - .02) for xi in range(nOutGroups - 1))

# WARNING: shared_params is not a complete parameter set! (some vals are None)

shared_params = dict(
    model_name='',  # name of simulation model
    n_sim=100,  # IMPORTANT!!! the number of replications for each parameter set (each value of epsilon)
    pop_size=2 * 10 ** 4,  
    pctHIVseed= 12, #1,
    burn_days=365 * 5,
    sim_days=365 * 30,
    out_interval=365,
    rndseed=7930881,  # From EHG code
    beta_M2F=None,
    beta_F2M=None,
    beta_FCSW2M=None,
    beta_M2F_ART=None,
    beta_F2M_ART=None,
    # durations during each stage: a sequence
    dur=None,
    # aggregate partnership formation paramter
    rho=None,
    sigma=None,
    p_nclients=None,
    p_nsexworkers=None,
    p_nM_ART=None,
    p_nF_ART=None
)

# get the shared params
ehg_staged01 = shared_params.copy()
beta_ehg  = (beta_p, beta_a, beta_s, beta_0)
beta_ART  = (beta_p_ART, beta_a_ART, beta_s_ART, beta_0)
beta_PREP = (beta_p_PREP,beta_a_PREP,beta_s_PREP,beta_0)

# update the sim specific params
ehg_staged01.update(
    sim_name='ehg-staged',
    rho=mk.rho,  
    sigma=mk.sigma,  # chk where is this partnership dissolution param used?
    # daily infection rate during each stage: a sequence
    beta_M2F=beta_ehg,
    beta_F2M=beta_ehg,
    beta_FCSW2M=beta_ehg,
    beta_M2F_ART=beta_ART,
    beta_F2M_ART=beta_ART,
    # durations during each stage: a sequence
    dur=(dur_p, dur_a, dur_s, dur_0),
    p_nclients=0,
    p_nsexworkers=0,
    p_nM_ART=0,
    p_nF_ART=0
)
