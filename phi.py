def phi_ehg(male, female, epsilon):
    """Return float, the probability of partnership formation,
    **given** that ``male`` and ``female`` are selected as potential partners.
    Based on individuals' current number of partnerships, and mixing parameter
    :note: epsilon is 1 if no concurrency and 0 if no resistance to concurrency
    :note: we will use closure to set the `epsilon` local to the `sim`
    """
    if (female.n_partners == 0 and male.n_partners == 0):
        return 1.0
    else:
        return (1.0 - epsilon)  # epsilon controls concurrency


def phi2(male, female, epsilon):
    """Return float, the probability of partnership formation,
    **given** that `male` and `female` are selected as potential partners.
    Based on each individual's current number of partnerships,
    and mixing parameter
    Sawers, Isaac, Stillwaggon - coital dilution example.
    :note:
      epsilon controls concurrency (1-> none, 0 -> no resistance),
      eF is female resistance to concurrency,
      eM is male resistance to concurrency,
      setting eF>eM (since epsilon in (0,1))
    """
    eF = 0 if female.n_partners == 0 else math.sqrt(epsilon)
    eM = 0 if male.n_partners == 0 else epsilon * epsilon
    return (1.0 - eF) * (1 - eM)  # epsilon controls concurrency


def phi3(male, female, epsilon):
    """Return float, the probability of partnership formation,
    **given** that `male` and `female` are selected as potential partners.
    Based on each individual's current number of partnerships,
    and mixing parameter ``epsilon``, which controls resistance
    to concurrency`
    note:
      epsilon controls concurrency (1-> none, 0 -> no resistance),
    jk: Here we add commercial sex workers with no
    it is the number of partners of the female that is controlling.
    Males will always accept a new partner who has no partners,
    and females without a partner will always accept a new partner.
    Otherwise partnerships form with prob 1-eps.
    """
    if (female.n_partners == 0):
        return 1.0
    else:
        return (1.0 - epsilon)  # epsilon controls concurrency


# jk: I do not understand why below averages are necessary?
def avtrans2asymtrans(beta, reltransfm):
    """Return dict, f2m->list and m2f->list.
    For average (stage-dependent) transmission rates beta,
    computes the f2m and m2f rates that give that average
    but have relative f2m transmission of reltransfm.

    Each transmission rate beta_i needs to be split into
    a M -> F and F -> M version so that two otherwise
    identical partnerships, one male infected and one
    female infected, will transmit at this rate on average.

    E.g., solving for reltransfm=0.7:
    b = (m + f)/2  and .7m = f
    yields
    b = 1.7m/2 or m= 2b/1.7

    beta : list of float
      the stage-dependent average transmission rates
    reltransfm : float
      female to male transmission rate relative to male to female

      #jk: fcsw to male transmission rate over different stages
      is the same as female to male transmission rate.
    """
    from copy import copy
    assert (reltransfm <= 1.0)

    mscale = 2. / (1 + reltransfm)
    m2f = tuple(mscale * bi for bi in beta)
    f2m = tuple(reltransfm * bi for bi in m2f)
    return {'m2f': m2f, 'f2m': tuple(reltransfm * bi for bi in m2f), 'fcsw2m': copy(f2m)}
