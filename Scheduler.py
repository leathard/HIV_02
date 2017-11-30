"""
Contains Scheduler and its helper routines.

Most of the code is by Alan G. Isaac. In Sim.py. 
Modified treatment and tracking of agents. 

"""
from collections import defaultdict
from constants import *
import logging
from Infection import Partnership
import numpy as np
from Persons import *

def random_pairings(n_pairs, males, females, params):
    """Return n_pairs random male-female pairs (**with** replacement).

    n_pairs : int
        the number of pairs to form
    males : list
        the potential male partners
    females : list
        the potential female partners
    params : dict
        params['prng'] is the random number generator;
        params['phi'] is a function that returns partnership formation probability
    """
    _phi = params['sim_phi']
    _prng = params['prng']
    _randint = _prng.randint
    _random_sample = _prng.random_sample
    if n_pairs < 0:
        raise ValueError('cannot form a negative number of pairs')
    nM = len(males)
    nF = len(females)
    assert n_pairs < nM and n_pairs < nF
    while n_pairs:
        # `n_pairs` random males and females (*with* replacement)
        males2pair = (males[idxM] for idxM in _randint(0, nM, n_pairs))
        females2pair = (females[idxF] for idxF in _randint(0, nF, n_pairs))
        draws = _random_sample(n_pairs)
        # must form pairs sequentially due to concurrency resistence
        for (male, female, draw) in zip(males2pair, females2pair, draws):
            # form partnership with probability phi, which varies with concurrency
            if test_random_pair(male, female, draw, _phi):
                n_pairs -= 1
                yield male, female


def test_random_pair(male, female, draw, phi):
    pairup = False
    if male.comsex and female.comsex:
        if female not in male.partners:
            pairup = True
        logging.debug("commercial partnership")
    elif (draw < phi(male, female)):  # ai: EHG use binomial; why?
        # ai: check that partnership doesn't already exist
        if female not in male.partners:
            assert male not in female.partners
            pairup = True
    return pairup



def seed_infections(males, females, day, schedule, params):
    assert sum(f.is_infected for f in females) == 0
    assert sum(m.is_infected for m in males) == 0
    logging.info('Seed infections.')
    Disease = params['Disease']
    nM = len(males)
    nF = len(females)
    nMaleSeed = int(round((nM // 100) * params['pctHIVseed']))
    nFemaleSeed = int(round((nF // 100) * params['pctHIVseed']))
    if not (nFemaleSeed < len(females) and nMaleSeed < len(males)):
        raise ValueError('choose smaller seeds')
    _prng = params['prng']
    idxF = _prng.permutation(nF)[:nFemaleSeed]  
    idxM = _prng.permutation(nM)[:nMaleSeed]
    diseases = list()

    for idx in idxF:
        diseases.append(females[idx].HIVSeedInfect(day, schedule))
    for idx in idxM:
        diseases.append(males[idx].HIVSeedInfect(day, schedule))

    assert nFemaleSeed == sum(f.is_infected for f in females)
    assert nMaleSeed == sum(m.is_infected for m in males)
    return diseases


def tally_transmission(day, transmitter, counters):
    """Tally the transmission (HIV stage for output).

    Returns: None.

    :note: EHG increment these in partnership.py; we do this instead
    """
    stage = transmitter.get_HIV_stage(day)
    if transmitter.sex == 'M':
        if (stage == 'primary'):
            counters['malePrimaryTransToday'] += 1
        elif (stage == 'asymptomatic'):
            counters['maleAsymptomaticTransToday'] += 1
        elif (stage == 'sympotmatic'):
            counters['maleSymptomaticTransToday'] += 1
    elif transmitter.sex == 'F':
        if (stage == 'primary'):
            counters['femalePrimaryTransToday'] += 1
        elif (stage == 'asymptomatic'):
            counters['femaleAsymptomaticTransToday'] += 1
        elif (stage == 'sympotmatic'):
            counters['femaleSymptomaticTransToday'] += 1
    else:
        raise ValueError('Unknown sex: ' + str(transmitter.sex))

def start_treatment(males, females, day, params):
    """Randomly pick agents and assign treatment.

    note: treated agents scheduled to die are not automatically replaced by
    another treated agent. 
    """
    p_nM_ART = params.get('p_nM_ART', 0.0)
    p_nF_ART = params.get('p_nF_ART', 0.0)
    p_PREP = params.get('p_PREP', 0.0)
    nF_ART = sum(1 for female in females if female.ART)
    nM_ART = sum(1 for male in males if male.ART)
    nART = nF_ART + nM_ART
    nF_PREP = sum(1 for female in females if female.PREP)
    nM_PREP = sum(1 for male in males if male.PREP)
    nPREP = nF_PREP + nM_PREP

    newARTm = int(round(p_nM_ART * len(males) - nM_ART))
    newARTf=int(round(p_nF_ART * len(females) - nF_ART))
    if newARTm > 0:
        nARTm = [idx for idx in range(len(males)) if males[idx].is_infected and not males[idx].ART and not males[idx].PrEP]
        #print(day)
        #print(len(nAM))
        #print(sART_M)
        if len(nAm)>0:
            sARTm = np.random.choice(min(len(nARTm),newARTm)) + 1 
            sARTmi = np.random.choice(nAm, sARTm, replace = False) 
            #start treatment
            for idx in sARTmi:
                males[idx].ART = True
                males[idx._Disease] = params['DiseaseART']

    if newARTf > 0:
        nAf = [idx for idx in range(len(females)) if females[idx].is_infected and not females[idx].ART and not females[idx].PrEP]
        if len(nAf)>0:
            sARTf = np.random.choice(min(len(nARTf),newARTf))+1
            sARTfi = np.random.choice(nAf, sARTf, replace = False)
            #start treatment 
            for idx in sARTfi: 
                females[idx].ART = True
                females[idx]._Disease = params['DiseaseART']

    newPREP = int(round((p_PREP * (len(males) + len(females)) - nPREP)/2))
    if newPREP > 0:
        nPREPm = [idx for idx in range(len(males)) if not males[idx].is_infected
                and not males[idx].ART and not males[idx].PREP]
        sPREPm= np.random.choice(newPREP) + 1
        if len(nPREPm)>0:
            sPREPmi = np.random.choice(nPREPm, sPREPm, replace = False)
            for idx in sPREPmi:
                males[idz].PREP = True
                males[idx]._Disease = params['DiseasePREP']

            nPREPf = [idx for idx in range(len(females)) if not
                    females[idx].is_infected and not females[idz].ART and not
                    females[idx].PREP]
            sPREPf = np.random.choice(newPREP) + 1
            sPREPfi = np.random.choice(nPREPf, sPREPf, replace = False)
            for idx in sPREPfi: 
                females[idx].PREP = True
                females[idx]._Disease = params['DiseasePREP']


class Scheduler(object):
    """Provides a schedular object that
    maintains a schedule of future events
    and also keeps a registry of other objects.
    """

    def __init__(self, params):
        self.params = params
        self.deaths = defaultdict(list)  # multimap for scheduling deaths: map date -> list of Person
        self.cures = defaultdict(list)  # multimap for scheduling cures: map date -> list of Person
        self.transmissions = defaultdict(list)  # map **upcoming** transmission date to list of Partnership
        self.dissolutions = defaultdict(list)  # map each **end** date list of Partnership
        # list of objects (switch to weakref set?)
        self.n_people = 0  # initialize partnership counter
        self.males = list()
        self.females = list()
        self.n_partnerships = 0  # initialize partnership counter
        self.partnerships = list()

    def coresim(self, males, females, day, params):
        """Return None; modify `params`.
        Run the core stages of the simulation (one iteration):
        form partnerships, do HIV transmissions,
        dissolve partnerships, do HIV deaths.
        These are run every iteration (including burn days).
        """
        self.form_partnerships(males, females, day)
        self.hiv_transmissions(day)
        self.dissolve_partnerships(day)
        self.hiv_deaths(params, day)

    def register_person(self, person):
        if person.sex == 'M':
            self.males.append(person)
            self.n_people += 1
        elif person.sex == 'F':
            self.females.append(person)
            self.n_people += 1
        else:
            raise ValueError('unknown sex {}'.format(person.sex))
        assert self.n_people == len(self.males) + len(self.females)

    def register_infection(self, disease):
        end = disease.end_date
        host = disease._host
        if disease.is_fatal:
            if host in self.deaths[end]:
                raise ValueError()
            self.deaths[end].append(host)
        else:
            raise ValueError('all currently fatal')
            # self.cures[end].append(host)

    def register_transmission(self, partnership):
        self.transmissions[partnership.transmission_date].append(partnership)

    def register_partnership(self, pship):
        # add partnership to partnership multimap with end date
        self.partnerships.append(pship)
        self.n_partnerships += 1
        assert self.n_partnerships == len(self.partnerships), \
            '{} != {}'.format(self.n_partnerships, len(self.partnerships))
        self.dissolutions[pship.end_date].append(pship)

    def delete_partnership(self, pship):
        """Return None; delete a partnership from ``partnerships`` multimap
        *before* scheduled end date (e.g., when partner dies from HIV).
        Also remove any scheduled transmission.
        """
        self.dissolutions[pship.end_date].remove(pship)
        self.partnerships.remove(pship)
        self.n_partnerships -= 1
        # destroy the partnership (match EHG's destructor for the partnership class)
        male, female = pship._mf
        male.remove_partnership(pship)  # removes from set *and* decrements this person's count
        female.remove_partnership(pship)
        # check if a future transmission is scheduled, and if so remove it
        if (pship.transmission_scheduled):
            self.transmissions[pship.transmission_date].remove(pship)

    def clear_partnerships(
            self):  # delete all partnerships (prepare for next simulation); compare EHG's DumpPartnerships
        # chk this may work very differently than the EHG code, but that shd't matter, it's just cleanup
        for lst in self.dissolutions.values():
            for pship in lst[:]:
                self.delete_partnership(pship)
        assert self.n_partnerships == 0
        self.dissolutions.clear()
        self.transmissions.clear()

    ######################  ai: convenience methods #######################################
    def count_partnerships(self):
        return sum(len(lst) for lst in self.dissolutions.values())

    def count_transmissions(self):
        return sum(len(lst) for lst in self.transmissions.values())

    def count_scheduled_deaths(self):
        return sum(len(lst) for lst in self.deaths.values())

    def show_one_day(self, day):
        """Return str, the schedule (so far) for `day`.
        (Used only for logging.)
        """
        males = self.males
        females = self.females
        report = """
        Schedule for day {day}
        ======================
        Number infected: {ninf} of {pop}
        Transmissions: {ntrans_total} total and {ntrans_today} today
        Deaths: {d_total} total and {d_today} today
        Dissolutions: {p_total} total and {p_today} today
        """.format(
            day=day,
            ninf=sum(p.is_infected for p in (males + females)),
            pop=len(males) + len(females),
            ntrans_total=sum(len(lst) for lst in self.transmissions.values()),
            d_total=sum(len(lst) for lst in self.deaths.values()),
            p_total=sum(len(lst) for lst in self.dissolutions.values()),
            # scheduled for today (BUT: avoid adding keys!)
            ntrans_today=len(self.transmissions.get(day, [])),
            d_today=len(self.deaths.get(day, [])),
            p_today=len(self.dissolutions.get(day, [])),
        )
        return report

    #################   MAIN SCHEDULE   ##################################
    def form_partnerships(self, males, females, day):  # chkchkchk schedule partnerships
        params = self.params
        ##### form partnerships
        logging.debug('coresim: BEGIN form partnerships')
        old_partnership_count = self.n_partnerships  # chkchkchk
        assert self.count_partnerships() == old_partnership_count
        # population sizes
        yearday = day % 365  # the day in the year
        filtered_males = filter(lambda x: x.is_available(yearDay=yearday), males)
        filtered_females = filter(lambda x: x.is_available(yearDay=yearday), females)
        nM = len(filtered_males);
        nF = len(filtered_females);
        mdiff = len(males) - nM
        fdiff = len(females) - nF

        if len(males) > len(filtered_males):
            msg = "day: {day}, year day: {yday} Filtered out {nminers} from partnership formation". \
                format(nminers=mdiff, day=day, yday=yearday)
            logging.info(msg)
        if len(females) > len(filtered_females):
            msg = "day: {day}, year day: {yday} Filtered out {nminers} from partnership formation". \
                format(nminers=fdiff, day=day, yday=yearday)
            logging.info(msg)

        # we follow EHG and impose a hard ceiling on the number of partnerships!
        # (there a lots of things not to like about this)
        max_new_partnerships = (nM + nF) // 2 - self.n_partnerships
        assert type(max_new_partnerships) is int
        # determine how many new partnerships to form (binomial draw with probability rho)
        prng = params['prng']
        rho = params['rho']  # aggregate partnership formation paramter
        assert type(rho) is float
        nFormPartnerships = 0  # needed to match EHG behavior (R's rbinom(0, rho)=0)
        if max_new_partnerships > 0:
            nFormPartnerships = prng.binomial(n=max_new_partnerships, p=rho)
        # form exactly nFormPartnerships new partnerships
        npairings = 0  # counter for testing that `random_pairings` returns the right number
        #:note: random pairings uses phi, which is in `params`
        for male, female in random_pairings(nFormPartnerships, filtered_males, filtered_females, params):
            newpartnership = Partnership(male, female, day, registry=self, params=params)  # form partnership
            npairings += 1
        assert npairings == nFormPartnerships
        logging.debug('coresim: END form partnerships')
        logging.info('\t{0} partnerships formed'.format(nFormPartnerships))
        assert self.count_partnerships() == old_partnership_count + nFormPartnerships

    def hiv_transmissions(self, day):
        params = self.params
        ##### do HIV transmissions
        logging.debug('coresim: begin HIV transmission')
        # recall: schedule.transmissions maps days to lists of partnerships
        Disease = params['Disease']
        transmissions4today = self.transmissions[day]  # list of partnerships
        ntrans_today = len(transmissions4today)
        ntrans_waiting = self.count_transmissions()
        logging.info('do {0} of {1} transmissions'.format(ntrans_today, ntrans_waiting))
        for pship in transmissions4today:
            transmitter, disease = pship.transmit(day)  # chk is there a better way?
            if transmitter is not None:
                assert disease is not None
                tally_transmission(day=day, transmitter=transmitter, counters=params['counters'])
        del self.transmissions[day]
        del transmissions4today
        # lose ntrans_today future transmissions but may have gotten some new via infections at transmission
        assert self.count_transmissions() >= ntrans_waiting - ntrans_today
        logging.debug('coresim: end HIV transmission')

    def dissolve_partnerships(self, day):
        ##### dissolve partnerships
        logging.debug('begin: do scheduled partnership dissolutions')
        partnerships4today = self.dissolutions[day]
        logging.info('Do {0} dissolutions'.format(len(partnerships4today)))
        for pship in list(partnerships4today):  # *termination* dates map to partnership lists
            self.delete_partnership(pship)
        assert not self.dissolutions[day]  # shd be empty
        del self.dissolutions[day]
        del partnerships4today
        logging.debug('end: do scheduled partnership dissolutions')

    def hiv_deaths(self, day):
        ##### do HIV deaths
        logging.info('begin: HIV deaths')
        deaths4today = self.deaths[day]
        logging.info('do {0} deaths'.format(len(deaths4today)))
        assert len(deaths4today) == len(set(deaths4today))
        for indiv in tuple(deaths4today):  # don't really need a copy (chk Person.die)
            self.kill(indiv)
        del self.deaths[day]
        del deaths4today
        logging.info('end: HIV deaths')

    def kill(self, indiv):  # compare EHG's `Kill`
        """Return None; "kills" the individual (i.e., reset, so population is constant).
        Terminate partnerships and reset individual to unpartnered susceptible.
        """
        logging.info('ENTER: kill indiv')
        # terminate all partnerships
        for pship in tuple(indiv.partnerships):
            self.delete_partnership(pship)
        assert len(indiv.partnerships) == 0
        indiv.reset()  # reborn anew, uninfected, same sex
        indiv.reset_treatment(params) #reborn not readily treated
        logging.info('EXIT: kill indiv')


#### add some helper functions


def record_output(males, females, params,day):
    """Return None. Write data to disk.
    :note: EHG store these values in arrays; we write them to disk instead
    :note: record same values (in same order) as EHG, then append additions
    """
    counters = params['counters']  # count transmissions by stage
    fout = params['fout']
    n_infected_males = sum(m.is_infected for m in males)
    n_infected_females = sum(f.is_infected for f in females)
    # build up this period's data as a list
    data = [n_infected_males, n_infected_females]
    data += [
        counters['malePrimaryTransToday'],
        counters['femalePrimaryTransToday'],
        counters['maleAsymptomaticTransToday'],
        counters['femaleAsymptomaticTransToday'],
        counters['maleSymptomaticTransToday'],
        counters['femaleSymptomaticTransToday'],
    ]
    maleDistPartnerships = np.bincount([male.n_partners for male in males])
    assert sum(maleDistPartnerships) == len(males)
    femaleDistPartnerships = np.bincount([female.n_partners for female in females])
    assert sum(femaleDistPartnerships) == len(females)
    if n_infected_males:
        maleDistHIVp = np.bincount([male.n_partners for male in males if male.is_infected])
    else:  # bincount balks at empty lists
        maleDistHIVp = [0] * nOutGroups
    if n_infected_females:
        femaleDistHIVp = np.bincount([female.n_partners for female in females if female.is_infected])
    else:
        femaleDistHIVp = [0] * nOutGroups
    # the following were provided as output arrays in the EHG code
    for dist in maleDistPartnerships, femaleDistPartnerships, maleDistHIVp, femaleDistHIVp:
        distlen = len(dist)
        data.extend(dist[:nOutGroups])
        if distlen < nOutGroups:
            data.extend([0] * (nOutGroups - distlen))
        elif distlen > nOutGroups:
            logging.warn('discarding high partnership counts')
    assert len(data) == 2 + 6 + 4 * nOutGroups
    assert all(isinstance(item, int) for item in data)
    # above shd match EHG's output; below are additions
    nMprimary = sum(1 for male in males if male.has_primary())
    nFprimary = sum(1 for female in females if female.has_primary())
    data += [nMprimary, nFprimary]

    nFSW   = sum(1 for female in females if female.is_infected and female.comsex)
    nCSW   = sum(1 for male in males if male.is_infected and male.comsex)
    nMin   = sum(1 for male in males if male.is_infected and isPerson(male,Person03))
    nF_ART = sum(1 for female in females if female.is_infected and female.ART)
    nMART  = sum(1 for male in males if male.is_infected and male.ART )
    nART   = nF_ART + nMART
    nFPREP = sum(1 for female in females if female.is_infected and female.PrEP)
    nMPREP = sum(1 for male in males if male.is_infected and male.PrEP)
    nPREP  = nFPREP + nMPREP
    data += [nFSW,nCSW,nMin,nART,nPREP]
    
    stages = ['primary','asymptomatic','sympotmatic']
    persons = [Person01, Person02, Person03, Person04] # what about Person?
    
    dataStages=[]

    for pers in persons:
        for stage in stages:
            countF = sum(1 for female in females if isPersonInfectedInStage(female,pers,stage,day))
            countM = sum(1 for male in males if isPersonInfectedInStage(male,pers,stage,day))
            count  = countF + countM
            if pers in [Person01]:
                dataStages += [countF,countM]                
            else:
                dataStages+= [count]
    #[ Person0x in stage primary,  Person0x in stage asymptomatic ..., Person0x+1 in stage primary, ...]
    data += dataStages
    
    nFSWART = sum(1 for female in females if female.comsex  and female.is_infected and female.ART)
    nCSWART = sum(1 for male in males if male.comsex  and male.is_infected and male.ART)
    nMinART = sum(1 for male in males if isPerson(male,Person03)  and male.is_infected and male.ART)
    
    nFSWPREP = sum(1 for female in females if female.comsex  and female.is_infected and female.PrEP)
    nCSWPREP = sum(1 for male in males if male.comsex  and male.is_infected and male.PrEP)
    nMinPREP = sum(1 for male in males if isPerson(male,Person03)  and male.is_infected and male.PrEP)
    
    data += [nFSWART,nCSWART,nMinART,nFSWPREP,nCSWPREP,nMinPREP]

    nFSWprim = sum(1 for female in females if female.comsex  and female.is_infected and female.ART and female.has_primary())
    nCSWprim = sum(1 for male in males if male.comsex  and male.is_infected and male.ART and male.has_primary())
    nF_ARTprim = sum(1 for female in females if female.is_infected and female.ART and female.has_primary())
    nMARTprim = sum(1 for male in males if male.is_infected and male.ART and male.has_primary()) 
    nARTprim = nF_ARTprim+nMARTprim
    nMinprim = sum(1 for male in males if isPerson(male,Person03)  and male.is_infected and male.has_primary())

    nFPREPprim = sum(1 for female in females if female.is_infected and female.PrEP and female.has_primary())
    nMPREPprim = sum(1 for male in males if male.is_infected and male.PrEP and male.has_primary())
    nPREPprim  = nFPREPprim+nMPREPprim

    data += [nFSWprim,nCSWprim,nMinprim,nARTprim,nPREPprim]
    
    # finally, write the data to file
    fout.write('\n')
    data = ','.join(str(d) for d in data)
    fout.write(data)
