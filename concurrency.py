"""Provides key classes for concurrency simulations.

The code in this module began as a replication of Eaton,
Hallet, and Garnett (2010) and is very heavily influenced by
their online C++ code (person.cpp and partnership.cpp).
http://www.springerlink.com/content/22t477k235477463/10461_2010_Article_9787_ESM.html

:note: EHG used rgeom(p)+1, which is the same as np.random.geometric(p)!
  (i.e., numpy uses the shifted geometric; R does not)
:note 2013-12-03: replace schedule with Schedule object, passed explicitly.
"""

import logging
from collections import defaultdict


def stagedHIVfactory(durations, transM2F, transF2M):
    """Return class for staged HIV infection.
    (We adopt this strategy to allow multiprocessing,
    since class variables must be assigned at runtime.)
    """
    class HIV(Infection):
        _durations = tuple(durations)
        _duration = sum(_durations)
        _beta_M2F = tuple(transM2F)
        _beta_F2M = tuple(transF2M)
        is_fatal = True
        assert len(_durations)==len(_beta_F2M)==len(_beta_M2F)
    return HIV

class Infection(object):
    """
    Provides a partially **ABSTRACT** infection class.
    DO NOT USE THIS CLASS.
    Instead, use the class factory `stagedHIVfactory`
    in order to add the following vital class variables:
    _durations : duration of each infection stage
    _duration : total duration of the infection
    _beta_M2F  : staged transmission rates for infected males
    _beta_F2M  : staged transmission rates for infected females
    """
    _durations = NotImplemented
    _beta_M2F = NotImplemented
    _beta_F2M = NotImplemented
    def __init__(self, host, day, registry):
        """Return None.
        Record start date and host.
        Schedule host death.
        """
        self._host = host
        self.start_date = day
        self.end_date = day + self.duration
        self.registry = registry
        if registry is not None:
            registry.register_infection(self)
    def stage(self, day):
        """Return str, the name of the current infection stage.
        (Utilized to tally transmissions.)
        Called by `Person.get_HIV_stage`,
        which used to hold this code.
        """
        assert self._host.is_infected
        dur_p, dur_a, dur_s, dur_0 = self._durations
        days_infected = day - self.start_date
        if(days_infected <= dur_p): # primary infection
            stage = 'primary'
        elif( days_infected <= (dur_p + dur_a) ): # asymptomatic infection
            stage = 'asymptomatic'
        elif( days_infected <= (dur_p + dur_a + dur_s) ): # sympotmatic infection
            stage = 'sympotmatic'
        else: # pre-death, no transmission
            stage = 'dying'
        return stage
    def expose(self, partnership, day): #called by Partnership::expose_transmission
        """Return None; expose a *partnership*;
        and possibly schedule a transmission event.
        
        We expose the partnership to infection
        by sampling a candidate transmission date during
        each stage based on the transmission rate during that stage,
        checking to see if the candidate date is before the end of that stage
        and before the end of the partnership.
        If so, transmission is scheduled for that date.
        If not, proceed to the next stage.
        If transmission is not scheduled during any stage,
        no transmission occurs during the partnership.

        Called by `Partnership.expose_transmission`,
        which is called upon partnership formation or upon infection.

        :comment:
            We schedule a transmission only if serodiscordant
            (exactly one of the two is infected),
            but with concurrency another partner might still infect
            the uninfected partner before this scheduled transmission happens.
        """
        male, female = partnership.partners
        assert male.is_infected != female.is_infected
        sexfreq = partnership.sexfreq
        if male.is_infected:
            beta = self._beta_M2F
        if female.is_infected:
            beta = self._beta_F2M
        beta_p, beta_a, beta_s, beta_0 = ( sexfreq * beta_i for beta_i in beta )
        dur_p, dur_a, dur_s, dur_0 = self._durations
        candidate = eligDate = day  #candidate transmission day
        if(eligDate >= partnership.end_date):
            #ai: partnership ends beforehand; I think this only happens if eligDate==partnership.end_date
            # (because o/w partnership shd have dissolved) chk
            assert eligDate==partnership.end_date
            return

        stageEndDate = self.start_date #initialize to date of infection
        #ai: small change from EHG: loop over common code
        #ai: primary, asymptomatic, symptomatic infection stage
        schedule = self.registry  #can we most this out of here? chkchk
        for (dur_i,beta_i) in [(dur_p,beta_p),(dur_a,beta_a),(dur_s,beta_s),]:
            stageEndDate += dur_i  # infection period (increment then compare)
            if(eligDate < stageEndDate):
                if(beta_i > 0):  #daily transmission probability during this stage
                    candidate = eligDate + partnership._prng.geometric(beta_i)  #chk
                    if(candidate <= stageEndDate and candidate <= partnership.end_date ):
                        #schedule a transmission
                        partnership.transmission_scheduled = True
                        partnership.transmission_date = candidate
                        schedule.register_transmission(partnership)
                        logging.info('Infection.expose: Transmission scheduled')
                        return
                eligDate = stageEndDate  #if beta_i=0 or candidate too far out
            if(eligDate >= partnership.end_date):
                return

        #ai: final days (could reuse the above code)
        stageEndDate += dur_0
        if(eligDate < stageEndDate and beta_0 > 0):
            candidate = eligDate + partnership._prng.geometric(beta_0)
            if(candidate <= stageEndDate and candidate <= partnership.end_date):
                partnership.transmission_scheduled = True
                partnership.transmission_date = candidate
                schedule.register_transmission(partnership)
                return
    @property
    def host(self):
        return self._host
    @property  #ai: really a class property
    def duration(self):
        return self._duration

class Person(object):
    """Provides a basic person class.
    A person has a sex and a set of partnerships.
    A person may have a disease and an associated date of death.
    :todo: add date of birth for age contingent outcomes
    """
    def __init__(self, sex, registry, params):
        """Return None;  initialize a Person.
        """
        self._sex = sex  #'M' or 'F' (vs EHG who use integer: 1=male, 2=female)
        self._prng = params['prng']
        self.comsex = False  #commercial sex sales (F) or purchase (M)
        self.reset()
        if registry is not None:
            registry.register_person(self)

    def reset(self):
        """Return None.
        Prepare a Person instance for reuse after death.
        """
        self.n_partners = 0
        self.partnerships = list() #use list (not set) to ensure order of iteration
        self.is_infected = False  #ai: rendundant? (Just use doi)
        self.disease = None
        self.iDOD = None   #ai: impose max (adult) lifespan? chk (need day!)

    def get_HIV_stage(self, day):
        """Return str, the current stage of infection
        """
        return self.disease.stage(day)

    def add_partnership(self, partnership):
        """Return None; add a new partnership to self.partnerships list.
        """
        self.partnerships.append(partnership)
        self.n_partners += 1
        assert (self.n_partners == len(self.partnerships) == len(set(self.partners))) 
        return

    def remove_partnership(self, partnership):
        """Return None; remove a partnership from self.partnerships list.
        """
        logging.debug('ENTER: Person.remove_partnership')
        self.partnerships.remove(partnership)
        self.n_partners -= 1
        assert (self.n_partners == len(self.partnerships) == len(set(self.partners)))
        logging.debug('EXIT: Person.remove_partnership')

    def infect(self, Disease, day, registry):  #compare EHG's HIVInfect
        """Return Disease instance,
        infect an individual with HIV, schedule death, and expose other partners to transmission 
        """
        logging.debug('ENTER: Person.infect')
        #change to infected state
        self.is_infected = True
        #fatal disease will schedule death (chkchk change to list?)
        disease = self.disease = Disease(self, day, registry=registry)
        # expose current partners to disease transmission
        for pship in self.partnerships:
            pship.expose_transmission(day)
        logging.debug('EXIT: Person.infect')
        return disease

    def HIVSeedInfect(self, Disease, day, registry):  #compare EHG's HIVSeedInfect in person.cpp
        """Return None.  seed infection
        :comment: There is a low probably of an offset of 0,
            which would mean the day of death is the day seeded.
            This is to match EHG and is only mildly annoying.
        """
        logging.debug('ENTER: Person.HIVSeedInfect')
        #change to infected state
        self.is_infected = True
        duration = Disease._duration   #chk works but ugly
        # schedule current date uniformly during the duration of infection
        offset = self._prng.randint(0, duration)  #same as random.randrange
        doi = day - (duration - offset)
        #registry will register the disease with the scheduler
        self.disease = disease = Disease(self, doi, registry)
        for pship in self.partnerships:
            pship.expose_transmission(day)  #odd to use day instead of DOI, but simpler, and matches EHG
        logging.debug('EXIT: Person.HIVSeedInfect')
        return disease
    
    def has_primary(self):
        return any(p.tag=='primary' for p in self.partnerships)



    ###################### ai: convenience additions ##############################
    @property
    def sex(self):
        return self._sex
    @property
    def date_of_infection(self):
        return getattr(self.disease, 'start_date', None)

    @property
    def partners(self):
        return (indiv for pship in self.partnerships for indiv in pship.partners if indiv is not self)

  
class Partnership(object):
    """Provides a paired heterosexual "partnership" (i.e., sexual relationship).
    A `Partnership` has a predetermined random duration,
    and it schedules HIV transmission events.

    A Partnership is initialized with a `start_date`;
    its duration is determined during intialization,
    using a geometric distribution parameterized by `sigma`
    (which must be provided in the `params` dict).
    The mean duration (in days) of the partnership is simply 1/sigma.
    """
    def __init__(self, male, female, day, registry, params):
        """Return None. Initialize partnership,
        schedule partnership dissolution,
        schedule disease transmission.
        """
        logging.debug('ENTER: Partnership.__init__')
        assert (male.sex=='M' and female.sex=='F')
        self._mf = male, female
        self._prng = params['prng']
        self.set_tag(params)  #must come before duration setting
        self.start_date = day
        self.end_date = day + self.pshipduration(params)
        self.transmission_scheduled = False    # T/F whether HIV transmission is scheduled for the partnership
        self.transmission_date = None      #date for scheduled transmission
        #register this partnership if possible
        self.registry = registry  #can we get rid of this?
        if registry is not None:
            registry.register_partnership(self)
        # add new partnership to each partner's partnership list
        male.add_partnership(self)
        female.add_partnership(self)
        # expose partnership to infection transmission (if exactly one is infected)
        self.expose_transmission(day)
        logging.debug('EXIT: Partnership.__init__')
    def pshipduration(self, params):
        #we follow EHG, but is geometric distribution the best choice? (constant risk of break up...)
        #add partnership duration (min of 1)
        tag = self.tag
        if (tag == 'primary'):
            sigma = params['sigma01']
        elif (tag == 'secondary'):
            sigma = params['sigma02']
        else:
            raise ValueError("unknown partnership tag: {}".format(tag))
        return params['prng'].geometric(sigma)
    def set_tag(self,params):
        #allow only one primary partner (check union of partnerships sets)
        #(note that this is a partnership characteristic; not an individual attitude)
        male, female = self._mf
        has_primary = any(p.tag=='primary' for p in (male.partnerships + female.partnerships))
        if has_primary:
            self.tag = 'secondary'
            #set the relative frequency for non 'primary' partners (default is 1)
            self.sexfreq = params.get('secondarySexfreq', 1)
        else:
            self.tag = 'primary'
            self.sexfreq = 1  #primary partnership normal frequency
        
    def expose_transmission(self, day):
        """Return None; expose and possibly schedule partnership for transmission event.
        Called by `Partnership.__init__`, so upon formation a partnership is
        checked for HIV exposure.

        :comment:
            We only schedule a transmission if exactly one of the two is infected,
            but with concurrency another partner might still infect the uninfected partner
            before this scheduled transmission happens.
        """
        logging.debug('ENTER: Partnership.expose_transmission')
        male, female = self._mf

        if male.is_infected and not female.is_infected:
            infected = male
        elif not male.is_infected and female.is_infected:
            infected = female
        else:                #partners concordant -> exit the function
            return

        infected.disease.expose(self, day)  #expose the partnership to infection

        logging.debug('EXIT: Partnership.expose_transmission')

    def transmit(self, Disease, day): #compare EHG's HIVTransmission in partnership.cpp
        """Return tuple of (Person,Disease);
        transmit the HIV infection as **scheduled**
        from the infected partner to the uninfected partner.
        (If both are infected, then the originally uninfected partner got
        infected by someone else since this partnership formed.)
        """
        logging.debug('ENTER: Partnership.transmit_HIV')
        assert self.transmission_scheduled, 'Only scheduled transmissions shd call this'
        assert self.transmission_date==day, 'Only scheduled transmissions shd call this'
        registry = self.registry
        male, female = self._mf
        transmitter = None
        disease = None
        if(male.is_infected and not female.is_infected): # M --> F transmission
            transmitter = male
            disease = female.infect(Disease, day, registry)   # infect the female
        elif(not male.is_infected and female.is_infected): # F --> M transmission
            transmitter = female
            disease = male.infect(Disease, day, registry)   # infect the male
        else:
            #partner who was uninfected when partnership formed subsequently became infected
            #(in this case there is no new transmission -> no reinfection)
            assert (male.is_infected and female.is_infected)
        self.transmission_scheduled = False
        logging.debug('EXIT: Partnership.transmit_HIV')
        return transmitter, disease


    #properties
    @property
    def partners(self):
        return self._mf
        


