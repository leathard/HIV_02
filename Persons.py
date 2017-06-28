"""
This module is the place where all types of persons-participants in the simulation are defined.
"""

import logging
import copy


class Individual(object):
    """
    The purpose of this class is to define the interface for any individual.
    Add new behaviors here and override in the specific class that implements it.
    Must define upfront all expected behaviors of the class, present and future.
    """

    def description(self):
        return self.__class__.__name__

    def __init__(self):
        self._sex = NotImplemented
        self._prng = NotImplemented
        self.comsex = False  # commercial sex sales (F) or purchase (M)
        self.ART = False  # person receiving ART
        self.PrEP = False # person receiving PrEP
        
    def init(self, sex=NotImplemented, params=NotImplemented):
        """
        convenience routine that will be inherited, therefore can be called within constructors of derived classes
        """
        self._sex = sex  # 'M' or 'F' (vs EHG who use integer: 1=male, 2=female)
        self._prng = params['prng']
        self._Disease = params['Disease']
        self.comsex = False  # commercial sex sales (F) or purchase (M)
        self.ART = False  # person receiving ART
        self.PrEP = False # person receiving PrEP
        self.sexfreq = 1        
        self.reset()

    def reset(self):
        """Return None.
        Prepare a Person instance for reuse after death.
        """
        self.n_partners = 0
        self.partnerships = list()  # use list (not set) to ensure order of iteration
        self.is_infected = False  # ai: rendundant? (Just use doi)
        self.iDOD = None  # ai: impose max (adult) lifespan? chk (need day!)
        self.disease = None

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

    def infect(self, day, registry):  # compare EHG's HIVInfect
        """Return Disease instance,
        infect an individual with HIV, schedule death, and expose other partners to transmission 
        """
        logging.debug('ENTER: Person.infect')
        # change to infected state
        self.is_infected = True
        # fatal disease will schedule death (chkchk change to list?)
        # rc:disease = self.disease = Disease(self, day, registry=registry)
        Disease = self.disease = self.Disease(self, day, registry=registry)
        # expose current partners to disease transmission
        for pship in self.partnerships:
            pship.expose_transmission(day)
        logging.debug('EXIT: Person.infect')
        return Disease

    def HIVSeedInfect(self, day, registry):  # compare EHG's HIVSeedInfect in person.cpp
        """Return None.  seed infection
        :comment: There is a low probably of an offset of 0,
        which would mean the day of death is the day seeded.
        This is to match EHG and is only mildly annoying.
        """
        logging.debug('ENTER: Person.HIVSeedInfect')
        # change to infected state
        self.is_infected = True
        duration = self.Disease._duration  # chk works but ugly
        # schedule current date uniformly during the duration of infection
        offset = self._prng.randint(0, duration)  # same as random.randrange
        doi = day - (duration - offset)
        # registry will register the disease with the scheduler
        # rc: self.disease = disease = Disease(self, doi, registry)
        disease = self.disease = self.Disease(self, doi, registry)
        for pship in self.partnerships:
            pship.expose_transmission(day)  # odd to use day instead of DOI, but simpler, and matches EHG
        logging.debug('EXIT: Person.HIVSeedInfect')
        return disease

    def has_primary(self):
        return any(p.tag == 'primary' for p in self.partnerships)

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

    def is_available(self, **kwargs):
        """
        Returns availability for partnership based on values of keyword arguments.
        """
        return True

    @property
    def Disease(self):
        self.log_beta()
        return self._Disease

    @Disease.setter
    def Disease(self, value):
        self.log_beta()
        self._disease = value

    def print_beta(self):
        if self.sex=='M':
            print("beta:%s"%str(self._Disease._beta_M2F))
        if self.sex=='F':
            print("beta:%s"%str(self._Disease._beta_F2M))

    def log_beta(self):
        if self.sex=='M':
            logging.debug("beta:%s"%str(self._Disease._beta_M2F))
        if self.sex=='F':
            logging.debug("beta:%s"%str(self._Disease._beta_F2M))


class Decorator(Individual):
    """
    The purpose of this class is to allow building more complex types and have access to its building blocks types behavior.
    """

    def __init__(self, individual, **kwargs):
        self.person = individual
        #self.xparams = kwargs

    def params(self):
        return self.xparams

    def description(self):
        return self.__class__.__name__ + '(' + self.person.description() + ')'

    def is_available(self, **kwargs):
        #looks for reasons the person might not be available (worst case)
        return True and self.person.is_available(**kwargs)


class Person(Decorator):
    """Provides a basic person class.
    A person has a sex and a set of partnerships.
    A person may have a disease and an associated date of death.
    :todo: add date of birth for age contingent outcomes
    """


    def __init__(self, individual, sex=None, params=None):
        """Return an initialized Person.
        """
        Decorator.__init__(self, individual, sex=sex, params=params)        
        self.init(sex=sex, params=params)
        self.comsex |= self.person.comsex
        self.ART    |= self.person.ART
        self.PrEP   |= self.person.PrEP
        
        if self.person.ART:
            logging.debug("Overriding ART "+self.__class__.__name__)
            self._Disease = copy.deepcopy(self.person._Disease)
        if self.person.PrEP:
            logging.debug("Overriding PrEP "+self.__class__.__name__)
            self._Disease = copy.deepcopy(self.person._Disease)
            self.logging.debug_beta()

class Person01(Decorator):
    """
    Person01 extends Person class for sexworkers(F). It is also used for sexworkers clients(M).
    """

    def __init__(self, individual, sex=None, params=None):
        """
        Looks inside the decorated object to check if it should receive specific treatment (thus overwriting _Disease)
        If so, self._Disease is overwritten with either beta_ART or beta_PREP, depending on the case.
        """
        Decorator.__init__(self, individual, sex=sex, params=params)
        self.init(sex=sex, params=params)

        self.comsex  = True
        self.sexfreq = 5 if sex=='F' else 1
        self.ART    |= self.person.ART
        self.PrEP   |= self.person.PrEP
        if self.person.ART:
            logging.debug("Overriding ART "+self.__class__.__name__)
            self._Disease = copy.deepcopy(self.person._Disease)
            self.log_beta()
        if self.person.PrEP:
            logging.debug("Overriding PrEP "+self.__class__.__name__)
            self._Disease = copy.deepcopy(self.person._Disease )
            self.log_beta()

class Person02(Decorator):
    """
	Person02 extends Person class
	new attribute: ART (meaning receiving ART treatment)
	"""
    
    def __init__(self, individual, sex=None, params=None):
        """
        ART and PREP are conflicting properties. beta_PREP overwrites beta_ART if both treatments are set.
        """
        Decorator.__init__(self, individual, sex=sex, params=params)
        self.init(sex=sex, params=params)
        self.comsex |= self.person.comsex
        self.ART     = True
        self.PrEP   |= self.person.PrEP
        if self.person.PrEP:
            logging.debug("Overriding PrEP "+self.__class__.__name__)
            self._Disease = copy.deepcopy(self.person._Disease )
            self.log_beta()
        
class Person03(Decorator):
    """
    This class represents the miner class, with higher sexfreq (3) only active ~ 200 days out of 365/year
    static data members:
    activeRange  = (0,200) : the range of days throughout the year when they can form partnerships
    sexfreq      = 3 : sex frequency parameter, native to the individual type, for all instances
    The type Person03 is used in scheduler to discriminate behavior in forming parnerships
    using a specific override is_available(yearDay=yearday)
    """

    def __init__(self, individual, sex=None, params=None):
        Decorator.__init__(self, individual, sex=sex, params=params)
        self.activeRange = (0, 200)
        self.init(sex=sex, params=params)
        self.sexfreq = 3
        self.comsex |= self.person.comsex
        self.ART    |= self.person.ART
        self.PrEP   |= self.person.PrEP
        if self.person.ART:
            logging.debug("Overriding ART "+self.__class__.__name__)
            self._Disease = copy.deepcopy(self.person._Disease)
            self.log_beta()
        if self.person.PrEP:
            logging.debug("Overriding PrEP "+self.__class__.__name__)
            self._Disease = copy.deepcopy(self.person._Disease)
            self.log_beta()
    def is_available(self, yearDay=None):
        return yearDay >= self.activeRange[0] and yearDay < self.activeRange[-1]


class Person04(Decorator):
    """
    This class models individuals taking PrEP drugs.
    Its only influence on the simulation is the beta parameters (set to 0,0,0,0).
    """

    def __init__(self, individual, sex=None, params=None):
        Decorator.__init__(self, individual, sex=sex, params=params)
        self.init(sex=sex, params=params)
        self.comsex |= self.person.comsex
        self.PrEP    = True
        self.ART    |= self.person.ART

def typesCombinations(origAttrsSet, numbIndividuals):
    """
    Engine to generate composed types from the size of the population and its percent-wise composition.
    input: a list of tuples, population size
    list of tuples: each containing (type "factory", construction arguments for type, {'fraction':percent in population}
    output: a population composed from subpopulations as prescibed by their fraction (or as close as possible).

    It generates all combinations of attributes in decreasing order of their complexity:
      - for [A,B,C,D] 
        - it generates (ABCD), then (ABC), (ABD), (BCD), (ABD), then (AB),(AC),(AD),...(CD)
           - then creates an mixed type which on the surface exhibits the characteritics of the outermost(rightmost) type
           - the decorator structure alows the inner types to be accessed
    classes with methods overriding Person methods should be placed close to the end of the list-of-tuples
    classes with parameters that needs overwriting towards encapsulating type should be placed near the beginning of the list.
    """
    import itertools
    population = []
    attrsSet = copy.deepcopy(origAttrsSet)
    for n in xrange(len(attrsSet), 0, -1):
        for combination in itertools.combinations(attrsSet, n):
            fraction = 1.0
            for Type, kwArgs, buildParams in combination:
                fraction *= buildParams['fraction']
            nComposedType=int(round(fraction * numbIndividuals))
            correction=fraction-nComposedType/float(numbIndividuals)
            for i in xrange(nComposedType):
                person = Individual()
                for Type, kwArgs, buildParams in combination:
                    person = Type(person, **copy.deepcopy(kwArgs))
                population.append(person)
            for Type, kwArgs, buildParams in combination:
                buildParams['fraction'] -= (fraction-correction)

    nleft = numbIndividuals - len(population)
    # and correcting for accumulated rounding errors
    # fill it with "Person" neutral type
    population += [attrsSet[0][0](Individual(), **attrsSet[0][1]) for i in xrange(nleft)]
    return population


def isPerson(instance,Type):
    if type(instance)==Individual:
        return False
    return isinstance(instance,Type) or isPerson(instance.person,Type)

def isPersonInfectedInStage(instance,Type,stage,day):
    return isPerson(instance,Type) and instance.is_infected and instance.get_HIV_stage(day)==stage
    
if __name__ == '__main__':
    """nM=100000
    p_miners=0.08;p_nM_ART=0.2;p_PREP=0.1;p_nclients=0.087
    class disease:
        def __init__(self,name):
            self._beta_M2F=1.2
            self._beta_F2M=1.0
            self.name=name
    params={'Disease':disease("disease"),'prng':'prng'}
    params_ART={'Disease':disease("disease_ART"),'prng':'prng'}
    prep_params={'Disease':disease("disease_prep"),'prng':'prng'}
    
    mPopComb = [(Person, {'sex': 'M', 'params': params}, {'fraction': 1 - p_PREP - p_miners - p_nM_ART - p_nclients})]
    mPopComb +=[(Person01,{'sex':'M','params':params},{'fraction':p_nclients})]
    mPopComb +=[(Person02, {'sex': 'M', 'params': params_ART}, {'fraction': p_nM_ART})]
    mPopComb +=[(Person04, {'sex': 'M', 'params': prep_params}, {'fraction': p_PREP})]
    mPopComb +=[(Person03, {'sex': 'M', 'params': params}, {'fraction': p_miners})]
    males = typesCombinations(copy.deepcopy(mPopComb), nM)
    n_ART=nM*p_nM_ART
    n_PREP=nM*p_PREP
    n_comsex=nM*p_nclients
    print
    c_ART=len(filter(lambda x:x.ART,males))
    c_PREP=len(filter(lambda x:x.PrEP,males))
    c_CS=len(filter(lambda x:x.comsex,males))
    print "comsex count known: %d, counted: %d"%(n_comsex,c_CS),(c_CS-n_comsex)/float(n_comsex)
    print "prep count  known: %d, counted: %d"%(n_PREP,c_PREP),(c_PREP-n_PREP)/float(n_PREP)
    print "art count  known: %d, counted: %d"%(n_ART,c_ART),(c_ART-n_ART)/float(n_ART)"""
    #Persons = [Person01, Person02, Person03, Person04] # what about Person?
    #print map(lambda x:x.__name__,Persons)
