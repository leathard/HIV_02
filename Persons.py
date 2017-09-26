""" This module is the place where all types of persons-participants in the
simulation are defined.
Individual object is refactored Person object from concurrency.py written by
Alan G. Isaac. 

"""

import logging
import copy
import itertools 
import concurrency 
import Scheduler 


class Persons(concurrency.Person):
    """ Person(concurrency.Person) builds on isaac.concurrency.Person class.

    This class is a basic person class used to build Person* subclasses.
    It holds common objects (attribues and methods) for all Person* subclasses.
    The __init__ method of Person class will be executed by each Person*
    subclass, neccesitating execution of only Person* changes (which are
    minimal) at each Person* __init__.

    Vars: 

    params: 
    """

    def __init__(self, ptype, params, person):
        """Return initialized Person.
        """
        concurrency.Person.__init__(self, sex, registry, params)        
        #self._sex = person._sex 
        #self._prng = person._prng
        self.comsex = person.comsex 
        #self._type = ptype
        #self._Disease = params['Disease']
        self.ART = False
        self.PrEP = False
        self.sexfreq = 1 
        #self.params = params 
        self.activeRange = (0, 365)
        if ptype == 'PersonCSW':
            self.comsex = True
            if self._sex=='F':
                self.sexfreq = 5 
        if ptype == 'PersonART':
            self.ART = True 
        if ptype == 'PersonMiner':
            self.activeRange = (0,200)
            self.sexfreq = 3
        if ptype == 'PersonPREP':
            self.PrEP = True
        #self.n_partners = person.n_partners
        #self.partnerships = person.partnerships
        #self.is_infected = person.is_infected
        #self.disease = person.disease
        #self.iDOD= person.iDOD

    def infect(self, day, registry): 
        """Return Disease instance,
        infect an individual with HIV, schedule death, and expose other partners to transmission 
        """
        logging.debug('ENTER: Person.infect')
        self.is_infected = True
        disease = self.disease = self.Disease(self, day, registry=registry)
        for pship in self.partnerships:
            pship.expose_transmission(day)
        logging.debug('EXIT: Person.infect')
        return disease

    def HIVSeedInfect(self, day, registry):  #compare EHG's HIVSeedInfect in person.cpp - tr: Disease,  removed
        """Return None.  seed infection
        :comment: There is a low probably of an offset of 0,
            which would mean the day of death is the day seeded.
            This is to match EHG and is only mildly annoying.
        """
        logging.debug('ENTER: Person.HIVSeedInfect')
        self.is_infected = True
        duration = self.Disease._duration   
        offset = self._prng.randint(0, duration)  
        doi = day - (duration - offset)
        disease = self.disease = self.Disease(self, doi, registry)  
        for pship in self.partnerships:
            pship.expose_transmission(day)
        logging.debug('EXIT: Person.HIVSeedInfect')
        return disease

    def is_available(self, yearDay=None):
        return yearDay >= self.activeRange[0] and yearDay < self.activeRange[-1]

    @property
    def Disease(self):
        return self._Disease 

def typesCombinations(origAttrsSet, numbIndividuals):
    """
    Engine to generate composed types from the size of the population and its
    percent-wise composition.  input: a list of tuples, population size list of
    tuples: each containing (type "factory", construction arguments for type,
    {'fraction':percent in population} output: a population composed from
    subpopulations as prescibed by their fraction (or as close as possible).
    """
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
                person = Person(**copy.deepcopy(kwArgs))
                for Type, kwArgs, buildParams in combination:
                    #person = Type(person, **copy.deepcopy(kwArgs))
                    person = PersonJ(Type, kwArgs['params'], person)    #**copy.deepcopy(kwArgs)
                population.append(person)
            for Type, kwArgs, buildParams in combination:
                buildParams['fraction'] -= (fraction-correction)

    nleft = numbIndividuals - len(population)
    population += [Persons(attrsSet[0][0], attrsSet[0][1]['params'], Person(**attrsSet[0][1])) for i in xrange(nleft)]    
    return population


def isPerson(instance, Type):
    """Add docstring
    """
    #if type(instance) == Person:
    #    return False
    #return isinstance(instance, Type) or isPerson(instance.person, Type)
    return instance._type == Type

def isPersonInfectedInStage(instance, Type, stage, day):
    """Add docstring
    """
    return isPerson(instance, Type) and instance.is_infected and instance.get_HIV_stage(day) == stage

if __name__ == "__main__":
    """just for convenience
    """


