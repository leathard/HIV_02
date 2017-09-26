from isaac import concurrency # isaac/concurrency.py


class PersonExample01(concurrency.Person):
    def __init__(self,
        sex,
        comsex,    #commercial sex worker?
        treatment, #treament 0 (none), 1 (art), 2 (prep)
        registry,
        params
        ):
        concurrency.Person.__init__(self, sex, registry, params)
        self.comsex = comsex
        self.treament = treatment
    def __str__(self):
        return """
        I am a {sex} {cls} who is {comsex} engaged in sex work.
        """.format(**dict(
            sex=self.sex,
            cls=self.__class__,
            comsex="not" if self.sex=="M" else ""
            ))

    @property
    def ART(self):
        return 1 == self.treatment

    @property
    def PrEP(self):
        return 2 == self.treatment

    def is_available(self, yearDay=None):
        #Probably better to just pass the day and compute the yearDay here ...
        return yearDay >= self.activeRange[0] and yearDay < self.activeRange[-1]


if __name__=="__main__":
    p = PersonExample('M', True, 0, None, dict(prng=None))
    print(p)

