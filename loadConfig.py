import ConfigParser
import StringIO

from phi import *

def loadConfig(ini_path,defsection=None):
    """
    loads a config without sections (if defsection=None) or `defsection` section otherwise.
    returns a params dictionary with key-values mirroring the ones in config/section.
    """
    section='__root__' if defsection is None else defsection
    ini_str = '[%s]\n'%section + open(ini_path, 'r').read()
    ini_fp = StringIO.StringIO(ini_str)
    config = ConfigParser.ConfigParser()
    config.optionxform=str
    config.readfp(ini_fp)
    rParams=dict(config.items(section))
    if len(rParams)==0:
        raise ConfigParser.NoSectionError("section %s does not exist or is empty"%section)
    item=None;value=None
    for k in rParams.keys():
        item=k; value=rParams[k]
        try:
            rParams[k]=eval(rParams[k])
        except Exception as inst:
            extra=" config item:%s value:%s "%(k,value)            
            inst.args =inst.args+(extra,)
            inst.message += extra
            raise
    return rParams

def updateComputableParts(params):
    """
    updates the params in received as arguments according to computations needed and also returns the updated value
    """
    import math
    nOutGroups=params['nOutGroups']
    params['epsilon'] = tuple(1 - 1.02 * (math.exp(3.9 * (xi / 10. - 1)) - .02) for xi in range(nOutGroups - 1))

    transRatesKeys=filter(lambda x:'transmissionRatesHIV' in x,params.keys())
    r=params['relativeTransmissionRatesHIV']
    l=len(r)
    for k in transRatesKeys:        
        keySuffix=k.replace('transmissionRatesHIV','')        
        beta_M2F=tuple(2*params[k][i]/(1.0+r[i]) for i in range(l))
        beta_F2M=tuple(2*params[k][i]*r[i]/(1.0+r[i]) for i in range(l))        
        params['beta_M2F'+keySuffix]=beta_M2F
        params['beta_F2M'+keySuffix]=beta_F2M
    params['model_phi']=params['sPhi']
    params['outfile_pattern'] = 'eps??sim??.out'
    
    return params

if __name__=='__main__':
    baseline = 'baseline.ini'
    
    rParams=loadConfig(baseline)
    experimentf='experiments.ini'
    experiment='FSWlow'
    rParams.update(loadConfig(experimentf,defsection=experiment))
    updateComputableParts(rParams)
    print rParams
