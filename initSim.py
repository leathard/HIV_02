"""
Global initialization of simulation.
Creates directories and initializes logging
"""

import logging
import datetime
import os


def tryMkDir(path):
    try:
        os.makedirs(path)
        return True
    except OSError:
        return False
        pass


def prepare(outputPath, logPath='out/temp.log'):
    import os
    if not tryMkDir(os.path.dirname(logPath)):
        print "couldn't create directories: %s"%os.path.dirname(logPath)
    logging.basicConfig(level=logging.ERROR, filename=logPath, filemode='w')
    logging.info(datetime.datetime.today().isoformat())
    tryMkDir(outputPath)
