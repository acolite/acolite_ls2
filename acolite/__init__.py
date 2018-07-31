from acolite.shared import *
from acolite import ac
from acolite import acolite
from acolite import aerlut

from acolite import landsat
from acolite import sentinel

from acolite import output
from acolite import plotting

##
import os
path = os.path.dirname(__file__)
config = import_config(path+'/../config/config.txt')

## test whether we can find the relative paths
for t in config:
    if os.path.exists(config[t]): continue
    tmp = path + '/../' + config[t]
    config[t] = os.path.abspath(tmp)
