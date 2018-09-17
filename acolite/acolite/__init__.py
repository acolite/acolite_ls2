from .acolite_ac import acolite_ac
from .acolite_settings import acolite_settings

from .acolite_run import acolite_run
from .acolite_gui import acolite_gui
from .acolite_cli import acolite_cli

from .acolite_toa_crop import acolite_toa_crop

from .acolite_l2w import acolite_l2w
from .acolite_l2w_qaa import acolite_l2w_qaa
from .l2w_required import l2w_required

from .acolite_map import acolite_map

from .pscale import pscale
from .settings_read import settings_read
from .settings_write import settings_write

from acolite.shared import *

import os
path=os.path.dirname(os.path.abspath(__file__))
for i in range(2): path = os.path.split(path)[0]

## check if binary distribution
if '{}dist{}acolite'.format(os.path.sep,os.path.sep) in path:
    ## two more levels for this file
    for i in range(2): path = os.path.split(path)[0]

cfile='{}{}config{}acolite_config.txt'.format(path,os.path.sep,os.path.sep)
config = import_config(cfile)

## test whether we can find the relative paths
for t in config:
    if t in ['version']: continue
    if os.path.exists(config[t]): continue
    tmp = path + os.path.sep + config[t]
    config[t] = os.path.abspath(tmp)
