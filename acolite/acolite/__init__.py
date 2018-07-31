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

from acolite import shared

import os
path=os.path.dirname(os.path.abspath(__file__))
config=shared.import_config(path+'/../../config/acolite_config.txt')

for t in config: config[t] = path + '/../../' + config[t]
