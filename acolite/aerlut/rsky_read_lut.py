## import sky reflectance lut
## QV 2020-03-17

def rsky_read_lut(model, lutbase='ACOLITE-RSKY-202003'):
    import os
    import numpy as np
    import scipy.interpolate
    from netCDF4 import Dataset
    import acolite as ac

    lutdir = '{}/LUT/RSKY/'.format(ac.config['pp_data_dir'])
    lutnc = '{}/{}-MOD{}.nc'.format(lutdir, lutbase, model)

    if os.path.isfile(lutnc):
            nc = Dataset(lutnc)
            meta = {}
            for attr in nc.ncattrs():
                attdata = getattr(nc,attr)
                if isinstance(attdata,str): attdata = attdata.split(',')
                meta[attr]=attdata
            lut = nc.variables['lut'][:]
            nc.close()

            dim = [meta['wave'], meta['azi'], meta['thv'], meta['ths'], meta['tau']]
            rgi = scipy.interpolate.RegularGridInterpolator(dim, lut, bounds_error=False, fill_value=np.nan)
            return(lut, meta, dim, rgi)
