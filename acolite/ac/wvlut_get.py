## wvlut_get
## imports WV transmittance LUT
##
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-10-17
## modifications:
##                2017-11-28 (QV) moved PP data directory
##                2018-07-18 (QV) changed acolite import name
def wvlut_get(config = "201710C"):
    import acolite as pp
    
    from time import time, strftime
    import os, sys
    
    #from source import config as pp_config
    pp_path = pp.config['pp_data_dir']

    lut_path = '{}/LUT/WV'.format(pp_path)
    lut_id = 'WV_{}'.format(config)
    lutnc = '{}/{}.nc'.format(lut_path,lut_id)

    ## read dataset from NetCDF
    try:
        #if os.path.isfile(lutnc):
            from netCDF4 import Dataset
            nc = Dataset(lutnc)
            meta=dict()
            for attr in nc.ncattrs():
                attdata = getattr(nc,attr)
                if isinstance(attdata,str): attdata = attdata.split(',')
                meta[attr]=attdata
            lut = nc.variables['lut'][:]
            nc.close()
    except:
        print(sys.exc_info()[0])
        print('Failed to open LUT {} data from NetCDF'.format(config))
        
    return(lut,meta)
