## o2lut_get
## imports o2 transmittance LUT
##
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2018-10-30
## modifications:

def o2lut_get(config = "201810C"):
    import acolite as ac
    
    from time import time, strftime
    import os, sys
    
    lut_path = '{}/LUT/O2'.format(ac.config['pp_data_dir'])
    lut_id = 'O2_{}'.format(config)
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
        print(lut_nc)

    return(lut,meta)
