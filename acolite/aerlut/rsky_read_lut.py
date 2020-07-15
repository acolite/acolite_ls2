## import sky reflectance lut
## QV 2020-03-17
## last updates: 2020-07-15 (QV) added sensor resampling

def rsky_read_lut(model, lutbase='ACOLITE-RSKY-202003', sensor=None, override=False):
    import os
    import numpy as np
    import scipy.interpolate
    from netCDF4 import Dataset
    import acolite as ac

    lutdir = '{}/LUT/RSKY/'.format(ac.config['pp_data_dir'])
    lutnc = '{}/{}-MOD{}.nc'.format(lutdir, lutbase, model)

    if os.path.isfile(lutnc):
            ## base lut
            if sensor is None:
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

            ## sensor specific lut
            else:
                lutnc_s = '{}/{}-MOD{}_{}.nc'.format(lutdir, lutbase, model, sensor)

                ## get sensor RSR
                pp_path = ac.config['pp_data_dir']
                lutdir=pp_path+'/LUT/'
                rsr_file = pp_path+'/RSR/'+sensor+'.txt'
                rsr, rsr_bands = ac.shared.rsr_read(file=rsr_file)


                ## make new sensor lutfile
                if (override) & (os.path.isfile(lutnc_s)): os.remove(lutnc_s)

                if not os.path.isfile(lutnc_s):
                    print('Resampling {} to {}'.format(model, sensor))

                    ## read lut
                    lut, meta, dim, rgi = ac.aerlut.rsky_read_lut(model)

                    ## resample to bands
                    lut_sensor = {}
                    for band in rsr_bands:
                        lut_sensor[band] = ac.shared.rsr_convolute_nd(lut, meta['wave'],rsr[band]['response'], rsr[band]['wave'], axis=0)

                    #return(lut_sensor, meta, dim)
                    ## save to new file
                    from netCDF4 import Dataset
                    nc = Dataset(lutnc_s, 'w', format='NETCDF4_CLASSIC')
                    for i in meta.keys():
                        attdata=meta[i]
                        if isinstance(attdata,list):
                            if isinstance(attdata[0],str):
                                attdata=','.join(attdata)
                        setattr(nc, i, attdata)

                    #return(dim[1])
                    nc.createDimension('azi', len(dim[1]))
                    nc.createDimension('thv', len(dim[2]))
                    nc.createDimension('ths', len(dim[3]))
                    nc.createDimension('tau', len(dim[4]))

                    ## write LUT
                    for band in rsr_bands:
                        var = nc.createVariable(band,float,('azi','thv','ths','tau'))
                        nc.variables[band][:] = lut_sensor[band]
                    nc.close()
                ## end resample lut

                ## lutfile was already resampled
                if os.path.isfile(lutnc_s):
                    nc = Dataset(lutnc_s)
                    meta = {}
                    for attr in nc.ncattrs():
                        attdata = getattr(nc,attr)
                        if isinstance(attdata,str): attdata = attdata.split(',')
                        meta[attr]=attdata

                    lut_sensor = {}
                    for band in rsr_bands:
                        lut_sensor[band] = nc.variables[band][:]
                    nc.close()


                    dim = [meta['azi'], meta['thv'], meta['ths'], meta['tau']]
                    rgi_sensor = {}
                    for band in rsr_bands:
                        rgi_sensor[band] = scipy.interpolate.RegularGridInterpolator(dim, lut_sensor[band], bounds_error=False, fill_value=np.nan)
                    return(lut_sensor, meta, dim, rgi_sensor)
