## def get_sensor_lut
## imports LUT, interpolates to hyperspectral and convolutes to sensor bands
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2016-07-13
## modifications: 2017-01-23 (QV) added generic lutdir/rsrfile locations
##                2018-06-07 (QV) reordered lutdir and rsr_file check
##                2018-07-18 (QV) changed acolite import name
def get_sensor_lut(sensor, rsr_file, override=0, lutdir=None, lutid='PONDER-LUT-201605-MOD2-1013mb'):
    import os, sys


    if True:
        import acolite as pp
        from acolite import rsr_convolute, rsr_read, aerlut
        ## get sensor RSR
        pp_path = pp.config['pp_data_dir']
        if lutdir is None:
            lutdir=pp_path+'/LUT/'
        if rsr_file is None:
            rsr_file = pp_path+'/RSR/'+sensor+'.txt'
        rsr, rsr_bands = rsr_read(file=rsr_file)

    ## sensor LUT NetCDF is stored here
    lutnc=lutdir+'/'+lutid+'/'+lutid+'_'+sensor+'.nc'
    if (os.path.isfile(lutnc)) & (override == 1): os.remove(lutnc)

    ## generate sensor LUT NetCDF if needed
    if (os.path.isfile(lutnc) == 0): 
        from numpy import linspace, interp, nan, zeros

        ## read LUT
        lut, meta = aerlut.import_lut(lutid,lutdir, override=0)
        lut_dims = lut.shape

        ## set up wavelength space
        wave_range=[0.2,2.4]
        wave_step=0.001
        wave_hyper = linspace(wave_range[0],wave_range[1],((wave_range[1]-wave_range[0])/wave_step)+2)

        ## interpolate RSR to same dimensions
        rsr_hyper = dict()
        for band in rsr: 
            band_wave_hyper = wave_hyper
            band_response_hyper = interp(wave_hyper, rsr[band]['wave'], rsr[band]['response'], left=0, right=0)
            band_response_sum = sum(band_response_hyper)
            rsr_hyper[band]={'wave':band_wave_hyper, 'response': band_response_hyper, 'sum':band_response_sum}

        ## interpolate lut to hyperspectral, and convolute to sensor bands        
        ### set up empty sensor LUT
        lut_sensor = dict()
        for band in rsr: lut_sensor[band]=zeros([len(meta['par']),len(meta['azi']),len(meta['thv']),len(meta['ths']), 1, len(meta['tau'])])

        ### only one wind setting, hard coded here
        i4=0
        v4=meta['wnd']

        ### run through all simulations
        data_wave = meta['wave']
        for i0, v0 in enumerate(meta['par']):
            for i1, v1 in enumerate(meta['azi']):
                for i2, v2 in enumerate(meta['thv']):
                    for i3, v3 in enumerate(meta['ths']):
                        #for i4, v4 in enumerate(meta['wnd']):
                            for i5, v5 in enumerate(meta['tau']):
                                data_hyper = interp(wave_hyper, data_wave, lut[i0,:,i1,i2,i3,i4,i5], left=0, right=0)
                                for band in rsr: 
                                    lut_sensor[band][i0,i1,i2,i3,i4,i5] = (sum(data_hyper*rsr_hyper[band]['response'])/rsr_hyper[band]['sum'])

        ## write nc file
        try:
            if os.path.isfile(lutnc) == 0:
                from netCDF4 import Dataset
                nc = Dataset(lutnc, 'w', format='NETCDF4_CLASSIC')
                ## write metadata
                for i in meta: 
                    attdata=meta[i]
                    if isinstance(attdata,list):
                        if isinstance(attdata[0],str): 
                            attdata=','.join(attdata)
                    setattr(nc, i, attdata)
                ## set up LUT dimension
                nc.createDimension('par', lut_dims[0])
                #nc.createDimension('wave', lut_dims[1]) # not used here
                nc.createDimension('azi', lut_dims[2])
                nc.createDimension('thv', lut_dims[3])
                nc.createDimension('ths', lut_dims[4])
                nc.createDimension('wnd', lut_dims[5])
                nc.createDimension('tau', lut_dims[6])
                ## write LUT
                for band in lut_sensor.keys():
                    var = nc.createVariable(band,float,('par','azi','thv','ths','wnd','tau'))
                    nc.variables[band][:] = lut_sensor[band]
                nc.close()
                nc = None
                arr = None
                meta = None
        except:
            if os.path.isfile(lutnc): os.remove(lutnc)
            print(sys.exc_info()[0])
            print('Failed to write LUT data to NetCDF (id='+lutid+')')
  
    ## read dataset from NetCDF
    if (os.path.isfile(lutnc) == 1): 
        try:
            if os.path.isfile(lutnc):
                from netCDF4 import Dataset
                nc = Dataset(lutnc)
                ## read in metadata
                meta=dict()
                for attr in nc.ncattrs():
                    attdata = getattr(nc,attr)
                    if isinstance(attdata,str): attdata = attdata.split(',')
                    meta[attr]=attdata
                ## read in LUT
                lut_sensor = dict()
                datasets = list(nc.variables.keys())
                for dataset in datasets: 
                    lut_sensor[dataset] = nc.variables[dataset][:]
                nc.close()
                nc = None
        except:
            print(sys.exc_info()[0])
            print('Failed to open LUT data from NetCDF (id='+lutid+')')

    return(lut_sensor, meta)
