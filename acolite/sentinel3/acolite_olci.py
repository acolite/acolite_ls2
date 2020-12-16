## def acolite_olci
## does DSF atmospheric correction for Sentinel-3/OLCI
## written by Quinten Vanhellemont, RBINS
## 2019-04-02
## last updates QV 2019-04-03 added pressure dependency, sky correction
##              QV 2019-04-04 added resolved geometry for LUT
##              QV 2019-04-30 new LUT interpolation
##              QV 2019-05-02 changed TPG resampling
##              QV 2020-07-23 new tiling and new LUTs and new LUT interpolation
##              QV 2020-07-24 added smile correction
##              QV 2020-08-03 added filtering of tiled tau
##              QV 2020-11-23 set suggested defaults (36x36 km tiles with darkest pixel)
##              QV 2020-12-15 added a simple (non-water) mask output
##              QV 2020-12-16 added support for user defined gains

def acolite_olci(bundle, output, limit=None,
                use_tpg=False,
                luts = ['ACOLITE-LUT-202003-MOD1', 'ACOLITE-LUT-202003-MOD2'],
                rhod_tgas_cutoff = 0.85, rhod_min_wave = 440.,
                map_rgb=True,

                tiled_processing = True,
                tile_dimensions = (36000, 36000), ## in m # (12000,12000)
                sensor_resolution = 300,
                dark_spectrum_option = 'darkest', # intercept
                dark_spectrum_percentile = 1,
                dif_tau_max = 0.02,

                verbosity=0,
                uoz_default=0.3, uwv_default=1.5,
                sky_correction=True, gas_transmittance=True,
                resolved_geometry=True, write_resolved_atmosphere = False,
                use_supplied_ancillary = True,

                rhot_mask_threshold = 0.4,
                rhot_mask_threshold_swir = 0.03,

                ##
                smile_correction = True,
                tiled_interpolation = 'tpg',

                ## allow for user specified gains
                use_gains = False,
                gains = [1.0]*21,

                ## these are not needed?
                use_supplied_pressure = True,
                use_supplied_altitude = True,

                ## these are not needed?
                elevation = None,
                pressure = 1013.25,
                dem_pressure_percentile = 25):

    import os, glob
    import scipy.interpolate
    import scipy.ndimage

    import numpy as np
    import acolite as ac
    import dateutil.parser

    import matplotlib
    #matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    from acolite.sentinel3 import olci_sub

    import datetime
    print('{} - Running ACOLITE/OLCI processing for {}'.format(datetime.datetime.now().isoformat()[0:19], bundle))

    dtor = np.pi/180

    ## identify sensor
    sensor = '{}_OLCI'.format(os.path.basename(bundle)[0:3])
    rsr_file = ac.config['pp_data_dir']+'/RSR/'+sensor+'.txt'
    rsr, rsr_bands = ac.shared.rsr_read(file=rsr_file)

    ## read luts
    lutdir= '{}/LUT/'.format(ac.config['pp_data_dir'])
    lutd = ac.aerlut.import_luts(base_luts=luts, sensor=sensor, add_rsky=True)

    ## set up
    new = True
    oname = os.path.splitext(os.path.basename(bundle))[0]
    ofile = '{}/{}.nc'.format(output, oname)

    dfiles = [os.path.basename(f) for f in glob.glob('{}/*.nc'.format(bundle))]
    dfiles.sort()

    mfile = [os.path.basename(f) for f in glob.glob('{}/*.xml'.format(bundle))]
    if len(mfile)==1: mfile = mfile[0]

    sub = None
    data_shape = None
    if limit is not None:
        sub, data_shape = olci_sub(bundle, limit, use_tpg=use_tpg)
        if any([s == 0 for s in sub[2:3]]):
            print('Limit {} out of scene {}'.format(limit, bundle))
            print(sub)
            return()

    ## read data
    dshape = None
    lfiles = {}
    data, meta = {}, {}
    for f in dfiles:
        fname = os.path.splitext(f)[0]
        file = '{}/{}'.format(bundle, f)
        datasets = ac.shared.nc_datasets(file)
        gatts = ac.shared.nc_gatts(file)
        start_time = gatts['start_time']
        stop_time = gatts['stop_time']

        ## tie point grids
        if ('tie' in fname) or (fname in ['removed_pixels',
                     'instrument_data', 'time_coordinates']):
            meta[fname]={}
            for ds in datasets:
                meta[fname][ds]= ac.shared.nc_data(file, ds)

                if ds == 'detector_index':
                    data[ds] = ac.shared.nc_data(file, ds, sub=sub)
        ## or full size data
        else:
            for ds in datasets:
                ## also save toa reflectance here?
                if '_radiance' in ds:
                    ## store full paths to read data later
                    lfiles[ds] = file
                    if smile_correction:
                        data[ds] = ac.shared.nc_data(file, ds, sub=sub)
                        if use_gains:
                            cg = 1.0
                            if len(gains) == 21:
                                gi = int(ds[2:4])-1
                                cg = float(gains[gi])
                            print('Applying gain {:.5f} for {}'.format(cg, ds))
                            data[ds]*=cg
                    continue
                ## end skip radiance

                data_ = ac.shared.nc_data(file, ds, sub=sub)
                if dshape is None: dshape = data_.shape
                if data_shape is None: data_shape = data_.shape

                ## save latitude and longitude already
                if ds in ['latitude', 'longitude']:
                    ## subset to lat/lon for compatibility with other scripts
                    ac.output.nc_write(ofile, ds[0:3], data_, new=new)
                else:
                    ac.output.nc_write(ofile, ds, data_, new=new)
                    del data_
                    new=False
                ## end lat lon
    ##
    full_scene = False
    if sub is None:
        #data_shape = data['latitude'].shape
        full_scene = True
        ## reversed 2019-04-03
        #subx = np.arange(0,data_shape[0])
        #suby = np.arange(0,data_shape[1])
        ## reversed again 2019-05-02
        ## added half a pixel offset - correspond to centre (and SNAP)
        subx = np.arange(0,data_shape[1])+0.5
        suby = np.arange(0,data_shape[0])+0.5
    else:
        subx = np.arange(sub[1],sub[1]+sub[3])
        suby = np.arange(sub[0],sub[0]+sub[2])
        ## added half a pixel offset - correspond to centre (and SNAP)
        subx = np.arange(sub[0],sub[0]+sub[2])+0.5
        suby = np.arange(sub[1],sub[1]+sub[3])+0.5

    ## interpolate tie point grids
    tpg_shape = meta['tie_geo_coordinates']['latitude'].shape

    tpx = np.linspace(0,data_shape[1]-1,tpg_shape[1])
    tpy = np.linspace(0,data_shape[0]-1,tpg_shape[0])
    #idx = np.linspace(0,data_shape[1]-1,data_shape[1])

    tpg = {}
    for k in meta.keys():
        if 'tie_' in k:
            for l in meta[k].keys():
                if l in ['atmospheric_temperature_profile',
                         'horizontal_wind', 'reference_pressure_level']: continue
                z = scipy.interpolate.interp2d(tpx, tpy, meta[k][l])
                tpg[l] = z(subx,suby)

    ## compute relative azimuth TPG
    tpg['RAA'] = abs(tpg['SAA']-tpg['OAA'])
    tpg['RAA'][tpg['RAA']>180]-=180

    ## get averages and make band names
    bands = ac.shared.sensor_wave(sensor)
    bands_data = ac.sentinel3.olci_bands_data()

    waves_ave = np.mean(meta['instrument_data']['lambda0'].data, axis=1)
    waves_mu = np.array(waves_ave)/1000

    ## take wavelengths and band names from external table
    ## the smile correction should bring things in line with these
    waves_names = ['{:.0f}'.format(bands_data[b]['wavelength']) for b in bands_data]
    dnames = ['{}_radiance'.format(b) for b in bands_data]
    bnames = [b for b in bands_data]

    ## average geometry
    sza = np.nanmean(tpg['SZA'])
    vza = np.nanmean(tpg['OZA'])
    raa = np.abs(np.nanmean(tpg['SAA'])-np.nanmean(tpg['OAA']))
    while raa >= 180: raa-=180

    ## compute gas transmittance
    if use_supplied_ancillary:
        ## convert ozone from kg.m-2 to cm.atm
        uoz = np.nanmean(tpg['total_ozone'])/0.02141419
        ## convert water from kg.m-2 to g.cm-2
        uwv = np.nanmean(tpg['total_columnar_water_vapour'])/10
    else:
        ## can get other ancillary here
        uoz = uoz_default
        uwv = uwv_default

    ##
    ttg = ac.ac.gas_transmittance(sza, vza, uoz=uoz, uwv=uwv, waves=waves_ave)

    if use_supplied_pressure:
        #pressure = np.nanmean(tpg['sea_level_pressure'])
        pressure = tpg['sea_level_pressure']
    else:
        ## get pressure from DEM
        if elevation is None:
            if use_supplied_altitude:
                ## read band from NetCDF
                dem = ac.shared.nc_data(ofile, 'altitude')
            else:
                ## interpolate DEM
                dem = ac.dem.hgt_lonlat(tpg['longitude'], tpg['latitude'])
                ## use percentile elevation
                dem = np.nanpercentile(dem, dem_pressure_percentile)
            ## compute DEM pressure
            pressure = ac.ac.pressure_elevation(dem, ratio=False)
        else:
            pressure = ac.ac.pressure_elevation(float(elevation), ratio=False)

    ## get per pixel detector index
    if sub is None:
        di = meta['instrument_data']['detector_index']
    else:
        di = meta['instrument_data']['detector_index'][sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2]]


    ## smile correction - from l2gen smile.c
    if smile_correction:
        smile = {}
        for band in bands_data:
            ## dataset name
            dname = '{}_radiance'.format(band)

            ## band index
            b_i = bands_data[band]['band']-1

            ## gas_correction:
            data['{}_radiance'.format(band)]/=ttg['tt_gas'][b_i]

            if verbosity > 1: print('{} - Smile correction for band {} {} nm'.format(datetime.datetime.now().isoformat()[0:19], band, bands_data[band]['wavelength'] ), end='\n')

            ## bounding bands
            b1_i = bands_data[band]['lower_water']-1
            b2_i = bands_data[band]['upper_water']-1
            band1 = 'Oa{}'.format(str(bands_data[band]['lower_water']).zfill(2))
            band2 = 'Oa{}'.format(str(bands_data[band]['upper_water']).zfill(2))

            ## compute reflectance using per detector F0
            r_ = (data['{}_radiance'.format(band)]) / meta['instrument_data']['solar_flux'][b_i][di]

            ## based on that reflectance, compute radiance for target F0
            r_ *= bands_data[band]['E0']

            ## difference in radiance
            smile[band] = r_-data['{}_radiance'.format(band)]
            del r_ ## free memory

            ## do additional correction based on two bounding bands
            ## currently applying water everywhere
            if bands_data[band]['switch_water'] > 0:
                if verbosity > 1: print('{} - Smile correction - bounding bands {}/{}'.format(datetime.datetime.now().isoformat()[0:19], band1, band2), end='\n')

                ## compute per pixel reflectance difference for bounding bands
                r21_diff = (data['{}_radiance'.format(band2)]) / meta['instrument_data']['solar_flux'][b2_i][di]-\
                           (data['{}_radiance'.format(band1)]) / meta['instrument_data']['solar_flux'][b1_i][di]

                ## wavelength difference ratio
                wdiff_ratio = (bands_data[band]['wavelength'] - meta['instrument_data']['lambda0'][b_i][di])/\
                              (meta['instrument_data']['lambda0'][b2_i][di] - meta['instrument_data']['lambda0'][b1_i][di])

                ## additional smile
                smile[band] += (r21_diff)*(wdiff_ratio)*(meta['instrument_data']['solar_flux'][b_i][di])
                del r21_diff, wdiff_ratio

        ## add smile effect to radiance data
        for band in smile: data['{}_radiance'.format(band)]+=smile[band]
        del smile
    ## end smile correction


    mu = np.cos(tpg['SZA']*dtor)
    time = dateutil.parser.parse(start_time)
    doy = time.strftime('%j')

    ## is the F0 already corrected for Sun - Earth distance?
    ## see OLCI_L2_ATBD_Instrumental_Correction.pdf
    #se = ac.shared.distance_se(doy)
    se = 1.
    se2 = se**2

    ## read TOA
    mask = None
    for iw, wave in enumerate(waves_names):
        if verbosity > 1: print('{} - Reading TOA data for {} nm'.format(datetime.datetime.now().isoformat()[0:19], wave), end='\n')
        # per pixel wavelength
        l = meta['instrument_data']['lambda0'][iw][di]
        # per pixel f0
        f0 = meta['instrument_data']['solar_flux'][iw][di]
        # per pixel fwhm
        fwhm = meta['instrument_data']['FWHM'][iw][di]

        #dname = 'Oa{}_radiance'.format(str(iw+1).zfill(2))
        dname = dnames[iw]
        #d = (np.pi * data[dname] * se2) / (f0*mu)

        ## smile corrected data is in data dict
        if smile_correction:
            d = (np.pi * data[dname] * se2) / (f0*mu)
        else:
            cg = 1.0
            if use_gains:
                if len(gains) == 21:
                    cg = float(gains[iw])
                    print('Applying gain {:.5f} for {}'.format(cg, wave))
            d = (np.pi * ac.shared.nc_data(lfiles[dname], dname, sub=sub) * se2) / (f0*mu)
            d *= cg
            
        ds_att  = {'wavelength':float(wave)}
        for key in ttg: ds_att[key]=ttg[key][iw]

        ac.output.nc_write(ofile, 'rhot_{}'.format(wave), d, dataset_attributes=ds_att)

        if mask is None:
            mask = (d > rhot_mask_threshold)
        else:
            mask[np.where(d > rhot_mask_threshold)] = True
        if float(wave) > 1000.:
            mask[np.where(d > rhot_mask_threshold_swir)] = True

    ac.output.nc_write(ofile, 'mask', mask.astype(np.int32))

    ## clear data
    data = None

    ## set up tiles
    tile_pix = (int(tile_dimensions[0]/sensor_resolution), int(tile_dimensions[1]/sensor_resolution))
    if tiled_processing:
        ntx = int(np.ceil(len(subx)/(tile_pix[0])))
        nty = int(np.ceil(len(suby)/(tile_pix[1])))

        if (ntx < 2) or (nty < 2):
            ntx, nty = 1, 1
            tiled_processing = False
    else:
        ntx, nty = 1, 1
    ntiles = ntx*nty
    tile_dict = {}

    #print('AOT map {} tiles'.format(ntiles))
    if verbosity > 1: print('{} - AOT map will be generated for {} tiles'.format(datetime.datetime.now().isoformat()[0:19], ntiles), end='\n')

    ## LUT dataset to fit to the dark spectrum
    if sky_correction:
        fit_tag = 'romix+rsky_t' ## combined romix rsky
    else:
        fit_tag = 'romix' ## just romix

    ## bt stores the retrieved tau for a given tile and band, and lut
    ## tiles_y, tiles_x, bands, luts
    bt = np.zeros((nty, ntx, len(waves_names), len(luts)))+np.nan

    ## bd stores rhod in each band
    bd = np.zeros((nty, ntx, len(waves_names)))+np.nan

    ## run through bands
    for iw, wave in enumerate(waves_names):
        if float(wave) < rhod_min_wave: continue
        if ttg['tt_gas'][iw] < rhod_tgas_cutoff: continue

        band = rsr_bands[iw]
        dname = dnames[iw]
        bname = bnames[iw]

        ## read band from NetCDF
        d, ds_att = ac.shared.nc_data(ofile, 'rhot_{}'.format(wave), attributes=True)

        ## run through tiles
        for xi in range(ntx):
            cropx = [0 + tile_pix[0]*xi, min(dshape[1], tile_pix[0] + tile_pix[0]*xi)]
            for yi in range(nty):
                cropy = [0 + tile_pix[1]*yi, min(dshape[0], tile_pix[1] + tile_pix[1]*yi)]

                ## subset data
                dsub = d[cropy[0]:cropy[1], cropx[0]:cropx[1]]

                ## get n darkest pixels for this tile
                if dark_spectrum_option == 'intercept':
                    npix=100
                    tmp = dsub.ravel()
                    tmp.sort()

                    x = np.where(np.isfinite(tmp[0:npix]))[0]
                    y = [tmp[0:npix][i] for i in x]

                    ## darkest pixel is the intercept of the regression line pixel, rhot
                    my, by, ry, smy, sby = ac.shared.regression.lsqfity(x, y)
                    rd = by*1

                    ## geometry is tile average
                    raa_ = np.nanmean(tpg['RAA'][cropy[0]:cropy[1], cropx[0]:cropx[1]])
                    vza_ = np.nanmean(tpg['OZA'][cropy[0]:cropy[1], cropx[0]:cropx[1]])
                    sza_ = np.nanmean(tpg['SZA'][cropy[0]:cropy[1], cropx[0]:cropx[1]])

                    if type(pressure) == np.ndarray:
                        pressure_ =  np.nanmean(pressure[cropy[0]:cropy[1], cropx[0]:cropx[1]])
                    else:
                        pressure_ = 1*pressure

                elif (dark_spectrum_option == 'darkest') | (dark_spectrum_option == 'percentile'):
                    rsub = dsub.ravel()

                    ## get minimum reflectance in current subset
                    if dark_spectrum_option == 'darkest':
                        subidx = np.argsort(rsub)
                        ## find pixel index of minimum
                        didx = int(subidx[0]/dsub.shape[1]),subidx[0]%dsub.shape[1]

                    ## get percentile and find brightest pixel darker or equal to percentile
                    if dark_spectrum_option == 'percentile':
                        rsub.sort()
                        rd = np.nanpercentile(rsub, dark_spectrum_percentile)
                        subidx = np.where(dsub <= rd)[0]
                        ## find pixel index of percentile
                        didx = int(subidx[-1]/dsub.shape[1]),subidx[-1]%dsub.shape[1]

                    ## darkest pixel
                    rd = dsub[didx[0], didx[1]]

                    ## geometry for that specific pixel
                    raa_ = tpg['RAA'][cropy[0]+didx[0], cropx[0]+didx[1]]
                    vza_ = tpg['OZA'][cropy[0]+didx[0], cropx[0]+didx[1]]
                    sza_ = tpg['SZA'][cropy[0]+didx[0], cropx[0]+didx[1]]

                    if type(pressure) == np.ndarray:
                        pressure_ =  pressure[cropy[0]+didx[0], cropx[0]+didx[1]]
                    else:
                        pressure_ = 1*pressure

                else:
                    print('Dark spectrum option {} not configured.'.format(dark_spectrum_option))
                    continue


                ## correct for gas transmittance
                if not smile_correction:
                    if gas_transmittance: rd/=ttg['tt_gas'][iw]

                ## store rhod for this band and tile
                bd[yi, xi, iw] = rd

                ## interpolate tau to this band observation for each LUT
                for li, lutid in enumerate(luts):
                    ret_ = lutd[lutid]['rgi'][band]((pressure_, lutd[lutid]['ipd'][fit_tag],
                                                        raa_, vza_, sza_,
                                                       lutd[lutid]['meta']['tau']))
                    bt[yi, xi, iw, li] = np.interp(rd, ret_, lutd[lutid]['meta']['tau'])
        if verbosity > 1: print('{} - Finished tiling and fitting for {} nm'.format(datetime.datetime.now().isoformat()[0:19], wave), end='\n')

        #print('Finished tiling and fitting {} nm'.format(wave))

    ## filter out minimum taus
    bt[bt == 0.001] = np.nan
    ## and fill nans in tiling with closest value
    for i in range(bt.shape[2]): ## bands
        for j in range(bt.shape[3]): ## aerosol models
            ind = scipy.ndimage.distance_transform_edt(np.isnan(bt[:,:,i,j]), return_distances=False, return_indices=True)
            bt[:,:,i,j] = bt[:,:,i,j][tuple(ind)]
    ## end fill tau

    ## get the two bands that give the lowest aot per model
    btau = np.zeros((nty, ntx, 2, len(luts))) + np.nan
    bidx = np.zeros((nty, ntx, 2, len(luts)))

    ## run through models and sort retrieved aot
    for li, lut in enumerate(luts):
        sort = np.argsort(bt[:,:,:,li], axis=2)
        for bi, band in enumerate(bnames):

            bsub = np.where(sort[:,:,0] == bi)
            if len(bsub[0]) > 0:
                btau[bsub[0], bsub[1], 0, li] = bt[bsub[0], bsub[1], bi, li]
                bidx[bsub[0], bsub[1], 0, li] = bi

            bsub = np.where(sort[:,:,1] == bi)
            if len(bsub[0]) > 0:
                btau[bsub[0], bsub[1], 1, li] = bt[bsub[0], bsub[1], bi, li]
                bidx[bsub[0], bsub[1], 1, li] = bi

    ## choose which model to use
    ## here by min tau difference between two bands giving min tau
    for li, lut in enumerate(luts):
        tdiff = btau[:,:,1, li]-btau[:,:,0, li]
        if li == 0:
            sel_mod = np.zeros((nty, ntx))
            sel_tau = btau[:,:,0, li]
            sel_dif = tdiff
        else:
            tsub = np.where(tdiff<sel_dif)
            if len(tsub[0])>0:
                sel_mod[tsub[0], tsub[1]] = li
                sel_tau[tsub[0], tsub[1]] = btau[tsub[0], tsub[1],0, li]
                sel_dif[tsub[0], tsub[1]] = tdiff[tsub[0], tsub[1]]

    ## make global attributes
    gatts={'sensor':sensor,
           'isodate':start_time,
           'obase':'{}'.format(oname),
           'tiles':ntiles}

    if limit is not None:
        gatts['limit']=limit
        gatts['sub']=sub
    update_attributes=True

    ## tile coordinates
    if sub is None:
        tlx = np.linspace(0,data_shape[1]-1,ntx)
        tly = np.linspace(0,data_shape[0]-1,nty)
    else:
        tlx = np.linspace(sub[0]+0.5,sub[0]+sub[2]+0.5,ntx)
        tly = np.linspace(sub[1]+0.5,sub[1]+sub[3]+0.5,nty)

    ## compute difference of each tile with the maximum of its surroundings
    #f1 = np.ones((3,3))
    #f1[1,1] = 0
    #dif_tau = np.abs(sel_tau - scipy.ndimage.maximum_filter(sel_tau, footprint=f1, mode='reflect'))

    ## compute difference of each tile with the median of its surroundings
    f1 = np.ones((5,5))
    f1[2,2] = 0
    #dif_tau = np.abs(scipy.ndimage.median_filter(sel_tau, footprint=f1, mode='reflect') - sel_tau)
    dif_tau = scipy.ndimage.median_filter(sel_tau, footprint=f1, mode='reflect') - sel_tau

    ## mask and replace nans
    #fil_tau = sel_tau * 1
    #fil_tau[dif_tau > dif_tau_max] = np.nan
    #ind = scipy.ndimage.distance_transform_edt(np.isnan(fil_tau), return_distances=False, return_indices=True)

    ## mask and replace outliers
    ind = scipy.ndimage.distance_transform_edt(dif_tau > dif_tau_max, return_distances=False, return_indices=True)
    fil_tau = sel_tau[tuple(ind)]
    fil_mod = sel_mod[tuple(ind)]

    if tiled_processing:
        ## function to nn resize an array to a given shape
        def nn_resize(arr, size):
            return(np.asarray([[ arr[int(arr.shape[0] * r / size[0])][int(arr.shape[1] * c / size[1])]
                     for c in range(size[1])] for r in range(size[0])]))

        ## interpolate selected tau
        #z = scipy.interpolate.interp2d(tlx, tly, fil_tau)
        #fil_tau_f = z(subx, suby)

        ## interpolate selected model (nearest neighbour)
        #fil_mod_f = nn_resize(fil_mod, dshape)

        ## interpolate selected tau
        if tiled_interpolation == 'tpg':
            z = scipy.interpolate.interp2d(tlx, tly, fil_tau)
            fil_tau_f = z(subx, suby)

            ## interpolate selected model (nearest neighbour)
            fil_mod_f = nn_resize(fil_mod, dshape)

        if tiled_interpolation == 'pixel':
            ## new method
            ## interpolate selected tau per_model
            for li, lutid in enumerate(luts):
                z = scipy.interpolate.interp2d(tlx, tly, btau[:,:,0,li])
                tmp = z(subx, suby)
                tmp = scipy.ndimage.gaussian_filter(tmp, 24, mode='reflect')

                if li == 0:
                    pp_tau = 1 * tmp
                else:
                    pp_tau = np.dstack((pp_tau, tmp))

            ## interpolate selected model
            z = scipy.interpolate.interp2d(tlx, tly, fil_mod)
            pp_mod = z(subx, suby)
            pp_mod = scipy.ndimage.gaussian_filter(pp_mod, 24, mode='reflect')
    else:
        fil_mod_f = np.squeeze(fil_mod)
        fil_tau_f = np.squeeze(fil_tau)

    if ntiles == 1:
        print('Fixed tau')
        gatts['ac_aot550'] = fil_tau_f
        gatts['ac_model'] = fil_mod_f
        #gatts['ac_band'] = fil_band_f

    ## make average tpg tiles
    if (tiled_processing) & (tiled_interpolation == 'tpg'):
        tpg_tiles = {}
        ## run through tiles
        for xi in range(ntx):
            cropx = [0 + tile_pix[0]*xi, min(dshape[1], tile_pix[0] + tile_pix[0]*xi)]
            for yi in range(nty):
                cropy = [0 + tile_pix[1]*yi, min(dshape[0], tile_pix[1] + tile_pix[1]*yi)]
                for t in tpg:
                    if t not in tpg_tiles:
                        tpg_tiles[t] = np.zeros((nty, ntx))
                    tpg_tiles[t][yi, xi] = np.nanmean(tpg[t][cropy[0]:cropy[1], cropx[0]:cropx[1]])

                ## add pressure
                if 'pressure' not in tpg_tiles:
                        tpg_tiles['pressure'] = np.zeros((nty, ntx))
                if type(pressure) == np.ndarray:
                    tpg_tiles['pressure'][yi, xi] = np.nanmean(pressure[cropy[0]:cropy[1], cropx[0]:cropx[1]])
                else:
                    tpg_tiles['pressure'][yi, xi] = pressure

    ## compute rhos
    pars = ['romix','dtott','utott','astot']
    for iw, wave in enumerate(waves_names):
        band = bnames[iw]

        if verbosity > 0: print('{} - Processing band {} nm to surface reflectance'.format(datetime.datetime.now().isoformat()[0:19], wave), end='\n')
        band_data, ds_att = ac.shared.nc_data(ofile, 'rhot_{}'.format(wave), attributes=True)
        ds_att = {d:ds_att[d] for d in ds_att if d not in ['_FillValue']}

        if not smile_correction:
            if gas_transmittance:
                band_data/= ds_att['tt_gas']
                if verbosity > 1: print('Gas transmittance {:.4f}'.format(ds_att['tt_gas']))

        ## get tile averaged parameters
        if (tiled_processing) & (tiled_interpolation == 'tpg'):
            ac_pars_tile = {par: np.zeros((nty, ntx)) for ip, par in enumerate(pars)}

            for li, lut in enumerate(luts):
                modidx = np.where(fil_mod!=li)
                if len(modidx[0]) == 0: continue

                for ip, par in enumerate(pars):
                    if par == 'romix':
                        par_ = '{}'.format(fit_tag)
                    else:
                        par_ = '{}'.format(par)

                    if (tiled_processing) & (tiled_interpolation == 'tpg'):
                        tmp = lutd[lut]['rgi'][band]((tpg_tiles['pressure'], lutd[lut]['ipd'][par_],
                                                    tpg_tiles['RAA'], tpg_tiles['OZA'], tpg_tiles['SZA'], fil_tau))
                    else:
                        tmp = lutd[lut]['rgi'][band]((tpg['sea_level_pressure'], lutd[lut]['ipd'][par_],
                                                    tpg['RAA'], tpg['OZA'], tpg['SZA'], fil_tau))

                    ac_pars_tile[par][modidx] = tmp[modidx]
                #print(ac_pars_tile.keys())

        ## make full size ac_pars
        ac_pars = {}
        for par in pars:
            if par == 'romix':
                par_ = '{}'.format(fit_tag)
            else:
                par_ = '{}'.format(par)

            ## tiled processing
            if tiled_processing:
                if tiled_interpolation == 'tpg':
                    #tmp = ac_pars_tile[par] * 1
                    tmp = scipy.ndimage.gaussian_filter(ac_pars_tile[par], 1.5, mode='reflect')

                    z = scipy.interpolate.interp2d(tlx, tly, tmp)
                    ac_pars[par] = z(subx, suby)
                    ac_pars[par] = scipy.ndimage.gaussian_filter(ac_pars[par], 24, mode='reflect')
                    #for li, lutid in enumerate(luts)
                    #    res = lutd[lut]['rgi'][band]((pressure, lutd[lut]['ipd'][par_],
                    #                                  tpg['RAA'], tpg['OZA'], tpg['SZA'], pp_tau[:,:,li]))

                    #ac_pars[par] = lutd[lut]['rgi'][band]((pressure, lutd[lut]['ipd'][par_],
                    #                                       tpg['RAA'], tpg['OZA'], tpg['SZA'], fil_tau))

                ### quite a bit slower
                if tiled_interpolation == 'pixel':
                    for li, lut in enumerate(luts):
                        ## get parameter for current model and per pixel tau
                        res = lutd[lut]['rgi'][band]((pressure, lutd[lut]['ipd'][par_],
                                                                          tpg['RAA'], tpg['OZA'], tpg['SZA'], pp_tau[:,:,li]))
                        ## weigh according to model
                        if li == 0:
                            ac_pars[par] = (res * (1-pp_mod))
                        elif li == 1:
                            ac_pars[par] += (res * pp_mod)

            else:
                ## fixed processing
                ac_pars[par] = lutd[lut]['rgi'][band]((tpg['sea_level_pressure'], lutd[lut]['ipd'][par_],
                                                       tpg['RAA'], tpg['OZA'], tpg['SZA'], fil_tau))

            #print(par_, ac_pars[par].shape)

            if write_resolved_atmosphere:
                ac.output.nc_write(ofile, '{}_{}'.format(par, wave), ac_pars[par])

        ## get rhos data
        rhos_data = (band_data) - ac_pars['romix']
        rhos_data = (rhos_data) / ((ac_pars['utott']*ac_pars['dtott']) + ac_pars['astot'] * rhos_data)
        band_data = None
        ac_pars=None

        rhos_data[np.isinf(rhos_data)]=np.nan
        ac.output.nc_write(ofile, 'rhos_{}'.format(wave), rhos_data, dataset_attributes=ds_att,
                            attributes=gatts, fillvalue=np.nan, update_attributes=update_attributes)
        update_attributes = False
        rhos_data = None

    print('{} - File written to {}'.format(datetime.datetime.now().isoformat()[0:19], ofile), end='\n')
    if map_rgb:
        ac.acolite.acolite_map(inputfile=ofile, output='{}/plotted'.format(output),
                               rgb_rhot = True, rgb_rhos = True, mapped=False)

    print('{} - Finished processing {}'.format(datetime.datetime.now().isoformat()[0:19], bundle))
