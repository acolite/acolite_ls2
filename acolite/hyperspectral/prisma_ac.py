## def prisma_ac
##
## QV 2020-01-14
##
## modifications: 2020-06-17 (QV) updated processing with V2020 OpticsExpress
##                2020-10-15 (QV) functionised
##                2020-10-20 (QV) added to acolite, made aot_average 1-1.8 micron default processing option
##                2020-10-29 (QV) fixed gas correction order in rhos computation, added aot and model to output NC file

def prisma_ac(file, output=None,
              uoz_default = 0.3,
              uwv_default = 1.5,
              tgas_cutoff = 0.95,
              ancillary_data = True,
              gas_correction = True,
              sky_correction = False,
              pressure=1013.25,

              aot_selection= 'average', #'minimum',
              aot_fixed=None,
              aerosol_model=None,
              correct_rhod=True,

              range_wave1 = 1000, #None
              range_wave2 = 1800, #None
              exclude_wave = [402, 750, 2050],
              exclude_window = [10, 10, 100],
              dark_prc_idx = 2,
              percentiles = (0, 0.1,1,5,10,25,50,75,90,95,99,99.9),

              make_plot=True, make_rgb=True):

    import acolite as ac

    if make_plot:
        import matplotlib.pyplot as plt

    import numpy as np
    import h5py
    import os
    import scipy.interpolate
    import dateutil

    ## import luts
    lutd = ac.aerlut.import_luts()
    ## import rsky luts
    rlutd = {}
    cluts = list(lutd.keys())
    mods = [int(l[-1]) for l in cluts]
    for aermod in mods:
        l, m, d, r = ac.aerlut.rsky_read_lut(aermod)
        rlutd[aermod] = {'lut':l, 'meta':m, 'dim':d, 'rgi':r}


    ## def get_cur_rhot
    ### get rhot with sky glint for current model and tau
    ## QV 2020-06-17
    def get_cur_rhot(ta, lutid, aermod):

            ipd = {par:ip for ip,par in enumerate(lutd[lutid]['meta']['par'])}
            cur_rsky = rlutd[aermod]['rgi']((lutd[lutid]['meta']['wave'], raa, vza, sza, ta))
            cur_rsky[np.isnan(cur_rsky)] = 0

            cur_romix = lutd[lutid]['rgi']((pressure, ipd['romix'], lutd[lutid]['meta']['wave'], raa, vza, sza, ta))
            cur_utott = lutd[lutid]['rgi']((pressure, ipd['utott'], lutd[lutid]['meta']['wave'], raa, vza, sza, ta))
            cur_dtott = lutd[lutid]['rgi']((pressure, ipd['dtott'], lutd[lutid]['meta']['wave'], raa, vza, sza, ta))
            cur_astot = lutd[lutid]['rgi']((pressure, ipd['astot'], lutd[lutid]['meta']['wave'], raa, vza, sza, ta))
            cur_romix[np.isnan(cur_romix)] = 0

            ## model toa reflectance for this tau
            cur_rhot = cur_romix + (cur_utott * cur_dtott * cur_rsky) / (1. - cur_rsky * cur_astot)
            return(cur_rhot, cur_romix)


    ## open HDF5 file
    f = h5py.File(file, mode='r')
    h5_gatts = {}
    for a in f.attrs.keys():
        h5_gatts[a] = f.attrs[a]


    ## read wavelengths and make rsr/f0
    waves_vnir = h5_gatts['List_Cw_Vnir']
    bands_vnir = ['{:.0f}'.format(w) for w in waves_vnir]
    fwhm_vnir = h5_gatts['List_Fwhm_Vnir']
    n_vnir = len(waves_vnir)

    waves_swir = h5_gatts['List_Cw_Swir']
    bands_swir = ['{:.0f}'.format(w) for w in waves_swir]
    fwhm_swir = h5_gatts['List_Fwhm_Swir']
    n_swir = len(waves_swir)

    f0_vnir_cw = [ac.shared.f0_wave(w) for w in waves_vnir]
    f0_swir_cw = [ac.shared.f0_wave(w) for w in waves_swir]

    rsr_vnir = [ac.shared.gauss_response(waves_vnir[b], fwhm_vnir[b], step=0.1) for b in range(0, n_vnir)]
    rsr_swir = [ac.shared.gauss_response(waves_swir[b], fwhm_swir[b], step=0.1) for b in range(0, n_swir)]

    f0_vnir = [ac.shared.f0_band(rsr_vnir[b][0],rsr_vnir[b][1]) for b in range(n_vnir)]
    f0_swir = [ac.shared.f0_band(rsr_swir[b][0],rsr_swir[b][1]) for b in range(n_swir)]

    f0_all = f0_vnir + f0_swir
    f0_all_centre = f0_vnir_cw + f0_swir_cw

    waves = [w for w in waves_vnir] + [w for w in waves_swir]
    band_rsr = [w for w in rsr_vnir] + [w for w in rsr_swir]

    waves_names = ['{:.0f}'.format(w) for w in waves]
    instrument = ['vnir']*n_vnir + ['swir']*n_swir
    band_index = [i for i in range(n_vnir)] + [i for i in range(n_swir)]

    ## order to store wavelengths
    idx = np.argsort(waves)
    bands = {}
    for i in idx:
        cwave = waves[i]
        if cwave == 0: continue

        swave = '{:.0f}'.format(cwave)
        bands[swave]= {'wave':cwave, 'wavelength':cwave, 'wave_mu':cwave/1000.,
                       'wave_name':waves_names[i],
                       'i':i, 'index':band_index[i],
                       'rsr': band_rsr[i],
                       'f0': f0_all[i], 'f0_centre': f0_all_centre[i],
                       'instrument':instrument[i],}


    ## output global attributes
    gatts = {}
    sza = h5_gatts['Sun_zenith_angle']
    saa = h5_gatts['Sun_azimuth_angle']
    cossza = np.cos(sza*(np.pi/180.))

    isotime = h5_gatts['Product_StartTime']
    time = dateutil.parser.parse(isotime)

    doy = int(time.strftime('%j'))
    d = ac.shared.distance_se(doy)

    ## compute view geometry from  LOS vector
    x_ = f['KDP_AUX']['LOS_Vnir'][:,0]
    y_ = f['KDP_AUX']['LOS_Vnir'][:,1]
    z_ = f['KDP_AUX']['LOS_Vnir'][:,2]
    dtor = np.pi/180
    vza = np.arctan2(y_,x_)/dtor
    vaa = np.arctan2(z_,np.sqrt(x_**2+y_**2))/dtor
    vza_ave = np.nanmean(np.abs(vza))
    vaa_ave = np.nanmean(np.abs(vaa))

    gatts['THS'] = sza
    gatts['THV'] = vza_ave

    gatts['PHIS'] = saa
    gatts['PHIV'] = vaa_ave
    print(gatts)

    mu0 = np.cos(gatts['THS']*(np.pi/180))
    muv = np.cos(gatts['THV']*(np.pi/180))

    ths = gatts['THS']
    thv = gatts['THV']
    azi = abs(gatts['PHIS'] - gatts['PHIV'])
    while azi >= 180:
        azi = abs(azi-360)


    ## read lat and lon
    src = 'HCO'
    lat = f['HDFEOS']['SWATHS']['PRS_L1_{}'.format(src)]['Geolocation Fields']['Latitude_SWIR'][:]
    lon = f['HDFEOS']['SWATHS']['PRS_L1_{}'.format(src)]['Geolocation Fields']['Longitude_SWIR'][:]


    ## get ancillary data
    ## set EarthData u/p from config
    if len(ac.acolite.config['EARTHDATA_u']) > 0:
        os.environ['EARTHDATA_u'] = ac.acolite.config['EARTHDATA_u']
    if len(ac.acolite.config['EARTHDATA_p']) > 0:
        os.environ['EARTHDATA_p'] = ac.acolite.config['EARTHDATA_p']


    ## get gas transmittances
    uoz=uoz_default*1.0
    uwv=uwv_default*1.0


    if ancillary_data:
        print('Getting ancillary data')
        pc_date = time.strftime('%Y-%m-%d')
        pc_time = time.hour + time.minute/60. + time.second/3600.
        pc_lon = np.nanpercentile(lon, 50)
        pc_lat = np.nanpercentile(lat, 50)


        ## get ancillary data
        pc_anc = ac.ac.ancillary.ancillary_get(pc_date, pc_lon, pc_lat, ftime=pc_time, kind='nearest')

        ## get pressure from ancillary data if not determined by user or by DEM
        pressure = pc_anc['press']['interp']

        if 'ozone' in pc_anc: uoz=pc_anc['ozone']['interp']
        else: print('No ancillary ozone found: using default {}.'.format(uoz))
        if 'p_water' in pc_anc: uwv=pc_anc['p_water']['interp']/10.
        else:print('No ancillary ozone found: using default {}.'.format(uwv))

    ## get ozone
    koz = ac.shared.ko3_get()
    koz['wave_nm'] = [w*1000 for w in koz['wave']]

    ## band ozone transmittance
    for bi, b in enumerate(bands):
        ## resampled to band
        #data, wave, rsr, rsrwave
        toz = ac.shared.rsr_convolute(koz['data'], koz['wave_nm'], bands[b]['rsr'][1], bands[b]['rsr'][0]) * uoz
        tt_oz = np.exp(-1.*(toz) / mu0) * np.exp(-1.*(toz) / muv)
        bands[b]['tt_oz']= tt_oz

        ## band centre
        c_band, c_wave = ac.shared.widx(koz['wave'], bands[b]['wave_mu'])
        toz = koz['data'][c_band] * uoz
        tt_oz = np.exp(-1.*(toz) / mu0) * np.exp(-1.*(toz) / muv)
        bands[b]['tt_oz_centre']= tt_oz

    ## get water vapour
    wave_wv, t_wv = ac.ac.wvlut_interp(gatts['THS'], gatts['THV'], uwv=uwv)
    wave_wv_nm = [w*1000 for w in wave_wv]

    ## band wv transmittance
    for bi, b in enumerate(bands):
        t_wv_b = ac.shared.rsr_convolute(t_wv, wave_wv_nm, bands[b]['rsr'][1], bands[b]['rsr'][0])
        bands[b]['tt_wv']= t_wv_b

        c_band, c_wave = ac.shared.widx(wave_wv, bands[b]['wave_mu'])
        bands[b]['tt_wv_centre'] = t_wv[c_band]

    ## get oxygen
    wave_o2, t_o2 = ac.ac.o2lut_interp(gatts['THS'], gatts['THV'])
    wave_o2_nm = [w*1000 for w in wave_o2]

    ## band o2 transmittance
    for bi, b in enumerate(bands):
        t_o2_b = ac.shared.rsr_convolute(t_o2, wave_o2_nm, bands[b]['rsr'][1], bands[b]['rsr'][0])
        bands[b]['tt_o2']= t_o2_b

        c_band, c_wave = ac.shared.widx(wave_o2, bands[b]['wave_mu'])
        bands[b]['tt_o2_centre']= t_o2[c_band]

    ## sky reflectance
    dtor = np.pi / 180.
    rsky = []
    rsky_centre =[]
    wave_hyp_nm = np.arange(250,1100, step=1)
    wave_hyp = wave_hyp_nm/1000.
    rsky_hyp = ac.ac.rayleigh.ray_refl_onlysky(wave_hyp, gatts['THS']*dtor,
                                                           max(0.01, gatts['THV'])*dtor,
                                                           gatts['PHIV']*dtor,
                                                           gatts['PHIS']*dtor, Patm=pressure)
    for bi, b in enumerate(bands):
        rsky_b = ac.shared.rsr_convolute(rsky_hyp, wave_hyp_nm, bands[b]['rsr'][1], bands[b]['rsr'][0])
        bands[b]['rsky_ray']= rsky_b

        only_sky = ac.ac.rayleigh.ray_refl_onlysky(bands[b]['wave_mu'], gatts['THS']*dtor,
                                                           max(0.01,gatts['THV'])*dtor,
                                                           gatts['PHIV']*dtor,
                                                           gatts['PHIS']*dtor, Patm=pressure)
        bands[b]['rsky_ray_centre']= only_sky

    ## total gas transmittance
    for bi, b in enumerate(bands):
        bands[b]['tt_gas'] = bands[b]['tt_wv'] * bands[b]['tt_oz'] * bands[b]['tt_o2']
        bands[b]['tt_gas_centre'] = bands[b]['tt_wv_centre'] * bands[b]['tt_oz_centre'] * bands[b]['tt_o2_centre']

    ## create output file
    obase  = '{}_{}_L2R'.format('PRISMA',  time.strftime('%Y_%m_%d_%H_%M_%S'))
    if not os.path.exists(output):
        os.makedirs(output)
    ofile = '{}/{}.nc'.format(output, obase)

    gatts['obase'] = obase
    gatts['sensor'] = 'PRISMA_PRISMA'
    gatts['isodate'] = time.isoformat()

    ac.output.nc_write(ofile, 'lat', np.flip(np.rot90(lat)), new=True, attributes=gatts)
    ac.output.nc_write(ofile, 'lon', np.flip(np.rot90(lon)))

    ## read bands in spectral order
    read_cube = True
    prc = []
    rdark = {}
    if read_cube:
        vnir_data = h5_gatts['Offset_Vnir'] + \
                    f['HDFEOS']['SWATHS']['PRS_L1_{}'.format(src)]['Data Fields']['VNIR_Cube'][:]/h5_gatts['ScaleFactor_Vnir']

        swir_data = h5_gatts['Offset_Swir'] + \
                    f['HDFEOS']['SWATHS']['PRS_L1_{}'.format(src)]['Data Fields']['SWIR_Cube'][:]/h5_gatts['ScaleFactor_Swir']

    for bi, b in enumerate(bands):
        wi = bands[b]['index']
        i = bands[b]['i']

        print('Reading rhot_{}'.format(bands[b]['wave_name']))

        if bands[b]['instrument'] == 'vnir':
            if read_cube:
                cdata_radiance = vnir_data[:,wi,:]
                cdata = cdata_radiance * (np.pi * d * d) / (bands[b]['f0'] * 10 * cossza)
            else:
                cdata_radiance = h5_gatts['Offset_Vnir'] + \
                        f['HDFEOS']['SWATHS']['PRS_L1_{}'.format(src)]['Data Fields']['VNIR_Cube'][:,i,:]/h5_gatts['ScaleFactor_Vnir']
                cdata = cdata_radiance * (np.pi * d * d) / (bands[b]['f0'] * 10 * cossza)

        if bands[b]['instrument'] == 'swir':
            if read_cube:
                cdata_radiance = swir_data[:,wi,:]
                cdata = cdata_radiance * (np.pi * d * d) / (bands[b]['f0'] * 10 * cossza)
            else:
                cdata_radiance = h5_gatts['Offset_Swir'] + \
                        f['HDFEOS']['SWATHS']['PRS_L1_{}'.format(src)]['Data Fields']['SWIR_Cube'][:,i,:]/h5_gatts['ScaleFactor_Swir']
                cdata = cdata_radiance * (np.pi * d * d) / (bands[b]['f0'] * 10 * cossza)

        ds_att = {k:bands[b][k] for k in bands[b] if k not in ['rsr']}

        if True:
            ## write toa radiance
            ac.output.nc_write(ofile, 'Lt_{}'.format(bands[b]['wave_name']), np.flip(np.rot90(cdata_radiance)),
                          dataset_attributes = ds_att)



        ## write toa reflectance
        ac.output.nc_write(ofile, 'rhot_{}'.format(bands[b]['wave_name']), np.flip(np.rot90(cdata)),
                          dataset_attributes = ds_att)

        rd = np.nanpercentile(cdata, percentiles)
        prc.append(rd[dark_prc_idx])
        rdark[b] = rd



    ## fit aot
    lut_results = {}

    raa = azi*1
    vza = thv*1
    sza = sza*1

    ## get aot
    if (aot_fixed is not None) & (aerosol_model is not None):
        sel_aot = float(aot_fixed)
    else:
        for lutid in lutd:
            lut_data = {}
            lut_data_rsky = {}

            taus = lutd[lutid]['meta']['tau']

            waves_nm = [w*1000 for w in lutd[lutid]['meta']['wave']]

            for tau550 in taus:
                ## get different parameters for this TAU550
                #cur = {}
                #cur_c = {}
                #for par in ['romix','rorayl','dtotr', 'utotr', 'dtott','utott', 'astot']:
                #for par in ['romix']:
                #    ip = [i for i,value in enumerate(lutd[lutid]['meta']['par']) if value == par]
                #    ret = lutd[lutid]['rgi']((pressure, ip, lutd[lutid]['meta']['wave'], azi, thv, ths, tau550))
                #    cur[par] = [ac.shared.rsr_convolute(ret, waves_nm, bands[b]['rsr'][1], bands[b]['rsr'][0]) for b in bands]
                #    #cur_c[par] = [lutd[lutid]['rgi']((pressure, ip, bands[b]['wave_mu'], azi, thv, ths, tau550)) for b in bands]
                #lut_data[tau550] = cur
                ##lut_data_centre[tau550] = cur_c

                cur_rhot, cur_romix = get_cur_rhot(tau550, lutid, int(lutid[-1]))
                lut_data_rsky[tau550] = {'rhot': [ac.shared.rsr_convolute(cur_rhot, waves_nm, bands[b]['rsr'][1], bands[b]['rsr'][0]) for b in bands],
                                         'romix': [ac.shared.rsr_convolute(cur_romix, waves_nm, bands[b]['rsr'][1], bands[b]['rsr'][0]) for b in bands]}


            ## get aot results from all bands
            aot = []
            aot_filt = {}
            for i,dk in enumerate(rdark.keys()):
                if bands[dk]['tt_gas'] < tgas_cutoff:
                    tau = np.nan
                else:
                    if True:
                        ## include sky reflectance
                        ratm = [lut_data_rsky[tau]['rhot'][i] for tau in taus]
                    else:
                        ## do not include sky reflectance
                        ratm = [lut_data[tau]['romix'][i] for tau in taus]

                    tau_id, tau_br = ac.aerlut.lutpos(ratm, rdark[dk][dark_prc_idx]/bands[dk]['tt_gas'])
                    tau = np.interp(tau_id, range(len(taus)),taus)
                    if tau == 0.001: tau = np.nan
                aot.append(tau)
                aot_filt[dk] = tau
            sel_aot = np.nanmin(aot)

            ## exclude wavelengths
            if not ((exclude_wave is None) & (exclude_window is None)):
                for iw, ew in enumerate(exclude_wave):
                    aot_filt = {w:aot_filt[w] for i, w in enumerate(aot_filt.keys()) if not ((bands[w]['wave'] > ew - exclude_window[iw]/2) & (bands[w]['wave'] < ew + exclude_window[iw]/2))}

            ## include range
            if not ((range_wave1 is None) & (range_wave2 is None)):
                aot_filt = {w:aot_filt[w] for i, w in enumerate(aot_filt.keys()) if ((bands[w]['wave'] > range_wave1) & (bands[w]['wave'] < range_wave2))}
                #print(len(aot_filt))


            ## select aot value
            if aot_selection == 'minimum':
                sel_aot = np.nanmin([aot_filt[k] for k in aot_filt.keys()])
                sel_band = [k for k in aot_filt.keys() if aot_filt[k] == sel_aot]
                if len(sel_band) == 1:
                    sel_band = sel_band[0]
                else:
                    print('Multiple bands with same aot: {}'.format(sel_band))
                    sel_band = sel_band[0]
                print('Selected wavelength: {:.1f} nm, aot_550: {:.3f}'.format(bands[sel_band]['wave'], sel_aot))


            elif aot_selection == 'average':
                sel_aot = np.nanmean([aot_filt[k] for k in aot_filt.keys()])
                if not ((range_wave1 is None) & (range_wave2 is None)):
                    print('Average aot_550 ({:.1f}-{:.1f} nm): {:.3f}'.format(range_wave1, range_wave2, sel_aot))
                else:
                    print('Average aot_550: {:.3f}'.format(sel_aot))


            ###
            pars = lutd[lutid]['meta']['par']
            sel = {}
            #sel_centre = {}
            for par in pars:
                ip = [i for i,value in enumerate(lutd[lutid]['meta']['par']) if value == par]
                ret = lutd[lutid]['rgi']((pressure, ip, lutd[lutid]['meta']['wave'], azi, thv, ths, sel_aot))
                sel[par] = [ac.shared.rsr_convolute(ret, waves_nm, bands[b]['rsr'][1], bands[b]['rsr'][0]) for b in bands]

            ret_rsky = rlutd[int(lutid[-1])]['rgi']((lutd[lutid]['meta']['wave'], azi, thv, ths, sel_aot))
            sel['rsky'] = [ac.shared.rsr_convolute(ret_rsky, waves_nm, bands[b]['rsr'][1], bands[b]['rsr'][0]) for b in bands]

            sel = {k: np.asarray(sel[k]) for k in sel}
            sel['rsky_toa'] = (sel['utott'] * sel['dtott'] * sel['rsky'])/(1. - sel['rsky'] * sel['astot'])


            lut_results[lutid] =  {'ret':sel, 'aot550':sel_aot, 'aot_all':aot} #, 'sel':sel, 'sel_centre': sel_centre}

    ## select model
    sel_rmsd = 10
    sel_lut = None
    x = [rdark[b][2]/bands[b]['tt_gas'] for b in bands if bands[b]['tt_gas']>tgas_cutoff]
    wvs = [bands[b]['wave'] for b in bands if bands[b]['tt_gas']>tgas_cutoff]
    for lutid in lut_results:
        y = [lut_results[lutid]['ret']['romix'][ib]+lut_results[lutid]['ret']['rsky_toa'][ib] for ib, b in enumerate(bands) if bands[b]['tt_gas']>tgas_cutoff]
        rmsd = ac.shared.rmsd(x, y)
        lut_results[lutid]['rmsd'] = rmsd
        if rmsd < sel_rmsd:
            sel_lut = lutid
            sel_rmsd = rmsd*1
    sel = lut_results[sel_lut]['ret']

    gatts['ac_model'] = sel_lut
    gatts['ac_aot550'] = lut_results[sel_lut]['aot550']

    ## compute surface reflectance
    for bi, b in enumerate(bands):
        wi = band_index[i]
        if bands[b]['wave'] > 2400: continue # skip this for now - not covered by Thuillier F0
        if bands[b]['tt_gas'] < 0.85: continue # skip bands with low gas transmittance

        print('Computing rhos_{}'.format(bands[b]['wave_name']))

        cur_data, cur_att = ac.shared.nc_data(ofile, 'rhot_{}'.format(b), attributes=True)

        if gas_correction:
            cur_data /= cur_att["tt_gas"]
        tt_gas_cur = 1.0

        if sky_correction:
            cur_att['rsky_toa'] = lut_results[sel_lut]['ret']['rsky_toa'][bi]
            cur_data -= cur_att['rsky_toa']


        ## compute surface reflectance
        cur_data = ac.shared.rtoa_to_rhos(cur_data, sel['romix'][bi],
                                        sel['utott'][bi], sel['dtott'][bi], sel['astot'][bi], tt_gas = tt_gas_cur)

        ## store parameters
        for ip, par in enumerate(['romix','rorayl','dtotr', 'utotr', 'dtott','utott', 'astot']):
            cur_att[par] = sel[par][bi]

        ## mask output data
        cur_data[cur_data > 1.5] = np.nan
        ac.output.nc_write(ofile, 'rhos_{}'.format(b),
                            cur_data, dataset_attributes=cur_att, attributes=gatts, update_attributes=(bi==0))
        new = False


    if make_plot:
        plt.plot([bands[b]['wave'] for b in bands], [rdark[b][2] for b in bands], color='Grey')

        x = [rdark[b][2]/bands[b]['tt_gas'] for b in bands if bands[b]['tt_gas']>tgas_cutoff]
        wvs = [bands[b]['wave'] for b in bands if bands[b]['tt_gas']>tgas_cutoff]

        plt.plot(wvs, x, '.:', color='Black', label=r'$\rho_d$ selected')
        for lutid in lut_results:
            y = [lut_results[lutid]['ret']['romix'][ib]+lut_results[lutid]['ret']['rsky_toa'][ib] for ib, b in enumerate(bands) if bands[b]['tt_gas']>tgas_cutoff]
            rmsd = ac.shared.rmsd(x, y)
            plt.plot(wvs, y, '.:', label=r'{} $\tau_a$: {:.3f}, RMSD: {:.1e}'.format(lutid.split('-')[-1],
                                                                               lut_results[lutid]['aot550'], rmsd))

        plt.legend()
        plt.xlabel('Wavelength (nm)')
        plt.ylabel(r'$\rho$ (-)')
        plt.savefig(ofile.replace('_L2R.nc', '_dark_spectrum.png'), dpi=300, bbox_inches='tight')
        plt.close();


    if make_rgb:
        ac.acolite.acolite_map(inputfile=ofile, output=output,
                       map_colorbar_edge=False,
                       mapped=False,
                       map_scalepos='UL',
                       rgb_rhos=True, rgb_rhot=True)
