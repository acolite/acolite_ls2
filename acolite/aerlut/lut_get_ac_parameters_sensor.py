## def lut_get_ac_parameters_sensor
## returns ac parameters from sensor LUT for given rtoa, azi, thv and ths (and LUT)
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2016-07-13
## modifications: 2016-12-13 (QV) avoid selecting tau equal to first LUT step
##                2017-10-19 (QV) added force_band option
##                2017-12-07 (QV) changed all tau to dicts
##                2018-07-18 (QV) changed acolite import name
def lut_get_ac_parameters_sensor(lut_sensor,meta,azi,thv,ths,rtoa, force_band=None):
    import acolite as pp

    ## calculate tau for each band observation
    tau_dict = pp.aerlut.lut_get_tau_sensor(lut_sensor,meta,azi,thv,ths,rtoa)

    ## get band wavelengths and order on wavelength
    band_wl = pp.aerlut.interplut_sensor(lut_sensor, meta, azi, thv, ths, 0., par='wl')
    bo = sorted(band_wl, key=band_wl.get)
    band_order = [b for i,b in enumerate(bo) if b in rtoa]

    ## extract from dict
    band_wave = [band_wl[b] for i,b in enumerate(band_order)]
    tau550_all_bands = [tau_dict[b] for i,b in enumerate(band_order)]
    band_rtoa = [rtoa[b] for b in band_order]

    #tau550_all_bands_rmsd = []
    tau550_all_bands_rmsd={}
    for band in band_order:
        if band in rtoa:
            ## use precomputed sensor specific LUT -- faster
            ratm_ = pp.aerlut.interplut_sensor(lut_sensor, meta, azi, thv, ths, tau_dict[band], par='romix')
            ## use sensor RSR and compute now -- slower
            ratm_1 = [ratm_[b] for b in band_order]
            tau550_all_bands_rmsd[band] = pp.rmsd(band_rtoa, ratm_1)

    if force_band == None:
        ## select lowest TAU550 to avoid <0 reflectances
        ## if tau550 is equal to minimum table tau then select next band
        tau550_inlut = [i for i in tau550_all_bands if (i > min(meta['tau']))]
        if len(tau550_inlut) > 0: 
             tau550 = min(tau550_inlut)
             ## find band to which this tau corresponds
             dark_idx = [i for i,value in enumerate(tau550_all_bands) if value == tau550][0]
        else:
             tau550=min(meta['tau'])
             print('Warning minimum LUT tau selected.')
             dark_idx = -1
    else:
         if force_band not in band_order:
             print('Band {} not recognised'.format(force_band))
             return()

         dark_idx = [i for i,band in enumerate(band_order) if band==force_band][0]
         tau550 = tau550_all_bands[dark_idx]

    dark_band = band_order[dark_idx]

    ## get different parameters for this TAU550
    ratm = pp.aerlut.interplut_sensor(lut_sensor, meta, azi, thv, ths, tau550, par='romix')
    rorayl = pp.aerlut.interplut_sensor(lut_sensor, meta, azi, thv, ths, tau550, par='rorayl')
    dtotr = pp.aerlut.interplut_sensor(lut_sensor, meta, azi, thv, ths, tau550, par='dtotr')
    utotr = pp.aerlut.interplut_sensor(lut_sensor, meta, azi, thv, ths, tau550, par='utotr')
    dtott = pp.aerlut.interplut_sensor(lut_sensor, meta, azi, thv, ths, tau550, par='dtott')
    utott = pp.aerlut.interplut_sensor(lut_sensor, meta, azi, thv, ths, tau550, par='utott')
    astot = pp.aerlut.interplut_sensor(lut_sensor, meta, azi, thv, ths, tau550, par='astot')

    #return ratm, rorayl, dtotr, utotr, dtott, utott, astot, tau550_all_bands, dark_idx, tau550_all_bands_rmsd
    return ratm, rorayl, dtotr, utotr, dtott, utott, astot, tau_dict, dark_band, tau550_all_bands_rmsd
