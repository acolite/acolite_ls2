## def lut_get_taufit_sensor
## returns tau/rmsd/band from sensor LUT for given rtoa list, azi, thv and ths
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-04-27
## modifications: 2017-10-19 QV added force_band option
##                2018-07-18 (QV) changed acolite import name
def lut_get_taufit_sensor(lut_sensor,meta,azi,thv,ths,rtoa, bestfit_bands=None, force_band=None):
    import acolite as pp
    from numpy import nanmin, nan, atleast_1d
    
    ## calculate tau for each band observation
    tau_dict = pp.aerlut.lut_get_tau_sensor(lut_sensor,meta,azi,thv,ths,rtoa)

    ## get band wavelengths and order on wavelength
    band_wl = pp.aerlut.interplut_sensor(lut_sensor, meta, azi, thv, ths, 0., par='wl')
    bo = sorted(band_wl, key=band_wl.get)
    band_order = [b for i,b in enumerate(bo) if b in rtoa]
    fit_bands = band_order if bestfit_bands is None else bestfit_bands
    
    ## extract from dict
    band_wave = [band_wl[b] for i,b in enumerate(band_order)]

    ## get fit for all bands
    tau_rmsd_dict = {}
    for band in band_order:
        if band in rtoa:
            ## use precomputed sensor specific LUT -- faster
            tau550band = []
            if len(atleast_1d(tau_dict[band])) == 1:
                ratm_ = pp.aerlut.interplut_sensor(lut_sensor, meta, azi, thv, ths, tau_dict[band], par='romix')
                tau550band.append(pp.rmsd([rtoa[b] for b in fit_bands], [ratm_[b] for b in fit_bands]))
            else:
                for ti,tau_i in enumerate(tau_dict[band]):
                    ratm_ = pp.aerlut.interplut_sensor(lut_sensor, meta, azi, thv, ths, tau_i, par='romix')
                    tau550band.append(pp.rmsd([rtoa[b][ti] for b in fit_bands], [ratm_[b] for b in fit_bands]))
            tau_rmsd_dict[band] = (tau550band)
    
    ndark = len(atleast_1d(tau_dict[band_order[0]]))
    
    if ndark > 1:
        ## select best band (lowest tau) for each pixel in rtoa (rdark) dictionary
        tau_band = []
        tau_550 = []
        tau_rmsd = []

        for i in range(0,ndark):
            tau_pix = [tau_dict[band][i] for band in band_order]
            tau_pix = [i if (i > min(meta['tau'])) else nan for i in tau_pix]
            rmsd_pix = [tau_rmsd_dict[band][i] for band in band_order]
            if force_band is None:
                min_tau = nanmin(tau_pix)
                dark_idx = [idx for idx,val in enumerate(tau_pix) if val == min_tau]
                dark_idx = 0 if len(dark_idx) == 0 else dark_idx[0]
            else:
                if force_band not in band_order:
                    print('Band {} not recognised'.format(force_band))
                    return()
                dark_idx = [i for i,band in enumerate(band_order) if band==force_band][0]

            tau_band.append(dark_idx)
            tau_550.append(tau_pix[dark_idx])
            tau_rmsd.append(rmsd_pix[dark_idx])
    else:
        tau_pix = [tau_dict[band] for band in band_order]
        tau_pix = [i if (i > min(meta['tau'])) else nan for i in tau_pix]
        rmsd_pix = [tau_rmsd_dict[band] for band in band_order]
        if force_band is None:
            min_tau = nanmin(tau_pix)
            dark_idx = [idx for idx,val in enumerate(tau_pix) if val == min_tau]#[0]
            dark_idx = 0 if len(dark_idx) == 0 else dark_idx[0]
        else:
            if force_band not in band_order:
                print('Band {} not recognised'.format(force_band))
                return()
            dark_idx = [i for i,band in enumerate(band_order) if band==force_band][0]

        tau_band = dark_idx
        tau_550 = tau_pix[dark_idx]
        tau_rmsd = rmsd_pix[dark_idx]
        
    return tau_550, tau_rmsd, tau_band
