## fit aerosol model to rho dark
##
## QV 2020-07-29 (rewrite)
##

def rhod_fit_model(rhod, raa, vza, sza,
                   pressure = 1013, fit_tag='romix',
                   lutd = None, sensor=None,
                   luts=['ACOLITE-LUT-202003-MOD1', 'ACOLITE-LUT-202004-MOD2']):

    import acolite as ac
    import numpy as np

    if lutd is None:
        lutd = ac.aerlut.import_luts(base_luts=luts, sensor=sensor, add_rsky=True)
    else:
        luts = list(lutd.keys())
        luts.sort()

    res = {}

    ## run through the luts
    for li, lutid in enumerate(luts):
        res[lutid] = {'tau_all': {},
                      'tau_sel':5.0,
                      'fit_rmsd':5.0,
                      'band_sel':'',
                      'band_sel2':''}

        ## compute aot for each band
        for ib, band in enumerate(rhod):
            ret_ = lutd[lutid]['rgi'][band]((pressure, lutd[lutid]['ipd'][fit_tag],
                                         raa, vza, sza, lutd[lutid]['meta']['tau']))
            res[lutid]['tau_all'][band] = np.interp(rhod[band], ret_, lutd[lutid]['meta']['tau'])

            ## track lowest tau
            if res[lutid]['tau_all'][band] < res[lutid]['tau_sel']:
                res[lutid]['tau_sel'] = res[lutid]['tau_all'][band]
                res[lutid]['band_sel'] = band

        ## compute spectrum for this tau
        res[lutid]['rpath'] = {band:lutd[lutid]['rgi'][band]((pressure, lutd[lutid]['ipd'][fit_tag],
                                    raa, vza, sza, res[lutid]['tau_sel'])) for band in rhod}

        ## compute fits to lowest tau band and the others
        for ib, band in enumerate(rhod):
            if band == res[lutid]['band_sel']: continue
            x = rhod[band], rhod[res[lutid]['band_sel']]
            y = res[lutid]['rpath'][band], res[lutid]['rpath'][res[lutid]['band_sel']]
            cur_rmsd = ac.shared.rmsd(x, y)
            if cur_rmsd < res[lutid]['fit_rmsd']:
                res[lutid]['fit_rmsd'] = cur_rmsd
                res[lutid]['band_sel2'] = band

        #print(lutid, res[lutid]['tau_sel'], res[lutid]['band_sel'], res[lutid]['band_sel2'], '{:.4f}'.format(res[lutid]['fit_rmsd']))
        
        if li == 0:
            sel_model = lutid
            sel_tau = res[lutid]['tau_sel']
        elif res[lutid]['tau_sel'] < sel_tau:
            sel_model = lutid
            sel_tau = res[lutid]['tau_sel']

    return(sel_model, sel_tau, res)
