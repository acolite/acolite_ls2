## fit model with sky reflectance to rho dark
## QV 2020-03-18
##                2020-06-22 (QV) added rsr keyword and removed band wavelength computation

def rhod_fit_model(raa, vza, sza, pressure = 1013, rhod = None,
                   sensor=None, resample=True, rsr=None,
                   lutd = None, luts=['PONDER-LUT-201704-MOD1', 'PONDER-LUT-201704-MOD2'],
                   rlutd = None):
    import numpy as np
    import acolite as ac

    ### get rhot with sky glint for current model and tau
    def get_cur_rhot(ta):
        cur_rsky = rlutd[aermod]['rgi']((waves, raa, vza, sza, ta))
        cur_rsky[np.isnan(cur_rsky)] = 0

        cur_romix = lutd[lutid]['rgi']((pressure, ipd['romix'], waves, raa, vza, sza, ta))
        cur_utott = lutd[lutid]['rgi']((pressure, ipd['utott'], waves, raa, vza, sza, ta))
        cur_dtott = lutd[lutid]['rgi']((pressure, ipd['dtott'], waves, raa, vza, sza, ta))
        cur_astot = lutd[lutid]['rgi']((pressure, ipd['astot'], waves, raa, vza, sza, ta))
        cur_romix[np.isnan(cur_romix)] = 0

        ## model toa reflectance for this tau
        cur_rhot = cur_romix + (cur_utott * cur_dtott * cur_rsky) / (1. - cur_rsky * cur_astot)
        return(cur_rhot, cur_romix)

    ## import aerosol luts
    if lutd is None:
        lutd = ac.aerlut.import_luts()
    luts = list(lutd.keys())
    mods = [int(l[-1]) for l in luts]

    ## import rsky luts
    if rlutd is None:
        rlutd = {}
        for aermod in mods:
            l, m, d, r = ac.aerlut.rsky_read_lut(aermod)
            rlutd[aermod] = {'lut':l, 'meta':m, 'dim':d, 'rgi':r}

    ## import sensor rsr
    if sensor is not None:
        resample = True
        if rsr is None:
            rsr_file = ac.config['pp_data_dir']+'/RSR/'+sensor+'.txt'
            rsr, rsr_bands = ac.shared.rsr_read(file=rsr_file)
    else:
        resample = False

    ## get results for current geometry
    results = {}
    for il, lutid in enumerate(luts):
        ## set up parameter index dict
        ipd = {par:ip for ip,par in enumerate(lutd[lutid]['meta']['par'])}
        taus = lutd[lutid]['meta']['tau']
        waves = lutd[lutid]['meta']['wave']

        ## model for rsky lut
        aermod=mods[il]

        ## model for each tau
        for it, ta in enumerate(taus):
            cur_rhot, cur_romix = get_cur_rhot(ta)

            if it == 0:
                romix_all = cur_romix
                rhot_all = cur_rhot
            else:
                romix_all = np.vstack((romix_all, cur_romix))
                rhot_all = np.vstack((rhot_all, cur_rhot))

        results[lutid] = {'rhot':rhot_all, 'romix':romix_all, 'taus':taus, 'waves':waves}

    ##
    if resample:
        results_resampled = {}
        for lutid in luts:
            ## run through results
            for it, tau in enumerate(results[lutid]['taus']):
                rs_rhot = ac.shared.rsr_convolute_dict(results[lutid]['waves'], results[lutid]['rhot'][it,:],  rsr)
                rs_romix = ac.shared.rsr_convolute_dict(results[lutid]['waves'], results[lutid]['romix'][it,:],  rsr)

                if it == 0:
                    res = {'tau':[], 'rhot':{b:[] for b in rsr}, 'romix':{b:[] for b in rsr}}
                res['tau'].append(tau)

                for b in rsr:
                    res['rhot'][b].append(rs_rhot[b])
                    res['romix'][b].append(rs_romix[b])

            ## fit here!
            if rhod is not None:
                tau_fits = {}
                tau_fits_nosky = {}
                for b in rhod:
                    tau_fits[b] = np.interp(rhod[b], res['rhot'][b], res['tau'])
                    tau_fits_nosky[b] = np.interp(rhod[b], res['romix'][b], res['tau'])

                min_tau = 5
                min_band = ''
                for b in rhod:
                    if (tau_fits[b] <= min_tau) or (min_band == ''):
                        min_tau = tau_fits[b]
                        min_band = b

                res['tau_fits'] = tau_fits
                res['tau_fit'] = min_tau
                res['tau_band'] = min_band

                 ## model toa reflectance for this tau
                cur_rhot, cur_romix = get_cur_rhot(res['tau_fit'])
                res['rhot_fit'] = cur_rhot
                res['rhot_fit_rs'] = ac.shared.rsr_convolute_dict(results[lutid]['waves'], cur_rhot,  rsr)
                res['romix_fit'] = cur_romix
                res['romix_fit_rs'] = ac.shared.rsr_convolute_dict(results[lutid]['waves'], cur_romix,  rsr)

                ## get best two-band fit
                min_rmsd = 10
                min_rmsd_band = ''
                for b in rhod:
                    if b == res['tau_band']: continue
                    if res['tau_band'] == '': continue
                    x = rhod[b], rhod[res['tau_band']]
                    y = res['rhot_fit_rs'][b], res['rhot_fit_rs'][res['tau_band']]
                    cur_rmsd = ac.shared.rmsd(x, y)
                    if (cur_rmsd < min_rmsd) or (min_rmsd_band == ''):
                        min_rmsd = cur_rmsd
                        min_rmsd_band = (b, res['tau_band'])
                    #print(b, res['tau_band'], cur_rmsd)

                res['rmsd_fit'] = min_rmsd
                res['rmsd_bands'] = min_rmsd_band

            ###
            results_resampled[lutid] = res

        if rhod is not None:
            mod_sel = None
            mod_rmsd = 10

            for lutid in luts:
                if (results_resampled[lutid]['rmsd_fit'] < mod_rmsd) or (mod_sel is None):
                    mod_rmsd=results_resampled[lutid]['rmsd_fit']
                    mod_sel = lutid
            results_resampled['mod_sel'] = mod_sel
            results_resampled['mod_rmsd'] = mod_rmsd

        return(results_resampled)
    else:
        return(results)
