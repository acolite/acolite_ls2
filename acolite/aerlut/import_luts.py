## read all luts and set up rgi
## QV 2020-01-15
def import_luts(pressures = [500, 1013, 1100], base_luts = ['PONDER-LUT-201704-MOD1', 'PONDER-LUT-201704-MOD2']):
    import scipy.interpolate
    import numpy as np
    import acolite as ac

    lut_dict = {}
    ## run through luts
    for lut in base_luts:
        ## run through pressures
        for ip, pr in enumerate(pressures):
            lutid = '{}-{}mb'.format(lut, '{}'.format(pr).zfill(4))

            lutdir = '{}/{}'.format(ac.config['pp_data_dir'], 'LUT')
            lut_data, lut_meta = ac.aerlut.import_lut(lutid,lutdir, override=0)

            ## set up lut dimensions
            if ip == 0:
                # lut data is par, wave, azi, thv, ths, wind, tau
                lut_dim = [[i for i in pressures]] + [[i for i in range(len(lut_meta['par']))]]
                lut_dim += [lut_meta[k] for k in ['wave', 'azi', 'thv', 'ths', 'tau']]
                lutd = []
                lut_dict[lut] = {'meta':lut_meta, 'dim':lut_dim}

            lutd.append(lut_data)
        lut_dict[lut]['lut'] = np.stack(lutd)

        ## set up LUT interpolator
        lut_dict[lut]['rgi'] = scipy.interpolate.RegularGridInterpolator(lut_dict[lut]['dim'],
                                                                         lut_dict[lut]['lut'][:,:,:,:,:,:,0,:],
                                                                         bounds_error=False, fill_value=np.nan)
    return(lut_dict)
