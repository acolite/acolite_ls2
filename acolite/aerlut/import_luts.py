## read all luts and set up rgi
## QV 2020-01-15
## Last modifications: 2020-07-14 (QV) added sensor option

def import_luts(pressures = [500, 1013, 1100],
                base_luts = ['PONDER-LUT-201704-MOD1', 'PONDER-LUT-201704-MOD2'],
                sensor=None):
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

            if sensor is None:
                ## LUT with 18 monochromatic wavelengths (0.39-2.4)
                lut_data, lut_meta = ac.aerlut.import_lut(lutid,lutdir, override=0)
            else:
                ## sensor specific lut
                lut_data_dict, lut_meta = ac.aerlut.get_sensor_lut(sensor, None, override=0, lutdir=lutdir, lutid=lutid)

                #bands = list(lut_data_dict.keys())
                # get bands from rsr_file as different systems may not keep dict keys in the same order
                rsr_file = ac.config['pp_data_dir']+'/RSR/'+sensor+'.txt'
                rsr, rsr_bands = ac.shared.rsr_read(file=rsr_file)

                ## make array from the individual band results
                for ib, band in enumerate(rsr_bands):
                    if ib == 0:
                        lut_meta['bands'] = []
                        lut_data = np.zeros(list(np.stack(lut_data_dict[band]).shape)+[len(rsr_bands)])
                    lut_meta['bands'].append(ib)
                    lut_data[:,:,:,:,:,:, ib] = lut_data_dict[band]

                del lut_data_dict

                ## make second axis correspond to "wave"
                ## in this case stricly "band"
                lut_data = np.moveaxis(lut_data, -1, 1)

            ## set up lut dimensions
            if ip == 0:
                # lut data is par, wave, azi, thv, ths, wind, tau
                lut_dim = [[i for i in pressures]] + [[i for i in range(len(lut_meta['par']))]]
                if sensor is None:
                    lut_dim += [lut_meta['wave']]
                else:
                    lut_dim += [lut_meta['bands']]
                lut_dim += [lut_meta[k] for k in ['azi', 'thv', 'ths', 'tau']]
                lutd = []
                lut_dict[lut] = {'meta':lut_meta, 'dim':lut_dim}

            lutd.append(lut_data)
        lut_dict[lut]['lut'] = np.stack(lutd)

        ## set up LUT interpolator
        lut_dict[lut]['rgi'] = scipy.interpolate.RegularGridInterpolator(lut_dict[lut]['dim'],
                                                                         lut_dict[lut]['lut'][:,:,:,:,:,:,0,:],
                                                                         bounds_error=False, fill_value=np.nan)
    return(lut_dict)
