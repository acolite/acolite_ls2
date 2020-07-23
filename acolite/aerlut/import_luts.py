## read all luts and set up rgi
## QV 2020-01-15
## Last modifications: 2020-07-14 (QV) added sensor option
##                     2020-07-14 (QV) changed sensor option to a dict per band
##                     2020-07-22 (QV) added the rsky lut to the lut and rgi - need to update rsky azimuths

def import_luts(pressures = [500, 1013, 1100],
                base_luts = ['PONDER-LUT-201704-MOD1', 'PONDER-LUT-201704-MOD2'],
                sensor=None, add_rsky=False):
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
            #print(lutdir, lutid)
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


            ## set up lut dimensions
            if ip == 0:
                # lut data is par, wave, azi, thv, ths, wind, tau
                lut_dim = [[i for i in pressures]] + [[i for i in range(len(lut_meta['par']))]]
                if sensor is None:
                    lutd = []
                    lut_dim += [lut_meta['wave']]

                lut_dim += [lut_meta[k] for k in ['azi', 'thv', 'ths', 'tau']]
                ipd = {par:ip for ip,par in enumerate(lut_meta['par'])}

                lut_dict[lut] = {'meta':lut_meta, 'dim':lut_dim, 'ipd':ipd}

                if sensor is not None:
                    lut_dict[lut]['lut'] = {band:[] for band in rsr_bands}

            if sensor is None:
                lutd.append(lut_data)
            else:
                for band in rsr_bands:
                    lut_dict[lut]['lut'][band].append(lut_data_dict[band])


        if sensor is None:
            if add_rsky:
                print('Adding rsky not yet supported for generic LUT')

            lut_dict[lut]['lut'] = np.stack(lutd)
            ## set up LUT interpolator
            lut_dict[lut]['rgi'] = scipy.interpolate.RegularGridInterpolator(lut_dict[lut]['dim'],
                                                                             lut_dict[lut]['lut'][:,:,:,:,:,:,0,:],
                                                                             bounds_error=False, fill_value=np.nan)
        else:
            ## make arrays
            for band in rsr_bands:
                lut_dict[lut]['lut'][band] = np.stack(lut_dict[lut]['lut'][band])

            ## add rsky if requested
            if add_rsky:
                print('Update of rsky azimuths needed!')
                rlut, rmeta, rdim, rrgi = ac.aerlut.rsky_read_lut(int(lut[-1]), sensor=sensor)

                ## current pars
                ipd = {p:i for i,p in enumerate(lut_dict[lut]['meta']['par'])}

                ## run through bands
                for band in rlut:
                    ## add 180 degree azimuth (==0 degree) to end of rsky lut
                    tmp = np.insert(rlut[band], (-1), rlut[band][ 0, :, :, :], axis=0)
                    ## and remove last two vza  of rsky lut
                    tmp = tmp[:,:-2, :,:]
                    ## add empty axis
                    tmp = tmp[:,:,:,None, :]

                    ## add to the LUT
                    ## model rsky at surface
                    ## rsky_s idx 22
                    tlut = np.insert(lut_dict[lut]['lut'][band], (22), tmp, axis=1)

                    ## model rsky at toa
                    ## rsky_t idx 23
                    ## (utott * dtott * rsky) / (1. - rsky * astot)
                    tmp = (tlut[:, ipd['utott'],:,:,:,:,:]*\
                           tlut[:, ipd['dtott'],:,:,:,:,:]*
                           tlut[:, 22,:,:,:,:,:]) /\
                          (1.-tlut[:, 22,:,:,:,:,:] *\
                           tlut[:, ipd['astot'],:,:,:,:,:])
                    tlut = np.insert(tlut, (23), tmp, axis=1)

                    ## add romix+rsky
                    ## idx 24
                    tmp = tlut[:, ipd['romix'],:,:,:,:,:] + tlut[:, 23,:,:,:,:,:]
                    tlut = np.insert(tlut, (24), tmp, axis=1)

                    ## replace in dict
                    lut_dict[lut]['lut'][band] = tlut

                ## add new pars
                lut_dict[lut]['meta']['par'] += ['rsky_s', 'rsky_t', 'romix+rsky_t']
                lut_dict[lut]['dim'][1]+= [22, 23, 24]
                ## end add rsky

            lut_dict[lut]['rgi'] = {}
            lut_dict[lut]['ipd'] = {p:i for i,p in enumerate(lut_dict[lut]['meta']['par'])}
            for band in rsr_bands:
                ## set up LUT interpolator per band
                lut_dict[lut]['rgi'][band] = scipy.interpolate.RegularGridInterpolator(lut_dict[lut]['dim'],
                                                                                       lut_dict[lut]['lut'][band][:,:,:,:,:,0,:],
                                                                                       bounds_error=False, fill_value=np.nan)

    return(lut_dict)
