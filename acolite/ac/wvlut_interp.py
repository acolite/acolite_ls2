## wvlut_interp
## gives WV transmittance for given sun and view zenith angles, wv and sensor
## uwv = water vapour in g/cm2
##
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-10-17
## modifications: 2017-10-18 (QV) added default uwv
##                2017-11-28 (QV) moved PP data directory
##                2018-01-31 (QV) fixed return when no sensor is given
##                2018-07-18 (QV) changed acolite import name

def wvlut_interp(ths, thv, uwv=1.5, sensor=None, config='201710C', par_id = 2):
    import os
    import acolite as pp
    
    wvlut, wvmeta = pp.ac.wvlut_get(config=config)
    
    ## find position of given angles in LUT
    ths_id, ths_br = pp.aerlut.lutpos(wvmeta['ths'], ths)
    thv_id, thv_br = pp.aerlut.lutpos(wvmeta['thv'], thv)
    wv_id, wv_br = pp.aerlut.lutpos(wvmeta['wv'], uwv)

    ## interpolate hyperspectral dataset
    iw = pp.aerlut.interp3d(wvlut[:,:,:,par_id,:], ths_id, thv_id, wv_id)

    if sensor is None: 
        ## return hyperspectral dataset for this geometry
        return(wvmeta['wave'], iw)
    else:
        # find RSR
        pp_path = pp.config['pp_data_dir']
        rsr_file = pp_path+'/RSR/{}.txt'.format(sensor)
        rsr,bands = pp.shared.rsr_read(file=rsr_file)

        ## make band averaged values
        band_averaged = pp.shared.rsr_convolute_dict(wvmeta['wave'], iw, rsr)
        return(band_averaged)
