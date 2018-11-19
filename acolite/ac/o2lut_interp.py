## o2lut_interp
## gives O2 transmittance for given sun and view zenith angles, sensor
##
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2018-10-30
## modifications:

def o2lut_interp(ths, thv, sensor=None, o2config='201810C', par_id = 2):
    import os
    import acolite as ac
    
    o2lut, o2meta = ac.ac.o2lut_get(config=o2config)
    
    ## find position of given angles in LUT
    ths_id, ths_br = ac.aerlut.lutpos(o2meta['ths'], ths)
    thv_id, thv_br = ac.aerlut.lutpos(o2meta['thv'], thv)
    
    ## interpolate hyperspectral dataset
    iw = ac.aerlut.interp2d(o2lut[:,:,par_id,:], ths_id, thv_id)

    if sensor is None: 
        ## return hyperspectral dataset for this geometry
        return(o2meta["wave"], iw)
    else:
        # find RSR
        pp_path = ac.config['pp_data_dir']
        rsr_file = pp_path+'/RSR/{}.txt'.format(sensor)
        rsr,bands = ac.shared.rsr_read(file=rsr_file)

        ## make band averaged values
        band_averaged = ac.shared.rsr_convolute_dict(o2meta['wave'], iw, rsr)
        return(band_averaged)
