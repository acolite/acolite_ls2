## computes analytically the sky reflectance contribution from Fresnel equations
## use band weighted Rayleigh optical thickness from 6SV LUTs
##
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-11-14
##
##                2017-11-28 (QV) moved PP data directory
##                2018-01-22 (QV) added satellite sensor check, and thv check
##                2018-03-14 (QV) added lut sensor check
##                2018-07-18 (QV) changed acolite import name

def toa_rsky(metadata, lut=None, lutdir=None, rsr_file=None, 
             sel_model_lut=None,sel_model_lut_meta=None, pressure=None, tau550=0.01):
    import acolite as pp
    import os

    if 'LUT_SENSOR' in metadata.keys():
        sensor = metadata['LUT_SENSOR']
    elif 'SATELLITE_SENSOR' in metadata.keys():
        sensor = metadata['SATELLITE_SENSOR']
    else:
        sensor = metadata['SENSOR']
    
    if sel_model_lut==None or sel_model_lut_meta==None:

        ## set LUT dir and rsr_file
        from acolite import config
        pp_path = config['pp_data_dir']
        if lutdir is None:
            lutdir=pp_path+'/LUT/'
        if rsr_file is None:
            rsr_file = pp_path+'/RSR/'+sensor+'.txt'

        if lut is None: lut = ['PONDER-LUT-201704-MOD1-1013mb', 'PONDER-LUT-201704-MOD2-1013mb', 'PONDER-LUT-201704-MOD3-1013mb'][1]
        if pressure is not None:
            ## interpolate LUTs to given pressure
            sel_model_lut, sel_model_lut_meta = pp.aerlut.aerlut_pressure(lut, lutdir, pressure, sensor, rsr_file)
        else:
            sel_model_lut, sel_model_lut_meta = pp.aerlut.get_sensor_lut(sensor, rsr_file, lutdir=lutdir, lutid=lut, override=0)

    thv = max(metadata['THV'],0.001)
    ## get Rayleigh optical thickness from LUT
    tray = pp.aerlut.interplut_sensor(sel_model_lut,sel_model_lut_meta,
                                      metadata['AZI'], thv, metadata['THS'], 
                                       tau550, par='tray')
    
    ## set wl to 550 nm (ignored anyway by Rayleigh computation if Rayleigh optical thickness is given)
    wl = {btag:0.55 for btag in tray.keys()}
    
    ## to convert angles in degrees in metadata to the radians used in the Rayleigh functions
    from numpy import pi
    dtor = pi / 180.
    
    rsky = {}
    for btag in tray.keys():
        if btag == 'Pan': continue       
        only_sky = pp.ac.rayleigh.ray_refl_onlysky(wl[btag], metadata['THS']*dtor, thv*dtor, 0, metadata['AZI']*dtor, 
                                           Patm=pressure, tau_ray=tray[btag])        
        rsky[btag] = only_sky
    return(rsky)
