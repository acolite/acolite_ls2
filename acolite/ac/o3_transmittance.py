# function o3_transmittance
# computes ozone transmittance for sensor rsr
# uoz = ozone in cm-1, 0.3 cm-1 = 300 DU
# QV 2017-04-14
##                2017-11-28 (QV) moved PP data directory
##                2018-07-18 (QV) changed acolite import name
def o3_transmittance(sensor, metadata, uoz=0.3, rsr_file=None):
    from numpy import cos, exp, pi
    from acolite.shared import rsr_read, ko3_band
    
    ## get sensor rsr
    import os
    ## set rsr_file
    from acolite import config
    pp_path = config['pp_data_dir']
    if rsr_file is None:
        rsr_file = pp_path+'/RSR/'+sensor+'.txt'
    rsr, rsr_bands = rsr_read(file=rsr_file)
    
    ## cosine of sun and sensor zenith angles
    mu0 = cos(metadata['THS']*(pi/180))
    muv = cos(metadata['THV']*(pi/180))
    
    ## compute ozone transmittance
    tt_oz={}
    for btag in rsr_bands:
        koz = ko3_band(rsr[btag]['wave'],rsr[btag]['response'])
        tau_oz = koz * uoz
        t0_ozone = exp(-1.*(tau_oz) / mu0)
        tv_ozone = exp(-1.*(tau_oz) / muv)
        tt_oz[btag] = t0_ozone * tv_ozone
    return(tt_oz)
