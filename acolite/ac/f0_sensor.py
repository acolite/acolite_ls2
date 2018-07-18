# function f0_sensor
# computes F0 for sensor rsr
# QV 2017-06-06
##                2017-11-28 (QV) moved PP data directory
##                2018-07-18 (QV) changed acolite import name
def f0_sensor(sensor, rsr_file=None):
    from numpy import cos, exp, pi
    from acolite.shared import rsr_read, f0_band
    import acolite as pp

    ## get sensor rsr
    import os
    
    ## set rsr_file
    from acolite import config
    pp_path = config['pp_data_dir']
    if rsr_file is None:
        rsr_file = pp_path+'/RSR/'+sensor+'.txt'
    rsr, rsr_bands = rsr_read(file=rsr_file)

    ## compute F0
    f0={}
    for btag in rsr_bands:
        wave = [w*1000 for w in rsr[btag]['wave'] ]
        f0[btag] = f0_band(wave,rsr[btag]['response'])*10.

    return(f0)
