## def lut_get_ac_parameters_fixed_tau_sensor
## returns ac parameters from sensor LUT for given tau, azi, thv and ths (and LUT)
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-01-30
## modifications: 
##                2018-07-18 (QV) changed acolite import name
def lut_get_ac_parameters_fixed_tau_sensor(lut_sensor,meta,azi,thv,ths,tau550):
    from acolite import aerlut

    ## get different parameters for this TAU550
    ratm = aerlut.interplut_sensor(lut_sensor, meta, azi, thv, ths, tau550, par='romix')
    rorayl = aerlut.interplut_sensor(lut_sensor, meta, azi, thv, ths, tau550, par='rorayl')
    dtotr = aerlut.interplut_sensor(lut_sensor, meta, azi, thv, ths, tau550, par='dtotr')
    utotr = aerlut.interplut_sensor(lut_sensor, meta, azi, thv, ths, tau550, par='utotr')
    dtott = aerlut.interplut_sensor(lut_sensor, meta, azi, thv, ths, tau550, par='dtott')
    utott = aerlut.interplut_sensor(lut_sensor, meta, azi, thv, ths, tau550, par='utott')
    astot = aerlut.interplut_sensor(lut_sensor, meta, azi, thv, ths, tau550, par='astot')
    
    return (ratm, rorayl, dtotr, utotr, dtott, utott, astot)
