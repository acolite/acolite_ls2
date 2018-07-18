## def lut_get_tau_sensor
## gets tau for observed rtoa (dict with same band names as LUT) and sensor LUT (and azi/thv/ths)
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2016-07-13
## modifications: 2017-04-27 (QV) added option for more than one rtoa element per band
##                2018-07-18 (QV) changed acolite import name
def lut_get_tau_sensor(lut_sensor,meta,azi,thv,ths,rtoa):
    from numpy import interp, atleast_1d
    from acolite.aerlut import interplut_sensor, lutpos
    
    tau = dict()
    for band in rtoa:
        if band in lut_sensor:
            ratm = [interplut_sensor(lut_sensor, meta, azi, thv, ths, tau550, band, par='romix') for tau550 in meta['tau']]
            if len(atleast_1d(rtoa[band])) == 1:
                 tau_id, tau_br = lutpos(ratm, rtoa[band]) 
                 tau[band] = interp(tau_id, range(len(meta['tau'])),meta['tau'])
            else:
                 tauband=[]
                 for i in rtoa[band]:
                     tau_id, tau_br = lutpos(ratm, i) 
                     tauband.append(interp(tau_id, range(len(meta['tau'])),meta['tau']))
                 tau[band] = tauband
    return tau
