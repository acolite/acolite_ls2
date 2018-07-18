## def interplut_sensor
## interpolates sensor LUT to given azi/thv/ths and tau (550)
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2016-07-13
## modifications:
##                2018-07-18 (QV) changed acolite import name

def interplut_sensor(lut, meta, azi, thv, ths, tau, band=None, par='romix'):
    from acolite.aerlut import lutpos, interp3d

    #nwave = len(meta['wave'])
    ntau = len(meta['tau'])
    nazi = len(meta['azi'])
    nthv = len(meta['thv'])
    nths = len(meta['ths'])

    azi_id, azi_br = lutpos(meta['azi'], azi)
    thv_id, thv_br = lutpos(meta['thv'], thv)
    ths_id, ths_br = lutpos(meta['ths'], ths)

    tau_id, tau_br = lutpos(meta['tau'], tau)
    tau_w = tau_id - tau_br[0]

    win_id = 0
    par_id = [i for i,value in enumerate(meta['par']) if value == par][0]
  
    ## adapted to sensor LUT
    blist = lut.keys()

    rinti = dict()
    for bd in blist:
        if bd in lut.keys():
            it1 = interp3d(lut[bd][par_id,:,:,:,win_id, tau_br[0]], azi_id, thv_id, ths_id)
            it2 = interp3d(lut[bd][par_id,:,:,:,win_id, tau_br[1]], azi_id, thv_id, ths_id)
            rinti[bd] = it2 * tau_w + it1 * (1.-tau_w)

    if band != None: 
        return rinti[band]
    else: return rinti
