## def rsr_convolute
## convolutes dataset to given rsr
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2016-07-11
## modifications: 2017-05-31 (QV) added isfinite check


def rsr_convolute(data, wave, rsr, rsrwave):
    from numpy import interp, where, isfinite
    
    dlen = len(wave)
    drange = [min(wave), max(wave)]
    dstep = (drange[1]-drange[0])/dlen
        
    rlen = len(rsrwave)
    rrange = [min(rsrwave), max(rsrwave)]
    rstep = (rrange[1]-rrange[0])/rlen    
    np = max(dlen, rlen)
    
    sub = where(isfinite(data))

    if len(sub[0]) is not 0:
        data = [data[s] for s in sub[0]]
        wave = [wave[s] for s in sub[0]]

    
    if np == dlen:
        rsr_to_data = interp(wave, rsrwave, rsr, left=0, right=0)
        stotal = sum([value*data[i] for i,value in enumerate(rsr_to_data)])
        ftotal = sum([value for i,value in enumerate(rsr_to_data)])
        
    if np == rlen:
        data_to_rsr = interp(rsrwave, wave, data)
        stotal = sum([value*data_to_rsr[i] for i,value in enumerate(rsr)])
        ftotal = sum([value for i,value in enumerate(rsr)])
    
    band_average = stotal/ftotal
    return band_average
