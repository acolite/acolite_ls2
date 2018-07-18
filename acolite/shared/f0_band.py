## def f0_band
## gets f0 for given wave and rsr
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-01-27
## modifications:
##                2018-07-18 (QV) changed acolite import name

def f0_band(wave, rsr, f0file=None):
    from acolite.shared import f0_get, rsr_convolute
    from numpy import linspace

    f0 = f0_get(f0file=f0file)
    
    band_f0 = rsr_convolute(f0['data'], f0['wave'], rsr, wave)
    return band_f0
