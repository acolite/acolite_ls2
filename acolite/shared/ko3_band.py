## def ko3_band
## gets ko3 for given wave and rsr
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-01-27
## modifications:
##                2018-07-18 (QV) changed acolite import name

def ko3_band(wave, rsr, ko3file=None):
    from acolite.shared import ko3_get, rsr_convolute
    from numpy import linspace

    ko3 = ko3_get(ko3file=ko3file)
    
    band_ko3 = rsr_convolute(ko3['data'], ko3['wave'], rsr, wave)
    return band_ko3
