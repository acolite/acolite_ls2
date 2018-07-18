## def f0_wave
## gets f0 for given wavelength and band width (optional)
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-01-24
## modifications:
##                2018-07-18 (QV) changed acolite import name

def f0_wave(wave, width=1, f0file=None):
    from acolite.shared import f0_get, rsr_convolute
    from numpy import linspace

    f0 = f0_get(f0file=f0file)
    
    off = width/2.
    band_wave=linspace(wave-off,wave+off,width+1)
    band_rsr=[1]*len(band_wave)

    band_f0 = rsr_convolute(f0['data'], f0['wave'], band_rsr, band_wave)
    return band_f0
