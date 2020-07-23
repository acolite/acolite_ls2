## def gas_transmittance
## gets hyperspectral gas transmittances
## written by Quinten Vanhellemont, RBINS
## 2019-04-02
## modifications: 

def gas_transmittance(sza, vza, waves=None, uoz=0.3, uwv=1.5):
    from numpy import cos, exp, pi, interp
    import acolite as ac

    ko3 = ac.shared.ko3_get()
    if waves is None: waves=ko3['wave']
    else: waves = [float(w/1000) for w in waves]
    
    ## compute water transmittance
    wv_wv_hs, tt_wv_hs = ac.ac.wvlut_interp(sza, vza, uwv=1.5)
    tt_wv = interp(waves, wv_wv_hs, tt_wv_hs)

    ## compute oxygen transmittance 
    wv_o2_hs, tt_o2_hs = ac.ac.o2lut_interp(sza, vza)
    tt_o2 = interp(waves, wv_o2_hs, tt_o2_hs)

    ## cosine of sun and sensor zenith angles
    mu0 = cos(sza*(pi/180))
    muv = cos(vza*(pi/180))

    ## compute ozone transmittance
    koz = interp(waves, ko3['wave'], ko3['data'])
    tau_oz = koz * uoz
    t0_ozone = exp(-1.*(tau_oz) / mu0)
    tv_ozone = exp(-1.*(tau_oz) / muv)
    tt_o3= t0_ozone * tv_ozone
    return({'wave':[1000*w for w in waves], 
            'tt_h2o':tt_wv,
            'tt_o3':tt_o3,
            'tt_o2':tt_o2, 'tt_gas':tt_wv*tt_o3*tt_o2})
