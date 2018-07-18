## def rtoa_to_rhos
## convert TOA reflectance to surface reflectance
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2016-07-05
## modifications:

def rtoa_to_rhos(rtoa, ratm, tu, td, sa, tt_gas = 1.):
    rtoa_noatm = (rtoa/tt_gas) - ratm
    rhos = (rtoa_noatm) / (tu * td + sa * rtoa_noatm)
    return rhos
