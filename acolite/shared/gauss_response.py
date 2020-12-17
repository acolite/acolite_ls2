## def gauss_response
## compute gaussian rsr for center wave and fwhm
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2018-01-30
## modifications: 2020-10-15 (QV) added int to numpy linspace

def gauss_response(center, fwhm, step=1):
    from numpy import linspace, sqrt, log, exp

    wrange = (center - 1.5*fwhm, center + 1.5*fwhm)
    sigma = fwhm / (2*sqrt(2*log(2)))

    x = linspace(wrange[0], wrange[1], int(1+(wrange[1]-wrange[0])/step))
    y = exp(-((x-center)/sigma)**2 )
    return(x,y)
