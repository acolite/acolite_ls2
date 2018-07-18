## def exponential_epsilon
## does exponential extrapolation of epsilon for given wavelengths
## idx1 and idx2 are the bands for which the ratio is given
## optional idxc is the desired wavelength (will be set automatically if a 2D epsilon is given)
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-06-20
## modifications: 

def exponential_epsilon(epsilon, waves, idx1=None, idx2=None, idxc=None):
    from numpy import asarray, power
    eps = asarray(epsilon)
    wavelengths=asarray(waves)
    
    if idx1 is None: idx1 = len(wavelengths)-2
    if idx2 is None: idx2 = len(wavelengths)-1
    
    if (eps.shape == ()) or (eps.shape == (1,)):
        if (eps.shape == (1,)): eps=eps[0]
    else:
        if idxc is None:
            idxc = len(wavelengths)-3
                    
    ldiff = (wavelengths[idx2] - wavelengths)
    if idxc is None:
        delta = ldiff / ldiff[idx1]
    else:
        delta = ldiff[idxc] / ldiff[idx1]
    epsall = power(eps,delta)
    return(epsall)
