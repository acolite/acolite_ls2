## to convolute an n-dimensional array (d1, d2, ..., dn, wave) to a band averaged array
## axis gives the position of the wavelength dimension
## QV 2020-07-14

def rsr_convolute_nd(data, wave, response, response_wave, axis=2):
    import numpy as np
    from scipy.interpolate import interp1d

    f = interp1d(wave, data, axis=axis)
    data_i = f(response_wave)

    # put the wavelength axis in front
    shape = np.swapaxes(data_i, data.ndim-1, axis).shape

    # broadcast the response to the same shape
    response_b = np.broadcast_to(response, shape)

    # put the wavelength axis back to the end
    response_b = np.swapaxes(response_b, data_i.ndim-1, axis)

    ## return band weighted data
    return(np.nansum(data_i * response_b, axis=axis) / np.nansum(response))
