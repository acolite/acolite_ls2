## make representative rhod, simplified from select_model.py
## QV 2020-07-15

def rhod_select(rdark, rdark_list_selection='intercept', lowess_frac = 0.5):

    import numpy as np
    import acolite as ac
    from statsmodels.nonparametric.smoothers_lowess import lowess

    ## empty dicts
    rdark_smooth = {}
    rdark_selected = {}

    ## find out what kind of rdark is given and which selection to use
    for band in rdark.keys():
        ## given rdark is a single spectrum
        if (type(rdark[band]) == np.float64):
            rdark_selected[band] = rdark[band]

        ## given rdark is a list of spectra
        elif (type(rdark[band]) == np.ndarray) or (type(rdark[band]) ==  np.ma.MaskedArray):
            pixel_range = [float(i) for i in range(0, len(rdark[band]))]

            ## select rdark by lowess smoothing the given spectra
            if rdark_list_selection == 'lowess':
                rdark_smooth[band] = lowess(rdark[band],pixel_range,frac=lowess_frac)[:,1]
                rdark_selected[band] = rdark_smooth[band][0]
                rdark_list = False

            ## select rdark by OLS intercept
            elif rdark_list_selection == 'intercept':
                for band in rdark.keys():
                    reg = ac.shared.regression.lsqfity(pixel_range, rdark[band]) # m, b, r, sm, sb
                    rdark_selected[band] = reg[1]
                    rdark_list = False

            ## fallback selection is darkest pixel in each band
            else:
                rdark_list_selection = 'darkest'
                rdark_selected[band] = np.nanmin(rdark[band])

        ## given dark is something else
        else:
            print('rdark type not recognised')
            rdark_selected = (rdark)

    return(rdark_selected)
