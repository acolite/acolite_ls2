## def select_rdark
## 
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-12-06
## modifications: QV  2017-12-11 added pixel range check
##                2018-07-18 (QV) changed acolite import name

def select_rdark(data, rdark_list_selection='intercept', pixel_range=[0,1000], 
                 lowess_frac = 0.5, pixel_idx=100):
    from numpy import isfinite, sort,linspace, nanmin
    from statsmodels.nonparametric.smoothers_lowess import lowess
    import acolite as pp
    
    dsorted = sort(data[isfinite(data)], axis=None)
    pixel_range[1] = nanmin((len(dsorted), pixel_range[1]))
    rdark_list = dsorted[pixel_range[0]:pixel_range[1]]

    if rdark_list_selection=='intercept':
        xi = linspace(pixel_range[0], pixel_range[1], num=pixel_range[1]-pixel_range[0])
        m, b, r, sm, sb = pp.shared.regression.lsqfity(xi,rdark_list)
        rdark_sel = b
    elif rdark_list_selection=='smooth':
        xi = linspace(pixel_range[0], pixel_range[1], num=pixel_range[1]-pixel_range[0])
        rdark_smooth = lowess(rdark_list,xi, frac=lowess_frac)[:,1]
        rdark_sel = rdark_smooth[0]
    elif rdark_list_selection=='absolute_index':
        rdark_sel = dsorted[pixel_idx]
    else:
        rdark_sel = dsorted[0]
    
    return(rdark_sel)
