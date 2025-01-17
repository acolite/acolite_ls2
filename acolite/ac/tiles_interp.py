## def tiles_interp
## interpolates tiled dataset to full scene extent
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-12-11
## modifications: 2020-11-17 (QV) added target mask for limiting interpolation extent and processing time
##                2020-11-18 (QV) added dtype to convert from griddata float64, by default float32
##                                this improves peak memory use when several datasets are kept in memory
##

def tiles_interp(data, xnew, ynew, kern_size=2, method='linear', mask=None,
                 target_mask=None, target_mask_full=False, dtype='float32'):

    import numpy as np
    from scipy.interpolate import griddata
    from scipy.ndimage import uniform_filter,percentile_filter, distance_transform_edt

    if mask is not None: data[mask] = np.nan

    ## fill nans with closest value
    ind = distance_transform_edt(np.isnan(data), return_distances=False, return_indices=True)
    cur_data = data[tuple(ind)]

    ## smooth dataset
    z = uniform_filter(cur_data, size=kern_size)
    zv=list(z.ravel())

    dim = data.shape

    ### tile centers
    #x = arange(0.5, dim[1], 1)
    #y = arange(0.5, dim[0], 1)

    ## tile edges
    x = np.arange(0., dim[1], 1)
    y = np.arange(0., dim[0], 1)

    xv, yv = np.meshgrid(x, y, sparse=False)
    ci = (list(xv.ravel()), list(yv.ravel()))

    ## interpolate
    if target_mask is None:
        ## full dataset
        znew = griddata(ci, zv, (xnew[None,:], ynew[:,None]), method=method)
    else:
        ## limit to target mask
        vd = np.where(target_mask)
        if target_mask_full:
            ## return a dataset with the proper dimensions
            znew = np.zeros((len(ynew), len(xnew)))+np.nan
            znew[vd] = griddata(ci, zv, (xnew[vd[1]], ynew[vd[0]]), method=method)
        else:
            ## return only target_mask data
            znew = griddata(ci, zv, (xnew[vd[1]], ynew[vd[0]]), method=method)

    if dtype is None:
        return(znew)
    else:
        return(znew.astype(np.dtype(dtype)))
