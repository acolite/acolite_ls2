## def tiles_interp
## interpolates tiled dataset to full scene extent
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-12-11
## modifications: 

def tiles_interp(data, xnew, ynew, kern_size=2, method='linear', mask=None):
    from numpy import nan, isnan, meshgrid, arange
    from scipy.interpolate import griddata
    from scipy.ndimage import uniform_filter,percentile_filter, distance_transform_edt

    if mask is not None: data[mask] = nan
    
    ## fill nans with closest value
    ind = distance_transform_edt(isnan(data), return_distances=False, return_indices=True)
    cur_data = data[tuple(ind)]

    ## smooth dataset
    z = uniform_filter(cur_data, size=kern_size)
    zv=list(z.ravel())

    dim = data.shape

    ### tile centers
    #x = arange(0.5, dim[1], 1)
    #y = arange(0.5, dim[0], 1)
                
    ## tile edges
    x = arange(0., dim[1], 1)
    y = arange(0., dim[0], 1)
    
    xv, yv = meshgrid(x, y, sparse=False)
    xv=list(xv.ravel())
    yv=list(yv.ravel())

    ## interpolate
    znew = griddata((xv, yv), zv, (xnew[None,:], ynew[:,None]), method=method)
    return(znew)
