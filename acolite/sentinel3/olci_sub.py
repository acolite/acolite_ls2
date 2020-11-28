## def olci_sub
## finds crop in OLCI image
## written by Quinten Vanhellemont, RBINS 2019-04-02
##              QV 2020-11-28 added more accurate limit subsetting

def olci_sub(bundle, limit, use_tpg=True):
    import acolite as ac
    from numpy import where

    if use_tpg:
        file = '{}/{}.nc'.format(bundle, 'geo_coordinates')
        lat = ac.shared.nc_data(file, 'latitude')
        data_shape = lat.shape
        lat = None

        file = '{}/{}.nc'.format(bundle, 'tie_geo_coordinates')
        lat = ac.shared.nc_data(file, 'latitude')
        lon = ac.shared.nc_data(file, 'longitude')
        tpg_shape = lat.shape
    else:
        file = '{}/{}.nc'.format(bundle, 'geo_coordinates')
        lat = ac.shared.nc_data(file, 'latitude')
        lon = ac.shared.nc_data(file, 'longitude')
        data_shape = lat.shape

    ## new version
    tmp = (lat >= limit[0]) & (lat <= limit[2]) & \
          (lon >= limit[1]) & (lon <= limit[3])
    lat = None
    lon = None

    roi = where(tmp)
    tmp = None
    x0, x1 = min(roi[0]),max(roi[0])
    y0, y1 = min(roi[1]),max(roi[1])
    roi = None

    if use_tpg:
        ## to recalculate to full scene!
        #y0 = max(0, y0-1)
        #y1 = min(tpg_shape[1]-1, y1+1)
        ys = int((data_shape[1]-1)/(tpg_shape[1]-1))
        y0s=max(0,y0-1)*ys
        y1s=min(tpg_shape[1],y1+1)*ys
        sub = [y0s, x0, y1s-y0s, x1-x0]
    else:
        sub = [y0, x0, y1-y0, x1-x0]

    return(sub, data_shape)
