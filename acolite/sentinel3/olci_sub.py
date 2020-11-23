## def olci_sub
## finds crop in OLCI image
## written by Quinten Vanhellemont, RBINS 2019-04-02

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
        
    ## get corners
    corner = {}
    for c in ['LL', 'UL', 'UR', 'LR']:
        if c[0]=='L':
            st_lat = limit[0]
        if c[0]=='U':
            st_lat = limit[2]

        if c[1]=='L':
            st_lon = limit[1]
        if c[1]=='R':
            st_lon = limit[3]

        londiff = abs(lon-st_lon)    
        latdiff = abs(lat-st_lat)

        diff = pow(pow(londiff,2)+pow(latdiff,2),0.5)
        x,y = where(diff == diff.min())
        corner[c]={'x':x[0], 'y':y[0]}
    diff = None
    londiff = None
    latdiff = None
    lat = None
    lon = None
    
    x0 = min(corner['LL']['x'], corner['UL']['x'])
    x1 = max(corner['LR']['x'], corner['UR']['x'])

    y0 = min(corner['LL']['y'], corner['LR']['y'])
    y1 = max(corner['UL']['y'], corner['UR']['y'])
    
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
