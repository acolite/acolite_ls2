## def get_ll
## gets lon and lat for given landsat metadata with optional limit
## written by Quinten Vanhellemont, RBINS
## 2017-04-13
## modifications:
##                2017-11-27 (QV) added extent_limit keyword: return larger lat/lon if requested limit extends over the scene
##                2018-06-06 (QV) changed xy output
##                2018-07-18 (QV) changed acolite import name
##                2018-10-01 (QV) added grid cell size option

def get_ll(metadata, limit=None, xy=False, extend_limit=False):
    from numpy import linspace, tile, flipud
    from acolite.landsat.geo import get_sub
    from acolite.landsat.geo import get_projection

    if 'GRID_CELL_SIZE_REFLECTIVE' in metadata:
        pixelsize = [float(metadata["GRID_CELL_SIZE_REFLECTIVE"])]*2
    elif 'GRID_CELL_SIZE_REF' in metadata:
        pixelsize = [float(metadata["GRID_CELL_SIZE_REF"])]*2
    else:
        return(1)

    if limit is not None:
        sub, p, (xrange, yrange, grid_region), proj4_string = get_sub(metadata, limit)
        if extend_limit is False:
            dims = [sub[2],sub[3]]
        else:
            dims = grid_region['dims']
            xrange, yrange = grid_region['xrange'], grid_region['yrange']
    else:
        p, (xrange,yrange), proj4_string = get_projection(metadata)
        dims = metadata["DIMS"]

    ## no mid pixel for Landsat?
    xdim =  linspace(xrange[0],xrange[1],dims[0]).reshape(1,dims[0])
    ydim =  linspace(yrange[0],yrange[1],dims[1]).reshape(dims[1],1)

    xdim = tile(xdim, (dims[1],1))
    ydim = flipud(tile(ydim, (1,dims[0])))

    if xy:
        return(xdim,ydim)
    else:
        lon,lat = p(xdim,ydim,inverse=True)
        xdim = None
        ydim = None
        return(lon,lat)
