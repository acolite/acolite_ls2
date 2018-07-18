## def get_ll
## gets lon and lat for given Sentinel granule metadata with optional limit (defaults to 10m)
## written by Quinten Vanhellemont, RBINS
## 2017-04-18
## modifications: 2017-04-19 (QV) fixed dims problem with full tile data
##                2017-11-27 (QV) added extent_limit keyword: return larger lat/lon if requested limit extends over the scene
##                2018-04-17 (QV) removed flipud for ydim in case of full tile
##                2018-05-15 (QV) added res string
##                2018-06-06 (QV) changed xy output
##                2018-07-18 (QV) changed acolite import name

def get_ll(metadata, limit=None, xy=False, resolution='10', extend_limit=False):
    from numpy import linspace, tile, flipud
    from acolite.sentinel.geo import get_sub
    from acolite.sentinel.geo import get_projection

    res = str(resolution)
    pixelsize = [float(resolution)]*2

    if limit is not None:
        grids, proj4_string = get_sub(metadata, limit)
        if extend_limit is False:
            sub, p, xrange, yrange = grids[res]['sub'], grids[res]['p'], grids[res]['xrange'], grids[res]['yrange'] 
            dims = [sub[2],sub[3]]
        else:
            ## compute full scene lat/lon
            sub, p, xrange, yrange = grids['grids_region'][res]['sub'], \
                                     grids['grids_region'][res]['p'], \
                                     grids['grids_region'][res]['xrange'], \
                                     grids['grids_region'][res]['yrange']
            dims = [sub[2],sub[3]]
    else:
        p,grids,proj4_string = get_projection(metadata)
        xrange, yrange = grids[res]['xrange'], grids[res]['yrange'] 
        dims = [metadata["GRIDS"][res]['NCOLS'],metadata["GRIDS"][res]['NROWS']]
        

    # testing 2018 06 06
    midpix = int(int(res)/2)
    xdim =  linspace(xrange[0]+midpix,xrange[1]-midpix,dims[0]).reshape(1,dims[0])
    ydim =  linspace(yrange[0]+midpix,yrange[1]-midpix,dims[1]).reshape(dims[1],1)

    xdim = tile(xdim, (dims[1],1))
    ydim = tile(ydim, (1,dims[0]))
    if limit is not None:
        ydim = flipud(ydim)
    
    if xy:
        return(xdim,ydim)
    else:
        lon,lat = p(xdim,ydim,inverse=True)
        xdim = None
        ydim = None
        return(lon,lat)
