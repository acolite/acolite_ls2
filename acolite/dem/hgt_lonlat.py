## def hgt_lonlat
## gets SRTM DEM for given lon/lat arrays
## written by Quinten Vanhellemont, RBINS
## QV 2017-07-17
## modifications:
## updated 

def hgt_lonlat(lon1, lat1, nearest=True, hgt_dir=None):
    if hgt_dir is None:
        import acolite as ac
        hgt_dir = ac.config['hgt_dir']

    from acolite.dem import hgt_find,hgt_read,hgt_geolocation
    from acolite.shared import reproject2

    ## find dem files
    limit=[0,0,0,0]
    if type(lon1) is float:
        limit[1]=lon1
        limit[3]=lon1
    else:
        limit[1]=lon1.min()
        limit[3]=lon1.max()

    if type(lat1) is float:
        limit[0]=lat1
        limit[2]=lat1
    else:
        limit[0]=lat1.min()
        limit[2]=lat1.max()

    hgt_files, hgt_required = hgt_find(limit, required=True, hgt_dir=hgt_dir)

    ## run through dem files and reproject data to target lat,lon
    for i, hgt_file in enumerate(hgt_files):
        ## read hgt data and geolocation
        hgt = hgt_read(hgt_file)

        if (type(lon1) is float) & (type(lat1) is float):
            lon0,lat0 = hgt_geolocation(hgt_file, grid=False)
            from scipy import interpolate
            hgtip = interpolate.RectBivariateSpline(lon0,lat0, hgt)
            result = hgtip(lon1,lat1)
        else:
            lon0,lat0 = hgt_geolocation(hgt_file, grid=True)
            ## reproject
            result = reproject2(hgt, lon0, lat0, lon1, lat1, nearest=nearest)
        
        ## make output
        if i == 0:
            dem = result
        else:
            ## copy new results where != 0 (pyresample fill does not seem to be working)
            dem[result != 0] = result[result != 0]
            
    return(dem)
