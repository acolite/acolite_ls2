## def get_projection
## get projection from Landsat metadata
##
## written by Quinten Vanhellemont, RBINS 
## 2017-04-13
## modifications:  2017-10-25 (QV) moved to Proj4 strings and added Polar Stereographic projection
##                 2018-06-06 (QV) added return of Proj4 string

def get_projection(metadata):
    #import osr
    from pyproj import Proj

    if 'GRID_CELL_SIZE_REFLECTIVE' in metadata:
        pixelsize = [float(metadata["GRID_CELL_SIZE_REFLECTIVE"])]*2
    elif 'GRID_CELL_SIZE_REF' in metadata:
        pixelsize = [float(metadata["GRID_CELL_SIZE_REF"])]*2
    else:
        return(1)

    proj = metadata['MAP_PROJECTION']
    if 'ELLIPSOID' in metadata:
        ellipsoid = metadata['ELLIPSOID']
    elif 'REFERENCE_ELLIPSOID' in metadata:
        ellipsoid = metadata['REFERENCE_ELLIPSOID']
    else:
        return(1)

    if 'DATUM' in metadata:
        datum = metadata['DATUM']
    elif 'REFERENCE_DATUM' in metadata:
        datum = metadata['REFERENCE_DATUM']
    else:
        return(1)

    if (proj == 'UTM') & (ellipsoid == 'WGS84') & (datum == 'WGS84'): is_utm = True
    else: is_utm = False
        
    if (proj == 'PS') & (ellipsoid == 'WGS84') & (datum == 'WGS84'): is_ps = True
    else: is_ps = False

    if (not is_utm) & (not is_ps): print('Projection not implemented')
    
    if is_utm:
        zone = int(metadata['UTM_ZONE'])

        #p = Proj(proj='utm',zone=zone,ellps=ellipsoid)
        proj4_list = ['+proj=utm',
                      '+zone={}'.format(zone),
                      '+datum={}'.format(datum),
                      '+units=m',
                      '+no_defs ']
        
    if is_ps:
        vertical_lon = float(metadata["VERTICAL_LON_FROM_POLE"])
        lat_ts = float(metadata["TRUE_SCALE_LAT"])
        false_e = int(metadata["FALSE_EASTING"])
        false_n = int(metadata["FALSE_NORTHING"])
         
        lat_0 = -90. if lat_ts < 0 else 90.
        proj4_list =['+proj=stere',
                     '+lat_0={}'.format(lat_0),
                     '+lat_ts={}'.format(lat_ts),
                     '+lon_0={}'.format(vertical_lon),
                     '+k=1',
                     '+x_0={}'.format(false_e),
                     '+y_0={}'.format(false_n),
                     '+datum={}'.format(datum),
                     '+units=m',
                     '+no_defs ']
    
    proj4_string = ' '.join(proj4_list)
    p = Proj(proj4_string)
    
    ## check corners of image
    x,y = [],[]
    for corner in ['LL','UL','UR','LR']:
            xtag = 'CORNER_{}_PROJECTION_X_PRODUCT'.format(corner)
            if xtag in metadata: 
                x.append(metadata[xtag])
            else:
                xtag = 'PRODUCT_{}_CORNER_MAPX'.format(corner)
                x.append(metadata[xtag])
            ytag = 'CORNER_{}_PROJECTION_Y_PRODUCT'.format(corner)
            if ytag in metadata: 
                y.append(metadata[ytag])
            else:
                ytag = 'PRODUCT_{}_CORNER_MAPY'.format(corner)
                y.append(metadata[ytag])

    xrange = [min(x),max(x)]
    yrange = [min(y),max(y)]

    return(p, (xrange,yrange), proj4_string)
