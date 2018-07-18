## def get_projection
## get projection from Sentinel granule metadata
##
## written by Quinten Vanhellemont, RBINS 
## 2017-04-18
## modifications: 2018-03-07 QV changed to Proj4 string to avoid _gdal import problem in PyInstaller
##                 2018-06-06 (QV) added return of Proj4 string

def get_projection(metadata):
    from pyproj import Proj

    proj = 'UTM'
    ellipsoid = 'WGS84'
    datum = 'WGS84'
    cs_code = metadata["HORIZONTAL_CS_CODE"]
    cs_name = metadata["HORIZONTAL_CS_NAME"]
    epsg = int(cs_code.split(':')[1])
    

    split = cs_name.split('/')
    datum = split[0].strip()
    zone_name = split[1].split()[-1]

    
    datum = 'WGS84'
    if 32600 < epsg <= 32660:
        zone = epsg - 32600
        proj4_list = ['+proj=utm',
                      '+zone={}'.format(zone),
                      '+datum={}'.format(datum),
                      '+units=m',
                      '+no_defs ']

    if 32700 < epsg <= 32760:
        zone = epsg - 32700
        proj4_list = ['+proj=utm',
                      '+zone={}'.format(zone),
                      '+south',
                      '+datum={}'.format(datum),
                      '+units=m',
                      '+no_defs ']

    proj4_string = ' '.join(proj4_list)
    p = Proj(proj4_string)
    
    grids = {}
    for i in metadata['GRIDS'].keys():
        grid=metadata['GRIDS'][i]
        x0 = float(grid['ULX'])
        y0 = float(grid['ULY'])
        xs = float(grid['XDIM'])
        ys = float(grid['YDIM'])
        nx = float(grid['NCOLS'])
        ny = float(grid['NROWS'])

        x1 = x0 + (xs * nx)
        y1 = y0 + (ys * ny)
        xrange = (x0,x1)
        yrange = (y0,y1)
        grids[i] = {'xrange':xrange,'yrange':yrange, 'nx':nx,'ny':ny,
                   'x0':x0, 'xs':xs, 'x1':x1, 'y0':y0, 'ys':ys, 'y1':y1}
    return(p,grids,proj4_string)
