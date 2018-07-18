## def get_rtoa
## reads Landsat TOA reflectance
## written by Quinten Vanhellemont, RBINS
## 2017-04-13
## modifications: (QV) 2017-04-14 added scene center sun zenith angle normalisation
##                (QV) 2018-01-23 added pan reading
##                (QV) 2018-04-17 changed to float32
##                2018-07-18 (QV) changed acolite import name
def get_rtoa(bundle, metadata, band, usgs_reflectance=True, radiance=False, sub=None, pan=None):
    from acolite.landsat import read_band

    from numpy import nan, pi, cos, float32
    
    file = metadata['B{}'.format(band)]
    
    if (pan is not None) and (sub is not None):
        pansub = [s*2 for s in sub]
        data = read_band(file, sub=pansub)
    else:
        data = read_band(file, sub=sub)

    if (data is None):
        print('Bad crop.')
        return(-1)

    nodata = data <= 0

    tmp = ['REFLECTANCE' in key for key in metadata.keys()]
    if usgs_reflectance:
        if True in tmp:
            usgs_reflectance=True
        else:
            usgs_reflectance=False
       
    if usgs_reflectance:
        ##print('Using USGS reflectance.')
        offset = metadata['REFLECTANCE_ADD_BAND_{}'.format(band)]
        slope = metadata['REFLECTANCE_MULT_BAND_{}'.format(band)]
        data=data*slope
        data+=offset
        ## normalise to sun zenith angle
        data /= cos(metadata['THS']*(pi/180.))
    else:
        #print('Using USGS radiance.')
        offset = metadata['RADIANCE_ADD_BAND_{}'.format(band)]
        slope = metadata['RADIANCE_MULT_BAND_{}'.format(band)]

        data=data*slope
        data+=offset

        if radiance is False:
            cossza = cos(metadata['THS']*(pi/180.))

            d = metadata['SE_DISTANCE']
            #print('need to fix F0!')
            f0 = metadata['B{}_F0'.format(band)]
            data *= (pi * d * d) / (f0 * cossza)
    
    data[nodata] = nan

    data = data.astype(float32)
    return data
