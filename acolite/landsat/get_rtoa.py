## def get_rtoa
## reads Landsat TOA reflectance
## written by Quinten Vanhellemont, RBINS
## 2017-04-13
## modifications: (QV) 2017-04-14 added scene center sun zenith angle normalisation
##                (QV) 2018-01-23 added pan reading
##                (QV) 2018-04-17 changed to float32
##                2018-07-18 (QV) changed acolite import name
##                2018-10-01 (QV) added ALI / old format

def get_rtoa(bundle, metadata, band, usgs_reflectance=True, radiance=False, sub=None, pan=None):
    from acolite.landsat import read_band

    from numpy import nan, pi, cos, float32
    
    if (not metadata['NEW_STYLE']):
        la_tag = 'B{}_OFFSET'
        lm_tag = 'B{}_SCALING_FACTOR'
        if metadata['SENSOR'] == 'ALI':
            ali_bands = {'Pan': 'BAND1', '1p': 'BAND2', '1': 'BAND3', '2': 'BAND4', '3': 'BAND5', 
                         '4': 'BAND6', '4p': 'BAND7', '5p': 'BAND8', '5': 'BAND9', '7': 'BAND10'}
            la_tag=la_tag.format(ali_bands[band].strip('B'))
            lm_tag=lm_tag.format(ali_bands[band].strip('B'))
            btmp = ali_bands[band].replace('BAND','').zfill(2)
            file_tag = 'B{}'.format(btmp)
            f0_tag = 'B{}_F0'.format(band.replace('BAND',''))
        else:
            file_tag = 'B{}'.format(band)
            f0_tag = 'B{}_F0'.format(band)
    else:
        ra_tag = 'REFLECTANCE_ADD_BAND_{}'.format(band)
        rm_tag = 'REFLECTANCE_MULT_BAND_{}'.format(band)
        la_tag = 'RADIANCE_ADD_BAND_{}'.format(band)
        lm_tag = 'RADIANCE_MULT_BAND_{}'.format(band)
        file_tag = 'B{}'.format(band)
        f0_tag = 'B{}_F0'.format(band)

    file = metadata[file_tag]
    
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
        offset = metadata[ra_tag]
        slope = metadata[rm_tag]
        data=data*slope
        data+=offset
        ## normalise to sun zenith angle
        data /= cos(metadata['THS']*(pi/180.))
    else:
        #print('Using USGS radiance.')
        offset = metadata[la_tag]
        slope = metadata[lm_tag]

        data=data*slope
        data+=offset

        if radiance is False:
            cossza = cos(metadata['THS']*(pi/180.))

            d = metadata['SE_DISTANCE']
            #print('need to fix F0!')
            f0 = metadata[f0_tag]
            data *= (pi * d * d) / (f0 * cossza)
    
    data[nodata] = nan

    data = data.astype(float32)
    return data
