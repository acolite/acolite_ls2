## def get_rtoa
## reads Sentinel-2 TOA reflectance
## written by Quinten Vanhellemont, RBINS
## 2017-04-18
## modifications: (QV) 2017-12-07 added zoom to target resolution here
##                (QV) 2018-04-17 changed to float32
##                2018-07-18 (QV) changed acolite import name

def get_rtoa(file, metadata, band_data, granule_files, band_name, radiance=False, sub=None, target_res=10):
    from acolite.sentinel import read_band
    from numpy import nan, pi, cos, float32
    from scipy.ndimage import zoom
    
    band_id = [i for i in band_data['BandNames'].keys() if band_data['BandNames'][i] == band_name][0]
    band_res = band_data['Resolution'][band_id]
    band_f0 = band_data['F0'][band_id]

    band_file = granule_files[band_name]['path']

    if sub is not None:
        sub_res = sub['{}'.format(band_res)]['sub']
        data = read_band(band_file, sub=sub_res)
    else:
        data = read_band(band_file, sub=sub)

    ## resample bands to target_res if needed
    if int(band_data['Resolution'][band_id]) != target_res:
        data = zoom(data,int(band_data['Resolution'][band_id])/float(target_res), order=0, mode='nearest')

    nodata = data == int(metadata['NODATA'])
    saturated = data == int(metadata['SATURATED'])

    quantisation = float(metadata['QUANTIFICATION_VALUE'])
    data=data/quantisation

    if radiance is True:
        cossza = cos(metadata['MEAN_SUN_ZENITH']*(pi/180.))
        d = metadata['SE_DISTANCE']
        data *= (pi * d * d) / (band_f0 * cossza)

    if len(nodata) > 0: data[nodata] = nan
    if len(saturated) > 0: data[saturated] = nan

    data = data.astype(float32)
    return data
