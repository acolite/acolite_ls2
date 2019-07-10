## def get_bt
## reads Landsat at sensor brightness temperatures
## written by Quinten Vanhellemont, RBINS
## 2017-12-05
## modifications: 
##                2018-07-18 (QV) changed acolite import name
##                2019-07-10 (QV) added LTOA return

def get_bt(bundle, metadata, band, return_radiance = False, sub=None):
    from acolite.landsat import read_band
    from numpy import log, nan
    
    btag = 'B{}'.format(band)
    if btag not in metadata:
        print('Band {} not found in {}.'.format(btag, bundle))
        return(None)

    file = metadata[btag]
    data = read_band(file, sub=sub)
    nodata = data <= 0

    ## compute LTOA
    offset = metadata['RADIANCE_ADD_BAND_{}'.format(band)]
    slope = metadata['RADIANCE_MULT_BAND_{}'.format(band)]
    data=data*slope
    data+=offset
    if return_radiance:
        data[nodata] = nan
        return(data)

    ## compute BT
    if 'K1_CONSTANT_BAND_{}'.format(band) in metadata.keys():
        K1=float(metadata['K1_CONSTANT_BAND_{}'.format(band)])
    else:
        ## from https://serc.carleton.edu/files/NAGTWorkshops/gis/activities2/student_handout_calculating_te.pdf
        if metadata['SATELLITE'] == 'LANDSAT_5': K1=607.76
        if metadata['SATELLITE'] == 'LANDSAT_7': K1=666.09

    if 'K2_CONSTANT_BAND_{}'.format(band) in metadata.keys():
        K2=float(metadata['K2_CONSTANT_BAND_{}'.format(band)])
    else:
        ## from https://serc.carleton.edu/files/NAGTWorkshops/gis/activities2/student_handout_calculating_te.pdf
        if metadata['SATELLITE'] == 'LANDSAT_5': K2=1282.71
        if metadata['SATELLITE'] == 'LANDSAT_7': K2=1260.56
    data = K2 / log((K1/data)+1.)

    data[nodata] = nan
    return data
