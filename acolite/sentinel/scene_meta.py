## def scene_meta
## imports S2 scene metadata
## written by Quinten Vanhellemont, RBINS
## 2017-04-18
## modifications: 2017-04-27 (QV) added bandnames fallback for older metadata
##                2017-06-06 (QV) added waves and rgb bands info
##                2017-11-22 (QV) added defaults for the model selection
##                2018-07-18 (QV) changed acolite import name
##                2018-10-01 (QV) removed obsolete bits
##                2020-06-02 (SV) ensure that linspace gets an integer as its third argument

def scene_meta(metafile):
    import dateutil.parser
    from xml.dom import minidom
        
    from acolite.shared import distance_se
    from numpy import linspace

    try: 
        xmldoc = minidom.parse(metafile)
    except: 
        print('Error opening metadata file.')
        sys.exit()
        
    xml_main = xmldoc.firstChild
    
    metadata = {}

    tags = ['PRODUCT_START_TIME','PRODUCT_STOP_TIME','PRODUCT_URI','PROCESSING_LEVEL',
            'PRODUCT_TYPE', 'PROCESSING_BASELINE', 'GENERATION_TIME','SPACECRAFT_NAME',
            'DATATAKE_SENSING_START', 'SENSING_ORBIT_NUMBER', 'SENSING_ORBIT_DIRECTION',
            'PRODUCT_FORMAT', 'QUANTIFICATION_VALUE', 'U']

    for tag in tags:
        tdom = xmldoc.getElementsByTagName(tag)
        if len(tdom) > 0: 
            if tdom[0].firstChild is not None:
                metadata[tag] = tdom[0].firstChild.nodeValue

    ##
    tdom = xmldoc.getElementsByTagName('Special_Values')
    for t in tdom:
        fill = (t.getElementsByTagName('SPECIAL_VALUE_TEXT')[0].firstChild.nodeValue)
        fill_value = (t.getElementsByTagName('SPECIAL_VALUE_INDEX')[0].firstChild.nodeValue)
        metadata[fill] = fill_value

    ## get information for sensor bands
    banddata = {}
    banddata['F0'] = {}
    tdom = xmldoc.getElementsByTagName('SOLAR_IRRADIANCE')
    for t in tdom:
        band = t.getAttribute('bandId')
        banddata['F0'][band] = float(t.firstChild.nodeValue) # 'unit':t.getAttribute('unit')

    banddata['PHYSICAL_GAINS'] = {}
    tdom = xmldoc.getElementsByTagName('PHYSICAL_GAINS')
    for t in tdom:
        band = t.getAttribute('bandId')
        banddata['PHYSICAL_GAINS'][band] = float(t.firstChild.nodeValue)

    banddata['BandNames'] = {}
    banddata['Resolution'] = {}
    banddata['Wavelength'] = {}
    banddata['RSR'] = {}

    tdom = xmldoc.getElementsByTagName('Spectral_Information')
    for t in tdom:
        band = t.getAttribute('bandId')
        banddata['BandNames'][band] = t.getAttribute('physicalBand')
        banddata['Resolution'][band] = t.getElementsByTagName('RESOLUTION')[0].firstChild.nodeValue
        banddata['Wavelength'][band] = {tag:float(t.getElementsByTagName(tag)[0].firstChild.nodeValue) for tag in ['CENTRAL','MIN','MAX']}
        tag = t.getElementsByTagName('Spectral_Response')
        if len(tag) > 0:
            step = float(tag[0].getElementsByTagName('STEP')[0].firstChild.nodeValue)
            rsr = [float(rs) for rs in tag[0].getElementsByTagName('VALUES')[0].firstChild.nodeValue.split(' ')]
            wave = linspace(banddata['Wavelength'][band]['MIN'],banddata['Wavelength'][band]['MAX'], int(ceil(((banddata['Wavelength'][band]['MAX']-banddata['Wavelength'][band]['MIN'])/step)+1)))                                                                                                                                                            
        banddata['RSR'][band] = {'response':rsr, 'wave':wave}
    #print(banddata['Wavelength'])

    if len(banddata['BandNames']) == 0:
        bandnames = ['B1','B2','B3','B4','B5','B6','B7','B8','B8A','B9','B10','B11','B12']
        bandresolutions=[60,10,10,10,20,20,20,10,20,60,60,20,20]
        bandwaves = [{'CENTRAL': 443.9, 'MAX': 457.0, 'MIN': 430.0},
                     {'CENTRAL': 496.6, 'MAX': 538.0, 'MIN': 440.0},
                     {'CENTRAL': 560.0, 'MAX': 582.0, 'MIN': 537.0},
                     {'CENTRAL': 664.5, 'MAX': 684.0, 'MIN': 646.0},
                     {'CENTRAL': 703.9, 'MAX': 713.0, 'MIN': 694.0},
                     {'CENTRAL': 740.2, 'MAX': 749.0, 'MIN': 731.0},
                     {'CENTRAL': 782.5, 'MAX': 797.0, 'MIN': 769.0},
                     {'CENTRAL': 835.1, 'MAX': 908.0, 'MIN': 760.0},
                     {'CENTRAL': 864.8, 'MAX': 881.0, 'MIN': 848.0},
                     {'CENTRAL': 945.0, 'MAX': 958.0, 'MIN': 932.0},
                     {'CENTRAL': 1373.5, 'MAX': 1412.0, 'MIN': 1337.0},
                     {'CENTRAL': 1613.7, 'MAX': 1682.0, 'MIN': 1539.0},
                     {'CENTRAL': 2202.4, 'MAX': 2320.0, 'MIN': 2078.0}]

        for b,band in enumerate(bandnames):
            banddata['BandNames'][str(b)] = band
            banddata['Wavelength'][str(b)] = bandwaves[b]
            banddata['Resolution'][str(b)] = str(bandresolutions[b])

    ## some interpretation
    metadata['TIME'] = dateutil.parser.parse(metadata['PRODUCT_STOP_TIME'])
    metadata["DOY"] = metadata["TIME"].strftime('%j')
    metadata["SE_DISTANCE"] = distance_se(metadata['DOY'])

    metadata["SATELLITE"] = metadata['SPACECRAFT_NAME'] # 'Sentinel-2A'
    metadata["SENSOR"] = 'MSI'

    if metadata['SATELLITE'] == 'Sentinel-2A':
        metadata['SATELLITE_SENSOR'] = 'S2A_MSI'
        metadata['RGB_BANDS']= [497,560,664]
        metadata['WAVES']= ['444', '497', '560', '664', '704', '740', '782', '835', '865',  '945', '1374', '1614', '2202'] #
    
    if metadata['SATELLITE'] == 'Sentinel-2B':
        metadata['SATELLITE_SENSOR'] = 'S2B_MSI'
        metadata['RGB_BANDS']= [497,560,664]
        metadata['WAVES']= ['444', '497', '560', '664', '704', '740', '782', '835', '865',  '945', '1374', '1614', '2202'] #

    bands = list(banddata['BandNames'].keys())
    metadata['BANDS'] = bands
    waves = ['{0:0.0f}'.format(banddata['Wavelength'][i]['CENTRAL']) for i in bands]
    metadata['WAVES'] = [int(w) for w in waves]
    band_names = [banddata['BandNames'][i] for i in bands]
    metadata['BAND_NAMES'] = band_names
    metadata['BANDS_ALL'] = ['1','2','3','4','5','6','7','8','8A','11','12']
    metadata['BAND_NAMES_ALL'] = band_names
    metadata['BANDS_BESTFIT'] = ['11','12']

    return(metadata,banddata)
