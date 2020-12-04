## def metadata_parse
## parses landsat metadata to get parameters for processing
## written by Quinten Vanhellemont, RBINS
## 2017-04-13
## modifications: (QV) 2017-06-15 changed band removal
##                2018-07-18 (QV) changed acolite import name
##                2018-10-01 (QV) removed obsolete bits
##                                added ALI
##                2019-02-28 (QV) improved support for LO8 scenes without TIRS
##                2019-03-26 (QV) fixed issue with time parsing
##                2020-12-04 (QV) added support for Coll2 data - whole function could fo with a rewrite

def metadata_parse(bundle):
    import dateutil.parser
    from acolite.shared import distance_se
    from acolite.landsat import metadata_read
    from acolite.landsat import bundle_test
    from acolite.ac import f0_sensor

    bundle_files = bundle_test(bundle)
    metafile = bundle_files['MTL']['path']
    mdata = metadata_read(metafile)


    metadata = {}
    if 'PRODUCT_METADATA' in mdata.keys():
        attkey = 'PRODUCT_METADATA'
    elif 'IMAGE_ATTRIBUTES' in mdata.keys():
        attkey = 'IMAGE_ATTRIBUTES'

    metadata["SATELLITE"] = mdata[attkey]['SPACECRAFT_ID'].strip('"')
    metadata["SENSOR"] = mdata[attkey]['SENSOR_ID'].strip('"')

    if 'IMAGE_ATTRIBUTES' in mdata.keys():
        metadata['NEW_STYLE']=True
        metadata['AZI'] = float(mdata['IMAGE_ATTRIBUTES']['SUN_AZIMUTH'])
        metadata['THS'] = 90. - float(mdata['IMAGE_ATTRIBUTES']['SUN_ELEVATION'])
        metadata['THV'] = 0.
        metadata['ISODATE'] = '{}T{}Z'.format(mdata[attkey]['DATE_ACQUIRED'].strip('"'),mdata[attkey]['SCENE_CENTER_TIME'].strip('"'))

        if 'METADATA_FILE_INFO' in mdata.keys():
            prokey = 'METADATA_FILE_INFO'
            metadata['PROCESSING_DATE'] = mdata[prokey]['FILE_DATE']
            prjkey = 'PRODUCT_METADATA'
            sclkey = 'RADIOMETRIC_RESCALING'
        elif 'LEVEL1_PROCESSING_RECORD' in mdata.keys():
            prokey = 'LEVEL1_PROCESSING_RECORD'
            metadata['PROCESSING_DATE'] = mdata[prokey]['DATE_PRODUCT_GENERATED']
            prjkey = 'PROJECTION_ATTRIBUTES'
            sclkey = 'LEVEL1_RADIOMETRIC_RESCALING'

        metadata['SCENE'] = mdata[prokey]['LANDSAT_SCENE_ID']
        if 'LANDSAT_PRODUCT_ID' in mdata[prokey].keys():
            metadata['PRODUCT'] = mdata[prokey]['LANDSAT_PRODUCT_ID']


        ### get projection info
        metadata["PATH"] = mdata[attkey]['WRS_PATH'].strip('"')
        metadata["ROW"] = mdata[attkey]['WRS_ROW'].strip('"')

        metadata["DIMS"] = [int(mdata[prjkey]['REFLECTIVE_SAMPLES'].strip('"')),
                            int(mdata[prjkey]['REFLECTIVE_LINES'].strip('"'))]

        tags = ['CORNER_LL_PROJECTION_X_PRODUCT', 'CORNER_LL_PROJECTION_Y_PRODUCT',
                'CORNER_LR_PROJECTION_X_PRODUCT', 'CORNER_LR_PROJECTION_Y_PRODUCT',
                'CORNER_UL_PROJECTION_X_PRODUCT', 'CORNER_UL_PROJECTION_Y_PRODUCT',
                'CORNER_UR_PROJECTION_X_PRODUCT', 'CORNER_UR_PROJECTION_Y_PRODUCT']
        for tag in tags: metadata[tag] = float(mdata[prjkey][tag])

        tags = ['CORNER_LL_LAT_PRODUCT', 'CORNER_LL_LON_PRODUCT',
                'CORNER_LR_LAT_PRODUCT', 'CORNER_LR_LON_PRODUCT',
                'CORNER_UL_LAT_PRODUCT', 'CORNER_UL_LON_PRODUCT',
                'CORNER_UR_LAT_PRODUCT', 'CORNER_UR_LON_PRODUCT']
        for tag in tags: metadata[tag] = float(mdata[prjkey][tag])
        ## get rescaling
        bands = set()
        for key in mdata[sclkey].keys():
            metadata[key] = float(mdata[sclkey][key])
            if 'REFLECTANCE' in key: bands.add(key.split('_')[3]) ## add optical bands to list
            if 'RADIANCE' in key: bands.add(key.split('_')[3]) ## add optical bands to list
    else:
        metadata['NEW_STYLE']=False
        metadata['SENSOR'] = mdata['PRODUCT_METADATA']['SENSOR_ID']
        spacecraft = mdata['PRODUCT_METADATA']['SPACECRAFT_ID']
        if spacecraft == 'Landsat5': metadata["SATELLITE"]='LANDSAT_5'
        if spacecraft == 'Landsat7': metadata["SATELLITE"]='LANDSAT_7'
        metadata['AZI'] = float(mdata['PRODUCT_PARAMETERS']['SUN_AZIMUTH'])
        metadata['THS'] = 90. - float(mdata['PRODUCT_PARAMETERS']['SUN_ELEVATION'])
        metadata['THV'] = 0.
        metadata['SCENE'] = mdata['PRODUCT_METADATA']['METADATA_L1_FILE_NAME']

        metadata['PROCESSING_DATE'] = mdata['METADATA_FILE_INFO']['PRODUCT_CREATION_TIME']

        ## ALI image
        if (metadata["SATELLITE"] == 'EO1') & (metadata["SENSOR"]  == 'ALI'):
            metadata['ISODATE'] = '{}T{}Z'.format(mdata['PRODUCT_METADATA']['ACQUISITION_DATE'].replace('-','/'),mdata['PRODUCT_METADATA']['START_TIME'].split()[-1])

            metadata["PATH"] = 'UNK'
            metadata["ROW"] = 'UNK'
            metadata["DIMS"] = [int(mdata['PRODUCT_METADATA']['PRODUCT_SAMPLES_REF'].strip('"')),
                                int(mdata['PRODUCT_METADATA']['PRODUCT_LINES_REF'].strip('"'))]
            ## get rescaling
            ali_bands={'BAND1':'Pan', 'BAND2':'1p', 'BAND3':'1', 'BAND4':'2',
               'BAND5':'3', 'BAND6':'4', 'BAND7':'4p', 'BAND8':'5p',
               'BAND9':'5', 'BAND10':'7'}

            bands = set()
            for key in mdata['RADIANCE_SCALING'].keys():
                metadata[key] = float(mdata['RADIANCE_SCALING'][key])
                band = key.split('_')[0]
                bands.add(ali_bands[band])
            bands=list(bands)
            bands.sort()
        else:
            metadata['ISODATE'] = '{}T{}Z'.format(mdata['PRODUCT_METADATA']['ACQUISITION_DATE'].strip('"'),mdata['PRODUCT_METADATA']['SCENE_CENTER_SCAN_TIME'].strip('"'))

            ### get projection info
            metadata["PATH"] = mdata['PRODUCT_METADATA']['WRS_PATH'].strip('"')
            metadata["ROW"] = mdata['PRODUCT_METADATA']['STARTING_ROW'].strip('"')
            metadata["DIMS"] = [int(mdata['PRODUCT_METADATA']['PRODUCT_SAMPLES_REF'].strip('"')),
                                int(mdata['PRODUCT_METADATA']['PRODUCT_LINES_REF'].strip('"'))]

            ## get rescaling
            bands = set()
            for key in mdata['MIN_MAX_RADIANCE'].keys():
                metadata[key] = float(mdata['MIN_MAX_RADIANCE'][key])
                bands.add(key[-1]) ## add optical bands to list
            for key in mdata['MIN_MAX_PIXEL_VALUE'].keys():
                metadata[key] = float(mdata['MIN_MAX_PIXEL_VALUE'][key])
            print(bands)

        ### get scene extent
        tags = ["PRODUCT_UL_CORNER_LAT","PRODUCT_UL_CORNER_LON",
                "PRODUCT_UR_CORNER_LAT","PRODUCT_UR_CORNER_LON",
                "PRODUCT_LL_CORNER_LAT","PRODUCT_LL_CORNER_LON",
                "PRODUCT_LR_CORNER_LAT","PRODUCT_LR_CORNER_LON",
                "PRODUCT_UL_CORNER_MAPX","PRODUCT_UL_CORNER_MAPY",
                "PRODUCT_UR_CORNER_MAPX","PRODUCT_UR_CORNER_MAPY",
                "PRODUCT_LL_CORNER_MAPX","PRODUCT_LL_CORNER_MAPY",
                "PRODUCT_LR_CORNER_MAPX","PRODUCT_LR_CORNER_MAPY"]
        for tag in tags: metadata[tag] = float(mdata['PRODUCT_METADATA'][tag])

    ### common things
    metadata["TIME"] = dateutil.parser.parse(''.join([metadata["ISODATE"]]).replace("Z", ''))
    metadata["DOY"] = metadata["TIME"].strftime('%j')
    metadata["SE_DISTANCE"] = distance_se(metadata['DOY'])

    if 'PROJECTION_PARAMETERS' in mdata.keys():
        for key in mdata['PROJECTION_PARAMETERS'].keys():
                metadata[key] = mdata['PROJECTION_PARAMETERS'][key]
    elif 'LEVEL1_PROJECTION_PARAMETERS' in mdata.keys():
        for key in mdata['LEVEL1_PROJECTION_PARAMETERS'].keys():
                metadata[key] = mdata['LEVEL1_PROJECTION_PARAMETERS'][key]
    if 'UTM_PARAMETERS' in mdata: metadata['UTM_ZONE']=mdata['UTM_PARAMETERS']['ZONE_NUMBER']

    ## set up bands
    metadata['BANDS']=list(bands)
    metadata['BANDS'].sort()
    metadata['BAND_NAMES'] = metadata['BANDS']

    ## some sensor specific config
    if metadata["SATELLITE"] == 'LANDSAT_8':
        ## get TIRS for L8
        if 'TIRS_THERMAL_CONSTANTS' in mdata:
            for key in mdata['TIRS_THERMAL_CONSTANTS'].keys():
                metadata[key] = float(mdata['TIRS_THERMAL_CONSTANTS'][key])
        ## bands to remove
        remove_bands = ['8','9','10','11'] ## pan cirrus thermal x2
        for rb in remove_bands:
            if rb in metadata['BANDS']: metadata['BANDS'].remove(rb)
        metadata['SATELLITE_SENSOR'] = 'L8_OLI'
        metadata['RGB_BANDS']= [483,561,655]
        metadata['WAVES']= [443,483,561,655,865,1609,2201]
        metadata['BANDS_ALL'] = ['1','2','3','4','5','6','7','9','10','11']
        metadata['BANDS_BESTFIT'] = ['6','7']
        metadata['BANDS_THERMAL'] = ['10','11']
        metadata['BAND_NAMES_ALL'] = ['1','2','3','4','5','6','7','9','10','11']
        metadata['WAVES_ALL']= [443,483,561,655,865,1609,2201,1375,10900,11500]

    if metadata["SATELLITE"] == 'LANDSAT_5':
        remove_bands = ['6'] ## thermal
        for rb in remove_bands:
            if rb in metadata['BANDS']: metadata['BANDS'].remove(rb)
        metadata['SATELLITE_SENSOR'] = 'L5_TM' if metadata['SENSOR'] == 'TM' else 'L5_MSS'
        metadata['RGB_BANDS']= [486,571,660]
        metadata['WAVES']= [486, 571, 660, 839, 1678, 2217]
        metadata['BANDS_ALL'] = ['1','2','3','4','5','6','7']
        metadata['BANDS_BESTFIT'] = ['5','7']
        metadata['BANDS_THERMAL'] = ['6']
        metadata['BAND_NAMES_ALL'] = ['1','2','3','4','5','6','7']
        metadata['WAVES_ALL']= [486, 571, 660, 839, 1678, 10000, 2217]

    if metadata["SATELLITE"] == 'LANDSAT_7':
        remove_bands = ['6','8'] ## thermal pan
        for rb in remove_bands:
            if rb in metadata['BANDS']: metadata['BANDS'].remove(rb)
        metadata['SATELLITE_SENSOR'] = 'L7_ETM'
        metadata['RGB_BANDS']= [479,561,661]
        metadata['WAVES']= [479, 561, 661, 835, 1650, 2208]
        metadata['BANDS_ALL'] = ['1','2','3','4','5','6','7']
        metadata['BANDS_BESTFIT'] = ['5','7']
        metadata['BANDS_THERMAL'] = ['6']
        metadata['BAND_NAMES_ALL'] = ['1','2','3','4','5','6','7']
        metadata['WAVES_ALL']= [479, 561, 661, 835, 1650, 10000, 2208]

    if metadata['SATELLITE'] == 'EO1':
        metadata['SATELLITE_SENSOR'] = 'EO1_ALI'
        metadata['BANDS_BESTFIT'] = ['5','7']
        metadata['BANDS_ALL'] = ['1p', '1', '2', '3', '4', '4p', '5p', '5', '7'] # 'Pan'
        metadata['WAVES_ALL']= [442, 485, 567, 660, 790, 866, 1244, 1640, 2226]

    ## make F0 for bands
    f0 = f0_sensor(metadata['SATELLITE_SENSOR'])
    for band in f0.keys():
        metadata['B{}_F0'.format(band)] = f0[band]

    ## get FILE names for bands and test if they exist in the bundle
    if 'PRODUCT_CONTENTS' in mdata:
        product_key = 'PRODUCT_CONTENTS'
    elif 'PRODUCT_METADATA' in mdata:
        product_key = 'PRODUCT_METADATA'

    for par in mdata[product_key]:
        if 'FILE_NAME' in par:
                fname = mdata[product_key][par]
                exists = False
                for file in bundle_files:
                    if bundle_files[file]['fname'] == fname:
                        exists = True
                        metadata[file] = bundle_files[file]['path']

    return(metadata)
