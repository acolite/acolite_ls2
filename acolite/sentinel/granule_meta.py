## def granule_meta
## imports S2 granule metadata
## written by Quinten Vanhellemont, RBINS
## 2017-04-18
## modifications:
##                2018-07-18 (QV) changed acolite import name
##                2020-10-28 (QV) fill nans in angles grids

def granule_meta(metafile):

    import dateutil.parser
    from xml.dom import minidom

    from acolite.shared import distance_se, fillnan
    from acolite.sentinel import safe_tile_grid

    from numpy import linspace, zeros, isfinite, isnan, where

    try:
        xmldoc = minidom.parse(metafile)
    except:
        print('Error opening metadata file.')
        sys.exit()

    xml_main = xmldoc.firstChild

    metadata = {}

    tags = ['TILE_ID','DATASTRIP_ID', 'SENSING_TIME']
    for tag in tags:
        tdom = xmldoc.getElementsByTagName(tag)
        if len(tdom) > 0: metadata[tag] = tdom[0].firstChild.nodeValue

    #Geometric_Info
    grids = {'10':{}, '20':{}, '60':{}}
    Geometric_Info = xmldoc.getElementsByTagName('n1:Geometric_Info')
    if len(Geometric_Info) > 0:
        for tg in Geometric_Info[0].getElementsByTagName('Tile_Geocoding'):
            tags = ['HORIZONTAL_CS_NAME','HORIZONTAL_CS_CODE']
            for tag in tags:
                tdom = tg.getElementsByTagName(tag)
                if len(tdom) > 0: metadata[tag] = tdom[0].firstChild.nodeValue

            for sub in  tg.getElementsByTagName('Size'):
                res = sub.getAttribute('resolution')
                grids[res]['RESOLUTION'] = float(res)
                tags = ['NROWS','NCOLS']
                for tag in tags:
                    tdom = sub.getElementsByTagName(tag)
                    if len(tdom) > 0: grids[res][tag] = int(tdom[0].firstChild.nodeValue)

            for sub in  tg.getElementsByTagName('Geoposition'):
                res = sub.getAttribute('resolution')
                tags = ['ULX','ULY','XDIM','YDIM']
                for tag in tags:
                    tdom = sub.getElementsByTagName(tag)
                    if len(tdom) > 0: grids[res][tag] = int(tdom[0].firstChild.nodeValue)

        for ta in Geometric_Info[0].getElementsByTagName('Tile_Angles'):
            ## sun angles
            sun_angles={}
            for tag in ['Zenith','Azimuth']:
                for sub in ta.getElementsByTagName('Sun_Angles_Grid')[0].getElementsByTagName(tag):
                    sun_angles[tag] = safe_tile_grid(sub)

            for sub in ta.getElementsByTagName('Mean_Sun_Angle'):
                sun_angles['Mean_Zenith'] = float(sub.getElementsByTagName('ZENITH_ANGLE')[0].firstChild.nodeValue)
                sun_angles['Mean_Azimuth'] = float(sub.getElementsByTagName('AZIMUTH_ANGLE')[0].firstChild.nodeValue)

            ## view angles (merge detectors)
            view_angles={}
            for sub in ta.getElementsByTagName('Viewing_Incidence_Angles_Grids'):
                band = sub.getAttribute('bandId')
                detector = sub.getAttribute('detectorId')

                band_view = {}
                for tag in ['Zenith','Azimuth']:
                    band_view[tag] = safe_tile_grid(sub.getElementsByTagName(tag)[0])

                if band not in view_angles.keys():
                    view_angles[band] = band_view
                else:
                    for tag in ['Zenith','Azimuth']:
                        mask = isfinite(band_view[tag]) & isnan(view_angles[band][tag])
                        view_angles[band][tag][mask] = band_view[tag][mask]

            for b,band in enumerate(view_angles.keys()):
                for tag in ['Zenith','Azimuth']:
                    view_angles[band][tag] = fillnan(view_angles[band][tag])

            ## average view angle grid
            ave = {}
            for b,band in enumerate(view_angles.keys()):
                for tag in ['Zenith','Azimuth']:
                    data = view_angles[band][tag]
                    count = isfinite(data)*1
                    if b == 0:
                        ave[tag] = data
                        ave['{}_Count'.format(tag)] = count
                    else:
                        ave[tag] += data
                        ave['{}_Count'.format(tag)] += count
            for tag in ['Zenith','Azimuth']: view_angles['Average_View_{}'.format(tag)] = ave[tag] / ave['{}_Count'.format(tag)]

    metadata["GRIDS"] = grids
    metadata["VIEW"] = view_angles
    metadata["SUN"] = sun_angles

    ## some interpretation
    metadata['TIME'] = dateutil.parser.parse(metadata['SENSING_TIME'])
    metadata["DOY"] = metadata["TIME"].strftime('%j')
    metadata["SE_DISTANCE"] = distance_se(metadata['DOY'])

    return(metadata)
