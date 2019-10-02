## def acolite_toa_crop
## crops (and merges) given Landsat or Sentinel tiles for a given limit
##
## note: geometry is taken from the first tile, e.g. in Landsat the scene center sun zenith angle could otherwise give a discontinuity between two tiles
##
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-11-27 (based on acolite_py)
## modifications: 2017-11-28 (QV) updated scene list and tile code naming
##                2017-12-05 (QV) tested old style sentinel 2 scenes, added all bands (cirrus+thermal) for L8
##                                renamed from toa_convert_ncdf
##                                updated nc write function
##                2018-01-11 (QV) Fixed multi scene offset issue
##                2018-01-23 (QV) Added utm zone check
##                2018-03-27 (QV) Added s2 target res
##                2018-04-19 (QV) added pan and pan_ms output options for Landsat 8
##                2018-05-07 (QV) added override keyword
##                2018-06-06 (QV) added support for xy outputs
##                2018-07-18 (QV) changed acolite import name
##                2018-10-24 (QV) fixed s2 grid name formatting
##                2018-11-19 (QV) fixed cirrus band output
##                2019-02-21 (QV) fixed S2 band name issue
##                2019-03-12 (QV) fixed issue with limits extending out of the scene and saving the L8 PAN/MS data
##                2019-03-26 (QV) added CF dataset names
##                2019-04-09 (QV) added pan global dimensions to correspond to 2xMS extended dimensions
##                2019-04-11 (QV) added check for valid data for cropped scenes (blackfill_skip)
##                2019-09-11 (QV) skipping Lt in thermal channels when merging S2 tiles
##                2019-10-02 (QV) added test for MSI L1C files, skip processing for L2A files

def acolite_toa_crop(scenes, odir, limit=None, 
                     ## skip cropped scenes that are in the "blackfill"
                     blackfill_skip=True,
                     blackfill_max=1.0,
                     blackfill_wave = 1600, 
                     nc_compression=True, chunking=True, tile_code=None, s2_target_res=10, 
                     nc_write_geo_xy = False, 
                     l8_output_pan=False, l8_output_pan_ms=False, override = True):
    import acolite as pp
    from numpy import nanmean, nan
    import numpy as np
    from scipy.ndimage import zoom
    import time, os

    new=True

    if limit is None:
        print('No limit given for tile merging/crop, nothing to do.')
        return(1)

    if type(scenes) is not list: scenes = [scenes]

    if tile_code is None:
        if len(scenes) > 1:
            tile_code = 'MERGED' 
        else:
             tile_code = 'CROP'

    for bundle_id, bundle in enumerate(scenes):
        try:
            metadata = pp.landsat.metadata_parse(bundle)
            data_type = "Landsat"
            granules = [bundle]
            if metadata['NEW_STYLE'] is False:
                 print('Old style Landsat not yet configured {}.'.format(bundle))
                 return(1)
        except:
            data_type = None

        if data_type is None:
            try:
                safe_files = pp.sentinel.safe_test(bundle)
                data_type = "Sentinel"
                metafile = safe_files['metadata']['path']
                granules = safe_files['granules']
                metadata,bdata= pp.sentinel.scene_meta(metafile)
                if 'PROCESSING_LEVEL' in metadata:
                    if metadata['PROCESSING_LEVEL'] != 'Level-1C':
                         print('Processing of {} files not supported, skipping {}'.format(metadata['PROCESSING_LEVEL'], bundle))
                         return(1)
            except:
                data_type = None

        if data_type is None:
            print('Inputfile {} not recognised.'.format(bundle))
            return(1)
            
        import os
        if os.path.exists(odir) is False: os.makedirs(odir)

        ## set up some variables
        bands = [b.strip('B') for b in metadata['BAND_NAMES_ALL']]

        ## get RSR wavelengths
        swaves = pp.shared.sensor_wave(metadata['SATELLITE_SENSOR'])
        swavesl = [swaves[b] if b in swaves else nan for b in bands]

        ## run through granules (will be one for Landsat, almost always one for S2 - old style S2 files not tested)
        for granule in granules:
            ptime = time.strftime('%Y-%m-%d %H:%M:%S %Z')
            print('{} - Reading {}'.format(ptime,granule))
            out_of_scene = False
            no_coverage = False

            if data_type == 'Landsat':
                pixel_size = (30,30)
                ## make output names
                metadata['TILE_CODE']=tile_code
                oname = '{}_{}_{}'.format(metadata['SATELLITE_SENSOR'],metadata['TIME'].strftime('%Y_%m_%d_%H_%M_%S'), tile_code)

                ## get scene extent for given limit
                if limit is None:
                    sub = None
                    p, (xrange,yrange), proj4_string = pp.landsat.geo.get_projection(metadata)
                else:
                    scene_extent = pp.landsat.geo.get_sub(metadata, limit)
                    if type(scene_extent) is int: out_of_scene = True
                    else: 
                        sub, p, (xrange,yrange,grid_region), proj4_string = scene_extent
                        ## get "extended" xyranges
                        xrange, yrange = grid_region['xrange'], grid_region['yrange']
                    if limit is None: sub=None
            
                ## calculate view azimuth and update metadata
                view_azi = pp.landsat.view_azimuth(bundle, metadata)
                azi = float(metadata['AZI'])-view_azi
                if azi > 180.: azi-= 180.
                if azi < 0.: azi+= 180.
                metadata['AZI']=azi

                ## get projection
                proj = metadata['MAP_PROJECTION']
                if (proj == 'UTM'):
                    zone = int(metadata['UTM_ZONE'])
                else: zone=None

            if data_type == 'Sentinel':
                pixel_size = (int(s2_target_res),int(s2_target_res))
                ## read granule metadata
                gr_metafile = safe_files[granule]['metadata']['path']
                grmeta = pp.sentinel.granule_meta(gr_metafile)

                ## get crop for given limit
                if limit is not None:
                    grids = pp.sentinel.geo.get_sub(grmeta, limit)
                    if type(grids) is int: out_of_scene = True
                    else: 
                        grids, proj4_string = grids
                        sub = grids['{}'.format(s2_target_res)]['sub']
                        ## get "extended" xyranges
                        xrange = grids['grids_region']['{}'.format(s2_target_res)]['xrange']
                        yrange = grids['grids_region']['{}'.format(s2_target_res)]['yrange']

                else:
                    p, (grids), proj4_string = pp.sentinel.geo.get_projection(grmeta)
                    xrange = grids['{}'.format(s2_target_res)]['xrange']
                    yrange = grids['{}'.format(s2_target_res)]['yrange']
                    grids = None
            
                ## make output names
                tile_id = grmeta['TILE_ID']

                metadata['TILE_CODE']=tile_code
                oname = '{}_{}_{}'.format(metadata['SATELLITE_SENSOR'],metadata['TIME'].strftime('%Y_%m_%d_%H_%M_%S'), tile_code)

                ## update metadata with sun and view geometry
                metadata['THS'] = grmeta['SUN']['Mean_Zenith']
                metadata['THV'] = nanmean(grmeta['VIEW']['Average_View_Zenith'])

                sun_azi = grmeta['SUN']['Mean_Azimuth']
                view_azi = nanmean(grmeta['VIEW']['Average_View_Azimuth'])

                azi = sun_azi - view_azi
                if azi > 180.: azi-= 180.
                if azi < 0.: azi+= 180.
                metadata['AZI']=azi
                metadata['ISODATE']=grmeta['SENSING_TIME']

                ## get zone
                cs_name = grmeta["HORIZONTAL_CS_NAME"]
                zone = cs_name.split(' ')[-1]

            ## exit if crop is out of scene
            if (out_of_scene):
                 print('Region {} out of scene {}'.format(limit,bundle))
                 continue

            ## skip cropped scenes that are in the "blackfill"
            if (sub is not None) & (blackfill_skip):
                bi, wave = pp.shared.closest_idx(swavesl, blackfill_wave)
                if data_type == 'Landsat':
                    band_name = metadata['BANDS_ALL'][bi]
                if data_type == 'Sentinel':
                    band_name = metadata['BAND_NAMES'][bi]
                parname_t = 'rhot_{}'.format(wave)
                if data_type == 'NetCDF':
                    band_data = pp.shared.nc_data(granule, parname_t)
                if data_type == 'Landsat':
                    band_data = pp.landsat.get_rtoa(bundle, metadata, band_name, sub=sub)
                if data_type == 'Sentinel':
                    band_data = pp.sentinel.get_rtoa(bundle, metadata, bdata, safe_files[granule], band_name, target_res=s2_target_res, sub=grids)
                npx = band_data.shape[0] * band_data.shape[1]
                nbf = npx - len(np.where(np.isfinite(band_data))[0])
                band_data = None
                if (nbf/npx) >= float(blackfill_max):
                    print('Skipping scene as crop is {:.0f}% blackfill'.format(100*nbf/npx))
                    continue
            
            ## set up new file and make lat/lon datasets
            if (bundle_id == 0) or (new):
                attributes = {tag:metadata[tag] for tag in metadata.keys()}
                attributes['proj4_string'] = proj4_string
                attributes['xrange'] = xrange
                attributes['yrange'] = yrange
                attributes['pixel_size'] = pixel_size

                attributes["file_type"] =  'Level 1 TOA Reflectance Merged'
                attributes["utm_zone"] = zone

                ## change attribute time to a string for writing to NetCDF
                attributes['TIME'] = attributes['TIME'].strftime('%Y-%m-%d %H:%M:%S')
                
                ## join band names to csv
                attributes['BANDS'] = ','.join(metadata['BANDS'])
                attributes['BAND_NAMES'] = ','.join(metadata['BAND_NAMES'])
                attributes['BANDS_ALL'] = ','.join(metadata['BANDS_ALL'])

                if limit is not None:
                    attributes['LIMIT'] = limit
                    attributes['SUB'] = sub

                ## make outputfiles     
                ncfile = '{}/{}_L1.nc'.format(odir,oname)
                if os.path.exists(ncfile) & (override == False):
                    return(ncfile)

                if l8_output_pan:
                    ncfile_pan = '{}/{}_L1_pan.nc'.format(odir,oname)
                    new_pan = True
                if l8_output_pan_ms:
                    ncfile_pan_ms = '{}/{}_L1_pan_ms.nc'.format(odir,oname)
                    new_pan_ms = True

                if data_type == 'Landsat':
                    lon, lat = pp.landsat.geo.get_ll(metadata, limit=limit, extend_limit=True)
                if data_type == 'Sentinel':
                    lon, lat = pp.sentinel.geo.get_ll(grmeta, limit=limit, extend_limit=True, resolution=s2_target_res)
                crop_dims = lon.shape

                pp.output.nc_write(ncfile, 'lon', lon, new=new, attributes=attributes,
                                    dataset_attributes={'standard_name':'longitude', 'units':'degree_east'},
                                    nc_compression=nc_compression, chunking=chunking)
                new = False
                pp.output.nc_write(ncfile, 'lat', lat, new=new, attributes=attributes, 
                                    dataset_attributes={'standard_name':'latitude', 'units':'degree_north'},
                                    nc_compression=nc_compression, chunking=chunking)

                ##########################
                ## write easting and northing datasets
                if nc_write_geo_xy:
                    print('Writing easting (x) and northing (y) to outputfiles.')
                    if data_type == 'Landsat':
                        x, y = pp.landsat.geo.get_ll(metadata, limit=limit, extend_limit=True, xy=True)
                    if data_type == 'Sentinel':
                        x, y = pp.sentinel.geo.get_ll(grmeta, limit=limit, extend_limit=True, resolution=s2_target_res, xy=True)
                    pp.output.nc_write(ncfile, 'x', x, new=new, attributes=attributes, 
                                        nc_compression=nc_compression, chunking=chunking)
                    pp.output.nc_write(ncfile, 'y', y, new=new, attributes=attributes, 
                                        nc_compression=nc_compression, chunking=chunking)
                ## end write geo_xy
                ##########################

            else:
                print(zone)
                if zone != attributes['utm_zone']:
                    print('UTM zones not compatible.')
                    continue

            ## read rhot and write for each band
            rhos_data={}
            ## moved these out of the band loop 2018-01-11
            offset = None
            replace_nan = False

            for b,band in enumerate(bands):
                if data_type == 'Landsat':
                    ## merge the metadata
                    if b == 0:
                        metadata['THS'] = attributes['THS'] ## to avoid scene edges in merged product
                        replace_nan = True ## tile edges on Landsat are filled with nans
                        offset = grid_region['off']

                    band_name = metadata['BANDS_ALL'][b]

                    if band_name in metadata['BANDS_THERMAL']:
                        band_data = pp.landsat.get_bt(bundle, metadata, band_name, sub=sub)
                        oname = 'bt{}'.format(band_name)
                        wave = None
                        ds_att = {'band_name':band_name}
                    else:
                        band_data = pp.landsat.get_rtoa(bundle, metadata, band_name, sub=sub)
                        wave = swavesl[b] #metadata['WAVES_ALL'][b]
                        oname = 'rhot_{}'.format(wave)
                        ds_att = {'wavelength':float(wave),'band_name':band_name}

                if data_type == 'Sentinel':
                    wave = swavesl[b]
                    oname = 'rhot_{}'.format(wave)
                    band_name = metadata['BAND_NAMES'][b]
                    
                    offset = grids['grids_region']['{}'.format(s2_target_res)]['off']
                    band_data = pp.sentinel.get_rtoa(bundle, metadata, bdata, safe_files[granule], band_name, sub=grids, target_res=s2_target_res)
                    bid = str(b)
                    ds_att = {'wavelength':float(wave), 'band_name':band_name}
                if band_data is None: continue

                ## add CF names
                if 'bt' in oname:
                    ds_att['standard_name']='toa_brightness_temperature'
                    ds_att['units']='K'
                else:
                    ds_att['standard_name']='toa_bidirectional_reflectance'
                    ds_att['units']=1

                ## write to NetCDF file
                pp.output.nc_write(ncfile, oname, band_data, dataset_attributes=ds_att, new=new, offset=offset, replace_nan=replace_nan,
                                           attributes=attributes, nc_compression=nc_compression, chunking=chunking)

                ## also write Lt for thermal bands
                if (data_type == 'Landsat'):
                    if band_name in metadata['BANDS_THERMAL']:
                        band_data = pp.landsat.get_bt(bundle, metadata, band_name, sub=sub, return_radiance=True)
                        oname = 'lt{}'.format(band_name)
                        wave = None
                        ds_att = {'band_name':band_name}
                        ds_att['standard_name']='toa_brightness_temperature'
                        ds_att['units']='W/m2srmum'
                        ## write to NetCDF file
                        pp.output.nc_write(ncfile, oname, band_data, dataset_attributes=ds_att, new=new, offset=offset, replace_nan=replace_nan,
                                           attributes=attributes, nc_compression=nc_compression, chunking=chunking)

            ## read pan band if requested
            if (data_type == 'Landsat') and (l8_output_pan or l8_output_pan_ms):
                if metadata["SATELLITE"] in ['LANDSAT_7','LANDSAT_8']:
                    pan_offset = [o*2 for o in offset]
                    pan_crop_dims = (crop_dims[0]*2, crop_dims[1]*2)

                    if metadata['SATELLITE'] == 'LANDSAT_7': panband = 8
                    if metadata['SATELLITE'] == 'LANDSAT_8': panband = 8
                    band_data = pp.landsat.get_rtoa(bundle, metadata, panband, sub=sub, pan=True)

                    ## write PAN to NetCDF file
                    if l8_output_pan:
                        oname = 'rhot_pan'
                        pp.output.nc_write(ncfile_pan, oname, band_data, new=new_pan, global_dims=pan_crop_dims, offset=pan_offset, replace_nan=replace_nan,
                                                   attributes=attributes, nc_compression=nc_compression, chunking=chunking)
                        new_pan=False

                    ## write PAN at MS resolution
                    if l8_output_pan_ms:
                        oname = 'rhot_pan_ms'
                        band_data = zoom(band_data, zoom=0.5, order=1)
                        pp.output.nc_write(ncfile_pan_ms, oname, band_data, new=new_pan_ms, global_dims=crop_dims, offset=offset, replace_nan=replace_nan,
                                                   attributes=attributes, nc_compression=nc_compression, chunking=chunking)
                        new_pan_ms=False
                    band_data=None

    if 'ncfile' in locals(): return(ncfile)
    else: return(1)
