## def acolite_ac
##
## written by Quinten Vanhellemont, RBINS for the PONDER project
## Apr-June 2017
## modifications: QV 2017-06-07 merged separate Landsat/Sentinel developments
##                QV 2017-06-20 renamed acolite_py - started adding exponential fit options (V&R 2015)
##                QV 2017-07-18 added pressure option to interpolate LUTs
##                QV 2017-10-24 added ancillary data and gas correction
##                QV 2017-11-27 added TOA NetCDF option (e.g. for merged TOA tiles)
##                QV 2017-11-28 updated Landsat "grids_region" and PP data position
##                QV 2017-12-05 updated metadata for landsat and ancillary fallback for old data
##                              renamed to acolite_ac
##                              updated nc_write function
##                QV 2018-01-22 added rsky
##                QV 2018-02-08 added rdark length check
##                QV 2018-03-05 integrated tiled processing
##                QV 2018-03-07 renamed some keyword arguments to include dsf/exp
##                QV 2018-03-08 added l8_output_bt keyword
##                QV 2018-03-12 fixed S2 strftime
##                QV 2018-04-17 fixed S2 azimuth, added exp_wave1, exp_wave2, exp_alpha keywords
##                2018-04-19 (QV) added pan and pan_ms output options for Landsat 8
##                2018-05-07 (QV) added THS_original to output attributes
##                2018-06-06 (QV) added support for xy outputs, added try & except clause for ancillary data
##                2018-06-07 (QV) added check for rhot < rray in DSF
##                2018-07-09 (QV) updated rho_rc 'new' file generation
##                2018-07-18 (QV) changed acolite import name
##                2018-07-24 (QV) added atmospheric correction parameters to metadata (for fixed DSF), renamed t_gas tag to tt_gas, changed flags type to int32
##                                added support for tiled DSF on merged scenes
##                2018-07-25 (QV) added orange band support for tiled DSF
##                                added check for ROI limits just at the edge of the scene
##                2018-07-30 (QV) added glint correction
##                2018-08-02 (QV) added threshold for glint correction (don't do GC over land)
##                2018-09-10 (QV) added glint in tiled mode
##                2018-10-01 (QV) added ALI support
##                2018-10-24 (QV) fixed rhorc output for EXP
##                2018-10-30 (QV) added oxygen transmittance, renamed band list, fixed cirrus band (toa) output
##                2018-11-19 (QV) fixed the glint correction transmittance scaling
##                                fixed tile size for processing s2 at 20 and 60 m
##                                fixed output of exp method when l1r file is present
##                2018-12-04 (QV) fixed DEM support, added elevation input options
##                2018-12-05 (QV) fixed cirrus and wv band output in tiled processing
##                                added resolved_angles for S2 processing
##                2019-03-12 (QV) global attributes xrange and yrange fixed for limits crossing the scene borders
##                2019-03-26 (QV) added some CF dataset names
##                2019-04-11 (QV) added check for valid data for cropped scenes (blackfill_skip)
##                                added check for bright scenes for cropped scenes (cropmask_skip)
##                2019-07-04 (QV) added l1r_nc_delete
##                2019-07-10 (QV) added l8_output_lt_tirs
##                2019-09-18 (QV) made glint correction slightly more RAM friendly for tiled dsf_path_reflectance option
##                2019-10-02 (QV) added test for MSI L1C files, skip processing for L2A files
##                2019-11-29 (QV) added output of extra ac parameters (for fixed DSF only at the moment)
##                2019-12-11 (QV) added output of extra ac parameters for tiled DSF
##                2020-06-22 (QV) added new sky reflectance correction taking the aerosol load into account (for fixed DSF)
##                                added option to skip bands from the DSF (e.g. SWIR bands on S2)
##                2020-06-23 (QV) added new sky reflectance correction for tiled processing (needs some speed improvement)
##                2020-07-14 (QV) added nan fillvalue to reflectance outputs
##                2020-07-15 (QV) commented out r_ray TOA masking, LUT reading changes, added some geometry for NetCDF processing
##                2020-07-29 (QV) massive cleanup started and new fitting procedure done for fixed DSF
##                2020-07-30 (QV) more cleaning and new fitting done for tiled DSF
##                2020-09-15 (QV) added CAMS grib ancillary data support (when the file is provided externally)
##                2020-10-28 (QV) fixed issue with tiled glint correction with resolved angles (S2)
##                2020-11-16 (QV) added garbage collection from Martin Boettcher
##                                reordered glint correction and added valid data slicing
##                2020-11-17 (QV) fixed reordered glint correction for non tiled processing
##                                improved performance of tiled processing for sparse tiles with new tiles_interp and slicing
##                2020-11-19 (QV) various fixes for slicing data and keywording of slicing option
##                2020-12-14 (QV) fix for exp options
##                2021-01-04 (QV) forced double precision for lat/lon writing
##                2021-02-01 (QV) fix for pan subsetting
##                2021-02-02 (QV) fix for masked arrays in sparse NetCDF data
##                2021-02-03 (QV) fix for rhorc outputs in tiled dsf

def acolite_ac(bundle, odir,
                scene_name=False,
                limit=None,
                aerosol_correction='dark_spectrum',

                ## skip cropped scenes that are in the "blackfill"
                blackfill_skip=True,
                blackfill_max=1.0,
                blackfill_wave=1600,

                ## skip cropped scenes that are largely masked
                cropmask_skip=False,
                cropmask_max=0.8,
                cropmask_wave=1600,
                cropmask_threshold=0.1,

                ## old options?
                pixel_idx=200,
                perc_idx=1,
                percentiles = [0,0.1,1,5,10,25,50,75,90,95,99,99.9,100],
                luts=['PONDER-LUT-201704-MOD1', 'PONDER-LUT-201704-MOD2'],#, 'PONDER-LUT-201704-MOD3-1013mb'],
                bestfit='bands',
                bestfit_bands=None,
                pixel_range_min=0,
                pixel_range_max=1000,
                #map_dark_pixels = False,

                ## ACOLITE dark_spectrum settings
                dsf_spectrum_option='dark_list', # 'absolute_pixel'
                dsf_full_scene=False,
                dsf_model_selection='min_drmsd',
                dsf_list_selection='intercept',
                dsf_path_reflectance='fixed', # 'tiled'
                dsf_plot_retrieved_tiles=True,
                dsf_plot_dark_spectrum=True,
                dsf_min_tile_cover = 0.10,
                dsf_min_tile_aot = 0.01,
                dsf_tile_dims = None, ## custom size of tiling grid
                dsf_write_tiled_parameters = False,
                dsf_force_band = None,
                dsf_wave_range = [400,2500],
                dsf_exclude_bands = [],
                extra_ac_parameters=False,

                dsf_fixed_aot=None,
                dsf_fixed_lut='PONDER-LUT-201704-MOD2',

                ## ACOLITE exponential settings
                exp_swir_threshold = 0.0215, ## swir non water mask
                exp_fixed_epsilon=True, ## use fixed (sub)scene epsilon
                exp_fixed_epsilon_percentile=50,
                exp_fixed_aerosol_reflectance=False, ## use fixed (sub)scene aerosol reflectance
                exp_fixed_aerosol_reflectance_percentile = 5,
                exp_wave1=1600,
                exp_wave2=2200,
                exp_alpha=None,
                exp_alpha_weighted=True,
                exp_epsilon=None,
                exp_initial_epsilon = 1.0,
                exp_gamma=None,

                ## Sentinel-2 target resolution
                s2_target_res = 10,

                ## L8/TIRS BT
                l8_output_bt = False,
                l8_output_lt_tirs = False,

                ## pan band use for Landsat 8
                l8_output_pan = False, ## output L1R_pan file
                l8_output_pan_ms = False, ## output L1R_pan_ms file at 30 m
                l8_output_orange = False, ## output L8 orange band
                ##pan_sharpen_rhow = False, ## sharpen VIS/NIR bands using TOA pan relationship - move to other script

                ## pressure
                lut_pressure = True,
                dem_pressure = False,
                dem_pressure_percentile = 25,
                pressure = None,
                elevation = None,

                ## apply gas corrections at TOA
                gas_transmittance = True,
                uoz_default = 0.3,
                uwv_default = 1.5,
                wvlut = '201710C',

                ## resolve the view/illumination angles during processing (for S2 only currently)
                resolved_angles=False,
                resolved_angles_write=False,

                ## use ancillary data for gas transmittances rather than defaults
                ancillary_data = True,
                cams_grib_ancillary = None,

                ## do sky glint correction
                sky_correction = False,
                sky_correction_option = 'all',

                ## glint correction
                glint_correction = False,
                glint_force_band = None,
                glint_mask_rhos_band = 1600,
                glint_mask_rhos_threshold = 0.05,
                glint_write_rhog_all = False,
                glint_write_rhog_ref = False,

                ## gains to be applied at TOA
                gains = False,
                gains_l5_tm=[1,1,1,1,1,1],
                gains_l7_etm=[1,1,1,1,1,1],
                gains_l8_oli=[1,1,1,1,1,1,1],
                gains_s2a_msi=[1,1,1,1,1,1,1,1,1,1,1,1,1],
                gains_s2b_msi=[1,1,1,1,1,1,1,1,1,1,1,1,1],

                ## max wavelength to check for l2_negatives
                neg_wave = 900,

                ## NetCDF outputs
                l1r_nc_compression = False,
                l1r_nc_override = True,
                l1r_nc_delete = False,
                l2r_nc_compression = False,
                nc_write = True,
                nc_write_rhot = True,
                nc_write_rhorc = False,
                nc_write_rhos = True,
                nc_write_rhow = False,
                nc_write_geo = True,
                nc_write_geo_xy = False,
                nc_delete=False,
                nc_compression=True,
                chunking=True,

                ## slicing for sparse mosaicked data
                slicing=True,

                ## debugging
                verbosity=0,
                ret_rdark=False):

    import acolite as pp
    from numpy import nanmean, nanpercentile, isfinite, count_nonzero, zeros, ceil, nan, linspace, min, max, where,float64, float32, int32, abs, minimum, nanmin, nanmax
    import numpy as np

    from netCDF4 import Dataset
    from scipy.ndimage import zoom
    import time, os
    import dateutil.parser
    from scipy.interpolate import interp2d
    import copy

    l2r_files = []
    sub = None

    #####################
    ## find out input file type
    data_type=None
    if '.nc' in bundle:
        try:
            metadata = pp.shared.nc_gatts(bundle)
            data_type = "NetCDF"
        except:
            data_type = None

    if data_type is None:
        try:
            metadata = pp.landsat.metadata_parse(bundle)
            data_type = "Landsat"
        except:
            data_type = None

    if data_type is None:
        safe_files = pp.sentinel.safe_test(bundle)
        try:
            safe_files = pp.sentinel.safe_test(bundle)
            metafile = safe_files['metadata']['path']
            granules = safe_files['granules']
            data_type = "Sentinel"
        except:
            data_type = None

    if data_type is None:
        print('Inputfile {} not recognised.'.format(bundle))
        return(1)
    ## end input file determination
    #####################

    bands_skip_corr = []
    pixel_range = [pixel_range_min, pixel_range_max]
    lowess_frac=0.5

    if data_type == 'NetCDF':
        sensor_family = ''
        metadata['TIME'] = dateutil.parser.parse(metadata['TIME'])

        if 'LANDSAT' in metadata['SATELLITE']:
            sensor_family = 'Landsat'
            bands = metadata['BANDS_ALL'].split(',')
            metadata['BANDS_ALL'] = bands
            band_names = bands
            waves = metadata['WAVES_ALL']
            resolution = 30
        else:
            bands = metadata['BANDS'].split(',')
            metadata['BANDS_ALL'] = bands
            band_names = metadata['BAND_NAMES'].split(',')
            sensor_family = 'Sentinel'
            waves = metadata['WAVES']
            #bands_skip_corr = ['B9','B10']
            resolution = 10

        ## get tiled dimensions
        pixel_size = metadata['pixel_size']
        if dsf_tile_dims is None:
            dsf_tile_dims = (int(6000/pixel_size[0]), int(6000/pixel_size[1]))
        else:
            dsf_tile_dims = int(dsf_tile_dims[0]),int(dsf_tile_dims[1])

        granules = [bundle]
        ## make band dict
        band_dict = {band_name:{'name':band_name, 'wave': waves[b], 'F0':[],
                                'resolution':resolution,
                                'index':b, 'index_name':bands[b], 'lut_name':'{}'.format(band_name.lstrip('B'))}
                                 for b, band_name in enumerate(band_names)}

    if data_type == 'Landsat':
        sensor_family = 'Landsat'
        bands = metadata['BANDS_ALL']
        band_names = metadata['BANDS_ALL']
        waves = metadata['WAVES_ALL']
        granules = [bundle]
        ## make band dict Landsat
        if metadata['SENSOR'] == 'ALI':
            band_dict = {band_name:{'name':band_name, 'wave': waves[b], 'F0':[],
                                    'resolution':30,
                                    'index':b, 'index_name':bands[b], 'lut_name':bands[b]}
                                     for b, band_name in enumerate(band_names)}

        else:
            band_dict = {band_name:{'name':band_name, 'wave': waves[b], 'F0':[],
                                'resolution':30,
                                'index':int(bands[b]), 'index_name':bands[b], 'lut_name':bands[b]}
                                 for b, band_name in enumerate(band_names)}

    if data_type == 'Sentinel':
        sensor_family = 'Sentinel'
        metafile = safe_files['metadata']['path']
        granules = safe_files['granules']
        metadata,bdata= pp.sentinel.scene_meta(metafile)

        if 'PROCESSING_LEVEL' in metadata:
            if metadata['PROCESSING_LEVEL'] != 'Level-1C':
                 print('Processing of {} files not supported, skipping {}'.format(metadata['PROCESSING_LEVEL'], bundle))
                 return(1)

        ftime = metadata['TIME'].strftime('%Y-%m-%dT%H:%M:%SZ')

        bands = list(bdata['BandNames'].keys())
        waves = ['{0:0.0f}'.format(bdata['Wavelength'][i]['CENTRAL']) for i in bands]
        band_names = [bdata['BandNames'][i] for i in bands]

        ## make band dict Sentinel
        band_dict = {band_name:{'name':band_name, 'wave': float(waves[b]), 'F0':bdata['F0'][bands[b]],
                                'resolution':bdata['Resolution'][bands[b]],
                                'index':int(bands[b]), 'index_name':bands[b], 'lut_name':'{}'.format(band_name.lstrip('B'))}
                                 for b, band_name in enumerate(band_names)}

    ## get RSR wavelengths
    swaves = pp.shared.sensor_wave(metadata['SATELLITE_SENSOR'])
    for b in swaves:
        for d in band_dict:
            if band_dict[d]['lut_name'] == b:
                band_dict[d]['wave'] = int(swaves[b])

    ## get band list ordered by wavelength
    ordered_waves = [band_dict[d]['wave'] for d in band_dict]
    ordered_waves.sort()
    ordered_bands = []
    for w in ordered_waves:
        band = [d for d in band_dict if band_dict[d]['wave'] == w]
        ordered_bands.append(band[0])

    ## dimensions for dark spectrum aot subtiles
    ## and Sentinel-2 target_resolution
    bands_skip_thermal = []
    bands_skip_corr = []

    if sensor_family == 'Landsat':
        if metadata['SATELLITE'] == 'LANDSAT_5' :
            bands_skip_corr = ['6']
            bands_skip_thermal = ['6']
            gains_dict = {bn:gains_l5_tm[bi] for bi,bn in enumerate(['1','2','3','4','5','7'])}
        if metadata['SATELLITE'] == 'LANDSAT_7' :
            bands_skip_corr = ['6','8']
            bands_skip_thermal = ['6']
            gains_dict = {bn:gains_l7_etm[bi] for bi,bn in enumerate(['1','2','3','4','5','7'])}
        if metadata['SATELLITE'] == 'LANDSAT_8' :
            ob_cfg = None
            bands_skip_corr = ['8','9','10','11']
            bands_skip_thermal = ['10','11']
            #gains_dict = {bn:gains_l8_oli[bi] for bi,bn in enumerate(ordered_bands) if bn not in bands_skip_corr}
            gains_dict = {bn:gains_l8_oli[bi] for bi,bn in enumerate(['1','2','3','4','5','6','7']) if bn not in bands_skip_corr}
        if metadata['SATELLITE'] == 'EO1' :
            bands_skip_corr = []
            bands_skip_thermal = []
            gains_eo1_tm = [1.0 for bi,bn in enumerate(ordered_bands) if bn not in bands_skip_corr]
            gains_dict = {bn:gains_eo1_tm[bi] for bi,bn in enumerate(ordered_bands) if bn not in bands_skip_corr}

        if dsf_tile_dims is None:
            dsf_tile_dims = [201,201] # 6x6 km
        else:
            dsf_tile_dims = (int(dsf_tile_dims[0]), int(dsf_tile_dims[1]))
    elif sensor_family == 'Sentinel':
        if "S2A" in metadata['SATELLITE_SENSOR']:
            gains_dict = {bn:gains_s2a_msi[bi] for bi,bn in enumerate(ordered_bands)}
        if "S2B" in metadata['SATELLITE_SENSOR']:
            gains_dict = {bn:gains_s2b_msi[bi] for bi,bn in enumerate(ordered_bands)}
        bands_skip_corr = ['B9','B10']
        if s2_target_res == None:
            s2_target_res = 10
        if dsf_tile_dims is None:
            dsf_tile_dims = [int(603 * 10/float(s2_target_res)),int(603 * 10/float(s2_target_res))] # 6x6 km
        else:
            dsf_tile_dims = (int(dsf_tile_dims[0]), int(dsf_tile_dims[1]))
    else:
        print('Sensor family {} not configured.'.format(sensor_family))

    ## make output directory
    if os.path.exists(odir) is False: os.makedirs(odir)

    ## see if a NetCDF file needs to be generated
    nc_write = (nc_write) & (nc_write_rhot | nc_write_rhorc | nc_write_rhos | nc_write_geo)

    ## set up some variables
    sensor = metadata['SATELLITE_SENSOR']
    sat, sen = sensor.split('_')

    ## get rsr
    rsr_file = pp.config['pp_data_dir']+'/RSR/'+sensor+'.txt'
    rsr, rsr_bands = pp.shared.rsr_read(file=rsr_file)

    ## make them float
    waves = [float(w) for w in waves]

    ## run through granules (will be one for Landsat, almost always one for S2)
    for granule in granules:
        l2r_nc_new = True
        ptime = time.strftime('%Y-%m-%d %H:%M:%S')
        print('{} - Processing {}'.format(ptime,granule))

        new=True
        out_of_scene = False
        no_coverage = False

        if data_type == 'NetCDF':
            ## make output names
            if scene_name:
                oname = metadata['SCENE']
            else:
                tile_code = metadata['TILE_CODE']
                oname = '{}_{}_{}'.format(sensor,metadata['TIME'].strftime('%Y_%m_%d_%H_%M_%S'), tile_code)
            ## add scene extent test
            xrange = metadata['xrange']
            yrange = metadata['yrange']
            proj4_string = metadata['proj4_string']
            pixel_size = metadata['pixel_size']

        if data_type == 'Landsat':
            pixel_size = (30,30)

            ## make output names
            if scene_name:
                oname = metadata['SCENE']
            else:
                tile_code = '{}{}'.format(metadata['PATH'].zfill(3),metadata['ROW'].zfill(3))
                oname = '{}_{}_{}'.format(sensor,metadata['TIME'].strftime('%Y_%m_%d_%H_%M_%S'), tile_code)

            ## get scene extent for given limit
            if limit is None:
                sub = None
                p, (xrange,yrange), proj4_string = pp.landsat.geo.get_projection(metadata)
            else:
                scene_extent = pp.landsat.geo.get_sub(metadata, limit)
                if type(scene_extent) is int: out_of_scene = True
                else:
                    sub, p, (xrange,yrange,grid_region), proj4_string = scene_extent
                if limit is None: sub=None

            ## calculate view azimuth and update metadata
            if metadata['SENSOR'] != 'ALI':
                view_azi = pp.landsat.view_azimuth(bundle, metadata)
                azi = float(metadata['AZI'])-view_azi
                if azi > 180.: azi-= 180.
                if azi < 0.: azi+= 180.
                metadata['AZI']=azi

        if data_type == 'Sentinel':
            ## read granule metadata
            gr_metafile = safe_files[granule]['metadata']['path']
            grmeta = pp.sentinel.granule_meta(gr_metafile)
            pixel_size = (int(s2_target_res),int(s2_target_res))
            ## get crop for given limit
            if limit is not None:
                grids = pp.sentinel.geo.get_sub(grmeta, limit)
                if type(grids) is int: out_of_scene = True
                else:
                    grids, proj4_string = grids
                    sub = grids['{}'.format(s2_target_res)]['sub']
                    xrange = grids['{}'.format(s2_target_res)]['xrange']
                    yrange = grids['{}'.format(s2_target_res)]['yrange']
            else:
                p, (grids), proj4_string = pp.sentinel.geo.get_projection(grmeta)
                xrange = grids['{}'.format(s2_target_res)]['xrange']
                yrange = grids['{}'.format(s2_target_res)]['yrange']
                grids = None
                sub = None

            ## make output names
            if scene_name:
                oname = grmeta['TILE_ID']
            else:
                tile_id = grmeta['TILE_ID']
                tile_code = tile_id.split('_')[-2]
                oname = '{}_{}_{}'.format(sensor,metadata['TIME'].strftime('%Y_%m_%d_%H_%M_%S'), tile_code)

            ## update metadata with sun and view geometry
            metadata['THS'] = grmeta['SUN']['Mean_Zenith']
            metadata['THV'] = nanmean(grmeta['VIEW']['Average_View_Zenith'])

            sun_azi = grmeta['SUN']['Mean_Azimuth']
            view_azi = nanmean(grmeta['VIEW']['Average_View_Azimuth'])

            azi = abs(sun_azi - view_azi)
            while(azi >= 180.):
                azi -= 180.

            metadata['AZI']=abs(azi)
            metadata['ISODATE']=grmeta['SENSING_TIME']

        ## exit if crop is too small
        if sub is not None:
            if (sub[2] <= 0) | (sub[3] <= 0):
                print('Insufficient coverange of region {} in scene {}'.format(limit,bundle))
                continue

        ## exit if crop is out of scene
        if (out_of_scene):
            print('Region {} out of scene {}'.format(limit,bundle))
            continue

        ##########################
        ## skip cropped scenes that are in the "blackfill"
        if (sub is not None) & (blackfill_skip):
            bi, bw = pp.shared.closest_idx(ordered_waves, blackfill_wave)
            band_name = ordered_bands[bi]
            wave = band_dict[band_name]['wave']
            parname_t = 'rhot_{:.0f}'.format(wave)
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
        #### end blackfill check
        ##########################

        ##########################
        ## skip cropped scenes that largely masked
        if (sub is not None) & (cropmask_skip):
            bi, bw = pp.shared.closest_idx(ordered_waves, cropmask_wave)
            band_name = ordered_bands[bi]
            wave = band_dict[band_name]['wave']
            parname_t = 'rhot_{:.0f}'.format(wave)
            if data_type == 'NetCDF':
                band_data = pp.shared.nc_data(granule, parname_t)
            if data_type == 'Landsat':
                band_data = pp.landsat.get_rtoa(bundle, metadata, band_name, sub=sub)
            if data_type == 'Sentinel':
                band_data = pp.sentinel.get_rtoa(bundle, metadata, bdata, safe_files[granule], band_name, target_res=s2_target_res, sub=grids)
            npx = band_data.shape[0] * band_data.shape[1]
            ncm = npx - len(np.where(band_data<cropmask_threshold)[0])
            band_data = None
            if (ncm/npx) >= float(cropmask_max):
                print('Skipping scene as crop is likely cloudy ({:.0f}% masked)'.format(100*ncm/npx))
                continue
        #### end mask check
        ##########################

        ## exit if THS is out of scope LUT
        metadata['THS-true'] = metadata['THS']

        ## override is not ideal!
        if (metadata['THS']>79.9):
             print('Warning sun zenith angle {} out of LUT \n overriding to 80 degrees'.format(metadata['THS']))
             metadata['THS']=79.9
             #continue

        ## use DEM to compute pressure at this location
        if (dem_pressure) & (pressure == None) & (lut_pressure):
            if data_type == 'NetCDF':
                lon = pp.shared.nc_data(granule, 'lon')
                lat = pp.shared.nc_data(granule, 'lat')
            if data_type == 'Landsat':
                lon, lat = pp.landsat.geo.get_ll(metadata, limit=limit)
            if data_type == 'Sentinel':
                lon, lat = pp.sentinel.geo.get_ll(grmeta, limit=limit, resolution=s2_target_res)
            #from numpy import nanpercentile
            ## interpolate DEM
            dem = pp.dem.hgt_lonlat(lon, lat)
            ## compute DEM pressure
            demp = pp.ac.pressure_elevation(dem, ratio=False)
            ## use percentile pressure
            pressure = nanpercentile(demp, dem_pressure_percentile)

        ## use given elevation
        if (pressure == None) & (lut_pressure) & (elevation is not None):
            pressure = pp.ac.pressure_elevation(float(elevation), ratio=False)

        ## get NCEP & TOAST ancillary data
        if ancillary_data:
            if ('lat' not in locals()) or ('lat' not in locals()):
                if data_type == 'NetCDF':
                    lon = pp.shared.nc_data(granule, 'lon')
                    lat = pp.shared.nc_data(granule, 'lat')
                if data_type == 'Landsat':
                    lon, lat = pp.landsat.geo.get_ll(metadata, limit=limit)
                if data_type == 'Sentinel':
                    lon, lat = pp.sentinel.geo.get_ll(grmeta, limit=limit, resolution=s2_target_res)

            ## use image/region mid-point
            pc_lon=lon[int(lon.shape[0]/2), int(lon.shape[1]/2)]
            pc_lat=lat[int(lat.shape[0]/2), int(lat.shape[1]/2)]
            pc_date = metadata['TIME'].strftime('%Y-%m-%d')
            pc_time=metadata['TIME'].hour + metadata['TIME'].minute/60. + metadata['TIME'].second/3600.

            ## use default NASA ancillary data
            if cams_grib_ancillary is None:
                try:
                    pc_anc = pp.ac.ancillary.ancillary_get(pc_date, pc_lon, pc_lat, ftime=pc_time, kind='nearest', verbosity=verbosity)
                except:
                    pc_anc = {}
                    print('Could not retrieve ancillary data, proceeding with default values.')
            else:
                print('Using CAMS grib data: {}'.format(cams_grib_ancillary))
                ## read cams grib ancillary
                pc_anc = pp.ac.ancillary.cams_read_grib(cams_grib_ancillary, pc_date, pc_lon, pc_lat, ftime=pc_time,
                                   spatial_selection = 'interpolated', verbosity = verbosity)

            ## get pressure from ancillary data if not determined by user or by DEM
            if (pressure == None) & (lut_pressure):
                if 'press' not in pc_anc:
                    print('No ancillary pressure found: using default.')
                    pressure=None
                else: pressure = pc_anc['press']['interp']

        ## get gas transmittances
        if gas_transmittance:
            uoz=uoz_default
            uwv=uwv_default

            ## use ancillary data if provided
            if ancillary_data:
                if 'ozone' in pc_anc: uoz=pc_anc['ozone']['interp']
                else: print('No ancillary ozone found: using default {}.'.format(uoz))
                if 'p_water' in pc_anc: uwv=pc_anc['p_water']['interp']/10.
                else:print('No ancillary ozone found: using default {}.'.format(uwv))
            ## compute transmittances
            tt_oz = pp.ac.o3_transmittance(sensor, metadata, uoz=uoz)
            tt_wv = pp.ac.wvlut_interp(metadata['THS'], metadata['THV'], uwv=uwv, sensor=sensor, config=wvlut)
            tt_o2 = pp.ac.o2lut_interp(metadata['THS'], metadata['THV'], sensor=sensor)
            tt_gas = {btag: tt_oz[btag] * tt_wv[btag] * tt_o2[btag] for btag in tt_oz.keys()}

        if pressure is None: pressure = 1013.
        pressure = float(pressure)

        ## Sky reflectance correction
        if (sky_correction) & (sky_correction_option=='all'):
            rsky = pp.ac.toa_rsky(metadata, pressure=pressure, lut = luts[0])

        ## make outputfiles
        metadata["OUTPUT_BASE"] =  '{}/{}'.format(odir,oname)
        l1r_ncfile = '{}/{}_L1R.nc'.format(odir,oname)
        l1r_ncfile_pan = '{}/{}_L1R_pan.nc'.format(odir,oname)
        l1r_ncfile_pan_ms = '{}/{}_L1R_pan_ms.nc'.format(odir,oname)
        l2r_ncfile = '{}/{}_L2R.nc'.format(odir,oname)

        ds_plot = '{}/{}_dark_spectrum.png'.format(odir,oname)
        dp_map = '{}/{}_dark_pixels.png'.format(odir,oname)

        ####################################
        ## set up attributes
        attributes = {'sensor':sensor, 'isodate':metadata['ISODATE'],
                              'THS':metadata['THS'],'THV':metadata['THV'], 'AZI':metadata['AZI'],
                              'pressure':pressure if pressure is not None else 1013.25,
                              'aerosol_correction':aerosol_correction}
        ## track projection
        attributes['proj4_string'] = proj4_string
        attributes['xrange'] = xrange
        attributes['yrange'] = yrange
        attributes['pixel_size'] = pixel_size

        ## add the used x and yranges to the attributes
        ## i.e. the ones inside the image
        ## crop dimensions = (x,y) (sub[2],sub[3])
        if limit is not None:
            #attributes['xrange_crop'] = (xrange[0], xrange[0]+(sub[2])*pixel_size[0])
            #attributes['yrange_crop'] = (yrange[1], yrange[1]-(sub[3])*pixel_size[1])
            attributes['xrange'] = (xrange[0], xrange[0]+(sub[2])*pixel_size[0])
            attributes['yrange'] = (yrange[1]-(sub[3])*pixel_size[1], yrange[1])

        if 'THS-true' in metadata: attributes['THS-true'] = metadata['THS-true']

        ## output attributes
        attributes["output_base"] =  metadata['OUTPUT_BASE']
        attributes["output_dir"] =  odir
        attributes["output_name"] =  oname
        attributes['l1_file'] = bundle
        attributes["l1r_file"] =  l1r_ncfile
        attributes["l2r_file"] =  l2r_ncfile
        attributes["file_type"] =  'Level 2 Reflectance Product'

        attributes['ancillary_data'] = str(ancillary_data)
        attributes['sky_correction'] = str(sky_correction)

        ## add gas transmittance information
        if gas_transmittance is True:
            attributes["uoz"] = uoz
            attributes["uwv"] = uwv
            for band in tt_gas.keys():
                attributes['{}_tt_gas'.format(band)] = tt_gas[band]

        ## add ancillary data information
        if pressure is not None: attributes["pressure"] = pressure
        if ancillary_data is True:
            for k in pc_anc.keys():
                if type(pc_anc[k]) is dict:
                    if 'series' in pc_anc[k].keys(): attributes['anc_{}_series'.format(k)] = pc_anc[k]['series']
                    if 'interp' in pc_anc[k].keys(): attributes['anc_{}'.format(k)] = pc_anc[k]['interp']
                else:
                    attributes['anc_{}'.format(k)] = pc_anc[k]

        ## add sky correction information
        if (sky_correction) & (sky_correction_option == 'all'):
            for band in rsky.keys():
                attributes['{}_r_sky'.format(band)] = rsky[band]
        ## end common attributes
        ####################################

        ## new PONDER atmospheric correction 'dark_spectrum'
        if aerosol_correction == 'dark_spectrum':
            attributes['dsf_spectrum_option']=dsf_spectrum_option
            attributes['dsf_full_scene']=1 if dsf_full_scene else 0
            attributes['dsf_path_reflectance']=dsf_path_reflectance

            ## set up l1r_ncfile
            ## read from existing l1r_ncfile
            l1r_read_nc = True
            l1r_read_nc = (l1r_read_nc) and (os.path.exists(l1r_ncfile)) and (limit is None)
            l1r_nc_write = False if l1r_read_nc else True
            if l1r_nc_override: l1r_nc_write=True
            ##
            l1r_datasets = []
            if os.path.exists(l1r_ncfile): l1r_datasets = pp.shared.nc_datasets(l1r_ncfile)
            l1r_nc_new = True

            l1r_read_nc = False

            ## read lut data
            t0 = time.time()
            print('Reading LUT data')
            lutd = pp.aerlut.import_luts(base_luts=luts, sensor=sensor, add_rsky = (sky_correction) & (sky_correction_option == 'rsky_new'))
            t1 = time.time()
            print('Reading LUTs took {:.1f} seconds'.format(t1-t0))

            ## LUT dataset to fit to the dark spectrum
            if (sky_correction) & (sky_correction_option=='rsky_new'):
                fit_tag = 'romix+rsky_t' ## combined romix rsk
            else:
                fit_tag = 'romix' ## just romix

            ## estimate Rayleigh for TOA masking  from new LUTS
            rray = {band:float(lutd[luts[0]]['rgi'][band]((pressure, lutd[luts[0]]['ipd'][fit_tag], \
                                                           attributes['AZI'], attributes['THV'], attributes['THS'], 0.001))) for band in rsr_bands}

            ## scene dimensions are X,Y
            if data_type == 'NetCDF':
                nc = Dataset(bundle)
                global_dims=(nc.dimensions['y'].size, nc.dimensions['x'].size)
                nc.close()
            else:
                if sensor_family == 'Landsat':
                    dims = metadata['DIMS']
                if sensor_family == 'Sentinel':
                    dims = grmeta['GRIDS']['{}'.format(s2_target_res)]['NCOLS'],grmeta['GRIDS']['{}'.format(s2_target_res)]['NROWS']

                ## python dims are Y,X
                global_dims = (dims[1], dims[0])

            ## python dims are Y,X
            tile_dims = (int(dsf_tile_dims[1]),int(dsf_tile_dims[0]))

            ## number of tiles
            tiles = [int(ceil(global_dims[0]/tile_dims[0])),int(ceil(global_dims[1]/tile_dims[1]))]
            ntiles = tiles[0]*tiles[1]

            if (dsf_path_reflectance == 'tiled') & ((tiles[0] < 2) | (tiles[1] < 2)):
                print('Scene too small to perform tiled DSF. Using fixed subscene DSF.')
                dsf_path_reflectance = 'fixed'

            ## set up resolved angles (S2 only at the moment)
            if resolved_angles:
                ## required ac parameters - use this approach also for fixed and tiled processing?
                ac_pars = ['romix','rorayl','dtotr','utotr',
                           'dtott','utott','astot', 'ttot']
                if extra_ac_parameters: ac_pars = lutd[luts[0]]['meta']['par']
                if (sky_correction) & (sky_correction_option == 'rsky_new'): ac_pars += ['rsky_t']

                if (data_type == 'Sentinel'):
                    tp_dim = grmeta['VIEW']['Average_View_Zenith'].shape
                    ## setup interpolation grid here
                    tp_xnew = linspace(0, tp_dim[1]-1, global_dims[1])
                    tp_ynew = linspace(0, tp_dim[0]-1, global_dims[0])
                    ## for tiled processing use the subtile centre interpolated values
                    ## perhaps it is better to match the tiled processing with the S2 grid?
                    if dsf_path_reflectance == 'tiled':
                        ## tie point coordinates
                        txg=linspace(0, tp_dim[1]-1,num=tp_dim[1])
                        tyg=linspace(0, tp_dim[0]-1,num=tp_dim[0])
                        ## interpolators for tiled processing
                        sunzi = interp2d(txg, tyg, grmeta['SUN']['Zenith'])
                        sunai = interp2d(txg, tyg, grmeta['SUN']['Azimuth'])
                        senzi = interp2d(txg, tyg, grmeta['VIEW']['Average_View_Zenith'])
                        senai = interp2d(txg, tyg, grmeta['VIEW']['Average_View_Azimuth'])

                elif (data_type == 'NetCDF'):
                    input_bundle_datasets = pp.shared.nc_datasets(bundle)
                    resolved_geometry_type = None
                    ## do average geometry
                    if ('raa' in input_bundle_datasets) & ('vza' in input_bundle_datasets) & ('sza' in input_bundle_datasets):
                        resolved_geometry_type = 'average'
                    ## or per band geometry
                    else:
                        tmp = True
                        for b in band_dict.keys():
                            if ('vaa_{}'.format(b) not in input_bundle_datasets): tmp = False
                            if ('vza_{}'.format(b) not in input_bundle_datasets): tmp = False
                        if tmp: resolved_geometry_type = 'per_band'
                    if resolved_geometry_type is None:
                        print('Could not find required geometry dataset.')
                        resolved_angles = False

                    else:
                        sza_ave = pp.shared.nc_data(bundle, 'sza')
                        if resolved_geometry_type == 'average':
                            raa_ave = pp.shared.nc_data(bundle, 'raa')
                            vza_ave = pp.shared.nc_data(bundle, 'vza')
                ## end NetCDF resolved angles
                else:
                    ## to add if e.g. L8 inputs from angles are given
                    print('Resolved angles processing not supported for {}'.format(data_type))
                    resolved_angles=False
            ## end resolved angles

            #########################
            ## do tiled DSF AOT retrieval
            if dsf_path_reflectance == 'tiled':
                print('Starting DSF aerosol correction with tiled path reflectance.')
                ## read data and save tiled rdark
                ## can also save L1T netcdf (recommended)

                ## this opens each band once and then does the tiling
                ## is about 10x faster than processing per tile (at least for the 6x6 km tiles)

                t0 = time.time()
                print('Reading TOA data and performing tiling')

                ## empty tile tracker
                tile_data = {band_name:{'rdark':np.zeros(tiles, dtype=np.float32),
                                        'cover':np.zeros(tiles, dtype=np.float32),
                                        'tau550':np.zeros(tiles, dtype=np.float32)} for band_name in ordered_bands}
                ## tracker for geometry
                tile_angles = {v: np.zeros(tiles, dtype=np.float32)+nan for v in ['VAA', 'SAA', 'VZA', 'SZA']}

                ## run through bands: (1) read, (2) write to NCDF, (3) get dark value
                for b,band_name in enumerate(ordered_bands):
                    ## read band data
                    if band_name in bands_skip_thermal: continue

                    ## set up band parameter
                    wave = band_dict[band_name]['wave']
                    parname = 'rhot_{:.0f}'.format(wave)
                    ds_att = {'wavelength':float(wave),'band_name':band_name}
                    ds_att['standard_name']='toa_bidirectional_reflectance'
                    ds_att['units']=1

                    ## read data and make full tile NC file
                    print('Reading band {}'.format(band_name))
                    if (l1r_read_nc) & (parname in l1r_datasets):
                        print('Reading band {} from L1R NetCDF: {}'.format(band_name,l1r_ncfile))
                        band_full = pp.shared.nc_data(l1r_ncfile, parname)
                    else:
                        print('Reading band {} from input bundle: {}'.format(band_name,bundle))
                        if data_type == 'NetCDF':
                            band_full = pp.shared.nc_data(bundle, parname)
                        if data_type == 'Landsat':
                            band_full = pp.landsat.get_rtoa(bundle, metadata, band_name)
                        if data_type == 'Sentinel':
                            band_full = pp.sentinel.get_rtoa(bundle, metadata, bdata, safe_files[granule], band_name, target_res=s2_target_res)

                        ## write to L1R file
                        if l1r_nc_write:
                            print('Writing band {} to {}'.format(band_name,l1r_ncfile))
                            if limit is not None:
                                band_sub = band_full[sub[1]:sub[1]+sub[3],sub[0]:sub[0]+sub[2]]
                                crop_dims = (sub[3],sub[2])
                                pp.output.nc_write(l1r_ncfile, parname, band_sub, dataset_attributes=ds_att, fillvalue=np.nan,
                                               new=l1r_nc_new, global_dims=crop_dims, nc_compression=l1r_nc_compression)
                            else:
                                pp.output.nc_write(l1r_ncfile, parname, band_full, dataset_attributes=ds_att, fillvalue=np.nan,
                                                   new=l1r_nc_new, global_dims=global_dims, nc_compression=l1r_nc_compression)
                            l1r_nc_new=False

                    ## skip here so the cirrus and wv bands are output to L1R
                    if band_name in bands_skip_corr: continue

                    if (gains) & (band_name in gains_dict):
                        band_full *= gains_dict[band_name]
                        print('Applied gain {} for band {}'.format(gains_dict[band_name],band_name))

                    ## mask toa < rray
                    band_full[band_full < rray[band_dict[band_name]['lut_name']]] = np.nan

                    ##
                    ## run through tiles
                    tileid = 0
                    for xi in range(tiles[1]):
                        x0 = xi*tile_dims[1]

                        xsub = [x0, min((x0+tile_dims[1], global_dims[1]))-x0]

                        #for yi in range(tiles[1]):
                        for yi in range(tiles[0]):

                            tileid += 1
                            print('Subsetting tile {} of {}'.format(tileid, ntiles), end='\r' if tileid < ntiles else '\n')

                            y0 = yi*tile_dims[0]
                            ysub = [y0, min((y0+tile_dims[0], global_dims[0]))-y0]

                            ## compute tile average geometry
                            if resolved_angles:
                                if (data_type == 'Sentinel'):
                                    ## get subtile average view and sun geometry
                                    tsub = [ysub[0], xsub[0], ysub[0]+ysub[1],xsub[0]+xsub[1]]
                                    tmidx = int((tsub[1]+tsub[3])/2)
                                    tmidy = int((tsub[0]+tsub[2])/2)
                                    tile_angles['VAA'][yi, xi] = senai(tp_xnew[tmidx], tp_ynew[tmidy])[0]
                                    tile_angles['SAA'][yi, xi] = sunai(tp_xnew[tmidx], tp_ynew[tmidy])[0]
                                    tile_angles['VZA'][yi, xi] = senzi(tp_xnew[tmidx], tp_ynew[tmidy])[0]
                                    tile_angles['SZA'][yi, xi] = sunzi(tp_xnew[tmidx], tp_ynew[tmidy])[0]
                                if (data_type == 'NetCDF'):
                                    tile_angles['VAA'][yi, xi] = np.nanmedian(raa_ave[ysub[0]:ysub[0]+ysub[1], xsub[0]:xsub[0]+xsub[1]])
                                    tile_angles['SAA'][yi, xi] = 0
                                    tile_angles['VZA'][yi, xi] = np.nanmedian(vza_ave[ysub[0]:ysub[0]+ysub[1], xsub[0]:xsub[0]+xsub[1]])
                                    tile_angles['SZA'][yi, xi] = np.nanmedian(sza_ave[ysub[0]:ysub[0]+ysub[1], xsub[0]:xsub[0]+xsub[1]])
                            ## end compute average geometry

                            ## offsets for writing data to NetCDF
                            offset=[xsub[0],ysub[0]]
                            band_data = band_full[ysub[0]:ysub[0]+ysub[1], xsub[0]:xsub[0]+xsub[1]]
                            ## fix for aot issue with masked sparse data
                            if type(band_data) == np.ma.core.MaskedArray:
                                band_data[band_data.mask] = np.nan
                                band_data = band_data.data

                            ## get tile coverage
                            elements = band_data.shape[0]*band_data.shape[1]
                            if elements == 0: continue
                            finite = count_nonzero(isfinite(band_data[:,:]))
                            cover = finite/elements
                            if elements < 1000:
                                mincover = 1
                            else:
                                mincover = 1000./elements

                            ## flipped yi, xi here to have proper tile orientation
                            tile_data[band_name]['cover'][yi, xi] = cover
                            if cover < mincover: continue

                            ## get rdark
                            rdark_band=pp.ac.select_rdark(band_data, rdark_list_selection=dsf_list_selection,
                                               pixel_range=pixel_range, lowess_frac = lowess_frac, pixel_idx=pixel_idx)

                            ## store rdark in tile data
                            ## flipped yi, xi here to have proper tile orientation
                            tile_data[band_name]['rdark'][yi, xi] = rdark_band

                            ## exclude bands from dsf
                            ## bit of a dirty solution
                            if (band_dict[band_name]['wave'] < dsf_wave_range[0]): tile_data[band_name]['rdark'][yi, xi] += 1
                            if (band_dict[band_name]['wave'] > dsf_wave_range[1]): tile_data[band_name]['rdark'][yi, xi] += 1
                            if (band_dict[band_name]['lut_name'] in dsf_exclude_bands): tile_data[band_name]['rdark'][yi, xi] += 1
                            ##

                t1 = time.time()
                print('Reading {} tiles took {:.1f} seconds'.format(ntiles, t1-t0))

                ## from now on read the l1r NCDF
                l1r_read_nc = (l1r_nc_write) and (os.path.exists(l1r_ncfile))

                ## compute tiled AOT map
                ## flipped yi,xi here, but it does not matter so much here
                t0 = time.time()
                print('Calculating tiled AOT map')
                ## set up empty tile_output
                tags = ['tau550', 'band', 'model']
                tile_output = {tag: np.zeros(tiles, dtype=np.float32)+np.nan for tag in tags}
                tags = ['ratm', 'rorayl','dtott', 'utott', 'astot', 'ttot']
                if (sky_correction) & (sky_correction_option == 'rsky_new'): tags.append('rsky')

                ## add extra ac parameters output
                if (extra_ac_parameters) & (dsf_write_tiled_parameters):
                    for p in lutd[luts[0]]['meta']['par']:
                        if p not in tags: tags.append(p)
                tile_output['atm'] = {band:{tag: np.zeros(tiles, dtype=np.float32)+np.nan for tag in tags} for band in ordered_bands}
                #if (sensor == 'L8_OLI') & (l8_output_pan | l8_output_pan_ms  | l8_output_orange):
                #    print(tile_output['atm'].keys())
                #    print(band_dict.keys())
                #    tile_output['atm']['8'] = {tag: np.zeros(tiles, dtype=np.float32)+np.nan for tag in tags}

                ## run through tiles
                ntiles = tiles[0]*tiles[1]
                tileid = 0

                for xi in range(tiles[1]):
                    for yi in range(tiles[0]):
                        tileid += 1
                        print('Processing tile {} of {}'.format(tileid, ntiles), end='\r' if tileid < ntiles else '\n')

                        ## skip sparse tiles
                        if tile_data[ordered_bands[0]]['cover'][yi,xi] < dsf_min_tile_cover: continue

                        ## make rdark for each tile
                        rdark = {band_dict[band_name]['lut_name']:tile_data[band_name]['rdark'][yi,xi]
                                     for band_name in ordered_bands if band_name not in bands_skip_corr}

                        ## make a copy of the metadata
                        metadata_tile = copy.copy(metadata)

                        ## compute tile average geometry
                        if resolved_angles:
                            if (data_type == 'Sentinel'):
                                ## add to metadata for select model (to improve)
                                metadata_tile['THS']=tile_angles['SZA'][yi, xi]
                                metadata_tile['THV']=tile_angles['VZA'][yi, xi]
                                azii = abs(tile_angles['SAA'][yi, xi]-tile_angles['VAA'][yi, xi])
                                while azii>180: azii=abs(azii-180)
                                metadata_tile['AZI']=azii
                            if (data_type == 'NetCDF'):
                                metadata_tile['THS']=tile_angles['SZA'][yi, xi]
                                metadata_tile['THV']=tile_angles['VZA'][yi, xi]
                                azii = abs(tile_angles['SAA'][yi, xi]-tile_angles['VAA'][yi, xi])
                                while azii>180: azii=abs(azii-180)
                                metadata_tile['AZI']=azii

                        ## updated rhod fitting
                        ## Jul 2020
                        rdark_sel = pp.ac.rhod_select(rdark, rdark_list_selection=dsf_list_selection)
                        raa = 1.0*metadata_tile['AZI']
                        vza = 1.0*metadata_tile['THV']
                        sza = 1.0*metadata_tile['THS']
                        fit_res = pp.ac.rhod_fit_model(rdark_sel, raa, vza, sza, fit_tag=fit_tag, pressure = pressure, lutd=lutd)
                        sel_mod, tau550, res = fit_res
                        sel_mod_number = int(sel_mod[-1])
                        dark_name = res[sel_mod]['band_sel']
                        ret_par = {}
                        for ap in ['romix','rorayl','dtotr','utotr','dtott', 'utott', 'astot', 'ttot']:
                            ret_par[ap] = {band:float(lutd[sel_mod]['rgi'][band]((pressure, lutd[sel_mod]['ipd'][ap], raa, vza, sza, tau550))) for band in rsr_bands}
                        ## end updated fitting

                        ## output extra ac parameters
                        if (extra_ac_parameters) & (dsf_write_tiled_parameters):
                            extra_ac = {}
                            for ap in lutd[luts[0]]['meta']['par']:
                                if ap in ['ratm', 'rorayl','dtott', 'utott', 'astot', 'ttot']: continue
                                extra_ac[ap] = {band:float(lutd[sel_mod]['rgi'][band]((pressure, lutd[sel_mod]['ipd'][ap],
                                                                                       raa, vza, sza, tau550))) for band in rsr_bands}

                        ## store retrievals per band
                        for b,band_name in enumerate(ordered_bands):
                            if band_name in bands_skip_corr: continue
                            ## tau computed based on rdark in this band
                            tile_data[band_name]['tau550'][yi,xi]=tau550#tau550_all_bands[band_dict[band_name]['lut_name']]

                            ## atmospheric parameters computed with lowest tau550 for this dark spectrum
                            for ap in ['romix', 'rorayl', 'dtott', 'utott', 'astot', 'ttot']:
                                apo = ap
                                if ap == 'romix': apo = 'ratm'
                                tile_output['atm'][band_name][apo][yi,xi] = ret_par[ap][band_dict[band_name]['lut_name']]

                            ## add ray cor option
                            if nc_write_rhorc:
                                for ap in ['dtotr', 'utotr']:
                                    if ap not in tile_output['atm'][band_name]:
                                        tile_output['atm'][band_name][ap] = np.zeros(tiles, dtype=np.float32)+np.nan
                                    tile_output['atm'][band_name][ap][yi,xi] = ret_par[ap][band_dict[band_name]['lut_name']]

                            ## add new sky correction
                            if fit_tag == 'romix+rsky_t':
                                tile_output['atm'][band_name]['rsky'][yi,xi] = \
                                    float(lutd[sel_mod]['rgi'][band_dict[band_name]['lut_name']]((pressure, lutd[sel_mod]['ipd']['rsky_t'], raa, vza, sza, tau550)))
                                #float(lutd[sel_mod]['rgi'][band]((pressure, lutd[sel_mod]['ipd']['rsky_t'], raa, vza, sza, tau550))) for band in rsr_bands}

                            #if (sky_correction) & (sky_correction_option == 'rsky_new'):
                            #    tile_output['atm'][band_name]['rsky'][yi,xi] = rsky_s[band_dict[band_name]['lut_name']]

                            ## if output extra ac parameters
                            if (extra_ac_parameters) & (dsf_write_tiled_parameters):
                                for ap in extra_ac:
                                    tile_output['atm'][band_name][ap][yi,xi] = extra_ac[ap][band_dict[band_name]['lut_name']]

                            ## glint in tiled mode
                            #if glint_correction:
                            #    #ttot_s = {band:float(lutd[sel_mod]['rgi'][band]((pressure, lutd[sel_mod]['ipd']['ttot'],
                            #    #                                   raa, vza, sza, tau550))) for band in rsr_bands}
                            #    #print(pressure, raa, vza, sza, tau550)
                            #    #ttot_s = {band:float(lutd[sel_mod]['rgi'][band_dict[band_name]['lut_name']]((pressure, lutd[sel_mod]['ipd']['ttot'], raa, vza, sza, tau550)))}
                            #    if 'ttot' not in tile_output['atm'][band_name]: tile_output['atm'][band_name]['ttot'] = zeros(tiles)+nan
                            #    print(band_name, tile_output['atm'][band_name]['ttot'][yi,xi])
                            #    # tile_output['atm'][band_name]['ttot'][yi,xi] = ttot_s[band_dict[band_name]['lut_name']]
                            ## end tiled glint

                        ## orange band in tiled mode
                        if (l8_output_orange) & (sensor_family == 'Landsat') & (metadata['SATELLITE'] == 'LANDSAT_8'):
                            if ob_cfg is None:
                                ## get orange band config
                                ob_cfg_file = pp.config['pp_data_dir']+'/Shared/oli_orange.cfg'
                                ob_cfg = pp.shared.import_config(ob_cfg_file)

                            if ob_cfg['combine'] == 'after':
                                ## save pan band results in current tile
                                band_name = '8'
                                if band_name not in tile_output['atm']: tile_output['atm'][band_name] = {tag: np.zeros(tiles, dtype=np.float32)+np.nan for tag in tags}
                                for ap in ['romix', 'rorayl', 'dtott', 'utott', 'astot']:
                                    apo = ap
                                    if ap == 'romix': apo = 'ratm'
                                    tile_output['atm'][band_name][apo][yi,xi] = ret_par[ap][band_name]

                                ## add new sky correction
                                if (sky_correction) & (sky_correction_option == 'rsky_new'):
                                    #tile_output['atm'][band_name]['rsky'][yi,xi] = rsky_s[band_name]
                                    if fit_tag == 'romix+rsky_t':
                                        tile_output['atm'][band_name]['rsky'][yi,xi] = \
                                            float(lutd[sel_mod]['rgi'][band_name]((pressure, lutd[sel_mod]['ipd']['rsky_t'], raa, vza, sza, tau550)))

                                ## if output extra ac parameters
                                if (extra_ac_parameters) & (dsf_write_tiled_parameters):
                                    for ap in extra_ac:
                                        tile_output['atm'][band_name][ap][yi,xi] = extra_ac[ap][band_name]
                            else:
                                ob_sensor = 'L8_OLI_ORANGE'
                                ac_model = sel_model_lut_meta['base']
                                if type(ac_model) == list: ac_model=ac_model[0]
                                ob_lut = ac_model.split('-')[0:4] + ['1013mb']
                                ob_lut = '-'.join(ob_lut)
                                ## get orange band LUT
                                ob_lutdir = pp.config['pp_data_dir']+'/LUT'
                                ob_rsr_file = pp.config['pp_data_dir']+'/RSR/'+ob_sensor+'.txt'
                                ob_lut_sensor, ob_meta_sensor = pp.aerlut.aerlut_pressure(ob_lut, ob_lutdir, pressure, ob_sensor, ob_rsr_file)
                                ## get atmospheric correction parameters
                                ratm_o,rorayl_o,dtotr_o,utotr_o,dtott_o,utott_o,astot_o=\
                                pp.aerlut.lut_get_ac_parameters_fixed_tau_sensor(ob_lut_sensor,ob_meta_sensor,\
                                                                                 cur_azi, cur_thv, cur_ths,tau550)
                                ## save orange band results in current tile
                                band_name = 'O'
                                if band_name not in tile_output['atm']: tile_output['atm'][band_name] = {tag: np.zeros(tiles, dtype=np.float32)+np.nan for tag in tags}

                                tile_output['atm'][band_name]['ratm'][yi,xi] = ratm_o[band_name]
                                tile_output['atm'][band_name]['rorayl'][yi,xi] = rorayl_o[band_name]
                                tile_output['atm'][band_name]['dtott'][yi,xi] = dtott_o[band_name]
                                tile_output['atm'][band_name]['utott'][yi,xi] = utott_o[band_name]
                                tile_output['atm'][band_name]['astot'][yi,xi] = astot_o[band_name]

                                ## add new sky correction
                                if (sky_correction) & (sky_correction_option == 'rsky_new'):
                                    rsr_o, rsr_o_bands = pp.shared.rsr_read(file=ob_rsr_file)
                                    rsky_o = pp.shared.rsr_convolute_dict(waves, ret_rsky, rsr_o)
                                    tile_output['atm'][band_name]['rsky'][yi,xi] = rsky_o[band_name]
                        ## end tiled ob

                        ## selected parameters per tile
                        tile_output['tau550'][yi,xi] = tau550
                        tile_output['model'][yi,xi] = sel_mod_number

                        #if len(dark_idx) == 1:
                        #     dark_name = dark_idx
                        #elif len(dark_idx) == 2:
                        #     dark_name = dark_idx[0]

                        if sensor_family == 'Sentinel': dark_name = 'B{}'.format(dark_name)
                        tile_output['band'][yi,xi] = band_dict[dark_name]['index']

                t1 = time.time()
                print('Computing AOT for {} tiles took {:.1f} seconds'.format(ntiles, t1-t0))

                #################
                ## set up tile grid for interpolation of path reflectance parameters
                ## output coordinates
                xnew = linspace(0, tiles[1]-1, global_dims[1])
                ynew = linspace(0, tiles[0]-1, global_dims[0])
                if limit is not None:
                    xnew=xnew[sub[0]:sub[0]+sub[2]]
                    ynew=ynew[sub[1]:sub[1]+sub[3]]
                    global_dims = (sub[3],sub[2])

                ## mask and add buffer to tiles
                tile_mask = tile_output['tau550']<dsf_min_tile_aot

                if dsf_plot_retrieved_tiles:
                    pp.plotting.plot_tiled_retrieval(tiles, tile_data, tile_output, metadata, rdark, odir)
            ## end tiled retrieval
            ##################################

            ##################################
            ## fixed path reflectance for the scene
            ## get a single path reflectance for the (sub) scene
            if dsf_path_reflectance == 'fixed':
                print('Starting DSF aerosol correction with fixed path reflectance.')

                ## read toa reflectance for dark spectrum
                rdark={}
                for b,band_name in enumerate(ordered_bands):
                    ## read band data
                    if band_name in bands_skip_thermal: continue
                    #if band_name in bands_skip_corr: continue

                    ## set up band parameter
                    wave = band_dict[band_name]['wave']
                    parname_t = 'rhot_{:.0f}'.format(wave)
                    ds_att = {'wavelength':float(wave),'band_name':band_name}
                    ds_att['standard_name']='toa_bidirectional_reflectance'
                    ds_att['units']=1

                    ## read data
                    if (l1r_read_nc) & (parname_t in l1r_datasets):
                        band_data = pp.shared.nc_data(l1r_ncfile, parname_t)
                    else:
                        if data_type == 'NetCDF':
                            band_data = pp.shared.nc_data(granule, parname_t)
                        if data_type == 'Landsat':
                            band_data = pp.landsat.get_rtoa(bundle, metadata, band_name, sub=sub)
                        if data_type == 'Sentinel':
                            band_data = pp.sentinel.get_rtoa(bundle, metadata, bdata, safe_files[granule], band_name, target_res=s2_target_res, sub=grids)

                        ## write to L1R file
                        if l1r_nc_write:
                            print('Writing band {} to {}'.format(band_name,l1r_ncfile))
                            if limit is not None:
                                crop_dims = (sub[3],sub[2])
                                pp.output.nc_write(l1r_ncfile, parname_t, band_data, dataset_attributes=ds_att, fillvalue=np.nan,
                                                   new=l1r_nc_new, global_dims=crop_dims, nc_compression=l1r_nc_compression)

                            else:
                                pp.output.nc_write(l1r_ncfile, parname_t, band_data, dataset_attributes=ds_att, fillvalue=np.nan,
                                                       new=l1r_nc_new, global_dims=global_dims, nc_compression=l1r_nc_compression)
                            l1r_nc_new=False

                    ## skip
                    if band_name in bands_skip_corr: continue

                    ## mask toa < rray
                    band_data[band_data < rray[band_dict[band_name]['lut_name']]] = np.nan

                    ## apply gains
                    if (gains) & (band_name in gains_dict):
                        band_data *= gains_dict[band_name]
                        print('Applied gain {} for band {}'.format(gains_dict[band_name],band_name))

                    ## apply gas correction
                    if gas_transmittance:
                        band_data/=tt_gas[band_dict[band_name]['lut_name']]
                        ds_att['tt_gas'] = tt_gas[band_dict[band_name]['lut_name']]

                    ## apply sky correction
                    if sky_correction:
                        if sky_correction_option == 'all':
                            band_data -= rsky[band_dict[band_name]['lut_name']]
                            ds_att['rsky'] = rsky[band_dict[band_name]['lut_name']]

                    data={}
                    data[band_dict[band_name]['lut_name']] = band_data
                    band_data = None

                    ## get dark spectrum
                    rtoa_dict_cur, perc_all_cur = pp.ac.get_dark_spectrum(data, option=dsf_spectrum_option,
                                                                              percentiles=percentiles, perc_idx=perc_idx,
                                                                              pixel_idx=pixel_idx)
                    data = None
                    rdark[band_dict[band_name]['lut_name']] = rtoa_dict_cur[band_dict[band_name]['lut_name']]

                    ## exclude bands from dsf
                    ## bit of a dirty solution
                    if (band_dict[band_name]['wave'] < dsf_wave_range[0]): rdark[band_dict[band_name]['lut_name']] += 1
                    if (band_dict[band_name]['wave'] > dsf_wave_range[1]): rdark[band_dict[band_name]['lut_name']] += 1
                    if (band_dict[band_name]['lut_name'] in dsf_exclude_bands): rdark[band_dict[band_name]['lut_name']] += 1
                    ##

                    ## check for empty rdark
                    if (dsf_spectrum_option=='absolute_pixel_list'):
                        if len(rdark[band_dict[band_name]['lut_name']]) == 0: no_coverage = True

                if no_coverage:
                    print('No valid pixels found for region {} in scene {}: region probably outside swath'.format(limit, granule))
                    continue

                ## check length of rdark
                if dsf_spectrum_option=='dark_list':
                    rdark_len = [len(rdark[b]) for b in rdark]
                    max_len = max(rdark_len)
                    min_len = min(rdark_len)
                    rdark = {b:rdark[b][0:min_len] for b in rdark}

                ## new fitting 20200729
                if True:
                    raa = attributes['AZI']
                    vza = attributes['THV']
                    sza = attributes['THS']
                    ## do rdark selection
                    rdark_sel = pp.ac.rhod_select(rdark, rdark_list_selection=dsf_list_selection)

                    if dsf_fixed_aot is None:
                        fit_res = pp.ac.rhod_fit_model(rdark_sel, raa, vza, sza, fit_tag=fit_tag, pressure = pressure, lutd=lutd)
                        sel_mod, tau550, res = fit_res
                    else:
                        sel_mod, tau550 = dsf_fixed_lut.strip(), float(dsf_fixed_aot)
                        res =  {sel_mod : {'fit_rmsd': np.nan, 'band_sel':'None', 'band_sel2':'None'}}

                    attributes['ac_model'] = sel_mod
                    attributes['ac_aot550'] = tau550
                    attributes['ac_rmsd'] = res[sel_mod]['fit_rmsd']
                    attributes['ac_band'] = res[sel_mod]['band_sel']+'/'+res[sel_mod]['band_sel2']

                    ## sel aermod
                    if 'MOD1' in sel_mod:
                        sel_aermod = 1
                        model_char = 'C'
                    if 'MOD2' in sel_mod:
                        sel_aermod = 2
                        model_char = 'M'
                    if 'MOD3' in sel_mod:
                        sel_aermod = 3
                        model_char = 'M'
                    attributes['ac_model_char'] = model_char
                    sel_mod_number = sel_aermod

                    print('Fitted model:{} band:{} aot={:.3f}'.format(attributes['ac_model'],attributes['ac_band'], attributes['ac_aot550']))

                    ## lut things
                    waves = lutd[sel_mod]['meta']['wave']

                    ### get resampled parameters
                    if True:
                    ## from the sensor lut
                        ratm_s = {band:float(lutd[sel_mod]['rgi'][band]((pressure, lutd[sel_mod]['ipd']['romix'], raa, vza, sza, tau550))) for band in rsr_bands}
                        rorayl_s = {band:float(lutd[sel_mod]['rgi'][band]((pressure, lutd[sel_mod]['ipd']['rorayl'], raa, vza, sza, tau550))) for band in rsr_bands}
                        dtotr_s = {band:float(lutd[sel_mod]['rgi'][band]((pressure, lutd[sel_mod]['ipd']['dtotr'], raa, vza, sza, tau550))) for band in rsr_bands}
                        utotr_s = {band:float(lutd[sel_mod]['rgi'][band]((pressure, lutd[sel_mod]['ipd']['utotr'], raa, vza, sza, tau550))) for band in rsr_bands}
                        dtott_s = {band:float(lutd[sel_mod]['rgi'][band]((pressure, lutd[sel_mod]['ipd']['dtott'], raa, vza, sza, tau550))) for band in rsr_bands}
                        utott_s = {band:float(lutd[sel_mod]['rgi'][band]((pressure, lutd[sel_mod]['ipd']['utott'], raa, vza, sza, tau550))) for band in rsr_bands}
                        astot_s = {band:float(lutd[sel_mod]['rgi'][band]((pressure, lutd[sel_mod]['ipd']['astot'], raa, vza, sza, tau550))) for band in rsr_bands}

                        ## add rsky results
                        if fit_tag == 'romix+rsky_t':
                            ## OSOAA LUT rsky
                            rsky_boa = {band:float(lutd[sel_mod]['rgi'][band]((pressure, lutd[sel_mod]['ipd']['rsky_s'], raa, vza, sza, tau550))) for band in rsr_bands}
                            rsky_toa = {band:float(lutd[sel_mod]['rgi'][band]((pressure, lutd[sel_mod]['ipd']['rsky_t'], raa, vza, sza, tau550))) for band in rsr_bands}
                            for band in rsky_toa.keys():
                                attributes['{}_r_sky'.format(band)] = rsky_toa[band]
                                attributes['{}_r_sky_boa'.format(band)] = rsky_boa[band]
                        else:
                            rsky_boa = None
                            rsky_toa = None
                            if sky_correction:
                                ## old Rayleigh only
                                rsky_boa = rsky
                                rsky_toa = rsky

                ## resolved angles for S2
                if resolved_angles:
                    if (data_type == 'Sentinel'):
                        ## run through Sentinel-2 angle grids
                        ## here the a/c parameters are computed at the S2 grid for the retrieved tau
                        ## tau is computed above from rdark for scene average geometry - which is hopefully acceptable
                        ac_data = {}
                        relazi = abs(grmeta['SUN']['Azimuth']-grmeta['VIEW']['Average_View_Azimuth'])
                        relazi[relazi>180] = abs(relazi[relazi>180]-180)
                        for par in ac_pars:
                            for ob in ordered_bands:
                                b = band_dict[ob]['lut_name']
                                if b not in ac_data: ac_data[b] = {}
                                ac_data[b][par] = lutd[sel_mod]['rgi'][b]((pressure, lutd[sel_mod]['ipd'][par],
                                                            relazi, grmeta['VIEW']['Average_View_Zenith'],
                                                            grmeta['SUN']['Zenith'], attributes['ac_aot550']))
                        ## final interpolation grid for given ROI
                        if limit is not None:
                            tp_xnew=tp_xnew[sub[0]:sub[0]+sub[2]]
                            tp_ynew=tp_ynew[sub[1]:sub[1]+sub[3]]
                            tp_sub_dims = (sub[3],sub[2])
                    ## end S2 resolved angles
                ## end resolved angles
                ### fixed retrieval done

                if ret_rdark:
                    return(metadata, rdark_sel, band_dict)

                ## set up attributes
                if (dsf_full_scene is False) and (limit != None): ds_origin = 'sub scene'
                else: ds_origin = 'full scene'
                attributes['dsf_origin']=ds_origin
                #attributes['dsf_spectrum_option']=dsf_spectrum_option
                attributes['dsf_pixel_idx']=pixel_idx
                attributes['dsf_percentile']=percentiles[perc_idx]
                attributes['dsf_bestfit']=bestfit
                attributes['dsf_model_selection']=dsf_model_selection

                #Output additional parameters from the ac
                if extra_ac_parameters:
                    extra_ac = {}
                    for ap in lutd[sel_mod]['meta']['par']:
                        extra_ac[ap] = {band:float(lutd[sel_mod]['rgi'][band]((pressure, lutd[sel_mod]['ipd'][ap],
                                                                               raa, vza, sza, attributes['ac_aot550']))) for band in rsr_bands}
                ## end extra ac parameters

                ## store the band attributes
                for band in ratm_s.keys():
                     attributes['{}_ratm'.format(band)] = ratm_s[band]
                     attributes['{}_rorayl'.format(band)] = rorayl_s[band]
                     attributes['{}_dtotr'.format(band)] = dtotr_s[band]
                     attributes['{}_utotr'.format(band)] = utotr_s[band]
                     attributes['{}_dtott'.format(band)] = dtott_s[band]
                     attributes['{}_utott'.format(band)] = utott_s[band]
                     attributes['{}_astot'.format(band)] = astot_s[band]
                     if extra_ac_parameters:
                         for ap in extra_ac:
                             atag = '{}_{}'.format(band,ap)
                             if atag not in attributes: attributes[atag] = extra_ac[ap][band]
                ### end band attributes

                if dsf_plot_dark_spectrum:
                    pp.plotting.plot_dark_spectrum(metadata, ds_plot, bands, ordered_bands, data_type, ordered_waves, ratm_s, rorayl_s, rdark, rdark_sel, dsf_spectrum_option, attributes['ac_band'], tau550,sel_mod_number, rsky_boa)

                ## from now on read the l1r NCDF
                l1r_read_nc = (l1r_nc_write) and (os.path.exists(l1r_ncfile))
        ## end dark spectrum fitting a/c
        ###########################

        ##########################
        ## old ACOLITE exponential a/c
        if aerosol_correction == 'exponential':
            l1r_read_nc = False

            ## get Rayleigh reflectance from only 1 LUT
            lutd = pp.aerlut.import_luts(base_luts=[luts[0]], sensor=sensor)
            raa = attributes['AZI']
            vza = attributes['THV']
            sza = attributes['THS']
            rorayl = {band:float(lutd[luts[0]]['rgi'][band]((pressure, lutd[luts[0]]['ipd']['rorayl'], raa, vza, sza, 0.001))) for band in rsr_bands}
            dtotr = {band:float(lutd[luts[0]]['rgi'][band]((pressure, lutd[luts[0]]['ipd']['dtotr'], raa, vza, sza, 0.001))) for band in rsr_bands}
            utotr = {band:float(lutd[luts[0]]['rgi'][band]((pressure, lutd[luts[0]]['ipd']['utotr'], raa, vza, sza, 0.001))) for band in rsr_bands}

            if metadata['SENSOR'] == 'MSI':
                #band_indices = [int(b) for i,b in enumerate(band_names)]
                band_indices = [int(i) for i,b in enumerate(ordered_bands)]
            else:
                #tags = list(bdict.keys())
                tags = []
                for b,band in enumerate(ordered_bands):
                    band_name = ordered_bands[b]
                    if band in bands_skip_corr: continue
                    if band in bands_skip_thermal: continue
                    btag = '{}'.format(band_name.lstrip('B'))
                    tags.append(btag)
                band_indices = [i for i,b in enumerate(ordered_bands) if b.lstrip('B') in tags]

            ## find SWIR bands
            ## find requested bands
            #exp_wave1 = 865
            #exp_wave2 = 1600
            swir1_idx, swir1_wv = pp.shared.closest_idx(ordered_waves, 1650.)
            short_idx, short_wv = pp.shared.closest_idx(ordered_waves, float(exp_wave1))
            long_idx, long_wv = pp.shared.closest_idx(ordered_waves, float(exp_wave2))

            short_wave = '{:.0f}'.format(short_wv)
            long_wave = '{:.0f}'.format(long_wv)
            short_tag = '{}'.format(ordered_bands[short_idx].lstrip('B'))
            long_tag = '{}'.format(ordered_bands[long_idx].lstrip('B'))

            print('Selected bands {}/{} ({}/{} nm)'.format(short_tag, long_tag, short_wave,long_wave))

            if data_type == 'NetCDF':
                short_data = pp.shared.nc_data(bundle, 'rhot_{}'.format(short_wave)).data
                long_data = pp.shared.nc_data(bundle, 'rhot_{}'.format(long_wave)).data
                if (short_idx != swir1_idx) & (long_idx != swir1_idx):
                    mask_data = pp.shared.nc_data(bundle, 'rhot_{:.0f}'.format(waves[swir1_idx])).data

            if data_type == 'Landsat':
                short_data = pp.landsat.get_rtoa(bundle, metadata, ordered_bands[short_idx], sub=sub)
                long_data = pp.landsat.get_rtoa(bundle, metadata, ordered_bands[long_idx], sub=sub)
                if (short_idx != swir1_idx) & (long_idx != swir1_idx):
                    mask_data = pp.landsat.get_rtoa(bundle, metadata, ordered_bands[swir1_idx], sub=sub)

            if data_type == 'Sentinel':
                short_data = pp.sentinel.get_rtoa(bundle, metadata, bdata, safe_files[granule], \
                                                  ordered_bands[short_idx], target_res=s2_target_res, sub=grids)
                long_data = pp.sentinel.get_rtoa(bundle, metadata, bdata, safe_files[granule], \
                                                  ordered_bands[long_idx], target_res=s2_target_res, sub=grids)
                if (short_idx != swir1_idx) & (long_idx != swir1_idx):
                    mask_data = pp.sentinel.get_rtoa(bundle, metadata, bdata, safe_files[granule], \
                                                      ordered_bands[swir1_idx], target_res=s2_target_res, sub=grids)
            ## apply gains
            if (gains) & (ordered_bands[short_idx] in gains_dict) & (ordered_bands[long_idx] in gains_dict):
                short_data *= gains_dict[ordered_bands[short_idx]]
                long_data *= gains_dict[ordered_bands[long_idx]]
                print('Applied gain {} for band {}'.format(gains_dict[ordered_bands[short_idx]],band_names[short_idx]))
                print('Applied gain {} for band {}'.format(gains_dict[ordered_bands[long_idx]],band_names[long_idx]))

            short_data -= rorayl[short_tag]
            long_data -= rorayl[long_tag]

            ## choose processing option
            if (short_wv < 900) & (long_wv < 900):
                exp_option = 'red/NIR'
            elif (short_wv < 900) & (long_wv > 1500):
                exp_option = 'NIR/SWIR'
            else:
                exp_option = 'SWIR'

            print('Using {} for exponential correction'.format(exp_option))

            ## compute mask
            if short_idx == swir1_idx:
                mask = short_data >= exp_swir_threshold
            elif long_idx == swir1_idx:
                mask = long_data >= exp_swir_threshold
            else:
                mask = mask_data >= exp_swir_threshold
                mask_data = None

            exp_cwlim = 0.005

            ## compute aerosol epsilon band ratio
            epsilon = short_data/long_data
            epsilon[np.where(mask)] = nan

            ## red/NIR option
            if exp_option == 'red/NIR':
                print('Using similarity spectrum for red/NIR')
                exp_fixed_epsilon = True

                ## transmittances in both bands
                tt_short = (dtotr[long_tag] * utotr[long_tag] * tt_gas[short_tag])
                tt_long = (dtotr[long_tag] * utotr[long_tag] * tt_gas[long_tag])

                ## compute gamma
                exp_gamma = tt_short/tt_long if exp_gamma == None else float(exp_gamma)
                print('Gamma: {:.2f}'.format(exp_gamma))

                ## compute alpha
                if exp_alpha == None:
                    if exp_alpha_weighted:
                        exp_alpha = pp.shared.similarity_ratio_sensor(sensor, short_tag, long_tag)
                    else:
                       exp_alpha = pp.shared.similarity_ratio_wave(float(short_wave), float(long_wave))
                    print('Using default alpha: {:0.2f}'.format(exp_alpha))
                else:
                    print('Using user supplied alpha: {:0.2f}'.format(exp_alpha))

                ## first estimate of rhow to find clear waters
                exp_c1 = (exp_alpha/tt_long)/(exp_alpha*exp_gamma-exp_initial_epsilon)
                exp_c2 = exp_initial_epsilon * exp_c1
                rhow_short = (exp_c1 * short_data) - (exp_c2 * long_data)

                ## additional masking for epsilon
                mask2 = (rhow_short < 0.) &\
                        (rhow_short > exp_cwlim)
                epsilon[mask2] = nan
                mask2 = None
                rhow_short = None
            elif exp_option == 'NIR/SWIR':
                exp_fixed_epsilon = True
                ## additional masking for epsilon
                mask2 = (long_data < ((short_data+0.005) * 1.5) ) &\
                        (long_data > ((short_data-0.005) * 0.8) ) &\
                        ((long_data + 0.005)/short_data > 0.8)
                epsilon[mask2] = nan
                mask2 = None
            elif exp_option == 'SWIR':
                tmp = None

            ## get fixed epsilon if requested
            if exp_fixed_aerosol_reflectance: exp_fixed_epsilon = True
            if exp_fixed_epsilon:
                if exp_epsilon is not None:
                    epsilon = float(exp_epsilon) * 1.0
                    print('User provided epsilon ({}/{} nm): {:.2f}'.format(short_wave,long_wave, epsilon))
                else:
                    epsilon = nanpercentile(epsilon,exp_fixed_epsilon_percentile)
                    print('{:.0f}th percentile epsilon ({}/{} nm): {:.2f}'.format(exp_fixed_epsilon_percentile, short_wave,long_wave, epsilon))
                epsall = pp.ac.exponential_epsilon(epsilon,waves,idx1=short_idx,idx2=long_idx)

            ## determination of rhoam in long wavelength
            if exp_option == 'red/NIR':
                rhoam = (exp_alpha * exp_gamma * long_data - short_data) / (exp_alpha * exp_gamma - epsilon)
            else:
                rhoam = long_data*1.0

            short_data = None
            long_data = None

            if exp_fixed_aerosol_reflectance:
                rhoam = nanpercentile(rhoam,exp_fixed_aerosol_reflectance_percentile)
                #if not exp_fixed_epsilon: exp_fixed_epsilon=True
                print('{:.0f}th percentile rhoam ({} nm): {:.5f}'.format(exp_fixed_aerosol_reflectance_percentile, long_wave, rhoam))
        ##### end old acolite a/c
        #########################

        #########################
        ## output L2R reflectances
        if nc_write:
            l2_negatives = None
            print('Computing surface reflectances.')
            t0 = time.time()

            ##########################
            ## read toa reflectance for writing files rhos
            for b,band_name in enumerate(ordered_bands):
                ## read band data
                if band_name in bands_skip_thermal: continue
                ## set up band parameter
                btag = band_dict[band_name]['lut_name']
                wave = band_dict[band_name]['wave']

                parname_t = 'rhot_{:.0f}'.format(wave)
                parname_s = 'rhos_{:.0f}'.format(wave)
                parname_w = 'rhos_{:.0f}'.format(wave)
                ds_att = {'wavelength':float(wave),'band_name':band_name}

                ## read data
                if (l1r_read_nc) & (os.path.exists(l1r_ncfile)):
                    datasets = pp.shared.nc_datasets(l1r_ncfile)
                    if parname_t not in datasets: continue
                    print('Reading band {} from L1R NetCDF: {}'.format(band_name,l1r_ncfile))
                    band_data = pp.shared.nc_data(l1r_ncfile, parname_t)
                else:
                    print('Reading band {} from input bundle: {}'.format(band_name,bundle))
                    if data_type == 'NetCDF':
                        band_data = pp.shared.nc_data(granule, parname_t)
                    if data_type == 'Landsat':
                        band_data = pp.landsat.get_rtoa(bundle, metadata, band_name, sub=sub)
                    if data_type == 'Sentinel':
                        band_data = pp.sentinel.get_rtoa(bundle, metadata, bdata, safe_files[granule], band_name, target_res=s2_target_res, sub=grids)
                ## get nans from toa product
                valid_mask = isfinite(band_data)
                rhot_dims = band_data.shape

                ## subset to valid data if slicing
                if slicing: band_data = band_data[valid_mask]

                ## apply gains
                if (gains) & (band_name in gains_dict):
                    band_data *= gains_dict[band_name]
                    print('Applied gain {} for band {}'.format(gains_dict[band_name],band_name))

                ## write toa reflectance
                if nc_write_rhot:
                    ds_att['standard_name']='toa_bidirectional_reflectance'
                    ds_att['units']=1
                    if slicing: # reconstitute full dataset
                        tmp = np.zeros(rhot_dims, dtype=np.float32)+np.nan
                        tmp[valid_mask] = band_data
                        pp.output.nc_write(l2r_ncfile, parname_t, tmp, dataset_attributes=ds_att, new=l2r_nc_new, fillvalue=np.nan,
                                           attributes=attributes, nc_compression=l2r_nc_compression, chunking=chunking)
                        tmp = None
                    else:
                        pp.output.nc_write(l2r_ncfile, parname_t, band_data, dataset_attributes=ds_att, new=l2r_nc_new, fillvalue=np.nan,
                                            attributes=attributes, nc_compression=l2r_nc_compression, chunking=chunking)
                    l2r_nc_new=False

                ## skip some bands
                if (band_name in bands_skip_corr): continue

                ## write Rayleigh corrected reflectance
                if nc_write_rhorc:
                    if aerosol_correction == 'dark_spectrum':
                        if dsf_path_reflectance == 'fixed':
                            rrc_cur = (band_data - rorayl_s[btag]) / (dtotr_s[btag]*utotr_s[btag])
                        else:
                            rorayl_s = pp.ac.tiles_interp(tile_output['atm'][band_name]['rorayl'], xnew, ynew, target_mask=(valid_mask if slicing else None))
                            dtotr_s = pp.ac.tiles_interp(tile_output['atm'][band_name]['dtotr'], xnew, ynew, target_mask=(valid_mask if slicing else None))
                            utotr_s = pp.ac.tiles_interp(tile_output['atm'][band_name]['utotr'], xnew, ynew, target_mask=(valid_mask if slicing else None))
                            rrc_cur = (band_data - rorayl_s) / (dtotr_s*utotr_s)
                            rorayl_s = None
                            dtotr_s = None
                            utotr_s = None
                    elif aerosol_correction == 'exponential':
                        rrc_cur = (band_data - rorayl[btag]) / (dtotr[btag]*utotr[btag])
                    ds_att['standard_name']='rayleigh_corrected_reflectance'
                    ds_att['units']=1
                    rrc_cur[valid_mask == 0] = nan
                    pp.output.nc_write(l2r_ncfile, 'rhorc_{}'.format(wave), rrc_cur, new=l2r_nc_new,  fillvalue=np.nan,
                                        attributes=attributes, dataset_attributes={'wavelength':float(wave),'band_name':band_name})
                    rrc_cur = None
                    l2r_nc_new=False

                ## apply gas correction
                if gas_transmittance:
                    band_data/=tt_gas[band_dict[band_name]['lut_name']]
                    ds_att['tt_gas'] = tt_gas[band_dict[band_name]['lut_name']]

                ## apply sky correction
                if sky_correction:
                    if sky_correction_option == 'all':
                        band_data -= rsky[band_dict[band_name]['lut_name']]
                        ds_att['rsky'] = rsky[band_dict[band_name]['lut_name']]
                    if (sky_correction_option == 'rsky_new') & (aerosol_correction == 'dark_spectrum'):
                        if dsf_path_reflectance == 'tiled':
                            rsky_t_cur = pp.ac.tiles_interp(tile_output['atm'][band_name]['rsky'], xnew, ynew, target_mask=(valid_mask if slicing else None))
                            #rsky_t_cur = pp.ac.tiles_interp(tile_output['atm'][band_name]['rsky'], xnew, ynew)
                            band_data -= rsky_t_cur
                            print('tiled rsky')
                        if dsf_path_reflectance == 'fixed':
                            if resolved_angles:
                                if data_type == 'Sentinel':
                                    print(ac_data[btag].keys())
                                    #rsky_t_cur = pp.ac.tiles_interp(ac_data[btag]['rsky_t'], tp_xnew, tp_ynew)
                                    rsky_t_cur = pp.ac.tiles_interp(ac_data[btag]['rsky_t'], tp_xnew, tp_ynew, target_mask=(valid_mask if slicing else None))

                                if data_type == 'NetCDF':
                                    if resolved_geometry_type == 'average':
                                        #print(lutd[sel_mod]['ipd'])
                                        #tau550_ = 1* tau550
                                        #dtott_s_cur = lutd[sel_mod]['rgi'][band_dict[band_name]['lut_name']]((pressure, lutd[sel_mod]['ipd']['dtott'], raa_ave, vza_ave, sza_ave, tau550_))
                                        #utott_s_cur = lutd[sel_mod]['rgi'][band_dict[band_name]['lut_name']]((pressure, lutd[sel_mod]['ipd']['utott'], raa_ave, vza_ave, sza_ave, tau550_))
                                        #astot_s_cur = lutd[sel_mod]['rgi'][band_dict[band_name]['lut_name']]((pressure, lutd[sel_mod]['ipd']['astot'], raa_ave, vza_ave, sza_ave, tau550_))
                                        #rsky_b_cur = rlutd[sel_aermod]['rgi'][band_dict[band_name]['lut_name']]((raa_ave, vza_ave, sza_ave, tau550_))
                                        #rsky_t_cur = (utott_s_cur * dtott_s_cur * rsky_b_cur)/(1. - rsky_b_cur * astot_s_cur)
                                        #print(band_name, band_dict[band_name]['lut_name'], dtott_s_cur.shape, utott_s_cur.shape, astot_s_cur.shape, rsky_b_cur.shape)
                                        tau550_ = 1* tau550
                                        rsky_t_cur = lutd[sel_mod]['rgi'][band_dict[band_name]['lut_name']]((pressure, lutd[sel_mod]['ipd']['romix+rsky_t'], raa_ave, vza_ave, sza_ave, tau550_))
                                        print(band_name, band_dict[band_name]['lut_name'], rsky_t_cur.shape)
                                band_data -= rsky_t_cur
                            else:
                                band_data -= rsky_toa[band_dict[band_name]['lut_name']]
                                ds_att['rsky'] = rsky_toa[band_dict[band_name]['lut_name']]

                            if (resolved_angles) & (dsf_write_tiled_parameters):
                                if slicing:
                                    tmp = np.zeros(rhot_dims, dtype=np.float32)+np.nan
                                    tmp[valid_mask] = rsky_t_cur
                                    print(rsky_t_cur.shape, tmp.shape)
                                    pp.output.nc_write(l2r_ncfile, 'rsky_{}'.format(wave), tmp, dataset_attributes={'wavelength':float(wave),'band_name':band_name}, new=l2r_nc_new, attributes=attributes, nc_compression=l2r_nc_compression, chunking=chunking)
                                    tmp = None
                                else:
                                    pp.output.nc_write(l2r_ncfile, 'rsky_{}'.format(wave), rsky_t_cur, dataset_attributes={'wavelength':float(wave),'band_name':band_name}, new=l2r_nc_new, attributes=attributes, nc_compression=l2r_nc_compression, chunking=chunking)
                            ## gc from MB
                            rsky_t_cur = None

                ## write surface reflectance
                if (nc_write_rhos | (nc_write_rhow)):
                                        ########################
                    ## dark spectrum fitting
                    if (aerosol_correction == 'dark_spectrum'):
                        ## with tiled aot
                        if dsf_path_reflectance == 'tiled':
                            print('Interpolating a/c tiles for band {}'.format(band_name))
                            ## interpolate tiles to full scene extent
                            #ratm_cur = pp.ac.tiles_interp(tile_output['atm'][band_name]['ratm'], xnew, ynew)
                            #astot_cur = pp.ac.tiles_interp(tile_output['atm'][band_name]['astot'], xnew, ynew)
                            #dutott_cur =  pp.ac.tiles_interp(tile_output['atm'][band_name]['dtott']*tile_output['atm'][band_name]['utott'], xnew, ynew)
                            ratm_cur = pp.ac.tiles_interp(tile_output['atm'][band_name]['ratm'], xnew, ynew, target_mask=(valid_mask if slicing else None))
                            astot_cur = pp.ac.tiles_interp(tile_output['atm'][band_name]['astot'], xnew, ynew, target_mask=(valid_mask if slicing else None))
                            dutott_cur =  pp.ac.tiles_interp(tile_output['atm'][band_name]['dtott']*tile_output['atm'][band_name]['utott'], xnew, ynew, target_mask=(valid_mask if slicing else None))

                            ##
                            ## write the resolved angle grid that was used
                            ## here interpolated for the subtile centres
                            if (b == 0) & (resolved_angles_write):
                                    #tmp =  pp.ac.tiles_interp(tile_angles['VAA'], xnew, ynew)
                                    tmp =  pp.ac.tiles_interp(tile_angles['VAA'], xnew, ynew, target_mask=(valid_mask if slicing else None), target_mask_full=slicing)
                                    pp.output.nc_write(l2r_ncfile, 'view_azimuth', tmp, new=l2r_nc_new, attributes=attributes, nc_compression=l2r_nc_compression, chunking=chunking)
                                    l2r_nc_new=False
                                    #tmp =  pp.ac.tiles_interp(tile_angles['SAA'], xnew, ynew)
                                    tmp =  pp.ac.tiles_interp(tile_angles['SAA'], xnew, ynew, target_mask=(valid_mask if slicing else None), target_mask_full=slicing)
                                    pp.output.nc_write(l2r_ncfile, 'sun_azimuth', tmp, new=l2r_nc_new)
                                    #tmp =  pp.ac.tiles_interp(tile_angles['SZA'], xnew, ynew)
                                    tmp =  pp.ac.tiles_interp(tile_angles['SZA'], xnew, ynew, target_mask=(valid_mask if slicing else None), target_mask_full=slicing)
                                    pp.output.nc_write(l2r_ncfile, 'sun_zenith', tmp, new=l2r_nc_new)
                                    #tmp =  pp.ac.tiles_interp(tile_angles['VZA'], xnew, ynew)
                                    tmp =  pp.ac.tiles_interp(tile_angles['VZA'], xnew, ynew, target_mask=(valid_mask if slicing else None), target_mask_full=slicing)
                                    pp.output.nc_write(l2r_ncfile, 'view_zenith', tmp, new=l2r_nc_new)
                            ## end resolved angle grid

                        ### with fixed path reflectance
                        if dsf_path_reflectance == 'fixed':
                            if resolved_angles:
                                print('Interpolating resolved angle tiles for {}'.format(band_name))
                                if data_type == 'Sentinel':
                                    ## interpolate tiles to full scene extent
                                    ratm_cur = pp.ac.tiles_interp(ac_data[btag]['romix'], tp_xnew, tp_ynew, target_mask=(valid_mask if slicing else None))
                                    astot_cur = pp.ac.tiles_interp(ac_data[btag]['astot'], tp_xnew, tp_ynew, target_mask=(valid_mask if slicing else None))
                                    dutott_cur =  pp.ac.tiles_interp(ac_data[btag]['dtott']*
                                                                     ac_data[btag]['utott'], tp_xnew, tp_ynew, target_mask=(valid_mask if slicing else None))
                                    ##
                                    ## write the resolved angle grid that was used
                                    ## here according to the S2 grid spacing
                                    if (b == 0) & (resolved_angles_write):
                                        tmp =  pp.ac.tiles_interp(grmeta['VIEW']['Average_View_Azimuth'], tp_xnew, tp_ynew, target_mask=(valid_mask if slicing else None), target_mask_full=slicing)
                                        pp.output.nc_write(l2r_ncfile, 'view_azimuth', tmp, new=l2r_nc_new, attributes=attributes, nc_compression=l2r_nc_compression, chunking=chunking)
                                        l2r_nc_new=False
                                        tmp =  pp.ac.tiles_interp(grmeta['VIEW']['Average_View_Zenith'], tp_xnew, tp_ynew, target_mask=(valid_mask if slicing else None), target_mask_full=slicing)
                                        pp.output.nc_write(l2r_ncfile, 'view_zenith', tmp, new=l2r_nc_new, attributes=attributes, nc_compression=l2r_nc_compression, chunking=chunking)
                                        tmp =  pp.ac.tiles_interp(grmeta['SUN']['Azimuth'], tp_xnew, tp_ynew, target_mask=(valid_mask if slicing else None), target_mask_full=slicing)
                                        pp.output.nc_write(l2r_ncfile, 'sun_azimuth', tmp, new=l2r_nc_new, attributes=attributes, nc_compression=l2r_nc_compression, chunking=chunking)
                                        tmp =  pp.ac.tiles_interp(grmeta['SUN']['Zenith'], tp_xnew, tp_ynew, target_mask=(valid_mask if slicing else None), target_mask_full=slicing)
                                        pp.output.nc_write(l2r_ncfile, 'sun_zenith', tmp, new=l2r_nc_new, attributes=attributes, nc_compression=l2r_nc_compression, chunking=chunking)
                                    ## end resolved angle grid

                                if data_type == 'NetCDF':
                                    if resolved_geometry_type == 'average':
                                        print(sel_mod)
                                        tau550_ = 1* tau550
                                        ratm_cur = lutd[sel_mod]['rgi'][band_dict[band_name]['lut_name']]((pressure, lutd[sel_mod]['ipd']['romix'], raa_ave, vza_ave, sza_ave, tau550_))
                                        dutott_cur = lutd[sel_mod]['rgi'][band_dict[band_name]['lut_name']]((pressure, lutd[sel_mod]['ipd']['dtott'], raa_ave, vza_ave, sza_ave, tau550_))*\
                                                     lutd[sel_mod]['rgi'][band_dict[band_name]['lut_name']]((pressure, lutd[sel_mod]['ipd']['utott'], raa_ave, vza_ave, sza_ave, tau550_))
                                        astot_cur = lutd[sel_mod]['rgi'][band_dict[band_name]['lut_name']]((pressure, lutd[sel_mod]['ipd']['astot'], raa_ave, vza_ave, sza_ave, tau550_))
                            else:
                                ratm_cur=ratm_s[btag]
                                astot_cur=astot_s[btag]
                                dutott_cur=utott_s[btag]*dtott_s[btag]

                                ## add atmosphere parameters to attributes
                                ds_att['ratm'] = ratm_s[btag]
                                ds_att['rorayl'] = rorayl_s[btag]
                                ds_att['dtotr'] = dtotr_s[btag]
                                ds_att['utotr'] = utotr_s[btag]
                                ds_att['dtott'] = dtott_s[btag]
                                ds_att['utott'] = utott_s[btag]
                                ds_att['astot'] = astot_s[btag]

                        if (dsf_write_tiled_parameters):
                            if (dsf_path_reflectance == 'tiled') | ((dsf_path_reflectance == 'fixed') & resolved_angles):
                                if slicing:
                                    tmp = np.zeros(rhot_dims, dtype=np.float32)+np.nan
                                    tmp[valid_mask] = ratm_cur
                                    pp.output.nc_write(l2r_ncfile, 'ratm_{}'.format(wave), tmp,
                                                       dataset_attributes={'wavelength':float(wave),'band_name':band_name}, new=l2r_nc_new, attributes=attributes, nc_compression=l2r_nc_compression, chunking=chunking)
                                    tmp[valid_mask] = dutott_cur
                                    pp.output.nc_write(l2r_ncfile, 'dutott_{}'.format(wave), tmp,
                                                       dataset_attributes={'wavelength':float(wave),'band_name':band_name}, new=l2r_nc_new, attributes=attributes, nc_compression=l2r_nc_compression, chunking=chunking)
                                    tmp[valid_mask] = astot_cur
                                    pp.output.nc_write(l2r_ncfile, 'astot_{}'.format(wave), tmp,
                                                       dataset_attributes={'wavelength':float(wave),'band_name':band_name}, new=l2r_nc_new, attributes=attributes, nc_compression=l2r_nc_compression, chunking=chunking)
                                    tmp = None
                                else:
                                    pp.output.nc_write(l2r_ncfile, 'ratm_{}'.format(wave), ratm_cur, dataset_attributes={'wavelength':float(wave),'band_name':band_name}, new=l2r_nc_new, attributes=attributes, nc_compression=l2r_nc_compression, chunking=chunking)
                                    pp.output.nc_write(l2r_ncfile, 'dutott_{}'.format(wave), dutott_cur, dataset_attributes={'wavelength':float(wave),'band_name':band_name}, new=l2r_nc_new, attributes=attributes, nc_compression=l2r_nc_compression, chunking=chunking)
                                    pp.output.nc_write(l2r_ncfile, 'astot_{}'.format(wave), astot_cur, dataset_attributes={'wavelength':float(wave),'band_name':band_name}, new=l2r_nc_new, attributes=attributes, nc_compression=l2r_nc_compression, chunking=chunking)

                                ## output extra ac parameters
                                if extra_ac_parameters:
                                    if (dsf_path_reflectance == 'tiled'):
                                        extra_ac_tags = tile_output['atm'][band_name].keys()
                                    else:
                                        extra_ac_tags = ac_data[btag].keys()

                                    for ap in extra_ac_tags:
                                        if ap in ['ratm', 'astot', 'dutott']: continue
                                        if (dsf_path_reflectance == 'tiled'):
                                            ap_cur =  pp.ac.tiles_interp(tile_output['atm'][band_name][ap], xnew, ynew, target_mask=(valid_mask if slicing else None), target_mask_full=slicing)
                                        else:
                                            ap_cur =  pp.ac.tiles_interp(ac_data[btag][ap], tp_xnew, tp_ynew, target_mask=(valid_mask if slicing else None), target_mask_full=slicing)
                                        if '+' in ap: ap = ap.replace('+', '_')
                                        pp.output.nc_write(l2r_ncfile, '{}_{}'.format(ap, wave), ap_cur, dataset_attributes={'wavelength':float(wave),'band_name':band_name}, new=l2r_nc_new, attributes=attributes, nc_compression=l2r_nc_compression, chunking=chunking)
                                l2r_nc_new = False

                        ## compute surface reflectance
                        ttg_cur = 1.0 ## gas correction should have been performed already
                        rhos_data = (band_data/ttg_cur) - ratm_cur
                        rhos_data = (rhos_data) / (dutott_cur + astot_cur * rhos_data)
                    ## end DSF
                    ########################

                    ########################
                    ## exponential
                    if (aerosol_correction == 'exponential'):
                        eps_cur = pp.ac.exponential_epsilon(epsilon,waves,idx1=short_idx,idx2=long_idx, idxc=b)
                        rhoam_cur = rhoam * eps_cur
                        rhos_data = (band_data - rorayl[btag] - rhoam_cur) / (utotr[btag]*dtotr[btag])
                        rhos_data[mask] = nan
                    ## end exponential
                    ########################

                    ## write surface reflectance
                    if (nc_write_rhos):
                            ds_att['standard_name']='surface_bidirectional_reflectance'
                            ds_att['units']=1
                            if slicing:
                                tmp = np.zeros(rhot_dims, dtype=np.float32)+np.nan
                                tmp[valid_mask] = rhos_data
                                rhos_data = tmp
                                tmp = None
                            else:
                                rhos_data[valid_mask == 0] = nan
                            pp.output.nc_write(l2r_ncfile, parname_s, rhos_data, dataset_attributes=ds_att, new=l2r_nc_new, fillvalue=np.nan,
                                           attributes=attributes, nc_compression=l2r_nc_compression, chunking=chunking)
                            l2r_nc_new = False

                    ## write water reflectance
                    ## add sky corr and mask
                    if (nc_write_rhow):
                            ds_att['standard_name']='water_bidirectional_reflectance'
                            ds_att['units']=1
                            rhos_data[valid_mask == 0] = nan
                            rhos_data[mask] = nan
                            pp.output.nc_write(l2r_ncfile, parname_w, rhos_data, dataset_attributes=ds_att, new=l2r_nc_new, fillvalue=np.nan,
                                               attributes=attributes, nc_compression=l2r_nc_compression, chunking=chunking)
                            l2r_nc_new=False

                ## track rhos negatives
                if l2_negatives is None:
                    l2_negatives = (rhos_data*0).astype(int32)

                if wave < neg_wave:
                    l2_negatives[rhos_data < 0] = 2**1

                rhos_data = None
                band_data = None
                ## end band loop
                ###########################
            ## end writing L2R
            ##############################

            ## write rhos negatives
            pp.output.nc_write(l2r_ncfile, 'l2_negatives', l2_negatives, new=l2r_nc_new, attributes=attributes,
                                   nc_compression=l2r_nc_compression, chunking=chunking)
            l2_negatives = None
            ## end write l2_negatives

            ##############################
            ## write EXP variable rhoam and epsilon
            if (aerosol_correction == 'exponential') & (exp_fixed_epsilon is False):
                rhoam[mask] = nan
                pp.output.nc_write(l2r_ncfile, 'rhoam', rhoam, new=l2r_nc_new, attributes=attributes, fillvalue=np.nan,
                                   nc_compression=l2r_nc_compression, chunking=chunking)
                l2r_nc_new = False
                epsilon[mask] = nan
                pp.output.nc_write(l2r_ncfile, 'epsilon', epsilon, new=l2r_nc_new, attributes=attributes, fillvalue=np.nan,
                                   nc_compression=l2r_nc_compression, chunking=chunking)
            ## end write EXP variable
            ##############################

            t1 = time.time()
            print('Computing surface reflectances took {:.1f} seconds'.format(t1-t0))

            ##########################
            ## write BT if requested datasets
            if (l8_output_bt) & (sensor_family == 'Landsat') & (metadata['SATELLITE'] == 'LANDSAT_8'):
                for band in ['10','11']:
                    parname = 'bt{}'.format(band)
                    ds_att = {'parameter':'at sensor brightness temperature B{}'.format(band)}
                    ds_att['standard_name']='toa_brightness_temperature'
                    ds_att['units']='K'
                    if data_type == 'NetCDF':
                        if parname in pp.nc_datasets(bundle):
                            data = pp.nc_data(bundle, parname)
                        else: continue
                    else:
                        data = pp.landsat.get_bt(bundle, metadata, band, sub=sub)
                    pp.output.nc_write(l2r_ncfile, parname, data, new=l2r_nc_new, attributes=attributes,  fillvalue=np.nan,
                                        dataset_attributes=ds_att, nc_compression=l2r_nc_compression, chunking=chunking)
            ## end write BT
            ##############################

            ##########################
            ## write Lt TIRS if requested
            if (l8_output_lt_tirs) & (sensor_family == 'Landsat') & (metadata['SATELLITE'] == 'LANDSAT_8'):
                for band in ['10','11']:
                    parname = 'lt{}'.format(band)
                    ds_att = {'parameter':'Lt B{}'.format(band)}
                    ds_att['standard_name']='toa_radiance'
                    ds_att['units']='W/(m^2 sr^1 µm^1)'
                    if data_type == 'NetCDF':
                        if parname in pp.nc_datasets(bundle):
                            data = pp.nc_data(bundle, parname)
                        else: continue
                    else:
                        data = pp.landsat.get_bt(bundle, metadata, band, sub=sub, return_radiance=True)
                    pp.output.nc_write(l2r_ncfile, parname, data, new=l2r_nc_new, attributes=attributes,  fillvalue=np.nan,
                                        dataset_attributes=ds_att, nc_compression=l2r_nc_compression, chunking=chunking)
            ## end write Lt TIRS
            ##############################


            ##########################
            ## write latitude and longitude datasets
            if nc_write_geo:
                print('Writing latitude and longitude to outputfiles.')
                if data_type == 'NetCDF':
                    lon = pp.shared.nc_data(granule, 'lon')
                    lat = pp.shared.nc_data(granule, 'lat')
                if data_type == 'Landsat':
                    lon, lat = pp.landsat.geo.get_ll(metadata, limit=limit)
                if data_type == 'Sentinel':
                    lon, lat = pp.sentinel.geo.get_ll(grmeta, limit=limit, resolution=s2_target_res)

                ### write to L1R NCDF
                if os.path.exists(l1r_ncfile) & l1r_read_nc:
                    pp.output.nc_write(l1r_ncfile, 'lon', lon, dataset_attributes={'standard_name':'longitude', 'units':'degree_east'}, double=True)
                    pp.output.nc_write(l1r_ncfile, 'lat', lat, dataset_attributes={'standard_name':'latitude', 'units':'degree_north'}, double=True)
                ### write to L2R NCDF
                if os.path.exists(l2r_ncfile):
                    pp.output.nc_write(l2r_ncfile, 'lon', lon, dataset_attributes={'standard_name':'longitude', 'units':'degree_east'}, double=True)
                    pp.output.nc_write(l2r_ncfile, 'lat', lat, dataset_attributes={'standard_name':'latitude', 'units':'degree_north'}, double=True)
                del lon
                del lat
            ## end write geo
            ##########################


            ##########################
            ## write easting and northing datasets
            if nc_write_geo_xy:
                print('Writing easting (x) and northing (y) to outputfiles.')
                if data_type == 'NetCDF':
                    x = pp.shared.nc_data(granule, 'x')
                    y = pp.shared.nc_data(granule, 'y')
                if data_type == 'Landsat':
                    x, y = pp.landsat.geo.get_ll(metadata, limit=limit, xy=True)
                if data_type == 'Sentinel':
                    x, y = pp.sentinel.geo.get_ll(grmeta, limit=limit, resolution=s2_target_res, xy=True)

                ### write to L1R NCDF
                if os.path.exists(l1r_ncfile):
                    pp.output.nc_write(l1r_ncfile, 'x', x)
                    pp.output.nc_write(l1r_ncfile, 'y', y)
                ### write to L2R NCDF
                if os.path.exists(l2r_ncfile):
                    pp.output.nc_write(l2r_ncfile, 'x', x)
                    pp.output.nc_write(l2r_ncfile, 'y', y)
                del x
                del y
            ## end write geo_xy
            ##########################

            ##########################
            ## write pan band if requested
            if (l8_output_pan | l8_output_pan_ms) & (sensor_family == 'Landsat') & (metadata['SATELLITE'] in ['LANDSAT_7','LANDSAT_8']):
                print('Writing pan band to outputfile.')
                parname = 'rhot_pan'
                parname_ms = 'rhot_pan_ms'
                if data_type == 'NetCDF':
                    tmp = os.path.splitext(bundle)
                    l1_pan_ncdf = '{}_pan{}'.format(tmp[0],tmp[1])
                    l1_pan_ms_ncdf = '{}_pan_ms{}'.format(tmp[0],tmp[1])

                    if os.path.exists(l1_pan_ncdf):
                        data = pp.nc_data(l1_pan_ncdf, parname)
                        if l8_output_pan:
                            pp.output.nc_write(l1r_ncfile_pan, parname, data, new=True, fillvalue=np.nan, nc_compression=l1r_nc_compression, chunking=chunking)
                    else:
                        print('Could not find L1 pan NetCDF.')

                    if os.path.exists(l1_pan_ms_ncdf):
                        data = pp.nc_data(l1_pan_ms_ncdf, parname_ms)
                        if l8_output_pan_ms:
                            pp.output.nc_write(l1r_ncfile_pan_ms, parname_ms, data, new=True, fillvalue=np.nan, nc_compression=l1r_nc_compression, chunking=chunking)
                    else:
                        print('Could not find L1 pan ms NetCDF.')
                else:
                    if metadata['SATELLITE'] == 'LANDSAT_7': panband = 8
                    if metadata['SATELLITE'] == 'LANDSAT_8': panband = 8
                    data = pp.landsat.get_rtoa(bundle, metadata, panband, sub=sub, pan=True)

                    ## fix scenes where number of pan pixels is off
                    if sub is None:
                        pan_ms_dims = global_dims
                        pan_dims = global_dims[0]*2, global_dims[1]*2
                    else:
                        pan_ms_dims = sub[3], sub[2]
                        pan_dims = sub[3]*2, sub[2]*2
                    if data.shape[0] <  pan_dims[0]:
                        data = np.vstack((data, np.zeros((pan_dims[0]-data.shape[0], data.shape[1]))))
                    elif data.shape[0] >  pan_dims[0]:
                        data = data[0:pan_dims[0], :]
                    if data.shape[1] < pan_dims[1]:
                        data = np.hstack((data, np.zeros((data.shape[0], pan_dims[1]-data.shape[1]))))
                    elif data.shape[1] > pan_dims[1]:
                        data = data[:, 0:pan_dims[1]]

                    if data is not None:
                        if l8_output_pan:
                            pp.output.nc_write(l1r_ncfile_pan, parname, data, new=True, fillvalue=np.nan, nc_compression=l1r_nc_compression, chunking=chunking)

                        if l8_output_pan_ms:
                            parname = 'rhot_pan_ms'
                            data = zoom(data, zoom=0.5, order=1)
                            pp.output.nc_write(l1r_ncfile_pan_ms, parname, data, new=True, fillvalue=np.nan, nc_compression=l1r_nc_compression, chunking=chunking)
                data = None
            ## end write pan
            ##############################

            ##########################
            ## write orange band if requested
            if (l8_output_orange) & (sensor_family == 'Landsat') & (metadata['SATELLITE'] == 'LANDSAT_8') & (aerosol_correction == 'dark_spectrum'):
                if os.path.exists(l1r_ncfile_pan_ms):
                    print('Calculating orange band.')
                    pan_ms = pp.shared.nc_data(l1r_ncfile_pan_ms, 'rhot_pan_ms')
                    pan_wave = swaves['8']

                    ## get nans from toa product
                    valid_mask = isfinite(pan_ms)

                    ## orange band config
                    ob_sensor = 'L8_OLI_ORANGE'
                    ob_rsr_file = pp.config['pp_data_dir']+'/RSR/'+ob_sensor+'.txt'
                    ob_rsr, ob_rsr_bands = pp.rsr_read(file=ob_rsr_file)
                    owave = pp.shared.rsr_convolute(ob_rsr['O']['wave'], ob_rsr['O']['wave'], ob_rsr['O']['response'], ob_rsr['O']['wave'])
                    ob_wave = '{:.0f}'.format(owave*1000.)
                    ds_att = {'wavelength': float(ob_wave), 'band_name':"O"}
                    btag = 'O'
                    wave_red,wave_green,wave_orange = 655,561,613

                    ## get orange band config
                    if ob_cfg is None:
                        ob_cfg_file = pp.config['pp_data_dir']+'/Shared/oli_orange.cfg'
                        ob_cfg = pp.shared.import_config(ob_cfg_file)

                    if ob_cfg['combine'] == 'after':
                        ob_ds = 'rhos_{}'
                        ptag = '8'
                        ds_att_pan = {'wavelength': float(pan_wave), 'band_name':"8",
                                      'tt_gas':tt_gas[ptag]}#, 'rsky':rsky_toa[ptag]}

                        # pan band gas transmittance
                        if gas_transmittance:
                            pan_ms/=tt_gas[ptag]

                        # pan band sky reflectance
                        if sky_correction:
                            if sky_correction_option == 'all':
                                pan_ms-=rsky[ptag]
                            if sky_correction_option == 'rsky_new':
                                if dsf_path_reflectance == 'tiled':
                                    rsky_t_cur = pp.ac.tiles_interp(tile_output['atm'][ptag]['rsky'], xnew, ynew, target_mask=(valid_mask if slicing else None), target_mask_full=slicing)
                                    pan_ms -= rsky_t_cur
                                if dsf_path_reflectance == 'fixed':
                                    pan_ms -= rsky_toa[ptag]
                                    ds_att_pan['rsky'] =rsky_toa[ptag]

                        # fixed or tiled DSF to get pan band rhos
                        if dsf_path_reflectance == 'fixed':
                            ## pan band
                            ds_att_pan['ratm'] = ratm_s[ptag]
                            ds_att_pan['rorayl'] = rorayl_s[ptag]
                            ds_att_pan['dtotr'] = dtotr_s[ptag]
                            ds_att_pan['utotr'] = utotr_s[ptag]
                            ds_att_pan['dtott'] = dtott_s[ptag]
                            ds_att_pan['utott'] = utott_s[ptag]
                            ds_att_pan['astot'] = astot_s[ptag]
                            pan_ms = pp.rtoa_to_rhos(pan_ms, ratm_s[ptag], utott_s[ptag], dtott_s[ptag], astot_s[ptag], tt_gas = 1)
                        else:
                            print('tiled orange band!')
                            ttg_cur = 1.0 ## gas correction done above
                            ## interpolate tiles to full scene extent
                            ratm_cur = pp.ac.tiles_interp(tile_output['atm'][ptag]['ratm'], xnew, ynew)
                            astot_cur = pp.ac.tiles_interp(tile_output['atm'][ptag]['astot'], xnew, ynew)
                            dutott_cur =  pp.ac.tiles_interp(tile_output['atm'][ptag]['dtott']*tile_output['atm'][ptag]['utott'], xnew, ynew)
                            pan_ms = (pan_ms/ttg_cur) - ratm_cur
                            pan_ms = (pan_ms) / (dutott_cur + astot_cur * pan_ms)

                        ## write surface reflectance
                        rhos_data = pan_ms * 1.0
                        rhos_data[valid_mask == 0] = nan
                        pp.output.nc_write(l2r_ncfile, 'rhos_{}'.format(pan_wave), rhos_data, dataset_attributes=ds_att_pan)
                        rhos_data = None
                    else:
                        ob_ds = 'rhot_{}'

                    ## compute orange band
                    green = pp.shared.nc_data(l2r_ncfile, ob_ds.format(wave_green))
                    red = pp.shared.nc_data(l2r_ncfile, ob_ds.format(wave_red))
                    valid_mask *= isfinite(green) * isfinite(red)

                    if ob_cfg['algorithm'] == 'A':
                        ds_att['algorithm'] = 'A'
                        ## green_p = (green * a_gf[0] + a_gf[1]) * a_gf[2]
                        ## red_p   = (red * a_rf[0] + a_rf[1]) * a_rf[2]
                        ## orange  = (pan_ms - (green_p + red_p)) / a_of
                        a_gf = [float(f.strip()) for f in ob_cfg['a_gf'].split(',')]
                        a_rf = [float(f.strip()) for f in ob_cfg['a_rf'].split(',')]
                        a_of = float(ob_cfg['a_of'])
                        green_p = (green * a_gf[0] + a_gf[1]) * a_gf[2]
                        red_p = (red * a_rf[0] + a_rf[1]) * a_rf[2]
                        band_data =  (pan_ms - (green_p + red_p)) / a_of
                    elif ob_cfg['algorithm'] == 'B':
                        ds_att['algorithm'] = 'B'
                        ## orange = b_pf * pan_ms + b_gf * green + b_rf * red
                        b_pf = float(ob_cfg['b_pf'])
                        b_gf = float(ob_cfg['b_gf'])
                        b_rf = float(ob_cfg['b_rf'])
                        band_data = b_pf * pan_ms + b_gf * green + b_rf * red

                    ## with fixed DSF
                    if ob_cfg['combine'] == 'before':
                        ## write TOA data
                        pp.output.nc_write(l2r_ncfile, 'rhot_{}'.format(ob_wave), band_data, new=l2r_nc_new, fillvalue=np.nan, attributes=attributes, dataset_attributes=ds_att)
                        new=False

                        if dsf_path_reflectance == 'fixed':
                            ac_model = attributes['ac_model']
                            if type(ac_model) == list: ac_model=ac_model[0]
                            ob_lut = ac_model.split('-')[0:4] + ['1013mb']
                            ob_lut = '-'.join(ob_lut)
                            ## get orange band LUT
                            ob_lutdir = pp.config['pp_data_dir']+'/LUT'
                            ob_lut_sensor, ob_meta_sensor = pp.aerlut.aerlut_pressure(ob_lut, ob_lutdir, attributes['pressure'], ob_sensor, ob_rsr_file)
                            ## get atmospheric correction parameters
                            ratm_o,rorayl_o,dtotr_o,utotr_o,dtott_o,utott_o,astot_o=\
                            pp.aerlut.lut_get_ac_parameters_fixed_tau_sensor(ob_lut_sensor,ob_meta_sensor,attributes['AZI'],attributes['THV'],attributes['THS'],attributes['ac_aot550'])
                            ## add atmosphere parameters to attributes
                            ds_att['ratm'] = ratm_o[btag]
                            ds_att['rorayl'] = rorayl_o[btag]
                            ds_att['dtotr'] = dtotr_o[btag]
                            ds_att['utotr'] = utotr_o[btag]
                            ds_att['dtott'] = dtott_o[btag]
                            ds_att['utott'] = utott_o[btag]
                            ds_att['astot'] = astot_o[btag]
                        elif dsf_path_reflectance == 'tiled':
                            ttg_cur = 1.0 ## gas correction done below
                            ## interpolate tiles to full scene extent
                            ratm_cur = pp.ac.tiles_interp(tile_output['atm'][btag]['ratm'], xnew, ynew)
                            astot_cur = pp.ac.tiles_interp(tile_output['atm'][btag]['astot'], xnew, ynew)
                            dutott_cur =  pp.ac.tiles_interp(tile_output['atm'][btag]['dtott']*tile_output['atm'][btag]['utott'], xnew, ynew)

                        ## write Rayleigh corrected reflectance
                        if nc_write_rhorc:
                            rrc_cur = (band_data - rorayl_o[btag]) / (dtotr_o[btag]*utotr_o[btag])
                            rrc_cur[valid_mask == 0] = nan
                            pp.output.nc_write(l2r_ncfile, 'rhorc_{}'.format(ob_wave), rrc_cur, new=new, attributes=attributes, fillvalue=np.nan, dataset_attributes=ds_att)
                            rrc_cur = None
                            new=False

                        if gas_transmittance:
                            ob_tt_oz = pp.ac.o3_transmittance(ob_sensor, metadata, uoz=uoz)
                            ob_tt_wv = pp.ac.wvlut_interp(attributes['THS'], attributes['THV'], uwv=uwv, sensor=ob_sensor, config=wvlut)
                            ob_tt_o2 = pp.ac.o2lut_interp(attributes['THS'], attributes['THV'], sensor=ob_sensor)
                            ob_tt_gas = {btag: ob_tt_oz[btag] * ob_tt_wv[btag] * ob_tt_o2[btag] for btag in ob_tt_oz.keys()}
                            band_data /= ob_tt_gas[btag]

                        if sky_correction:
                            metadata_O =metadata.copy()
                            metadata_O["SATELLITE_SENSOR"] = ob_sensor
                            rsky_O = pp.ac.toa_rsky(metadata_O, pressure=pressure, lut=luts[0])
                            if sky_correction_option == 'all':
                                band_data -= rsky_O[btag]
                            if sky_correction_option == 'rsky_new':
                                rsky_ret_O = pp.shared.rsr_convolute_dict(waves, cur_rsky, ob_rsr)
                                band_data -= rsky_ret_O[btag]
                                ds_att['rsky'] =rsky_ret_O[btag]

                        ## compute surface reflectance
                        if dsf_path_reflectance == 'fixed':
                            rhos_data = pp.rtoa_to_rhos(band_data, ratm_o[btag], utott_o[btag], dtott_o[btag], astot_o[btag], tt_gas = 1)
                        elif dsf_path_reflectance == 'tiled':
                            rhos_data = (band_data/ttg_cur) - ratm_cur
                            rhos_data = (rhos_data) / (dutott_cur + astot_cur * rhos_data)
                            if dsf_write_tiled_parameters:
                                pp.output.nc_write(l2r_ncfile, 'ratm_{}'.format(ob_wave), ratm_cur, dataset_attributes=ds_att, new=l2r_nc_new, attributes=attributes, nc_compression=l2r_nc_compression, chunking=chunking)
                                pp.output.nc_write(l2r_ncfile, 'dutott_{}'.format(ob_wave), dutott_cur, dataset_attributes=ds_att, new=l2r_nc_new, attributes=attributes, nc_compression=l2r_nc_compression, chunking=chunking)
                                pp.output.nc_write(l2r_ncfile, 'astot_{}'.format(ob_wave), astot_cur, dataset_attributes=ds_att, new=l2r_nc_new, attributes=attributes, nc_compression=l2r_nc_compression, chunking=chunking)

                                ## if output extra ac parameters
                                if extra_ac_parameters:
                                    for ap in tile_output['atm'][btag]:
                                        if ap in ['ratm', 'astot', 'dutott']: continue
                                        ap_cur =  pp.ac.tiles_interp(tile_output['atm'][btag][ap], xnew, ynew)
                                        if '+' in ap: ap = ap.replace('+', '_')
                                        pp.output.nc_write(l2r_ncfile, '{}_{}'.format(ap, ob_wave), ap_cur, dataset_attributes=ds_att, new=l2r_nc_new, attributes=attributes, nc_compression=l2r_nc_compression, chunking=chunking)
                                l2r_nc_new = False
                        band_data = None
                    else:
                        rhos_data = band_data * 1.0
                        band_data = None

                    ## write surface reflectance
                    rhos_data[valid_mask == 0] = nan
                    pp.output.nc_write(l2r_ncfile, 'rhos_{}'.format(ob_wave), rhos_data, fillvalue=np.nan, dataset_attributes=ds_att)
                    rhos_data = None
                else:
                    print('L1 pan ms NetCDF file not found')
            ## end write orange band
            ##############################

        ## gc from MB
        ratm_cur = None
        astot_cur = None
        dutott_cur = None

        ## end nc writing
        ####################################

        ####################################
        ## glint correction
        if (aerosol_correction == 'dark_spectrum') & glint_correction:
            t0 = time.time()
            print('Starting glint correction')

            ## compute scattering angle
            from numpy import arccos, sqrt, cos, sin, pi, exp

            dtor = pi / 180.
            ths = attributes['THS'] * dtor
            thv = attributes['THV'] * dtor
            azi = attributes['AZI'] * dtor

            muv = cos(thv)
            mus = cos(ths)

            cos2omega = mus*muv + sin(ths)*sin(thv)*cos(azi)
            omega = arccos(sqrt(cos2omega))
            omega = arccos(cos2omega)/2

            ## read refractive index
            refri = pp.shared.read_refri()

            ## compute fresnel reflectance for each n
            Rf = [pp.ac.sky_refl(omega, n_w=n) for n in refri['n']]

            ## convolute to sensor bands
            rsr_file = pp.config['pp_data_dir']+'/RSR/'+attributes['sensor']+'.txt'
            rsr, rsr_bands = pp.shared.rsr_read(file=rsr_file)
            wave = [w/1000. for w in refri['wave']]
            Rf_sen = pp.shared.rsr_convolute_dict(wave, Rf, rsr)

            ## get SWIR waves
            gc_waves = [band_dict[b]['wave'] for b in ordered_bands]
            gc_swir1_idx, gc_swir1_wv = pp.shared.closest_idx(gc_waves, 1650.)
            gc_swir2_idx, gc_swir2_wv = pp.shared.closest_idx(gc_waves, 2200.)
            gc_swir1_band = band_dict[ordered_bands[gc_swir1_idx]]['lut_name']
            gc_swir2_band = band_dict[ordered_bands[gc_swir2_idx]]['lut_name']
            gc_t_keep = [gc_swir1_band, gc_swir2_band]

            gc_user_band = None
            if glint_force_band is not None:
                gc_user_idx, gc_user_wv = pp.shared.closest_idx(gc_waves, float(glint_force_band))
                gc_user_band = band_dict[ordered_bands[gc_user_idx]]['lut_name']
                gc_t_keep.append(gc_user_band)

            ## get total atmosphere optical thickness
            if dsf_path_reflectance == 'fixed':
                if not resolved_angles:
                    ttot = {band:float(lutd[sel_mod]['rgi'][band]((pressure, lutd[sel_mod]['ipd']['ttot'],
                                                                   raa, vza, sza, attributes['ac_aot550']))) for band in rsr_bands}


            ## QV 2020-11-16
            ## updated subsetting

            ## get swir threshold for glint correction
            gc_mask_idx, gc_mask_wave = gc_swir1_idx, gc_swir1_wv = pp.shared.closest_idx(gc_waves, glint_mask_rhos_band)
            glint_ref_rhos = pp.shared.nc_data(l2r_ncfile, 'rhos_{}'.format('{:.0f}'.format(gc_waves[gc_mask_idx])))

            ## sub_gc has the idx for non masked data with rhos_ref below the masking threshold
            sub_gc = np.where((~glint_ref_rhos.mask) & (glint_ref_rhos.data<=glint_mask_rhos_threshold))
            glint_ref_rhos = None

            ## glint mask for tile interpolation
            if (dsf_path_reflectance == 'tiled') | (resolved_angles & (data_type == 'Sentinel')):
                glint_mask = valid_mask * False
                glint_mask[sub_gc] = True

            ## compute T for glint correction reference bands
            for b,band_name in enumerate(ordered_bands):
                if band_name in bands_skip_thermal: continue
                if band_name in bands_skip_corr: continue
                btag = band_dict[band_name]['lut_name']

                if glint_force_band is not None:
                    if btag not in [gc_user_band]: continue
                else:
                    if btag not in [gc_swir1_band, gc_swir2_band]: continue

                ## direct up and down transmittance
                if dsf_path_reflectance == 'fixed':
                    if resolved_angles:
                        if data_type == 'Sentinel':
                            #ttot_cur = pp.ac.tiles_interp(ac_data[btag]['ttot'], tp_xnew, tp_ynew)[sub_gc]
                            ttot_cur = pp.ac.tiles_interp(ac_data[btag]['ttot'], tp_xnew, tp_ynew, target_mask=glint_mask)
                        if data_type == 'NetCDF':
                            if resolved_geometry_type == 'average':
                                tau550_ = 1* tau550
                                ttot_cur = lutd[sel_mod]['rgi'][band_dict[band_name]['lut_name']]((pressure, lutd[sel_mod]['ipd']['ttot'], raa_ave, vza_ave, sza_ave, tau550_))
                    else:
                        ttot_cur = ttot[btag]
                else: ## tiled ttot
                    #ttot_cur = pp.ac.tiles_interp(tile_output['atm'][band_name]['ttot'], xnew, ynew)[sub_gc]
                    ttot_cur = pp.ac.tiles_interp(tile_output['atm'][band_name]['ttot'], xnew, ynew, target_mask=glint_mask)

                ## two way direct transmittance
                T_cur  = exp(-1.*(ttot_cur/muv)) * exp(-1.*(ttot_cur/mus))
                ttot_cur = None

                if glint_force_band is not None:
                    T_USER = T_cur * 1.0
                else:
                    if btag == gc_swir1_band:
                        T_SWIR1 = T_cur * 1.0
                    if btag == gc_swir2_band:
                        T_SWIR2 = T_cur * 1.0
                T_cur = None


            ## compute glint correction factors and correct per band
            for b,band_name in enumerate(ordered_bands):
                if band_name in bands_skip_thermal: continue
                if band_name in bands_skip_corr: continue
                ## set up band parameter
                btag = band_dict[band_name]['lut_name']
                wave = band_dict[band_name]['wave']

                print('Performing glint correction for band {} ({} nm)'.format(band_name, wave))

                ## read rhos
                cur_rhos = pp.shared.nc_data(l2r_ncfile, 'rhos_{}'.format(wave))

                ## compute required parameters for this band
                ## direct up and down transmittance
                if (glint_force_band is None) & (btag == gc_swir1_band):
                    T_cur = T_SWIR1 * 1.0
                elif (glint_force_band is None) & (btag == gc_swir2_band):
                    T_cur = T_SWIR2 * 1.0
                elif (glint_force_band is not None) & (btag == gc_user_band):
                    T_cur = T_USER * 1.0
                else:
                    print('Computing transmittance for band {}'.format(band_name))
                    if dsf_path_reflectance == 'fixed':
                        if resolved_angles:
                            if data_type == 'Sentinel':
                                #ttot_cur = pp.ac.tiles_interp(ac_data[btag]['ttot'], tp_xnew, tp_ynew)[sub_gc]
                                ttot_cur = pp.ac.tiles_interp(ac_data[btag]['ttot'], tp_xnew, tp_ynew, target_mask=glint_mask)
                            if data_type == 'NetCDF':
                                if resolved_geometry_type == 'average':
                                    tau550_ = 1* tau550
                                    ttot_cur = lutd[sel_mod]['rgi'][band_dict[band_name]['lut_name']]((pressure, lutd[sel_mod]['ipd']['ttot'], raa_ave, vza_ave, sza_ave, tau550_))#[sub_gc]
                        else:
                            ttot_cur = ttot[btag]
                    else: ## tiled ttot
                        #ttot_cur = pp.ac.tiles_interp(tile_output['atm'][band_name]['ttot'], xnew, ynew)[sub_gc]
                        ttot_cur = pp.ac.tiles_interp(tile_output['atm'][band_name]['ttot'], xnew, ynew, target_mask=glint_mask)

                    ## two way direct transmittance
                    T_cur  = exp(-1.*(ttot_cur/muv)) * exp(-1.*(ttot_cur/mus))

                    ttot_cur = None

                if glint_force_band is None:
                    gc_SWIR1 = (T_cur/T_SWIR1) * (Rf_sen[btag]/Rf_sen[gc_swir1_band])
                    gc_SWIR2 = (T_cur/T_SWIR2) * (Rf_sen[btag]/Rf_sen[gc_swir2_band])
                    #print(band_name,gc_SWIR1,gc_SWIR2)
                else:
                    gc_USER = (T_cur/T_USER) * (Rf_sen[btag]/Rf_sen[gc_user_band])
                T_cur = None

                ## choose SWIR band with lowest values
                ## use blue band results
                if wave == min(ordered_waves):
                    ## read in rhos
                    ## from SWIR bands
                    if glint_force_band is None:
                        swir1_rhos = pp.shared.nc_data(l2r_ncfile, 'rhos_{}'.format('{:.0f}'.format(gc_waves[gc_swir1_idx])))[sub_gc]
                        swir2_rhos = pp.shared.nc_data(l2r_ncfile, 'rhos_{}'.format('{:.0f}'.format(gc_waves[gc_swir2_idx])))[sub_gc]
                        ## set negatives to 0
                        swir1_rhos[swir1_rhos<0] = 0
                        swir2_rhos[swir2_rhos<0] = 0

                        ## estimate glint correction in the blue band
                        g1_blue = gc_SWIR1 * swir1_rhos
                        g2_blue = gc_SWIR2 * swir2_rhos

                        ## use SWIR1 or SWIR2 based glint correction
                        use_swir1 = np.where(g1_blue<g2_blue)

                        rhog_ref = swir2_rhos
                        rhog_ref[use_swir1] = swir1_rhos[use_swir1]
                        swir1_rhos, swir2_rhos = None, None
                        g1_blue, g2_blue = None, None
                    ## otherwise from user defined band
                    else:
                        rhog_ref = pp.shared.nc_data(l2r_ncfile, 'rhos_{}'.format('{:.0f}'.format(gc_waves[gc_user_idx])))[sub_gc]
                        ## set negatives to 0
                        rhog_ref[rhog_ref<0] = 0

                    ## write reference glint
                    if glint_write_rhog_ref:
                        rhog_ref_full = np.zeros(rhot_dims, dtype=np.float32) + np.nan
                        rhog_ref_full[sub_gc] = rhog_ref
                        pp.output.nc_write(l2r_ncfile, 'rhog_ref', rhog_ref_full, fillvalue=np.nan)
                        rhog_ref_full = None
                ### end choose band for GC

                ## estimate current band glint
                if glint_force_band is None:
                    cur_rhog = gc_SWIR2 * rhog_ref
                    try:
                        cur_rhog[use_swir1] = gc_SWIR1[use_swir1] * rhog_ref[use_swir1]
                    except:
                        cur_rhog[use_swir1] = gc_SWIR1 * rhog_ref[use_swir1]
                else:
                    cur_rhog = gc_USER * rhog_ref

                ## write full scene glint dataset
                if glint_write_rhog_all:
                    cur_rhog_full = np.zeros(rhot_dims, dtype=np.float32) + np.nan
                    cur_rhog_full[sub_gc] = cur_rhog
                    pp.output.nc_write(l2r_ncfile, 'rhog_{}'.format(wave), cur_rhog_full, fillvalue=np.nan)
                    cur_rhog_full = None

                cur_rhos[sub_gc]-=cur_rhog
                pp.output.nc_write(l2r_ncfile, 'rhos_{}'.format(wave), cur_rhos, fillvalue=np.nan)
                cur_rhos, cur_rhog = None, None
            t1 = time.time()
            print('Glint correction took {:.1f} seconds'.format(t1-t0))
       ## end glint correction
       ####################################

        ## remove l1r_nc_files
        if l1r_nc_delete:
            for f in [l1r_ncfile, l1r_ncfile_pan, l1r_ncfile_pan_ms]:
                if os.path.exists(f): os.remove(f)

        ## remove nc file
        if (nc_delete) & (os.path.exists(l2r_ncfile)):
            os.remove(l2r_ncfile)
        else:
            l2r_files.append(l2r_ncfile)

    lutd = None

    return(l2r_files)
