## def select_model
## initialised interpolators for geolocation of Pleiades image
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2016-07
## modifications: QV 2017-01-18 converted to function
##                QV 2017-01-27 added vis bestfit option
##                QV 2017-04-13 added landsat 8 config
##                QV 2017-04-18 added Sentinel-2A config and updated to 201704 LUTs
##                QV 2017-04-19 addedd bestfit 'bands' option for user specified bands
##                QV 2017-07-18 added pressure keyword and linear weighting of LUTs according to pressure
##                2017-10-19 (QV) added force_band option
##                2017-10-24 (QV) switched Landsat default bands to SWIR
##                2017-11-13 (QV) added Pléiades correct view azimuth
##                2017-11-21 (QV) rewrite of geometry and defaults extraction (should now be standardised in metadata)
##                                merged in rtoa 'list' and 'list_smoothed' methods (with some keywords)
##                                added model_selection keyword: 'min_tau' (choose model giving lowest tau, default)
##                                                               'min_rmsd' (choose model giving lowest rmsd for selected bands)
##                                added rdark_list_selection keyword (when len(rdark) != 1): 'intercept' (OLS intercept, default), 
##                                                                                      'smooth' (lowess smoothing)
##                                                                                      'list' (min rmsd in list - not implemented here)
##                2017-11-22 (QV) replace ratm by rorayl when ratm == nan
##                2017-11-28 (QV) moved PP data directory
##                2018-01-25 (QV) added MaskedArray as rdark type, and check rmsd return for being masked
##                2018-03-06 (QV) added new delta AOT method for best model selection (model_selection='min_dtau')
##                2018-03-12 (QV) changed pixel_range to floats for Windows processing
##                2018-03-14 (QV) added lut sensor metadata check
##                2018-07-18 (QV) changed acolite import name
##                2018-10-01 (QV) removed obsolete bits

def select_model(metadata, rdark, rsr_file=None, lutdir=None, bestfit='bands', bestfit_bands=None, 
                 force_band=None, pressure=None, model_selection='min_tau', rdark_list_selection='intercept',
                 lowess_frac = 0.5, lut_data_dict=None, 
                 luts=['PONDER-LUT-201704-MOD1-1013mb', 'PONDER-LUT-201704-MOD2-1013mb', 'PONDER-LUT-201704-MOD3-1013mb']):
    
    import acolite as pp
    import os
    from numpy import float64, ndarray, arange, nanmin, where, nan, isnan, min
    from numpy.ma import MaskedArray,is_masked
    from statsmodels.nonparametric.smoothers_lowess import lowess

    ## get scene geometry and default bands from metadata
    try:
        if 'SE_DISTANCE' in metadata.keys():
            se_distance = metadata['SE_DISTANCE']
        else:
            se_distance = pp.distance_se(metadata['DOY'])
        ths = metadata['THS']
        thv = metadata['THV']
        azi = metadata['AZI']
        if 'LUT_SENSOR' in metadata.keys():
            sensor = metadata['LUT_SENSOR']
        elif 'SATELLITE_SENSOR' in metadata.keys():
            sensor = metadata['SATELLITE_SENSOR']
        else:
            sensor = metadata['SENSOR']
        bestfit_bands_defaults = metadata['BANDS_BESTFIT']
        bands_sorted = metadata['BANDS_ALL']
    except:
        print('Could not get appropriate metadata for model selection for satellite {}'.format(metadata['SATELLITE']))
        print(metadata.keys())
        return(1)
    
    bestfit_bands_all = rdark.keys()
    if bestfit_bands is None: bestfit_bands = bestfit_bands_defaults
    if bestfit_bands == 'default': bestfit_bands = bestfit_bands_defaults
    if bestfit_bands == 'all': bestfit_bands = bestfit_bands_all

    ## set LUT dir and rsr_file
    pp_path = pp.config['pp_data_dir']
    if lutdir is None:
        lutdir=pp_path+'/LUT/'
    if rsr_file is None:
        rsr_file = pp_path+'/RSR/'+sensor+'.txt'
    rsr, rsr_bands = pp.shared.rsr_read(file=rsr_file)

    ## make lut data dictionary that can be reused in next runs
    if lut_data_dict is None:
        lut_data_dict = {}
        for li,lut in enumerate(luts):
                ## get sensor LUT
                lut_sensor, meta_sensor = pp.aerlut.get_sensor_lut(sensor, rsr_file, lutdir=lutdir, lutid=lut, override=0)
                lut_data_dict[lut] = {'lut':lut_sensor, 'meta':meta_sensor}

                ## read luts at other pressures if needed
                if pressure is not None:
                    lut_split = lut.split('-')
                    lut0 = '-'.join(lut_split[0:-1]+['0500mb'])
                    lut_sensor, meta_sensor = pp.aerlut.get_sensor_lut(sensor, rsr_file, lutdir=lutdir, lutid=lut0, override=0)
                    lut_data_dict[lut0] = {'lut':lut_sensor, 'meta':meta_sensor}

                    lut1 = '-'.join(lut_split[0:-1]+['1100mb'])
                    lut_sensor, meta_sensor = pp.aerlut.get_sensor_lut(sensor, rsr_file, lutdir=lutdir, lutid=lut1, override=0)
                    lut_data_dict[lut1] = {'lut':lut_sensor, 'meta':meta_sensor}

    ## empty dicts
    rdark_smooth = {}
    rdark_selected = {}
    
    ## default is not to use the 'list' type model/rdark selection
    ## will only be used if the rdark is a list and the rdark_list_selection=="list"
    rdark_list = False

    ## find out what kind of rdark is given and which selection to use
    for band in rdark.keys():
        ## given rdark is a single spectrum
        if (type(rdark[band]) == float64):
            rdark_selected[band] = rdark[band]
            
        ## given rdark is a list of spectra
        elif (type(rdark[band]) == ndarray) or (type(rdark[band]) == MaskedArray):
            #pixel_range = range(0, len(rdark[band]))
            pixel_range = [float(i) for i in range(0, len(rdark[band]))]
            #print(pixel_range)
            rdark_list = True

            ## select rdark by lowess smoothing the given spectra
            if rdark_list_selection == 'smooth':
                rdark_smooth[band] = lowess(rdark[band],pixel_range,frac=lowess_frac)[:,1]
                rdark_selected[band] = rdark_smooth[band][0]
                rdark_list = False
                
            ## select rdark by OLS intercept
            elif rdark_list_selection == 'intercept':
                for band in rdark.keys():
                    reg = pp.shared.regression.lsqfity(pixel_range, rdark[band]) # m, b, r, sm, sb
                    rdark_selected[band] = reg[1]
                    rdark_list = False
            
            ## select rdark from all given pixels according to min rmsd
            elif rdark_list_selection == 'list':
                rdark_selected[band] = rdark[band]
                rdark_list = True
                
            ## select rdark from all given pixels according to min rmsd, but first smooth the given rdark
            elif rdark_list_selection == 'list_smooth':
                rdark_selected[band] = lowess(rdark[band],pixel_range,frac=lowess_frac)[:,1]
                rdark_list = True

            ## fallback selection is darkest pixel in each band
            else:
                rdark_list_selection = 'darkest'
                rdark_selected[band] = nanmin(rdark[band])
                
        ## given dark is something else
        else:
            print('rdark type not recognised')
            rdark_selected = (rdark)

    if rdark_list:
                from numpy import min,max
                rdark_len = [len(rdark[b]) for b in rdark]
                max_len = max(rdark_len)
                min_len = min(rdark_len)
                rdark = {b:rdark[b][0:min_len] for b in rdark}

    ## find best fitting model
    sel_rmsd = 1.
    sel_tau = 5.
    taufits = []
    
    daot_minimum = 5

    ## run through luts
    for li,lut in enumerate(luts):
        
        ## get current sensor LUT
        if pressure is not None: ## interpolate LUTs to given pressure
            lut_sensor, meta_sensor = pp.aerlut.aerlut_pressure(lut, lutdir, pressure, sensor, rsr_file, lut_data_dict=lut_data_dict)
        else: ## just the default LUTs
            lut_sensor, meta_sensor = (lut_data_dict[lut]["lut"], lut_data_dict[lut]["meta"])
        
        ## get tau for selected dark spectrum
        tau_550, tau_rmsd, tau_band = pp.aerlut.lut_get_taufit_sensor(lut_sensor, meta_sensor,azi,thv,ths,rdark_selected, 
                                                                      bestfit_bands=bestfit_bands, force_band=force_band)
        
        ## if rdark list method do not select tau and model here, just append to 'taufits' list for later processing
        if (rdark_list):
            taufits.append((tau_550, tau_rmsd, tau_band))
            
        ## if single spectra selected before, select tau and model here
        else:
            ## get tau and ac parameters for this model
            tmp = pp.aerlut.lut_get_ac_parameters_sensor(lut_sensor,meta_sensor,azi,thv,ths,
                                                        rdark_selected, force_band=force_band)
            

            ## band index for this model (i.e. band giving lowest tau)
            idx = tmp[8]
            ## selected tau for this model
            tau550_cur = tmp[7][idx]
            
            ## probably obsolete
            if bestfit == 'bands':
                rmsd_y = [rdark_selected[band] for band in bestfit_bands]
                rmsd_x = [tmp[0][band] for band in bestfit_bands]
                rmsd = pp.rmsd(rmsd_x,rmsd_y)
            else:
                rmsd = tmp[9][idx]
                
            ## find best spectral fit
            if model_selection == 'min_rmsd':
                if is_masked(rmsd):
                    rmsd = 1.0

                if (rmsd <= sel_rmsd) | (li == 0):
                    sel_rmsd = rmsd
                    sel_idx = idx
                    sel_ac_par = tmp
                    sel_model_lut, sel_model_lut_meta = (lut_sensor, meta_sensor)
                    pixel_idx = -1
            
            ## find lowest tau fit
            if model_selection == 'min_tau':
                if (tau550_cur <= sel_tau) | (li == 0):
                    sel_tau = tau550_cur
                    sel_rmsd = rmsd
                    sel_idx = idx
                    sel_ac_par = tmp
                    sel_model_lut, sel_model_lut_meta = (lut_sensor, meta_sensor)
                    pixel_idx = -1
                    
            ## find lowest tau and second band with most similar tau
            if model_selection == 'min_dtau':
                allt = [tmp[7][b] for i,b in enumerate(tmp[7].keys())]
                min_tau = min([t for t in allt if t > 0.001])
                allt_diff_sorted = [t-min_tau for t in allt if t > min_tau]
                allt_diff_sorted.sort()
                if allt_diff_sorted[0] < daot_minimum:
                    daot_minimum = allt_diff_sorted[0]
                    daot_second_band = [b for i,b in enumerate(tmp[7].keys()) \
                                        if allt[i] == allt_diff_sorted[0] + min_tau][0]
                    
                    sel_tau = tau550_cur
                    sel_rmsd = rmsd
                    sel_idx = (idx,daot_second_band)
                    sel_ac_par = tmp
                    sel_model_lut, sel_model_lut_meta = (lut_sensor, meta_sensor)
                    pixel_idx = -1

            ## find lowest tau and minimum rmsd with any other band
            if model_selection == 'min_drmsd':
                rmsdi = {}
                rmsdi_band = ''
                rmsdi_min = 5
                for i,b in enumerate(tmp[7].keys()):
                    if b == idx: continue

                    rmsd_yi = [rdark_selected[b], rdark_selected[idx]]
                    rmsd_xi = [tmp[0][b],tmp[0][idx]]
                    rmsdi[b] = pp.rmsd(rmsd_xi,rmsd_yi)
                    if is_masked(rmsdi[b]): rmsdi[b] = 5
                    if (rmsdi[b] < rmsdi_min) & (b != idx):
                        rmsdi_band = b
                        rmsdi_min = rmsdi[b]

                # this happens if no band fits better (e.g. for rhopath==rhorayleigh)
                if rmsdi_band == '': 
                    rmsdi_band = b
                    sel_rmsd = -1
                    sel_idx = (idx,idx)
                    sel_ac_par = tmp
                    sel_model_lut, sel_model_lut_meta = (lut_sensor, meta_sensor)
                    pixel_idx = -1
                elif (rmsdi[rmsdi_band] <= sel_rmsd) | (li == 0):
                    sel_rmsd = rmsdi[rmsdi_band]
                    sel_idx = (idx,rmsdi_band)
                    sel_ac_par = tmp
                    sel_model_lut, sel_model_lut_meta = (lut_sensor, meta_sensor)
                    pixel_idx = -1

    ## get rdark/tau+model according to min rmsd/tau in the given pixel lists
    ## when rdark_list_selection did not select a single spectrum
    if (rdark_list):
        ## choose best lut for each pixel
        tau_550_sel=[]
        tau_rmsd_sel=[]
        tau_band_sel=[]
        tau_model_sel=[]
        
        ## get for each pixel the minimum rmsd [1] according to the models
        ## default select on min tau
        list_min_idx = 0 
        
        if model_selection == "min_tau":
            list_min_idx = 0
        elif model_selection == "min_rmsd":
            list_min_idx = 1
        
        for i in pixel_range:
            val_c = [taufits[j][list_min_idx][i] for j in range(0,len(taufits))] ## selected tau/rmsd per model for this pixel
            idx = [i for i,t in enumerate(val_c) if (t == nanmin(val_c))]
            if len(idx) == 0: 
                tau_550_sel.append(nan)
                tau_rmsd_sel.append(nan)
                tau_band_sel.append(nan)
                tau_model_sel.append(nan)
            else:
                idx = idx[0]
                tau_550_sel.append(taufits[idx][0][i])
                tau_rmsd_sel.append(taufits[idx][1][i])
                tau_band_sel.append(taufits[idx][2][i])
                tau_model_sel.append(idx+1)
        
        ## select pixel with min rmsd
        pixel_idx = where(tau_rmsd_sel == nanmin(tau_rmsd_sel))[0][0]
    
        ## select lut and band
        lut_sel = luts[tau_model_sel[pixel_idx]-1]
        sel_idx = tau_band_sel[pixel_idx]
        sel_rmsd = tau_rmsd_sel[pixel_idx]
        
        ## construct the selected dark pixel
        for band in rdark_selected.keys(): rdark_selected[band] = rdark_selected[band][pixel_idx]
            
        ## get selected LUT
        if pressure is not None:
            sel_model_lut, sel_model_lut_meta = pp.aerlut.aerlut_pressure(lut_sel, lutdir, pressure, sensor, rsr_file, lut_data_dict=lut_data_dict)
        else:
            sel_model_lut, sel_model_lut_meta = (lut_data_dict[lut_sel]["lut"], lut_data_dict[lut_sel]["meta"])
        
        ## get ac parameters for the selected rdark  - probably better to use the TAU to retrieve parameters ?
        sel_ac_par = pp.aerlut.lut_get_ac_parameters_sensor(sel_model_lut,sel_model_lut_meta,azi,thv,ths,
                                                            rdark_selected, force_band=force_band)


    ## get ac parameters from previously selected model
    (ratm_s,rorayl_s,dtotr_s,utotr_s,dtott_s,utott_s,astot_s) = sel_ac_par[0:7]
    tau550_all_bands = sel_ac_par[7]
    if type(sel_idx) == str:
        tau550 =tau550_all_bands[sel_idx]
    else:
        tau550 =tau550_all_bands[sel_idx[0]]
    dark_idx = sel_idx

    ## if tau is the minimum in the model ratm may be nans
    ## in that case replace by Rayleigh reflectance
    from numpy import isnan
    for band in ratm_s.keys():
        if isnan(ratm_s[band]): ratm_s[band] = rorayl_s[band]

    return (ratm_s,rorayl_s,dtotr_s,utotr_s,dtott_s,utott_s,astot_s, tau550),\
           (bands_sorted, tau550_all_bands, dark_idx, sel_rmsd, rdark_selected, pixel_idx), (sel_model_lut, sel_model_lut_meta)

