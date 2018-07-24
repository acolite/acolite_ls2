## def acolite_l2w
## computes l2w parameters from l2r ACOLITE output file
##
## Written by Quinten Vanhellemont 2017-12-05
## Last modifications: 2018-03-07 (QV) added l2w masking and several parameters:
##                                         t_nechad, t_dogliotti and chl_oc
##                     2018-03-08 (QV) added more parameters: chl_re, fai, ndvi, qaa
##                                     and improved band selection
##                                     added bt and parameter aliases
##                     2018-03-22 (QV) fixed nechad type algorithms
##                     2018-03-28 (QV) added rhow and Rrs
##                     2018-04-17 (QV) tracking dataset attributes from inputfile, added hue angle
##                     2018-04-18 (QV) added hue angle coeff file, updated chl_oc wavelength requirement
##                     2018-06-06 (QV) added xy outputs
##                     2018-06-07 (QV) added l2_flags output and l2w water masking option, added OLH, changed S2A wavelengths for T
##                     2018-07-18 (QV) changed acolite import name
##                     2018-07-24 (QV) fixed problem with shortest blue band after chl_oc3 computation
##                                     added l2w_mask_negative_rhow keyword

def acolite_l2w(inputfile, output, parameters=None, output_map=False, retain_data_read=False,
                l2w_mask=True, l2w_mask_wave=1600, l2w_mask_threshold=0.0215, l2w_mask_water_parameters=True, l2w_mask_negative_rhow=True,
                nc_compression=True, chunking=True):
    import os
    import datetime, time

    from acolite.shared import datascl,nc_data,nc_datasets,nc_gatts,closest_idx,\
                               coef_hue_angle,coef_nechad_hs, coef_chl_re_gons, coef_nechad2016, coef_nechad_spm_hs, coef_chl_oc
    from acolite.acolite import l2w_required, acolite_l2w_qaa
    from acolite.output import nc_write

    from numpy import pi, nan, where, log10, isfinite, power

    if not os.path.exists(inputfile):
        print('File {} not found.'.format(inputfile))
        return(False)

    ## read attributes from L2R inputfile
    gatts = nc_gatts(inputfile)
    l2r_datasets = nc_datasets(inputfile)
    print(l2r_datasets)

    ## get wavelengths in datafile
    rhos_waves = [ds.split('_')[1] for ds in l2r_datasets if 'rhos_' in ds]
    rhot_waves = [ds.split('_')[1] for ds in l2r_datasets if 'rhot_' in ds]

    if ("l2r_file" not in gatts) or (len(rhos_waves)==0):
        print('File {} is probably not ACOLITE L2R file.'.format(inputfile))
                
    if parameters is not None:

        ## make mask dataset
        if l2w_mask:
            mask_wave = 0
            mask_dataset = ''
            mind = l2w_mask_wave*1.0

            for ds in l2r_datasets:
                if 'rhot' not in ds: continue
                dst, dsw = ds.split('_')
                if abs(float(dsw)-float(l2w_mask_wave)) < mind:
                    mind = abs(float(dsw)-float(l2w_mask_wave))
                    mask_wave = dsw
                    mask_dataset = ds
            if mask_dataset in l2r_datasets:
                 mask = nc_data(inputfile,mask_dataset) > float(l2w_mask_threshold)
            else:
                 mask = None

            ## make l2_flags dataset
            l2_flags = None
            if ("l2_negatives" in l2r_datasets) & (l2w_mask_negative_rhow):
                l2_flags = nc_data(inputfile, "l2_negatives")
            if l2_flags is None:
                l2_flags = mask.astype(int)*(2**0)
            else:
                l2_flags += mask.astype(int)*(2**0)
            mask = l2_flags != 0

        ## make outputfile
        if output is not None:
            gatts['output_dir'] = output
            
        ncfile = '{}/{}_L2W.nc'.format(gatts['output_dir'],gatts['output_name'])
        nc_new = True
        gatts["file_type"] =  'Level 2 Water Product'
        gatts["l2w_file"] =  ncfile

        ## make sure parameters is a list
        if type(parameters) is not list: parameters=[parameters]

        aliases = ['qaa_all']
        added_par, removed_par = [],[]
        for par in parameters:
            if par[-1]=='*':
                print(par[0:-1])
                if par[0:-2] in ['Rrs','rhow']:
                    added_par += ['{}_{}'.format(par[0:-2], ds.split('_')[1]) for ds in l2r_datasets if 'rhos_' in ds]
                else:
                    added_par += [ds for ds in l2r_datasets if par[0:-1] in ds]
                removed_par.append(par)
            if par == 'qaa_all':
                added_par+=['a443_qaasw','bbp443_qaasw','kd443_qaasw',
                            'a490_qaasw','bbp490_qaasw','kd490_qaasw',
                            'a560_qaasw','bbp560_qaasw','kd560_qaasw',
                            'a665_qaasw','bbp665_qaasw','kd665_qaasw',
                            'kdpar_qaasw']
                removed_par.append(par)
        if len(added_par)>0: 
            parameters+=added_par
        if len(removed_par)>0: 
            for par in removed_par: parameters.remove(par)

    
        ## track whether we ran qaa already
        qaa_computed = False
        qaa_data = {}
        
        ## store read data
        data_read = {}
        att_read = {}

        for par in parameters:
            ## clearing the read data from memory might be needed for large datasets 
            if not retain_data_read: 
                data_read = {}
                att_read = {}

            ## reset variables
            mask_data = True
            
            req = None
            par_computed = False
            par_exists = False
            par_data = None
            par_attributes = None
            par_split = None
            required_datasets = []

            ## fold parameter case and split
            par_name = par.casefold().strip()
            par_split = par_name.split('_')

            #################################
            ## Nechad SPM
            if par_name[0:10] == 'spm_nechad':
                par_exists = True
                mask_data = True
                par_split = par_name.split('_')
                
                if par_split[1] == 'nechad':
                    par_attributes = {'algorithm':'Nechad et al. 2010', 'title':'Nechad SPM'}
                    ### get required datasets
                    if len(par_split) > 2:
                        required_datasets = ['rhos_{}'.format(par_split[2])]
                        par_attributes['wave']=par_split[2]
                    else: ## defaults
                        if gatts['sensor'] == 'L5_TM':
                            required_datasets = ['rhos_660']
                        elif gatts['sensor'] == 'L7_ETM':
                            required_datasets = ['rhos_661']
                        elif gatts['sensor'] == 'L8_OLI':
                            required_datasets = ['rhos_655']
                        elif gatts['sensor'] == 'S2A_MSI':
                            required_datasets = ['rhos_665']
                        elif gatts['sensor'] == 'S2B_MSI':
                            required_datasets = ['rhos_665']
                        else:
                            print('Parameter {} not configured for {}.'.format(par_name,gatts['sensor']))
                            continue

                        par_attributes['wave']=required_datasets[0].split('_')[1]

                    ## test if required datasets are present
                    if len(required_datasets) > 0: 
                        req = l2w_required(inputfile, required_datasets, data_read, att_read)
                        ## compute dataset
                        if req:
                            cwave = float(required_datasets[0].split('_')[1])/1000.
                            nechad = coef_nechad_spm_hs()
                            idx,xret = min(enumerate(nechad['wave']), key=lambda x: abs(x[1]-cwave))
                            par_attributes['wave']=nechad['wave'][idx]
                            par_attributes['A_SPM']=nechad['A'][idx]
                            par_attributes['C_SPM']=nechad['C'][idx]
                            par_data = (par_attributes['A_SPM'] * data_read[required_datasets[0]]) /\
                                       (1.-data_read[required_datasets[0]]/par_attributes['C_SPM'])

                if par_split[1] == 'nechad2016':
                    par_attributes = {'algorithm':'Nechad et al. 2010, 2016 calibration', 'title':'Nechad SPM'}
                    npar = 'SPM'
                    nsen = gatts['sensor'].split('_')[1]
                    ### get required datasets
                    if len(par_split) > 2:
                        required_datasets = ['rhos_{}'.format(par_split[2])]
                        par_attributes['wave']=par_split[2]
                    else: ## defaults
                        if gatts['sensor'] == 'L8_OLI':
                            required_datasets = ['rhos_655']
                        elif gatts['sensor'] == 'S2A_MSI':
                            required_datasets = ['rhos_665']
                        elif gatts['sensor'] == 'S2B_MSI':
                            required_datasets = ['rhos_665']
                        else:
                            print('Parameter {} not configured for {}.'.format(par_name,gatts['sensor']))
                            continue
                        par_attributes['wave']=required_datasets[0].split('_')[1]

                    ## test if required datasets are present
                    if len(required_datasets) > 0: 
                        req = l2w_required(inputfile, required_datasets, data_read, att_read)
                        ## compute dataset
                        if req:
                            nechad = coef_nechad2016()
                            for n in nechad:
                                if n['par'] != npar: continue
                                if n['sensor'] != nsen: continue
                                if abs(n['wave'] - float(par_attributes['wave'])) > 10: continue
                                par_attributes['wave']=n['wave']
                                par_attributes['A_SPM']=n['A']
                                par_attributes['C_SPM']=n['C']

                            if ('A_SPM' in par_attributes) & ('C_SPM' in par_attributes):
                                par_data = (par_attributes['A_SPM'] * data_read[required_datasets[0]]) /\
                                           (1.-data_read[required_datasets[0]]/par_attributes['C_SPM'])
            ##### END Nechad SPM
            #################################
            
            #################################
            ## Nechad turbidity
            if par_name[0:8] == 't_nechad':
                par_exists = True
                mask_data = True
                par_split = par_name.split('_')
                par_attributes = {'algorithm':'Nechad et al. 2009'}
                
                if par_split[1] == 'nechad':
                    if len(par_split) > 2:
                        if par_split[2].lower() in ['red', 'nir']:
                            if par_split[2].lower() == 'red':
                                twave = 660
                            if par_split[2].lower() == 'nir':
                                twave = 845

                            ### get required datasets
                            if gatts['sensor'] == 'L5_TM':
                                if twave == 660:
                                    A_TUR = 249.61
                                    C_TUR = 0.1701
                                if twave == 845:
                                    A_TUR = 1382.10
                                    C_TUR = 0.1782
                            elif gatts['sensor'] == 'L7_ETM':
                                if twave == 660:
                                    A_TUR = 251.58
                                    C_TUR = 0.1710
                                if twave == 845:
                                    A_TUR = 1467.57
                                    C_TUR = 0.1891
                            elif gatts['sensor'] == 'L8_OLI':
                                if twave == 660:
                                    A_TUR = 242.27
                                    C_TUR = 0.1682
                                if twave == 845:
                                    A_TUR = 2106.14
                                    C_TUR = 0.2113
                            elif gatts['sensor'] == 'S2A_MSI':
                                if twave == 660:
                                    A_TUR = 268.52
                                    C_TUR = 0.1725
                                if twave == 845:
                                    A_TUR = 1441.19
                                    C_TUR = 0.1928
                            elif gatts['sensor'] == 'S2B_MSI':
                                if twave == 660:
                                    A_TUR = 270.17
                                    C_TUR = 0.1726
                                if twave == 845:
                                    A_TUR = 1450.73
                                    C_TUR = 0.1933

                            widx,selwave = closest_idx(rhos_waves, twave)
                            required_datasets = ['rhos_{}'.format(selwave)]
                            algwave = 'ave'
                            print(algwave, A_TUR, C_TUR)
                        else:
                            widx,selwave = closest_idx(rhos_waves, float(par_split[2]))
                            required_datasets = ['rhos_{}'.format(selwave)]
                            par_attributes['wave']=float(selwave)
                            ndt = coef_nechad_hs('T')
                            didx,algwave = closest_idx(ndt['wave'], selwave)
                            A_TUR = ndt['A'][didx]
                            C_TUR = ndt['C'][didx]
                            print(algwave, A_TUR, C_TUR)
                    else: ## defaults
                        A_TUR = 228.1
                        C_TUR = 0.1641
                        algwave = 645
                        ### get required datasets
                        if gatts['sensor'] == 'L5_TM':
                             required_datasets = ['rhos_660']
                        elif gatts['sensor'] == 'L7_ETM':
                             required_datasets = ['rhos_661']
                        elif gatts['sensor'] == 'L8_OLI':
                             required_datasets = ['rhos_655']
                        elif gatts['sensor'] == 'S2A_MSI':
                             required_datasets = ['rhos_665']
                        elif gatts['sensor'] == 'S2B_MSI':
                             required_datasets = ['rhos_665']
                        else:
                            print('Parameter {} not configured for {}.'.format(par_name,gatts['sensor']))
                            continue

                    ## test if required datasets are present
                    if len(required_datasets) > 0: 
                        req = l2w_required(inputfile, required_datasets, data_read, att_read)
                        ## compute dataset
                        if req:
                            par_data = None
                            par_attributes['algo_wave']=algwave
                            par_attributes['wave']=required_datasets[0].split('_')[1]
                            par_attributes['A_TUR']=A_TUR
                            par_attributes['C_TUR']=C_TUR
                            ## compute turbidity
                            par_data=(par_attributes['A_TUR'] * data_read[required_datasets[0]]) \
                                   / (1.-data_read[required_datasets[0]]/par_attributes['C_TUR'])
                                
                if par_split[1] == 'nechad2016':
                    par_attributes = {'algorithm':'Nechad et al. 2009, 2016 calibration', 'title':'Nechad TUR'}
                    npar = 'TUR'
                    nsen = gatts['sensor'].split('_')[1]
                    ### get required datasets
                    if len(par_split) > 2:
                        required_datasets = ['rhos_{}'.format(par_split[2])]
                        par_attributes['wave']=par_split[2]
                    else: ## defaults
                        if gatts['sensor'] == 'L8_OLI':
                            required_datasets = ['rhos_655']
                        elif gatts['sensor'] == 'S2A_MSI':
                            required_datasets = ['rhos_665']
                        elif gatts['sensor'] == 'S2B_MSI':
                            required_datasets = ['rhos_665']
                        else:
                            print('Parameter {} not configured for {}.'.format(par_name,gatts['sensor']))
                            continue
                        par_attributes['wave']=required_datasets[0].split('_')[1]

                    ## test if required datasets are present
                    if len(required_datasets) > 0: 
                        req = l2w_required(inputfile, required_datasets, data_read, att_read)
                        ## compute dataset
                        if req:
                            nechad = coef_nechad2016()
                            for n in nechad:
                                if n['par'] != npar: continue
                                if n['sensor'] != nsen: continue
                                if abs(n['wave'] - float(par_attributes['wave'])) > 10: continue
                                par_attributes['wave']=n['wave']
                                par_attributes['A_TUR']=n['A']
                                par_attributes['C_TUR']=n['C']
                                
                            par_data = (par_attributes['A_TUR'] * data_read[required_datasets[0]]) /\
                                       (1.-data_read[required_datasets[0]]/par_attributes['C_TUR'])
            ##### END Nechad T
            #################################

            #################################
            ## Dogliotti turbidity
            if par_name[0:11] == 't_dogliotti':
                par_exists = True
                mask_data = True
                par_split = par_name.split('_')
                par_attributes = {'algorithm':'Dogliotti et al. 2015'}
                
                ### get required datasets
                if gatts['sensor'] == 'L5_TM':
                     required_datasets = ['rhos_660','rhos_839']
                elif gatts['sensor'] == 'L7_ETM':
                     required_datasets = ['rhos_661', 'rhos_835']
                elif gatts['sensor'] == 'L8_OLI':
                     required_datasets = ['rhos_655','rhos_865']
                elif gatts['sensor'] == 'S2A_MSI':
                     required_datasets = ['rhos_665','rhos_833']
                elif gatts['sensor'] == 'S2B_MSI':
                     required_datasets = ['rhos_665','rhos_833']
                else:
                     print('Parameter {} not configured for {}.'.format(par_name,gatts['sensor']))
                     continue

                ## test if required datasets are present
                if len(required_datasets) > 0: 
                    req = l2w_required(inputfile, required_datasets, data_read, att_read)
                    ## compute dataset
                    if req:
                        par_data = None
                        par_attributes['algo_wave_red']=645
                        par_attributes['algo_wave_nir']=859
                        par_attributes['wave_red']=required_datasets[0].split('_')[1]
                        par_attributes['wave_nir']=required_datasets[1].split('_')[1]
                        par_attributes['A_T_red']=228.1
                        par_attributes['C_T_red']=0.1641
                        par_attributes['A_T_nir']=3078.9
                        par_attributes['C_T_nir']=0.2112
                        par_attributes['lower_lim']=0.05
                        par_attributes['upper_lim']=0.07

                        par_data=(par_attributes['A_T_red'] * data_read[required_datasets[0]]) \
                               / (1.-data_read[required_datasets[0]]/par_attributes['C_T_red'])

                        par_data_nir=(par_attributes['A_T_nir'] * data_read[required_datasets[1]]) \
                               / (1.-data_read[required_datasets[1]]/par_attributes['C_T_nir'])

                        ## blending algorithm
                        if len(par_split) == 2:
                            ## replace most turbid with nir band
                            sub = where(data_read[required_datasets[0]] >= par_attributes['upper_lim'])
                            if len(sub[0]) > 0: par_data[sub] = par_data_nir[sub]

                            ## blend in between
                            sub = where((data_read[required_datasets[0]] < par_attributes['upper_lim']) & \
                                        (data_read[required_datasets[0]] >= par_attributes['lower_lim']))
                            if len(sub[0]) > 0: 
                                w=(data_read[required_datasets[0]]  - par_attributes['lower_lim']) \
                                         / (par_attributes['upper_lim']-par_attributes['lower_lim'])
                                par_data[sub] = ((1.-w[sub]) * par_data[sub]) + (w[sub]*par_data_nir[sub])
                        
                        ## individual algorithms
                        if len(par_split) == 3:
                            if par_split[2].lower() == 'red':
                                par_data = par_data*1.0
                            elif par_split[2].lower() == 'nir':
                                par_data = par_data_nir*1.0
                            else:
                                par_data = None
                        par_data_nir = None
            ##### END Dogliotti T
            #################################
            
            #################################
            ## CHL_OC
            if par_name[0:6] == 'chl_oc':
                chl_oc_wl_diff = 20
                par_exists = True
                mask_data = True
                par_split = par_name.split('_')
                par_attributes = {'algorithm':'Chlorophyll a blue/green ratio', 'dataset':'rhos'}
                ds_waves = [w for w in rhos_waves]

                chl_oc = coef_chl_oc()
                
                if par_name in chl_oc:
                    blue_wave = chl_oc[par_name]['blue']
                    green_wave = chl_oc[par_name]['green']
                    chl_coef = chl_oc[par_name]['chl_coef']
                    chl_ref = chl_oc[par_name]['reference']
                    blue_wave_sel = []
                    green_wave_sel = []
                    
                    for i, reqw in enumerate(blue_wave+green_wave):
                        widx,selwave = closest_idx(ds_waves, reqw)
                        if abs(float(selwave)-float(reqw)) > chl_oc_wl_diff: continue
                        selds='{}_{}'.format(par_attributes['dataset'],selwave)

                        if reqw in blue_wave: 
                            blue_wave_sel.append(selwave)
                            required_datasets.append(selds)
                        if reqw in green_wave: 
                            green_wave_sel.append(selwave)
                            required_datasets.append(selds)

                    ## store selections
                    par_attributes['blue']=blue_wave_sel
                    par_attributes['green']=green_wave_sel
                    par_attributes['reference']=chl_ref
                    par_attributes['coefficients']=chl_coef

                if (len(blue_wave) != len(blue_wave_sel)) or (len(green_wave) != len(green_wave_sel)):
                    print('Error selecting blue and green bands.')
                    par_data = None
                else:
                    ## test if required datasets are present
                    if len(required_datasets) > 0: 
                        req = l2w_required(inputfile, required_datasets, data_read, att_read)
                        
                        ## max blue
                        for i,wv in enumerate(blue_wave_sel):
                            cur_tag = 'rhos_{}'.format(wv)
                            if i == 0:
                                blue_rhow = data_read[cur_tag]*1.0
                            else:
                                sub = where((data_read[cur_tag]>blue_rhow) & isfinite(data_read[cur_tag]))
                                if len(sub[0]) > 0:
                                    blue_rhow[sub]=data_read[cur_tag][sub]
                                
                        ## max green
                        for i,wv in enumerate(green_wave_sel):
                            cur_tag = 'rhos_{}'.format(wv)
                            if i == 0:
                                green_rhow = data_read[cur_tag]*1.0
                            else:
                                sub = where((data_read[cur_tag]>green_rhow) & isfinite(data_read[cur_tag]))
                                if len(sub[0]) > 0:
                                    green_rhow[sub]=data_read[cur_tag][sub]
                                
                        ## compute dataset
                        if req:
                            par_data = None
                            ratio = log10(blue_rhow/green_rhow)
                            blue_rhow = None
                            green_rhow = None
                            par_data = chl_coef[0] + chl_coef[1] * ratio + \
                                                     chl_coef[2] * ratio * ratio + \
                                                     chl_coef[3] * ratio * ratio * ratio + \
                                                     chl_coef[4] * ratio * ratio * ratio * ratio
                            par_data = power(10, par_data)
            ##### END CHL_OC
            #################################
            
            #################################
            ## CHL RE
            if par_name[0:6] == 'chl_re':
                par_exists = True
                mask_data = True
                par_split = par_name.split('_')
                par_attributes = {'algorithm':'Red-edge chlorophyll', 'dataset':'rhos'}
                req_waves,req_waves_selected = [],[]
                ds_waves = [w for w in rhos_waves]

                if len(par_split) >= 3:
                    if par_split[2][0:4] == 'gons': 
                        par_attributes['algorithm']='Gons et al. 3 band'
                        gons = coef_chl_re_gons()

                        ### get required datasets
                        if gatts['sensor'] in ['S2A_MSI', 'S2B_MSI']:
                            if par_split[2] == 'gons':
                                req_waves = [670,705,780]
                                gons_name = 'chl_re_gons'
                                req_waves = [gons[gons_name][tag] for tag in ['red_band', 'rededge_band', 'nir_band']]
                            if par_split[2] == 'gons740':
                                req_waves = [670,705,740]    
                                gons_name = 'chl_re_gons740'
                                req_waves = [gons[gons_name][tag] for tag in ['red_band', 'rededge_band', 'nir_band']]
                        else:
                            print('Parameter {} not configured for {}.'.format(par_name,gatts['sensor']))
                            continue

                        if len(req_waves) > 0:
                            for i, reqw in enumerate(req_waves):
                                widx,selwave = closest_idx(ds_waves, reqw)
                                if abs(float(selwave)-float(reqw)) > 10: continue
                                selds='{}_{}'.format(par_attributes['dataset'],selwave)
                                required_datasets.append(selds)
                                req_waves_selected.append(selwave)
                            par_attributes['waves']=req_waves_selected

                        ## test if required datasets are present
                        if len(required_datasets) > 0: 
                            req = l2w_required(inputfile, required_datasets, data_read, att_read)
                            ## compute dataset
                            ## put coefficients in external file
                            if req:
                                #bb = (1.61 * data_read[required_datasets[2]]) \
                                #      / (0.082 - 0.6 * data_read[required_datasets[2]])
                                #rm = data_read[required_datasets[1]]/data_read[required_datasets[0]]
                                #par_data = ((rm * (0.7 + bb)) - 0.40 - power(bb,1.05))
                                #par_data /= 0.015

                                gc = gons[gons_name]['chl_coef']
                                bb = (gc[0] * data_read[required_datasets[2]]) \
                                      / (gc[1] - gc[2] * data_read[required_datasets[2]])
                                rm = data_read[required_datasets[1]]/data_read[required_datasets[0]]
                                par_data = ((rm * (gc[3] + bb)) - gc[4] - power(bb,gc[5]))
                                par_data /= gons[gons_name]['astar_chl']

                                par_data[par_data<0]=nan

                                gm = gons[gons_name]['validity']
                                par_data[((data_read[required_datasets[0]] <= gm[0]) | (data_read[required_datasets[1]]/data_read[required_datasets[0]] <= gm[1]))]=nan
                                
                    if par_split[2][0:5] == 'moses':
                        par_attributes['algorithm']='Moses et al. 3 band'
                        par_attributes['a']=(232.29,23.173) ## put coefficients in external file
                                                    
                        ### get required datasets
                        if gatts['sensor'] in ['S2A_MSI', 'S2B_MSI']:
                            if par_split[2] == 'moses3b':
                                req_waves = [670,705,780]
                            if par_split[2] == 'moses3b740':
                                req_waves = [670,705,740]    
                        else:
                            print('Parameter {} not configured for {}.'.format(par_name,gatts['sensor']))
                            continue

                        if len(req_waves) > 0:
                            for i, reqw in enumerate(req_waves):
                                widx,selwave = closest_idx(ds_waves, reqw)
                                if abs(float(selwave)-float(reqw)) > 10: continue
                                selds='{}_{}'.format(par_attributes['dataset'],selwave)
                                required_datasets.append(selds)
                                req_waves_selected.append(selwave)
                            par_attributes['waves']=req_waves_selected

                        ## test if required datasets are present
                        if len(required_datasets) > 0: 
                            req = l2w_required(inputfile, required_datasets, data_read, att_read)
                            ## compute dataset
                            if req:
                                par_data = par_attributes['a'][0]* \
                                        ((power(data_read[required_datasets[0]],-1)-power(data_read[required_datasets[1]],-1)) * \
                                        data_read[required_datasets[2]]) + par_attributes['a'][1]
                                par_data[par_data<0]=nan  
                                
                    if par_split[2][0:6] == 'mishra': 
                        par_attributes['algorithm']='Mishra et al. 2014, NDCI'
                        par_attributes['a']=(14.039, 86.115, 194.325) ## put coefficients in external file
                        
                        ### get required datasets
                        if gatts['sensor'] in ['S2A_MSI', 'S2B_MSI']:
                            req_waves = [670,705]
                        else:
                            print('Parameter {} not configured for {}.'.format(par_name,gatts['sensor']))
                            continue

                        if len(req_waves) > 0:
                            for i, reqw in enumerate(req_waves):
                                widx,selwave = closest_idx(ds_waves, reqw)
                                if abs(float(selwave)-float(reqw)) > 10: continue
                                selds='{}_{}'.format(par_attributes['dataset'],selwave)
                                required_datasets.append(selds)
                                req_waves_selected.append(selwave)
                            par_attributes['waves']=req_waves_selected

                        ## test if required datasets are present
                        if len(required_datasets) > 0: 
                            req = l2w_required(inputfile, required_datasets, data_read, att_read)
                            ## compute dataset
                            if req:
                                ndci = (data_read[required_datasets[1]]-data_read[required_datasets[0]]) / \
                                       (data_read[required_datasets[1]]+data_read[required_datasets[0]])
                                par_data = par_attributes['a'][0] + par_attributes['a'][1]*ndci + par_attributes['a'][2]*ndci*ndci
                                ndci = None
                                par_data[par_data<0]=nan  
            ##### END CHL_RE
            #################################
            
            #################################
            ## NDCI
            if par_name[0:4] == 'ndci':
                par_exists = True
                par_split = par_name.split('_')
                par_attributes = {'algorithm':'Mishra et al. 2014, NDCI', 'dataset':'rhos'}
                req_waves,req_waves_selected = [],[]
                ds_waves = [w for w in rhos_waves]
                
                ### get required datasets
                if gatts['sensor'] in ['S2A_MSI', 'S2B_MSI']:
                    req_waves = [670,705]
                else:
                    print('Parameter {} not configured for {}.'.format(par_name,gatts['sensor']))
                    continue

                if len(req_waves) > 0:
                    for i, reqw in enumerate(req_waves):
                        widx,selwave = closest_idx(ds_waves, reqw)
                        if abs(float(selwave)-float(reqw)) > 10: continue
                        selds='{}_{}'.format(par_attributes['dataset'],selwave)
                        required_datasets.append(selds)
                        req_waves_selected.append(selwave)
                    par_attributes['waves']=req_waves_selected
                    
                ## test if required datasets are present
                if len(required_datasets) > 0: 
                    req = l2w_required(inputfile, required_datasets, data_read, att_read)
                    ## compute dataset
                    if req:
                        par_data = (data_read[required_datasets[1]]-data_read[required_datasets[0]]) / \
                                   (data_read[required_datasets[1]]+data_read[required_datasets[0]])
            ##### END NDCI
            #################################

            #################################
            ## SLH
            if par_name[0:3] == 'slh':
                par_exists = True
                par_split = par_name.split('_')
                par_attributes = {'algorithm':'Kudela et al. 2015, SLH', 'dataset':'rhos'}
                req_waves,req_waves_selected = [],[]
                ds_waves = [w for w in rhos_waves] 

                ### get required datasets
                if gatts['sensor'] in ['S2A_MSI', 'S2B_MSI']:
                    req_waves = [670,705,780]
                else:
                    print('Parameter {} not configured for {}.'.format(par_name,gatts['sensor']))
                    continue

                if len(req_waves) > 0:
                    for i, reqw in enumerate(req_waves):
                        widx,selwave = closest_idx(ds_waves, reqw)
                        if abs(float(selwave)-float(reqw)) > 10: continue
                        selds='{}_{}'.format(par_attributes['dataset'],selwave)
                        required_datasets.append(selds)
                        req_waves_selected.append(selwave)
                    par_attributes['waves']=req_waves_selected
                    
                ## test if required datasets are present
                if len(required_datasets) > 0: 
                    req = l2w_required(inputfile, required_datasets, data_read, att_read)
                    ## compute dataset
                    if req:
                        slh_waves = [float(ds.split('_')[1]) for ds in required_datasets]
                        ratio = (data_read[required_datasets[2]]-data_read[required_datasets[0]]) / \
                                (slh_waves[2]+slh_waves[0])
                                
                        par_data = data_read[required_datasets[1]] - \
                                    (data_read[required_datasets[0]] + (ratio)*(slh_waves[1]+slh_waves[0]))
                        ratio = None
            ##### END SLH
            #################################

            
            #################################
            ## QAA
            if 'qaa' in par_name:
                par_data = None
                par_exists = True
                par_split = par_name.split('_')
                par_attributes = {'algorithm':'QAA', 'dataset':'rhos'}

                if not qaa_computed:
                    req_waves,req_waves_selected = [],[]
                    ds_waves = [w for w in rhos_waves] 

                    ### get required datasets
                    if gatts['sensor'] in ['L8_OLI','S2A_MSI', 'S2B_MSI']:
                        req_waves = [443, 490, 560, 665]
                    else:
                        print('Parameter {} not configured for {}.'.format(par_name,gatts['sensor']))
                        continue

                    for i, reqw in enumerate(req_waves):
                        widx,selwave = closest_idx(ds_waves, reqw)
                        if abs(float(selwave)-float(reqw)) > 10: continue
                        selds='{}_{}'.format(par_attributes['dataset'],selwave)
                        required_datasets.append(selds)
                        req_waves_selected.append(selwave)
                    par_attributes['waves']=req_waves_selected
                    
                    if len(required_datasets) > 0: 
                        req = l2w_required(inputfile, required_datasets, data_read, att_read)
                        #return(data_read, par_attributes)
                        if req:
                            qaa_data = acolite_l2w_qaa(data_read,par_attributes,satellite=gatts['sensor'],ths=gatts['THS'])
                            qaa_computed = True
            
                #print(qaa_data.keys())
                ## get parameter from qaa outputs
                for k in qaa_data.keys():
                    if k.lower() == par_name: 
                        par_data = qaa_data[k]
            ##### END QAA
            #################################

            #################################
            ## FAI
            if par_name[0:3] == 'fai':
                par_exists = True
                mask_data = False
                par_split = par_name.split('_')
                par_attributes = {'algorithm':'Floating Algal Index, Hu et al. 2009', 'dataset':'rhos'}
                req_waves,req_waves_selected = [],[]
                ds_waves = [w for w in rhos_waves] 

                if par_name=='fai_rhot': 
                    par_attributes['dataset']='rhot'
                    ds_waves = [w for w in rhot_waves]
                    
                fai_diff = [10, 30, 80]
                req_waves = [660,865,1610]                
                for i, reqw in enumerate(req_waves):
                    widx,selwave = closest_idx(ds_waves, reqw)
                    if abs(float(selwave)-float(reqw)) > fai_diff[i]: continue
                    selds='{}_{}'.format(par_attributes['dataset'],selwave)
                    required_datasets.append(selds)
                    req_waves_selected.append(selwave)
                par_attributes['waves']=req_waves_selected
                print(req_waves_selected)
                print(required_datasets)

                if len(required_datasets) == len(req_waves): 
                    req = l2w_required(inputfile, required_datasets, data_read, att_read)                            
                    ## compute dataset
                    if req:
                        par_data = None
                        fai_sc = (float(par_attributes['waves'][1])-float(par_attributes['waves'][0]))/\
                                 (float(par_attributes['waves'][2])-float(par_attributes['waves'][0]))
                        nir_prime = data_read[required_datasets[0]] + \
                                   (data_read[required_datasets[2]]-data_read[required_datasets[0]]) * fai_sc
                        par_data = data_read[required_datasets[1]] - nir_prime
                        nir_prime = None
            ##### END FAI
            #################################

            #################################
            ## NDVI
            if par_name[0:4] == 'ndvi':
                par_exists = True
                mask_data = False
                par_split = par_name.split('_')
                par_attributes = {'algorithm':'NDVI', 'dataset':'rhos'}
                
                req_waves,req_waves_selected = [],[]
                ds_waves = [w for w in rhos_waves] 

                if par_name=='ndvi_rhot': 
                    par_attributes['dataset']='rhot'
                    ds_waves = [w for w in rhot_waves]

                ### get required datasets
                if gatts['sensor'] in ['L5_TM','L7_ETM','S2A_MSI', 'S2B_MSI']:
                    req_waves = [660,840]
                else:
                    req_waves = [660,865]
                req_waves_selected = []
                
                for i, reqw in enumerate(req_waves):
                    widx,selwave = closest_idx(ds_waves, reqw)
                    if abs(float(selwave)-float(reqw)) > 20: continue
                    selds='{}_{}'.format(par_attributes['dataset'],selwave)
                    required_datasets.append(selds)
                    req_waves_selected.append(selwave)
                par_attributes['waves']=req_waves_selected

                if len(required_datasets) > 0: 
                    req = l2w_required(inputfile, required_datasets, data_read, att_read)                     
                    ## compute dataset
                    if req:
                        par_data = None
                        par_data = (data_read[required_datasets[1]]-data_read[required_datasets[0]])/\
                                   (data_read[required_datasets[1]]+data_read[required_datasets[0]])
            ##### END NDVI
            #################################
            
            #################################
            ## rhow
            if par_name[0:4] == 'rhow':
                par_exists = True
                mask_data = True
                par_split = par_name.split('_')
                par_attributes = {'algorithm':'Water reflectance', 'dataset':'rhos'}
                required_datasets = ['rhos_{}'.format(par_split[1])]
                if len(required_datasets) > 0: 
                        req = l2w_required(inputfile, required_datasets, data_read, att_read)
                        ## compute dataset
                        if req:
                            par_data = None
                            par_data = data_read[required_datasets[0]] * 1.0
                            for tag in att_read[required_datasets[0]]:
                                par_attributes[tag] = att_read[required_datasets[0]][tag]
            ##### END rhow
            #################################

            #################################
            ## Rrs
            if par_name[0:3] == 'rrs':
                par_exists = True
                mask_data = True
                par_split = par_name.split('_')
                par_name = 'Rrs_{}'.format(par_split[1])
                par_attributes = {'algorithm':'Remote sensing reflectance', 'dataset':'rhos', 'parname':par_name}
                required_datasets = ['rhos_{}'.format(par_split[1])]
                if len(required_datasets) > 0: 
                        req = l2w_required(inputfile, required_datasets, data_read, att_read)
                        ## compute dataset
                        if req:
                            par_data = None
                            par_data = data_read[required_datasets[0]]/pi
                            for tag in att_read[required_datasets[0]]:
                                par_attributes[tag] = att_read[required_datasets[0]][tag]
            ##### END Rrs
            #################################

            #################################
            ## Hue Angle
            if par_name == 'hue_angle':
                par_data = None
                par_exists = True
                mask_data = True
                par_split = par_name.split('_')
                par_attributes = {'algorithm':'Hue Angle', 'dataset':'rhos'}

                if not qaa_computed:
                    req_waves,req_waves_selected = [],[]
                    ds_waves = [w for w in rhos_waves] 
                    
                    hue_coeff = coef_hue_angle()

                    ### get required datasets
                    if gatts['sensor'] in hue_coeff:
                        req_waves = hue_coeff[gatts['sensor']]['req_waves']
                        hac = hue_coeff[gatts['sensor']]
                    else:
                        print('Parameter {} not configured for {}.'.format(par_name,gatts['sensor']))
                        continue

                    for i, reqw in enumerate(req_waves):
                        widx,selwave = closest_idx(ds_waves, reqw)
                        if abs(float(selwave)-float(reqw)) > 10: continue
                        selds='{}_{}'.format(par_attributes['dataset'],selwave)
                        required_datasets.append(selds)
                        req_waves_selected.append(selwave)
                    par_attributes['waves']=req_waves_selected
                    
                    if len(required_datasets) > 0: 
                        req = l2w_required(inputfile, required_datasets, data_read, att_read)
                        if req:
                            from numpy import where, nan, mod, pi, arctan2, nanmedian, power
                            yw = 1/3.
                            xw = 1/3.
                            print(req_waves_selected)
                            for iw, w in enumerate(req_waves_selected):
                                idx, w_ = closest_idx(hac['lambda'], req_waves_selected[iw])
                                print(idx, w, w_)
                                if iw == 0:
                                    X = data_read[required_datasets[iw]] * hac['X'][idx]
                                    Y = data_read[required_datasets[iw]] * hac['Y'][idx]
                                    Z = data_read[required_datasets[iw]] * hac['Z'][idx]
                                else:
                                    X += data_read[required_datasets[iw]] * hac['X'][idx]
                                    Y += data_read[required_datasets[iw]] * hac['Y'][idx]
                                    Z += data_read[required_datasets[iw]] * hac['Z'][idx]
                            X[where(mask)] = nan
                            Y[where(mask)] = nan
                            Z[where(mask)] = nan
                            den = (X+Y+Z)
                            x = X/den
                            y = Y/den
                            den = None
    
                            ## calculate alpha
                            alpha = mod(arctan2(y-yw, x-xw),2* pi)
                            x,y = None, None
                            alpha*=(180/pi)
                            hues_100 = alpha/100.
                            corr = (hac['coef'][0] * power(hues_100,5)) + \
                                   (hac['coef'][1] * power(hues_100,4)) + \
                                   (hac['coef'][2] * power(hues_100,3)) + \
                                   (hac['coef'][3] * power(hues_100,2)) + \
                                   (hac['coef'][4] * hues_100) + hac['coef'][5]
                            hues_100 = None
                            alpha += corr
                            corr = None
                            par_data = alpha
                            alpha = None
            ##### END hue angle
            #################################

            #################################
            ## OLH
            if par_name[0:3] == 'olh':
                par_exists = True
                mask_data = True
                par_split = par_name.split('_')
                par_attributes = {'algorithm':'Castagna et al. in prep'}

                ### get required datasets
                if gatts['sensor'] == 'L8_OLI':
                    req_waves = [561,613,655]
                    required_datasets = ['rhos_{}'.format(w) for w in req_waves]

                if len(required_datasets) > 0: 
                        req = l2w_required(inputfile, required_datasets, data_read, att_read)
                                                        
                        ## compute dataset
                        if req:
                            ow = (float(req_waves[2])-req_waves[1])/(float(req_waves[2])-float(req_waves[0]))
                            par_data = data_read[required_datasets[0]]*ow + data_read[required_datasets[2]]*(1-ow)
                            par_data = data_read[required_datasets[1]]-par_data
                            #par_data = None

            ##### END OLH
            #################################

            #################################
            ## sample
            if par_name[0:3] == 'sample':
                par_exists = True
                mask_data = False
                par_split = par_name.split('_')
                par_attributes = {'algorithm':'sample'}
                

                if len(required_datasets) > 0: 
                        req = l2w_required(inputfile, required_datasets, data_read, att_read)
                                                        
                        ## compute dataset
                        if req:
                            par_data = None
            ##### END sample
            #################################
            
            ## copy requested L2R variables
            if (par_data is None) and (par_name in l2r_datasets):
                if par_name in data_read: 
                    par_data = data_read[par_name]
                    par_attributes = att_read[par_name]
                else: 
                    par_data,par_attributes = nc_data(inputfile, par_name, attributes=True)
                if (par_name in ['bt10','bt11']) | ('rhot_' in par_name) | ('rhos_' in par_name) | ('rhorc_' in par_name) : 
                    mask_data = False
            ##
            if par_data is not None:
                par_computed = True

            ## if this parameter does not exist
            if not par_computed:
                if not par_exists: 
                    print('Parameter {} not configured.'.format(par))
                else: 
                    print('Parameter {} not computed for {}.'.format(par, inputfile))
                    if req is False: print('Required datasets for parameter {} are not in file'.format(par_name))
                continue

            ## create NetCDF
            if nc_new:
                nc_write(ncfile, 'lon', nc_data(inputfile, 'lon'), attributes=gatts, new=True,
                                   nc_compression=nc_compression, chunking=chunking)
                nc_write(ncfile, 'lat', nc_data(inputfile, 'lat'))

                ## write easting and northing if requested
                if 'x' in l2r_datasets: nc_write(ncfile, 'x', nc_data(inputfile, 'x'))
                if 'y' in l2r_datasets: nc_write(ncfile, 'y', nc_data(inputfile, 'y'))
                
                ## write flags
                if l2_flags is not None:
                    nc_write(ncfile, "l2_flags", l2_flags)
                    l2_flags = None

            ## mask output
            if (l2w_mask) & (l2w_mask_water_parameters):
                if (mask is not None) & (mask_data):
                    par_data[where(mask)] = nan

            ## write dataset to NetCDF
            nc_write(ncfile, par_name, par_data, dataset_attributes=par_attributes)
            nc_new=False
    return(ncfile)

