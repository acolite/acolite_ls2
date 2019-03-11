## def acolite_map
## outputs maps from given L2R/W NetCDF file
##
## Written by Quinten Vanhellemont 2017-12-05
## Last modifications: 2018-03-08 (QV) updated output directory options
##                     2018-04-17 (QV) added autoscale for * in the config files, added 'mapped' option using qmap
##                     2018-04-18 (QV) integrated RGB mapping, added dataset rescaling
##                     2018-04-19 (QV) added RGB pan sharpening option
##                     2018-07-18 (QV) changed acolite import name
##                     2018-10-24 (QV) changed matplotlib figure call
##                     2019-03-06 (QV) added plt.close
##                     2019-03-11 (QV) added strip to parameters

def acolite_map(inputfile=None, output=None, parameters=None, 
                dpi=300, ext='png', mapped=True, max_dim = 1000, limit=None,
                auto_range=False, range_percentiles=(5,95), dataset_rescale=False,
                map_title=True, 
                map_colorbar=False, map_colorbar_orientation='vertical',#'horizontal', 
                rgb_rhot = False, rgb_rhos = False, 
                red_wl = 660, green_wl = 560, blue_wl = 480, rgb_min = [0.0]*3, rgb_max = [0.15]*3, rgb_pan_sharpen = False, map_parameters_pan=True,
                map_fillcolor='White',
                map_scalepos = 'LR', map_scalebar = True, map_scalecolor='Black', map_scalecolor_rgb='White', map_scalelen=None, map_projection='tmerc',
                map_colorbar_edge=True, map_points=None, return_image=False, map_raster=False):

    import os, copy
    import datetime, time, dateutil.parser

    from acolite.shared import datascl,nc_data,nc_datasets,nc_gatts,qmap,closest_idx
    from acolite.acolite import pscale
    import acolite as ac

    from numpy import nanpercentile, log10, isnan, dstack
    from scipy.ndimage import zoom

    import matplotlib
    import matplotlib.pyplot as plt


    if not os.path.exists(inputfile):
        print('File {} not found.'.format(inputfile))
        return(False)

    ## run through maps
    maps = {'rhot':rgb_rhot,'rhos':rgb_rhos, 'parameters':parameters != None}
    if all([maps[m] == False for m in maps]): return()
        
    ## get parameter scaling
    psc = pscale()

    ## read netcdf info    
    l2w_datasets = nc_datasets(inputfile)
    print(l2w_datasets)

    gatts = nc_gatts(inputfile)
    if 'MISSION_INDEX' in gatts:
        sat, sen = gatts['MISSION'], gatts['MISSION_INDEX']
        stime = dateutil.parser.parse(gatts['IMAGING_DATE']+' '+gatts['IMAGING_TIME']) 
        obase = '{}_{}_{}'.format(sat, sen, stime.strftime('%Y_%m_%d_%H_%M_%S'))
    else:
        sp = gatts['sensor'].split('_') if 'sensor' in gatts else gatts['SATELLITE_SENSOR'].split('_')  
        sat, sen = sp[0], sp[1]
        stime = dateutil.parser.parse(gatts['isodate'] if 'isodate' in gatts else gatts['ISODATE']) 
        obase = gatts['output_name'] if 'output_name' in gatts else gatts['obase']

    ## find pan sharpening dataset
    if rgb_pan_sharpen:
        if sat not in ['L7','L8']: rgb_pan_sharpen = False
        tmp = os.path.splitext(inputfile)
        l1_pan_ncdf = '{}L1R_pan{}'.format(tmp[0][0:-3],tmp[1])
        if os.path.exists(l1_pan_ncdf):
            pan_data = nc_data(l1_pan_ncdf, 'rhot_pan')
        else:
            print('L1 pan NetCDF file not found')
            rgb_pan_sharpen=False

    if output is not None:
        odir = output
    else:
        odir = gatts['output_dir']
        
    if not os.path.exists(odir): os.makedirs(odir)

    scf= 1.
    rescale = 1.0

    #if dataset_rescale or mapped:
    lon = nc_data(inputfile, 'lon')
    if mapped: 
        lat = nc_data(inputfile, 'lat')
        if rgb_pan_sharpen:
            lon_pan = zoom(lon, zoom=2, order=1)
            lat_pan = zoom(lat, zoom=2, order=1)

    ## set up mapping info
    if True:
        from numpy import linspace, tile, ceil, isnan, nan
        mask_val = -9999.9999
        from scipy.ndimage.interpolation import map_coordinates

        ## rescale to save memory
        dims = lon.shape
        dsc = (dims[0]/max_dim, dims[1]/max_dim)
        scf/=max(dsc)

        if rgb_pan_sharpen: scf = 1.0

        if (scf < 1.) and dataset_rescale:
            sc_dims = (int(ceil(dims[0] * scf)), int(ceil(dims[1] * scf)))
            xdim =  linspace(0,dims[1],sc_dims[1]).reshape(1,sc_dims[1])
            ydim =  linspace(0,dims[0],sc_dims[0]).reshape(sc_dims[0],1)
            xdim = tile(xdim, (sc_dims[0],1))
            ydim = tile(ydim, (1,sc_dims[1]))

            resc = [ydim,xdim]
            xdim, ydim = None, None
            lon = map_coordinates(lon, resc, mode='nearest')
            lat = map_coordinates(lat, resc, mode='nearest')
        else:
            rescale = scf

    ## run through parameters
    for mi in maps:
        if not maps[mi]: continue

        if mi == 'parameters':
            if rgb_pan_sharpen:
                if map_parameters_pan & mapped:
                    lon = lon_pan * 1.0
                    lon_pan = None
                    lat = lat_pan * 1.0
                    lat_pan = None
                pan_data, lon_pan, lat_pan = None, None, None

            print('Mapping {}'.format(mi))
            if type(parameters) is not list: parameters=[parameters]
            for pid, par in enumerate(parameters):
                par = par.strip()
                pard = None

                ## check if this parameter exists
                if par not in l2w_datasets:
                    print('Parameter {} not in file {}.'.format(par, inputfile))
                    continue
                    
                print('Mapping {}'.format(par))
                ## read data
                data = nc_data(inputfile, par)
                if (rgb_pan_sharpen) & (map_parameters_pan):
                    data = zoom(data, zoom=2, order=1)

                ## rescale data
                if (scf != 1.0) and dataset_rescale:
                    data[isnan(data)] = mask_val
                    data = map_coordinates(data, resc, cval=mask_val)
                    data[data <= int(mask_val)] = nan
                    data[data <= 0] = nan

                data_range = nanpercentile(data, range_percentiles)

                ## get parameter mapping configuration
                if par in psc:
                    pard = copy.deepcopy(psc[par])
                else:
                     tmp = par.split('_')
                     par_generic = '_'.join((tmp[0:-1]+['*']))
                     if par_generic in psc: 
                         pard = copy.deepcopy(psc[par_generic])
                         try: ## add wavelength to generic name
                             wave = int(tmp[len(tmp)-1])
                             pard['name'] = '{} ({} nm)'.format(pard['name'], wave)
                         except:
                             pass
                     else: pard= {'color table':'default', 'min':data_range[0], 'max':data_range[1],
                                  'log': False, 'name':par, 'unit':'', 'parameter':par}

                if pard['color table'] == 'default': pard['color table']='viridis'
                ctfile = "{}/{}/{}.txt".format(ac.config['pp_data_dir'], 'Shared/ColourTables', pard['color table'])

                if os.path.exists(ctfile):
                    from matplotlib.colors import ListedColormap
                    from numpy import loadtxt
                    pard['color table'] = ListedColormap(loadtxt(ctfile)/255.)

                if 'title' not in pard: pard['title']='{} [{}]'.format(pard['name'],pard['unit'])
                if auto_range:
                    pard['min']=data_range[0]
                    pard['max']=data_range[1]

                if isnan(pard['min']): pard['min']=data_range[0]
                if isnan(pard['max']): pard['max']=data_range[1]

                ## outputfile
                outputfile = '{}/{}_{}.png'.format(odir,obase,par)

                if map_title:
                    title = '{} {}/{} {}'.format(pard['name'], sat, sen, stime.strftime('%Y-%m-%d (%H:%M UTC)'))
                else:
                    title = None

                ## use qmap option
                if mapped:
                    range = (pard['min'], pard['max'])
                    if 'limit' in gatts:
                        limit = gatts['limit']

                    if ('xx' not in locals()):
                        xx, yy, m = qmap(data, lon, lat, outputfile=outputfile, title=title, rescale=rescale,
                                           colorbar=map_colorbar_orientation, colorbar_edge=map_colorbar_edge, cmap=pard['color table'],
                                           label=pard['title'], range=range, log = pard['log'], map_fillcolor=map_fillcolor,
                                           limit=limit, dpi=dpi, points=map_points, projection=map_projection,
                                           scalebar=map_scalebar, scalepos=map_scalepos, 
                                           scalecolor=map_scalecolor, scalelen=map_scalelen)                
                    else:
                        xx, yy, m = qmap(data, lon, lat, outputfile=outputfile, title=title, rescale=rescale,
                                           colorbar=map_colorbar_orientation, colorbar_edge=map_colorbar_edge, cmap=pard['color table'],
                                           label=pard['title'], range=range, log = pard['log'], map_fillcolor=map_fillcolor, 
                                           limit=limit, dpi=dpi, points=map_points, projection=map_projection,
                                           scalebar=map_scalebar, scalepos=map_scalepos, 
                                           scalecolor=map_scalecolor, scalelen=map_scalelen, xx=xx, yy=yy, m=m)

                else:
                    import matplotlib.cm as cm
                    from matplotlib.colors import ListedColormap
                    cmap = cm.get_cmap(pard['color table'])
                    cmap.set_bad(map_fillcolor)
                    cmap.set_under(map_fillcolor)

                    if not map_raster:
                        ## set up plot
                        fig = plt.figure()
                        canvas = matplotlib.backends.backend_agg.FigureCanvasAgg(fig)
                        ax = fig.add_subplot(111)

                        print(pard['min'], pard['max'])

                        if pard['log']:
                            from matplotlib.colors import LogNorm
                            cax = ax.imshow(data, vmin=pard['min'], vmax=pard['max'], cmap=cmap,
                                               norm=LogNorm(vmin=pard['min'], vmax=pard['max']))
                        else:
                            cax = ax.imshow(data, vmin=pard['min'], vmax=pard['max'], cmap=cmap)

                        if map_colorbar:
                            if map_colorbar_orientation == 'vertical':
                                cbar = fig.colorbar(cax, orientation='vertical')
                                cbar.ax.set_ylabel(pard['title'])
                            else:
                                cbar = fig.colorbar(cax, orientation='horizontal')
                                cbar.ax.set_xlabel(pard['title'])

                            if map_title: ax.set_title(title)
                            ax.axis('off')
                            canvas.print_figure(outputfile, dpi=dpi, bbox_inches='tight')
                        plt.close()
                    else:
                        from PIL import Image
                        ## rescale for mapping
                        if pard['log']:
                            from numpy import log10
                            datasc = datascl(log10(data), dmin=log10(pard['min']), dmax=log10(pard['max']))
                        else:
                            datasc = datascl(data, dmin=pard['min'], dmax=pard['max'])

                        d = cmap(datasc)
                        for wi in (0,1,2):
                            ## convert back to 8 bit channels (not ideal)
                            d_ = datascl(d[:,:,wi], dmin=0, dmax=1)
                            if wi == 0: im = d_
                            else: im = dstack((im,d_))

                        img = Image.fromarray(im)

                        ## output image    
                        img.save(outputfile)

                print('Wrote {}'.format(outputfile))
        else:
            print('Mapping RGB {}'.format(mi))
            ## RGBs
            waves = [float(ds.split('_')[1]) for ds in l2w_datasets if ds[0:4] == mi]
            if len(waves) == 0:
                print('No appropriate datasets found for RGB {} in {}'.format(mi, inputfile))
                continue

            ## read datasets
            for wi, wl in enumerate([red_wl, green_wl, blue_wl]):
                idx, wave = closest_idx(waves, wl)
                cpar = '{}_{}'.format(mi, int(wave))
                ## read data
                data = nc_data(inputfile, cpar)

                if rgb_pan_sharpen:
                    data = zoom(data, zoom=2, order=1)
                    if wi == 0: vis_i = data * 1.0
                    else: vis_i += data
                    if wi == 2:
                        vis_i /= 3
                        pan_i = vis_i/pan_data
                        vis_i = None

                ## rescale data
                if (scf != 1.0) and dataset_rescale:
                    data[isnan(data)] = mask_val
                    data = map_coordinates(data, resc, cval=mask_val)
                    data[data <= int(mask_val)] = nan
                    data[data <= 0] = nan

                ## stack image
                if wi == 0:
                    image = data
                else:
                    image = dstack((image,data))
                
            ## rescale data between 0 and 1
            for wi in (2,1,0):
                if rgb_pan_sharpen: image[:,:,wi] /= pan_i
                image[:,:,wi] = datascl(image[:,:,wi], dmin=rgb_min[wi], dmax=rgb_max[wi])/255.

            par = r'$\rho_{}$'.format(mi[3]) + ' RGB'
            if map_title:
                title = '{} {}/{} {}'.format(par, sat, sen, stime.strftime('%Y-%m-%d (%H:%M UTC)'))
            else:
                title = None

            ## outputfile
            if rgb_pan_sharpen: 
                outputfile = '{}/{}_rgb_{}_pan.png'.format(odir,obase,mi)
            else:
                outputfile = '{}/{}_rgb_{}.png'.format(odir,obase,mi)

            # use qmap option
            if mapped:
                if 'limit' in gatts:
                    limit = gatts['limit']

                if rgb_pan_sharpen:
                    ret = qmap(image, lon_pan, lat_pan, outputfile=outputfile, title=title, rescale=rescale,
                                               colorbar=map_colorbar_orientation, colorbar_edge=map_colorbar_edge,
                                               limit=limit, dpi=dpi, points=map_points, projection=map_projection,
                                               scalebar=map_scalebar, scalepos=map_scalepos, 
                                               scalecolor=map_scalecolor_rgb, scalelen=map_scalelen)      
                    ret = None     
                else:
                    if ('xx' not in locals()):
                        xx, yy, m = qmap(image, lon, lat, outputfile=outputfile, title=title, rescale=rescale,
                                               colorbar=map_colorbar_orientation, colorbar_edge=map_colorbar_edge,
                                               limit=limit, dpi=dpi, points=map_points, projection=map_projection,
                                               scalebar=map_scalebar, scalepos=map_scalepos, 
                                               scalecolor=map_scalecolor_rgb, scalelen=map_scalelen)                
                    else:
                        xx, yy, m = qmap(image, lon, lat, outputfile=outputfile, title=title, rescale=rescale,
                                               colorbar=map_colorbar_orientation, colorbar_edge=map_colorbar_edge,
                                               limit=limit, dpi=dpi, points=map_points, projection=map_projection,
                                               scalebar=map_scalebar, scalepos=map_scalepos, 
                                               scalecolor=map_scalecolor_rgb, scalelen=map_scalelen, xx=xx, yy=yy, m=m)


            else:
                if not map_raster:
                    ## set up plot
                    fig = plt.figure()
                    canvas = matplotlib.backends.backend_agg.FigureCanvasAgg(fig)
                    ax = fig.add_subplot(111)
                    ax.imshow(image)
                    image = None
                    
                    if map_title: ax.set_title(title)
                    ax.axis('off')
                    canvas.print_figure(outputfile, dpi=dpi, bbox_inches='tight')
                    plt.close()
                else:
                    from PIL import Image
                    for wi in (0,1,2):
                        # convert again to 8 bit channels (not ideal)
                        data = datascl(image[:,:,wi], dmin=0, dmax=1)
                        if wi == 0:
                            im = data
                        else:
                            im = dstack((im,data))

                    img = Image.fromarray(im)
                    img.save(outputfile)

            print('Wrote {}'.format(outputfile))
