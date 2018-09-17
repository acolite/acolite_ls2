## def qmap
## output map with lat/lon annotation, scale bar, colour bar
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-03
## modifications: 2017-03-09 (QV) added mask and projection options
##                2018-03-29 (QV) major changes, improved scalebar
##                2018-03-30 (QV) new annotation of meridian/parallel lines, added check for proximity to plot border
##                2018-04-06 (QV) new determination of plot window
##                2018-04-17 (QV) functionized, added cmap keyword
##                2018-04-18 (QV) new lat/lon step determination, "ideal" number of ticks
##                2018-07-18 (QV) changed acolite import name

def qmap(data,lon,lat, mask=None, outputfile=None, mscale=None, 
         colorbar_edge=True, colorbar = "horizontal", cbsize='3%', cbpad=0.5, cmap=None,
         limit=None, 
         scalebar=False, scalepos='UR', scalelen=None, scalecolor='Black', map_fillcolor='White',
         points=None, axes_linewidth=2, 
         range=None, label='', rescale=1., log=False,
         max_edge_tick=(0.33,0.25), max_scale_frac=0.33, verbose=False, geo_minticks=2, geo_maxticks=5,
         title=None, projection='tmerc', dpi=72, **kwargs):
    
    import time
    from acolite.shared import distance_in_ll, read_points, scale_dist
    
    def getstep(drange):
        dstep = 0.005
        if drange > 0.02: dstep=0.01
        if drange > 0.04: dstep=0.02
        if drange > 0.15: dstep=0.05
        if drange > 0.25: dstep=0.1
        if drange > 1.: dstep=0.5
        if drange > 5.: dstep=1.
        return(dstep)

    def getstep2(drange, mint=2, maxt=5):
        from acolite.shared import closest_idx
        steps = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2., 5., 10., 20., 50.]
        numb = [int(drange/s) for s in steps]
        idmax, ntmax = closest_idx(numb, maxt)
        idmin, ntmin = closest_idx(numb, mint)
        return(steps[idmin],steps[idmax])
   
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from mpl_toolkits.basemap import Basemap
    import matplotlib.text
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
 
    import pyproj
    from numpy.ma import masked_where
    from numpy import loadtxt, insert, arange
    
    if len(data.shape) is 2:
        y_size, x_size = data.shape 
    if len(data.shape) is 3:
        y_size, x_size, z_size = data.shape 
        
    x_size_r = x_size * rescale
    y_size_r = y_size * rescale
    
    #x0, y0, x1, y1
    plot_window = [0.,0.,1.,1.]
    
    ## these should be the approximate borders to the plot area
    border_y = [70, 120]
    border_x = [150, 100]
    
    ## add some space for the colorbar
    ## add to border_x and _y
    if (colorbar_edge):
        if (colorbar == 'vertical'):
            border_x[0] += 100
            border_x[1] += 400
            border_y[1] += 25


        if (colorbar == 'horizontal'):
            border_y[0] += 100
            border_y[1] += 150
            border_x[0] += 0
            border_x[1] += 25

        
    ## testing output shape
    if True:
        figsize = [x_size_r / dpi, y_size_r / dpi]
                
        border_y0 = border_y[0]/dpi
        border_y1 = border_y[1]/dpi
        figsize[1] += border_y0 + border_y1
        plot_window[1] = border_y0/figsize[1]
        plot_window[3] = 1-(border_y1/figsize[1])
        
        border_x0 = border_x[0]/dpi
        border_x1 = border_x[1]/dpi
        figsize[0] += border_x0 + border_x1
        plot_window[0] = border_x0/figsize[0]
        plot_window[2] = 1-(border_x1/figsize[0])
                
    if limit is not None:            
        if len(limit) is 4:
            lonw, lone = limit[1],limit[3]
            lats, latn = limit[0],limit[2]
    else:
        lone = lon[0,x_size-1]
        lonw = lon[y_size-1,0]
        latn = lat[0,x_size-1]
        lats = lat[y_size-1,0]

    lonm = (lonw+lone)/2
    latm = (lats+latn)/2
    
    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes(plot_window,zorder=20)
    
    ## set axis thickness
    for axis in ['top','bottom','left','right']:
      ax.spines[axis].set_linewidth(axes_linewidth)

    fig.patch.set_facecolor('white')
    fig.patch.set_alpha(1.0)
    
    lonrange = abs(lonw-lone)
    latrange = abs(latn-lats)
    t0 = time.time()
    if ('m' not in kwargs):
        if projection == 'wgs':
            m = Basemap(llcrnrlon=lon[x_size-1,0],llcrnrlat=lat[x_size-1,0],
                        urcrnrlon=lon[0,y_size-1],urcrnrlat=lat[0,y_size-1],
                    resolution = 'c', epsg=4326, suppress_ticks=True)
        if projection == 'tmerc':
            m = Basemap(llcrnrlon=lonw,llcrnrlat=lats,urcrnrlon=lone,urcrnrlat=latn,
                resolution = 'c', projection='tmerc', suppress_ticks=True, lat_0=latm,lon_0=lonm)
        if projection == 'cyl':
            m = Basemap(llcrnrlon=lonw,llcrnrlat=lats,urcrnrlon=lone,urcrnrlat=latn,
                resolution = 'c', projection=projection, suppress_ticks=True, lat_0=latm,lon_0=lonm)
        if projection == 'stere':

            proj = 'npstere' if latm > 0 else 'spstere'

            prj = Basemap(projection=proj, lon_0=lonm, lat_0=latm,
                          boundinglat=0, resolution='c')

            lonl = [lonw, lone, lonw, lone, lonm, lonm]
            latl = [lats, lats, latn, latn, lats, latn]
            xm, ym = prj(lonl, latl)
            ll_lon, ll_lat = prj(min(xm), min(ym), inverse=True)
            ur_lon, ur_lat = prj(max(xm), max(ym), inverse=True)
            m = Basemap(projection=projection, lat_ts=latm, lat_0=latm, lon_0=lonm, boundinglat=0,
                        llcrnrlon=lonw,llcrnrlat=lats,urcrnrlon=lone,urcrnrlat=latn, suppress_ticks=True, resolution = 'c')
    else:
        m = kwargs['m']
    t1 = time.time()
    if verbose: print('set up basemap {:.2f} sec'.format(t1-t0))
    
    # transform coordinates from lat/lon to map coordinates
    if ('xx' not in kwargs) | ('yy' not in kwargs):
        xx, yy = m(lon, lat)
    else:
        xx, yy = kwargs['xx'], kwargs['yy']
    t2 = time.time()
    if verbose: print('set up xx yy {:.2f} sec'.format(t2-t1))
    
    ## get LL, UR corners
    xll, yll = m(lonw,lats)
    xur, yur = m(lone,latn)

    xd = xur-xll
    yd = yur-yll

    ## different treatment of 2D and 3D arrays (maps and rgbs)
    if len(data.shape) is 2:
        # load colormap
        from matplotlib.colors import ListedColormap
        if cmap == None:
            cmap = cm.get_cmap('viridis') #loadtxt("Planck_Parchment_RGB.txt")/255.
        else:
            if type(cmap) == str:
                 cmap = cm.get_cmap(cmap) #loadtxt("Planck_Parchment_RGB.txt")/255.

        cmap.set_bad(map_fillcolor)
        cmap.set_under(map_fillcolor)

        if len(range) != 2: 
            from numpy import nanpercentile
            range = nanpercentile(data, (10,90))

        if log:
            from matplotlib.colors import LogNorm
            norm = LogNorm(vmin=range[0], vmax=range[1])
        else:
            norm = None

        if mask is None:
            img=m.pcolormesh(xx, yy, data, cmap=cmap, norm=norm, vmin=range[0], vmax=range[1]) #zorder=15, 
        else:
            img=m.pcolormesh(xx, yy, masked_where(mask,data), cmap=cmap, norm=norm, vmin=range[0], vmax=range[1])#zorder=15, 

    ## rgb colorTuple takes a bit longer
    if len(data.shape) is 3:
        mesh_rgb = data[:, :-1, :]
        colorTuple = mesh_rgb.reshape((mesh_rgb.shape[0] * mesh_rgb.shape[1]), 3)
        colorTuple = insert(colorTuple,3,1.0,axis=1)
        img = m.pcolormesh(xx, yy, data[:,:,0], latlon=False,color=colorTuple)

    t3 = time.time()
    if verbose: print('set up colormesh {:.2f} sec'.format(t3-t2))

    lonsteps = getstep2(lonrange, mint=geo_minticks, maxt=geo_maxticks)
    latsteps = getstep2(latrange, mint=geo_minticks, maxt=geo_maxticks)

    for ri in (0,1):
        lonstep = lonsteps[ri]
        lonedge=lonstep*max_edge_tick[0]
        meridians = arange(-180, 180, lonstep)
        meridians = meridians[(meridians >= lonw) & (meridians <= lone)]
        meridians = [mr for mr in meridians if (abs(mr-lone) > lonedge) & (abs(mr-lonw) > lonedge)]
        if len(meridians) >= geo_minticks: break

    for ri in (0,1):
        latstep = latsteps[ri]
        latedge=latstep*max_edge_tick[1]
        parallels = arange(-90, 90, latstep)
        parallels = parallels[(parallels >= lats) & (parallels <= latn)]
        parallels = [pl for pl in parallels if (abs(pl-lats) > latedge) & (abs(pl-latn) > latedge)]
        if len(parallels) >= geo_minticks: break
            

    flon = (len(str(lonstep))-2)
    xtl,xtv=[],[]
    for mr in meridians:
        xs, ys = m(mr,lats)
        xtl.append(('{'+':.{}'.format(flon)+'f}'+'°{}').format(abs(mr), 'E' if mr >=0 else 'W'))
        xtv.append(xs)
    ax.set_xticklabels(xtl)
    ax.set_xticks(xtv)

    flat = (len(str(latstep))-2)
    ytl,ytv=[],[]
    for pl in parallels:
        xs, ys = m(lonw,pl)
        ytl.append(('{'+':.{}'.format(flat)+'f}'+'°{}').format(abs(pl), 'N' if pl >=0 else 'S'))
        ytv.append(ys)
    ax.set_yticklabels(ytl)
    ax.set_yticks(ytv)

    t4 = time.time()
    if verbose: print('draw grid {:.2f} sec'.format(t4-t3))

    ## add scale bar
    if scalebar:
            if scalepos not in ['UR','UL','LL','LR']:
                print('Scalepos {} not recognised. Using default.'.format(scalepos))
                scalepos = 'UR'
                
            xlab_pix, ylab_pix = 100, 125
            ylab_off = (20/(figsize[1]*dpi))

            pos = {}
            
            pos['Upper'] = 1. - (ylab_pix/(figsize[1]*dpi))
            pos['Lower'] = (ylab_pix/(figsize[1]*dpi))
            
            pos['Right'] = 1. - (xlab_pix/(figsize[0]*dpi))
            pos['Left'] = (xlab_pix/(figsize[0]*dpi))
            

            if scalepos[0]=='U':
                latsc = lats+abs(latn-lats)*pos['Upper'] #0.87
            if scalepos[0]=='L':
                latsc = lats+abs(latn-lats)*pos['Lower'] #0.08
            if scalepos[1]=='R':
                lonsc = lonw+abs(lone-lonw)*pos['Right'] #0.92
                scale_sign=-1.
            if scalepos[1]=='L':
                lonsc = lonw+abs(lone-lonw)*pos['Left'] #0.08
                scale_sign=1.
            
            ## length of the scalebar in degrees
            dist = distance_in_ll(latm)[0]
            if scalelen is None:
                dd = dist*abs(lone-lonw)
                scalelen = dd * max_scale_frac #* rescale
                scalelen, unit = scale_dist(scalelen)
                sf = 1
                if unit == 'm':
                    scaleline = (scalelen / 1000) / dist
                else:
                    scaleline = scalelen / dist
            else:
                scalelen = float(scalelen)
                if scalelen < 1:
                    unit = 'm'
                    sf = 1000
                else:
                    unit = 'km'
                    sf = 1
                scaleline = scalelen / dist

            ## find start and end points in the scene
            xs, ys = m(lonsc,latsc)
            xse, ys = m(lonsc+scale_sign*scaleline,latsc)

            ## plot the line
            plt.plot((xs, xse), (ys,ys), '-', color=scalecolor, zorder=17, linewidth=2)

            ## add the label
            fz=14
            sclabel = '{} {}'.format(scalelen*sf, unit)
            if abs(scalelen-3.14)<0.01:
                sclabel = '{} {}'.format(r'$\pi$', unit)
            # label offset from scale bar
            sc_l_off = yd * ylab_off
            plt.text(xs+(xse-xs)/2, ys+sc_l_off, sclabel, color=scalecolor, zorder=17, 
                                                 horizontalalignment='center', fontsize=fz)
    ##### end scale bar
    
            
    ## plot provided points
    if points is not None:
        if type(points) is str:
            points = read_points(points)

        ypts_off = (50/(figsize[1]*dpi))
        xpts_off = (30/(figsize[0]*dpi))

        for pt in points:
            p=points[pt]
            sm = 'o' if 'sym' not in p else p['sym']
            cl = 'Black' if 'color' not in p else p['color']
            lb = None if 'label' not in p else p['label']
            fz = 14 if 'label_size' not in p else p['label_size']
            
            xp, yp = m(p['lon'], p['lat'])
            plt.plot(xp, yp, sm, color=cl, zorder=17)
        
            label_side = 'right' if 'label_side' not in p else p['label_side']
            va = 'center' if 'va' not in p else p['va']
            
            ## offset for labels
            h_l_off = xd * 0.025
            v_l_off = yd * 0.035

            h_l_off = xd * xpts_off
            v_l_off = yd * ypts_off

            if label_side == 'right':
                xpt, ypt = xp+h_l_off, yp*1
                ha = 'left' if 'ha' not in p else p['ha']
            if label_side == 'left':
                xpt, ypt = xp-h_l_off, yp*1
                ha = 'right' if 'ha' not in p else p['ha']
            if label_side == 'top':
                xpt, ypt = xp*1, yp+v_l_off
                ha = 'center' if 'ha' not in p else p['ha']
            if label_side == 'bottom':
                xpt, ypt = xp*1, yp-v_l_off
                ha = 'center' if 'ha' not in p else p['ha']

            plt.text(xpt, ypt, lb, color=cl, zorder=17, 
                     horizontalalignment=ha, 
                     verticalalignment=va, fontsize=fz)

    if title is not None:
        plt.title(title)

    ## make room for colorbar and plot if 2D dataset    
    divider = make_axes_locatable(ax)
    if colorbar_edge:
        if colorbar == "horizontal":
            cax = divider.append_axes("bottom", size=cbsize, pad=cbpad)
            orientation='horizontal'
        else:
            cax = divider.append_axes("right", size=cbsize, pad=cbpad)
            orientation='vertical'


        if len(data.shape) is 2:
            cbar = plt.colorbar(img, cax=cax, orientation=orientation)
            cbar.set_label(label=label)
        else: 
            cax.axis('off')
    else:
        pad = [0.1,0.1,0.1,0.3]
        pad = [0.15,0.15,0.15,0.15]
        for si, side in enumerate(['bottom','top','left', 'right']):
            cax = divider.append_axes(side, size='1%', pad=pad[si])
            cax.axis('off')
            
    if outputfile is None:
        plt.show()
    else:
        t5 = time.time()
        plt.savefig(outputfile, dpi=dpi) # , bbox_inches='tight'
        t6 = time.time()
        if verbose: print('write file {:.2f} sec'.format(t6-t5))
        if verbose: print('total {:.2f} sec'.format(t6-t0))

    plt.close()
    return xx, yy, m
