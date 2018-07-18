## def get_sub
## gets sub in Sentinel-2 image grids for given metadata and limit
## written by Quinten Vanhellemont, RBINS
## 2017-04-18
## modifications:  2017-05-10 (QV) update for granules not starting on XY divisible by 60
##                                 now limits sub to part inside granule
##                 2017-11-27 (QV) added extended xrange and yrange for limits extending the scene
##                 2018-01-24 (QV) added scmod offsets to xrange & yrange region
##                2018-07-18 (QV) changed acolite import name

def get_sub(metadata, limit):
    from acolite.sentinel.geo import get_projection
    p, grids, proj4_string = get_projection(metadata)
    xscene, yscene = grids['60']['xrange'],grids['60']['yrange']
    
    ## compute x and y limits, round to pixel size of coarsest grid
    pixelsize_coarse = [60.,-60.]
    xrange_raw, yrange_raw = p([limit[1],limit[3]],[limit[0],limit[2]])

    ## offsets in case the scene extent is not divisible by coarse resolution
    xscmod = (xscene[0] % pixelsize_coarse[0])
    yscmod = (yscene[0] % pixelsize_coarse[1])

    ## new region xrange
    xrange_coarse = [xrange_raw[0] - (xrange_raw[0] % pixelsize_coarse[0]) + xscmod, 
                      xrange_raw[1] - (xrange_raw[1] % pixelsize_coarse[0]) + xscmod]
    ## new region yrange
    yrange_coarse = [yrange_raw[0] - (yrange_raw[0] % pixelsize_coarse[1]) + yscmod, 
                      yrange_raw[1] - (yrange_raw[1] % pixelsize_coarse[1]) + yscmod]

    ## keep requested limit information (e.g. if outside granule)
    ## 2018-01-24 added scmod offsets
    xrange_region = [xrange_raw[0] - (xrange_raw[0] % pixelsize_coarse[0]) + xscmod, 
                      xrange_raw[1] - (xrange_raw[1] % pixelsize_coarse[0]) + xscmod]
    yrange_region = [yrange_raw[0] - (yrange_raw[0] % pixelsize_coarse[1])+ yscmod, 
                      yrange_raw[1] - (yrange_raw[1] % pixelsize_coarse[1])+ yscmod]
    grids_region = {'xrange':xrange_region, 'yrange':yrange_region}

    ## limit sub to granule (scene) extent
    if xrange_coarse[0] < xscene[0]: xrange_coarse[0] = xscene[0]
    if xrange_coarse[1] > xscene[1]: xrange_coarse[1] = xscene[1]
    if yrange_coarse[0] < yscene[1]: yrange_coarse[0] = yscene[1]
    if yrange_coarse[1] > yscene[0]: yrange_coarse[1] = yscene[0]

    grids_out = {}
    grids_out = {'granule_p':p, 'granule_grids':grids}
    
    for grid in grids.keys():
        dims = (grids[grid]['nx'],grids[grid]['ny'])
        pixelsize = (grids[grid]['xs'],grids[grid]['ys'])

        xrange=[i for i in xrange_coarse]
        yrange=[i for i in yrange_coarse]
        
        ## crop size
        x_size = int((xrange[1]-xrange[0])/pixelsize[0])
        y_size = int((yrange[1]-yrange[0])/abs(pixelsize[1]))
    
        if (xrange[1] < min(xscene)) or (xrange[0] > max(xscene)):
            print('Limits out of scene longitude')
            return(1)
        elif (yrange[1] < min(yscene)) or (yrange[0] > max(yscene)):
            print('Limits out of scene latitude')
            return(1)
        else:
            xoff = [(i - min(xscene))/pixelsize[0] for i in xrange]
            yoff = [(i - min(yscene))/abs(pixelsize[1]) for i in yrange]

            ## is this needed still?
            if xoff[0] < 0: xoff[0] = 0
            if yoff[0] < 0: yoff[0] = 0
            if xoff[1] >= dims[0]: xoff[1] = dims[0]#-1
            if yoff[1] >= dims[1]: yoff[1] = dims[1]#-1

            ## 2017-05-10
            yoff = [dims[1]-yoff[1], dims[1]-yoff[0]]
            
            ## xoff, yoff, nx, ny
            sub = [xoff[0], yoff[0], xoff[1]-xoff[0], yoff[1]-yoff[0]]
            sub = [int(s) for s in sub]

            grids_out[grid] = {'sub':sub, 'p':p, 'xrange':xrange, 'yrange':yrange, 
                                                 'xscene':xscene, 'yscene':yscene}
    ## compute grids for full region
    for grid in grids.keys():
        dims = (grids[grid]['nx'],grids[grid]['ny'])
        pixelsize = (grids[grid]['xs'],grids[grid]['ys'])

        xrange=[i for i in grids_region["xrange"]]
        yrange=[i for i in grids_region["yrange"]]
        
        ## crop size
        x_size = int((xrange[1]-xrange[0])/pixelsize[0])
        y_size = int((yrange[1]-yrange[0])/abs(pixelsize[1]))
        
        xoff = [(i - min(xscene))/pixelsize[0] for i in xrange]
        yoff = [(i - min(yscene))/abs(pixelsize[1]) for i in yrange]

        ## xoff, yoff, nx, ny
        yoff = [dims[1]-yoff[1], dims[1]-yoff[0]]
        sub = [xoff[0], yoff[0], xoff[1]-xoff[0], yoff[1]-yoff[0]]
        sub = [int(s) for s in sub]
        
        ## compute offset in region for this scene
        sub_scene = grids_out[grid]['sub']
        off = [sub_scene[0]-sub[0], sub_scene[1]-sub[1]]

        grids_region[grid] = {'sub':sub, 'off':off, 'p':p, 'xrange':xrange, 'yrange':yrange, 'dims':(x_size,y_size)}
            
            
    grids_out['grids_region']=grids_region
    return(grids_out, proj4_string)
