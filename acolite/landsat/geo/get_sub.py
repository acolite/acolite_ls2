## def get_sub
## gets sub in landsat image for given metadata and limit
## written by Quinten Vanhellemont, RBINS
## 2017-04-13
## modifications: 2017-06-22 (QV) fixed xoff/yoff bug
##                2017-10-25 (QV) added polar stereographic projection crop (min and max x and y ranges are selected)
##                                Maybe also better for UTM?
##                2018-02-08 (QV) tried fixing the cropping southern Y edge error
##                 2018-06-06 (QV) added return of Proj4 string
##                2018-07-18 (QV) changed acolite import name
##                2018-10-01 (QV) added grid cell size option
##                2019-03-12 (QV) changed the Y extent crop
##                2019-04-09 (QV) changed the region x/y sizes

def get_sub(metadata, limit):
    from acolite.landsat.geo import get_projection

    dims = metadata["DIMS"]

    if 'GRID_CELL_SIZE_REFLECTIVE' in metadata:
        pixelsize = [float(metadata["GRID_CELL_SIZE_REFLECTIVE"])]*2
    elif 'GRID_CELL_SIZE_REF' in metadata:
        pixelsize = [float(metadata["GRID_CELL_SIZE_REF"])]*2
    else:
        return(1)

    p, (xscene,yscene), proj4_string = get_projection(metadata)
        
    ## compute x and y limits, round to pixel sizes
    if metadata['MAP_PROJECTION'] == 'UTM':
        xrange_raw, yrange_raw = p([limit[1],limit[3]],[limit[0],limit[2]])
        xrange = [xrange_raw[0] - (xrange_raw[0] % pixelsize[0]), xrange_raw[1]+pixelsize[0]-(xrange_raw[1] % pixelsize[0])]
        yrange = [yrange_raw[0] - (yrange_raw[0] % pixelsize[1]), yrange_raw[1]+pixelsize[1]-(yrange_raw[1] % pixelsize[1])]

    if metadata['MAP_PROJECTION'] == 'PS':     
        # West, West, East, East
        # South, North, North, South
        xrange_raw, yrange_raw = p((limit[1],limit[1],limit[3],limit[3]),
                                   (limit[0],limit[2],limit[2],limit[0]))
        
        #
        xrange_raw = (min(xrange_raw), max(xrange_raw))
        yrange_raw = (min(yrange_raw), max(yrange_raw))
        
        xrange = [xrange_raw[0] - (xrange_raw[0] % pixelsize[0]), xrange_raw[1]+pixelsize[0]-(xrange_raw[1] % pixelsize[0])]
        yrange = [yrange_raw[0] - (yrange_raw[0] % pixelsize[1]), yrange_raw[1]+pixelsize[1]-(yrange_raw[1] % pixelsize[1])]

        

    ## crop size
    #x_size = int((xrange[1]-xrange[0])/pixelsize[0])
    #y_size = int((yrange[1]-yrange[0])/pixelsize[1])-1
    
    ## 7 june 2018
    x_size = int((xrange[1]-xrange[0])/pixelsize[0])+1
    y_size = int((yrange[1]-yrange[0])/pixelsize[1])+1

    ## 9 april 2019
    x_size = int((xrange[1]-xrange[0])/pixelsize[0])+1
    y_size = int((yrange[1]-yrange[0])/pixelsize[1])-1#+1

    grid_region = {'dims':(x_size,y_size), 'xrange':xrange, 'yrange':yrange}

    if (xrange[1] < min(xscene)) or (xrange[0] > max(xscene)):
        print('Limits out of scene longitude')
        return(1)
    elif (yrange[1] < min(yscene)) or (yrange[0] > max(yscene)):
        print('Limits out of scene latitude')
        return(1)
    else:
        xoff = [(i - min(xscene))/pixelsize[0] for i in xrange]
        #yoff = [(i - min(yscene))/pixelsize[1] for i in yrange]
        ## 6 june 2018 add 1 pixel for the reversal later
        yoff = [((i+pixelsize[1]) - min(yscene))/pixelsize[1] for i in yrange]

        xoff_region = [(i - min(xscene))/pixelsize[0] for i in xrange]
        #yoff_region = [(i - min(yscene))/pixelsize[1] for i in yrange]
        ## 6 june 2018  add 1 pixel for the reversal later
        yoff_region = [((i+pixelsize[1]) - min(yscene))/pixelsize[1] for i in yrange]

        if xoff[0] < 0: xoff[0] = 0
        if yoff[0] < 0: yoff[0] = 0
        if xoff[1] >= dims[0]: xoff[1] = dims[0]-1
        if yoff[1] >= dims[1]: yoff[1] = dims[1]-1

        ## flip the y subset
        yoff = [dims[1]-yoff[1], dims[1]-yoff[0]]
        yoff_region = [dims[1]-yoff_region[1], dims[1]-yoff_region[0]]

        ## 12 Mar 2019
        sub = [xoff[0], yoff[0], xoff[1]-xoff[0]+1, yoff[1]-yoff[0]-1]
        sub = [int(s) for s in sub]
        sub_region = [xoff_region[0], yoff_region[0], xoff_region[1]-xoff_region[0]+1, yoff_region[1]-yoff_region[0]-1]
        sub_region = [int(s) for s in sub_region]

        off = [sub[0]-sub_region[0], sub[1]-sub_region[1]]
        ## 7 june
        grid_region['off'] = off
        grid_region['sub'] = sub_region

        return sub, p, (xrange,yrange,grid_region), proj4_string
