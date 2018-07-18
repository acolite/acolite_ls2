## def image_corners
## find the approximate image corners (not LL/UL RR/UR scene corners!) in the landsat scene 
## written by Quinten Vanhellemont, RBINS
## 2017-04-14
## modifications: 2017-06-06 (QV) use_band='5' instead of '6'
##                2018-07-18 (QV) changed acolite import name
def image_corners(bundle, metadata, use_band='5', use_rtoa=True, use_mean=True, return_geo=True):
    from numpy import isnan, where, mean
    from acolite.landsat import get_rtoa, read_band
    from acolite.landsat.geo import get_projection
        
    if use_rtoa:
        band_image =  get_rtoa(bundle, metadata, use_band)
        def test_sub(data):
            return((where(isnan(data) == False))[0])
    else:
        band_image =  read_band(metadata['B{}'.format(use_band)])
        def test_sub(data):
            return((where(data == 0))[0])
    dim = band_image.shape

    ## get top corner
    top = None
    row=0
    while top == None:
        sub = test_sub(band_image[row,:])
        row+=1
        if len(sub) == 0: continue
        if use_mean: sub=mean(sub)
        else: sub = min(sub)
        top = (row,sub)

    ## get bottom corner
    bottom = None
    row=dim[0]-1
    while bottom == None:
        sub = test_sub(band_image[row,:])
        row-=1
        if len(sub) == 0: continue
        if use_mean: sub=mean(sub)
        else: sub = max(sub)
        bottom = (row,sub)

    ## get left corner
    left = None
    col=0
    while left == None:
        sub = test_sub(band_image[:,col])
        col+=1
        if len(sub) == 0: continue
        if use_mean: sub=mean(sub)
        else: sub = min(sub)
        left = (sub,col)

    ## get right corner
    right = None
    col=dim[1]-1
    while right == None:
        sub = test_sub(band_image[:,col])
        col-=1
        if len(sub) == 0: continue
        if use_mean: sub=mean(sub)
        else: sub = max(sub)
        right = (sub,col)

    ## get intersects of nadir line with top and bottom edge
    nadir_top = ((top[0]+right[0])/2.,(top[1]+right[1])/2.)
    nadir_bottom = ((left[0]+bottom[0])/2.,(left[1]+bottom[1])/2.)
    nadir_middle = ((nadir_top[0]+nadir_bottom[0])/2.,(nadir_top[1]+nadir_bottom[1])/2.)
    
    if return_geo:
        ## projection info
        p, (xextent, yextent), proj4_string = get_projection(metadata)
        cell = float(metadata['GRID_CELL_SIZE_REFLECTIVE'])

        right_ll = p(xextent[0]+(right[1]*cell),yextent[1]-(right[0]*cell), inverse=True)
        left_ll = p(xextent[0]+(left[1]*cell),yextent[1]-(left[0]*cell), inverse=True)
        top_ll = p(xextent[0]+(top[1]*cell),yextent[1]-(top[0]*cell), inverse=True)
        bottom_ll = p(xextent[0]+(bottom[1]*cell),yextent[1]-(bottom[0]*cell), inverse=True)

        nadir_top_ll = p(xextent[0]+(nadir_top[1]*cell),yextent[1]-(nadir_top[0]*cell), inverse=True)
        nadir_bottom_ll = p(xextent[0]+(nadir_bottom[1]*cell),yextent[1]-(nadir_bottom[0]*cell), inverse=True)
        nadir_middle_ll = p(xextent[0]+(nadir_middle[1]*cell),yextent[1]-(nadir_middle[0]*cell), inverse=True)
        
        return(right_ll,left_ll,top_ll,bottom_ll, nadir_top_ll, nadir_bottom_ll,nadir_middle_ll)
    else:
        return(right,left,top,bottom,nadir_top,nadir_bottom,nadir_middle)

