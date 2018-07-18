## def view_azimuth
## computes approximate azimuth for nadir line for landsat scene
## written by Quinten Vanhellemont, RBINS
## 2017-04-14
## modifications:
##                2018-07-18 (QV) changed acolite import name
def view_azimuth(bundle, metadata):
    from acolite.landsat import image_corners
    from acolite.shared import azimuth_two_points

    ## get geolocation for image (not scene) corners 
    r,l,t,b,nadir_top, nadir_bottom,nadir_middle = image_corners(bundle, metadata)
    
    ## get azimuth between top and bottom nadir line
    azi = azimuth_two_points(nadir_top[0],nadir_top[1],nadir_bottom[0],nadir_bottom[1])
    return(azi)
