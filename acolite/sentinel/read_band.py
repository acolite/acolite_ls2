## def read_band
## simple image reading for Sentinel band files
## sub keyword (xoff, yoff, xcount, ycount)
##
## written by Quinten Vanhellemont, RBINS 
## 2017-04-18
## modifications: 

def read_band(file, sub=None):
    import os, sys, fnmatch
    if not os.path.isfile(file):
        print('File '+file+' not found.')
        sys.exit()
       
    if fnmatch.fnmatch(file,'*.jp2'):
        from osgeo import gdal
        gdal.UseExceptions()
        band = gdal.Open(file)
        nrows=band.RasterYSize
        ncols=band.RasterXSize

        if sub is None:
            data = band.ReadAsArray()
        else:
             data = band.ReadAsArray(sub[0],sub[1],sub[2],sub[3])
        ds = None

    return(data)
