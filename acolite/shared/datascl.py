## def datascl
## linearly rescales data from the range [dmin,dmax] to [tmin,tmax], defaults to scaling to byte range and type
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2016-06-29
## modifications:

def datascl(data, dmin=None, dmax=None, tmin=0, tmax=255, dtype='uint8'):
    from numpy import interp
    if dmin == None: dmin = data.min()
    if dmax == None: dmax = data.max()
    data=interp(data, [dmin,dmax],[tmin,tmax])
    if dtype != None: data=data.astype(dtype)
    return data
