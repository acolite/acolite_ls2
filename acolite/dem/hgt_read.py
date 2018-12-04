## def hgt_read
## reads DEM HGT SRTM files
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-07-17

def hgt_read(file):
    import struct
    from numpy import asarray
    
    if '.gz' in file:
        import gzip
        with gzip.open(file,'rb') as f:
            data_read = f.read()
    else:
        with open(file,'rb') as f:
            data_read = f.read()

    dim = (1201,1201)

    ## big endian, unsigned shorts
    data = struct.unpack('>{}'.format('H'*dim[0]*dim[1]),data_read)
    data = asarray(data).reshape(dim)
    
    data[data > 32768] -= 65535
    return(data)
