## def import_lut
## imports LUT made with 6SV and converts to NetCDF
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2016-07-05
## modifications:

def import_lut(lutid,lutdir,override=0):
    class Structure(object):
        pass

    import os, sys

    lutdat=lutdir+'/'+lutid+'/'+lutid+'.dat'
    lutdim=lutdir+'/'+lutid+'/'+lutid+'.dim'
    lutnc=lutdir+'/'+lutid+'/'+lutid+'.nc'

    if (os.path.isfile(lutnc)) & (override == 1): os.remove(lutnc)
        
    if (os.path.isfile(lutnc) == 0) & (os.path.isfile(lutdat)) & (os.path.isfile(lutdim)): 
        import csv, re
        import numpy as np

        ## read metadata
        try:
            file_open = open(lutdim, newline='')
            data = csv.reader(file_open, delimiter=';', quotechar='|')
            waves="0.39 0.41 0.44 0.47 0.51 0.55 0.61 0.67 0.75 0.865 1.04 1.24 1.61 2.25"
            ## for April 2017 simulations
            waves="0.39 0.41 0.44 0.47 0.51 0.55 0.61 0.67 0.75 0.865 1.04 1.24 1.55 1.61 1.66 2.10 2.25 2.40"
            wave = [float(i) for i in waves.split(' ')] 
            meta=Structure()
            meta.wave = wave
            for line in data:
                split = line[0].split('=')
                if len(split) < 2: continue
                if len(split[1]) == 0: continue
                if re.match('base=',line[0]): meta.base = split[1]
                if re.match('template=',line[0]): meta.template = split[1]
                if re.match('aermod=',line[0]): meta.aermod = split[1]
                if re.match('par=',line[0]): meta.par = [i for i in split[1].split('  ')]
                if re.match('tau=',line[0]): meta.tau = [float(i) for i in split[1].split(' ')]
                if re.match('ths=',line[0]): meta.ths = [float(i) for i in split[1].split(' ')]
                if re.match('thv=',line[0]): meta.thv = [float(i) for i in split[1].split(' ')]
                if re.match('azi=',line[0]): meta.azi = [float(i) for i in split[1].split(' ')]
                if re.match('wnd=',line[0]): meta.wnd = [float(i) for i in split[1].split(' ')]
                if re.match('elev=',line[0]): meta.elev = [float(i) for i in split[1].split(' ')]
                if re.match('press=',line[0]): meta.press = [float(i) for i in split[1].split(' ')]
                if re.match('MOD1=',line[0]): meta.mod1 = [float(i) for i in split[1].split(' ')]
                if re.match('MOD2=',line[0]): meta.mod2 = [float(i) for i in split[1].split(' ')]
                if re.match('MOD3=',line[0]): meta.mod3 = [float(i) for i in split[1].split(' ')]
                if re.match('MOD4=',line[0]): meta.mod4 = [float(i) for i in split[1].split(' ')]
            #print(meta.keys())
            #if 'wnd' in meta.keys(): 
            #    print('wnd')
            #    lut_dims = (len(meta.par), len(meta.wave), len(meta.azi), len(meta.thv), len(meta.ths), len(meta.wnd), len(meta.tau))
            #else: 
            #lut_dims = (len(meta.par), len(meta.wave), len(meta.azi), len(meta.thv), len(meta.ths), 1, len(meta.tau))
            meta.wnd=[-1]
            lut_dims = (len(meta.par), len(meta.wave), len(meta.azi), len(meta.thv), len(meta.ths), len(meta.wnd), len(meta.tau))
        except:
            print(sys.exc_info()[0])
            print('Failed to import LUT metadata (id='+lutid+')')
            
        ## read data
        try:
            file_open = open(lutdat, newline='')
            data = csv.reader(file_open)
            ### list of records
            arr = np.array([float(i) for line in data for i in line[0].split()])
            print(arr.shape)
            arr = np.reshape(arr,lut_dims,order='F') # reshape works
            ### record per simulation / parameters
            #arr = np.array([np.array(line[0].split()).astype(float) for line in data],order='F')
            #print(arr.shape)
            #arr = np.reshape(arr,lut_dims) # not correct
        except:
            print(sys.exc_info()[0])
            print('Failed to import LUT data (id='+lutid+')')

        ## write nc file
        try:
            if os.path.isfile(lutnc) == 0:
                from netCDF4 import Dataset
                nc = Dataset(lutnc, 'w', format='NETCDF4_CLASSIC')
                for i in meta.__dict__.keys(): 
                    attdata=meta.__dict__[i]
                    if isinstance(attdata,list):
                        if isinstance(attdata[0],str): 
                            attdata=','.join(attdata)
                    setattr(nc, i, attdata)
                #for i in meta.__dict__.keys(): print(meta.__dict__[i])

                nc.createDimension('par', lut_dims[0])
                nc.createDimension('wave', lut_dims[1])
                nc.createDimension('azi', lut_dims[2])
                nc.createDimension('thv', lut_dims[3])
                nc.createDimension('ths', lut_dims[4])
                nc.createDimension('wnd', lut_dims[5])
                nc.createDimension('tau', lut_dims[6])
                var = nc.createVariable('lut',float,('par','wave','azi','thv','ths','wnd','tau'))
                var[:] = arr
                nc.close()
                arr = None
                meta = None
        except:
            if os.path.isfile(lutnc): os.remove(lutnc)
            print(sys.exc_info()[0])
            print('Failed to write LUT data to NetCDF (id='+lutid+')')
            
    ## read dataset from NetCDF
    try:
        if os.path.isfile(lutnc):
            from netCDF4 import Dataset
            nc = Dataset(lutnc)
            meta=dict()
            for attr in nc.ncattrs():
                attdata = getattr(nc,attr)
                if isinstance(attdata,str): attdata = attdata.split(',')
                meta[attr]=attdata
            #meta = {attr : getattr(nc,attr) for attr in nc.ncattrs()}
            #for attr in meta.__dict__.keys():
            #    print(attr) #if isinstance(meta[attr],str): print(meta[attr].split(','))
            lut = nc.variables['lut'][:]
            nc.close()
    except:
        print(sys.exc_info()[0])
        print('Failed to open LUT data from NetCDF (id='+lutid+')')

        
    return lut, meta
