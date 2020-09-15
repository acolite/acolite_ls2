## cams_read_grib
## reads ancillary data from CAMS NRT grib files
##
## written by Quinten Vanhellemont, RBINS
## 2020-09-15
## modifications:
##

def cams_read_grib(file, date, lon, lat, ftime=12,
                   spatial_selection = 'interpolated', verbosity = 0):
    import pygrib, os
    import numpy as np
    import scipy.interpolate

    anc = {}
    if os.path.exists(file):
        ## convert isodate into integer
        dateint = [int(s) for s in date.split('-')]
        dateint = (dateint[0]*10000) + (dateint[1]*100) + (dateint[2])

        ## open grib file
        grbs = pygrib.open(file)
        grbs.seek(0)

        for grb in grbs:
            if grb.date != dateint: continue

            if grb.name not in ['Total column water vapour',
                                'Mean sea level pressure',
                                'GEMS Total column ozone']:
                continue

            opar = None
            unit = grb.units
            par = grb.name
            tmp = grb.data()

            data = tmp[0]
            dlat = tmp[1]
            dlon = tmp[2]

            if par == 'Mean sea level pressure':
                ## Pa to hPa
                data/=100
                opar = 'press'
            elif par == 'GEMS Total column ozone':
                ## 1 DU = 2.1415E-5 kg m-2
                ## https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview
                ## to atmosphere cm
                data/=2.144e-2
                opar = 'ozone'
            elif par == 'Total column water vapour':
                ## to g/cm2
                data/=10
                opar = 'p_water'
            else:
                print('{} not configured'.format(par))
                continue

            ## find pixel containing the lat/lon position
            tmp = np.power(np.power(dlat-lat, 2) + np.power(dlon-lon, 2), 0.5)
            y, x = np.where(tmp == np.nanmin(tmp))
            if len(x)>1:x=x[0]
            if len(y)>1:y=y[0]
            y = y.squeeze()
            x = x.squeeze()

            if opar not in anc: anc[opar] = {'value':[], 'time':[]}

            ## get model data for position
            if spatial_selection == 'nearest':
                dv = data[y,x]
            elif spatial_selection == 'interpolated':
                ip = scipy.interpolate.interp2d(dlon[y-1:y+2,x-1:x+2],
                                                dlat[y-1:y+2,x-1:x+2],
                                                data[y-1:y+2,x-1:x+2])
                dv = ip(lon, lat)[0]
                ip = None

            ## get model time
            dt = grb.hour + grb.minute/60
            anc[opar]['time'].append(dt)
            anc[opar]['value'].append(dv)
        grbs.close()

        if len(anc) == 0:
            if verbosity > 0:
                print('Ancillary data for {} not in {}'.format(date, file))
        else:
            ## interpolate in time
            for opar in anc:
                anc[opar]['interp'] = np.interp(ftime, anc[opar]['time'],anc[opar]['value'])

    return(anc)
