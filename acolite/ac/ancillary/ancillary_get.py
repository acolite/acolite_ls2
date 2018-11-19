## ancillary_get
## downloads and interpolates ancillary data from the ocean data server
##
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-10-18
## modifications: 2017-10-24 (QV) added kind keyword
##                2017-12-05 (QV) added test if data is missing
##                2018-07-18 (QV) changed acolite import name
##                2018-11-19 (QV) added verbosity option
def ancillary_get(date, lon, lat, ftime=12.0, local_dir=None, quiet=True, kind='linear', verbosity=0):
    import acolite as pp
    
    if local_dir == None: local_dir=pp.config['met_dir']
    
    ## list and download files
    anc_files = pp.ac.ancillary.ancillary_list(date)
    anc_local = pp.ac.ancillary.ancillary_download(ancillary_files=anc_files, local_dir = local_dir, verbosity=verbosity)
    
    ## get ozone file
    auraomi_file = [anc_local[i] for i, j in enumerate(anc_local) if "AURAOMI" in j]
    eptoms_file = [anc_local[i] for i, j in enumerate(anc_local) if "EPTOMS" in j]
    toast_file = [anc_local[i] for i, j in enumerate(anc_local) if "TOAST" in j]
    
    ## use toast as fallback
    ozone_file=None
    if len(auraomi_file) == 1: ozone_file = auraomi_file[0]
    elif len(eptoms_file) == 1: ozone_file = eptoms_file[0]
    elif len(toast_file) == 1: ozone_file = toast_file[0]
    
    ## get ncep MET files
    ncep_files = [anc_local[i] for i, j in enumerate(anc_local) if "NCEPR2" in j]

    anc = {'date':date, 'lon':lon, 'lat': lat, 'ftime':ftime}

    ## interpolate ozone
    if ozone_file is None:
        print('No ozone file found for {}'.format(date))
    else:
        anc_ozone = pp.ac.ancillary.ancillary_interp_ozone(ozone_file, lon, lat, kind=kind)
        for k in anc_ozone.keys(): anc[k] = anc_ozone[k]

    ## interpolate MET
    if len(ncep_files) == 0:
        print('No NCEP files found for {}'.format(date))
    else:
        anc_met = pp.ac.ancillary.ancillary_interp_met(ncep_files,  lon, lat, ftime, kind=kind)
        for k in anc_met.keys(): anc[k] = anc_met[k]
    
    return(anc)
