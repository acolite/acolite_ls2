## def bundle_test
## lists files in given directory and returns dict with band and file path
## written by Quinten Vanhellemont, RBINS
## 2017-04-13
## modifications:QV 2018-10-01 added test for level name

def bundle_test(bundle):
    import os
    files = os.listdir(bundle)
    datafiles = {}
    for i, fname in enumerate(files):
        split = fname.split('_')
        tmp = split[-1].split('.')
        if len(tmp) is 2: band,ext = tmp
        else: continue
        if band in ['L1GST', 'L1T']: band = split[-2]
        if ext not in ['TIF', 'txt', 'TXT']: continue
        file = '{}/{}'.format(bundle,fname)
        datafiles[band] = {"path":file, "fname":fname}
    return(datafiles)
