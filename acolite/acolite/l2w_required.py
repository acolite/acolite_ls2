## def l2w_required
## tests whether a netcdf file contains the required datasets, and reads them if all are present
##
## Written by Quinten Vanhellemont 2017-12-05
## Last modifications: 2018-04-17 (QV) tracking attributes from input NCDF
##                2018-07-18 (QV) changed acolite import name
def l2w_required(inputfile, required, data={}, att={}):
    from acolite.shared import nc_data,nc_datasets, nc_atts
    datasets = nc_datasets(inputfile)

    if all([rd in datasets for rd in required]):
        for rd in required:
            if rd not in data: 
                data[rd],att[rd] = nc_data(inputfile,rd, attributes=True)
                #att[rd] = nc_atts(inputfile, rd)
        return(True)
    else: return(False)
