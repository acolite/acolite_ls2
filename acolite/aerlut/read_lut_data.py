## def read_lut_data
## imports given luts for reuse over tiled processing
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-12-06
## modifications:
##                2018-07-18 (QV) changed acolite import name

def read_lut_data(sensor, lutdir=None, rsr_file=None, lut_data_dict=None, 
                         luts=['PONDER-LUT-201704-MOD1-1013mb', 
                               'PONDER-LUT-201704-MOD2-1013mb', 
                               'PONDER-LUT-201704-MOD3-1013mb'], pressure=True):
    import acolite as pp

    pp_path = pp.config['pp_data_dir']
    if lutdir is None: lutdir=pp_path+'/LUT/'
    if rsr_file is None: rsr_file = pp_path+'/RSR/'+sensor+'.txt'
    rsr, rsr_bands = pp.shared.rsr_read(file=rsr_file)
    
    ## make lut data dictionary that can be reused in next runs
    if lut_data_dict is None:
        lut_data_dict = {}
        for li,lut in enumerate(luts):
                ## get sensor LUT
                lut_sensor, meta_sensor = pp.aerlut.get_sensor_lut(sensor, rsr_file, lutdir=lutdir, lutid=lut, override=0)
                lut_data_dict[lut] = {'lut':lut_sensor, 'meta':meta_sensor}

                ## read luts at other pressures if needed
                if pressure:
                    lut_split = lut.split('-')
                    lut0 = '-'.join(lut_split[0:-1]+['0500mb'])
                    lut_sensor, meta_sensor = pp.aerlut.get_sensor_lut(sensor, rsr_file, lutdir=lutdir, lutid=lut0, override=0)
                    lut_data_dict[lut0] = {'lut':lut_sensor, 'meta':meta_sensor}

                    lut1 = '-'.join(lut_split[0:-1]+['1100mb'])
                    lut_sensor, meta_sensor = pp.aerlut.get_sensor_lut(sensor, rsr_file, lutdir=lutdir, lutid=lut1, override=0)
                    lut_data_dict[lut1] = {'lut':lut_sensor, 'meta':meta_sensor}
    return(lut_data_dict)
