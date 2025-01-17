## def coef_nechad2016
## reads Nechad 2016 calibration data
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2018-03-07
## modifications:
##                2018-07-18 (QV) changed acolite import name

def coef_nechad2016():
    import os,sys
    from acolite import config
    pp_path = config['pp_data_dir']
    nfile = pp_path+'/Shared/REMSEM/Nechad_calibration_201609.cfg'

    data = [] 
    with open(nfile, 'r') as f:
        for line in f.readlines():
            if line[0] in [';','!', '/']: continue
            line = line.strip()
            split = line.split('\t')
            if len(split) != 6: continue
            cd = {'par':split[0], 'sensor':split[1], 'band':split[2],
                  'wave':float(split[3]), 'A':float(split[4]), 'C':float(split[5])}
    
            data.append(cd)

    return(data)
