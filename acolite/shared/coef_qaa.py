## def coef_qaa
## reads qaa coefficients
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2018-03-08
## modifications:
##                2018-07-18 (QV) changed acolite import name

def coef_qaa():
    import os,sys
    from acolite import config
    pp_path = config['pp_data_dir']
    nfile = pp_path+'/Shared/qaa_settings.cfg'

    data = {} 
    with open(nfile, 'r') as f:
        for line in f.readlines():
            if line[0] in [';','!', '/']: continue
            line = line.strip()
            if len(line)==0: continue
            split = line.split('=')
            if len(split)==2:
                if ',' in split[1]:
                    data[split[0]]=[float(d) for d in split[1].split(',')]
                else:
                    data[split[0]]=split[1]
    if 'useconfig' in data:
        for tag in ['g','h','k','l','m']:
            data[tag] = data['{}_{}'.format(tag,data['useconfig'])]
    return(data)
