## def coef_nechad_hs
## reads Nechad hyperspectral data
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2018-05-15
## modifications:
##                2018-07-18 (QV) changed acolite import name

def coef_nechad_hs(par):
    import os,sys
    from acolite import config
    pp_path = config['pp_data_dir']
    
    if par == 'SPM':
        file = pp_path+'/Shared/REMSEM/SPM_N2010_Published.txt'
        
    if par == 'T':
        file = pp_path+'/Shared/REMSEM/Turbidity_N2009_Published.txt'

    keys = ['wave','A','B','Rsq','C']
    data = {k:[] for k in keys}
    with open(file,'r') as f:
        for line in f.readlines():
            if line[0] in ['#',';', '%']: continue
            sp = line.strip().split(',')
            if len(sp) != 5: continue
            for i,k in enumerate(keys):
                data[k].append(float(sp[i]))
            

    return(data)
