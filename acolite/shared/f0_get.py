## def f0_get
## reads f0 data
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-01-24
## modifications:
##                2017-11-28 (QV) moved PP data directory
##                2018-07-18 (QV) changed acolite import name
##                2018-09-19 (QV) added encoding

def f0_get(f0file=None):

    if f0file is None:
        import os,sys
        from acolite import config
        pp_path = config['pp_data_dir']
        f0file = pp_path+'/Shared/f0.txt'

    f0data=[]
    f0wave=[]
    with open(f0file, 'r', encoding='utf-8') as f:
        for line in f:
            if line[0] == '!': continue
            if line[0] == '/': continue
            split = line.split(' ')
            if len(split) != 2: continue
            f0data.append(float(split[1]))
            f0wave.append(float(split[0]))
    f0={"wave":f0wave, "data":f0data}
    return f0
