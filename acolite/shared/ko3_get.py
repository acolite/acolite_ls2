## def ko3_get
## reads ko3 data
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-04-14
## modifications:
##                2017-11-28 (QV) moved PP data directory
##                2018-07-18 (QV) changed acolite import name

def ko3_get(ko3file=None):
    if ko3file is None:
        import os,sys
        from acolite import config
        pp_path = config['pp_data_dir']
        ko3file = pp_path+'/Shared/k_o3_anderson.txt'
        
    ko3data=[]
    ko3wave=[]
    with open(ko3file, 'r') as f:
        for line in f:
            if line[0] == '!': continue
            if line[0] == '/': continue
            split = line.split(' ')
            if len(split) != 2: continue
            ko3data.append(float(split[1]))
            ko3wave.append(float(split[0])/1000.)
    ko3={"wave":ko3wave, "data":ko3data}
    return ko3
