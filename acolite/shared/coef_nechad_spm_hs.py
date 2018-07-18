## def coef_nechad_spm_hs
## reads Nechad hyperspectral data
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-12-05
## modifications: 2018-03-07 (QV) renamed from nechad_spm_hs
##                2018-07-18 (QV) changed acolite import name

def coef_nechad_spm_hs():
    import os,sys
    from acolite import config
    pp_path = config['pp_data_dir']
    abfile = pp_path+'/Shared//REMSEM/Nechad_SPM_AB_2010.cfg'
    cfile = pp_path+'/Shared//REMSEM/Nechad_SPM_C_2010.cfg'

    data = {'wave':[], 'A':[], 'B':[], 'C':[], 'Rsq':[]}
    with open(abfile, 'r') as f:
        for line in f.readlines():
            if line[0] in [';','!', '/']: continue
            line = line.strip()
            split = line.split('\t')
            if len(split) != 4: continue
            data['wave'].append(float(split[0])/1000.)
            data['A'].append(float(split[1]))
            data['B'].append(float(split[2]))
            data['Rsq'].append(float(split[3]))
            
    with open(cfile, 'r') as f:
        i = 0
        for line in f.readlines():
            if line[0] in [';','!', '/']: continue
            line = line.strip()
            split = line.split('\t')
            if len(split) != 2: continue
            if data['wave'][i] != float(split[0])/1000.:
                print('Error')
            data['wave'].append(float(split[0])/1000.)
            data['C'].append(float(split[1])/100.)
            i+=1

    return(data)
