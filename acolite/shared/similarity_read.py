## def similarity_read
## reads similarity spectrum
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-12
## modifications: 2018-04-17 (QV) added to acolite function
##                2018-07-18 (QV) changed acolite import name

def similarity_read(file=None):
    if file == None:
        import acolite as ac 
        file = '{}/Shared/REMSEM/similarityspectrumtable.txt'.format(ac.config['pp_data_dir'])
    
    from numpy import argsort, ndarray, append

    ss_data = {'wave':ndarray(0), 'ave':ndarray(0), 'std':ndarray(0), 'cv':ndarray(0)}

    with open(file, 'r') as f:
        for i,line in enumerate(f.readlines()):
            line = line.strip()
            sp = line.split()
            if i == 0: 
                continue
            else:
                j = 0
                while j<len(sp):
                    ss_data['wave']=append(ss_data['wave'],float(sp[j])/1000.)
                    ss_data['ave']=append(ss_data['ave'],float(sp[j+1]))
                    ss_data['std']=append(ss_data['std'],float(sp[j+2]))
                    ss_data['cv']=append(ss_data['cv'],float(sp[j+3]))
                    j+=4

    ## sort indices
    idx = argsort(ss_data['wave'])
    for k in ss_data.keys(): ss_data[k]=ss_data[k][idx]
    return(ss_data)
