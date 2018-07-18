## def similarity_ratio_sensor
## returns similarity spectrum ratio for given bands
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2018-04-17
## modifications: 
##                2018-07-18 (QV) changed acolite import name

def similarity_ratio_sensor(sensor, b1, b2):
    import acolite as ac
    
    ## get similarity spectrum table
    ssd = ac.shared.similarity_read()

    ## get rsr
    pp_path = ac.config['pp_data_dir']
    rsr_file = pp_path+'/RSR/{}.txt'.format(sensor)
    rsr,bands = ac.shared.rsr_read(file=rsr_file)

    ## make band averaged values
    band_averaged = ac.shared.rsr_convolute_dict(ssd['wave'], ssd['ave'], rsr)
    
    return(band_averaged[b1]/band_averaged[b2])
