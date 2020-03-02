## def sensor_wave
## gets band wavelengths based on RSR
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2018-06-07
## modifications: 
##                2018-07-18 (QV) changed acolite import name
##                2020-03-02 (QV) added ceil to number of elements in linspace

def sensor_wave(satsen, wave_range = (0.39,2.4), wave_step = 0.001):
    from numpy import linspace, ceil
    from acolite import config, rsr_convolute_dict, rsr_read

    pp_path = config['pp_data_dir']
    rsr_file = pp_path+'/RSR/'+satsen+'.txt'
    rsr, rsr_bands = rsr_read(file=rsr_file)
    
    wave_hyper = linspace(wave_range[0],wave_range[1],ceil(((wave_range[1]-wave_range[0])/wave_step)+2))

    rsr_wave = rsr_convolute_dict(wave_hyper, wave_hyper, rsr)
    waves = {rband:'{:.0f}'.format(round(rsr_wave[rband]*1000.)) for rband in rsr_wave}
    return(waves)
