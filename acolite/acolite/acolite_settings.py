## def acolite_settings
## 
##
## Written by Quinten Vanhellemont 2017-11-30
## Last modifications:
##                2018-07-18 (QV) changed acolite import name
def acolite_settings(settings):
    import os 
    from acolite import acolite
    
    ## read defaults
    default_settings = acolite.config['acolite_defaults']
    setd = acolite.settings_read(default_settings)
    
    ## read settings file
    if settings is not None:
        ## path to settings file given
        if type(settings) is str:
            if os.path.exists(settings):
                setu = acolite.settings_read(settings)
            else:
                print('Settings file {} not found.'.format(settings))
                return(1)
        elif type(settings) is dict:
            setu = settings
        else:
            print('Settings not recognised.')
            return(1)
    else: setu={}
        
    ## set defaults for options not specified
    for key in setd.keys():
        if key in setu.keys(): continue
        else: setu[key] = setd[key]
    return(setu)
