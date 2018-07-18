## def settings_write
## write "settings" file for ACOLITE based on settings dict
## 
## Written by Quinten Vanhellemont 2017-11-30
## Last modifications:

def settings_write(file, settings):
    import os
    if os.path.exists(os.path.dirname(file)) is False: os.makedirs(os.path.dirname(file)) 
        
    import datetime
    
    comf='## {}\n'
    valf='{}={}\n'
    with open(file,'w') as f:
        f.write(comf.format('ACOLITE Python settings'))
        f.write(comf.format('Written at {}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))))
        for key in settings.keys():
            vals = settings[key]
            if type(vals) is list: 
                vals = ','.join([str(i) for i in settings[key]])
            f.write(valf.format(key,vals))
