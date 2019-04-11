## def acolite_cli
## 
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-12-05
## modifications: 2018-03-07 (QV) added *args
##                2018-03-12 (QV) added path check
##                2018-03-22 (QV) added option for multiple images
##                2018-07-18 (QV) changed acolite import name
##                2018-09-12 (QV) added --nogfx option

def acolite_cli(*args):
    import os

    ## object for logging stdout to log file when processing
    class LogTee(object):
        def __init__(self, name):
            self.name=name
            ## make new file
            if os.path.exists(os.path.dirname(self.name)) is False:
                os.makedirs(os.path.dirname(self.name))
            self.file = open(self.name, 'w')
            self.file.close()
            self.mode='a'
            self.stdout = sys.stdout
            sys.stdout = self
        def __del__(self):
            sys.stdout = self.stdout
        def write(self, data):
            self.stdout.write(data)
            data = data.strip()
            if len(data) > 0:
                self.file = open(self.name, self.mode)
                self.file.write('{}: {}\n'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),data))
                self.file.close()
        def flush(self):
            pass

    import argparse
    parser = argparse.ArgumentParser(description='ACOLITE CLI')
    parser.add_argument('--settings', help='settings file')
    parser.add_argument('--images', help='list of images')
    args, unknown = parser.parse_known_args()

    if args.settings is None:
       print('No settings file given.')
       return(1)

    import sys, datetime
    from acolite.acolite import acolite_run, settings_read

    acolite_settings = settings_read(args.settings)
    if 'runid' not in acolite_settings:
        acolite_settings['runid'] = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    if 'output' not in acolite_settings:
        print('No output specified in settings file.')
        return(1)

    logfile = '{}/{}'.format(acolite_settings['output'],'acolite_run_{}_log.txt'.format(acolite_settings['runid']))
    log = LogTee(logfile)

    if '--nogfx' in unknown:
        print('Disabling matplotlib and graphical outputs.')
        acolite_settings['rgb_rhot']=False
        acolite_settings['rgb_rhos']=False
        acolite_settings['map_l2w']=False
        acolite_settings['dsf_plot_retrieved_tiles']=False
        acolite_settings['dsf_plot_dark_spectrum']=False

    print('Running ACOLITE')
    if args.images is None:
        acolite_run(settings=acolite_settings)
    else:
        images = args.images
        if ',' in images:
            images = images.split(',')

        if type(images) is not list: 
            if not os.path.isdir(images):
                with open(images, 'r') as f:
                    images = [line.strip() for line in f.readlines()]
            else:
                images = [images]

        acolite_settings['inputfile']=images
        acolite_run(settings=acolite_settings)

if __name__ == '__main__':
    acolite_cli()




