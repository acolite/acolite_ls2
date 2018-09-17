## wrapper to launch ACOLITE GUI/CLI
## QV 2018
## last modifications QV 2018-09-12 renamed from acolite.py and added main test

def launch_acolite():
    ## import sys to parse arguments
    import sys

    ## fix matplotlib backend to Agg
    ## skip import if --nogfx is given
    if not (('--cli' in sys.argv) & ('--nogfx' in sys.argv)):
        import matplotlib
        matplotlib.use("Agg")

    ## import acolite source
    import acolite as ac

    ## run command line if --cli provided, otherwise use gui
    if '--cli' in sys.argv:
        ac.acolite.acolite_cli(sys.argv)
    else:
        ret = ac.acolite.acolite_gui(sys.argv, version=ac.acolite.config['version'])

if __name__ == '__main__':
    launch_acolite()
